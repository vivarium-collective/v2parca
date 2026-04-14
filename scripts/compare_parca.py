#!/usr/bin/env python
"""
compare_parca.py — per-step side-by-side report comparing v2parca (port-
first composite) against the original ``fitSimData_1`` in vivarium-ecoli.

Walks every available step checkpoint in both directories and emits a
single self-contained HTML report with a sticky left navigation.  Each
Step gets its own section with:

  * runtime (v2parca vs vEcoli, ratio)
  * Input / Output port manifest (documents the step's declared data flow)
  * State comparison — scalars, distributions, and cell_specs entries
    that the step produced, each with max |Δ|, max rel Δ, and KS p-value
  * overlaid histograms for array-valued outputs

Inputs
------
  --v2parca-outdir       DIR  contains checkpoint_step_N.pkl + runtimes.json
                             (produced by scripts/parca_bigraph.py)
  --original-intermediates DIR  contains sim_data_<step>.cPickle +
                             cell_specs_<step>.cPickle (produced by
                             vivarium-ecoli's runscripts/parca.py
                             --save-intermediates)
  -o                    PATH output HTML file

Missing checkpoints on either side render as "not compared" with the
reason listed; partial pipelines don't break the report.
"""

from __future__ import annotations

import argparse
import base64
import importlib
import io
import json
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from scipy import stats as _scipy_stats
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


# ---------------------------------------------------------------------------
# Per-step metadata — name, v2parca checkpoint filename, vEcoli filenames,
# and the port module we pull INPUT_PORTS / OUTPUT_PORTS from.
# ---------------------------------------------------------------------------

STEPS: List[Dict[str, str]] = [
    dict(n=1, name='initialize',        long='Initialize sim_data + scatter',
         v2parca='checkpoint_step_1.pkl',  vecoli_stub='initialize',
         module='v2parca.steps.step_01_initialize'),
    dict(n=2, name='input_adjustments', long='Pre-fitted adjustments to expression + deg rates',
         v2parca='checkpoint_step_2.pkl',  vecoli_stub='input_adjustments',
         module='v2parca.steps.step_02_input_adjustments'),
    dict(n=3, name='basal_specs',       long='Build basal cell specifications',
         v2parca='checkpoint_step_3.pkl',  vecoli_stub='basal_specs',
         module='v2parca.steps.step_03_basal_specs'),
    dict(n=4, name='tf_condition_specs',long='Per-TF + combined condition cell specs',
         v2parca='checkpoint_step_4.pkl',  vecoli_stub='tf_condition_specs',
         module='v2parca.steps.step_04_tf_condition_specs'),
    dict(n=5, name='fit_condition',     long='Bulk distributions + translation supply rates',
         v2parca='checkpoint_step_5.pkl',  vecoli_stub='fit_condition',
         module='v2parca.steps.step_05_fit_condition'),
    dict(n=6, name='promoter_binding',  long='TF-promoter binding probabilities (CVXPY)',
         v2parca='checkpoint_step_6.pkl',  vecoli_stub='promoter_binding',
         module='v2parca.steps.step_06_promoter_binding'),
    dict(n=7, name='adjust_promoters',  long='Ligand concentrations + RNAP recruitment',
         v2parca='checkpoint_step_7.pkl',  vecoli_stub='adjust_promoters',
         module='v2parca.steps.step_07_adjust_promoters'),
    dict(n=8, name='set_conditions',    long='Per-nutrient dicts + mass rescaling',
         v2parca='checkpoint_step_8.pkl',  vecoli_stub='set_conditions',
         module='v2parca.steps.step_08_set_conditions'),
    dict(n=9, name='final_adjustments', long='ppGpp kinetics + amino-acid supply constants',
         v2parca='checkpoint_step_9.pkl',  vecoli_stub='final_adjustments',
         module='v2parca.steps.step_09_final_adjustments'),
]


# Fields to compare as scalars (sim_data attr paths).
SCALARS: List[Tuple[str, Tuple[str, ...]]] = [
    ('mass.avg_cell_dry_mass_init',       ('mass', 'avg_cell_dry_mass_init')),
    ('mass.avg_cell_dry_mass',            ('mass', 'avg_cell_dry_mass')),
    ('mass.avg_cell_water_mass_init',     ('mass', 'avg_cell_water_mass_init')),
    ('mass.fitAvgSolubleTargetMolMass',   ('mass', 'fitAvgSolubleTargetMolMass')),
    ('constants.darkATP',                 ('constants', 'darkATP')),
]

# Array-valued distributions.
DISTRIBUTIONS: List[Tuple[str, Tuple[str, ...]]] = [
    ('RNA expression — basal',    ('process', 'transcription', 'rna_expression', 'basal')),
    ('RNA synthesis prob — basal',('process', 'transcription', 'rna_synth_prob', 'basal')),
    ('RNA deg rates',              ('process', 'transcription', 'rna_data', 'deg_rate')),
    ('Cistron deg rates',          ('process', 'transcription', 'cistron_data', 'deg_rate')),
    ('Protein deg rates',          ('process', 'translation', 'monomer_data', 'deg_rate')),
    ('Translation efficiencies',   ('process', 'translation', 'translation_efficiencies_by_monomer')),
    ('Km endoRNase (transcribed)', ('process', 'transcription', 'rna_data', 'Km_endoRNase')),
    ('Km endoRNase (mature)',      ('process', 'transcription', 'mature_rna_data', 'Km_endoRNase')),
]

CELL_SPECS_FIELDS = [
    'expression', 'synthProb', 'fit_cistron_expression',
    'doubling_time', 'avgCellDryMassInit', 'fitAvgSolubleTargetMolMass',
    'bulkContainer',
]


# ---------------------------------------------------------------------------
# vEcoli pickle compatibility
# ---------------------------------------------------------------------------

def _alias_vivarium_ecoli_modules() -> None:
    """Register vendored ``v2parca.*`` modules under two legacy aliases:
    (1) the top-level vEcoli names (``reconstruction.ecoli.*``, etc.)
        so vivarium-ecoli pickles unpickle.
    (2) the pre-rename ``vparca.*`` names so old v2parca checkpoints
        pickled before the package was renamed still unpickle."""
    for modpath in (
        'v2parca',
        'v2parca.reconstruction.ecoli.simulation_data',
        'v2parca.reconstruction.ecoli.dataclasses',
        'v2parca.wholecell.utils.units',
        'v2parca.ecoli.library.schema',
    ):
        try:
            importlib.import_module(modpath)
        except Exception:
            pass
    for name, mod in list(sys.modules.items()):
        for top in ('v2parca.reconstruction', 'v2parca.wholecell', 'v2parca.ecoli'):
            if name == top or name.startswith(top + '.'):
                alias = name[len('v2parca.'):]
                sys.modules.setdefault(alias, mod)
        if name == 'v2parca' or name.startswith('v2parca.'):
            legacy = 'vparca' + name[len('v2parca'):]
            sys.modules.setdefault(legacy, mod)


def _load_pickle(path: Optional[str]) -> Any:
    if path is None or not os.path.exists(path):
        return None
    _alias_vivarium_ecoli_modules()
    with open(path, 'rb') as f:
        return pickle.load(f)


# ---------------------------------------------------------------------------
# Data-shape navigation (flat v2parca state vs nested SimulationDataEcoli)
# ---------------------------------------------------------------------------

def _get(obj, attr):
    if obj is None:
        return None
    if isinstance(obj, dict):
        return obj.get(attr)
    val = getattr(obj, attr, None)
    if val is not None:
        return val
    try:
        return obj[attr]
    except (KeyError, IndexError, TypeError):
        return None


def _reach(obj, path: Tuple[str, ...]):
    if obj is None:
        return None
    if isinstance(obj, dict) and 'transcription' in obj and 'process' not in obj:
        if path and path[0] == 'process' and len(path) > 1:
            path = path[1:]
        if path and path[0] == 'internal_state' and len(path) > 1:
            if path[1] == 'bulk_molecules':
                path = ('bulk_molecules',) + tuple(path[2:])
    for p in path:
        if obj is None:
            return None
        obj = _get(obj, p)
    return obj


# ---------------------------------------------------------------------------
# Numerical helpers
# ---------------------------------------------------------------------------

def _as_array(x) -> Optional[np.ndarray]:
    if x is None:
        return None
    try:
        if hasattr(x, 'asNumber'):
            x = x.asNumber()
    except Exception:
        pass
    try:
        arr = np.asarray(x)
    except Exception:
        return None
    if arr.dtype.kind not in ('i', 'f', 'u'):
        if arr.dtype.names:
            for name in ('count', 'counts', 'deg_rate'):
                if name in arr.dtype.names:
                    return np.asarray(arr[name], dtype=float)
            for name in arr.dtype.names:
                sub = arr[name]
                if sub.dtype.kind in ('i', 'f', 'u'):
                    return np.asarray(sub, dtype=float)
        return None
    return arr.astype(float, copy=False)


def _safe_rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return float('nan')
    eps = 1e-30
    denom = np.maximum(np.abs(a) + np.abs(b) + eps, eps)
    return float(np.max(np.abs(a - b) / denom))


def _safe_max_abs(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return float('nan')
    return float(np.max(np.abs(a - b)))


def _b64fig(fig, dpi: int = 110) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('ascii')


def _hist_overlay(a, b, title, log=False) -> str:
    aa = _as_array(a); bb = _as_array(b)
    fig, ax = plt.subplots(figsize=(6, 2.8))
    if aa is not None and aa.size:
        v = aa[np.isfinite(aa)]
        if log and (v > 0).any():
            v = np.log10(v[v > 0])
        ax.hist(v, bins=60, alpha=0.55, label='v2parca', color='#2563eb', density=True)
    if bb is not None and bb.size:
        v = bb[np.isfinite(bb)]
        if log and (v > 0).any():
            v = np.log10(v[v > 0])
        ax.hist(v, bins=60, alpha=0.55, label='vEcoli',  color='#dc2626', density=True)
    ax.set_title(title + (' (log10)' if log else ''))
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)
    return _b64fig(fig)


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

def _row_class(rel_diff: Optional[float],
               tol_pass: float = 1e-6, tol_warn: float = 1e-3) -> str:
    if rel_diff is None or not np.isfinite(rel_diff):
        return 'warn'
    if rel_diff < tol_pass:
        return 'pass'
    if rel_diff < tol_warn:
        return 'warn'
    return 'fail'


def _fmt(x):
    if x is None:
        return ''
    if isinstance(x, float):
        if not np.isfinite(x):
            return 'n/a'
        if abs(x) >= 1e4 or (abs(x) > 0 and abs(x) < 1e-3):
            return f'{x:.3e}'
        return f'{x:.4g}'
    return str(x)


# ---------------------------------------------------------------------------
# Per-step section builders
# ---------------------------------------------------------------------------

_STORE_PATH_CACHE: Optional[Dict[str, List[str]]] = None

def _store_paths() -> Dict[str, List[str]]:
    global _STORE_PATH_CACHE
    if _STORE_PATH_CACHE is None:
        try:
            from v2parca.composite import STORE_PATH as _SP
            _STORE_PATH_CACHE = dict(_SP)
        except Exception:
            _STORE_PATH_CACHE = {}
    return _STORE_PATH_CACHE


def _port_rows(ports: Dict[str, Any], store_paths: Dict[str, List[str]]) -> str:
    """Render a port dict as a table: port name, store path, description."""
    visible = [(k, v) for k, v in ports.items() if not k.startswith('tick_')]
    tick_ports = [k for k in ports if k.startswith('tick_')]
    if not visible:
        return ('<p class="meta">(none — this step has no upstream '
                'subsystem reads; it only consumes ordering tokens: '
                f'<code>{", ".join(tick_ports) or "—"}</code>)</p>')
    rows = ['<table class="ports"><tr>'
            '<th>port</th><th>store path</th><th>role</th></tr>']
    for k, v in visible:
        path = store_paths.get(k)
        path_str = '/' + '/'.join(path) if path else '<span class="meta">(not in STORE_PATH)</span>'
        rows.append(
            f'<tr><td><code>{k}</code></td>'
            f'<td><code>{path_str}</code></td>'
            f'<td class="meta">{v}</td></tr>'
        )
    rows.append('</table>')
    if tick_ports:
        rows.append(f'<p class="meta">Plus ordering tokens: '
                    f'<code>{", ".join(tick_ports)}</code> — these enforce '
                    'Step execution order in the bigraph but carry no data.</p>')
    return '\n'.join(rows)


def _port_table(module_name: str) -> str:
    """Render a step's declared data flow: prose header + per-direction
    port tables showing port name, store path, and role."""
    try:
        m = importlib.import_module(module_name)
    except Exception as e:
        return f'<p class="meta">(port manifest unavailable: {e})</p>'
    ins  = getattr(m, 'INPUT_PORTS',  {}) or {}
    outs = getattr(m, 'OUTPUT_PORTS', {}) or {}
    doc = (m.__doc__ or '').strip()
    # Keep the first paragraph of the module docstring (if any) as context.
    doc_para = doc.split('\n\n', 1)[0].replace('\n', ' ').strip() if doc else ''

    intro = [
        '<p class="meta">Each Step declares <code>INPUT_PORTS</code> '
        '(stores it reads) and <code>OUTPUT_PORTS</code> (stores it writes). '
        'Port names are resolved to absolute store paths via '
        '<code>STORE_PATH</code> in <code>v2parca/composite.py</code>, and '
        'the composite wires each port to the corresponding location in the '
        'nested bigraph. The <em>role</em> column shows the description '
        'given in the step module (e.g. <code>sim_data.transcription</code> '
        'means "mirrors the <code>transcription</code> subsystem on '
        '<code>SimulationDataEcoli</code>"; <code>overwrite</code> means '
        'the port replaces whatever value lives at that store).</p>'
    ]
    if doc_para:
        intro.append(f'<p><strong>Step purpose:</strong> {doc_para}</p>')

    sp = _store_paths()
    parts = list(intro)
    parts.append(f'<h4>Inputs — reads ({sum(1 for k in ins if not k.startswith("tick_"))})</h4>')
    parts.append(_port_rows(ins, sp))
    parts.append(f'<h4>Outputs — writes ({sum(1 for k in outs if not k.startswith("tick_"))})</h4>')
    parts.append(_port_rows(outs, sp))
    return '\n'.join(parts)


def _section_runtime(vt: Optional[float], ot: Optional[float]) -> str:
    ratio = (vt / ot) if (vt and ot and ot > 0) else None
    return (
        '<h3>Runtime</h3>'
        '<table><tr><th>v2parca</th><th>vEcoli</th><th>ratio (v2parca/vEcoli)</th></tr>'
        f'<tr><td>{_fmt(vt)} s</td><td>{_fmt(ot)} s</td><td>{_fmt(ratio)}</td></tr>'
        '</table>'
    )


def _section_scalars(v2parca, original) -> str:
    out = ['<h3>Scalar state</h3>',
           '<table><tr><th>path</th><th>v2parca</th><th>vEcoli</th>'
           '<th>rel Δ</th></tr>']
    any_row = False
    for label, path in SCALARS:
        a = _reach(v2parca, path); b = _reach(original, path)
        a_num = float(a.asNumber()) if a is not None and hasattr(a, 'asNumber') else (
            float(a) if isinstance(a, (int, float, np.number)) else None)
        b_num = float(b.asNumber()) if b is not None and hasattr(b, 'asNumber') else (
            float(b) if isinstance(b, (int, float, np.number)) else None)
        if a_num is None and b_num is None:
            continue
        any_row = True
        rd = (_safe_rel_diff(np.array([a_num]), np.array([b_num]))
              if (a_num is not None and b_num is not None) else float('nan'))
        out.append(
            f'<tr class="{_row_class(rd)}">'
            f'<td><code>{label}</code></td>'
            f'<td>{_fmt(a_num)}</td><td>{_fmt(b_num)}</td>'
            f'<td>{_fmt(rd)}</td></tr>')
    out.append('</table>')
    return '\n'.join(out) if any_row else ''


def _section_distributions(v2parca, original) -> str:
    figs_html = []
    tbl = ['<h3>Distribution numerical summary</h3>',
           '<table><tr><th>distribution</th><th>shape</th>'
           '<th>max |Δ|</th><th>max rel Δ</th><th>KS p-value</th></tr>']
    any_dist = False
    for label, path in DISTRIBUTIONS:
        va = _reach(v2parca, path); oa = _reach(original, path)
        a = _as_array(va); b = _as_array(oa)
        if a is None and b is None:
            continue
        any_dist = True
        shape = (a.shape if a is not None else b.shape)
        if a is not None and b is not None and a.shape == b.shape:
            rd = _safe_rel_diff(a, b); ma = _safe_max_abs(a, b)
            ks = (_scipy_stats.ks_2samp(a.ravel(), b.ravel()).pvalue
                  if HAVE_SCIPY else None)
        else:
            rd, ma, ks = float('nan'), float('nan'), None
        log = (a is not None and a.size > 0 and (np.asarray(a) > 0).any()
               and (float(a.max()) / max(float(a[a > 0].min()) if (a > 0).any() else 1, 1e-30) > 100))
        img = _hist_overlay(va, oa, label, log=log)
        figs_html.append(f'<div><strong>{label}</strong><br>'
                         f'<img src="data:image/png;base64,{img}"/></div>')
        tbl.append(
            f'<tr class="{_row_class(rd)}">'
            f'<td><code>{label}</code></td>'
            f'<td>{shape}</td>'
            f'<td>{_fmt(ma)}</td><td>{_fmt(rd)}</td><td>{_fmt(ks)}</td></tr>')
    tbl.append('</table>')
    if not any_dist:
        return ''
    return ('<h3>Distributions</h3>'
            '<div class="grid">' + '\n'.join(figs_html) + '</div>' + '\n' + '\n'.join(tbl))


def _section_cell_specs(v2parca, original_cell_specs) -> str:
    cs_v = _reach(v2parca, ('cell_specs',)) or {}
    cs_o = original_cell_specs or {}
    if not cs_v and not cs_o:
        return ''
    common = sorted(set(cs_v.keys()) & set(cs_o.keys()))
    only_v = sorted(set(cs_v.keys()) - set(cs_o.keys()))
    only_o = sorted(set(cs_o.keys()) - set(cs_v.keys()))
    out = ['<h3>cell_specs (per-condition max rel Δ)</h3>']
    if only_v or only_o:
        out.append(
            f'<p class="meta">Only v2parca: <code>{only_v}</code>. '
            f'Only vEcoli: <code>{only_o}</code>.</p>')
    if not common:
        return '\n'.join(out) + '<p class="meta">(no shared conditions)</p>'
    out.append('<table><tr><th>condition</th>'
               + ''.join(f'<th>{f}</th>' for f in CELL_SPECS_FIELDS)
               + '</tr>')
    for cond in common:
        v_spec = cs_v[cond]; o_spec = cs_o[cond]
        cells = []
        for field in CELL_SPECS_FIELDS:
            va = v_spec.get(field) if isinstance(v_spec, dict) else None
            oa = o_spec.get(field) if isinstance(o_spec, dict) else None
            if va is None or oa is None:
                cells.append('<td class="meta">—</td>'); continue
            aa = _as_array(va); ob = _as_array(oa)
            if aa is None or ob is None:
                try:
                    av = float(va.asNumber()) if hasattr(va, 'asNumber') else float(va)
                    bv = float(oa.asNumber()) if hasattr(oa, 'asNumber') else float(oa)
                    rd = _safe_rel_diff(np.array([av]), np.array([bv]))
                    cells.append(f'<td class="{_row_class(rd)}">{_fmt(rd)}</td>')
                except Exception:
                    cells.append('<td class="meta">?</td>')
                continue
            if aa.shape != ob.shape:
                cells.append(f'<td class="fail">shape {aa.shape}≠{ob.shape}</td>')
                continue
            rd = _safe_rel_diff(aa, ob)
            cells.append(f'<td class="{_row_class(rd)}">{_fmt(rd)}</td>')
        out.append(f'<tr><td><code>{cond}</code></td>{"".join(cells)}</tr>')
    out.append('</table>')
    return '\n'.join(out)


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = r"""<!doctype html>
<html><head><meta charset='utf-8'>
<title>v2parca vs ParCa — per-step comparison</title>
<style>
  :root {{ --accent: #2563eb; --pass: #ecfdf5; --warn: #fffbeb; --fail: #fef2f2; }}
  * {{ box-sizing: border-box; }}
  body  {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
           margin: 0; color: #111; background: #fafafa; }}
  #layout {{ display: grid; grid-template-columns: 220px 1fr; min-height: 100vh; }}
  nav    {{ position: sticky; top: 0; height: 100vh; overflow-y: auto;
            background: #111; color: #eee; padding: 18px 12px; }}
  nav h2  {{ color: #fff; font-size: 14px; margin: 0 0 12px 4px;
             letter-spacing: .5px; text-transform: uppercase; }}
  nav a   {{ display: block; color: #ccc; text-decoration: none;
             padding: 6px 10px; border-radius: 4px; font-size: 13px;
             margin-bottom: 2px; }}
  nav a:hover {{ background: #1f2937; color: #fff; }}
  nav a.active {{ background: var(--accent); color: #fff; }}
  nav .status-pass {{ color: #34d399; }}
  nav .status-warn {{ color: #fbbf24; }}
  nav .status-fail {{ color: #f87171; }}
  nav .status-na   {{ color: #6b7280; }}
  main   {{ padding: 24px 40px; max-width: 1100px; }}
  h1     {{ border-bottom: 3px solid var(--accent); padding-bottom: 6px; }}
  section {{ margin-bottom: 48px; padding: 24px 28px;
             background: #fff; border-radius: 8px;
             box-shadow: 0 1px 3px rgba(0,0,0,.06); }}
  section h2 {{ margin-top: 0; border-bottom: 1px solid #eee; padding-bottom: 6px; }}
  section h3 {{ margin-top: 24px; color: #374151; }}
  table  {{ border-collapse: collapse; width: 100%; font-size: 13px;
            margin-top: 8px; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 5px 8px; text-align: left; }}
  th     {{ background: #f3f4f6; }}
  tr.pass {{ background: var(--pass); }}
  tr.warn {{ background: var(--warn); }}
  tr.fail {{ background: var(--fail); }}
  td.pass {{ background: var(--pass); }}
  td.warn {{ background: var(--warn); }}
  td.fail {{ background: var(--fail); }}
  img    {{ max-width: 100%; border: 1px solid #eee; margin: 6px 0;
            background: #fff; }}
  .grid  {{ display: grid; grid-template-columns: 1fr 1fr; gap: 10px 16px;
            margin-top: 10px; }}
  code   {{ background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }}
  .meta  {{ color: #6b7280; font-size: 12px; }}
  .banner-na {{ background: #f9fafb; border-left: 4px solid #9ca3af;
                padding: 10px 14px; color: #6b7280; }}
</style>
<script>
  // Highlight the current section in the left nav while scrolling.
  document.addEventListener('DOMContentLoaded', () => {{
    const sections = document.querySelectorAll('section[id]');
    const links    = document.querySelectorAll('nav a[href^="#"]');
    const byId = {{}};
    links.forEach(a => byId[a.getAttribute('href').slice(1)] = a);
    const io = new IntersectionObserver(entries => {{
      entries.forEach(e => {{
        if (e.isIntersecting) {{
          links.forEach(a => a.classList.remove('active'));
          const a = byId[e.target.id]; if (a) a.classList.add('active');
        }}
      }});
    }}, {{ rootMargin: '-45% 0px -55% 0px' }});
    sections.forEach(s => io.observe(s));
  }});
</script>
</head><body>
<div id='layout'>
<nav>
  <h2>v2parca report</h2>
  <a href="#overview">Overview</a>
  {nav_links}
</nav>
<main>
<h1>v2parca vs vEcoli ParCa</h1>
<p class='meta'>{meta}</p>
{sections}
</main>
</div></body></html>"""


def _runtime_bar(v2parca_times: Dict[str, float],
                 original_times: Dict[str, float]) -> str:
    steps = [f'step_{n}' for n in range(1, 10)]
    v = [v2parca_times.get(s, 0.0)  for s in steps]
    o = [original_times.get(s, 0.0) for s in steps]
    fig, ax = plt.subplots(figsize=(8, 3.2))
    x = np.arange(len(steps)); w = 0.38
    ax.bar(x - w/2, v, w, label='v2parca', color='#2563eb')
    ax.bar(x + w/2, o, w, label='vEcoli',  color='#dc2626')
    ax.set_xticks(x); ax.set_xticklabels([f'step {n}' for n in range(1, 10)],
                                          rotation=30)
    ax.set_ylabel('seconds')
    ax.set_title('Per-step runtime')
    ax.set_yscale('symlog', linthresh=1)
    ax.grid(True, axis='y', alpha=0.3); ax.legend(fontsize=9)
    return _b64fig(fig)


def build_report(v2parca_outdir: str, vecoli_dir: Optional[str],
                 output_path: str) -> None:
    # Load runtimes.
    v2parca_rt_path = os.path.join(v2parca_outdir, 'runtimes.json')
    v2parca_rt = json.load(open(v2parca_rt_path)) if os.path.exists(v2parca_rt_path) else {}
    # vEcoli's --save-intermediates doesn't write a runtimes.json; we
    # parse its per-stage "Ran X in Ys" prints if a log is available.
    original_rt = _maybe_parse_vecoli_runtimes(vecoli_dir)

    # Per-step section data — status ('pass' / 'warn' / 'fail' / 'na')
    # is a quick summary for the left-nav colored dots.
    section_items = []
    for step in STEPS:
        n = step['n']
        v2parca_pkl = os.path.join(v2parca_outdir, step['v2parca'])
        vecoli_sd  = (os.path.join(vecoli_dir, f"sim_data_{step['vecoli_stub']}.cPickle")
                      if vecoli_dir else None)
        vecoli_cs  = (os.path.join(vecoli_dir, f"cell_specs_{step['vecoli_stub']}.cPickle")
                      if vecoli_dir else None)

        v2parca_state = _load_pickle(v2parca_pkl)
        original     = _load_pickle(vecoli_sd)
        original_cs  = _load_pickle(vecoli_cs)

        # Attach cell_specs to the original sim_data if present.
        if original is not None and original_cs is not None:
            try:
                original.cell_specs = original_cs
            except Exception:
                pass

        vt = v2parca_rt.get(f'step_{n}')
        ot = original_rt.get(step['vecoli_stub'])

        section_id = f'step_{n}'
        parts = []
        parts.append(f'<section id="{section_id}">')
        parts.append(f'<h2>Step {n} — {step["name"]}</h2>')
        parts.append(f'<p class="meta">{step["long"]}</p>')

        # Availability banner.
        have_vp = v2parca_state is not None
        have_ve = original is not None
        if not have_vp and not have_ve:
            parts.append('<div class="banner-na">No checkpoint on either side. '
                         'Run `scripts/parca_bigraph.py` for v2parca and '
                         '`runscripts/parca.py --save-intermediates` for vEcoli.</div>')
            status = 'na'
        elif not have_ve:
            parts.append('<div class="banner-na">vEcoli reference pickle unavailable '
                         '— showing v2parca side only.</div>')
            status = 'warn'
        elif not have_vp:
            parts.append('<div class="banner-na">v2parca checkpoint unavailable '
                         '— showing vEcoli side only.</div>')
            status = 'warn'
        else:
            status = 'pass'

        parts.append(_section_runtime(vt, ot))
        parts.append('<h3>Declared data flow</h3>')
        parts.append(_port_table(step['module']))

        if have_vp and have_ve:
            parts.append(_section_scalars(v2parca_state, original))
            parts.append(_section_distributions(v2parca_state, original))
            parts.append(_section_cell_specs(v2parca_state, original_cs))

        parts.append('</section>')
        section_items.append({
            'n': n, 'name': step['name'], 'status': status,
            'html': '\n'.join(parts),
        })

    # Overview section.
    overview = ['<section id="overview">', '<h2>Overview</h2>']
    overview.append(f'<p>Generated {time.strftime("%Y-%m-%d %H:%M:%S")}</p>')
    overview.append(f'<p class="meta">v2parca checkpoints: <code>{v2parca_outdir}</code>. '
                    f'vEcoli intermediates: <code>{vecoli_dir}</code>.</p>')
    overview.append('<h3>Per-step runtime</h3>')
    overview.append(f'<img src="data:image/png;base64,{_runtime_bar(v2parca_rt, original_rt)}"/>')
    overview.append('<h3>Step-by-step availability</h3>')
    overview.append('<table><tr><th>step</th><th>v2parca</th><th>vEcoli</th>'
                    '<th>compared</th></tr>')
    for step, item in zip(STEPS, section_items):
        n = step['n']
        vp_ok = os.path.exists(os.path.join(v2parca_outdir, step['v2parca']))
        ve_ok = (vecoli_dir is not None and
                 os.path.exists(os.path.join(vecoli_dir,
                                             f"sim_data_{step['vecoli_stub']}.cPickle")))
        overview.append(
            f'<tr><td>step {n} — {step["name"]}</td>'
            f'<td>{"✓" if vp_ok else "—"}</td>'
            f'<td>{"✓" if ve_ok else "—"}</td>'
            f'<td>{"yes" if vp_ok and ve_ok else "no"}</td></tr>')
    overview.append('</table>')
    overview.append('</section>')

    nav_links = '\n'.join(
        f'<a href="#step_{it["n"]}"><span class="status-{it["status"]}">●</span> '
        f'Step {it["n"]} — {it["name"]}</a>'
        for it in section_items
    )

    meta = (f'Generated {time.strftime("%Y-%m-%d %H:%M:%S")} — '
            f'tolerance: pass &lt; 1e-6, warn &lt; 1e-3')
    sections_html = '\n'.join(overview) + '\n' + \
                    '\n'.join(it['html'] for it in section_items)

    html = _HTML_TEMPLATE.format(
        nav_links=nav_links, meta=meta, sections=sections_html)

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(html)
    size_kb = os.path.getsize(output_path) / 1024
    print(f'wrote {output_path} ({size_kb:.1f} KB)')


def _maybe_parse_vecoli_runtimes(vecoli_dir: Optional[str]) -> Dict[str, float]:
    """Best-effort: look for a 'Ran X in Ys' log next to the pickles."""
    if not vecoli_dir:
        return {}
    for candidate in (os.path.join(vecoli_dir, 'runtimes.json'),
                      os.path.join(os.path.dirname(vecoli_dir), 'runtimes.json')):
        if os.path.exists(candidate):
            try:
                return json.load(open(candidate))
            except Exception:
                pass
    # Parse from vecoli_parca logs in /tmp if present.
    for log in ('/tmp/vecoli_intermediates.log',
                '/tmp/vecoli_solo.log', '/tmp/vecoli_parca_run3.log'):
        if os.path.exists(log):
            import re
            out = {}
            for m in re.finditer(r'Ran (\S+) in (\d+) s', open(log).read()):
                out[m.group(1)] = float(m.group(2))
            if out:
                return out
    return {}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--v2parca-outdir', default='out/sim_data',
                   help='dir with checkpoint_step_N.pkl + runtimes.json')
    p.add_argument('--original-intermediates', default='out/original_intermediates',
                   help='dir with sim_data_<step>.cPickle + cell_specs_<step>.cPickle')
    p.add_argument('-o', '--output', default='out/compare/report.html')
    args = p.parse_args()

    vecoli_dir = args.original_intermediates if os.path.isdir(
        args.original_intermediates) else None
    build_report(args.v2parca_outdir, vecoli_dir, args.output)


if __name__ == '__main__':
    main()
