#!/usr/bin/env python
"""
compare_parca.py — side-by-side report comparing the original
``fitSimData_1`` (vivarium-ecoli) with vParCa's port-first
``build_parca_composite``.

The ParCa is a one-shot fitting operation, not a dynamic simulation, so
the useful comparison axes are:

  1.  **Runtimes** per pipeline step and overall
  2.  **Fitted distributions** that the ParCa emits — RNA expression,
      synthesis probability, RNA / cistron / protein degradation rates,
      translation efficiencies, bulk counts
  3.  **Per-condition initial conditions** in ``cell_specs`` — expression,
      synthProb, doubling_time, avgCellDryMassInit, bulkContainer, …
  4.  **Scalar constants** that the ParCa computes (e.g. ``darkATP``,
      ``avg_cell_dry_mass_init``)

Output is a single self-contained HTML report with base64-embedded
matplotlib figures, modeled after v2ecoli's ``compare_report.py``.

Usage
-----

    # Compare two pre-computed outputs (typical flow):
    python scripts/compare_parca.py \\
        --vparca-state      out/sim_data/parca_state.pkl \\
        --original-sim-data out/orig/sim_data.cPickle \\
        --original-runtimes out/orig/runtimes.json \\
        -o out/compare/report.html

    # Run both engines from scratch and compare (slow — ~60+ min):
    python scripts/compare_parca.py --run --mode fast -o out/compare/report.html

The ``--original-sim-data`` / ``--original-runtimes`` inputs are produced
by running ``vivarium-ecoli``'s ``runscripts/parca.py`` once; point the
script at the resulting ``sim_data.cPickle`` and a JSON you write with
``{"step_1": secs, "step_2": secs, …}`` if the original's per-step timing
is available (otherwise omit and the report shows only total).
"""

from __future__ import annotations

import argparse
import base64
import io
import json
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

# Make vparca/ imports work when run from anywhere.
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
# Utilities
# ---------------------------------------------------------------------------

def _b64fig(fig, dpi: int = 120) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('ascii')


def _as_array(x) -> np.ndarray:
    """Coerce Unum / structured-array / list-like to a float numpy array."""
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
        # Try to pull out a 'count' or first numeric column from structured.
        if arr.dtype.names:
            for name in ('count', 'counts', 'deg_rate'):
                if name in arr.dtype.names:
                    return np.asarray(arr[name], dtype=float)
            # fall back: first numeric field
            for name in arr.dtype.names:
                sub = arr[name]
                if sub.dtype.kind in ('i', 'f', 'u'):
                    return np.asarray(sub, dtype=float)
        return None
    return arr.astype(float, copy=False)


def _safe_rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    """Max |a - b| / (|a| + |b| + eps) — symmetric, avoids div-by-zero."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return np.nan
    eps = 1e-30
    denom = np.maximum(np.abs(a) + np.abs(b) + eps, eps)
    return float(np.max(np.abs(a - b) / denom))


def _safe_max_abs(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return np.nan
    return float(np.max(np.abs(a - b)))


# ---------------------------------------------------------------------------
# Data access — walk both sim_data shapes into a uniform (path, array) map
# ---------------------------------------------------------------------------

def _get(obj, attr):
    """Attribute access that tolerates dicts, namespaces, and None."""
    if obj is None:
        return None
    if isinstance(obj, dict):
        return obj.get(attr)
    return getattr(obj, attr, None)


def _walk(obj, path=None):
    """Yield (dot-path, attr-object) tuples for a nested sim_data-like thing."""
    if path is None:
        path = []
    # Skip Python built-ins on containers
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _walk(v, path + [str(k)])
    elif hasattr(obj, '__dict__'):
        for k, v in obj.__dict__.items():
            if k.startswith('_'):
                continue
            yield from _walk(v, path + [k])
    else:
        yield ('.'.join(path), obj)


# Distributions to extract.  Each entry: (label, (accessor, …))
# accessors are applied in sequence to reach the leaf.
DISTRIBUTIONS: List[Tuple[str, Tuple[str, ...]]] = [
    ('RNA expression — basal',             ('process', 'transcription', 'rna_expression', 'basal')),
    ('RNA synthesis prob — basal',         ('process', 'transcription', 'rna_synth_prob', 'basal')),
    ('RNA deg rates',                       ('process', 'transcription', 'rna_data', 'deg_rate')),
    ('Cistron deg rates',                   ('process', 'transcription', 'cistron_data', 'deg_rate')),
    ('Protein deg rates',                   ('process', 'translation', 'monomer_data', 'deg_rate')),
    ('Translation efficiencies',            ('process', 'translation', 'translation_efficiencies_by_monomer')),
    ('Km endoRNase (transcribed)',          ('process', 'transcription', 'rna_data', 'Km_endoRNase')),
    ('Km endoRNase (mature)',               ('process', 'transcription', 'mature_rna_data', 'Km_endoRNase')),
    ('cistron_data deg_rate (via struct)',  ('process', 'transcription', 'cistron_data')),
]


SCALARS: List[Tuple[str, Tuple[str, ...]]] = [
    ('mass.avg_cell_dry_mass_init',       ('mass', 'avg_cell_dry_mass_init')),
    ('mass.avg_cell_dry_mass',            ('mass', 'avg_cell_dry_mass')),
    ('mass.avg_cell_water_mass_init',     ('mass', 'avg_cell_water_mass_init')),
    ('mass.fitAvgSolubleTargetMolMass',   ('mass', 'fitAvgSolubleTargetMolMass')),
    ('constants.darkATP',                 ('constants', 'darkATP')),
]


# cell_specs fields to compare per condition.
CELL_SPECS_FIELDS = [
    'expression', 'synthProb', 'fit_cistron_expression',
    'doubling_time', 'avgCellDryMassInit', 'fitAvgSolubleTargetMolMass',
    'bulkContainer',
]


def _reach(obj, path: Tuple[str, ...]):
    """Apply path sequentially; tolerates dict or attr access."""
    for p in path:
        if obj is None:
            return None
        obj = _get(obj, p)
    return obj


def _sim_data_like_from_vparca_state(state: Dict[str, Any]) -> Any:
    """Wrap a vParCa composite.state into a sim_data-shaped namespace so
    ``_reach(path)`` works uniformly for both engines.

    The vParCa state already stores each sim_data subsystem as a live
    object at its canonical path (``state['process']['transcription']``,
    ``state['mass']``, …); we just hand it back as-is since it's
    dict-shaped matching sim_data's nested-attribute layout.
    """
    return state


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _hist_overlay(arr_a, arr_b, label_a, label_b, title, log=False) -> str:
    a = _as_array(arr_a); b = _as_array(arr_b)
    if a is None and b is None:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, 'data unavailable', ha='center', va='center')
        ax.set_axis_off(); return _b64fig(fig)

    fig, ax = plt.subplots(figsize=(6, 3))
    if a is not None and a.size:
        av = a[np.isfinite(a)]
        if log and (av > 0).any():
            av = np.log10(av[av > 0])
        ax.hist(av, bins=60, alpha=0.55, label=label_a, color='#dc2626', density=True)
    if b is not None and b.size:
        bv = b[np.isfinite(b)]
        if log and (bv > 0).any():
            bv = np.log10(bv[bv > 0])
        ax.hist(bv, bins=60, alpha=0.55, label=label_b, color='#2563eb', density=True)
    ax.set_title(title + (' (log10)' if log else ''))
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)
    return _b64fig(fig)


def _runtime_bar(vparca_times: Dict[str, float],
                 original_total: Optional[float]) -> str:
    steps = [f'step_{n}' for n in range(1, 10)]
    vparca_vals = [vparca_times.get(s, 0.0) for s in steps]
    fig, ax = plt.subplots(figsize=(8, 3.5))
    x = np.arange(len(steps))
    ax.bar(x, vparca_vals, width=0.6, color='#2563eb', label='vParCa (per step)')
    if original_total is not None:
        ax.axhline(original_total, color='#dc2626', lw=1.5, ls='--',
                   label=f'original total = {original_total:.0f}s')
    ax.set_xticks(x); ax.set_xticklabels(steps, rotation=30)
    ax.set_ylabel('seconds'); ax.set_title('Pipeline runtime')
    ax.grid(True, axis='y', alpha=0.3); ax.legend(fontsize=8)
    return _b64fig(fig)


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """<!doctype html>
<html><head><meta charset='utf-8'>
<title>vParCa vs ParCa — fitting comparison</title>
<style>
  body  {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
           max-width: 1100px; margin: 20px auto; padding: 0 24px; color: #111; }}
  h1    {{ border-bottom: 3px solid #2563eb; padding-bottom: 6px; }}
  h2    {{ margin-top: 36px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }}
  h3    {{ margin-top: 20px; color: #333; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  th, td {{ border: 1px solid #ccc; padding: 5px 8px; text-align: left; }}
  th    {{ background: #f3f4f6; }}
  tr.pass {{ background: #ecfdf5; }}
  tr.warn {{ background: #fffbeb; }}
  tr.fail {{ background: #fef2f2; }}
  img   {{ max-width: 100%; border: 1px solid #eee; margin: 6px 0; }}
  .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 10px 16px; }}
  code  {{ background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }}
  .meta {{ color: #555; font-size: 13px; }}
</style></head><body>
<h1>vParCa vs ParCa — fitting comparison</h1>
<p class='meta'>{meta}</p>
{body}
</body></html>"""


def _row_class(rel_diff: float, tol_pass: float = 1e-6, tol_warn: float = 1e-3) -> str:
    if rel_diff is None or not np.isfinite(rel_diff):
        return 'warn'
    if rel_diff < tol_pass:
        return 'pass'
    if rel_diff < tol_warn:
        return 'warn'
    return 'fail'


def _fmt(x):
    if x is None: return ''
    if isinstance(x, float):
        if not np.isfinite(x): return 'n/a'
        if abs(x) >= 1e4 or (abs(x) > 0 and abs(x) < 1e-3):
            return f'{x:.3e}'
        return f'{x:.4g}'
    return str(x)


# ---------------------------------------------------------------------------
# Report builder
# ---------------------------------------------------------------------------

def build_report(vparca, original, vparca_times, original_times,
                 output_path: str, tol_pass: float = 1e-6,
                 tol_warn: float = 1e-3) -> None:
    """Write the HTML report.  vparca/original should be dict-like or
    namespace-like containers that respond to ``_reach(path)`` uniformly."""
    body_parts: List[str] = []

    # -- 1. Runtimes ------------------------------------------------------
    body_parts.append('<h2>1. Runtimes</h2>')
    orig_total = sum(original_times.values()) if original_times else None
    body_parts.append(f"<img src='data:image/png;base64,{_runtime_bar(vparca_times, orig_total)}'/>")
    body_parts.append('<table><tr><th>step</th><th>vParCa (s)</th>'
                      '<th>original (s)</th><th>ratio</th></tr>')
    for n in range(1, 10):
        k = f'step_{n}'
        vt = vparca_times.get(k); ot = original_times.get(k) if original_times else None
        ratio = (vt / ot) if (vt and ot) else None
        body_parts.append(
            f'<tr><td>{k}</td><td>{_fmt(vt)}</td>'
            f'<td>{_fmt(ot)}</td><td>{_fmt(ratio)}</td></tr>')
    v_tot = sum(vparca_times.values()) if vparca_times else None
    body_parts.append(
        f'<tr><th>TOTAL</th><th>{_fmt(v_tot)}</th>'
        f'<th>{_fmt(orig_total)}</th>'
        f'<th>{_fmt((v_tot / orig_total) if (v_tot and orig_total) else None)}</th></tr>')
    body_parts.append('</table>')

    # -- 2. Scalar constants ---------------------------------------------
    body_parts.append('<h2>2. Scalar constants</h2>')
    body_parts.append('<table><tr><th>path</th><th>vParCa</th>'
                      '<th>original</th><th>rel diff</th></tr>')
    for label, path in SCALARS:
        a = _reach(vparca, path); b = _reach(original, path)
        a_num = float(a.asNumber()) if a is not None and hasattr(a, 'asNumber') else (float(a) if a is not None else None)
        b_num = float(b.asNumber()) if b is not None and hasattr(b, 'asNumber') else (float(b) if b is not None else None)
        if a_num is not None and b_num is not None:
            rd = _safe_rel_diff(np.array([a_num]), np.array([b_num]))
        else:
            rd = np.nan
        body_parts.append(
            f'<tr class="{_row_class(rd, tol_pass, tol_warn)}">'
            f'<td><code>{label}</code></td>'
            f'<td>{_fmt(a_num)}</td><td>{_fmt(b_num)}</td>'
            f'<td>{_fmt(rd)}</td></tr>')
    body_parts.append('</table>')

    # -- 3. Fitted distributions -----------------------------------------
    body_parts.append('<h2>3. Fitted distributions</h2>')
    body_parts.append('<p class="meta">Histograms are density-normalized; log10 '
                      'where the underlying quantity is strictly positive with '
                      'a large dynamic range.</p>')
    body_parts.append('<div class="grid">')
    dist_table_rows = []
    for label, path in DISTRIBUTIONS:
        va = _reach(vparca, path); oa = _reach(original, path)
        a = _as_array(va); b = _as_array(oa)
        shape_a = tuple(a.shape) if a is not None else None
        shape_b = tuple(b.shape) if b is not None else None
        if a is not None and b is not None and shape_a == shape_b:
            rd = _safe_rel_diff(a, b); ma = _safe_max_abs(a, b)
            ks = (_scipy_stats.ks_2samp(a.ravel(), b.ravel()).pvalue
                  if HAVE_SCIPY else None)
        else:
            rd, ma, ks = np.nan, np.nan, None
        log_scale = (a is not None and a.size and np.all(a[np.isfinite(a)] >= 0)
                     and (a.max() / max(a.min(), 1e-30) > 100 if a.size else False))
        img = _hist_overlay(va, oa, 'vParCa', 'original', label, log=log_scale)
        body_parts.append(
            f'<div><strong>{label}</strong><br>'
            f'<img src="data:image/png;base64,{img}"/></div>')
        dist_table_rows.append((label, shape_a, shape_b, rd, ma, ks))
    body_parts.append('</div>')

    body_parts.append('<h3>Per-distribution numerical summary</h3>')
    body_parts.append('<table><tr><th>distribution</th><th>shape (vParCa)</th>'
                      '<th>shape (orig)</th><th>max |Δ|</th>'
                      '<th>max rel Δ</th><th>KS p-value</th></tr>')
    for label, sa, sb, rd, ma, ks in dist_table_rows:
        body_parts.append(
            f'<tr class="{_row_class(rd, tol_pass, tol_warn)}">'
            f'<td><code>{label}</code></td>'
            f'<td>{sa}</td><td>{sb}</td>'
            f'<td>{_fmt(ma)}</td><td>{_fmt(rd)}</td><td>{_fmt(ks)}</td></tr>')
    body_parts.append('</table>')

    # -- 4. cell_specs per-condition -------------------------------------
    body_parts.append('<h2>4. Initial conditions (cell_specs)</h2>')
    cs_v = _reach(vparca, ('cell_specs',)) or {}
    cs_o = _reach(original, ('cell_specs',)) or {}
    common_conditions = sorted(set(cs_v.keys()) & set(cs_o.keys()))
    only_v = sorted(set(cs_v.keys()) - set(cs_o.keys()))
    only_o = sorted(set(cs_o.keys()) - set(cs_v.keys()))
    body_parts.append(f'<p class="meta">Common conditions: <code>{common_conditions}</code>. '
                      f'Only vParCa: <code>{only_v}</code>. '
                      f'Only original: <code>{only_o}</code>.</p>')
    body_parts.append('<table><tr><th>condition</th>'
                      + ''.join(f'<th>{f}</th>' for f in CELL_SPECS_FIELDS)
                      + '</tr>')
    for cond in common_conditions:
        v_spec = cs_v[cond]; o_spec = cs_o[cond]
        cells = []
        for field in CELL_SPECS_FIELDS:
            va = v_spec.get(field) if isinstance(v_spec, dict) else None
            oa = o_spec.get(field) if isinstance(o_spec, dict) else None
            if va is None or oa is None:
                cells.append('<td class="meta">—</td>'); continue
            aa = _as_array(va); ob = _as_array(oa)
            if aa is None or ob is None:
                # Scalar with units?
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
        body_parts.append(f'<tr><td><code>{cond}</code></td>{"".join(cells)}</tr>')
    body_parts.append('</table>')

    # -- 5. Top divergences ----------------------------------------------
    body_parts.append('<h2>5. Top divergences (distributions)</h2>')
    worst = sorted(
        [(label, rd) for (label, _, _, rd, _, _) in dist_table_rows
         if rd is not None and np.isfinite(rd)],
        key=lambda x: -x[1])[:10]
    if worst:
        body_parts.append('<table><tr><th>rank</th><th>distribution</th>'
                          '<th>max rel Δ</th></tr>')
        for i, (label, rd) in enumerate(worst, 1):
            body_parts.append(
                f'<tr class="{_row_class(rd, tol_pass, tol_warn)}">'
                f'<td>{i}</td><td><code>{label}</code></td>'
                f'<td>{_fmt(rd)}</td></tr>')
        body_parts.append('</table>')
    else:
        body_parts.append('<p class="meta">no comparable distributions loaded.</p>')

    meta = (f'generated {time.strftime("%Y-%m-%d %H:%M:%S")} — '
            f'tolerances: pass &lt; {tol_pass:.0e}, warn &lt; {tol_warn:.0e}, '
            f'else fail')
    html = _HTML_TEMPLATE.format(meta=meta, body='\n'.join(body_parts))

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(html)
    print(f'wrote {output_path} ({os.path.getsize(output_path) / 1024:.1f} KB)')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _load_pickle(path: str) -> Any:
    with open(path, 'rb') as f:
        return pickle.load(f)


def _load_json(path: str) -> Dict:
    with open(path) as f:
        return json.load(f)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vparca-state', type=str,
                   help='pickle of composite.state from scripts/parca_bigraph.py')
    p.add_argument('--original-sim-data', type=str,
                   help='pickle of a fitted SimulationDataEcoli from vivarium-ecoli')
    p.add_argument('--vparca-runtimes', type=str,
                   help='JSON {step_N: seconds, …} for vParCa (optional)')
    p.add_argument('--original-runtimes', type=str,
                   help='JSON {step_N: seconds, …} for the original ParCa')
    p.add_argument('--run', action='store_true',
                   help='run both engines from scratch (not yet implemented)')
    p.add_argument('--mode', choices=['fast', 'full'], default='fast',
                   help='(--run only) debug (fast) vs full ParCa')
    p.add_argument('-o', '--output', type=str, default='out/compare/report.html')
    p.add_argument('--tol-pass', type=float, default=1e-6)
    p.add_argument('--tol-warn', type=float, default=1e-3)
    args = p.parse_args()

    if args.run:
        raise SystemExit(
            "--run is a stub; run scripts/parca_bigraph.py and vivarium-ecoli's "
            "runscripts/parca.py separately, then pass --vparca-state and "
            "--original-sim-data here.")

    if not args.vparca_state or not args.original_sim_data:
        raise SystemExit(
            "must provide --vparca-state and --original-sim-data "
            "(or use --run when that path is implemented)")

    vparca_state = _load_pickle(args.vparca_state)
    original     = _load_pickle(args.original_sim_data)
    vparca       = _sim_data_like_from_vparca_state(vparca_state)

    vt = _load_json(args.vparca_runtimes)   if args.vparca_runtimes   else {}
    ot = _load_json(args.original_runtimes) if args.original_runtimes else {}

    build_report(vparca, original, vt, ot, args.output,
                 tol_pass=args.tol_pass, tol_warn=args.tol_warn)


if __name__ == '__main__':
    main()
