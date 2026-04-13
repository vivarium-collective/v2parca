# vParCa

A **process-bigraph migration of the E. coli Parameter Calculator (ParCa)**.

The ParCa transforms raw knowledge-base flat files into a fitted `SimulationDataEcoli`
object that parameterizes the whole-cell model. Traditionally it runs as a single
monolithic `fitSimData_1()` that mutates a deeply nested `sim_data` object across
nine sequential stages. vParCa decomposes each stage into a `process_bigraph.Step`
with **explicit, named input and output ports**, so the data flowing between stages
is no longer hidden inside `sim_data` — it's wired through the Composite.

This is the ParCa analogue of [v2ecoli](https://github.com/eagmon/v2ecoli), which
did the same decomposition for the simulation processes.

## Pipeline overview

```
raw_data ─▶ [1 initialize]        ─▶ state_1
state_1  ─▶ [2 input_adjustments] ─▶ state_2        (PURE — Extract / Compute / Merge)
state_2  ─▶ [3 basal_specs]       ─▶ state_3 + named outs
state_3  ─▶ [4 tf_condition_specs]─▶ state_4 + named outs
state_4  ─▶ [5 fit_condition]     ─▶ state_5 + named outs
state_5  ─▶ [6 promoter_binding]  ─▶ state_6 + named outs
state_6  ─▶ [7 adjust_promoters]  ─▶ state_7 + named outs
state_7  ─▶ [8 set_conditions]    ─▶ state_8        (PURE — Extract / Compute / Merge)
state_8  ─▶ [9 final_adjustments] ─▶ state_9
```

Each stage in `vparca/stage_0N_*.py` exposes three functions:

| function | role |
| --- | --- |
| `extract_input(sim_data, cell_specs, **kwargs) -> StageInput` | pull only the fields the stage reads |
| `compute_*(inp: StageInput) -> StageOutput` | the core math — pure where possible |
| `merge_output(sim_data, cell_specs, out: StageOutput)` | write back only the fields the stage produces |

The dataclasses in `vparca/_types.py` make the real I/O of every stage explicit
and inspectable without reading through 4400 lines of ParCa code.

Purity legend (see `vparca/__init__.py` and `vparca/DATA_FLOW.md`):

- **PURE** — compute has no `sim_data` access; fully testable with synthetic data
- **READ-ONLY** — compute reads `sim_data` via ref but does not mutate
- **COUPLED** — compute still mutates `sim_data` via ref (future refactor target)

Pure stages (2 and 8) are decomposed in the Composite into an
`Extract → Compute → Merge` triplet so the Compute step has *only* explicit
typed ports.

## Layout

```
vparca/
  __init__.py              stage registry + purity table
  DATA_FLOW.md             annotated data-flow documentation
  _types.py                Input/Output dataclasses for stages 2–9
  _math.py                 pure math helpers
  _fitting.py              sim_data-reading fitting helpers
  _shared.py               backward-compatible re-export shim
  parca_types.py           process-bigraph ParcaState type + registration
  steps.py                 Step classes (one per stage + Extract/Merge for pure stages)
  composite.py             build_parca_composite() / run_parca()
  parca_promoter_fitting.py   matrix builders + CVXPY optimization for stages 6/7
  stage_02_input_adjustments.py
  stage_03_basal_specs.py
  stage_04_tf_condition_specs.py
  stage_05_fit_condition.py
  stage_06_promoter_binding.py
  stage_07_adjust_promoters.py
  stage_08_set_conditions.py
  stage_09_final_adjustments.py
tests/
  test_parca_stage_02.py ... test_parca_stage_09.py
scripts/
  parca_bigraph.py         CLI: run ParCa as a process-bigraph Composite
  parca_workflow.py        CLI: run ParCa stages as plain function calls
```

## Dependency on vEcoli / vivarium-ecoli

vParCa relies on vEcoli for the non-ParCa substrate: `SimulationDataEcoli`,
the raw knowledge base, `wholecell.utils`, and `ecoli.library`. These are not
vendored — install vEcoli (or vivarium-ecoli) and make its root importable.

```bash
# one-time — install vEcoli so its modules are on PYTHONPATH
cd /path/to/vEcoli           # or vivarium-ecoli
pip install -e .

# then install vParCa
cd /path/to/vParCa
pip install -e .[test]
```

vParCa imports `reconstruction.ecoli.simulation_data`, `wholecell.*`, and
`ecoli.library.*` directly from the installed vEcoli; everything else resolves
inside the `vparca` package.

## Running the pipeline

```bash
# via process-bigraph Composite
python scripts/parca_bigraph.py --mode fast --cpus 4

# as plain function calls (for comparison)
python scripts/parca_workflow.py --mode fast --cpus 4
```

`run_parca(raw_data, **kwargs)` in `vparca.composite` returns the fitted
`SimulationDataEcoli` after the Composite's Steps execute in DAG order.

## Tests

```bash
pytest tests/
```

Each `tests/test_parca_stage_0N.py` exercises `extract_input`, `compute_*`, and
`merge_output` for one stage with synthetic fixtures so the pure slice of each
stage is validated without requiring a full ParCa run.

## Origin

This repo was extracted from
`vivarium-ecoli`, branch `refactor-parca-again` — the preliminary work that
split `fitSimData_1` into modular, typed stages. Intra-parca imports were
rewritten from `reconstruction.ecoli.parca.*` to `vparca.*`; all other imports
into vEcoli internals are preserved.
