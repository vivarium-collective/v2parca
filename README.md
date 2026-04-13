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

Each stage module in `vparca/stages/` exposes three functions plus a Step class:

| function / class | role |
| --- | --- |
| `extract_input(sim_data, cell_specs, **kwargs) -> StageInput` | pull only the fields the stage reads |
| `compute_*(inp: StageInput) -> StageOutput` | the core math — pure where possible |
| `merge_output(sim_data, cell_specs, out: StageOutput)` | write back only the fields the stage produces |
| `<Stage>Step(Step)` | the process-bigraph wrapper that runs in the Composite |

The dataclasses in `vparca/types.py` make the real I/O of every stage explicit
and inspectable without reading through 4400 lines of ParCa code.

Purity legend (see `vparca/__init__.py` and `docs/DATA_FLOW.md`):

- **PURE** — compute has no `sim_data` access; fully testable with synthetic data
- **READ-ONLY** — compute reads `sim_data` via ref but does not mutate
- **COUPLED** — compute still mutates `sim_data` via ref (future refactor target)

Pure stages (2 and 8) are decomposed in the Composite into an
`Extract → Compute → Merge` triplet so the Compute step has *only* explicit
typed ports.

## Layout

```
vparca/                          # migration code (process-bigraph Steps)
  composite.py                   # build_parca_composite() / run_parca()
  state.py                       # ParcaState bigraph-schema type
  types.py                       # Input/Output dataclasses per stage
  fitting.py                     # pure math + sim_data-reading helpers
  promoter_fitting.py            # matrix builders + CVXPY for stages 6/7
  stages/
    __init__.py                  # ALL_STEP_CLASSES registry
    stage_01_initialize.py
    stage_02_input_adjustments.py         (PURE — Extract/Compute/Merge)
    stage_03_basal_specs.py
    stage_04_tf_condition_specs.py
    stage_05_fit_condition.py               (READ-ONLY)
    stage_06_promoter_binding.py
    stage_07_adjust_promoters.py
    stage_08_set_conditions.py            (PURE — Extract/Compute/Merge)
    stage_09_final_adjustments.py

reconstruction/                  # vendored vEcoli substrate
  ecoli/
    knowledge_base_raw.py        # KnowledgeBaseEcoli — loads flat/
    simulation_data.py           # SimulationDataEcoli
    dataclasses/                 # sim_data subclasses
    flat/                        # RAW DATA — all KB TSVs + media/mass_fractions/etc.
    scripts/update_biocyc_files.py
  spreadsheets.py
wholecell/                       # vendored — units, fitting, parallelization, etc.
ecoli/                           # vendored (trimmed)
  library/                       # schema + initial_conditions
  processes/                     # only the 5 modules needed by initial_conditions

tests/                           # test_parca_stage_02.py ... test_parca_stage_09.py
scripts/                         # parca_bigraph.py, parca_workflow.py
docs/                            # DATA_FLOW.md
```

## Running from raw data

vParCa is self-contained: the raw KB (`reconstruction/ecoli/flat/`) and all vEcoli
modules the ParCa needs are vendored inside this repo. You do **not** need
vEcoli or vivarium-ecoli installed to run the ParCa.

```bash
pip install -e .[test]

# full pipeline as a process-bigraph Composite
python scripts/parca_bigraph.py --mode fast --cpus 4

# same pipeline as plain function calls (side-by-side reference)
python scripts/parca_workflow.py --mode fast --cpus 4
```

Programmatic entry point:

```python
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from vparca.composite import run_parca

raw = KnowledgeBaseEcoli(
    operons_on=True,
    remove_rrna_operons=False,
    remove_rrff=False,
    stable_rrna=False,
)
sim_data = run_parca(raw, cpus=4, debug=True)
```

## Tests

```bash
pytest tests/
```

Each `tests/test_parca_stage_0N.py` exercises `extract_input`, `compute_*`, and
`merge_output` for one stage with synthetic fixtures, so the pure slice of each
stage is validated without requiring a full ParCa run.

## Comparison with the original ParCa

The comparison harness — which runs the original `fitSimData_1` from
`vivarium-ecoli` and diffs `sim_data` after each stage — is the **only** part of
vParCa that imports vEcoli. It's optional; install it with:

```bash
pip install -e .[compare]
# and make vivarium-ecoli importable (editable install, or PYTHONPATH)
```

## Vendoring notes

- `reconstruction/`, `wholecell/`, `ecoli/` are lifted from
  `vivarium-ecoli@refactor-parca-again` and left at their original import paths,
  so `from reconstruction.ecoli.simulation_data import SimulationDataEcoli` works
  the same as upstream.
- `ecoli/__init__.py` and `ecoli/processes/__init__.py` were **replaced** with
  empty stubs; upstream those files registered vivarium emitters, dividers,
  serializers, and process classes that belong to the simulation side and have
  no role in the ParCa.
- `ecoli/processes/metabolism.py` had its `test_metabolism_listener` function
  stripped because it imported `ecoli.experiments.ecoli_master_sim`, which is
  not vendored.
- The raw knowledge base under `reconstruction/ecoli/flat/` is included in full
  (~13 MB of TSVs + nested `adjustments/`, `condition/media/`, `mass_fractions/`,
  `rna_seq_data/`, `trna_data/`, `rrna_options/`, `new_gene_data/`, etc.).
- `reconstruction/ecoli/scripts/update_biocyc_files.py` is vendored so raw data
  can be regenerated from BioCyc within this repo.

## Origin

Extracted from `vivarium-ecoli`, branch `refactor-parca-again`, which split
`fitSimData_1` into modular, typed stages. Intra-parca imports were rewritten
from `reconstruction.ecoli.parca.*` to `vparca.*`; everything else imports the
vendored substrate.
