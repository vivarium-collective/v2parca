# v2parca

A **process-bigraph migration of the E. coli Parameter Calculator (ParCa)**.

The ParCa transforms raw knowledge-base flat files into a fitted `SimulationDataEcoli`
object that parameterizes the whole-cell model. Traditionally it runs as a single
monolithic `fitSimData_1()` that mutates a deeply nested `sim_data` object across
nine sequential stages. v2parca decomposes each stage into a `process_bigraph.Step`
with **explicit, named input and output ports**, so the data flowing between stages
is no longer hidden inside `sim_data` — it's wired through the Composite.

This is the ParCa analogue of [v2ecoli](https://github.com/eagmon/v2ecoli), which
did the same decomposition for the simulation processes.

## Architecture

v2parca is a **port-first, nested-store** pipeline.  The composite's state
is a bigraph tree that mirrors `SimulationDataEcoli`'s own structure —
subsystem objects (`Transcription`, `Mass`, `Constants`, …) live at
natural paths like `process/transcription` or `constants`, and pure-data
top-level dicts (`conditions`, `tf_to_active_inactive_conditions`,
`cell_specs`, …) live at sibling paths.  Each Step declares an explicit
port manifest — one port per subsystem or leaf it touches — and the
composite wires each port to its store path.  No `sim_data` or
`cell_specs` blob travels.

Pipeline overview:

```
             raw_data (config)
                  │
                  ▼
 ┌────────────────────────────────────────┐
 │ Step 1  initialize + scatter           │
 └────────────────────────────────────────┘
       │   │   │   │   (18 subsystem + 9 top-level leaves, + tick_1)
       ▼   ▼   ▼   ▼
 ┌────────────────────────────────────────┐
 │  nested bigraph store                  │
 │    process/transcription               │
 │    process/translation                 │
 │    process/metabolism ...              │
 │    mass / constants / ...              │
 │    conditions / cell_specs / ...       │
 └────────────────────────────────────────┘
       ▲   ▲   ▲   ▲
       │   │   │   │   (each Step wires port → store path)
 ┌────────────────────────────────────────┐
 │ Step 2  input_adjustments              │  tick_1 ──▶ tick_2
 │ Step 3  basal_specs                    │  tick_2 ──▶ tick_3
 │ Step 4  tf_condition_specs             │  tick_3 ──▶ tick_4
 │ Step 5  fit_condition                  │  tick_4 ──▶ tick_5
 │ Step 6  promoter_binding               │  tick_5 ──▶ tick_6
 │ Step 7  adjust_promoters               │  tick_6 ──▶ tick_7
 │ Step 8  set_conditions                 │  tick_7 ──▶ tick_8
 │ Step 9  final_adjustments              │  tick_8 ──▶ tick_9
 └────────────────────────────────────────┘
```

Each Step is a thin `process_bigraph.Step` subclass in
`v2parca/steps/step_NN_*.py` with two module-level dicts — `INPUT_PORTS`
and `OUTPUT_PORTS` — and an `update(state)` method that reads port values
from `state`, calls the ParCa sub-functions through a `SimpleNamespace`
facade built by `v2parca/steps/_facade.make_sim_data_facade`, and returns
a dict keyed by output-port name.  Subsystem objects carry their own
mutations back out via their output ports.

**Tick ordering.** Because several Steps read and write overlapping
subsystem objects, process-bigraph can't infer a total order from the
data-level wires alone.  Opaque `tick_0..tick_9` leaves are wired as
`tick_{N-1}` input → `tick_N` output per step, forcing a strict serial
execution Step 1 → Step 9 as part of the Composite's initial DAG fire.

**Bigraph-schema types** (`v2parca/schema.py`) register Overwrite
subclasses named `sim_data.transcription`, `sim_data.mass`, etc., so the
port manifests document which kind of Python object lives at each
subsystem leaf.  Today they're behaviorally identical to `overwrite`;
future work can hook dispatch (serialize, diff) per type.

## Layout

Everything lives under a single Python namespace: `v2parca`.

```
v2parca/
  __init__.py
  composite.py                   # build_parca_composite() / run_parca()
                                 # + STORE_PATH (port-name → nested-store path)
  schema.py                      # bigraph-schema type registry (subsystem leaves)
  fitting.py                     # sim_data-reading fitting helpers (pure math
                                 #   + expressionConverge / Km / mass rescaling)
  promoter_fitting.py            # matrix builders + CVXPY for steps 6/7
  trna_charging.py               # calculate_trna_charging + constants
  steps/
    __init__.py                  # ALL_STEP_CLASSES registry
    _facade.py                   # make_sim_data_facade(ports) → SimpleNamespace
    step_01_initialize.py        # scatter: sim_data.initialize → 18 subsystems
    step_02_input_adjustments.py
    step_03_basal_specs.py
    step_04_tf_condition_specs.py
    step_05_fit_condition.py
    step_06_promoter_binding.py
    step_07_adjust_promoters.py
    step_08_set_conditions.py
    step_09_final_adjustments.py

  # Vendored vEcoli substrate — all under v2parca/ so there's one namespace
  reconstruction/
    spreadsheets.py
    ecoli/
      knowledge_base_raw.py      # KnowledgeBaseEcoli — loads flat/
      simulation_data.py         # SimulationDataEcoli
      dataclasses/               # sim_data subclasses
      flat/                      # RAW DATA — KB TSVs + media/mass_fractions/etc.
      scripts/update_biocyc_files.py
  wholecell/                     # units, fitting, parallelization, …
  ecoli/
    library/                     # schema.py + initial_conditions.py (only)

tests/                           # test_ports_and_wiring.py  (fast static checks)
scripts/                         # parca_bigraph.py, parca_workflow.py
docs/                            # PORT_MAP.md, DATA_FLOW.md
```

Imports: `v2parca.reconstruction.ecoli.knowledge_base_raw`,
`v2parca.wholecell.utils`, `v2parca.ecoli.library.schema`, etc.  The substrate
is not a separate installable package — it's part of v2parca.

## Running from raw data

v2parca is self-contained: the raw KB (`reconstruction/ecoli/flat/`) and all vEcoli
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
from v2parca.reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from v2parca.composite import run_parca

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

`tests/test_ports_and_wiring.py` (~16 tests, ~5s) validates the *static*
port-first architecture without running the ParCa:

- every port a Step declares is covered by `STORE_PATH`
- Step 1 scatters every subsystem / leaf any downstream Step reads
- tick chain `tick_0 → tick_1 → … → tick_9` is serial and complete
- every `sim_data.*` schema type is registered
- `build_parca_composite(...)` constructs without error

These catch port-manifest drift immediately rather than surfacing as an
`AttributeError` 30 minutes into a real pipeline run.

## Comparison with the original ParCa

The comparison harness — which runs the original `fitSimData_1` from
`vivarium-ecoli` and diffs `sim_data` after each stage — is the **only** part of
v2parca that imports vEcoli. It's optional; install it with:

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
from `reconstruction.ecoli.parca.*` to `v2parca.*`; everything else imports the
vendored substrate.
