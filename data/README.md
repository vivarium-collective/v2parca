# Pre-computed ParCa state

`parca_state.pkl.gz` is the final composite state of the full v2parca
pipeline in `--mode fast` (debug=True, reduced TF condition set). Produced
by `scripts/parca_bigraph.py` and gzipped for the repo. Consumers can use
this as a drop-in `sim_data`-shaped dict without re-running the 70-minute
step 5.

## Contents

- 24 top-level stores: the nine `process/*` subsystem objects
  (transcription, translation, metabolism, rna_decay, complexation,
  equilibrium, two_component_system, transcription_regulation, replication)
  plus top-level dataclasses (mass, constants, growth_rate_parameters,
  adjustments, molecule_groups, molecule_ids, relation, getter,
  external_state), top-level data dicts (conditions,
  condition_to_doubling_time, condition_active_tfs, condition_inactive_tfs,
  expected_dry_mass_increase_dict, pPromoterBound), plus
  `internal_state/bulk_molecules` and `cell_specs`.
- 7 conditions in `cell_specs`: `basal`, `with_aa`, `acetate`, `succinate`,
  `no_oxygen`, `CPLX-125__active`, `CPLX-125__inactive`.
- `runtimes.json`: per-step wall-clock timings from the run that
  produced this pickle.

## Loading

```python
import gzip
import pickle

with gzip.open('data/parca_state.pkl.gz', 'rb') as f:
    state = pickle.load(f)

# Access any store:
transcription = state['process']['transcription']
cell_specs    = state['cell_specs']
mass          = state['mass']
```

Or via the convenience loader:

```python
from v2parca.data_loader import load_parca_state
state = load_parca_state()
```

## Provenance

- Pipeline: `scripts/parca_bigraph.py --mode fast --cpus 2`
- Duration: 71.6 min end-to-end (step 5: 70 min; steps 1–4 + 6–9: ~1.5 min)
- debug=True, operons_on=True, remove_rrna_operons=False, stable_rrna=False

To regenerate from scratch (takes ~70 minutes):

```bash
python scripts/parca_bigraph.py --mode fast --cpus 2
gzip -k out/sim_data/parca_state.pkl && mv out/sim_data/parca_state.pkl.gz data/
cp out/sim_data/runtimes.json data/
```
