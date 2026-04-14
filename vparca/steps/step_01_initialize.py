"""Stage 1 — Initialize sim_data from raw_data.

This stage has no extract/compute/merge split: InitializeStep
directly invokes ``sim_data.initialize(raw_data=...)`` which populates
the whole nested sim_data from the KB flat files.
"""

import time

from process_bigraph import Step

from vparca.state import ParcaState


# ---------------------------------------------------------------------------
# Stage 1: Initialize (COUPLED)
# ---------------------------------------------------------------------------

class InitializeStep(Step):
    """Stage 1: Initialize sim_data from raw_data."""

    config_schema = {
        'basal_expression_condition': {
            '_type': 'string',
            '_default': 'M9 Glucose minus AAs',
        },
    }

    def inputs(self):
        return {
            'state': 'parca_state',
            'raw_data': 'overwrite',
        }

    def outputs(self):
        return {
            'state': 'parca_state',
        }

    def update(self, state):
        t0 = time.time()
        parca_state = state['state']
        sim_data = parca_state.sim_data
        sim_data.initialize(
            raw_data=state['raw_data'],
            basal_expression_condition=self.config.get(
                'basal_expression_condition', 'M9 Glucose minus AAs'),
        )
        print(f"  Stage 1 (initialize) completed in {time.time() - t0:.1f}s")
        return {
            'state': ParcaState(sim_data=sim_data,
                                cell_specs=parca_state.cell_specs),
        }


