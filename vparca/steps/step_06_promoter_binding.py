"""
Step 6 — promoter_binding.  Fit transcription-factor binding probabilities
and their effect on RNA synthesis.

Uses CVXPY/ECOS to fit recruitment parameters ``r`` and binding
probabilities ``P`` such that computed RNA synthesis matches the measured
condition-specific synth probabilities.  Mutates the ``transcription`` and
``transcription_regulation`` subsystems via the shared sub-function
``fitPromoterBoundProbability``; writes ``r_vector`` / ``r_columns`` into
``cell_specs["basal"]``.

Store paths wired by the composite
----------------------------------

READS (subsystems):
  transcription, transcription_regulation, equilibrium, replication,
  mass, constants, molecule_ids, molecule_groups, bulk_molecules
READS (data leaves): cell_specs, conditions, condition_to_doubling_time

WRITES: transcription, transcription_regulation (mutated), cell_specs
"""

import time

from process_bigraph import Step

from vparca.promoter_fitting import fitPromoterBoundProbability
from vparca.steps._facade import make_sim_data_facade


INPUT_PORTS = {
    'tick_5'                            : 'overwrite',
    'transcription':            'sim_data.transcription',
    'transcription_regulation': 'sim_data.transcription_regulation',
    'equilibrium':              'sim_data.equilibrium',
    'two_component_system':     'sim_data.two_component_system',
    'replication':              'sim_data.replication',
    'mass':                     'sim_data.mass',
    'constants':                'sim_data.constants',
    'molecule_ids':             'overwrite',
    'molecule_groups':          'overwrite',
    'relation':                 'overwrite',
    'getter':                   'overwrite',
    'bulk_molecules':           'overwrite',
    'conditions':               'overwrite',
    'condition_to_doubling_time': 'overwrite',
    'tf_to_active_inactive_conditions': 'overwrite',
    'tf_to_fold_change':        'overwrite',
    'tf_to_direction':          'overwrite',
    'condition_active_tfs':     'overwrite',
    'condition_inactive_tfs':   'overwrite',
    'cell_specs':               'overwrite',
}

OUTPUT_PORTS = {
    'tick_6'                            : 'overwrite',
    'transcription':            'sim_data.transcription',
    'transcription_regulation': 'sim_data.transcription_regulation',
    'cell_specs':               'overwrite',
    # fitPromoterBoundProbability sets sim_data.pPromoterBound as a
    # dynamic top-level attr; propagate it so step 7 can read it.
    'pPromoterBound':           'overwrite',
}


class PromoterBindingStep(Step):
    """Step 6 — promoter_binding.  See module docstring."""

    def inputs(self):
        return dict(INPUT_PORTS)

    def outputs(self):
        return dict(OUTPUT_PORTS)

    def update(self, state):
        t0 = time.time()

        sd = make_sim_data_facade(state)
        cell_specs = dict(state['cell_specs'])

        print("Fitting promoter binding")
        r_vector, r_columns = fitPromoterBoundProbability(sd, cell_specs)

        cell_specs.setdefault("basal", {})
        cell_specs["basal"]["r_vector"]  = r_vector
        cell_specs["basal"]["r_columns"] = r_columns

        print(f"  Step 6 (promoter_binding) completed in {time.time() - t0:.1f}s")
        return {
            'transcription':            sd.process.transcription,
            'transcription_regulation': sd.process.transcription_regulation,
            'cell_specs':               cell_specs,
            'pPromoterBound':           getattr(sd, 'pPromoterBound', {}),
            'tick_6': True,
        }
