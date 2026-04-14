"""
Step 7 — adjust_promoters.  Tune ligand concentrations and build the RNAP
recruitment parameters from the fitted promoter binding probabilities.

Store paths wired by the composite
----------------------------------

READS (subsystems):
  transcription, transcription_regulation, equilibrium, metabolism,
  replication, mass, constants, molecule_ids, molecule_groups, bulk_molecules
READS (data leaves): cell_specs, conditions, condition_to_doubling_time

WRITES: transcription_regulation (basal_prob + delta_prob, mutated),
        metabolism (molecule_set_amounts, mutated),
        equilibrium (reverse_rates, mutated).
"""

import time

from process_bigraph import Step

from v2parca.promoter_fitting import (
    fitLigandConcentrations,
    calculateRnapRecruitment,
)
from v2parca.steps._facade import make_sim_data_facade


INPUT_PORTS = {
    'tick_6'                            : 'overwrite',
    'transcription':            'sim_data.transcription',
    'transcription_regulation': 'sim_data.transcription_regulation',
    'equilibrium':              'sim_data.equilibrium',
    'two_component_system':     'sim_data.two_component_system',
    'metabolism':               'sim_data.metabolism',
    'replication':              'sim_data.replication',
    'mass':                     'sim_data.mass',
    'constants':                'sim_data.constants',
    'molecule_ids':             'overwrite',
    'molecule_groups':          'overwrite',
    'relation':                 'overwrite',
    'getter':                   'overwrite',
    'bulk_molecules':           'overwrite',
    'sim_data_root':            'overwrite',
    'conditions':               'overwrite',
    'condition_to_doubling_time': 'overwrite',
    'tf_to_active_inactive_conditions': 'overwrite',
    'tf_to_fold_change':        'overwrite',
    'tf_to_direction':          'overwrite',
    'condition_active_tfs':     'overwrite',
    'condition_inactive_tfs':   'overwrite',
    'cell_specs':               'overwrite',
    # From step 6; fitLigandConcentrations reads sim_data.pPromoterBound.
    'pPromoterBound':           'overwrite',
}

OUTPUT_PORTS = {
    'tick_7'                            : 'overwrite',
    'transcription_regulation': 'sim_data.transcription_regulation',
    'metabolism':               'sim_data.metabolism',
    'equilibrium':              'sim_data.equilibrium',
}


class AdjustPromotersStep(Step):
    """Step 7 — adjust_promoters.  See module docstring."""

    def inputs(self):
        return dict(INPUT_PORTS)

    def outputs(self):
        return dict(OUTPUT_PORTS)

    def update(self, state):
        t0 = time.time()

        sd = make_sim_data_facade(state)
        cell_specs = state['cell_specs']

        fitLigandConcentrations(sd, cell_specs)
        basal_prob, delta_prob = calculateRnapRecruitment(sd, cell_specs)

        sd.process.transcription_regulation.basal_prob = basal_prob
        sd.process.transcription_regulation.delta_prob = delta_prob

        print(f"  Step 7 (adjust_promoters) completed in {time.time() - t0:.1f}s")
        return {
            'transcription_regulation': sd.process.transcription_regulation,
            'metabolism':               sd.process.metabolism,
            'equilibrium':              sd.process.equilibrium,
        
            'tick_7': True,}
