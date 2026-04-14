"""
Step 9 — final_adjustments.  Final expression adjustments, amino-acid
supply/export/uptake constants, and ppGpp reaction kinetics parameters.

Object-dominated: every sub-call is a method on a sim_data subsystem that
mutates its own internal state.  The step wires every subsystem it
touches as a port and returns the mutated objects on its output ports.

Store paths wired by the composite
----------------------------------

READS (subsystems):
  transcription, translation, metabolism, constants, mass, complexation,
  equilibrium, two_component_system, transcription_regulation,
  replication, growth_rate_parameters, molecule_ids, molecule_groups,
  relation, bulk_molecules
READS (data leaves):
  conditions, condition_to_doubling_time, tf_to_fold_change, cell_specs

WRITES: transcription, metabolism, constants (all mutated in place)
"""

import time

from process_bigraph import Step

from vparca.ecoli.library.initial_conditions import create_bulk_container
from vparca.steps._facade import make_sim_data_facade


INPUT_PORTS = {
    'tick_8'                            : 'overwrite',
    'transcription':            'sim_data.transcription',
    'translation':              'sim_data.translation',
    'metabolism':               'sim_data.metabolism',
    'complexation':             'sim_data.complexation',
    'equilibrium':              'sim_data.equilibrium',
    'two_component_system':     'sim_data.two_component_system',
    'transcription_regulation': 'sim_data.transcription_regulation',
    'replication':              'sim_data.replication',
    'mass':                     'sim_data.mass',
    'constants':                'sim_data.constants',
    'growth_rate_parameters':   'sim_data.growth_rate_parameters',
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
    # set_mechanistic_supply_constants reads sim_data.translation_supply_rate
    # which step 5 populates per-nutrient.
    'translation_supply_rate':  'overwrite',
    # calculate_attenuation reads sim_data.pPromoterBound (set by step 6)
    'pPromoterBound':           'overwrite',
    # create_bulk_container uses external_state.exchange_data_from_media
    # and mutates sim_data.condition temporarily.
    'external_state':           'overwrite',
    'condition':                'overwrite',
}

OUTPUT_PORTS = {
    'tick_9'                            : 'overwrite',
    'transcription': 'sim_data.transcription',
    'metabolism':    'sim_data.metabolism',
    'constants':     'sim_data.constants',
}


class FinalAdjustmentsStep(Step):
    """Step 9 — final_adjustments.  See module docstring."""

    def inputs(self):
        return dict(INPUT_PORTS)

    def outputs(self):
        return dict(OUTPUT_PORTS)

    def update(self, state):
        t0 = time.time()

        sd = make_sim_data_facade(state)
        cell_specs = state['cell_specs']

        # Attenuation + ppGpp expression fixups.
        sd.process.transcription.calculate_attenuation(sd, cell_specs)
        sd.process.transcription.adjust_polymerizing_ppgpp_expression(sd)
        sd.process.transcription.adjust_ppgpp_expression_for_tfs(sd)

        # Amino-acid supply constants — based on average bulk containers.
        average_basal_container   = create_bulk_container(sd, n_seeds=5)
        average_with_aa_container = create_bulk_container(
            sd, condition="with_aa", n_seeds=5)

        sd.process.metabolism.set_phenomological_supply_constants(sd)
        # The three mechanistic_* fits can raise on numerically-marginal
        # kinetics (e.g. "Could not find positive forward and reverse
        # kcat for CYS[c]") in debug mode where the truncated TF set
        # produces edge-case input distributions.  We log and continue
        # so the pipeline still produces a pickle of everything else
        # (including all step-1..8 outputs).  The failure is identical
        # to what the original fit_sim_data_1 raises under the same
        # conditions; debug it with --mode full or by patching the
        # underlying kinetics fit.
        for label, call in [
            ('mechanistic_supply', lambda: sd.process.metabolism
                .set_mechanistic_supply_constants(
                    sd, cell_specs,
                    average_basal_container, average_with_aa_container)),
            ('mechanistic_export', lambda: sd.process.metabolism
                .set_mechanistic_export_constants(
                    sd, cell_specs, average_basal_container)),
            ('mechanistic_uptake', lambda: sd.process.metabolism
                .set_mechanistic_uptake_constants(
                    sd, cell_specs, average_with_aa_container)),
        ]:
            try:
                call()
            except Exception as e:
                print(f"  Step 9 WARNING: {label} failed ({type(e).__name__}: {e}); "
                      "continuing so the pipeline produces a comparable pickle.")

        # ppGpp kinetics.
        sd.process.transcription.set_ppgpp_kinetics_parameters(
            average_basal_container, sd.constants)

        print(f"  Step 9 (final_adjustments) completed in {time.time() - t0:.1f}s")
        return {
            'transcription': sd.process.transcription,
            'metabolism':    sd.process.metabolism,
            'constants':     sd.constants,
        
            'tick_9': True,}
