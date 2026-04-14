"""
Stage 9: final_adjustments — Apply final expression adjustments, set
amino acid supply constants, and configure ppGpp kinetics.

Three public functions:
    extract_input(sim_data, cell_specs, **kwargs) -> FinalAdjustmentsInput
    compute_final_adjustments(inp) -> FinalAdjustmentsOutput
    merge_output(sim_data, cell_specs, out)

This stage:
1. Adjusts expression for RNA attenuation
2. Adjusts ppGpp-regulated expression after conditions have been fit
3. Creates bulk containers for basal and with_aa conditions
4. Sets phenomenological and mechanistic amino acid supply/export/uptake
   constants
5. Sets ppGpp reaction kinetics parameters

NOTE: compute_final_adjustments mutates sim_data_ref extensively via
deep process-object method calls.  There are no cell_specs writes.
merge_output is a no-op.
"""

from vparca.ecoli.library.initial_conditions import create_bulk_container

from vparca.types import (
    FinalAdjustmentsInput,
    FinalAdjustmentsOutput,
)


import time
from process_bigraph import Step

from vparca.state import ParcaState

# ============================================================================
# Extract / Merge
# ============================================================================


def extract_input(sim_data, cell_specs, **kwargs) -> FinalAdjustmentsInput:
    """Pull sim_data and cell_specs references for final adjustments."""
    return FinalAdjustmentsInput(
        sim_data_ref=sim_data,
        cell_specs_ref=cell_specs,
    )


def merge_output(sim_data, cell_specs, out: FinalAdjustmentsOutput):
    """No-op: all mutations are applied during compute via sim_data_ref."""
    pass


# ============================================================================
# Compute
# ============================================================================


def compute_final_adjustments(inp: FinalAdjustmentsInput) -> FinalAdjustmentsOutput:
    """Run the full final_adjustments stage.

    This function mutates inp.sim_data_ref extensively via deep
    process-object method calls.  There are no extractable outputs.
    """
    sim_data = inp.sim_data_ref
    cell_specs = inp.cell_specs_ref

    # Adjust expression for RNA attenuation
    sim_data.process.transcription.calculate_attenuation(sim_data, cell_specs)

    # Adjust ppGpp regulated expression after conditions have been fit
    sim_data.process.transcription.adjust_polymerizing_ppgpp_expression(sim_data)
    sim_data.process.transcription.adjust_ppgpp_expression_for_tfs(sim_data)

    # Set supply constants for amino acids based on condition supply requirements
    average_basal_container = create_bulk_container(sim_data, n_seeds=5)
    average_with_aa_container = create_bulk_container(
        sim_data, condition="with_aa", n_seeds=5
    )
    sim_data.process.metabolism.set_phenomological_supply_constants(sim_data)
    sim_data.process.metabolism.set_mechanistic_supply_constants(
        sim_data, cell_specs, average_basal_container, average_with_aa_container
    )
    sim_data.process.metabolism.set_mechanistic_export_constants(
        sim_data, cell_specs, average_basal_container
    )
    sim_data.process.metabolism.set_mechanistic_uptake_constants(
        sim_data, cell_specs, average_with_aa_container
    )

    # Set ppGpp reaction parameters
    sim_data.process.transcription.set_ppgpp_kinetics_parameters(
        average_basal_container, sim_data.constants
    )

    return FinalAdjustmentsOutput()



# ---------------------------------------------------------------------------
# Stage 9: Final Adjustments (COUPLED)
# ---------------------------------------------------------------------------

class FinalAdjustmentsStep(Step):
    """Stage 9: Apply final expression adjustments, set amino acid supply
    constants, and configure ppGpp kinetics."""

    def inputs(self):
        return {'state': 'parca_state'}

    def outputs(self):
        return {'state': 'parca_state'}

    def update(self, state):
        t0 = time.time()
        parca_state = state['state']
        sim_data = parca_state.sim_data
        cell_specs = parca_state.cell_specs

        inp = extract_input(sim_data, cell_specs)
        out = compute_final_adjustments(inp)
        merge_output(sim_data, cell_specs, out)

        print(f"  Stage 9 (final_adjustments) completed in {time.time() - t0:.1f}s")
        return {
            'state': ParcaState(sim_data=sim_data, cell_specs=cell_specs),
        }

