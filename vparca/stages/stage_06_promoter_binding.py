"""
Stage 6: promoter_binding — Fit transcription factor binding probabilities
and their effects on RNA synthesis.

Three public functions:
    extract_input(sim_data, cell_specs, **kwargs) -> PromoterBindingInput
    compute_promoter_binding(inp) -> PromoterBindingOutput
    merge_output(sim_data, cell_specs, out)

This stage:
1. Calculates initial TF-promoter binding probabilities from bulk average
   counts of TFs and ligands for each condition
2. Uses convex optimization (CVXPY/ECOS) to fit parameters alpha and r
   such that computed RNA synthesis probabilities match measured values
3. Iteratively optimizes both the binding probabilities (P) and the
   recruitment parameters (r) until convergence

NOTE: compute_promoter_binding mutates sim_data_ref because
fitPromoterBoundProbability writes pPromoterBound and rna_synth_prob
to sim_data.  merge_output only writes cell_specs.
"""

from vparca.types import (
    PromoterBindingInput,
    PromoterBindingOutput,
)
from vparca.promoter_fitting import (

    fitPromoterBoundProbability,
)


import time
from process_bigraph import Step

from vparca.state import ParcaState

# ============================================================================
# Extract / Merge
# ============================================================================


def extract_input(sim_data, cell_specs, **kwargs) -> PromoterBindingInput:
    """Pull sim_data and cell_specs references for promoter fitting."""
    return PromoterBindingInput(
        sim_data_ref=sim_data,
        cell_specs_ref=cell_specs,
    )


def merge_output(sim_data, cell_specs, out: PromoterBindingOutput):
    """Write computed results into cell_specs.

    sim_data mutations (pPromoterBound, rna_synth_prob) are already
    applied by compute_promoter_binding via sim_data_ref.
    """
    cell_specs["basal"]["r_vector"] = out.r_vector
    cell_specs["basal"]["r_columns"] = out.r_columns


# ============================================================================
# Compute
# ============================================================================


def compute_promoter_binding(inp: PromoterBindingInput) -> PromoterBindingOutput:
    """Run the full promoter_binding stage.

    This function mutates inp.sim_data_ref as a side effect because
    fitPromoterBoundProbability sets pPromoterBound and updates
    rna_synth_prob on sim_data.
    """
    sim_data = inp.sim_data_ref
    cell_specs = inp.cell_specs_ref

    print("Fitting promoter binding")

    r_vector, r_columns = fitPromoterBoundProbability(sim_data, cell_specs)

    return PromoterBindingOutput(
        r_vector=r_vector,
        r_columns=r_columns,
    )



# ---------------------------------------------------------------------------
# Stage 6: Promoter Binding (COUPLED)
# ---------------------------------------------------------------------------

class PromoterBindingStep(Step):
    """Stage 6: Fit transcription factor binding probabilities."""

    def inputs(self):
        return {'state': 'parca_state'}

    def outputs(self):
        return {
            'state': 'parca_state',
            'r_vector': 'overwrite',
            'r_columns': 'overwrite',
        }

    def update(self, state):
        t0 = time.time()
        parca_state = state['state']
        sim_data = parca_state.sim_data
        cell_specs = parca_state.cell_specs

        inp = extract_input(sim_data, cell_specs)
        out = compute_promoter_binding(inp)
        merge_output(sim_data, cell_specs, out)

        print(f"  Stage 6 (promoter_binding) completed in {time.time() - t0:.1f}s")
        return {
            'state': ParcaState(sim_data=sim_data, cell_specs=cell_specs),
            'r_vector': out.r_vector,
            'r_columns': out.r_columns,
        }


