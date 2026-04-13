"""
Stage 2: input_adjustments — Apply pre-fitted adjustments to translation
efficiencies, RNA expression, and degradation rates.

Three public functions:
    extract_input(sim_data, cell_specs, **kwargs) -> InputAdjustmentsInput
    compute_input_adjustments(inp) -> InputAdjustmentsOutput
    merge_output(sim_data, cell_specs, out)

Pure sub-functions (no sim_data dependency):
    adjust_translation_efficiencies
    balance_translation_efficiencies
    adjust_rna_expression
    adjust_rna_deg_rates
    adjust_protein_deg_rates
"""

import numpy as np

from vparca.types import (
    InputAdjustmentsInput,
    InputAdjustmentsOutput,
)


import time
from process_bigraph import Step

from vparca.state import ParcaState

# ============================================================================
# Extract / Merge
# ============================================================================


def extract_input(sim_data, cell_specs, **kwargs) -> InputAdjustmentsInput:
    """Pull specific data from sim_data into a typed dataclass."""
    transcription = sim_data.process.transcription
    translation = sim_data.process.translation

    # Pre-compute the cistron-to-RNA mapping as a plain dict
    cistron_ids = transcription.cistron_data["id"]
    cistron_id_to_rna_indexes = {}
    for cid in cistron_ids:
        cistron_id_to_rna_indexes[cid] = transcription.cistron_id_to_rna_indexes(cid)

    return InputAdjustmentsInput(
        debug=kwargs.get("debug", False),
        # Translation efficiencies
        monomer_ids=translation.monomer_data["id"],
        translation_efficiencies=translation.translation_efficiencies_by_monomer.copy(),
        translation_eff_adjustments=dict(
            sim_data.adjustments.translation_efficiencies_adjustments
        ),
        balanced_translation_groups=list(
            sim_data.adjustments.balanced_translation_efficiencies
        ),
        # RNA expression
        rna_ids=transcription.rna_data["id"],
        cistron_ids=cistron_ids,
        basal_rna_expression=transcription.rna_expression["basal"].copy(),
        rna_expression_adjustments=dict(
            sim_data.adjustments.rna_expression_adjustments
        ),
        cistron_id_to_rna_indexes=cistron_id_to_rna_indexes,
        # Degradation rates
        rna_deg_rates=transcription.rna_data.struct_array["deg_rate"].copy(),
        cistron_deg_rates=transcription.cistron_data.struct_array["deg_rate"].copy(),
        rna_deg_rate_adjustments=dict(
            sim_data.adjustments.rna_deg_rates_adjustments
        ),
        protein_deg_rates=translation.monomer_data.struct_array["deg_rate"].copy(),
        protein_deg_rate_adjustments=dict(
            sim_data.adjustments.protein_deg_rates_adjustments
        ),
        # TF conditions
        tf_to_active_inactive_conditions=dict(
            sim_data.tf_to_active_inactive_conditions
        ),
    )


def merge_output(sim_data, cell_specs, out: InputAdjustmentsOutput):
    """Write computed results back into sim_data."""
    transcription = sim_data.process.transcription
    translation = sim_data.process.translation

    translation.translation_efficiencies_by_monomer[:] = out.translation_efficiencies
    transcription.rna_expression["basal"][:] = out.basal_rna_expression
    transcription.rna_data.struct_array["deg_rate"][:] = out.rna_deg_rates
    transcription.cistron_data.struct_array["deg_rate"][:] = out.cistron_deg_rates
    translation.monomer_data.struct_array["deg_rate"][:] = out.protein_deg_rates

    if out.tf_to_active_inactive_conditions is not None:
        sim_data.tf_to_active_inactive_conditions = (
            out.tf_to_active_inactive_conditions
        )


# ============================================================================
# Pure sub-functions
# ============================================================================


def adjust_translation_efficiencies(monomer_ids, efficiencies, adjustments):
    """
    Multiply translation efficiencies by specified adjustment factors.

    Args:
        monomer_ids: array of monomer ID strings (with "[c]" suffix)
        efficiencies: array of translation efficiencies (modified in-place copy)
        adjustments: {protein_id: multiplier}

    Returns:
        Modified efficiencies array.
    """
    result = efficiencies.copy()
    for protein, multiplier in adjustments.items():
        idx = np.where(monomer_ids == protein)[0]
        result[idx] *= multiplier
    return result


def balance_translation_efficiencies(monomer_ids, efficiencies, groups):
    """
    Set translation efficiencies within each group to the group mean.

    Args:
        monomer_ids: array of monomer ID strings (with "[c]" suffix)
        efficiencies: array of translation efficiencies
        groups: list of lists of protein IDs (without "[c]" suffix)

    Returns:
        Modified efficiencies array.
    """
    result = efficiencies.copy()
    monomer_id_to_index = {
        mid[:-3]: i for i, mid in enumerate(monomer_ids)
    }
    for proteins in groups:
        protein_indexes = np.array([monomer_id_to_index[m] for m in proteins])
        mean_eff = result[protein_indexes].mean()
        result[protein_indexes] = mean_eff
    return result


def adjust_rna_expression(
    rna_ids, cistron_ids, expression, adjustments, cistron_id_to_rna_indexes
):
    """
    Adjust basal RNA expression levels by specified factors, then normalize.

    If a mol_id is a cistron, all RNAs containing that cistron are adjusted.
    If multiple adjustments affect the same RNA, the maximum factor is used.

    Args:
        rna_ids: array of RNA ID strings
        cistron_ids: array of cistron ID strings
        expression: array of basal expression values
        adjustments: {mol_id: multiplier}
        cistron_id_to_rna_indexes: {cistron_id: array of RNA indexes}

    Returns:
        Normalized expression array.
    """
    result = expression.copy()
    cistron_id_set = set(cistron_ids)
    rna_id_to_index = {rna_id[:-3]: i for i, rna_id in enumerate(rna_ids)}

    rna_index_to_adjustment = {}

    for mol_id, adj_factor in adjustments.items():
        if mol_id in cistron_id_set:
            rna_indexes = cistron_id_to_rna_indexes[mol_id]
        elif mol_id in rna_id_to_index:
            rna_indexes = rna_id_to_index[mol_id]
        else:
            raise ValueError(
                f"Molecule ID {mol_id} not found in list of cistrons or"
                " transcription units."
            )

        # If multiple adjustments hit the same RNA, take the maximum
        for rna_index in np.atleast_1d(rna_indexes):
            rna_index_to_adjustment[rna_index] = max(
                rna_index_to_adjustment.get(rna_index, 0), adj_factor
            )

    for rna_index, adj_factor in rna_index_to_adjustment.items():
        result[rna_index] *= adj_factor

    result /= result.sum()
    return result


def adjust_rna_deg_rates(
    rna_ids, cistron_ids, rna_rates, cistron_rates, adjustments,
    cistron_id_to_rna_indexes
):
    """
    Adjust RNA and cistron degradation rates by specified factors.

    If a mol_id is a cistron, both the cistron rate and the rates of all
    RNAs containing that cistron are adjusted. If multiple adjustments hit
    the same RNA, the maximum factor is used.

    Args:
        rna_ids: array of RNA ID strings
        cistron_ids: array of cistron ID strings
        rna_rates: array of RNA degradation rates
        cistron_rates: array of cistron degradation rates
        adjustments: {mol_id: multiplier}
        cistron_id_to_rna_indexes: {cistron_id: array of RNA indexes}

    Returns:
        (rna_rates, cistron_rates) — both modified copies.
    """
    rna_result = rna_rates.copy()
    cistron_result = cistron_rates.copy()

    cistron_id_to_index = {cid: i for i, cid in enumerate(cistron_ids)}
    rna_id_to_index = {rna_id[:-3]: i for i, rna_id in enumerate(rna_ids)}

    rna_index_to_adjustment = {}

    for mol_id, adj_factor in adjustments.items():
        if mol_id in cistron_id_to_index:
            # Adjust the cistron degradation rate
            cistron_index = cistron_id_to_index[mol_id]
            cistron_result[cistron_index] *= adj_factor

            # Find all RNAs containing this cistron
            rna_indexes = cistron_id_to_rna_indexes[mol_id]
        elif mol_id in rna_id_to_index:
            rna_indexes = rna_id_to_index[mol_id]
        else:
            raise ValueError(
                f"Molecule ID {mol_id} not found in list of cistrons or"
                " transcription units."
            )

        for rna_index in np.atleast_1d(rna_indexes):
            rna_index_to_adjustment[rna_index] = max(
                rna_index_to_adjustment.get(rna_index, 0), adj_factor
            )

    for rna_index, adj_factor in rna_index_to_adjustment.items():
        rna_result[rna_index] *= adj_factor

    return rna_result, cistron_result


def adjust_protein_deg_rates(monomer_ids, rates, adjustments):
    """
    Multiply protein degradation rates by specified adjustment factors.

    Args:
        monomer_ids: array of monomer ID strings
        rates: array of protein degradation rates
        adjustments: {protein_id: multiplier}

    Returns:
        Modified rates array.
    """
    result = rates.copy()
    for protein, multiplier in adjustments.items():
        idx = np.where(monomer_ids == protein)[0]
        result[idx] *= multiplier
    return result


# ============================================================================
# Main compute function
# ============================================================================


def compute_input_adjustments(inp: InputAdjustmentsInput) -> InputAdjustmentsOutput:
    """
    Pure function: apply all input adjustments.

    No sim_data, no cell_specs, no side effects.
    """
    # Debug mode: limit TF conditions
    tf_conditions = None
    if inp.debug:
        print(
            "Warning: Running the Parca in debug mode"
            " - not all conditions will be fit"
        )
        key = list(inp.tf_to_active_inactive_conditions.keys())[0]
        tf_conditions = {key: inp.tf_to_active_inactive_conditions[key]}

    # Translation efficiencies
    eff = adjust_translation_efficiencies(
        inp.monomer_ids, inp.translation_efficiencies, inp.translation_eff_adjustments
    )
    eff = balance_translation_efficiencies(
        inp.monomer_ids, eff, inp.balanced_translation_groups
    )

    # RNA expression
    expr = adjust_rna_expression(
        inp.rna_ids,
        inp.cistron_ids,
        inp.basal_rna_expression,
        inp.rna_expression_adjustments,
        inp.cistron_id_to_rna_indexes,
    )

    # RNA degradation rates
    rna_deg, cistron_deg = adjust_rna_deg_rates(
        inp.rna_ids,
        inp.cistron_ids,
        inp.rna_deg_rates,
        inp.cistron_deg_rates,
        inp.rna_deg_rate_adjustments,
        inp.cistron_id_to_rna_indexes,
    )

    # Protein degradation rates
    prot_deg = adjust_protein_deg_rates(
        inp.monomer_ids, inp.protein_deg_rates, inp.protein_deg_rate_adjustments
    )

    return InputAdjustmentsOutput(
        translation_efficiencies=eff,
        basal_rna_expression=expr,
        rna_deg_rates=rna_deg,
        cistron_deg_rates=cistron_deg,
        protein_deg_rates=prot_deg,
        tf_to_active_inactive_conditions=tf_conditions,
    )



# ---------------------------------------------------------------------------
# Stage 2: Input Adjustments — Extract / Pure / Merge
# ---------------------------------------------------------------------------

# Port schema shared by Extract (outputs) and InputAdjustmentsStep (inputs)
_STAGE_02_INPUT_PORTS = {
    'monomer_ids': 'overwrite',
    'translation_efficiencies': 'overwrite',
    'translation_eff_adjustments': 'overwrite',
    'balanced_translation_groups': 'overwrite',
    'rna_ids': 'overwrite',
    'cistron_ids': 'overwrite',
    'basal_rna_expression': 'overwrite',
    'rna_expression_adjustments': 'overwrite',
    'cistron_id_to_rna_indexes': 'overwrite',
    'rna_deg_rates': 'overwrite',
    'cistron_deg_rates': 'overwrite',
    'rna_deg_rate_adjustments': 'overwrite',
    'protein_deg_rates': 'overwrite',
    'protein_deg_rate_adjustments': 'overwrite',
    'tf_to_active_inactive_conditions': 'overwrite',
}

# Port schema shared by InputAdjustmentsStep (outputs) and Merge (inputs)
_STAGE_02_OUTPUT_PORTS = {
    'translation_efficiencies': 'overwrite',
    'basal_rna_expression': 'overwrite',
    'rna_deg_rates': 'overwrite',
    'cistron_deg_rates': 'overwrite',
    'protein_deg_rates': 'overwrite',
    'tf_to_active_inactive_conditions': 'overwrite',
}


class ExtractForStage2Step(Step):
    """Helper: extract fields from parca_state for the pure Stage 2 step."""

    config_schema = {
        'debug': {'_type': 'boolean', '_default': False},
    }

    def inputs(self):
        return {'state': 'parca_state'}

    def outputs(self):
        return dict(_STAGE_02_INPUT_PORTS)

    def update(self, state):
        parca_state = state['state']
        inp = extract_input(parca_state.sim_data, parca_state.cell_specs,
                          debug=self.config.get('debug', False))
        return {
            'monomer_ids': inp.monomer_ids,
            'translation_efficiencies': inp.translation_efficiencies,
            'translation_eff_adjustments': inp.translation_eff_adjustments,
            'balanced_translation_groups': inp.balanced_translation_groups,
            'rna_ids': inp.rna_ids,
            'cistron_ids': inp.cistron_ids,
            'basal_rna_expression': inp.basal_rna_expression,
            'rna_expression_adjustments': inp.rna_expression_adjustments,
            'cistron_id_to_rna_indexes': inp.cistron_id_to_rna_indexes,
            'rna_deg_rates': inp.rna_deg_rates,
            'cistron_deg_rates': inp.cistron_deg_rates,
            'rna_deg_rate_adjustments': inp.rna_deg_rate_adjustments,
            'protein_deg_rates': inp.protein_deg_rates,
            'protein_deg_rate_adjustments': inp.protein_deg_rate_adjustments,
            'tf_to_active_inactive_conditions': inp.tf_to_active_inactive_conditions,
        }


class InputAdjustmentsStep(Step):
    """Stage 2: PURE — every field is an explicit typed port.

    No parca_state in inputs or outputs.  Operates entirely on
    individually-named data ports extracted by ExtractForStage2Step.
    """

    config_schema = {
        'debug': {'_type': 'boolean', '_default': False},
    }

    def inputs(self):
        return dict(_STAGE_02_INPUT_PORTS)

    def outputs(self):
        return dict(_STAGE_02_OUTPUT_PORTS)

    def update(self, state):
        t0 = time.time()
        inp = InputAdjustmentsInput(
            debug=self.config.get('debug', False),
            monomer_ids=state['monomer_ids'],
            translation_efficiencies=state['translation_efficiencies'],
            translation_eff_adjustments=state['translation_eff_adjustments'],
            balanced_translation_groups=state['balanced_translation_groups'],
            rna_ids=state['rna_ids'],
            cistron_ids=state['cistron_ids'],
            basal_rna_expression=state['basal_rna_expression'],
            rna_expression_adjustments=state['rna_expression_adjustments'],
            cistron_id_to_rna_indexes=state['cistron_id_to_rna_indexes'],
            rna_deg_rates=state['rna_deg_rates'],
            cistron_deg_rates=state['cistron_deg_rates'],
            rna_deg_rate_adjustments=state['rna_deg_rate_adjustments'],
            protein_deg_rates=state['protein_deg_rates'],
            protein_deg_rate_adjustments=state['protein_deg_rate_adjustments'],
            tf_to_active_inactive_conditions=state['tf_to_active_inactive_conditions'],
        )
        out = compute_input_adjustments(inp)
        print(f"  Stage 2 (input_adjustments) completed in {time.time() - t0:.1f}s")
        return {
            'translation_efficiencies': out.translation_efficiencies,
            'basal_rna_expression': out.basal_rna_expression,
            'rna_deg_rates': out.rna_deg_rates,
            'cistron_deg_rates': out.cistron_deg_rates,
            'protein_deg_rates': out.protein_deg_rates,
            'tf_to_active_inactive_conditions': out.tf_to_active_inactive_conditions,
        }


class MergeAfterStage2Step(Step):
    """Helper: merge Stage 2 outputs back into parca_state."""

    def inputs(self):
        return {
            'state': 'parca_state',
            **_STAGE_02_OUTPUT_PORTS,
        }

    def outputs(self):
        return {'state': 'parca_state'}

    def update(self, state):
        parca_state = state['state']
        sim_data = parca_state.sim_data
        cell_specs = parca_state.cell_specs

        out = InputAdjustmentsOutput(
            translation_efficiencies=state['translation_efficiencies'],
            basal_rna_expression=state['basal_rna_expression'],
            rna_deg_rates=state['rna_deg_rates'],
            cistron_deg_rates=state['cistron_deg_rates'],
            protein_deg_rates=state['protein_deg_rates'],
            tf_to_active_inactive_conditions=state['tf_to_active_inactive_conditions'],
        )
        merge_output(sim_data, cell_specs, out)

        return {
            'state': ParcaState(sim_data=sim_data, cell_specs=cell_specs),
        }


