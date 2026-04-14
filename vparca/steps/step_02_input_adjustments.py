"""
Step 2 — input_adjustments.  Apply pre-fitted adjustments to translation
efficiencies, RNA expression, and degradation rates.

Port-first design: the Step declares each sim_data leaf it reads or writes
as an explicit port.  The composite wires each port to a path inside a
nested bigraph store that mirrors sim_data's structure — no ``sim_data`` or
``cell_specs`` object is passed between Steps.

Store paths wired by the composite (see ``vparca/composite.py``):

    READS
      monomer_ids                        process / translation / monomer_data / id
      translation_efficiencies           process / translation / translation_efficiencies_by_monomer
      translation_eff_adjustments        adjustments / translation_efficiencies_adjustments
      balanced_translation_groups        adjustments / balanced_translation_efficiencies
      rna_ids                            process / transcription / rna_data / id
      cistron_ids                        process / transcription / cistron_data / id
      basal_rna_expression               process / transcription / rna_expression / basal
      rna_expression_adjustments         adjustments / rna_expression_adjustments
      cistron_id_to_rna_indexes          process / transcription / cistron_id_to_rna_indexes_map
      rna_deg_rates                      process / transcription / rna_data / deg_rate
      cistron_deg_rates                  process / transcription / cistron_data / deg_rate
      rna_deg_rate_adjustments           adjustments / rna_deg_rates_adjustments
      protein_deg_rates                  process / translation / monomer_data / deg_rate
      protein_deg_rate_adjustments       adjustments / protein_deg_rates_adjustments
      tf_to_active_inactive_conditions   tf_to_active_inactive_conditions

    WRITES (same paths, round-trip for the 5 mutated arrays + optional TF dict)
      translation_efficiencies           process / translation / translation_efficiencies_by_monomer
      basal_rna_expression               process / transcription / rna_expression / basal
      rna_deg_rates                      process / transcription / rna_data / deg_rate
      cistron_deg_rates                  process / transcription / cistron_data / deg_rate
      protein_deg_rates                  process / translation / monomer_data / deg_rate
      tf_to_active_inactive_conditions   tf_to_active_inactive_conditions (conditional)

All sub-functions below are **pure**: numpy-only, no sim_data, no side
effects.  The Step simply assembles its inputs, calls them, and emits its
outputs by port name.
"""

import time

import numpy as np

from process_bigraph import Step


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


# ============================================================================
# Step class — one port per leaf
# ============================================================================

# Port schema — every entry is an explicit sim_data (or cell_specs) leaf
# the Step reads or writes.  Names here are the port identifiers used by
# the Composite to wire into store paths.  Value ``'overwrite'`` is the
# bigraph-schema type used for replace-semantics on every leaf.
INPUT_PORTS = {
    'monomer_ids':                      'overwrite',
    'translation_efficiencies':         'overwrite',
    'translation_eff_adjustments':      'overwrite',
    'balanced_translation_groups':      'overwrite',
    'rna_ids':                          'overwrite',
    'cistron_ids':                      'overwrite',
    'basal_rna_expression':             'overwrite',
    'rna_expression_adjustments':       'overwrite',
    'cistron_id_to_rna_indexes':        'overwrite',
    'rna_deg_rates':                    'overwrite',
    'cistron_deg_rates':                'overwrite',
    'rna_deg_rate_adjustments':         'overwrite',
    'protein_deg_rates':                'overwrite',
    'protein_deg_rate_adjustments':     'overwrite',
    'tf_to_active_inactive_conditions': 'overwrite',
}

OUTPUT_PORTS = {
    'translation_efficiencies':         'overwrite',
    'basal_rna_expression':             'overwrite',
    'rna_deg_rates':                    'overwrite',
    'cistron_deg_rates':                'overwrite',
    'protein_deg_rates':                'overwrite',
    'tf_to_active_inactive_conditions': 'overwrite',
}


class InputAdjustmentsStep(Step):
    """Step 2 — input_adjustments.

    Declares a port for every sim_data leaf read or written.  The composite
    wires each port directly to the corresponding store path (see this
    module's docstring for the read/write path table).  No ``sim_data`` or
    ``cell_specs`` blob is passed through any port.
    """

    config_schema = {
        'debug': {'_type': 'boolean', '_default': False},
    }

    def inputs(self):
        return dict(INPUT_PORTS)

    def outputs(self):
        return dict(OUTPUT_PORTS)

    def update(self, state):
        t0 = time.time()

        # ---- debug switch: optionally trim TF conditions -------------------
        tf_conditions_out = None
        if self.config.get('debug', False):
            print(
                "  Step 2: debug mode — reducing tf_to_active_inactive_conditions"
                " to a single key"
            )
            tf_cond = state['tf_to_active_inactive_conditions']
            first_key = next(iter(tf_cond))
            tf_conditions_out = {first_key: tf_cond[first_key]}

        # ---- translation efficiencies -------------------------------------
        eff = adjust_translation_efficiencies(
            state['monomer_ids'],
            # defensively copy the live array — the store holds the live ref
            np.asarray(state['translation_efficiencies']).copy(),
            state['translation_eff_adjustments'],
        )
        eff = balance_translation_efficiencies(
            state['monomer_ids'],
            eff,
            state['balanced_translation_groups'],
        )

        # ---- RNA expression -----------------------------------------------
        expr = adjust_rna_expression(
            state['rna_ids'],
            state['cistron_ids'],
            np.asarray(state['basal_rna_expression']).copy(),
            state['rna_expression_adjustments'],
            state['cistron_id_to_rna_indexes'],
        )

        # ---- RNA + cistron degradation rates ------------------------------
        rna_deg, cistron_deg = adjust_rna_deg_rates(
            state['rna_ids'],
            state['cistron_ids'],
            np.asarray(state['rna_deg_rates']).copy(),
            np.asarray(state['cistron_deg_rates']).copy(),
            state['rna_deg_rate_adjustments'],
            state['cistron_id_to_rna_indexes'],
        )

        # ---- protein degradation rates ------------------------------------
        prot_deg = adjust_protein_deg_rates(
            state['monomer_ids'],
            np.asarray(state['protein_deg_rates']).copy(),
            state['protein_deg_rate_adjustments'],
        )

        print(f"  Step 2 (input_adjustments) completed in {time.time() - t0:.1f}s")

        out = {
            'translation_efficiencies': eff,
            'basal_rna_expression':     expr,
            'rna_deg_rates':            rna_deg,
            'cistron_deg_rates':        cistron_deg,
            'protein_deg_rates':        prot_deg,
        }
        if tf_conditions_out is not None:
            out['tf_to_active_inactive_conditions'] = tf_conditions_out
        return out


