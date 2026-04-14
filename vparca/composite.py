"""
Composite builder for the ParCa pipeline.

Port-first design: every step declares a port per sim_data leaf it reads
or writes, and the composite wires those ports directly to paths inside a
nested bigraph store that mirrors sim_data's structure.  No ``sim_data``
or ``cell_specs`` blob is carried through the pipeline — the store *is*
the pipeline state.

**This is the PoC composite** covering Step 1 (scatter) and Step 2
(input_adjustments).  Steps 3–9 will be added after they are converted
to the leaf-port style.

Store tree currently wired::

    adjustments/
        balanced_translation_efficiencies
        protein_deg_rates_adjustments
        rna_deg_rates_adjustments
        rna_expression_adjustments
        translation_efficiencies_adjustments
    process/
        transcription/
            cistron_data/id
            cistron_data/deg_rate
            cistron_id_to_rna_indexes_map
            rna_data/id
            rna_data/deg_rate
            rna_expression/basal
        translation/
            monomer_data/id
            monomer_data/deg_rate
            translation_efficiencies_by_monomer
    tf_to_active_inactive_conditions
    raw_data                        (input to Step 1)

Step 1 writes every leaf above (raw_data stays as-is).  Step 2 reads them
all, and writes back the 5 arrays + the TF dict that it mutates.  Dataflow
is a pure DAG: Step 1 → store → Step 2 → store.
"""

from process_bigraph import Composite, allocate_core

from vparca.schema import register_parca_schema
from vparca.steps import ALL_STEP_CLASSES


# ---------------------------------------------------------------------------
# Store paths — single source of truth for port-to-path wiring.
# ---------------------------------------------------------------------------

# Keyed by port name (stable across Step 1's outputs and Step 2's inputs
# where they overlap).  Value is the path list into the store tree.
STORE_PATH = {
    'monomer_ids':                      ['process', 'translation', 'monomer_data', 'id'],
    'translation_efficiencies':         ['process', 'translation', 'translation_efficiencies_by_monomer'],
    'translation_eff_adjustments':      ['adjustments', 'translation_efficiencies_adjustments'],
    'balanced_translation_groups':      ['adjustments', 'balanced_translation_efficiencies'],
    'rna_ids':                          ['process', 'transcription', 'rna_data', 'id'],
    'cistron_ids':                      ['process', 'transcription', 'cistron_data', 'id'],
    'basal_rna_expression':             ['process', 'transcription', 'rna_expression', 'basal'],
    'rna_expression_adjustments':       ['adjustments', 'rna_expression_adjustments'],
    'cistron_id_to_rna_indexes':        ['process', 'transcription', 'cistron_id_to_rna_indexes_map'],
    'rna_deg_rates':                    ['process', 'transcription', 'rna_data', 'deg_rate'],
    'cistron_deg_rates':                ['process', 'transcription', 'cistron_data', 'deg_rate'],
    'rna_deg_rate_adjustments':         ['adjustments', 'rna_deg_rates_adjustments'],
    'protein_deg_rates':                ['process', 'translation', 'monomer_data', 'deg_rate'],
    'protein_deg_rate_adjustments':     ['adjustments', 'protein_deg_rates_adjustments'],
    'tf_to_active_inactive_conditions': ['tf_to_active_inactive_conditions'],
}


def _wires(port_names):
    """Produce a composite ``wires`` dict for the given ports."""
    return {name: STORE_PATH[name] for name in port_names}


def build_parca_composite(raw_data, debug=False, core=None):
    """Build a Composite that runs the PoC ParCa pipeline (steps 1+2).

    Args:
        raw_data: a ``KnowledgeBaseEcoli`` instance.
        debug:    if True, Step 2 reduces tf_to_active_inactive_conditions
                  to a single key.
        core:     optional pre-built core; if omitted one is allocated and
                  schema types + Step classes are registered on it.

    Returns:
        The ``Composite`` instance with the pipeline already executed.
        The final store state is at ``composite.state``.
    """
    if core is None:
        core = allocate_core(top=ALL_STEP_CLASSES)
        register_parca_schema(core)

    # Pull port manifests from the Steps themselves so we can't drift.
    from vparca.steps.step_01_initialize import OUTPUT_PORTS as _step1_out
    from vparca.steps.step_02_input_adjustments import (
        INPUT_PORTS  as _step2_in,
        OUTPUT_PORTS as _step2_out,
    )

    spec = {
        # Fire all Steps once as part of Composite construction — we have
        # no time-advancing processes here, only a DAG of Steps, so this
        # is effectively "run the pipeline".
        'run_steps_on_init': True,
        'state': {
            'initialize': {
                '_type':   'step',
                'address': 'local:InitializeStep',
                # raw_data is carried in config, not a store port, so
                # bigraph-schema doesn't introspect the KB's internals.
                'config':  {'raw_data': raw_data},
                'inputs':  {},
                'outputs': _wires(_step1_out.keys()),
            },

            'input_adjustments': {
                '_type':   'step',
                'address': 'local:InputAdjustmentsStep',
                'config':  {'debug': debug},
                'inputs':  _wires(_step2_in.keys()),
                'outputs': _wires(_step2_out.keys()),
            },
        },
    }

    return Composite(spec, core=core)


def run_parca(raw_data, debug=False):
    """Build the PoC composite, let the Step DAG execute, return store state."""
    composite = build_parca_composite(raw_data, debug=debug)
    return composite.state
