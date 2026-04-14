"""
Composite builder for the ParCa pipeline.

Port-first, nested-store design.  The composite's state is a bigraph tree
that mirrors ``SimulationDataEcoli``'s own structure: subsystem objects
(``Transcription``, ``Mass``, ``Constants``, …) sit at natural paths like
``process/transcription`` or ``constants``, and pure-data top-level dicts
(``tf_to_active_inactive_conditions``, ``conditions``, …) sit at sibling
paths.  Each Step's ports wire to the subsystems and leaves it actually
touches — nothing travels through a monolithic ``sim_data`` blob.

**This PoC composite** runs Step 1 (scatter) and Step 2 (input_adjustments).
Steps 3–9 will be wired in after they're converted to the object/leaf
port style.
"""

from process_bigraph import Composite, allocate_core

from vparca.schema import register_parca_schema
from vparca.steps import ALL_STEP_CLASSES


# ---------------------------------------------------------------------------
# Port-to-store-path table — single source of truth for every wire in the
# composite.  Port names are kept globally unique across Steps so this one
# table suffices for all wiring.
# ---------------------------------------------------------------------------

STORE_PATH = {
    # subsystem object leaves
    'transcription':                    ['process', 'transcription'],
    'translation':                      ['process', 'translation'],
    'metabolism':                       ['process', 'metabolism'],
    'rna_decay':                        ['process', 'rna_decay'],
    'complexation':                     ['process', 'complexation'],
    'equilibrium':                      ['process', 'equilibrium'],
    'two_component_system':             ['process', 'two_component_system'],
    'transcription_regulation':         ['process', 'transcription_regulation'],
    'replication':                      ['process', 'replication'],
    'mass':                             ['mass'],
    'constants':                        ['constants'],
    'growth_rate_parameters':           ['growth_rate_parameters'],
    'adjustments':                      ['adjustments'],
    'molecule_groups':                  ['molecule_groups'],
    'bulk_molecules':                   ['internal_state', 'bulk_molecules'],

    # pure-data top-level dicts
    'tf_to_active_inactive_conditions': ['tf_to_active_inactive_conditions'],
    'conditions':                       ['conditions'],
    'condition_to_doubling_time':       ['condition_to_doubling_time'],
    'tf_to_fold_change':                ['tf_to_fold_change'],
    'condition_active_tfs':             ['condition_active_tfs'],
    'condition_inactive_tfs':           ['condition_inactive_tfs'],
}


def _wires(port_names):
    """Produce a composite ``wires`` dict for the given ports."""
    return {name: STORE_PATH[name] for name in port_names}


def build_parca_composite(raw_data, debug=False, core=None):
    """Build a Composite that runs the PoC ParCa pipeline (steps 1+2).

    Args:
        raw_data: a ``KnowledgeBaseEcoli`` instance.  Passed through
            InitializeStep's config to keep bigraph-schema from walking
            its nested KB internals at composite construction time.
        debug:    if True, Step 2 reduces tf_to_active_inactive_conditions
                  to a single key.
        core:     optional pre-built core; if omitted one is allocated
                  and schema types + Step classes are registered on it.
    Returns:
        The ``Composite`` instance with the pipeline already executed.
        The final store state is at ``composite.state``.
    """
    if core is None:
        core = allocate_core(top=ALL_STEP_CLASSES)
        register_parca_schema(core)

    from vparca.steps.step_01_initialize import OUTPUT_PORTS as _step1_out
    from vparca.steps.step_02_input_adjustments import (
        INPUT_PORTS  as _step2_in,
        OUTPUT_PORTS as _step2_out,
    )

    spec = {
        'run_steps_on_init': True,
        'state': {
            'initialize': {
                '_type':   'step',
                'address': 'local:InitializeStep',
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
