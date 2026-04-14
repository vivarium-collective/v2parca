"""
Helper for Steps 3–9: build a ``sim_data``-shaped ``SimpleNamespace``
facade from port values so the existing sub-functions
(``expressionConverge``, ``setKmCooperative…``, ``fitPromoterBoundProbability``,
…) — which were written to take a whole ``SimulationDataEcoli`` — keep
working unchanged.

The facade holds *live references* to the subsystem objects handed in
through the Step's input ports.  Mutations inside sub-functions go
through those references, so when the Step returns the subsystem objects
on its output ports, the store receives the mutated state.
"""

from types import SimpleNamespace


def make_sim_data_facade(ports):
    """Build a minimal sim_data-shaped namespace from a port-value dict.

    Any port listed below that is present in ``ports`` is installed at
    its canonical sim_data path.  Missing ports are skipped — each Step
    only wires in what it actually uses.

    Args:
        ports: dict mapping port names (from any Step's inputs() keys) to
            their live values.
    Returns:
        a ``SimpleNamespace`` that quacks like ``SimulationDataEcoli``
        for the attributes needed by ParCa sub-functions.
    """
    proc = SimpleNamespace()
    internal = SimpleNamespace()
    sd = SimpleNamespace(process=proc, internal_state=internal)

    # Subsystem objects — installed onto sim_data.process.*
    process_attrs = {
        'transcription':            'transcription',
        'translation':              'translation',
        'metabolism':               'metabolism',
        'rna_decay':                'rna_decay',
        'complexation':             'complexation',
        'equilibrium':              'equilibrium',
        'two_component_system':     'two_component_system',
        'transcription_regulation': 'transcription_regulation',
        'replication':              'replication',
    }
    for port_name, attr_name in process_attrs.items():
        if port_name in ports:
            setattr(proc, attr_name, ports[port_name])

    # Top-level subsystems on sim_data directly.
    top_level = [
        'mass', 'constants', 'growth_rate_parameters',
        'adjustments', 'molecule_groups', 'molecule_ids', 'relation',
        'getter',
    ]
    for name in top_level:
        if name in ports:
            setattr(sd, name, ports[name])

    # internal_state.bulk_molecules
    if 'bulk_molecules' in ports:
        internal.bulk_molecules = ports['bulk_molecules']

    # Pure-data top-level dicts.
    data_leaf_attrs = [
        'tf_to_active_inactive_conditions', 'conditions',
        'condition_to_doubling_time', 'tf_to_fold_change',
        'condition_active_tfs', 'condition_inactive_tfs',
    ]
    for name in data_leaf_attrs:
        if name in ports:
            setattr(sd, name, ports[name])

    return sd
