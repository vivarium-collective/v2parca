"""
Helper for Steps 3–9: build a ``sim_data``-shaped object from port values
so the existing ParCa sub-functions (``expressionConverge``,
``setKmCooperative…``, ``fitPromoterBoundProbability``, …) — which were
written to take a whole ``SimulationDataEcoli`` — keep working unchanged.

Two-layer design:

1.  A ``SimpleNamespace`` we populate from the port values wired into
    the Step's inputs.  This is the primary, explicit surface: every
    subsystem and top-level dict the Step declares gets installed at
    its canonical ``sim_data`` attribute path.

2.  An **attribute-fallback proxy** that routes misses to the live
    ``SimulationDataEcoli`` instance handed in via the ``sim_data_root``
    port.  This covers the handful of rare paths a Step's callees reach
    that aren't worth hand-wiring — e.g. ``sim_data.calculate_ppgpp_expression``
    (a method defined on sim_data itself, not a subsystem) or
    ``sim_data.genetic_perturbations`` (a top-level dict used only by
    ``create_bulk_container``).  Because the root holds references to
    the same subsystem objects that are in the ports, mutations made
    through either route stay coherent.
"""

from types import SimpleNamespace


class _FacadeProxy:
    """Attribute-lookup that first consults the ``SimpleNamespace`` of
    wired ports, then falls back to the root SimulationDataEcoli."""

    __slots__ = ('_ns', '_root')

    def __init__(self, ns, root):
        object.__setattr__(self, '_ns', ns)
        object.__setattr__(self, '_root', root)

    def __getattr__(self, name):
        # Ports take priority.
        if hasattr(self._ns, name):
            return getattr(self._ns, name)
        if self._root is not None:
            return getattr(self._root, name)
        raise AttributeError(name)

    def __setattr__(self, name, value):
        # Set on the namespace *and* the root so any sub-function that
        # mutates a top-level attr (e.g. fitPromoterBoundProbability does
        # ``sim_data.pPromoterBound = …``) is visible on both surfaces.
        setattr(self._ns, name, value)
        if self._root is not None:
            setattr(self._root, name, value)


def make_sim_data_facade(ports):
    """Build a sim_data-shaped object from a port-value dict.

    Any port listed below that is present in ``ports`` is installed at
    its canonical sim_data path on a namespace.  Unwired attribute
    accesses fall through to ``ports['sim_data_root']`` when present.

    Args:
        ports: dict mapping port names (from any Step's ``inputs()`` keys)
            to their live values.
    Returns:
        a proxy that quacks like ``SimulationDataEcoli`` across the union
        of the wired ports and the root.
    """
    proc = SimpleNamespace()
    internal = SimpleNamespace()
    ns = SimpleNamespace(process=proc, internal_state=internal)

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

    top_level = [
        'mass', 'constants', 'growth_rate_parameters',
        'adjustments', 'molecule_groups', 'molecule_ids', 'relation',
        'getter',
    ]
    for name in top_level:
        if name in ports:
            setattr(ns, name, ports[name])

    if 'bulk_molecules' in ports:
        internal.bulk_molecules = ports['bulk_molecules']

    if 'external_state' in ports:
        ns.external_state = ports['external_state']

    data_leaf_attrs = [
        'tf_to_active_inactive_conditions', 'conditions',
        'condition_to_doubling_time', 'tf_to_fold_change',
        'tf_to_direction',
        'condition_active_tfs', 'condition_inactive_tfs',
        'translation_supply_rate',
        'expected_dry_mass_increase_dict',
        'pPromoterBound',
        'condition',
    ]
    for name in data_leaf_attrs:
        if name in ports:
            setattr(ns, name, ports[name])

    if 'expected_dry_mass_increase_dict' in ports:
        ns.expectedDryMassIncreaseDict = ports['expected_dry_mass_increase_dict']

    return _FacadeProxy(ns, ports.get('sim_data_root'))
