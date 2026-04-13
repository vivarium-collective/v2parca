"""Custom bigraph-schema type for the ParCa pipeline state.

ParcaState is a value object that wraps sim_data + cell_specs, replacing
opaque 'overwrite'-typed ports named 'sim_data'/'cell_specs' with a single
registered type that carries both.
"""

from dataclasses import dataclass

from bigraph_schema.schema import Overwrite


@dataclass(kw_only=True)
class ParcaStateSchema(Overwrite):
    """Bigraph-schema type for the ParCa pipeline context.

    Behaves identically to 'overwrite' (wholesale replacement on update)
    but provides a distinct type name for documentation and introspection.
    """
    pass


@dataclass
class ParcaState:
    """Value object wrapping sim_data + cell_specs for the ParCa pipeline."""
    sim_data: object = None
    cell_specs: dict = None


def register_parca_types(core):
    """Register the 'parca_state' type with the given bigraph-schema core."""
    core.register_type('parca_state', ParcaStateSchema)
    return core
