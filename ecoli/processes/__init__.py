"""Trimmed `ecoli.processes` package for vParCa.

The original `vivarium-ecoli` `ecoli/processes/__init__.py` registered every
process class in the vivarium process_registry. vParCa only needs a handful of
modules reachable from `ecoli.library.initial_conditions.create_bulk_container`:
`registries`, `partition`, `metabolism`, and `polypeptide_elongation`. None of
the registrations in the original are required for the ParCa, so this is left
empty.
"""
