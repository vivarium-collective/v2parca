"""
vParCa — process-bigraph migration of the E. coli Parameter Calculator.

Layout
======

    vparca/
        composite.py          build_parca_composite / run_parca / register_parca_steps
        state.py              ParcaState bigraph-schema type + register_parca_types
        types.py              Input/Output dataclasses for every stage
        fitting.py            pure math + sim_data-reading fitting helpers
        promoter_fitting.py   matrix builders + CVXPY optimization (stages 6/7)
        stages/
            __init__.py       ALL_STEP_CLASSES registry
            stage_01_initialize.py
            stage_02_input_adjustments.py      (Extract → Compute → Merge, PURE)
            stage_03_basal_specs.py            (COUPLED)
            stage_04_tf_condition_specs.py     (COUPLED)
            stage_05_fit_condition.py          (READ-ONLY)
            stage_06_promoter_binding.py       (COUPLED)
            stage_07_adjust_promoters.py       (COUPLED)
            stage_08_set_conditions.py         (Extract → Compute → Merge, PURE)
            stage_09_final_adjustments.py      (COUPLED)

Each stage module (2–9) provides:

    extract_input(sim_data, cell_specs, **kwargs) -> StageInput
    compute_<name>(inp: StageInput)              -> StageOutput
    merge_output(sim_data, cell_specs, out)      -> None
    <StageName>Step(Step)                         — process-bigraph wrapper

Pure stages (2, 8) additionally expose ExtractFor... and MergeAfter...
Step classes so the compute Step has only explicit typed ports.

Purity legend
=============
  PURE       compute has no sim_data/cell_specs access
  READ-ONLY  compute reads sim_data via ref but does not mutate it
  COUPLED    compute still mutates sim_data via ref (future refactor target)
"""
