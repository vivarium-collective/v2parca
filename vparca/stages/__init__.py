"""
Process-bigraph Step classes for the ParCa pipeline, one module per stage.

Each stage module exposes its Step class(es) alongside the
``extract_input`` / ``compute_*`` / ``merge_output`` functions that
do the actual work.  Pure stages (2, 8) are decomposed into an
``Extract → Compute → Merge`` triplet of Steps so the compute Step
has only explicit typed ports.

``ALL_STEP_CLASSES`` is the flat registry used by
``process_bigraph.allocate_core(top=...)``.
"""

from vparca.stages.stage_01_initialize import InitializeStep
from vparca.stages.stage_02_input_adjustments import (
    ExtractForStage2Step,
    InputAdjustmentsStep,
    MergeAfterStage2Step,
)
from vparca.stages.stage_03_basal_specs import BasalSpecsStep
from vparca.stages.stage_04_tf_condition_specs import TfConditionSpecsStep
from vparca.stages.stage_05_fit_condition import FitConditionStep
from vparca.stages.stage_06_promoter_binding import PromoterBindingStep
from vparca.stages.stage_07_adjust_promoters import AdjustPromotersStep
from vparca.stages.stage_08_set_conditions import (
    ExtractForStage8Step,
    SetConditionsStep,
    MergeAfterStage8Step,
)
from vparca.stages.stage_09_final_adjustments import FinalAdjustmentsStep


ALL_STEP_CLASSES = {
    'InitializeStep': InitializeStep,
    'ExtractForStage2Step': ExtractForStage2Step,
    'InputAdjustmentsStep': InputAdjustmentsStep,
    'MergeAfterStage2Step': MergeAfterStage2Step,
    'BasalSpecsStep': BasalSpecsStep,
    'TfConditionSpecsStep': TfConditionSpecsStep,
    'FitConditionStep': FitConditionStep,
    'PromoterBindingStep': PromoterBindingStep,
    'AdjustPromotersStep': AdjustPromotersStep,
    'ExtractForStage8Step': ExtractForStage8Step,
    'SetConditionsStep': SetConditionsStep,
    'MergeAfterStage8Step': MergeAfterStage8Step,
    'FinalAdjustmentsStep': FinalAdjustmentsStep,
}


__all__ = [
    'ALL_STEP_CLASSES',
    'InitializeStep',
    'ExtractForStage2Step',
    'InputAdjustmentsStep',
    'MergeAfterStage2Step',
    'BasalSpecsStep',
    'TfConditionSpecsStep',
    'FitConditionStep',
    'PromoterBindingStep',
    'AdjustPromotersStep',
    'ExtractForStage8Step',
    'SetConditionsStep',
    'MergeAfterStage8Step',
    'FinalAdjustmentsStep',
]
