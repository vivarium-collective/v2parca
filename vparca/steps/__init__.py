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

from vparca.steps.step_01_initialize import InitializeStep
from vparca.steps.step_02_input_adjustments import (
    ExtractForStep2Step,
    InputAdjustmentsStep,
    MergeAfterStep2Step,
)
from vparca.steps.step_03_basal_specs import BasalSpecsStep
from vparca.steps.step_04_tf_condition_specs import TfConditionSpecsStep
from vparca.steps.step_05_fit_condition import FitConditionStep
from vparca.steps.step_06_promoter_binding import PromoterBindingStep
from vparca.steps.step_07_adjust_promoters import AdjustPromotersStep
from vparca.steps.step_08_set_conditions import (
    ExtractForStep8Step,
    SetConditionsStep,
    MergeAfterStep8Step,
)
from vparca.steps.step_09_final_adjustments import FinalAdjustmentsStep


ALL_STEP_CLASSES = {
    'InitializeStep': InitializeStep,
    'ExtractForStep2Step': ExtractForStep2Step,
    'InputAdjustmentsStep': InputAdjustmentsStep,
    'MergeAfterStep2Step': MergeAfterStep2Step,
    'BasalSpecsStep': BasalSpecsStep,
    'TfConditionSpecsStep': TfConditionSpecsStep,
    'FitConditionStep': FitConditionStep,
    'PromoterBindingStep': PromoterBindingStep,
    'AdjustPromotersStep': AdjustPromotersStep,
    'ExtractForStep8Step': ExtractForStep8Step,
    'SetConditionsStep': SetConditionsStep,
    'MergeAfterStep8Step': MergeAfterStep8Step,
    'FinalAdjustmentsStep': FinalAdjustmentsStep,
}


__all__ = [
    'ALL_STEP_CLASSES',
    'InitializeStep',
    'ExtractForStep2Step',
    'InputAdjustmentsStep',
    'MergeAfterStep2Step',
    'BasalSpecsStep',
    'TfConditionSpecsStep',
    'FitConditionStep',
    'PromoterBindingStep',
    'AdjustPromotersStep',
    'ExtractForStep8Step',
    'SetConditionsStep',
    'MergeAfterStep8Step',
    'FinalAdjustmentsStep',
]
