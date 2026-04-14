"""
Process-bigraph Step classes for the ParCa pipeline, one module per step.

Each step module exposes its Step class and the pure helpers it delegates
to.  Every Step declares **leaf-level ports**: one port per sim_data (or
cell_specs) field it reads or writes.  The composite wires port names to
store paths in a nested bigraph that mirrors sim_data's structure.

NOTE: this is mid-redesign.  Step 2 (``InputAdjustmentsStep``) has been
converted to the leaf-port model.  Steps 1, 3–9 still use the older
``extract_input`` / ``merge_output`` wrappers around a ``ParcaState`` blob;
they will be converted next.  ``ALL_STEP_CLASSES`` currently reflects the
transitional state.
"""

from vparca.steps.step_01_initialize import InitializeStep
from vparca.steps.step_02_input_adjustments import InputAdjustmentsStep
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
    'InputAdjustmentsStep': InputAdjustmentsStep,
    'BasalSpecsStep': BasalSpecsStep,
    'TfConditionSpecsStep': TfConditionSpecsStep,
    'FitConditionStep': FitConditionStep,
    'PromoterBindingStep': PromoterBindingStep,
    'AdjustPromotersStep': AdjustPromotersStep,
    # step 8 is still on the old extract/compute/merge triplet
    'ExtractForStep8Step': ExtractForStep8Step,
    'SetConditionsStep': SetConditionsStep,
    'MergeAfterStep8Step': MergeAfterStep8Step,
    'FinalAdjustmentsStep': FinalAdjustmentsStep,
}


__all__ = [
    'ALL_STEP_CLASSES',
    'InitializeStep',
    'InputAdjustmentsStep',
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
