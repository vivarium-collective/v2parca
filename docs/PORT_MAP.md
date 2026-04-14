# vParCa Port Map — per-step sim_data leaves

Single source of truth for the port-first redesign. Derived by reading every
`extract_input()` and `merge_output()` in `vparca/steps/step_NN_*.py` as of
the last extract/merge-paradigm commit. Every leaf is given its **dotted
path into sim_data** (or `cell_specs`), which becomes the store path the
Step's port is wired to.

Path notation:
- `process.transcription.rna_data.deg_rate` = attribute of `rna_data`'s
  structured array, field `"deg_rate"`.
- `process.transcription.rna_expression["basal"]` = entry at key `"basal"`
  of the `rna_expression` dict.
- `cell_specs["basal"]["bulkContainer"]` = `cell_specs` tree entry.

Leaves that are Python objects with *methods the step invokes* (e.g.
`process.transcription` as a `Transcription` instance) are marked
**OBJECT** — those remain opaque leaves in the store because we don't
have a tractable way to expose their method surface as data ports. Every
other field is a value leaf.

## Step 1 — initialize

Scatters the outputs of `SimulationDataEcoli.initialize(raw_data=...)` into
every store leaf that any of steps 2–9 reads (see aggregate at the bottom).
Step 1 runs the vendored `sim_data.initialize()` internally as machinery
but does not expose `sim_data` as a port.

## Step 2 — input_adjustments  (PURE)

READS:
- `process.transcription.cistron_data.id`                                  — array[str]
- `process.transcription.cistron_id_to_rna_indexes_map`                    — dict[str, array[int]] *(derived; computed at init)*
- `process.translation.monomer_data.id`                                    — array[str]
- `process.translation.translation_efficiencies_by_monomer`                — array[float]
- `adjustments.translation_efficiencies_adjustments`                       — dict[str, float]
- `adjustments.balanced_translation_efficiencies`                          — list[list[str]]
- `process.transcription.rna_data.id`                                      — array[str]
- `process.transcription.rna_expression["basal"]`                          — array[float]
- `adjustments.rna_expression_adjustments`                                 — dict[str, float]
- `process.transcription.rna_data.deg_rate`                                — array[float]
- `process.transcription.cistron_data.deg_rate`                            — array[float]
- `adjustments.rna_deg_rates_adjustments`                                  — dict[str, float]
- `process.translation.monomer_data.deg_rate`                              — array[float]
- `adjustments.protein_deg_rates_adjustments`                              — dict[str, float]
- `tf_to_active_inactive_conditions`                                       — dict[str, dict]
- config `debug`                                                           — bool

WRITES:
- `process.translation.translation_efficiencies_by_monomer`                — array[float]
- `process.transcription.rna_expression["basal"]`                          — array[float]
- `process.transcription.rna_data.deg_rate`                                — array[float]
- `process.transcription.cistron_data.deg_rate`                            — array[float]
- `process.translation.monomer_data.deg_rate`                              — array[float]
- `tf_to_active_inactive_conditions`                                       — dict (only if debug reduced it)

## Step 3 — basal_specs  (COUPLED)

READS (data leaves — `OBJECT` ports elided for brevity):
- `process.metabolism`                                                     — **OBJECT** (used for `concentration_updates.concentrations_based_on_nutrients("minimal")`)
- `process.transcription`                                                  — **OBJECT** (used for `set_ppgpp_expression` + `rna_data` struct access)
- `process.rna_decay`                                                      — **OBJECT** (used for endoRNase Km cooperative fitting)
- `process.transcription.rna_expression["basal"]`                          — array[float] (copied)
- `condition_to_doubling_time["basal"]`                                    — scalar⋅time
- `mass.avg_cell_dry_mass_init`                                            — scalar⋅mass
- `mass.avg_cell_to_initial_cell_conversion_factor`                        — float
- `mass.cell_dry_mass_fraction`                                            — float
- `mass.cell_water_mass_fraction`                                          — float
- `constants.cell_density`                                                 — scalar with units
- `constants.n_avogadro`                                                   — scalar with units
- `constants.growth_associated_maintenance`                                — scalar with units
- `constants.sensitivity_analysis_kcat_endo`                               — bool
- `molecule_groups.endoRNase_rnas`                                         — array
- `internal_state.bulk_molecules.bulk_data.id`                             — array[str]
- kwargs `variable_elongation_transcription`, `variable_elongation_translation`,
  `disable_ribosome_capacity_fitting`, `disable_rnapoly_capacity_fitting`,
  `cache_dir`

WRITES:
- `cell_specs["basal"].concDict`                                           — dict
- `cell_specs["basal"].expression`                                         — array
- `cell_specs["basal"].synthProb`                                          — array
- `cell_specs["basal"].fit_cistron_expression`                             — array
- `cell_specs["basal"].doubling_time`                                      — scalar⋅time
- `cell_specs["basal"].avgCellDryMassInit`                                 — scalar⋅mass
- `cell_specs["basal"].fitAvgSolubleTargetMolMass`                         — scalar⋅mass
- `cell_specs["basal"].bulkContainer`                                      — structured-array
- `mass.avg_cell_dry_mass_init`                                            — scalar⋅mass (updated)
- `mass.avg_cell_dry_mass`                                                 — scalar⋅mass
- `mass.avg_cell_water_mass_init`                                          — scalar⋅mass
- `mass.fitAvgSolubleTargetMolMass`                                        — scalar⋅mass
- `process.transcription.rna_expression["basal"]`                          — array
- `process.transcription.rna_synth_prob["basal"]`                          — array
- `process.transcription.fit_cistron_expression["basal"]`                  — array
- `process.transcription.rna_data.Km_endoRNase`                            — array
- `process.transcription.mature_rna_data.Km_endoRNase`                     — array
- `process.rna_decay.Km_first_order_decay`                                 — array
- `process.rna_decay.sensitivity_analysis_alpha_residual`                  — dict
- `process.rna_decay.sensitivity_analysis_kcat_res_ini`                    — dict
- `process.rna_decay.sensitivity_analysis_kcat_res_opt`                    — dict
- `process.rna_decay.stats_fit`                                            — dict
- `constants.darkATP`                                                      — scalar

## Step 4 — tf_condition_specs  (COUPLED; per-condition fan-out)

READS:
- `tf_to_active_inactive_conditions`                                       — dict
- `conditions[<label>]`                                                    — dict (per condition)
- `conditions[<label>].perturbations`                                      — list
- `conditions[<label>].nutrients`                                          — str
- `condition_active_tfs`                                                   — dict
- `condition_inactive_tfs`                                                 — dict
- `tf_to_fold_change[<tf>]`                                                — dict
- `process.transcription`                                                  — **OBJECT**
- `process.metabolism`                                                     — **OBJECT**
- `mass`                                                                   — **OBJECT**
- `condition_to_doubling_time[<label>]`                                    — scalar⋅time
- `process.transcription.rna_data.Km_endoRNase`                            — array
- kwargs `cpus`, elongation + capacity-fitting flags

WRITES (per condition label ℓ):
- `cell_specs[ℓ].concDict`                                                 — dict
- `cell_specs[ℓ].expression`                                               — array
- `cell_specs[ℓ].synthProb`                                                — array
- `cell_specs[ℓ].fit_cistron_expression`                                   — array
- `cell_specs[ℓ].doubling_time`                                            — scalar⋅time
- `cell_specs[ℓ].avgCellDryMassInit`                                       — scalar⋅mass
- `cell_specs[ℓ].fitAvgSolubleTargetMolMass`                               — scalar⋅mass
- `cell_specs[ℓ].bulkContainer`                                            — structured-array
- `cell_specs[ℓ].cistron_expression`                                       — array (non-basal only)
- `process.transcription.rna_expression[ℓ]`                                — array
- `process.transcription.rna_synth_prob[ℓ]`                                — array
- `process.transcription.cistron_expression[ℓ]`                            — array
- `process.transcription.fit_cistron_expression[ℓ]`                        — array

## Step 5 — fit_condition  (READ-ONLY on sim_data)

READS:
- `cell_specs[<label>].expression`                                         — array
- `cell_specs[<label>].concDict`                                           — dict
- `cell_specs[<label>].avgCellDryMassInit`                                 — scalar⋅mass
- `cell_specs[<label>].doubling_time`                                      — scalar⋅time
- `conditions[<label>].nutrients`                                          — str
- `process.complexation`                                                   — **OBJECT**
- `process.equilibrium`                                                    — **OBJECT**
- `process.two_component_system`                                           — **OBJECT**
- `process.translation.monomer_data.id`                                    — array
- `process.translation.monomer_data.aa_counts`                             — array
- `internal_state.bulk_molecules.bulk_data.id`                             — array
- `constants.cell_density`, `constants.n_avogadro`                         — scalars
- `mass.cell_dry_mass_fraction`                                            — float
- kwargs `cpus`

WRITES:
- `cell_specs[ℓ].bulkAverageContainer`                                     — structured-array
- `cell_specs[ℓ].bulkDeviationContainer`                                   — structured-array
- `cell_specs[ℓ].proteinMonomerAverageContainer`                           — structured-array
- `cell_specs[ℓ].proteinMonomerDeviationContainer`                         — structured-array
- `cell_specs[ℓ].translation_aa_supply`                                    — array
- `translation_supply_rate[<nutrients>]`                                   — array

## Step 6 — promoter_binding  (COUPLED)

READS:
- `process.equilibrium`                                                    — **OBJECT**
- `process.transcription_regulation`                                       — **OBJECT**
- `process.replication`                                                    — **OBJECT**
- `cell_specs[<condition>].bulkAverageContainer`                           — structured-array
- `cell_specs[<condition>].doubling_time`                                  — scalar⋅time
- `cell_specs[<condition>].avgCellDryMassInit`                             — scalar⋅mass

WRITES:
- `cell_specs["basal"].r_vector`                                           — array
- `cell_specs["basal"].r_columns`                                          — dict
- `process.transcription_regulation.pPromoterBound[…]`                     — array assignments
- `process.transcription.rna_synth_prob[…]`                                — array assignments

## Step 7 — adjust_promoters  (COUPLED)

READS:
- `process.equilibrium`                                                    — **OBJECT**
- `process.metabolism`                                                     — **OBJECT**
- `process.transcription_regulation`                                       — **OBJECT**
- `cell_specs[<condition>]` — `r_vector`, `r_columns`, `bulkAverageContainer`, `avgCellDryMassInit`

WRITES:
- `process.transcription_regulation.basal_prob`                            — array
- `process.transcription_regulation.delta_prob`                            — dict
- `process.metabolism.molecule_set_amounts[…]`                             — assignments
- `process.equilibrium.reverse_rates[…]`                                   — assignments

## Step 8 — set_conditions  (PURE)

READS:
- `process.transcription.rna_data.is_mRNA`, `is_tRNA`, `is_rRNA`           — bool masks
- `process.transcription.rna_data.includes_ribosomal_protein`, `includes_RNAP` — bool masks
- per condition ℓ:
  - `conditions[ℓ].nutrients`, `conditions[ℓ].perturbations`
  - `condition_to_doubling_time[ℓ]`
  - `cell_specs[ℓ].bulkContainer`, `.avgCellDryMassInit`
  - `process.metabolism.concentration_updates.concentrations_based_on_nutrients(<nutrients>)` — **OBJECT** method
  - `mass.getBiomassAsConcentrations(dt)`, `mass.get_component_masses(dt)` — **OBJECT** methods
  - `mass.avg_cell_to_initial_cell_conversion_factor`                      — float
  - `constants.cell_density`, `constants.n_avogadro`                       — scalars
  - `growth_rate_parameters`                                               — **OBJECT**
  - `process.transcription.rna_synth_prob[ℓ]`                              — array

WRITES:
- `process.transcription.rnaSynthProbFraction`                             — dict
- `process.transcription.rnapFractionActiveDict`                           — dict
- `process.transcription.rnaSynthProbRProtein`                             — dict
- `process.transcription.rnaSynthProbRnaPolymerase`                        — dict
- `process.transcription.rnaPolymeraseElongationRateDict`                  — dict
- `expectedDryMassIncreaseDict`                                            — dict
- `process.translation.ribosomeElongationRateDict`                         — dict
- `process.translation.ribosomeFractionActiveDict`                         — dict
- per condition ℓ:
  - `cell_specs[ℓ].avgCellDryMassInit`, `.fitAvgSolublePoolMass`, `.bulkContainer`

## Step 9 — final_adjustments  (COUPLED, OBJECT-dominated)

READS:
- `process.transcription`, `process.metabolism`, `constants`               — **OBJECT**s
- `cell_specs`                                                             — dict-of-dicts

WRITES (all via method calls on the objects above):
- wide surface in `process.transcription.*`, `process.metabolism.*`, `constants.*`

---

## Step 1 aggregate — leaves to populate at init

Union of the leaf paths any of steps 2–9 reads (OBJECT leaves implicit).
This is the complete set of named ports on `InitializeStep`'s outputs.

_(see per-step sections above; this section is the union by construction.)_
