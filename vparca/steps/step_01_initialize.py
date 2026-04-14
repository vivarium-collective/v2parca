"""
Step 1 — initialize (scatter).

Runs ``SimulationDataEcoli.initialize(raw_data=...)`` internally and scatters
the resulting leaves out through explicit output ports.  No ``sim_data`` or
``cell_specs`` blob leaves this Step; from here on, every downstream Step
reads and writes leaves directly through the nested bigraph store.

The scatter's output surface is exactly the union of leaves any downstream
Step reads.  **This is the PoC version that covers Step 2 only**; as
subsequent steps are converted to leaf-ports the scatter grows.  Store
paths wired by the composite:

  raw_data                           (input)    raw_data

  OUTPUTS (port → store path)
    monomer_ids                       process / translation / monomer_data / id
    translation_efficiencies          process / translation / translation_efficiencies_by_monomer
    translation_eff_adjustments       adjustments / translation_efficiencies_adjustments
    balanced_translation_groups       adjustments / balanced_translation_efficiencies
    rna_ids                           process / transcription / rna_data / id
    cistron_ids                       process / transcription / cistron_data / id
    basal_rna_expression              process / transcription / rna_expression / basal
    rna_expression_adjustments        adjustments / rna_expression_adjustments
    cistron_id_to_rna_indexes         process / transcription / cistron_id_to_rna_indexes_map
    rna_deg_rates                     process / transcription / rna_data / deg_rate
    cistron_deg_rates                 process / transcription / cistron_data / deg_rate
    rna_deg_rate_adjustments          adjustments / rna_deg_rates_adjustments
    protein_deg_rates                 process / translation / monomer_data / deg_rate
    protein_deg_rate_adjustments      adjustments / protein_deg_rates_adjustments
    tf_to_active_inactive_conditions  tf_to_active_inactive_conditions

It is idiomatic for a ``Step`` with no inputs (besides ``raw_data``) to
execute once at the head of the composite DAG.
"""

import time

from process_bigraph import Step

from vparca.reconstruction.ecoli.simulation_data import SimulationDataEcoli


# Every scatter output — keep in sync with what downstream steps read.
OUTPUT_PORTS = {
    'monomer_ids':                      'overwrite',
    'translation_efficiencies':         'overwrite',
    'translation_eff_adjustments':      'overwrite',
    'balanced_translation_groups':      'overwrite',
    'rna_ids':                          'overwrite',
    'cistron_ids':                      'overwrite',
    'basal_rna_expression':             'overwrite',
    'rna_expression_adjustments':       'overwrite',
    'cistron_id_to_rna_indexes':        'overwrite',
    'rna_deg_rates':                    'overwrite',
    'cistron_deg_rates':                'overwrite',
    'rna_deg_rate_adjustments':         'overwrite',
    'protein_deg_rates':                'overwrite',
    'protein_deg_rate_adjustments':     'overwrite',
    'tf_to_active_inactive_conditions': 'overwrite',
}


class InitializeStep(Step):
    """Run ``sim_data.initialize(raw_data=...)`` and scatter its leaves.

    ``raw_data`` (a ``KnowledgeBaseEcoli`` instance) is delivered via
    ``config`` rather than through a store port, so bigraph-schema doesn't
    try to type-infer the nested KB structure at composite construction
    time.  The Step has **no inputs** — it's the head of the DAG.
    """

    config_schema = {
        'raw_data':                   'overwrite',
        'basal_expression_condition': {
            '_type': 'string',
            '_default': 'M9 Glucose minus AAs',
        },
    }

    def inputs(self):
        return {}

    def outputs(self):
        return dict(OUTPUT_PORTS)

    def update(self, state):
        t0 = time.time()
        raw_data = self.config['raw_data']

        # Build and initialize sim_data.  This is the one place sim_data
        # exists; from here, we scatter it into store leaves.
        sim_data = SimulationDataEcoli()
        sim_data.initialize(
            raw_data=raw_data,
            basal_expression_condition=self.config.get(
                'basal_expression_condition', 'M9 Glucose minus AAs'),
        )

        transcription = sim_data.process.transcription
        translation   = sim_data.process.translation
        adjustments   = sim_data.adjustments

        # Precompute the cistron→rna-indexes mapping — it's derived from
        # transcription's internal structure and feeding it through the
        # store as a plain dict means step 2 doesn't need the live
        # Transcription object just for this lookup.
        cistron_ids = transcription.cistron_data['id']
        cistron_id_to_rna_indexes_map = {
            cid: transcription.cistron_id_to_rna_indexes(cid)
            for cid in cistron_ids
        }

        out = {
            # ---- translation ----
            'monomer_ids': translation.monomer_data['id'],
            'translation_efficiencies':
                translation.translation_efficiencies_by_monomer.copy(),
            'protein_deg_rates':
                translation.monomer_data.struct_array['deg_rate'].copy(),

            # ---- transcription ----
            'rna_ids':                 transcription.rna_data['id'],
            'cistron_ids':             cistron_ids,
            'basal_rna_expression':    transcription.rna_expression['basal'].copy(),
            'rna_deg_rates':
                transcription.rna_data.struct_array['deg_rate'].copy(),
            'cistron_deg_rates':
                transcription.cistron_data.struct_array['deg_rate'].copy(),
            'cistron_id_to_rna_indexes': cistron_id_to_rna_indexes_map,

            # ---- adjustments (copied to decouple from live sim_data) ----
            'translation_eff_adjustments':
                dict(adjustments.translation_efficiencies_adjustments),
            'balanced_translation_groups':
                list(adjustments.balanced_translation_efficiencies),
            'rna_expression_adjustments':
                dict(adjustments.rna_expression_adjustments),
            'rna_deg_rate_adjustments':
                dict(adjustments.rna_deg_rates_adjustments),
            'protein_deg_rate_adjustments':
                dict(adjustments.protein_deg_rates_adjustments),

            # ---- top-level ----
            'tf_to_active_inactive_conditions':
                dict(sim_data.tf_to_active_inactive_conditions),
        }

        print(f"  Step 1 (initialize + scatter) completed in {time.time() - t0:.1f}s")
        return out
