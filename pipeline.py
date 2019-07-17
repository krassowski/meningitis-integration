from nbpipeline.rules import NotebookRule, Group

Group(
    'Proteomics',
    color='#e9eee6'
)

Group(
    'RNA',
    color='#f5e8e8'
)

Group(
    'Clinical',
    color='#cdedf6'
)

Group('Integration', color='#e9e8f6')


NotebookRule(
    'Extract SOMAScan protein data',
    input={'csf_soma_path': 'data/raw/Protein/CSF_SOMA_Hyb_RawData.xlsx'},
    output={'output_path': 'data/clean/protein/levels.csv'},
    notebook='exploration/protein/Data_extraction.ipynb',
    group='Proteomics'
)


NotebookRule(
    'Explore SOMAScan technology & check quality',
    input={'protein_levels_path': 'data/clean/protein/levels.csv'},
    output=dict(aptamers_path = 'data/other/relevant_aptamers.csv'),
    notebook='exploration/protein/Exploration_and_quality_control.ipynb',
    group='Proteomics'
)

transformed_proteins = NotebookRule(
    'Transform and normalize (double z-score)',
    input=dict(
        protein_levels_path = 'data/clean/protein/levels.csv'
    ),
    output=dict(
        indexed_by_target_path = 'data/clean/protein/indexed_by_target.csv',
        zz_log_path = 'data/clean/protein/zz_log_10.csv',
        log_matrix_path = 'data/clean/protein/log_10.csv',
        patients_variance_at_one_path = 'data/clean/protein/z_log_10-patients_variance_at_one.csv'
    ),
    notebook='exploration/protein/Transformation_and_normalization.ipynb',
    group='Proteomics'
)

NotebookRule(
    'Unsupervised analysis: Clustering, PCA',
    input=dict(
        protein_path='data/clean/protein/indexed_by_target.csv',
        # TODO: maybe this (below) is even better way of representing the pipeline?
        zz_log_path = transformed_proteins.outputs['zz_log_path'],
        # 'data/clean/protein/zz_log_10.csv',
        log_matrix_path = 'data/clean/protein/log_10.csv',
        
        clinical_path='data/clean/clinical/data.csv',
        aptamers_path='data/other/relevant_aptamers.csv',
        
        uniprot_to_go_path = 'data/ontologies/goa_human.gaf',
        gene_ontology_path = 'data/ontologies/go-basic.obo',
    ),
    notebook='exploration/protein/Unsupervised_analysis.ipynb',
    group='Proteomics'
)


proteins_mapped_to_genes = NotebookRule(
    'Map proteins to genes',
    notebook='exploration/protein/Gene_level_mapping.ipynb',
    input=dict(
        uniprot_to_entrez_path = 'data/uniprot_to_gene_id.tsv',
        aptamers_path = 'data/other/relevant_aptamers.csv',
        raw_protein_levels_path = 'data/clean/protein/indexed_by_target.csv'
    ),
    output=dict(
        raw_gene_path = 'data/clean/protein/gene_levels_by_entrez.csv',
        target_to_entrez_path = 'data/clean/protein/protein_to_entrez.csv'
    ),
    group='Proteomics'
)

NotebookRule(
    'Extract samples list and verify samples presence',
    notebook='exploration/Samples_extraction.ipynb',
    input=dict(
        raw_sample_list_path='data/raw/SampleListHead.txt',
        protein_levels_path='data/clean/protein/levels.csv',
        rna_seq_path='data/clean/rna/all_samples.csv'
    ),
    output={'clean_sample_list_path': 'data/clean/samples_list.csv'}
)

NotebookRule(
    'Extract, munge and explore clinical variables',
    notebook='exploration/clinical/Clinical_data_first_look.ipynb',
    input=dict(
        raw_clinical_path = 'data/raw/PatientClinicalData.xlsx',
        samples_path = 'data/clean/samples_list.csv'
    ),
    output=dict(
        clean_clinical_path = 'data/clean/clinical/data.csv',
        glossary_path = 'data/clean/clinical/glossary.csv'
    ),
    group='Clinical'
)


NotebookRule(
    'Survival analysis - clinical variables alone',
    notebook='analyses/Clinical_survival.ipynb',
    input=dict(
        clinical_path = 'data/clean/clinical/data_with_derived_variables.csv'
    ),
    group='Clinical'
)




NotebookRule(
    'Survival analysis - protein levels',
    notebook='analyses/protein_vs_clinical/Survival.ipynb',
    input=dict(
        clinical_path = 'data/clean/clinical/data_with_derived_variables.csv',
        zz_log_path = 'data/clean/protein/zz_log_10.csv'
    ),
    group='Proteomics'
)


NotebookRule(
    'Explain and explore PCA implementation details',
    notebook='exploration/Notes_on_PCA_with_prcomp_and_factoextra.ipynb',
    input={'path': 'data/clean/protein/levels.csv'}
)

NotebookRule(
    'Extract RNASeq and the preliminary differential expression data (DSeq2).',
    notebook='exploration/rna/Data_extraction.ipynb',
    input=dict(
        definite_tbm_rna_path = 'data/raw/RNA-Seq/DefiniteTBM_CM_VM/DefTBM_CM_VM.xlsx',
        all_samples_path = 'data/raw/RNA-Seq/AllSamples/NormalisedCounts_AllSamples_ConditionFiltered.txt',
        metadata_path = 'data/raw/RNA-Seq/AllSamples/ColData.txt'
    ),
    output=dict(
        tbm_subset_clean_path = 'data/clean/rna/definite_tbm_against_all.csv',
        duplicates_path = 'data/other/duplicates_rna_definite_tbm_subset.csv',
        definite_tbm_cm_deg_path = 'data/preliminary_analyses/deg/definite_tbm-cm.csv',
        definite_tbm_vm_deg_path = 'data/preliminary_analyses/deg/definite_tbm-vm.csv',
        cm_vm_deg_path = 'data/preliminary_analyses/deg/cm-vm.csv',
        all_samples_clean_path = 'data/clean/rna/all_samples.csv',
        all_samples_duplicates_path = 'data/other/duplicates_rna_all_samples.csv',
    ),
    group='RNA'
)

NotebookRule(
    'Match transcripts, proteins and patients',
    notebook='analyses/integration/Transcript-protein_matching.ipynb',
    group='Integration'
)

### ANALYSES

NotebookRule(
    'Analyse distributions and derive additional variables',
    input={
        'clinical_data_path': 'data/clean/clinical/data.csv'
    },
    output={
        'output_path': 'data/clean/clinical/data_with_derived_variables.csv'
    },
    notebook='analyses/Clinical_data.ipynb',
    group='Clinical'
)


#cdecf6

clinical_subset_for_proteomics = NotebookRule(
    'Subset clinical data for proteomics analyses',
    notebook='analyses/protein_vs_clinical/Subset_clinical_data_and_compare_transformations.ipynb',
    input=dict(
        indexed_by_target_path = 'data/clean/protein/indexed_by_target.csv',
        patients_variance_at_one_path ='data/clean/protein/z_log_10-patients_variance_at_one.csv',
        zz_log_path = 'data/clean/protein/zz_log_10.csv',
        log_matrix_path = 'data/clean/protein/log_10.csv',

        clinical_path = 'data/clean/clinical/data_with_derived_variables.csv'
    ),
    output=dict(
        subset_path = 'data/clean/protein/clinical_data_ordered_to_match_proteins_matrix.csv'
    ),
    group='Proteomics'
)

NotebookRule(
    'Total CSF protein',
    notebook='analyses/protein_vs_clinical/Total_CSF_protein.ipynb',
    input={
        **clinical_subset_for_proteomics.outputs,
        **dict(
            indexed_by_target_path = 'data/clean/protein/indexed_by_target.csv',
            log_matrix_path = 'data/clean/protein/log_10.csv'
        )
    },
    group='Proteomics'
)


NotebookRule(
    'Compare effects of transformations',
    notebook='analyses/protein_vs_clinical/Compare_effects_of_transformations.ipynb',
    input=dict(
        indexed_by_target_path = 'data/clean/protein/indexed_by_target.csv',
        patients_variance_at_one_path ='data/clean/protein/z_log_10-patients_variance_at_one.csv',
        zz_log_path = 'data/clean/protein/zz_log_10.csv',
        log_matrix_path = 'data/clean/protein/log_10.csv',
        
        clinical_path = 'data/clean/protein/clinical_data_ordered_to_match_proteins_matrix.csv',
    ),
    group='Proteomics'
)


NotebookRule(
    'Differential abundance of proteins',
    notebook='analyses/protein_vs_clinical/Differential_levels_and_ORA.ipynb',
    input=dict(
        indexed_by_target_path = 'data/clean/protein/indexed_by_target.csv',
        clinical_path = 'data/clean/protein/clinical_data_ordered_to_match_proteins_matrix.csv',
        log_matrix_path = 'data/clean/protein/log_10.csv',
        **proteins_mapped_to_genes.outputs
    ),
    output={
        **{
            f'out_{subset}': f'data/preliminary_analyses/differential_protein_levels/{subset}.csv'
            for subset in ['tbm', 'crypto', 'viral']
        },
        **{
            'out_log_matrix_filtered_path': 'data/clean/protein/log_10_filtered.csv'
        },
        **dict(
            out_tmm_normalized_counts_path = 'data/preliminary_analyses/differential_protein_levels/normalized_counts/tmm_for_subsets.csv',
            out_tmm_normalized_counts_gene_level_path = 'data/preliminary_analyses/differential_protein_levels/normalized_counts/gene_level_tmm_for_subsets.csv',
            out_rle_normalized_counts_gene_level_path = 'data/preliminary_analyses/differential_protein_levels/normalized_counts/gene_level_rle_for_subsets.csv'
        )
    },
    group='Proteomics'
)


NotebookRule(
    'GSEA for proteins',
    notebook='analyses/protein_vs_clinical/GSEA.ipynb',
    input=dict(
        indexed_by_target_path = 'data/clean/protein/indexed_by_target.csv',
        clinical_path = 'data/clean/protein/clinical_data_ordered_to_match_proteins_matrix.csv',
        log_matrix_path = 'data/clean/protein/log_10.csv',
        out_log_matrix_filtered_path = 'data/clean/protein/log_10_filtered.csv',
        **proteins_mapped_to_genes.outputs
    ),
    group='Proteomics'
)



NotebookRule(
    'Notes: Limma expects log-transformed data ',
    notebook='analyses/notes/Limma_expects_log_transformed_data.ipynb',
    input=dict(
        indexed_by_target_path = 'data/clean/protein/indexed_by_target.csv',
clinical_path = 'data/clean/protein/clinical_data_ordered_to_match_proteins_matrix.csv'
    ),
    group='Proteomics'
)

# TODO:
NotebookRule(
    'Differential expression',
    notebook='analyses/rnaseq_vs_clinical/Differential_expression.ipynb',
    group='RNA'
)

NotebookRule(
    'Subset clinical data for transcriptomic analyses',
    notebook='analyses/rnaseq_vs_clinical/Subset_clinical_data.ipynb',
    group='RNA'
)

NotebookRule(
    'RNAseq sanity checks',
    notebook='analyses/rnaseq_vs_clinical/RNAseq_sanity_checks.ipynb',
    group='RNA'
)


NotebookRule(
    'Limma_vs_DESeq2.ipynb',
    notebook='analyses/rnaseq_vs_clinical/_Limma_vs_DESeq2_-_comparison.ipynb',
    group='RNA'
)
