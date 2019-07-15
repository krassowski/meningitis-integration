from pandas import concat, Categorical, DataFrame


def prepare_opposite_correlations(
    cm_tm, rna_shared, protein_shared, protein_hc, clinical,
    n=10, additional_columns=None
):
    top_differences = cm_tm.head(n)
    protein_hc_indexed = protein_hc.set_index('index')
    df = concat([
        DataFrame(dict(
            RNA=rna_shared.loc[k][rna_shared.loc[k]>0],
            Protein=protein_shared.loc[k][rna_shared.loc[k]>0],
            hcp_q1=protein_hc_indexed.loc[k].value.quantile(.25),
            hcp_q2=protein_hc_indexed.loc[k].value.quantile(.5),
            hcp_q3=protein_hc_indexed.loc[k].value.quantile(.75),
            gene=k,
            diff=top_differences.loc[k]['diff'],
            **{
                column: top_differences.loc[k][column]
                for column in additional_columns or []
            }
        )).rename_axis('patient').reset_index()
        for k in top_differences.index
    ])
    df['Meningitis'] = df.patient.map(clinical['Meningitis'])
    df['Tuberculosis'] = df.patient.map(clinical['Tuberculosis'])
    df.gene = Categorical(df.gene, categories=top_differences.index, ordered=True)
    return df