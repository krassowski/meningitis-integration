from itertools import combinations

from pandas import Categorical, concat
from scipy.stats import spearmanr


de_labels = {
    'deseq2_Rachel': 'DESeq2 from dr Rachel',
    **{
        method.lower().replace('-', '_') + id_suffix: method + suffix
        for method in ['DESeq2', 'voom-TMM', 'voom-RLE', 'voom-EEE', 'voom-qtotal']
        for id_suffix, suffix in {
            '_no_filter': ' without filtering',
            '_filtered': ' filtered',
            '_reproduction': ' reproduced',
            '_weighted': ' weighted',
            '_shrinkage_normal': ' "normal" shrinkage',
            '_shrinkage_apeglm': ' "apeglm" shrinkage',
            '_shrinkage_ashr': ' "ashr" shrinkage'
        }.items()
    }
}

deseq_cols = ['pvalue', 'padj']
voom_cols = {
    'P.Value': 'pvalue',
    'adj.P.Val': 'padj'
}



def significant_set(p_values, threshold=0.05):
    return set(p_values.query(f'padj < {threshold}').index)


def significant_union(p_values, a, b, threshold=0.05):
    return (
        significant_set(p_values[a], threshold)
        |
        significant_set(p_values[b], threshold)
    )


def get_ranks(p_values, method, subset, labels):
    return (
        p_values[method]
        .dropna()
        .reindex(subset)
        .rename_axis('gene').reset_index()
        .sort_values(['padj', 'pvalue', 'gene'])
        .set_index('gene')
        .rank(method='dense', numeric_only=True)
        .assign(method=labels[method])
    ).rename(columns={'pvalue': 'pvalue_rank'})


def corr_label(res):
    return f'corr = {res.correlation:.2f}, p = {res.pvalue:.2f}'


def generate_comparison(p_values, methods, threshold=0.05, labels=de_labels, method=spearmanr, nan_policy='omit', rename=None):
    df = concat([
        concat([ranks_a, ranks_b]).assign(
            contrast=f'{labels[a]} - {labels[b]}',
            corr=corr_label(method(
                ranks_a.loc[significanct_in_either]['pvalue_rank'],
                ranks_b.loc[significanct_in_either]['pvalue_rank'],
                nan_policy=nan_policy
            ))
        )
        for a, b in combinations(sorted(methods), 2)
        for a, b in [sorted([a, b])]
        for significanct_in_either in [significant_union(p_values, a, b, threshold)]
        for ranks_a, ranks_b in [[
            get_ranks(p_values, a, significanct_in_either, labels).assign(side=-1),
            get_ranks(p_values, b, significanct_in_either, labels).assign(side=+1)
        ]]
    ]).reset_index()
    if rename is not None:
        df['name'] = rename.loc[df.gene].values
    df = df.sort_values('method')
    df['method'] = Categorical(df['method'], ordered=True)
    df = df.sort_values('contrast')
    df['contrast'] = Categorical(df['contrast'], ordered=True)
    return df
