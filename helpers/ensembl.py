from numpy import where
from pandas import DataFrame, Series
from pyensembl import EnsemblRelease

from helpers.utilities import get_or_dummy


class Ensembl(EnsemblRelease):

    def get_gene(self, gene_id):
        return get_or_dummy(self.gene_by_id, gene_id)

    def merge_gene_data(self, df, include=['gene_name', 'biotype', 'contig', 'strand'], sort=None):
        temp = DataFrame(
            Series(
                df.index
                .map(self.get_gene)
                .map(lambda g: g.__dict__)
            )
            .tolist()
        ).set_index('gene_id')
        df = df.merge(temp[include], left_index=True, right_index=True)
        if sort:
            return df.sort_values(sort)
        return df

    def map_index(self, matrix: DataFrame, new_index='gene_name'):
        ensembl_to_gene_name = self.merge_gene_data(matrix)[[new_index]]
        return (
            matrix.index
            .map(lambda x: ensembl_to_gene_name[new_index].get(x, x))
        )

    def collapse_and_reindex_to(self, matrix: DataFrame, to: str):
        """for transposed (profiles) integration"""
        rna_gene_index = self.map_index(matrix, new_index=to)

        # NOTE: for transposed integration I could reject some isoforms
        #  which are not translated or otherwise malfucioning
        rna_matrix_isoforms_collapsed = matrix.copy()
        rna_matrix_isoforms_collapsed.index = rna_gene_index
        rna_matrix_isoforms_collapsed = rna_matrix_isoforms_collapsed.groupby(rna_gene_index).sum()

        return rna_matrix_isoforms_collapsed

    def reindex_to(self, matrix: DataFrame, to: str):

        rna_ensembl_index = matrix.index
        rna_gene_index = self.map_index(matrix, new_index=to)

        matrix.index = where(
            # if gene symbol is duplicated
            rna_gene_index.duplicated(keep=False),
            # then append the ensembl id:
            rna_gene_index + ' (' + rna_ensembl_index + ')',
            # otherwise just use the gene symbol:
            rna_gene_index
        )

        return matrix
