from pandas import DataFrame, Series
from pyensembl import EnsemblRelease

from helpers.utilities import get_or_dummy


class Ensembl(EnsemblRelease):
    
    def get_gene(self, gene_id):
        return get_or_dummy(self.gene_by_id, gene_id)

    def merge_gene_data(self, df, include=['gene_name', 'biotype', 'contig', 'strand'], sort=None):
        temp = DataFrame(Series(df.index.map(self.get_gene).map(lambda g: g.__dict__)).tolist()).set_index('gene_id')
        df = df.merge(temp[include], left_index=True, right_index=True)
        if sort:
            return df.sort_values(sort)
        return df