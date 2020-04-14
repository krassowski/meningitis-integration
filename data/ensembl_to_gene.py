from pyensembl import EnsemblRelease
from pandas import DataFrame

ensembl = EnsemblRelease(95)

df = DataFrame(
    [
        [gene.gene_id, gene.gene_name]
        for gene in ensembl.genes()
    ],
    columns=['gene_id', 'gene_name']
)

df.to_csv('ensembl_to_gene_symbol.csv', index=False)
