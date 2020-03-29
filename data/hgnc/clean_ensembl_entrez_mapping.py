from pandas import read_table

in_ensembl_to_entrez_path = 'ensembl_gene_to_entrez.tsv'
ensembl_to_entrez_path = '../ensembl_to_entrez.csv'

ensembl_to_entrez = read_table(in_ensembl_to_entrez_path).dropna().set_index('Ensembl gene ID').astype(int).astype(str)
ensembl_to_entrez.to_csv(ensembl_to_entrez_path)
print(ensembl_to_entrez.head())
