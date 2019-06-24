from rpy2.robjects import StrVector, ListVector
from rpy2.robjects import globalenv


def nice_reactome(name):
    return name.replace('REACTOME_', '').replace('_', ' ').title()


def calculate_overlap(gmt, matrix):
    all_identifiers = gmt.all_identifiers

    soma = set(map(str, matrix.index))
    precentage = len(all_identifiers & soma) / len(soma) * 100
    print(f'{precentage:.2f}%')


def collection_to_R(collection, trim_to, min=5, max=500):
    gene_ids = trim_to
    filtered = {
        gene_set.name: StrVector(list(gene_set.genes))
        for gene_set in (
            collection
            # limma::cameraPR goes crazy without this
            # limma::mroast seems to work fine (and be more aware of the limitted statistical support)
            .subset(gene_ids)
            .gene_sets
        )
        if max > len(gene_set.genes) > min
    }
    gene_sets_r = ListVector(filtered)

    globalenv[collection.name] = gene_sets_r
    
    return filtered

