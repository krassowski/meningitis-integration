import re

from gsea_api.molecular_signatures_db import GeneMatrixTransposed
from rpy2.robjects import StrVector, ListVector
from rpy2.robjects import globalenv


def nice_reactome(name):
    return name.replace('REACTOME_', '').replace('_', ' ').title()


def nice_kegg(name):
    return name.replace('KEGG_', '').replace('_', ' ')


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


def simplify_to_match_msigdb(name):
    return name.lower().replace('(', '').replace(')', '').replace('/', ' ').replace('-', ' ').replace(':', ' ').replace('  ', ' ')


TERM_LETTER_CASE_MAP = {
    # just fix the most obvious and annoying ones
    # which would obsure readability if not fixed
    r'\bRna\b': 'RNA',
    r'\bMrna\b': 'mRNA'
}


def formatter_to_fix_letter_case(reference: GeneMatrixTransposed, custom_words_map=TERM_LETTER_CASE_MAP):
    """Returns formatter designed to improve the all-uppper-case-and-no-special-characters
    pathway names as provided by MSigDB."""
    reactome_map = {
        simplify_to_match_msigdb(pathway.name): pathway.name
        for pathway in reference.gene_sets
    }

    def try_to_correct_letter_case(pathway):
        nice = nice_reactome(pathway)
        simple = simplify_to_match_msigdb(nice)
        if simple in reactome_map:
            return reactome_map[simple]
        else:
            for pattern, substitute in custom_words_map.items():
                nice = re.sub(pattern, substitute, nice)
        return nice

    return try_to_correct_letter_case
