from itertools import starmap
from types import FunctionType
from typing import NamedTuple, Dict, Type

import numpy as np
from pandas import DataFrame, concat
from multiprocess import Pool


from .permutations import (
    CorrelationsPermuter, PermuterKeepingOmicsAndPhenotypes,
    PermuterKeepingOmicsAndObservationPairs,
    PermuterKeepingOmics,
)


class ComparativeCorrelations(NamedTuple):

    rna: DataFrame
    protein: DataFrame
    tb_patients: DataFrame
    cm_patients: DataFrame
    method: FunctionType
    genes_in_both: list

    permutation_constraints: Dict[str, Type[CorrelationsPermuter]] = {
        # a) keep phenotypes together
        'omics-and-phenotypes': PermuterKeepingOmicsAndPhenotypes,
        # b) keep observation pairs together
        'omics-and-observation-pairs': PermuterKeepingOmicsAndObservationPairs,
        # c) reshuffle all
        'omics': PermuterKeepingOmics
    }

    def _fast_correlation(self, rna, protein, filter_out_low_rna=True):
        assert all(rna.index == protein.index)

        # order of columns is not guaranteed to be the same (due to permutations),
        # here a map from rna columns to protein columns is created
        rna_to_protein_order = np.array([
            protein.columns.get_loc(rna_patient)
            for rna_patient in rna.columns
        ])

        method = self.method

        def _fast_gene_correlation(r, p, patients_with_enough_rna):
            if not all(patients_with_enough_rna):
                p = p[rna_to_protein_order[patients_with_enough_rna]]
                r = r[patients_with_enough_rna]
            else:
                # reorder to rna patients (all are true)
                p = p[rna_to_protein_order[patients_with_enough_rna]]
            cor = method(r, p)
            return cor[0], cor[1], len(r)

        return DataFrame(
            starmap(
                _fast_gene_correlation,
                zip(
                    rna.values, protein.values,
                    (rna.values > 0) if filter_out_low_rna else ([False])
                )
            ),
            columns=['correlation', 'correlation_pvalue', 'n'],
            index=rna.index
        )

    def random_permutation(self, i, permuter, metric):
        permutation = permuter.permute()

        correlations_tb = self._fast_correlation(permutation.tb_rna, permutation.tb_protein)
        correlations_cm = self._fast_correlation(permutation.cm_rna, permutation.cm_protein)

        return metric(correlations_cm, correlations_tb)

    def compute_null_distribution(self, metric, constraint, n):

        assert constraint in self.permutation_constraints

        permuter_type = self.permutation_constraints[constraint]

        permuter = permuter_type(
            rna=self.rna,
            protein=self.protein,
            tb_patients=self.tb_patients,
            cm_patients=self.cm_patients,
            genes_in_both=self.genes_in_both
        )

        pool = Pool()
        null = pool.imap(self.random_permutation, range(n), shared_args=(permuter, metric), total=n)

        return concat(null, axis=1)
