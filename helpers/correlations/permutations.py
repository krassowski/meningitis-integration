import random
from abc import ABC, abstractmethod
from typing import NamedTuple

from pandas import DataFrame, concat


class PermutationResult(NamedTuple):
    tb_rna: DataFrame
    cm_rna: DataFrame
    tb_protein: DataFrame
    cm_protein: DataFrame


class CorrelationsPermuter(ABC):

    def __init__(self, rna, protein, tb_patients, cm_patients, genes_in_both):
        self.tb_patients = list(tb_patients)
        self.cm_patients = list(cm_patients)

        # reorder RNA and protein frames to match genes and patients order
        self._tb_rna = rna[tb_patients].loc[genes_in_both]
        self._cm_rna = rna[cm_patients].loc[genes_in_both]

        self._tb_protein = protein[tb_patients].loc[genes_in_both]
        self._cm_protein = protein[cm_patients].loc[genes_in_both]

    @abstractmethod
    def permute(self) -> PermutationResult:
        pass


class PermuterKeepingOmicsAndPhenotypes(CorrelationsPermuter):
    """Fast, in-place permutation of the labels of RNA and Protein frames.

    Important: because the permutation is done in-place, the previously returned result changes with each permutation.
    """

    def permute(self):
        # it is sufficient to shuffle patient identifiers for one of (rna, protein) frames only
        random.shuffle(self.tb_patients)
        random.shuffle(self.cm_patients)

        self._tb_protein.columns = self.tb_patients
        self._cm_protein.columns = self.cm_patients

        return PermutationResult(
            tb_rna=self._tb_rna,
            cm_rna=self._cm_rna,
            tb_protein=self._tb_protein,
            cm_protein=self._cm_protein
        )


class PhenotypePermuter(CorrelationsPermuter):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._all_rna = concat([self._tb_rna, self._cm_rna], axis='columns')
        self._all_protein = concat([self._tb_protein, self._cm_protein], axis='columns')
        self.all_patients = self.tb_patients + self.cm_patients

        # keep memory footprint low
        del self._tb_rna, self._cm_rna, self._tb_protein, self._cm_protein

    @property
    @abstractmethod
    def keep_observation_pairs(self) -> int:
        pass

    def permute(self):
        all_patients = self.all_patients
        tb_patients = self.tb_patients

        random.shuffle(all_patients)
        rna_random_tb = all_patients[:len(tb_patients)]    # select patients for TB
        rna_random_cm = all_patients[len(tb_patients):]    # all other patients go to CM

        if self.keep_observation_pairs:
            protein_random_tb = rna_random_tb
            protein_random_cm = rna_random_cm
        else:
            protein_random_tb = rna_random_tb[:]
            protein_random_cm = rna_random_cm[:]
            random.shuffle(protein_random_tb)
            random.shuffle(protein_random_cm)

        return PermutationResult(
            tb_rna=self._all_rna[rna_random_tb],
            cm_rna=self._all_rna[rna_random_cm],
            tb_protein=self._all_protein[protein_random_tb],
            cm_protein=self._all_protein[protein_random_cm]
        )


class PermuterKeepingOmicsAndObservationPairs(PhenotypePermuter):

    keep_observation_pairs = True


class PermuterKeepingOmics(PhenotypePermuter):

    keep_observation_pairs = False
