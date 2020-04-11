from functools import partial
from pandas import DataFrame


class ProteinOntologyClassifier:
    
    def __init__(self, ontology, associations):
        self.ontology = ontology
        self.associations = associations

    def classify_by(self, parent_term):
        def classifier(uniprot):
            processes = set()
            associations = (
                self.associations[uniprot]
                if uniprot in self.associations else
                []
            )
            for go in associations:
                if go not in self.ontology:
                    print('skipping ' + go)
                    continue
                all_parents = self.ontology[go].get_all_parents()
                for parent in all_parents:
                    parent = self.ontology[parent]
                    if any(
                        grandparent.name == parent_term
                        for grandparent in parent.parents
                    ):
                        processes.add(parent.name)
            return processes
        return classifier
    
    def classify_all(self, protein_ids, parent_term, multi=False):
        return classes_as_columns(
            protein_ids.apply(self.classify_by(parent_term)),
            parent_term if multi else None
        )
    
    def classify(self, index, proteins, parent_term):
        return DataFrame({
            'target': index,
            **self.classify_all(
                proteins,
                parent_term
            )
        }).set_index('target')


def classes_as_columns(processes, name=None):
    unique_processes = set.union(*processes.to_list())
    return {
        (name, process) if name else process: [
            process in target_processes
            for target_processes in processes
        ]
        for process in unique_processes
    }


def transform_patient_proteins(protein_values_of_patient, ordered_processes):
    v = ordered_processes.apply(
        lambda process_presence: process_presence * protein_values_of_patient
    ).sum()
    
    return v


def transform_to_classes(data, classification):
    # spread the signal from gene across processes
    relative_importance = classification.div(classification.sum(axis=1), axis=0)
    signal_lost = relative_importance.isna().any(axis=1).sum()
    percentage = signal_lost / len(relative_importance) * 100
    print(f'Going to loose signal from {percentage:.2f}% ({signal_lost}) proteins')
    relative_importance = relative_importance.fillna(0)
    
    # reorder proteins to match the protein_values_of_patient
    ordered_processes = relative_importance.reindex(data.index)
    transformation = partial(transform_patient_proteins, ordered_processes=ordered_processes)
    reduced_to_processes = data.apply(transformation)
    
    return reduced_to_processes
