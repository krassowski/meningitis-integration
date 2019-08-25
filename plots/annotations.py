from pandas import DataFrame

conditions_names = {
    'CM': 'Cryptococcal',
    'HC': 'Healthy control',
    'TMD': 'Tuberculosis',
    'TMR': 'Tuberculosis',
    'TMS': 'Tuberculosis',
    'VM': 'Viral',
    'BM': 'Bacterial',  # no such patients for protein data, but there are some in general
}

tuberculosis_status = {
    'TMD': 'Definite',
    'TMR': 'Probable',
    'TMS': 'Possible',
}


def generate_patient_annotations(clinical, hiv=True, fillna='?', **kwargs):
    if hiv:
        kwargs['HIV status'] = 'HIVResult'
    return DataFrame({
        **{
            'Meningitis': clinical.condition.replace(conditions_names),
            'Tuberculosis status': clinical.condition.map(tuberculosis_status).fillna('-')
        },
        **{
            name: getattr(clinical, column)
            for name, column in kwargs.items()
        }
    }).set_index(clinical.index).fillna(fillna)
