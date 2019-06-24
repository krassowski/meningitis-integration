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


def generate_patient_annotations(clinical, hiv=True, fillna='?'):
    return DataFrame({
        **{
            'Meningitis': clinical.condition.replace(conditions_names),
            'Tuberculosis status': clinical.condition.map(tuberculosis_status).fillna('-')
        },
        **(
            {'HIV status': clinical.HIVResult}
            if hiv else
            {}
        )
    }).set_index(clinical.index).fillna(fillna)
