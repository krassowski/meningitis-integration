from pandas import concat


def cm_vs_tmdr_or_tms_roc_curves(model, tms_test_set):

    tms_result = model.validate(tms_test_set)
    tms_result_sub = model.cross_validation.validate(tms_test_set)

    return concat([
        model.cross_validation_results.roc_auc.assign(group='CV train: TMD+TMR vs CM', models='Cross-Validation', set='Train: TMD+TMR vs CM'),
        model.sub_sampling_test_results.roc_auc.assign(group='CV test: TMD+TMR vs CM', models='Cross-Validation', set='Test: TMD+TMR vs CM'),
        tms_result_sub.roc_auc.assign(group='CV test: TMS vs CM', models='Cross-Validation', set='Test: TMS vs CM'),
        tms_result.roc_auc.assign(group='Full test: TMS vs CM', models='Full model', set='Test: TMS vs CM'),
        model.test_result.roc_auc.assign(group='Full test: TMD+TMR vs CM', models='Full model', set='Test: TMD+TMR vs CM')
    ])