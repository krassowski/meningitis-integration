from functools import partial
from statistics import mean

from pandas import DataFrame, Series
import numpy as np
from scipy import interp
from sklearn import metrics
from helpers.r import r_function, r
from tqdm.auto import tqdm


r('library(pROC)')
r('library(cvAUC)')

test_roc = partial(r_function, 'roc.test')
auc = partial(r_function, 'auc')
ci_se = partial(r_function, 'ci.se')
ci_auc = partial(r_function, 'ci.auc')
roc_object = partial(r_function, 'roc')
ci_cv_auc = partial(r_function, 'ci.cvAUC')
cv_auc = partial(r_function, 'cvAUC')
power_roc_test = partial(r_function, 'power.roc.test')


def compare_roc_curves(
    a_response, a_prediction,
    b_response, b_prediction,
    paired=False,
    compute_power=False,
    power_level=0.8
):
    roc1 = roc_object(a_response, a_prediction, quiet=True)
    roc2 = roc_object(b_response, b_prediction, quiet=True)

    p = test_roc(roc1, roc2, paired=paired).rx2('p.value')[0]

    if compute_power:
        assert paired
        # get power at significance level 0.05 for two-sided alternative
        power = power_roc_test(roc1, roc2)
        # get number of patients required to get a=0.05 and B=0.8
        n_at_power = power_roc_test(roc1, roc2, power=power_level)

        power_report = {
            f'cases_needed_for_{power_level}': n_at_power.rx2('ncases')[0],
            f'controls_needed_for_{power_level}': n_at_power.rx2('ncontrols')[0],
            **{
                k: v[0]
                for k, v in power.items()
            }
        }
    else:
        power_report = None

    return p, power_report


def roc_auc_plot_data(probabilities, true_responses):

    aucs, mean_fpr, random_tprs, tprs = compute_cv_statistics(probabilities, true_responses)
    mean_tpr = np.mean(tprs, axis=0)

    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    
    ci_cv = ci_cv_auc(probabilities, true_responses)
    cv_mean_auc = ci_cv.rx2('cvAUC')
    if cv_mean_auc != mean(aucs):
        print(f'Warning: AUC computed by cvAUC differs from mean: {cv_mean_auc} vs {mean(aucs)}')

    std_auc = np.std(aucs)
    std_tpr = np.std(tprs, axis=0)

    # correct the visuals
    mean_tpr[-1] = 1.0
    for i in range(len(tprs)):
        tprs[i][0] = 0.0
    mean_tpr = np.mean(tprs, axis=0)

    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    
    random_std_tpr = np.std(tprs, axis=0)
    random_mean_tpr = np.mean(random_tprs, axis=0)

    return DataFrame(dict(
        x_linear_space=mean_fpr,
        true_positive_mean=mean_tpr,
        random_expected=random_mean_tpr,
        pooled_auc=mean_auc,
        average_auc=mean(aucs),
        std_auc=std_auc,
        average_auc_ci_lower=ci_cv.rx2('ci')[0],
        average_auc_ci_upper=ci_cv.rx2('ci')[1],
        average_auc_se=ci_cv.rx2('se'),
        random_expected_lower_ci=np.maximum(random_mean_tpr - random_std_tpr, 0),
        random_expected_upper_ci=np.minimum(random_mean_tpr + random_std_tpr, 1),
        true_positive_upper_ci=tprs_upper,
        true_positive_lower_ci=tprs_lower
        #true_positive_lower_ci=np.flip(np.mean(sensitivity_ci_lower, axis=0)),
        #true_positive_upper_ci=np.flip(np.mean(sensitivity_ci_upper, axis=0)),
    ))


def compute_cv_statistics(probabilities, true_responses, boot_n=200):
    """Based on the official plot_roc_crossval example
    which is © 2007 - 2019, scikit-learn developers,
    and distributed under the terms of BSD License.

    scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    Licence: github.com/scikit-learn/scikit-learn/blob/master/COPYING
    """
    tprs = []
    aucs = []
    # proc_aucs = []
    # sensitivity_ci_lower = []
    # sensitivity_ci_upper = []
    mean_fpr = np.linspace(0, 1, max(max(len(t) for t in true_responses), 100))
    random_tprs = []
    random_aucs = []
    # for p, t in tqdm(zip(probabilities, true_responses), total=len(true_responses)):
    for p, t in zip(probabilities, true_responses):

        # roc = roc_object(Series(t), Series(p), quiet=True)

        # this takes too long...
        # sensitivity_ci_estimate = ci_se(
        #    roc, specificities=mean_fpr,
        #    **{'boot.n': boot_n}
        # )

        # TODO: enable parallel backend for ci_se?

        # sensitivity_ci_lower.append(sensitivity_ci_estimate[:,0])
        # sensitivity_ci_upper.append(sensitivity_ci_estimate[:,2])
        # sensitivity_ci_median.append(sensitivity_ci_estimate[:,1])

        fpr, tpr, threshold = metrics.roc_curve(t, p)
        interpolated_tpr = interp(mean_fpr, fpr, tpr)
        tprs.append(interpolated_tpr)
        aucs.append(metrics.auc(fpr, tpr))

        random_fpr, random_tpr, random_threshold = metrics.roc_curve(
            t, np.repeat(1, len(t))
        )

        random_tprs.append(interp(mean_fpr, random_fpr, random_tpr))
        random_tprs[-1][0] = 0.0
        random_aucs.append(metrics.auc(random_fpr, random_tpr))
    return aucs, mean_fpr, random_tprs, tprs
