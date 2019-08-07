from functools import partial
from pandas import DataFrame
import numpy as np
from scipy import interp
from sklearn import metrics
from helpers.r import r_function


test_roc = partial(r_function, 'roc.test')
auc_ci = partial(r_function, 'ci')
ci_se = partial(r_function, 'ci.se')
roc_object = partial(r_function, 'roc')
power_roc_test = partial(r_function, 'power.roc.test')


def compare_roc_curves(
    a_response, a_prediction,
    b_response, b_prediction,
    paired=False,
    compute_power=False,
    power=0.8
):
    roc1 = roc_object(a_response, a_prediction, quiet=True)
    roc2 = roc_object(b_response, b_prediction, quiet=True)

    p = test_roc(roc1, roc2, paired=paired).rx2('p.value')[0]

    if compute_power:
        assert paired
        # get power at significance level 0.05 for two-sided alternative
        power = power_roc_test(roc1, roc2)
        # get number of patients required to get a=0.05 and B=0.8
        n_at_power = power_roc_test(roc1, roc2, power=power)

        power_report = {
            f'cases_needed_for_{power}': n_at_power.rx2('ncases')[0],
            f'controls_needed_for_{power}': n_at_power.rx2('ncontrols')[0],
            **{
                k: v[0]
                for k, v in power.items()
            }
        }
    else:
        power_report = None

    return p, power_report


def roc_auc_plot_data(probabilities, true_responses):
    """Based on the official plot_roc_crossval example
    which is © 2007 - 2019, scikit-learn developers,
    and distributed under the terms of BSD License.
    
    scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    Licence: github.com/scikit-learn/scikit-learn/blob/master/COPYING
    """

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, len(true_responses))

    random_tprs = []
    random_aucs = []
    
    for p, t in zip(probabilities, true_responses):
        fpr, tpr, threshold = metrics.roc_curve(
            t, p
        )
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        aucs.append(metrics.auc(fpr, tpr))
        
        # random
        #np.mean(t)
        #class_impalance
        fpr, tpr, threshold = metrics.roc_curve(
            t, np.repeat(1, len(t))
        )
        random_tprs.append(interp(mean_fpr, fpr, tpr))
        random_tprs[-1][0] = 0.0
        random_aucs.append(metrics.auc(fpr, tpr))
        
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    std_tpr = np.std(tprs, axis=0)

    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    
    random_std_tpr = np.std(tprs, axis=0)
    random_mean_tpr = np.mean(random_tprs, axis=0)
    
    return DataFrame(dict(
        x_linear_space=mean_fpr,
        true_positive_mean=mean_tpr,
        random_expected=random_mean_tpr,
        mean_auc=mean_auc,
        auc=aucs,
        std_auc=std_auc,
        random_expected_lower_ci=np.maximum(random_mean_tpr - random_std_tpr, 0),
        random_expected_upper_ci=np.minimum(random_mean_tpr + random_std_tpr, 1),
        true_positive_upper_ci=tprs_upper,
        true_positive_lower_ci=tprs_lower
    ))