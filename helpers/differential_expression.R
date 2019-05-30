import::here(space_to_dot, dot_to_space, .from='utilities.R')


design_from_conditions = function(conditions, intercept=F) {

    # "The levels must by syntactically valid names in R"
    conditions = space_to_dot(conditions)

    groups = as.factor(conditions)
    if (intercept != F) {
        groups = relevel(groups, ref=intercept)
        design = model.matrix(~ groups)
    } else {
        design = model.matrix(~ 0 + groups)
    }
    colnames(design) = gsub('groups', '', colnames(design))
    design
}

calculate_means = function(data, design) {
    fit = limma::lmFit(data, design)
    fit
}

limma_fit = function(data, conditions_vector, a, b, use_all=T, workaround_for_non_transformed_data=F) {

    if(workaround_for_non_transformed_data == T && use_all != T)
        stop('the workaround is only supported when using the full data')

    if(use_all) {
        # in this scenario I benefit from the additional information about the distribution of protein abundances from other conditions;
        # this might or might not be desirable, however https://support.bioconductor.org/p/73107/ indicates that this is usually beneficial, unless:
        # """
        # The most obvious motivation for partitioning would be if some groups had unequal variances,
        # such that you wouldn't want to analyze them together for fear of getting an inappropriate variance estimate.
        # However, this effect can probably be ignored -see https://support.bioconductor.org/p/67984/#67987 for a good explanation.
        # There may well be other scenarios where partitioning is appropriate... [by Aaron Lun]
        # """
        # I do not have the experience to decide if in our case the variances are to unequal between condition
        # (although my intuition is that yes, they are quite diferent as we analyze quite an extereme situtation
        # where certain cells come up in the CSF fluid where there are usually none or little of them), sp
        # I will just look at the results of both scenarios (and conclude that further reading is needed if the
        # results differ or that the difference is negligible if the difference is negligible) 
        
        design <- design_from_conditions(conditions_vector)
        fit <- calculate_means(data, design)

        contrast_specification <- paste(
            space_to_dot(a),
            space_to_dot(b),
            sep='-'
        )
        contrast.matrix <- limma::makeContrasts(contrasts=contrast_specification, levels=design)
        #fit <- limma::contrasts.fit(fit, contrast.matrix)
        fit = contrasts_fit(fit, contrast.matrix, divide=workaround_for_non_transformed_data)
    } else {
        considered <- cbind(
            data[conditions_vector == a],
            data[conditions_vector == b]
        )
        a_cnt <- sum(conditions_vector == a)
        b_cnt <- sum(conditions_vector == b)
        groups <- c(rep(1, a_cnt), rep(0, b_cnt))
        design <- cbind(Intercept=1, Group=groups)

        fit <- limma::lmFit(considered, design)
    }
    fit
}


limma_diff_ebayes <- function(a, b, ...){
    fit = limma_fit(a=a, b=b, ...)
    limma::eBayes(fit, trend=T, robust=T)
}

full_table <- function(e_bayes_result, coef=1){    
    table = limma::topTable(e_bayes_result, coef=coef, number=Inf)
    table$protein = rownames(table)
    table
}


contrasts_fit = function(fit, contrast.matrix, divide, ...) {
    if (divide == F) {
        diff_fit <- limma::contrasts.fit(fit, contrast.matrix)
    } else {
        print("
            see notes/Limma_expects_log_transformed_data.ipynb
            NOTE: this is just to indicate how using non log-transformed data works
            and that it is wrong; the results of divide=T were not validated
            and are not expected to be full correct but rather to be used as
            a demenstration of the log-transformation impact and where the issue lies.
        ")
        diff_fit <- limma::contrasts.fit(fit, contrast.matrix)
        a = rownames(contrast.matrix[contrast.matrix == 1,, drop=F])
        b = rownames(contrast.matrix[contrast.matrix == -1,, drop=F])
        diff_fit$coefficients = fit$coefficients[,a] / fit$coefficients[,b]
    }
    diff_fit
}