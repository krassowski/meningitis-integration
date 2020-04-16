import::here(space_to_dot, dot_to_space, .from='utilities.R')


design_from_conditions = function(conditions, intercept=T) {

    # "The levels must by syntactically valid names in R"
    conditions = space_to_dot(conditions)

    groups = as.factor(conditions)
    if (intercept != F) {
        groups = relevel(groups, ref=intercept)
        design = model.matrix(~ groups)
    } else {
        # note: use of no-intercept mode matrix was causing issues for voom, see:
        # https://support.bioconductor.org/p/57631/
        design = model.matrix(~ groups)
    }
    colnames(design) = gsub('groups', '', colnames(design))
    design
}

calculate_means = function(data, design) {
    fit = limma::lmFit(data, design)
    fit
}

limma_fit = function(data, conditions_vector, a, b, use_all=T, workaround_for_non_transformed_data=F, intercept=T) {

    if (workaround_for_non_transformed_data == T && use_all != T)
        stop('the workaround is only supported when using the full data')

    if (use_all) {
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
        
        design <- design_from_conditions(conditions_vector, intercept=intercept)
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
        # TODO test change below against protein diff expr notebook

        considered <- cbind(
            data[,conditions_vector == a],
            data[,conditions_vector == b]
        )
        a_cnt <- sum(conditions_vector == a)
        b_cnt <- sum(conditions_vector == b)
        groups <- c(rep(1, a_cnt), rep(-1, b_cnt))
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


extract_counts = function(counts) {
    # TODO: warn about unmatched class?
    #  or maybe rename to "extract_counts_if_needed"
    if (class(counts) == 'DESeqDataSet')
        counts = as.data.frame(DESeq2::counts(counts, normalized=TRUE))

    if (class(counts) == 'DGEList')
        counts = edgeR::cpm(counts)

    counts
}


select_by_expression_level = function(dge, conditions_vector) {
    design <- design_from_conditions(conditions_vector)
    counts_condition <- length(rownames(dge))
    keep <- edgeR::filterByExpr(dge, design)
    ratio <- sum(keep) / length(rownames(dge)) * 100
    print(paste0('Retaining: ', round(ratio, 2), '%'))
    keep
}


filter_out_low_expression = function(dge, conditions_vector) {
    keep = select_by_expression_level(dge, conditions_vector)
    dge[keep,,keep.lib.sizes=FALSE]
}


# Move from RNA differential expression notebook:

import::here(space_to_dot, .from='utilities.R')
import::here(replace_ids, .from='identifiers_mapping.R')


choose_ref = function(a, b, conditions) {
    space_to_dot(setdiff(unique(conditions), c(a, b))[1])
}

calc_dge = function(a, b, data, conditions, voom=F, quality_weights=F, ref=F, as_is=F) {
    design <- design_from_conditions(conditions, intercept=ref)

    if (voom) {
        if (quality_weights)
            # see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4551905/
            # NB: it might possible to skip removal of the outliers when using "quality_weights"
            dge = limma::voomWithQualityWeights(data, design, normalize.method="none")
        else
            # TODO: normalize="quantile"
            dge = limma::voom(data, design)
    }
    else if (as_is) {
        dge = data
    }
    else {
        # expecting data to be fitered counts with pre-computed normalization factors
        dge <- edgeR::cpm(data, log=TRUE, prior.count=0.25)
    }
    dge
}

calc_fit = function(a, b, data, conditions, voom=F, robust=T, use_all=T, ...) {
    ref = choose_ref(a, b, conditions)
    if (is.na(ref)) {
        ref = F
        use_all = F
        print('No reference group, switching use_all off')
    }

    dge = calc_dge(a, b, data, conditions, voom=voom, ref=ref, ...)

    raw_fit <- limma_fit(data=dge, conditions_vector=conditions, a=a, b=b, intercept=ref, use_all=use_all)

    # when using voom, the trend correction is redundant
    limma::eBayes(raw_fit, trend=!voom, robust=robust)
}


calc_de = function(a, b, data, conditions, confint=F, ...) {
    fit = calc_fit(a, b, data, conditions, ...)
    limma::topTable(fit, n=Inf, confint=confint)
}


camera_with_statistic = function(table, statistic, collection, ...) {
    statistic = table[,statistic]
    statistic = replace_ids(statistic, table, ...)
    limma::cameraPR(statistic, collection)
}


calc_camera = function(a, b, data, conditions, statistic=NULL, collection=reactome_new_symbol, convert_to='symbol', ...) {
    set.seed(0)
    if(!is.null(statistic)) {
        table = calc_de(a, b, data, conditions, ...)

        camera_with_statistic(table, statistic, collection=collection, convert_to=convert_to)

        # fit = calc_fit(a, b, data, conditions, ...)
        # rownames(fit$t) = ensembl_to_entrez[rownames(fit),]
        # limma::cameraPR(fit$t[,1], collection)
    }
    else {
        ref = choose_ref(a, b, conditions)
        if(is.na(ref)) {
            ref = F
            print('No reference group, setting to FALSE')
        }

        dge = calc_dge(a, b, data, conditions, ref=ref, ...)
        design <- design_from_conditions(conditions, intercept=ref)

        if(ref != F){
            contrast_specification <- paste(
                space_to_dot(a),
                space_to_dot(b),
                sep='-'
            )
            contrast.matrix <- limma::makeContrasts(contrasts=contrast_specification, levels=design)
        }
        if(!is.null(convert_to))
            dge = replace_ids(dge, dge, convert_to=convert_to)
        if (ref != F)
            limma::camera(dge, collection, design, contrast=contrast.matrix)
        else
            limma::camera(dge, collection, design)
    }
}


get_deseq_definite_tb_cm = function(dds, shrinkage=NULL, ref=NULL, ...) {
    if (is.null(shrinkage)) {
        result = DESeq2::results(
            dds,
            contrast=c('conditions_for_deseq', 'Definite.tuberculosis', 'Cryptococcal'),
            parallel=T,
            ...
        )
    } else {
        dds$conditions_for_deseq = relevel(dds$conditions_for_deseq, 'Cryptococcal')
        dds = nbinomWaldTest(dds)

        result = lfcShrink(
            dds,
            coef='conditions_for_deseq_Definite.tuberculosis_vs_Cryptococcal',
            type=shrinkage,
            parallel=T,
            ...
        )
    }
    as.data.frame(result)
}