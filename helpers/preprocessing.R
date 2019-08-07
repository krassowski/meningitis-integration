import::here(
    filter_out_low_expression, select_by_expression_level,
    extract_counts, design_from_conditions,
    .from = 'differential_expression.R'
)
import::here(remove_leading_X, .from = 'utilities.R')

voom_cpm = function(matrix, dge) {
    # Based upon limma source code which is licenced GPL (>=2)
    # and copyright of prof. Gordon Smyth

    # NOTE: using edgeR::cpm produces very different results, due
    # to subtle implementation differences:
    #   - scaling of prior.count (log damping factor)
    #   - prior count
    # which can lead to completely different mean-variance trend,
    # see https://support.bioconductor.org/p/59846/ for details.

    lib_size <- dge$samples$lib.size * dge$samples$norm.factors
    t(log2(t(matrix + 0.5) / (lib_size + 1) * 1e6))
}

mean_and_residual_variance = function(y, design) {
    # this is not "blind" for the design!

    # this is identical to the voom method which does not look
    # at stdev but residual stdev, i.e fits linear model for each
    # gene and looks at stdev of residuals
    fit <- limma::lmFit(y, design)
    list(mean=fit$Amean, variance=sqrt(fit$sigma))
}

mean_and_variance = function(y, ...) {
    # "blind" for the design

    # the mean should be identical to fit$Amean above;
    # also, note that this is much faster
    list(
        mean=rowMeans(y),
        variance=sqrt(apply(y, FUN=sd, MARGIN=1))
    )
}

mean_variance_plot = function(stats, color='red', ...) {
    plot(stats$mean, stats$variance, ...)
    tryCatch( {
            lines(lowess(stats$variance ~ stats$mean), col=color)
        },
        error=function(cond){
            print('failed to compute loess (delta problem)')
        }
    )
}

fit_trend = function(stats) {
    loess(stats$variance ~ stats$mean, span=2/3, degree=1, family="symmetric", iterations=3)
}


choose_regions_above_the_mean = function(
    trend, corrected_so_far, only_monotonic=TRUE,
    verbose=FALSE
) {
    distances = c(Inf)
    for (i in 2:length(trend)) {
        distances = c(
            distances,
            abs(trend[i] - mean(trend[1:i]))
        )
    }

    if(only_monotonic) {
        chosen_minimum = min(distances)
        minima = sort(distances)
        candidate_breakpoint = which(distances == chosen_minimum)
        nth = 1
        while(any(trend[1:(candidate_breakpoint-1)] < trend[candidate_breakpoint])) {
            
            nth = nth + 1
            chosen_minimum = minima[nth]
            candidate_breakpoint = which(distances == chosen_minimum)
        }
        if(verbose && nth != 1) {
            print(paste(
                'Skipped', nth, 'non-monotonic breakpoint candidates'
            ))
        }
            
    }
    breakpoint = which(distances == chosen_minimum)

    breakpoint_value = trend[breakpoint]

    if (verbose) {
        print(paste0(
            'Chosen mean breakpoint at ',
            round(100 * (breakpoint / length(trend)), 2),
            '%, for variance=',
            round(breakpoint_value, 3),
            ' with the correction error of ',
            formatC(min(distances), format='e', digits=2)
        ))
    }
    to_correct = !corrected_so_far & (trend > breakpoint_value)

    # once we can no longer select regions above mean, we shrink all that remained
    # (the regions which are very close to the mean).
    if(!any(to_correct) & any(!corrected_so_far)) {
        if (verbose)
            print('No unseen observations above selected breakpoint; correcting remaining observations')
        list(to_correct=!corrected_so_far, swap_to=mean(trend[!corrected_so_far]))
    }
    else
        list(to_correct=to_correct, swap_to=breakpoint_value)
}


choose_all_regions = function(predicted_variance, corrected_so_far) {
    list(
        to_correct=!is.null(predicted_variance),
        swap_to=mean(predicted_variance)
    )
}

loess_trend_correction = function(
    matrix, dge, design, estimate=mean_and_residual_variance,
    epsilon=0.0001, ratio=1, max_iter=20, max_steps=3,
    shirink_relative_to_diff=FALSE,
    plot_iterations=FALSE,
    # given trend estimate, return TRUE for regions that should be corrected
    # and FALSE for regions which do not require correction.
    # by default always return true, but one useful thing to do is
    # to only correct the initial, decreasing part of the trend, or the
    # part of the trend which is above the mean or median, e.g.:
    #   predicted_variance > mean(predicted_variance)
    # one coud also use a condition on the mean abundance, to only correct
    # the lowely expressed genes; There are two predefined function for that:
    # - 'choose_all_regions'
    # - 'choose_regions_above_the_mean'
    choose_regions_to_correct=choose_regions_above_the_mean,
    # only return the results which were produced after all observations were corrected
    # (avoid returning partial results from iterations in which the algorithm did not converge yet)
    require_all_corrected=TRUE
) {
    sh = ifelse(matrix > 0, 1, 1)

    delta = Inf
    result = NULL

    y = voom_cpm(matrix, dge)
    stats = estimate(y, design)
    predicted_variance <- predict(fit_trend(stats), stats$mean)

    last_mean = mean(predicted_variance)
    i = 0
    steps = 0
    corrected_so_far = is.null(predicted_variance)

    stop_condition = function() {
        # if require_all_corrected is set, require that we got at least one result
        if (require_all_corrected) {
            at_least_one_result = !is.null(result)
            if (at_least_one_result)
                # got a result, let's see if we can get a better one
                # within the limitations we got
                delta < epsilon
            else
                # not a single ressult, continue
                FALSE
        }
        else
            delta < epsilon
    }
    
    while(
        # max_iter condition is always applied to prevent too long loops
        !stop_condition() && i < max_iter && steps < max_steps
    ) {

        if(plot_iterations)
            mean_variance_plot(stats, color="yellow", main=paste('Iteration', i))
        
        centered = t(scale(t(matrix), scale=F))

        # sh multiplication is shape broadcasting,
        # not sure how to do that better in R
        regions = choose_regions_to_correct(predicted_variance, corrected_so_far)
        to_correct = sh * regions$to_correct
        corrected_so_far = corrected_so_far | regions$to_correct

        
        max_abs_std_local = max(ifelse(to_correct, abs(stats$variance - predicted_variance), -Inf))
        diff = ifelse(
            to_correct,
            (
                centered
                /
                (median(predicted_variance) / predicted_variance) * ratio
                # preserve the stdev structure
                * (
                    if (shirink_relative_to_diff)
                        1 - (abs(stats$variance - predicted_variance) / (max_abs_std_local))
                    else
                        (1 - abs(stats$variance - predicted_variance) / regions$swap_to)
                )
            ),
            0
        )
        matrix = matrix - diff
        
        y = voom_cpm(matrix, dge)
        stats = estimate(y, design)

        predicted_variance <- predict(fit_trend(stats), stats$mean)

        new_mean = mean(predicted_variance, na.rm=TRUE)

        delta = abs(last_mean - new_mean)
        last_mean = new_mean

        i = i + 1
            
        ## we got full coverage, let's reset the correction mask and store the result
        if(all(corrected_so_far)) {
            corrected_so_far = is.null(predicted_variance)
            result = matrix
            steps = steps + 1
        }
        else if(!require_all_corrected) {
            result = matrix
            steps = steps + 1
        }
    }
            
    if(plot_iterations)
        mean_variance_plot(stats, color="yellow", main=paste('Iteration', i))

    print(paste(
        'Loess mean-variance correction done in',
        steps, 'steps', i, 'iterations'
    ))
    result
}


deseq_transfroms = list(
    vst=DESeq2::vst,
    rlog=DESeq2::rlog,
    ntd=function(dds, blind){
        if(!blind)
            stop('log-norm transform is always blind, pass blind=T to silence this error')
        DESeq2::normTransform(dds)
    }
)
    

correct_trend = function(
    filtered_dge, matrix, dge, design, method,
    to_plot='before-after', blind=TRUE,
    ...
) {
    estimate = ifelse(blind, mean_and_variance, mean_and_residual_variance)

    # full_matrix = dge$counts
    matrix = filtered_dge$counts

    # use voom algorithm for cpm calculation:
    y = voom_cpm(matrix, filtered_dge)

    if (to_plot == 'before-after') {
        par(mfrow=c(1, 2))

        stats = estimate(y, design)
        mean_variance_plot(stats, color="red", main='Before correction')
    }

    stopifnot(method %in% c('voom', 'loess', names(deseq_transfroms)))
    
    if(method %in% names(deseq_transfroms)) {
        transform = deseq_transfroms[[method]]
        coldata = data.frame(group=filtered_dge$samples$group)

        dds = DESeq2::DESeqDataSetFromMatrix(
            countData = round(matrix),
            colData = coldata,
            design = ~ group
        )
        vsd = transform(dds, blind=blind)
        transformed = SummarizedExperiment::assay(vsd)
    }
    if(method == 'voom') {
        # unfortunately, this does not work well :(
        v <- limma::voom(filtered_dge, design, plot=F)

        centered = t(scale(t(matrix), scale=F))

        lib_size <- filtered_dge$samples$lib.size * filtered_dge$samples$norm.factors
        weights = (t(2^(t(v$weights) / 1e6) * (lib_size + 1)) - 0.5)
        
        matrix = matrix - centered * (1-
            (weights - min(weights))
            /
            (max(weights) - min(weights))
        )

        matrix[matrix < 0] = 0

        transformed = matrix
    }
    if(method == 'loess') {
        transformed = loess_trend_correction(
            matrix, filtered_dge, design,
            estimate=estimate, plot_iterations=(to_plot != 'before-after'),
            ...
        )
    }

    if(to_plot == 'before-after') {
        y = voom_cpm(transformed, filtered_dge)
        stats = estimate(y, design)
        mean_variance_plot(stats, color="green", main='After correction')
    }

    transformed
}


calc_normalization_factors = function(dge, method, test_method=NULL, iterations=1, rescale_lib=T, rescale_factor=T) {
    if(is.null(test_method)) {
        
        lib_sizes = colSums(dge$counts)
        stopifnot(all(dge$samples$lib.size == lib_sizes))
        
        if(iterations != 1)
            stop('iterations do not make sense without test_method')
        if(method == 'qtotal') {

            # note from the paper's code it seems that they always scaled the size
            # factors, dividing by the mean -  even for TMM (which is not done in the edgeR::cpm)
            # see https://github.com/wtaoyang/RNASeqComparison/blob/master/Rscript/NormMethods.r
            # this may change things a bit..
            size_factors = ABSSeq::qtotalNormalized(dge$counts)

            # silly mistake or a nice way to demonstrate overfitting:
            # size_factors = size_factors * mean(size_factors)

            if(rescale_factor)
                size_factors = size_factors / mean(size_factors)

            dge$samples$norm.factors = size_factors
            # so this on one hand follows what they did for the paper more closely,
            # but on the other changes the order of magnitude of the size
            # factors by 6-9 orders, thus causing a lot of lowly expressed
            # genes to appear significant in voom - so it seem that this
            # is not the right way to go. Also the paper uses design
            # without intercept which is known to cause problems for voom
            # while I could get great results, I am afraid that this is overfitting.
            # dge$samples$lib.size = lib_sizes / lib_sizes # a vector of 1s
            
            # but this one does not suffer from the above mentioned issues
            dge$samples$lib.size = lib_sizes / lib_sizes * ifelse(rescale_lib, mean(lib_sizes), 1)

        }
        else {
            dge = edgeR::calcNormFactors(dge, method = method)
            
            # unscale the lib-sizes if equested
            dge$samples$lib.size = lib_sizes / ifelse(!rescale_lib, mean(lib_sizes), 1)
        }
    }
    else {
        library(TCC)

        tcc = new("TCC", dge$counts, dge$samples$group)
        tcc = TCC::calcNormFactors(
            tcc,
            norm.method=method,
            test.method=test_method,
            iteration=iterations
        )
        dge$samples$norm.factors = tcc$norm.factors
    }
    dge
}


filter_out_low_expression_by_n = function(matrix, smallest_n, return_frame=TRUE) {
    dummy_condition = c(
        rep(1, smallest_n),
        rep(0, ncol(matrix) - smallest_n)
    )
    dge = edgeR::DGEList(counts=matrix, group=dummy_condition)
    if (return_frame) {
        filtered_dge = filter_out_low_expression(dge, dummy_condition)
        as.data.frame(filtered_dge$counts)
    } else {
        select_by_expression_level(dge, dummy_condition)
    }
}


normalize_abundance = function(
    matrix, by_condition, normalization_method='TMM',
    trend_correction=NULL, log=TRUE, prior.count=2,
    filter=TRUE, subset_rows=FALSE,
    trend_args=list(), 
    ...
) {
    dge = edgeR::DGEList(counts = matrix, group = by_condition)
    
    
    if(filter) {
        stopifnot(subset_rows == FALSE)
        filtered_dge = filter_out_low_expression(dge, by_condition)
    }
    else {
        if(subset_rows != FALSE) {
            filtered_dge = dge[subset_rows,,keep.lib.sizes=FALSE]
        } else {
            filtered_dge = dge
        }
    }
    filtered_dge = calc_normalization_factors(
        filtered_dge, method=normalization_method, ...
    )

    if(!is.null(trend_correction)) {
        design = design_from_conditions(by_condition)

        filtered_dge$counts = do.call(
            correct_trend, 
            c(
                list(filtered_dge, matrix, dge, design),
                trend_args,
                method=trend_correction
            )
        )
        # TODO: should I recompute norm factors now?
    }
    # do not transfrom data corrected by DESeq2, those are already transformed 
    if(
        !is.null(trend_correction) && (!trend_correction %in% deseq_transfroms)
    )
        transformed = filtered_dge$counts
    else
        transformed = edgeR::cpm(filtered_dge, log=log, prior.count=prior.count)


    colnames(transformed) = remove_leading_X(colnames(transformed))
    as.data.frame(transformed)
}