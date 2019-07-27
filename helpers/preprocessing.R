import::here(filter_out_low_expression, extract_counts, design_from_conditions, .from = 'differential_expression.R')
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
    lines(lowess(stats$variance ~ stats$mean), col=color)
}

fit_trend = function(stats) {
    loess(stats$variance ~ stats$mean, span=2/3, degree=1, family="symmetric", iterations=3)
}

correct_trend = function(
    filtered_dge, matrix, dge, design, method,
    to_plot='before-after', blind=TRUE
) {
    keep <- edgeR::filterByExpr(dge, design)

    estimate = ifelse(blind, mean_and_variance, mean_and_residual_variance)

    matrix = matrix[keep,]

    # use voom algorithm for cpm calculation:
    y = voom_cpm(matrix, dge)

    if (to_plot == 'before-after') {
        par(mfrow=c(1, 2))

        stats = estimate(y, design)
        mean_variance_plot(stats, color="red", main='Before correction')
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

    if(method %in% names(deseq_transfroms)) {
        transform = deseq_transfroms[[method]]
        coldata = data.frame(group=dge$samples$group)

        dds = DESeq2::DESeqDataSetFromMatrix(
            countData = round(matrix),
            colData = coldata,
            design = ~ group
        )
        vsd = transform(dds, blind=blind)
        transformed = SummarizedExperiment::assay(vsd)
    }
    if(method == 'voom') {
        # unfortunately, this does not work well
        v <- limma::voom(filtered_dge, design, plot=F)

        centered = t(scale(t(matrix), scale=F))

        lib_size <- dge$samples$lib.size * dge$samples$norm.factors
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
        sh = ifelse(matrix > 0, 1, 1)

        delta = Inf
        epsilon = 0.1
        ratio = 1
        max_iter = 10

        # y = voom_cpm(matrix, dge)
        stats = estimate(y, design)
        predicted_variance <- predict(fit_trend(stats), stats$mean)

        last_mean = mean(predicted_variance)
        print(last_mean)
        i = 0

        # optimize_what = 'all'
        # TODO lowly_expressed needs to check monotonicity of the trend!
        # fun = switch(optimize_what, all=min, lowly_expressed=mean)

        while(delta > epsilon && i < max_iter) {

            centered = t(scale(t(matrix), scale=F))

            to_correct = (sh * predicted_variance) > min(predicted_variance)
            max_abs_std_local = max(ifelse(to_correct, abs(stats$variance - predicted_variance), -Inf))
            diff = ifelse(
                to_correct,
                (
                    centered
                    /
                    (median(predicted_variance) / predicted_variance) * ratio
                    # preserve the stdev structure
                    * (
                        1 - (abs(stats$variance - predicted_variance) / (max_abs_std_local))
                    )
                    #* ifelse(predicted_variance > mean(predicted_variance), 1, -1)
                ),
                0
            )
            matrix = matrix - diff

            y = voom_cpm(matrix, dge)
            stats = estimate(y, design)

            if(to_plot != 'before-after')
                mean_variance_plot(stats, color="yellow", main=paste('Step', i))

            predicted_variance <- predict(fit_trend(stats), stats$mean)

            new_mean = mean(predicted_variance, na.rm=TRUE)
            delta = abs(last_mean - new_mean)
            last_mean = new_mean

            i = i + 1
        }

        # TODO: this is missing transformation step!
        transformed = matrix
    }

    if(to_plot == 'before-after') {
        y = voom_cpm(transformed, dge)
        stats = estimate(y, design)
        mean_variance_plot(stats, color="green", main='After correction')
    }

    transformed
}


calc_normalization_factors = function(dge, method, test_method=NULL, iterations=1, rescale_lib=T, rescale_factor=T) {

    if(is.null(test_method)) {
        if(iterations != 1)
            stop('iterations do not make sense without test_method')
        if(method == 'qtotal') {
            lib_sizes = colSums(dge$counts)
            stopifnot(all(dge$samples$lib.size == lib_sizes))

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
            #dge$samples$lib.size = lib_sizes / lib_sizes # a vector of 1s
            
            # but this one does not suffer from the above mentioned issues
            dge$samples$lib.size = lib_sizes / lib_sizes * ifelse(rescale_lib, mean(lib_sizes), 1)
        }
        else {
            dge = edgeR::calcNormFactors(dge, method = method)
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


normalize_abundance = function(
    matrix, by_condition, normalization_method='TMM',
    trend_correction=NULL, ...
) {
    dge = edgeR::DGEList(counts = matrix, group = by_condition)

    filtered_dge = filter_out_low_expression(dge, by_condition)
    dge = calc_normalization_factors(dge, method=normalization_method)
    design = design_from_conditions(by_condition)

    if(!is.null(trend_correction)) {
        filtered_dge$counts = correct_trend(
            filtered_dge, matrix, dge, design, method=trend_correction, ...
        )
    }

    transformed = edgeR::cpm(filtered_dge, log=T, 0.5)
    colnames(transformed) = remove_leading_X(colnames(transformed))

    if(!is.null(file))
        write.csv(transformed, file = out_tmm_normalized_counts_path)
    else
        as.data.frame(transformed)
}


normalize_abundance_for_pls = function(
    matrix, by_condition, file=NULL, ...
) {
    transformed = edgeR::cpm(filtered_dge, log=T, 0.5)
    colnames(transformed) = remove_leading_X(colnames(transformed))

    if(!is.null(file))
        write.csv(transformed, file = out_tmm_normalized_counts_path)
    else
        as.data.frame(transformed)
}
