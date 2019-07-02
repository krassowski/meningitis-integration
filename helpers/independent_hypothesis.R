import::here(deseq2_prevent_aggressive_filtering, .from='thirdparty.R')
library('genefilter')

# Note: Deseq2 uses alpha=0.1 by default
filter_by_mean_expression = function(table, fdr_threshold=0.05, step=0.005, expression_col='AveExpr', pvalue_col='P.Value', show_plot=T, use_deseq2_method=F) {
    import::from(genefilter, filtered_p)
    
    theta = seq(from=0.1, to=1, by=step)
    
    if(is.character(expression_col))
        col = table[,expression_col]
    else
        col = expression_col
    filtPadj <- filtered_p(filter=col, test=table[,pvalue_col], theta=theta, method="BH")
    sigGenes <- colSums(filtPadj < fdr_threshold, na.rm = TRUE)

    if(show_plot)
        plot(theta, sigGenes, type="b", xlab="Quantile filtered out", ylab="Significant genes")
    else {
        if(use_deseq2_method)
            j = deseq2_prevent_aggressive_filtering(sigGenes, theta)
        else
            j = which.max(sigGenes)
        filtPadj[,j]
    }
}

add_p_value_for_filtered_subset = function(data, ...) {
    data$unfiltered.adj.P.Val = data$adj.P.Val
    data$adj.P.Val = filter_by_mean_expression(data, ..., show_plot=F)
    data
}


library("IHW")
# https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-6-hypothesis-weighting/introduction_to_ihw.html
# https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html

weight_by_mean_expression = function(table, fdr_threshold=0.05, expression_col='AveExpr', pvalue_col='P.Value') {
    formula = as.formula(paste(pvalue_col, '~', expression_col))
    ihw_res <- IHW::ihw(formula, data=table, alpha=fdr_threshold)
    ihw_res
}

add_p_value_for_weighted_hypothesis = function(data, ...) {
    data$unweighted.adj.P.Val = data$adj.P.Val
    data$adj.P.Val = IHW::adj_pvalues(weight_by_mean_expression(data, ...))
    data
}
