signif_thresholds = c(0, 0.001, 0.01, 0.05, 0.1, 1)
signif_codes = c("\U2042", "\U2051", "\U204E", "\U26AC", " ")

convert_signif_codes = function(pvalues) {
    # extracted from base R, see https://stackoverflow.com/a/41263378/6646912
    symnum(
        pvalues, corr=FALSE, na=FALSE,
        cutpoints=signif_thresholds,
        symbols=signif_codes
    )
}