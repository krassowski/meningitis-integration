strip_strata_prefix = function(fit) {
    gsub('.*=', '', names(fit$strata))
}
