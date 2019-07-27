is_outlier = function(y) {
    y < quantile(y, 0.25) - IQR(y) * 1.5 | y > quantile(y, 0.75) + IQR(y) * 1.5
}

label_outliers = function(label, value, group, outlier_test=is_outlier) {
    ids = 1:length(label)
    group = list(group)
    outliers = unlist(aggregate(ids, group, identity)$x)[
        unlist(aggregate(value, group, outlier_test)$x)
    ]
    ifelse(ids %in% outliers, as.character(label), NA)
}
