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


# For a plot which has randomly placed elements
# (for example based on the force directed layout)
# iterate over seeds and find the starting point which maximizes
# space usage, assuming that the background is white  
optimize_space_usage = function(plotting_func, end, start=1) {
    library(imager)

    whites_count = 0
    count_whites = function(r,g,b) {
        whites_count <<- sum(r== 1 & r == g & g == b)    
        rgb(r,g,b)
    }

    whites = list()

    for (i in start:end) {
        set.seed(i)

        p = plotting_func()
        suppressMessages(ggsave('.temp.png', p, dpi=80, w=10, h=7, scale=1.2))

        img = load.image('.temp.png')
        img = as.raster(img, colourscale=count_whites, rescale=F)

        whites[i] = whites_count
    }

    whites
}
