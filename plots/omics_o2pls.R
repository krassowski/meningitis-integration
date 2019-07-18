# Note: many functions in this file are GPL-3,
# due to an extraction function data_for_plot which is based on GPL-3 derived code

data_for_plot = function (
    x, loading_name = c("Xjoint", "Yjoint", "Xorth", "Yorth"), 
    i = 1, j = NULL, label = c("number", "colnames"),
    ...
) {
    # Based on plot.o2p by Said el Bouhaddani et al., licence: GPL-3
    # Note that the plotting functions usign this function are also GPL-3
    # due to the licence virulence
    loading_name = match.arg(loading_name)
    which_load = switch(
        loading_name,
        Xjoint = "W.", Yjoint = "C.", 
        Xorth = "P_Yosc.", Yorth = "P_Xosc."
    )
    load = as.matrix(x[which_load][[1]])
    if (ncol(load) < max(i, j)) 
        stop("i and j cannot exceed #components = ", ncol(load))
    load = load[, c(i, j)]
    p = nrow(as.matrix(load))
    if (is.null(j)) {
        load = cbind(1:p, load)
        colnames(load) = c("index", paste(loading_name, "loadings", i))
    }
    else {
        stopifnot(j == round(j))
        colnames(load) = c(
            paste(loading_name, "loadings", i),
            paste(loading_name, "loadings", j)
        )
    }
    label = match.arg(label)
    if (label == "colnames" && !is.null(rownames(x[which_load][[1]]))) {
        label = rownames(x[which_load][[1]])
    }
    load = as.data.frame(load)
    load$label = rownames(load)
    load
}

plot_most_extreme_loadings = function(fit, n=25, loading_name='Xjoint', i=1, x_axis='index', x_values=NA) {
    data = data_for_plot(fit, i=i, loading_name=loading_name)
    loading = paste(loading_name, 'loadings', i) 
    if(!is.na(x_values))
        data[[x_axis]] = x_values
    most_extreme = data[head(order(-abs(data[[loading]])), n=n),]
    
    p = (
        ggplot(data, aes_string(y=paste0('`', loading, '`'), x=x_axis))
        + geom_point(alpha=0.5)
        + ggrepel::geom_label_repel(data=most_extreme, aes(label=label))
    )
    p
}

get_scores = function(fit, matrix) {
    id = switch(matrix, X='Tt', Y='U')
    fit[[id]]
}

matrix_scores_plot = function(fit, matrix, x=1, y=2, color='black', ...) {
    name = switch(matrix, X='T', Y='U')
    m = get_scores(fit, matrix)
    if(length(colnames(m)) < 2) {
        y_data = 1
        y_lab = '-'
        y = '-'
    } else {
        y_data = m[,y]
        y_lab = paste0(matrix, '-Score ', name, y)
        y_aes = 'y'
    }
    data = data.frame(x=m[,x], y=y_data, label=rownames(m), color=color)
    (
        ggplot(data, aes_string(x='x', y='y', label='label', color='color'))
        + geom_point()
        + ggrepel::geom_label_repel(...)
        + xlab(paste0(matrix, '-Score ', name, x))
        + ylab(y_lab)
        + stat_ellipse()
    )
}


single_scores_bar = function(fit, matrix, color='black') {
    (
        matrix_scores_plot(fit, matrix, color=color, nudge_y=0.15, hjust=0, direction='y', y=1)
        + nice_theme + coord_flip() + theme(legend.position='none') + ylim(1, 1.35)
        + theme(
            axis.line.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank()
          )
        + ylab(paste0(matrix, '-Score'))
    )
}