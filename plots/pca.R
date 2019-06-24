library(ggplot2)
source('plots/colors.R')

small_pca = function(matrix, name='', axes=c(1, 2)) {
    colnames(matrix) <- sub('^X', '', colnames(matrix))
    pca = prcomp(t(matrix), scale=T)
    (
        fviz_pca_biplot(
            axes=axes,
            pca, repel=T,
            geom.ind="point",
            gradient.cols=c("#c49c94", "#2CA02C"),
            select.var=list(contrib=25),
            col.ind='white',
            # TODO: generalise
            fill.ind=patient_annotations$Meningitis,
            palette=patient_colors$Meningitis,
            pointshape=21, pointsize=3,
            #addEllipses=T,
            #alpha.var="contrib",
            col.var="contrib",
            alpha.ind=.3,
            alpha.var=.9,
            title=paste('Biplot of', axes[1], 'and', axes[2], 'PCs of', name)
        )
        + labs(fill='Meningitis', color='Contribution')
        + theme(legend.position='bottom')
        + guides(color=guide_colourbar(barheight=0.5, title.vjust=1))
    )
}

first_four_pcs = function(matrix, name) {
    # TODO: rewrite
    a = small_pca(matrix, axes=c(1, 2), name=name)
    b = small_pca(matrix, axes=c(3, 4), name=name)
    cowplot::plot_grid(
        cowplot::plot_grid(
            a + theme(legend.position = "none"),
            b + theme(legend.position = "none"),
            labels = "AUTO"
        ),
        cowplot::plot_grid(
            cowplot::get_legend(a + guides(fill=F)),
            cowplot::get_legend(b + guides(fill=F))
        ),
        cowplot::get_legend(a + guides(color=F)),
        nrow=3,
        rel_heights=c(1, .08, .06)
    )
}

cumulative_variance_explained <- function(data, threshold){
    # note: preferred way to select number of PC would be to use CV,
    # though this gives a nice overwiew of the PCs importance anyway
    (
        ggplot(data, aes(x=PC, y=cumulative.variance.percent))
        + geom_bar(stat='identity', fill='grey70')
        + annotate('segment', x=0, xend=Inf, y=threshold, yend=threshold, colour="red")
        + theme(legend.position='none', axis.text.x=element_text(angle=90))
        + xlab('Prinicpal component') + ylab('Cumulative variance explained [%]')
    )
}


plot_pca_with = function(pca, key, func=fviz_pca_ind, data=normalized_form, ...){
    # TODO: rewrite using grids::ggplot_columns_grid
    plot = fviz_pca_ind(pca, geom='point', col.ind=data[[key]])
    legend = cowplot::get_legend(plot + theme(legend.position='top') + guides(color=guide_legend(title=key)))
    plots = list()
    for(i in seq(1, 9))
       plots[[i]] = (
           # TODO: geom?
           func(pca, geom='point', axes=c(i, i+1), col.ind=data[[key]], ...)
           + theme(legend.position='none')
       )
    gridExtra::grid.arrange(
        legend,
        gridExtra::arrangeGrob(grobs = plots, ncol=3),
        heights=c(1, 10), nrow=2
    )
}
