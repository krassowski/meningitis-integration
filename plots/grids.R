library(ggplot2)


ggplot_columns_grid = function(data, ggplot, ncol=8, prepare_legend=function(plot){plot}) {
    first_column_name = colnames(data)[1]
    plot = ggplot(data, first_column_name)
    legend = cowplot::get_legend(prepare_legend(plot + theme(legend.position='top')))

    plots = list()

    for (column_name in colnames(data)) {
        column = data[[column_name]]
        if(is.numeric(column))
        {
            data$val = column_name

            plots[[length(plots) + 1]] = (
                ggplot(data, column_name) +
                facet_wrap('val') + theme(legend.position='none', axis.title.x = element_blank(), axis.title.y = element_blank())

            )
        }
    }

    gridExtra::grid.arrange(
        legend,
        gridExtra::arrangeGrob(grobs=plots, ncol=ncol),
        heights=c(1, 10), nrow=2
    )
}
