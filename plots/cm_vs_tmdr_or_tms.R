import::here(plot_roc_auc, .from='roc_auc.R')

default_colors = c('Test: TMD+TMR vs CM'='#E69F00', 'Test: TMS vs CM'='#56B4E9', 'Train: TMD+TMR vs CM'='#999999')


roc_auc_tmdr_cm = function(roc_auc, colors=default_colors, annotation=0.05, table_xmin=0.45, table_xmax=0.75, ...) {
    (
        plot_roc_auc(
            roc_auc, annotate_lines=T, table='full',
            color='set', linetype='models',
            random_color='grey50', random_linetype=3,
            table_row_colors=colors, annotation_at=annotation,
            annotation_x=annotation, nudge_x=1,
            nudge_y=-1,
            direction='y',
            ylim=c(0.35, 0.65),
            table_xmin=table_xmin, table_xmax=table_xmax, ...
        )
        + scale_color_manual(name='Set', values=colors)
        + scale_fill_manual(name='Set', values=colors)
        + scale_linetype_manual(
            name='Tested on',
            values=c('Cross-Validation'='dashed', 'Full model'='solid')
        )
        + guides(
            linetype=guide_legend(title.position='top'),
            fill=guide_legend(title.position='top'),
            color=guide_legend(title.position='top')
        )
    )
}