plot_roc_auc = function(
    roc_auc, random_ci=FALSE, annotate_lines=FALSE,
    annotation_x=0.5, ci=TRUE, mean_auc=FALSE,
    annotation_at=0.05, annotation_nudge_x=NULL,
    label=TRUE, table=FALSE, table_xmin=0.4, table_xmax=0.75,
    color='group', linetype=NULL,
    random_linetype=2, random_color='red', table_row_colors=NULL,
    nudge_x='auto', simple=FALSE, ...
) {
    if (!(table %in% c(FALSE, 'minimal', 'full'))) {
        stop('Unrecognized table value')
    }

    if (is.null(roc_auc$group))
        roc_auc$group = 'AUC'

    grouped = aggregate(roc_auc, list(group_name=roc_auc$group, color=roc_auc[[color]]), mean)
    grouped$group = grouped$group_name
    grouped[[color]] = grouped$color

    p = (
        ggplot(roc_auc, aes(x=x_linear_space))
        + geom_line(aes(y=random_expected), color=random_color, linetype=random_linetype)
    )

    if (random_ci)
        p = p + geom_ribbon(
            aes(ymin=random_expected_lower_ci, ymax=random_expected_upper_ci),
            fill='red', alpha=0.1
        )

    if (ci)
        p = (p
            + geom_ribbon(
                aes_string(
                    ymin='true_positive_lower_ci',
                    ymax='true_positive_upper_ci',
                    fill=color,
                    group='group'
                ),
                alpha=0.1
            )
        )

    x_steps = round(
        nrow(roc_auc) / length(unique(roc_auc$group))
    )

    lowest_fpr = aggregate(
        roc_auc$true_positive_mean,
        list(g=roc_auc$group),
        function(x) x[[ceiling(annotation_at * x_steps)]]
    )
    lowest_fpr = lowest_fpr[match(grouped$group, lowest_fpr$g), ]$x
    annotate_line = if (label) ggrepel::geom_label_repel else ggrepel::geom_text_repel

    if (annotate_lines) {
        if (nudge_x == 'auto') {
            nudge_x = (
                if(!is.null(annotation_nudge_x))
                    annotation_nudge_x
                else
                    min(lowest_fpr) - 0.05
            )
        }
        p = (p
            + annotate_line(
                data=grouped,
                aes_string(
                    x='annotation_x',
                    group='group',
                    y='lowest_fpr',
                    label=(
                        if(simple)
                        "paste0(
                            group
                        )"
                        else
                        "paste0(
                            group,
                            ', AUC: ',
                            round(pooled_auc, 3),
                            ' \U00B1 ',
                            round(std_auc, 2)
                        )"
                    ),
                    color=color
                ),
                hjust=1,
                segment.alpha=0.4,
                nudge_x=nudge_x,
                key_glyph='blank',
                ...
            )
        )
    }


    p = (p
        + geom_line(
            aes_string(
                y='true_positive_mean', group='group',
                color=color, linetype=linetype
            ),
            key_glyph='timeseries'
        )
    )
    if (mean_auc) {
        p = (
            p + annotate(
                'text',
                y=mean(roc_auc$true_positive_mean), x=0.5,
                label=paste('AUC =', round(roc_auc$pooled_auc, 3), '\U00B1', round(roc_auc$std_auc, 2)),
                check_overlap=T
            )
        )
    }

    if (table != FALSE) {
        table_theme = gridExtra::ttheme_minimal(
            core = list(
                fg_params=list(cex = 0.9),
                bg_params=list(fill='grey97')
            ),
            colhead = list(fg_params=list(cex = 0.9)),
            rowhead = list(fg_params=list(cex = 0.9, fontface=1))
        )

        if (table == 'full') {
            table = full_auc_table(grouped, table_row_colors, table_theme)
        } else if (table == 'minimal') {
            table = minimal_auc_table(grouped, table_row_colors, table_theme)
        }

        p = p + annotation_custom(
            table,
            xmin=table_xmin, xmax=table_xmax, ymin=0.05 * nrow(grouped) + 0.07, ymax=0
        )
    }

    p = (p
        + xlab('False positive rate')
        + ylab('True positive rate')
        + coord_fixed()
        + nice_theme
    )
    p
}


add_linetype_to_ggrepel = function(g, linetypes, modify='segment') {
    library(grid)
    g = grid.force(g)

    labels <- grid.ls(
        getGrob(g, c("GRID.labelrepelgrob"), grep=TRUE, global=TRUE),
        print=F
    )$name

    for (i in 1:(length(labels) / 4)) {
        g = editGrob(
            grid.force(g),
            c(labels[[4 * (i - 1) + 1]], modify), grep = TRUE,
            gp = gpar(lty=linetypes[[i]]),
            global=T
        )
    }
    g
}


minimal_auc_table = function(grouped, table_row_colors, table_theme) {

    df = grouped[, c('average_auc'), FALSE]
    df = round(df, 3)
    df$ci = paste(
        round(grouped$average_auc_ci_lower, 3), '-', round(grouped$average_auc_ci_upper, 3)
    )
    colnames(df) = c('AUC', 'AUC 95% CI')
    rownames(df) = grouped$group

    table_grid = gridExtra::tableGrob(df, theme=table_theme)
    table_grid = add_colors_to_rows(table_grid, grouped, table_row_colors)

    header <- gridExtra::tableGrob(
        df[, 1], theme=table_theme,
        cols=c('Averaged')
    )
}


full_auc_table = function(grouped, table_row_colors, table_theme) {

    df = grouped[, c('pooled_auc', 'average_auc')]
    df = round(df, 3)

    df$ci = paste(
        round(grouped$average_auc_ci_lower, 3), '-', round(grouped$average_auc_ci_upper, 3)
    )
    colnames(df) = c('AUC', 'AUC', 'AUC 95% CI')
    rownames(df) = grouped$group

    table_grid = gridExtra::tableGrob(df, theme=table_theme)
    table_grid = add_colors_to_rows(table_grid, grouped, table_row_colors)

    header <- gridExtra::tableGrob(
        df[, 1:2], theme=table_theme,
        cols=c('Pooled', 'Averaged')
    )

    merged_table <- gridExtra::combine(header[1, ], table_grid, along=2)
    merged_table$layout[1:4, c("l", "r")] <- list(c(2, 3), c(2, 4))

    merged_table
}


add_colors_to_rows = function(table_grid, grouped, table_row_colors) {
    if (!is.null(table_row_colors)) {
        # note this very muh does not work due to some kind of bug in tableGrob,
        # a workaround is defined below
        col = table_row_colors[grouped$color]
        col = as.vector(col)
    }
    else {
        col='black'
    }
    if (!is.null(table_row_colors)) {
        i = 2
        for (row_col in col) {
            # row column, for our purpose starts at 2
            ind <- find_cell(table_grid, i, 1, "rowhead-fg")
            table_grid$grobs[ind][[1]][["gp"]] <- grid::gpar(col=row_col, cex=0.9)
            i = i + 1
        }
    }
    table_grid
}


# function from documentation of tableGrob
# which is authored by Baptiste Auguie and GPL (≥ 2) licences
# https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html 
find_cell = function(table, row, col, name="core-fg") {
    l <- table$layout
    which(l$t==row & l$l==col & l$name==name)
}