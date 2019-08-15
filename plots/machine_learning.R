import::here(convert_signif_codes, signif_codes, signif_thresholds, .from='significance_codes.R')
library(scales)


select_coeffs = function(coeffs, n, by='p_value', descending=FALSE, abs_value=FALSE) {
    transform = ifelse(abs_value, abs, identity)
    head(
        coeffs[
            order(ifelse(descending, 1, -1) * transform(coeffs[[by]])), 
        ],
        n
    )
}


plot_coefficients = function(
    coeffs, limit_to_n_most_extreme=NULL, label_margin=0.05,
    format='e', digits=0,
    select='p_value', label='p_value', fill='selected_in',
    select_descending=TRUE, select_abs=FALSE,
    label_name=NULL, fill_name=NULL
) {
    if(!is.null(limit_to_n_most_extreme)) {
        coeffs = select_coeffs(
            coeffs, n=limit_to_n_most_extreme,
            by=select, descending=select_descending,
            abs_value=select_abs
        )
    }
    p = (
        ggplot(coeffs, aes(x=gene, y=mean))
        + geom_bar(aes_string(fill=fill), stat='identity')
        + geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci))
    )
    
    if (label != 'fdr') {
        p = (p
            + geom_text(aes(
                label=formatC(coeffs[[label]], format=format, digits=digits),
                color=select,
                y=ifelse(mean>0, mean + ci + label_margin, mean - ci - label_margin),
                vjust=ifelse(mean>0, 0, 1)
            ))
            + scale_color_manual(
                name='', values='black',
                labels=ifelse(is.null(label_name), label, label_name)
            )
        )
    } else {
        names(signif_codes) = signif_codes
        names(signif_thresholds) = c('', signif_codes)
        p = (p
            + geom_point(aes(
                    y=ifelse(mean>0, mean + ci + label_margin, mean - ci - label_margin),
                    shape=as.character(convert_signif_codes(coeffs[[label]]))
                ),
                size=5, color='red'
            )
            + scale_shape_manual(
                name=paste0('Significance level (', 'FDR', ')'),
                values=signif_codes, breaks=signif_codes,
                labels=signif_thresholds[signif_codes]
            )
        )
    }
    
    (p
        + nice_theme
        + theme(axis.text.x=element_text(angle=90, hjust=1))
       
        + ylab('Mean contribution (coefficient)')
        #+ ggforce::facet_zoom(x=value < -0.1)
    )
}



plot_highest_magnitutude = function(coeffs, n, fill_name='Ratio of models including this feature', ...) {
    (
        plot_coefficients(
            coeffs, limit_to_n_most_extreme=n, select='mean',
            label='fdr', label_name='FDR', select_abs=TRUE, select_descending=F,,
            fill='selected_in', fill_name=fill_name, ...
        )
        + xlab('Feature')
        + ggtitle('Features with the highest magnitude of contributions to the logistic regression models')
    )
}


plot_most_significant = function(coeffs, n, fill_name='Ratio of models including this feature', ...) {
    fill_scale=c(-10, -1, 0, 1, 10)
    (
        plot_coefficients(
            coeffs, limit_to_n_most_extreme=n, select='p_value',
            label='fdr', label_name='FDR',
            fill='selected_in', ...
        )
        + scale_fill_gradientn(
            name=ifelse(is.null(fill_name), fill, fill_name),
            #midpoint=0.5,
            colors=RColorBrewer::brewer.pal(n = 5, name = "RdBu"),
            values = scales::rescale(fill_scale)
        )
        + xlab('Feature')
        + ggtitle('Features with the most significant contributions to the logistic regression models')
    )
}


plots_most_frequently_included = function(coeffs, n=20, label_name='Ratio of models including this feature', ...) {
    (
        plot_coefficients(
            coeffs, limit_to_n_most_extreme=n, format='f', digits=2,
            select_descending=F,
            select='selected_in', label='selected_in', label_name=label_name,
            fill='-log10(p_value)', fill_name='P-value',
            ...
        )
        + scale_fill_gradient(
            #name=ifelse(is.null(fill_name), fill, fill_name),
        )
    )
}

how_volatile = function(volatile) {
    ifelse(
        volatile > 0.01,
        ifelse(
            volatile > 0.1,
            'volatile',
            'unstable'
        ),
        'stable'
    )

}

how_frequent = function(selected_in) {
    ifelse(
        selected_in > 0.5,
        ifelse(selected_in > 0.75, '>0.75', '>0.50'),
        '\U2264 0.50'
    )
}

frequency_fill = list(
    '\U2264 0.50'='grey',
    '>0.50'='orange',
    '>0.75'='red'
)

volatility_shape = list(
    stable=21,
    volatile=1,
    unstable=4
)



mean_vs_coefficients = function(coeffs, n=10, transformation=identity, fdr_threshold=0.01) {
    coeffs$is_frequent = coeffs$selected_in > 0.5
    top_coeffs = select_coeffs(coeffs, n, 'mean', descending=0, abs_value=1)
    
    (
        ggplot(coeffs, aes(x=transformation(mean_abundance), y=mean))
        + geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), color='grey40')
        + geom_point(
            aes(
                fill=how_frequent(selected_in),
                shape=how_volatile(volatile)
            ),
            size=2.5, alpha=0.8, color=ifelse(coeffs$is_frequent, 'grey10', 'grey30')
        )
        + scale_shape_manual(name='sign', values=volatility_shape)
        + ggrepel::geom_label_repel(
            data=top_coeffs,
            aes(label=gene, y=mean, color=ifelse(is.na(fdr), 'black', fdr < fdr_threshold)),
            box.padding=0.5,
            segment.color='grey10'
        )
        + nice_theme
        + scale_fill_manual(name='Selected in', values=frequency_fill)
        + scale_color_manual(name=paste('FDR < ', fdr_threshold), values=c('grey10', 'darkblue'))
        # see bug https://github.com/tidyverse/ggplot2/issues/2322
        + guides(
            fill=guide_legend(override.aes=list(shape=21)),
            shape=guide_legend(override.aes=list(fill='brown'))
        )
    )
}

mean_vs_coefficients = function(coeffs, n=10, transformation=identity, fdr_threshold=0.01) {
    coeffs$is_frequent = coeffs$selected_in > 0.5
    top_coeffs = select_coeffs(coeffs, n, 'mean', descending=0, abs_value=1)
    
    p = (
        ggplot(coeffs, aes(x=transformation(mean_abundance), y=mean))
        + geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), color='grey40')
        + geom_point(
            aes(
                fill=how_frequent(selected_in),
                shape=how_volatile(volatile)
            ),
            size=2.5, alpha=0.8, color=ifelse(coeffs$is_frequent, 'grey10', 'grey30')
        )
        + scale_shape_manual(name='sign', values=volatility_shape)
        + ggrepel::geom_label_repel(
            data=top_coeffs,
            aes(label=gene, y=mean, color=ifelse(is.na(fdr_threshold), 'none', fdr < fdr_threshold)),
            box.padding=0.5,
            segment.color='grey10'
        )
        + nice_theme
        + scale_fill_manual(name='Selected in', values=frequency_fill)
    )
    if(!is.na(fdr_threshold))
        p = p + scale_color_manual(name=paste('FDR < ', fdr_threshold), values=c('grey10', 'darkblue'))
    else
        p = p + scale_color_manual(name='', values=c('grey10'), labels='gene')
       # see bug https://github.com/tidyverse/ggplot2/issues/2322
    p = (p + guides(
            fill=guide_legend(override.aes=list(shape=21)),
            shape=guide_legend(override.aes=list(fill='brown'))
        
        )
    )
    p
}


coefficients_volcano_plot = function(coeffs, n=10, transformation=identity, p_value='p_value', neg=TRUE) {
    coeffs$is_frequent = coeffs$selected_in > 0.5
    coeffs$gene = rownames(coeffs)
    top_coeffs = select_coeffs(coeffs, n, by=p_value, descending=neg)
    
    (
        ggplot(coeffs, aes_string(y='mean', x=paste0((if (neg) '-' else ''), 'log10(', p_value, ')')))
        + geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), color='grey40')
        + geom_point(
            aes(
                fill=how_frequent(selected_in),
                shape=how_volatile(volatile)
            ),
            size=2.5, alpha=0.8, color=ifelse(coeffs$is_frequent, 'grey10', 'grey30')
        )
        + scale_shape_manual(name='sign', values=volatility_shape)
        + ggrepel::geom_text_repel(
            data=top_coeffs,
            aes(label=gene),
            #box.padding=0.5,
            segment.color='grey10'
        )
        + coord_flip()
        + nice_theme
        + scale_color_manual(values=c('grey', 'red'))
        + scale_fill_manual(name='Selected in', values=frequency_fill)
        + guides(
            fill=guide_legend(override.aes=list(shape=21)),
            shape=guide_legend(override.aes=list(fill='brown'))
        
        )
    )
}

                                                     
# function from documentation of tableGrob
# which is authored by Baptiste Auguie and GPL (≥ 2) licences
# https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html 
find_cell <- function(table, row, col, name="core-fg") {
    l <- table$layout
    which(l$t==row & l$l==col & l$name==name)
}


plot_roc_auc = function(
    roc_auc, random_ci=FALSE, annotate_lines=FALSE,
    annotation_x=0.5, ci=TRUE, mean_auc=FALSE,
    annotation_at=0.05, annotation_nudge_x=NULL,
    label=TRUE, table=FALSE, table_xmin=0.4, table_xmax=0.75,
    color='group', linetype=NULL,
    random_linetype=2, random_color='red', table_row_colors=NULL,
    nudge_x='auto', ...
) {
    
    if(is.null(roc_auc$group))
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
        function(x) {x[[ceiling(annotation_at * x_steps)]]}
    )
    lowest_fpr = lowest_fpr[match(grouped$group, lowest_fpr$g),]$x
    annotate_line = if (label) ggrepel::geom_label_repel else ggrepel::geom_text_repel

    if (annotate_lines) {
        if(nudge_x == 'auto') {
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
                    label="paste0(
                        group,
                        ', AUC: ',
                        round(pooled_auc, 3),
                        ' \U00B1 ',
                        round(std_auc, 2)
                    )",
                    color=color
                ),
                #direction='x',
                hjust=1,
                segment.alpha=0.4,
                nudge_x=nudge_x,
                key_glyph='blank',
                ...
                #nudge_y=c(0,)
                ## - subset(dat, wt > 3)$wt
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
    
    if(table) {
        mytheme <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = 2.0)),
        colhead = list(fg_params=list(cex = 1.0)),
        rowhead = list(fg_params=list(cex = 1.0)))
        
        df = grouped[,c('pooled_auc', 'average_auc')]
        df = round(df, 3)
        df$ci = paste(
            round(grouped$average_auc_ci_lower, 3), '-', round(grouped$average_auc_ci_upper, 3)
        )
        colnames(df) = c('AUC', 'AUC', 'AUC 95% CI')
        rownames(df) = grouped$group
        
        if(!is.null(table_row_colors)) {
            # note this very muh does not work due to some kind of bug in tableGrob,
            # a workaround is defined below
            col = table_row_colors[grouped$color]
            col = as.vector(col)
        }
        else {
            col='black'
        }
        
        table_theme = gridExtra::ttheme_minimal(
            core = list(
                fg_params=list(cex = 0.9),
                bg_params=list(fill='grey97')
            ),
            colhead = list(fg_params=list(cex = 0.9)),
            rowhead = list(fg_params=list(cex = 0.9, col=col, fontface=1))
        )
        table_grid = gridExtra::tableGrob(
            df, theme=table_theme
        )
        
        if(!is.null(table_row_colors)) {
            i = 2
            for(row_col in col) {
                # row column, for our purpose starts at 2
                ind <- find_cell(table_grid, i, 1, "rowhead-fg")
                table_grid$grobs[ind][[1]][["gp"]] <- grid::gpar(col=row_col, cex=0.9)
                i = i + 1
            }
        }
        #print(table_grid)
        #print(df[1:2, 1:2])
        header <- gridExtra::tableGrob(
            df[, 1:2], theme=table_theme,
            #rows=rownames(df),
            cols=c('Pooled', 'Averaged')
        ) 

        merged_table <- gridExtra::combine(header[1,], table_grid, along=2)
        merged_table$layout[1:4 , c("l", "r")] <- list(c(2, 3), c(2, 4))

        p = p + annotation_custom(
            merged_table,
            xmin=table_xmin, xmax=table_xmax, ymin=0.05 * nrow(df) + 0.07, ymax=0
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


                                                     

frequency_line = list(
    # greater equal
    '\U2265 0.5'='dotted',
    '\U2265 0.6'='dotdash',
    '\U2265 0.7'='dashed',
    '\U2265 0.8'='longdash',
    '\U2265 0.9'='solid'
)

select_edges = function(cn, key, q, names) {
    top = cn[cn[[key]] >= quantile(cn[[key]], 1 - q), ]
    top$group = names[1]
    bottom = cn[cn[[key]] <= quantile(cn[[key]], q), ]
    bottom$group = names[2]
    df = rbind(top, bottom)
    df$edge_label = df[[key]]
    df$edge_label_group = tools::toTitleCase(gsub('_', ' ', key))
    df
}

plot_contributions_network = function(
    cn, cn_contributions, q=0.05, conditions=c('Tuberculosis', 'Cryptococcal')
) {

    library(ggnetwork)
    library(network)
    
    cn = rbind(
        select_edges(
            cn, 'frequency', q,
            c('Most frequently co-selected', 'Least frequently co-selected')
        ),
        select_edges(
            cn, 'log2_frequency_ratio', q,
            c('Co-selected more frequently than expected', 'Co-selected less frequently than expected')
        )
    )
    
    cn_contributions$condition = as.character(ifelse(
        cn_contributions$mean > 0,
        conditions[1], conditions[2]
    ))
    cn_contributions$gene = as.character(cn_contributions$gene)

    net = ggnetwork(
        network(
            cn,
            vertex.attr=cn_contributions,
            matrix.type='edgelist',
            ignore.eval=F,
            directed=F
        ),
        by='group'
    )
    
    new_net = c()
    for (group in unique(net$group)) {
        n = net[net$group == group,]
        
        edges = n[!is.na(n$frequency),]
        genes_with_edges = c(edges$a, edges$b)
        selected = n[n$gene %in% genes_with_edges,]
        new_net = rbind(new_net, selected)
    }
    
    (
        ggplot(
            new_net,
            aes(x=x, y=y, xend=xend, yend=yend)
        )
        + facet_wrap(~ group, scale='free', ncol=2)
        + geom_edges(color="grey50", curvature=0.01)
        + geom_edgetext(
            aes(
                label=round(edge_label, 3),
                fill=edge_label_group
            ),
            color='black',
            key_glyph='rect'
        )
        + scale_fill_manual(name='Edge shows', values=c('#FDD4D0', 'grey90'))
        + geom_nodelabel_repel(
            aes(
                color=condition,
                label=gene
            ),
            fill='grey99',
            fontface='bold',
            box.padding = unit(0.7, "lines"),
            direction='y',
            segment.color='grey80',
            key_glyph='blank'
        )
        + geom_nodes(aes(size=selected_in, color=condition))
        + scale_size_area(name='Selected in')
        + scale_color_manual(values=patient_colors$Meningitis, name='High in')
        + nice_theme
        + theme(
            axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()
        )
        + guides(
            size=guide_legend(override.aes=list(color='brown')),
            color=guide_legend(override.aes=list(size=4))
        )
    )
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
            c(labels[[4*(i-1)+1]], modify), grep = TRUE,
            gp = gpar(lty=linetypes[[i]]),
            global=T
        )
    }
    g
}