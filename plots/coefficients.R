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
    label_pos=NULL,
    select='p_value', label='p_value', fill='selected_in',
    select_descending=TRUE, select_abs=FALSE,
    label_name=NULL, fill_name=NULL, type='contribution',
    order_by='mean'
) {
    if (!is.null(limit_to_n_most_extreme)) {
        coeffs = select_coeffs(
            coeffs, n=limit_to_n_most_extreme,
            by=select, descending=select_descending,
            abs_value=select_abs
        )
    }
    coeffs = coeffs[order(coeffs[order_by]), ]
    coeffs$gene = factor(coeffs$gene, levels=coeffs$gene)

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
                y=(
                    if (is.null(label_pos))
                    ifelse(mean>0, mean + ci + label_margin, mean - ci - label_margin)
                    else
                    label_pos
                ),
                vjust=(
                    if (is.null(label_pos))
                    ifelse(mean>0, 0, 1)
                    else
                    0.5
                )
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
        + ylab(paste('Mean', type))
        #+ scale_fill_manual(=fill_name)
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


plots_most_frequently_included = function(
    coeffs, n=20,
    frequency='selected_in',
    fill='-log10(p_value)', fill_name='P-value',
    label_name='Ratio of models including this feature',
    ...
) {
    (
        plot_coefficients(
            coeffs, limit_to_n_most_extreme=n, format='f', digits=2,
            select_descending=F,
            select=frequency, label=frequency, label_name=label_name,
            fill=fill, fill_name=fill_name,
            ...
        )
        + scale_fill_gradient(
            name=ifelse(is.null(fill_name), fill, fill_name),
            #labels = function(x) format(x, scientific = TRUE)),
            #guide=guide_legend(label.theme = element_text(angle = 90))
        )
    )
}

how_volatile = function(volatile) {
    ifelse(
        volatile > 0.01,
        ifelse(
            volatile > 0.1,
            'volatile',
            'quite stable'
        ),
        'stable'
    )

}

how_frequent = function(selected_in) {
    ifelse(
        selected_in > 0.5,
        ifelse(
            selected_in > 0.75,
            ifelse(
                selected_in == 1,
                '=1',
                '>0.75'
            ),
            '>0.50'
        ),
        '\U2264 0.50'
    )
}

frequency_fill = list(
    '\U2264 0.50'='grey',
    '>0.50'='orange',
    '>0.75'='red',
    '=1'='darkred'
)

volatility_shape = list(
    stable=21,
    volatile=4,
    'quite stable'=1
)



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


coefficients_volcano_plot = function(
    coeffs, n=10, transformation=identity, p_value='p_value', neg=TRUE,
    quite_stable_alpha=0.8, volatile_alpha=0.8, highlight_frequent=FALSE,
    safeguard=100, frequency='selected_in', frequency_title='Selected in',
    highlight_top=TRUE, highlight_volatile=FALSE, label_seed=1, box_padding=0.25,
    errorbar_alpha='default', point_color='default', label_color='default'
) {
    coeffs$is_frequent = coeffs[frequency] > 0.5
    coeffs$gene = rownames(coeffs)
    coeffs_how_volatile = how_volatile(coeffs$volatile)
    coeffs = coeffs[order(-coeffs[frequency], coeffs$volatile), ]

    if (highlight_top) {
        top_coeffs = select_coeffs(coeffs, n, by=p_value, descending=neg)
        top_coeffs = top_coeffs[
            (highlight_volatile | how_volatile(top_coeffs$volatile) != 'volatile')
            ,
        ]
    }
    else
        top_coeffs = data.frame()

    frequency_fill_aes = paste('how_frequent(', frequency, ')')

    p = (
        ggplot(
            coeffs,
            aes_string(y='mean', x=paste0((if (neg) '-' else ''), 'log10(', p_value, ')'))
        )
    )
    if (errorbar_alpha == 'default') {
        p = (
            p
            + geom_errorbar(
                aes_string(
                    ymin='mean-ci',
                    ymax='mean+ci',
                    alpha='how_volatile(volatile)',
                    color=frequency_fill_aes
                )
            )
        )
    } else {
        p = (
            p
            + geom_errorbar(
                aes_string(
                    ymin='mean-ci',
                    ymax='mean+ci',
                    color=frequency_fill_aes
                ),
                alpha=errorbar_alpha
            )
        )
    }
    point_aes = aes_string(
        fill=frequency_fill_aes,
        color=frequency_fill_aes,
        shape='how_volatile(volatile)',
        alpha='how_volatile(volatile)'
    )
    point_args = list(point_aes)
    if (point_color != 'default') {
        point_args[['color']] = point_color
    }
    point_args[['size']] = 2.5

    p = p + do.call('geom_point', point_args)

    if (highlight_frequent != FALSE) {
        frequent = coeffs[
            (
                coeffs[frequency] >= highlight_frequent
                &
                (coeffs_how_volatile != 'volatile' | highlight_volatile)
            )
            ,
        ]
        if (nrow(frequent) == 0) {
            print('none selected, change frequency threshold')
            return(ggplot())
        }
        if (nrow(frequent) > safeguard) {
            print(paste0(
                'Too many (', nrow(frequent), ') selected for highlighting.',
                ' Increase safeguard or change the criterion'
            ))
            return(ggplot())
        }
    }
    else {
        frequent = data.frame()
    }
    highlight_data = rbind(frequent, top_coeffs)
    highlight_data = highlight_data[!duplicated(highlight_data), ]
    p = (
        p
        + ggrepel::geom_label_repel(
            data=highlight_data,
            aes_string(label='gene', color=frequency_fill_aes),
            segment.color='grey40',
            box.padding=box_padding,
            force=2,
            max.overlaps=Inf,
            seed=label_seed
        )
        + ggrepel::geom_label_repel(
            data=highlight_data,
            aes_string(label='gene', color=frequency_fill_aes),
            segment.color='grey40',
            segment.alpha=0,
            box.padding=box_padding,
            force=2,
            max.overlaps=Inf,
            seed=label_seed
        )
        + scale_shape_manual(name='sign', values=volatility_shape)
        #+ ggrepel::geom_label_repel(
        #    data=top_coeffs,
        #    aes(label=gene),
        #    #box.padding=0.5,
        #    segment.color='grey20'
        #)
        + coord_flip()
        + nice_theme
        #+ scale_color_manual(values=c('grey', 'red'))
        + scale_color_manual(name=frequency_title, values=frequency_fill)
        + scale_fill_manual(name=frequency_title, values=frequency_fill)
        + scale_alpha_manual(name='sign', values=c(
            stable=1,
            'quite stable'=quite_stable_alpha,
            volatile=volatile_alpha
        ))
        + guides(
            fill=guide_legend(override.aes=list(shape=21)),
            shape=guide_legend(override.aes=list(fill='blue'))
        )
    )
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
        n = net[net$group == group, ]
        edges = n[!is.na(n$frequency), ]
        genes_with_edges = c(as.character(edges$a), as.character(edges$b))
        selected = n[n$gene %in% genes_with_edges, ]
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
