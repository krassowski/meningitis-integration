
select_n_most_significant = function(data, n) {
    head(data[order(data$adj.P.Val),], n)
}

select_n_most_extreme = function(data, n) {
    head(data[order(-abs(data$logFC)),], n)
}

select_n_most_extreme_two_sided = function(data, n) {
    reduced_abundance = data[data$logFC < 0,]
    increased_abundance = data[data$logFC > 0,]
    rbind(
        head(reduced_abundance[order(reduced_abundance$logFC),], n),
        head(increased_abundance[order(-increased_abundance$logFC),], n)
    )
}

annotate_significance_thresholds = function(alpha, alpha_suggestive, transform=NULL) {
    if (!is.null(transform)) {
        alpha = transform(alpha)
        alpha_suggestive = transform(alpha_suggestive)
    }
    list(
        annotate('segment', x=-Inf, xend=Inf, y=alpha, yend=alpha, color='red', lty='dashed'),
        annotate('segment', x=-Inf, xend=Inf, y=alpha_suggestive, yend=alpha_suggestive, color='grey', lty='dashed')
    )
}

volcano_plot = function(data, facet=F, n=15, select=select_n_most_significant, significant_only=T, scale='fixed') {
    # set seed to increase labels positions reproducibility
    set.seed(0)
    data$is_significant = ifelse(data$adj.P.Val < 0.05, 'significant', 'non-significant')
    significant = data[data$is_significant=='significant',]
    
    select_n_top <- function(x) { select(x, n=n) }

    highlight_data = if(significant_only) significant else data
    
    if(facet != F){
        x = by(highlight_data, highlight_data[[facet]], select_n_top)
        selected_proteins = do.call('rbind', x)
    } else {
        selected_proteins = select_n_top(highlight_data)
    }
    label_color = sapply(selected_proteins$adj.P.Val, function(x){
        if (x < 0.05) color = 'black'
        else if (x < 0.10) color = 'grey30'
        else color = 'grey60'
        color
    })
    
    g = (
        ggplot(data, aes(x=logFC, y=-log10(adj.P.Val), color=is_significant))
        + theme_bw()
        + geom_point(alpha=0.4)
        + scale_color_manual(
            name='',
            values=c('significant'='red', 'non-significant'='black')
        )
        + ggrepel::geom_text_repel(
            data=selected_proteins,
            aes(label=protein),
            show.legend=F,
            color=label_color,
            segment.color='grey50',
            max.iter=50000, force=1, point.padding=0.3, box.padding=0.8, min.segment.length=0.2
        )
        + theme(legend.position='bottom')
    )
    if(facet != F) {
        g = g + facet_wrap(facet, scale=scale)
    }
    # disable sed
    rm(.Random.seed, envir=globalenv())
    g
}