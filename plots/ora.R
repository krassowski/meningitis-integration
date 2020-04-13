penalise_for_chance_in_background_selection = function(with_bg, without_bg) {
    bg_and_no_bg = merge(with_bg$result, without_bg$result, by='term_name', suffixes = c(".bg", ".all"))
    bg_and_no_bg$ratio_studied = bg_and_no_bg$term_size.bg / bg_and_no_bg$term_size.all
    bg_and_no_bg$chance_of_inclusion_in_SOMAScan.p_value = unlist(apply(bg_and_no_bg, 1, function(x) {
        binom.test(as.numeric(x['term_size.bg']), as.numeric(x['term_size.all']), alternative='greater')$p.value
    }))
    # this is no longer easily interpretable p-value (!)
    bg_and_no_bg$weigted_p = bg_and_no_bg$p_value.bg / (1 - bg_and_no_bg$chance_of_inclusion_in_SOMAScan.p_value)
    bg_and_no_bg = bg_and_no_bg[order(bg_and_no_bg$weigted_p), ]
    # high_confidence_subset = bg_and_no_bg[bg_and_no_bg$significant.all == T,]

    bg_and_no_bg
}

abbreviate_terms = function(df, abbreviations) {
    if (!is.null(abbreviations)) {
        df$term_name = stringr::str_replace_all(df$term_name, abbreviations)
    }
    df
}

highlight_keywords = function(df, keywords) {
    if (!is.null(keywords)) {
        for (keyword_group in keywords) {
            for (keyword in keyword_group) {
                k = paste0(c('\\\\textbf{', keyword, '}'), collapse='')
                df$term_name = gsub(keyword, k, df$term_name)
            }
        }
    }
    df
}


explain_abbreviations = function(abbreviations, fontface='bold', ...) {
    annotate(
        'text', fontface=fontface,
        label=paste(abbreviations, names(abbreviations), sep=' = ', collapse='\n'),
        ...
    )
}

annotated_manhattan = function(gprofiler_ora_result, n=15, abbreviations=NULL, keywords=NULL, seed=0) {
    p = gprofiler2::gostplot(gprofiler_ora_result, capped=F, interactive=F)
    set.seed(seed)
    selected = head(tb_go_result$result, n=n)$term_id
    df <- p$data
    df <- base::subset(df, term_id %in% selected)

    df = abbreviate_terms(df, abbreviations)
    df = highlight_keywords(df, keywords)

    ymax = -log10(min(df$p_value / 2))

    p = (
        p
        + ggrepel::geom_text_repel(
            data=df, size=3, aes(label=latex2exp::TeX(term_name, output="character")),
            box.padding=.6, min.segment.length=0.2, max.iter=50000, direction='y',
            parse=T
        )
        + scale_y_continuous(expand=c(0, 0), limits=c(-1, ymax))
        + explain_abbreviations(abbreviations, x=1, y=Inf, hjust=0, vjust=1)
    )
    p
}

chance_plot = function(with_bg, without_bg, n=15, abbreviations=NULL, keywords=NULL, seed=0) {
    set.seed(seed)

    weighted_by_likelihood = penalise_for_chance_in_background_selection(with_bg=with_bg, without_bg=without_bg)

    df = head(weighted_by_likelihood, n=n)

    df = abbreviate_terms(df, abbreviations)
    df = highlight_keywords(df, keywords)

    (
        ggplot(weighted_by_likelihood, aes(y=-log10(p_value.bg), x=-log10(chance_of_inclusion_in_SOMAScan.p_value)))
        + theme_bw()
        + geom_point(aes(color=source.all), size=3)
        + ggrepel::geom_text_repel(
            data=df,
            size=4, aes(label=latex2exp::TeX(df$term_name, output="character")),
            box.padding=.6, min.segment.length=0.2, max.iter=50000, direction='both',
            parse=T
        )
        + explain_abbreviations(abbreviations, x=Inf, y=Inf, hjust=1, vjust=1)
    )
}
