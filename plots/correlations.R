plot_opposite_correlations = function(df, fdr_column=NA, p_column=NA, rows=2) {
    df = df[df$Meningitis!='Viral',]
    p = (
        ggplot(df, aes(x=RNA, y=Protein, color=Meningitis))
        + facet_wrap('gene', rows, scale='free')
        + geom_point(aes(shape=Tuberculosis))
        + geom_smooth(method='lm')
        + geom_text(aes(label=paste('d = ', round(diff, 2))), x=-Inf, y=Inf, check_overlap=T, vjust=1.5, hjust=-0.15, color='grey30')
        + geom_rect(aes(ymin=hcp_q1, ymax=hcp_q3, xmin=-Inf, xmax=Inf), alpha=0.002, fill='grey', color=NA)
        + geom_line(aes(y=hcp_q2, ymin=hcp_q1, ymax=hcp_q3, xmin=-Inf, xmax=Inf, color='Healthy control'))
        + color_meningitis + nice_theme
        + scale_shape_manual(values=c('-'=19, 'Definite'=19, 'Probable'=10, 'Possible'=1))
    )
    if(!is.na(fdr_column)) {
        p = p + geom_text(aes(label=paste('FDR = ', round(df[[fdr_column]], 2))), x=-Inf, y=Inf, check_overlap=T, vjust=3.25, hjust=-0.09, color='grey30')
    }
    if(!is.na(p_column)) {
        p = p + geom_text(aes(label=paste('p = ', formatC(df[[p_column]], format = "e", digits = 2))), x=-Inf, y=Inf, check_overlap=T, vjust=5.0, hjust=-0.1, color='grey30')
    }
    p
}
