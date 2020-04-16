import::here(extract_counts, .from='../helpers/differential_expression.R')
import::here(counts_to_pathways_space, .from='../helpers/pathways.R')
import::here(pvclust_heatmap, simple_clinical_annotation, .from='complex_heatmap.R')
import::here(space_to_dot, dot_to_space, remove_leading_X, .from='../helpers/utilities.R')
 

significant.limma = function(dds_result, alpha=0.05) {
    g = dds_result[complete.cases(dds_result),]
    g[g$adj.P.Val < alpha, ]
}


significant.deseq = function(dds_result, alpha=0.05) {
    g = dds_result[complete.cases(dds_result),]
    g[g$padj < alpha, ]
}


# TODO: outliers
# gene subset may be a character vector of gene ids, or a subset of (limma or DESeq2) results data frame
differential_expression_heatmap = function(genes_subset, counts, skip_cols=outliers, id_to_gene_name=NULL, ...) {

    if (nrow(genes_subset) == 0)
        return(NULL)

    padj = NA

    # DESeq2
    if ('padj' %in% colnames(genes_subset)) {
        padj = genes_subset$padj
        genes_subset = rownames(genes_subset)
    }
    # limma
    if ('adj.P.Val' %in% colnames(genes_subset)) {
        padj = genes_subset$adj.P.Val
        genes_subset = rownames(genes_subset)
    }

    counts = log2(extract_counts(counts) + 0.25)

    #if (class(counts) == 'DESeqDataSet')
    #    counts = log2(as.data.frame(DESeq2::counts(counts, normalized=TRUE)) + 0.25)

    #if (class(counts) == 'DGEList')
    #    counts = edgeR::cpm(counts, log=TRUE, prior.count=0.25)

    counts = clean_and_subset_counts(counts, genes_subset, skip_cols=skip_cols)

    if (!is.null(id_to_gene_name)) {
        rownames(counts) <- id_to_gene_name[rownames(counts), 'gene_name']
    }

    if (!is.null(padj)) {
        rownames(counts) = paste(formatC(padj, format="e", digits=1), rownames(counts), sep='\t')
    }

    set.seed(0)

    pheatmap::pheatmap(
        counts,
        show_colnames=T,
        show_rownames=T,
        annotation_col=counts_patient_annotations,
        annotation_colors=patient_colors,
        clustering_method='ward.D2',
        clustering_distance_cols='correlation',
        cluster_rows=T,
        cluster_cols=T,
        scale='row',
        clustering_distance_rows='correlation',
        ...
    )
}


clean_and_subset_counts = function(counts, subset, skip_cols=NA) {
    colnames(counts) <- remove_leading_X(colnames(counts))
    if (!is.null(skip_cols))
        counts = counts[,!colnames(counts) %in% skip_cols]
    counts[subset,]
}


gene_set_heatmap = function(
    pathways_subset, counts, collection, id_type,
    skip_cols=outliers, id_to_gene_name=NA, trim=35, ...
) {

    if (nrow(pathways_subset) == 0)
        return(NULL)

    padj = NA

    # camera
    if ('FDR' %in% colnames(pathways_subset)) {
        padj = pathways_subset$FDR
        pathways_subset = rownames(pathways_subset)
    }

    counts = extract_counts(counts)

    counts = counts_to_pathways_space(counts, collection, id_type=id_type)
    counts = log2(counts + 0.25)

    counts = clean_and_subset_counts(counts, pathways_subset, skip_cols=skip_cols)

    if (!is.null(padj)) {
        rownames(counts) = paste(formatC(padj, format="e", digits=1), rownames(counts), sep='\t')
    }

    rownames(counts) = strtrim(rownames(counts), trim)
    set.seed(0)

    pheatmap::pheatmap(
        counts,
        show_colnames=T,
        show_rownames=T,
        annotation_col=counts_patient_annotations,
        annotation_colors=patient_colors,
        clustering_method='ward.D2',
        clustering_distance_cols='correlation',
        cluster_rows=T,
        cluster_cols=T,
        scale='row',
        clustering_distance_rows='correlation',
        ...
    )
}


# TODO:
#advanced_gene_set_heatmap = function(pathways_subset, counts, collection, id_type, skip_cols=outliers, id_to_gene_name=NA, trim=35, ...) {
#    
#}

advanced_differential_expression_heatmap = function(
    genes_subset, counts, skip_cols=outliers, id_to_gene_name=NULL,
    main='Differential expression', patients_clustering=NULL,
    ...
) {

    if (nrow(genes_subset) == 0)
        return(NULL)

    annotation_row = NA
    padj = NA

    # DESeq2
    if ('padj' %in% colnames(genes_subset)) {
        padj = genes_subset$padj
        lfc = genes_subset$log2FoldChange
        average_expression = genes_subset$baseMean
        genes_subset = rownames(genes_subset)
    }
    # limma
    if ('adj.P.Val' %in% colnames(genes_subset)) {
        padj = genes_subset$adj.P.Val
        lfc = genes_subset$logFC
        average_expression = genes_subset$AveExpr
        genes_subset = rownames(genes_subset)
    }

    counts = log2(extract_counts(counts) + 0.25)
    counts = clean_and_subset_counts(counts, genes_subset, skip_cols=skip_cols)

    if (!is.null(id_to_gene_name)) {
        original = rownames(counts)
        rownames(counts) <- id_to_gene_name[original]
        is_missing = is.na(rownames(counts))
        rownames(counts)[is_missing] = original[is_missing]
    }

    if (is.null(patients_clustering)) {
        if (nrow(counts) > 2)
            patients_clustering <- pvclust::pvclust(
                scale(counts), parallel=T,
                method.hclust="ward.D2", method.dist="correlation",
                quiet=T, ...
            )
        else {
            print('Too few observations to assess the clustering p-values')
            patients_clustering = NULL
        }
    }

    lfc_abs_max = max(abs(lfc))

    pvclust_heatmap(
        counts, patients_clustering, title=main,
        fill_title="Normalized RNA-seq z-score",
        right_annotation=ComplexHeatmap::rowAnnotation(
            'mean' = ComplexHeatmap::anno_simple(
                average_expression, col=circlize::colorRamp2(
                    c(min(average_expression), max(average_expression)),
                    c('white', '#56B870')
                )
            ),
            'log(FC)'=ComplexHeatmap::anno_simple(
                lfc, col=circlize::colorRamp2(
                    c(-lfc_abs_max, 0, +lfc_abs_max),
                    c('cornflowerblue', 'white', 'darkgoldenrod2')
                )
            ),
            '-log(FDR)'=ComplexHeatmap::anno_barplot(
                -log(padj),
                axis_param=list(side ='bottom'),
                width = unit(2, "cm")
            )
        )
    )
}

barplot = function(values, gene_sets, collection) {
    df = as.data.frame(values)
    df$gene = rownames(df)

    sets = collection[gene_sets]

    l = list()
    i = 0
    for (set_name in names(sets)) {
        genes = sets[[set_name]]
        for (gene in genes) {
            i = i + 1
            x = df[rownames(df) == gene, 'values']
            l[[i]] = data.frame(
                gene=gene,
                gene_set=set_name,
                values=x
            )
        }
    }
    df$gene_set = 'All genes'
    df2 = rbind(df, do.call(rbind, l))
    df2$mean = ave(df2$values, df2$gene_set)

    df2 <- transform(df2, variable=reorder(gene, values))
    (
        ggplot(df2, aes(x=reorder(gene, values), y=reorder(gene_set, mean), fill=values))
        + theme_bw()
        + geom_tile(colour='white')
        + theme(
            legend.position = 'bottom',
            axis.text.x = element_text(angle = 90),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
        + scale_x_discrete(breaks=df2[df2$gene_set != 'All genes', 'gene'][0:0])
        + scale_fill_gradient2(low='navy', mid='white', high='red', midpoint=0)
    )
}


de_comparison_pairplot = function(comparisons, label_top_n=5) {
    # paired plot
    df = na.omit(comparisons[comparisons$pvalue_rank < label_top_n,])
    (
        ggplot(na.omit(comparisons), aes(x=method, y=(1/pvalue_rank), color=method))
        + facet_wrap('contrast', scales='free_x')
        + geom_boxplot(alpha=0.2)
        + geom_point()
        + scale_y_log10()
        + geom_line(aes(group=gene), color='grey', alpha=0.8)
        + ggrepel::geom_label_repel(data=df, aes(label=name), direction='y', nudge_x=df$side)
        + nice_theme
        + scale_color_manual(values=ggthemes::tableau_color_pal('Tableau 10')(3))
        + geom_text(aes(label=corr), x=-Inf, y=-Inf, vjust=-1, hjust=-0.1, check_overlap=T)
    )
}
