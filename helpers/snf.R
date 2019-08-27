library(SNFtool)

heatmap_list = function(
    heatmaps, method='spearmanr', cluster=F, order=NULL, clusters=NULL,
    legend_for=NULL, transform=identity,
    ...
) {
    heatmap_list = NULL
    for(name in names(heatmaps)) {
        heatmap = heatmaps[[name]]
        if(!is.null(order)) {
            heatmap = heatmap[order, order]
        }
        heatmap = transform(heatmap)
        heatmap_plot = omic_profiles_corr_heatmap(
            heatmap, name,
            heatmap_legend_param=if(name %in% legend_for) list() else no_legend,
            cluster=cluster, show_annotation_name=(name %in% legend_for),
            method=method, clusters=clusters,
            ...
        )
        if (name %in% legend_for)
            heatmap_list = heatmap_list + ComplexHeatmap::rowAnnotation(
                arrow=ComplexHeatmap::anno_text(c(
                    rep('', round(nrow(heatmap)/2)), '\U2192', rep('', round(nrow(heatmap)/2))
                ), gp = grid::gpar(fontsize = 15))
            )
        if (is.null(heatmap_list))
            heatmap_list = heatmap_plot
        else
            heatmap_list = heatmap_list + heatmap_plot
    }
    heatmap_list
}

run_snf = function(input_matrices, K, alpha, T) {

    distance_matrices = list()
    affinity_matrices = list()

    for (name in names(input_matrices)) {
        matrix = input_matrices[[name]]

        distance_matrices[[name]] = (dist2(as.matrix(matrix), as.matrix(matrix)))^(1/2)
        affinity_matrices[[name]] = affinityMatrix(distance_matrices[[name]], K, alpha)
    }

    result = SNF(affinity_matrices, K, T)
    result_name = do.call('paste', c(as.list(names(input_matrices)), sep=' + '))

    list(
        result=result,
        affinity_matrices=affinity_matrices,
        distance_matrices=distance_matrices,
        name=result_name
    )
}


compare_silhouette = function(
    consensus,
    distance_matrices,
    clustering_function=spectralClustering, max_n=5
) {
    df = list()
    for (C in 2:4) {
        group_CSF = clustering_function(consensus, C)
        scores = c()
        for (name in names(distance_matrices)) {
            distance_matrix = distance_matrices[name]
            sil = cluster::silhouette(group_CSF, as.data.frame(distance_matrix))
            scores[[name]] = mean(sil[,'sil_width'])
        }
        df[[length(df)+1]] = c(
            n=C,
            scores
        )
    }
    df = as.data.frame(do.call('rbind', df))
    df[['Mean']] = rowMeans(df[,names(distance_matrices)])
    df
}


get_order = function(result, clusters) {
    patients = rownames(result)
    order = c()
    for (cluster in unique(clusters)) {
        order = c(order, patients[clusters==cluster])
    }
    order
}

# multidimensional scailing
mds_plots = function(
    distance_matrices, clusters,
    patient=identity, truth=NULL, possible_outlier_zscore=1.75,
    definite_outlier_zscore=3
) {
    mds_all = data.frame()
    for (name in names(distance_matrices)) {
        matrix = distance_matrices[[name]]
        mds = cmdscale(matrix)
        colnames(mds) <- c('Dim.1', 'Dim.2')
        mds = as.data.frame(mds)
        mds$cluster = as.factor(clusters)
        mds$patient = patient(rownames(mds))
        mds$group = name
        if (!is.null(truth))
            mds$truth = truth[rownames(mds)]
        else
            mds$truth = '?'
        mds_1 = cmdscale(matrix, k=1)
        mds_1 = as.data.frame(mds_1)
        for (cluster in unique(clusters)) {
            cluster_patients = clusters == cluster
            worst_z_score = pmax(
                abs(scale(mds_1[cluster_patients, 'V1'])),
                abs(scale(mds[cluster_patients, 'Dim.1'])),
                abs(scale(mds[cluster_patients, 'Dim.2']))
            )
            mds[cluster_patients, 'outlier'] = ifelse(
                (worst_z_score > definite_outlier_zscore),
                'definite',
                ifelse(
                    (worst_z_score > possible_outlier_zscore),
                    'possible',
                    '-'
                )
            )
        }
        mds_all = rbind(mds_all, mds)
    }
    shapes = c(23, 24, 21)
    used_shapes = shapes[1:length(distance_matrices)]
    (
        ggplot(mds_all, aes(x=Dim.1, y=Dim.2, fill=truth, shape=cluster, label=patient))
        + facet_wrap('~ group', scale='free')
        + stat_ellipse(aes(group=cluster, color=cluster))
        + geom_point(size=4)
        + ggrepel::geom_label_repel(
            data=mds_all[mds_all$outlier=='definite',],
            color='darkred',
            fill='white', alpha=0.5, box.padding=0.75
        )
        + ggrepel::geom_label_repel(
            data=mds_all[mds_all$outlier=='possible',],
            color='purple',
            fill='white', alpha=0.45, box.padding=0.75
        )
        + fill_meningitis
        + nice_theme
        + scale_shape_manual(values=used_shapes)
        + guides(
            fill=guide_legend(override.aes=list(shape=21))
        )
    )
}
strip_id = function(id) {substr(id, 5, 10)}


silhouette_plots = function(clusters, distance_matrices, index) {
    par(mfrow=c(1, length(distance_matrices)))
    for (name in names(distance_matrices)) {
        distance_matrix = as.data.frame(distance_matrices[[name]])
        sil = cluster::silhouette(clusters, distance_matrix)
        rownames(sil) = index
        plot(sil, main=name, col = heat.colors(2))
    }
}

affinities_and_result = function(snf_results_list) {
    affinites_with_result = snf_results_list$affinity_matrices
    affinites_with_result[[snf_results_list$name]] = snf_results_list$result
    affinites_with_result
}

plot_matrices = function(matrices, ordering) {
    par(mfrow=c(1, length(matrices)))
    omics = names(matrices)
    for(name in omics) {
        W = matrices[[name]]
        displayClusters(W, ordering)
        title(main=name)
    }
}

plot_mean_scores_by_parameter = function (scores_grid) {
    scores_grid$alpha = scores_grid$alpha * 10
    cols = names(scores_grid)
    df = reshape::melt(scores_grid, id=cols[!cols %in% c('alpha', 'K')], variable_name='parameter')
    df$parameter = ifelse(df$parameter=='alpha', '10 * alpha', 'K')
    scores_grid$alpha = scores_grid$alpha / 10
    (
        ggplot(df, aes(x=value, y=Mean, group=value))
        + geom_boxplot()
        + facet_wrap('parameter', scale='free', shrink=FALSE)
        + ylab('Mean score')
        + xlab('')
    )
}

plot_scores_grid = function(scores_grid) {
    scores_grid$K_label = paste('K =',scores_grid$K)
    scores_grid$a_label = paste('\U03B1 =',scores_grid$alpha)
    scores_grid$n = as.factor(scores_grid$n)
    ggplot(scores_grid, aes(x=n, y=Mean, group=n)) + geom_boxplot() + facet_grid(a_label ~ K_label)
}
