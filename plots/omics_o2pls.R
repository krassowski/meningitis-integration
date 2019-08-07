loadings = function(fit, which=c("Xjoint", "Yjoint", "Xorth", "Yorth"), subset=NULL) {
    which = match.arg(which)
    loading_id = switch(
        which,
        Xjoint = "W.", Yjoint = "C.",
        Xorth = "P_Yosc.", Yorth = "P_Xosc."
    )
    loadings = as.data.frame(fit[loading_id][[1]])

    colnames(loadings) = paste(which, 'loadings', 1:ncol(loadings))
    rownames(loadings) = rownames(fit[loading_id][[1]])

    if (!is.null(subset))
        loadings = loadings[, paste(which, 'loadings', subset), drop = FALSE]

    loadings
}


plot_most_extreme_loadings = function(fit, n=25, loading_name='Xjoint', i=1, x_axis='index', x_values=NA) {
    data = loadings(fit, which=loading_name, subset=i)
    loading = paste(loading_name, 'loadings', i) 
    if(!is.na(x_values))
        data[[x_axis]] = x_values
    most_extreme = data[head(order(-abs(data[[loading]])), n=n),]
    
    p = (
        ggplot(data, aes_string(y=paste0('`', loading, '`'), x=x_axis))
        + geom_point(alpha=0.5)
        + ggrepel::geom_label_repel(data=most_extreme, aes(label=label))
    )
    p
}


get_scores = function(fit, matrix) {
    id = switch(matrix, X='Tt', Y='U')
    fit[[id]]
}


matrix_scores_plot = function(fit, matrix, x=1, y=2, color='black', ...) {
    name = switch(matrix, X='T', Y='U')
    m = get_scores(fit, matrix)
    if(length(colnames(m)) < 2) {
        y_data = 1
        y_lab = '-'
        y = '-'
    } else {
        y_data = m[,y]
        y_lab = paste0(matrix, '-Score ', name, y)
        y_aes = 'y'
    }
    data = data.frame(x=m[,x], y=y_data, label=rownames(m), color=color)
    (
        ggplot(data, aes_string(x='x', y='y', label='label', color='color'))
        + geom_point()
        + ggrepel::geom_label_repel(...)
        + xlab(paste0(matrix, '-Score ', name, x))
        + ylab(y_lab)
        + stat_ellipse()
    )
}


single_scores_bar = function(fit, matrix, color='black') {
    (
        matrix_scores_plot(fit, matrix, color=color, nudge_y=0.15, hjust=0, direction='y', y=1)
        + nice_theme + coord_flip() + theme(legend.position='none') + ylim(1, 1.35)
        + theme(
            axis.line.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank()
          )
        + ylab(paste0(matrix, '-Score'))
    )
}



biplot_1D = function(fit, i=1, color, theme=color_meningitis) {
    gridExtra::grid.arrange(
        plot_most_extreme_loadings(fit, x_values=colMeans(log2(t(raw_rna_matrix[colnames(ran),rownames(ran)]))), x_axis='Mean.RNA', loading_name='Xjoint'),
        single_scores_bar(fit, 'X', color=color) + theme,
        plot_most_extreme_loadings(fit, x_values=colMeans(log2(t(raw_protein_matrix[colnames(pan),rownames(pan)]))), x_axis='Mean.Protein', loading_name='Yjoint'),
        single_scores_bar(fit, 'Y', color=color) + theme,
        ncol=4,
        widths=c(0.4, 0.1, 0.4, 0.1)
    )
}


scores.o2m <- function(x, which_part = c("Xjoint", "Yjoint", "Xorth", "Yorth"), 
                         subset = 0, ...) {
  # OmicsPLS::scores but working for orth, to be removed once 1.2.1 version is released
  # by Said el Bouhaddani et al., licence: GPL-3
  if(any(subset != abs(round(subset)))) stop("subset must be a vector of non-negative integers")
  
  which_part = match.arg(which_part)
  which_scores = switch(which_part, Xjoint = "Tt", Yjoint = "U", Xorth = "T_Yosc", Yorth = "U_Xosc")
  scores_matrix = x[[which_scores]]
  dim_names = dimnames(scores_matrix)
  if(length(subset) == 1 && subset == 0) subset = 1:ncol(scores_matrix)
  if(max(subset) > ncol(scores_matrix)) stop("Elements in subset exceed #components")
  scores_matrix = as.matrix(scores_matrix[,subset])
  dimnames(scores_matrix) <- dim_names
  
  return(scores_matrix)
}


scores_plot = function(fit, t=1, u=1, color='black', i='Xjoint', j='Yjoint') {
    scores_t = scores.o2m(fit, i)
    scores_u = scores.o2m(fit, j)

    data = data.frame(t=scores_t[,t], u=scores_u[,u], label=rownames(scores_t), color=color)
    (
        ggplot(data, aes(x=t, y=u, label=label, color=color))
        + geom_point()
        + ggrepel::geom_label_repel()
        + xlab(paste0(i, ' score ', t))
        + ylab(paste0(j, ' score ', u))
        + stat_ellipse(level=0.95)
    )
}

plus_minus_one_scale = function(x) {
    2 * (x - min(x)) / (max(x) - min(x)) - 1
}


# block: joint or orth
single_biplot = function(
    fit, m_i, m_j, i, j, block='joint',
    color='black', highlight='circle', n=25
) {
    i_block_id = paste0(m_i, block)
    j_block_id = paste0(m_j, block)

    scores_i = scores.o2m(fit, i_block_id)
    scores_j = scores.o2m(fit, j_block_id)

    loading_i_id = paste0(i_block_id, ' loadings ', i)
    loading_j_id = paste0(j_block_id, ' loadings ', j)

    loading_i_df = loadings(fit, i_block_id, subset=i)
    loading_j_df = loadings(fit, j_block_id, subset=j)

    loadings_df = merge(loading_i_df, loading_j_df, by="row.names")
    loadings_df$label = loadings_df$Row.names

    if(nrow(loadings_df) != nrow(loading_i_df) || nrow(loadings_df) != nrow(loading_j_df))
        warning('Loadings differ, using shared subset')

    lambda_i = max(abs(scores_i[,i]))
    lambda_j = max(abs(scores_j[,j]))

    if(!all(rownames(scores_i) == rownames(scores_j)))
        stop('Scores rownames do not match')

    scores = data.frame(
        t=scores_i[, i],# / lambda_i,
        u=scores_j[, j],# / lambda_j,
        label=rownames(scores_i),
        color=color
    )

    # impotant: loadings were scaled to the scores
    # print(mean(loadings[[loading_i]]) / mean(loadings[[loading_j]]))
    # print(lambda_i / lambda_j)

    groups = unique(scores$color)

    centroids = aggregate(. ~ color, scores, mean)

    a = scores[scores$color == groups[[1]],]
    b = scores[scores$color == groups[[2]],]

    loadings_df[[loading_i_id]] = plus_minus_one_scale(loadings_df[[loading_i_id]]) * lambda_i
    loadings_df[[loading_j_id]] = plus_minus_one_scale(loadings_df[[loading_j_id]]) * lambda_j

    loading_i = loadings_df[[loading_i_id]]
    loading_j = loadings_df[[loading_j_id]]

    m1x = mean(a$t)
    m2x = mean(b$t)
    m1y = mean(a$u)
    m2y = mean(b$u)
    x = m1x - m2x
    y = m1y - m2y
    alpha = atan2(y, x)

    intercept = (m1y - y/x * m1x + m2y - y/x * m2x)/2
    beta = (alpha - pi/2)

    intercept2 = (m1y + m2y) / 2 - tan(beta) * (m1x + m2x) / 2

    if (highlight == 'circle') {
        selection = loading_i^2 + loading_j^2
    }
    else if (highlight == 'DA') {
        ty = loading_j
        tx = loading_i
        w = (ty - intercept2) / tan(beta) - tx
        h = ty - (tx * tan(beta) + intercept2)
        gamma = atan2(w, h)
        distance = sin(gamma) * h
        selection = distance
    }
    most_extreme = loadings_df[
        head(order(-abs(selection)), n=n),
    ]

    g = ggplot(scores, aes(x=t, y=u))

    if(highlight == 'DA') {
        g = (
            g
            # DA line
            + geom_abline(intercept=intercept, slope=tan(alpha), color='grey')
            # separation line
            + geom_abline(intercept=intercept2, slope=tan(beta), color='grey')
        )
    }

    g = (
        g
        + geom_point(aes(color=color), size=2, shape=lapply(match(color, groups), function(i){
            c(2,17,0,15)[i]
        }))
        + geom_point(
            data=loadings_df,
            aes_string(
                x=paste0('`', loading_i_id, '`'),
                y=paste0('`', loading_j_id, '`')
            ),
            color='grey70'
            # alpha=0.2
        )
        + geom_point(
            data=most_extreme,
            aes_string(
                x=paste0('`', loading_i_id, '`'),
                y=paste0('`', loading_j_id, '`')
            )
        )
        # centroid
        + geom_point(data=centroids, aes(color=color), size=5, shape=18)
        + stat_ellipse(aes(color=color))

        + ggrepel::geom_label_repel(
            data=most_extreme,
            aes_string(
                x=paste0('`', loading_i_id, '`'),
                y=paste0('`', loading_j_id, '`'),
                label='label'
            )
        )
        + xlab(paste0(block, ' ', m_i, '-score ', i, ' and ', block, ' ', m_i, '-loading W', i))
        + ylab(paste0(block, ' ', m_j, '-score ', j, ' and ', block, ' ', m_j, '-loading W', j))
        #+ coord_fixed()
    )
    g
}


biplot_2D = function(fit, i=1, j=2, color, theme, ...) {
    gridExtra::grid.arrange(
        single_biplot(fit, 'X', 'X', i=i, j=j, color=color, ...) + theme,
        single_biplot(fit, 'Y', 'Y', i=i, j=j, color=color, ...) + theme,
        #single_scores_bar(fit, , ) + color_meningitis,
        #plot_most_extreme_loadings(fit, x_values=colMeans(log2(t(raw_protein_matrix[colnames(pan),rownames(pan)]))), x_axis='Mean.Protein', loading_name='Yjoint'),
        #single_scores_bar(fit, 'Y', color=color) + color_meningitis,
        ncol=2
    )
}


s_plot = function(data_matrix, fit, m, i, block='joint', n=20) {
    # References:
    # - http://metabolomics.se/Courses/MVA/MVA%20in%20Omics_Handouts_Exercises_Solutions_Thu-Fri.pdf
    # note - data_matrix should be in the same form as provided for fitting (transformed/normalized)
    # x axis represents variable magnitude (irrelevant for scaled data)
    # y axis represents "reliability"
    i_block_id = paste0(m, block)

    scores_i = scores.o2m(fit, i_block_id)[,i]

    loading_i_id = paste0(i_block_id, ' loadings ', i)
    loading_i = loadings(fit, i_block_id, subset=i)[[loading_i_id]]

    df = data.frame(
        # TODO centering?
        cov=cov(scores_i, data_matrix)[1,],
        cor=cor(scores_i, data_matrix)[1,],
        cor_pvalue=apply(
            data_matrix, 2, function(x){
                cor.test(x, scores_i)$p.value
            }
        )
    )
    df$fdr = p.adjust(df$cor_pvalue, 'BH')

    print(cor(loading_i, df$cov))
    print(cor(loading_i, df$cor))

    df$label = rownames(df)
    df$significant = df$fdr < 0.05

    low_risk = df[
        df$significant &
        abs(df$cor) > 0.7
        &
        (
            df$cov > quantile(df$cov, 0.8)
            |
            df$cov < quantile(df$cov, 0.2)
        ),
    ]

    if (nrow(low_risk) > n) {
        print(paste('Trimming', nrow(low_risk) - n))
        low_risk = head(low_risk[order(
            -abs(low_risk$cor)
            -abs(rank(low_risk$cov) - median(rank(low_risk$cov))) / nrow(low_risk)
        ), ], n)
    }

    (
        ggplot(df, aes(x=cov, y=cor))
        + geom_point(aes(color=significant))
        + ggrepel::geom_label_repel(data=low_risk, aes(label=label))
        + ggtitle(paste('S-plot of', block, m, 'component', i))
    )
}


grid_search_plot = function(grid, ncol=3) {
    (
        ggplot(grid, aes(
            y=value,
            x=param_joint_components,
            group=param_joint_components
        ))
        + facet_wrap(
            'scoring_label', scale='free_y', ncol=ncol,
            labeller=as_labeller(latex2exp::TeX, default = label_parsed)
        )
        + geom_boxplot()
        + nice_theme
    )
}


grid_orthogonal_components_plot = function(grid, value='average_cv_predictions', value_legend='mean of CV (Q^2X + Q^2Y)/2') {
    (
        ggplot(grid[grid$scoring == value,], aes(
            fill = ortho_gain_percent,
            x = param_x_ortho_components,
            y = param_y_ortho_components
        ))
        + facet_wrap(
            'param_joint_components',
            scale='free_y',
            labeller=as_labeller(function(x) paste('Joint components: ', x))
        )
        + geom_tile()
        + shadowtext::geom_shadowtext(
            aes(
                label=round(mean, 3),
                color=value_legend
            ),
            vjust=0, check_overlap=T,
            bg.r=0.2, bg.color='grey65',
        )
        + shadowtext::geom_shadowtext(
            aes(
                label=paste0(
                    sprintf(
                        "%+.2f",
                        ortho_gain_percent
                    ),
                    '%'
                ),
                color='gain over the joint component average'
            ),
            bg.r=0.2, bg.color='grey40',
            vjust=1.5, check_overlap=T
        )
        + nice_theme
        + xlab('Number of X-orthogonal components')
        + ylab('Number of Y-orthogonal components')
        + scale_color_manual(values=c('chartreuse', 'black'))
        + ggthemes::scale_fill_gradient2_tableau('Red-Green Diverging') 
        #+ scale_fill_gradient2(low='darkred', mid='grey80', high='darkgreen')
    )
}