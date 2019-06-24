import::here(patient_colors, .from='colors.R')
#library(ComplexHeatmap)

symbol_for_tb_status = function (statuses) {
    unlist(lapply(
        statuses,
        function(status) {
            if (status == '-') {
                NA
            } else if (status == 'Probable') {
                'R'
            } else
                substr(status, 1, 1)
        }
    ))
}

simple_clinical_annotation = function(annotations, limit_to=NULL) {
    if (!is.null(limit_to))
        annotations = annotations[limit_to,]
    
    include_tb_status = nrow(unique(annotations['Tuberculosis status'])) > 1
    
    if(include_tb_status) {
        clinical_annotation = ComplexHeatmap::HeatmapAnnotation(
            'HIV status'=annotations[,'HIV status'],
            'Tuberculosis status'=ComplexHeatmap::anno_simple(
                annotations[,'Tuberculosis status'],
                pch=symbol_for_tb_status(annotations[,'Tuberculosis status']),
                col=patient_colors[['Tuberculosis status']],
                pt_gp = grid::gpar(col="white"),

            ),
            'Meningitis'=annotations[,'Meningitis'],
            col=patient_colors,
            annotation_legend_param = list(
                'HIV status' = list(nrow=1, title_position="topleft"),
                'Meningitis' = list(nrow=1, title_position="topleft")
            )
        )
    } else {
        clinical_annotation = ComplexHeatmap::HeatmapAnnotation(
            'HIV status'=annotations[,'HIV status'],
            'Meningitis'=annotations[,'Meningitis'],
            col=patient_colors,
            annotation_legend_param = list(
                'HIV status' = list(nrow=1, title_position="topleft"),
                'Meningitis' = list(nrow=1, title_position="topleft")
            )
        )
    }
    
    legends = list()
    if(include_tb_status) {
        status_legend = ComplexHeatmap::Legend(
            title="Tuberculosis status",
            legend_gp = grid::gpar(
                col ='white',
                fill=patient_colors[['Tuberculosis status']]
            ),
            # pch=symbol_for_tb_status(names(patient_colors[['Tuberculosis status']])),
            labels=names(patient_colors[['Tuberculosis status']]),
            direction = "horizontal",
            nrow = 1
        )
        legends = c(legends, status_legend)
    }
    
    list(annotation=clinical_annotation, legends=legends)
}



signif_thresholds = c(0, 0.001, 0.01, 0.05, 0.1, 1)
signif_codes = c("***", "**", "*", "'", " ") # using ' instead of a dot as this is in the same line as *
        
convert_signif_codes = function(pvalues) {
    # extracted from base R, see https://stackoverflow.com/a/41263378/6646912
    symnum(
        pvalues, corr=FALSE, na=FALSE,
        cutpoints=signif_thresholds,
        symbols=signif_codes
    )
}

# annotations for use in ComplexHeatmap
annotate_pvclust = function(x, n, mapper=convert_signif_codes, col='red') {
    
    n = nobs(as.dendrogram(x$hclust))
    # Based on text.pvclust
    # Licence: GPL (≥ 2)
    # Author: Ryota Suzuki and Hidetoshi Shimodaira
    axes <- pvclust:::hc2axes(x$hclust)
    usr  <- par()$usr; wid <- usr[4] - usr[3]
    mapped = mapper(1 - x$edges[,"au"])
    au <- as.character(mapped)
    a <- grid.text(
        x=(axes[,1]-0.5)/n,
        y=axes[,2]/max(axes[,2]) + 0.01,
        au,
        gp=gpar(col=col, fontsize=20)
    )
    mapped
}

pheatmap_palette = function(mat) {
    
    max_abs = max(abs(mat))
    
    circlize::colorRamp2(
        # having symmetric 
        seq(-max_abs, max_abs, length=100),
        # TODO this is from  - add formal attribution?
        colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100)
    )
}

pheatmap_palette_with_zero = function(mat) {
    
    max_abs = max(abs(mat))
    
    circlize::colorRamp2(
        c(-max_abs, 0, max_abs),
        # TODO this is from  - add formal attribution?
        rev(c("#D73027", "#FFFFBF", "#4575B4"))
    )
}

pvclust_heatmap = function(counts_collapsed, samples_clustering, title, ...) {
    
    patients_with_rnaseq = colnames(counts_collapsed)

    clinical_annotations = simple_clinical_annotation(
        counts_patient_annotations, limit_to=patients_with_rnaseq
    )

    # first scale columns then rows
    mat = t(scale(t(scale(counts_collapsed))))

    name = 'RNA-seq heatmap'
    
    if(!is.null(samples_clustering))
        cluster_columns = samples_clustering$hclust
    else
        cluster_columns = FALSE
    
    ht_list = ComplexHeatmap::Heatmap(
        mat,
        name=name,
        column_title = title,
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        cluster_columns=cluster_columns,
        clustering_distance_rows='pearson',
        column_dend_reorder = FALSE,
        top_annotation=clinical_annotations$annotation,
        column_dend_height = unit(3, "cm"), 
        col=pheatmap_palette(mat),
        row_names_max_width = max_text_width(
            rownames(mat)
        ),
        heatmap_legend_param=list(
            title = "Normalized RNA-seq z-score",
            title_position = "topleft",
            direction = "horizontal"
            # maybe TODO, limit at=min, max
        ),
        ...
    )

    p_legend = ComplexHeatmap::Legend(
        title="Clustering p-value approximation",
        # maybe TODO: 
        # l = convert_signif_codes(1 - patients_clustering$edges$au)
        # and filter the thresholds using values in l
        labels=signif_thresholds[3:length(signif_thresholds) - 1],
        pch=signif_codes[1:length(signif_codes) - 1],
        legend_gp = grid::gpar(
            col ='red'
        ),
        type = "points",
        nrow = 1
    )

    ht_lista = ComplexHeatmap::draw(
        ht_list,
        merge_legend = TRUE,
        heatmap_legend_side = "top", 
        annotation_legend_side = "top",
        annotation_legend_list=append(clinical_annotations$legends, p_legend),
    )
    
    if(!is.null(samples_clustering))
        ComplexHeatmap::decorate_dend(name, {
            # see https://support.bioconductor.org/p/95294/#95318
            tree = ComplexHeatmap::column_dend(ht_lista)

            mapped = annotate_pvclust(samples_clustering)

            h_samples = samples_clustering$hclust
            h_heatmap = as.hclust(tree)

            if (h_samples$order != h_heatmap$order | h_samples$height != h_heatmap$height) {
                dendextend::tanglegram(dendextend::untangle(dendextend::dendlist(as.dendrogram(samples_clustering$hclust), tree), method = "step1side"))
                stop('Internal errror: passed and processed dendrograms do not match')
            }

            # maybe TODO: add rects for clusters
        }, envir = as.environment(1L))
}


compose_title = function(main, major, minor, latex=T) {
    if(is.null(main) & (!is.null(major))) {
        main = paste0(
            '\\overset{',
            major,
            '}{\\small{',
            minor,
            '}}'
        )
    }
    main = latex2exp::TeX(main)
    main
}
