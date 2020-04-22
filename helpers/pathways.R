import::here(extract_counts, .from='differential_expression.R')
import::here(replace_ids, .from='identifiers_mapping.R')


counts_to_pathways_space = function(counts, collection, id_type) {

    counts = extract_counts(counts)

    counts = replace_ids(counts, counts, convert_to=id_type)
    counts = as.matrix(scale(counts))

    new_counts = list()

    for (pathway in names(collection)){
        genes_in_pathway = collection[[pathway]]
        values_for_pathway = counts[rownames(counts) %in% genes_in_pathway,]
        mean_count_for_pathway = colMeans(values_for_pathway) / length(genes_in_pathway)
        new_counts[[pathway]] = mean_count_for_pathway
    }
    counts_in_pathways_space = as.data.frame(t(as.data.frame(new_counts)))

    rownames(counts_in_pathways_space) = names(collection)
    counts_in_pathways_space
}


average_pathways_expression = function(counts_with_normalization_factors, collection, id_type, pathways_subset=NULL, patients_subset=NULL) {
    normalized_gene_space_counts = extract_counts(counts_with_normalization_factors)

    pathway_space_counts = counts_to_pathways_space(normalized_gene_space_counts, collection, id_type)
    if (!is.null(pathways_subset))
        pathway_space_counts = pathway_space_counts[pathways_subset,]
    if (!is.null(patients_subset))
        pathway_space_counts = pathway_space_counts[,patients_subset]
    rowMeans(pathway_space_counts)
}



get_direct_parents = function(parents, child, parenthood) {
    direct_parents = list()

    for (parent in parents$parent) {
        # there is no other person who is my descendant and a parent (or an ancestor) of the child
        is_direct = T

        for (other in parents$parent) {
            # start walking up from another child's parent, see if can traverse to tested parent
            if (parent == other)
                next

            if (any(parenthood$parent == parent & parenthood$child == other)) {
                is_direct = F
                break
            }
        }

        if (is_direct) {
            direct_parents[length(direct_parents) + 1] = parent
        }
    }
    direct_parents
}


establish_parenthood = function(pathways, pathways_data) {
    parenthood = list()

    for (pathway in pathways) {

        for (other in pathways) {
            if (pathway == other)
                next

            other_genes = pathways_data[[other]]
            candidate_genes = pathways_data[[pathway]]
            is_super_pathway = length(setdiff(other_genes, candidate_genes)) == 0

            if (is_super_pathway) {
                common_genes = unique(c(other_genes, candidate_genes))

                parenthood[[length(parenthood) + 1]] = data.frame(
                    parent=pathway,
                    child=other,
                    distance=length(setdiff(candidate_genes, other_genes)),
                    stringsAsFactors=FALSE,
                    common_genes=paste(common_genes[order(common_genes)], collapse=',')
                )
            }

        }
    }
    as.data.frame(do.call('rbind', parenthood), stringsAsFactors=FALSE)
}


create_cluster_hierarchies = function(leaves, parenthood, collapse_parents_larger_than) {
    hierarchies = list()

    for (leaf in leaves) {
        breadcrumbs = c(leaf)
        node = leaf
        parents = parenthood[parenthood$child == node,]

        while (nrow(parents)) {

            if(nrow(parents) > 1) {
                direct_parents = get_direct_parents(parents, node, parenthood)
                if (length(direct_parents) > 1) {

                    non_direct_ancestors = setdiff(parents$parent, direct_parents)
                    is_ancestor_to_all_parents = F

                    if (length(non_direct_ancestors) == 1) {
                        ancestor = non_direct_ancestors[[1]]
                        is_ancestor_to_all_parents = T
                        for (parent in direct_parents) {
                            if (!any(parenthood$parent == ancestor & parenthood$child == parent))
                                is_ancestor_to_all_parents = F
                        }

                        if (is_ancestor_to_all_parents)
                            skipped = c()
                        else
                            skipped = non_direct_ancestors

                    } else {
                        skipped = non_direct_ancestors
                    }

                    breadcrumbs = c(
                        paste0('[', paste0(direct_parents, collapse=' | '), ']'),
                        breadcrumbs
                    )

                    if (is_ancestor_to_all_parents) {
                        breadcrumbs = c(ancestor, breadcrumbs)
                    }

                    if (length(skipped) > 0)
                        print(paste('Warning: parent node', skipped, 'will not be displayed'))

                    if (length(skipped) > 0)
                        print(breadcrumbs)
                    break
                } else {
                    node = direct_parents[1]
                }
            } else {
                node = parents$parent
            }

            breadcrumbs = c(node, breadcrumbs)
            parents = parenthood[parenthood$child == node,]
        }

        hierarchies[[length(hierarchies) + 1]] = breadcrumbs
    }

    hierarchies
}


names_for_pathway_clusters = function(
    clusters, pathways_data, pathways_ranking, character_limit=150,
    abbreviations=NULL
) {
    # abbreviations: data.frame with columns "full" and "abbreviation"
    # to replace long words with corresponding abbreviations

    cluster_names = list()

    l = character_limit

    conditions = data.frame(
        # if the generate name is longer than
        nchar=c(NA, l, l, l, l, l, l, Inf),
        # then try to take only x most significant pathways
        head=c(Inf, 10, 5, 4, 3, 2, 1, 1),
        # unless there are fewer pathways than:
        min_length=c(-Inf, 0, 5, 4, 3, 0, 0, NA)
    )
    for (pathways in clusters) {
        previous_condition = conditions[1, ]

        ranking = pathways_ranking[pathways, ]
        ranking = ranking[order(ranking$FDR), ]

        is_uncertain_cluster = FALSE
        consistency = 1

        if (length(unique(ranking$Direction)) > 1) {
            up = ranking[ranking$Direction == 'Up', ]
            down = ranking[ranking$Direction == 'Down', ]

            for (i in 1:nrow(ranking)) {
                if (length(unique(head(ranking, i)$Direction)) > 1)
                    break
            }
            consistent_until = i

            if (head(ranking, 1)$Direction == 'Up') {
                consistency = mean(ranking$Direction == 'Up')
                inconsistent_pathways = ranking[ranking$Direction == 'Down',]
            } else {
                consistency = mean(ranking$Direction == 'Down')
                inconsistent_pathways =  ranking[ranking$Direction == 'Up',]
            }

            print(
                paste0(
                    'Cluster not fully consistent; consistent until ',
                    consistent_until,
                    ' position in the FDR-based ranking (',
                    round(100 * consistent_until / nrow(ranking), 2),
                    '%); consistency = ', round(100 * consistency, 2), '%'))

            print('Top 3 pathways in the inconsistent cluster:')
            print(head(ranking, 3))
            print('Top 3 inconsistent pathways:')
            print(head(inconsistent_pathways, 3))

            is_uncertain_cluster = TRUE
        }

        for (i in 2:nrow(conditions)) {

            next_condition = conditions[i,]

            if (length(pathways) > previous_condition$min_length) {

                most_significant_subset = rownames(head(ranking, previous_condition$head))
                parenthood = establish_parenthood(most_significant_subset, pathways_data)

                # some pathways overlap in 100% percent, for example
                # "Diseases of Immune System"
                # and 
                # "Diseases associated with the TLR signaling cascade"
                full_overlap = parenthood[parenthood$distance == 0, ]

                representatives = list()

                # to resolve this, an additional annotations file would be needed;
                # instead I simply assign them as equivalent, taking the shortest one first:
                if (nrow(full_overlap) > 0) {

                    # remove the ones that overlap fully
                    most_significant_subset = most_significant_subset[!most_significant_subset %in% full_overlap$child]

                    groups = list()

                    # but keep the representatives
                    for (genes in unique(full_overlap$common_genes)) {
                        group = full_overlap[full_overlap$common_genes == genes, ]
                        group = group[order(nchar(group$child)),]
                        representative = head(group[, 'child'], 1)

                        most_significant_subset = c(most_significant_subset, representative)

                        groups[[length(groups) + 1]] = group
                        representatives[[length(representatives) + 1]] = representative
                    }

                    names(groups) = representatives

                    # recalculate parenthood, without duplicates
                    parenthood = establish_parenthood(most_significant_subset, pathways_data)
                }

                # leaf = anyone without children
                leaves = setdiff(most_significant_subset, parenthood$parent)

                if (length(leaves) == 0) {
                    print(most_significant_subset)
                    print(parenthood)
                    stop('At least one leaf required!')
                }

                hierarchies = create_cluster_hierarchies(leaves, parenthood, collapse_parents_larger_than)

                for (i in 1:length(hierarchies)) {
                    hierarchies[[i]] = lapply(hierarchies[[i]], function(pathway) {
                        if (pathway %in% representatives) {
                            paste(groups[[pathway]]$child, collapse=' >= ')
                        }
                        pathway
                    })
                }

                # abbreviate longer terms
                if (!is.null(abbreviations)) {
                    for (i in 1:length(hierarchies)) {
                        hierarchies[[i]] = lapply(hierarchies[[i]], function(pathway) {
                            for (j in rownames(abbreviations)) {
                                full = abbreviations[j, 'full']
                                abbrev = abbreviations[j, 'abbreviation']
                                pathway = gsub(full, abbrev, pathway, ignore.case=T)
                            }
                            pathway
                        })
                    }
                }

                cluster_name = cluster_name_from_hierarchies(hierarchies)

                excluded = length(pathways) - min(previous_condition$head, length(pathways))
                if (excluded != 0)
                    cluster_name = paste(cluster_name, '(+', excluded, 'more)')

                if (nchar(cluster_name) < next_condition$nchar) {
                    break
                }

            }
            previous_condition = next_condition

        }

        if (is_uncertain_cluster)
            cluster_name = paste0(cluster_name, ' [up/down consistency ', round(100 * consistency, 2), '%]')

        # TODO maybe show the most common words occuring?
        cluster_names[[length(cluster_names) + 1]] = cluster_name

    }
    cluster_names
}


collapse_pathways_groups = function(hierarchies, depth=0) {
    if (length(hierarchies) > 1) {
        top_levels = unique(unlist(lapply(hierarchies, function(h){h[[1]]})))

        if (length(top_levels) != length(hierarchies)) {
            new_hierarchies = list()

            for (top in top_levels) {

                children = Filter(function(h){h[[1]] == top & length(h) > 1}, hierarchies)
                children = lapply(children, function(c){c[2:length(c)]})

                n = 0
                if (length(children) > 1) {
                    n = length(collapse_pathways_groups(children, depth=depth+1))
                    children = cluster_name_from_hierarchies(children, sep=' | ', depth=depth+1)

                    if (n > 1)
                        children = paste0('[', children, ']')
                }

                if (depth < 1)
                    hierarchy = top
                else
                    hierarchy = '...'

                if (n > 0) {
                    hierarchy = c(hierarchy, children)
                }

                new_hierarchies[[length(new_hierarchies) + 1]] = hierarchy
            }
            hierarchies = new_hierarchies
        }
    }


    collapsed = lapply(hierarchies, function(breadcrumbs){ 
        if (length(breadcrumbs) >= 3) {
            list(breadcrumbs[[1]], '...', breadcrumbs[[length(breadcrumbs)]])
        }
        else {
            breadcrumbs
        }
    })
    pathway_groups = sapply(collapsed, paste, collapse=" > ")

    pathway_groups
}


cluster_name_from_hierarchies = function(hierarchies, sep='; ', depth=0) {
    # TODO: highlight the one with the lowest FDR

    pathway_groups = collapse_pathways_groups(hierarchies, depth=depth)

    name = paste(
        pathway_groups,
        collapse=sep
    )

    name
}


collapse_count_clusters = function(counts, clusters, names, average=colMeans) {

    collapsed_counts = data.frame()
    collapsed_pathways = c()
    i = 1
    for (cluster in clusters) {

        cluster_counts = average(counts[cluster, ])
        collapsed_counts[names[[i]],names(counts)] = cluster_counts
        i = i + 1
        collapsed_pathways = c(collapsed_pathways, cluster)
    }

    rownames(collapsed_counts) = names

    standalone_pathways = setdiff(rownames(counts), collapsed_pathways)

    standalone_counts = counts[standalone_pathways, ]
    rbind(collapsed_counts, standalone_counts)
}


collapse_ranking_clusters = function(pathways_ranking, clusters, names, collapse=colMeans) {
    pathways_ranking = pathways_ranking[, unlist(lapply(pathways_ranking, is.numeric))]
    collapsed_counts = data.frame()
    collapsed_pathways = c()
    i = 1
    for (cluster in clusters) {

        cluster_counts = collapse(pathways_ranking[cluster, ])
        collapsed_counts[names[[i]], names(pathways_ranking)] = cluster_counts
        i = i + 1
        collapsed_pathways = c(collapsed_pathways, cluster)
    }
    rownames(collapsed_counts) = names

    standalone_pathways = setdiff(rownames(pathways_ranking), collapsed_pathways)

    standalone_counts = pathways_ranking[standalone_pathways, ]
    rbind(collapsed_counts, standalone_counts)
}


get_statistic = function(coeffs, statistic, na_policy='raise') {
    if (length(statistic) != 1) {
        stop('Statistic should be a name of a single column')
    }
    values = as.numeric(coeffs[, statistic])
    names(values) = rownames(coeffs)
    if (any(is.na(values))) {
        if (na_policy == 'raise') {
            stop('NA in statistic')
        } else if (na_policy == 'drop') {
            before = length(values)
            values = na.omit(values)
            print(paste('Removed', before - length(values), 'NAs'))
        } else {
            stop('Unknown na_policy')
        }
    }
    values
}


camera_pr = function(data, statistic='mean', collection, na_policy='drop') {
    data = get_statistic(data, statistic, na_policy=na_policy)
    result = limma::cameraPR(
        data,
        collection
    )
    result[order(result$PValue), ]
}
