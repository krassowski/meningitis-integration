ensembl_to_entrez = read.csv('../data/ensembl_to_entrez.csv', row.names=1)
ensembl_to_gene_name = read.csv('../data/ensembl_to_gene_name.csv', row.names=1)


names_maybe_row = function(df, new_names, row) {
    if (row)
        rownames(df) = new_names
    else
        names(df) = new_names
    df
}


# assumes that rownames of source are Ensembl Gene ID, will assign new ids to target;
# target and source have to be aligned!
replace_ids = function(target, source, convert_to='symbol') {
    allowed_values = c('entrez', 'symbol')
    
    if (!convert_to %in% allowed_values)
        stop(paste('Unknown convert_to value:', convert_to, '. Choose from: ', allowed_values))

    if(convert_to == 'entrez')
        new_names = ensembl_to_entrez[rownames(source),]
    if(convert_to == 'symbol')
        new_names = ensembl_to_gene_name[rownames(source),]

    use_row = class(target) != 'numeric'
    names_maybe_row(target, new_names, row=use_row)
    
}