import::here(filter_out_low_expression, extract_counts, .from='differential_expression.R')
import::here(remove_leading_X, .from='utilities.R')


normalize_abundance = function(matrix, by_condition, file=NULL, normalization_method='TMM') {
    dge = edgeR::DGEList(counts=raw_protein_matrix, group=by_condition)
    dge = edgeR::calcNormFactors(dge, method=normalization_method)
    filtered_dge = filter_out_low_expression(dge, by_condition)
    transformed = edgeR::cpm(filtered_dge, log=T)
    colnames(transformed) = remove_leading_X(colnames(transformed))
    if(!is.null(file))
        write.csv(transformed, file=out_tmm_normalized_counts_path)
    else
        transformed
}
