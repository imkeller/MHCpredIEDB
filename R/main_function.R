#' Predict MHC binding peptides using IEDB analysis tools
#'
#' @param protein_sequence Protein sequence used for binding prediciton.
#' @param hla_alleles List of HLA alleles.
#' @param prediction_method Name of the prediction method to be used.
#' @param peptide_length List of peptide lengths (one HLA allele - one peptide length)
#' @param mhc_class MHC class. Currently only "MHC-I" is supported.
#'
#' @return dataframe
#' @export
#'
#' @examples
#' 
#' 
#' 
predict_mhc_peptides <- function(protein_input,
                                  hla_alleles,
                                  prediction_method = "consensus",
                                  peptide_length,
                                  mhc_class = "MHC-I") {
    
    # check if protein input is in AAStringSet format or convert into such
    protein_sequences <- check_protein_sequence(protein_input)
    
    # check if input parameters are valid (defined in input_control.R)
    if(!check_input_parameters(hla_alleles,
                          prediction_method,
                          peptide_length,
                          mhc_class)) {
        stop("Something is wrong with the input.")
        } else {
    
    # query IEDB database (defined in query_IEDB.R)
    res_file <- query_database(protein_sequences,
                   hla_alleles,
                   prediction_method,
                   peptide_length,
                   mhc_class)
    
    # return the output of the query
    return_output(res_file)
        }
}