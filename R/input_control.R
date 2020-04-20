#' Title
#'
#' @param protein_input 
#'
#' @return
#' @import Biostrings
#'
#' @examples
check_protein_sequence <- function(protein_input) {
    
    # If input is character string, convert to AAStringSet
    if (is.character(protein_input)) {
        protein_sequences <- AAStringSet(protein_input)
    }
    # If AAStringSet, then it's ok
    else if (class(protein_input) == "AAStringSet") {
        protein_sequences <- protein_input
    } else {stop("Protein input is neither character string nor AAStringSet")}
    
    # The MHC tools can only handle protein sequences of standard AA alphabet 
    # no point mutations are allowed
    # is protein sequence functional and contains only protein
    letter_set <- unique(unlist(strsplit(as.character(protein_sequences), split="")))
    # return TRUE or FALSE and print error if FALSE
    protein_valid <- all(letter_set %in% AA_STANDARD)
    if (!protein_valid) {
        stop("The input protein sequence does not only contain valid AA letters.")
        # return protein sequences if letters are part of standard alphabet
    } else {protein_sequences}
}

check_hla_alleles <- function(hla_alleles, mhc_class) {
    # check if hla_alleles are in the list that is supported by IEDB
    # reference allele list
    reference_list <- tryCatch(
        # try to download reference list and report error in case not possible
        read.table('https://help.iedb.org/hc/en-us/article_attachments/114094079071/hla_ref_set.class_i.txt'),
        error = function(e){
            print("Could not download the reference HLA list.")
            print(e)
    })
    
    # data frame needs some splitting
    reference_hlas <- sapply(as.character(reference_list$V1), 
           function(X) {strsplit(X, ",")[[1]]})[1,]
    
    # check if inut alleles are contained in reference list
    bool_contained <- hla_alleles %in% reference_hlas
    
    # which is the missing allele?
    missing_allele <- hla_alleles[!bool_contained]
    # is there a missing allele?
    if(all(bool_contained)) {
        TRUE
    } else {
            stop(sprintf(
                "Not all your input HLA alleles are present in the reference list. Following input is missing: %s",
                paste(missing_allele, collapse = ", ")))
    }
}

check_prediction_method <- function(prediction_method) {
    # library of valid prediciton methods
    valid_methods <- c("ann",
        "comblib_sidney2008",
        "consensus",
        "netmhccons",
        "netmhcpan_ba",
        "netmhcpan_el",
        "netmhcstabpan",
        "pickpocket",
        "recommended",
        "smm")
    # return boolean and error message
    prediction_valid <- prediction_method %in% valid_methods
    if (!prediction_valid) {
        stop(sprintf("The input prediction method is not contained in the list of supported methods, which are: %s",
                     paste(valid_methods, collapse = ", ")))
    } else {TRUE}
}

check_peptide_length <- function(peptide_length, hla_alleles, mhc_class) {
    # peptide_length and hla_alleles must have the same lengths!
    # currently only support MHC-I
    if (mhc_class == "MHC-I") {
        # check if lengths of lists are the same
        if (length(peptide_length) == length(hla_alleles)) {
            # supported length for MHC-I: 8, 9, 10, 11, 12, 13, 14
            valid_length = c(8, 9, 10, 11, 12, 13, 14)
            # check if all input length are valid
            pep_length_valid <- all(peptide_length %in% valid_length)
            if (!pep_length_valid) {
                stop(sprintf("Input peptide lengths are not valid. Only supported lengths: %s",
                     paste(valid_length, collapse = ", ")))
            } else {TRUE}
        } else {
            stop("The list of HLA alleles and the list of peptide lengths do not have the same length.")
        }
    } else {stop("Only MHC-I!")}
}



check_input_parameters <- function(hla_alleles,
                                    prediction_method = "consensus",
                                    peptide_length,
                                    mhc_class) {
    if(mhc_class != "MHC-I") {
        stop("Only MHC-I currently supported as mhc_class.")
    } else {
    all(
        check_hla_alleles(hla_alleles, mhc_class),
        check_prediction_method(prediction_method),
        check_peptide_length(peptide_length, hla_alleles, mhc_class)
    )
    }
}