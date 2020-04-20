#' Title
#'
#' @param protein_sequence 
#' @param hla_alleles 
#' @param prediction_method 
#' @param peptide_length 
#'
#' @return
#' @import Biostrings
#'
#' @examples
assemble_curlopts <- function(protein_sequence,
                              hla_alleles,
                              prediction_method,
                              peptide_length) {
    
    hla_string <- paste(hla_alleles, collapse=",")
    length_string <- paste(as.character(peptide_length), collapse=",")
    
    # if multiple sequences should be used, they can be provided in 
    # fasta format with URI-codes
    # e.g. curl --data "method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*03:01&length=9,10" 
    # http://tools-cluster-interface.iedb.org/tools_api/mhci/ 
    # > is %3E
    # \n is %0A
    
    # define a protein string that can also handle Biostring input with multiple sequences
    if (class(protein_sequence) == "AAStringSet") {
        seq_list <- as.character(protein_sequence)
        # paste together the URI formated fasta string
        protein_string <- paste0(
            paste0("%3Eseq", seq(1,length(protein_sequence)), "%0A"),
            seq_list,collapse="%0A"
            )
    } else {stop("The sequence input is not a Biostring::AAStringSet class and couls also not be converted into one.")}
    
    sprintf("method=%s&sequence_text=%s&allele=%s&length=%s",
        prediction_method, protein_string,
        hla_string, length_string)
    
}

#' Title
#'
#' @param opts_string 
#'
#' @return
#' @import curl

create_handle <- function(opts_string) {
    h <- new_handle()
    handle_setopt(h, copypostfields = opts_string)
    h
}


#' Title
#'
#' @param handle
#' @param mhc_class
#'
#' @return
#' @import curl
create_tmp_download <- function(handle, mhc_class) {
    tmp_file <- tempfile()
    if(mhc_class == "MHC-I") {
        curl_download("http://tools-cluster-interface.iedb.org/tools_api/mhci/", 
                      tmp_file, handle = handle)
        # return tmp file with output
        tmp_file
    } else {
        stop("this MHC class is not currently supported")
    }
}

query_database <- function(protein_sequence,
               hla_alleles,
               prediction_method = "consensus",
               peptide_length,
               mhc_class) {
    
    # assemble string with curl options
    opts_string <- assemble_curlopts(protein_sequence,
                      hla_alleles,
                      prediction_method,
                      peptide_length)
    
    h <- create_handle(opts_string)
    
    # return tmp file with output
    create_tmp_download(h, mhc_class= "MHC-I")
}