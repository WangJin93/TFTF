# register symbols used via non-standard evaluation (dplyr) and bundled data
# objects, so that R CMD check does not report spurious "no visible binding"
# notes.
globalVariables(c(
  ".", "TF", "Target", "Log2FC", "P_value", "cor",
  "r", "logP", "sigmark", "label", "geneA", "geneB",
  "tf_list", "tissue", "data_info", "knocktf_data",
  "ccle_info", "abrr_full"
))

.onAttach <- function(lib,pkg) {
  msg <- paste0("=========================================================================================
", pkg,"
Package URL: https://github.com/WangJin93/TFTF

If you use it in published research, please cite:
  Wang J. TFTF: An R-Based Integrative Tool for Decoding Human Transcription Factor-Target Interactions. Biomolecules. 2024; 14(7):749. https://doi.org/10.3390/biom14070749.
========================================================================================
                              --Enjoy it--")
  base::packageStartupMessage(msg)
}
