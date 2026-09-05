#' Transcription factors and their coverage in the nine TF-target datasets
#'
#' The 1575 transcription factors included in the package. Only TFs present in
#' at least 2 of the 9 datasets (hTFtarget, KnockTF, TRRUST, ENCODE, FIMO_JASPAR,
#' PWMEnrich_JASPAR, CHEA, GTRD, ChIP_Atlas) are kept. Each column (besides TF)
#' is "T"/"F" indicating whether the TF is covered by that dataset.
#'
#' @name tf_list
#' @docType data
#' @format A data.frame with 1575 rows and 10 columns.
#' @source Compiled from the 9 TF-target sources listed in
#' \url{https://github.com/WangJin93/TFTF}.
#' @keywords datasets
"tf_list"

#' Sample metadata used for pan-tissue correlation analysis
#'
#' Sample-level annotation (sample id, tissue, sample type) mapping expression
#' samples of TCGA and GTEx used by \code{pantissue_cor_analysis()}.
#'
#' @name tcga_gtex
#' @docType data
#' @format A data.frame with 18102 rows and 5 columns (sample, tissue, type,
#' type2, type1).
#' @source TCGA and GTEx sample annotations.
#' @keywords datasets
"tcga_gtex"

#' CCLE cell line annotation
#'
#' Annotation of the Cancer Cell Line Encyclopedia (CCLE) cell lines, used to
#' map expression samples to cell-line primary sites.
#'
#' @name ccle_info
#' @docType data
#' @format A data.frame with 1046 rows and 14 columns.
#' @source \url{https://sites.broadinstitute.org/ccle/}
#' @keywords datasets
"ccle_info"

#' TCGA cancer type abbreviations and full names
#'
#' Mapping between TCGA cancer-type abbreviations (e.g. COAD) and their full
#' names, used when reporting pan-tissue correlation results.
#'
#' @name abrr_full
#' @docType data
#' @format A data.frame with 33 rows and 2 columns (tissue, full.name).
#' @source TCGA.
#' @keywords datasets
"abrr_full"

#' Tissue type lists for TCGA and GTEx
#'
#' Lists of available TCGA cancer types (33) and GTEx tissue types (30) used as
#' choices for correlation analyses, see \code{tissue_type()}.
#'
#' @name tissue
#' @docType data
#' @format A list of two character vectors named "TCGA" and "GTEx".
#' @source TCGA and GTEx.
#' @keywords datasets
"tissue"

#' Description of the TF-target datasets
#'
#' Links and descriptions of the databases integrated in the package
#' (hTFtarget, KnockTF, ENCODE, CHEA, TRRUST, GTRD, ChIP_Atlas and JASPAR).
#'
#' @name data_info
#' @docType data
#' @format A list of 8 entries, each with link and description.
#' @source \url{https://github.com/WangJin93/TFTF}
#' @keywords datasets
"data_info"

#' KnockTF experimental annotation
#'
#' Annotation of the KnockTF knockdown/knockout experiments (sample id, TF,
#' molecular type, knock method, biosample, platform, PubMed id ...), merged
#' with KnockTF prediction results in \code{predict_target()} and
#' \code{predict_TF()}.
#'
#' @name knocktf_data
#' @docType data
#' @format A data.frame with 1086 rows and 8 columns.
#' @source \url{https://bio.liclab.net/KnockTF/}
#' @keywords datasets
"knocktf_data"

#' Gene symbol to RefSeq/Ensembl identifier mapping
#'
#' Identifier map used to translate between gene symbols and sequence
#' identifiers when querying external TF databases.
#'
#' @name idmap
#' @docType data
#' @format A data.frame with 61852 rows and 2 columns (ID, gene).
#' @source Compiled for the package.
#' @keywords datasets
"idmap"

#' RefSeq transcript to gene symbol map
#'
#' RefSeq transcript accessions mapped to gene symbols.
#'
#' @name refgene
#' @docType data
#' @format A data.frame with 53318 rows and 2 columns (Target, symbol).
#' @source RefSeq.
#' @keywords datasets
"refgene"
