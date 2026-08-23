#' Convert a raw count matrix to CPM, FPKM, or TPM
#'
#' Import exon coordinates from a GTF annotation, calculate non-overlapping
#' exon length for each gene, retain genes shared by the count matrix and the
#' annotation, and calculate CPM, FPKM, and TPM expression values.
#'
#' @param df_count A numeric gene-by-sample count matrix or data frame. Gene
#'   identifiers must be stored in row names.
#' @param gtf_path A path to a GTF or gzip-compressed GTF annotation readable
#'   by [rtracklayer::import()]. If `NULL`, a built-in human or mouse path is
#'   selected from `species`.
#' @param core_number A positive integer giving the number of parallel workers
#'   used to calculate gene lengths.
#' @param trans_to_obj The output type: `"cpm"`, `"fpkm"`, or `"tpm"`.
#'   Any other value returns all three matrices in a list.
#' @param species An NCBI taxonomy identifier used only to select a fallback
#'   GTF path when `gtf_path = NULL`. Supported values are `9606` for human and
#'   `10090` for mouse.
#' @param gene_type The count-row identifier type. Use `"Symbol"` to match GTF
#'   `gene_name` values or `"Ensembl"` to match version-stripped GTF `gene_id`
#'   values. Matching is case-sensitive. The historical default `"Sybmol"`
#'   does not select either identifier branch, so callers should specify this
#'   argument explicitly.
#'
#' @return A numeric matrix containing the requested CPM, FPKM, or TPM values.
#'   If `trans_to_obj` is not one of the three recognized lowercase values, a
#'   named list containing `cpm`, `fpkm`, and `tpm` matrices is returned. Only
#'   genes shared by `df_count` and the selected GTF identifier are retained.
#'
#' @details
#' Gene length is calculated as the number of unique genomic positions covered
#' by all exons assigned to a gene, so overlapping exon segments are counted
#' once. Ensembl version suffixes following a period are removed before
#' matching. CPM is calculated with [edgeR::cpm()]. FPKM and TPM are calculated
#' from library totals and gene lengths in kilobases.
#'
#' The function creates a PSOCK cluster with [parallel::makeCluster()] for the
#' gene-length calculation. The current implementation does not explicitly
#' stop that cluster before returning.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(
#'   c(100, 120, 80, 95, 200, 240),
#'   nrow = 2,
#'   dimnames = list(c("TP53", "EGFR"), c("Sample1", "Sample2", "Sample3"))
#' )
#'
#' tpm <- count_to_cpm_fpkm_tpm(
#'   df_count = counts,
#'   gtf_path = "annotation/Homo_sapiens.GRCh38.gtf.gz",
#'   core_number = 2,
#'   trans_to_obj = "tpm",
#'   gene_type = "Symbol"
#' )
#' }
#'
#' @seealso [edgeR::cpm()], [rtracklayer::import()],
#'   [parallel::makeCluster()]
#' @export
count_to_cpm_fpkm_tpm <- function(df_count,
                                  gtf_path = "/data/home/zuoanjian/bioinfo_documents/annotation/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz",
                                  core_number = 10,
                                  trans_to_obj = "fpkm",  # fpkm/cpm/tpm
                                  species = 9606,  # 9606/10090
                                  gene_type = "Sybmol"  # Symbol/Ensembl
                                 ){
    library(dplyr)
    
    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    # 检查 df_count
    if (!is.matrix(df_count) && !is.data.frame(df_count)) {
        stop("输入 df_count 必须是 matrix 或 data.frame。")
    }
    if (is.null(rownames(df_count))) {
        stop("df_count 必须包含基因名作为行名。")
    }
    counts <- df_count
    

    # gtf path
    if(is.null(gtf_path)){
        if(species == 9606){
            gtf_path = "/data/home/zuoanjian/bioinfo_documents/annotation/Homo_sapiens.GRCh38.110.gtf.gz"
        }else if(species == 10090){
            gtf_path = "/data/home/zuoanjian/bioinfo_documents/annotation/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz"
        }else{
            cat("\n",yellow,"No species Obtain, Use 9606.",reset,"\n")
            gtf_path = "/data/home/zuoanjian/bioinfo_documents/annotation/Homo_sapiens.GRCh38.110.gtf.gz"
        }
            
    }
    
    gtf <- as.data.frame(x = rtracklayer::import(con = gtf_path))
    
    exon <- gtf[gtf$type == "exon",c("start", "end", "gene_name","gene_id")]
    exon$gene_id <-  exon$gene_id %>% strsplit(.,"[.]") %>% sapply(.,"[",1)
    if(gene_type == "Symbol"){
        exon_by_gene_name <- split(exon, exon$gene_name)
    }else if (gene_type == "Ensembl"){
        exon_by_gene_name <- split(exon, exon$gene_id)
    }
    
    cl <- parallel::makeCluster(core_number)
    gene_length <- parallel::parLapply(cl = cl, X = exon_by_gene_name,
                                       fun = function(x) {
                                           tmp <- apply(x, 1, function(y) { y[1]:y[2]  })
                                           length(unique(unlist(tmp)))
                                       }
                                      )
    gene_length <- data.frame(gene_name = names(gene_length),length = as.numeric(gene_length))
    intersect_features <- intersect(rownames(counts),gene_length$gene_name)
    counts <- counts[intersect_features,]

    gene_length <- as.matrix(gene_length %>% dplyr::filter(gene_name %in% intersect_features) %>% tibble::column_to_rownames("gene_name") )
    
    y <- edgeR::DGEList(counts = counts)
    cpm <- edgeR::cpm(y = y, log=FALSE)
    
    total_counts <- colSums(counts)
    gene_length_kb <- as.vector(gene_length / 1000)
    fpkm <- sweep(counts, 2, total_counts, FUN = "/") * 1e6 / gene_length_kb
    
    counts_per_kb <- sweep(counts, 1, gene_length_kb, FUN = "/")
    tpm <- sweep(counts_per_kb, 2, colSums(x = counts_per_kb), FUN = "/") * 1e6

    
    if(trans_to_obj == "cpm"){
        return(cpm)
    }else if(trans_to_obj == "fpkm"){
        return(fpkm)
    }else if(trans_to_obj == "tpm"){
        return(tpm)
    }else{
        return(list(cpm = cpm,fpkm = fpkm,tpm = tpm))
    }

}
