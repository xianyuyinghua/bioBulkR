#' Convert gene identifiers to gene symbols
#'
#' Convert UCSC, Ensembl, Entrez, or RefSeq gene identifiers to gene symbols.
#' Gene symbols can also be supplied directly to obtain a standardized
#' two-column mapping table.
#'
#' @param gene_ID_type A single character string specifying the input
#'   identifier type. Supported values are `"UCSC"`, `"ENSEMBL"`,
#'   `"ENTREZID"`, `"REFSEQ"`, and `"SYMBOL"`.
#' @param useMart_dataset A character string naming the Ensembl BioMart data
#'   set, for example, `"hsapiens_gene_ensembl"`. This argument is used only
#'   when `gene_ID_type = "UCSC"`.
#' @param gene_id A vector of gene identifiers to convert.
#' @param species_OrgDb A character string naming the organism annotation
#'   package used by [clusterProfiler::bitr()], such as `"org.Hs.eg.db"` or
#'   `"org.Mm.eg.db"`. For human Ensembl identifiers, the value
#'   `"org.Hs.eg.db"` selects the local GENCODE annotation branch.
#' @param gencode_annotation A character string giving the path to a
#'   comma-separated GENCODE annotation file. The file must contain
#'   `ensembl_id`, `gene_name`, and `gene_type` columns. It is used only for
#'   human Ensembl identifiers when `species_OrgDb = "org.Hs.eg.db"`.
#' @param gene_type A character string specifying the GENCODE gene biotype to
#'   retain. The default is `"protein_coding"`.
#'
#' @return A two-column object with columns `Gene_id` and `Symbol`. Identifiers
#'   without a matching annotation may be omitted. One input identifier can
#'   produce multiple rows when multiple mappings are available.
#'
#' @details Conversion is performed differently for each identifier type:
#'
#'   * `UCSC` uses [biomaRt::useMart()] and [biomaRt::getBM()] and therefore
#'     requires an internet connection.
#'   * Human `ENSEMBL` identifiers use the local file supplied through
#'     `gencode_annotation` when `species_OrgDb = "org.Hs.eg.db"`.
#'   * Other `ENSEMBL`, `ENTREZID`, and `REFSEQ` identifiers use
#'     [clusterProfiler::bitr()] with `species_OrgDb`.
#'   * `SYMBOL` identifiers are returned unchanged in the standardized table.
#'
#' @seealso [identify_gene_ID_type()], [clusterProfiler::bitr()],
#'   [biomaRt::getBM()]
#'
#' @examples
#' \dontrun{
#' trans_gene_ID_type(
#'   gene_ID_type = "ENTREZID",
#'   gene_id = c("7157", "672"),
#'   species_OrgDb = "org.Hs.eg.db"
#' )
#'
#' trans_gene_ID_type(
#'   gene_ID_type = "UCSC",
#'   useMart_dataset = "hsapiens_gene_ensembl",
#'   gene_id = c("uc001gzh.2", "uc003tqk.3"),
#'   species_OrgDb = "org.Hs.eg.db"
#' )
#' }
#'
#' @export
trans_gene_ID_type <- function(gene_ID_type = gene_ID_type, useMart_dataset = useMart_dataset, gene_id = gene_id, species_OrgDb = species_OrgDb,
                               gencode_annotation = "/mnt/g/BioinforMationData/Annotation/gencode.v36.annotation.csv",gene_type = "protein_coding"
                              ){
    # 转换 geneID 为 Symbol
    if( gene_ID_type == "UCSC" ){
        mart <- biomaRt::useMart("ensembl", dataset = useMart_dataset)
        df <- biomaRt::getBM(attributes = c("ucsc", "hgnc_symbol"),filters = "ucsc",values = gene_id, mart = mart)
    }else if(gene_ID_type %in% c("ENSEMBL")){
        if(species_OrgDb == "org.Hs.eg.db"){
            df_anno <- data.table::fread(gencode_annotation) %>% dplyr::filter(gene_type == gene_type) %>% dplyr::select(ensembl_id,gene_name)
            df <- df_anno %>% dplyr::filter(ensembl_id %in% gene_id)
        }else{
            df <- bitr(gene_id, fromType = gene_ID_type, toType = "SYMBOL", OrgDb = species_OrgDb)
        }
    }else if (gene_ID_type %in% c("ENTREZID","REFSEQ")){ 
        df <- bitr(gene_id, fromType = gene_ID_type, toType = "SYMBOL", OrgDb = species_OrgDb)
    }else if(gene_ID_type == "SYMBOL"){
        df <- data.frame(Gene_id = gene_id, Symbol = gene_id )
    }else{
        stop("\n",yellow,"Gene ID type is unknown. Please provide a valid gene ID type.",reset,"\n")
    }
      
    colnames(df) <- c("Gene_id","Symbol")
    
    return(df)
}
