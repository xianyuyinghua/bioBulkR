#' Retrieve protein-protein interactions from STRING
#'
#' Map a vector of gene identifiers to preferred STRING names, query the
#' STRING network endpoint, and retain interactions for which both preferred
#' endpoint names occur in the original input vector.
#'
#' @param genes A character vector of gene or protein identifiers accepted by
#'   the STRING `get_string_ids` API endpoint, typically gene symbols.
#' @param required_score A numeric STRING combined-score threshold, usually
#'   between 0 and 1000. Only interactions meeting the requested confidence
#'   threshold are returned by STRING.
#' @param species An NCBI taxonomy identifier passed to STRING, such as `9606`
#'   for human or `10090` for mouse.
#'
#' @return A tibble containing the STRING network interactions. Columns are
#'   supplied by the STRING TSV network endpoint and typically include STRING
#'   identifiers, `preferredName_A`, `preferredName_B`, the taxonomy ID,
#'   `score`, and evidence-specific scores. Only rows whose two preferred names
#'   both occur in `genes` are retained.
#'
#' @details
#' The function first calls the STRING `get_string_ids` endpoint and writes its
#' TSV response to a temporary file. Preferred names from that response are
#' then submitted to the STRING `network` endpoint with `required_score` and
#' `species`. The network response is downloaded into memory and parsed with
#' [readr::read_tsv()]. The temporary identifier-mapping file is removed after
#' it is read.
#'
#' Because both API requests use `https://string-db.org`, an active internet
#' connection is required. STRING synonyms that map to a preferred name not
#' present in the original `genes` vector can be removed by the final endpoint
#' filter.
#'
#' @examples
#' \dontrun{
#' interactions <- string_ppi_get_data(
#'   genes = c("TP53", "EGFR", "MYC", "AKT1"),
#'   required_score = 400,
#'   species = 9606
#' )
#' head(interactions)
#' }
#'
#' @seealso
#' \url{https://string-db.org/help/api/}
#' @export
string_ppi_get_data <- function(genes, required_score = 400, species = 9606){
    
        ids <- genes
        
        # set the http
        string_api_url = "https://string-db.org/api"
        
        identifiers <- paste0(ids,collapse = '%0d')
    
        url <- paste0(c(string_api_url, "tsv", paste0("get_string_ids?species=",species,"&identifiers=")),collapse = '/')
        url <- paste0(url,identifiers)
        tmp <- tempfile()
        curl::curl_download(url,tmp)
        identifiers <- read.delim2(tmp) %>% pull(preferredName) %>% paste0(collapse = '%0d')
        unlink(tmp)
        
    
        # Getting the STRING network interactions-----
        method = paste0("network?species=",species,"&required_score=",required_score,"&identifiers=")
        url <- paste0(c(string_api_url, "tsv", method),collapse = '/')
        url <- paste0(url,identifiers)
    
        # 直接拉取数据到内存
        res <- curl::curl_fetch_memory(url)
        
        # 将内容解码为字符并读入为 data.frame
        ppi_df <- readr::read_tsv(rawToChar(res$content), show_col_types = FALSE)
        ppi_df <- ppi_df %>%  dplyr::filter(preferredName_A %in% ids & preferredName_B %in% ids)
    
        return(ppi_df)
}
