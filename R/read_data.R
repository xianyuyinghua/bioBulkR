#' Read tabular data based on its file extension
#'
#' Read a CSV, tab-delimited text, Excel, or gzip-compressed tabular file using
#' a reader selected from the file extension.
#'
#' @param file_path A single character string giving the path to the input
#'   file. Supported extensions are `.csv`, `.tsv`, `.txt`, `.xls`, `.xlsx`,
#'   and `.gz`. Extension matching is case-sensitive.
#'
#' @return A tabular object containing the imported data. The exact class
#'   depends on the input format: CSV and tab-delimited text files generally
#'   return a `data.frame`, Excel files return a tibble, and gzip-compressed
#'   files return the object produced by [data.table::fread()] or, for the
#'   fallback parser, [readr::read_delim()].
#'
#' @details Files are handled as follows:
#'
#'   * `.csv` files are read with [utils::read.csv()], with the first column
#'     used as row names and syntactic column-name conversion disabled.
#'   * `.tsv` and `.txt` files are read as header-containing, tab-delimited
#'     files with [utils::read.table()].
#'   * `.xls` and `.xlsx` files are read with [readxl::read_excel()].
#'   * `.gz` files are decompressed beside the input file using
#'     [R.utils::gunzip()] and then read with [data.table::fread()]. If the
#'     resulting first column is named `V1` and the final column contains only
#'     missing values, the uncompressed file is read again as tab-delimited
#'     data with [readr::read_delim()].
#'
#'   The `.gz` branch writes an uncompressed file with the `.gz` suffix removed
#'   and overwrites an existing file at that path. Under the default behavior
#'   of [R.utils::gunzip()], the original compressed file is removed after
#'   successful decompression.
#'
#' @examples
#' \dontrun{
#' expression_data <- read_data("data/expression.csv")
#' metadata <- read_data("data/metadata.xlsx")
#' counts <- read_data("data/counts.tsv.gz")
#' }
#'
#' @export
read_data <- function(file_path) {
    # 读取不同格式的文件
    file_ext <- tools::file_ext(file_path)
    
    if (file_ext == "csv") {
        data <- read.csv(file_path, check.names = F,row.names = 1)
    } else if (file_ext == "tsv" || file_ext == "txt") {
        data <- read.table(file_path, sep = "\t", header = TRUE)
    } else if (file_ext == "xls" || file_ext == "xlsx") {
        data <- readxl::read_excel(file_path)
    } else if (file_ext == "gz") {
        # 解压并读取
        uncompressed_file <- gsub(".gz$", "", file_path)
        R.utils::gunzip(file_path, destname = uncompressed_file, overwrite = TRUE)
        data <- data.table::fread(file = uncompressed_file)
        if( grepl("^V1$",colnames(data)[1]) & all(is.na(data[[ncol(data)]])) ){
          data <- readr::read_delim(file = uncompressed_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
        }
    } else{
        stop("Unsupported file format: ", file_ext)
    }
    return(data)
}
