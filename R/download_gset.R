#' Download and load a GEO data set
#'
#' Download a Gene Expression Omnibus (GEO) Series (`GSE`) or Platform
#' (`GPL`) record with [GEOquery::getGEO()]. If the automatic download fails,
#' the function prompts the user to download the required file manually and
#' then attempts to load it from `destdir`.
#'
#' For a `GSE` accession, the manually downloaded Series Matrix file should be
#' placed in `destdir`. For a `GPL` accession, place the uncompressed platform
#' annotation text file (a file whose name starts with `GPL` and ends in
#' `.txt`) in `destdir`. Because the fallback uses an interactive prompt, this
#' function is intended primarily for interactive use.
#'
#' @param data_id A single GEO accession supplied as a character string,
#'   beginning with `"GSE"` or `"GPL"` (for example, `"GSE1009"`).
#' @param destdir A character string giving the directory in which downloaded
#'   files are stored and from which manually downloaded files are read.
#'
#' @return On a successful automatic download, the object returned by
#'   [GEOquery::getGEO()] (typically a list of `ExpressionSet` objects for a
#'   `GSE` accession or a `GPL` object for a platform accession). When a
#'   manually downloaded GPL annotation file is loaded, a `data.frame` is
#'   returned. If the data cannot be downloaded or loaded, `NULL` is returned.
#'
#' @details The GEO download timeout is set to 120 seconds for the duration of
#'   the call. Setting `AnnotGPL = FALSE` and `getGPL = FALSE` means that GSE
#'   downloads do not request GPL-supplied gene annotation or a separate GPL
#'   record. `destdir` must already exist and be writable.
#'
#' @seealso [GEOquery::getGEO()]
#'
#' @examples
#' \dontrun{
#' gset <- download_gset("GSE1009", destdir = tempdir())
#' }
#'
#' @export
download_gset <- function(data_id, destdir) {
    if (!is.character(data_id) || length(data_id) != 1L ||
        is.na(data_id) || !grepl("^(GSE|GPL)[0-9]+$", data_id)) {
        stop(
            "`data_id` must be a single GEO accession in the form 'GSE' or 'GPL' followed by digits.",
            call. = FALSE
        )
    }
    if (!is.character(destdir) || length(destdir) != 1L ||
        is.na(destdir) || !nzchar(destdir)) {
        stop("`destdir` must be a single non-empty character string.", call. = FALSE)
    }
    if (!dir.exists(destdir)) {
        stop("`destdir` must be an existing directory.", call. = FALSE)
    }
    if (file.access(destdir, mode = 2L) != 0L) {
        stop("`destdir` must be writable.", call. = FALSE)
    }

    # ANSI 颜色代码
    yellow <- "\033[33m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    # 定义自动下载读取GEO表达数据及转人工下载
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = 120)
    tryCatch(
    {
        # 尝试下载数据
        gset <- GEOquery::getGEO(data_id, destdir = destdir, AnnotGPL = FALSE, getGPL = FALSE)
        if (!is.null(gset)) {
            message("Download successful!")
            return(gset)  # 成功返回 gset 对象
        }
    },
    error = function(e) {
    cat(paste0(cyan,"Attempt failed: ", e$message,reset,"\n"))
    if(grepl("^GSE",data_id)){strs <- " Series Matrix File(s)."}
    if(grepl("^GPL",data_id)){strs <- " annotation_txt."}
    cat(paste0("\n",yellow,"Please manually download the ",data_id,strs,reset))
    flush.console()  # 强制刷新控制台输出
    completed <- readline(prompt = "Have you manually downloaded the file(s)? (yes/no): ")
    if (tolower(completed) == "yes") {
        if(grepl("^GSE",data_id)){
            gset <- GEOquery::getGEO(data_id, destdir = destdir, AnnotGPL = FALSE, getGPL = FALSE)
        }else if(grepl("^GPL",data_id)){
            file_name <- list.files(path = destdir, pattern = "^GPL.*\\.txt$", full.names = FALSE)
            gset <- utils::read.delim(file = file.path(destdir,file_name), comment.char="#")
        }else{
            gset <- NULL
        }
        if (!is.null(gset)) {
          message("Manual download successful and file loaded!")
          return(gset)
        }
    }  
    }
    )
}
