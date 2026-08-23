#' Download and load raw Affymetrix data from GEO
#'
#' Download the supplementary files associated with a Gene Expression Omnibus
#' (GEO) Series and load its extracted Affymetrix CEL files into an
#' [affy::AffyBatch] object.
#'
#' The function first calls [GEOquery::getGEOSuppFiles()] to download all
#' supplementary files for `data_id` into `destdir`. It then asks the user to
#' extract the downloaded archive. If the automatic download fails, the user
#' is prompted to download and extract the raw files manually before loading
#' them with [affy::ReadAffy()].
#'
#' @param data_id A single GEO Series accession supplied as a character string,
#'   for example, `"GSE1009"`.
#' @param destdir A character string giving the existing directory in which
#'   supplementary files are downloaded and extracted CEL files are located.
#'
#' @return An [affy::AffyBatch] object when CEL files are successfully loaded;
#'   otherwise `NULL`.
#'
#' @details The extracted CEL files (including files ending in `.CEL.gz`) must
#'   be located directly inside `destdir`, rather than in a nested directory.
#'   This function is intended for Affymetrix data supported by the `affy`
#'   package and is interactive because it asks the user to confirm that file
#'   extraction or manual downloading is complete. The download timeout is set
#'   to 240 seconds for the duration of the call.
#'
#' @seealso [GEOquery::getGEOSuppFiles()], [affy::ReadAffy()],
#'   [download_gset()]
#'
#' @examples
#' \dontrun{
#' raw_gset <- download_array_raw_gset("GSE1009", destdir = tempdir())
#' }
#'
#' @export
download_array_raw_gset <- function(data_id, destdir) {
    # ANSI 颜色代码
    yellow <- "\033[33m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    # 定义自动下载读取array原始GEO表达数据及人工下载
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = 240)
    gset <- NULL
    tryCatch(
        {
        # 尝试下载原始数据
        GEOquery::getGEOSuppFiles(GEO = data_id, makeDirectory = FALSE, baseDir = destdir,fetch_files = TRUE)
        message("Download successful!")
        cat(paste0("\n",yellow,"Note: The downloaded raw files need to be extracted, and the folder should directly contain the corresponding CEL.gz files for each sample.",reset))
        flush.console()  # 强制刷新控制台输出
        extracted <- readline(prompt = "Have the downloaded raw files been extracted? (yes/no): ")  
        if(tolower(extracted) == "yes"){
            gset <- affy::ReadAffy(celfile.path = destdir )
        }
        if (!is.null(gset)) {
            message("Manual download successful and file loaded!")
            return(gset)
        }
        },
        error = function(e) {
            cat(paste0(cyan,"Attempt failed: ", e$message, reset, "\n"))
            cat(paste0("\n",yellow,"Please manually download the ", data_id," Raw data.",reset))
            cat(paste0("\n",yellow,"Note: The downloaded raw files need to be extracted, and the folder should directly contain the corresponding CEL.gz files for each sample.",reset))
            flush.console()  # 强制刷新控制台输出
            completed <- readline(prompt = "Have you manually downloaded the file(s)? (yes/no): ")
            if (tolower(completed) == "yes") {
                gset <- affy::ReadAffy(celfile.path = destdir)
            }
            if (!is.null(gset)) {
                message("Manual download successful and file loaded!")
                return(gset)
            }
        }
    )
}
