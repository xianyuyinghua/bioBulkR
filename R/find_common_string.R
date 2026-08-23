#' Find a longest common substring
#'
#' Find a longest contiguous character sequence that occurs in every element
#' of a character vector.
#'
#' @param strings A character vector whose elements are searched for a common
#'   substring. Missing values are not supported.
#'
#' @return A single character string. If `strings` is empty, or no common
#'   substring is found, an empty string (`""`) is returned.
#'
#' @details Candidate substrings are generated from the first element of
#'   `strings`, starting with the longest possible length. Each candidate is
#'   tested against every input element. The first candidate found at the
#'   maximum matching length is returned, so when several longest common
#'   substrings exist, the one occurring earliest in the first string is
#'   selected.
#'
#'   Matching is case-sensitive and uses [base::grepl()]. Consequently,
#'   regular-expression metacharacters in a candidate substring are interpreted
#'   as regular-expression syntax rather than matched literally.
#'
#' @examples
#' find_common_string(c("sample_A", "sample_B", "sample_C"))
#' find_common_string(c("GSM123_counts", "GSM456_counts"))
#' find_common_string(c("abc", "xyz"))
#' find_common_string(character())
#'
#' @export
find_common_string <- function(strings) {
    # 找出多个字符串的共同字符
    # 如果字符串向量为空，返回空字符串
    if (length(strings) == 0) return("")
    
    # 从第一个字符串开始，逐步减少子串长度
    base_string <- strings[1]
    common_substr <- ""
    
    # 从最长的子串开始，依次减少长度
    for (len in nchar(base_string):1) {
        for (start in 1:(nchar(base_string) - len + 1)) {
            # 提取子串
            substr_candidate <- substr(base_string, start, start + len - 1)
            
            # 检查这个子串是否在所有字符串中都出现
            if (all(sapply(strings, function(x) grepl(substr_candidate, x)))) {
                return(substr_candidate)  # 找到第一个最长的公共子串时返回
            }
        }
    }
    
    return(common_substr)  # 如果没有找到公共子串，返回空字符串
}
