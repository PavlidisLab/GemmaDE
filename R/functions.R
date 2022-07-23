#' JSONify
#'
#' Assuming the input string is entirely unquoted, make it suitable for JSON.
jsonify <- function(str) {
  if (is.null(str) || length(grep("(\\[|\\])", str, value = F)) == 0) {
    return(str)
  }
  
  gsub("\\[", '["', gsub("\\]", '"]', gsub(",", '","', gsub(", ", ",", str)))) %>% parse_json(simplifyVector = T)
}

#' Parse List Entry
#'
#' Given an entry (that may or may not be NA) of semicolon delimited (; ) list (in which, any may be NA),
#' return a vector of those values, coercing NA strings into true NAs.
#'
#' @param entry NA or a character of semicolon delimited (; ) values.
#'
#' @return A character vector of values in the input character.
parseListEntry <- function(entry) {
  if (length(entry) == 0) {
    return(NA)
  }
  
  if (is.na(entry) || !grepl("; ", entry, fixed = T)) {
    return(entry)
  }
  
  cleaned <- strsplit(entry %>% as.character(), "; ", fixed = T) %>% unlist()
  ifelse("NA" == cleaned, NA, cleaned)
}
