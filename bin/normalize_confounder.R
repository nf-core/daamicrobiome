normalize_confounders <- function(confounder) {
  if (is.null(confounder)) return(character(0))

  # if it already came as a vector, keep it
  if (length(confounder) > 1) {
    x <- as.character(confounder)
    x <- trimws(x)
    x <- x[nzchar(x)]
    return(x)
  }

  # length == 1: could be "Run" or "Run,Sex" or "[Run, Sex]"
  x <- as.character(confounder[[1]])
  x <- trimws(x)
  if (!nzchar(x) || tolower(x) %in% c("null", "na")) return(character(0))

  # remove brackets/quotes if present
  x <- gsub("^\\[", "", x)
  x <- gsub("\\]$", "", x)
  x <- gsub("[\"']", "", x)

  # split on commas
  parts <- unlist(strsplit(x, "\\s*,\\s*"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  parts
}