
# Helper functions: Detect input type

isFile <- function(x){
  !is.null(x) && is.character(x) && length(x) == 1 &&
    tryCatch(file.exists(x), error = function(e) FALSE)
}

isURL <- function(x){
  !is.null(x) && is.character(x) && length(x) == 1 &&
    any(sapply(c("^https:", "^http:", "^www\\."), grepl, x, ignore.case = TRUE))
}

isValue <- function(x){
  !is.null(x) && is.vector(x) && length(x) == 1 &&
    !isFile(x) && !isURL(x)
}



