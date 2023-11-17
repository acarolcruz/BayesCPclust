#' Auxiliary functions
#' @param x A vector with numbers
#' @export
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



logsumexp <- function(x, min_x = Inf){
  min_x <- min(min_x, min(x))
  const <- min_x / (-740) #-740
  return(list(x = x/const, min_x = min_x))
}
