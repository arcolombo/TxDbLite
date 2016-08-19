#' shove new element(s) x onto a (possibly splittable on y) string x after [z]
#' unlike Perl or C versions, this function does not modify w behind the scenes 
#' 
#' @param w   a string
#' @param x   new element(s) to shove onto w
#' @param y   a split character (default is " ")
#' @param z   after which element is x appended (default is 0, i.e. prepend x)
#' 
#' @return    the string produced by shoving x onto str2vec(w, y) after [z]
#' @importFrom BiocGenerics paste
#' @export
strunshift <- function(w, x, y=" ", z=0) {
  ww<-append(str2vec(w, y), x, z) 
  return(paste(ww,collapse=y))
}
