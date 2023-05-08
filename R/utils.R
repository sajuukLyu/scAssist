
`%ni%` <- function(a, b) return(! a %in% b)

#' Get mild material color
#'
#' @param pal *string* describing the palette used.
#' @param num *integer* describing the number of colors.
#'
#' @return *string vector* of generated colors.
#' 
#' @importFrom ggsci pal_material
#' 
#' @export
#'
getMatColor <- function(pal = "blue", num = 3L) {
  pal_material(palette = pal, n = 2 * num + 1)(2 * num + 1)[2 * 1:num]
}
