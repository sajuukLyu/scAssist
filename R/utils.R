
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

uniqueID <- function(name, ID) {
  stopifnot("ID must be unique!" = all(!duplicated(ID)))
  stopifnot("ID and name must have the same length!" = length(name) == length(ID))
  
  dupName <- unique(name[duplicated(name)])
  name[name %in% dupName] <- paste0(name[name %in% dupName], "_", ID[name %in% dupName])
  name
}

