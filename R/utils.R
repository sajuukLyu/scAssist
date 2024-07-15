
#' 
#' @export
#'
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

#' Replace duplicated gene name with ID
#' 
#' @export
#'
uniqueID <- function(name, ID) {
  stopifnot("ID must be unique!" = all(!duplicated(ID)))
  stopifnot("ID and name must have the same length!" = length(name) == length(ID))
  
  dupName <- unique(name[duplicated(name)])
  name[name %in% dupName] <- paste0(name[name %in% dupName], "_", ID[name %in% dupName])
  name
}

#' Calculate the percentage of all counts that belong to given sets of genes using Seurat function
#'
#' @importFrom Seurat PercentageFeatureSet
#' 
#' @export
#'
calcPercentageGeneSets <- function(obj, geneSets, assay = "RNA") {
  
  allGene <- rownames(obj[[assay]]$counts)
  for(i in names(geneSets)) {
    shareGene <- intersect(allGene, geneSets[[i]])
    obj[[i]] <- Seurat::PercentageFeatureSet(obj, features = shareGene, assay = assay)
  }
  
  obj
}

#' Find variable features while excluding given set of features using Seurat function
#'
#' @importFrom Seurat FindVariableFeatures
#' 
#' @export
#'
excluFindVariableFeatures <- function(obj, exFeat, assay = "RNA", ...) {
  
  rawCountMtx <- obj[[assay]]$counts
  clearCountMtx <- rawCountMtx[setdiff(rownames(rawCountMtx), exFeat), ]
  
  obj[[assay]]$counts <- clearCountMtx
  obj <- Seurat::FindVariableFeatures(obj, assay = assay, ...)
  obj[[assay]]$counts <- rawCountMtx
  
  obj
}

#' Calculate module scores for feature expression programs in single cells using Seurat function
#'
#' @importFrom Seurat AddModuleScore
#' 
#' @export
#'
renameAddModuleScore <- function(obj, geneList, ...) {
  
  obj <- Seurat::AddModuleScore(obj, features = geneList, ...)
  scoreMat <- obj@meta.data[, tail(1:ncol(obj@meta.data), length(geneList)), drop = F]
  
  colnames(scoreMat) <- names(geneList)
  obj@meta.data <- obj@meta.data[, 1:(ncol(obj@meta.data) - length(geneList))] %>% cbind(scoreMat)
  
  obj
}