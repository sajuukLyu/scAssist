
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

#' Transfer given human genes to mouse homologous genes
#' 
#' Mouse genes annotation from Cellranger mm10 v.2024-A.
#' Human genes annotation from Cellranger hg38 v.2024-A.
#' biomaRt v.2.60.1
#'
#' @import magrittr
#' @import data.table
#' @importFrom homologene human2mouse
#' 
#' @export
#'
hs2mmGene <- function(genes) {
  
  mid2gene <- set_names(mmAnno$id, mmAnno$Geneid)
  geneMeta <- set_colnames(hsAnno[id %in% genes, .(id, Geneid)], c("gene", "ens"))
  
  homo_biomart <- m2h[hsapiens_homolog_ensembl_gene %in% geneMeta$ens]
  homo_biomart[, mgene := mid2gene[ensembl_gene_id]]
  
  geneMeta$mouse <- tapply(homo_biomart$mgene, homo_biomart$hsapiens_homolog_ensembl_gene, c)[geneMeta$ens]
  
  homo_ncbi <- homologene::human2mouse(genes = geneMeta$gene) %>% as.data.table()
  geneMeta$mouse_ncbi <- tapply(homo_ncbi$mouseGene, homo_ncbi$humanGene, c)[geneMeta$gene]
  
  return(geneMeta)
}

#' Transfer given mouse genes to human homologous genes
#' 
#' Mouse genes annotation from Cellranger mm10 v.2024-A.
#' Human genes annotation from Cellranger hg38 v.2024-A.
#' biomaRt v.2.60.1
#'
#' @import magrittr
#' @import data.table
#' @importFrom homologene mouse2human
#' 
#' @export
#'
mm2hsGene <- function(genes) {
  
  hid2gene <- set_names(hsAnno$id, hsAnno$Geneid)
  geneMeta <- set_colnames(mmAnno[id %in% genes, .(id, Geneid)], c("gene", "ens"))
  
  homo_biomart <- m2h[ensembl_gene_id %in% geneMeta$ens]
  homo_biomart[, hgene := hid2gene[hsapiens_homolog_ensembl_gene]]
  
  geneMeta$human <- tapply(homo_biomart$hgene, homo_biomart$ensembl_gene_id, c)[geneMeta$ens]
  
  homo_ncbi <- homologene::mouse2human(genes = geneMeta$gene) %>% as.data.table()
  geneMeta$human_ncbi <- tapply(homo_ncbi$humanGene, homo_ncbi$mouseGene, c)[geneMeta$gene]
  
  return(geneMeta)
}

#' Transfer given human gene matrix to mouse homologous gene matrix
#'
#' @import magrittr
#' @import data.table
#' 
#' @export
#'
hs2mmMat <- function(mtx) {
  
  homoData <- hs2mmGene(rownames(mtx))
  
  mmList <- map2(homoData$mouse, homoData$mouse_ncbi, ~ unique(c(.x, .y) %>% na.omit))
  mmList <- set_names(mmList, homoData$gene)
  mmLen <- lengths(mmList)
  
  monoMtx <- mtx[mmLen == 1, ] %>% set_rownames(mmList[mmLen == 1] %>% unlist %>% unname)
  multiMtx <- mtx[mmLen > 1, ]
  multiMtx <- multiMtx[mmList[mmLen > 1] %>% {rep(names(.), lengths(.))}, ] %>% set_rownames(mmList[mmLen > 1] %>% unlist %>% unname)
  
  combMtx <- rbind(monoMtx, multiMtx)
  combLen <- rownames(combMtx) %>% table
  (combLen > 1) %>% table
  
  combMono <- combMtx[rownames(combMtx) %in% names(combLen)[combLen == 1], ]
  combMulti <- combMtx[rownames(combMtx) %in% names(combLen)[combLen > 1], ] %>% {apply(., 2, function(x) tapply(x, rownames(.), mean))}
  
  return(rbind(combMono, combMulti))
}

#' Transfer given mouse gene matrix to human homologous gene matrix
#'
#' @import magrittr
#' @import data.table
#' 
#' @export
#'
mm2hsMat <- function(mtx) {
  
  homoData <- mm2hsGene(rownames(mtx))
  
  hsList <- map2(homoData$human, homoData$human_ncbi, ~ unique(c(.x, .y) %>% na.omit))
  hsList <- set_names(hsList, homoData$gene)
  hsLen <- lengths(hsList)
  
  monoMtx <- mtx[hsLen == 1, ] %>% set_rownames(hsList[hsLen == 1] %>% unlist %>% unname)
  multiMtx <- mtx[hsLen > 1, ]
  multiMtx <- multiMtx[hsList[hsLen > 1] %>% {rep(names(.), lengths(.))}, ] %>% set_rownames(hsList[hsLen > 1] %>% unlist %>% unname)
  
  combMtx <- rbind(monoMtx, multiMtx)
  combLen <- rownames(combMtx) %>% table
  (combLen > 1) %>% table
  
  combMono <- combMtx[rownames(combMtx) %in% names(combLen)[combLen == 1], ]
  combMulti <- combMtx[rownames(combMtx) %in% names(combLen)[combLen > 1], ] %>% {apply(., 2, function(x) tapply(x, rownames(.), mean))}
  
  return(rbind(combMono, combMulti))
}


