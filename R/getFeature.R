#' Get feature for 2D dimensional reduction plot
#' 
#' Get 2D dimensional reduction coordinate and the values of the feature
#' from given object (i.e. gene expression, score value, metadata, etc.)
#' 
#' @param obj *object* containing data needed for visualization. (Seurat object, etc.)
#' @param assayMtx *matrix* of imputed expression value or other assay.
#' @param feat *string* describing the feature to plot. Can be a gene or a metadata.
#' @param feat_type *string* describing the type of feature. Can be one of:
#'  * "gene" : The feature is the name of a gene.
#'  * "meta" : The feature is the name of a numeric metadata.
#' @param assay *string* describing which assay is used to get gene expression.
#'  * For <u>Seurat</u> object, "RNA" is used by default.
#' @param slot *string* describing which slot of the assay is used to get gene expression.
#'  * For <u>Seurat</u> object, "data" is used by default.
#' @param reduc *string* describing which dimension reduction is used to plot cells.
#' @param dims *integer vector* describing which 2 dimensions are used to plot cells.
#' @param label *logical* describing whether to label cell clusters.
#' @param label_meta *string* describing which metadata is used to label cell clusters.
#' @param label_max *integer* describing the maximum number of different labels.
#' @param vln *logical* describing whether to add violin plot.
#' @param vln_meta *string* describing which metadata is used for violin plot.
#' @param vln_max *integer* describing the maximum number of categories of violin plot.
#' @param split *logical* describing whether to split subplots.
#' @param split_meta *string* describing which metadata is used for split subplots.
#' @param split_max *integer* describing the maximum number of categories of split subplots.
#' @param split_level *string* describing the order of categories of split subplots.
#' @param split_nrow *integer* describing the number of rows of split subplots if vln = F.
#' @param ... Other parameters
#'
#' @return *data.table* of the feature and needed metadata
#' 
#' @import purrr
#' @import magrittr
#' @import data.table
#' 
#' @rdname getFeature
#' @export
#'
getFeature <- function(
    obj,
    feat = "TOP2A",
    assayMtx = NULL,
    feat_type = c("gene", "meta"),
    assay = "RNA",
    slot = "data",
    reduc = "umap", dims = c(1L, 2L),
    label = T, label_meta = "orig.ident", label_max = 50L,
    vln = T, vln_meta = label_meta, vln_max = 50L,
    split = F, split_meta = "orig.ident", split_max = 10L, split_level = NULL,
    ...
) {
  
  stopifnot("Parameter 'dims' must be 2 different whole numbers!" = length(dims) == 2 && dims[1] != dims[2] && all.equal(dims, as.integer(dims)))
  
  if(class(obj) == "Seurat") {
    plotData <- as.data.table(obj[[reduc]]@cell.embeddings[, dims]) %>% set_colnames(c("dim1", "dim2"))
    metaData <- obj@meta.data
  } else if(class(obj) == "ArchRProject") {
    plotData <- as.data.table(obj@embeddings[[reduc]]$df[, dims]) %>% set_colnames(c("dim1", "dim2"))
    metaData <- obj@cellColData
  } else {
    stop("Sorry, object not supported yet.")
  }
  
  # get feature
  feat_type <- feat_type[1]
  stopifnot("Parameter 'feat' must be 1 single string." = length(feat) == 1)
  
  if(feat_type == "gene") {
    if(class(obj) == "Seurat") {
      obj_version <- Version(obj)$major
      if(obj_version >= 5) {
        assayMtx <- obj[[assay]][slot]
      } else {
        assayMtx <- obj[[assay]][[slot]]
      }
    }
    stopifnot("No assay matrix!" = !is.null(assayMtx))
    stopifnot("Not such gene in this assay!" = feat %in% rownames(assayMtx))
    plotData$feat <- assayMtx[feat, ]
    
  } else if(feat_type == "meta") {
    stopifnot("No such metadata in this object!" = feat %in% colnames(metaData))
    plotData$feat <- metaData[, feat] %>% as.vector()
  }
  
  # split
  if(split) {
    if(!split_meta %in% colnames(metaData)) {
      stop("Parameter 'split_meta' should be one of metadata names.", call. = F)
    } else if(length(unique(metaData[[split_meta]])) == 1) {
      warning("Only 1 categories to split, should find a better 'split_meta'.", call. = F)
    } else if(length(unique(metaData[[split_meta]])) > split_max) {
      warning("Too many categories to split, should find a simpler 'split_meta'.", call. = F)
    } else {
      plotData$split <- metaData[, split_meta] %>% as.vector()
      
      if(!is.null(split_level)) {
        splitOrder <- split_level
      } else {
        splitOrder <- unique(plotData$split)
      }
      plotData$split %<>% factor(levels = splitOrder)
    }
  }
  
  # add label
  if(label) {
    if(!label_meta %in% colnames(metaData)) {
      stop("Parameter 'label_meta' should be one of metadata names.", call. = F)
    } else if(length(unique(metaData[[label_meta]])) > label_max) {
      warning("Too many labels, should find a simpler 'label_meta'.", call. = F)
    } else {
      plotData$label <- metaData[, label_meta] %>% as.vector()
    }
  }
  
  # add violin plot
  if(vln) {
    if(!vln_meta %in% colnames(metaData)) {
      stop("Parameter 'vln_meta' should be one of metadata names.", call. = F)
    } else if(length(unique(metaData[[vln_meta]])) > vln_max) {
      warning("Too many categories, should find a simpler 'vln_meta'.", call. = F)
    } else {
      plotData$meta <- metaData[, vln_meta] %>% as.vector()
    }
  }
  
  return(plotData)
}
