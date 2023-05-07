
#' Better visualize features on 2D dimensional reduction plot
#' 
#' Place cells based on given 2D dimensional reduction coordinate in a specific order, and color
#' them using the values of a given feature (i.e. gene expression, some scores, etc.)
#'
#' @param obj *object* containing data needed for visualization. (Seurat object, etc.)
#' @param feat *string* describing the feature to plot. Can be a gene or a numeric metadata.
#' @param ... Other parameters
#'
#' @return *ggplot object* of the plot
#' 
#' @importFrom stats median
#' @importFrom methods is
#' @import ggplot2
#' @import ggrastr
#' @import RColorBrewer
#' @import data.table
#' 
#' @rdname plotFeature
#' @export
#'
setGeneric(
  "plotFeature",
  function(obj, feat = "TOP2A", ...) standardGeneric("plotFeature"),
  signature = "obj"
)

#'
#' @param feat_type *string* describing the type of feature. Can be one of:
#'  * "gene" : The feature is the name of a gene.
#'  * "meta" : The feature is the name of a numeric metadata.
#' @param assay *string* describing which assay is used to get gene expression.
#'  * For <u>Seurat</u> object, "RNA" is used by default.
#' @param slot *string* describing which slot of the assay is used to get gene expression.
#'  * For <u>Seurat</u> object, "data" is used by default.
#' @param reduc *string* describing which dimension reduction is used to plot cells.
#' @param dims *integer vector* describing which 2 dimensions are used to plot cells.
#' @param col *string vector* describing the colors used to setup the continuous color bar.
#' @param disp *numeric vector* describing the distribution pattern of the colors,
#' shoud have the same length as the 'col' parameter.
#' @param col_na *string* describing the color used for NA values.
#' @param cbar *string* describing the type of color bar. Can be one of:
#'  * "normal" : Normal color bar showing quartitles of the expression range.
#'  * "min-max" : Simplified color bar with only minimum and maximum labels.
#'  * "0-1" : Color bar forced to show values ranged from 0 to 1.
#'  * "pn" : Color bar showing values with both positive and negative values using
#'  different colors.
#' @param order *string* describing the way to order cells. Can be one of:
#'  * "value" : Cells will be ordered based on values of the feature.
#'  * "abs" : Cells will be ordered based on absolute values of the feature.
#'  * "given" : Cells will be ordered based on a given order.
#'  * "random" : Cells will be ordered randomly.
#' @param order_by *integer vector* describing the order of cells when the parameter
#' 'order' is set to "given".
#' @param size *numeric* describing the size of cell points.
#' @param alpha *numeric* describing the transparency of cell points.
#' @param do_raster *logical* describing whether to perform rasterization to plots.
#' @param dpi *integer* describing the quality of rasterization (dots per inch) when
#' the parameter 'do_raster' is set to TRUE.
#' @param x_name *string* describing the name of x axis.
#' @param y_name *string* describing the name of y axis.
#' @param label *logical* describing whether
#' @param label_meta *string* describing which metadata is used to label cell clusters.
#' @param label_max *integer* describing the maximum number of different labels.
#' 
#' @rdname plotFeature
#' @export
#'
setMethod(
  "plotFeature", "Seurat",
  function(
    obj,
    feat = "TOP2A", feat_type = c("gene", "meta"),
    assay = "RNA", slot = "data",
    reduc = "umap", dims = c(1L, 2L),
    col = c("gray90", brewer.pal(9, "YlOrRd")), disp = seq(0, 1, len = length(col)), col_na = "gray90",
    cbar = c("normal", "min-max", "0-1", "pn"),
    order =  c("value", "abs", "given", "random"), order_by = 1:ncol(obj),
    size = 0.4, alpha = 1,
    do_raster = F, dpi = 300L,
    x_name = paste0(toupper(reduc), "_", dims[1]), y_name = paste0(toupper(reduc), "_", dims[2]),
    label = T, label_meta = "orig.ident", label_max = 50L
  ) {
    
    stopifnot("Parameter 'dims' must be 2 different whole numbers!" = length(dims) == 2 && dims[1] != dims[2] && all.equal(dims, as.integer(dims)))
    plotData <- as.data.table(obj[[reduc]]@cell.embeddings[, dims])
    colnames(plotData) <- c("dim1", "dim2")
    
    # get feature
    feat_type <- match.arg(feat_type)
    stopifnot("Parameter 'feat' must be 1 single string." = length(feat) == 1)
    if(feat_type == "gene") {
      stopifnot("No such gene in this assay!" = feat %in% rownames(obj[[assay]][[slot]]))
      plotData$feat <- obj[[assay]][[slot]][feat, ]
    } else if(feat_type == "meta") {
      stopifnot("No such metadata in this object!" = feat %in% colnames(obj@meta.data))
      plotData$feat <- obj@meta.data[, feat]
    }
    stopifnot("Feature must be numeric!" = is.numeric(plotData$feat))
    
    # setup color bar
    cbar <- match.arg(cbar)
    if(cbar == "normal") {
      x <- plotData$feat
      b <- min(x) + diff(range(x))*seq(0, 1, len = 5)
      plotColorBar <- scale_color_gradientn(
        colors = col, values = disp,
        na.value = col_na,
        breaks = b,
        label = round(b, 1)
      )
    } else if (cbar == "min-max") {
      plotColorBar <- scale_color_gradientn(
        colors = col, values = disp,
        na.value = col_na,
        breaks = range(plotData$feat),
        label = c("min", "max")
      )
    } else if (cbar == "0-1") {
      plotColorBar <- scale_color_gradientn(
        colors = col, values = disp,
        na.value = col_na,
        limits = c(0, 1)
      )
    } else if(cbar == "pn") {
      pn <- range(plotData$feat)
      zpos <- -pn[1] / (pn[2] - pn[1])
      col <- c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF", "gray90", "#EAD397", "#FDB31A", "#E42A2A", "#A31D1D")
      disp <- c(seq(0, zpos, len = 5)[1:4], seq(zpos, 1, len = 5))
      plotColorBar <- scale_color_gradientn(
        colors = col, values = disp,
        na.value = col_na
      )
    }
    
    # setup theme
    plotTheme <- theme(
      aspect.ratio = 1,
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(hjust = 0.5),
      plot.title = element_text(face = "bold.italic", hjust = .5)
    )
    
    # setup order
    order <- match.arg(order)
    if(order == "value") {
      plotData <- plotData[order(feat)]
    } else if(order == "abs") {
      plotData <- plotData[order(abs(feat))]
    } else if(order == "given") {
      plotData <- plotData[order_by]
    } else if(order == "random") {
      plotData <- plotData[sample(1:.N)]
    }
    
    # main plot
    if(do_raster) {
      p <- ggplot(plotData, aes(x = dim1, y = dim2, color = feat)) +
        geom_point_rast(shape = 16, size = size, alpha = alpha, raster.dpi = dpi) +
        labs(x = x_name, y = y_name, title = feat, color = "") +
        plotColorBar +
        plotTheme
    } else {
      p <- ggplot(plotData, aes(x = dim1, y = dim2, color = feat)) +
        geom_point(shape = 16, size = size, alpha = alpha) +
        labs(x = x_name, y = y_name, title = feat, color = "") +
        plotColorBar +
        plotTheme
    }
    
    # add label
    if(label) {
      # check label
      if(!label_meta %in% colnames(obj@meta.data)) {
        warning("Parameter 'label_meta' should be one of metadata names.", call. = F)
      } else if(length(unique(obj@meta.data[[label_meta]])) > label_max) {
        warning("Too many labels, should find a simpler 'label_meta'.", call. = F)
      } else {
        label_x = tapply(obj[[reduc]]@cell.embeddings[, dims[1]], obj@meta.data[[label_meta]], median)
        label_y = tapply(obj[[reduc]]@cell.embeddings[, dims[2]], obj@meta.data[[label_meta]], median)
        
        labelData <- data.table(
          label = names(label_x),
          x = label_x,
          y = label_y
        )
        
        p <- p +
          geom_label(
            data = labelData, aes(x = x, y = y, label = label),
            fill = "#FFFFFF", size = 3, color = NA, alpha = .6, fontface = "bold",
            label.size = 0.5, label.r = unit(0.25, "lines"), label.padding = unit(0.15, "lines")
          ) +
          geom_text(
            data = labelData, aes(x = x, y = y, label = label),
            size = 3, color = "black", fontface = "bold"
          )
      }
    }
    
    p
  }
)
