#' Better visualize categorized metadata on 2D dimensional reduction plot
#' 
#' Place cells based on given 2D dimensional reduction coordinate in a specific order, and color
#' them using the values of a given categorized metadata (i.e. sample, celltype, etc.)
#'
#' @param obj *object* containing data needed for visualization. (Seurat object, etc.)
#' @param meta *string* describing the metadata to plot. Can be a categorized metadata.
#' @param ... Other parameters
#'
#' @return *ggplot object* of the plot
#' 
#' @importFrom stats median
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#' @import ggrastr
#' @import ggsci
#' @import data.table
#' 
#' @rdname plotMeta
#' @export
#'
setGeneric(
  "plotMeta",
  function(
    obj,
    meta = "orig.ident", meta_max = 100L,
    label = T, label_meta = meta, label_max = 50L,
    split = F, split_meta = meta, split_max = 10L,
    reduc = "umap", dims = c(1L, 2L),
    col = ggsci::pal_d3("category20")(20), legend_ncol = 1L,
    order =  c("random", "given"), order_by = 1:ncol(obj),
    size = if(ncol(obj) < 5000) {1} else {0.4},
    alpha = 1,
    do_raster = F, dpi = 300L,
    x_name = paste0(toupper(reduc), "_", dims[1]), y_name = paste0(toupper(reduc), "_", dims[2]),
    axis_tick = F,
    ...
  ) standardGeneric("plotMeta"),
  signature = "obj"
)

#'
#' @param meta_max *integer* describing the maximum number of different categories.
#' @param label *logical* describing whether to label cell clusters.
#' @param label_meta *string* describing which metadata is used to label cell clusters.
#' @param label_max *integer* describing the maximum number of different labels.
#' @param split *logical* describing whether to split subplots.
#' @param split_meta *string* describing which metadata is used for split subplots.
#' @param split_max *integer* describing the maximum number of categories of split subplots.
#' @param reduc *string* describing which dimension reduction is used to plot cells.
#' @param dims *integer vector* describing which 2 dimensions are used to plot cells.
#' @param col *string vector* describing the colors used to distinguish categories.
#' @param legend_ncol *integer* describing the column number of legends.
#' @param order *string* describing the way to order cells. Can be one of:
#'  * "random" : Cells will be ordered randomly.
#'  * "given" : Cells will be ordered based on a given order.
#' @param order_by *integer vector* describing the order of cells when the parameter
#' 'order' is set to "given".
#' @param size *numeric* describing the size of cell points.
#' @param alpha *numeric* describing the transparency of cell points.
#' @param do_raster *logical* describing whether to perform rasterization to plot.
#' @param dpi *integer* describing the quality of rasterization (dots per inch) when
#' the parameter 'do_raster' is set to TRUE.
#' @param x_name *string* describing the name of x axis.
#' @param y_name *string* describing the name of y axis.
#' @param axis_tick *logical* describing whether to show axis ticks of plot.
#' 
#' @rdname plotMeta
#' @export
#'
setMethod(
  "plotMeta", "Seurat",
  function(
    obj,
    meta = "orig.ident", meta_max = 100L,
    label = T, label_meta = meta, label_max = 50L,
    split = F, split_meta = meta, split_max = 10L,
    reduc = "umap", dims = c(1L, 2L),
    col = pal_d3("category20")(20), legend_ncol = 1L,
    order =  c("random", "given"), order_by = 1:ncol(obj),
    size = if(ncol(obj) < 5000) {1} else {0.4}, alpha = 1,
    do_raster = F, dpi = 300L,
    x_name = paste0(toupper(reduc), "_", dims[1]), y_name = paste0(toupper(reduc), "_", dims[2]),
    axis_tick = F
  ) {
    
    stopifnot("Parameter 'dims' must be 2 different whole numbers!" = length(dims) == 2 && dims[1] != dims[2] && all.equal(dims, as.integer(dims)))
    plotData <- as.data.table(obj[[reduc]]@cell.embeddings[, dims])
    colnames(plotData) <- c("dim1", "dim2")
    
    # get metadata
    stopifnot("Parameter 'meta' must be 1 single string." = length(meta) == 1)
    stopifnot("No such metadata in this object!" = meta %in% colnames(obj@meta.data))
    plotData$meta <- obj@meta.data[, meta]
    
    # get enough colors
    n <- length(unique(plotData$meta))
    stopifnot("Too many categories, must find a simpler 'meta'." = n <= meta_max)
    if(length(col) < n) {col <- colorRampPalette(col)(n)}
    
    # setup theme
    plotTheme <- theme(
      aspect.ratio = 1,
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      legend.key = element_blank(),
      axis.line = element_blank(),
      axis.title = element_text(hjust = 0.5),
      plot.title = element_text(face = "bold", hjust = .5),
      strip.background = element_blank()
    )
    if(!axis_tick) {
      plotTheme <- plotTheme + theme(
        axis.ticks = element_blank(),
        axis.text = element_blank()
      )
    }
    
    # split
    if(split) {
      if(!split_meta %in% colnames(obj@meta.data)) {
        warning("Parameter 'split_meta' should be one of metadata names.", call. = F)
      } else if(length(unique(obj@meta.data[[split_meta]])) == 1) {
        warning("Only 1 categories to split, should find a better 'split_meta'.", call. = F)
      } else if(length(unique(obj@meta.data[[split_meta]])) > split_max) {
        warning("Too many categories to split, should find a simpler 'split_meta'.", call. = F)
      } else {
        plotData$split <- obj@meta.data[, split_meta]
      }
    }
    
    # setup order
    order <- match.arg(order)
    if(order == "random") {
      plotDataOrdered <- plotData[sample(1:.N)]
    } else if(order == "given") {
      plotDataOrdered <- plotData[order_by]
    }
    
    # main plot
    if(do_raster) {
      p <- ggplot(plotDataOrdered, aes(x = dim1, y = dim2, color = meta)) +
        geom_point_rast(shape = 16, size = size, alpha = alpha, raster.dpi = dpi) +
        scale_color_manual(values = col) +
        guides(color = guide_legend(override.aes = list(size = 5, shape = 15, alpha = 1), ncol = legend_ncol)) +
        labs(x = x_name, y = y_name, title = meta, color = "") +
        plotTheme
      
      if("split" %in% colnames(plotData)) {
        p <- p + facet_wrap(~ split, nrow = 1)
      }
    } else {
      p <- ggplot(plotDataOrdered, aes(x = dim1, y = dim2, color = meta)) +
        geom_point(shape = 16, size = size, alpha = alpha) +
        scale_color_manual(values = col) +
        guides(color = guide_legend(override.aes = list(size = 5, shape = 15, alpha = 1), ncol = legend_ncol)) +
        labs(x = x_name, y = y_name, title = meta, color = "") +
        plotTheme
      
      if("split" %in% colnames(plotData)) {
        p <- p + facet_wrap(~ split, nrow = 1)
      }
    }
    
    # add label
    if(label) {
      # check label
      if(!label_meta %in% colnames(obj@meta.data)) {
        warning("Parameter 'label_meta' should be one of metadata names.", call. = F)
      } else if(length(unique(obj@meta.data[[label_meta]])) > label_max) {
        warning("Too many labels, should find a simpler 'label_meta'.", call. = F)
      } else {
        plotData$label <- obj@meta.data[[label_meta]]
        
        if(!"split" %in% colnames(plotData)) {
          labelData <- plotData[, .(x = median(dim1), y = median(dim2)), by = label]
        } else {
          labelData <- plotData[, .(x = median(dim1), y = median(dim2)), by = .(label, split)]
        }
        
        p <- p +
          geom_text(
            data = labelData, aes(x = x, y = y, label = label),
            size = 3, color = "black", fontface = "plain"
          )
      }
    }
    
    p
  }
)

#'
#' @param meta_max *integer* describing the maximum number of different categories.
#' @param label *logical* describing whether to label cell clusters.
#' @param label_meta *string* describing which metadata is used to label cell clusters.
#' @param label_max *integer* describing the maximum number of different labels.
#' @param split *logical* describing whether to split subplots.
#' @param split_meta *string* describing which metadata is used for split subplots.
#' @param split_max *integer* describing the maximum number of categories of split subplots.
#' @param reduc *string* describing which dimension reduction is used to plot cells.
#' @param dims *integer vector* describing which 2 dimensions are used to plot cells.
#' @param col *string vector* describing the colors used to distinguish categories.
#' @param legend_ncol *integer* describing the column number of legends.
#' @param order *string* describing the way to order cells. Can be one of:
#'  * "random" : Cells will be ordered randomly.
#'  * "given" : Cells will be ordered based on a given order.
#' @param order_by *integer vector* describing the order of cells when the parameter
#' 'order' is set to "given".
#' @param size *numeric* describing the size of cell points.
#' @param alpha *numeric* describing the transparency of cell points.
#' @param do_raster *logical* describing whether to perform rasterization to plot.
#' @param dpi *integer* describing the quality of rasterization (dots per inch) when
#' the parameter 'do_raster' is set to TRUE.
#' @param x_name *string* describing the name of x axis.
#' @param y_name *string* describing the name of y axis.
#' @param axis_tick *logical* describing whether to show axis ticks of plot.
#' 
#' @rdname plotMeta
#' @export
#'
setMethod(
  "plotMeta", "ArchRProject",
  function(
    obj,
    meta = "Sample", meta_max = 100L,
    label = T, label_meta = meta, label_max = 50L,
    split = F, split_meta = meta, split_max = 10L,
    reduc = "UMAP", dims = c(1L, 2L),
    col = pal_d3("category20")(20), legend_ncol = 1L,
    order =  c("random", "given"), order_by = 1:ncol(obj),
    size = if(length(obj$cellNames) < 5000) {1} else {0.4}, alpha = 1,
    do_raster = F, dpi = 300L,
    x_name = paste0(toupper(reduc), "_", dims[1]), y_name = paste0(toupper(reduc), "_", dims[2]),
    axis_tick = F
  ) {
    
    stopifnot("Parameter 'dims' must be 2 different whole numbers!" = length(dims) == 2 && dims[1] != dims[2] && all.equal(dims, as.integer(dims)))
    plotData <- as.data.table(obj@embeddings[[reduc]]$df[, dims])
    colnames(plotData) <- c("dim1", "dim2")
    
    # get metadata
    stopifnot("Parameter 'meta' must be 1 single string." = length(meta) == 1)
    stopifnot("No such metadata in this object!" = meta %in% colnames(obj@cellColData))
    plotData$meta <- as.vector(obj@cellColData[, meta])
    
    # get enough colors
    n <- length(unique(plotData$meta))
    stopifnot("Too many categories, must find a simpler 'meta'." = n <= meta_max)
    if(length(col) < n) {col <- colorRampPalette(col)(n)}
    
    # setup theme
    plotTheme <- theme(
      aspect.ratio = 1,
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      legend.key = element_blank(),
      axis.line = element_blank(),
      axis.title = element_text(hjust = 0.5),
      plot.title = element_text(face = "bold", hjust = .5),
      strip.background = element_blank()
    )
    if(!axis_tick) {
      plotTheme <- plotTheme + theme(
        axis.ticks = element_blank(),
        axis.text = element_blank()
      )
    }
    
    # split
    if(split) {
      if(!split_meta %in% colnames(obj@cellColData)) {
        warning("Parameter 'split_meta' should be one of metadata names.", call. = F)
      } else if(length(unique(as.vector(obj@cellColData[[split_meta]]))) == 1) {
        warning("Only 1 categories to split, should find a better 'split_meta'.", call. = F)
      } else if(length(unique(as.vector(obj@cellColData[[split_meta]]))) > split_max) {
        warning("Too many categories to split, should find a simpler 'split_meta'.", call. = F)
      } else {
        plotData$split <- as.vector(obj@cellColData[, split_meta])
      }
    }
    
    # setup order
    order <- match.arg(order)
    if(order == "random") {
      plotDataOrdered <- plotData[sample(1:.N)]
    } else if(order == "given") {
      plotDataOrdered <- plotData[order_by]
    }
    
    # main plot
    if(do_raster) {
      p <- ggplot(plotDataOrdered, aes(x = dim1, y = dim2, color = meta)) +
        geom_point_rast(shape = 16, size = size, alpha = alpha, raster.dpi = dpi) +
        scale_color_manual(values = col) +
        guides(color = guide_legend(override.aes = list(size = 5, shape = 15, alpha = 1), ncol = legend_ncol)) +
        labs(x = x_name, y = y_name, title = meta, color = "") +
        plotTheme
      
      if("split" %in% colnames(plotData)) {
        p <- p + facet_wrap(~ split, nrow = 1)
      }
    } else {
      p <- ggplot(plotDataOrdered, aes(x = dim1, y = dim2, color = meta)) +
        geom_point(shape = 16, size = size, alpha = alpha) +
        scale_color_manual(values = col) +
        guides(color = guide_legend(override.aes = list(size = 5, shape = 15, alpha = 1), ncol = legend_ncol)) +
        labs(x = x_name, y = y_name, title = meta, color = "") +
        plotTheme
      
      if("split" %in% colnames(plotData)) {
        p <- p + facet_wrap(~ split, nrow = 1)
      }
    }
    
    # add label
    if(label) {
      # check label
      if(!label_meta %in% colnames(obj@cellColData)) {
        warning("Parameter 'label_meta' should be one of metadata names.", call. = F)
      } else if(length(unique(as.vector(obj@cellColData[[label_meta]]))) > label_max) {
        warning("Too many labels, should find a simpler 'label_meta'.", call. = F)
      } else {
        plotData$label <- as.vector(obj@cellColData[[label_meta]])
        
        if(!"split" %in% colnames(plotData)) {
          labelData <- plotData[, .(x = median(dim1), y = median(dim2)), by = label]
        } else {
          labelData <- plotData[, .(x = median(dim1), y = median(dim2)), by = .(label, split)]
        }
        
        p <- p +
          geom_text(
            data = labelData, aes(x = x, y = y, label = label),
            size = 3, color = "black", fontface = "plain"
          )
      }
    }
    
    p
  }
)

