
#' Better visualize feature on 2D dimensional reduction plot
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
#' @importFrom circlize colorRamp2
#' @importFrom forcats fct_reorder
#' @importFrom purrr set_names
#' @import ggplot2
#' @import ggrastr
#' @import RColorBrewer
#' @import data.table
#' @import patchwork
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
#' @param label *logical* describing whether to label cell clusters.
#' @param label_meta *string* describing which metadata is used to label cell clusters.
#' @param label_max *integer* describing the maximum number of different labels.
#' @param vln *logical* describing whether to add violin plot.
#' @param vln_meta *string* describing which metadata is used for violin plot.
#' @param vln_max *integer* describing the maximum number of categories of violin plot.
#' @param split *logical* describing whether to split subplots.
#' @param split_meta *string* describing which metadata is used for split subplots.
#' @param split_max *integer* describing the maximum number of categories of split subplots.
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
    size = if(ncol(obj) < 5000) {1} else {0.4}, alpha = 1,
    do_raster = F, dpi = 300L,
    x_name = paste0(toupper(reduc), "_", dims[1]), y_name = paste0(toupper(reduc), "_", dims[2]),
    label = T, label_meta = "orig.ident", label_max = 50L,
    vln = T, vln_meta = label_meta, vln_max = 50L,
    split = F, split_meta = "orig.ident", split_max = 10L
  ) {
    
    stopifnot("Parameter 'dims' must be 2 different whole numbers!" = length(dims) == 2 && dims[1] != dims[2] && all.equal(dims, as.integer(dims)))
    plotData <- as.data.table(obj[[reduc]]@cell.embeddings[, dims])
    colnames(plotData) <- c("dim1", "dim2")
    
    # get feature
    obj_version <- Version(obj)$major
    feat_type <- match.arg(feat_type)
    stopifnot("Parameter 'feat' must be 1 single string." = length(feat) == 1)
    if(feat_type == "gene") {
      if(obj_version >= 5) {
        mat <- obj[[assay]][slot]
      } else {
        mat <- obj[[assay]][[slot]]
      }
      stopifnot("No such gene in this assay!" = feat %in% rownames(mat))
      plotData$feat <- mat[feat, ]
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
        limits = range(x),
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
        limits = pn,
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
      plot.title = element_text(face = "bold.italic", hjust = .5),
      strip.background = element_blank()
    )
    
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
    if(order == "value") {
      plotDataOrdered <- plotData[order(feat)]
    } else if(order == "abs") {
      plotDataOrdered <- plotData[order(abs(feat))]
    } else if(order == "given") {
      plotDataOrdered <- plotData[order_by]
    } else if(order == "random") {
      plotDataOrdered <- plotData[sample(1:.N)]
    }
    
    # main plot
    if(do_raster) {
      if(!"split" %in% colnames(plotData)) {
        p <- ggplot(plotDataOrdered, aes(x = dim1, y = dim2, color = feat)) +
          geom_point_rast(shape = 16, size = size, alpha = alpha, raster.dpi = dpi) +
          labs(x = x_name, y = y_name, title = feat, color = "") +
          plotColorBar +
          plotTheme
      } else {
        pList <- list()
        for(i in unique(plotDataOrdered$split)) {
          pList[[i]] <- ggplot(plotDataOrdered[split == i], aes(x = dim1, y = dim2, color = feat)) +
            geom_point_rast(shape = 16, size = size, alpha = alpha, raster.dpi = dpi) +
            labs(x = x_name, y = y_name, title = feat, color = "") +
            scale_x_continuous(limits = range(plotDataOrdered$dim1)) +
            scale_y_continuous(limits = range(plotDataOrdered$dim2)) +
            plotColorBar +
            plotTheme
        }
      }
    } else {
      if(!"split" %in% colnames(plotData)) {
        p <- ggplot(plotDataOrdered, aes(x = dim1, y = dim2, color = feat)) +
          geom_point(shape = 16, size = size, alpha = alpha) +
          labs(x = x_name, y = y_name, title = feat, color = "") +
          plotColorBar +
          plotTheme
      } else {
        pList <- list()
        for(i in unique(plotDataOrdered$split)) {
          pList[[i]] <- ggplot(plotDataOrdered[split == i], aes(x = dim1, y = dim2, color = feat)) +
            geom_point(shape = 16, size = size, alpha = alpha) +
            labs(x = x_name, y = y_name, title = i, color = "") +
            scale_x_continuous(limits = range(plotDataOrdered$dim1)) +
            scale_y_continuous(limits = range(plotDataOrdered$dim2)) +
            plotColorBar +
            plotTheme
        }
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
        } else {
          labelData <- plotData[, .(x = median(dim1), y = median(dim2)), by = .(label, split)]
          
          for(i in names(pList)) {
            pList[[i]] <- pList[[i]] +
              geom_label(
                data = labelData[split == i], aes(x = x, y = y, label = label),
                fill = "#FFFFFF", size = 3, color = NA, alpha = .6, fontface = "bold",
                label.size = 0.5, label.r = unit(0.25, "lines"), label.padding = unit(0.15, "lines")
              ) +
              geom_text(
                data = labelData[split == i], aes(x = x, y = y, label = label),
                size = 3, color = "black", fontface = "bold"
              )
          }
        }
      }
    }
    
    # add violin plot
    if(vln) {
      # check violin metadata
      if(!vln_meta %in% colnames(obj@meta.data)) {
        warning("Parameter 'vln_meta' should be one of metadata names.", call. = F)
      } else if(length(unique(obj@meta.data[[vln_meta]])) > vln_max) {
        warning("Too many categories, should find a simpler 'vln_meta'.", call. = F)
      } else {
        plotData$meta <- obj@meta.data[, vln_meta]
        
        if(!"split" %in% colnames(plotData)) {
          vlnData <- plotData[, .(.N, avg = mean(feat)), by = meta]
          vlnData <- vlnData[N >= 5]
          
          x <- range(plotData$feat)
          colFun <- colorRamp2(x[1] + diff(x)*disp, col)
          vlnCol <- set_names(colFun(vlnData$avg), vlnData$meta)
          
          v <- ggplot(plotData[meta %in% vlnData$meta], aes(x = fct_reorder(meta, feat, .desc = T), y = feat)) +
            geom_violin(aes(fill = meta), scale = "width", adjust = 1.5, show.legend = F) +
            scale_fill_manual(values = vlnCol) +
            scale_y_continuous(expand = expansion(c(0, 0.02))) +
            scale_x_discrete(position = "top") +
            theme(
              aspect.ratio = 1/3,
              panel.background = element_blank(),
              panel.border = element_rect(fill = NA, color = "black"),
              plot.background = element_blank(),
              panel.grid = element_blank(),
              panel.grid.major.x = element_line(color = "gray80"),
              axis.text.x = element_text(angle = 45, hjust = 0),
              axis.title = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank()
            )
          
          v + p + theme(plot.title = element_blank()) +
            plot_layout(ncol = 1, guides = "collect") +
            plot_annotation(
              title = feat,
              theme = theme(plot.title = element_text(face = "bold.italic", hjust = .5))
            )
        } else {
          vlnData <- plotData[, .(.N, avg = mean(feat)), by = .(meta, split)]
          vlnData <- vlnData[N >= 5]
          
          x <- range(plotData$feat)
          colFun <- colorRamp2(x[1] + diff(x)*disp, col)
          
          vList <- list()
          for(i in names(pList)) {
            vlnData_i <- vlnData[split == i]
            vlnCol_i <- set_names(colFun(vlnData_i$avg), vlnData_i$meta)
            
            vList[[i]] <- ggplot(plotData[meta %in% vlnData_i$meta], aes(x = fct_reorder(meta, feat, .desc = T), y = feat)) +
              geom_violin(aes(fill = meta), scale = "width", adjust = 1.5, show.legend = F) +
              scale_fill_manual(values = vlnCol_i) +
              scale_y_continuous(expand = expansion(c(0, 0.02))) +
              scale_x_discrete(position = "top") +
              theme(
                aspect.ratio = 1/3,
                panel.background = element_blank(),
                panel.border = element_rect(fill = NA, color = "black"),
                plot.background = element_blank(),
                panel.grid = element_blank(),
                panel.grid.major.x = element_line(color = "gray80"),
                axis.text.x = element_text(angle = 45, hjust = 0),
                axis.title = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank()
              )
          }
          
          z <- vList[[1]]
          for(i in 2:length(vList)) {
            z <- z + vList[[i]]
          }
          for(i in 1:length(pList)) {
            z <- z + pList[[i]]
          }
          z + plot_layout(nrow = 2, guides = "collect") +
            plot_annotation(
              title = feat,
              theme = theme(plot.title = element_text(face = "bold.italic", hjust = .5))
            )
        }
      }
    } else {
      if(!"split" %in% colnames(plotData)) {
        p
      } else {
        z <- pList[[1]]
        for(i in 2:length(pList)) {
          z <- z + pList[[i]]
        }
        z + plot_layout(nrow = 1, guides = "collect") +
          plot_annotation(
            title = feat,
            theme = theme(plot.title = element_text(face = "bold.italic", hjust = .5))
          )
      }
    }
  }
)
