library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)

library(miniUI)
library(grDevices)
library(RColorBrewer)
library(grid)

# sce <- Load10X_Spatial("/Volumes/neurogen/ASAP/ST/ftp.broadinstitute.org/novaseq_210727/a/A/outs")
sce <- Load10X_Spatial("/Volumes/NEUROGEN/ASAP/ST/ftp.broadinstitute.org/novaseq_230114/d/D/outs")
sample_name <- "novaseq_230114_D"
outfile_dir <- "/Volumes/neurogen/ASAP/ST/data/fold_removal/"
sce@meta.data["orig.ident"] = sample_name
sce@meta.data["selected_spot"] = FALSE

#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#'
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

# Feature plot palettes
#
FeaturePalettes <- list(
  'Spatial' = SpatialColors(n = 100),
  'Seurat' = c('lightgrey', 'blue')
)



# Update slots in an object
#
# @param object An object to update
#
# @return \code{object} with the latest slot definitions
#
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(mode = class(x = xobj), length = 1L)
    }
  }
  return(object)
}

DefaultImage <- function(object) {
  object <- UpdateSlots(object = object)
  images <- Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) < 1) {
    images <- Images(object = object)
  }
  return(images[[1]])
}





# Base plotting function for all Spatial plots
#
# @param data Data.frame with info to be plotted
# @param image SpatialImage object to be plotted
# @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
# or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
# By default, ggplot2 assigns colors
# @param image.alpha Adjust the opacity of the background images. Set to 0 to
# remove.
# @param pt.alpha Adjust the opacity of the points if plotting a SpatialDimPlot
# @param crop Crop the plot in to focus on points plotted. Set to FALSE to show
# entire background image.
# @param pt.size.factor Sets the size of the points relative to spot.radius
# @param stroke Control the width of the border around the spots
# @param col.by Mapping variable for the point color
# @param alpha.by Mapping variable for the point alpha value
# @param cells.highlight A list of character or numeric vectors of cells to
# highlight. If only one group of cells desired, can simply pass a vector
# instead of a list. If set, colors selected cells to the color(s) in
# cols.highlight
# @param cols.highlight A vector of colors to highlight the cells as; ordered
# the same as the groups in cells.highlight; last color corresponds to
# unselected cells.
# @param geom Switch between normal spatial geom and geom to enable hover
# functionality
# @param na.value Color for spots with NA values

#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes_string coord_fixed geom_point xlim ylim
#' coord_cartesian labs theme_void theme scale_fill_brewer
#'
SingleSpatialPlot <- function(
  data,
  image,
  cols = NULL,
  image.alpha = 1,
  pt.alpha = NULL,
  crop = TRUE,
  pt.size.factor = NULL,
  stroke = 0.25,
  col.by = NULL,
  alpha.by = NULL,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  geom = c('spatial', 'interactive', 'poly'),
  na.value = 'grey50'
) {
  geom <- match.arg(arg = geom)
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size.factor,
      cols.highlight = cols.highlight[1],
      col.base = cols.highlight[2]
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    levels(x = data$ident) <- c(order, setdiff(x = levels(x = data$ident), y = order))
    data <- data[order(data$ident), ]
  }
  plot <- ggplot(data = data, aes_string(
    x = colnames(x = data)[2],
    y = colnames(x = data)[1],
    fill = col.by,
    alpha = alpha.by
  ))
  plot <- switch(
    EXPR = geom,
    'spatial' = {
      if (is.null(x = pt.alpha)) {
        plot <- plot + geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          crop = crop,
          stroke = stroke,
        )
      } else {
        plot <- plot + geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          crop = crop,
          stroke = stroke,
          alpha = pt.alpha
        )
      }
      plot + coord_fixed() + theme(aspect.ratio = 1)
    },
    'interactive' = {
      plot + geom_spatial_interactive(
        data = tibble(grob = list(GetImage(object = image, mode = 'grob'))),
        mapping = aes_string(grob = 'grob'),
        x = 0.5,
        y = 0.5
      ) +
        geom_point(mapping = aes_string(color = col.by)) +
        xlim(0, ncol(x = image)) +
        ylim(nrow(x = image), 0) +
        coord_cartesian(expand = FALSE)
    },
    'poly' = {
      data$cell <- rownames(x = data)
      data[, c('x', 'y')] <- NULL
      data <- merge(
        x = data,
        y = GetTissueCoordinates(object = image, qhulls = TRUE),
        by = "cell"
      )
      plot + geom_polygon(
        data = data,
        mapping = aes_string(fill = col.by, group = 'cell')
      ) + coord_fixed() + theme_cowplot()

    },
    stop("Unknown geom, choose from 'spatial' or 'interactive'", call. = FALSE)
  )
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  if (!is.null(x = cols) && is.null(x = cells.highlight)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_fill_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_fill_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}





# For plotting the tissue image
#' @importFrom ggplot2 ggproto Geom aes ggproto_parent alpha draw_key_point
#' @importFrom grid unit gpar editGrob pointsGrob viewport gTree addGrob grobName
#'
GeomSpatial <- ggproto(
  "GeomSpatial",
  Geom,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "image", "image.alpha", "crop"),
  default_aes = aes(
    shape = 21,
    colour = "black",
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = 0.25
  ),
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    # We need to flip the image as the Y coordinates are reversed
    data$y = max(data$y) - data$y + min(data$y)
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, image.alpha, crop) {
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    if (!crop) {
      y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, ncol(x = image))
      panel_scales$y$continuous_range <- c(0, nrow(x = image))
      panel_scales$y.range <- c(0, nrow(x = image))
      panel_scales$x.range <- c(0, ncol(x = image))
    }
    z <- coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    # Flip Y axis for image
    z$y <- -rev(z$y) + 1
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)

    img <- editGrob(grob = img.grob, vp = vp)
    # spot.size <- slot(object = image, name = "spot.radius")
    spot.size <- Radius(object = image)
    coords <- coord$transform(data, panel_scales)
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = data$shape,
      size = unit(spot.size, "npc") * data$point.size.factor,
      gp = gpar(
        col = alpha(colour = coords$colour, alpha = coords$alpha),
        fill = alpha(colour = coords$fill, alpha = coords$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    if (image.alpha > 0) {
      if (image.alpha != 1) {
        img$raster = as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(x = img$raster),
            ncol = ncol(x = img$raster),
            byrow = TRUE)
        )
      }
      gt <- addGrob(gTree = gt, child = img)
    }
    gt <- addGrob(gTree = gt, child = pts)
    # Replacement for ggname
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
    return(gt)
    # ggplot2:::ggname("geom_spatial", gt)
  }
)


# influenced by: https://stackoverflow.com/questions/49475201/adding-tables-to-ggplot2-with-facet-wrap-in-r
# https://ggplot2.tidyverse.org/articles/extending-ggplot2.html
#' @importFrom ggplot2 layer
#'
#'
geom_spatial <-  function(
  mapping = NULL,
  data = NULL,
  image = image,
  image.alpha = image.alpha,
  crop = crop,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
) {
  layer(
    geom = GeomSpatial,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, image.alpha = image.alpha, crop = crop, ...)
  )
}





ISpatialFeaturePlot <- function(
  object,
  feature,
  image = NULL,
  slot = 'data',
  # alpha = c(0.1, 1)
  alpha = c(1.0, 1)
) {
  # Set inital data values
  assay.keys <- Key(object = object)[Assays(object = object)]
  keyed <- sapply(X = assay.keys, FUN = grepl, x = feature)
  assay <- if (any(keyed)) {
    names(x = which(x = keyed))[1]
  } else {
    DefaultAssay(object = object)
  }
  
  # which features to show? currently, only the loaded one
  # features <- sort(x = rownames(x = GetAssayData(
  #   object = object,
  #   slot = slot,
  #   assay = assay
  # )))
  features <- c(feature)
  feature.label <- 'Feature to visualize'
  assays.use <- vapply(
    X = Assays(object = object),
    FUN = function(x) {
      return(!IsMatrixEmpty(x = GetAssayData(
        object = object,
        slot = slot,
        assay = x
      )))
    },
    FUN.VALUE = logical(length = 1L)
  )
  assays.use <- sort(x = Assays(object = object)[assays.use])
  # Setup gadget UI
  ui <- miniPage(
    miniButtonBlock(miniTitleBarButton(
      inputId = 'done',
      label = 'Done',
      primary = TRUE
    )),
    miniContentPanel(
      fillRow(
        sidebarPanel(
          sliderInput(
            inputId = 'alpha',
            label = 'Alpha intensity',
            min = 0,
            max = max(alpha),
            value = min(alpha),
            step = 0.01,
            width = '100%'
          ),
          sliderInput(
            inputId = 'pt.size',
            label = 'Point size',
            min = 0,
            max = 5,
            value = 1.6,
            step = 0.1,
            width = '100%'
          ),
          selectInput(
            inputId = 'assay',
            label = 'Assay',
            choices = assays.use,
            selected = assay,
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'feature',
            label = feature.label,
            choices = features,
            selected = feature,
            selectize = FALSE,
            width = '100%'
          ),
          selectInput(
            inputId = 'palette',
            label = 'Color scheme',
            choices = names(x = FeaturePalettes),
            selected = 'Spatial',
            selectize = FALSE,
            width = '100%'
          ),
          actionButton("plot_reset", "Reset Spot Filter"),
          actionButton("do", "Download Processed Data"),
          width = '100%'
        ),
        # plotOutput(outputId = 'plot', height = '100%', click="plot_click",
        #            brush = "plot_brush", dblclick = "plot_reset"),
        plotOutput(outputId = 'plot', height = '100%', dblclick="plot_click",
                   brush = "plot_brush"),
        flex = c(1, 4)
      )
    ),
    verbatimTextOutput("info"),
    # tableOutput("highlight_spots")
  )
  # Prepare plotting data
  image <- image %||% DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  coords <- GetTissueCoordinates(object = object[[image]])
  feature.data <- FetchData(
    object = object,
    vars = feature,
    cells = cells.use,
    slot = slot
  )
  plot.data <- cbind(coords, feature.data)
  server <- function(input, output, session) {




    plot.env <- reactiveValues(
      data = plot.data,
      feature = feature,
      palette = 'Spatial'
    )
    # Observe events
    observeEvent(
      eventExpr = input$done,
      handlerExpr = stopApp(returnValue = plot.env$plot)
    )
    observe(x = {
      assay <- input$assay
      feature.use <- input$feature
      features.assay <- c(feature)
      # features.assay <- sort(x = rownames(x = GetAssayData(
      #   object = object,
      #   slot = slot,
      #   assay = assay
      # )))
      feature.use <- ifelse(
        test = feature.use %in% features.assay,
        yes = feature.use,
        no = features.assay[1]
      )
      updateSelectInput(
        session = session,
        inputId = 'assay',
        label = 'Assay',
        choices = assays.use,
        selected = assay
      )
      updateSelectInput(
        session = session,
        inputId = 'feature',
        label = feature.label,
        choices = features.assay,
        selected = feature.use
      )
    })
    observe(x = {
      feature.use <- input$feature
      try(
        expr = {
          feature.data <- FetchData(
            object = object,
            vars = paste0(Key(object = object[[input$assay]]), feature.use),
            cells = cells.use,
            slot = slot
          )
          colnames(x = feature.data) <- feature.use
          plot.env$data <- cbind(coords, feature.data)
          plot.env$feature <- feature.use
        },
        silent = TRUE
      )
    })
    observe(x = {
      plot.env$palette <- input$palette
    })
    # Create plot

    # Added for brush selection
    # for now - create a copy of plot.data and then use this for annotation. Is there a better way to do this?
    plot.data_interact <- data.frame(plot.data)
    plot.data_interact$imagerow <- max(plot.data_interact$imagerow) -
                                    plot.data_interact$imagerow +
                                    min(plot.data_interact$imagerow)
    showPlot <- reactiveVal(1)

    # selected <- reactiveVal(rep(FALSE, nrow(plot.data)))
    selected <- reactiveVal(rep(FALSE, nrow(plot.data_interact)))
    observeEvent(input$plot_brush, {
      # brushed <- brushedPoints(plot.data, input$plot_brush, allRows = TRUE)$selected_
      brushed <- brushedPoints(plot.data_interact, input$plot_brush, allRows = TRUE)$selected_
      # print the selected columns
      print("---")
      print(plot.data_interact[brushed,])
      # print(head(sce@meta.data))
      selected(brushed | selected())
    })
    observeEvent(input$plot_reset, {
      selected(rep(FALSE, nrow(plot.data_interact)))
    })


    output$plot <- renderPlot(expr = {
      # plot.data$imagerow_flip <- max(plot.data$imagerow) - plot.data$imagerow + min(plot.data$imagerow)
      plot.data_interact$selected_spot <- selected()
      plot.env$plot <- SingleSpatialPlot(
        data = plot.env$data,
        image = object[[image]],
        # col.by = plot.env$feature,
        col.by = NULL,
        pt.size.factor = input$pt.size,
        crop = TRUE,
        alpha.by = plot.env$feature
      ) + geom_point(plot.data_interact,
                               mapping=aes(colour = selected_spot, alpha=showPlot()), size=2.0) +
        # scale_fill_gradientn(name = plot.env$feature, colours = cols) +
        scale_fill_gradientn(name = plot.env$feature, colours = FeaturePalettes[[plot.env$palette]]) +
        theme(legend.position = 'top') +
        scale_alpha(range = c(input$alpha, 1)) +
        # guides(alpha = FALSE)
        guides(alpha = "none")
      plot.env$plot
    })

    # brush points in table
    # output$highlight_spots <- renderTable({
    # brushedPoints(plot.data_interact, input$plot_brush)
    # })

    # click input
    # output$info <- renderPrint({
    #   req(input$plot_click)
    #   x <- round(input$plot_click$x, 2)
    #   y <- round(input$plot_click$y, 2)
    #   cat("[", x, ", ", y, "]", sep = "")
    #   print(x)
    #   print(y)
    # })
    observeEvent(input$plot_click, {
      # brushed <- brushedPoints(plot.data, input$plot_brush, allRows = TRUE)$selected_
      nearpts <- nearPoints(plot.data_interact, input$plot_click, allRows = TRUE)$selected_
      # print the selected columns
      print(plot.data_interact[nearpts,])
      # selected(nearpts | selected())
      selected(xor(nearpts, selected()))
    })

    # button click
    # this downloads a new Seurat object with the fold removed
    observeEvent(input$do, {
      print("download button pressed")
      # print(showPlot)
      # newshowPlot <- 1 - showPlot()
      # showPlot(newshowPlot)
      # print(showPlot)
      sce@meta.data["selected_spot"] = selected()
      print(head(sce@meta.data))
      
      # write the new Seurat object
      print(paste("writing to ", outfile_dir, " ...", sep=""))
      saveRDS(sce, paste(outfile_dir, "st_seurat_rmfold_", sample_name, ".rds", sep=""))
      write.csv(sce@meta.data, paste(outfile_dir, "st_seurat_rmfold_meta_", sample_name, ".csv", sep=""), row.names = TRUE)
      print("Done!")
    })
  }
  runGadget(app = ui, server = server)
}



ISpatialFeaturePlot(sce, feature="nCount_Spatial")
# ISpatialFeaturePlot(sce, feature="selected_spot")
