# Spatial Transcriptomics Fold Removal App
#
# Interactive Shiny app for manually selecting and removing tissue fold
# artifacts from 10x Visium spatial transcriptomics data. Users can click spots
# overlapping folds or click and drag a box, and the app saves a Seurat object
# object and metadata with the annotation to disk.
#
# Usage:
#   Rscript app.R --datapath <10x_outs_dir> --sample <sample_name> --outdir <output_dir>
#
# Author: Jie Yuan
#
# Attribution:
# This code was adapted originally from the Seurat R package v4.4.0
# (https://github.com/satijalab/seurat), specifically from its spatial
# transcriptomics visualization module. The original source is licensed under
# MIT License (c) Satija Lab. Modifications were made to support interactive
# spot selection via Shiny for the purpose of tissue fold removal.

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(miniUI)
library(grDevices)
library(RColorBrewer)
library(grid)
library(optparse)

option_list <- list(
  make_option(c("-d", "--datapath"), type = "character", default = NULL,
              help = "Path to 10x Visium 'outs' directory (required)",
              metavar = "character"),
  make_option(c("-s", "--sample"), type = "character", default = NULL,
              help = "Sample name, used for labeling and output filenames (required)",
              metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "./",
              help = "Output directory for processed Seurat object and metadata CSV [default: ./]",
              metavar = "character"),
  make_option(c("-f", "--feature"), type = "character", default = "nCount_Spatial",
              help = "Seurat feature to display in the app [default: nCount_Spatial]",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$datapath)) {
  print_help(opt_parser)
  stop("--datapath is required.", call. = FALSE)
}
if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("--sample is required.", call. = FALSE)
}

data_path <- opt$datapath
sample_name <- opt$sample
outfile_dir <- opt$outdir
feature_name <- opt$feature

# Ensure output directory exists
if (!dir.exists(outfile_dir)) {
  dir.create(outfile_dir, recursive = TRUE)
}


# Data loading
message("Loading 10x Spatial data from: ", data_path)
sce <- Load10X_Spatial(data_path)

# Tag spots with sample identity; initialize selection column to FALSE
sce@meta.data["orig.ident"]     <- sample_name
sce@meta.data["selected_spot"]  <- FALSE

# Spectral palette used for spatial feature plots (reversed so low = blue, high = red)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

FeaturePalettes <- list(
  'Spatial' = SpatialColors(n = 100),
  'Seurat'  = c('lightgrey', 'blue')
)


# Helper functions

# Refresh all slots of a Seurat object to the current class definition.
# Needed when loading objects created by an older Seurat version.
UpdateSlots <- function(object) {
  object.list <- sapply(
    X   = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr  = slot(object = object, name = x),
        error = function(...) { return(NULL) }
      ))
    },
    simplify  = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)

  # Fill any newly added slots with zero-length defaults
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(mode = class(x = xobj), length = 1L)
    }
  }
  return(object)
}

# Return the name of the default image associated with the active assay.
DefaultImage <- function(object) {
  object <- UpdateSlots(object = object)
  images <- Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) < 1) {
    images <- Images(object = object)
  }
  return(images[[1]])
}

# ggproto object that renders the tissue image as a background grob and
# overlays Visium spots as pch-21 points at their correct spatial positions.
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
    # Visium imagerow increases downward; flip so ggplot y-axis is correct
    data$y <- max(data$y) - data$y + min(data$y)
    data
  },
  draw_key   = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, image.alpha, crop) {

    # When not cropping, expand panel range to the full image dimensions
    if (!crop) {
      y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, ncol(x = image))
      panel_scales$y$continuous_range <- c(0, nrow(x = image))
      panel_scales$y.range <- c(0, nrow(x = image))
      panel_scales$x.range <- c(0, ncol(x = image))
    }

    # Map image corners to normalized panel coordinates
    z <- coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    z$y <- -rev(z$y) + 1   # flip y for image raster orientation
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]

    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth,   units = "npc"),
      height = unit(x = hgth,   units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)
    img <- editGrob(grob = img.grob, vp = vp)

    # Draw spots as points scaled to the Visium spot radius
    spot.size <- Radius(object = image)
    coords <- coord$transform(data, panel_scales)
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = data$shape,
      size = unit(spot.size, "npc") * data$point.size.factor,
      gp = gpar(
        col = alpha(colour = coords$colour, alpha = coords$alpha),
        fill = alpha(colour = coords$fill,   alpha = coords$alpha),
        lwd = coords$stroke
      )
    )

    gt <- gTree(vp = viewport())
    if (image.alpha > 0) {
      # Apply transparency to the tissue image raster if requested
      if (image.alpha != 1) {
        img$raster <- as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(x = img$raster),
            ncol = ncol(x = img$raster),
            byrow = TRUE
          )
        )
      }
      gt <- addGrob(gTree = gt, child = img)
    }
    gt <- addGrob(gTree = gt, child = pts)
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
    return(gt)
  }
)

# ggplot2 layer constructor for GeomSpatial
geom_spatial <- function(
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
    params = list(na.rm = na.rm, image = image,
                  image.alpha = image.alpha, crop = crop, ...)
  )
}

# Core spatial ggplot builder. Handles color/alpha mappings, spot highlighting,
# and palette selection; returns a ggplot object without axes.
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
    warning("Cannot find '", col.by, "' in data, not coloring",
            call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }

  # Wrap column names in backticks to handle names with special characters
  col.by <- col.by   %iff% paste0("`", col.by,   "`")
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
          stroke = stroke
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

  # Apply highlight or custom color scale
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  if (!is.null(x = cols) && is.null(x = cells.highlight)) {
    if (length(x = cols) == 1 &&
        (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_fill_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 &&
               cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped')) {
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


# Shiny gadget
#
# Interactive spatial feature plot with fold-removal selection
#
# Launches a Shiny gadget that displays a spatial feature map and allows the
# user to brush-select or click spots to mark them as belonging to a tissue
# fold. On clicking "Download Processed Data", saves a Seurat RDS and metadata
# CSV with a `selected_spot` column indicating fold membership.
#
# @param object  A Seurat object with spatial data loaded.
# @param feature Character. Feature to plot (must be present in the object).
# @param image   Character. Name of the image slot to use; defaults to the
#                active assay image.
# @param slot    Character. Assay slot to pull feature data from.
# @param alpha   Numeric vector of length 2. Min/max alpha range for the slider.
ISpatialFeaturePlot <- function(
  object,
  feature,
  image = NULL,
  slot  = 'data',
  alpha = c(1.0, 1)
) {

  # Assay / feature setup
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

  # Only list assays that actually contain data in the requested slot
  assays.use <- vapply(
    X = Assays(object = object),
    FUN = function(x) {
      return(!IsMatrixEmpty(x = GetAssayData(object = object, slot = slot, assay = x)))
    },
    FUN.VALUE = logical(length = 1L)
  )
  assays.use <- sort(x = Assays(object = object)[assays.use])


  # UI layout
  ui <- miniPage(
    miniButtonBlock(miniTitleBarButton(
      inputId = 'done',
      label   = 'Done',
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
          actionButton("do",         "Download Processed Data"),
          width = '100%'
        ),
        # Double-click selects/deselects individual spots; brush selects regions
        plotOutput(outputId = 'plot', height = '100%',
                   dblclick = "plot_click",
                   brush = "plot_brush"),
        flex = c(1, 4)
      )
    ),
    verbatimTextOutput("info")
  )

  # Prepare plotting data
  image <- image %||% DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  coords <- GetTissueCoordinates(object = object[[image]])

  feature.data <- FetchData(object = object, vars = feature,
                            cells = cells.use, slot = slot)
  plot.data <- cbind(coords, feature.data)

  server <- function(input, output, session) {

    plot.env <- reactiveValues(
      data = plot.data,
      feature = feature,
      palette = 'Spatial'
    )

    # Close gadget and return the last rendered plot
    observeEvent(input$done, stopApp(returnValue = plot.env$plot))

    # Sync assay/feature dropdowns
    observe({
      assay <- input$assay
      feature.use <- input$feature
      features.assay <- c(feature)
      feature.use <- ifelse(
        test = feature.use %in% features.assay,
        yes = feature.use,
        no = features.assay[1]
      )
      updateSelectInput(
        session,
        'assay',
        choices = assays.use,
        selected = assay
      )
      updateSelectInput(
        session,
        'feature',
        choices = features.assay,
        selected = feature.use
      )
    })

    # Re-fetch feature data when assay/feature selection changes
    observe({
      feature.use <- input$feature
      try({
        feature.data <- FetchData(
          object = object,
          vars = paste0(Key(object = object[[input$assay]]), feature.use),
          cells = cells.use,
          slot = slot
        )
        colnames(x = feature.data) <- feature.use
        plot.env$data <- cbind(coords, feature.data)
        plot.env$feature <- feature.use
      }, silent = TRUE)
    })

    observe({ plot.env$palette <- input$palette })


    # Create plot

    # Added for brush selection
    # for now - create a copy of plot.data and then use this for annotation. Is there a better way to do this?
    plot.data_interact <- data.frame(plot.data)
    plot.data_interact$imagerow <- max(plot.data_interact$imagerow) -
                                   plot.data_interact$imagerow +
                                   min(plot.data_interact$imagerow)
    showPlot <- reactiveVal(1)

    # Boolean vector tracking which spots are currently selected (fold members)
    selected <- reactiveVal(rep(FALSE, nrow(plot.data_interact)))

    # Brush: union of currently selected and newly brushed spots
    observeEvent(input$plot_brush, {
      brushed <- brushedPoints(plot.data_interact, input$plot_brush, allRows = TRUE)$selected_
      print(plot.data_interact[brushed, ])
      selected(brushed | selected())
    })

    # Reset button: deselect all spots
    observeEvent(input$plot_reset, {
      selected(rep(FALSE, nrow(plot.data_interact)))
    })

    # Double-click: toggle individual spot selection
    observeEvent(input$plot_click, {
      nearpts <- nearPoints(plot.data_interact, input$plot_click, allRows = TRUE)$selected_
      print(plot.data_interact[nearpts, ])
      selected(xor(nearpts, selected()))
    })


    #  Render plot

    output$plot <- renderPlot({
      plot.data_interact$selected_spot <- selected()

      plot.env$plot <- SingleSpatialPlot(
        data = plot.env$data,
        image = object[[image]],
        col.by = NULL,          # color driven by feature alpha, not discrete fill
        pt.size.factor = input$pt.size,
        crop = TRUE,
        alpha.by = plot.env$feature
      ) +
        # Overlay a second point layer colored by selection status
        geom_point(plot.data_interact,
                   mapping = aes(colour = selected_spot, alpha = showPlot()),
                   size = 2.0) +
        scale_fill_gradientn(name = plot.env$feature,
                             colours = FeaturePalettes[[plot.env$palette]]) +
        theme(legend.position = 'top') +
        scale_alpha(range = c(input$alpha, 1)) +
        guides(alpha = "none")

      plot.env$plot
    })


    #  Save output

    # Write a Seurat RDS and metadata CSV with selected_spot column populated
    observeEvent(input$do, {
      message("Saving processed data to: ", outfile_dir)
      sce@meta.data["selected_spot"] <- selected()
      print(head(sce@meta.data))

      rds_path <- file.path(outfile_dir,
                            paste0("st_seurat_rmfold_", sample_name, ".rds"))
      csv_path <- file.path(outfile_dir,
                            paste0("st_seurat_rmfold_meta_", sample_name, ".csv"))

      saveRDS(sce, rds_path)
      write.csv(sce@meta.data, csv_path, row.names = TRUE)
      message("Done! Files written:\n  ", rds_path, "\n  ", csv_path)
    })
  }

  runGadget(app = ui, server = server)
}


# Run code
# ISpatialFeaturePlot(sce, feature="nCount_Spatial")
ISpatialFeaturePlot(sce, feature = feature_name)
