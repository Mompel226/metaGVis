#' @title 3D plots
#'
#' @description
#' This function is used to create a MultiDimensional Scaling ordination plot in 3D with plotly.
#'
#' @param physeq (Required) A phyloseq object.
#' @param explanatory_variable (Required) A character string of the metadata variable that is going to be used to group the samples.
#' @param distance A character string defining the distance used to calculate the distance matrix:
#' It can be any distance accepted by the distance function of the phyloseq package and also "PhilR".
#' To use the Aitchison distance, first transform the otu_table of the phyloseq object with one of the CLR functions
#' (clr_min_percent, clr_0_plus_1 or clr_runif) and then use the "euclidean" distance.
#' Default: "bray"
#' @param method A character string defining the ordination method: "PCA", "PCoA", "RDA", "CAP", "CA", "CCA" or "NMDS"
#' Default: "PCoA"
#' @param colors (Optional) A vector with the hexadecimal colors associated to the levels of the explanatory variable.
#' @param constrained A logical value. If TRUE a constrained ordination is carried out.
#' Default: FALSE
#' @param formula (Optional) The right hand side of the model formula -- the constraining variables and conditioning variables. e.g: "~ constraining variable"
#' Only used if constrained = TRUE
#'
#' @return A plotly object.
#'
#' @examples
#' Unconstrained ordination plot:
#' plot_3D <- MDSin3D(physeq = phyloseq.object, explanatory_variable = "Sex", colors = c("#F8766D", "#00C064"))
#'
#' Constrained ordination plot:
#' plot_3D <- MDSin3D(physeq = phyloseq.object, explanatory_variable = "Sex", constrained = TRUE, formula = "~ Country")
#' @export
MDSin3D <- function(physeq,
                    explanatory_variable,
                    distance = "bray",
                    method = "PCoA",
                    ordination.object = NULL,
                    colors,
                    constrained = FALSE,
                    display.ellipses = TRUE,
                    formula){

  # Get level names of the explanatory variable and remove spaces and parenthesis
  level.names <- levels(sample_data(physeq)[[explanatory_variable]])
  clean_levels <- gsub(" ", "_", gsub("\\(", "", gsub("\\)", "", level.names)))

  if(missing(colors)){
    colors <- unname(Polychrome::createPalette(length(level.names), c("#ff0000", "#00ff00", "#0000ff")))
  }

  ord_labels <- function(ord){
    if (class(ord)[1] == "pcoa"){
      ev <- ord$values[,1]
    } else {
      ev <- vegan::eigenvals(ord)
    }

      if (!is.na(ev)[1]) {
        tol <- -(1e-07)*ev[1]
        ord.labels <- rep("", length(ev))
        if ((any(is.na(ev))) | (any(ev < tol))) {
          for ( i in 1:length(ev)) {
            ord.labels[i] <- paste("PCo", i, sep = "")
          }
        }
        else {
          ev.pc <- round(100*(ev/sum(ev)), 2)
          axis.names <- names(ev)
          if (is.null(axis.names)) {
            for ( i in 1:length(ev.pc)) {
              ord.labels[i] <- paste("PCo", i, " ", sprintf(ev.pc[i], fmt = '%#.1f'), "%", sep="")
            }
          } else {
            for (i in 1:length(ev.pc)){
              ord.labels[i] <- paste(axis.names[i], " ", ev.pc[i],"%", sep="")
            }
          }
        }
      } else {
        ord.labels <- colnames(vegan::scores(ord))
      }

      return(ord.labels)
    }

  # Calculate dissimilarity matrix
  if (is.null(ordination.object)){
    if (distance == "PhilR"){
      physeq@otu_table <- physeq@otu_table + 1
      philr <- philr::philr(t(as(physeq@otu_table, "matrix")), physeq@phy_tree,
                            part.weights='enorm.x.gm.counts',
                            ilr.weights='blw.sqrt')
      dist <- dist(philr, method = "euclidean")
    } else {
      dist <- phyloseq::distance(physeq = physeq, method = distance)
    }

    # Ordination analysis
    if (isTRUE(constrained) && !method %in% c("PCA", "NMDS")){
      cap_ord <- ordinate(
        physeq = physeq,
        method = method,
        distance = dist,
        formula = as.formula(formula))
    } else if (!isTRUE(constrained) && !method %in% c("PCA", "NMDS")){
      cap_ord <- ordinate(
        physeq = physeq,
        method = method,
        distance = dist,
        formula = ~ 1)
    } else if (method == "PCA"){
      # With euclidean distance visual results are the same as when choosing methods:
      # "RDA" -- vegan::rda() or "MDS" -- ape:pcoa() OR stats::prcomp() and labdsv::pca()
      cap_ord <- vegan::rda(t(as(physeq@otu_table, "matrix")))
    } else if (method == "NMDS"){
      cap_ord <- vegan::metaMDS(dist, k = 3)
    }
  } else {
    cap_ord <- ordination.object
  }

  # Create a data frame with coordinates and related metadata
  if (method %in% c("MDS", "PCoA")){
    df <- as.data.frame(cap_ord$vectors[,1:3])
  } else {
    df <- as.data.frame(vegan::scores(cap_ord, display = "sites", choices = c(1, 2, 3)))
  }
  colnames(df) <- c("AXISx", "AXISy", "AXISz")
  df$Group <- sample_data(physeq)[, explanatory_variable]
  df$Group <- as.factor(as.matrix(df$Group))
  df$Group <- factor(df$Group, levels = level.names)

  # Calculate ellipse centroids
  mean.data <- df
  mean.data$Group <- sapply(mean.data[, "Group"], as, "numeric") #Coerce Group column to numeric
  mean.data <- aggregate( mean.data, list(df$Group), mean) #Axis mean per group
  rownames(mean.data) <- mean.data$Group.1 #Rename rownames
  mean.data <- mean.data[, -1, drop=FALSE] #"pop" the first column

  # Extract the axis information from the ordination analysis
  labels <- ord_labels(cap_ord)

  # Create 3D ellipses using the ordination analysis data and add the information to the data frame
  k <- 1
  ellip_df.list <- list()
  for (i in level.names){
    group <- df %>% filter(Group == i)
    covariance <- cov(group[1:3])
    covariance <- lqmm::make.positive.definite(covariance)
    ellip <- rgl::ellipse3d(covariance, centre = as.vector((mean.data[,1:3])[match(i, level.names),], mode = "numeric"), subdivide = 3)
    ellip_df <- ellip$vb[1:3,]
    ellip_df <- cbind(data.frame(Coord = c(paste0(clean_levels[k],"_x"),paste0(clean_levels[k],"_y"),paste0(clean_levels[k],"_z"))), ellip_df)
    ellip_df <- cbind(data.frame(Group = rep(level.names[k], 3)), ellip_df)
    ellip_df <- cbind(data.frame(Color= rep(colors[k], 3)), ellip_df)
    ellip_df.list[[i]] <- ellip_df
    k <- k + 1
  }
  df$SampleID <- rownames(df)
  ellip_df <- bind_rows(ellip_df.list)
  DF <- merge(df, ellip_df, by = "Group")


  # Create 3D graph with plot_ly()
  if (isTRUE(constrained)){
    title = paste0(method, " (constrained) ordination plot")
  } else {
    title = paste0(method, " (unconstrained) ordination plot")
  }

  fig <- plotly::plot_ly() %>%
         plotly::add_trace(data = DF, type = "scatter3d", mode = 'markers',
                           x = ~AXISx, y = ~AXISy, z = ~AXISz,
                           text = ~SampleID,
                           color = ~Group,
                           hovertemplate = "<b>%{text}</b>",
                           colors = colors,
                           marker = list(symbol = 'circle',
                                         size = 8,
                                         opacity = 0.7,
                                         line = list(color = 'rgba(255,255,255, 0.9)', width = 1))) %>%
        plotly::layout(title = title,
                       scene = list(
                         xaxis = list(title = labels[1],
                                     showline= T,
                                     linewidth=2,
                                     linecolor='black',
                                     mirror = T,
                                     gridcolor = '#E2E2E2',
                                     zerolinecolor = 'black',
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwidth = 2),
                        yaxis = list(title = labels[2],
                                     showline= T,
                                     linewidth=2,
                                     linecolor='black',
                                     mirror = T,
                                     gridcolor = '#E2E2E2',
                                     zerolinecolor = 'black',
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwith = 2),
                        zaxis = list(title = labels[3],
                                     showline= T,
                                     linewidth=2,
                                     linecolor='black',
                                     mirror = T,
                                     gridcolor = '#E2E2E2',
                                     zerolinecolor = 'black',
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwith = 2),
                        paper_bgcolor = 'rgb(243, 243, 243)',
                        plot_bgcolor = 'rgb(243, 243, 243)')
                        )

  if (display.ellipses){
    # Add the 3D ellipses
    for (i in clean_levels){

      x <- unname(unique(DF[DF$Coord %like% paste0(i,"_x"), 8:ncol(DF)]))
      y <- unname(unique(DF[DF$Coord %like% paste0(i,"_y"), 8:ncol(DF)]))
      z <- unname(unique(DF[DF$Coord %like% paste0(i,"_z"), 8:ncol(DF)]))

      fig <- fig %>% plotly::add_trace(type = 'mesh3d', alphahull = 0,
                                       x = x,
                                       y = y,
                                       z = z,
                                       opacity = 0.1,
                                       showlegend = FALSE)
    }

    # Match sample colors with ellipses colors
    decomp <- plotly::plotly_build(fig)
    k <- 1
    for (i in seq(from = length(level.names) + 1, to = length(level.names)*2, by = 1)){

      decomp$x$data[[i]]$color <- paste0("rgba(", paste0(as.vector(col2rgb(colors[k])), collapse = ","), ",1)")
      k <- k + 1
    }
  } else {
    decomp <- plotly::plotly_build(fig)
  }

  return(decomp)
}

####################################################################################################################################

#' @title 3D plots
#'
#' @description
#' This function is used to create an animation from a 3D ordination plot created with plotly.
#'
#' @param plotly_object A phyloseq object.
#' @param title A character string defining the title of the plot.
#' @param path_to_conda A character string defining the PATH to the bin folder of the Conda package manager
#'
#' @return A gif file.
#'
#' @examples
#' Animated_MDS(plotly_object = plot_3D, title = "PCoA_bray", path_to_conda = "/Users/John/anaconda3/bin/")
#'
#' @export
Animated_MDS <- function(plotly_object, title, path_to_conda = "/Users/daniel/anaconda3/bin/"){

  Sys.setenv("PATH" = paste(Sys.getenv("PATH"), path_to_conda, sep = .Platform$path.sep))

  dir.create(paste0("Animations/Fotograms/", title), recursive = TRUE)

  # Horizontal movement to the right
  for(i in seq(1.6,3,by=0.05)){
    outfile <- paste0("Animations/Fotograms/", title, "/", title, ".1.", format(round(i, digits=2), nsmall = 2), ".png")
    cam.zoom = 2.7
    ver.angle = 0
    graph <- plotly_object %>% plotly::layout(scene = list(camera = list(eye = list(x = cos(i)*cam.zoom,y = sin(i)*cam.zoom, z=0.2),
                                                                                 center = list(x = 0,
                                                                                               y = 0,
                                                                                               z = 0))))
    plotly::orca(graph, file = outfile)
  }
  # Vertical movement to the top
  for(i in seq(0.2,1.9,by=0.05)){
    outfile <- paste0("Animations/Fotograms/", title, "/", title, ".2.", format(round(i, digits=3), nsmall = 3), ".png")
    cam.zoom = 2.7
    ver.angle = 0
    graph <- plotly_object %>% plotly::layout(scene = list(camera = list(eye = list(x = cos(3)*cam.zoom,y = sin(3)*cam.zoom, z=i),
                                                                         center = list(x = 0,
                                                                                       y = 0,
                                                                                       z = 0))))
    plotly::orca(graph, file = outfile)
  }
  # Horizontal movement to the left
  k <- 1
  for(i in rev(seq(1.6, 3, by = 0.05))){
    outfile <- paste0("Animations/Fotograms/", title, "/", title, ".3.", sprintf("%02d", k), ".png")
    cam.zoom = 2.7
    ver.angle = 0
    graph <- plotly_object %>% plotly::layout(scene = list(camera = list(eye = list(x = cos(i)*cam.zoom,y = sin(i)*cam.zoom, z=2),
                                                                         center = list(x = 0,
                                                                                       y = 0,
                                                                                       z = 0))))
    plotly::orca(graph, file = outfile)
    k <- k + 1
  }
  k <- 1
  # Vertical movement to the bottom
  for(i in rev(seq(0.2, 1.9, by=0.05))){
    outfile <- paste0("Animations/Fotograms/", title, "/", title, ".4.", sprintf("%02d", k), ".png")
    cam.zoom = 2.7
    ver.angle = 0
    graph <- plotly_object %>% plotly::layout(scene = list(camera = list(eye = list(x = cos(1.6)*cam.zoom,y = sin(1.6)*cam.zoom, z=i),
                                                                         center = list(x = 0,
                                                                                       y = 0,
                                                                                       z = 0))))
    plotly::orca(graph, file = outfile)
    k <- k + 1
  }

  ## list file names and read in
  imgs <- list.files(paste0("Animations/Fotograms/", title, "/"), full.names = TRUE)
  img_list <- lapply(imgs, magick::image_read)

  ## join the images together
  img_joined <- magick::image_join(img_list)

  ## animate at 2 frames per second
  img_animated <- magick::image_animate(img_joined, fps = 10)

  ## save to disk
  magick::image_write(image = img_animated, path = paste0("Animations/", title,".gif"))
}
