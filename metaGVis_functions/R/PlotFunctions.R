#' @title Rarefaction Curves
#'
#' @description
#' This function is used to create a plot that shows the rarefaction curves of all samples.
#' Curves can be colored or have different linetypes depending on the grouping variables specified.
#'
#' @param physeq (Required) A phyloseq object.
#' @param vegandata (Optional) Output from the rarecurve function of the Vegan package.
#' Providing this data reduces the computational time, therefore changing graphical parameters is quickly performed. DEFAULT: NULL
#' @param step (Required) Step size for sample sizes in rarefaction curves. DEFAULT: 500
#' @param col_var (Optional) A character string indicating the metadata variable used to set the color of the lines. DEFAULT: NULL
#' @param colors (Optional) A vector with the hexadecimal colors associated to the levels of the color variable. DEFAULT: NULL
#' @param line_var (Optional) A character string indicating the metadata variable used to set the linetype of the lines. DEFAULT: NULL
#' @param linetypes (Optional) A vector with the line type names associated to the levels of the linetype variable. DEFAULT: NULL
#' @param labels (Optional) A logical value. If TRUE sample names are shown at the end of each curve. DEFAULT: TRUE
#' @param title (Optional) A character string indicating the plot title. DEFAULT: NULL
#' @param ylab (Optional) A character string indicating the title of the y axis. DEFAULT: "Number of different ASVs"
#' @param file_name (Required) A character string indicating the name of the png plot file. DEFAULT: "plot"
#' @param legend.box.col (Optional) A vector indicating the coordinates (position) of the color legend box. DEFAULT: c(0.86, 0.9)
#' @param legend.box.line (Optional) A vector indicating the coordinates (position) of the linetype legend box. DEFAULT: c(0.825, 0.4)
#' @param legend.box1.width (Optional) A numerical value indicating the width of the legend box 1. DEFAULT: 0.3
#' @param legend.box2.width (Optional) A numerical value indicating the width of the legend box 2. DEFAULT: 0.3
#' @param legend.box1.height (Optional) A numerical value indicating the height of the legend box 1. DEFAULT: 0
#' @param legend.box2.height (Optional) A numerical value indicating the height of the legend box 2. DEFAULT: 0
#' @param legend.titlecol (Optional) A character string indicating the hexadecimal color used for the legend title. DEFAULT: "#00008B"
#' @param legend.back (Optional) A character string indicating the hexadecimal color used for the legend background. DEFAULT: "#D3D3D3"
#' @param str2rm (Optional) A vector of character strings indicating the patterns to be removed from the names of the samples. DEFAULT: ""
#' @param return.vegandata (Optional) A logical value. If TRUE a plot_grid object is returned. If FALSE plot is saved in the specified directory. DEFAULT: TRUE
#' @param return.plot (Optional) A logical value. If TRUE the output from the rarecurve function of the Vegan package is returned. DEFAULT: FALSE
#'
#' @return
#' A filtered phyloseq object. It also writes summary files in the specified directory.
#' Or if out.table = TRUE a data frame showing the prevalence of each feature (OTU/ASV/Species) across all samples.
#'
#' @examples
#' data("GlobalPatterns")
#' setwd("~/PATHtoDirectory")
#' rare.curve(GlobalPatterns, file_name = "Rarefaction_curves", col_var = "SampleType")
#'
#' @export
rare.curve <- function(physeq,
                       vegandata = NULL,
                       step = 500,
                       rank = "Species",
                       col_var = NULL,
                       line_var = NULL,
                       colors = NULL,
                       linetypes = NULL,
                       labels = TRUE,
                       title = NULL,
                       y_lab = NULL,
                       legend.text.size = 1,
                       file_name = "plot",
                       legend.box.col = c(0.9, 0.97),
                       legend.box.line = c(0.9, 0.60),
                       legend.box1.width = 0.3,
                       legend.box2.width = 0.3,
                       legend.box1.height = 0,
                       legend.box2.height = 0,
                       legend.titlecol = "#00008B",
                       legend.back = "#D3D3D3",
                       str2rm = "",
                       return.vegandata = FALSE,
                       return.plot = TRUE){

  # To remove warning messages
  oldw <- getOption("warn")
  options(warn = -1)
  ifelse(is.null(y_lab), {y_lab <- "Number of different ASVs"}, "")
  ifelse(is.null(title), {title <- ""}, "")
  raremax <- min(colSums(physeq@otu_table))
  raresample <- sample_names(physeq@otu_table)[colSums(physeq@otu_table) == raremax]

  # Create original plot data
  if (is.null(vegandata)){
    # Agglomerate to specific rank
    if(rank != "Species"){
      agg <- phyloseq::tax_glom(physeq, taxrank = rank, NArm=FALSE)
      if (!is.null(col_var)){
        for (i in seq_along(levels(get_variable(physeq, col_var)))){
          agg@sam_data[,col_var][agg@sam_data[,col_var] == i] <- levels(get_variable(physeq, col_var))[i]
        }
      }
      if (!is.null(col_var)){
        for (i in seq_along(levels(get_variable(physeq, line_var)))){
          agg@sam_data[,line_var][agg@sam_data[,line_var] == i] <- levels(get_variable(physeq, line_var))[i]
        }
      }
      physeq <- agg
    }
    # Create rarefaction curves
    rarefaction.curve <- vegan::rarecurve(t(otu_table(physeq)), step = step, sample = raremax, cex = 1)
    if (isTRUE(return.vegandata)){
      return(rarefaction.curve)
    }
  } else {
    rarefaction.curve <- vegandata
  }
  # Add group color code
  if (!is.null(col_var)){
    var1 <- get_variable(physeq, col_var)
    if(is.null(colors)){
      colors <- Polychrome::createPalette(length(levels(var1)), c("#ff0000", "#00ff00", "#0000ff"))
    }
    cols <- colors[var1]
  }
  # Add group line type code
  if (!is.null(line_var)){
    var2 <- get_variable(physeq, line_var)
    if(is.null(linetypes)){
      linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
      linetypes[1:length(levels(var2))]
      ltys <- linetypes[var2]
    } else {
      ltys <- linetypes[var2]
    }
  }

  # Plot rarefaction curves
  plot.new()
  pdf(NULL)
  dev.control(displaylist = "enable")
  Nmax <- sapply(rarefaction.curve, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(rarefaction.curve, max)
  par(mar = c(5.1, 5.1, 4.1, 2.1))
  plot(c(1, max(Nmax)+max(Nmax)*0.2), c(1, max(Smax)), xlab = "Number of reads", ylab = y_lab, type = "n", main = title, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  abline(v = raremax, lty = "dashed")
  text((raremax+40),-1.6, paste0("Minimum read depth (", raresample, " = ", format(raremax, big.mark=",", scientific = FALSE), " reads)"), col = "gray60",
       adj = c(0, -.1), cex = 1.2)

  for (i in seq_along(rarefaction.curve)) {
    N <- attr(rarefaction.curve[[i]], "Subsample")
    if (!is.null(col_var) & !is.null(line_var)){
      lines(N, rarefaction.curve[[i]], col = cols[i], lty = ltys[i])
    } else if (is.null(col_var) & !is.null(line_var)){
      lines(N, rarefaction.curve[[i]], lty = ltys[i])
    } else if (!is.null(col_var) & is.null(line_var)){
      lines(N, rarefaction.curve[[i]], col = cols[i])
    } else {
      lines(N, rarefaction.curve[[i]])
    }
  }
  # Simplify sample labels and add them
  lab <- colnames(physeq@otu_table)
  for (i in str2rm){
    lab <- sub(i, "", lab)
  }
  if(isTRUE(labels)){
    vegan::ordilabel(cbind(Nmax, Smax), labels = lab, fill = NULL, lty = 0, cex = 1, adj = 0)
  }
  if (!is.null(col_var)){
    n1 <- length(levels(get_variable(physeq, col_var)))
    n1 <- n1 + legend.box1.height
    legend(c( (max(Nmax)*legend.box.col[1]), (max(Nmax)*legend.box.col[1]) + (max(Nmax)*legend.box.col[1])*legend.box1.width),
           c( (max(Smax)*legend.box.col[2]), (max(Smax)*legend.box.col[2] - max(Smax) * as.numeric(paste0("0.0", gsub("\\.","", 6-(n1*0.1)))) * (n1+1))),
           legend = levels(get_variable(physeq, col_var)), lty = 1, lwd = 2, col = colors,
           cex = legend.text.size, bg = legend.back, title = paste0(col_var, ":"), title.col = legend.titlecol, title.adj = 0.1, xpd = TRUE)
  }
  if (!is.null(line_var)){
    n2 <- length(levels(get_variable(physeq, line_var)))
    n2 <- n2 + legend.box2.height
    legend(c((max(Nmax)*legend.box.line[1]), (max(Nmax)*legend.box.line[1])+(max(Nmax)*legend.box.line[1])*legend.box2.width),
           c((max(Smax)*legend.box.line[2]), (max(Smax)*legend.box.line[2] - max(Smax) * as.numeric(paste0("0.0", gsub("\\.","", 6-(n2*0.1)))) * (n2+1))),
           legend = levels(get_variable(physeq, line_var)), lty = linetypes, lwd = 2,
           cex = legend.text.size, bg = legend.back, title = paste0(line_var, ":"), title.col = legend.titlecol, title.adj = 0.1, xpd = TRUE)
  }
  plot <- recordPlot()
  invisible(dev.off())
  if (isTRUE(return.plot)){
    plot.new()
    return(plot)
  } else {
    # Create directory
    ifelse(!dir.exists("Rarefaction_curves"), dir.create("Rarefaction_curves"), "")
    png(paste0("Rarefaction_curves/", file_name ,".png"), height = 7, width = 18, units = "in", res = 300)
    print({plot})
    dev.off()
  }
  options(warn = oldw)
}
###############################################################################################################################################################################
#' @title # Add colors to ggplot wrap strips
#'
#' @description
#' This function adds random or specified colors to the strips of panels created by facet_wrap.
#' Colors can be added to the Top and/or Bottom strips. They also can be assigned to the levels of the variables selected.
#'
#' @param df (Required) Data frame used to create the existing plot.
#' @param plot (Required) Existing plot created using the facet_wrap option.
#' @param Top.var (Required) A character string indicating the data frame variable used in facet_wrap. DEFAULT: NULL
#' @param Bot.var (Optional) A character string indicating a second data frame variable used in facet_wrap. DEFAULT: NULL
#' @param Top.strip.fill (Optional) A vector with the hexadecimal colors associated to the levels of the Top variable.
#' OR a character string ("by.Top.var" or "by.Bot.var" or "blank") indicating if colors are assigned to the levels of the Top or Bottom variables. DEFAULT: "by.Top.var"
#' @param Bot.strip.fill (Optional) A vector with the hexadecimal colors associated to the levels of the Bottom variable.
#' OR a character string ("by.Top.var", "by.Bot.var" or "blank") indicating if colors are assigned to the levels of the Top or Bottom variables. DEFAULT: "by.Bot.var"
#' @param col.palette (Optional) A character vector containing the hexadecimal representations of one or more colors (seedcolors). DEFAULT: c("#ff0000", "#00ff00", "#0000ff").
#' @param display (Optional) A logical value. If TRUE the function returns a ggplot object, otherwise, it returns a gtable() object. DEFAULT: TRUE
#'
#' @return
#' A ggplot object or if display = FALSE a gtable() object.
#'
#' @examples
#' data("GlobalPatterns")
#' colnames(tax_table(GlobalPatterns)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' filt.data <- FilterAndSummarise(GlobalPatterns)
#' mdf = psmelt(filt.data)
#' p = ggplot(mdf, aes(x=SampleType, y=Abundance, fill=Genus))
#' p = p + geom_bar(color="black", stat="identity", position="stack") + facet_wrap(~SampleType) + theme(legend.position = "none")
#' print(p)
#' p2 <- cowplot::plot_grid(wrap.strip.colors(mdf, p, Top.var = "SampleType", Top.strip.fill = "by.Top.var"))
#' print(p2)
#'
#' @export
wrap.strip.colors <- function(df,
                              plot,
                              Top.var = NULL,
                              Bot.var = NULL,
                              Top.strip.fill = "by.Top.var",
                              Bot.strip.fill = "by.Bot.var",
                              col.palette = c("#ff0000", "#00ff00", "#0000ff"),
                              display = TRUE){

  Table.of.Plot <- ggplot_gtable(ggplot_build(plot)) # Create ggplot table

  stript <- which(grepl('strip-t', Table.of.Plot$layout$name)) # Find where are the strips within the plot layout

  if (!is.null(Top.var)){
    # Add filling to the Top strip
    if (Top.strip.fill == "blank"){
      colT <- rep("#ffffff", length(levels(df[,Top.var])))
    } else if (Top.strip.fill == "by.Top.var"){
      colT <- unname(Polychrome::createPalette(length(levels(df[,Top.var])), col.palette))
    } else if (Bot.strip.fill == "by.Bot.var"){
      colT <- unname(Polychrome::createPalette(length(levels(df[,Bot.var])), col.palette))
    } else {
      colT <- Top.strip.fill
    }

    Top.fill <- list()
    if (Top.strip.fill == "by.Bot.var"){
      df$Col <- ""
      for (i in c(1:length(levels(df[,Bot.var])))){
        level <- levels(df[,Bot.var])[i]
        df <- df %>% mutate(Col = ifelse(df[,Bot.var] == level, colT[i], Col))
      }
      df$Col <- factor(df$Col, levels = colT)
      for (i in levels(df[,Top.var])){
        x <- df %>% filter(df[,Top.var] == i)
        x <- with(x, x[order(x[,Bot.var]),]) # Really important to maintain the levels order
        x[,Bot.var] <- factor(x[,Bot.var], levels = unique(x[,Bot.var]))
        x$Col <- factor(x$Col, levels = unique(x$Col))
        Top.fill[[i]] <- levels(x$Col)
      }
    } else {
      n <- 1
      for (i in levels(df[,Top.var])){
        x <- df %>% dplyr::filter(df[,Top.var] == i)
        if (!is.null(Bot.var)){
          x[,Bot.var] <- factor(x[,Bot.var], levels = unique(x[,Bot.var]))
          Top.fill[[n]] <- rep(colT[n], length(levels(x[,Bot.var])))
        } else {
          Top.fill[[n]] <- colT[n]
        }
        n <- n + 1
      }
    }
    Top.strip <- as.vector(unlist(Top.fill))
  }

  if (!is.null(Bot.var)){
    # Add filling to the Bottom strip
    if (Bot.strip.fill == "blank"){
      colB <- rep("#ffffff", length(levels(df[,Bot.var])))
    } else if (Bot.strip.fill == "by.Bot.var" & Top.strip.fill == "by.Bot.var"){
      colB <- colT
    } else if (Bot.strip.fill == "by.Bot.var"){
      colB <- unname(Polychrome::createPalette(length(levels(df[,Bot.var])), col.palette))
    } else if (Bot.strip.fill == "by.Top.var" & Top.strip.fill == "by.Top.var"){
      colB <- colT
    } else if (Bot.strip.fill == "by.Top.var"){
      colB <- unname(Polychrome::createPalette(length(levels(df[,Top.var])), col.palette))
    } else {
      colB <- Bot.strip.fill
    }

    Bottom.fill <- list()
    if (Bot.strip.fill == "by.Top.var"){
      n <- 1
      for (i in levels(df[,Top.var])){
        x <- df %>% filter(df[,Top.var] == i)
        x[,Bot.var] <- factor(x[,Bot.var], levels = unique(x[,Bot.var]))
        Bottom.fill[[n]] <- rep(colB[n], length(levels(x[,Bot.var])))
        n <- n + 1
      }
    } else {
      df$Col <- ""
      for (i in c(1:length(levels(df[,Bot.var])))){
        level <- levels(df[,Bot.var])[i]
        df <- df %>% mutate(Col = ifelse(df[,Bot.var] == level, colB[i], Col))
      }
      for (i in levels(df[,Top.var])){
        x <- df %>% filter(df[,Top.var] == i)
        x <- with(x, x[order(x[,Bot.var]),]) # Really important to maintain the levels order
        x[,Bot.var] <- factor(x[,Bot.var], levels = unique(x[,Bot.var]))
        x$Col <- factor(x$Col, levels = unique(x$Col))
        if (nlevels(x$Col) != nlevels(x[,Bot.var])){
          ncols <- nlevels(x[,Bot.var]) - nlevels(x$Col)
          Bottom.fill[[i]] <- c(levels(x$Col), rep("#ffffff", ncols))
        } else {
          Bottom.fill[[i]] <- levels(x$Col)
        }
      }
    }
    Bot.strip <- as.vector(unlist(Bottom.fill))
  }

  k <- 1
  for (i in stript) {
    if (!is.null(Top.var)){
      j1 <- which(grepl('rect', Table.of.Plot$grobs[[i]]$grobs[[1]]$childrenOrder))
      j3 <- which(grepl('title', Table.of.Plot$grobs[[i]]$grobs[[1]]$childrenOrder))
      Table.of.Plot$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- Top.strip[k]
      Table.of.Plot$grobs[[i]]$grobs[[1]]$children[[j3]]$children[[1]]$gp$font <- as.integer(2)
    }
    if (!is.null(Bot.var)){
      j2 <- which(grepl('rect', Table.of.Plot$grobs[[i]]$grobs[[2]]$childrenOrder))
      j3 <- which(grepl('title', Table.of.Plot$grobs[[i]]$grobs[[2]]$childrenOrder))
      Table.of.Plot$grobs[[i]]$grobs[[2]]$children[[j2]]$gp$fill <- Bot.strip[k]
      Table.of.Plot$grobs[[i]]$grobs[[2]]$children[[j3]]$children[[1]]$gp$font <- as.integer(1)
    }
    k <- k + 1
  }

  ifelse(isTRUE(display), {grid.object<-ggplotify::as.ggplot(Table.of.Plot)}, {grid.object<-Table.of.Plot})
  return(grid.object)
}
###############################################################################################################################################################################
#' @title Taxonomic profile, count data and sample clustering tree plots.
#'
#' @description
#' This function can create 3 different plots: a taxonomic profile (barplot or poinplot), a stacked bar plot with count data per domain
#' and a clustering tree diagram of the samples. The 3 plots can be visualised individually or together.
#'
#' @param physeq (Required) A phyloseq object.
#' @param abundance.plot (Optional) A character string indicating how to display the taxonomic profile ("barplot", "poinplot" or "blank"). DEFAULT: "pointplot"
#' @param count.plot (Optional) A logical value. If TRUE a stacked bar plot with count data per domain is displayed. DEFAULT: FALSE
#' @param data.count.plot (Optional, only required if count.plot = TRUE) A phyloseq object with the data used to create the count.plot.
#' Since this plot is made to provide library depth information unfiltered data should be provided. DEFAULT: NULL
#' @param count.plot.rel (Optional)
#' @param domains (Optional)
#' @param dendogram (Optional) A logical value. If TRUE a clustering tree diagram of the samples is displayed. DEFAULT: FALSE
#' @param rank (Optional) A character string indicating the taxonomic rank shown in the abundance.plot. DEFAULT: "Phylum"
#' @param ntaxa (Optional) A numerical value indicating the number of taxa shown in the abundance.plot. DEFAULT: 15
#' @param taxa.txt.size (Optional) A numerical value indicating the text size of the taxonomic names.
#' Only required if an abundance.plot is created. DEFAULT: 12
#' @param sample.txt.size (Optional) A numerical value indicating the text size of the sample names. DEFAULT: 8
#' @param show.unidentified (Optional) A logical value. If TRUE the unclassified/n.a taxa are renamed using the information their higher taxonomic ranks. DEFAULT: TRUE
#' @param name.format (Optional) A character string indicating the name format of the renamed unclassified taxa ("short" or "long", for more info: \link[metaGVis]{unclassified_ranks}). DEFAULT: "short"
#' @param ASV.count (Optional) A logical value. If TRUE the the number of species agglomerated under a higher taxonomic rank is appended at the end of the taxonomic name. DEFAULT: TRUE
#' @param col.palette (Optional) A character vector containing the hexadecimal representations of one or more colors (seedcolors). DEFAULT: c("#ff0000", "#00ff00", "#0000ff")
#' @param title (Optional) Main plot title. DEFAULT: ""
#' @param return.plot (Optional) A logical value. If TRUE a plot_grid object is returned. If FALSE plot is saved in the specified directory. DEFAULT: TRUE
#' @param folder (Optional) A character string indicating the PATH to the directory where the plot is saved. DEFAULT: NULL
#' @param file.name (Optional) File name. DEFAULT: "plot.png"
#' @param str2rm (Optional) A vector of character strings defining the patterns to be removed from the names of the samples. DEFAULT: NULL
#' @param Top.var (Optional) A character string indicating the data frame variable used in facet_wrap.
#' For more info: \link[metaGVis]{wrap.strip.colors}). DEFAULT: NULL
#' @param Top.strip.fill (Optional) A vector with the hexadecimal colors associated to the levels of the Top variable.
#' OR a character string ("by.Top.var" or "by.Bot.var" or "blank") indicating if colors are assigned to the levels of the Top or Bottom variables.
#' For more info: \link[metaGVis]{wrap.strip.colors}). DEFAULT: "by.Top.var"
#' @param Bot.var (Optional) A character string indicating a second data frame variable used in facet_wrap. DEFAULT: NULL
#' @param Bot.strip.fill (Optional) A vector with the hexadecimal colors associated to the levels of the Bottom variable.
#' OR a character string ("by.Top.var", "by.Bot.var" or "blank") indicating if colors are assigned to the levels of the Top or Bottom variables. DEFAULT: "by.Bot.var"
#' @param nboot (Optional) A numeric value indicating the number of bootstrap replications. DEFAULT: 500
#' @param aglom.clust (Optional) A logical value. If TRUE the clustering tree diagram of the samples is performed using the agglomerated data at the specified rank. DEFAULT: FALSE
#' @param transf (Optional) A character string indicating the transformation performed on the data before clustering
#' (Hellinger, logChord_0plus1, logChord_1p, CLR_0plus1, CLR_min or CLR_runif). DEFAULT: "CLR_runif"
#' @param distance (Optional) A character string indicating the clustering distance (Bray, wUniFrac, UniFrac, Euclidean). DEFAULT: "Euclidean"
#' @param clust.method (Optional) A character string indicating the hierarchical clustering algorithm used (single, complete, average, ward.D). DEFAULT: "ward.D"
#' @param n.clusters (Optional) A numeric value indicating the number of clusters to color code. DEFAULT: 2
#' @param clust.data.rel.abund (Optional) A logical value. If TRUE the abundance data is transformed into relative abundance before clustering (only applicable for Bray and wUniFrac). DEFAULT: FALSE
#' @param nudge.bs.label (Optional) A numeric value (see y.axis range) indicating how much the bootstrap values are moved vertically. DEFAULT: NULL
#' @param bs.label.size (Optional) A numeric value indicating the size of the bootstrap values. DEFAULT: 2.5
#' @param icons.dendogram (Optional) A logical value. If TRUE specified icons representing samples are located at the end of the dendogram branches. DEFAULT: FALSE
#' @param icons.variable (Optional) A character string indicating the variable use to assign icons. DEFAULT: NULL
#' @param icons.folder (Optional) A character string indicating the PATH to the directory where the icon files are located.
#' PNG files must have the same names as the icons.variable levels. DEFAULT: NULL
#' @param icons.col (Optional) A vector with the hexadecimal colors associated to the levels of the icons.variable. This option will make the icons color-coded. DEFAULT: "blank"
#' @param tile.col (Optional) A vector with the hexadecimal colors associated to the levels of the Bottom variable. The selected variable is used to color code the tile found at the end of the dendogram branches.
#' OR a character string ("by.Top.var", "by.Bot.var" or "blank") indicating if colors are assigned to the levels of the Top or Bottom variables. DEFAULT: "by.Bot.var".
#' The selected variable (Top or Bottom) is used to color code a tile found at the end of the dendogram branches.
#'
#' @return
#' The plot is saved as a png file in the specified directory. If return.plot = TRUE, it returns a plot_grid() object.
#'
#' @examples
#' data("GlobalPatterns")
#' # Since all metaGVis functions use "Domain" rank instead of Kingdom, let's change it:
#' colnames(tax_table(GlobalPatterns)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' # To make sure the first metadata column is named "SampleID":
#' colnames(sample_data(GlobalPatterns))[1] <- "SampleID"
#' filt.data <- FilterAndSummarise(GlobalPatterns)
#'
#' stackedbar.plot <- Taxonomic.analyis(filt.data, abundance.plot = "barplot", rank = "Order",
#' title = "Relative abundance of the 15 most predominant bacterial orders")
#' print(stackedbar.plot)
#'
#' stackedbar.plot.with.facets.and.1.variable <- Taxonomic.analyis(filt.data, abundance.plot = "barplot", rank = "Order",
#' Top.var = "SampleType", title = "Relative abundance of the 15 most predominant bacterial orders")
#' print(stackedbar.plot.with.facets.and.1.variable)
#'
#' stackedbar.plot.with.facets.and.2.variables <- Taxonomic.analyis(filt.data, abundance.plot = "barplot", rank = "Order",
#' Top.var = "SampleType", Top.strip.fill = "by.Bot.var", Bot.var = "Primer", Bot.strip.fill = "by.Bot.var",
#' title = "Relative abundance of the 15 most predominant bacterial orders")
#' print(stackedbar.plot.with.facets.and.2.variables)
#'
#' point.plot <- Taxonomic.analyis(filt.data, abundance.plot = "pointplot", rank = "Order",
#' title = "Relative abundance of the 15 most predominant bacterial orders")
#' print(point.plot)
#'
#' count.plot <- Taxonomic.analyis(data.count.plot = GlobalPatterns, count.plot = TRUE, abundance.plot = "blank",
#' title = "Relative abundance of the 15 most predominant bacterial orders")
#' print(count.plot)
#'
#' combined.pointplot <- Taxonomic.analyis(filt.data, data.count.plot = GlobalPatterns, abundance.plot = "pointplot", rank = "Order", count.plot = TRUE, dendogram = TRUE, nboot = 1,
#' title = "Relative abundance of the 15 most predominant bacterial orders")
#' print(combined.plot)
#'
#' combined.barplot.with.1.variable <- Taxonomic.analyis(filt.data, data.count.plot = GlobalPatterns, abundance.plot = "barplot", rank = "Phylum",
#' count.plot = TRUE, dendogram = T, nboot = 1, Top.var = "SampleType", tile.col = "by.Top.var",
#' title = "Relative abundance of the 15 most predominant bacterial orders")
#' print(combined.barplot.with.variables)
#'
#' combined.pointplot.with.1.variable <- Taxonomic.analyis(filt.data, data.count.plot = GlobalPatterns, abundance.plot = "pointplot", rank = "Phylum",
#' count.plot = TRUE, dendogram = T, nboot = 1, Top.var = "SampleType", tile.col = "by.Top.var",
#' title = "Relative abundance of the 15 most predominant bacterial orders")
#'
#' combined.pointplot.with.2.variables <- Taxonomic.analyis(filt.data, data.count.plot = GlobalPatterns, abundance.plot = "pointplot", rank = "Phylum",
#' count.plot = TRUE, dendogram = F, nboot = 1, Top.var = "SampleType", Bot.var = "Primer", Top.strip.fill = "blank", Bot.strip.fill = "by.Bot.var",
#' title = "Relative abundance of the 15 most predominant bacterial orders")
#'
#' @export
Taxonomic.analyis <- function(physeq,
                              merge.by = "None",
                              abundance.plot = "pointplot",
                              count.plot = FALSE,
                              data.count.plot = NULL,
                              count.plot.rel = FALSE,
                              domains = c("Viruses", "Archaea", "Bacteria", "Eukaryota", "Unclassified"),
                              dendogram = FALSE,
                              rank = "Phylum",
                              ntaxa = 15,
                              taxa.txt.size = 12,
                              sample.txt.size = 8,
                              strip.text.size = 8,
                              show.unidentified = TRUE,
                              name.format = "short",
                              ASV.count = TRUE,
                              col.palette = c("#ff0000", "#00ff00", "#0000ff"),
                              title = "Bacterial phyla",
                              return.plot = TRUE,
                              folder = NULL,
                              file.name = "plot.png",
                              str2rm = NULL,
                              Top.var = NULL,
                              Top.strip.fill = "by.Top.var",
                              Bot.var = NULL,
                              Bot.strip.fill = "by.Bot.var",
                              nboot = 500,
                              aglom.clust = FALSE,
                              transf = "CLR_runif",
                              distance = "Euclidean",
                              clust.method = "ward.D",
                              n.clusters = 2,
                              clust.data.rel.abund = FALSE,
                              nudge.bs.label = NULL,
                              bs.label.size = 2.5,
                              icons.dendogram = FALSE,
                              icons.variable = NULL,
                              icons.folder = NULL,
                              icons.col = "blank",
                              tile.col = "blank"){

  oldw <- getOption("warn")
  options(warn = -1)

  if (isTRUE(dendogram) & isTRUE(count.plot) & abundance.plot == "blank"){
    cat("This configuration is not possible. \nPlease specify an abundance plot or set the dendogram or count.plot arguments to FALSE.")
    on.exit(options(show.error.messages = FALSE))
    stop()
  }
  if (isTRUE(dendogram) & merge.by != "None"){
    cat("This configuration is not possible. \nThe dendogram and count plots cannot be created using a merged phyloseq object.")
    on.exit(options(show.error.messages = FALSE))
    stop()
  } else if (isTRUE(count.plot) & merge.by != "None"){
    cat("This configuration is not possible. \nThe dendogram and count plots cannot be created using a merged phyloseq object.")
    on.exit(options(show.error.messages = FALSE))
  }
  ################################################################################
  # CODE and TEXT adapted from https://rdrr.io/bioc/phyloseq/src/R/ordination-methods.R
  # Define an internal function for accessing and orienting the OTU table
  # in a fashion suitable for vegan functions
  # @keywords internal
  veganifyOTU <- function(physeq){
    if(taxa_are_rows(physeq@otu_table)){physeq <- t(physeq@otu_table)}
    return(as(otu_table(physeq), "matrix"))
  }
  ################################################################################

  if (!missing(physeq)){

    if (isTRUE(show.unidentified)){
      agg <- unclassified_ranks(physeq, rank, format = name.format, agglomerate = TRUE, count = ASV.count, output = "physeq")
      n.of.taxa <- nlevels(as.factor(tax_table(agg)[, rank]))
    } else {
      n.of.taxa <- nlevels(as.factor(tax_table(physeq)[, rank]))
      agg <- phyloseq::tax_glom(physeq, taxrank = rank, NArm=FALSE, bad_empty=c()) # agglomerate data
    }
    if(merge.by != "None"){
      agg <- merge_samples(agg, merge.by, fun = sum)
      # Re-create factor levels
      agg@sam_data[,"SampleID"] <- agg@sam_data[,merge.by] <- factor(rownames(agg@sam_data),
                                                                     levels = levels(get_variable(physeq, merge.by)))
      Top.var <- merge.by

      for (i in seq_along(levels(get_variable(physeq, Bot.var)))){
        agg@sam_data[,Bot.var][agg@sam_data[,Bot.var] == i] <- levels(get_variable(physeq, Bot.var))[i]
      }
      agg@sam_data[,Bot.var] <- factor(get_variable(agg, Bot.var),
                                       levels = levels(get_variable(physeq, Bot.var)))
    }
    data = transform_sample_counts(agg, function(x) x/sum(x)) # transform to relative abundance

    if(n.of.taxa > ntaxa){
      # Top taxa names
      Top.seq = names(sort(taxa_sums(data), TRUE)[1:(ntaxa-1)])
      Top.pruned <- prune_taxa(Top.seq, data)
      # Rename NA to Unclassified
      for (i in colnames(Top.pruned@tax_table)){
        Top.pruned@tax_table[,i][is.na(Top.pruned@tax_table[,i])] <- "Unclassified"
      }
      # Make all unclassified names unique
      Top.pruned@tax_table[,rank] <- make.unique(unname(Top.pruned@tax_table[,rank]))
      # Remove trees to merge phyloseq objects
      Phylo4tree <- Top.pruned
      data@phy_tree = NULL
      Top.pruned@phy_tree = NULL
      # Bottom taxa names
      ASV.seq <- rownames(data@tax_table)
      Bottom.seq = subset(ASV.seq, !(ASV.seq %in% Top.seq))
      Bottom.pruned <- prune_taxa(Bottom.seq, data)
      Bottom.pruned@tax_table[,rank] <- "Other"
      # Reassemble
      Final.object <- merge_phyloseq(Top.pruned, Bottom.pruned)
      # Melt phyloseq object
      Melted.final.object <- psmelt(Final.object)
    } else {
      Top.seq = names(sort(taxa_sums(data), TRUE))
      # Re-subsample
      Top.pruned <- prune_taxa(Top.seq, data)
      # Rename NA to Unclassified
      for (i in colnames(Top.pruned@tax_table)){
        Top.pruned@tax_table[,i][is.na(Top.pruned@tax_table[,i])] <- "Unclassified"
      }
      # Make all unclassified names unique
      Top.pruned@tax_table[,rank] <- make.unique(unname(Top.pruned@tax_table[,rank]))
      # Melt phyloseq object
      Melted.final.object <- psmelt(Top.pruned)
    }

    # Add "Other" at the beginning
    if(any(grepl("Other", Melted.final.object[,rank]))){
      taxnames <- c("Other", Top.pruned@tax_table[,rank][names(sort(taxa_sums(Top.pruned)))])
    } else {
      taxnames <- Top.pruned@tax_table[,rank][names(sort(taxa_sums(Top.pruned)))]
    }

    # Convert variables to factors
    Melted.final.object[,rank] <- factor(Melted.final.object[,rank], levels = taxnames)

    # Make sure a "SampleID" variable exists (to change sample names without errors)
    if (!"SampleID" %in% colnames(Melted.final.object)){
      Melted.final.object[,"SampleID"] <- Melted.final.object[,"Sample"]
    }

    # Make Sample names shorter
    if (!is.null(str2rm)){
      for (i in str2rm){
        Melted.final.object[,"SampleID"] <- sub(i, "", Melted.final.object[,"SampleID"])
      }
    }

    # To omit relative abundance dot when bacteria is absent in sample
    Melted.final.object <- Melted.final.object[Melted.final.object$Abundance > 0, ]

    # Number of colors
    if(any(taxnames == "Other")){
      mycolors <- Polychrome::createPalette((length(taxnames)-1), col.palette)
      names(mycolors) <- NULL
      mycolors <- rev(mycolors)
      mycolors <- c("#DDDDDD", mycolors)
    } else {
      mycolors <- Polychrome::createPalette(length(taxnames), col.palette)
      names(mycolors) <- NULL
      mycolors <- rev(mycolors)
    }
  }

  # To create dendogram plot
  if (isTRUE(dendogram)){

    if (!isTRUE(aglom.clust)){
      data.dendo = physeq
    } else {
      data.dendo = data
    }

    if (clust.data.rel.abund == TRUE){
      data.dendo = transform_sample_counts(data.dendo, function(x) x/sum(x)) # transform to relative abundance
    }

    # Create dendogram plot
    if (distance == "Bray" | distance == "bray"){
      dist <- sqrt(vegan::vegdist(veganifyOTU(data.dendo), method = "bray"))
    } else if (distance == "unifrac"| distance == "UniFrac"){
      dist <- phyloseq::distance(data.dendo, method = "unifrac")
    } else if (distance == "wunifrac"| distance == "wUniFrac"){
      dist <- phyloseq::distance(data.dendo, method = "wunifrac")
    } else if (distance == "Euclidean" & transf == "logChord_0plus1"){
      input_data <- logchord_0_plus_1(data.dendo)
      dist <- vegan::vegdist(veganifyOTU(input_data), method = "euclidean")
    } else if (distance == "Euclidean" & transf == "logChord_1p"){
      input_data <- logchord_1p(data.dendo)
      dist <- vegan::vegdist(veganifyOTU(input_data), method = "euclidean")
    } else if (distance == "Euclidean" & transf == "Hellinger"){
      input_data <- hellinger(data.dendo)
      dist <- vegan::vegdist(veganifyOTU(input_data), method = "euclidean")
    } else if (distance == "Euclidean" & transf == "CLR_min"){
      input_data <- clr_min_percent(data.dendo)
      dist <- vegan::vegdist(veganifyOTU(input_data), method = "euclidean")
    } else if (distance == "Euclidean" & transf == "CLR_0plus1"){
      input_data <- clr_0_plus_1(data.dendo)
      dist <- vegan::vegdist(veganifyOTU(input_data), method = "euclidean")
    } else if (distance == "Euclidean" & transf == "CLR_runif"){
      input_data <- clr_runif(data.dendo)
      dist <- vegan::vegdist(veganifyOTU(input_data), method = "euclidean")
    }

    hc <- hclust(dist, method = clust.method) # heirarchal clustering
    dendr <- ggdendro::dendro_data(hc, type="rectangle") # convert to ggplot
    dendr$labels$label <- factor(dendr$labels$label,
                                 levels = levels(as.factor(dendr$labels$label)))
    clust <- dendextend::cutree(hc, k = n.clusters, order_clusters_as_data = FALSE) # find k clusters
    clust.df <- data.frame(label = names(clust), cluster = factor(clust))
    # dendr[["labels"]] has the labels, merge with clust.df based on label column
    dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label", sort = FALSE)

    # Split dendrogram into upper grey section and lower coloured section
    height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
    cut.height <- mean(c(height[n.clusters], height[n.clusters - 1]))
    dendr$segments$line <- ifelse(dendr$segments$y == dendr$segments$yend & dendr$segments$y > cut.height, 1, 2)
    dendr$segments$line <- ifelse(dendr$segments$yend  > cut.height, 1, dendr$segments$line)

    # Number the clusters
    dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
    change <- which(dendr$segments$cluster == 1)
    for (i in 1:n.clusters){
      dendr$segments$cluster[change[i]] = i + 1
    }
    dendr$segments$cluster <-  ifelse(dendr$segments$line == 1, 1,
                                      ifelse(dendr$segments$cluster == 0, NA,
                                             dendr$segments$cluster))
    dendr$segments$cluster <- zoo::na.locf(dendr$segments$cluster)

    # Establish horizontal line height cut (3.5/4 of last shortest colored segment)
    col.seg <- dendr$segments %>% filter(cluster != 1) %>% group_by(cluster) %>% filter(yend == max(yend)) %>% arrange(-y)
    last.col.seg <- col.seg[1:n.clusters,] # last colored segments of each cluster
    seg.col.y <- min(last.col.seg$y) # smallest "y value"
    seg.col.yend <- max(last.col.seg$yend[last.col.seg$y == seg.col.y ]) # largest "yend value" after filtering
    h.line = seg.col.yend + ((seg.col.y - seg.col.yend)*3.5)/4
    y.tile = paste0(-max(height)*0.33) # y value where the plot tile starts
    height.tile = paste0(max(height)*0.11) # height of the plot tile
    y.image = max(height)*0.15 # y value where the plot images start (y - y.image)

    # Shorten labels
    dendr$labels$shortlabel <- dendr$labels$label
    if (!is.null(str2rm)){
      for (i in str2rm){
        dendr$labels$shortlabel <- sub(i, "", dendr$labels$shortlabel)
      }
    }
    # Consistent numbering between segment$cluster and label$cluster
    clust.df$label <- factor(clust.df$label, levels = levels(dendr$labels$label))
    clust.df <- dplyr::arrange(clust.df, label)
    clust.df$cluster <- factor((clust.df$cluster), levels = unique(clust.df$cluster),
                               labels = (1:n.clusters) + 1)
    dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")

    # Positions for cluster labels
    n.rle <- rle(dendr$segments$cluster)
    N <- cumsum(n.rle$lengths)
    N <- N[seq(1, length(N), 2)] + 1
    N.df <- dendr$segments[N, ]
    N.df$cluster <- N.df$cluster - 1

    if (isTRUE(icons.dendogram) & !is.null(icons.variable)){
      # Add metadata variable to dendr
      var <- data.frame(data.dendo@sam_data[,icons.variable])
      var$label <- rownames(var)
      dendr$labels <- merge(dendr$labels, var, by = "label")
      # Add icons as labels and save icons for legend
      dendr$labels$IMAGE <- NA
      icons <- c()
      for (i in levels(dendr$labels[,icons.variable])){
        dendr$labels[dendr$labels[,icons.variable] == i, ]$IMAGE <- paste0(icons.folder, "/", i, ".png")
        icons[i] <- paste0("<p><img src='", icons.folder, "/", i, ".png' width='16'> ",i,"</p>")
      }
    }

    # Re-order "Melted.final.object" to match "dendr" labels
    sample_pos_table <- with(dendr$labels, data.frame(x = x, Sample = as.character(label)))
    Final_data <- Melted.final.object %>% dplyr::left_join(sample_pos_table, by = "Sample")

    # Limits for horizontal panel
    axis_limits <- with(sample_pos_table, c(min(x - 1), max(x + 1)))

    ## Add information for nodes (see: https://rpubs.com/Roeland-KINDT/702542)
    # Calculate bootstrap p-values
    if (distance == "Bray" & clust.method == "ward.D"){
      # Ward.D requires data that can be run with the euclidean distance, thus Bray distance values are square rooted
      #pvclust requires a function, not the distance
      pvclust.dist = function(data){sqrt(vegan::vegdist(veganifyOTU(data), method = "bray"))}
    } else {
      pvclust.dist = "euc"
    }
    boot <- pvclust::pvclust(input_data@otu_table,
                             method.hclust = clust.method,
                             method.dist = pvclust.dist,
                             nboot = nboot,
                             parallel = FALSE,
                             iseed = 711)
    # Make sure that all heights are unique
    if (anyDuplicated(boot$hclust$height) != 0) {
      .Machine$double.eps
      small.add <- seq_along(boot$hclust$height) * 1e-15
      boot$hclust$height <- boot$hclust$height + small.add
    }

    # Obtain data sets for nodes and edges
    data.graph <- as.list(tidygraph::as_tbl_graph(boot$hclust))
    data.nodes <- data.frame(data.graph$nodes)
    data.edges <- data.frame(data.graph$edges)

    # Add information on nodes about the merging process
    merge.data <- data.frame(cbind(boot$hclust$merge, boot$hclust$height))
    names(merge.data) <- c("m1", "m2", "height")
    merge.data$ID <- c(1:nrow(merge.data))
    data.nodes2 <- dplyr::left_join(data.nodes, merge.data, by="height")
    edges <- data.frame(boot$edges)
    edges$ID <- c(1:nrow(edges))
    data.nodes3 <- dplyr::left_join(data.nodes2, edges, by="ID")

    # Add information to edges
    data.nodes3$SEQ <- c(1:nrow(data.nodes3))
    data.edges2 <- dplyr::left_join(data.edges,
                                    data.nodes3[, c("SEQ", "height", "au")],
                                    by = c("from" = "SEQ"))

    # Generate the dendrogram with ggraph
    data.graph <- tidygraph::tbl_graph(nodes = data.nodes3, edges = data.edges2, directed = TRUE)
    # Plot
    Dendo <- ggraph::ggraph(data.graph, layout = "dendrogram", circular = FALSE, height = height) +
      geom_segment(data = ggdendro::segment(dendr), aes(x = x, y = y, xend = xend, yend = yend,
                                                        size = factor(line),
                                                        colour = factor(cluster)),
                   lineend = "square", show.legend = FALSE) +
      geom_hline(yintercept = h.line, linetype = "dashed") +
      scale_colour_manual(values = c("grey60", rainbow(n.clusters))) +
      scale_size_manual(values = c(.1, 1)) +
      geom_text(data = N.df,
                aes(x = x, y = (y + min(y)*0.15), label = factor(cluster), colour = factor(cluster + 1)),
                hjust = 1.5,
                show.legend = FALSE) +
      labs(title = title,
           y = paste0(distance, " distance"),
           x = "Samples") +
      scale_x_continuous(breaks = sample_pos_table$x,
                         limits = axis_limits,
                         expand = c(0, 0),
                         labels = sample_pos_table$Sample) +
      scale_y_continuous(breaks = signif(seq.int(0, max(height), length.out = 6),
                                         digits = ifelse(nchar(round(max(height))) == 1,1,
                                                         ifelse(nchar(round(max(height))) == 2,1,
                                                                ifelse(nchar(round(max(height))) == 3,2,
                                                                       ifelse(nchar(round(max(height))) == 4,3,4)))))) +
      theme(legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.background=element_rect(fill="white"),
            plot.title = element_text(face = "bold",
                                      margin = margin(t = 10, r = 0, b = 10, l = 0)),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 10, colour = 'black'),
            axis.line.y = element_line(),
            axis.ticks.y = element_line(),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(fill = NA, colour = "grey50"),
            plot.margin = unit(c(0, 0, 0, 0), "cm"))

    if (tile.col[1] != "blank"){
      # Create colors for tile
      if (tile.col[1]  == "by.Bot.var"){
        col.tile <- unname(Polychrome::createPalette(length(levels(Final_data[,Bot.var])), col.palette))
        var.tile <- Bot.var
      } else if (tile.col[1]  == "by.Top.var"){
        col.tile <- unname(Polychrome::createPalette(length(levels(Final_data[,Top.var])), col.palette))
        var.tile <- Top.var
      } else {
        col.tile <- tile.col
        ifelse(is.null(Bot.var), {var.tile<-Top.var}, {var.tile<-Bot.var})
      }
      Dendo <- Dendo + geom_tile(data = Final_data,
                                 aes_string(x = "x",
                                            y = y.tile,
                                            height = height.tile,
                                            fill = var.tile),
                                 stat = "identity") +
        scale_fill_manual(values = col.tile) +
        guides(fill = guide_legend(ncol = ifelse(length(levels(Final_data[,var.tile]))>5, 2, 1)))
    } else {
      var.tile <- NULL
    }
    if (isTRUE(icons.dendogram) & !is.null(icons.variable)){
      # Detect number of icons
      icons.lev <- levels(Final_data[,icons.variable])
      num.icon.var <- length(icons.lev)
      # Assign icons to legend keys ## Reference: https://stackoverflow.com/questions/61327081/how-to-use-an-image-as-a-legend-key-glyph
      options(icons.lev = icons.lev)
      options(icons.folder = icons.folder)
      n.icon = 1 # establish global counter
      draw_key_image_col <- function(data, params, size){
        icons.folder <- getOption("icons.folder")
        icons.lev <- getOption("icons.lev")
        icons.col <- getOption("icons.col")
        grobs <- lapply(seq_along(data$colour),
                        function(i) {
                          img <- magick::image_read(paste0(icons.folder, "/",  icons.lev[n.icon], ".png"))
                          img <- ggimage:::color_image(img, data$colour[i], data$alpha[i])
                          n.icon <<- n.icon + 1
                          grid::rasterGrob(0.5, 0.5, image = img, width = 1, height = 1)
                        })
        class(grobs) <- "gList"
        ggplot2:::ggname("image_key", grid::gTree(children = grobs))
      }
      draw_key_image <- function(data, params, size){
        icons.folder <- getOption("icons.folder")
        icons.lev <- getOption("icons.lev")
        icons.col <- getOption("icons.col")
        grobs <- lapply(seq_along(data$colour),
                        function(i) {
                          img <- magick::image_read(paste0(icons.folder, "/",  icons.lev[n.icon], ".png"))
                          n.icon <<- n.icon + 1
                          grid::rasterGrob(0.5, 0.5, image = img, width = 1, height = 1)
                        })
        class(grobs) <- "gList"
        ggplot2:::ggname("image_key", grid::gTree(children = grobs))
      }

      if(all(grepl("^#[[:xdigit:]]{6}$", icons.col)) | all(grepl("^#[[:xdigit:]]{8}$", icons.col))){
        options(ggimage.keytype = "image")
        Dendo <- Dendo + ggnewscale::new_scale_colour() +
          ggimage::geom_image(data = ggdendro::label(dendr),
                              aes(x = x, y = y - y.image,
                                  image = IMAGE, color = get(icons.variable)),
                              size = 0.02,
                              by = "height",
                              key_glyph = draw_key_image_col) +
          scale_color_manual(values = icons.col) +
          guides(color = guide_legend(title=icons.variable,
                                      ncol = ifelse(length(icons)>5, 2, 1),
                                      order = 1))
      } else {
        Dendo <- Dendo + ggnewscale::new_scale_colour() +
          ggimage::geom_image(data = ggdendro::label(dendr),
                              aes(x = x, y = y - y.image,
                                  image = IMAGE),
                              size = 0.02,
                              by = "height")  +
          geom_point(data = ggdendro::label(dendr),
                     aes_string(x = "x", y = "y", color = icons.variable),
                     size = 0, key_glyph = draw_key_image) +
          guides(color = guide_legend(ncol = ifelse(length(icons)>5, 2, 1), order = 1))
      }
      if(all(icons.col == tile.col)){
        Dendo <- Dendo + guides(fill = "none")
      }

    } else if (isTRUE(icons.dendogram) & is.null(icons.variable)){
      cat("The icons.variable argument is missing, icons can't be associated to their respective samples.\n")
    }

    if (isTRUE(count.plot)){
      Dendo <- Dendo + theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 2), size = 12))
    } else {
      Dendo <- Dendo + theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 2), size = 12))
    }

    ifelse(is.null(nudge.bs.label), {nudge.bs.label<-(-max(height)*0.08)}, "")

    Dendo2 <- Dendo +
      ggraph::geom_node_point(aes(x = x + 1), shape = 21, show.legend=FALSE) +
      ggraph::geom_node_text(aes(x = x + 1, filter=!leaf,
                                 label = as.character(round(au, 2))),
                             nudge_y = nudge.bs.label,
                             size = bs.label.size)

    if (abundance.plot == "blank"){
      Dendo2 <- Dendo2 +
        theme(
          axis.text.x  = element_text(angle = 40, vjust = 1, hjust =1, size = sample.txt.size, colour = 'black'),
          axis.title.x = element_text(size = sample.txt.size + 2,
                                      margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.ticks.x = element_line(),
          plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
    }

    Dendogram <- ggplotGrob(Dendo2)

    # Icons height
    if (isTRUE(icons.dendogram) & !is.null(icons.variable)){
      if (abundance.plot == "blank"){
        for (i in 1:length(Dendogram$grobs[[6]]$children[[7]]$children)) {
          Dendogram$grobs[[6]]$children[[7]]$children[[i]]$height <- unit(0.04,"native")
        }
      } else if (!is.null(var.tile)){
        for (i in 1:length(Dendogram$grobs[[6]]$children[[7]]$children)) {
          Dendogram$grobs[[6]]$children[[7]]$children[[i]]$height <- unit(0.1,"native")
        }
      } else {
        for (i in 1:length(Dendogram$grobs[[6]]$children[[6]]$children)) {
          Dendogram$grobs[[6]]$children[[6]]$children[[i]]$height <- unit(0.1,"native")
        }
      }
    }
  }

  # To create read count plot
  if (isTRUE(count.plot)){
    if (isTRUE(count.plot.rel)){
      data.count.plot = transform_sample_counts(data.count.plot, function(x) x/sum(x))
    }
    Melted.data <- psmelt(data.count.plot)

    # Rename NA or unclassified Domain to Unclassified
    Melted.data$Domain[Melted.data$Domain == NA] <- "Unclassified"
    Melted.data$Domain[Melted.data$Domain == "unclassified"] <- "Unclassified"

    # Make sure a "SampleID" variable exists (to change sample names without errors)
    if (!"SampleID" %in% colnames(Melted.data)){
      Melted.data[,"SampleID"] <- Melted.data[,"Sample"]
    }

    # Make Sample names shorter -- IF NEEDED --
    if (!is.null(str2rm)){
      for (i in str2rm){
        Melted.data[,"SampleID"] <- sub(i, "", Melted.data[,"SampleID"])
      }
    }

    # Specified the levels of the Domain categorical variable
    Melted.data$Domain <- factor(Melted.data$Domain, levels = domains)

    if (isTRUE(dendogram)){
      # Re-order "Melted.data" to match "dendr" labels
      dendr2 <- ggdendro::dendro_data(hc, type = "rectangle") # convert to ggplot
      sample_pos_table <- with(dendr2$labels, data.frame(x = x, Sample = as.character(label)))
      Final.count.data <- Melted.data %>% dplyr::left_join(sample_pos_table, by = "Sample")

      # Plot
      Plot <- ggplot(data = Final.count.data, aes(x = x, y = Abundance, fill = Domain)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(y = "Number of reads") +
        scale_x_continuous(breaks = sample_pos_table$x,
                           limits = axis_limits,
                           expand = c(0, 0))
    } else if (!isTRUE(dendogram)){
      Plot <- ggplot(data = Melted.data, aes(x = SampleID, y = Abundance, fill = Domain)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(y = "Number of reads") +
        ggtitle(title)

      if (!is.null(Top.var)){
        ifelse(!is.null(Bot.var), {f<-paste0("~", Top.var, "+", Bot.var)}, {f<-paste0("~", Top.var)})
        Plot <- Plot + facet_wrap(as.formula(f),
                                  scales = "free_x",
                                  nrow = 1) +
          theme(strip.background.x = element_rect(color = "black", linetype = "solid"),
                strip.text.x = element_text(size = strip.text.size)
          )
      }
    }
    if (isTRUE(count.plot.rel)){
      label.format = scales::label_percent(accuracy = 1L)
    } else {
      label.format = function(x) formatC(x, format = "e", digits = 1)
    }
    Count.Plot <- Plot +
      scale_y_continuous(labels = label.format,
                         expand = expansion(mult = c(0, .1))) +
      scale_fill_manual(breaks = domains,
                        values = c("Bacteria" = "#e6194B", "Archaea" = "#3cb44b",
                                   "Eukaryota" = "#4363d8", "Fungi" = "#f7a240",
                                   "Viruses" = "#e8dd07",
                                   "Unclassified" = "#DDDDDD"),
                        name = "Domain",
                        limits = levels(Melted.data$Domain)) +
      theme (
        plot.title = element_text(face = "bold.italic"),
        panel.background = element_rect(fill = "white"),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 2),
                                    size = 12),
        axis.text.y = element_text(size = 10, colour = 'black'),
        panel.border = element_rect(fill = NA, colour = "grey50"),
        panel.grid.major.y = element_line(colour = "grey",  size = 0.5,
                                          linetype = 3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

    if (abundance.plot == "blank"){
      Count.Plot <- Count.Plot +
        theme(
          axis.title.x = element_text(size = sample.txt.size + 2),
          axis.text.x  = element_text(angle = 40, vjust = 1, hjust = 1, size = sample.txt.size, colour = 'black'),
          axis.ticks.x = element_line())
    }
    if (!isTRUE(dendogram)){
      Count.Plot <- wrap.strip.colors(Melted.data, Count.Plot, Top.var = Top.var,
                                      Bot.var = Bot.var, Top.strip.fill = Top.strip.fill,
                                      Bot.strip.fill = Bot.strip.fill, display = FALSE)
    }
  }

  # To create abundance plot
  if (abundance.plot == "pointplot"){
    if (!isTRUE(dendogram)){
      Top <- ggplot(data = Melted.final.object, aes_string(x = "SampleID", y = rank,
                                                           size = "Abundance", col = rank)) +
        geom_point() +
        labs(title = title, x = "Samples", y = rank, size = "Relative Abundance") +
        theme(plot.title = element_text(face = "bold"),
              axis.text.x  = element_text(angle = 40, vjust = 1, hjust = 1, size = sample.txt.size, colour = 'black'),
              axis.title.x = element_text(size = sample.txt.size + 2,
                                          margin = margin(t = 10, r = 0, b = 0, l = 0))
        )

      if (!is.null(Top.var)){
        ifelse(!is.null(Bot.var), {f<-paste0("~", Top.var, "+", Bot.var)}, {f<-paste0("~", Top.var)})
        Top <- Top + facet_wrap(as.formula(f),
                                scales = "free_x",
                                nrow = 1) +
          theme(strip.background.x = element_rect(color = "black", linetype = "solid"),
                strip.text.x = element_text(size = strip.text.size)
          )
      }

      if (isTRUE(count.plot)){
        Top <- Top + theme(strip.background.x = element_blank(),
                           strip.text.x = element_blank(),
                           plot.title = element_blank())
      }
    } else {
      Top <- ggplot(data = Final_data, aes_string(x = "x", y = rank, size = "Abundance", col = rank)) +
        geom_point() +
        scale_x_continuous(breaks = unique(sort(dendr$labels$x)),
                           limits = axis_limits, expand = c(0, 0),
                           labels = as.vector(dendr$labels$shortlabel[order(dendr$labels$x)])) +
        labs(x = "Samples", y = rank, size = "Relative Abundance") +
        theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
              axis.title.y = element_text(size = 12),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(fill = NA, colour = "grey50"))
    }
  } else if (abundance.plot == "barplot"){
    if (!isTRUE(dendogram)){
      Top <- ggplot(data = Melted.final.object, aes_string(x = "SampleID", y = "Abundance", fill = rank)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(title = title, x = "Samples", y = "Relative abundance") +
        theme(plot.title = element_text(face = "bold"),
              axis.text.x  = element_text(angle = 40, vjust = 1, hjust = 1, size = sample.txt.size, colour = 'black'),
              axis.title.x = element_text(size = sample.txt.size + 2,
                                          margin = margin(t = 10, r = 0, b = 0, l = 0)),
              axis.title.y = element_text(size = 12),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(fill = NA, colour = "grey50")
        )
      if (!is.null(Top.var)){
        ifelse(!is.null(Bot.var), {f<-paste0("~", Top.var, "+", Bot.var)}, {f<-paste0("~", Top.var)})
        Top <- Top + facet_wrap(as.formula(f),
                                scales = "free_x",
                                nrow = 1) +
          theme(strip.background.x = element_rect(color = "black", linetype = "solid"),
                strip.text.x = element_text(size = strip.text.size)
          )
      }
      if (isTRUE(count.plot)){
        Top <- Top + theme(strip.background.x = element_blank(),
                           strip.text.x = element_blank(),
                           plot.title = element_blank())
      }
    } else {
      Top <- ggplot(data = Final_data, aes(x = x, y = Abundance, fill = Final_data[,rank])) +
        geom_bar(stat = "identity", position = "stack") +
        scale_x_continuous(breaks = unique(sort(dendr$labels$x)),
                           limits = axis_limits, expand = c(0, 0),
                           labels = as.vector(dendr$labels$shortlabel[order(dendr$labels$x)])) +
        labs(x = "Samples", y = "Relative abundance") +
        theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
              axis.title.y = element_text(size = 12),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(fill = NA, colour = "grey50"))
    }
  }

  if (abundance.plot == "pointplot"){
    Abundance.Plot <- Top +
      scale_size_continuous(labels = scales::label_percent(accuracy = 1L),
                            limits = c(0,1), range = c(0,8),
                            breaks = c(0.01, 0.1, 0.5, 0.75)) +
      scale_color_manual(values = mycolors, guide = FALSE) +
      theme(
        axis.title.x = element_text(size = sample.txt.size + 2),
        axis.text.x  = element_text(angle = 40, vjust = 1, hjust = 1, size = sample.txt.size,
                                    colour = 'black'),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 2),
                                    size = taxa.txt.size + 2),
        axis.text.y = element_text(size = taxa.txt.size, colour = 'black'),
        panel.border = element_rect(fill = NA, colour = "grey50"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text())
  } else if (abundance.plot == "barplot"){
    Abundance.Plot <- Top +
      scale_y_continuous(labels = scales::percent,
                         breaks = c(0, .25, .50, .75, 1),
                         expand = c(0, 0), limits = c(0, 1)) +
      scale_fill_manual(values = mycolors, name = rank,
                        limits = levels(Melted.final.object[,rank])) +
      theme(
        axis.title.x = element_text(size = sample.txt.size + 2),
        axis.title.y = element_text(size = taxa.txt.size + 2),
        axis.text.x  = element_text(angle = 40, vjust = 1, hjust = 1, size = sample.txt.size, colour = 'black'),
        axis.text.y = element_text(size = taxa.txt.size, colour = 'black'),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = taxa.txt.size + 2),
        legend.text = element_text(size = taxa.txt.size)) +
      guides(fill = guide_legend(nrow = length(levels(Melted.final.object[,rank])), ncol = 1))

    if (isTRUE(count.plot)){
      Abundance.Plot <- Abundance.Plot + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 2), size = 12))
    } else {
      Abundance.Plot <- Abundance.Plot + theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 2),
                                                                           size = taxa.txt.size + 2))
    }

  }

  if (merge.by != "None"){
    Abundance.Plot <- Abundance.Plot +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }

  if (!isTRUE(dendogram) & !isTRUE(count.plot) & abundance.plot != "blank"){
    grid.object <- wrap.strip.colors(Melted.final.object, Abundance.Plot, Top.var = Top.var, Bot.var = Bot.var, Top.strip.fill = Top.strip.fill, Bot.strip.fill = Bot.strip.fill)
  } else if (isTRUE(dendogram) & isTRUE(count.plot) & abundance.plot != "blank"){
    grid.object <- cowplot::plot_grid(Dendogram, Count.Plot, Abundance.Plot, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.4,0.3,1))
  } else if (isTRUE(dendogram) & !isTRUE(count.plot) & abundance.plot != "blank"){
    grid.object <- cowplot::plot_grid(Dendogram, Abundance.Plot, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.5,1))
  } else if (!isTRUE(dendogram) & isTRUE(count.plot) & abundance.plot != "blank"){
    grid.object <- cowplot::plot_grid(Count.Plot, Abundance.Plot, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.3,1))
  } else if (!isTRUE(dendogram) & isTRUE(count.plot) & abundance.plot == "blank"){
    grid.object <- cowplot::plot_grid(Count.Plot)
  } else if (isTRUE(dendogram) & !isTRUE(count.plot) & abundance.plot == "blank"){
    grid.object <- cowplot::plot_grid(Dendogram)
  }

  if(!is.null(folder)){
    ifelse(!dir.exists(folder), dir.create(folder, recursive = TRUE), "")
  }

  if (isTRUE(return.plot)){
    return(grid.object)
  } else {
    ggsave(filename = paste0(folder, file.name),
           plot = grid.object, height = 9, width = 17, unit = "in", dpi = 600)
  }
  options(warn = oldw)
}
