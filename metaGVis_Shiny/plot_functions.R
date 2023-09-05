# # Globally defined plot functions
################################################################################
# Count histogram plot code
################################################################################
count.histogram <- function(data, subtitle){
  Table.bact <- as.data.frame(data)
  Table.bact$SUM <- rowSums(Table.bact)
  ggplot(Table.bact, aes(x = SUM)) + 
    geom_histogram(bins = 100) + 
    scale_x_log10(breaks = c(0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(
      title = "Histogram of log10 feature counts",
      subtitle = subtitle,
      y = "Taxa count", x = "Sum of read counts per taxa (log10)") +
    theme (
      plot.title = element_text(face = "bold.italic"),
      axis.text.x  = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
}
################################################################################
# Features per sample bar plot code
################################################################################
filtering_QC <- function(data, subtitle, max.y = 420){
  
  df <- otu_table(data)
  ASVs <- colSums(df != 0)
  ASVs <- as.data.frame(ASVs)
  colnames(ASVs) <- "Total_ASVs"
  ASVs$SampleID <- rownames(ASVs)
  
  # Create bar plot
  plot <- ggplot(data = ASVs, aes(x = SampleID, y = Total_ASVs)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = "Features (OTUs, ASVs, species,...) per sample",
      subtitle = subtitle,
      y = "Features", x = "Samples") +
    scale_y_continuous(breaks = signif(seq.int(0, max.y, length.out = 6),
                                       digits = ifelse(nchar(max.y) == 1,1,
                                                       ifelse(nchar(max.y) == 2,1,
                                                              ifelse(nchar(max.y) == 3,2, 3)))), 
                       expand = c(0, 0), 
                       limits = c(0, max.y)) +
    theme (
      plot.title = element_text(face = "bold.italic"),
      strip.background.x = element_rect(color = "black", linetype = "solid"),
      strip.text.x = element_text(face = "plain", size = 6),
      axis.text.x  = element_text(angle = 40, vjust = 1, hjust = 1, size = 6),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
      panel.grid.major.y = element_line(colour = "grey",  size = 0.5,  linetype = 3),
      panel.grid.major.x = element_blank(), 
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 7.5)
    )
  
  histo <- count.histogram(data@otu_table, subtitle)
  
  Final.plot <- cowplot::plot_grid(plot, histo, nrow = 2, ncol = 1, align = "hv", rel_heights = c(1, 0.7) )
}
################################################################################
# Functional Profile plot code
################################################################################
Alpha_Diversity_plot <- function(physeq,
                                 color_variable = NULL, 
                                 shape_variable = NULL, 
                                 group_variable = NULL, 
                                 boxplot = TRUE,
                                 measures = c("Observed", "Shannon", "Simpson"), 
                                 colors = NULL, 
                                 pair_compare = TRUE, 
                                 hide.ns = FALSE, 
                                 remove.ns = TRUE, 
                                 pair_compare_method = "wilcox.test",
                                 paired = FALSE,
                                 p.adjust = "BH",
                                 title = "",
                                 size.text.samples = 16,
                                 tip.length = 0.01,
                                 ggtheme = theme_bw(),
                                 shown.legend = "colors",
                                 return.p.values = TRUE){
    
    plot.list <- list()
    p.values.list <- list()
    for (i in 1:length(measures)){
      richness <- plot_richness(physeq, x = group_variable, color = color_variable, 
                                shape = shape_variable, measures = measures[i]) +
        labs(title = title, x = "Type of sample") +
        expand_limits(y = 0) +
        ggtheme +
        theme (
          plot.title = element_text(face = "bold.italic", hjust = 0.5, 
                                    margin = margin(t = 0, r = 0, b = 20, l = 0), size = 18),
          axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = size.text.samples),
          axis.text.y  = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 18),
          strip.text.x = element_text(size = 18),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.position = "bottom")
        
      # Apply custom color
      if (!is.null(color_variable)){
        richness <- richness + 
          scale_color_manual(values = colors, name = color_variable) +
          guides(colour = guide_legend(override.aes = list(size=6)))
      }
      
      # Extract legend
      if (!is.null(color_variable) | !is.null(shape_variable)){
        legend <- cowplot::get_legend(richness)
      } else {
        legend <- NULL
      }
      # Extract title
      ext.title <- cowplot::get_title(richness) 
      
      richness$layers <- richness$layers[-1] # Remove points (latter added by geom_jitter)
      
      # Apply custom color
      if (!is.null(color_variable) & isTRUE(boxplot)){
        richness <- richness + 
          geom_boxplot(data = richness$data, alpha = 0.1, show.legend = FALSE,
                       aes(x = richness$data[,group_variable], y = value, 
                           color = richness$data[,color_variable])) 
      } else if (isTRUE(boxplot)){
        richness <- richness + 
          geom_boxplot(data = richness$data, alpha = 0.1, show.legend = FALSE,
                       aes(x = richness$data[,group_variable], y = value)) 
      }
      richness <- richness +
        geom_jitter(position = position_jitter(width = .2)) +
        theme(axis.title.x = element_blank(),
              plot.title = element_blank(),
              legend.position = "none")
      
      if (shown.legend == "colors"){
        richness <- richness +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
      }
      
      if (i != 1){
        richness <- richness + theme(axis.title.y = element_blank())
      }

      if(pair_compare == TRUE){
        # construct and filter the paired comparisons list
        comparisons_list <- unique(as.character(richness$data[,group_variable])) %>% combn(., 2) %>% 
          {lapply(seq_len(ncol(.)), function(x) .[, x])}
        
        if (remove.ns == TRUE){
          # Calculate pvalues
          comparison.formula <- paste0("value", "~", group_variable) %>% as.formula()
          p.values.list[[i]] <- ggpubr::compare_means(comparison.formula,  data = richness$data, 
                                                      method = pair_compare_method, group.by = "variable", 
                                                      paired = paired, p.adjust.method = p.adjust)
          non.ns.pvalues <- p.values.list[[i]] %>% filter(p.adj < 0.05)
          
          if (nrow(non.ns.pvalues) != 0){
            comparisons_list <- list()
            for (k in c(1:nrow(non.ns.pvalues))){
              comparisons_list[[k]] <- c(non.ns.pvalues$group1[k], non.ns.pvalues$group2[k])
            }
            comparisons_list <- unique(comparisons_list)
          } else {
            comparisons_list <- NULL
          }
          
        }
        
        if (!is.null(comparisons_list)){
          richness <- richness + ggpubr::stat_compare_means(
            comparisons = comparisons_list,
            method = pair_compare_method, 
            paired = paired,
            tip.length = tip.length,
            label = "p.signif",
            hide.ns = hide.ns,
            step.increase = 0.1,
            symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                               symbols = c("****", "***", "**", "*", "ns")))
        }
      }
      plot.list[[i]] <- richness
    }    
    
    plot <- ggpubr::ggarrange(plotlist = plot.list, ncol = length(plot.list), nrow = 1, 
                              align = "hv")
    plot <- ggpubr::ggarrange(ext.title, plot, ncol = 1, nrow = 2, heights = c(0.1,0.8))
    if (shown.legend %in% c("both", "colors")){
      if (!is.null(color_variable) | !is.null(shape_variable)){
        plot <- ggpubr::ggarrange(plot, legend, ncol = 1, nrow = 2, heights = c(1, 0.1))
      } else {
        plot <- ggpubr::annotate_figure(plot, 
                                        bottom = grid::textGrob(group_variable, gp = grid::gpar(cex = 1.5)))
      }
    } else {
      plot <- ggpubr::annotate_figure(plot, 
                                      bottom = grid::textGrob(group_variable, gp = grid::gpar(cex = 1.5)))
    }

    
    if(isTRUE(return.p.values)){
      ls <- list(plot, p.values.list)
      return(ls)
    } else {
      ls <- list(plot, NULL)
      return(ls)
    }
}
################################################################################
# Whittaker plot code
################################################################################
whittaker <- function(physeq, 
                      group = NULL,
                      max.x = 1000,
                      ncol = 3, 
                      title = "Boxplot-based rank abundance curve",
                      ggtheme = theme_bw()){
  # Convert taxa counts in each sample to relative abundance
  rel.physeq = transform_sample_counts(physeq, function(x) 1e+02 * x/sum(x))
  melted.data = psmelt(rel.physeq)
  melted.data = filter(melted.data, Abundance > 0)
  # Calculate the mean of Abundance grouping by group and taxa
  df = melted.data %>% group_by(OTU, melted.data[,group]) %>% mutate(Mean.abundance = mean(Abundance))
  # Change colnames
  col.names <- colnames(df)
  col.names[(length(col.names)-1)] <- "Group"
  colnames(df) <- col.names
  # Order data frame
  df = df[with(df, order(Group, -Mean.abundance)), ]
  # Add Rank
  df = df %>% group_by(Group) %>% mutate(Rank = as.factor(dense_rank(-Mean.abundance)))
  # Plot
  plot <- ggplot(df) + 
    geom_boxplot(aes(x=Rank, y=Abundance), size = 1) +
    facet_wrap(as.formula(paste0("~", group)), scales = "free_x", ncol = ncol) +
    ggtheme +
    labs(
      title = title, 
      subtitle = "For each ASV rank a boxplot helps visualising the within group variance",
      x = "Species Abundance Rank",
      y = "Relative abundance (log10)") +
    scale_x_discrete(breaks = signif(seq.int(1, max.x, length.out = 6),
                                     digits = ifelse(nchar(max.x) == 1,1,
                                                     ifelse(nchar(max.x) == 2,1,
                                                            ifelse(nchar(max.x) == 3,2, 3))))) +
    scale_y_log10() 
  return(plot)
}
################################################################################
# Functional Profile plot code
################################################################################
FunctionalProfile <- function(physeq, 
                              taxa.rank, 
                              func.rank, 
                              variable, 
                              show.unclass.taxa = TRUE, 
                              show.unclass.func = FALSE, 
                              taxa.color.palette,
                              function.color.palette,
                              plot.mat,
                              labelsize,
                              numbersize,
                              return.legend = FALSE){
  set.seed(711)
  # Create the data frame used for the chordiagram to know the number of colors to generate
  m <- as.data.frame(physeq@otu_table)
  t <- as.data.frame(physeq@tax_table)
  t[,taxa.rank][t[,taxa.rank] == "unclassified"] <- "Unclassified taxa"
  t[,func.rank][t[,func.rank] == "unclassified"] <- "Unclassified functions"
  physeq@tax_table <- tax_table(as.matrix(t))
  chord.df = data.frame(from = gsub("[._].*","", gsub(".__", "", t[,taxa.rank])), 
                        to = t[,func.rank],
                        value = rowMeans(m))
  chord.df %<>% group_by(to) %>% mutate(Mean = sum(value)) %>% ungroup() %>% arrange(desc(Mean))
  taxa <- unique(chord.df$from)
  functions <- unique(chord.df$to)
  
  # Create color palette for specific taxa and function names
  if (isTRUE(show.unclass.taxa)){
    grid.col1 = unname(Polychrome::createPalette((length(taxa)-1), taxa.color.palette))
    grid.col1 = c("#DDDDDD", grid.col1)
    names(grid.col1) <- taxa
  } else {
    grid.col1 = unname(Polychrome::createPalette((length(taxa)-1), taxa.color.palette))
    taxa <- taxa[taxa != "Unclassified taxa"]
    names(grid.col1) <- taxa
  }
  if (isTRUE(show.unclass.func)){
    grid.col2 = unname(Polychrome::createPalette((length(functions)-1), function.color.palette))
    grid.col2 = c("#B0B0B0", grid.col2)
    names(grid.col2) <- functions
  } else {
    grid.col2 = unname(Polychrome::createPalette((length(functions)-1), function.color.palette))
    functions <- functions[functions != "Unclassified functions"]
    names(grid.col2) <- functions
  }
  grid.col = c(grid.col1, grid.col2)
  
  if (!isTRUE(return.legend)){
    # Plot chordiagrams
    pdf(NULL)
    dev.control(displaylist = "enable")
    circlize::circos.clear()
    par(xpd = TRUE, mar = c(1, 0, 1, 0), mfrow = plot.mat) 
    # Set xpd=NA and expand the top and bottom margins
    matrices <- list()
    lev <- levels(get_variable(physeq, variable))
    for (i in lev){
      oldDF <- as(sample_data(physeq), "data.frame")
      newDF <- subset(oldDF, eval(str2expression(variable)) == i)
      sub.phylo <- physeq
      sample_data(sub.phylo) <- sample_data(newDF)
      chord.df = data.frame(from = gsub("[._].*","",gsub(".__", "", tax_table(sub.phylo)[,taxa.rank])),
                            to = tax_table(sub.phylo)[,func.rank],
                            value = rowMeans(sub.phylo@otu_table))
      colnames(chord.df) <- c("from", "to", "value")
      if (!isTRUE(show.unclass.func)){
        chord.df <- chord.df[chord.df$to != "Unclassified functions", ]
      } 
      if (!isTRUE(show.unclass.taxa)){
        chord.df <- chord.df[chord.df$from != "Unclassified taxa", ]
      }
      chord.df <- chord.df%>%mutate(value = value / sum(value) * 100)
      
      # Create matrix
      mat <- reshape2::dcast(chord.df, from ~ to, sum)
      rownames(mat) <- mat$from
      mat <- as.matrix(mat[-1]) # convert to matrix and remove "from" column
      if (nrow(mat) > 1){
        # remove taxa with rowSums equal to 0
        mat %<>% .[rowSums(.) != 0, ] 
        # Sort taxa (rows) by total abundance (decrasing order)
        mat %<>% .[names(sort(rowSums(.), T)),]
        # Sort functions (columns) by total abundance 
        matrices[[i]] <- mat[,names(sort(colSums(mat), T))]
      } else {
        m <- mat[,names(sort(colSums(mat), T))]
        m2 <- matrix(m, dimnames=list(NULL, rownames(mat))) # Colname lost in previous step
        rownames(m2) <- names(sort(colSums(mat), T)) # Rownames lost in previous step
        matrices[[i]] <- t(m2) # translocation occurs in previous steps when only one column is present
      }
      
      ### Basic circos graphic parameters
      circlize::circos.par(cell.padding = c(0, 0, 0, 0))
      circlize::chordDiagram(matrices[[i]], annotationTrack = c("grid"), grid.col = grid.col, big.gap = 20)
      for(si in circlize::get.all.sector.index()) {
        xlim = circlize::get.cell.meta.data("xlim", sector.index = si, track.index = 1)
        ylim = circlize::get.cell.meta.data("ylim", sector.index = si, track.index = 1)
        circlize::circos.text(mean(xlim), 1.8, ifelse(round(xlim[2],1) >= 1, round(xlim[2],1),"."),
                              sector.index = si, cex = numbersize, track.index = 1, facing = "clockwise",
                              niceFacing = TRUE, col = "black")
      }
      abline(h = 0, lty = 2, col = "#00000080")
      title(main = i, line = -6, cex.main = labelsize)
      circlize::circos.clear()
    }
    plot <- grDevices::recordPlot()
    invisible(dev.off())
    return(plot)
  } else {
    # Create legends
    lgd_taxa = ComplexHeatmap::Legend(at = taxa, legend_gp =  grid::gpar(fill = grid.col1),
                                      title_position = "topleft", title = paste(taxa.rank, "(Bottom)"), ncol = 5, 
                                      labels_gp = grid::gpar(fontsize = 6))
    
    lgd_functions = ComplexHeatmap::Legend(at = functions, legend_gp = grid::gpar(fill = grid.col2),
                                           title_position = "topleft", title = "Functions (Top)", ncol = 3, 
                                           labels_gp = grid::gpar(fontsize = 6))
    
    lgd <- ComplexHeatmap::packLegend(lgd_functions, lgd_taxa)
    
    # Plot legends
    tmp.legend <- tempfile("Legends", fileext = ".png")
    plot.new()
    png(tmp.legend, height = 4, width = 7, units = "in", res = 600)
    par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1, 1)) 
    ComplexHeatmap::draw(lgd)
    dev.off()
    rl <- png::readPNG(tmp.legend)
    gl <- grid::rasterGrob(rl)
    legend <- cowplot::plot_grid(gl)
    return(legend)
  }
}
################################################################################
# PERMANOVA report code
################################################################################
PERMANOVA_report <- function(input_data, 
                             distance, 
                             formula, 
                             p.adj.method = "BH", 
                             sqrt.dist = F,
                             add = F){
  var <- unlist(strsplit(formula, " "))[1]
  df <- data.frame(sample_data(input_data))
  meta_dist <- merge(df, as.matrix(distance), by = "row.names")[-1] # merge metadata with distance matrix
  rownames(meta_dist) <- meta_dist$SampleID
  beta <- vegan::betadisper(distance, df[,var]) # Homogeneity of dispersion test (first step)
  homogeneity_of_dispersion <- vegan::permutest(beta) # Homogeneity of dispersion - Overall test
  permanova <- vegan::adonis2(as.formula(paste("distance ~ ", formula)), 
                              data = df, 
                              sqrt.dist = sqrt.dist,
                              add = add) # PERMANOVA test
  group.levels <- levels(as.factor(df[,var])) # extract levels
  # Multiple Pairwise Comparison to identify the centroid groups which are different -- if any --
  if (permanova[, "Pr(>F)"][1] < 0.05){
    other.levels <- group.levels
    pairs = c()
    permutations = c()
    F.value =c()
    R2 = c()
    p.value = c()
    for (l in group.levels){
      other.levels <- other.levels[other.levels != l]
      for (i in other.levels){
        pair.levels <- c(l, i)
        pairs <- c(pairs, paste(l,'vs',i))
        mod.df <- meta_dist %>% filter(meta_dist[,var] == l | meta_dist[,var] == i) # filter 2 levels of interest
        mod.dist <- mod.df %>% select(all_of(.[["SampleID"]])) %>% as.dist() # convert to distance
        p.pair <- vegan::adonis2(as.formula(paste("mod.dist ~ ", formula)), 
                                 data = mod.df,
                                 sqrt.dist = sqrt.dist,
                                 add = add) # PERMANOVA test
        permutations <- c(permutations, strsplit(capture.output(p.pair)[4], " ")[[1]][4]) # extract nº of permutations
        F.value <- c(F.value, p.pair$F[1]) # extract F value
        R2 <- c(R2, p.pair$R2[1]) # extract R2 value
        p.value <- c(p.value, p.pair$`Pr(>F)`[1]) # extract P-value

      }
    }
    p.adj <- p.adjust(p.value, method = p.adj.method) # adjust P-value
    sig = c(rep('', length(p.adj)))
    sig[p.adj <= 0.05] <-'*'
    sig[p.adj <= 0.01] <-'**'
    sig[p.adj <= 0.001] <-'***'
    permanova.pair <- data.frame(pairs, permutations, F.value, R2, p.value, p.adj, sig)
    rownames(permanova.pair) <- NULL
  }
  # Multiple Pairwise Comparison to identify which groups have difference in dispersion between groups -- if any --
  if (homogeneity_of_dispersion$tab[, "Pr(>F)"][1] < 0.05){
    pairwise.HSD <- TukeyHSD(beta) # pairwise comparison
    df <- as.data.frame(pairwise.HSD$group)
    if (!is.null(pairs)){
      rownames(df) <- pairs
    }
    df$sig <- rep("", nrow(df))
    df$sig[df$`p adj` <= 0.05] <-'*'
    df$sig[df$`p adj` <= 0.01] <-'**'
    df$sig[df$`p adj` <= 0.001] <-'***'
    pairwise.HSD$group <- df
  }
  stat.list <- list()
  stat.list[[1]] <- homogeneity_of_dispersion
  if (homogeneity_of_dispersion$tab[, "Pr(>F)"][1] < 0.05){
    stat.list[[2]] <- pairwise.HSD
  } else {
    stat.list[[2]] <- "No overall significant differences observed."
  }
  stat.list[[3]] <- permanova
  if (permanova[, "Pr(>F)"][1] < 0.05){
    stat.list[[4]] <- permanova.pair
  } else {
    stat.list[[4]] <- "No overall significant differences observed."
  }
  return(stat.list)
}
################################################################################
# 2D ordination plot code
################################################################################
ordination2D <- function(physeq, 
                 norm = "TSS",
                 transf = "CLR_runif",
                 method = "PCA",
                 distance = "euclidean",
                 color_var, 
                 shape_var = NULL, 
                 formula = "",
                 colors = "Set2", 
                 display.names = TRUE, 
                 display.ellipses = TRUE, 
                 group.ellipses, 
                 ellipse.fill = "Pastel2",
                 title = "",
                 axes.size = 12,
                 leg.size = 12,
                 str2rm = "",
                 point.size = 3,
                 names.size = 2.5,
                 legend.position = "bottom",
                 p.adj.method = "BH",
                 sqrt.dist = F,
                 add = F){

  if(isTRUE(display.ellipses) & missing(group.ellipses)){
    group.ellipses <- color_var
  }
  input_data <- normalisation(physeq, norm, color_var)
  input_data <- transformation(input_data, transf)
  if (title == ""){
    if (norm == "TSS" | norm == "Rarefy"){
      title = paste0("Normalisation: ", norm, " (computed on relative abundances)")
    } else {
      title = paste0("Normalisation: ", norm, " (computed on absolute abundances)")
    }
  }
  # Calculate dissimilarity matrices
  if (distance == "philr"){
    input_data@otu_table <- input_data@otu_table + 1
    philr <- philr::philr(t(as(input_data@otu_table, "matrix")), input_data@phy_tree,
                          part.weights='enorm.x.gm.counts',
                          ilr.weights='blw.sqrt')
    dist.m <- dist(philr, method = "euclidean")
  } else {
    dist.m <- phyloseq::distance(input_data, method = distance)
  }

  # Statistics
  if(formula != ""){
    var = formula
  } else {
    var = color_var
  }
  report <- PERMANOVA_report(input_data, distance = dist.m, 
                             formula = var, p.adj.method = p.adj.method,
                             sqrt.dist = sqrt.dist,
                             add = add)
  
  # PCA/PCoA analysis for 6 components (to visualise the components check plot)
  eig <- NULL
  ifelse(method == "NMDS",dim <- 3, dim <- 6) 
  if (method %in% c("PCA", "PCoA")){
    ord <- cmdscale(dist.m, k = dim, eig = TRUE)
    for (i in 1:dim) {
      eig[i] <- (ord$eig[i] / sum(ord$eig))
    }
  } else if (method == "NMDS"){
    for (i in 1:dim) {
      ord <- vegan::metaMDS(dist.m, k = i, trymax = 30)
      eig[i] <- ord$stress
    }
  } else if (method %in% c("CCA", "RDA", "CAP") & formula != ""){
    ord <- ordinate(physeq, method, distance, as.formula(paste("1 ~ ", formula)))
    for (i in 1:dim) {
      ifelse(method == "CAP", eig[i] <- ord$CCA$eig, eig[i] <- ord$CA$eig)
    }
  }
  
  # Stress plot -> to find the breakpoint" that can instruct selection of a minimum number of dimensions
  if(method == "NMDS"){
    stress.plot <- ggplotify::as.grob(function(){
      plot.new()
      par(mar=c(6.1,5.1,4.1,2.1), mgp = c(2, 1, 0))
      plot(seq(1, dim, 1), eig, 
           main = paste0("Stress value in tested dimensions ", "(", norm, ")"), 
           xlab = "Dimension", ylab = "Stress", ylim = c(0, max(eig) + .05), type = "n")
      lines(seq(1, dim,1), eig)
      points(seq(1, dim,1), eig, pch=21, col = "black", bg = "red", cex = 1.2)
      abline(0.2, 0, col = "red", lty = 2)
      abline(0.1, 0, col = "orange", lty = 2)
      abline(0.05, 0, col = "green", lty = 2)
    })
  } else {
    if (method == "PCoA"){xlab="Principal coordinates"} else {xlab="Principal components"}
    stress.plot <- ggplotify::as.grob(function(){
      plot.new()
      par(
        mar = c(6.1,5.1,4.1,2.1),
        mgp = c(3, 1, 0))
      bp <- barplot(eig, 
                    main = paste0("Scree plot ", "(", norm, ")"), 
                    xlab = xlab, ylab="Relative eigenvalues", ylim = c(0, max(eig) + .05), xaxt = "n")
      bp
      axis(1, at = bp[1:dim], labels = seq(1, dim, by = 1))
      abline(0.1, 0, col = "black", lty = 2)
    })
  }
  if(method == "NMDS"){
    # Shepard Diagram -> to show the relationship between the actual dissimilarities between objects (from the original dissimilarity matrix) and the ordination distances (i.e. the distances on the final plot).
    shepard.plot <- ggplotify::as.grob(function() {
      plot.new()
      par(mar=c(3,3,2,1), mgp = c(2, 1, 0))
      vegan::stressplot(ord, 
                        main = paste0("Shepard Diagram (", norm, " transformation)"))
    })
  }

  # Simplify sample labels and add them -- IF NEEDED --
  labels <- colnames(input_data@otu_table)
  for (i in str2rm){
    labels <- sub(i, "", labels)
  }
  
  # Color palette
  ifelse(!length(colors) > 1, {col = RColorBrewer::brewer.pal(nrow(unique(input_data@sam_data[,color_var])), colors)}, 
         {col = colors})
  ifelse(!length(ellipse.fill) > 1, {fill = RColorBrewer::brewer.pal(nrow(unique(input_data@sam_data[,group.ellipses])), ellipse.fill)}, 
         {fill = ellipse.fill})

  # MDS plot
  if (method == "PCoA"){
    ord = ape::pcoa(dist.m)
    p = phyloseq::plot_ordination(input_data, ord, type = "samples", 
                        color = color_var, 
                        shape = shape_variable) +
      labs(x = paste0("PCo1 [",round(eig[1]*100,1),"%]"), 
           y = paste0("PCo2 [",round(eig[2]*100,1),"%]"))
  } else if (method == "PCA"){
    p = phyloseq::plot_ordination(input_data, ord, type = "samples", 
                        color = color_var, 
                        shape = shape_var) +
      labs(x = paste0("PC1 [",round(eig[1]*100,1),"%]"), 
           y = paste0("PC2 [",round(eig[2]*100,1),"%]"))
  } else if (method == "NMDS"){
    p = plot_ordination(input_data, ord, type = "samples", 
                        color = color_var, 
                        shape = shape_var) +
      labs(x = paste0("NMDS1 (stress: ", round(eig[2],2),")"))
  }

  p = p + geom_point(size = point.size) + scale_colour_manual(values = col) + ggtitle(title) +
    theme(legend.position = legend.position)
  if (isTRUE(display.names)){
    p <- p + ggrepel::geom_text_repel(size = names.size, segment.alpha = 0.5, aes(label=labels))
  }
  
  if (isTRUE(display.ellipses)){
    meta <- as.matrix(input_data@sam_data)[,group.ellipses]
    counts <- dplyr::count(data.frame(meta), meta) 
    if (any(counts <= 3)){
      p <- p + geom_mark_ellipse(aes_string(color = group.ellipses,
                                            fill = group.ellipses, 
                                            group = group.ellipses))  +
        scale_fill_manual(values = fill)
    } else {
      p <- p + stat_ellipse(aes_string(group = group.ellipses, fill = group.ellipses),
                            geom = "polygon",type = "t", linetype = 2, alpha = 0.2) +
        scale_fill_manual(values = fill)
    }
  }
  plot.list <- list()
  plot.list[[1]] <- p + theme(axis.title = element_text(size = axes.size + 2),
                              axis.text = element_text(size = axes.size),
                              legend.title = element_text(size = leg.size + 2),
                              legend.text = element_text(size = leg.size),
                              legend.key.size = unit(axes.size * 0.06, 'cm'))
  if (method == "NMDS"){
    plot.list[[2]] <- ggpubr::ggarrange(stress.plot, shepard.plot, ncol = 2, align = "hv")
  } else {
    plot.list[[2]] <- stress.plot
  }
  plot.list[[3]] <- report
  plot.list[[4]] <- ord
  return(plot.list)
}
################################################################################
# Functions for Venn/Euler diagram
################################################################################
# Function from the eulerr package
# For more info: https://rdrr.io/github/jolars/eulerr/src/R/utils.R#sym-mix_colors
mix_colors <- function(rcol_in) {
  rgb_in <- t(grDevices::col2rgb(rcol_in))
  lab_in <- grDevices::convertColor(rgb_in, "sRGB", "Lab", scale.in = 255)
  mean_col <- colMeans(lab_in)
  rgb_out <- grDevices::convertColor(mean_col, "Lab", "sRGB", scale.out = 1)
  grDevices::rgb(rgb_out)
}
# Code adapted from the internal function plot_venn of the ggVennDiagram packages
# For more info: https://github.com/gaospecial/ggVennDiagram/blob/master/R/ggVennDiagram.R
plot_venn <- function(physeq,
                      group,
                      fraction = 0,
                      show_intersect = FALSE,
                      relative_abund = TRUE,
                      heatmap = TRUE,
                      edge_color = "black",
                      fill_color = NULL,
                      gradient2colors = c("#F8C9BE", "#FF0000"), 
                      title = NULL,
                      label_color = "black",
                      label_size = NA,
                      label = "both",
                      number_geom = "text",
                      number_alpha = 0.5,
                      number_color = "black",
                      number_size = NA,
                      number_digit = 0,
                      number_txtWidth = 60,
                      edge_lty = "solid",
                      edge_size = 1,
                      leg.txt.size = 10,
                      leg.key.height = 1,
                      leg.key.width = 1,
                      ...){

  ps_melted <- phyloseq::psmelt(physeq)
  # Aggregate sample abundances (mean) per sample type (group/selected variable)
  ps_agg <- stats::aggregate(as.formula(paste("Abundance ~ OTU +", group)), data = ps_melted,
                             function(x) (sum(x > 0)/length(x) >= fraction) * mean(x))
  # Reshape the matrix to show sample type as columns and taxaIDs/OTUs as rows
  ps_mat <- reshape2::dcast(as.formula(paste("OTU ~ ", group)), data = ps_agg, 
                            value.var = "Abundance")
  rownames(ps_mat) <- ps_mat[, 1] # Add rownames
  ps_mat <- ps_mat[, -1] # Remove 1st column
  ps_mat <- ps_mat[rowSums(ps_mat) > 0, ] # Remove rows with a sum of 0
  ps_mat_bin <- (ps_mat>0)*1 # Get a binary matrix
  
  # Create a list of sets that works with the function RVenn::Venn
  mdf <- merge(ps_mat_bin, tax_table(physeq), by = 0)
  ranks <- colnames(tax_table(physeq))
  mdf$taxa <- make.unique(apply( mdf[ , ranks] , 1 , paste, collapse = ";"))
  mdf$taxa <- gsub(" ", "_", mdf$taxa)
  # Create final matrix
  ldf <- ps_mat_bin
  # Replace 1s by the taxa rank selected
  ldf[] <- replace(mdf$taxa[row(ldf)], !ldf, 0) 
  k <- list()
  for (i in colnames(ldf)){
    k[[i]]  <- ldf[,i][ldf[,i] != 0]
  }
  venn <- RVenn::Venn(k)
  data <- ggVennDiagram::process_data(venn)
  if(label == "abundance"){
    # Aggregate by group the abundance of the OTUs found in each intersection
    weight <- vector()
    for (i in data@region$name){
      lev <- unlist(strsplit(i,"[.][.]"))
      df <- intersect.df(physeq, group, levels = lev, 
                         relative = relative_abund, weighted = T)
      sum <- sum(df[,1:length(lev)])
      weight <- c(weight, sum)
    }
    data@region$abundance <- round(weight, number_digit)
  }
  
  if(!is.null(fill_color) & !isTRUE(heatmap)){
    n_fills <- length(fill_color)
    n_id <- length(data@region$id)
    values <- strsplit(data@region$id,"")
    mix.fills <- vector()
    if (n_fills < n_id) {
      for (i in 1:n_id) {
        mix.fills[i] <- mix_colors(fill_color[as.numeric(unlist(values[[i]]))])
      }
    }
    p <- ggplot() + geom_sf(fill=mix.fills, data = data@region)
  } else if(label == "abundance"){
    p <- ggplot() + geom_sf(aes_string(fill="abundance"), data = data@region)
  } else {
    p <- ggplot() + geom_sf(aes_string(fill="count"), data = data@region)
  }
    p <- p +
      geom_sf(color = edge_color, data = data@setEdge, show.legend = F,
              lty = edge_lty, size = edge_size) +
      geom_sf_text(aes_string(label = "name"), data = data@setLabel,
                   size = label_size,
                   color = label_color) +
      theme_void() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust=0.5, size=18, margin = margin(0,0,30,0))) +
      coord_sf(clip = 'off')
    
    if (label != "none" & show_intersect == FALSE){
      region_label <- data@region %>%
        dplyr::filter(.data$component == "region") %>%
        dplyr::mutate(percent = paste(round(.data$count*100/sum(.data$count),
                                            digits = number_digit),"%", sep=""),
                      both = paste(.data$count,paste0("(",.data$percent,")"),sep = "\n")
        )
      if (label == "abundance"){
        region_label <- region_label %>%
          dplyr::mutate(abundance = paste0(data@region$abundance, "%")
          )
      }
      if (number_geom == "label"){
        p <- p + geom_sf_label(aes_string(label=label),
                               data = region_label,
                               alpha=number_alpha,
                               color = number_color,
                               size = number_size,
                               lineheight = 0.85,
                               label.size = NA)
      }
      if (number_geom == "text"){
        p <- p + geom_sf_text(aes_string(label=label),
                              data = region_label,
                              alpha=number_alpha,
                              color = number_color,
                              size = number_size,
                              lineheight = 0.85)
      }
    }
    if(heatmap){
      p <- p +
        scale_fill_gradient(low=gradient2colors[1],high = gradient2colors[2]) +
        theme(legend.key.height = unit(leg.key.height, 'cm'),
              legend.key.width = unit(leg.key.width, 'cm'),
              legend.text = element_text(size=leg.txt.size))
    }
    if (show_intersect == TRUE){
      items <- data@region %>%
        dplyr::rowwise() %>%
        dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, collapse = " "),
                                               width = number_txtWidth)) %>%
        sf::st_as_sf()
      label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
      if(!is.null(fill_color) & !isTRUE(heatmap)){
        p <- ggplot(items) +
          geom_sf(fill=mix.fills)
      } else{
        p <- ggplot(items) +
          geom_sf(aes_string(fill="count"))
      }
      p <- p +
        geom_sf_text(aes_string(label = "name"),
                     data = data@setLabel,
                     inherit.aes = F)
      if(label[1] == "abundance"){
        p <- p + geom_text(aes_string(label = "abundance", text = "text"),
                           x = label_coord[,1],
                           y = label_coord[,2],
                           show.legend = FALSE) +
          theme_void()
      } else {
        p <- p + geom_text(aes_string(label = "count", text = "text"),
                           x = label_coord[,1],
                           y = label_coord[,2],
                           show.legend = FALSE) +
          theme_void()     
      }
      if(heatmap){
        p <- p +
          scale_fill_gradient(low=gradient2colors[1],high = gradient2colors[2])
      }
      p <- p +
        ggtitle(title) +
        theme(plot.title = element_text(hjust=0.5, size=18, margin = margin(0,0,30,0)))
      ax <- list(showline = FALSE)
      p <- plotly::ggplotly(p, tooltip = c("text")) %>%
        plotly::layout(xaxis = ax, yaxis = ax)
    }
    p
  }
# Function to output multiple venn/euler diagrams together as well as diagrams
# from two different packages (eulerr and ggVennDiagram)
multi_venn <- function(physeq, 
                       group, 
                       euler = FALSE, 
                       relative_abund = TRUE, 
                       package = "eulerr", 
                       fraction = 0, 
                       label = "percent", 
                       show2plots = FALSE, 
                       heatmap = TRUE,
                       label_color ="black", 
                       label_size = 1.2, 
                       number_digit = 1,
                       number_color = "black", 
                       number_size = 1.3, 
                       fill_color = "Set2", 
                       fill_alpha = 0.7,
                       edge_color = "black", 
                       edge_size = 1, 
                       edge_lty = "solid", 
                       gradient2colors = c("#F8C9BE", "#FF0000"), 
                       plotly = FALSE,
                       plot.align = "horizontal",
                       title1 = "default", 
                       title2 = "default",
                       leg.txt.size = 10,
                       leg.key.height = 1,
                       leg.key.width = 1,...){
  
  if (!isTRUE(show2plots)){
    if (label!="abundance" & title1 == "default"){
      title1 = "Overlaps representing the number of ASVs"
    } else if(label=="abundance" & title1 == "default"){
      title1 = "Overlaps weighted by abundance"
    } 
  } else {
    if (label[1]!="abundance" & title1 == "default"){
      title1 = "Overlaps representing the number of ASVs"
    } else {
      title1 = "Overlaps weighted by abundance"
    } 
    if (label[2]!="abundance" & title2 == "default"){
      title2 = "Overlaps representing the number of ASVs"
    } else {
      title2 = "Overlaps weighted by abundance"
    } 
  }
  
  ps_melted <- phyloseq::psmelt(physeq)
  # Aggregate sample abundances (mean) per sample type (group/selected variable)
  ps_agg <- stats::aggregate(as.formula(paste("Abundance ~ OTU +", group)), data = ps_melted,
                             function(x) (sum(x > 0)/length(x) >= fraction) * mean(x))
  # Reshape the matrix to show sample type as columns and taxaIDs/OTUs as rows
  ps_mat <- reshape2::dcast(as.formula(paste("OTU ~ ", group)), data = ps_agg, 
                            value.var = "Abundance")
  
  rownames(ps_mat) <- ps_mat[, 1] # Add rownames
  ps_mat <- ps_mat[, -1] # Remove 1st column
  ps_mat <- ps_mat[rowSums(ps_mat) > 0, ] # Remove rows with a sum of 0
  ps_mat_bin <- (ps_mat>0)*1 # Get a binary matrix
  
  
  if(isTRUE(show2plots) & label[1] == "abundance"){
    weights.1 = rowMeans(ps_mat)
  } else if (isTRUE(show2plots) & label[1] != "abundance") {
    weights.1 = NULL
  }
  if(isTRUE(show2plots) & label[2] == "abundance"){
    weights.2 = rowMeans(ps_mat)
  } else if (isTRUE(show2plots) & label[2] != "abundance") {
    weights.2 = NULL
  }
  if(!isTRUE(show2plots) & label[1] == "abundance"){
    weights = rowMeans(ps_mat)
  } else {
    weights = NULL
  }
  gg.label <- gsub("counts", "count", label)
  eulerr.label <- list()
  ifelse(label[1]=="both", eulerr.label[[1]] <- c("counts", "percent"),
         ifelse((label[1]=="abundance" & isTRUE(relative_abund)), eulerr.label[[1]] <- "percent", 
                ifelse((label[1]=="abundance" & !isTRUE(relative_abund)), eulerr.label[[1]] <- "counts",
                       eulerr.label[[1]] <- label[1])))
  if (isTRUE(show2plots)){
    ifelse(label[2]=="both", eulerr.label[[2]] <- c("counts", "percent"),
           ifelse((label[2]=="abundance" & isTRUE(relative_abund)), eulerr.label[[2]] <- "percent", 
                  ifelse((label[2]=="abundance" & !isTRUE(relative_abund)), eulerr.label[[2]] <- "counts",
                         eulerr.label[[2]] <- label[2])))
  }
  
  ifelse(length(fill_color) == 1, {col=RColorBrewer::brewer.pal(ncol(ps_mat), fill_color)}, {col = fill_color})
  
  args <- list(physeq, group = group, fraction = fraction, show_intersect = plotly, 
               relative_abund = relative_abund, edge_color = edge_color, fill_color = col, heatmap = heatmap,
               number_color = number_color, number_size = number_size + 3,
               number_geom = "text", number_alpha = 1, label_color = label_color,
               label_size = label_size + 3, number_digit = number_digit, 
               number_txtWidth = 40, edge_lty = edge_lty, edge_size = edge_size, title = title1,
               leg.txt.size = leg.txt.size, leg.key.height = leg.key.height, leg.key.width = leg.key.width,
               gradient2colors = gradient2colors)
  
  if (!isTRUE(euler) & isTRUE(show2plots)){
    if (package == "ggVennDiagram"){
      plot.1 <- do.call(plot_venn, c(args, label = gg.label[1]))
      plot.2 <- do.call(plot_venn, c(args, label = gg.label[2]))
    } else if (package == "eulerr"){
      df.1 <- eulerr::venn(ps_mat_bin, weights = weights.1)
      df.1$original.values <- round(df.1$original.values, number_digit)
      df.2 <- eulerr::venn(ps_mat_bin, weights = weights.2)
      df.2$original.values <- round(df.2$original.values, number_digit)
    }
  } else if (!isTRUE(euler) & !isTRUE(show2plots)){
    if (package == "ggVennDiagram"){
      plot <- do.call(plot_venn, c(args, label = gg.label))
    } else if (package == "eulerr"){
      df <- eulerr::venn(ps_mat_bin, weights = weights)
      df$original.values <- round(df$original.values, number_digit)
    }
  } else if (isTRUE(euler) & isTRUE(show2plots)){
    df.1 <- eulerr::euler(ps_mat_bin, weights = weights.1)
    df.1$original.values <- round(df.1$original.values, number_digit)
    df.2 <- eulerr::euler(ps_mat_bin, weights = weights.2)
    df.2$original.values <- round(df.2$original.values, number_digit)
  } else if (isTRUE(euler) & !isTRUE(show2plots)){
    df <- eulerr::euler(ps_mat_bin, weights = weights)
    df$original.values <- round(df$original.values, number_digit)
  } 
  if (show2plots){
    if (package == "eulerr"){
      class(df.1)="euler" # To avoid "Error in drawVennDiagram"
      plot.1 <- plot(df.1, quantities = list(col=number_color, type=eulerr.label[[1]], cex=number_size), 
                     edges = list(col=edge_color, lty=edge_lty, lwd=edge_size),
                     labels = list(cex=label_size, col=label_color), fills = col, alpha = fill_alpha,
                     main = title1,...)
      class(df.2)="euler"
      plot.2 <- plot(df.2, quantities = list(col=number_color, type=eulerr.label[[2]], cex=number_size), 
                     edges = list(col=edge_color, lty=edge_lty, lwd=edge_size),
                     labels = list(cex=label_size, col=label_color), fills = col, alpha = fill_alpha,
                     main = title2,...)
    } 
    if (plot.align == "horizontal"){
      Final.plot <- cowplot::plot_grid(NULL, plot.1, NULL, plot.2, NULL, ncol = 5, 
                                       rel_widths = c(0.1, 2, 0.3, 2, 0.1), 
                                       labels = c("","(a)","","(b)",""), label_y = 0.8, label_size = 18)
    } else {
      Final.plot <- cowplot::plot_grid(NULL, plot.1, NULL, plot.2, ncol = 1, nrow = 4,
                                       rel_heights = c(0.1, 1, 0.1, 1), 
                                       labels = c("","(a)","","(b)"), 
                                       label_y = 0.95, label_x = 0.25, label_size = 18)
    }
    
    return(Final.plot)
  } else if (package == "ggVennDiagram"){
    plot
  } else if (package == "eulerr"){
    class(df)="euler" # To avoid "Error in drawVennDiagram"
    plot(df, quantities = list(col=number_color, type=eulerr.label[[1]], cex=number_size), 
         edges = list(col=edge_color, lty=edge_lty, lwd=edge_size), 
         labels = list(cex=label_size, col=label_color), fills = col, alpha = fill_alpha,
         main = title1,...)
  }
}
# Function to output a data.frame with the mean abundances of each taxa per group of samples
intersect.df <- function(physeq, group, levels, rank = "None", relative = T, 
                         weighted = T, fraction = 0){
  
  if(relative){ physeq <- transform_sample_counts(physeq, function(x) x / sum(x)) }
  ps_melted <- psmelt(physeq)
  # Aggregate sample abundances (mean) per sample type (group/selected variable)
  if(relative){
    ps_agg <- stats::aggregate(as.formula(paste("Abundance ~ OTU +", group)), data = ps_melted,
                               function(x) (sum(x > 0)/length(x) >= fraction) * mean(x))
  } else {
    ps_agg <- stats::aggregate(as.formula(paste("Abundance ~ OTU +", group)), data = ps_melted,
                               function(x) (sum(x > 0)/length(x) >= fraction) * sum(x))
  }

  # Reshape the matrix to show sample type as columns and taxaIDs/OTUs as rows
  ps_mat <- reshape2::dcast(as.formula(paste("OTU ~ ", group)), data = ps_agg, 
                            value.var = "Abundance")
  
  rownames(ps_mat) <- ps_mat[, 1] # Add rownames
  ps_mat <- ps_mat[, -1] # Remove 1st column
  ps_mat <- ps_mat[rowSums(ps_mat) > 0, ] # Remove rows with a sum of 0
  ps_mat_bin <- (ps_mat>0)*1 # Get a binary matrix
  
  # data frame providing the common taxa (non-null) for the selected levels
  if (isTRUE(relative) & isTRUE(weighted)){
    main <- ps_mat[,levels, drop=FALSE]
    main <- apply(main, 2, function(x) x/sum(ps_mat)*100)
    digits <- 3
  } else if (isTRUE(weighted)){
    main <- ps_mat[,levels, drop=FALSE]
    digits <- 0
  } else {
    main <- ps_mat_bin[,levels, drop=FALSE]
    digits <- 0
  }
  main.clean = main[apply(main, 1, function(x) !any(x==0)), , drop=FALSE]
  extra = ps_mat[ , !(names(ps_mat) %in% levels), drop=FALSE]
  extra.clean = extra[apply(extra, 1, function(x) !all(x==0)), , drop=FALSE]
  common = main.clean[!(rownames(main.clean) %in% rownames(extra.clean)), , drop=FALSE]
  final <- cbind(common, Total = rowSums(common))
  tax=as.data.frame(tax_table(physeq))
  all=merge(round(final, digits), tax, by=0, all=FALSE)
  colnames(all)[1] <- "TaxaID"
  all %<>% dplyr::relocate("TaxaID", .after = last_col())
  all <- all[order(-all$Total), ]
  # Aglomerate taxa per rank
  if(rank != "None"){
    tax <- unique(all[,rank])
    all[is.na(all)] <- "N/A"
    for (i in tax) {
      if(is.na(i)){
        all <- all
      } else {
        x <- all[all[rank] == i, , drop=FALSE]
        x[1, c(levels, "Total")] = colSums(x[,c(levels, "Total")])
        x = x[1, , drop = FALSE]
        n.rank.col <- which(colnames(x)==rank) + 1
        x[, n.rank.col:length(colnames(x))] <- "NA"
        all[all[rank] == i,] <- x
        all <- all[!duplicated(all), ]
      }
    }
  }
  # Add % to colnames
  if (relative){
    colnames(all)[colnames(all) %in% c(levels, "Total")] <- paste0(c(levels, "Total"),"(%)")
  }
  return(all)
}
