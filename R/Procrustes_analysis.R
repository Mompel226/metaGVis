#' @title 3D plots
#'
#' @description
#' This function is used to create a MultiDimensional Scaling ordination plot in 3D with plotly.
#'
#' @param physeq A phyloseq object.
#' @param explanatory_variable A character string of the metadata variable that is going to be used to group the samples.
#' @param distance A character string defining the distance used to calculate the distance matrix:
#' It can be any distance accepted by the distance function of the phyloseq package and also "PhilR".
#' To use the Aitchison distance, first transform the otu_table of the phyloseq object with one of the CLR functions
#' (clr_min_percent, clr_0_plus_1 or clr_runif) and then use the "euclidean" distance.
#' Default: "bray"
#' @param method A character string defining the ordination method: "PCA", "PCoA", "RDA", "CAP", "CA", "CCA" or "NMDS"
#' Default: "PCoA"
#' @param colors A vector with the hexadecimal colors associated to the levels of the explanatory variable. (OPTIONAL)
#' @param constrained A logical value. If TRUE a constrained ordination is carried out.
#' Default: FALSE
#' @param formula The right hand side of the model formula -- the constraining variables and conditioning variables. e.g: "~ constraining variable"
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
Procruste_analysis <- function(list_of_distances){
  library(harrietr)
  library(tidyr)
  library(plotly)
  library(xtable)

  DistList2Sel <- list_of_distances
  #print(names(DistList2Sel[])) and check if "bray_loguq" is the 14th

  # Remove "bray_clr" and change name of bray_loguq"

  DistList2Sel$PCoA_bray_CLR_min_percent <- NULL
  DistList2Sel$PCoA_bray_CLR_runif <- NULL
  DistList2Sel$PCoA_bray_CLR_0_plus_1 <- NULL
  DistList2Sel <- DistList2Sel[!grepl("NMDS", names(DistList2Sel))]

  #To create procrustes PCoA/PCA
  ProCCol<-c("None" = "#000000", "TSS" = "#A65628", "Rarefy" = "#ebc83f", "UQ" = "#ff961f","CSS" = "#c203fc",
             "DESeq" = "#E41A1C", "DESeq-VS" = "#f593ca", "TMM"="#265880", "TMM-wsp" = "#34ebd8", "CLR" = "#4DAF4A")

  # Subset colors depending on the input data
  norms <- strsplit(names(dist.list), "_")
  norms2 <- unique(as.vector(sapply(norms, function(x){tail(x, n = 2)})))
  ProCCol <- subset(ProCCol, names(ProCCol) %in% norms2)

  ScreeList <- list()
  StressList <- list()
  PCoAList <- list()

  k<-1
  l<-1

  #Create dataframe to hold the correlations
  proCpairAllmethods<-setNames(data.frame(matrix(ncol=length(DistList2Sel), nrow=length(DistList2Sel))), c(names(DistList2Sel)))
  rownames(proCpairAllmethods)<-c(names(DistList2Sel))
  proCpairAllmethodsCor<-proCpairAllmethods
  proCpairAllmethodsSS<-proCpairAllmethods
  for (i in 1:length(DistList2Sel)) {
    for (j in 1:length(DistList2Sel)) {
      j <- DistList2Sel[[k]]
      i <- DistList2Sel[[l]]

      # Perform Distance-based Redundancy Analysis using capscale and then perform protest.
      # HOWEVER, since capscale has ~1 on the LHS of the formula, capscale is just performing PCoA.
      # Function protest tests the non-randomness (significance) between two configurations.
      prot<-protest(capscale(i~1), capscale(j~1))

      proCpairAllmethods[k,l]<-(1-prot$t0) #Can try with different measures correlations (prot$t0), sum of squares (prot$ss). PCoA/PCA title have to be changed accordingly below. Changed this to always being correlations, but other metrics could be provided.
      proCpairAllmethodsCor[k,l]<-(prot$t0) #Dist raw correlations
      proCpairAllmethodsSS[k,l]<-(prot$ss) #Have made a seperate object to store sum of squares
      k <- k+1
    }
    #Reasign 1 to k and iterate l
    print((l/length(DistList2Sel))*100)
    l <- l+1
    k <- 1
  }
  proCpairAllmethods<-as.dist(proCpairAllmethods) #I'm getting the right number of observations corresponding to the lower diagonal the as.dist function default is diag=FALSE, upper=FALSE, auto_convert_data_frames=TRUE.
  proCpairAllmethodsCor<-as.dist(proCpairAllmethodsCor)
  proCpairAllmethodsSS<-as.dist(proCpairAllmethodsSS)

  ### Create PCoA with Distance-based redundancy analysis of pairwise procrustes correlations and sum of squares.

  PCoAcsObject<-capscale(proCpairAllmethods~1)
  PCoAcsObjectSS<-capscale(proCpairAllmethodsSS~1)

  #Make stressplot
  #Extract ordination distances and merge with observed dissimilarity
  #Correlations
  stress<-stressplot(PCoAcsObject)
  df <- reshape2::melt(as.matrix(stress))
  names(df)<-c("rowOrd", "colOrd", "OrdDist")
  df<-filter(df, OrdDist>0)
  df2 <- reshape2::melt(as.matrix(proCpairAllmethods))
  names(df2)<-c("rowObs", "colObs", "ObsDism")
  df2<-filter(df2, ObsDism>0)
  df<-unite(df, mergecol, c(rowOrd, colOrd), remove=FALSE)
  df2<-unite(df2, mergecol, c(rowObs, colObs), remove=FALSE)
  ggstress<-merge(df, df2, by="mergecol")
  #Create plot name correlations
  pltName <- paste( 'stress', 'procrustesCor', sep = '' )
  #Create Shepard diagram
  StressList[[ pltName ]]<-ggplot(ggstress) +
    geom_point(aes(ObsDism, OrdDist), size=1) + #Can change to 0.01 when running all
    ggtitle(paste("Shepard diagram", sep=" ")) +
    labs(x = "Observed dissimilarity", y = "Ordination distance") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12))
  #Sum of squares
  stressSS<-stressplot(PCoAcsObjectSS)
  dfSS <- reshape2::melt(as.matrix(stressSS))
  names(dfSS)<-c("rowOrd", "colOrd", "OrdDist")
  dfSS<-filter(dfSS, OrdDist>0)
  df2SS <- reshape2::melt(as.matrix(proCpairAllmethodsSS))
  names(df2SS)<-c("rowObs", "colObs", "ObsDism")
  df2SS<-filter(df2SS, ObsDism>0)
  dfSS<-unite(dfSS, mergecol, c(rowOrd, colOrd), remove=FALSE)
  df2SS<-unite(df2SS, mergecol, c(rowObs, colObs), remove=FALSE)
  ggstressSS<-merge(dfSS, df2SS, by="mergecol")
  #Create plot name sum of squares
  pltName <- paste( 'stress', 'procrustesSS', sep = '' )
  #Create stressplot sum of squares
  StressList[[ pltName ]]<-ggplot(ggstressSS) +
    geom_point(aes(ObsDism, OrdDist), size=1) + #Can change to 0.01 when running all
    ggtitle(paste("Shepard diagram", sep=" ")) +
    labs(x = "Observed dissimilarity", y = "Ordination distance") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12))
  ##PCoA correlations
  ##Add eig to plot axes. with cmdscale there are negative values not with capscale
  eig <- PCoAcsObject$CA$eig
  # Calculate the variation explained by PCoA1, 2, 3 and 4
  # and use it to generate axis labels
  eig_1_2 <- eig[1:4] / sum(eig) * 100
  eig_1 <- paste("PCoA1", round(eig_1_2[1], digits = 2), "% variance")
  eig_2 <- paste("PCoA2", round(eig_1_2[2], digits = 2), "% variance")
  eig_3 <- paste("PCoA3", round(eig_1_2[3], digits = 2), "% variance")
  eig_4 <- paste("PCoA4", round(eig_1_2[4], digits = 2), "% variance")
  ##Pull out coordinates for plotting from the ca object
  #Structuring to add to Metadata2
  PCoACA<-PCoAcsObject$CA #The ca object contains the actual ordination results: u ((Weighted) orthonormal site        scores), v ((Weighted)       orthonormal species scores) all na in mine, Xbar (The standardized data matrix after previous stages of analysis), and imaginary.u.eig ???. Info http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/cca.object.html
  PCoA<-as.data.frame(PCoACA$u)
  #Change colnames. Now add dis and trans info to names
  if (ncol(PCoA) < 4){PCoA$MDS4 <- rep(0, nrow(PCoA))}
  colnames(PCoA) <- c("MDS1","MDS2", "MDS3", "MDS4")
  #Add row names to df
  PCoA$Sample <- row.names(PCoA)

  # Subset shapes depending on the input data
  norms2 <- unique(sapply(norms, "[", -1))

  is.integer0 <- function(x)
  {is.integer(x) && length(x) == 0L}

  x <- vector()
  if(any(!is.integer0(grep("bray*", names(dist.list))), x <- append(x, c("Bray-Curtis")), ""))
  ifelse(!is.integer0(grep("^wunifrac*", names(dist.list))), x <- append(x, c("UniFrac (weighted)")),
         ifelse(!is.integer0(grep("^unifrac*", names(dist.list))), x <- append(x, c("UniFrac (non-weighted)")), ""))
  ifelse(!is.integer0(grep("*euclidean_CLR_0plus1*", names(dist.list))) & !is.integer0(grep("*euclidean_logChord_0plus1*", names(dist.list))),
         x <- append(x, c("Aitchison or log-Chord (0+1)")),
         ifelse(!is.integer0(grep("*euclidean_CLR_0plus1*", names(dist.list))), x <- append(x, c("Aitchison (0+1)")),
                ifelse(!is.integer0(grep("*euclidean_logChord_0plus1*", names(dist.list))), x <- append(x, c("log-Chord (0+1)")), "")))
  ifelse(!is.integer0(grep("*euclidean_CLR_min*", names(dist.list))) & !is.integer0(grep("*euclidean_logChord_min*", names(dist.list))),
         x <- append(x, c("Aitchison or log-Chord (min/2)")),
         ifelse(!is.integer0(grep("*euclidean_CLR_min*", names(dist.list))), x <- append(x, c("Aitchison (min/2)")),
                ifelse(!is.integer0(grep("*euclidean_logChord_min*", names(dist.list))), x <- append(x, c("log-Chord (min/2)")), "")))
  if (any(!is.integer0(grep("*euclidean_Hellinger*", names(dist.list))))){ x <- append(x, c("Hellinger")) }
  if (any(!is.integer0(grep("*euclidean_logChord_1p*", names(dist.list))))){ x <- append(x, c("log-Chord (All+1)")) }
  if (any(!is.integer0(grep("*euclidean_CLR_runif*", names(dist.list))))){ x <- append(x, c("Aitchison (runif)")) }
  ProCShape <- c(1:length(x))
  names(ProCShape) <- x
  if (any(grep("Aitchison or log-Chord \\(0\\+1\\)", names(ProCShape)))){
    value1 <- "Aitchison or log-Chord (0+1)"
    value2 <- "Aitchison or log-Chord (0+1)"
  } else {
    value1 <- "log-Chord (0+1)"
    value2 <- "Aitchison (0+1)"
  }
  if (any(grep("Aitchison or log-Chord \\(min\\/2\\)", names(ProCShape)))){
    value3 <- "Aitchison or log-Chord (min/2)"
    value4 <- "Aitchison or log-Chord (min/2)"
  } else {
    value3 <- "log-Chord (min/2)"
    value4 <- "Aitchison (min/2)"
  }
  #Create metadata
  MetadataProC<-data.frame(Sample=rownames(PCoA))
  MetadataProC$Trans <-
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*None", MetadataProC$Sample), "None",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*TSS", MetadataProC$Sample), "TSS",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*Rarefy", MetadataProC$Sample), "Rarefy",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*UQ", MetadataProC$Sample), "UQ",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*CSS", MetadataProC$Sample), "CSS",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*DESeq_VS", MetadataProC$Sample), "DESeq-VS",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*DESeq", MetadataProC$Sample), "DESeq",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*TMM_wsp", MetadataProC$Sample), "TMM-wsp",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*MM", MetadataProC$Sample), "TMM",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*CLR", MetadataProC$Sample), "CLR", NA))))))))))
  MetadataProC$Dist <-
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*euclidean_Hellinger*", MetadataProC$Sample), "Hellinger",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*euclidean_logChord_0plus1*", MetadataProC$Sample), value1,
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*euclidean_logChord_min*", MetadataProC$Sample), value3,
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*euclidean_CLR_min*", MetadataProC$Sample), value4,
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*euclidean_CLR_0plus1*", MetadataProC$Sample), value2,
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*euclidean_CLR_runif*", MetadataProC$Sample), "Aitchison (runif)",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*euclidean_logChord_1p*", MetadataProC$Sample), "log-Chord (All+1)",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*wunifrac*", MetadataProC$Sample), "UniFrac (weighted)",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*unifrac*", MetadataProC$Sample), "UniFrac (non-weighted)",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*bray*", MetadataProC$Sample), "Bray-Curtis", NA))))))))))

  #Merge according to Sample"Euclidean (Log runif)"
  MetadataProC2<-merge(MetadataProC, PCoA, by="Sample")
  #Create plot name
  pltName <- paste( 'PCoA', 'procrustesCor', sep = '' )
  #create PCoA
  PCoAList[[ pltName ]]<-ggplot(MetadataProC2) +
    #geom_line(aes(x=MDS1, y=MDS2, group=Trans), size=0.1, linetype="dotted") +
    geom_jitter(aes(MDS1, MDS2, col=Trans, shape=Dist), width=0.004, height=0.004, alpha=0.8, size=3, stroke=1.5) +
    #geom_point(aes(MDS1, MDS2, color = Sample_LPSX, group = Sample, shape = Temperature), size=5) +
    scale_color_manual(values=ProCCol) +
    scale_shape_manual(values=ProCShape) +
    ggtitle(paste("PCoA ", 'procrustes rotation 1-correlations')) +
    labs(colour="Transformation", shape="Distance", x = eig_1, y = eig_2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12), legend.position="bottom") +
    guides(colour=guide_legend(title.vjust = 0.8, nrow = 3), shape=guide_legend(title.vjust = 0.8, nrow = 3))

  ##PCoA sum of squares overwriting the objects created above for correlations using the alreade created metadata
  ##Add eig to plot axes. with cmdscale there are negative values not with capscale
  eig <- PCoAcsObjectSS$CA$eig
  # Calculate the variation explained by PCoA1, 2, 3 and 4
  # and use it to generate axis labels
  eig_1_2 <- eig[1:4] / sum(eig) * 100
  eig_1 <- paste("PCoA1", round(eig_1_2[1], digits = 2), "% variance")
  eig_2 <- paste("PCoA2", round(eig_1_2[2], digits = 2), "% variance")
  eig_3 <- paste("PCoA3", round(eig_1_2[3], digits = 2), "% variance")
  eig_4 <- paste("PCoA4", round(eig_1_2[4], digits = 2), "% variance")
  ##Pull out coordinates for plotting from the ca object
  #Structuring to add to Metadata2
  PCoACA<-PCoAcsObjectSS$CA #The ca object contains the actual ordination results: u ((Weighted) orthonormal site        scores), v ((Weighted)       orthonormal species scores) all na in mine, Xbar (The standardized data matrix after previous stages of analysis), and imaginary.u.eig ???. Info http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/cca.object.html
  PCoA<-as.data.frame(PCoACA$u)
  #Change colnames. Now add dis and trans info to names
  if (ncol(PCoA) < 4){PCoA$MDS4 <- rep(0, nrow(PCoA))}
  colnames(PCoA) <- c("MDS1","MDS2", "MDS3","MDS4")
  #Add row names to df
  PCoA$Sample <- row.names(PCoA)

  #Merge according to Sample
  MetadataProC2<-merge(MetadataProC, PCoA, by="Sample")
  #Create plot name
  pltName <- paste( 'PCoA', 'procrustesSS', sep = '' )
  #create PCoA
  PCoAList[[ pltName ]]<-ggplot(MetadataProC2) +
    geom_jitter(aes(MDS1, MDS2, col=Trans, shape=Dist), width=0.00, height=0.00, alpha=0.8, size=3, stroke=1.5) +
    scale_color_manual(values=ProCCol) +
    scale_shape_manual(values=ProCShape) +
    ggtitle(paste("PCoA ", 'procrustes rotation sum of squares')) +
    labs(colour="Transformation", shape="Distance", x = eig_1, y = eig_2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12), legend.position="bottom") +
    guides(colour=guide_legend(title.vjust = 0.8, nrow = 3), shape=guide_legend(title.vjust = 0.8, nrow = 3))

  #Screeplot correlations
  screeplot<-data.frame(PCoAcsObject$CA$eig)
  colnames(screeplot)<-c("eig")
  screeplot$eig <- screeplot$eig[1:length(screeplot$eig)] / sum(screeplot$eig) * 100
  screeplot<-add_rownames(screeplot, "MDS")
  screeplot$MDS <- factor(screeplot$MDS, levels=c(sprintf("MDS%d", 1:length(screeplot$eig))))
  #Create plot name
  pltName <- paste( 'scree', 'procrustesCor', sep = '' )
  #create screeplot
  ScreeList[[ pltName ]]<-ggplot(screeplot, aes(x=MDS, y=eig)) +
    geom_bar(stat="identity") +
    labs(x ="Principal coordinates", y ="Eigenvalues (%)") +
    ggtitle(paste("Scree plot ")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  #Screeplot sum of squares
  screeplot<-data.frame(PCoAcsObjectSS$CA$eig)
  colnames(screeplot)<-c("eig")
  screeplot$eig <- screeplot$eig[1:length(screeplot$eig)] / sum(screeplot$eig) * 100
  screeplot<-add_rownames(screeplot, "MDS")
  screeplot$MDS <- factor(screeplot$MDS, levels=c(sprintf("MDS%d", 1:length(screeplot$eig))))
  #Create plot name
  pltName <- paste( 'scree', 'procrustesSS', sep = '' )
  #create screeplot
  ScreeList[[ pltName ]]<-ggplot(screeplot, aes(x=MDS, y=eig)) +
    geom_bar(stat="identity") +
    labs(x ="Principal coordinates", y ="Eigenvalues (%)") +
    ggtitle(paste("Scree plot ")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12), axis.text.x=element_blank(), axis.ticks.x=element_blank())

  Correls<-melt_dist(matrix(proCpairAllmethodsCor))

  PCoAList$densityplotCor<-ggplot(Correls, aes(x=dist)) +
    geom_histogram(aes(y=..density..), binwidth=0.01, colour="black", fill="white") + #..density.. is to have it scale with geom_density
    geom_density(alpha=.25, fill="#FF6666") +
    labs(x ="Correlation", y ="Density") +
    ggtitle(paste("Density plot ")) +
    xlim(0,1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12))
  Sumssqs<-melt_dist(matrix(proCpairAllmethodsSS))
  PCoAList$densityplotSS<-ggplot(Sumssqs, aes(x=dist)) +
    geom_histogram(aes(y=..density..), binwidth=0.01, colour="black", fill="white") + #..density.. is to have it scale with geom_density
    geom_density(alpha=.25, fill="#FF6666") +
    labs(x ="Sum of squares", y ="Density") +
    ggtitle(paste("Density plot ")) +
    xlim(0,1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12))

  lay <- rbind(c(1,1),
               c(1,1),
               c(2,3))

  ## Create figure with correlations
  dir.create("Ordination_plots/")
  Top <- ggarrange(PCoAList$PCoAprocrustesCor)
  Bottom <- ggarrange(PCoAList$densityplotCor, StressList$stressprocrustesCor, ScreeList$screeprocrustesCor, ncol = 3)
  Final <- ggarrange(Top, Bottom, nrow = 2)
  ggsave(filename = paste("Ordination_plots/Procrustes_Correlations", ".png", sep=""),
         plot = Final, height = 10, width = 12.5, unit = "in", dpi = 300)

  ## Create figure with sum of squares
  Top <- ggarrange(PCoAList$PCoAprocrustesSS)
  Bottom <- ggarrange(PCoAList$densityplotSS, StressList$stressprocrustesSS, ScreeList$screeprocrustesSS, ncol = 3)
  Final <- ggarrange(Top, Bottom, nrow = 2)
  ggsave(filename = paste("Ordination_plots/Procrustes_SumOfSquares", ".png", sep=""),
         plot = Final, height = 10, width = 12.5, unit = "in", dpi = 300)
}
