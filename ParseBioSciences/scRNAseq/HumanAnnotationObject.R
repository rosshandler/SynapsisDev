library(scran)
library(Matrix)
library(ggplot2)

path2data <- "/data2/hanna/axonoutgrowth/analysis/DimRed/"

sce  <- readRDS(paste0(path2data, "sce_decontXDimRed.rds"))

sce  <- logNormCounts(sce, assay.type = "decontXcounts", name = "decontXlogcounts")
gexp_decont <- logcounts(sce)
sce  <- logNormCounts(sce, assay.type = "counts", name = "logcounts")
gexp <- logcounts(sce)

rownames(gexp) <- rowData(sce)$gene_name
rownames(gexp_decont) <- rowData(sce)$gene_name

#sce_filt <- sce_a[calculateAverage(sce_a)>0.1,]
#sce_filt <- logNormCounts(sce_filt,assay.type = "decontXcounts")

#decomp  <- modelGeneVar(sce_filt)
#hvgs    <- rownames(decomp)[decomp$FDR < 0.1]
#pca      <- irlba::prcomp_irlba(t(logcounts(gexp_decont[hvgs,])), n = 30)
#graph <- buildSNNGraph(pca$x, d = NA, transposed = TRUE)
#set.seed(42)
#clusters <- leiden(graph, resolution_parameter = 2)
#names(clusters) <-  colData(sce_a)$cell
#colData(sce_a)$leidenClustersdecontX <- clusters

df_plot <- data.frame(colData(sce))

day_colours <- rev(wesanderson::wes_palette("Zissou1", 13, type = "continuous"))
days <- sort(as.numeric(unique(df_plot$day_H)))

level_order <- days
df_plot$day <- factor(df_plot$day_H,levels=level_order)

umap <- reducedDim(sce,"UMAP")
colnames(umap) <- c("UMAP1", "UMAP2")
umap_decontX <- reducedDim(sce,"UMAP_decontXlogcounts")
colnames(umap_decontX) <- c("UMAP1_decontX", "UMAP2_decontX")

df_plot <- cbind(df_plot, umap, umap_decontX)

leiden_colours <- c(
"#7166d9",
"#58c655",
"#af57c6",
"#acbb37",
"#d149ac",
"#5ea136",
"#d9407f",
"#6fc480",
"#d83b52",
"#4ebfa9",
"#bf3c24",
"#46aed7",
"#e0752d",
"#7293dd",
"#d29e37",
"#5d64ac",
"#838c26",
"#c88ed9",
"#3c9153",
"#964d88",
"#53772f",
"#de80ad",
"#337d5c",
"#e06953",
"#b3b269",
"#a04456",
"#7b702d",
"#db7b7d",
"#9b5e2d",
"#db9869"
)

leiden_colours_decontX <- c(
"#dd862f",
"#915bc6",
"#5ec456",
"#d2489a",
"#9ab635",
"#637dc9",
"#d3b246",
"#ca88ca",
"#45953f",
"#d14558",
"#58c39b",
"#ce4f2e",
"#48adcf",
"#a78832",
"#ac5676",
"#5a792b",
"#e39175",
"#37835c",
"#9c5d31",
"#9fb26a",
"#736c2c"
)

plotLayoutExpression <- function(gene="CLIC6", layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
      if (layout=="UMAP"){ 
        df_tmp      <- data.frame(cbind(df_plot, logcounts))
        plot.index  <- order(df_tmp$logcounts)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1, y = UMAP2, colour = logcounts)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="darkgreen") +
          labs(color = paste0(gene,"\nlog(counts)")) +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMAP_decontX"){
        df_tmp      <- data.frame(cbind(df_plot, logcounts))
        plot.index  <- order(df_tmp$logcounts)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1_decontX, y = UMAP2_decontX, colour = logcounts)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="darkgreen") +
          labs(color = paste0(gene,"\nlog(counts)")) +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")    
    }else{
    message("Layout name not found")
    }
 }else{
    message(gene," was not detected in the expression matrix")
 }
}

plotLayoutExpressionDecont <- function(gene="CLIC6", layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp_decont[gene,]))
    if (sum(logcounts)>0){
      if (layout=="UMAP"){ 
        df_tmp      <- data.frame(cbind(df_plot, logcounts))
        plot.index  <- order(df_tmp$logcounts)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1, y = UMAP2, colour = logcounts)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="darkgreen") +
          labs(color = paste0(gene,"\nlog(counts)")) +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMAP_decontX"){
        df_tmp      <- data.frame(cbind(df_plot, logcounts))
        plot.index  <- order(df_tmp$logcounts)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1_decontX, y = UMAP2_decontX, colour = logcounts)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="darkgreen") +
          labs(color = paste0(gene,"\nlog(counts)")) +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")    
    }else{
    message("Layout name not found")
    }
 }else{
    message(gene," was not detected in the expression matrix")
 }
}

# Not activated yet
plotLayoutDenoisedExpression <- function(gene="CLIC6", layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    logcounts <- as.vector(as.matrix(dgexp[gene,]))
    if (sum(logcounts)>0){
      if (layout=="FA"){ 
        df_tmp      <- data.frame(cbind(df_plot, logcounts))
        plot.index  <- order(df_tmp$logcounts)
        ggplot(df_tmp[plot.index,], aes(x = Fa1, y = Fa2, colour = logcounts)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="darkgreen") +
          labs(color = paste0(gene,"\nlog(counts)")) +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMAP"){
        df_tmp      <- data.frame(cbind(df_plot, logcounts))
        plot.index  <- order(df_tmp$logcounts)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1, y = UMAP2, colour = logcounts)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="darkgreen") +
          labs(color = paste0(gene,"\nlog(counts)")) +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")    
    }else{
    message("Layout name not found")
    }
 }else{
    message(gene," was not detected in the expression matrix")
 }
}

plotLayoutMTfraction <- function(layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    mt.fraction <- df_plot$mt.fraction
      if (layout=="UMAP"){ 
        df_tmp      <- data.frame(cbind(df_plot, mt.fraction))
        plot.index  <- order(df_tmp$mt.fraction)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1, y = UMAP2, colour = mt.fraction)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="red") +
          labs(color = "MT fraction") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMA_decontXP"){
        df_tmp      <- data.frame(cbind(df_plot, mt.fraction))
        plot.index  <- order(df_tmp$mt.fraction)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1_decontX, y = UMAP2_decontX, colour = mt.fraction)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="red") +
          labs(color = "MT fraction") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")    
    }else{
    message("Layout name not found")
    }
}

plotLayoutContamination <- function(layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    contamination <- df_plot$decontX_contamination
      if (layout=="UMAP"){ 
        df_tmp      <- data.frame(cbind(df_plot, contamination))
        plot.index  <- order(df_tmp$contamination)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1, y = UMAP2, colour = contamination)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="red") +
          labs(color = "Contamination") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMAP_decontX"){
        df_tmp      <- data.frame(cbind(df_plot, contamination))
        plot.index  <- order(df_tmp$contamination)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1_decontX, y = UMAP2_decontX, colour = contamination)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="red") +
          labs(color = "Contamination") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")    
    }else{
    message("Layout name not found")
    }
}

plotLayoutDoubletScore <- function(layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    doublet.score <- df_plot$scDblFinder.score
      if (layout=="UMAP"){ 
        df_tmp      <- data.frame(cbind(df_plot, doublet.score))
        plot.index  <- order(df_tmp$doublet.score)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1, y = UMAP2, colour = doublet.score)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="black") +
          labs(color = "Doublet Score") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMAP_decontX"){
        df_tmp      <- data.frame(cbind(df_plot, doublet.score))
        plot.index  <- order(df_tmp$doublet.score)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1_decontX, y = UMAP2_decontX, colour = doublet.score)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="black") +
          labs(color = "Doublet Score") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")    
    }else{
    message("Layout name not found")
    }
}

plotLayoutLeiden <- function(layout="UMAP"){
  require(ggplot2)
  if (layout=="UMAP"){ 
    ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leidenclusters_before_decontX))) +
      geom_point(size = 1) +        
      scale_color_manual(values=leiden_colours, name = "leiden") +
      theme_minimal() + 
      labs(col="Leiden") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP_decontX"){
    ggplot(df_plot, aes(x = UMAP1_decontX, y = UMAP2_decontX, col = factor(leidenclusters_after_decontX))) +
      geom_point(size = 1) +
      scale_color_manual(values=leiden_colours_decontX, name = "leiden_decontX") +
      theme_minimal() + 
      labs(col="Leiden") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
  }else{
    message("Layout name not found")
  }
}

plotLayoutDay <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="UMAP"){ 
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(day))) +
      geom_point(size = 1) +        
      scale_color_manual(values=day_colours, name = "Day") +
      theme_minimal() + 
      labs(col="Cell type") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP_decontX"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1_decontX, y = UMAP2_decontX, col = factor(day))) +
      geom_point(size = 1) +
      scale_color_manual(values=day_colours, name = "Day") +
      theme_minimal() + 
      labs(col="Cell type") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
  }else{
    message("Layout name not found")
  }
}

plotLayoutBatch <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = batch_H)) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Batch") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP_decontX"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1_decontX, y = UMAP2_decontX, col = batch_H)) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Batch") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else{
    message("Layout name not found")
  }
}

plotViolinExpressionLeidenDecont <- function(gene="CLIC6"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp_decont[gene,]))
    if (sum(logcounts)>0){
        leiden <- df_plot$leidenclusters_after_decontX
        ggplot(mapping =  aes(x = leiden, 
                              y = logcounts, 
                              fill = factor(leiden))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=leiden_colours, name = "Leiden") +
        labs(y = "Log2 normalised count", x = "Leiden clusters") + 
        ggtitle(paste0(gene," expression across leiden clusters")) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_blank(),
              axis.title.x = element_blank()) 
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

plotViolinExpressionDayDecont <- function(gene="CLIC6"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp_decont[gene,]))
    if (sum(logcounts)>0){
        day <- df_plot$day
        ggplot(mapping =  aes(x = day, 
                              y = logcounts, 
                              fill = factor(day))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=day_colours, name = "Day") +
        labs(y = "Log2 normalised count") + 
        ggtitle(paste0(gene," expression across days")) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
              legend.position = "none",
              axis.title.x = element_blank())
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

plotViolinExpressionLeiden <- function(gene="CLIC6"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        leiden <- df_plot$leidenclusters_before_decontX
        ggplot(mapping =  aes(x = leiden, 
                              y = logcounts, 
                              fill = factor(leiden))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=leiden_colours, name = "Leiden") +
        labs(y = "Log2 normalised count", x = "Leiden clusters") + 
        ggtitle(paste0(gene," expression across leiden clusters")) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_blank(),
              axis.title.x = element_blank()) 
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

plotViolinExpressionDay <- function(gene="CLIC6"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        day <- df_plot$day
        ggplot(mapping =  aes(x = day, 
                              y = logcounts, 
                              fill = factor(day))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=day_colours, name = "Day") +
        labs(y = "Log2 normalised count") + 
        ggtitle(paste0(gene," expression across days")) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
              legend.position = "none",
              axis.title.x = element_blank())
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

save.image(file=paste0(path2data, 'plots_feline_annotationObjectMay2024.RData'))
