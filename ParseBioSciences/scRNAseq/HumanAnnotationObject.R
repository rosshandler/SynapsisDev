library(scran)
library(Matrix)
library(ggplot2)

path2data <- "/data2/hanna/synaptogenesis/newvolume/analysis/Jan2024/25-3-24/"

sce  <- readRDS(paste0(path2data, "sce.delta_h_t_26-3-24.rds"))
sce  <- logNormCounts(sce, exprs_values = "decontXcounts", name = "decontXlogcounts")
gexp <- logcounts(sce)
gexp_decont <- logcounts(sce, assay.type = "decontXcounts")

rownames(gexp) <- rowData(sce)$gene_name
rownames(gexp_decont) <- rowData(sce)$gene_name

df_plot <- data.frame(colData(sce))

day_colours <- rev(wesanderson::wes_palette("Zissou1", 8, type = "continuous"))
days <- sort(as.numeric(unique(df_plot$day_H)))

level_order <- days
df_plot$day <- factor(df_plot$day_H,levels=level_order)

umap <- reducedDim(sce)
colnames(umap) <- c("UMAP1", "UMAP2")

df_plot <- cbind(df_plot, umap)

leiden_colours <- c(
"#655f91",
"#4bba30",
"#ef49c3",
"#cdbe00",
"#9b85ff",
"#da8d00",
"#3949ab",
"#efbf61",
"#b1007a",
"#73d8bb",
"#df2831",
"#028ac3",
"#a3190c",
"#7eb6ff",
"#9a6700",
"#ff8ae4",
"#4f6600",
"#ff789b",
"#8c3832",
"#ffaa8b"
)

plotLayoutExpression <- function(gene="CLIC6", layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
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

plotLayoutExpressionDecont <- function(gene="CLIC6", layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp_decont[gene,]))
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
      if (layout=="FA"){ 
        df_tmp      <- data.frame(cbind(df_plot, mt.fraction))
        plot.index  <- order(df_tmp$mt.fraction)
        ggplot(df_tmp[plot.index,], aes(x = Fa1, y = Fa2, colour = mt.fraction)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="red") +
          labs(color = "MT fraction") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMAP"){
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
    }else{
    message("Layout name not found")
    }
}

plotLayoutContamination <- function(layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    contamination <- df_plot$decontX_contamination
      if (layout=="FA"){ 
        df_tmp      <- data.frame(cbind(df_plot, contamination))
        plot.index  <- order(df_tmp$contamination)
        ggplot(df_tmp[plot.index,], aes(x = Fa1, y = Fa2, colour = contamination)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="red") +
          labs(color = "Contamination") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMAP"){
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
    }else{
    message("Layout name not found")
    }
}

plotLayoutDoubletScore <- function(layout="UMAP"){
  require(Matrix)
  require(ggplot2)
    doublet.score <- df_plot$scDblFinder.score
      if (layout=="FA"){ 
        df_tmp      <- data.frame(cbind(df_plot, doublet.score))
        plot.index  <- order(df_tmp$doublet.score)
        ggplot(df_tmp[plot.index,], aes(x = Fa1, y = Fa2, colour = doublet.score)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="black") +
          labs(color = "Doublet Score") +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else if(layout=="UMAP"){
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
    }else{
    message("Layout name not found")
    }
}

plotLayoutLeiden <- function(layout="UMAP"){
  require(ggplot2)
  if (layout=="UMAP"){ 
    ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leiden_clusters))) +
      geom_point(size = 1) +        
      scale_color_manual(values=leiden_colours, name = "leiden") +
      theme_minimal() + 
      labs(col="Leiden") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="FA"){
    ggplot(df_plot, aes(x = Fa1, y = Fa2, col = factor(leiden))) +
      geom_point(size = 1) +
      scale_color_manual(values=leiden_colours, name = "leiden_clusters") +
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
  }else if(layout =="FA"){
    ggplot(df_plot[plot.index,], aes(x = Fa1, y = Fa2, col = factor(day))) +
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
  if (layout=="FA"){
    ggplot(df_plot[plot.index,], aes(x = Fa1, y = Fa2, col = batch_H)) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Batch") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = batch_H)) +
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
        leiden <- df_plot$leiden
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
        leiden <- df_plot$leiden
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
