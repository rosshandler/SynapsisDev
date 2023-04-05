library(scuttle)
library(scran)
library(irlba)
library(Rtsne)
library(Matrix)
library(ggplot2)
library(biomaRt)
library(viridisLite)
library(viridis)
library(scDblFinder)

library(umap)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

umap = import('umap')

path2data   <- '/data2/ivanir/Feline2023/ParseBS/newvolume/analysis/sCell/all-well/DGE_unfiltered/'
sample_info <- read.table('/data2/ivanir/Feline2023/ParseBS/newvolume/analysis/sample_info.tab',
  sep = "\t", header = TRUE)

#read the spare matrix into counts
#read the geneIDs, names and genome-of-origine into genes 
counts    <- t(readMM(paste0(path2data, "/DGE.mtx")))
genes     <- read.csv(paste0(path2data, "/all_genes.csv"))
metadata  <- read.csv(paste0(path2data, "/cell_metadata.csv"))

lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

genes_human <- genes[genes$genome == "hg38",]
genes_mouse <- genes[genes$genome == "mm10",]

dim(counts)
#[1] 119684 2660423
dim(counts[,ngenes > 400 & lib.sizes > 500])
#[1] 119684  70566

counts   <- counts[,ngenes > 400 & lib.sizes > 500]
metadata <- metadata[ngenes > 400 & lib.sizes > 500,]
lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

pdf("hist_ngenes_libsize_ratio.pdf")
hist(ngenes/lib.sizes)
dev.off()

counts   <- counts[,ngenes/lib.sizes < 0.9]
metadata <- metadata[ngenes/lib.sizes < 0.9,]
lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

#create a vector with nrow(metadata) many NAs for smaple_bc1_well, smaple_nuber and smaple_name_human/mouse 
sample_bc1_well <- rep(NA, nrow(metadata))        
sample_number   <- rep(NA, nrow(metadata))
sample_name_human     <- rep(NA, nrow(metadata))
sample_name_mouse     <- rep(NA, nrow(metadata))

#changing the sample_info.tab for Feline's data 
sample_info$H_day <- sample_info$H_Timepoint 
sample_info$H_day <- gsub("55\\+","",sample_info$H_day)
sample_info$H_day <- as.integer(sample_info$H_day)
sample_info$H_day <-  sample_info$H_day +55
sample_info$Sample_name_H <- paste(sample_info$H_Batch, sample_info$H_day, sample_info$H_Replicate, sep="_")

sample_info$M_day <- sample_info$M_Timepoint
sample_info$M_day <- gsub("8\\+","",sample_info$M_day)
sample_info$M_day <- as.integer(sample_info$M_day)
sample_info$M_day <-  sample_info$M_day +8
sample_info$Sample_name_M <- paste(sample_info$M_Batch, sample_info$M_day, sample_info$M_Replicate, sep="_")

#write.table(sample_info, file = "/data2/ivanir/Feline2023/ParseBS/newvolume/analysis/sample_info_alt.tab"
            , sep = "\t", row.names=FALSE, quote=FALSE)
#write.csv(metadata, "/data2/hanna/synaptogenesis/newvolume/analysis/metadata_alt.csv", row.names=FALSE, quote = FALSE)

#creating a vector with information of which well the cell is in
samples <- unique(sample_info$Sample_well)
for (i in 1:length(samples)){
  sample_bc1_well[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))] <- sample_info$Sample_well[i]
  sample_number[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]   <- sample_info$Sample_Number[i]
  sample_name_human[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]     <- sample_info$Sample_name_H[i]
  sample_name_mouse[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]     <- sample_info$Sample_name_M[i]
}
sample_name_human <- gsub(" ","_",sample_name_human)
sample_name_mouse <- gsub(" ","_",sample_name_mouse)


#creating submeta columns for human and mouse samples 
submeta_human <- data.frame(rlist::list.rbind(strsplit(sample_name_human, split="_")))
colnames(submeta_human) <- c("batch", "day", "replicate")
submeta_human$day <- gsub("d","",submeta_human$day)

submeta_mouse <- data.frame(rlist::list.rbind(strsplit(sample_name_mouse, split="_")))
colnames(submeta_mouse) <- c("batch", "day", "replicate")
submeta_mouse$day <- gsub("d","",submeta_mouse$day)

#adding new columns to metadata
metadata <- data.frame(cbind(metadata, lib.sizes, sample_number, sample_bc1_well, sample_name_human, submeta_human, sample_name_mouse, submeta_mouse))

plot_df <- metadata

#setwd('/data2/ivanir/Feline2023/ParseBS/newvolume/analysis/sCell/QC')
setwd('/data2/hanna/synaptogenesis/newvolume/analysis/QC')

#visualising metadata for human sample 
ggplot(plot_df, aes (x = factor(sample_name_human), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("UMIsBySample_beforeQC_H.pdf")

#visualising metadata for mouse sample 
ggplot(plot_df, aes (x = factor(sample_name_mouse), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("UMIsBySample_beforeQC_M.pdf")
 
pdf("cell_complexity.pdf")
qplot(lib.sizes, ngenes, col = ifelse(ngenes > 400 & lib.sizes > 500 , "drop", "keep")) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  labs(x = "UMI count", y = "Number of expressed genes") +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()
print(ggp)


#connecting to BioMart datasets
ensembl_human <- useEnsembl(biomart = "ensembl",  dataset = "hsapiens_gene_ensembl",mirror = "useast")
ensembl_mouse <- useEnsembl(biomart = "ensembl",  dataset = "mmusculus_gene_ensembl",mirror = "useast")

#retriving user specified atributes form the BioMart dataset of interest
gene_map_H  <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  filters = "hgnc_symbol", values = genes$gene_name, mart = ensembl_human)
gene_map_M  <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  filters = "hgnc_symbol", values = genes$gene_name, mart = ensembl_mouse)

#retiving mitochodnral genes to find the mitochondrial fraction  
mt.index_H    <- gene_map_H$chromosome_name == "MT"
mt.index_M    <- gene_map_M$chromosome_name == "MT"
mt.counts   <- counts[which(genes$gene_name %in% gene_map_H$hgnc_symbol[mt.index_H]), ]
dim(mt.counts)
mt.counts   <- rbind(mt.counts,counts[which(genes$gene_name %in% gene_map_M$hgnc_symbol[mt.index_M]), ])
dim(mt.counts)
mt.fraction <- colSums(mt.counts/lib.sizes)

mt.p  <- pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
mt.lim<- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.05)])

mt.lim
#[1] 0.04396735

mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.001)])

#Threhdold
mt.lim
#[1] 0.06210699

metadata <- data.frame(cbind(metadata,mt.fraction))

#mtreadfraction
pdf("Mtreadfraction.pdf")
qplot(lib.sizes, mt.fraction, col = ifelse(mt.fraction>mt.lim, "drop", "keep")) +
  scale_x_log10() +
  labs(x = "UMI count", y = "MT read fraction") +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()
print(plot)

dim(counts[,mt.fraction < mt.lim])
#[1] 112627  15897
dim(counts[,mt.fraction < mt.lim])
#[1] 112627      0

dim(counts[,mt.fraction < 0.2])
#[1] 112627  16069

mt.lim <- 0.2

#chooseing the counts with mt. fraction lover than treshold for sce.
sce <- SingleCellExperiment(list(counts=counts[,mt.fraction < mt.lim]),
  colData=DataFrame(metadata[mt.fraction < mt.lim,]))
rownames(sce) <- genes$gene_id

rownames(genes) <- rownames(sce)
rowData(sce) <- DataFrame(genes)

colnames(sce) <- metadata$bc_wells[mt.fraction  < mt.lim]
colData(sce)  <- DataFrame(metadata[mt.fraction < mt.lim,])


lib.sizes <- colSums(counts(sce))
sce_filt  <- sce[calculateAverage(sce)>0.05,]

clusts <- as.numeric(quickCluster(sce_filt, method = "igraph", min.size = 100))

min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce_filt <- computeSumFactors(sce_filt, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

sizeFactors(sce) <- sizeFactors(sce_filt)

pdf("sizefactors.pdf")
ggplot(data = data.frame(X = lib.sizes, Y = sizeFactors(sce)), mapping = aes(x = X, y = Y)) +
  geom_point() +
  scale_x_log10(breaks = c(500, 2000, 5000, 10000, 30000), labels = c("5,00", "2,000", "5,000", "10,000", "30,000") ) +
  scale_y_log10(breaks = c(0.2, 1, 5)) +
  theme_minimal() +
  theme(text = element_text(size=20))  +
  labs(x = "Number of UMIs", y = "Size Factor")
dev.off()

ggplot(data.frame(colData(sce)), aes (x = factor(sample_name), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("UMIsBySample_afterQC.pdf")

library(BiocParallel)

bp <- MulticoreParam(12, RNGseed=1234)
bpstart(bp)
#detection and evaluation of doublets/multiplets
sce <- scDblFinder(sce, samples="bc1_well", dbr=.03, dims=30, BPPARAM=bp)
bpstop(bp)
table(sce_human$scDblFinder.class)
#singlet doublet 
#  49995    1915

#normalisation
sce_filt <- sce[calculateAverage(sce)>0.05,]
sce_filt <- logNormCounts(sce_filt)

decomp <- modelGeneVar(sce_filt)
hvgs   <- rownames(decomp)[decomp$FDR < 0.5]
pca    <- prcomp_irlba(t(logcounts(sce_filt[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_filt)
tsne <- Rtsne(pca$x, pca = FALSE, check_duplicates = FALSE, num_threads=30)









layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=30, min_dist=.25)

df_plot <- data.frame(
 colData(sce),
 doublet  = colData(sce)$scDblFinder.class,
 tSNE1    = tsne$Y[, 1],
 tSNE2    = tsne$Y[, 2], 
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

plot.index <- order(df_plot$doublet)
ggplot(df_plot[plot.index,], aes(x = tSNE1, y = tSNE2, col = factor(doublet))) +
  geom_point(size = 0.4) +
  scale_color_manual(values=c("gray","#0169c1"), name = "") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("tsne_doublets.pdf")

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(doublet))) +
  geom_point(size = 0.4) +
  scale_color_manual(values=c("gray","#0169c1"), name = "") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))

colData(sce) <- DataFrame(df_plot)

saveRDS(sce,paste0(path2data,"sce.rds"))





