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
sample_info <- read.table('/data2/hanna/synaptogenesis/newvolume/analysis/sample_info_alt.tab',
  sep = "\t", header = TRUE)

counts    <- t(readMM(paste0(path2data, "DGE.mtx")))
genes     <- read.csv(paste0(path2data, "all_genes.csv"))
metadata  <- read.csv(paste0(path2data, "cell_metadata.csv"))

lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

dim(counts)
#[1] 112627 381176
dim(counts[,ngenes > 500])
#[1] 112627  16073

sample_bc1_well     <- rep(NA, nrow(metadata))        
sample_number       <- rep(NA, nrow(metadata))
sample_name_human   <- rep(NA, nrow(metadata))
sample_name_mouse   <- rep(NA, nrow(metadata))

samples <- unique(sample_info$Sample_well)
for (i in 1:length(samples)){
  sample_bc1_well[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]    <- sample_info$Sample_well[i]
  sample_number[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]      <- sample_info$Sample_Number[i]
  sample_name_human[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]  <- sample_info$Sample_name_H[i]
  sample_name_mouse[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]  <- sample_info$Sample_name_M[i]
}
sample_name_human <- gsub(" ","_",sample_name_human)
sample_name_mouse <- gsub(" ","_",sample_name_mouse)

submeta_human <- data.frame(rlist::list.rbind(strsplit(sample_name_human, split="_")))
colnames(submeta_human) <- c("batch", "day", "replicate")

submeta_mouse <- data.frame(rlist::list.rbind(strsplit(sample_name_mouse, split="_")))
colnames(submeta_mouse) <- c("batch", "day", "replicate")

metadata <- data.frame(cbind(metadata, lib.sizes, sample_number, sample_bc1_well, sample_name_human, submeta_human, sample_name_mouse, submeta_mouse))

plot_df <- metadata

ggplot(plot_df, aes (x = factor(sample_name_human), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("/data2/hanna/synaptogenesis/newvolume/analysis/UMIsBySample_beforeQC_H.pdf")

ggplot(plot_df, aes (x = factor(sample_name_mouse), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("/data2/hanna/synaptogenesis/newvolume/analysis/UMIsBySample_beforeQC_M.pdf")
 
pdf("cell_complexity.pdf")
ggp <- qplot(lib.sizes, ngenes, col = ifelse(ngenes < 500, "drop", "keep")) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  labs(x = "UMI count", y = "Number of expressed genes") +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()
print(ggp)

dim(counts[,ngenes > 500])
#[1] 62703 52719

counts   <- counts[,ngenes > 500]
metadata <- metadata[ngenes > 500,]
lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

ensembl_human <- useEnsembl(biomart = "ensembl",  dataset = "hsapiens_gene_ensembl",mirror="useast")
ensembl_mouse <- useEnsembl(biomart = "ensembl",  dataset = "mmusculus_gene_ensembl",mirror="useast")

gene_map_H  <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  filters = "hgnc_symbol", values = genes$gene_name, mart = ensembl_human)
gene_map_M  <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "chromosome_name"),
  filters = "mgi_symbol", values = genes$gene_name, mart = ensembl_mouse)

mt.index_H    <- gene_map_H$chromosome_name == "MT"
mt.counts_H   <- counts[which(genes$gene_name %in% gene_map_H$hgnc_symbol[mt.index_H]), ]
mt.fraction_H <- colSums(mt.counts_H)/lib.sizes

mt.p_H   <- pnorm(mt.fraction_H, mean = median(mt.fraction_H), sd = mad(mt.fraction_H), lower.tail = FALSE)
mt.lim_H <- min(mt.fraction_H[which(p.adjust(mt.p_H, method = "fdr") < 0.05)])

mt.index_M    <- gene_map_M$chromosome_name == "MT"
mt.counts_M   <- counts[which(genes$gene_name %in% gene_map_M$mgi_symbol[mt.index_M]), ]
mt.fraction_M <- colSums(mt.counts_M)/lib.sizes

mt.p_M   <- pnorm(mt.fraction_M, mean = median(mt.fraction_M), sd = mad(mt.fraction_M), lower.tail = FALSE)
mt.lim_M <- min(mt.fraction_M[which(p.adjust(mt.p_M, method = "fdr") < 0.05)])

#Threhdold
mt.lim_H
mt.lim_M
#[1] 0.08196721

mt.lim_H <- min(mt.fraction_H[which(p.adjust(mt.p_H, method = "fdr") < 0.001)])
mt.lim_M <- min(mt.fraction_M[which(p.adjust(mt.p_M, method = "fdr") < 0.001)])

#Threhdold
mt.lim_H
mt.lim_M
#[1] 0.1046875

metadata <- data.frame(cbind(metadata,mt.fraction_H))
metadata <- data.frame(cbind(metadata,mt.fraction_M))

pdf("H_mtreadfraction1.pdf")
hplot <- qplot(lib.sizes, mt.fraction_H, col = ifelse(mt.fraction_H>mt.lim_H, "drop", "keep")) +
  scale_x_log10() +
  labs(x = "UMI count", y = "MT read fraction") +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()
print(hplot)

pdf("M_mtreadfraction1.pdf")
mplot <- qplot(lib.sizes, mt.fraction_M, col = ifelse(mt.fraction_M>mt.lim_M, "drop", "keep")) +
  scale_x_log10() +
  labs(x = "UMI count", y = "MT read fraction") +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()
print(mplot)

dim(counts[,mt.fraction_H < mt.lim_H])
#[1] 112627  15897
dim(counts[,mt.fraction_M < mt.lim_M])
#[1] 112627      0

dim(counts[,mt.fraction_H < 0.2])
#[1] 112627  16069
dim(counts[,mt.fraction_M < 0.2])
#[1] 112627  16073

mtlim <- 0.2

sce_human <- SingleCellExperiment(list(counts=counts[,mt.fraction_H < mt.lim_H]),
  colData=DataFrame(metadata[mt.fraction_H < mt.lim_H,]))
rownames(sce_human) <- genes$gene_id
sce_mouse <- SingleCellExperiment(list(counts=counts[,mt.fraction_M < mt.lim_M]),
  colData=DataFrame(metadata[mt.fraction_M < mt.lim_M,]))
rownames(sce_mouse) <- genes$gene_id

rownames(genes) <- rownames(sce_human)
rowData(sce_human) <- DataFrame(genes)
rownames(genes) <- rownames(sce_mouse)
rowData(sce_mouse) <- DataFrame(genes)

colnames(sce_human) <- metadata$bc_wells[mt.fraction_H  < mt.lim_H]
colData(sce_human)  <- DataFrame(metadata[mt.fraction_H < mt.lim_H,])
colnames(sce_mouse) <- metadata$bc_wells[mt.fraction_M  < mt.lim_M]
colData(sce_mouse)  <- DataFrame(metadata[mt.fraction_M < mt.lim_M,])

lib.sizes_human <- colSums(counts(sce_human))
sce_filt_human  <- sce_human[calculateAverage(sce_human)>0.05,]
lib.sizes_mouse <- colSums(counts(sce_mouse))
sce_filt_mouse  <- sce_mouse[calculateAverage(sce_mouse)>0.05,]

clusts_human <- as.numeric(quickCluster(sce_filt_human, method = "igraph", min.size = 100))
clusts_mouse <- as.numeric(quickCluster(sce_filt_mouse, method = "igraph", min.size = 100))

min.clust_human <- min(table(clusts_human))/2
new_sizes_human <- c(floor(min.clust_human/3), floor(min.clust_human/2), floor(min.clust_human))
sce_filt_human <- computeSumFactors(sce_filt_human, clusters = clusts_human, sizes = new_sizes_human, max.cluster.size = 3000)
min.clust_mouse <- min(table(clusts_mouse))/2
new_sizes_mouse <- c(floor(min.clust_mouse/3), floor(min.clust_mouse/2), floor(min.clust_mouse))
sce_filt_mouse <- computeSumFactors(sce_filt_mouse, clusters = clusts_mouse, sizes = new_sizes_mouse, max.cluster.size = 3000)

sizeFactors(sce_human) <- sizeFactors(sce_filt_human)
sizeFactors(sce_mouse) <- sizeFactors(sce_filt_mouse)

pdf("sizefactors_H.pdf")
ggplot(data = data.frame(X = lib.sizes_human, Y = sizeFactors(sce_human)), mapping = aes(x = X, y = Y)) +
  geom_point() +
  scale_x_log10(breaks = c(500, 2000, 5000, 10000, 30000), labels = c("5,00", "2,000", "5,000", "10,000", "30,000") ) +
  scale_y_log10(breaks = c(0.2, 1, 5)) +
  theme_minimal() +
  theme(text = element_text(size=20))  +
  labs(x = "Number of UMIs", y = "Size Factor")
dev.off()

pdf("sizefactors_M.pdf")
ggplot(data = data.frame(X = lib.sizes_mouse, Y = sizeFactors(sce_mouse)), mapping = aes(x = X, y = Y)) +
  geom_point() +
  scale_x_log10(breaks = c(500, 2000, 5000, 10000, 30000), labels = c("5,00", "2,000", "5,000", "10,000", "30,000") ) +
  scale_y_log10(breaks = c(0.2, 1, 5)) +
  theme_minimal() +
  theme(text = element_text(size=20))  +
  labs(x = "Number of UMIs", y = "Size Factor")
dev.off()

ggplot(data.frame(colData(sce_human)), aes (x = factor(sample_name_human), y = as.numeric(lib.sizes_human))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("UMIsBySample_afterQC_H.pdf")

ggplot(data.frame(colData(sce_mouse)), aes (x = factor(sample_name_mouse), y = as.numeric(lib.sizes_mouse))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("UMIsBySample_afterQC_M.pdf")


library(BiocParallel)
bp <- MulticoreParam(12, RNGseed=1234)
bpstart(bp)
sce_human <- scDblFinder(sce_human, samples="bc1_well", dbr=.03, dims=30, BPPARAM=bp)
bpstop(bp)
table(sce_human$scDblFinder.class)
#singlet doublet 
# 15445     817 

bp <- MulticoreParam(12, RNGseed=1234)
bpstart(bp)
sce_mouse <- scDblFinder(sce_mouse, samples="bc1_well", dbr=.03, dims=30, BPPARAM=bp)
bpstop(bp)
table(sce_mouse$scDblFinder.class)
#singlet doublet 
#  49995    1915

sce_filt_human <- sce_human[calculateAverage(sce_human)>0.05,]
sce_filt_mouse <- sce_mouse[calculateAverage(sce_mouse)>0.05,]

sce_filt_human <- logNormCounts(sce_filt_human)
sce_filt_mouse <- logNormCounts(sce_filt_mouse)

decomp_human  <- modelGeneVar(sce_filt_human)
hvgs_human    <- rownames(decomp_human)[decomp_human$FDR < 0.5]
pca_human     <- prcomp_irlba(t(logcounts(sce_filt_human[hvgs_human,])), n = 30)
rownames(pca_human$x) <- colnames(sce_filt_human)
tsne <- Rtsne(pca_human$x, pca = FALSE, check_duplicates = FALSE, num_threads=30)

decomp_mouse  <- modelGeneVar(sce_filt_mouse)
hvgs_mouse    <- rownames(decomp_mouse)[decomp_mouse$FDR < 0.5]
pca_mouse     <- prcomp_irlba(t(logcounts(sce_filt_mouse[hvgs_mouse,])), n = 30)
rownames(pca_mouse$x) <- colnames(sce_filt_mouse)
tsne <- Rtsne(pca_mouse$x, pca = FALSE, check_duplicates = FALSE, num_threads=30)

######################################################################
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






