library(tximport)
library(readr)
library(dplyr)
setwd ("C:/Users/Alicia/Desktop/NKp46")

#find and download tabular file for canine from NCBI (feature_table.txt.gz)
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_feature_table.txt.gz",
              destfile = "inputs/GCF_000002285.3_CanFam3.1_feature_table.txt.gz")
feat_table <- read_tsv('inputs/GCF_000002285.3_CanFam3.1_feature_table.txt.gz')

#select relevant columns: product accession and GeneID (in that order)
feat_table <- select(feat_table, product_accession, symbol)
feat_table <- feat_table%>%
  unique()

#make sure table is capturing all salmon counts, should be 82430 obs in canine
quant <- read_tsv("outputs/Canine_NK_Day_0_HVE5_2_15_19_quant.sf")

table(quant$Name %in% feat_table$product_accession)
write.table(feat_table, "outputs/tx2gene/dogtx2gene.tsv", quote = F, row.names = F, sep = "\t")

# read in file names 
files <- read_csv("samples.csv")

# read in dog counts
tx2gene <- read_tsv("outputs/tx2gene/dogtx2gene.tsv")
canine <- tximport(files = files$files, type = "salmon", tx2gene = tx2gene)
caninecounts <- canine$counts

# change rownames for provenance
colnames(caninecounts) <- files$sample

# write to file
write.csv(caninecounts, "outputs/counts/counts.csv", quote = F, row.names = T)

library(DESeq2)

counts <- read.csv("outputs/counts/counts.csv",row.names = 1)
counts <- apply(counts, 1:2, round)
samples <- read.csv("samples.csv")


dds<- DESeqDataSetFromMatrix(counts,
                             colData = samples,
                             design = ~ combined)
dds <- DESeq(dds)

res1 <- results(dds, contrast=c("combined","CD5resting","NKp46resting"))
res2 <- results(dds, contrast=c("combined","NKp46resting","NKp46coculture"))
res3 <- results(dds, contrast=c("combined","NKp46coculture","Unbiasedcoculture"))
res4 <- results(dds, contrast=c("combined", "CD5resting", "Unbiasedcoculture"))
res5 <- results(dds, contrast=c("combined", "CD5resting", "NKp46coculture"))
res6 <- results(dds, contrast=c("combined", "NKp46resting", "Unbiasedcoculture"))

resOrdered1 <- res1[order(res1$pvalue),]
resOrdered2 <- res2[order(res2$pvalue),]
resOrdered3 <- res3[order(res3$pvalue),]
resOrdered4 <- res4[order(res4$pvalue),]
resOrdered5 <- res5[order(res5$pvalue),]
resOrdered6 <- res6[order(res6$pvalue),]

plotMA(res1, ylim=c(-10,10), alpha=0.05)
plotMA(res2, ylim=c(-10,10), alpha=0.05)
plotMA(res3, ylim=c(-10,10), alpha=0.05)
plotMA(res4, ylim=c(-10,10), alpha=0.05)
plotMA(res5, ylim=c(-10,10), alpha=0.05)
plotMA(res6, ylim=c(-10,10), alpha=0.05)

idx <- identify(res3$baseMean, res3$log2FoldChange)
rownames(res3)[idx]

plotCounts(dds, gene=which.min(res3$padj), intgroup="combined")
plotCounts(dds, gene=which(rownames(res1)=="KLRK1"), intgroup="combined")

write.csv(as.data.frame(resOrdered1), 
          file="res6_treated_results.csv")

vsd <- vst(dds, blind=FALSE)


library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("marker","treatment")])
heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                    cluster_cols=TRUE, annotation_col=df, legend = TRUE, annotation_colors = list(
                      marker = c(CD5dim = "#33CC33", NKp46 = "#CC0000", unbiased = "#6600FF"),
                      treatment = c(resting = "#CCFFFF", coculture = "#3399FF")
                    ))
library(ggplot2)

pcaData <- plotPCA(vsd, intgroup=c("marker", "treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=marker, shape=treatment)) +
  geom_point(size=3) +
  scale_color_manual(values = c(CD5dim = "#33CC33", NKp46 = "#CC0000", unbiased = "#6600FF"))+
  scale_shape_manual(values = c(resting = 16, coculture = 17))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_light()
