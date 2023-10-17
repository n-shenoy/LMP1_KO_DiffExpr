# Author: Navami Shenoy
# script for differential gene expression analysis of control vs. 
# LMP1 knockout lymphoblastoid cell lines
setwd("F:/bulkRNASeq/GM12878")

#load libraries
library(dplyr) 
library(tidyverse)
library(tibble)
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(ggridges)
library(clusterProfiler)
library(EnhancedVolcano)
library("EnsDb.Hsapiens.v86")
library("org.Hs.eg.db")


#-------Read in the counts matrix--------#

counts_matrix <- read.csv("counts/counts_matrix_5_samples.csv", header = TRUE, 
                          row.names = 1, sep = "\t", skip = 1)

#change column names to reflect sample IDs
colnames(counts_matrix) <- c("Control_Rep2", "Control_Rep3",
                             "LMP1_KO_Rep1", "LMP1_KO_Rep2", "LMP1_KO_Rep3")
head(counts_matrix)


#------------Pre-filtering--------------#

#keep rows with more than 10 reads total
counts_matrix <- counts_matrix[which(rowSums(counts_matrix) > 10), ]
#reduced the number of rows by more than half!


#---------------DESeq2------------------#

#create dataframe containing sample information
condition <- factor(c("control", "control",
                      "knockout", "knockout", "knockout"))
coldata <- data.frame(row.names = colnames(counts_matrix), condition)
head(coldata)

#construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata,
                              design = ~condition)

#set reference level
dds$condition <- relevel(dds$condition, ref = "control")

#run DESeq2
dds <- DESeq(dds)

#change p-value
res_0.05 <- results(dds, alpha = 0.05, contrast = c("condition", "knockout", "control"))

#drop NAs
sigs <- na.omit(res_0.05)

#get a summary of the dds object
summary(res_0.05)


#------------Collect significant genes-------------#


#get differentially expressed genes with adjusted p-value below 0.05
sigs <- sigs[sigs$padj < 0.05, ]
sigs

#map gene IDs to gene names 
sigs.df <- as.data.frame(sigs)
sigs.df$symbol <- ensembldb::mapIds(EnsDb.Hsapiens.v86, keys = row.names(sigs.df), 
                                    keytype = "GENEID", column = "SYMBOL")


#use gene ID as gene name for genes that did not 
#map to a symbol in the previous step
for(i in seq_along(sigs.df$symbol)){
    if(is.na(sigs.df$symbol[i])){
        sigs.df$symbol[i] = row.names(sigs.df)[i]
    }
}


#make a copy so the original df is available if 
#we need it later
sigs.df.copy <- sigs.df



#-----------Visualization------------#

#PCA plot
vsdata <- vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "condition")

#Dispersion plot
plotDispEsts(dds)

#MA plot 
plotMA(res_0.05)

#Volcano plot
select.genes <-  c("EBI3", "CCL22", "CD40", "MDM2","BATF", 
                   "NFKB1", "NFKB2", "IRF4", "CBS", "BCAT1")


keyvals <- ifelse(
    sigs.df.copy$log2FoldChange < -1, '#3ebdc8',
    ifelse(sigs.df.copy$log2FoldChange > 1, '#F6BE00',
           'darkgrey'))

#reorder df
x <- sigs.df.copy$symbol %in% select.genes
sigs.reordered <- rbind(sigs.df.copy[!x,], sigs.df.copy[x,])

#generate key-value pairs to change color and label schemes
keyvals <- ifelse(
    sigs.reordered$symbol %in% select.genes, 
    ifelse(sigs.reordered$log2FoldChange < 0, '#3ebdc8',
           ifelse(sigs.reordered$log2FoldChange > 0, '#F6BE00', 
                  'darkgrey')),
    'darkgrey'
)

keyvals[is.na(keyvals)] <- 'darkgrey'
names(keyvals)[keyvals == '#F6BE00'] <- 'Upregulated'
names(keyvals)[keyvals == 'darkgrey'] <- 'Not significant'
names(keyvals)[keyvals == '#3ebdc8'] <- 'Downregulated'

EnhancedVolcano(sigs.reordered, x = "log2FoldChange", y = "padj", lab = sigs.reordered$symbol,
                pCutoff = 0.005, FCcutoff = 1, 
                pointSize = (ifelse(sigs.reordered$symbol %in% select.genes, 3,1)),
                labSize = 4, 
                selectLab = select.genes,
                drawConnectors = TRUE,
                title = NULL, 
                subtitle = NULL, 
                caption = NULL,
                colCustom = keyvals,
                colAlpha = (ifelse(sigs.reordered$symbol %in% select.genes, 15, 0.8)),
                axisLabSize = 10,
                legendPosition = 'none',    
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                xlim = c(-27,27))


#Heatmap
#get normalized counts for significant genes from the dds object
mat <- counts(dds, normalized = TRUE)[rownames(sigs.df.copy), ]

#get the z-score for each row
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- colnames(mat)


Heatmap(mat.z, cluster_rows = TRUE, cluster_columns = TRUE, column_labels = colnames(mat.z), 
        name = "Z-score", 
        column_names_gp = grid::gpar(fontsize = 6),
        column_names_side = "top", column_names_rot = 0, 
        column_names_centered = TRUE,
        col = c("#0056b3","white","#ff0055"),
        show_row_names = FALSE,
        show_row_dend = FALSE,
        km = 2)



#-------------Gene Ontology------------#
#over-representated pathway analysis

#get ENSEMBL IDs of genes with highly altered expression
#upregulated genes
upreg.genes <- row.names(sigs.df.copy[sigs.df.copy$log2FoldChange > 1, ])
upreg.pathways <- enrichGO(gene = upreg.genes, OrgDb = "org.Hs.eg.db", 
                           keyType = "ENSEMBL", ont = "BP")

head(as.data.frame(upreg.pathways))

upreg <- upreg.pathways %>% dplyr::arrange(desc(Count))
head(upreg)

#plot upregulated pathways
fit_upreg <- plot(barplot(upreg.pathways, showCategory = 20))

#plot only top 10 
ggplot(upreg[1:10,], aes(y = reorder(Description, +Count), x = Count, fill = p.adjust)) + 
    geom_col() +
    scale_fill_gradient(low="blue", high="red") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

#downregulated genes
downreg.genes <- rownames(sigs.df.copy[sigs.df.copy$log2FoldChange < -1, ])
downreg.pathways <- enrichGO(gene = downreg.genes, OrgDb = "org.Hs.eg.db", 
                             keyType = "ENSEMBL", ont = "BP")

head(as.data.frame(downreg.pathways))

downreg <- downreg.pathways %>% dplyr::arrange(desc(Count))
head(downreg)

#plot downregulated pathways
fit_downreg <- plot(barplot(downreg.pathways, showCategory = 20))

#plot only top 10 
ggplot(downreg[1:10,], aes(y = reorder(Description, +Count), x = Count, fill = p.adjust)) + 
    geom_col() +
    scale_fill_gradient(low="blue", high="red") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())




#-----------Gene set enrichment analysis----------#

#--------GSEA GO---------#

#order genes by log2 fold change
res_gse <- sigs.df.copy[order(-sigs.df.copy$log2FoldChange), ] 
genes_list <- res_gse$log2FoldChange
names(genes_list) <- rownames(sigs.df.copy)


#run GSEA
gsea.go <- gseGO(genes_list, ont = "BP", keyType = "ENSEMBL", 
                 OrgDb = "org.Hs.eg.db")

gsea.df <- as.data.frame(gsea.go)
gsea.df.desc <- gsea.df[order(-gsea.df$NES), ]

#GSEA plot
#Downregulated top 3
#(none were upregulated)
n <- nrow(gsea.df)
gseaplot(gsea.go, geneSetID = which(gsea.df$Description == gsea.df.desc$Description[n]), 
         title = gsea.df.desc$Description[n])
gseaplot(gsea.go, geneSetID = which(gsea.df$Description == gsea.df.desc$Description[n-1]), 
         title = gsea.df.desc$Description[n-1])
gseaplot(gsea.go, geneSetID = which(gsea.df$Description == gsea.df.desc$Description[n-2]),
         title = gsea.df.desc$Description[n-2])


#Ridge plot
ridgeplot(gsea.go) + labs(x = "enrichment distribution")

#--------------GSEA KEGG-----------------#

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
sigs.entrez <- sigs.df.copy
sigs.entrez$entrezId <- ensembldb::mapIds(EnsDb.Hsapiens.v86, keys = row.names(sigs.df.copy), 
                                          keytype = "GENEID", column = "ENTREZID")

# omit any NA values 
sigs.entrez <- na.omit(sigs.entrez)

#-----deal with duplicates------#

#check for duplicate gene names
#they need to be dealt with for GSEA
nrow(sigs.entrez[duplicated(sigs.entrez$entrezId), ])
#10 ENTREZ IDs with duplicates

#extract ENTREZ IDs of duplicates
dupes <- sigs.entrez[duplicated(sigs.entrez$entrezId), ]$entrezId 

#filter out rows with duplicated ENTREZ IDs to view them
y <- sigs.entrez[sigs.entrez$entrezId %in% dupes, ]
y[order(y$entrezId, decreasing = TRUE), ]


#for ENSEMBL IDs that map to the same gene symbol,
#we'll keep the ones with the largest base mean
for(entrez in dupes){
    to_drop = min(sigs.entrez$baseMean[c(which(sigs.entrez$entrezId == entrez))])
    sigs.entrez <- sigs.entrez[-c(which(sigs.entrez$baseMean == to_drop)), ] 
}



# Create a vector of the log2 fold changes
kegg_gene_list <- sigs.entrez$log2FoldChange

# Name vector with ENTREZ IDs
names(kegg_gene_list) <- sigs.entrez$entrezId

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

#create KEGG object
gsea.kegg <- gseKEGG(kegg_gene_list, organism = "hsa", 
                     minGSSize = 1, maxGSSize = 20000, pvalueCutoff = 0.05,
                     pAdjustMethod = "BH", keyType = "ncbi-geneid")

kegg.df <- as.data.frame(gsea.kegg) 

kegg.df <- kegg.df[order(kegg.df$NES, decreasing = TRUE), ]

#Bar plot 
ggplot(kegg.df, aes(y = reorder(Description, -NES), x = NES, fill = p.adjust)) + 
    geom_col() +
    scale_fill_gradient(low="blue", high="red") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

#cnetplot
gsea.kegg2 <- setReadable(gsea.kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
kegg2.df <- as.data.frame(gsea.kegg2)
options(ggrepel.max.overlaps = Inf)
cnetplot(gsea.kegg2, categorySize = "p.adjust", color.params = list(foldChange = kegg_gene_list)) # categorySize can be either 'pvalue' or 'geneNum'

#Ridge plot
ridgeplot(gsea.kegg2) + labs(x = "enrichment distribution")

#GSEA plot
#top 3 most negative (downregulated) NES 
gseaplot(gsea.kegg, by = "all", title = gsea.kegg$Description[1], geneSetID = 1)
gseaplot(gsea.kegg, by = "all", title = gsea.kegg$Description[2], geneSetID = 2)
gseaplot(gsea.kegg, by = "all", title = gsea.kegg$Description[3], geneSetID = 3)
