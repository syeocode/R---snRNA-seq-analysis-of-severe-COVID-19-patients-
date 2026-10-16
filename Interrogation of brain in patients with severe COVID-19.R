set.seed(1)
library(data.table)
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(reticulate)
library(harmony)
library(DoubletFinder)

scale_colour_discrete=function(...) {scale_colour_brewer(..., palette="Set2")}
scale_fill_discrete=function(...) {scale_fill_brewer(...,palette="Set2")}
color=c(RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(6,"Set3"),RColorBrewer::brewer.pal(10,"RdYlGn"),RColorBrewer::brewer.pal(10,"BrBG"),RColorBrewer::brewer.pal(10,"PiYG"),RColorBrewer::brewer.pal(10,"PuOr"),RColorBrewer::brewer.pal(10,"RdBu"),RColorBrewer::brewer.pal(10,"PRGn"))

# read --------------------------------------------
# Data can be found at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164485
gc()
file = 'C:/Users/palan/Desktop/stat430/Finalproject/GSE164485_gene_expression_20210104.csv.gz'
exprf = fread(file)
expr = as.matrix(exprf[,-1])
rownames(expr) = exprf[,V1]
colnames(expr) = colnames(exprf)[-1]
dim(expr)
expr = t(expr)

file = 'C:/Users/palan/Desktop/stat430/Finalproject/GSE164485_meta_data_20210104.csv.gz'
meta = fread(file)
meta = as.data.frame(meta)
row.names(meta) = meta[,'V1']
meta = meta[,-1]




# make it into smaller  --------------------------------------------
# # samidx = sample(x = 1:dim(expr)[2],size=5000, replace=FALSE)
# geneidx = sample(x = 1:dim(expr)[1],size=5000, replace=FALSE)
rn = colnames(expr)
rn = rn[samidx]
sexpr = expr[sort(geneidx),rn]
# sexpr = t(sexpr)
smeta = meta[rn,]

# Trouble shooting --- sample not matching-- DONE!
# rownames(smeta[1:10,])
# colnames(sexpr[1:3,1:10])
# setdiff(colnames(sexpr), rownames(smeta))


# Testing --- 
# df <- data.frame(umap1= smeta$UMAP1, umap2= smeta$UMAP2)
# tmp = data.frame()
#   ggplot(
#     df, 
#     aes(x=umap1,
#         y=umap2,
#         color = 'k' )
#   ) +
#   geom_point() + 
#   xlab('UMAP1') + ylab('UMAP2')


# create seurat obj  --------------------------------------------
sexpr= expr
smeta= meta
gc()
sobj = CreateSeuratObject(counts = sexpr, meta.data= smeta)
abrain = sobj 
abrain[["percent.mito"]] = PercentageFeatureSet(object = abrain, pattern = "^MT-")
abrain = subset(abrain, subset = nFeature_RNA > 200 & percent.mito < 5)
abrain=NormalizeData(object = abrain, normalization.method = "LogNormalize", scale.factor = 10000)
abrain=FindVariableFeatures(object = abrain, selection.method = "vst", nfeatures = 2000)

abrain=ScaleData(object = abrain, vars.to.regress = c("percent.mito"))
abrain=RunPCA(object = abrain, features = VariableFeatures(object = abrain), npcs = 100, verbose = FALSE)
# abrain=ProjectDim(object = abrain, verbose = FALSE)
npc=20
abrain=FindNeighbors(object = abrain, dims = 1:npc, k.param = 100) 
# abrain=FindClusters(object = abrain, resolution = 0.2, graph.name = 'RNA_snn')
abrain=FindClusters(object = abrain, resolution = 0.2)
# abrain=RunUMAP(object = abrain, min.dist = 1.5, graph = 'RNA_snn') 
abrain=RunUMAP(object = abrain, dims=1:20) 
# gc()


# visualize  --------------------------------------------

pdf(paste0("2.seurat.analysis.pdf"))
DimPlot(object = abrain, reduction='pca', cols=color, group.by = 'clinical')
ElbowPlot(object = abrain, ndims = 100)
DimPlot(abrain, cols=color)
DimPlot(abrain, cols=color, label=T)
DimPlot(abrain, cols=color, group.by = 'subject_id')
DimPlot(abrain, cols=color, group.by = 'clinical')
DimPlot(abrain, cols=color, group.by = 'is_covid')
DimPlot(abrain, cols=color, group.by = 'tissue')
DimPlot(abrain, cols=color, group.by = 'celltype')
FeaturePlot(abrain, features = c("TYROBP"))
FeaturePlot(object = abrain, reduction='umap', features = c('nCount_RNA'))
FeaturePlot(object = abrain, reduction='umap', features = c('nFeature_RNA'))
dev.off()


saveRDS(abrain, file = "seurat.object.rds")

# Differentially expressed gene --- 
# abrain.markers <- FindAllMarkers(abrain, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
abrain.markers <- FindAllMarkers(abrain, min.pct = 0.3, logfc.threshold = 0.25)



saveRDS(abrain.markers, file = paste0("seurat.marker.list.rds"))
topgene=abrain.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(topgene, paste0('seurat_topgene.csv'))
write.csv(as.data.frame(abrain.markers), paste0('seurat_markers.csv') )

#

nsrt = abrain 
Idents(nsrt) <- 'is_covid'
VlnPlot(nsrt, features = "TYROBP")

covid.de.markers <- FindMarkers(nsrt, ident.1 = "1", ident.2 = "0")
covid.de.markers


write.csv(as.data.frame(covid.de.markers), paste0('seurat_covid_markers.csv') )
#




