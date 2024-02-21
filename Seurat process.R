# Seurat process

sample1.data <-Read10X(data.dir ="~/sample1")
sample1 <- CreateSeuratObject(counts = sample1.data, project = "sample1", min.cells = 3, min.features = 200)
sample1@meta.data$sample <- 'sample1'

sample2.data <-Read10X(data.dir ="~/sample2")
sample2 <- CreateSeuratObject(counts = sample2.data, project = "sample2", min.cells = 3, min.features = 200)
sample2@meta.data$sample <- 'sample2'

sample3.data <-Read10X(data.dir ="~/sample3")
sample3 <- CreateSeuratObject(counts = sample3.data, project = "sample3", min.cells = 3, min.features = 200)
sample3@meta.data$sample <- 'sample3'

mergelist <-  list(sample1,sample2,sample3)
for (i in 1:length(mergelist)) {
  mergelist[[i]] <- NormalizeData(mergelist[[i]])
}
mergeALL <- merge(sample1,y=c(sample2,sample3))
mergeALL[["percent.mt"]] <- PercentageFeatureSet(object = mergeALL, pattern = "^MT-")

VlnPlot(object = mergeALL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)
VlnPlot(object = mergeALL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

mergeALL <- NormalizeData(mergeALL, normalization.method = "LogNormalize", scale.factor = 10000)

quantile(mergeALL@meta.data$nFeature_RNA,c(0.90,0.95,0.99))
quantile(mergeALL@meta.data$nCount_RNA,c(0.90,0.95,0.99))
quantile(mergeALL@meta.data$percent.mt,c(0.90,0.95,0.99))

mergeALL_filter <- subset(x = mergeALL, subset = 200<nFeature<6000 & 500<nCount<34000 & percent.mt < 13)    

VlnPlot(object = mergeALL_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)
VlnPlot(object = mergeALL_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

mergeALL_filter <- NormalizeData(mergeALL_filter, normalization.method = "LogNormalize", scale.factor = 10000)



# UMAP process
targetcell<-readRDS("./data/targetcell.rds")
levels(targetcell)

targetcell <- FindVariableFeatures(object = targetcell, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(targetcell)
targetcell <- ScaleData(targetcell, features = all.genes)

targetcell=RunPCA(object= targetcell,npcs = 40,pc.genes=VariableFeatures(object = targetcell))   

DimHeatmap(object = targetcell, dims = 1:40, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
targetcell <- JackStraw(object = targetcell, num.replicate = 100)
targetcell <- ScoreJackStraw(object = targetcell, dims = 1:40)
JackStrawPlot(object = targetcell, dims = 1:20)
ElbowPlot(object = targetcell)

targetcell <- targetcell %>% 
  RunHarmony("Batch", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(targetcell, 'harmony')

pcSelect=30
targetcell <- FindNeighbors(object = targetcell, dims = 1:pcSelect,reduction = "harmony")                
targetcell <- FindClusters(object = targetcell, resolution = 1.2) 
targetcell <- RunUMAP(object = targetcell, dims = 1:pcSelect,reduction = "harmony") 
targetcell <- RunTSNE(object = targetcell, dims = 1:pcSelect,reduction = "harmony") 

Idents(object=targetcell) <- targetcell@meta.data$celltype
levels(targetcell)

UMAPPlot(object = targetcell, pt.size = 0.05, label = F,raster=F,group.by="celltype",
         cols=c("#7F7F7F","#B54D86","#0DA17B","#007857","#65B7E5","#0D7AB5","#DB82A0","#E6A52F","#FC9999","#FF0000","#B983FF","#EEE44C","#D46932","#FCE61C"))

DotPlot(targetcell,features=c("CD68","MNDA","CD14","FCGR3A","CD1C","CLEC10A","CLEC4C","IL3RA","CLC","CD79A","MS4A1","IGHD","CD3D","CD4","CD8A","NCAM1","CD34","ITGA2B"),
        scale.min = 0,scale.max = 70,group.by="celltype",cols = c("#2A778E","#FCE61C"))

DotPlot(targetcell,features=c("CD68","MNDA","CD14","FCGR3A","CD1C","CLEC10A","CLEC4C","IL3RA","CLC","CD79A","MS4A1","IGHD","JCHAIN","CD27","CD3D","CCR7","CD4","CD8A","CD8B","NCAM1","KLRC1","CD34","ITGA2B"),
        scale.min = 0,scale.max = 70,group.by="celltype",cols = c("#2A778E","#FCE61C"))

prop <- prop.table(table(Idents(targetcell), targetcell@meta.data$sample), margin = 2)
write.csv(prop,'targetcell_proportion.csv')

logFCfilter=0.5
adjPvalFilter=0.05
merge1.markers <- FindAllMarkers(object = targetcell,
                                 group.by="celltype",
                                 only.pos = F,
                                 min.pct = 0.1,
                                 logfc.threshold = logFCfilter)
sig.markers=merge1.markers[(abs(as.numeric(as.vector(merge1.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(merge1.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,'cell_markers_celltype.csv')