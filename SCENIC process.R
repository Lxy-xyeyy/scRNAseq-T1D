library(Seurat)
library(AUCell)
library(GENIE3)
library(SCENIC)
library(R2HTML)
library(rbokeh)
library(doParallel)

targetcell<-readRDS("./data/targetcell.rds")

scenicOptions <- initializeScenic(org = "hgnc", dbDir = "cisTarget_databases", 
                                  dbs=c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"),
                                  nCores = 8)

exprMat123 <- as.matrix(targetcell@assays$RNA@counts)

exprMat123[1:4,1:4]
cellInfo <-targetcell@meta.data[,c("seurat_clusters","nCount_RNA","nFeature_RNA","group","sample")]
colnames(cellInfo) <- c("seurat_clusters",  'nGene' ,'nUMI',"group","sample")
head(cellInfo)
table(cellInfo$celltype)
exprMat_log <- log2(exprMat123+1)


library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log) 

aucellApp <- plotTsne_AUCellApp(scenicOptions,exprMat_log)
savedSelections <- shiny::runApp(aucellApp)

newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions,exprMat_log,skipBoxplot=F,skipHeatmaps=F,skipTsne=F)


tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat123)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC

targetcellauc <- AddMetaData(targetcell, AUCmatrix)
targetcellauc@assays$integrated <- NULL
saveRDS(targetcellauc,'targetcellauc.rds')


BINmatrix <- readRDS("int/4.2_binaryRegulonActivity_nonDupl.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN

targetcellbin <- AddMetaData(targetcell, BINmatrix)
targetcellbin@assays$integrated <- NULL
saveRDS(targetcellbin, 'targetcellbin.rds')

write.csv(names(AUCmatrix),"AUCmatrix_TF.csv")
write.csv(names(BINmatrix),"BINmatrix_TF.csv")

a <- targetcellauc@meta.data
write.csv(a,"targetcell_AUC.csv")
a <- targetcellbin@meta.data
write.csv(a,"targetcell_BIN2.csv")



targetcellbin <- readRDS("targetcellbin.rds")
targetcellauc <- readRDS("targetcellauc.rds")

#FeaturePlot
p1 = FeaturePlot(targetcellauc, features='GATA3_extended_52g', label=T, reduction = 'umap')
p2 = FeaturePlot(targetcellbin, features='GATA3_extended_52g', label=T, reduction = 'umap')
p3 = DimPlot(targetcell, reduction = 'umap', group.by = "celltype", label=T)
plotc = p1|p2|p3
plotc
library(ggplot2)


p1 =FeaturePlot(targetcellbin, features='GATA3_extended_52g', label=T, reduction = 'umap')
p2 =FeaturePlot(targetcellbin, features='TBX21_extended_92g', label=T, reduction = 'umap')
plotc = p1|p2
plotc
ggsave('GATA3_extended_52g.pdf', plotc, width=14 ,height=4)

colss2<-c("#DE77AE","#8E0152","#FB9A99","#B2DF8A", "#006D2C")
pdf(file="GATA3.2.pdf",width=5,height=3.5)
VlnPlot(object = targetcellauc, features = c("GATA3_extended_52g"),group.by="sample",pt.size = 0,cols=colss2) +geom_boxplot(width=0.2,col="black",fill="white")
dev.off()



cellInfo <- data.frame(targetcell@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="sample")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
saveRDS(cellInfo, file="cellInfo.Rds")

library(pheatmap)
cellInfo <- readRDS("cellInfo.Rds")
celltype = subset(cellInfo,select = 'celltype')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

my.regulons <- c("TCF7_extended_101g","TBX21_extended_92g","GATA3_extended_52g")
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]

pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype,
         width = 6, height = 5,cluster_rows=T,cluster_cols=T)

pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         color = colorRampPalette(colors = c("white","black"))(100),
         width = 6, height = 5,cluster_rows=T,cluster_cols=T)