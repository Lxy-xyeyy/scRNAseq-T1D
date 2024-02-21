# Find doublets in scRNA-seq dataset using scDblFinder
library(scDblFinder)
library(Seurat)

targetcell<-readRDS("./data/targetcell.rds")
levels(targetcell)
quantile(targetcell@meta.data$nFeature_RNA,c(0.00,0.50,1.00))
quantile(targetcell@meta.data$nCount_RNA,c(0.00,0.50,1.00))
quantile(targetcell@meta.data$percent.mt,c(0.00,0.50,1.00))

library(BiocParallel)
set.seed(123)
sce <- as.SingleCellExperiment(targetcell)

#dbrset=0.15
#dbrset=0.10
#dbrset=0.05
#dbrset=0.01
sce <- scDblFinder(sce,samples="sample", dbr=dbrset)
table(sce$scDblFinder.class)
sce$doublet_logic <- ifelse(sce$scDblFinder.class == "doublet", TRUE, FALSE)
#plotDoubletMap(sce)
table(sce$scDblFinder.class)
scDblObj <- colData(sce)
scDblObj <- as.data.frame(scDblObj)
identical(rownames(scDblObj), rownames(targetcell@meta.data)) #TRUE
targetcell@meta.data$scDblFinder.class <- scDblObj$scDblFinder.class
targetcell1=subset(targetcell, subset = scDblFinder.class == "singlet")
targetcell2=subset(targetcell, subset = scDblFinder.class == "doublet")

VlnPlot(object = targetcell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0,group.by = "sample",split.by = "scDblFinder.class")