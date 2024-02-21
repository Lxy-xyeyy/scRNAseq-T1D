library(monocle3)
library(Seurat)
library(ggplot2)

targetcell<-readRDS("./data/targetcell.rds")

data <- GetAssayData(targetcell, assay ='RNA', slot = 'counts')
cell_metadata <- targetcell@meta.data
gene_annotation <-data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <-rownames(data)
cds <- new_cell_data_set(data,cell_metadata =cell_metadata, gene_metadata =gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = "PCA") 
plot_cells(cds)
cds <- cluster_cells(cds) 
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('cds.umap')

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(targetcell, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('int.umap')

p = p1|p2
p
ciliated_genes <- c("CD4","CD52","JUN")
plot_cells(cds,genes=ciliated_genes,label_cell_groups=FALSE,show_trajectory_graph=FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype2", label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
plot_cells(cds, color_cells_by = "celltype2", label_groups_by_cluster=FALSE,label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=3)

cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", 
           label_cell_groups = FALSE, label_leaves = F,label_branch_points = F,
           cell_size = 0.2,alpha = 2,trajectory_graph_segment_size = 0.9,
           trajectory_graph_color = "grey28",rasterize = F)

Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=8)
Track_genes <-Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes,"Trajectory_genes.csv", row.names = F)

Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
Track_genes_sig <- c("TIGIT","CCR7","FOXP3","IL2RA","TGFB1","CTLA4")

colss1=c("#E74A32","#00A188","#8592B5","#8ED0C1","#52BED7","#395389")

P2=plot_genes_in_pseudotime(cds["TIGIT",],cell_size = 2.5,color_cells_by="celltype2",min_expr=0.5, ncol= 2) + scale_color_manual(values=colss1)
P1=plot_genes_in_pseudotime(cds["CCR7",],cell_size = 2.5,color_cells_by="celltype2",min_expr=0.5, ncol= 2)  + scale_color_manual(values=colss1)
P3=plot_genes_in_pseudotime(cds["FOXP3",],cell_size = 2.5,color_cells_by="celltype2",min_expr=0.5, ncol= 2)  + scale_color_manual(values=colss1)
P4=plot_genes_in_pseudotime(cds["IL2RA",],cell_size = 2.5,color_cells_by="celltype2",min_expr=0.5, ncol= 2)  + scale_color_manual(values=colss1)
P5=plot_genes_in_pseudotime(cds["TGFB1",],cell_size = 2.5,color_cells_by="celltype2",min_expr=0.5, ncol= 2)  + scale_color_manual(values=colss1)
P6=plot_genes_in_pseudotime(cds["CTLA4",],cell_size = 2.5,color_cells_by="celltype2",min_expr=0.5, ncol= 2)  + scale_color_manual(values=colss1)

P1+P2+P3+P4+P5+P6+ plot_layout(ncol = 6, guides = "collect")
ggsave("tmp2-6GENE3.pdf",width = 34,height =8,units = "cm")


colData(cds)$monocle3_partitions <- as.character(partitions(cds))
colData(cds)$monocle3_clusters <- as.character(clusters(cds))
colData(cds)$monocle3_pseudotime <- as.numeric(pseudotime(cds))