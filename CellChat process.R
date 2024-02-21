library(CellChat)
library(ggplot2)
library(Seurat)
options(stringsAsFactors = FALSE)

celldata.final <-readRDS("./data/targetcell.rds")

Idents(object=celldata.final)<-celldata.final@meta.data$celltype
levels(celldata.final)

data.input  <- celldata.final@assays$RNA@data
identity = data.frame(group =celldata.final$celltype, row.names = names(celldata.final$celltype)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels

cellchat <- createCellChat(data.input)
cellchat
summary(cellchat)

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
unique(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell Contact","ECM-Receptor")) # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

cellchat@netP$pathways
head(cellchat@LR$LRsig)

pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
#vertex.receiver = c(1,11)
dir.create("all_pathways_com_circle") 
setwd("all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
            vertex.receiver = vertex.receiver, layout = "circle") #????????ͼ
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
         plot=gg, width = 10, height = 5, units = 'in', dpi = 300)
}
setwd("../")

levels(cellchat@idents) # show factor levels of the cell labels


p1 = netVisual_bubble(cellchat, sources.use = c(5,6,7), signaling=c("IL10","TGFb"),
                      targets.use = c(2,3,4), remove.isolate = F)+  theme(
                        axis.text.x = element_text(angle = 35,hjust=1,vjust = 1))
p1

ggsave("TGFb-IL10-2.pdf", p1, width = 6.5, height = 3.1) 

vertex.receiver = c(2,3,4,5,6,7) # a numeric vector

###TGFb##################
netAnalysis_contribution(cellchat, signaling = "TGFb")
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = "TGFb", geneLR.return = FALSE) #??ȡ??TGFb?й??׵????????????? 
LR.show <- pairLR.TGFb[1,] 

netVisual_chord_cell(cellchat,signaling ="TGFb",sources.use = c(5,6,7),targets.use = c(2,3,4))
netVisual_aggregate(cellchat, signaling = "TGFb",  vertex.receiver = c(2,3,4),sources.use = c(5,6,7), layout = c( "hierarchy"))

netVisual_aggregate(cellchat, signaling = "TGFb",  vertex.receiver = vertex.receiver,  sources.use = c(5,6,7,8), targets.use = c(2,3,4))

netVisual_individual(cellchat, signaling ="TGFb",  vertex.receiver = vertex.receiver,  sources.use = c(5,6,7), targets.use = c(2,3,4))

netVisual_aggregate(cellchat, signaling = "TGFb",  vertex.receiver = vertex.receiver,  sources.use = c(5,6,7), targets.use = c(2,3,4))
netVisual_individual(cellchat, signaling ="TGFb",  pairLR.use = pairLR.TGFb[1,]  ,vertex.receiver = vertex.receiver,  sources.use = c(5,6,7), targets.use = c(2,3,4))
netVisual_individual(cellchat, signaling ="TGFb",  pairLR.use = pairLR.TGFb[2,]  ,vertex.receiver = vertex.receiver,  sources.use = c(5,6,7), targets.use = c(2,3,4))
netVisual_individual(cellchat, signaling ="TGFb",  pairLR.use = pairLR.TGFb[3,]  ,vertex.receiver = vertex.receiver,  sources.use = c(5,6,7), targets.use = c(2,3,4))


###IL10##################
netAnalysis_contribution(cellchat, signaling = "IL10")
pairLR.IL10 <- extractEnrichedLR(cellchat, signaling = "IL10", geneLR.return = FALSE) #??ȡ??IL10?й??׵????????????? 
LR.show <- pairLR.IL10[1,] 
netVisual_chord_cell(cellchat,signaling ="IL10",sources.use = c(6,7,8),targets.use = c(2,3,4,5))

netVisual_aggregate(cellchat, signaling = "IL10",  vertex.receiver = c(2,3,4,5),sources.use = c(6,7,8), layout = c( "hierarchy"))


netVisual_aggregate(cellchat, signaling = "IL10",  vertex.receiver = vertex.receiver,  sources.use = c(6,7,8), targets.use = c(2,3,4,5))
netVisual_individual(cellchat, signaling ="IL10",  vertex.receiver = vertex.receiver,  sources.use = c(6,7,8), targets.use = c(2,3,4,5))
dev.off()
netVisual_individual(cellchat, signaling ="IL10",  pairLR.use = pairLR.IL10[1,]  ,vertex.receiver = vertex.receiver,  sources.use = c(6,7,8), targets.use = c(2,3,4,5))

netVisual_aggregate(cellchat,idents.use =c(3,4,5,6,7,8,9), signaling = "IL10",  sources.use = c(7,8,9), targets.use = c(3,4,5,6))
netAnalysis_contribution(cellchat, signaling =  "IL10")


netVisual_aggregate(cellchat, signaling = "IL10",  idents.use =c(3,4,5,6,7,8,9), vertex.receiver = vertex.receiver,sources.use = c(7,8,9), targets.use = c(3,4,5,6),layout = c( "hierarchy"))
netVisual_aggregate(cellchat, signaling = "IL10",  idents.use =c(3,4,5,6,7,8,9), vertex.receiver = c(3,4,5,6),layout = c( "hierarchy"))
netVisual_aggregate(cellchat, signaling = "TGFb",  
                    LR.show <- pairLR.CXCL[1,], layout = c( "chord"))

netVisual_individual(cellchat, signaling =  "TGFb", pairLR.use = "TGFb", layout = "chord")
plotGeneExpression(cellchat, signaling = "IL10", enriched.only = FALSE)

netVisual_aggregate(cellchat, signaling = "IL10",layout = c( "hierarchy"),sources.use = c(7,8,9), targets.use = c(3,4,5,6))

netVisual_chord_cell(cellchat,   vertex.receiver = vertex.receiver,signaling ="TGFb",sources.use = c(7,8,9),
                     targets.use = c(3,4,5,6))