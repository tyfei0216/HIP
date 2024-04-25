###Figure 3A
seuratobj <- readRDS("E:/单细胞数据/海马数据/20231015_filter_data/macaque1/spatial_rds_filter/T25.rds")
Idents(seuratobj) <- seuratobj$subtype_rename
table(seuratobj$subtype_rename)
seuratobj <- subset(seuratobj, idents = c("Glu CA3/4 1","Glu CA3/4 2"))


seurat_list <- list()

# 使用 for 循环读入三个 rds 文件
for (i in c(25,27,28,29,30,31,32,33,34,36,37,39,40,41,42,43,44,45,47,48,49,50,51,52,53,54,55,56,57,59)) {
  filename <- paste("E:/单细胞数据/海马数据/20231015_filter_data/macaque1/spatial_rds_filter/T", i, ".rds", sep="")
  seuratobj <- readRDS(filename)
  Idents(seuratobj) <- seuratobj$subtype_rename
  # 将 "Glu CA2/3 1" 的细胞存储到列表中
  seurat_list[[i]] <- subset(seuratobj, idents = c("Glu CA3/4 1","Glu CA3/4 2"))
}

# 合并列表中的 Seurat 对象
merged_seurat <- merge(x=seurat_list[[25]],y=c(seurat_list[[27]],seurat_list[[28]],seurat_list[[29]],seurat_list[[30]],seurat_list[[31]],seurat_list[[32]],seurat_list[[33]],seurat_list[[34]],seurat_list[[36]],seurat_list[[37]],seurat_list[[39]],seurat_list[[40]],seurat_list[[41]],seurat_list[[42]],seurat_list[[43]],seurat_list[[44]],seurat_list[[45]],seurat_list[[47]],seurat_list[[48]],seurat_list[[49]],seurat_list[[50]],seurat_list[[51]],seurat_list[[52]],seurat_list[[53]],seurat_list[[54]],seurat_list[[55]],seurat_list[[56]],seurat_list[[57]],seurat_list[[59]]))



for (i in c(25,27,28,29,30,31,32,33,34,36,37,39,40,41,42,43,44,45,47,48,49,50,51,52,53,54,55,57,59)) {
  filename <- paste("E:/单细胞数据/海马数据/20231015_filter_data/macaque1/spatial_rds_filter/T", i, ".rds", sep="")
  seuratobj <- readRDS(filename)
  Idents(seuratobj) <- seuratobj$subtype_rename
  # 将 "Glu CA2/3 1" 的细胞存储到列表中
  seurat_list[[i]] <- subset(seuratobj, idents = c("Glu CA3/4 1","Glu CA3/4 2"))
}

# 合并列表中的 Seurat 对象
merged_seurat <- merge(x=seurat_list[[25]],y=c(seurat_list[[27]],seurat_list[[28]],seurat_list[[29]],seurat_list[[30]],seurat_list[[31]],seurat_list[[32]],seurat_list[[33]],seurat_list[[34]],seurat_list[[36]],seurat_list[[37]],seurat_list[[39]],seurat_list[[40]],seurat_list[[41]],seurat_list[[42]],seurat_list[[43]],seurat_list[[44]],seurat_list[[45]],seurat_list[[47]],seurat_list[[48]],seurat_list[[49]],seurat_list[[50]],seurat_list[[51]],seurat_list[[52]],seurat_list[[53]],seurat_list[[54]],seurat_list[[55]],seurat_list[[57]],seurat_list[[59]]))

saveRDS(merged_seurat,file = 'spa_macaque_CA3_4.rds')


# 使用 for 循环读入三个 rds 文件
for (i in c(447,448,451,452,454,455,456,457,458,460,461,462,463,464,465,466,467,470,472,473,474,475,476,478,479,480,482)) {
  filename <- paste("E:/单细胞数据/海马数据/20231015_filter_data/marmoset1/spatial_rds_filter/T", i, ".rds", sep="")
  seuratobj <- readRDS(filename)
  Idents(seuratobj) <- seuratobj$subtype_rename
  # 将 "Glu CA2/3 1" 的细胞存储到列表中
  seurat_list[[i]] <- subset(seuratobj, idents = c("Glu CA2/3 1"))
}

# 合并列表中的 Seurat 对象
merged_seurat <- merge(x=seurat_list[[447]],y=c(seurat_list[[448]],seurat_list[[451]],seurat_list[[452]],seurat_list[[454]],seurat_list[[455]],seurat_list[[456]],seurat_list[[457]],seurat_list[[458]],seurat_list[[460]],seurat_list[[461]],seurat_list[[462]],seurat_list[[463]],seurat_list[[464]],seurat_list[[465]],seurat_list[[466]],seurat_list[[467]],seurat_list[[470]],seurat_list[[472]],seurat_list[[473]],seurat_list[[474]],seurat_list[[475]],seurat_list[[476]],seurat_list[[478]],seurat_list[[479]],seurat_list[[480]],seurat_list[[482]]))

saveRDS(merged_seurat,file = 'spa_marmoset_CA2_3.rds')

for (i in c(447,448,451,452,454,455,456,457,458,460,461,462,463,464,465,466,467,470,472,473,474,475,476,478,479,480,482)) {
  filename <- paste("E:/单细胞数据/海马数据/20231015_filter_data/marmoset1/spatial_rds_filter/T", i, ".rds", sep="")
  seuratobj <- readRDS(filename)
  Idents(seuratobj) <- seuratobj$subtype_rename
  # 将 "Glu CA2/3 1" 的细胞存储到列表中
  seurat_list[[i]] <- subset(seuratobj, idents = c("Glu CA3/4 1","Glu CA3/4 2"))
}

# 合并列表中的 Seurat 对象
merged_seurat <- merge(x=seurat_list[[447]],y=c(seurat_list[[448]],seurat_list[[451]],seurat_list[[452]],seurat_list[[454]],seurat_list[[455]],seurat_list[[456]],seurat_list[[457]],seurat_list[[458]],seurat_list[[460]],seurat_list[[461]],seurat_list[[462]],seurat_list[[463]],seurat_list[[464]],seurat_list[[465]],seurat_list[[466]],seurat_list[[467]],seurat_list[[470]],seurat_list[[472]],seurat_list[[473]],seurat_list[[474]],seurat_list[[475]],seurat_list[[476]],seurat_list[[478]],seurat_list[[479]],seurat_list[[480]],seurat_list[[482]]))

saveRDS(merged_seurat,file = 'spa_marmoset_CA3_4.rds')

for (i in c(461,470)) {
  filename <- paste("E:/单细胞数据/海马数据/20231015_filter_data/marmoset1/spatial_rds_filter/T", i, ".rds", sep="")
  seuratobj <- readRDS(filename)
  Idents(seuratobj) <- seuratobj$subtype_rename
  # 将 "Glu CA2/3 1" 的细胞存储到列表中
  seurat_list[[i]] <- subset(seuratobj, idents = c("Glu CA3/4 1","Glu CA3/4 2"))
}

# 合并列表中的 Seurat 对象
merged_seurat <- merge(x=seurat_list[[461]],y=seurat_list[[470]])

saveRDS(merged_seurat,file = 'spa_marmoset_CA3_4.rds')


seuratobj=pancreas.integrated

#define spatial coefficient and num_pc to use


library(Seurat)

setwd("E:/单细胞数据/海马数据/20231015_filter_data")

##macaque-CA3/4
seuratobj <- readRDS("E:/单细胞数据/海马数据/20231015_filter_data/macaque1/spatial_rds_filter/T25.rds")
Idents(seuratobj) <- seuratobj$subtype_rename
table(seuratobj$subtype_rename)
seuratobj <- subset(seuratobj, idents = c("Glu CA3/4 1","Glu CA3/4 2"))
Idents(seuratobj) <- seuratobj$mapped
seuratobj <- subset(seuratobj, idents = c("CA3","CA4"))



spa <- 0.02
num_pc <- 20

# normalize and calculate pca
seuratobj <- SCTransform(object=seuratobj)
seuratobj  <- RunPCA(seuratobj)

# fix pca
cor <- seuratobj@meta.data

seuratobj@reductions$pca@cell.embeddings[,num_pc+1] <- cor[,10]*spa
seuratobj@reductions$pca@cell.embeddings[,num_pc+2] <- cor[,11]*spa

# cluster 
r <- 1
seuratobj <- FindNeighbors(seuratobj,dims = 1:num_pc+2)
seuratobj  <- FindClusters(seuratobj, resolution = r, algorithm = 1)
seuratobj  <- RunUMAP(object = seuratobj, dims = 1:num_pc+2) 

DimPlot(seuratobj)
#Idents(seuratobj)=seuratobj$mapped
#DimPlot(seuratobj)

seuratobj@reductions$spatial <- seuratobj@reductions$umap
seuratobj@reductions$spatial@cell.embeddings[,1] <- cor[,10]
seuratobj@reductions$spatial@cell.embeddings[,2] <- cor[,11]

p1=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p1

Idents(seuratobj)=seuratobj$mapped
p2=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p2

# 重命名所有类群的身份

Idents(seuratobj) <- seuratobj$SCT_snn_res.1

seuratobj <- RenameIdents(object = seuratobj, 
                          
                          "0" = "CA3",
                          "1" = "CA3",
                          "2" = "CA3",
                          "3" = "CA4",
                          "4" = "CA3",
                          "5" = "CA4",
                          "6" = "CA4"
)  
seuratobj$celltype <- Idents(seuratobj)

p3=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p3

p=p2|p1|p3


ggsave(p,filename = 'T25_CA3_4.pdf',width = 15,height = 5)


saveRDS(seuratobj,file = 'T25-new-CA3_4.rds')





##macaque-CA2/3

seuratobj <- readRDS("E:/单细胞数据/海马数据/20231015_filter_data2/macaque1/spatial_rds_filter/T59.rds")
Idents(seuratobj) <- seuratobj$subtype_rename
table(seuratobj$subtype_rename)
seuratobj <- subset(seuratobj, idents = c("Glu CA2/3 1"))
Idents(seuratobj) <- seuratobj$mapped
seuratobj <- subset(seuratobj, idents = c("CA2","CA3"))


spa <- 0.02
num_pc <- 20

# normalize and calculate pca
seuratobj <- SCTransform(object=seuratobj)
seuratobj  <- RunPCA(seuratobj)

# fix pca
cor <- seuratobj@meta.data

seuratobj@reductions$pca@cell.embeddings[,num_pc+1] <- cor[,10]*spa
seuratobj@reductions$pca@cell.embeddings[,num_pc+2] <- cor[,11]*spa

# cluster 
r <- 1
seuratobj <- FindNeighbors(seuratobj,dims = 1:num_pc+2)
seuratobj  <- FindClusters(seuratobj, resolution = r, algorithm = 1)
seuratobj  <- RunUMAP(object = seuratobj, dims = 1:num_pc+2) 

DimPlot(seuratobj)
#Idents(seuratobj)=seuratobj$mapped
#DimPlot(seuratobj)

seuratobj@reductions$spatial <- seuratobj@reductions$umap
seuratobj@reductions$spatial@cell.embeddings[,1] <- cor[,10]
seuratobj@reductions$spatial@cell.embeddings[,2] <- cor[,11]

p1=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p1

Idents(seuratobj)=seuratobj$mapped
p2=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p2

# 重命名所有类群的身份

Idents(seuratobj) <- seuratobj$SCT_snn_res.1

seuratobj <- RenameIdents(object = seuratobj, 
                          
                          "0" = "CA2",
                          "1" = "CA3",
                          "2" = "CA2",
                          "3"="CA3",
                          "4"="CA3",
                          "5"="CA3",
                          "6"="CA3",
                          "7"="CA2")
seuratobj$celltype <- Idents(seuratobj)

p3=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p3

p=p2|p1|p3


ggsave(p,filename = 'T49_CA2_3.pdf',width = 15,height = 5)


saveRDS(seuratobj,file = 'T49-new-CA2_3.rds')



##marmoset-CA3/4

seuratobj <- readRDS("E:/单细胞数据/海马数据/20231015_filter_data/marmoset1/spatial_rds_filter/T470.rds")
Idents(seuratobj) <- seuratobj$subtype_rename
table(seuratobj$subtype_rename)
seuratobj <- subset(seuratobj, idents = c("Glu CA3/4 1","Glu CA3/4 2"))
Idents(seuratobj) <- seuratobj$mapped
seuratobj <- subset(seuratobj, idents = c("CA3","CA4"))



spa <- 0.02
num_pc <- 20

# normalize and calculate pca
seuratobj <- SCTransform(object=seuratobj)
seuratobj  <- RunPCA(seuratobj)

# fix pca
cor <- seuratobj@meta.data

seuratobj@reductions$pca@cell.embeddings[,num_pc+1] <- cor[,37]*spa
seuratobj@reductions$pca@cell.embeddings[,num_pc+2] <- cor[,38]*spa

# cluster 
r <- 1
seuratobj <- FindNeighbors(seuratobj,dims = 1:num_pc+2)
seuratobj  <- FindClusters(seuratobj, resolution = r, algorithm = 1)
seuratobj  <- RunUMAP(object = seuratobj, dims = 1:num_pc+2) 

DimPlot(seuratobj)
#Idents(seuratobj)=seuratobj$mapped
#DimPlot(seuratobj)

seuratobj@reductions$spatial <- seuratobj@reductions$umap
seuratobj@reductions$spatial@cell.embeddings[,1] <- cor[,37]
seuratobj@reductions$spatial@cell.embeddings[,2] <- cor[,38]

p1=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p1

Idents(seuratobj)=seuratobj$mapped
p2=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p2

# 重命名所有类群的身份

Idents(seuratobj) <- seuratobj$SCT_snn_res.1

seuratobj <- RenameIdents(object = seuratobj, 
                          
                          "0" = "CA3",
                          "1" = "CA3",
                          "2" = "CA3",
                          "3" = "CA4",
                          "4" = "CA4",
                          "5" = "CA4",
                          "6" = "CA3",
                          "7" = "CA4"
)  
seuratobj$celltype <- Idents(seuratobj)

p3=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p3

p=p2|p1|p3


ggsave(p,filename = 'T470_CA3_4.pdf',width = 15,height = 5)


saveRDS(seuratobj,file = 'T470-new-CA3_4.rds')



##marmoset-CA2/3

seuratobj <- readRDS("E:/单细胞数据/海马数据/20231015_filter_data2/marmoset1/spatial_rds_filter/T482.rds")
Idents(seuratobj) <- seuratobj$subtype_rename
table(seuratobj$subtype_rename)
seuratobj <- subset(seuratobj, idents = c("Glu CA2/3 1"))
Idents(seuratobj) <- seuratobj$mapped
seuratobj <- subset(seuratobj, idents = c("CA2","CA3"))



spa <- 0.02
num_pc <- 20

# normalize and calculate pca
seuratobj <- SCTransform(object=seuratobj)
seuratobj  <- RunPCA(seuratobj)

# fix pca
cor <- seuratobj@meta.data

seuratobj@reductions$pca@cell.embeddings[,num_pc+1] <- cor[,37]*spa
seuratobj@reductions$pca@cell.embeddings[,num_pc+2] <- cor[,38]*spa

# cluster 
r <- 1
seuratobj <- FindNeighbors(seuratobj,dims = 1:num_pc+2)
seuratobj  <- FindClusters(seuratobj, resolution = r, algorithm = 1)
seuratobj  <- RunUMAP(object = seuratobj, dims = 1:num_pc+2) 

DimPlot(seuratobj)
#Idents(seuratobj)=seuratobj$mapped
#DimPlot(seuratobj)

seuratobj@reductions$spatial <- seuratobj@reductions$umap
seuratobj@reductions$spatial@cell.embeddings[,1] <- cor[,37]
seuratobj@reductions$spatial@cell.embeddings[,2] <- cor[,38]

p1=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p1

Idents(seuratobj)=seuratobj$mapped
p2=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p2

# 重命名所有类群的身份

Idents(seuratobj) <- seuratobj$SCT_snn_res.1

seuratobj <- RenameIdents(object = seuratobj, 
                          
                          "0" = "CA2",
                          "1" = "CA3",
                          "2" = "CA2")  
seuratobj$celltype <- Idents(seuratobj)

p3=DimPlot(object = seuratobj, 
           reduction = "spatial",
           label = TRUE,
           label.size = 4,
           repel = TRUE)
p3

p=p2|p1|p3


ggsave(p,filename = 'T451_CA2_3.pdf',width = 15,height = 5)


saveRDS(seuratobj,file = 'T451-new-CA2_3.rds')



##macaque
seurat_integrated <- AddModuleScore(seurat_integrated,
                                    features = c("ZNF507","CEMIP","MTHFR","DNAH8","CRTAC1"),
                                    ctrl = 100,
                                    name = "CA4_macaque")
VlnPlot(object = seurat_integrated, features=c("CA4_macaque"),pt.size = 0)

ggsave(p,filename = 'macaque_CA4_module.pdf',width = 15,height = 5)

seurat_integrated <- AddModuleScore(seurat_integrated,
                                    features = c("CAMK2N1","COPG2","HAPLN4","NNAT","PCP4"),
                                    ctrl = 100,
                                    name = "CA2_macaque")
VlnPlot(object = seurat_integrated, features=c("CA2_macaque"),pt.size = 0)

ggsave(p,filename = 'macaque_CA2_module.pdf',width = 15,height = 5)


##marmoset

seurat_integrated <- AddModuleScore(seurat_integrated,
                                    features = c("CARTPT","GPR162","SYN1","GRIN1","DGCR2","RPL21","OCRL","B3GALNT1","PCSK1N"),
                                    ctrl = 100,
                                    name = "CA4_marmoset")
VlnPlot(object = seurat_integrated, features=c("CA4_marmoset"),pt.size = 0)

ggsave(p,filename = 'marmoset_CA4_module.pdf',width = 15,height = 5)

seurat_integrated <- AddModuleScore(seurat_integrated,
                                    features = c("CABP1", "COTL1", 'WIPF3'),
                                    ctrl = 100,
                                    name = "CA2_marmoset")

VlnPlot(object = seurat_integrated, features=c("CA2_marmoset"),pt.size = 0)

ggsave(p,filename = 'macaque_CA2_module.pdf',width = 15,height = 5)


##Figure 3B
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)




my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
library(stringr)
library(ggplot2)
library(ggforce)
library(ggalluvial)
library(ggVolcano)


#CA2_CA3_macaque
results<- FindMarkers(seurat_integrated, ident.1 = "CA4", ident.2 = "DG", logfc.threshold = 0)


colnames(results)
results["gene"]=rownames(results)
data <- add_regulate(results, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05) 

write.csv(data,file = "CA4_DG_marmoset-2.csv")


p <- ggvolcano(data, x = "log2FoldChange", y = "padj",pointSize = 1.5,
               label ="gene" ,legend_position	="DR",
               #fills = c("brown","grey","steelblue"),
               #colors = c("brown","grey","steelblue"),
               label_number = 10, output = F,  log2FC_cut = 0.5,
               FDR_cut = 0.05)
p
ggsave(p,filename = 'CA4_DG_marmoset.pdf',width = 6,height = 5)


#CA2_CA3_marmoset
results<- FindMarkers(seurat_integrated, ident.1 = "CA4", ident.2 = "DG", logfc.threshold = 0)


colnames(results)
results["gene"]=rownames(results)
data <- add_regulate(results, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05) 

write.csv(data,file = "CA4_DG_mouse.csv")

p <- ggvolcano(data, x = "log2FoldChange", y = "padj",pointSize = 1.5,
               label ="gene" ,legend_position	="DR",
               #fills = c("brown","grey","steelblue"),
               #colors = c("brown","grey","steelblue"),
               label_number = 10, output = F,  log2FC_cut = 0.5,
               FDR_cut = 0.05)
p
ggsave(p,filename = 'CA4_DG_mouse.pdf',width = 6,height = 5)



#CA2_CA3_mouse
results<- FindMarkers(seurat_integrated, ident.1 = "CA2_mouse", ident.2 = "CA3_mouse", logfc.threshold = 0)


colnames(results)
results["gene"]=rownames(results)
data <- add_regulate(results, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05) 

write.csv(data,file = "CA2_CA3_mouse.csv")

data <- read.csv(file = 'CA2_CA3_mouse.csv',row.names = 1,sep = ',',header = T)

p <- ggvolcano(data, x = "log2FoldChange", y = "padj",pointSize = 1.5,
               label ="gene" ,legend_position	="DR",
               #fills = c("brown","grey","steelblue"),
               #colors = c("brown","grey","steelblue"),
               label_number = 20, output = F,  log2FC_cut = 0.5,
               FDR_cut = 0.05)

ggsave(p,filename = 'CA2_CA3_mouse.pdf',width = 6,height = 5)


##Figure 3C

library(UpSetR)         #Upset图（upset 包，适用样本数 2-7）
library(VennDiagram) 
library(ggsci)
##define the color

cors<-pal_npg("nrc", alpha = 0.6)(10)
# 读取数据文件

marmoset <- read.csv('E:/单细胞数据/海马数据/CA2_CA3_marmoset.csv')
macaque <- read.csv('E:/单细胞数据/海马数据/CA2_CA3_macaque.csv')
mouse <- read.csv('E:/单细胞数据/海马数据/CA2_CA3_mouse.csv')


up1 = marmoset$regulate == 'Up'
data_sgni1= marmoset[up1,]

up2 = macaque$regulate == 'Up'
data_sgni2= macaque[up2,]

up3 = mouse$regulate == 'Up'
data_sgni3= mouse[up3,]


v1 <- data_sgni1$X
v2 <- data_sgni2$X
v3 <- data_sgni3$X


upset_list <- list(v1,v2,v3)   # 制作Upset图搜所需要的列表文件
names(upset_list) <- c('marmoset','macaque','mouse')    # 把列名赋值给列表的key值

#作图
pp1 <- upset(fromList(upset_list),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
             matrix.color = '#E64B3599',
             main.bar.color = '#4DBBD599',
             sets.bar.color = "#00A08799",
             nsets = 100,     # 绘制的最大集合个数
             nintersects = 40, #绘制的最大交集个数，NA则全部绘制
             order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
             keep.order = F, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
             mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
             text.scale = 2 # 文字标签的大小
)

pp1


##Figure 3D
DimPlot(object = seurat_integrated, 
        reduction = "umap",
        label = TRUE,
        label.size = 2,
        pt.size = 0.02,
        repel = TRUE,raster=FALSE
)


##Figure 3E


#CA4_CA3_macaque
Idents(seurat_integrated) <- seurat_integrated$predicted.id


results<- FindMarkers(seurat_integrated, ident.1 = "CA4", ident.2 = "CA3", logfc.threshold = 0)


colnames(results)
results["gene"]=rownames(results)
data <- add_regulate(results, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05) 

write.csv(data,file = "CA4_CA3_macaque.csv")


p <- ggvolcano(data, x = "log2FoldChange", y = "padj",pointSize = 1.5,
               label ="gene" ,legend_position	="DR",
               #fills = c("brown","grey","steelblue"),
               #colors = c("brown","grey","steelblue"),
               label_number = 10, output = F,  log2FC_cut = 0.5,
               FDR_cut = 0.05)

ggsave(p,filename = 'CA4_CA3_macaque.pdf',width = 6,height = 5)


#CA4_CA3_marmoset
results<- FindMarkers(seurat_integrated, ident.1 = "CA4", ident.2 = "CA3", logfc.threshold = 0)


colnames(results)
results["gene"]=rownames(results)
data <- add_regulate(results, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05) 

write.csv(data,file = "CA4_CA3_marmoset.csv")


p <- ggvolcano(data, x = "log2FoldChange", y = "padj",pointSize = 1.5,
               label ="gene" ,legend_position	="DR",
               #fills = c("brown","grey","steelblue"),
               #colors = c("brown","grey","steelblue"),
               label_number = 20, output = F,  log2FC_cut = 0.5,
               FDR_cut = 0.05)

ggsave(p,filename = 'CA4_CA3_marmoset.pdf',width = 5,height = 5)



#CA4_CA3_mouse
results<- FindMarkers(seurat_integrated, ident.1 = "CA4_mouse", ident.2 = "CA3_mouse", logfc.threshold = 0)


colnames(results)
results["gene"]=rownames(results)
data <- add_regulate(results, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05) 

write.csv(data,file = "CA4_CA3_mouse.csv")


p <- ggvolcano(data, x = "log2FoldChange", y = "padj",pointSize = 1.5,
               label ="gene" ,legend_position	="DR",
               #fills = c("brown","grey","steelblue"),
               #colors = c("brown","grey","steelblue"),
               label_number = 10, output = F,  log2FC_cut = 0.5,
               FDR_cut = 0.05)

ggsave(p,filename = 'CA4_CA3_mouse.pdf',width = 5,height = 5)



##Figure 3F
p=VlnPlot(seurat_integrated,features = c("SHISA9","BMPR1B","SCG2","CARTPT","NR2F2","DCN","ARHGAP36","N4BP2","CSF2RB2","EPB41L2","MYO1B","NECAB1","CNTN6","GRM8"), pt.size = 0, adjust = 2,split.by = "celltype",cols = my36colors,stack = T)
ggsave(p,filename = 'new.pdf',width =12,height = 6 )

##Figure 3H
library(Seurat)
library(ggplot2)
library(Seurat)
library(ggplot2)
help(FetchData)
Key(seurat)

Seurat=seurat

a1=seurat@assays$integrated@scale.data
a1=t(a1)
exprs=a1

#exprs <- data.frame(FetchData(object = seurat, vars = a1))#vars 输入基因
colnames(exprs)
exprs$Barcod<-rownames(exprs)
ident<-data.frame()
seurat$celltype_exp
#barcode与聚类信息提取
ident<-data.frame(Barcod=names(seurat$celltype),orig.ident=seurat$celltype)
ident
#通过merge函数，将表达量与聚类号对应起来
c<-merge(exprs,ident,by='Barcod')
#对其进行排序
colnames(c)
c
#c$orig.ident<-factor(c$orig.ident,levels=c(sort(unique(seurat$sample))))
c$orig.ident<-factor(c$orig.ident,levels=c("CA4","CA3"))
c$orig.ident
Idents(seurat)
levels(seurat)
library(ggplot2)
my_comparisons <- list(c("CA4","CA3"))
#noise <- rnorm(n = length(x = data[,c('Cdh17')])) / 100000
#data[,c('Cdh17')] <- data[, c('Cdh17')] + noise
#ggplot(data = c,mapping = aes(x = factor(x = orig.ident),y = Junb)) +geom_violin(scale = "width",adjust =1,trim = TRUE,mapping = aes(fill = factor(x = orig.ident)))+
library(ggviolin)
library(ggpubr)

e=ggviolin(c,x = "orig.ident", y = c(colnames(c)[1]), #基因所在的行
           combine = T,color = "orig.ident",alpha = 0.5,error.plot="none",
           palette = "npg",#颜色设定
           trim=F,
           ylab="Normalized Expression",size=0,
           add = "boxplot", add.params = list(fill = "white"))+
  NoLegend()+
  stat_summary(fun=median, geom="line",  linetype = 2,aes(group=0), color="#E63863", size=0.5)+ 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75))

p=e+stat_compare_means(method = "t.test",
                       #label = "p.signif",##星号设置
                       comparisons = my_comparisons)

p


##Figure 3I
library(UpSetR)         #Upset图（upset 包，适用样本数 2-7）
library(VennDiagram) 
library(ggsci)
##define the color

cors<-pal_npg("nrc", alpha = 0.6)(10)
# 读取数据文件

marmoset <- read.csv('E:/单细胞数据/海马数据/CA4_CA3_marmoset.csv')
macaque <- read.csv('E:/单细胞数据/海马数据/CA4_CA3_macaque.csv')
mouse <- read.csv('E:/单细胞数据/海马数据/CA4_CA3_mouse.csv')


up1 = marmoset$regulate == 'Up'
data_sgni1= marmoset[up1,]

up2 = macaque$regulate == 'Up'
data_sgni2= macaque[up2,]

up3 = mouse$regulate == 'Up'
data_sgni3= mouse[up3,]


v1 <- data_sgni1$X
v2 <- data_sgni2$X
v3 <- data_sgni3$X


upset_list <- list(v1,v2,v3)   # 制作Upset图搜所需要的列表文件
names(upset_list) <- c('marmoset','macaque','mouse')    # 把列名赋值给列表的key值

#作图
pp1 <- upset(fromList(upset_list),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
             matrix.color = '#E64B3599',
             main.bar.color = '#4DBBD599',
             sets.bar.color = "#00A08799",
             nsets = 100,     # 绘制的最大集合个数
             nintersects = 40, #绘制的最大交集个数，NA则全部绘制
             order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
             keep.order = F, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
             mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
             text.scale = 2 # 文字标签的大小
)

pp1



##Figure 3J-K-TF factor
rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)



sce_SCENIC <- open_loom("sce_SCENIC.loom")
sce_SCENIC <- open_loom("sce_SCENIC_macaque.loom")
sce_SCENIC <- open_loom("sce_SCENIC_marmoset.loom")
#eprMat <- get_dgem(sce_SCENIC)#从sce_SCENIC文件提取表达矩阵
#exprMat_log <- log2(eprMat+1) # log处理
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)
regulons
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)



#sink("va_result.txt") 
regulons 
#sink()
#RSS分析，查看细胞类型特异性转录因子
library("Seurat")
Seurat <- readRDS("E:/单细胞数据/海马数据/scRNA/转录因子/output/marmoset_zhengehe.rds")
sub_cells=Seurat
table(Seurat$celltype)
#sub_cells@meta.data$celltype_exp=paste(sub_cells@meta.data$sample,sub_cells@meta.data$celltype,sep="_")
sub_cells@meta.data$cell = sub_cells@meta.data$sample
table(sub_cells$cell)
cellinfo <- sub_cells@meta.data[,c("cell","nFeature_RNA","nCount_RNA","virus_level")]#细胞meta信息
colnames(cellinfo)=c( "cell",'nGene' ,'nUMI',"virus_level")
######计算细胞特异性TF

######计算细胞特异性TF
#cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
cellTypes <- as.data.frame(cellinfo)

selectedResolution <- "cell"
sub_regulonAUC <- regulonAUC
sub_regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])

rss=na.omit(rss) 





selectedResolution <- "cell" # select resolution
# Split the cells by cluster:
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution]) 
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] # 去除extened regulons
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression:
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
pheatmap::pheatmap(regulonActivity_byGroup_Scaled,cluster_rows = T, cluster_cols = T)


ggsave(p,filename = 'mouse-pheatmap.pdf',width = 6,height = 25)
ggsave(p,filename = 'macaque-pheatmap.pdf',width = 6,height = 25)
ggsave(p,filename = 'marmoset-pheatmap.pdf',width = 6,height = 25)
ggsave(p,filename = 'mouse-pheatmap2.pdf',width = 6,height = 10)
ggsave(p,filename = 'macaque-pheatmap2.pdf',width = 6,height = 10)
ggsave(p,filename = 'marmoset-pheatmap2.pdf',width = 6,height = 10)

library(dplyr)

#rss1=DataFrame(rss)
#new_cols = c("W3","W9","W15","W33")
#rss1=rss1[,new_cols]
#rss = as.matrix(rss1)
rss
rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    zThreshold = 1,
    #col.low = '#330066',
    #col.mid = '#66CC66',
    #col.high = '#FFCC33'
  )
p = rssPlot$plot
p
ggsave(p,filename = 'marmoset-dotplot.pdf',width=6.18,height = 10)




#rss特异性TF结果
library(ggrepel)
B_rss<- as.data.frame(rss)
B_rss
#需要作图的细胞类型
celltype <- c("CA1","CA2","CA3","CA4","DG")
rssRanklist <-list()
B_rss
write.csv(B_rss,"mouse-rss.CSV")

for(i in 1:length(celltype)){
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))#提取数据
  colnames(data_rank_plot)<- c("TF","celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  p<- ggplot(data_rank_plot, aes(x=rank, y=celltype))+
    geom_point(size=3, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:6,],
               size=3,color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour ='black',size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    geom_text_repel(data= data_rank_plot[1:6,],
                    aes(label=TF), color="black", size=3,fontface="italic",
                    arrow = arrow(ends="first",length = unit(0.01,"npc")),
                    point.padding = 0.3, segment.color = 'black',
                    segment.size = 0.3, force = 1, max.iter = 3e3)
  rssRanklist[[i]] <- p
}


library(cowplot)
pdf("mouse.pdf",width=12, height=4)
pdf("macaque.pdf",width=12, height=4)
pdf("marmoset.pdf",width=12, height=4)
plot_grid(rssRanklist[[1]],rssRanklist[[2]],rssRanklist[[3]],rssRanklist[[4]],rssRanklist[[5]],ncol=12)
dev.off()


#第二个可视化：将转录因子分析结果与seurat对象结合，可视化类似于seurat！

next_regulonAUC <- regulonAUC[,match(colnames(sub_cells),colnames(regulonAUC))]
dim(next_regulonAUC)

regulon_AUC <- regulonAUC@NAMES
regulon_AUC
sub_cells@meta.data = cbind(sub_cells@meta.data ,t(assay(next_regulonAUC[regulon_AUC,])))
sub_cells@meta.data
#自己选定感兴趣的或者比较重要的转录因子，这里我是随机的
next_regulonAUC[regulon_AUC,]
TF_plot <- c("CUX1(+)","MYCN(+)",
             "NR1I2(+)","TEAD4(−)","USF2(+)","ZMZ1(+)","ZSCAN31(−)")

pdf("TFdotplot-macaque.pdf",width=10, height=4)

DotPlot(sub_cells, features = TF_plot,group.by = "celltype")+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
dev.off()


#展示转录因子平均活性！
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 

regulonActivity_byGroup_Scaled 
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

pdf("TF平均.pdf")
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byGroup_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6),
                                   show_row_names = F)) 


exp <- cor(rss, method= "spearman")

#行列注释
annotation_col = data.frame(
  celltype = c("CA3_macaque","CA2_macaque", "CA4_macaque","CA3_marmoset","CA2_marmoset", "CA4_marmoset","CA3_mouse","CA2_mouse", "CA4_mouse"))

row.names(annotation_col) <- colnames(exp)

annotation_row = data.frame(
  celltype = c("CA3_macaque","CA2_macaque", "CA4_macaque","CA3_marmoset","CA2_marmoset", "CA4_marmoset","CA3_mouse","CA2_mouse", "CA4_mouse"))

row.names(annotation_row) <- colnames(exp)

#做热图

p=pheatmap::pheatmap(exp, annotation_col = annotation_col,
                     annotation_row = annotation_row, 
                     #color = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")),
                     cluster_cols = F, cluster_rows = F)
p
ggsave(p,filename = 'correlation.pdf',width = 6,height = 5)
write.csv(exp,file = 'correlation.csv')


