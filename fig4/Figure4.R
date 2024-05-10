#====A-B
meta_infor_mouse=readRDS('/work/dy/project/temp/singlecell/hippocampus/mouse_MJ/hippocampus/analysis_d/mouse_seurat_GABA_GLu_other_metainfor_new.rds.gz')
meta_infor_moneky=readRDS('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/macaque_seurat_GABA_GLu_other_metainfor_new.rds.gz')
meta_infor_moneky=meta_infor_moneky[,colnames(meta_infor_mouse)]
meta_infor_Marmaset=readRDS('/work/dy/project/temp/singlecell/hippocampus/Marmaset_MJ/new_MT_data08/Marmaset_seurat_GABA_GLu_other_metainfor_new.rds.gz')

meta_infor_human=readRDS('/work/dy/project/temp/singlecell/hippocampus/merge/multispecice/macaque_tansfer_human_all.data.format.rds.gz')
meta_infor_human=meta_infor_human[,colnames(meta_infor_moneky)]
meta_infor_human$field[which(meta_infor_human$field=='LAMP')]='LAMP5'
				
anno.comb <- Reduce("rbind", list(meta_infor_mouse,meta_infor_moneky,meta_infor_Marmaset[,colnames(meta_infor_moneky)],
					meta_infor_human))					
anno.comb$field=toupper(anno.comb$field)
anno.comb$cluster_gene=toupper(anno.comb$cluster_gene)
anno.comb$donor_id=anno.comb$orig.ident


Neuron=anno.comb %>% 
  dplyr::count(species, Neuronal_or_not) 
Neuron$freq=unlist(lapply(split(Neuron[,"n"],Neuron$species),function(x){x/sum(x)*100}))
P=1
Neuron_plot=ggplot(Neuron,aes(x=species,y=freq))+
    geom_bar(aes(fill=Neuronal_or_not),stat="identity",width=0.7)+
    scale_y_continuous(expand=c(0,0))+
    #scale_fill_manual(limit=c("Benefit","NonBenefit"),values=RColorBrewer::brewer.pal(8,"Set1")[1:2],guide=F)+
    labs(title=paste("\np = ",P,sep=""),y="Percentage (%)")+
    theme(panel.background = element_blank(),axis.line = element_line(color="black"),axis.title.x = element_blank(),axis.text.x = element_text(angle = 60))
ggsave(Neuron_plot, file = paste0(fig.dir, "Neuron.prop.subset.bar.pdf"), width = 4, height = 3)

#============================= GABA
GABA_Glu=anno.comb %>% 
  filter(MainType %in% c('GABA','Glu')) %>% #filter(species %in% c('macaca'))%>% 
  dplyr::count(species, MainType) 
GABA_Glu$freq=unlist(lapply(split(GABA_Glu[,"n"],GABA_Glu$species),function(x){x/sum(x)*100}))

ggplot(subset(GABA_Glu,MainType=='GABA'),aes(x=species,y=freq))+
    geom_bar(aes(fill=MainType),stat="identity",width=0.7)+
    scale_y_continuous(expand=c(0,0))+
    labs(y="Percentage (%)")+
    theme(panel.background = element_blank(),axis.line = element_line(color="black"),axis.title.x = element_blank())

#=subtype
subset_GABA=c('CCK','LAMP5','PVALB','SST','VIP')

GABA_sub=anno.comb %>% 
  filter(MainType %in% c('GABA')) %>% #filter(species %in% c('macaca'))%>% 
  dplyr::count(species, field)
GABA_sub$freq=unlist(lapply(split(GABA_sub[,"n"],GABA_sub$species),function(x){x/sum(x)*100}))
GABA_sub=subset(GABA_sub,field %in% subset_GABA)

ggplot(data = subset(subset(GABA_sub,field %in% c('VIP','PVALB','SST')),species %in% c('human','macaque','marmaset','mouse')
				), mapping = aes(x = species, y = freq, group = field ,color=field)) + geom_line(color='grey')+ geom_point(aes(color=field,shape=field))+ xlab(unique(GABA_sub$field)) +
         theme(axis.text.x = element_text(angle = 60), axis.title.x = element_blank())+ geom_smooth(method='lm')+
			scale_color_manual(values=c('#F18A15','#8964A8','#7B4195'))+
			scale_x_discrete(limits =rev(unique(GABA_sub$species)))

#===C
library(Seurat)
library(cowplot)
#library(feather)
library(dplyr)
library(Matrix)
library(Hmisc)
library(tidyr)
library(stringr)

library(matrixStats)
library(ggplot2)
library(igraph)
set.seed(42)
.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

sample.combined=readRDS('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/new_SST/GABA/seurat_Macaca_Marmaset_Mouse_SCT_intergrat_WuZhongJian_GABA.rds.gz')
DefaultAssay(sample.combined)="SCT"

sample.combined$cell_clusters <- paste(sample.combined$orig.ident_species,sample.combined$field,sample.combined$cluster_gene,sample.combined$cluster_id, sep='_')

png("/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/new_SST/GABA/intergrat_unmap_no_other.png",width = 500,height = 400)
plots1 <- DimPlot(sample.combined, group.by = c("orig.ident_species"), combine = T, cols =.cluster_cols) #两套数据 方式可视化
print(plots1)
dev.off()

#重新设置颜色 mouse 841F5E marmoset E37644 macaque 486BB3  human 96D3D8

#重新设置颜色2 macaque #CC3399 marmoset #FFC000 mouse  #5C74B0 human #57F3D1

.cluster_cols2=c('#96D3D8','#841F5E','#E37644','#486BB3')
#不画人的
.cluster_cols3=c('#841F5E','#E37644','#486BB3')
.cluster_cols4=c('#CC3399','#FFC000','#5C74B0')

Idents(sample.combined)='orig.ident_species'

plots1 <- DimPlot(sample.combined, group.by = c("orig.ident_species"), combine = T, cols =.cluster_cols4,pt.size=0.1,
			order = c('macaque','marmoset','mouse')) #两套数据 方式可视化
print(plots1)

plots3 <- DimPlot(sample.combined, group.by = c("field"), split.by = c("orig.ident_species"), combine = T,cols =.cluster_cols) 
print(plots3)


#===D-E
library(Seurat)
library(dplyr)
library(matrixStats)
library(Matrix)
library(ggplot2)
library(scrattch.hicat)
sample.combined=readRDS('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/GABA/seurat_Macaca_Marmaset_Mouse_SCT_intergrat_WuZhongJian_GABA.rds.gz')

sample.combined <- RunUMAP(sample.combined, dims = 1:100)

Idents(sample.combined) <- sample.combined$orig.ident_species
macaque_data <- subset(sample.combined, idents = "macaque")
marmoset_data <- subset(sample.combined, idents = "marmoset")
mouse_data <- subset(sample.combined, idents = "mouse")

marmoset_data$subclass_label=marmoset_data$field
mouse_data$subclass_label=mouse_data$field
macaque_data$subclass_label=macaque_data$field

Idents(marmoset_data) <- marmoset_data$subclass_label
Idents(mouse_data) <- mouse_data$subclass_label
Idents(macaque_data) <- macaque_data$subclass_label

subclasses <- unique(macaque_data$subclass_label)
marmoset_data <- subset(marmoset_data, idents = subclasses)
mouse_data <- subset(mouse_data, idents = subclasses)

marmoset_data <- SCTransform(marmoset_data,  vars.to.regress = c("percent.mt","nFeature_RNA","nCount_RNA"),  verbose = TRUE)
mouse_data <- SCTransform(mouse_data,  vars.to.regress = c("percent.mt","nFeature_RNA","nCount_RNA"),  verbose = TRUE)
macaque_data <- SCTransform(macaque_data,  vars.to.regress = c("percent.mt","nFeature_RNA","nCount_RNA"),  verbose = TRUE)

mouse_cells_markers <- FindAllMarkers(mouse_data, assay = "SCT", slot = "data")#, test.use = "roc")
marmoset_cells_markers <- FindAllMarkers(marmoset_data, assay = "SCT", slot = "data")#,min.pct = 0.1, logfc.threshold = 0.1)#, test.use = "roc")
macaque_cells_markers <- FindAllMarkers(macaque_data, assay = "SCT", slot = "data")#, test.use = "roc")


subclasses 
tmp <- 1  #change this
marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
macaque_genes <- macaque_cells_markers[grep(subclasses[tmp], macaque_cells_markers$cluster), ]

marmoset_genes <- marmoset_genes[which(marmoset_genes$avg_log2FC > 0), ]
mouse_genes <- mouse_genes[which(mouse_genes$avg_log2FC > 0), ]
macaque_genes <- macaque_genes[which(macaque_genes$avg_log2FC > 0), ]

all.genes <- data.frame(genes = unique(c(macaque_genes$gene, marmoset_genes$gene, mouse_genes$gene)))
all.genes$macaque <- as.character(match(all.genes$genes, macaque_genes$gene))
all.genes$Marmoset <- as.character(match(all.genes$genes, marmoset_genes$gene))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
all.genes$macaque[which(is.na(all.genes$macaque))] <- FALSE 
all.genes$macaque[which(all.genes$macaque != FALSE)] <- TRUE
all.genes$Marmoset[which(is.na(all.genes$Marmoset))] <- FALSE 
all.genes$Marmoset[which(all.genes$Marmoset != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$macaque <- as.logical(all.genes$macaque)
all.genes$Marmoset <- as.logical(all.genes$Marmoset)
all.genes$Mouse <- as.logical(all.genes$Mouse)

apply(all.genes[,-1],2,sum)

library(eulerr)
plot(euler(
  all.genes[ ,2:4]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0(subclasses[tmp], " vs. All Inh"),
  fills = c("royalblue1", "maroon4", "sienna2")
)

subclasses

#heatmap - conserved genes  
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  macaque_genes <- macaque_cells_markers[grep(subclasses[tmp], macaque_cells_markers$cluster), ]
  marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  macaque_genes <- macaque_genes[which(macaque_genes$avg_log2FC > 0), ]
  marmoset_genes <- marmoset_genes[which(marmoset_genes$avg_log2FC > 0), ]
  mouse_genes <- mouse_genes[which(mouse_genes$avg_log2FC > 0), ]
  all.genes <- data.frame(genes = unique(c(macaque_genes$gene, marmoset_genes$gene, mouse_genes$gene)))
  all.genes$macaque <- as.character(match(all.genes$genes, macaque_genes$gene))
  all.genes$Marmoset <- as.character(match(all.genes$genes, marmoset_genes$gene))
  all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
  all.genes$macaque[which(is.na(all.genes$macaque))] <- FALSE 
  all.genes$macaque[which(all.genes$macaque != FALSE)] <- TRUE
  all.genes$Marmoset[which(is.na(all.genes$Marmoset))] <- FALSE 
  all.genes$Marmoset[which(all.genes$Marmoset != FALSE)] <- TRUE
  all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
  all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
  all.genes$macaque <- as.logical(all.genes$macaque)
  all.genes$Marmoset <- as.logical(all.genes$Marmoset)
  all.genes$Mouse <- as.logical(all.genes$Mouse)
  all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,3] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,4] == TRUE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot1 <- genes_to_plot[-1]

data_to_plot <- subset(macaque_data, downsample = 50)  
data_to_plot <- ScaleData(data_to_plot, assay = "SCT")
levels(data_to_plot) <- rev(subclasses)
library(viridis)
converve_macaque_plot=DoHeatmap(data_to_plot, assay = "SCT", slot = "scale.data", features = unique(genes_to_plot)) +
  #scale_fill_viridis(option = "viridis") + NoLegend()
  scale_fill_gradientn(colors = c("white", "white", "sienna2")) + NoLegend() +
  theme(axis.text.y=element_blank()) #+ NoLegend()

data_to_plot <- subset(mouse_data, downsample = 50)  
data_to_plot <- ScaleData(data_to_plot, assay = "SCT")
levels(data_to_plot) <- rev(subclasses)
library(viridis)
converve_mouse_plot=DoHeatmap(data_to_plot, assay = "SCT", slot = "scale.data", features = unique(genes_to_plot)) +
  #scale_fill_viridis(option = "viridis") + NoLegend()
  scale_fill_gradientn(colors = c("white", "white", "sienna2")) + NoLegend() +
  theme(axis.text.y=element_blank()) #+ NoLegend()
  
data_to_plot <- subset(marmoset_data, downsample = 50)  
data_to_plot <- ScaleData(data_to_plot, assay = "SCT")
levels(data_to_plot) <- rev(subclasses)
library(viridis)
converve_marmoset_plot=DoHeatmap(data_to_plot, assay = "SCT", slot = "scale.data", features = unique(genes_to_plot)) +
  #scale_fill_viridis(option = "viridis") + NoLegend()
  scale_fill_gradientn(colors = c("white", "white", "sienna2")) + NoLegend() +
  theme(axis.text.y=element_blank()) #+ NoLegend()
  
library(cowplot)
converve_all_plot = plot_grid(plotlist = list(converve_macaque_plot,converve_marmoset_plot,converve_mouse_plot),nrow = round(sqrt(length(3))))


#heatmap - macaque_genes_specific genes
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  macaque_genes <- macaque_cells_markers[grep(subclasses[tmp], macaque_cells_markers$cluster), ]
  marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  macaque_genes <- macaque_genes[which(macaque_genes$avg_log2FC  > 0), ] #猕猴特异 基因 让他 > 0
  all.genes <- data.frame(genes = unique(c(macaque_genes$gene, marmoset_genes$gene, mouse_genes$gene)))
  all.genes$macaque <- as.character(match(all.genes$genes, macaque_genes$gene))
  all.genes$Marmoset <- as.character(match(all.genes$genes, marmoset_genes$gene))
  all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
  all.genes$macaque[which(is.na(all.genes$macaque))] <- FALSE 
  all.genes$macaque[which(all.genes$macaque != FALSE)] <- TRUE
  all.genes$Marmoset[which(is.na(all.genes$Marmoset))] <- FALSE 
  all.genes$Marmoset[which(all.genes$Marmoset != FALSE)] <- TRUE
  all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
  all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
  all.genes$macaque <- as.logical(all.genes$macaque)
  all.genes$Marmoset <- as.logical(all.genes$Marmoset)
  all.genes$Mouse <- as.logical(all.genes$Mouse)
  all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,3] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,4] == FALSE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot2 <- genes_to_plot[-1]

#heatmap - marmoset_specific genes
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  macaque_genes <- macaque_cells_markers[grep(subclasses[tmp], macaque_cells_markers$cluster), ]
  marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  marmoset_genes <- marmoset_genes[which(marmoset_genes$avg_log2FC > 0), ]
  all.genes <- data.frame(genes = unique(c(macaque_genes$gene, marmoset_genes$gene, mouse_genes$gene)))
  all.genes$macaque <- as.character(match(all.genes$genes, macaque_genes$gene))
  all.genes$Marmoset <- as.character(match(all.genes$genes, marmoset_genes$gene))
  all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
  all.genes$macaque[which(is.na(all.genes$macaque))] <- FALSE 
  all.genes$macaque[which(all.genes$macaque != FALSE)] <- TRUE
  all.genes$Marmoset[which(is.na(all.genes$Marmoset))] <- FALSE 
  all.genes$Marmoset[which(all.genes$Marmoset != FALSE)] <- TRUE
  all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
  all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
  all.genes$macaque <- as.logical(all.genes$macaque)
  all.genes$Marmoset <- as.logical(all.genes$Marmoset)
  all.genes$Mouse <- as.logical(all.genes$Mouse)
  all.genes <- all.genes[which(all.genes[,2] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,3] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,4] == FALSE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot3 <- genes_to_plot[-1]

#heatmap - mouse_specific genes
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  macaque_genes <- macaque_cells_markers[grep(subclasses[tmp], macaque_cells_markers$cluster), ]
  marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  mouse_genes <- mouse_genes[which(mouse_genes$avg_log2FC > 0), ]
  all.genes <- data.frame(genes = unique(c(macaque_genes$gene, marmoset_genes$gene, mouse_genes$gene)))
  all.genes$macaque <- as.character(match(all.genes$genes, macaque_genes$gene))
  all.genes$Marmoset <- as.character(match(all.genes$genes, marmoset_genes$gene))
  all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
  all.genes$macaque[which(is.na(all.genes$macaque))] <- FALSE 
  all.genes$macaque[which(all.genes$macaque != FALSE)] <- TRUE
  all.genes$Marmoset[which(is.na(all.genes$Marmoset))] <- FALSE 
  all.genes$Marmoset[which(all.genes$Marmoset != FALSE)] <- TRUE
  all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
  all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
  all.genes$macaque <- as.logical(all.genes$macaque)
  all.genes$Marmoset <- as.logical(all.genes$Marmoset)
  all.genes$Mouse <- as.logical(all.genes$Mouse)
  all.genes <- all.genes[which(all.genes[,2] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,3] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,4] == TRUE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot4 <- genes_to_plot[-1]

genes_to_plot <- c(genes_to_plot1, genes_to_plot2, genes_to_plot3, genes_to_plot4)
genes_to_plot <- c( genes_to_plot2, genes_to_plot3, genes_to_plot4) #1957=865 458 634

#macaque
data_to_plot <- subset(macaque_data, downsample = 50) #抽样 但是抽样后 热图差异不明显
data_to_plot <- ScaleData(macaque_data, assay = "SCT")
levels(data_to_plot) <- rev(subclasses)

specific_macaque_plot=DoHeatmap(data_to_plot, assay = "SCT", slot = "scale.data", features = unique(genes_to_plot)) +
  scale_fill_gradientn(colors = c("white", "white", "royalblue1")) +
	 theme(axis.text.y=element_blank()) + NoLegend()
	 
#marmoset
data_to_plot <- subset(marmoset_data, downsample = 50)
data_to_plot <- ScaleData(marmoset_data, assay = "SCT")
levels(data_to_plot) <- rev(subclasses)

specific_marmoset_plot=DoHeatmap(data_to_plot, assay = "SCT", slot = "scale.data", features = unique(genes_to_plot)) +
  scale_fill_gradientn(colors = c("white", "white", "sienna2")) +
 theme(axis.text.y=element_blank()) + NoLegend()
 
#mouse
data_to_plot <- subset(mouse_data, downsample = 50)
data_to_plot <- ScaleData(mouse_data, assay = "SCT")
levels(data_to_plot) <- rev(subclasses)

specific_mouse_plot=DoHeatmap(data_to_plot, assay = "SCT", slot = "scale.data", features = unique(genes_to_plot)) +
  scale_fill_gradientn(colors = c("white", "white", "#ae5a8c")) +
 theme(axis.text.y=element_blank()) + NoLegend()

specific_all_plot = plot_grid(plotlist = list(specific_macaque_plot,specific_marmoset_plot,specific_mouse_plot),nrow = round(sqrt(length(3))))


#===F-G
library(Seurat)
library(cowplot)
library(dplyr)
library(Matrix)
library(Hmisc)
library(tidyr)
library(stringr)

library(matrixStats)
library(ggplot2)
library(igraph)

reorder_matrix <- function(matrix1, by.rows = TRUE) {
  if (by.rows == TRUE) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
}

compare_cl <- function(cl, ref.cl,
                       plot.title = NA, plot.silent = TRUE,
                       heat.colors = colorRampPalette(c("white", "grey70", "black"))(100),
                       row.cl.num = min(length(unique(cl)),
                                        length(unique(ref.cl)))) {
  library(grid)
  library(pheatmap)
  
  conf1 <- table(cl, ref.cl) 
  conf1 <- sweep(conf1, 1, rowSums(conf1), "/") 
  conf2 <- reorder_matrix(conf1)
  
  # Cluster co-occurence
  cl.prop.cocl <- apply(conf1, 2, function(x) {
    grid1 <- expand.grid(x, x)
    min.prop <- apply(grid1, 1, min)
  })
  
  cl.prop.cocl.total <- apply(cl.prop.cocl, 1, sum)
  cl.prop.cocl.m <- matrix(cl.prop.cocl.total, nrow(conf1), nrow(conf1),
                           dimnames = list(rownames(conf1), rownames(conf1)))
  
  ph1 <- pheatmap(conf2, cutree_rows = row.cl.num, clustering_method = "ward.D2",
                  # annotation_row = ref.cl.anno[, -grep("cluster_label", colnames(ref.cl.anno))],
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}

our.combined=readRDS('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/new_SST/GABA/seurat_Macaca_Marmaset_Mouse_SCT_intergrat_WuZhongJian_GABA.rds.gz')
our.combined=our.combined@meta.data

our.combined$label_for_heatmap=paste(our.combined$orig.ident_species, our.combined$field,our.combined$cluster_gene,our.combined$cluster_id, sep = "_")
our.combined$label_for_heatmap=paste(our.combined$orig.ident_species, our.combined$field,our.combined$cluster_id, sep = "_")

heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)

ref.cl <- our.combined$label_for_heatmap  
cca.cl <- our.combined$seurat_clusters.new 
compare.species <- unique(our.combined$orig.ident_species)

cl.conf <- compare_cl(our.combined$label_for_heatmap, our.combined$seurat_clusters.new)
cocl <- cl.conf$cocl
grid.draw(cl.conf[["ph"]]$gtable)

library(fpc)
clus.method <- "single"

clus.num <- pamk(cocl, 1:(min(nrow(cocl), ncol(cocl)) - 1))$nc

cocl.subset <- cocl[grepl(compare.species[1], row.names(cocl)),
                    grepl(compare.species[3], row.names(cocl))] 

row.names(cocl.subset) <- sub(paste0(compare.species[1], "_"), "",
                              row.names(cocl.subset))

colnames(cocl.subset) <- sub(paste0(compare.species[3], "_"), "", 
                             colnames(cocl.subset))


cocl.subset2=cocl.subset
temp_macaque_mouse=cocl.subset2

cocl.subset_plot1=pheatmap(cocl.subset2, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat.colors,
         fontsize = 12) 			


cocl.subset <- cocl[grepl(compare.species[1], row.names(cocl)),
                    grepl(compare.species[2], row.names(cocl))] 

row.names(cocl.subset) <- sub(paste0(compare.species[1], "_"), "",
                              row.names(cocl.subset))

colnames(cocl.subset) <- sub(paste0(compare.species[2], "_"), "", 
                             colnames(cocl.subset))

cocl.subset2=cocl.subset
temp_macaque_marmost=cocl.subset2							 

cocl.subset_plot2=pheatmap(cocl.subset2[,c(1,4,3,2,6,7,8,5,9,10:13)], cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat.colors,
         fontsize = 12) 

		 set.seed(42)
library(magrittr)
library(Seurat)
library(reshape)
library(ggplot2)
library(pheatmap)
set.seed(42)
.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

macaque_transfer_big=read.csv('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/macaque1/macaque1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice.csv',header=T)		
macaque_transfer_big=macaque_transfer_big[-which(macaque_transfer_big$region=='NAN'),]
macaque_transfer_big=macaque_transfer_big[grep('GABA',macaque_transfer_big$Cell_Type),]
macaque_transfer_big=table(macaque_transfer_big$Cell_Type,macaque_transfer_big$region)
rownames(macaque_transfer_big)=sub(' ','.',rownames(macaque_transfer_big))
#mouse_transfer_big=mouse_transfer_big[,-1]			
macaque_transfer_big=t(macaque_transfer_big)

macaque_transfer_big_GABA=macaque_transfer_big[sort(rownames(macaque_transfer_big)),]#c("GABA.Cck","GABA.Lamp5","GABA.Pvalb","GABA.Sst","GABA.Vip")


#marmoset_transfer_big=read.csv('/home/share/share_gouxj/transfer_0818_result/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_filter_region.csv',header=T)		
marmoset_transfer_big=read.csv('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice.csv',header=T)		

marmoset_transfer_big=marmoset_transfer_big[-which(marmoset_transfer_big$region=='NAN'),]
marmoset_transfer_big=marmoset_transfer_big[grep('GABA',marmoset_transfer_big$Cell_Type),]

marmoset_transfer_big=table(marmoset_transfer_big$Cell_Type,marmoset_transfer_big$region)
rownames(marmoset_transfer_big)=sub(' ','.',rownames(marmoset_transfer_big))
#mouse_transfer_big=mouse_transfer_big[,-1]			
marmoset_transfer_big=t(marmoset_transfer_big)

marmoset_transfer_big_GABA=marmoset_transfer_big[sort(rownames(marmoset_transfer_big)),]#c("GABA.Cck","GABA.Lamp5","GABA.Pvalb","GABA.Sst","GABA.Vip")


mouse_transfer_big=read.csv('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/mouse1/mouse1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice.csv',header=T)		

mouse_transfer_big=mouse_transfer_big[-which(mouse_transfer_big$region=='NAN'),]

mouse_transfer_big=mouse_transfer_big[grep('GABA',mouse_transfer_big$Cell_Type),]

mouse_transfer_big=table(mouse_transfer_big$Cell_Type,mouse_transfer_big$region)
rownames(mouse_transfer_big)=sub(' ','.',rownames(mouse_transfer_big))
mouse_transfer_big=t(mouse_transfer_big)

mouse_transfer_big_GABA=mouse_transfer_big[sort(rownames(mouse_transfer_big)),]#c("GABA.Cck","GABA.Lamp5","GABA.Pvalb","GABA.Sst","GABA.Vip")

zijulei_color=list('CA2/3 pyr'='#ffe648',
"CA3/4 pyr"='#cbb83a',
"CA L.Mol"='#7f3a46',
"CA1 Orien & Rad"='#994554',
 "CA1 pyr"='#ff748d',
"DG Mol"='#246d7f',
"DG Gr"='#48daff',
"SUB"='#4affb1',
"SUB sup."='#7f4c99',
"SUB int."='#d57eff',
"SUB deep"='#aa65cc',
"other"='#000000',
"HIP alveus"='#808080')

rownames(mouse_transfer_big_GABA)
mouse_color=list(
 "CA2 pyr"='#ffe648', 
"CA3 pyr"='#cbb83a', 
"CA L.Mol"= '#7f3a46', 
# "CA mid":'#994554',
"CA1 pyr"= '#ff748d',
"DG Mol"='#246d7f',
"DG Gr"='#48daff', 
 "DG Po"= '#2b8399',
 "SUB"='#4affb1', 
"CA3 Orien & Rad"="#dad44f",
 "CA1 Orien & Rad"= "#da7184", 
 # "SUB inner":'#7f4c99', 
 # "SUB mid":'#d57eff', 
 "FC"='#33b07a', 
 "other"='#000000', 
 "HIP alveus"='#808080'
)


library(reshape)
macaque_transfer_big_GABA_bar=macaque_transfer_big_GABA[-which(rownames(macaque_transfer_big_GABA)=='other'),]
macaque_transfer_big_GABA_bar=apply(macaque_transfer_big_GABA_bar,2,function(x){x/sum(x)})
macaque_transfer_big_GABA_bar=melt(macaque_transfer_big_GABA_bar)
macaque_transfer_big_GABA_bar$wuzhong='macaque'
colnames(macaque_transfer_big_GABA_bar)=colnames(macaque_transfer_big_GABA_bar)
	
p_macaque <- ggplot(macaque_transfer_big_GABA_bar, aes(x=Var.2 ,y= value))
p3_macaque=p_macaque+
		geom_bar(aes(fill=Var.1),stat = "identity",color="black")+  #stat_smooth(method='lm')+#aes(shape = Var1)
		scale_fill_manual(values = unlist(zijulei_color), limits=names(zijulei_color))+ #facet_grid(wuzhong~.,scale="free_y") +
		theme(axis.text.x = element_text(angle = 60))#+, vjust = 0.5, hjust=1
options(repr.plot.width=8, repr.plot.height=7)
p3_macaque

marmoset_transfer_big_GABA_bar=marmoset_transfer_big_GABA[-which(rownames(marmoset_transfer_big_GABA)=='other'),]
marmoset_transfer_big_GABA_bar=apply(marmoset_transfer_big_GABA_bar,2,function(x){x/sum(x)})
marmoset_transfer_big_GABA_bar=melt(marmoset_transfer_big_GABA_bar)
marmoset_transfer_big_GABA_bar$wuzhong='marmoset'
colnames(marmoset_transfer_big_GABA_bar)=colnames(marmoset_transfer_big_GABA_bar)

p_marmoset <- ggplot(marmoset_transfer_big_GABA_bar, aes(x=Var.2 ,y= value))
p3_marmoset=p_marmoset+
		#geom_point(lwd = 0.8)+  #stat_smooth(method='lm')+#aes(shape = Var1)
		geom_bar(aes(fill=Var.1),stat = "identity",color="black")+  #stat_smooth(method='lm')+#aes(shape = Var1)
		scale_fill_manual(values = unlist(zijulei_color), limits=names(zijulei_color))+ #facet_grid(wuzhong~.,scale="free_y") +
		theme(axis.text.x = element_text(angle = 60))#+, vjust = 0.5, hjust=1
p3_marmoset

mouse_transfer_big_GABA_bar=mouse_transfer_big_GABA[-which(rownames(mouse_transfer_big_GABA)=='other'),]
mouse_transfer_big_GABA_bar=apply(mouse_transfer_big_GABA_bar,2,function(x){x/sum(x)})
mouse_transfer_big_GABA_bar=melt(mouse_transfer_big_GABA_bar)
mouse_transfer_big_GABA_bar$wuzhong='mouse'
mouse_transfer_big_GABA_bar$Var.1=as.character(mouse_transfer_big_GABA_bar$Var.1)

mouse_transfer_big_GABA_bar=subset(mouse_transfer_big_GABA_bar, !(Var.2 %in% c('GABA.1 Cdh7 (Sst)','GABA.8 Plxna4 (Pvalb)')))

p_mouse <- ggplot(mouse_transfer_big_GABA_bar, aes(x=Var.2 ,y= value))
p3_mouse=p_mouse+
		geom_bar(aes(fill=Var.1),stat = "identity",color="black")+  #stat_smooth(method='lm')+#aes(shape = Var1)
		scale_fill_manual(values = unlist(mouse_color), limits=names(mouse_color))+ #facet_grid(wuzhong~.,scale="free_y") +
		theme(axis.text.x = element_text(angle = 60))#+, vjust = 0.5, hjust=1
options(repr.plot.width=8, repr.plot.height=7)
p3_mouse


library(cowplot)
p3_macaque_marmoset_mouse = plot_grid(plotlist = list(p3_macaque,p3_marmoset,p3_mouse),nrow=1)#,nrow = round(sqrt(length(3))))
p3_macaque_marmoset_mouse

#==H-I
library(Seurat)
set.seed(42)

library(magrittr)
library(ggplot2)
#=====================================================macaque
cellinfor=read.csv('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/macaque1/macaque1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell.csv',header=T)		

cellinfor_type=lapply(split(cellinfor,cellinfor$Slice),function(xx){
    length(unique(xx$Cell_Type))
   
})%>% unlist
cellinfor_type=cellinfor_type[cellinfor_type>=30]

#
slice_use=which(names(cellinfor_type)=='T36')

mytype = 'GABA new_SST_1 (SST)'
for(ii1 in c(slice_use)){
    xx=names(cellinfor_type)[ii1]
	T59 <- readRDS(paste('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/macaque1/spatial_rds_filter/',xx,'.rds',sep=''))

			  temp_infor=subset(cellinfor,Slice==xx)
  T59$Cell_Type=temp_infor$Cell_Type[match(T59$cell,temp_infor$cell)]


	T59$all_cell='all_cell_group'
    umap <- new(
		Class = 'DimReduc',
		cell.embeddings = as.matrix(T59@meta.data[,c('x','y')])
        )
    seurat_obj_umap = T59
    seurat_obj_umap@reductions[["umap"]]=umap 
        seurat_obj_umap@reductions$umap@key="UMAP_"
	colnames(seurat_obj_umap@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")

       seurat_obj_change=seurat_obj_umap
    Idents(seurat_obj_change) = "Cell_Type"

  temp.seurat = seurat_obj_change
  metadata = temp.seurat@meta.data
  metadata[,"tempType"] = ifelse(metadata$Cell_Type==mytype,mytype,"Other")
  
  temp.seurat@meta.data = metadata
  Idents(temp.seurat) = "tempType"
  
  colors = c('yellow','navy')
  names(colors) = c(mytype,"Other")
  
  SliceName = xx
  
  p = DimPlot(temp.seurat,
              reduction = "umap",group.by='tempType',
              cols = colors,
        order = c(mytype,'Other') ,
              pt.size = 1,
              label = F) + 
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())+  
    ggtitle(SliceName)+ coord_fixed()+
    NoLegend()+ geom_point(data=data.frame(x=seq(min(temp.seurat@meta.data$x+10),min(temp.seurat@meta.data$x+10+2000),1),y=min(temp.seurat@meta.data$y)+10), 
                  aes(x=x,y=y),color='red',size=1)
  print(p)
}


slice_use=which(names(cellinfor_type)=='T30')

mytype = 'GABA new_SST_0 (SST)'

for(ii1 in c(slice_use)){
    xx=names(cellinfor_type)[ii1]
	T59 <- readRDS(paste('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/macaque1/spatial_rds_filter/',xx,'.rds',sep=''))
	
			  temp_infor=subset(cellinfor,Slice==xx)
  T59$Cell_Type=temp_infor$Cell_Type[match(T59$cell,temp_infor$cell)]

	T59$all_cell='all_cell_group'
    umap <- new(
		Class = 'DimReduc',
		cell.embeddings = as.matrix(T59@meta.data[,c('x','y')])
        )
    seurat_obj_umap = T59
    seurat_obj_umap@reductions[["umap"]]=umap 
        seurat_obj_umap@reductions$umap@key="UMAP_"
	colnames(seurat_obj_umap@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")

       seurat_obj_change=seurat_obj_umap
    Idents(seurat_obj_change) = "Cell_Type"

  temp.seurat = seurat_obj_change
  metadata = temp.seurat@meta.data
  metadata[,"tempType"] = ifelse(metadata$Cell_Type==mytype,mytype,"Other")
  
  temp.seurat@meta.data = metadata
  Idents(temp.seurat) = "tempType"
  
  colors = c('yellow','navy')
  names(colors) = c(mytype,"Other")
  
  SliceName = xx
  
  p = DimPlot(temp.seurat,
              reduction = "umap",group.by='tempType',
              cols = colors,
        order = c(mytype,'Other') ,
              pt.size = 1,
              label = F) + 
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())+  
    ggtitle(SliceName)+ coord_fixed()+
    NoLegend()+ geom_point(data=data.frame(x=seq(min(temp.seurat@meta.data$x+10),min(temp.seurat@meta.data$x+10+2000),1),y=min(temp.seurat@meta.data$y)+10), 
                  aes(x=x,y=y),color='red',size=1)
  print(p)
}


#=====================================================marmoset
cellinfor=read.csv('/home/dongyu/project/hippo/spaceResult/transfer_from_xiaojuan20230926/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell.csv',header=T)
cellinfor$gene_area_old=cellinfor$gene_area
cellinfor$gene_area=cellinfor$mapped

cellinfor_type=lapply(split(cellinfor,cellinfor$Slice),function(xx){
    length(unique(xx$Cell_Type))
   
})%>% unlist

mytype = 'GABA 11 11 (SST)'
 
slice_use=which(names(cellinfor_type)=='T457')

for(ii1 in c(slice_use)){
    xx=names(cellinfor_type)[ii1]
	T59 <- readRDS(paste('/home/dongyu/project/hippo/spaceResult/transfer_from_xiaojuan20230926/marmoset1/spatial_rds_filter/',xx,'.rds',sep=''))
	temp_infor=subset(cellinfor,Slice==xx)
  T59$Cell_Type=temp_infor$Cell_Type[match(T59$cell,temp_infor$cell)]

	T59$all_cell='all_cell_group'
    umap <- new(
		Class = 'DimReduc',
		cell.embeddings = as.matrix(T59@meta.data[,c('x','y')])
        )
    seurat_obj_umap = T59
    seurat_obj_umap@reductions[["umap"]]=umap
        seurat_obj_umap@reductions$umap@key="UMAP_"
	colnames(seurat_obj_umap@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")

       seurat_obj_change=seurat_obj_umap
    Idents(seurat_obj_change) = "Cell_Type"

  temp.seurat = seurat_obj_change
  metadata = temp.seurat@meta.data
  metadata[,"tempType"] = ifelse(metadata$Cell_Type==mytype,mytype,"Other")
  
  temp.seurat@meta.data = metadata
  Idents(temp.seurat) = "tempType"
  
  colors = c('yellow','navy')
  names(colors) = c(mytype,"Other")
  
  SliceName = xx
  
  p = DimPlot(temp.seurat, 
              reduction = "umap",group.by='tempType',
              cols = colors,
	        order = c(mytype,'Other') ,
              pt.size = 2,
              label = F) +
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())+ 
    ggtitle(SliceName)+ coord_fixed()+
    NoLegend()+ 
    geom_point(data=data.frame(x=seq(min(temp.seurat@meta.data$x+10),min(temp.seurat@meta.data$x+10+2000),1),y=min(temp.seurat@meta.data$y)+10), 
                  aes(x=x,y=y),color='red',size=1)
  print(p)
}

mytype = 'GABA 1 1 (SST)'
slice_use=which(names(cellinfor_type)=='T470')

for(ii1 in c(slice_use)){
    xx=names(cellinfor_type)[ii1]
	T59 <- readRDS(paste('/home/dongyu/project/hippo/spaceResult/transfer_from_xiaojuan20230926/marmoset1/spatial_rds_filter/',xx,'.rds',sep=''))
			  temp_infor=subset(cellinfor,Slice==xx)
  T59$Cell_Type=temp_infor$Cell_Type[match(T59$cell,temp_infor$cell)]

	T59$all_cell='all_cell_group'
    umap <- new(
		Class = 'DimReduc',
		cell.embeddings = as.matrix(T59@meta.data[,c('x','y')])
        )
    seurat_obj_umap = T59
    seurat_obj_umap@reductions[["umap"]]=umap 
        seurat_obj_umap@reductions$umap@key="UMAP_"
	colnames(seurat_obj_umap@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")

       seurat_obj_change=seurat_obj_umap
    Idents(seurat_obj_change) = "Cell_Type"

  temp.seurat = seurat_obj_change
  metadata = temp.seurat@meta.data
  metadata[,"tempType"] = ifelse(metadata$Cell_Type==mytype,mytype,"Other")
  
  temp.seurat@meta.data = metadata
  Idents(temp.seurat) = "tempType"
  
  colors = c('yellow','navy')
  names(colors) = c(mytype,"Other")
  
  SliceName = xx
  
  p = DimPlot(temp.seurat, 
              reduction = "umap",group.by='tempType',
              cols = colors,
	        order = c(mytype,'Other') ,
              pt.size = 2,
              label = F) +
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())+ 
    ggtitle(SliceName)+ coord_fixed()+
    NoLegend()+ 
    geom_point(data=data.frame(x=seq(min(temp.seurat@meta.data$x+10),min(temp.seurat@meta.data$x+10+2000),1),y=min(temp.seurat@meta.data$y)+10), 
                  aes(x=x,y=y),color='red',size=1)
  print(p)
}


#=====================================================mouse
library(magrittr)
library(Seurat)
cellinfor=read.csv('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/mouse1/mouse1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell.csv',header=T)

cellinfor_type=lapply(split(cellinfor,cellinfor$Slice),function(xx){
    length(unique(xx$Cell_Type))
   
})%>% unlist

all_rds=list.files('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/mouse1/spatial_rds_filter')
cellinfor_type=cellinfor_type[sub('.rds','',all_rds)]


mytype = 'GABA 11 Cbln4 (Sst)'

slice_use=which(names(cellinfor_type)=='T330')

all_plot=list(slice_use)
for(ii1 in c(slice_use)){
    xx=names(cellinfor_type)[ii1]
	T59 <- readRDS(paste('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/mouse1/spatial_rds_filter/',xx,'.rds',sep=''))

  			  temp_infor=subset(cellinfor,Slice==xx)
  T59$Cell_Type=temp_infor$Cell_Type[match(T59$cell,temp_infor$cell)]

	T59$all_cell='all_cell_group'
    umap <- new(
		Class = 'DimReduc',
		cell.embeddings = as.matrix(T59@meta.data[,c('x','y')])
        )
    seurat_obj_umap = T59
    seurat_obj_umap@reductions[["umap"]]=umap
        seurat_obj_umap@reductions$umap@key="UMAP_"
	colnames(seurat_obj_umap@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")

       seurat_obj_change=seurat_obj_umap
    Idents(seurat_obj_change) = "Cell_Type"

  temp.seurat = seurat_obj_change
  metadata = temp.seurat@meta.data
  metadata[,"tempType"] = ifelse(metadata$Cell_Type==mytype,mytype,"Other")
  
  temp.seurat@meta.data = metadata
  Idents(temp.seurat) = "tempType"
  
   colors = c('yellow','navy')
 names(colors) = c(mytype,"Other")
  
  SliceName = xx
  
  p = DimPlot(temp.seurat, 
              reduction = "umap",group.by='tempType',
              cols = colors,
		        order = c(mytype,'Other') ,
              pt.size = 2.5,
              label = F) +
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())+  
    ggtitle(SliceName)+ coord_fixed()+
    NoLegend()+ geom_point(data=data.frame(x=seq(min(temp.seurat@meta.data$x+10),min(temp.seurat@meta.data$x+10+1000),1),y=min(temp.seurat@meta.data$y)+10), 
                  aes(x=x,y=y),color='red',size=1)
  print(p)

}

slice_use=which(names(cellinfor_type)=='T330')

mytype = 'GABA 10 Pde1c (Sst)'
for(ii1 in c(slice_use)){
    xx=names(cellinfor_type)[ii1]
	T59 <- readRDS(paste('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/mouse1/spatial_rds_filter/',xx,'.rds',sep=''))

  			  temp_infor=subset(cellinfor,Slice==xx)
  T59$Cell_Type=temp_infor$Cell_Type[match(T59$cell,temp_infor$cell)]

	T59$all_cell='all_cell_group'
    umap <- new(
		Class = 'DimReduc',
		cell.embeddings = as.matrix(T59@meta.data[,c('x','y')])
        )
    seurat_obj_umap = T59
    seurat_obj_umap@reductions[["umap"]]=umap
        seurat_obj_umap@reductions$umap@key="UMAP_"
	colnames(seurat_obj_umap@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")

       seurat_obj_change=seurat_obj_umap
    Idents(seurat_obj_change) = "Cell_Type"

  temp.seurat = seurat_obj_change
  metadata = temp.seurat@meta.data
  metadata[,"tempType"] = ifelse(metadata$Cell_Type==mytype,mytype,"Other")
  
  temp.seurat@meta.data = metadata
  Idents(temp.seurat) = "tempType"
  
   colors = c('yellow','navy')
 names(colors) = c(mytype,"Other")
  
  SliceName = xx
  
  p = DimPlot(temp.seurat, 
              reduction = "umap",group.by='tempType',
              cols = colors,
		        order = c(mytype,'Other') ,
              pt.size = 2.5,
              label = F) +
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())+ 
    ggtitle(SliceName)+ coord_fixed()+
    NoLegend()+ geom_point(data=data.frame(x=seq(min(temp.seurat@meta.data$x+10),min(temp.seurat@meta.data$x+10+1000),1),y=min(temp.seurat@meta.data$y)+10), 
                  aes(x=x,y=y),color='red',size=1)
  print(p)

}

library(magrittr)
library(Seurat)

#================================================将 Glu GABA 其他类 合并为一个 seurat
setwd('/home/dongyu/project/hippo/monkey_two_pi/analysis_d/harmoney')
library(ggplot2)
library(Seurat)
library(magrittr)
library(dplyr)

set.seed(42)
.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

#猕猴的 seurat
seurat_harmony_filter_refine_GABA_GLU=readRDS('/home/dongyu/project/hippo/monkey_two_pi/analysis_d/harmoney/Macaque_seurat_GABA_GLu_other_seurat.newSST01.rds.gz')
#Glu 24275
#seurat_harmony_Glu_filter=subset(seurat_harmony_filter_refine_GABA_GLU, subset= MainType=='Glu')


#GABA 4641  
seurat_harmony_GABA_filter=subset(seurat_harmony_filter_refine_GABA_GLU, subset= MainType=='GABA')


#读入 猕猴的 GABA的差异基因
seurat_harmony_GABA_filter_deg=readRDS('/home/dongyu/project/hippo/monkey_two_pi/analysis_d/harmoney/Macaque_seurat_GABA_GLu_other_seurat.newSST01.GABA.marker.gene.rds.gz')
seurat_harmony_GABA_filter_deg_SSt1=subset(seurat_harmony_GABA_filter_deg,cluster=='new_SST_1')
head(seurat_harmony_GABA_filter_deg_SSt1)

#狨猴 的 seurat
marmoset_seurat_harmony_filter_refine_GABA_GLU=readRDS('/home/dongyu/project/hippo/merge/plasticity/marmoset/marmoset_MT/analysis_d/harmoney/Marmaset_seurat_GABA_GLu_other_filter_seurat.rds.gz')
marmoset_seurat_harmony_filter_refine_GABA_GLU_GABA=subset(marmoset_seurat_harmony_filter_refine_GABA_GLU, subset= MainType=='GABA')


#狨猴 的 seurat的 GABA的差异基因
marmoset_seurat_harmony_GABA_filter_deg=readRDS('/home/dongyu/project/hippo/merge/plasticity/marmoset/marmoset_MT/analysis_d/harmoney/GABA/marmoset_GABA_cell_res1_markers_eachCluster-final.rds.gz')
marmoset_seurat_harmony_GABA_filter_deg_SSt1=subset(marmoset_seurat_harmony_GABA_filter_deg,cluster=='11')
head(marmoset_seurat_harmony_GABA_filter_deg_SSt1)

#==猕猴和狨猴 sst1 的差异基因overlap  top 50
seurat_harmony_GABA_filter_deg_SSt1_sort=seurat_harmony_GABA_filter_deg_SSt1[order(seurat_harmony_GABA_filter_deg_SSt1$avg_log2FC,decreasing =T),]
marmoset_seurat_harmony_GABA_filter_deg_SSt1_sort=marmoset_seurat_harmony_GABA_filter_deg_SSt1[order(marmoset_seurat_harmony_GABA_filter_deg_SSt1$avg_log2FC,decreasing =T),]

inter_SST1_overlap=intersect(seurat_harmony_GABA_filter_deg_SSt1_sort[c(1:50),]$gene,
    marmoset_seurat_harmony_GABA_filter_deg_SSt1_sort[c(1:50),]$gene)
inter_SST1_overlap


#分成 2类  画 sst1的marker
marmoset_seurat_harmony_filter_refine_GABA_GLU_GABA$cluster_id_SSt1=ifelse(marmoset_seurat_harmony_filter_refine_GABA_GLU_GABA$cluster_id=='11',
			'new_SST_1','other')
DotPlot(marmoset_seurat_harmony_filter_refine_GABA_GLU_GABA,features=inter_SST1_overlap[c(1:10)],cols =c("#fcf8f2","red"),
			scale=F,group.by='cluster_id_SSt1') +#,
			theme(axis.text.x = element_text(angle = 45, hjust=1))+#group.by='SCT_snn_res.0.2')+ 
			#scale_y_discrete(limits=candidate_maintyp_order)	+
			RotatedAxis()



seurat_harmony_GABA_filter$cluster_id_SSt1=ifelse(seurat_harmony_GABA_filter$cluster_id=='new_SST_1',
			'new_SST_1','other')
DotPlot(seurat_harmony_GABA_filter,features=inter_SST1_overlapc(1:10)],cols =c("#fcf8f2","red"),
			scale=F,group.by='cluster_id_SSt1') +#,
			theme(axis.text.x = element_text(angle = 45, hjust=1))+#group.by='SCT_snn_res.0.2')+ 
			#scale_y_discrete(limits=candidate_maintyp_order)	+
			RotatedAxis()	





















			
			
			
			