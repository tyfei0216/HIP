library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(scales)
library(ggsci)
library(pheatmap)
library(dplyr)
library(openxlsx)
library(reshape2)
library(gridExtra)

#Define drawing function
merge_plot <- function(df, group.by=NULL, group.cols=NULL, pt.size=0.5, title='',
                       text.size=20, legend.size=5){
  
  if (is.null(group.by)){
    group.by <- 'group.by'
  }else{
    group.by <- df[,group.by]
  }
  
  if (is.null(group.cols)){
    group.cols <- hue_pal()(length(table(group.by)))
  }else{
    group.cols <- group.cols[1:length(table(group.by))]
  }
  
  p <- ggplot(data = df) +
    geom_point(mapping = aes(x=fx, y=fy, color=group.by), size=pt.size) +
    scale_color_manual(values=group.cols) +
    labs(x=' ', y=' ') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.title=element_blank(), 
          text = element_text(size=text.size))+
    guides(colour = guide_legend(override.aes = list(size=legend.size)))+
    xlim(-4355, 4355)+ylim(-3775, 3775)+ coord_fixed()+
    theme_classic() +
    guides(color= "none") +
    theme_void() 
  
  return(p)
}

Bior_Highlight_point <- function(df, group.by='subtype_rename', highlight_group, highlight_color, 
                                 pt.size=1){
  df$color <- 'no'
  df$color[which(df[,group.by] == highlight_group)] <- 'yes'
  df$color <- factor(df$color, levels = c('no', 'yes'))
  df <- df[order(df$color),]
  p <- ggplot(data = df) +
    geom_point(mapping = aes(x = fx, y = fy, color = color), size = pt.size) +
    scale_color_manual(values = c('navy', highlight_color)) +
    labs(x = NULL, y = NULL, title = NULL) + 
    theme_void() +
    theme(legend.position = 'none',
          text = element_text(size = 20),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"))+ coord_fixed()
  return(p)
}

### Figure A : shows the UMAP diagram of all groups of three species
three_species <- readRDS('E:/seurat_spatial/transfer_0818/singcell/rename_data/SCTIntegration/seurat_Macaca_Marmaset_Mouse_SCT_intergrat.rds.gz')
color <- c('Glu CA1'="#1F77B4FF","Glu CA2"="#FFBB78FF",'Glu CA2/3'="#FFBB78FF","Glu CA3"="#FF7F0EFF",'Glu CA3/4'="#FF7F0EFF", 'Glu DG'="#D62728FF",
           'Glu SUB'="#17BECFFF",'Glu pSUB-deep'="#2CA02CFF", 'Glu pSUB-int'="#9EDAE5FF","Glu HIP"="#9467BDFF",'GABA CCK'="#AEC7E8FF",'GABA LAMP5'="#F7B6D2FF",'GABA PVALB'="#98DF8AFF",
           'GABA SST'="#FF9896FF", 'GABA VIP'="#C49C94FF",'astrocyte'="#C5B0D5FF", 'olig'="#DBDB8DFF", 'OPC'="#8C564BFF", 'VLMC'="#E377C2FF","microglia"="#BCBD22FF")
table(three_species$bigtype_copy)
celltype_levels <- names(color)
Idents(three_species) <- factor(three_species$bigtype_copy,levels=celltype_levels)
table(Idents(three_species))
p <- DimPlot(three_species,cols = color,label = T,label.color = "white",raster=FALSE)+theme_void()+
    theme(legend.position = 'right',
          text = element_text(size = 20))

p <- DimPlot(three_species,cols = color,label = F,label.color = "white")+theme_void()+
    theme(plot.title = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          panel.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white",size = 25),
          legend.title = element_text(color = "white"))
ggsave(p,filename = "./three_species/three_UMAP_legend_black.pdf", width =10,height =8)

### Figure B: Zoom in on the marker gene of the group
new_names <- character(length(names(color)) * 2)
new_names[seq(1, length(new_names), by = 2)] <- names(color)
new_names <- c(new_names,rep("",9))
markers.to.plot <- c("MAN1A1","AMIGO2","RGS14","AKAP13","TRPS1","RFX3","RXFP1","TLE4","RASGEF1B","TSHZ2","CCK","LAMP5","POSTN","SST","VIP","GFAP","MOBP","PDGFRA","IGFBP7","APBB1IP")
p <- DotPlot(three_species, features = rev(markers.to.plot), cols = c("royalblue1", "sienna2", "maroon4"), split.by = "orig.ident_species", assay = "SCT")  +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15,angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank())+
  scale_y_discrete(labels = new_names) + coord_flip()
p


ggsave(p,filename = "./three_species/three_dotplot_legend.pdf", width =10,height =6)

### Figure C: The picture after the slice transfer of the enlarged group
#### mouse
mouse_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/mouse1/mouse1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)

plot_dir <- './mouse/plot_png/'
dir.create(plot_dir, recursive = TRUE)
summary(mouse_Merge_df$fx)
summary(mouse_Merge_df$fy)
# Allcelltype for each slice
Merge_df <- mouse_Merge_df
celltype_levels <- names(mouse_color)
Merge_df$celltype_new <- factor(Merge_df$subtype_rename, levels = celltype_levels)
Slices <- names(table(Merge_df$Slice))
for (i in 1:length(Slices)){
  plot_df <- Merge_df[which(Merge_df$Slice==Slices[i]),]
  # Maintype_merge
  p1 <- merge_plot(plot_df, group.by='celltype_new', group.cols=mouse_color, pt.size=0.1, title=Slices[i])
  p1_file <- paste(plot_dir, Slices[i], '_Allcelltype.png', sep = '')
  ggsave(p1, filename = p1_file, width = 3, height = 3, bg = 'transparent')
}

#marmoset
marmoset_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
plot_dir <- './marmoset/plot_png/'
dir.create(plot_dir, recursive = TRUE)
summary(marmoset_Merge_df$fx)
summary(marmoset_Merge_df$fy)
# Allcelltype for each slice
Merge_df <- marmoset_Merge_df
celltype_levels <- names(marmoset_color)
Merge_df$celltype_new <- factor(Merge_df$subtype_rename, levels = celltype_levels)
Slices <- names(table(Merge_df$Slice))
for (i in 1:length(Slices)){
  plot_df <- Merge_df[which(Merge_df$Slice==Slices[i]),]
  # Maintype_merge
  p1 <- merge_plot(plot_df, group.by='celltype_new', group.cols=marmoset_color, pt.size=0.1, title=Slices[i])
  p1_file <- paste(plot_dir, Slices[i], '_Allcelltype.png', sep = '')
  ggsave(p1, filename = p1_file, width = 3, height = 3, bg = 'transparent')
}

#macaque
macaque_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque1/macaque1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
plot_dir <- './macaque/plot_png/'
dir.create(plot_dir, recursive = TRUE)
summary(macaque_Merge_df$fx)
summary(macaque_Merge_df$fy)
# Allcelltype for each slice
Merge_df <- macaque_Merge_df
celltype_levels <- names(macaque_color)
Merge_df$celltype_new <- factor(Merge_df$subtype_rename, levels = celltype_levels)
Slices <- names(table(Merge_df$Slice))
for (i in 1:length(Slices)){
  plot_df <- Merge_df[which(Merge_df$Slice==Slices[i]),]
  # Maintype_merge
  p1 <- merge_plot(plot_df, group.by='celltype_new', group.cols=macaque_color, pt.size=0.1, title=Slices[i])
  p1_file <- paste(plot_dir, Slices[i], '_Allcelltype.png', sep = '')
  ggsave(p1, filename = p1_file, width = 3, height = 3, bg = 'transparent')
}

#Figure D Stereo-Seq spatial distribution diagram of three marmoset genes
#Image path: E:/seurat_spatial/newtransfer/gene_fty/marmoset_3genes_T461_TRPS1.pdf
#Figure E: FISH results of three genes of marmoset
#Picture path: D:/FISH/FISH result picture

### Figure F: Distribution of small groups in spatial self-clustering
#### mouse
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/mouse1/mouse1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)

all_Merge_df$subtype_rename2 <- gsub("^(\\S+)\\s(\\S+)\\s(\\S+)", "\\1 \\2-\\3", all_Merge_df$subtype_rename)
Slices <- names(table(all_Merge_df$Slice))
new_mouse_color <- mouse_color
names(new_mouse_color) <- gsub("^(\\S+)\\s(\\S+)\\s(\\S+)", "\\1 \\2-\\3", names(mouse_color))
mouse_new_celltype_levels <- names(new_mouse_color)
glu_mouse_new_celltype_levels <- mouse_new_celltype_levels[grep("Glu", mouse_new_celltype_levels)]
all_Merge_df$subtype_rename2 <- factor(all_Merge_df$subtype_rename2,levels =glu_mouse_new_celltype_levels)
all_Merge_df <- all_Merge_df[!all_Merge_df$region %in% c("HIP alveus","other","CA L.Mol","NAN","FC"),]
mouse_df <- table(all_Merge_df$region,all_Merge_df$subtype_rename2)
for (j in 1:ncol(mouse_df)){
  mouse_df[,j] <- mouse_df[,j]/colSums(mouse_df)[j]
}
mouse_df <- mouse_df[,complete.cases(t(mouse_df))]
colnames(mouse_df) <- gsub("Glu", "", colnames(mouse_df))
breaks <- seq(0.1,1,0.01)
p1 <- pheatmap(mouse_df, cluster_rows = F, cluster_cols = F, 
               breaks = breaks, legend_breaks = seq(0,1,0.2),
               display_numbers=F,number_format = "%.2f",
               fontsize = 15,border_color = "#696969",angle_col = 45,
               color = colorRampPalette(c("navy","white", "red"))(length(breaks)),fontsize_row=18,fontsize_col=18,legend = F)
p1

#### marmoset
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)

all_Merge_df$subtype_rename2 <- gsub("^(\\S+)\\s(\\S+)\\s(\\S+)", "\\1 \\2-\\3", all_Merge_df$subtype_rename)
Slices <- names(table(all_Merge_df$Slice))
new_marmoset_color <- marmoset_color
marmoset_new_celltype_levels <- names(new_marmoset_color)
glu_marmoset_new_celltype_levels <- marmoset_new_celltype_levels[grep("Glu", marmoset_new_celltype_levels)]
all_Merge_df$subtype_rename2 <- factor(all_Merge_df$subtype_rename2,levels =glu_marmoset_new_celltype_levels)
all_Merge_df <- all_Merge_df[!all_Merge_df$region %in% c("HIP alveus","other","CA L.Mol","NAN"),]
all_Merge_df$region <- factor(all_Merge_df$region,levels = c("CA1-ori/rad","CA1-pyr","CA2/3-pyr","CA3/4-pyr","DG-gr","DG mol","SUB","pSUB-deep","pSUB-int","pSUB-sup"))
marmoset_df <- table(all_Merge_df$region,all_Merge_df$subtype_rename2)
for (j in 1:ncol(marmoset_df)){
  marmoset_df[,j] <- marmoset_df[,j]/colSums(marmoset_df)[j]
}
marmoset_df <- marmoset_df[,complete.cases(t(marmoset_df))]
colnames(marmoset_df) <- gsub("Glu", "", colnames(marmoset_df))
breaks <- seq(0.1,1,0.01)
p2 <- pheatmap(marmoset_df, cluster_rows = F, cluster_cols = F, 
               breaks = breaks, legend_breaks = seq(0,1,0.2),
               display_numbers=F,number_format = "%.2f",
               fontsize = 15,border_color = "#696969",angle_col = 45,
               color = colorRampPalette(c("navy","white", "red"))(length(breaks)),fontsize_row=18,fontsize_col=18,legend = F)
p2

#### macaque
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque1/macaque1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)

all_Merge_df$subtype_rename2 <- gsub("^(\\S+)\\s(\\S+)\\s(\\S+)", "\\1 \\2-\\3", all_Merge_df$subtype_rename)
table(all_Merge_df$subtype_rename2)

Slices <- names(table(all_Merge_df$Slice))
new_macaque_color <- macaque_color
macaque_new_celltype_levels <- names(new_macaque_color)
glu_macaque_new_celltype_levels <- macaque_new_celltype_levels[grep("Glu", macaque_new_celltype_levels)]
all_Merge_df$subtype_rename2 <- factor(all_Merge_df$subtype_rename2,levels =glu_macaque_new_celltype_levels)
all_Merge_df <- all_Merge_df[!all_Merge_df$region %in% c("HIP alveus","other","CA L.Mol","NAN"),]
all_Merge_df$region <- factor(all_Merge_df$region,levels = c("CA1-ori/rad","CA1-pyr","CA2/3-pyr","CA3/4-pyr","DG-gr","DG mol","SUB","pSUB-deep","pSUB-int","pSUB-sup"))
all_Merge_df$region <- gsub("-", "\u002D", all_Merge_df$region)
all_Merge_df$subtype_rename2 <- gsub("-", "\u002D", all_Merge_df$subtype_rename2)
macaque_df <- table(all_Merge_df$region,all_Merge_df$subtype_rename2)
for (j in 1:ncol(macaque_df)){
  macaque_df[,j] <- macaque_df[,j]/colSums(macaque_df)[j]
}
macaque_df <- macaque_df[,complete.cases(t(macaque_df))]
colnames(macaque_df) <- gsub("Glu", "", colnames(macaque_df))
breaks <- seq(0.1,1,0.01)
p3 <- pheatmap(macaque_df, cluster_rows = F, cluster_cols = F, 
               breaks = breaks, legend_breaks = seq(0,1,0.2),
               display_numbers=F,number_format = "%.2f",
               fontsize = 15,border_color = "#696969",angle_col = 45,
               color = colorRampPalette(c("navy","white", "red"))(length(breaks)),fontsize_row=18,fontsize_col=18,legend = F)
p3
ggsave(p3, filename = "./macaque/macaque_spatialdomian_percent.pdf",width =8, height = 5)
p <- cowplot::plot_grid(p3$gtable, p2$gtable, p1$gtable,ncol= 3, align = "h") 
p <- p + theme(plot.margin = margin(25, 0, 25, 0))
p
ggsave(p, filename = "./three_species/three_spatialdomian_percent.pdf",width =24, height = 5)

### Figure G：Distribution of specific SUB small groups and gene expression
#### marmoset
plot_dir <- './marmoset/plot_celltype_split_T470/' 
dir.create(plot_dir, recursive = TRUE)
Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
Merge_df <- Merge_df[Merge_df$Slice == "T470",]

group <- names(table(Merge_df$subtype_rename))
celltype_levels <- names(marmoset_color)
# Allcelltype for each slice
Merge_df$subtype_rename <- factor(Merge_df$subtype_rename, levels = celltype_levels)
for (i in 1:length(group)){
  # Maintype_split
  p2 <- Bior_Highlight_point(Merge_df, highlight_group=group[i],highlight_color="yellow")
  if (group[i]=="Glu CA2/3 1"){
    p2_file <- paste(plot_dir, "Glu CA2_3_1", '.pdf', sep = '')
  }else if (group[i]=="Glu CA3/4 1") {
    p2_file<- paste(plot_dir, "Glu CA3_4_1", '.pdf', sep = '')
  }else if (group[i]=="Glu CA3/4 2") {
    p2_file<- paste(plot_dir, "Glu CA3_4_2", '.pdf', sep = '')
  }else{
    p2_file <- paste(plot_dir, group[i], ' Allslice.pdf', sep = '')
  }
  ggsave(p2, filename = p2_file, width =5,height =4)
}
plot_dir <- './marmoset/plot_gene_ex_T470/' 
dir.create(plot_dir, recursive = TRUE)

celltype_marker <- c("ADAMTSL1","HS3ST4","DPP10","LRRTM4","TSHZ2","RORB")
file <- paste0("E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset1/spatial_rds_filter/T470.rds")
stereoobj <- function(file){
  stereo_obj <- readRDS(file)
  mat_loc <- stereo_obj@meta.data[,c('fx','fy')]
  rownames(mat_loc) <- stereo_obj@meta.data$cell
  umap <- new(
    Class = 'DimReduc',
    cell.embeddings = as.matrix(mat_loc))
  stereo_obj@reductions[["umap"]]=umap
  stereo_obj@reductions$umap@key="UMAP_"
  colnames(stereo_obj@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  return(stereo_obj)
}
T470_stereo <- stereoobj(file)
for (i in 1:length(celltype_marker)){
  p <-FeaturePlot(T470_stereo,features=celltype_marker[i],pt.size=1,slot="data",cols = c("navy","red","yellow"),order=T) +
    labs(x = NULL, y = NULL, title = NULL) + 
    theme_void() + 
    theme(legend.position = 'none',
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"))+ coord_fixed()
  ggsave(p,filename = paste("./marmoset/plot_gene_ex_T470/",celltype_marker[i],".pdf",sep = ""), width =5,height =4)
}
  p <-FeaturePlot(T470_stereo,features=celltype_marker[1],pt.size=1,slot="data",cols = c("navy","red","yellow"),order=T) +
    labs(x = NULL, y = NULL, title = NULL) + 
    theme_void() + 
    theme(legend.position = 'right',
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white"))+ coord_fixed()
  ggsave(p,filename = paste("./marmoset/plot_gene_ex_T470/",celltype_marker[1],"_legend.pdf",sep = ""), width =5,height =4)
#### macaque
plot_dir <- './macaque/plot_celltype_split_T33/' 
dir.create(plot_dir, recursive = TRUE)
Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque1/macaque1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
Merge_df <- Merge_df[Merge_df$Slice == "T33",]

group <- names(table(Merge_df$subtype_rename))
celltype_levels <- names(macaque_color)
# Allcelltype for each slice
Merge_df$subtype_rename <- factor(Merge_df$subtype_rename, levels = celltype_levels)
for (i in 1:length(group)){
  # Maintype_split
  p2 <- Bior_Highlight_point(Merge_df, highlight_group=group[i],highlight_color="yellow")
  if (group[i]=="Glu CA2/3 1"){
    p2_file <- paste(plot_dir, "Glu CA2_3_1", '.pdf', sep = '')
  }else if (group[i]=="Glu CA3/4 1") {
    p2_file<- paste(plot_dir, "Glu CA3_4_1", '.pdf', sep = '')
  }else if (group[i]=="Glu CA3/4 2") {
    p2_file<- paste(plot_dir, "Glu CA3_4_2", '.pdf', sep = '')
  }else{
    p2_file <- paste(plot_dir, group[i], ' Allslice.pdf', sep = '')
  }
  ggsave(p2, filename = p2_file, width =5,height =4)
}
plot_dir <- './macaque/plot_gene_ex_T33/' 
dir.create(plot_dir, recursive = TRUE)
celltype_marker <- c("NDST4","HTR2C","RXFP1","KCNIP4","TSHZ2","ABI3BP")
file <- paste0("E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque1/spatial_rds_filter/T33.rds")
stereoobj <- function(file){
  stereo_obj <- readRDS(file)
  mat_loc <- stereo_obj@meta.data[,c('fx','fy')]
  rownames(mat_loc) <- stereo_obj@meta.data$cell
  umap <- new(
    Class = 'DimReduc',
    cell.embeddings = as.matrix(mat_loc))
  stereo_obj@reductions[["umap"]]=umap
  stereo_obj@reductions$umap@key="UMAP_"
  colnames(stereo_obj@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  return(stereo_obj)
}
T33_stereo <- stereoobj(file)
for (i in 1:length(celltype_marker)){
  p <-FeaturePlot(T33_stereo,features=celltype_marker[i],pt.size=1,slot="data",cols = c("navy","red","yellow"),order=T) +
    labs(x = NULL, y = NULL, title = NULL) + 
    theme_void() + 
    theme(legend.position = 'none',
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"))+ coord_fixed()
  ggsave(p,filename = paste("./macaque/plot_gene_ex_T33/",celltype_marker[i],".pdf",sep = ""), width =5,height =4)
}
  p <-FeaturePlot(T33_stereo,features=celltype_marker[1],pt.size=1,slot="data",cols = c("navy","red","yellow"),order=T) +
    labs(x = NULL, y = NULL, title = NULL) + 
    theme_void() + 
    theme(legend.position = 'right',
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white"))+ coord_fixed()
  ggsave(p,filename = paste("./macaque/plot_gene_ex_T33/",celltype_marker[1],"_legend.pdf",sep = ""), width =5,height =4)

## Attached Picture

### Figure A：UAMP diagram of three small species groups
#### mouse

mouse_celltype_levels <- mouse_new_celltype_levels
mouse_data1$subtype_rename_clean <- mouse_data1$subtype_rename
mouse_data1$subtype_rename_clean <- gsub("Glu ", "", mouse_data1$subtype_rename_clean)
mouse_data1$subtype_rename_clean <- gsub("GABA ", "", mouse_data1$subtype_rename_clean)
mouse_celltype_levels_clean <- mouse_celltype_levels
mouse_celltype_levels_clean <- gsub("Glu ", "", mouse_celltype_levels_clean)
mouse_celltype_levels_clean <- gsub("GABA ", "", mouse_celltype_levels_clean)

Idents(mouse_data1) <- factor(mouse_data1$subtype_rename_clean,levels=mouse_celltype_levels_clean)
generate_gradient <- function(base_color, n_steps) {
  generate_gradient_func <- colorRampPalette(c(base_color, "white"))
  return(generate_gradient_func(n_steps))
}
mouse_list <- lapply(mouse_color, function(base_color) {
  generate_gradient(base_color, n_steps = 10)
})

mouse_gradient_color <- c("Glu CA1 1"=mouse_list$`Glu CA1 1`[1], "Glu CA1 2"=mouse_list$`Glu CA1 1`[3],"Glu CA1 3"=mouse_list$`Glu CA1 1`[5], "Glu CA1 4"=mouse_list$`Glu CA1 1`[7], 
                 "Glu CA2 1"="#FFBB78FF", 
                 "Glu CA3 1"=mouse_list$`Glu CA3 1`[1], "Glu CA3 2"=mouse_list$`Glu CA3 1`[3], "Glu CA3 3"=mouse_list$`Glu CA3 1`[5],"Glu CA3 4"=mouse_list$`Glu CA3 1`[7],
                 "Glu DG 1"=mouse_list$`Glu DG 1`[1],"Glu DG 2"=mouse_list$`Glu DG 1`[3],"Glu DG 3"=mouse_list$`Glu DG 1`[5],"Glu DG 4"=mouse_list$`Glu DG 1`[7],       
                 "Glu HIP 1"=mouse_list$`Glu HIP 1`[1], "Glu HIP 2" =mouse_list$`Glu HIP 1`[2], "Glu HIP 3"=mouse_list$`Glu HIP 1`[3], "Glu HIP 4"=mouse_list$`Glu HIP 1`[4], "Glu HIP 5" =mouse_list$`Glu HIP 1`[5], "Glu HIP 6"=mouse_list$`Glu HIP 1`[6],
                 "Glu SUB 1"="#17BECFFF",
                 "GABA Cck 1" =mouse_list$`GABA Cck 1`[1], "GABA Cck 2"=mouse_list$`GABA Cck 1`[3],
                 "GABA Lamp5 1"=mouse_list$`GABA Lamp5 1`[1], "GABA Lamp5 2"=mouse_list$`GABA Lamp5 1`[3], "GABA Lamp5 3"=mouse_list$`GABA Lamp5 1`[5], 
                 "GABA Pvalb 1"=mouse_list$`GABA Pvalb 1`[1], "GABA Pvalb 2" =mouse_list$`GABA Pvalb 1`[3],"GABA Pvalb 3"=mouse_list$`GABA Pvalb 1`[5],
                 "GABA Sst 1"=mouse_list$`GABA Sst 1`[1],"GABA Sst 2"=mouse_list$`GABA Sst 1`[3],"GABA Sst 3"=mouse_list$`GABA Sst 1`[5],"GABA Sst 4" =mouse_list$`GABA Sst 1`[7],
                 "GABA Vip 1" =mouse_list$`GABA Vip 1`[1],"GABA Vip 2"=mouse_list$`GABA Vip 1`[3],
                 "astrocyte"="#C5B0D5FF", "OPC"="#8C564BFF", "VLMC"="#E377C2FF","microglia"="#BCBD22FF")
names(mouse_gradient_color) <- gsub("^(\\S+)\\s(\\S+)\\s(\\S+)", "\\1 \\2-\\3", names(mouse_gradient_color))
cleaned_labels <- gsub("Glu ", "", names(mouse_gradient_color))
cleaned_labels <- gsub("GABA ", "", cleaned_labels)
cleaned_mouse_gradient_color <- setNames(mouse_gradient_color, cleaned_labels)
DimPlot(mouse_data1, cols = cleaned_mouse_gradient_color, label = FALSE, label.color = "white") +
  theme_void() +
  theme(legend.text = element_text(size = 25))


ggsave("./mouse/mouse_UMAP_fu1a.pdf",width = 10.5, height = 7)

#### marmoset
marmoset_celltype_levels <- marmoset_new_celltype_levels
marmoset_data1$subtype_rename_clean <- marmoset_data1$subtype_rename
marmoset_data1$subtype_rename_clean <- gsub("Glu ", "", marmoset_data1$subtype_rename_clean)
marmoset_data1$subtype_rename_clean <- gsub("GABA ", "", marmoset_data1$subtype_rename_clean)

marmoset_celltype_levels_clean <- marmoset_celltype_levels
marmoset_celltype_levels_clean <- gsub("Glu ", "", marmoset_celltype_levels_clean)
marmoset_celltype_levels_clean <- gsub("GABA ", "", marmoset_celltype_levels_clean)

Idents(marmoset_data1) <- factor(marmoset_data1$subtype_rename_clean,levels=marmoset_celltype_levels_clean)
marmoset_list <- lapply(marmoset_color, function(base_color) {
  generate_gradient(base_color, n_steps = 10)
})

marmoset_gradient_color <- c('Glu CA1 1'="#1F77B4FF", 
                    'Glu CA2/3 1'="#FFBB78FF",
                    'Glu CA3/4 1'=marmoset_list$`Glu CA3/4 1`[1], 'Glu CA3/4 2'=marmoset_list$`Glu CA3/4 1`[3], 
                    'Glu DG 1'="#D62728FF",
                    'Glu SUB 1'="#17BECFFF",
                    'Glu SUB deep 1'=marmoset_list$`Glu SUB deep 1`[1],  'Glu SUB deep 2'=marmoset_list$`Glu SUB deep 1`[3],'Glu SUB deep 3'=marmoset_list$`Glu SUB deep 1`[5], 
                    'Glu SUB int. 1'=marmoset_list$`Glu SUB int. 1`[1],'Glu SUB int. 2'=marmoset_list$`Glu SUB int. 1`[2], 'Glu SUB int. 3'=marmoset_list$`Glu SUB int. 1`[3], 'Glu SUB int. 4'=marmoset_list$`Glu SUB int. 1`[4], 'Glu SUB int. 5'=marmoset_list$`Glu SUB int. 1`[5], 'Glu SUB int. 6'=marmoset_list$`Glu SUB int. 1`[6], 'Glu SUB int. 7'=marmoset_list$`Glu SUB int. 1`[7],'Glu SUB int. 8'=marmoset_list$`Glu SUB int. 1`[8], 
                    'GABA CCK 1'="#AEC7E8FF",
                    'GABA LAMP5 1'=marmoset_list$`GABA LAMP5 1`[1], 'GABA LAMP5 2'=marmoset_list$`GABA LAMP5 1`[1], 'GABA LAMP5 3'=marmoset_list$`GABA LAMP5 1`[1], 
                    'GABA PVALB 1'=marmoset_list$`GABA PVALB 1`[1],'GABA PVALB 2'=marmoset_list$`GABA PVALB 1`[3], 
                    'GABA SST 1'=marmoset_list$`GABA SST 1`[1], 'GABA SST 2'=marmoset_list$`GABA SST 1`[3], 'GABA SST 3'=marmoset_list$`GABA SST 1`[5], 'GABA SST 4'=marmoset_list$`GABA SST 1`[7], 
                    'GABA Vip 1'=marmoset_list$`GABA Vip 1`[1],'GABA Vip 2'=marmoset_list$`GABA Vip 1`[3],'GABA Vip 3'=marmoset_list$`GABA Vip 1`[5],
                    'astrocyte'="#C5B0D5FF", 'olig'="#DBDB8DFF", 'OPC'="#8C564BFF", 'VLMC'="#E377C2FF","microglia"="#BCBD22FF")

names(marmoset_gradient_color) <- gsub("^(\\S+)\\s(\\S+)\\s(\\S+)", "\\1 \\2-\\3", names(marmoset_color))
# Remove substrings from labels
cleaned_labels <- gsub("Glu ", "", names(marmoset_gradient_color))
cleaned_labels <- gsub("GABA ", "", cleaned_labels)
# Create a new color vector with cleaned labels
cleaned_marmoset_gradient_color <- setNames(marmoset_gradient_color, cleaned_labels)
DimPlot(marmoset_data1, cols = cleaned_marmoset_gradient_color, label = FALSE, label.color = "white") +
  theme_void() +
  theme(legend.text = element_text(size = 25))

ggsave("./marmoset/marmoset_UMAP_fu1a.pdf",width = 10.5, height = 7)

#### macaque
macaque_celltype_levels <- macaque_new_celltype_levels
macaque_data1$subtype_rename_clean <- macaque_data1$subtype_rename
macaque_data1$subtype_rename_clean <- gsub("Glu ", "", macaque_data1$subtype_rename_clean)
macaque_data1$subtype_rename_clean <- gsub("GABA ", "", macaque_data1$subtype_rename_clean)

macaque_celltype_levels_clean <- macaque_celltype_levels
macaque_celltype_levels_clean <- gsub("Glu ", "", macaque_celltype_levels_clean)
macaque_celltype_levels_clean <- gsub("GABA ", "", macaque_celltype_levels_clean)

Idents(macaque_data1) <- factor(macaque_data1$subtype_rename_clean,levels=macaque_celltype_levels_clean)
macaque_list <- lapply(macaque_color, function(base_color) {
  generate_gradient(base_color, n_steps = 10)
})

macaque_gradient_color <- c('Glu CA1 1'=macaque_list$`Glu CA1 1`[1], 'Glu CA1 2'=macaque_list$`Glu CA1 1`[3],
                   'Glu CA2/3 1'="#FFBB78FF",
                   'Glu CA3/4 1'=macaque_list$`Glu CA3/4 1`[1], 'Glu CA3/4 2'=macaque_list$`Glu CA3/4 1`[3], 
                   'Glu DG 1'=macaque_list$`Glu DG 1`[1],'Glu DG 2'=macaque_list$`Glu DG 1`[3],'Glu DG 3'=macaque_list$`Glu DG 1`[5],'Glu DG 4'=macaque_list$`Glu DG 1`[7],
                   'Glu SUB 1'="#17BECFFF",
                   'Glu SUB deep 1'=macaque_list$`Glu SUB deep 1`[1],  'Glu SUB deep 2'=macaque_list$`Glu SUB deep 1`[3],'Glu SUB deep 3'=macaque_list$`Glu SUB deep 1`[5], 
                   'Glu SUB int. 1'=macaque_list$`Glu SUB int. 1`[1],'Glu SUB int. 2'=macaque_list$`Glu SUB int. 1`[2], 'Glu SUB int. 3'=macaque_list$`Glu SUB int. 1`[3], 'Glu SUB int. 4'=macaque_list$`Glu SUB int. 1`[4], 'Glu SUB int. 5'=macaque_list$`Glu SUB int. 1`[5], 'Glu SUB int. 6'=macaque_list$`Glu SUB int. 1`[6], 'Glu SUB int. 7'=macaque_list$`Glu SUB int. 1`[7],
                   'GABA CCK 1'=macaque_list$`GABA CCK 1`[1],'GABA CCK 2'=macaque_list$`GABA CCK 1`[3],'GABA CCK 3'=macaque_list$`GABA CCK 1`[5],'GABA CCK 4'=macaque_list$`GABA CCK 1`[7],
                   'GABA LAMP5 1'=macaque_list$`GABA LAMP5 1`[1], 'GABA LAMP5 2'=macaque_list$`GABA LAMP5 1`[3], 'GABA LAMP5 3'=macaque_list$`GABA LAMP5 1`[5], 'GABA LAMP5 4'=macaque_list$`GABA LAMP5 1`[7],
                   'GABA PVALB 1'=macaque_list$`GABA PVALB 1`[1],'GABA PVALB 2'=macaque_list$`GABA PVALB 1`[3], 'GABA PVALB 3'=macaque_list$`GABA PVALB 1`[5],
                   'GABA SST 1'=macaque_list$`GABA SST 1`[1], 'GABA SST 2'=macaque_list$`GABA SST 1`[3],
                   'GABA Vip 1'="#C49C94FF",
                   'astrocyte'="#C5B0D5FF", 'olig'="#DBDB8DFF", 'OPC'="#8C564BFF", 'VLMC'="#E377C2FF","microglia"="#BCBD22FF")
# Remove substrings from labels
cleaned_labels <- gsub("Glu ", "", names(macaque_gradient_color))
cleaned_labels <- gsub("GABA ", "", cleaned_labels)
# Create a new color vector with cleaned labels
cleaned_macaque_gradient_color <- setNames(macaque_gradient_color, cleaned_labels)
DimPlot(macaque_data1, cols = cleaned_macaque_gradient_color, label = FALSE, label.color = "white") +
  theme_void() +
  theme(legend.text = element_text(size = 25))

ggsave("./macaque/macaque_UMAP_fu1a.pdf",width = 10.5, height = 7)

### Figure B：Quality control situation of single cell small group
#### mouse
Idents(mouse_data1) <- factor(mouse_data1$subtype_rename,levels=mouse_new_celltype_levels)
p1 <- VlnPlot(mouse_data1, features = c("nCount_RNA"), pt.size = 0, cols = new_mouse_color) + NoLegend() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 32),
        plot.title = element_blank()) +
  ylab("nCount") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 40000, by = 20000))

p2 <- VlnPlot(mouse_data1, features = c("nFeature_RNA"), pt.size = 0, cols = new_mouse_color) + NoLegend() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 32),
        plot.title = element_blank()) +
  ylab("nFeature") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 8000, by = 4000))

p3 <- VlnPlot(mouse_data1, features = c("percent.mt"), pt.size = 0, cols = new_mouse_color) + NoLegend() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1,size = 35),
        plot.title = element_blank(),
        axis.text.y = element_text(size = 32))+
  ylab("percent.mt") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.2))
p3 <- p3 + scale_x_discrete(labels = function(x) gsub("(Glu|GABA)\\s*", "", x))
#p3 <- p3 + scale_y_continuous(limits = c(0, 5))

combined_plot <- p1 / p2 / p3
combined_plot

ggsave(combined_plot,filename = "./mouse/mouse_qc_fu1b.pdf",width =25, height = 9)

#### marmoset
Idents(marmoset_data1) <- factor(marmoset_data1$subtype_rename,levels=marmoset_new_celltype_levels)
p1 <- VlnPlot(marmoset_data1, features = c("nCount_RNA"), pt.size = 0, cols = new_marmoset_color) + NoLegend() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 32),
        plot.title = element_blank()) +
  ylab("nCount") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 40000, by = 20000))

p2 <- VlnPlot(marmoset_data1, features = c("nFeature_RNA"), pt.size = 0, cols = new_marmoset_color) + NoLegend() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 32),
        plot.title = element_blank()) +
  ylab("nFeature") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 8000, by = 4000))

p3 <- VlnPlot(marmoset_data1, features = c("percent.mt"), pt.size = 0, cols = new_marmoset_color) + NoLegend() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1,size = 35),
        plot.title = element_blank(),
        axis.text.y = element_text(size = 32))+
  ylab("percent.mt") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 5, by = 2))
p3 <- p3 + scale_x_discrete(labels = function(x) gsub("(Glu|GABA)\\s*", "", x))

combined_plot <- p1 / p2 / p3
combined_plot

ggsave(combined_plot,filename = "./marmoset/marmoset_qc_fu1b.pdf",width =25, height = 9)

#### macaque
Idents(macaque_data1) <- factor(macaque_data1$subtype_rename,levels=macaque_new_celltype_levels)
p1 <- VlnPlot(macaque_data1, features = c("nCount_RNA"), pt.size = 0, cols = new_macaque_color) + NoLegend() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 32),
        plot.title = element_blank()) +
  ylab("nCount") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 100000, by = 50000), labels = function(x) format(x, scientific = FALSE))

p2 <- VlnPlot(macaque_data1, features = c("nFeature_RNA"), pt.size = 0, cols = new_macaque_color) + NoLegend() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 32),
        plot.title = element_blank()) +
  ylab("nFeature") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 8000, by = 4000))

p3 <- VlnPlot(macaque_data1, features = c("percent.mt"), pt.size = 0, cols = new_macaque_color) + NoLegend() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1,size = 35),
        plot.title = element_blank(),
        axis.text.y = element_text(size = 32))+
  ylab("percent.mt") +
  theme(axis.title.y = element_text(size = 35)) +
  xlab("")+
  scale_y_continuous(breaks = seq(0, 5, by = 2))
p3 <- p3 + scale_x_discrete(labels = function(x) gsub("(Glu|GABA)\\s*", "", x))
combined_plot <- p1 / p2 / p3
combined_plot

ggsave(combined_plot,filename = "./macaque/macaque_qc_fu1b.pdf",width =25, height = 9)

### Figure CD：Comparison of similarities between each slice
#### mouse
#Two adjacent slices
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/mouse1/mouse1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
Merge_df <- all_Merge_df[all_Merge_df$Slice %in% c("T313","T312") ,]
a <- as.matrix(table(Merge_df$Cell_Type,Merge_df$Slice))
df <- as.data.frame(a)
reshaped_df <- dcast(df, Var1 ~ Var2, value.var = "Freq")
reshaped_df$T313 <- reshaped_df$T313 / sum(reshaped_df$T313)
reshaped_df$T312 <- reshaped_df$T312 / sum(reshaped_df$T312)
cor_result <- cor.test(reshaped_df$T313, reshaped_df$T312, method = 'pearson')
correlation <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", correlation))
print(paste("p-value:", p_value))
cor <- round(cor(reshaped_df$T313, reshaped_df$T312, method = 'pearson')[1],2)
p <- ggplot(reshaped_df, aes(x = T313, y = T312)) +
  geom_point(col = "grey", alpha = 0.5) +
  geom_smooth(method = lm, level = 0.90, color = 'red', fill = 'lightblue') +
  theme(panel.grid = element_blank(), text = element_text(size = 10)) +
  labs(x = '+2.61 (#1)', y = '+2.71 (#1)') +
  theme_classic() +
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -10),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 

p1 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))

#Similarity of slices of two mice at almost the same position
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/mouse1/mouse1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
mouse1_Merge_df <- all_Merge_df[all_Merge_df$Slice %in% c("T313") ,]
mouse1_Merge_df <- mouse1_Merge_df[,c("Cell_Type","Slice")]
mouse2_all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/mouse2/mouse2_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
mouse2_Merge_df <- mouse2_all_Merge_df[mouse2_all_Merge_df$Slice %in% c("T232") ,]
mouse2_Merge_df <- mouse2_Merge_df[,c("Cell_Type","Slice")]
Merge_df <- rbind(mouse1_Merge_df,mouse2_Merge_df)
a <- as.matrix(table(Merge_df$Cell_Type,Merge_df$Slice))
df <- as.data.frame(a)

reshaped_df <- dcast(df, Var1 ~ Var2, value.var = "Freq")
reshaped_df$T313 <- reshaped_df$T313 / sum(reshaped_df$T313)
reshaped_df$T232 <- reshaped_df$T232 / sum(reshaped_df$T232)
cor_result <- cor.test(reshaped_df$T313, reshaped_df$T232, method = 'pearson')
correlation <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", correlation))
print(paste("p-value:", p_value))
cor <- round(cor(reshaped_df$T313, reshaped_df$T232, method = 'pearson')[1],2)
p <- ggplot(reshaped_df, aes(x=T313, y=T232)) +
  geom_point(col="grey",alpha=0.5) +
  geom_smooth(method=lm, level=0.90,color='red', fill='lightblue') +
  theme(panel.grid = element_blank(), text = element_text(size = 15)) +
  labs(x='+2.61 (#1)', y='+2.66 (#2)',
       title = paste('R = ',cor,'  p<0.001',sep=''))+theme_classic() +
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -10),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 

p2 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))


#mouse
#The similarity of the composition of the slice ratio and the correlation of the distance
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/mouse1/mouse1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
Merge_df_matrix <- as.matrix(table(all_Merge_df$Cell_Type,all_Merge_df$ebz))
df1new <- Merge_df_matrix
library(reshape)
cor1=cor(df1new)
cor1=as.data.frame(cor1)
cor1$slide1=rownames(cor1)
cor1db=melt(cor1)
cor1db$slide1=factor(cor1db$slide1,levels=unique(cor1db$slide1))
cor1db$loc1=cor1db$slide1
cor1db$loc2=cor1db$variable
cor1db$loc1=as.numeric(as.character(cor1db$loc1))
cor1db$loc2=as.numeric(as.character(cor1db$loc2))
cor1db$distance=cor1db$loc2-cor1db$loc1
cor1db1=cor1db[cor1db$distance>=0,]

cor_result <- cor.test(cor1db1$value, cor1db1$distance, method = 'pearson')
cor <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", cor))
print(paste("p-value:", p_value))

library(ggplot2)
p <- ggplot(data=cor1db1, aes(x=distance, y=value)) +
  geom_point(col="grey",alpha=0.5) +
  geom_smooth(method=lm, level=0.90,color='red', fill='lightblue')+theme_classic()+
  labs(x='Distance between sections (mm)', y='Correlation coefficient',
       title = paste('R = ',cor,'  p<0.001',sep='')) +
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -35),  
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 
p3 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))
combined_plot <- grid.arrange(p1, p2, p3, ncol = 3)
combined_plot

ggsave(combined_plot,filename = "./mouse/mouse_correction_fu1c_20230304.pdf", width = 30, height = 10)

#### marmoset
#Two adjacent slices
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
Merge_df <- all_Merge_df[all_Merge_df$Slice %in% c("T474","T473") ,]
a <- as.matrix(table(Merge_df$Cell_Type,Merge_df$Slice))
df <- as.data.frame(a)
library(reshape2)
library(dplyr)
reshaped_df <- dcast(df, Var1 ~ Var2, value.var = "Freq")
reshaped_df$T474 <- reshaped_df$T474 / sum(reshaped_df$T474)
reshaped_df$T473 <- reshaped_df$T473 / sum(reshaped_df$T473)
cor_result <- cor.test(reshaped_df$T474, reshaped_df$T473, method = 'pearson')
correlation <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", correlation))
print(paste("p-value:", p_value))
cor <- round(cor(reshaped_df$T474, reshaped_df$T473, method = 'pearson')[1],2)
p <- ggplot(reshaped_df, aes(x=T474, y=T473)) +
  geom_point(col="grey",alpha=0.5) +
  geom_smooth(method=lm, level=0.90,color='red', fill='lightblue') +
  theme(panel.grid = element_blank(), text = element_text(size = 15)) +
  labs(x='+1.77 (#1)', y='+2.02 (#1)',
       title = paste('R = ',cor,'  p<0.001',sep=''))+theme_classic() +
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -10),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 

p1 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))

#The similarity of two marmoset slices at almost the same position
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
marmoset1_Merge_df <- all_Merge_df[all_Merge_df$Slice %in% c("T465") ,]
marmoset1_Merge_df <- marmoset1_Merge_df[,c("Cell_Type","Slice")]
marmoset2_all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset2/marmoset2_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_20240304.csv', row.names = NULL)
marmoset2_Merge_df <- marmoset2_all_Merge_df[marmoset2_all_Merge_df$Slice %in% c("T658") ,]
marmoset2_Merge_df <- marmoset2_Merge_df[,c("Cell_Type","Slice")]
Merge_df <- rbind(marmoset1_Merge_df,marmoset2_Merge_df)
a <- as.matrix(table(Merge_df$Cell_Type,Merge_df$Slice))
df <- as.data.frame(a)

reshaped_df <- dcast(df, Var1 ~ Var2, value.var = "Freq")
reshaped_df$T465 <- reshaped_df$T465 / sum(reshaped_df$T465)
reshaped_df$T658 <- reshaped_df$T658 / sum(reshaped_df$T658)
cor_result <- cor.test(reshaped_df$T465, reshaped_df$T658, method = 'pearson')
correlation <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", correlation))
print(paste("p-value:", p_value))
cor <- round(cor(reshaped_df$T465, reshaped_df$T658, method = 'pearson')[1],2)
p <- ggplot(reshaped_df, aes(x=T465, y=T658)) +
  geom_point(col="grey",alpha=0.5) +
  geom_smooth(method=lm, level=0.90,color='red', fill='lightblue') +
  theme(panel.grid = element_blank(), text = element_text(size = 15)) +
  labs(x='+3.97 (#1)', y='+4.0 (#2)',
       title = paste('R = ',cor,'  p<0.001',sep=''))+theme_classic()+
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -10),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 

p2 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))

#Draw the similarity of the composition of the marmoset slice ratio and the correlation of the distance
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset1/marmoset1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
Merge_df_matrix <- as.matrix(table(all_Merge_df$Cell_Type,all_Merge_df$ebz))
df1new <- Merge_df_matrix
library(reshape)
cor1=cor(df1new)
cor1=as.data.frame(cor1)
cor1$slide1=rownames(cor1)
cor1db=melt(cor1)
cor1db$slide1=factor(cor1db$slide1,levels=unique(cor1db$slide1))
cor1db$loc1=cor1db$slide1
cor1db$loc2=cor1db$variable
cor1db$loc1=as.numeric(as.character(cor1db$loc1))
cor1db$loc2=as.numeric(as.character(cor1db$loc2))
cor1db$distance=cor1db$loc2-cor1db$loc1
cor1db1=cor1db[cor1db$distance>=0,]

cor_result <- cor.test(cor1db1$value, cor1db1$distance, method = 'pearson')
cor <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", correlation))
print(paste("p-value:", p_value))

library(ggplot2)
p <- ggplot(data=cor1db1, aes(x=distance, y=value)) +
  geom_point(col="grey",alpha=0.5) +
  geom_smooth(method=lm, level=0.90,color='red', fill='lightblue')+theme_classic()+
  labs(x='Distance between sections (mm)', y='Correlation coefficient',
       title = paste('R = ',cor,'  p<0.001',sep=''))+
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -35),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 
p3 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))
combined_plot <- grid.arrange(p1, p2, p3, ncol = 3)
combined_plot

ggsave(combined_plot,filename = "./marmoset/marmoset_correction_fu1c_20230304.pdf.pdf", width = 30, height = 10)

#### macaque

#Two adjacent slices
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque1/macaque1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
Merge_df <- all_Merge_df[all_Merge_df$Slice %in% c("T33","T31") ,]
a <- as.matrix(table(Merge_df$Cell_Type,Merge_df$Slice))
df <- as.data.frame(a)
library(reshape2)
library(dplyr)
reshaped_df <- dcast(df, Var1 ~ Var2, value.var = "Freq")
reshaped_df$T33 <- reshaped_df$T33 / sum(reshaped_df$T33)
reshaped_df$T31 <- reshaped_df$T31 / sum(reshaped_df$T31)
cor_result <- cor.test(reshaped_df$T33, reshaped_df$T31, method = 'pearson')
correlation <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", correlation))
print(paste("p-value:", p_value))
cor <- round(cor(reshaped_df$T33, reshaped_df$T31, method = 'pearson')[1],2)
p <- ggplot(reshaped_df, aes(x=T33, y=T31)) +
  geom_point(col="grey",alpha=0.5) +
  geom_smooth(method=lm, level=0.90,color='red', fill='lightblue') +
  theme(panel.grid = element_blank(), text = element_text(size = 15)) +
  labs(x='+10 (#1)', y='+9.51 (#1)',
       title = paste('R = ',cor,'  p<0.001',sep=''))+theme_classic()+
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -10),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 

p
p1 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))

#Similarity of two macaque slices at almost the same position
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque1/macaque1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
macaque1_Merge_df <- all_Merge_df[all_Merge_df$Slice %in% c("T33") ,]
macaque1_Merge_df <- macaque1_Merge_df[,c("Cell_Type","Slice")]
macaque2_all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/macaque2/merge/macaque2_maintype_GABA_Glu_subtype_rxry_filter.csv', row.names = NULL)
macaque2_Merge_df <- macaque2_all_Merge_df[macaque2_all_Merge_df$Slice %in% c("T573") ,]
macaque2_Merge_df <- macaque2_Merge_df[,c("Cell_Type","Slice")]
Merge_df <- rbind(macaque1_Merge_df,macaque2_Merge_df)
a <- as.matrix(table(Merge_df$Cell_Type,Merge_df$Slice))
df <- as.data.frame(a)

reshaped_df <- dcast(df, Var1 ~ Var2, value.var = "Freq")
reshaped_df$T33 <- reshaped_df$T33 / sum(reshaped_df$T33)
reshaped_df$T573 <- reshaped_df$T573 / sum(reshaped_df$T573)
cor_result <- cor.test(reshaped_df$T33, reshaped_df$T573, method = 'pearson')
correlation <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", correlation))
print(paste("p-value:", p_value))
cor <- round(cor(reshaped_df$T33, reshaped_df$T573, method = 'pearson')[1],2)
p <- ggplot(reshaped_df, aes(x=T33, y=T573)) +
  geom_point(col="grey",alpha=0.5) +
  geom_smooth(method=lm, level=0.90,color='red', fill='lightblue') +
  theme(panel.grid = element_blank(), text = element_text(size = 15)) +
  labs(x='+10 (#1)', y='+9.06 (#2)',
       title = paste('R = ',cor,'  p<0.001',sep=''))+theme_classic()+
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -10),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 

p

p2 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))

#Draw the similarity of the composition of the macaque slice ratio and the correlation of the distance
all_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque1/macaque1_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = NULL)
Merge_df_matrix <- as.matrix(table(all_Merge_df$Cell_Type,all_Merge_df$ebz))
df1new <- Merge_df_matrix
library(reshape)
cor1=cor(df1new)
cor1=as.data.frame(cor1)
cor1$slide1=rownames(cor1)
cor1db=melt(cor1)
cor1db$slide1=factor(cor1db$slide1,levels=unique(cor1db$slide1))
cor1db$loc1=cor1db$slide1
cor1db$loc2=cor1db$variable
cor1db$loc1=as.numeric(as.character(cor1db$loc1))
cor1db$loc2=as.numeric(as.character(cor1db$loc2))
cor1db$distance=cor1db$loc2-cor1db$loc1
cor1db1=cor1db[cor1db$distance>=0,]

cor_result <- cor.test(cor1db1$value, cor1db1$distance, method = 'pearson')
cor <- round(cor_result$estimate, 2)
p_value <- cor_result$p.value
print(paste("Correlation:", correlation))
print(paste("p-value:", p_value))

library(ggplot2)
p <- ggplot(data=cor1db1, aes(x=distance, y=value)) +
  geom_point(col="grey",alpha=0.5) +
  geom_smooth(method=lm, level=0.90,color='red', fill='lightblue')+theme_classic()+
  labs(x='Distance between sections (mm)', y='Correlation coefficient',
       title = paste('R = ',cor,'  p<0.001',sep=''))+
  theme(
    axis.title = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.1, vjust = -35),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45)
  ) 
p3 <- p + ggtitle(paste("R =", cor, "\n", "p < 0.001"))
combined_plot <- grid.arrange(p1, p2, p3, ncol = 3)
combined_plot

ggsave(combined_plot,filename = "./macaque/macaque_correction_fu1c_20230304.pdf.pdf", width = 30, height = 10)

### Figure E：Correlation between single cell data and spatial group data
#### mouse
celltype_NV1 <- read.csv("E:/seurat_spatial/transfer_0818/metaneighbor/mouse1/celltype_NV1_20230304.csv",check.names = FALSE,row.names=1)
colnames(celltype_NV1) <- gsub("(Glu|GABA)\\s*", "", colnames(celltype_NV1))
rownames(celltype_NV1) <- gsub("(Glu|GABA)\\s*", "", rownames(celltype_NV1))
breaks <- seq(0.1,1,0.01)
p1 <- pheatmap(celltype_NV1, cluster_rows = F, cluster_cols = F, 
              breaks = breaks, legend_breaks = seq(0,1,0.2),
              display_numbers=F,number_format = "%.2f",
              fontsize = 15,border_color = "#696969",angle_col = 90,
              color = colorRampPalette(c("navy","white", "red"))(length(breaks)),fontsize_row=15,fontsize_col=15,legend = F)

p1

#### marmoset
celltype_NV1 <- read.csv("E:/seurat_spatial/transfer_0818/metaneighbor/marmoset1/celltype_NV1_20230304.csv",check.names = FALSE,row.names=1)
colnames(celltype_NV1) <- gsub("(Glu|GABA)\\s*", "", colnames(celltype_NV1))
rownames(celltype_NV1) <- gsub("(Glu|GABA)\\s*", "", rownames(celltype_NV1))
breaks <- seq(0.1,1,0.01)
p2 <- pheatmap(celltype_NV1, cluster_rows = F, cluster_cols = F, 
              breaks = breaks, legend_breaks = seq(0,1,0.2),
              display_numbers=F,number_format = "%.2f",
              fontsize = 15,border_color = "#696969",angle_col = 90,
              color = colorRampPalette(c("navy","white", "red"))(length(breaks)),fontsize_row=15,fontsize_col=15,legend = F)
p2

#### macaque
celltype_NV1 <- read.csv("E:/seurat_spatial/transfer_0818/metaneighbor/macaque1/celltype_NV1_20230304.csv",check.names = FALSE,row.names=1)
colnames(celltype_NV1) <- gsub("(Glu|GABA)\\s*", "", colnames(celltype_NV1))
rownames(celltype_NV1) <- gsub("(Glu|GABA)\\s*", "", rownames(celltype_NV1))
breaks <- seq(0.1,1,0.01)
p3 <- pheatmap(celltype_NV1, cluster_rows = F, cluster_cols = F, 
              breaks = breaks, legend_breaks = seq(0,1,0.2),
              display_numbers=F,number_format = "%.2f",
              fontsize = 15,border_color = "#696969",angle_col = 90,
              color = colorRampPalette(c("navy","white", "red"))(length(breaks)),fontsize_row=15,fontsize_col=15,legend = F)
p3
p <- cowplot::plot_grid(p3$gtable, p2$gtable, p1$gtable,ncol= 3, align = "h") 
p <- p + theme(plot.margin = margin(25, 0, 25, 0))
p

ggsave(p, filename = "./three_species/three_cor_heatmap_20240304.pdf",width =27, height = 10)

### Picture F：The slices of the second animal and the picture after the slice transfer of the large group
#### mouse2
mouse_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/mouse2/mouse2_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = 1)
mouse_Merge_df <- mouse_Merge_df[mouse_Merge_df$Slice %in% c("T230","T229","T228","T227","T226","T225"),]
plot_dir <- './mouse2/plot_png/'
summary(mouse_Merge_df$fx)
summary(mouse_Merge_df$fy)
# Allcelltype for each slice
Merge_df <- mouse_Merge_df
celltype_levels <- names(mouse_color)
Merge_df$celltype_new <- factor(Merge_df$subtype_rename, levels = celltype_levels)
Slices <- names(table(Merge_df$Slice))
for (i in 1:length(Slices)){
  plot_df <- Merge_df[which(Merge_df$Slice==Slices[i]),]
  # Maintype_merge
  p1 <- merge_plot(plot_df, group.by='celltype_new', group.cols=mouse_color, pt.size=0.1, title=Slices[i])
  p1_file <- paste(plot_dir, Slices[i], '_Allcelltype.png', sep = '')
  ggsave(p1, filename = p1_file, width = 3, height = 3, bg = 'transparent')
}


#### marmoset2
marmoset_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/marmoset2/marmoset2_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_20240304.csv', row.names = 1)

marmoset_Merge_df <- marmoset_Merge_df[marmoset_Merge_df$Slice %in% c("T652","T654","T658","T661"),]
plot_dir <- './marmoset2/plot_png/'
summary(marmoset_Merge_df$fx)
summary(marmoset_Merge_df$fy)
# Allcelltype for each slice
Merge_df <- marmoset_Merge_df
celltype_levels <- names(marmoset_color)
Merge_df$celltype_new <- factor(Merge_df$subtype_rename, levels = celltype_levels)
Slices <- names(table(Merge_df$Slice))
for (i in 1:length(Slices)){
  plot_df <- Merge_df[which(Merge_df$Slice==Slices[i]),]
  # Maintype_merge
  p1 <- merge_plot(plot_df, group.by='celltype_new', group.cols=marmoset_color, pt.size=1, title=Slices[i])
  p1_file <- paste(plot_dir, Slices[i], '_Allcelltype.png', sep = '')
  ggsave(p1, filename = p1_file, width = 3, height = 3, bg = 'transparent')
}

#### macaque2
macaque_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque2/macaque2_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_20240304.csv', row.names = 1)
macaque_Merge_df <- macaque_Merge_df[macaque_Merge_df$Slice %in% c("T573","T574","T575","T576"),]
plot_dir <- './macaque2/plot_png/'
summary(macaque_Merge_df$fx)
summary(macaque_Merge_df$fy)

# Allcelltype for each slice
Merge_df <- macaque_Merge_df
celltype_levels <- names(macaque_color)
Merge_df$celltype_new <- factor(Merge_df$subtype_rename, levels = celltype_levels)
Slices <- names(table(Merge_df$Slice))
for (i in 1:length(Slices)){
  plot_df <- Merge_df[which(Merge_df$Slice==Slices[i]),]
  # Maintype_merge
  p1 <- merge_plot(plot_df, group.by='celltype_new', group.cols=macaque_color, pt.size=0.1, title=Slices[i])
  p1_file <- paste(plot_dir, Slices[i], '_Allcelltype.png', sep = '')
  ggsave(p1, filename = p1_file, width = 3, height = 3, bg = 'transparent')
}

#### macaque3
macaque_Merge_df <- read.csv('E:/seurat_spatial/transfer_0818/results_0919/filter_data/macaque3/macaque3_maintype_GABA_Glu_subtype_rxry_rename_self_rotate_manual_filter_cell_slice_20240304.csv', row.names = 1)
plot_dir <- './macaque3/plot_png/'
summary(macaque_Merge_df$fx)
summary(macaque_Merge_df$fy)
# Allcelltype for each slice
Merge_df <- macaque_Merge_df
celltype_levels <- names(macaque_color)
Merge_df$celltype_new <- factor(Merge_df$subtype_rename, levels = celltype_levels)
Slices <- names(table(Merge_df$Slice))
for (i in 1:length(Slices)){
  plot_df <- Merge_df[which(Merge_df$Slice==Slices[i]),]
  # Maintype_merge
  p1 <- merge_plot(plot_df, group.by='celltype_new', group.cols=macaque_color, pt.size=0.1, title=Slices[i])
  p1_file <- paste(plot_dir, Slices[i], '_Allcelltype.png', sep = '')
  ggsave(p1, filename = p1_file, width = 3, height = 3, bg = 'transparent')
}
                            
