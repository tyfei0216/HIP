#============================A
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
set.seed(42)
.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")


sample.combined=readRDS('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/GLu/seurat_Macaca_Marmaset_Mouse_SCT_intergrat_WuZhongJian_GLu.100.1500.1.rds.gz')
DefaultAssay(sample.combined)="SCT"

.cluster_cols3=c('#841F5E','#E37644','#486BB3')

plots1=DimPlot(sample.combined, split.by = c("orig.ident_species"), combine = T, 
			cols =.cluster_cols,label=T,order = c('macaque','marmoset','mouse'))+NoLegend()
print(plots1) # group.by = c("name_transfer2"),

			
.cluster_cols3=c('#841F5E','#E37644','#486BB3')
plots1=DimPlot(sample.combined, group.by = c("orig.ident_species"), combine = T, cols =.cluster_cols3,pt.size=0.1,
			order = c('macaque','marmoset','mouse'))+NoLegend() 
print(plots1)


sample.combined$lingzhang_GRIA4=paste(sample.combined$orig.ident_species,
	sample.combined$cluster_gene,sep='_')

sample.combined$name_transfer_lingzhang=sample.combined$lingzhang_GRIA4
sample.combined$name_transfer_lingzhang=ifelse(sample.combined$name_transfer_lingzhang %in% c('macaque_19','marmoset_GRIA4'),
						sample.combined$name_transfer_lingzhang,'other')
						
		
DimPlot(sample.combined, group.by = c("name_transfer_lingzhang"), split.by = c("orig.ident_species"), combine = T, 
				cols = c('grey',.cluster_cols),
			label=T,order = c("macaque_19","marmoset_GRIA4","other"))+NoLegend() 

#===============B			
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
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}

sample.combined=readRDS('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/GLu/seurat_Macaca_Marmaset_Mouse_SCT_intergrat_WuZhongJian_GLu.100.1500.1.second.rds.gz')
our.combined=sample.combined
our.combined=our.combined@meta.data

our.combined$label_for_heatmap=paste(our.combined$orig.ident_species, our.combined$field,our.combined$cluster_id,our.combined$cluster_gene, sep = "_")

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
temp_mac_mou=cocl.subset2

heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)
#heat.colors <- colorRampPalette(c("white", "grey70", "black"))(5)

cocl.subset_plot1=pheatmap(cocl.subset2, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat.colors,
         fontsize = 12) 			


cocl.subset <- cocl[grepl(compare.species[1], row.names(cocl)),
                    grepl(compare.species[2], row.names(cocl))] #change 3 to 2 for marmoset

row.names(cocl.subset) <- sub(paste0(compare.species[1], "_"), "",
                              row.names(cocl.subset))

colnames(cocl.subset) <- sub(paste0(compare.species[2], "_"), "", #change 3 to 2 for marmoset
                             colnames(cocl.subset))


heat.colors <- colorRampPalette(c("white", "grey70", "black"))(5)
heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)
#heat.colors <- =rev(RColorBrewer::brewer.pal(7,'RdYlBu'))
#heat.colors <- colorRampPalette(c("white", "grey70", "black"))(50)

cocl.subset2=cocl.subset
temp_mac_mar=cocl.subset2

cocl.subset_plot2=pheatmap(cocl.subset2, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat.colors,
         fontsize = 12) 

#=====D
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
slice_use=which(names(cellinfor_type)=='T30')

mytype = 'Glu 19 (CA)'
for(ii1 in c(slice_use)){
    xx=names(cellinfor_type)[ii1]
	#============================每张 空间的片子
	#T59 <- readRDS(paste('/work/dy/project/temp/singlecell/hippocampus/spcace_transfer/data/macaque/',xx,'.rds',sep=''))
	#T59 <- readRDS(paste('/home/dongyu/project/hippo/spaceResult/transfer_from_xiaojuan20230926/macaque1/spatial_rds_filter/',xx,'.rds',sep=''))
	T59 <- readRDS(paste('/home/share/share_gouxj/transfer_0818_result/results_0919/filter_data/macaque1/spatial_rds_filter/',xx,'.rds',sep=''))
	
  #T59$Cell_Type=subset(cellinfor,Slice==xx)$Cell_Type
	#T59$gene_area=subset(cellinfor,Slice==xx)$gene_area
	#所有细胞放一起
	#T59$all_cell='all_cell_group'
			  temp_infor=subset(cellinfor,Slice==xx)
  #加入Cell_Type
  T59$Cell_Type=temp_infor$Cell_Type[match(T59$cell,temp_infor$cell)]


    #所有细胞放一起
	T59$all_cell='all_cell_group'
    umap <- new(
		Class = 'DimReduc',
		cell.embeddings = as.matrix(T59@meta.data[,c('x','y')])
        )
    seurat_obj_umap = T59
    seurat_obj_umap@reductions[["umap"]]=umap #将 对象中 的unmap 坐标变成 真实的空间 坐标
        seurat_obj_umap@reductions$umap@key="UMAP_"
	colnames(seurat_obj_umap@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")

       seurat_obj_change=seurat_obj_umap
    #每个细胞类型 单独画 原位图
    Idents(seurat_obj_change) = "Cell_Type"
    #设置每个要画的celltype
    #mytype = as.character(unique(seurat_obj_change@meta.data$Cellmarker))

  temp.seurat = seurat_obj_change
  metadata = temp.seurat@meta.data
  metadata[,"tempType"] = ifelse(metadata$Cell_Type==mytype,mytype,"Other")
  
  temp.seurat@meta.data = metadata
  Idents(temp.seurat) = "tempType"
  
  #colors = c(unique(as.character(celltypecolor))[i],adjustcolor("gray95",0.4))
  colors = c('yellow','navy')
  names(colors) = c(mytype,"Other")
  
  #filename = gsub("-|/",".",mytype[i])
  SliceName = xx
  
  p = DimPlot(temp.seurat,
              reduction = "umap",group.by='tempType',
              cols = colors,
			  #cols = c('navy','yellow'),
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


slice_use=which(names(cellinfor_type)=='T31')

mytype = 'Glu 19 (CA)'

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

mytype = 'Glu 0 GRIA4 (CA)'

slice_use=which(names(cellinfor_type)=='T463')

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

mytype = 'Glu 0 GRIA4 (CA)'
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

#======================E
library(Seurat)
set.seed(42)
.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

#==========macaque

seurat_Glu.markers.subtypes=readRDS(paste("/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/GLu/Glumonkey_Glu_cell_res_markers_eachCluster_filter.negative.positive.rds.gz",sep=''))

SubDataSig <- subset(seurat_Glu.markers.subtypes,cluster=='19')
SubDataSig$Sig = ifelse(SubDataSig$avg_log2FC > 0.5 & SubDataSig$p_val < 0.05, "Up",ifelse(SubDataSig$avg_log2FC < -0.5 & SubDataSig$p_val < 0.05,"Down","None"))
SubDataSig$cluster= as.character(SubDataSig$cluster)

top50Genes_up=SubDataSig %>%  dplyr::top_n(5, avg_log2FC )
top50Genes_down=SubDataSig %>%  dplyr::top_n(-5, avg_log2FC )

SubDataSigTop10 <- rbind(top50Genes_up,top50Genes_down)

p <- ggplot(SubDataSig,aes(x= avg_log2FC,y= -log10(p_val)))+
    geom_point(aes(shape=Sig,color=Sig))+labs(x="log2(Foldchange)",y="-log10(p-value)",title='GRIA4 vs other')+
    scale_shape_manual(values=c("Up"=17,"Down"=15,"None"=20),guide=F)+
    scale_color_manual(values=c("Up"="#e41a1c","Down"="#3182bd","None"="lightgray"),guide=F)+
   theme(panel.background = element_rect(fill=NA,color="black"),axis.title= element_text(color="black",size=12),panel.grid.major = element_blank(),
          axis.text= element_text(color="black",size=10),legend.position = c(0.9,0.1))+
    geom_vline(xintercept = c(-0.5,0.5),linetype="dashed")+
    geom_hline(yintercept = -log10(0.05),linetype="dashed")+
    ggrepel::geom_text_repel(data=SubDataSigTop10,aes(x= avg_log2FC,y= -log10(p_val),
	label=gene),force=1,max.overlaps =30)

print(p)

library(org.Hs.eg.db)
seurat_Purkinje.markers_enrichment=lapply(split(subset(SubDataSig,avg_log2FC>=0.5),subset(SubDataSig,avg_log2FC>=0.5)$cluster),function(xx){
seurat_Purkinje.markers.markers_top100 <- xx %>% dplyr::group_by(cluster) %>% dplyr::top_n(100, avg_log2FC) %>% unique()


		Tissue_GO <- clusterProfiler::enrichGO(gene= seurat_Purkinje.markers.markers_top100$gene %>%  unique(), 
                                         keyType = "SYMBOL",
                                         OrgDb         = "org.Hs.eg.db",
                                         ont           = "BP",
                                         pvalueCutoff = 0.05,
										 pAdjustMethod = "fdr",
                                         minGSSize = 2,
                                         qvalueCutoff  = 1) 
		Tissue_GO <- data.frame(Tissue_GO)
		list(Tissue_GO)
		}
)

	xx=1
	x=seurat_Purkinje.markers_enrichment[[xx]]
	temp_KEGG=x[[1]]
	temp_KEGG$log10p=-log10(temp_KEGG$pvalue)
	sortXKEGG=temp_KEGG[order(temp_KEGG$log10p,decreasing = T),]#,temp_KEGG[order(temp_KEGG$pvalue,decreasing = T),])
	#top 10
	if(nrow(sortXKEGG)>=10){
	sortXKEGG=sortXKEGG[c(10:1),]	
	}else{
	sortXKEGG
	}

	
	ggbarplot(sortXKEGG, x="Description", y='log10p',fill='log10p',col=NA,title='Function',xlab='KEGG',ylab='-log10(P-value)',rotate=TRUE)+
				scale_fill_gradient2(low='white',high='red')+RotatedAxis()
								

#==========marmoset
ii=1
cut=0.1
seurat_Glu.markers.subtypes=readRDS(paste("/work/dy/project/temp/singlecell/hippocampus/Marmaset_MJ/analysis_d/Glu/plasticity/",ii,'_',cut,"_markers_Glutamate_eachCluster_final_annotation_pos_negative_remove13.rds.gz",sep=''))

SubDataSig <- subset(seurat_Glu.markers.subtypes,cluster=='0')
SubDataSig$Sig = ifelse(SubDataSig$avg_log2FC > 0.5 & SubDataSig$p_val < 0.05, "Up",ifelse(SubDataSig$avg_log2FC < -0.5 & SubDataSig$p_val < 0.05,"Down","None"))
SubDataSig$cluster= as.character(SubDataSig$cluster)

top50Genes_up=SubDataSig %>%  dplyr::top_n(10, avg_log2FC )
top50Genes_down=SubDataSig %>%  dplyr::top_n(-10, avg_log2FC )

SubDataSigTop10 <- rbind(top50Genes_up,top50Genes_down)

p <- ggplot(SubDataSig,aes(x= avg_log2FC,y= -log10(p_val)))+
    geom_point(aes(shape=Sig,color=Sig))+labs(x="log2(Foldchange)",y="-log10(p-value)",title='GRIA4 vs other')+
    scale_shape_manual(values=c("Up"=17,"Down"=15,"None"=20),guide=F)+
    scale_color_manual(values=c("Up"="#e41a1c","Down"="#3182bd","None"="lightgray"),guide=F)+
   theme(panel.background = element_rect(fill=NA,color="black"),axis.title= element_text(color="black",size=12),panel.grid.major = element_blank(),
          axis.text= element_text(color="black",size=10),legend.position = c(0.9,0.1))+
    geom_vline(xintercept = c(-0.5,0.5),linetype="dashed")+
    geom_hline(yintercept = -log10(0.05),linetype="dashed")+
    ggrepel::geom_text_repel(data=SubDataSigTop10,aes(x= avg_log2FC,y= -log10(p_val),
			label=gene),force=1,max.overlaps =30)
print(p)


library(org.Hs.eg.db)

seurat_Purkinje.markers_enrichment=lapply(split(subset(SubDataSig,avg_log2FC>=0.5),subset(SubDataSig,avg_log2FC>=0.5)$cluster),function(xx){

seurat_Purkinje.markers.markers_top100 <- xx %>% dplyr::group_by(cluster) %>% dplyr::top_n(100, avg_log2FC) %>% unique()

		Tissue_GO <- clusterProfiler::enrichGO(gene= seurat_Purkinje.markers.markers_top100$gene %>%  unique(), #靶基因,
                                         #universe      = names(geneList),
                                         keyType = "SYMBOL",
                                         OrgDb         = "org.Hs.eg.db",
                                         ont           = "BP",
                                         pvalueCutoff = 0.05,
										 pAdjustMethod = "fdr",
                                         minGSSize = 2,
                                         qvalueCutoff  = 1) 
		##Tissue_GO <- clusterProfiler::simplify(Tissue_GO,cutoff=0.7)										 
		Tissue_GO <- data.frame(Tissue_GO)
		list(Tissue_GO)
		}
)

	
	xx=1
	print(xx)
	x=seurat_Purkinje.markers_enrichment[[xx]]
	print('go')
	temp_GO=x[[1]]
	temp_GO$log10p=-log10(temp_GO$pvalue)
	sortX=temp_GO[order(temp_GO$log10p,decreasing = T),]#,temp_KEGG[order(temp_KEGG$pvalue,decreasing = T),])
	#top 10
	if(nrow(sortX)>=10){
	sortX=sortX[c(1:10),]	
	}else{
	sortX
	}
	ggbarplot(sortXKEGG, x="Description", y='log10p',fill='log10p',col=NA,title='Function',xlab='KEGG',ylab='-log10(P-value)',rotate=TRUE)+
				scale_fill_gradient2(low='white',high='red')+RotatedAxis()
				
#==========G-H

#===============================macaque
seurat_Macaca_GLu=readRDS('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/Macaque_seurat_GABA_GLu_other_seurat.glu.rds.gz')

filename='macaque_GRIA4_vs_other_ampa1.2.3.4'

#saveRDS(seurat_Macaca_GLu,'/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/Macaque_seurat_GABA_GLu_other_seurat.glu.rds.gz')


#exp
seurat_Macaca_GLu_exp=as.matrix(seurat_Macaca_GLu@assays$SCT@data)[c('GRIA4','GRIA3','GRIA2','GRIA1','CACNG2','CACNG3','CACNG4'),]
#cell type
meta=seurat_Macaca_GLu@meta.data
meta$subclass_label=ifelse(meta$cluster_id=='19','GRIA4_cell','other')#seurat_Macaca_GLu$field

dbmat=as.data.frame(t(seurat_Macaca_GLu_exp))
dbmat$subcluster=meta$subclass_label
dbmat$region=meta$MainType
dbm=melt(dbmat)
dbm$subcluster=factor(dbm$subcluster,levels = unique(dbm$subcluster))#[c(12,19,17,9,8,24,15,6,3,4,5,13,1,7,21,10,22,20,16,18,14,11,31,30,29,28,26,25,27,23,2,32,34,37,35,36,33)])


dbmat$subcluster=factor(dbmat$subcluster,levels = levels(dbm$subcluster))
dbco1=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:4){
    dbco1[((i-1)*4+j),1]=levels(dbmat$subcluster)[i]
    dbco1[((i-1)*4+j),2]=colnames(dbtmp[j])
    dbco1[((i-1)*4+j),3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,j]==dbtmp$sum,])
    dbco1[((i-1)*4+j),4]=nrow(dbtmp)
  }
}

dbco2=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:3){
    for (k in (j+1):4){
      dbtmp$sum1=dbtmp[,j]+dbtmp[,k]
      n=n+1
      dbco2[n,1]=levels(dbmat$subcluster)[i]
      dbco2[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k])
      dbco2[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp$sum1==dbtmp$sum,])
      dbco2[n,4]=nrow(dbtmp)
    }
  }
}

dbco3=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:2){
    for (k in (j+1):3){
      for (m in (k+1):4){
        dbtmp$sum1=dbtmp[,j]+dbtmp[,k]+dbtmp[,m]
        n=n+1
        dbco3[n,1]=levels(dbmat$subcluster)[i]
        dbco3[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k],colnames(dbtmp)[m])
        dbco3[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp[,m]>0 & dbtmp$sum1==dbtmp$sum,])
        dbco3[n,4]=nrow(dbtmp)
      }
    }
  }
}

dbco4=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:1){
    for (k in (j+1):2){
      for (m in (k+1):3){
        for (p in (m+1):4){
          dbtmp$sum1=dbtmp[,j]+dbtmp[,k]+dbtmp[,m]+dbtmp[,p]
          n=n+1
          dbco4[n,1]=levels(dbmat$subcluster)[i]
          dbco4[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k],colnames(dbtmp)[m],colnames(dbtmp)[p])
          dbco4[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp[,m]>0 & dbtmp[,p]>0 & dbtmp$sum1==dbtmp$sum,])
          dbco4[n,4]=nrow(dbtmp)
        }
      }
    }
  }
}


dbco=rbind(dbco1,dbco2,dbco3,dbco4)
dbco$ratio=dbco$num/dbco$totalnum
dbco$genes=factor(dbco$genes,levels = unique(dbco$genes))
dbco$subcluster=factor(dbco$subcluster,levels = levels(dbmat$subcluster))

library(ggplot2)

p2=ggplot(dbco, aes(genes,subcluster)) +
  geom_tile(aes(fill = ratio)) + 
  scale_fill_gradient(low = "white", high = "red") +theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2

#==两两基因可视化	
dbco2=read.table(text = "",col.names = c("subcluster","gene1","gene2","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:7]
  for (j in 1:6){
    for (k in (j+1):7){
      n=n+1
      dbco2[n,1]=levels(dbmat$subcluster)[i]
      dbco2[n,2]=colnames(dbtmp)[j]
      dbco2[n,3]=colnames(dbtmp)[k]
      dbco2[n,4]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0,])
      dbco2[n,5]=nrow(dbtmp[dbtmp[,j]>0 | dbtmp[,k]>0,])
    }
  }
}
dbco2$ratio=dbco2$num/dbco2$totalnum
#dbco2$gene1=apply(dbco2,1,function(x){strsplit(x[2]," ")[[1]][1]})
#dbco2$gene2=apply(dbco2,1,function(x){strsplit(x[2]," ")[[1]][2]})
dbco2$subcluster=factor(dbco2$subcluster,levels = levels(dbmat$subcluster))


		#================red
p3=ggplot(dbco2, aes(gene1, gene2)) +
  geom_tile(aes(fill = ratio),colour="black") + 
  geom_text(aes(label = round(ratio, 2)),size=3)+
  scale_fill_gradient(low = "white", high = "red") +theme_minimal()+theme(axis.text.x = element_text(angle = 45,hjust = 1),panel.grid.major = element_blank())
p3=p3+facet_wrap(~subcluster)
p3



 macaque_name=openxlsx::read.xlsx('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/name_for_clsuter/rename2 的副本.xlsx',
	sheet='macaque')
 macaque_name=macaque_name[grep('Glu', macaque_name$subtype),]
 macaque_name$subtype2=sub(' \\(CA)| \\(DG)','',macaque_name$subtype)
  macaque_name$subtype2=gsub(' ','_',macaque_name$subtype2)
   macaque_name$subtype2=gsub('Glu_','',macaque_name$subtype2)
macaque_name$new=paste(macaque_name[,2],macaque_name[,3],macaque_name[,4],sep='_')


filename='macaque_GRIA4_vs_other_ampa1.2.3.4-new-names-sub'

seurat_Macaca_GLu_exp=as.matrix(seurat_Macaca_GLu@assays$SCT@data)[c('GRIA4','GRIA3','GRIA2','GRIA1','CACNG2','CACNG3','CACNG4'),]
#cell type
meta=seurat_Macaca_GLu@meta.data

#新名字
#meta$subtype =paste(meta$MainType,meta$cluster_id,meta$cluster_gene,sep='_')
meta=inner_join(meta,macaque_name,by=c('cluster_id'='subtype2'))

meta$subclass_label= meta$new #ifelse(meta$cluster_id=='4','GRIA4_cell','other')#seurat_Macaca_GLu$field


#meta$subclass_label=ifelse(meta$cluster_id=='19','GRIA4_cell','other')#seurat_Macaca_GLu$field

dbmat=as.data.frame(t(seurat_Macaca_GLu_exp))
dbmat$subcluster=meta$subclass_label
dbmat$region=meta$MainType

#==提取 sub 类群
dbmat=dbmat[grep('SUB',dbmat$subcluster),]

dbm=melt(dbmat)
dbm$subcluster=factor(dbm$subcluster,levels = unique(dbm$subcluster))#[c(12,19,17,9,8,24,15,6,3,4,5,13,1,7,21,10,22,20,16,18,14,11,31,30,29,28,26,25,27,23,2,32,34,37,35,36,33)])

if(F){
p=ggviolin(dbm,x = "subcluster",y = "value",fill = "subcluster")
p1=p+facet_wrap(~variable,ncol = 1)+
	theme(axis.text.x = element_text(angle = 45,hjust = 1))+#ylim(0,6)+
	 theme(legend.position = "right") 
ggsave(paste0("/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/GLu/AMAP_zuhe_exp_ratio/",
		filename,".violin.pdf"),p1,height=8,width=5)
}
	 
if(F){
p=ggviolin(dbm,x = "region",y = "value",fill = "region")
p+facet_wrap(~variable,ncol = 1)+theme(axis.text.x = element_text(angle = 45,hjust = 1))+ylim(0,4)+
	 theme(legend.position = "right") 
}

dbmat$subcluster=factor(dbmat$subcluster,levels = levels(dbm$subcluster))
dbco1=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:4){
    dbco1[((i-1)*4+j),1]=levels(dbmat$subcluster)[i]
    dbco1[((i-1)*4+j),2]=colnames(dbtmp[j])
    dbco1[((i-1)*4+j),3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,j]==dbtmp$sum,])
    dbco1[((i-1)*4+j),4]=nrow(dbtmp)
  }
}

dbco2=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:3){
    for (k in (j+1):4){
      dbtmp$sum1=dbtmp[,j]+dbtmp[,k]
      n=n+1
      dbco2[n,1]=levels(dbmat$subcluster)[i]
      dbco2[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k])
      dbco2[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp$sum1==dbtmp$sum,])
      dbco2[n,4]=nrow(dbtmp)
    }
  }
}

dbco3=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:2){
    for (k in (j+1):3){
      for (m in (k+1):4){
        dbtmp$sum1=dbtmp[,j]+dbtmp[,k]+dbtmp[,m]
        n=n+1
        dbco3[n,1]=levels(dbmat$subcluster)[i]
        dbco3[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k],colnames(dbtmp)[m])
        dbco3[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp[,m]>0 & dbtmp$sum1==dbtmp$sum,])
        dbco3[n,4]=nrow(dbtmp)
      }
    }
  }
}

dbco4=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:1){
    for (k in (j+1):2){
      for (m in (k+1):3){
        for (p in (m+1):4){
          dbtmp$sum1=dbtmp[,j]+dbtmp[,k]+dbtmp[,m]+dbtmp[,p]
          n=n+1
          dbco4[n,1]=levels(dbmat$subcluster)[i]
          dbco4[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k],colnames(dbtmp)[m],colnames(dbtmp)[p])
          dbco4[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp[,m]>0 & dbtmp[,p]>0 & dbtmp$sum1==dbtmp$sum,])
          dbco4[n,4]=nrow(dbtmp)
        }
      }
    }
  }
}


dbco=rbind(dbco1,dbco2,dbco3,dbco4)
dbco$ratio=dbco$num/dbco$totalnum
dbco$genes=factor(dbco$genes,levels = unique(dbco$genes))
dbco$subcluster=factor(dbco$subcluster,levels = levels(dbmat$subcluster))

library(ggplot2)
p2=ggplot(dbco, aes(genes,subcluster)) +
  geom_tile(aes(fill = ratio)) + 
  scale_fill_gradient(low = "white", high = "#DC6018") +theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2



#=======================================marmoset
filename='Marmaset_GRIA4_vs_other_ampa1.2.3.4'

#去掉14 类
seurat_Marmaset_GLu=subset(seurat_Marmaset_GLu, cluster_id != '13' )

#exp
seurat_Marmaset_GLu_exp=as.matrix(seurat_Marmaset_GLu@assays$SCT@data)[c('GRIA4','GRIA3','GRIA2','GRIA1','CACNG2','CACNG3','CACNG4'),]
#cell type
meta=seurat_Marmaset_GLu@meta.data
meta$subclass_label=ifelse(meta$cluster_id=='0','GRIA4_cell','other')#seurat_Marmaset_GLu$field

dbmat=as.data.frame(t(seurat_Marmaset_GLu_exp))
dbmat$subcluster=meta$subclass_label
dbmat$region=meta$MainType
dbm=melt(dbmat)
dbm$subcluster=factor(dbm$subcluster,levels = unique(dbm$subcluster))#[c(12,19,17,9,8,24,15,6,3,4,5,13,1,7,21,10,22,20,16,18,14,11,31,30,29,28,26,25,27,23,2,32,34,37,35,36,33)])


dbmat$subcluster=factor(dbmat$subcluster,levels = levels(dbm$subcluster))
dbco1=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:4){
    dbco1[((i-1)*4+j),1]=levels(dbmat$subcluster)[i]
    dbco1[((i-1)*4+j),2]=colnames(dbtmp[j])
    dbco1[((i-1)*4+j),3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,j]==dbtmp$sum,])
    dbco1[((i-1)*4+j),4]=nrow(dbtmp)
  }
}

dbco2=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:3){
    for (k in (j+1):4){
      dbtmp$sum1=dbtmp[,j]+dbtmp[,k]
      n=n+1
      dbco2[n,1]=levels(dbmat$subcluster)[i]
      dbco2[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k])
      dbco2[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp$sum1==dbtmp$sum,])
      dbco2[n,4]=nrow(dbtmp)
    }
  }
}

dbco3=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:2){
    for (k in (j+1):3){
      for (m in (k+1):4){
        dbtmp$sum1=dbtmp[,j]+dbtmp[,k]+dbtmp[,m]
        n=n+1
        dbco3[n,1]=levels(dbmat$subcluster)[i]
        dbco3[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k],colnames(dbtmp)[m])
        dbco3[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp[,m]>0 & dbtmp$sum1==dbtmp$sum,])
        dbco3[n,4]=nrow(dbtmp)
      }
    }
  }
}

dbco4=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:1){
    for (k in (j+1):2){
      for (m in (k+1):3){
        for (p in (m+1):4){
          dbtmp$sum1=dbtmp[,j]+dbtmp[,k]+dbtmp[,m]+dbtmp[,p]
          n=n+1
          dbco4[n,1]=levels(dbmat$subcluster)[i]
          dbco4[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k],colnames(dbtmp)[m],colnames(dbtmp)[p])
          dbco4[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp[,m]>0 & dbtmp[,p]>0 & dbtmp$sum1==dbtmp$sum,])
          dbco4[n,4]=nrow(dbtmp)
        }
      }
    }
  }
}

dbco=rbind(dbco1,dbco2,dbco3,dbco4)
dbco$ratio=dbco$num/dbco$totalnum
dbco$genes=factor(dbco$genes,levels = unique(dbco$genes))
dbco$subcluster=factor(dbco$subcluster,levels = levels(dbmat$subcluster))

library(ggplot2)

p2=ggplot(dbco, aes(genes,subcluster)) +
  geom_tile(aes(fill = ratio)) + 
  scale_fill_gradient(low = "white", high = "red") +theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2


dbco2=read.table(text = "",col.names = c("subcluster","gene1","gene2","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:7]
  for (j in 1:6){
    for (k in (j+1):7){
      n=n+1
      dbco2[n,1]=levels(dbmat$subcluster)[i]
      dbco2[n,2]=colnames(dbtmp)[j]
      dbco2[n,3]=colnames(dbtmp)[k]
      dbco2[n,4]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0,])
      dbco2[n,5]=nrow(dbtmp[dbtmp[,j]>0 | dbtmp[,k]>0,])
    }
  }
}
dbco2$ratio=dbco2$num/dbco2$totalnum
#dbco2$gene1=apply(dbco2,1,function(x){strsplit(x[2]," ")[[1]][1]})
#dbco2$gene2=apply(dbco2,1,function(x){strsplit(x[2]," ")[[1]][2]})
dbco2$subcluster=factor(dbco2$subcluster,levels = levels(dbmat$subcluster))

p3=ggplot(dbco2, aes(gene1, gene2)) +
  geom_tile(aes(fill = ratio),colour="black") + 
  geom_text(aes(label = round(ratio, 2)),size=3)+
  scale_fill_gradient(low = "white", high = "red") +theme_minimal()+theme(axis.text.x = element_text(angle = 45,hjust = 1),panel.grid.major = element_blank())
p3=p3+facet_wrap(~subcluster)
p3



maar_name=openxlsx::read.xlsx('/work/dy/project/temp/singlecell/hippocampus/merge_new_macaque_marmoset/name_for_clsuter/rename2 的副本.xlsx',
                              sheet='marmoset')
maar_name=maar_name[grep('Glu', maar_name$subtype),]
maar_name$subtype2=sub(' \\(CA)| \\(DG)','',maar_name$subtype)
maar_name$subtype2=gsub(' ','_',maar_name$subtype2)
maar_name$subtype2=gsub('Glu_','',maar_name$subtype2)
maar_name$new=paste(maar_name[,2],maar_name[,3],maar_name[,4],sep='_')


filename='marmoset_GRIA4_vs_other_ampa1.2.3.4-new-names-sub'

seurat_Marmaset_GLu_exp=as.matrix(seurat_Marmaset_GLu@assays$SCT@data)[c('GRIA4','GRIA3','GRIA2','GRIA1','CACNG2','CACNG3','CACNG4'),]
#cell type
meta=seurat_Marmaset_GLu@meta.data

#新名字
meta$subtype =paste(meta$cluster_id,meta$cluster_gene,sep='_')
meta=inner_join(meta,maar_name,by=c('subtype'='subtype2'))

meta$subclass_label= meta$new #ifelse(meta$cluster_id=='4','GRIA4_cell','other')#seurat_Macaca_GLu$field
 #ifelse(meta$cluster_id=='4','GRIA4_cell','other')#seurat_Macaca_GLu$field

 
dbmat=as.data.frame(t(seurat_Marmaset_GLu_exp))
dbmat$subcluster=meta$subclass_label
dbmat$region=meta$MainType

#==提取 sub 类群
dbmat=dbmat[grep('SUB',dbmat$subcluster),]

dbm=melt(dbmat)
dbm$subcluster=factor(dbm$subcluster,levels = unique(dbm$subcluster))


dbmat$subcluster=factor(dbmat$subcluster,levels = levels(dbm$subcluster))
dbco1=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:4){
    dbco1[((i-1)*4+j),1]=levels(dbmat$subcluster)[i]
    dbco1[((i-1)*4+j),2]=colnames(dbtmp[j])
    dbco1[((i-1)*4+j),3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,j]==dbtmp$sum,])
    dbco1[((i-1)*4+j),4]=nrow(dbtmp)
  }
}

dbco2=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:3){
    for (k in (j+1):4){
      dbtmp$sum1=dbtmp[,j]+dbtmp[,k]
      n=n+1
      dbco2[n,1]=levels(dbmat$subcluster)[i]
      dbco2[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k])
      dbco2[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp$sum1==dbtmp$sum,])
      dbco2[n,4]=nrow(dbtmp)
    }
  }
}

dbco3=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:2){
    for (k in (j+1):3){
      for (m in (k+1):4){
        dbtmp$sum1=dbtmp[,j]+dbtmp[,k]+dbtmp[,m]
        n=n+1
        dbco3[n,1]=levels(dbmat$subcluster)[i]
        dbco3[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k],colnames(dbtmp)[m])
        dbco3[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp[,m]>0 & dbtmp$sum1==dbtmp$sum,])
        dbco3[n,4]=nrow(dbtmp)
      }
    }
  }
}

dbco4=read.table(text = "",col.names = c("subcluster","genes","num","totalnum"))
n=0
for (i in 1:length(levels(dbm$subcluster))){
  dbtmp=dbmat[dbmat$subcluster==levels(dbmat$subcluster)[i],1:4]
  dbtmp$sum=apply(dbtmp,1,sum)
  for (j in 1:1){
    for (k in (j+1):2){
      for (m in (k+1):3){
        for (p in (m+1):4){
          dbtmp$sum1=dbtmp[,j]+dbtmp[,k]+dbtmp[,m]+dbtmp[,p]
          n=n+1
          dbco4[n,1]=levels(dbmat$subcluster)[i]
          dbco4[n,2]=paste(colnames(dbtmp)[j],colnames(dbtmp)[k],colnames(dbtmp)[m],colnames(dbtmp)[p])
          dbco4[n,3]=nrow(dbtmp[dbtmp[,j]>0 & dbtmp[,k]>0 & dbtmp[,m]>0 & dbtmp[,p]>0 & dbtmp$sum1==dbtmp$sum,])
          dbco4[n,4]=nrow(dbtmp)
        }
      }
    }
  }
}

dbco=rbind(dbco1,dbco2,dbco3,dbco4)
dbco$ratio=dbco$num/dbco$totalnum
dbco$genes=factor(dbco$genes,levels = unique(dbco$genes))
dbco$subcluster=factor(dbco$subcluster,levels = levels(dbmat$subcluster))

library(ggplot2)

p2=ggplot(dbco, aes(genes,subcluster)) +
  geom_tile(aes(fill = ratio)) + 
  scale_fill_gradient(low = "white", high = "red") +theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2	














			
			
			
			










