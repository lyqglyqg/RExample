step1: 导入数据

# install.packages('remotes')
#remotes::install_github(repo = 'genecell/COSGR')


rm(list=ls())
getwd()
options(stringsAsFactors = F) 
library(COSG)
library(harmony)
library(ggsci)
library(dplyr) 
library(future)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggsci)
#BiocManager::install("ggsci",ask = F,update = F)


setwd("/media/lyqg/WD_BLACK/000000_YanmingLi_Work_Z820/000000复现/非肿瘤之scRNA-seq和bulk RNA结合思路")

###### step1:导入数据 ######   

## 解压缩到 GSE181279_raw/， 名字如“GSM5494107_AD1_barcodes.tsv.gz”， 改为“GSM5494107-AD1_barcodes.tsv.gz”
library(stringr)
fs = list.files('./GSE181279_raw/',pattern = '^GSM')
#执行这一步需要解压tar -xvf
samples=str_split(fs,'_',simplify = T)[,1]
samples
 [1] "GSM5494107-AD1" "GSM5494107-AD1" "GSM5494107-AD1" "GSM5494110-AD2"
 [5] "GSM5494110-AD2" "GSM5494110-AD2" "GSM5494113-AD3" "GSM5494113-AD3"
 [9] "GSM5494113-AD3" "GSM5494116-NC1" "GSM5494116-NC1" "GSM5494116-NC1"
[13] "GSM5494119-NC2" "GSM5494119-NC2" "GSM5494119-NC2"



lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste0("GSE181279_raw/", str_split(y[1],'_',simplify = T)[,1])
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("GSE181279_raw/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("GSE181279_raw/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE181279_raw/",y[3]),file.path(folder,"matrix.mtx.gz"))
})

dir='./GSE181279_raw/'

samples=list.files( dir )
samples
# samples = head(samples,10) 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  sce =CreateSeuratObject(counts =  Read10X(file.path(dir,pro )) ,
                          project =  gsub('^GSM[0-9]*_','',pro)  ,
                          min.cells = 5,
                          min.features = 500 )
  return(sce)
})
#gsub函数： gsub (匹配内容，替换内容，操作对象）
samples
a=gsub('^GSM[0-9]*_','',samples)
a
b=gsub('_gene_cell_exprs_table.txt.gz','',gsub('^GSM[0-9]*_','',samples))
b

names(sceList) 
# gsub('^GSM[0-9]*','',samples)
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =  gsub('_gene_cell_exprs_table.txt.gz','',gsub('^GSM[0-9]*_','',samples) )     )

as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 



step2: QC质控
###### step2:QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")


input_sce<-sce.all
  #计算线粒体基因比例
  mito_genes=rownames(input_sce)[grep("^MT-", rownames(input_sce),ignore.case = T)] 
  print(mito_genes) #可能是13个线粒体基因
  #input_sce=PercentageFeatureSet(input_sce, "^MT-", col.name = "percent_mito")
  input_sce=PercentageFeatureSet(input_sce, features = mito_genes, col.name = "percent_mito")
  fivenum(input_sce@meta.data$percent_mito)#fivenum函数：返回五个数据：最小值、下四分位数、中位数、上四分位数、最大值。
  
  #计算核糖体基因比例
  ribo_genes=rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce),ignore.case = T)]
  print(ribo_genes)#可能是100个核糖体基因
  input_sce=PercentageFeatureSet(input_sce,  features = ribo_genes, col.name = "percent_ribo")
  fivenum(input_sce@meta.data$percent_ribo)
  
  #计算红血细胞基因比例
  Hb_genes=rownames(input_sce)[grep("^Hb[^(p)]", rownames(input_sce),ignore.case = T)]
  print(Hb_genes)#可能是7个红血细胞基因
  input_sce=PercentageFeatureSet(input_sce,  features = Hb_genes,col.name = "percent_hb")
  fivenum(input_sce@meta.data$percent_hb)
  
  #可视化细胞的上述比例情况
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  p1 
  w=length(unique(input_sce$orig.ident))/3+5;w
  ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
    scale_y_continuous(breaks=seq(0, 100, 5)) +
    NoLegend()
  p2  
  w=length(unique(input_sce$orig.ident))/2+5;w
  ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)
  
  p3=FeatureScatter(input_sce, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
  ggsave(filename="Scatterplot.pdf",plot=p3)



#根据上述指标，过滤低质量细胞/基因
  #过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
  # 一般来说，在CreateSeuratObject的时候已经是进行了这个过滤操作
  # 如果后期看到了自己的单细胞降维聚类分群结果很诡异，就可以回过头来看质量控制环节
  # 先走默认流程即可
#  if(F){
    selected_c <- WhichCells(input_sce, expression = (nFeature_RNA > 200 & nFeature_RNA < 2500))
    selected_f <- rownames(input_sce)[Matrix::rowSums(input_sce@assays$RNA@counts > 0 ) > 3]
    input_sce.filt <- subset(input_sce, features = selected_f, cells = selected_c)
    dim(input_sce) 
    [1] 16382 36252
    dim(input_sce.filt) 
    [1] 16382 35265
#  }
#  input_sce.filt =  input_sce




## par(mar = c(4, 8, 2, 1))
#  # 这里的C 这个矩阵，有一点大，可以考虑随抽样 
#  C=subset(input_sce.filt,downsample=100)@assays$RNA@counts # downSample将随机抽样一个数据集
#  dim(C)
#  C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
#
#  most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
#  
# pdf("TOP50_most_expressed_gene.pdf",width=21,height = 21)
#  boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
#          cex = 0.1, las = 1, 
#          xlab = "% total count per cell", 
#          col = (scales::hue_pal())(50)[50:1], 
#          horizontal = TRUE)
#  dev.off()
#  rm(C)




#过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
  selected_mito <- WhichCells(input_sce.filt, expression = percent_mito < 5)
#  selected_ribo <- WhichCells(input_sce.filt, expression = percent_ribo > 3)
#  selected_hb <- WhichCells(input_sce.filt, expression = percent_hb < 1 )
#  length(selected_hb)
#  length(selected_ribo)
  length(selected_mito)
  [1] 24026
  
  input_sce.filt <- subset(input_sce.filt, cells = selected_mito)
#  input_sce.filt <- subset(input_sce.filt, cells = selected_ribo)
#  input_sce.filt <- subset(input_sce.filt, cells = selected_hb)
  dim(input_sce.filt)#[1] 16382 24026

  table(input_sce.filt$orig.ident) 
GSM5494107-AD1 GSM5494110-AD2 GSM5494113-AD3 GSM5494116-NC1 GSM5494119-NC2 
          5317           4964           3152           6074           4519 




#可视化过滤后的情况
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  p1_filtered
  w=length(unique(input_sce.filt$orig.ident))/3+5;w 
  ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)
  
  feats <- c("percent_mito")
  p2_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 1) + 
    NoLegend()
  p2_filtered
  w=length(unique(input_sce.filt$orig.ident))/2+5;w 
  ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 

 setwd('../')





step3: harmony整合多个单细胞样品
###### step3: harmony整合多个单细胞样品 ######
dir.create("2-harmony")
setwd("2-harmony")
getwd()
#source('../scRNA_scripts/harmony.R')
# 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"


input_sce<-input_sce.filt
  print(dim(input_sce))
  
  input_sce <- NormalizeData(input_sce, 
                             normalization.method = "LogNormalize",
                             scale.factor = 1e4) 
 
  input_sce <- FindVariableFeatures(input_sce)
  input_sce <- ScaleData(input_sce)
  input_sce <- RunPCA(input_sce, features = VariableFeatures(object = input_sce))
  seuratObj <- RunHarmony(input_sce, "orig.ident")
  names(seuratObj@reductions)
  seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                       reduction = "harmony") #umap
  #seuratObj=RunTSNE(seuratObj,  dims = 1:15, reduction = "harmony") #tsne
  p = DimPlot(seuratObj,reduction = "umap",label=T ) 
  p #看整合效果
#  ggsave(filename='umap-by-orig.ident-after-harmony.pdf',plot = p)





#####
input_sce=seuratObj
  input_sce <- FindNeighbors(input_sce, reduction = "harmony",
                             dims = 1:15) 
  input_sce.all=input_sce
  
  #设置不同的分辨率，观察分群效果(选择哪一个？)
  for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
    input_sce.all=FindClusters(input_sce.all, #graph.name = "CCA_snn", 
                               resolution = res, algorithm = 1)
  }
  colnames(input_sce.all@meta.data)
  apply(input_sce.all@meta.data[,grep("RNA_snn",colnames(input_sce.all@meta.data))],2,table)



DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2")








step4: 细胞亚群生物学命名
# 需要自行看图,定细胞亚群：
celltype=data.frame(ClusterID=0:7,
                    celltype= 0:7) 

#定义细胞亚群 
celltype[celltype$ClusterID %in% c( 0 ),2]='Memory CD4 T1'  
celltype[celltype$ClusterID %in% c( 1),2]='Memory CD4 T2'   
celltype[celltype$ClusterID %in% c( 2),2]='NKT' 
celltype[celltype$ClusterID %in% c( 3),2]='NK'  
 
celltype[celltype$ClusterID %in% c( 4),2]='B Cell'  
celltype[celltype$ClusterID %in% c( 5),2]='DC'
celltype[celltype$ClusterID %in% c( 6),2]='CD8 T' 
celltype[celltype$ClusterID %in% c( 7),2]='Platelet'


head(celltype)
celltype
table(celltype$celltype)

sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(input_sce.all@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)

       B Cell         CD8 T            DC Memory CD4 T1 Memory CD4 T2 
         1870           183           249          8421          6312 
           NA            NK           NKT      Platelet 
        12226          2047          4889            55 

       B Cell         CD8 T            DC Memory CD4 T1 Memory CD4 T2 
         1870           183           249          8421          6312 
           NK           NKT      Platelet 
         2047          4889            55 


## #去除干扰亚群
table(sce.all$celltype)
       B Cell         CD8 T            DC Memory CD4 T1 Memory CD4 T2 
         1870           183           249          8421          6312 
           NK           NKT      Platelet 
         2047          4889            55 

#sce=sce.all
##sce <- sce[, !(sce$celltype %in% c("24","47"))]
#table(sce$celltype)





可视化
#sce.all=input_sce.all

#th=theme(axis.text.x = element_text(angle = 45, 
#                                    vjust = 0.5, hjust=0.5)) 
#library(patchwork)
#p_all_markers=DotPlot(sce.all, features = unique(genes_to_check),
#                      assay='RNA', group.by = 'celltype' )  + coord_flip()+th
p_umap=DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T,label.box = T)
p_umap
ggsave('umap_by_celltype.pdf',width = 8,height = 8)
ggsave('umap_by_celltype.png',width = 8,height = 8)

#p_all_markers+p_umap
#ggsave('markers_umap_by_celltype.pdf',width = 18,height = 8)



DimPlot(sce.all, reduction = "umap",split.by = 'orig.ident',
        group.by = "celltype",label = T) 
ggsave('umap_by_celltype.each_subj.pdf',width = 40, height = 8)
ggsave('umap_by_celltype.each_subj.png',width = 40, height = 8)


table(input_sce.filt$orig.ident) 

#GSM5494107-AD1 GSM5494110-AD2 GSM5494113-AD3 GSM5494116-NC1 GSM5494119-NC2 
#          5317           4964           3152           6074           4519 

sce.all@meta.data$disease = "NA"

sce.all@meta.data[which(sce.all@meta.data$orig.ident == "GSM5494107-AD1"),'disease'] <- "AD"
sce.all@meta.data[which(sce.all@meta.data$orig.ident == "GSM5494110-AD2"),'disease'] <- "AD"
sce.all@meta.data[which(sce.all@meta.data$orig.ident == "GSM5494113-AD3"),'disease'] <- "AD"
sce.all@meta.data[which(sce.all@meta.data$orig.ident == "GSM5494116-NC1"),'disease'] <- "NC"
sce.all@meta.data[which(sce.all@meta.data$orig.ident == "GSM5494119-NC2"),'disease'] <- "NC"

  
table(sce.all@meta.data$disease)
   AD    NC 
13433 10593 


DimPlot(sce.all, reduction = "umap",split.by = 'disease',
        group.by = "celltype",label = T) 

ggsave('umap_by_celltype.ADvsNC.pdf',width = 40, height = 8)





在健康滑膜和PTOA滑膜中，成纤维细胞、髓系细胞和内皮细胞在数量上占主导地位，其中成纤维细胞和髓系细胞在ACLR 7天的丰度增加最急剧。还检测到罕见的细胞群，包括T细胞和Schwann cells，这些细胞在滑膜中没有很好的特征。



TSNEPlot(sce.all, reduction = "tsne",split.by = 'disease',
        group.by = "celltype",label = T) 


  seuratObj.1 <- RunTSNE(seuratObj,  dims = 1:15, 
                       reduction = "harmony") #tsne
  #seuratObj=RunTSNE(seuratObj,  dims = 1:15, reduction = "harmony") #tsne
  p = TSNEPlot(seuratObj.1,reduction = "tsne",label=T ) 



#####
input_sce=seuratObj.1
  input_sce <- FindNeighbors(input_sce, reduction = "harmony",
                             dims = 1:15) 
  input_sce.all=input_sce
  
  #设置不同的分辨率，观察分群效果(选择哪一个？)
  for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
    input_sce.all=FindClusters(input_sce.all, #graph.name = "CCA_snn", 
                               resolution = res, algorithm = 1)
  }
  colnames(input_sce.all@meta.data)
  apply(input_sce.all@meta.data[,grep("RNA_snn",colnames(input_sce.all@meta.data))],2,table)



TSNEPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.2")
ggsave('TSNE_by_celltype.nonanno.pdf')



step4: 细胞亚群生物学命名
# 需要自行看图,定细胞亚群：
celltype=data.frame(ClusterID=0:7,
                    celltype= 0:7) 

#定义细胞亚群 
celltype[celltype$ClusterID %in% c( 0 ),2]='Memory CD4 T1'  
celltype[celltype$ClusterID %in% c( 1),2]='Memory CD4 T2'   
celltype[celltype$ClusterID %in% c( 2),2]='NKT' 
celltype[celltype$ClusterID %in% c( 3),2]='NK'  
 
celltype[celltype$ClusterID %in% c( 4),2]='B Cell'  
celltype[celltype$ClusterID %in% c( 5),2]='DC'
celltype[celltype$ClusterID %in% c( 6),2]='CD8 T' 
celltype[celltype$ClusterID %in% c( 7),2]='Platelet'


head(celltype)
celltype
table(celltype$celltype)

input_sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  input_sce.all@meta.data[which(input_sce.all@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(input_sce.all@meta.data$celltype)

      B Cell         CD8 T            DC Memory CD4 T1 Memory CD4 T2 
         1870           183           249          8421          6312 
           NK           NKT      Platelet 
         2047          4889            55 


p_tsne=TSNEPlot(input_sce.all, reduction = "tsne", group.by = "celltype",label = T,label.box = T)
p_tsne
ggsave('TSNE_by_celltype.pdf',width = 8,height = 8)
ggsave('TSNE_by_celltype.png',width = 8,height = 8)



input_sce.all@meta.data$disease = "NA"

input_sce.all@meta.data[which(input_sce.all@meta.data$orig.ident == "GSM5494107-AD1"),'disease'] <- "AD"
input_sce.all@meta.data[which(input_sce.all@meta.data$orig.ident == "GSM5494110-AD2"),'disease'] <- "AD"
input_sce.all@meta.data[which(input_sce.all@meta.data$orig.ident == "GSM5494113-AD3"),'disease'] <- "AD"
input_sce.all@meta.data[which(input_sce.all@meta.data$orig.ident == "GSM5494116-NC1"),'disease'] <- "NC"
input_sce.all@meta.data[which(input_sce.all@meta.data$orig.ident == "GSM5494119-NC2"),'disease'] <- "NC"

  
table(input_sce.all@meta.data$disease)
   AD    NC 
13433 10593 


TSNEPlot(input_sce.all, reduction = "tsne",split.by = 'disease',
        group.by = "celltype",label = T) 

ggsave('TSNE_by_celltype.ADvsNC.pdf',width = 10, height = 5)

###############################################################################
###############################################################################

C=subset(input_sce.all,downsample=100)@assays$RNA@counts # downSample将随机抽样一个数据集
  dim(C)
  [1] 16382  1873
  C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100

  most_expressed <- order(apply(C, 1, median), decreasing = T)[sample(c(1:16382), 50)]

most_expressed
 [1] 13378  8646 10607 10127 15033  9436  8967 12791  7658  2399 15308  6949
[13]  1184  8283 13407 10178  3403 12027   957  8230 15315  3850 15165  1982
[25]  4065  9768 13484  9876  4389   223  7892 15931 11689 15502 10220  3063
[37]  7445 11090 15677  8192 12390 11161  8428 12191 13681 14944   555 10964
[49]  4692 11150


gene_to_check <- rownames(C)[most_expressed]
> gene_to_check
 [1] "FRG1B"         "PTPRJ"         "MED6"          "ALOX5AP"      
 [5] "Z83844.1"      "RP11-996F15.2" "ARRB1"         "BRIP1"        
 [9] "WDR34"         "AC005037.3"    "RSPH1"         "RPL7"         
[13] "UAP1"          "RP11-348N5.7"  "EIF2S2"        "LCP1"         
[17] "RTP4"          "RNF166"        "MRPS21"        "NPM3"         
[21] "AP001046.5"    "ANXA5"         "CERK"          "RNF181"       
[25] "RAI14"         "METTL25"       "TTPAL"         "UNG"          
[29] "SEPT8"         "NBPF3"         "PLXDC2"        "AC005786.7"   
[33] "ZNF768"        "AC092620.3"    "THSD1"         "RP11-221J22.2"
[37] "NFIL3"         "PDCD7"         "FAM86B2"       "ZFYVE27"      
[41] "TRAF4"         "COX5A"         "TOLLIP-AS1"    "ASGR2"        
[45] "RPS15"         "CHEK2"         "GPX7"          "SERF2"        
[49] "BPHL"          "UBL7-AS1"     




th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
#library(patchwork)
p_all_markers=DotPlot(sce.all, features = unique(gene_to_check),
                      assay='RNA', group.by = 'celltype' )  + coord_flip()+th

ggsave('dotplot.pdf',width = 10, height = 5)




