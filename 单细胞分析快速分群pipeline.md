# 单细胞分析快速分群pipeline

> 杨晓龙 2023-11-2
>
> 该流程封装了大量seurat代码为数个函数，能极大幅度精简分析速度，未经允许请勿外传。

## 数据读取

``` R
library(SeuratDisk)
library(Seurat)
library(writexl) #excel处理包
library(ggplot2)
#LoadLoom路径不能有中文，不然报错
##直接h5ad转换成h5seurat 
##overwrite参数：覆盖源文件
#Convert("./peiwai/res1.5/", dest = "h5seurat", overwrite = F)
load("./data/CS7.blood.rds")
seurat_obj <- readRDS()
seurat_obj <- LoadH5Seurat("./cs7_res_tunning/res2/bak.h5seurat") #本流程推荐的读取方式

#如果数据里本来就有分群数据，则可选启用下面这段代码
#seurat_obj@meta.data$clusters <- clusters$clusters
#seurat_obj@meta.data$name1 <- seurat_obj@meta.data$name
#originalname <- 'name1'
```

## 高变基因检索

``` R
huoshan <- function(seurat_obj, dir){
  # 如果目录不存在，创建目录
  if(!dir.exists(dir)){
    dir.create(dir)
  }
  
  #标准化；寻找在不同细胞中高变异的基因绘制火山图；归一化。
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj,nfeatures = 3000)
  
  top10_ss2 <- head(VariableFeatures(seurat_obj),10) #前面十个高变基因
  #归一化
  #标准化是对列进行的，归一化是按照行进行的
  seurat_obj <- ScaleData(seurat_obj)
  #默认的时候scale只做2000个高表达变异基因，有需要也可以用参数换成所有基因
  
  #画火山图看高变基因
  plot10x1 <- VariableFeaturePlot(seurat_obj)
  plot10x2 <- LabelPoints(plot = plot10x1, points = top10_ss2, repel = F)
  pdf(paste0(dir, "/huoshan.pdf"), width = 16, height = 9) #输出16*10的pdf
  print(plot10x1 + plot10x2)
  dev.off()
  return(seurat_obj)
}
seurat_obj <- huoshan(seurat_obj, 'huoshan') #如果不赋值则不会保存
```

## PCA,UMAP

``` R
#归一化，PCA，UMAP画图
seurat_obj <- ScaleData(seurat_obj, verbose = T)
seurat_obj <- RunPCA(seurat_obj,npcs = 30,verbose = T)
# Plot PCA
DimPlot(seurat_obj, reduction = "pca")
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
ElbowPlot(seurat_obj, ndims = 30)
#从图中看出前30个PC就够了
seurat_obj <- RunUMAP(seurat_obj,reduction = "pca",dims = 1:20 ) #，min.dist= ,spread=这几个参数是循环迭代的核心，可以写个循环函数来调整UMAP的。
seurat_obj <- FindNeighbors(seurat_obj,reduction = "pca",dims = 1:20)#默认20
seurat_obj <- FindClusters(seurat_obj,resolution = 0.1)#resolution是分辨率。
#记得要保存
# saveRDS(seurat_obj, file = "seurat_obj.rds")
# seurat_obj <- readRDS("seurat_obj.rds")
```

## 变量定义

``` R
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') 
#如果数据里本来就有分群数据，则可选启用下面这段代码
# resname <- 'seurat_clusters'
```

## 自选基因dotplot图绘制函数

> 这个函数可以输入任意一组gene，它可以实现全自动绘制并将基因按照每个cluster的表达情况呈阶梯排序，不需要进行人工排序。

``` R
#画dotplot并且自动梯度排序gene                           
plot_dot_sl <- function(seurat_obj, genes, save.path){
  avg.expression <- AverageExpression(seurat_obj, features = genes)
  avg.expression <- as.data.frame(avg.expression)
  
  get_sort_genes <- function(data){
    # 计算每个基因的最高表达值对应的cluster
    max_exp_cluster <- apply(data, 1, function(x) names(which.max(x)))
    
    # 按照cluster归类基因
    cluster_genes <- split(names(max_exp_cluster), max_exp_cluster)
    
    return(cluster_genes)
  }
  
  combine_unique_genes <- function(cluster_genes){
    # 对cluster名进行排序
    sorted_names <- order(as.numeric(gsub("RNA.", "", names(cluster_genes))))
    
    # 按照排序结果将基因拼接成向量并去重
    unique_genes <- rev(unique(unlist(cluster_genes[sorted_names])))
    
    return(unique_genes)
  }
  
  unique_genes <- combine_unique_genes(get_sort_genes(avg.expression))
  
  dplot <- DotPlot(object = seurat_obj, features =  unique_genes, dot.scale = 4,scale.min=0 ,scale.max = 100) + theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5)) +
    scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7F7F7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")) +
    labs(x=NULL,y=NULL) +
    guides(size=guide_legend(order=3)) +
    theme(legend.position = "bottom",legend.direction =  "horizontal") +scale_size_continuous (range = c (1.5,5.5))
    
  ggsave(save.path, plot = dplot, width = 10, height = 7,limitsize = FALSE)
  
  cat(paste0("dotplot保存路径:", save.path, "\n"))
}

#使用方法
 plot_dot_sl(seurat_obj= , genes= unique(genes), save.path= '')
```

## 分群核心函数

> 生成大量数据

``` R
plot_others_save <- function(seurat_obj, resultdir) {
  # 创建保存图片的文件夹
  if (!dir.exists(resultdir)) {
    dir.create(resultdir)
  }
  
  # 保存umap图
  print("Saving UMAP plot...")
  #UMAPplot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 10, cols = c("#532C8A","#F7901D","#B51D8D","#F397C0","#C594BF", "#DFCDE4","#A0CC47","#3F84AA","#B3793B","#683612", "#C72228","#EF4E22", "#989898","#333333","#7F6874", "#7253A2","#65A83E","#EF5A9D","#FBBE92","#139992", "#C9EBFB","#8DB5CE","#CE4E82", "#354E23","#77783C", "#8EC792","#0F4A9C","#FACB12","#BBDCA8","#1A1A1A", "#C3C388","#DABE99","#005579", "#CDE088","#FFF574","#F6BFCB", "gray"),pt.size = 1.5,group.by = c(resname,originalname)) #resname来自下一个函数，originalname自己赋值一下
      UMAPplot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 6, cols = my36colors,pt.size = 1.5,group.by = c(resname,originalname)) #group.by 参数可选
  save.name <- file.path(resultdir, "UMAP_plot.pdf")
  ggsave(save.name, UMAPplot, width = 12, height = 6)
  cat(paste0("UMAP图片保存路径为：", save.name, "\n"))
  
  # # 画肘图并保存
  # print("Saving Elbow Plot...")
  # ElbowPlot(seurat_obj, ndims = 30)
  # save.name <- file.path(resultdir, "ElbowPlot.pdf")
  # ggsave(save.name, width = 10, height = 10)
  # cat(paste0("肘图图片保存路径为：", save.name, "\n"))
  
  # 寻找差异基因并保存csv
  print("Saving significant markers...")
  # 平均表达量
  avg.expr <- AverageExpression(seurat_obj)
  # 找到每个聚类的显著基因
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # 将markers保存为CSV文件
  save.name <- file.path(resultdir, "significant_markers.csv")
  write.csv(markers, save.name, row.names = FALSE)
  cat(paste0("差异基因csv保存路径为：", save.name, "\n"))
  library(readxl)
  save.name <- file.path(resultdir, "significant_markers.xlsx")
  write_xlsx(markers,path = save.name)
    
  library(dplyr)
markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top10
  hmap <- DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
  save.name <- file.path(resultdir, "top20_heatmap.pdf")
  ggsave(save.name, plot = hmap, width = 18, height = 18)
  cat(paste0("heatmap保存路径:", save.name, "\n"))
  
  save.name <- file.path(resultdir, "top20_dotplot.pdf")
  unique_genes <- unique(top10$gene)
dplot <- DotPlot(object = seurat_obj, features = unique_genes, dot.scale = 4,scale.min=0 ,scale.max = 100) + theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7F7F7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")) +
  labs(x=NULL,y=NULL) +
  guides(size=guide_legend(order=3)) +
  theme(legend.position = "bottom",legend.direction =  "horizontal") +scale_size_continuous (range = c (1.5,5.5))
ggsave(save.name, plot = dplot, width = 25, height = 8,limitsize = FALSE)
  cat(paste0("dotplot保存路径:", save.name, "\n"))
  save.name <- file.path(resultdir, "sl_dotplot.pdf")
   #调用自选基因点图函数(可选)
#  plot_dot_sl(seurat_obj= seurat_obj, genes = c('PCNA', 'ITGA4', 'PCLAF', 'MYB', 'GAD1', 'CD52', 'LIN28A', 'RGS16', 'THY1', 'CD34', 'HAND1', 'BMP4', 'KDR', 'CDH5', 'KIT', 'CD44', 'PROCR', 'ITGA2B', 'TEK', 'MEF2C', 'PECAM1', 'GJA4', 'GJA5', 'TAL1', 'GYPB', 'ALAS2', 'LMO2', 'KLF1', 'GP1BB', 'TREM2', 'NFE2', 'GATA1', 'PF4', 'GP1BB', 'TREM2', 'NFE2', 'GATA1', 'PF4', 'PPBP', 'GP9', 'CMTM5','VWF', 'PLEK', 'RUNX1', 'SLC2A3', 'DDIT4', 'CD14', 'CD68', 'PKIB', 'PTPRC', 'ADGRG1', 'SPI1', 'HBZ', 'HBA1', 'HBA2', 'GYPA', 'HBG2', 'HBE1', 'HBM', 'HBG1', 'TTR', 'APOE', 'APOA1', 'FGB','LYVE1','FOLR2'), save.path = save.name)
    
  # 保存每个聚类的细胞数目为csv
  print("Saving cluster cell count...")
  table_seurat <- table(seurat_obj@meta.data$seurat_clusters)
  save.name <- file.path(resultdir, "cluster_cell_count.csv")
  write.csv(table_seurat, save.name, row.names = FALSE)
  # write.csv(table_seurat, paste(save.name, '.xlsx', seq= ''), row.names = FALSE)
  cat(paste0("聚类细胞数目csv保存路径：", save.name, "\n"))
  #library(readxl)
  table_seurat <- as.data.frame(table_seurat)
  save.name <- file.path(resultdir, "cluster_cell_count.xlsx")
  write_xlsx(table_seurat,path = save.name)
  
  print("开始存储数据备份")
  save.name <- file.path(resultdir, "bak.h5seurat")
  SaveH5Seurat(seurat_obj,filename=save.name, overwrite = TRUE)
  #数据转为最终h5ad格式
  Convert(save.name, dest = "h5ad", overwrite = TRUE)
  print("执行结束")
}
plot_others_save(seurat_obj, 'UMAP')

resultdir <- './peiwai_HQ1/res2'
## 对指定范围内的分辨率，分别执行函数算出所有所需的图片

resolution_tunning <- function(seurat_obj, minres = 2.2, maxres = 3.2, dirpath = 'cs7_res_tunning') {
  if (!dir.exists(dirpath)) {
    dir.create(dirpath)
  }
  
  for (i in seq(minres, maxres, by = 0.2)) {
    print(paste("当前分辨率为：", i))
    seurat_obj <- FindClusters(seurat_obj, resolution = i) # 修改了分辨率参数
    subdir <- file.path(dirpath, paste0("res", i)) # 新建子文件夹
    resname <- paste0("RNA_snn_res.", i, sep='')
    if (!dir.exists(subdir)) {
      dir.create(subdir)
    }
    plot_others_save(seurat_obj, file.path(subdir)) # 修改存储路径和文件名
  }
}
#使用方法
plot_others_save(seurat_obj, 'UMAP')
```

## 分辨率调整函数

> 调用核心函数，生成不同分辨率下的大量图片与表格

``` R
resolution_tunning <- function(seurat_obj, minres = 2.2, maxres = 3.2, dirpath = 'cs7_res_tunning') {
  if (!dir.exists(dirpath)) {
    dir.create(dirpath)
  }
  
  for (i in seq(minres, maxres, by = 0.2)) {
    print(paste("当前分辨率为：", i))
    seurat_obj <- FindClusters(seurat_obj, resolution = i) # 修改了分辨率参数
    subdir <- file.path(dirpath, paste0("res", i)) # 新建子文件夹
    resname <- paste0("RNA_snn_res.", i, sep='')
    if (!dir.exists(subdir)) {
      dir.create(subdir)
    }
    plot_others_save(seurat_obj, file.path(subdir)) # 修改存储路径和文件名
  }
}
#使用方法
resolution_tunning(seurat_obj, minres = 2.2, maxres = 3.2, dirpath = 'cs7_res_tunning')
```

## DEGS表格优化函数

需要自己构建一个my_markers.csv,如图：这个表格和DEGs是算法实现的基础

<img src="https://picgo-temporary-storage-yxl.oss-cn-hangzhou.aliyuncs.com/Yxl-Picgo-temporary/image-20231102104253137.png" alt="image-20231102104253137" style="zoom:33%;" />

``` R
# 加载需要的包
library(tidyverse)
#library(readxl)
#library(readr)
library(writexl)

remove_spaces_from_genes <- function(gene_string) {
  genes <- unlist(strsplit(gene_string, split = ","))
  cleaned_genes <- trimws(genes, which = "both") # 去除每个gene名称的前后空格
  return(paste(cleaned_genes, collapse = ","))
}
# 创建一个函数
my_gene_in_DEGs <- function(marker_filepath, DEG_filepath, output_dir) {

  # 读取数据
  my_markers <- read.csv(marker_filepath, stringsAsFactors = FALSE)
  DEGs <- read.csv(DEG_filepath, stringsAsFactors = FALSE)
#如果没有提供output_dir参数，使用DEG_filepath的路径
  if(output_dir == "") {
      output_dir <- dirname(DEG_filepath) 
  } else { 
      output_dir <- paste0(output_dir, "/") 
  }
  # 创建一个空的列表来存储结果
  result_list <- list()

  # 创建循环遍历my_markers中的每个组织类型
  for(i in colnames(my_markers)[-1]){
    # 提取特定的marker基因
    markers <- my_markers %>% pull(!!sym(i))
    # 删除空值
    markers <- markers[markers != ""]

    # 创建一个空的数据框用于存储结果
    res <- data.frame(cluster = integer(), num = integer(), gene = character())

    # 遍历DEGs中的每个cluster
    for(j in unique(DEGs$cluster)){
      # 提取特定cluster中的基因
      genes <- DEGs %>% filter(cluster == j) %>% pull(gene)

      # 计算交集基因数量及其名字
      common <- intersect(markers, genes)
      num <- length(common)

      # 添加到结果数据框中
      res <- rbind(res, data.frame(cluster = j, num = num, gene = toString(common)))
    }

    # 按照num降序排列，并只取前x个
    res <- res %>% arrange(desc(num)) %>% head(40) #灵活调整

    # 保存到结果列表
    result_list[[i]] <- res
  }

  # 输出结果到指定文件夹内
  for(i in names(result_list)){
    write_xlsx(result_list[[i]], paste0(output_dir,'/', i, "_my_gene_in_DEGs.xlsx"))
    # result_list[[i]]$gene <- sapply(result_list[[i]]$gene, remove_spaces_from_genes)
    # write.csv(result_list[[i]], paste0(output_dir,'/', i, "_my_gene_in_DEGs.csv"))
  }
}

# 调用函数
#my_gene_in_DEGs("./my_markers.csv", ".//res1/significant_markers.csv", "")


##搜索路径下需要的文件并执行功能的函数
find_and_run <- function(dir = '', name = '') {
  #先判断路径是否存在
  if (!file.exists(dir)) {
    cat('Directory does not exist. Please check your path.\n')
    return()
  }
  
  #列出所有在dir路径下的文件
  files <- list.files(path = dir, recursive = TRUE, full.names = TRUE)

  #使用正则匹配文件名包含name字符串的文件
  matched_files <- grep(name, files, value = TRUE)
  
  #打印匹配到的文件路径
  if (length(matched_files) > 0) {
    for (i in (1:length(matched_files))) {
    cat('Running Cammand...Using file:\n')
    print(matched_files[i])
    my_gene_in_DEGs("my_markers.csv", matched_files[i], "") }
  } else {
    cat('No files matched.\n')
  }
}

find_and_run(dir = './', name = 'significant_markers.csv')


```

