library(Seurat) #3.2.3
library(viridis)
library(magrittr)
library(dplyr)
library(ggplot2)


library(hrbrthemes)

library(patchwork)


library(loomR)

library(ComplexHeatmap)

sample = "PLA-8w-RNA"



##reload original object and cluster.df.add
placenta <- readRDS( "PLA-8w-RNA.final.rds")
cluster.df.add <- readRDS('cluster.df.add.rds')



#modified for CTB with dark red colors
color_snap_mod1 = c('1'='#777711','2'='#E31A1C','3'='#68228B','4'='#771122','5'='grey','6'='#1F78B4','7'='#FFD700','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')

color_gradient_my <- c(
    rgb(5,48,97,maxColorValue = 255),
    rgb(42,113,178,maxColorValue = 255),
    rgb(147,198,222,maxColorValue = 255),
    rgb(239,243,245,maxColorValue = 255),
    rgb(253,219,199,maxColorValue = 255),
    rgb(214,96,77,maxColorValue = 255),
    rgb(121,5,34,maxColorValue = 255)

)

color_cellranger <-c('#820610','#C50F1E','#F42428','#F86D30','#FBB33D','#FCFB4E','#C0FB61','#87FB8F','#41FEF9','#2BAED7','#155CB1','#08238D') #for depth


color_pyreds <- c(
"#fff5f0","#fef4ef","#fef3ee","#fef3ed","#fef2ec","#fef1eb","#fef1ea","#fef0e9","#feefe8","#feefe7","#feeee6","#feede5","#feede4","#feece3","#feebe2","#feebe1","#feeae0","#fee9e0","#fee9df","#fee8de"
)


####read in color palette list###
color_set_yellowbrick.flat <- readRDS('color_set_yellowbrick.flat.rds')
color_list_archr <- readRDS('ArchR.color_list.rds')


options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_list_archr)){ 
  len = length(color_list_archr[[name]])
  color = color_list_archr[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}

options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_set_yellowbrick.flat)){ 
  len = length(color_set_yellowbrick.flat[[name]])
  color = color_set_yellowbrick.flat[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}



#####get cluster.df again and compare
# cluster <- Idents(placenta)
# umap <- Embeddings(placenta,reduction = 'umap') #placenta_filter@reductions$umap@cell.embeddings
# all.equal(names(cluster),rownames(umap)) #TRUE


# cluster.df <- data.frame(cluster=cluster,umap)

# metadata <- placenta@meta.data
# all.equal(rownames(cluster.df),rownames(metadata))#TRUE

# cluster.df.add <- cbind(cluster.df, metadata)

# all.equal(cluster.df.add,readRDS('cluster.df.add.rds') ) #TRUE



###quick plot embedding and check
DimPlot(placenta,reduction = 'umap',label = TRUE, label.size = 10)


##filter cluster, rotate umap coordinate, label cluster , set cluster color

table(cluster.df.add$cluster)

   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74 


###filter STR and cluster 7

cluster.df.add <- subset(cluster.df.add, ! cluster %in% c('12','10','13','14','15','16') ) #10202 #9392
table(cluster.df.add$cluster)


   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
1720 1428 1225 1093 1083  963  893  853  656    0  288    0    0    0    0    0 

# 0    1    2    3    4    5    6    7    8    9   10   11 
# 2372 1757 1609 1580 1026  497  287    0  165   99    0    0 


#####################the final filtering object###################
idx <- which (rownames(placenta@meta.data) %in% rownames(cluster.df.add)     )
table(rownames(placenta@meta.data) %in% rownames(cluster.df.add)  )
FALSE  TRUE 
 1004 10202 

FALSE  TRUE 
  352  9392 


placenta <- subset(x = placenta, cells = rownames(cluster.df.add) )
placenta
24307 features across 10202 samples within 1 assay 
Active assay: RNA (24307 features, 2000 variable features)
 3 dimensional reductions calculated: pca, umap, tsne

# 22355 features across 9392 samples within 1 assay 
# Active assay: RNA (22355 features, 2000 variable features)
#  3 dimensional reductions calculated: pca, umap, tsne

all.equal(rownames(placenta@meta.data),rownames(cluster.df.add) ) #TRUE
all.equal(colnames(placenta@assays$RNA@data),rownames(cluster.df.add) ) #TRUE

#drop levels

table(Idents(placenta))

   1    2    3    4    5    6    7    8    9   11 
1720 1428 1225 1093 1083  963  893  853  656  288 


DimPlot(placenta,reduction = 'umap',label = TRUE, label.size = 10)+NoLegend()


#rename cluster id
mapid <- list('1'='1','2'='2','3'='3', '4'='4','5'='5','6'='6','7'='7',
'8'='8',
'9'='9',
'11'='10'
 )

names(mapid)

cluster <- as.character(Idents(placenta) )
names(cluster) <- names(Idents(placenta))
for(i in 1:length(cluster)){ cluster[i] = mapid[[ cluster[i] ]] }

cluster <- factor(cluster,levels=gtools::mixedsort( names(table(cluster)) ))

table(cluster)
cluster
   1    2    3    4    5    6    7    8    9   10 
1720 1428 1225 1093 1083  963  893  853  656  288 


placenta <- SetIdent(object = placenta,value =  cluster )
table(Idents(placenta))


##quick look the filtered umap embedding and cluster
DimPlot(placenta,reduction = 'umap',label = TRUE, label.size = 10)+NoLegend()


##add cell name and cell color
map_cellname <- list(
    '1'='STB1',
    '2'='STB5', 
    '3'='STB2',
    '4'='STB3',
    '5'='CTB2',
    '6'='STB naive',
    '7'='STB4',
    '8'='CTB1',
    '9'='CTB3',
    '10'='Fusion component'
 )





reds <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf')
blues <- c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','darkblue')




oranges <- c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6')
purples <- c('#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')

map_cellcolor <- list(
    '1'=purples[2],#'STB1',
    '2'=purples[3],#'STB5', 
    '3'=purples[5],#'STB2',
    '4'=purples[2],#'STB3',
    '5'=blues[5],#'CTB-2',
    '6'=purples[1],#naive STB',
    '7'=purples[4],#'STB4',
    '8'=blues[3],#'lightblue',#'CTB-1',
    '9'=blues[4],#'CTB-3',
    '10'='darkgreen'#'Fusion component'
    
#     '1'=purples[4], #STB4
#     '2'=purples[5],  #STB5
#     '3'=purples[2],  #STB2
#     '4'=purples[1], #STB1
#     '5'=purples[3], #STB3
#     '6'='#8B0000', #Syncytial knot
#     '7'='#d8daeb', #naive STB
#     '8'='#7f3b08', #STB-new
#     '9'='darkgreen', #CTB
#     '10'=''
 )


##do some filtering for UMAP and add the centered UMAP

cluster.df.add <- subset(cluster.df.add,UMAP_1 < 8) #10198

placenta <- subset(x = placenta,cells = rownames(cluster.df.add) )
placenta
#24307 features across 10198 samples

all.equal(rownames(cluster.df.add),rownames(placenta@meta.data) ) #TRUE





##get tne final cluster.df.add
umap <- Embeddings(placenta,reduction = 'umap')
cluster <- Idents(placenta)
cluster.df = data.frame(cluster=cluster,UMAP_1=umap[,1],UMAP_2=umap[,2]) 
rownames(cluster.df) = rownames(placenta@meta.data)

all.equal (rownames(cluster.df),rownames(placenta@meta.data) )#TRUE
cluster.df.add <- cbind(cluster.df,placenta@meta.data)


all.equal(Idents(placenta),cluster.df.add$cluster,check.attributes = FALSE)#TRUE

cellname <- as.character(Idents(placenta))
cellcolor<- as.character(Idents(placenta))
for(i in 1:length(cellname)){ cellname[i] = map_cellname[[ cellname[i] ]] }
for(i in 1:length(cellcolor)){ cellcolor[i] = map_cellcolor[[ cellcolor[i] ]] }

cluster.df.add[,'cellname'] <- factor(cellname)
cluster.df.add[,'cellcolor'] <- factor(cellcolor)




all.equal(Idents(placenta),cluster.df.add$cluster,check.attributes = FALSE) #TRUE


#####save/reload the final object#######
#saveRDS(placenta,"placenta.final.final.rds")     
placenta <- readRDS("placenta.final.final.rds") #TE10198

#write.table(cluster.df,file='snapATAC.PLA-8w-RNA.umap.bin5k.final.txt',sep='\t',quote = FALSE, row.names = TRUE, col.names = TRUE)

#saveRDS(cluster.df.add,"cluster.df.add.final.rds")

cluster.df.add <- readRDS("cluster.df.add.final.rds") #TE10198


# ##transform to loom format (will save to file), better create loom in a standalone r script, or this file will locked

# lf <- as.loom(x = placenta,verbose=TRUE,filename = 'seurat2loom/PLA-8w-RNA.TE10198.new1.loom',overwrite = TRUE)#,filename = '',assay='')

# lf$close_all #?


###how many CTB nuclei in 8w data?### 
subset(cluster.df.add, cluster %in% c('8','5','9','10'))
##2878



######plot the final fixed cluster, color and label with customized plot function
# centers <- cluster.df.add %>% dplyr::group_by(cellname) %>% dplyr::summarize(x = median(x = UMAP_1), 
#         y = median(x = UMAP_2))

centers <- cluster.df.add %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
       y = median(x = UMAP_2))

centers_shift = data.frame()
##add little shift for text x and y, to plot text halo, borrow from snapATAC
theta= seq(0, 2*pi, length.out=50)
r=0.1
strwidth = 0.5 
strheight = 0.5
xo <- r*strwidth # r*strwidth('A')
yo <- r*strheight #r*strheight('A')
for (i in seq_len(nrow(centers))){
  for (j in theta) {
        centers_shift = rbind(centers_shift,
                              data.frame(
                                  ##cluster=as.character(unlist(centers[i,'cellname'])),
                                  cluster=as.character(unlist(centers[i,'cluster'])),
                                  x=centers[i,'x'] + cos(j)*xo, 
                                  y=centers[i,'y'] + sin(j)*yo
                                 )
                       )
      }
}

####the UMAP plot with annotation


##the UMAP range
UMAP = Embeddings(placenta,reduction = 'umap')
all.equal(as.data.frame(UMAP),cluster.df.add[,c('UMAP_1','UMAP_2')],check.attributes = FALSE) #TRUE

left <- 2.1*min(UMAP[,1])
right <- 2.1*max(UMAP[,1])
top <- 1*max(UMAP[,2])
bottom <- 1*min(UMAP[,2])


##label on cluster
options(repr.plot.height=5,repr.plot.width=5.5,repr.plot.res = 150)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
##ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cellname  )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = color_good)  +
  scale_colour_manual(values = unlist(map_cellcolor) )  +
  ##scale_colour_manual(values = color_snap_mod1)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  theme(
        legend.position = 'none',
        #axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        #axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "two donors, total cells:",nrow(cluster.df.add),  sep=" ") ) +
  geom_text(data = centers_shift, #the halo
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "white", 
            ##size = 4.5) +
            size = 6.5) +
  geom_text(data = centers, 
            ##mapping = aes(x=x,y=y,label = cellname), 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            ##size = 4.5) +
            size = 6.5) +
  ##guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  #ylim(0,1) + 
  xlim(left,right) + ylim(bottom,top) + ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/PLA-8w-RNA-UMAP.labelon.pdf",height=5,width=5.5,useDingbats=FALSE)
##ggsave(filename = "pdfs/PLA-8w-RNA-UMAP.labelon.clusterid.pdf",height=5,width=5.5,useDingbats=FALSE)
##ggsave(filename = "pdfs/PLA-8w-RNA-UMAP.labeloff.pdf",height=5,width=5.5,useDingbats=FALSE)
###############



##by donors
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= sample  )) +
  geom_point(size = .1,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = c('red','pink'))  + #'#94C6DD', '#1273AE'; 'red','navy'
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Source of donor ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  xlim(left,right) + ylim(bottom,top) +
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/PLA-8w-RNA-source-of-donor.pdf",height=5.5,width=5,useDingbats=FALSE)




###by depth
#options(repr.plot.height=5.5,repr.plot.width=5)
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= log(nCount_RNA) )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = rev(color_cellranger))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'D'))  + 
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Sequence depth (log)",  sep=" ") ) +
  #ggtitle(paste(sample, "Sequence depth ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  xlim(left,right) + ylim(bottom,top) +
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/PLA-8w-RNA-depth.pdf",height=5.5,width=5,useDingbats=FALSE)


###by gene captured
#options(repr.plot.height=5.5,repr.plot.width=5)
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=  nFeature_RNA)) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = rev(color_cellranger))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'D'))  + 
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Gene Captured",  sep=" ") ) +
  #ggtitle(paste(sample, "Sequence depth ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  xlim(left,right) + ylim(bottom,top) +
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/PLA-8w-RNA-gene-captured.pdf",height=5.5,width=5,useDingbats=FALSE)



#####by cell cycling stage: plot layers of Phase score ##########
cluster.df.add.G1 <- subset(cluster.df.add,Phase == 'G1') #6246
cluster.df.add.G2M <- subset(cluster.df.add,Phase == 'G2M') #1950
cluster.df.add.S <- subset(cluster.df.add,Phase == 'S') #3010

cluster.df.add.phase <- cluster.df.add.G1
stage_name <- 'G1'
color_stage <- '#282F76'
outpdf <- "pdfs/PLA-8w-RNA-cellcycling-G1.pdf"

cluster.df.add.phase <- cluster.df.add.G2M
stage_name <- 'G2M'
color_stage <- 'brown'
outpdf <- "pdfs/PLA-8w-RNA-cellcycling-G2M.pdf"

cluster.df.add.phase <- cluster.df.add.S
stage_name <- 'S'
color_stage <- 'darkgreen'
outpdf <- "pdfs/PLA-8w-RNA-cellcycling-S.pdf"


options(repr.plot.height=5,repr.plot.width=5)
ggplot(cluster.df.add.phase,aes(x=UMAP_1,y=UMAP_2,col= Phase  )) +
  #geom_point(size = 0.1,alpha= 0.5,col='#282F76' ) +
  ##geom_point(data=cluster.df.add.G1,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col='#282F76' ) +
  geom_point(data=cluster.df.add.phase,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col= color_stage ) +
  ##geom_point(data=cluster.df.add.S,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col='darkgreen' ) +
  #scale_colour_manual(values = c('#282F76','brown','darkgreen'))  +
  #scale_colour_manual(values = c('red','navy','orange'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "cell cycling phase ", stage_name , sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = outpdf,height=5.5,width=5)


##statistic for cellcycling with barplot
res.cellcycling <- table(cluster.df.add$Phase,cluster.df.add$cluster)
res.cellcycling <- as.matrix(res.cellcycling )
res.cellcycling <- res.cellcycling[,c(8,5,9,10,6,4,7,1,3)]

res.cellcycling.perc <- 100*t(t(res.cellcycling)/colSums(res.cellcycling))

pal.cellcycling <- c('G1'='#282F76', 'G2M'='brown', 'S'='darkgreen')

pdf('pdfs/PLA-8w-RNA-cellcycling.barplot.pdf',height=5.5,width=5)
options(repr.plot.width = 5, repr.plot.height = 5.5)
#par(xpd=TRUE)
barplot( res.cellcycling.perc, main = 'Statitics of cell cycling state', ylab = 'Percentage (%)' ,xlab = 'Clusters',col=pal.cellcycling )
legend(10,120,legend = c('G1','G2M','S'), fill = pal.cellcycling, xpd = TRUE )
dev.off()

####cellcycling done




# ## Heretical clustering ###
# # calculate the ensemble signals for each cluster
# ensemble.ls = lapply(split(seq(length(cluster.df.add$cluster)), cluster.df.add$cluster), function(x){
# 	Matrix::colMeans(GetAssayData(placenta,slot = 'data')[,x]);
# 	})
# # cluster using 1-cor as distance  
# hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");

# pdf( "pdfs/hclust.cluster.pdf",height=5.5,width=4,useDingbats = FALSE)
# options(repr.plot.height=5.5,repr.plot.width=4)
# plot(hc, hang=-1, xlab="");
# dev.off()


# ####bin ave coverage(depth) of cluster (bmat)
# ensemble.ls.df <- do.call(cbind,ensemble.ls) #620094  
# #all.equal(t(do.call(rbind, ensemble.ls)),do.call(cbind,ensemble.ls)) #TRUE
# ensemble.ls.df <- ensemble.ls.df[rowSums(ensemble.ls.df) != 0,] #568573 #575002
# ensemble.ls.df

# options(repr.plot.height=5,repr.plot.width=5)
# boxplot(ensemble.ls.df,outline = FALSE,col = color_snap_mod1,las=2)


# ##bin sum coverage(depth) of cluster (bmat)
# ensembleSum.ls = lapply(split(seq(length(placenta@cluster)), placenta@cluster), function(x){
# 	SnapATAC::colSums(placenta[x,], mat="bmat");
# 	})
# ensembleSum.ls.df <- do.call(cbind,ensembleSum.ls) #620094

# ensembleSum.ls.df <- ensembleSum.ls.df[rowSums(ensembleSum.ls.df) != 0,] #568573 #574062
# ensembleSum.ls.df

# options(repr.plot.height=5,repr.plot.width=5)
# boxplot(ensembleSum.ls.df,outline = FALSE,col = color_snap_mod1,las=2) #quick

# ##bin sum coverage(depth) of cluster (bmat) by donor
# #all.equal(cluster.df.add[,-c(1,2,3)],placenta@metaData) #TRUE
# ensembleSum.ls = lapply(split(seq(length(placenta@cluster)),  #group by donor and cluster
#                               list(cluster=placenta@cluster,
#                                    lib=placenta@metaData$library)
#                              ), 
#                         function(x){SnapATAC::colSums(placenta[x,], mat="bmat")}
#                        )
# ensembleSum.ls.df <- do.call(cbind,ensembleSum.ls) #620094
# saveRDS(ensembleSum.ls.df,'ensembleSum.ls.df.rds')


# ensembleSum.ls.df.d1 <- ensembleSum.ls.df[,grep(pattern = 'D1',x=colnames(ensembleSum.ls.df))]
# ensembleSum.ls.df.d2 <- ensembleSum.ls.df[,grep(pattern = 'D2',x=colnames(ensembleSum.ls.df))]



# ensembleSum.ls.df.d1 <- ensembleSum.ls.df.d1[rowSums(ensembleSum.ls.df.d1) != 0,] #566361 #572203
# ensembleSum.ls.df.d1
# ensembleSum.ls.df.d2 <- ensembleSum.ls.df.d2[rowSums(ensembleSum.ls.df.d2) != 0,] #565418 #568099
# ensembleSum.ls.df.d2


# pdf( "pdfs/boxplot_depth_by_donor.pdf",height=7.5,width=10,useDingbats = FALSE)
# par(mfrow = c(2,1),mar=c(3,3,3,1))
# #options(repr.plot.height=10,repr.plot.width=5)
# boxplot(ensembleSum.ls.df.d1,outline = FALSE,col = color_snap_mod1,las=2,ylim=c(0,40),main='donor1') #quick
# boxplot(ensembleSum.ls.df.d2,outline = FALSE,col = color_snap_mod1,las=2,ylim=c(0,40),main='donor2') #quick
# dev.off()



######stat d1 d2 and cluster cell number correlation######
res.stat <- table(cluster.df.add$cluster,cluster.df.add$sample) 
    
      D1  D2
  1  992 728
  2  827 601
  3  712 513
  4  669 424
  5  645 438
  6  529 434
  7  473 420
  8  544 309
  9  406 250
  10 185 103



as.data.frame(res.stat) %>% dplyr::group_by(Var2) %>% summarise(sum = sum(Freq))

  D1   D2 
5982 4216



res.cor <- cor(res.stat) #0.924 #0.977 #0.983 #0.988
res.stat[,1] = -1 * res.stat[,1]

res.stat.df <- as.data.frame(res.stat)
#res.stat.df <- reshape2::melt(res.stat)
colnames(res.stat.df) <- c('cluster','sample','count')

##plot horizonal barplot  #  '#94C6DD', '#1273AE'  vs '#F3A585','#C80927'
options(repr.plot.height=5.5,repr.plot.width=4)
ggplot(res.stat.df, aes(fill=sample, y=count, x=cluster )) + 
    #geom_hline(yintercept = c(0,25,50,75),linetype='solid',size=.3,col='black') +
    #geom_hline(yintercept = c(0),linetype='solid',size=1,col='black') +
    #geom_bar(position="stack", stat="identity",alpha=1,width = 0.5) +  
    geom_bar(position=position_stack(reverse=TRUE), stat="identity",width = 0.5) +  
    #xlim(100,0) +
    #scale_x_continuous(breaks = seq(100,0,-25),labels=seq(100,0,-25)) +
    #scale_y_reverse() +
    #annotate('text',x = 16, y = 1000,label = paste('spearman cor=',cor(res.stat)[1,2]) ) +
    coord_flip() +
    #scale_fill_viridis(discrete = T,option = "E") +
    scale_fill_manual(values = c('red', 'pink'),labels=c('D1','D2'),name='sample' ) +
    #ggtitle("cell number count") +
    labs(title = "cell number count", subtitle=paste('spearman cor=',round(res.cor[1,2],digits=2)) ) +
    theme_ipsum(base_family = 'sans') + #to avoid Arial narrow problem in ggsave
    ylab("count")

ggsave(filename="pdfs/donor1.vs.donor2.cluster.cor.pdf",width = 4, height = 5.5,useDingbats=FALSE)





#######




##### annatation with know genes #######
count <- GetAssayData(object = placenta,slot = 'counts')
data <- GetAssayData(object = placenta,slot = 'data') #24307 x 10198
scale.data <- GetAssayData(object = placenta,slot = 'scale.data') 

#saveRDS(count,'exprMat.count.rds')
#saveRDS(data,'exprMat.data.rds')
#saveRDS(scale.data,'exprMat.scaledata.rds')


#data <- GetAssayData(placenta,slot = 'data')
data <- as(data,'matrix')

data <- as(scale.data,'matrix')


marker.genes = c(
    "DNMT1", "CDH1", "MKI67",
    "FLT1", "CSHL1", "PSG8", 
    "ERVFRD-1", "LAIR2", "PLAC8",
    'PECAM1','CD14'
  );


marker.genes = c( #try to distinguish STB terminals 
    "FLT1", 'LEP','INSIG2', #FLT1-enriched early stage
    "CSHL1",'CSH1','PAPPA',  #PAPPA-enriched late stage
    'PSG1','CGA','GCM1' #general

  );

##general STB: PSG1,PSG8,CGA,;naive STB: ERVFRD-1,GCM1; hormone-bias STB:PAPPA,CSHL1;FLT1, ENG,LEP, INSIG2
#CTB DNMT1, CDH1
#Syncytial knot, CDKN1A SPATA18, CROT, PTCHD4
#STB-new:  DDX60,LIFR,MAPKAPK3,VAC14
#total 9 or 12 subplot
marker.genes = c(
    'PSG1','PSG8','CGA','CGB3',
     'PAPPA','CSHL1', 'CSH1','CSH2',
     'FLT1', 'ENG','LEP','INSIG2',
    'ERVFRD-1','GCM1','DNMT1', 'CDH1',
     'TEAD3','TEAD4','TEAD1','GATA3',
    'CDKN1A', 'SPATA18', 'CROT', 'PTCHD4',
    'DDX60','LIFR','MAPKAPK3','VAC14'
)

marker.genes <- c('PAPPA', 'FLT1', 'LEP', 'PSG8', 'PTCHD4', 'SH3TC2')
p
marker.genes <- c('DNMT1','ERVFRD-1', 'PSG8', 'CGA', 'LEP','CSHL1' , 'PAPPA','FLT1')


marker.genes <- marker.genes[marker.genes %in% rownames(data)]

marker.genes <- c('SLC19A2','SLC46A1','FOLR1','FOLR2','FOLR3','LRP2','SLC19A1','SLC25A32')

options(repr.plot.height=15,repr.plot.width = 15)
FeaturePlot(placenta, features = marker.genes, reduction = "umap",slot = 'scale.data',cols = c("lightgrey", "darkred") )

######customized plot marker gene scatter plot#########


marker.gmat <- t(as.data.frame (data[marker.genes,]) )
#marker.gmat <- as.data.frame(as(placenta@gmat[, marker.genes],'matrix'))
#rownames(marker.gmat) <- rownames(placenta@metaData)
all.equal (rownames(cluster.df.add),rownames(marker.gmat) ) #TRUE
marker.gmat.df <- cbind(cluster.df.add[,c('cluster','UMAP_1','UMAP_2')],marker.gmat )


##cutoff quantiles
#low.q <- 0.1
#high.q <- 0.9

# low.q <- c('DNMT1'=0.1,'CDH1'=0.1,'PAGE4'=0.1,'FLT1'=0,'CSHL1'=0.1,'PSG8'=0,'ERVFRD-1'=0,'LAIR2'=0,'PLAC8'=0.1,'VIM'=0,'PECAM1'=0,'MKI67'=0,'CD14'=0)
# high.q <- c('DNMT1'=1,'CDH1'=1,'PAGE4'=1,'FLT1'=1,'CSHL1'=1,'PSG8'=1,'ERVFRD-1'=1,'LAIR2'=1,'PLAC8'=1,'VIM'=1,'PECAM1'=1,'MKI67'=1,'CD14'=1)

# low.q <- c('FLT1'=0,'LEP'=0,'INSIG2'=0,'CSHL1'=0,'CSH1'=0,'PAPPA'=0,'PSG1'=0,'CGA'=0,'GCM1'=0)
# high.q <- c('FLT1'=1,'LEP'=1,'INSIG2'=1,'CSHL1'=1,'CSH1'=1,'PAPPA'=1,'PSG1'=1,'CGA'=1,'GCM1'=1)


# ##used for viridis D 256 color
# low.q <- c('PSG1'=0,'PSG8'=0,'CGA'=0,'CGB3'=0,'PAPPA'=0.1,'CSHL1'=0.1,'CSH1'=0.1,'CSH2'=0.1,'FLT1'=0.1,'ENG'=0.1,'LEP'=0.1,'INSIG2'=0.05,'ERVFRD-1'=0.1,'GCM1'=0.05,'DNMT1'=0,'CDH1'=0,'TEAD3'=0.1,'TEAD4'=0,'TEAD1'=0,'GATA3'=0.1,'CDKN1A'=0.1,'SPATA18'=0,'CROT'=0,'PTCHD4'=0,'DDX60'=0.1,'LIFR'=0.05,'MAPKAPK3'=0.1,'VAC14'=0.1)
# high.q <- c('PSG1'=1,'PSG8'=1,'CGA'=1,'CGB3'=1,'PAPPA'=1,'CSHL1'=1,'CSH1'=1,'CSH2'=1,'FLT1'=1,'ENG'=1,'LEP'=1,'INSIG2'=1,'ERVFRD-1'=1,'GCM1'=1,'DNMT1'=1,'CDH1'=1,'TEAD3'=1,'TEAD4'=1,'TEAD1'=1,'GATA3'=1,'CDKN1A'=1,'SPATA18'=1,'CROT'=1,'PTCHD4'=1,'DDX60'=1,'LIFR'=1,'MAPKAPK3'=1,'VAC14'=1)


# ##used for color_gradient_my
# low.q <- c('PSG1'=0,'PSG8'=0,'CGA'=0,'CGB3'=0,'PAPPA'=0.1,'CSHL1'=0.1,'CSH1'=0.1,'CSH2'=0.1,'FLT1'=0.1,'ENG'=0.1,'LEP'=0.1,'INSIG2'=0.05,'ERVFRD-1'=0.1,'GCM1'=0.05,'DNMT1'=0,'CDH1'=0,'TEAD3'=0.1,'TEAD4'=0,'TEAD1'=0,'GATA3'=0.1,'CDKN1A'=0.1,'SPATA18'=0,'CROT'=0,'PTCHD4'=0,'DDX60'=0.06,'LIFR'=0.02,'MAPKAPK3'=0.02,'VAC14'=0.08)
# high.q <- c('PSG1'=1,'PSG8'=1,'CGA'=1,'CGB3'=1,'PAPPA'=1,'CSHL1'=1,'CSH1'=1,'CSH2'=1,'FLT1'=1,'ENG'=1,'LEP'=1,'INSIG2'=1,'ERVFRD-1'=1,'GCM1'=1,'DNMT1'=1,'CDH1'=1,'TEAD3'=1,'TEAD4'=1,'TEAD1'=1,'GATA3'=1,'CDKN1A'=1,'SPATA18'=1,'CROT'=1,'PTCHD4'=1,'DDX60'=1,'LIFR'=1,'MAPKAPK3'=1,'VAC14'=1)

##use for grey-red
low.q <- c('PSG1'=0,'PSG8'=0,'CGA'=0,'CGB3'=0,'PAPPA'=0.1,'CSHL1'=0.1,'CSH1'=0.1,'CSH2'=0.1,'FLT1'=0.1,'ENG'=0,'LEP'=0.1,'INSIG2'=0.05,'ERVFRD-1'=0.1,'GCM1'=0.05,'DNMT1'=0,'CDH1'=0,'TEAD3'=0.1,'TEAD4'=0,'TEAD1'=0,'GATA3'=0.1,'CDKN1A'=0.1,'SPATA18'=0,'CROT'=0,'PTCHD4'=0,'DDX60'=0,'LIFR'=0.5,'MAPKAPK3'=0,'VAC14'=0)
high.q <- c('PSG1'=1,'PSG8'=1,'CGA'=1,'CGB3'=1,'PAPPA'=1,'CSHL1'=1,'CSH1'=1,'CSH2'=1,'FLT1'=1,'ENG'=1,'LEP'=1,'INSIG2'=1,'ERVFRD-1'=1,'GCM1'=1,'DNMT1'=1,'CDH1'=1,'TEAD3'=1,'TEAD4'=1,'TEAD1'=1,'GATA3'=1,'CDKN1A'=1,'SPATA18'=1,'CROT'=1,'PTCHD4'=1,'DDX60'=1,'LIFR'=1,'MAPKAPK3'=1,'VAC14'=1)

##for data
low.q <- c('PAPPA'=0, 'FLT1'=0.5, 'LEP'=0.1, 'PSG8'=0.1, 'PTCHD4'=0, 'SH3TC2'=0.6)
high.q <- c('PAPPA'=1, 'FLT1'=1, 'LEP'=1, 'PSG8'=1, 'PTCHD4'=1, 'SH3TC2'=1)

##for scale.data
low.q <- c('PAPPA'=0, 'FLT1'=0.5, 'LEP'=0.1, 'PSG8'=0.2, 'PTCHD4'=0, 'SH3TC2'=0.5)
high.q <- c('PAPPA'=1, 'FLT1'=1, 'LEP'=1, 'PSG8'=1, 'PTCHD4'=1, 'SH3TC2'=1)

low.q <- c('DNMT1' = 0,'ERVFRD-1'= 0, 'PSG8'= 0.2, 'CGA'= 0.3, 'LEP'= 0.2,'CSHL1'= 0 , 'PAPPA'= 0.3,'FLT1'= 0.5)
high.q <- c('DNMT1' = 1,'ERVFRD-1'= 1, 'PSG8'= 1, 'CGA'= 1, 'LEP'= 1,'CSHL1'= 1 , 'PAPPA'= 1,'FLT1'= 1)




marker.gmat.df.cutoff <- marker.gmat.df[,1:3]
for(gene in marker.genes){
    #cat('gene ',gene)
    feature.value <- marker.gmat.df[,gene]
    #cutoff <- quantile(feature.value,probs = c(0,1) )
    cutoff <- quantile(feature.value,probs = c(low.q[gene],high.q[gene])) 
    cutoff.low <- cutoff[1]
    cutoff.high <- cutoff[2]
    feature.value[feature.value < cutoff.low] <- cutoff.low
    feature.value[feature.value > cutoff.high] <- cutoff.high
    marker.gmat.df.cutoff <- cbind(marker.gmat.df.cutoff,feature.value)
}
colnames(marker.gmat.df.cutoff)[4:ncol(marker.gmat.df.cutoff)] <- marker.genes


######


#color_set_yellowbrick.flat 
#color_list_archr

color_marker_rna <- color_list_archr[['horizonExtra']]
color_marker_rna <- color_set_yellowbrick.flat [['GnBu.6']]
color_marker_rna <- color_set_yellowbrick.flat [['Reds.9']] #use this


#options(repr.plot.height=5.8,repr.plot.width=5)
res.marker <- list()
for(gene in marker.genes){
    p <- ggplot(marker.gmat.df.cutoff,aes_string(x="UMAP_1",y="UMAP_2",col= paste0("`",gene,"`") )) +
    #p <- ggplot(marker.gmat.df,aes_string(x="UMAP_1",y="UMAP_2",col= paste0("`",gene,"`") )) +
      geom_point(size = .2,show.legend = TRUE,alpha= .5 ) +
      #scale_colour_manual(values = c('red','navy'))  +
      #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
      ##scale_colour_gradientn(colors = color_gradient_my )  +
      ##scale_colour_gradientn(colors = color_pyreds )  +
      scale_colour_gradientn(colors = color_marker_rna)  +
      #scale_colour_gradientn(colors = c("lightgrey", "darkred") )  +
      #scale_colour_gradientn(colors = viridis(256, option = "D")) + #good
      ##scale_colour_gradientn(colors = color_tfdev )  +
      #theme_classic() +
      #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
      theme(
            #legend.position = 'none',
            legend.position = 'right',
            axis.text=element_blank(), 
            axis.title = element_text(size = 15, face = "bold"),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(color="black", fill = NA,size=1),
            plot.title = element_text(size = 15, face = "bold"),
            #complete = TRUE
            plot.margin = unit(c(1,1,1,1), "lines")
           )+
      ggtitle(gene) +
      #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
      #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
      #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
      #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
      xlim(left,right) + ylim(bottom,top) +
      labs(x = "UMAP1", y = "UMAP2")
      #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
     ## ggsave(filename = paste('pdfs/marker.expression.plot.',gene,'.pdf',sep=''),height=5,width=5,useDingbats=FALSE)
    ggsave(filename = paste('pdfs/marker.expression.plot.',gene,'.with_legend.pdf',sep=''),height=5,width=6,useDingbats=FALSE)
    res.marker[[gene]] <- p
    #print(p)
}

# ##arrange plot by patchwork
# options(repr.plot.height=15,repr.plot.width=15)
# res.marker[['DNMT1']] + res.marker[['CDH1']] + res.marker[['ERVFRD-1']] + res.marker[['FLT1']] +
# res.marker[['CSHL1']] + res.marker[['PSG8']] + res.marker[['VIM']] + res.marker[['PECAM1']] +
# res.marker[['CD14']] #+ res.marker[['MKI67']] + res.marker[['PLAC8']] + res.marker[['CD14']] +
# plot_layout(ncol=3,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ##ggsave(filename = 'pdfs/marker.gene.common.umap.col.grad.my.pdf',height=13.5,width=18,useDingbats=FALSE) #can save grid multiple plots
# ggsave(filename = 'pdfs/marker.gene.common.umap.col.SnapATAC_like.pdf',height=13.5,width=18,useDingbats=FALSE)

# ##arrange plot by patchwork
# options(repr.plot.height=15,repr.plot.width=15)
# res.marker[['FLT1']] + res.marker[['LEP']] + res.marker[['INSIG2']] + res.marker[['CSHL1']] +
# res.marker[['CSH1']] + res.marker[['PAPPA']] + res.marker[['PSG1']] + res.marker[['CGA']] +
# res.marker[['GCM1']] 
# plot_layout(ncol=3,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ggsave(filename = 'pdfs/marker.gene.hormone.bias.umap.pdf',height=13.5,width=16,useDingbats=FALSE) #can save grid multiple plots

##arrange plot by patchwork
options(repr.plot.height=15,repr.plot.width=15)
#     'PSG1','PSG8','CGA','CGB3',
#      'PAPPA','CSHL1', 'CSH1','CSH2',
#      'FLT1', 'ENG','LEP','INSIG2',
#     'ERVFRD-1','GCM1','DNMT1', 'CDH1',
#      'TEAD3','TEAD4','TEAD1','GATA3',
#     'CDKN1A', 'SPATA18', 'CROT', 'PTCHD4',
#     'DDX60','LIFR','MAPKAPK3','VAC14'

res.marker[['PSG1']] + res.marker[['PSG8']] + res.marker[['CGA']] + res.marker[['CGB3']] +
res.marker[['PAPPA']] + res.marker[['CSHL1']] + res.marker[['CSH1']] + res.marker[['CSH2']] +
res.marker[['FLT1']] + res.marker[['ENG']] + res.marker[['LEP']] + res.marker[['INSIG2']] +
res.marker[['ERVFRD-1']] + res.marker[['GCM1']] + res.marker[['DNMT1']] + res.marker[['CDH1']] 
plot_layout(ncol=4,nrow=4)
#print(res.marker[['FLT1']], vp=viewport(angle=-185))
##ggsave(filename = 'pdfs/marker.gene.hormone.bias.umap.pdf',height=15,width=15,useDingbats=FALSE) #can save grid multiple plots

options(repr.plot.height=12,repr.plot.width=15)
res.marker[['TEAD3']] + res.marker[['TEAD4']] + res.marker[['TEAD1']] + res.marker[['GATA3']] +
res.marker[['CDKN1A']] + res.marker[['SPATA18']] + res.marker[['CROT']] + res.marker[['PTCHD4']] +
res.marker[['DDX60']] + res.marker[['LIFR']] + res.marker[['MAPKAPK3']] + res.marker[['VAC14']] 
plot_layout(ncol=4,nrow=3)



##pick 4 x 3
options(repr.plot.height=12,repr.plot.width=15)
res.marker[['PSG8']] + res.marker[['CGA']] + res.marker[['PAPPA']] + res.marker[['CSHL1']] +
res.marker[['FLT1']] + res.marker[['ENG']] + res.marker[['ERVFRD-1']] + res.marker[['DNMT1']] +
res.marker[['CDKN1A']] + res.marker[['SPATA18']] + res.marker[['MAPKAPK3']] + res.marker[['LIFR']] 
plot_layout(ncol=4,nrow=3)

ggsave(filename = 'pdfs/marker.gene.picked.umap.grey-red.pdf',height=12,width=15,useDingbats=FALSE)


###do for probe gene
pdf(file = "pdfs/marker.expression.combined.pdf",useDingbats = FALSE,height = 10, width = 15,compress = FALSE)
#options(repr.plot.height=10,repr.plot.width=15) #fix 5 x 12 for each row with 3 col
#patchwork::wrap_plots(res.marker, nrow = 2, ncol = 3)
options(repr.plot.height=10,repr.plot.width=20)
patchwork::wrap_plots(res.marker, nrow = 2, ncol = 4)
dev.off()



# ###########customized FeaturePlot for marker genes ####

# res.marker <- list()
# for (gene in marker.genes){
#     options(repr.plot.height=7.5,repr.plot.width=8)
#     p<- FeaturePlot(placenta, features = gene, reduction = "umap",slot = 'data',cols = c("lightgrey", "darkred") )+ #c("lightgrey", "#ff0000", "#00ff00") c("lightgrey", "darkred")
#       #scale_color_gradientn(colours = color_ga)+
#       #scale_color_gradientn(colours = color_peak)+
    
#       #scale_color_gradientn(colours = color_gradient_my)+
#       #scale_color_gradientn(colours = color_pyreds)+
#       theme(
#             legend.position = 'right',
#             axis.text=element_blank(), 
#             axis.title = element_text(size = 15, face = "bold"),
#             axis.ticks = element_blank(),
#             panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             panel.background = element_rect(color="black", fill = NA,size=0.8),
#     #         panel.background = element_rect(fill = "white", colour = "white", 
#     #                 size = rel(1)),
#             #panel.border = element_blank(),
#             plot.title = element_text(size = 15, face = "bold"),
#             #complete = TRUE
#             plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
#            ) 
#     print(p)
#     #ggsave(filename = paste('marker_gene.red_blue',gene,'.pdf',sep=''),height=5,width=5.5,useDingbats=FALSE  )
#      res.marker[[gene]] <- p

# }


# ######control low_threshold and high_threshold for certain gene
# quantile(placenta@assays$RNA@data['PAPPA',],probs = seq(from = 0,to=1,by=0.1) )
# FeaturePlot(placenta, features = 'PAPPA', reduction = "umap",slot = 'data',cols = c("lightgrey", "darkred"),min.cutoff = 4.83577655714442, max.cutoff = 5.63957108331984 )+ #c("lightgrey", "#ff0000", "#00ff00") c("lightgrey", "darkred")
#       #scale_color_gradientn(colours = color_ga)+
#       #scale_color_gradientn(colours = color_peak)+
    
#       #scale_color_gradientn(colours = color_gradient_my)+
#       #scale_color_gradientn(colours = color_pyreds)+
#       theme(
#             legend.position = 'right',
#             axis.text=element_blank(), 
#             axis.title = element_text(size = 15, face = "bold"),
#             axis.ticks = element_blank(),
#             panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             panel.background = element_rect(color="black", fill = NA,size=0.8),
#     #         panel.background = element_rect(fill = "white", colour = "white", 
#     #                 size = rel(1)),
#             #panel.border = element_blank(),
#             plot.title = element_text(size = 15, face = "bold"),
#             #complete = TRUE
#             plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
#            ) 

# ##arrange plot by patchwork
# options(repr.plot.height=13.5,repr.plot.width=18)
# res.marker[['PSG8']] + res.marker[['DNMT1']] + res.marker[['ERVFRD-1']] + res.marker[['PECAM1']] +
# res.marker[['CSHL1']] + res.marker[['CDH1']] + res.marker[['LAIR2']] + res.marker[['VIM']] +
# res.marker[['FLT1']] + res.marker[['MKI67']] + res.marker[['PLAC8']] + res.marker[['CD14']] +
# plot_layout(ncol=4,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ggsave(filename = 'marker.gene.umap.pdf',height=13.5,width=18,useDingbats=FALSE) #can save grid multiple plots








######do findAllMarkers() and visualization######

#check the clusters
options(repr.plot.height = 7.5, repr.plot.width = 7.5)
DimPlot(placenta,label = TRUE,label.size = 15,reduction = 'umap') + NoLegend()

table(Idents(placenta))
  1    2    3    4    5    6    7    8    9   10 
1720 1428 1225 1093 1083  962  892  852  655  288 

##do DEG identification for specific cluster 

##do all DEG identification
placenta.markers <- FindAllMarkers(placenta, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
de.markers.all <- placenta.markers

#saveRDS(object=de.markers.all,file='DEGs/de.markers.all.rds')
de.markers.all <-readRDS(file='DEGs/de.markers.all.rds') #5317 #1342, 
table(de.markers.all$cluster)

#   1    2    3    4    5    6    7    8    9   10 
#  145  132  231  121 1246  137  237 1543  863  662 

#   1   2   3   4   5   6   7   8   9 
#   3 173  14  56  75  62 280 165 514 


de.markers.all$cluster <- factor(de.markers.all$cluster,levels = c('8','5','9','10','6','1','3','2','4','7')  )
table(de.markers.all$cluster)
 8    5    9   10    6    1    3    2    4    7 
1543 1246  863  662  137  145  231  132  121  237 

#    8    5    9   10    6    1    3    4    7    2 
# 1543 1246  863  662  137  145  231  121  237  132 


de.markers.all <- de.markers.all[order(de.markers.all$cluster),]

saveRDS(object=de.markers.all,file='DEGs/de.markers.all.order.rds')




###############get top25 DEG for each cluster##################

#top10$cluster <- factor(top10$cluster,levels = c('8','5','9','10','6','1','3','4','7','2')  )
##top10 <- de.markers.all %>% group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_logFC)
top25 <- de.markers.all %>% group_by(cluster) %>% dplyr::top_n(n = 25, wt = avg_logFC)
top50 <- de.markers.all %>% group_by(cluster) %>% dplyr::top_n(n = 50, wt = avg_logFC)





cluster_order <- factor(Idents(placenta),levels = c('8','5','9','10','6','1','3','2','4','7')  )
placenta <- AddMetaData(placenta,metadata = cluster_order, col.name = 'cluster_order')

all.equal(placenta@meta.data$cluster_order, cluster_order,check.attributes = FALSE ) #TRUE

##save top25
##saveRDS(object=top25,file='DEGs/top25.rds')

##saveRDS(object=top50,file='DEGs/top50.rds')

top25 <- readRDS(file='DEGs/top25.rds')
top50 <- readRDS(file='DEGs/top50.rds')


for(i in c('8','5','9','10','6','1','3','2','4','7') ){
    top50.sel <- subset(top50,cluster == i)
    saveRDS(top50.sel,paste('DEGs/top50_split/top50.c',i,'.rds',sep="") )
    write.table(x = unique(top50.sel$gene),paste('DEGs/top50_split/top50.c',i,'.txt',sep="") ,quote = FALSE, row.names = FALSE, col.name=FALSE )
}

for(i in c('8','5','9','10','6','1','3','2','4','7') ){
    top25.sel <- subset(top25,cluster == i)
    saveRDS(top25.sel,paste('DEGs/top25_split/top25.c',i,'.rds',sep="") )
    write.table(x = unique(top25.sel$gene),paste('DEGs/top25_split/top25.c',i,'.txt',sep=""),quote = FALSE,row.names = FALSE, col.name=FALSE )
}
    
            
# top50.c6 <- subset(top50,cluster == '6')$gene #48 #50
# top50.c8 <- subset(top50,cluster == '8')$gene #41 #50
# top50.c2 <- subset(top50,cluster == '2')$gene #50 #50
# top50.c7 <- subset(top50,cluster == '7')$gene #29 #50
# top50.c4 <- subset(top50,cluster == '4')$gene #31 #37
# dtop50.c5 <- subset(top50,cluster == '5')$gene #44 #50
# top50.c9 <- subset(top50,cluster == '9')$gene #46 #46

# saveRDS(de.markers.c6,'DEGs/de.markers.c6.rds')
# saveRDS(de.markers.c8,'DEGs/de.markers.c8.rds')
# saveRDS(de.markers.c2,'DEGs/de.markers.c2.rds')
# saveRDS(de.markers.c7,'DEGs/de.markers.c7.rds')
# saveRDS(de.markers.c4,'DEGs/de.markers.c4.rds')
# saveRDS(de.markers.c5,'DEGs/de.markers.c5.rds')
# saveRDS(de.markers.c9,'DEGs/de.markers.c9.rds')




##plot DAG heatmap (DoHeatmap, use scale.data by default) in a long heatmap
options(repr.plot.height = 20, repr.plot.width = 12)
DoHeatmap(placenta, 
          #features = top10$gene,
          #features = top25$gene,
          ##features = top25.dedup$gene,
          features = top50$gene,
          slot='scale.data',
          group.by = 'cluster_order',
          group.colors = unlist(map_cellcolor)#,
          #disp.min = -2.5,disp.max=1.5,
         ) +
  scale_fill_gradientn(colors =  color_gradient_my)+
  #scale_fill_gradientn(colors =  c('blue','white','red'))#+
  #scale_fill_gradientn(colors =  viridis(256, option = "D") )#+
  #scale_fill_gradientn(colors =  mapal)#+
  #scale_fill_gradientn(colors =  color_pyreds)#+ #not this
  NoLegend()  

#ggsave(filename = 'DEGs/de.gene.top25.heatmap.pdf',height=20,width=12,useDingbats=FALSE)
ggsave(filename = 'DEGs/de.gene.top25.dedup.heatmap.pdf',height=20,width=12,useDingbats=FALSE)

###




##########################plot aggregate bulk heatmap#############################


###does de gene in more than one cluster?
table(top25$gene %in% rownames(scale.data) ) #250 TRUE
table(duplicated(top25$gene))
FALSE  TRUE 
  220    30 

top25.dedup <- top25[!duplicated(top25$gene ),]

#########how many know markers in these DEG?
marker_genes <- c('DNMT1','ERVFRD-1','LEP','CSHL1','PAPPA','FLT1')
table(marker_genes %in% top25.dedup$gene   ) #all six

subset(top25.dedup,gene %in% marker_genes)[,c('cluster','gene')]
5	DNMT1
10	ERVFRD-1
1	CSHL1
3	PAPPA
4	FLT1
7	LEP

subset(top25,gene %in% marker_genes)[,c('cluster','gene')]
5	DNMT1
9	DNMT1
10	ERVFRD-1
1	CSHL1
3	PAPPA
3	CSHL1
4	FLT1
7	FLT1
7	LEP
2	PAPPA
########


###overlap of fstb de gene list and whole de gene list##

top25.fstb <- readRDS("trajectory_infer/monocle2/DEGs/seurat.top25.fstb.rds") #used for trajectory construction
#175 x 7

table(top25.fstb$gene %in% top25.dedup$gene)
FALSE  TRUE 
   36   139 

top25.fstb$gene[!top25.fstb$gene %in% top25.dedup$gene] #detected in fstb but not in de_gene_all_cluster
'AC011287.1''AZIN1''CDYL2''FRMD6-AS2''CSH1''CYTH3''C1QTNF6''THSD7A''AGAP1''PAG1''SLC27A6''IQGAP2''XDH''DDB1''P4HA1''UBC''AFAP1''SASH1''PTPRM''AC087857.1''EGFR''SLC38A9''DSP''LGR4''ATP10D''UBASH3B''PCDH11X''TIMP3''PCED1B''SLC7A1''GRK3''FLNB''CCND3''KIAA1217''PBX1''TCAF1'


subset(top25.fstb,gene %in% marker_genes)[,c('cluster','gene')]
1	CSHL1
3	PAPPA
3	CSHL1
4	FLT1
7	LEP
10	ERVFRD-1


########

top25.df <- scale.data[top25$gene,]
sum(duplicated( rownames(top25.df) )) #30 duplicated line

sum(duplicated(make.names(rownames(top25.df),unique=TRUE) ) ) #0
rownames(top25.df) <- make.names(rownames(top25.df),unique=TRUE)

aggregateClusters <- function(idents,ids,data){
   data.res <- data.frame(row.names = row.names(data))
   if(ids == "all"){ids = levels(idents)}
   for(id in ids){
     barcodes.sel <- names(idents[idents == id])
     #data.sel <- Matrix::rowSums(data[,barcodes.sel])
     data.sel <- Matrix::rowMeans(data[,barcodes.sel])
     data.res <- cbind(data.res,data.sel)

   }
   colnames(data.res) <- ids
   return(data.res)

}


all.equal(rownames(cluster.df.add),colnames(top25.df) ) #TRUE
cluster <- cluster.df.add$cluster
names(cluster) <- rownames(cluster.df.add)

cluster <- factor(cluster,levels = c('8','5','9','10','6','1','3','2','4','7') )

top25.df.aggre <- aggregateClusters(cluster,'all',top25.df) #already zscored



#quick look the zscore range
options(repr.plot.height = 4.5,repr.plot.width = 4.5)
boxplot(top25.df.aggre,las=2,main="Ave expression aggregated by cluster",xlab="clusters",ylab="peak accessibility",col=unlist(map_cellcolor) ) #cellcolor
abline(h= c(-2,2),lty=2)


for (i in  quantile(unlist(top25.df.aggre),probs = seq(0,1,by = 0.1) ) ){cat(names(i),i,'\n') } ##already zscored
 -1.342801 
 -0.5065482 
 -0.4089049 
 -0.3241194 
 -0.2380356 
 -0.1382739 
 0.0007044859 
 0.1753885 
 0.4381418 
 0.7735507 
 3.119016 


# top25.df.aggre.z <- scale(top25.df.aggre,scale = TRUE)

# for (i in  quantile(unlist(top25.df.aggre.z),probs = seq(0,1,by = 0.1) ) ){cat(names(i),i,'\n') }
# -0.7047654 
#  -0.6066142 
#  -0.5635115 
#  -0.5087242 
#  -0.4445581 
#  -0.362609 
#  -0.2400061 
#  -0.03976284 
#  0.2960676 
#  1.017162 
#  12.54854

#plot aggregate heatmap
options(repr.plot.height = 7.5,repr.plot.width = 5)
#options(repr.plot.height = 5,repr.plot.width = 3)
hp = Heatmap(top25.df.aggre, name = "Z-score", 
             cluster_rows = FALSE, cluster_columns = FALSE, 
             show_row_names = FALSE,show_column_names = TRUE,
             show_column_dend = FALSE,show_row_dend = FALSE,
             use_raster = FALSE,#will use raster if >2000 row or cols, however rstudio do not support raster
             #col = circlize::colorRamp2(seq(-0.6,1.2,by=1.8/10), viridis(n = 11,option = "C")),
             #col = circlize::colorRamp2(seq(-0.2,0.3,by=0.5/10), viridis(n = 11,option = "C")), 
             #col = circlize::colorRamp2(seq(-0.2,0.3,by=0.5/7), color_set_yellowbrick.flat[['RdPu.8']]),
             #col = circlize::colorRamp2(seq(0,0.09,by=0.08/7), color_set_yellowbrick.flat[['RdPu.8']]),
             col = circlize::colorRamp2(seq(-2,1.8,by=3.8/6), color_gradient_my),
             column_title = "DE gene with aggregated expression ",
             #column_split = pdata.order$Cluster,
             #cluster_column_slices = TRUE,
             #show_parent_dend_line = FALSE,
             #column_title = "%s",
             #clustering_method_columns = "complete", ##ward.D,complete
             #clustering_distance_columns  = function (m) dist(m,method="euclidean"), #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #column_gap = unit(1, "mm"),
             border = TRUE,
             heatmap_legend_param = list(legend_direction = 'horizontal',title_position = 'lefttop',border = FALSE)#,grid_height = unit(4, "mm") )
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha,
             #top_annotation = hca
             #heatmap_width=unit(8, "cm")
             
) #%v% HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color_good[1:8]),  #colors   ### + horizonal, %v% : vertical, append annotations independent of heatmap 
                                      ##equal to HeatmapAnnotation()
 #                                     labels = 1:8,#paste0("cluster",1:8),
 #                                     labels_gp = gpar(col = "white", fontsize = 12),
  #                                    which = "column"
                                      
  #        )
  #     )  
#width = max(grobWidth(textGrob(labels))))
pdf('DEGs/DEG_heatmap.aggregated.pdf',height=7.5,width=5,useDingbats=FALSE)
draw(hp, heatmap_legend_side = "bottom",gap = unit(0.1, "cm"))

dev.off()



##########





#####save de gene table#####
write.table(top25,file = "DEGs/top25.de.gene.txt",col.names = TRUE, row.names = FALSE,sep='\t',quote=FALSE)

de.gene.c8 <- subset(top25,cluster == '8')$gene
de.gene.c5 <- subset(top25,cluster == '5')$gene
de.gene.c9 <- subset(top25,cluster == '9')$gene
de.gene.c10 <- subset(top25,cluster == '10')$gene
de.gene.c6 <- subset(top25,cluster == '6')$gene
de.gene.c1 <- subset(top25,cluster == '1')$gene
de.gene.c3 <- subset(top25,cluster == '3')$gene
de.gene.c4 <- subset(top25,cluster == '4')$gene
de.gene.c7 <- subset(top25,cluster == '7')$gene
de.gene.c2 <- subset(top25,cluster == '2')$gene

write.table(de.gene.c8,file = "DEGs/de.gene.c8.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c5,file = "DEGs/de.gene.c5.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c9,file = "DEGs/de.gene.c9.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c10,file = "DEGs/de.gene.c10.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c6,file = "DEGs/de.gene.c6.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c1,file = "DEGs/de.gene.c1.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c3,file = "DEGs/de.gene.c3.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c4,file = "DEGs/de.gene.c4.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c7,file = "DEGs/de.gene.c7.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)
write.table(de.gene.c2,file = "DEGs/de.gene.c2.txt",col.names = FALSE, row.names = FALSE,sep='\t',quote=FALSE)

##are de.gene.c8 genes related to cell cycle?
table(de.gene.c8 %in% cc.genes$s.genes)
table(de.gene.c8 %in% cc.genes$g2m.genes)

de.gene.c8[de.gene.c8 %in% cc.genes$s.genes] #'BRIP1''DTL''ATAD2'
de.gene.c8[de.gene.c8 %in% cc.genes$g2m.genes] #'SMC4''TOP2A''CENPE'

###



write.table(top50,file = "DEGs/top50.de.gene.txt",col.names = TRUE, row.names = FALSE,sep='\t',quote=FALSE)







# #####use magic smoothed data/scale.data slot, no improvement####
# placenta_magic <- placenta

# data_magic <- readRDS('data_magic.rds')
# placenta_magic@assays$RNA@data <- as(data_magic,'dgCMatrix')
# scale.data_magic <- readRDS('scale.data_magic.rds') #dgCMatrix
# placenta_magic@assays$RNA@scale.data <-  scale.data_magic


# options(repr.plot.height = 10, repr.plot.width = 7.5)
# DoHeatmap(placenta_magic, 
#           features = top10$gene,
#           slot= 'data', #'scale.data',
#           group.colors = unlist(map_cellcolor)
#           #disp.min = 0,disp.max=2,
#          ) +
#   scale_fill_gradientn(colors =  color_gradient_my)#+
#   #scale_fill_gradientn(colors =  c('blue','white','red'))#+
#   #scale_fill_gradientn(colors =  viridis(256, option = "D") )#+
#   #scale_fill_gradientn(colors =  mapal)#+
#   #scale_fill_gradientn(colors =  color_pyreds)#+ #not this
#   #NoLegend()  



# #####zoom in for 5 6 7 8 and plot again
# top10_c56789 <- top10[top10$cluster %in% c('5','6','7','8','9'),]
# #saveRDS(object=top10_c56789,file='DEGs/top10_c56789.rds',)
# top10_c56789 <- readRDS(file='DEGs/top10_c56789.rds',)

# #options(repr.plot.height = 7.5, repr.plot.width = 7.5)
# #DoHeatmap(placenta, features = top10_c56789$gene) + NoLegend()

# cluster <- Idents(placenta)# placenta@meta.data
# cluster.filter <- cluster[cluster %in%  c('5','6','7','8','9')]

# #pdf(file = 'DEGs/de.marker.c56789.heatmap.pdf',width = 6,height = 5,useDingbats = FALSE)
# pdf(file = 'DEGs/de.marker.c56789.heatmap.withlegend.pdf',width = 6,height = 5,useDingbats = FALSE)
# options(repr.plot.height = 5, repr.plot.width = 6)
# DoHeatmap(placenta, 
#           features = top10_c56789$gene,
#           cells = names(cluster.filter),
#           slot='scale.data',#'data',
#           group.colors = unlist(map_cellcolor)#,
#           #disp.min = -3,disp.max=3.5,
#          ) +
#   scale_fill_gradientn(colors =  color_gradient_my)#+
#   #scale_fill_gradientn(colors =  c('blue','white','red'))#+
#   #scale_fill_gradientn(colors =  viridis(256, option = "D") )#+
#   #scale_fill_gradientn(colors =  mapal)#+
#   #scale_fill_gradientn(colors =  color_pyreds)#+ #not this
#   #NoLegend()  
# dev.off()




###de.markers list from DEG heatmap##

#all.equal(subset(top10_c56789,cluster == '5'),subset(top10,cluster == '5') ) #TRUE

de.markers.c6 <- subset(top10_c56789,cluster == '6')$gene
de.markers.c8 <- subset(top10_c56789,cluster == '8')$gene
de.markers.c2 <- subset(top10,cluster == '2')$gene
de.markers.c7 <- subset(top10,cluster == '7')$gene
de.markers.c4 <- subset(top10,cluster == '4')$gene
de.markers.c5 <- subset(top10,cluster == '5')$gene

saveRDS(de.markers.c6,'DEGs/de.markers.c6.rds')
saveRDS(de.markers.c8,'DEGs/de.markers.c8.rds')
saveRDS(de.markers.c2,'DEGs/de.markers.c2.rds')
saveRDS(de.markers.c7,'DEGs/de.markers.c7.rds')
saveRDS(de.markers.c4,'DEGs/de.markers.c4.rds')
saveRDS(de.markers.c5,'DEGs/de.markers.c5.rds')


#######plot DEG with umap-dotheatmap(FeaturePlot, use data slot by default), violinPlot dotmatrix######

###quick look with FeaturePlot for all de.markers.cx
options(repr.plot.width=15,repr.plot.height=5)
FeaturePlot(placenta, features = de.markers.c5,ncol = 5,slot='scale.data',reduction = 'umapmirror',cols = color_gradient_my ) 
FeaturePlot(placenta_magic, features = de.markers.c5,ncol = 5,slot='scale.data',reduction = 'umapmirror',cols = color_gradient_my ) 

##quick look with VlnPlot (use data slot by default) for all de.markers.sel
options(repr.plot.width=15,repr.plot.height=7.5)
VlnPlot(placenta, features = de.markers.c5,pt.size = 0.1,slot = 'scale.data',cols = map_cellcolor)




####visualization and arrange single gene plot for identified top10 DEG


#####the plot function###
featurePlot_single <- function(gene = NULL,
                               min.cutoff = "q10", 
                               max.cutoff = "q90",
                               cols = c('lightgrey','darkred')
                              ){

    options(repr.plot.width=5,repr.plot.height=5)
    res.p <- FeaturePlot(placenta, 
                features = gene,
                ncol = 1,
                reduction = 'umap',
                #reduction = 'umapmirror',
                slot = 'scale.data',
                min.cutoff = min.cutoff , 
                max.cutoff = max.cutoff ,
                cols = cols,
                #cols =  viridis(256, option = "D"),
                #cols = c('lightgrey','darkred'),
                #cols = c('blue','white','red')
                #cols = rev(mapal),
                #cols = color_gradient_my
                pt.size = .1,
                label = FALSE #TRUE
               ) +

        theme(
            #legend.position = 'none',
            legend.position = 'top',
            axis.text=element_blank(), 
            axis.title = element_text(size = 15, face = "bold"),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(color="black", fill = NA,size=.3),
            plot.title = element_text(size = 15, face = "bold"),
            #complete = TRUE
            plot.margin = unit(c(1,1,1,1), "lines")
           )+
      ggtitle(gene) +
      xlim(left,right) + ylim(bottom,top) +
      labs(x = "UMAP1", y = "UMAP2")
    
      #print(res.p)
      return(res.p)

}
###

res.p <- list()
####cluster 2
##use for grey-red
##for(i in de.markers.c2){cat("'",i,"'","=","'","q1","',", sep = "") }
##for(i in de.markers.c2){cat("'",i,"'","=","'","q99","',", sep = "") }
low.q <- c('LINC02291'='q60','ADAMTS6'='q60','KLRD1'='q60','ADGRG6'='q60','LINC01483'='q60','EXT1'='q60','LINC00882'='q60','DTNB'='q60','MAP3K13'='q60','AC004784.1'='q10')
high.q <- c('LINC02291'='q90','ADAMTS6'='q90','KLRD1'='q90','ADGRG6'='q90','LINC01483'='q90','EXT1'='q90','LINC00882'='q90','DTNB'='q90','MAP3K13'='q90','AC004784.1'='q90')



for(gene in de.markers.c2){
  cat("gene:",gene,"cutoff:",low.q[[gene]],'-',high.q[[gene]])
  res.p[[gene]] <- featurePlot_single(gene = gene,min.cutoff = low.q[[gene]], max.cutoff = high.q[[gene]])
}

##pick 4 x 3
#for(i in de.markers.c2){cat("res.p[['",i,"']] + ", sep = "") }
options(repr.plot.height=12,repr.plot.width=15)
res.p[['LINC02291']] + res.p[['ADAMTS6']] + res.p[['KLRD1']] + res.p[['ADGRG6']] + res.p[['LINC01483']] + res.p[['EXT1']] + res.p[['LINC00882']] + res.p[['DTNB']] + res.p[['MAP3K13']] + res.p[['AC004784.1']] +
plot_layout(ncol=4,nrow=3)

ggsave(filename = 'DEGs/de.gene.picked.umap.grey-red.cluster2.pdf',height=12,width=15,useDingbats=FALSE)
####


##cluster 8
##for(i in de.markers.c8){cat("'",i,"'","=","'","q1","',", sep = "") }
##for(i in de.markers.c8){cat("'",i,"'","=","'","q99","',", sep = "") }

low.q <- c('DDX60'='q60','PPP1R9A'='q50','NMNAT2'='q60','CDK14'='q1','NEAT1'='q60','LVRN'='q60','SASH1'='q60','LIFR'='q60','NEDD4L'='q10','TAF4B'='q10')
high.q <- c('DDX60'='q99','PPP1R9A'='q99','NMNAT2'='q99','CDK14'='q99','NEAT1'='q99','LVRN'='q99','SASH1'='q99','LIFR'='q99','NEDD4L'='q99','TAF4B'='q99')


for(gene in de.markers.c8){
  cat("gene:",gene,"cutoff:",low.q[[gene]],'-',high.q[[gene]])
  res.p[[gene]] <- featurePlot_single(gene = gene,min.cutoff = low.q[[gene]], max.cutoff = high.q[[gene]])
}

##pick 4 x 3
#for(i in de.markers.c8){cat("res.p[['",i,"']] + ", sep = "") }
options(repr.plot.height=12,repr.plot.width=15)
res.p[['DDX60']] + res.p[['PPP1R9A']] + res.p[['NMNAT2']] + res.p[['CDK14']] + res.p[['NEAT1']] + res.p[['LVRN']] + res.p[['SASH1']] + res.p[['LIFR']] + res.p[['NEDD4L']] + res.p[['TAF4B']] +
plot_layout(ncol=4,nrow=3)

ggsave(filename = 'DEGs/de.gene.picked.umap.grey-red.cluster8.pdf',height=12,width=15,useDingbats=FALSE)
####

##cluster 6
##for(i in de.markers.c6){cat("'",i,"'","=","'","q1","',", sep = "") }
##for(i in de.markers.c6){cat("'",i,"'","=","'","q99","',", sep = "") }

low.q <- c('PTCHD4'='q1','AC097478.1'='q50','INPP5D'='q1','MIR34AHG'='q1','AC092167.1'='q1','PTPRM'='q50','KIAA1217'='q50','GAS7'='q60','BBS9'='q50','MDM2'='q60')
high.q <- c('PTCHD4'='q99','AC097478.1'='q99','INPP5D'='q99','MIR34AHG'='q99','AC092167.1'='q99','PTPRM'='q99','KIAA1217'='q99','GAS7'='q99','BBS9'='q99','MDM2'='q99')


for(gene in de.markers.c6){
  cat("gene:",gene,"cutoff:",low.q[[gene]],'-',high.q[[gene]])
  res.p[[gene]] <- featurePlot_single(gene = gene,min.cutoff = low.q[[gene]], max.cutoff = high.q[[gene]])
}

##pick 4 x 3
#for(i in de.markers.c6){cat("res.p[['",i,"']] + ", sep = "") }
options(repr.plot.height=12,repr.plot.width=15)
res.p[['PTCHD4']] + res.p[['AC097478.1']] + res.p[['INPP5D']] + res.p[['MIR34AHG']] + res.p[['AC092167.1']] + res.p[['PTPRM']] + res.p[['KIAA1217']] + res.p[['GAS7']] + res.p[['BBS9']] + res.p[['MDM2']] +
plot_layout(ncol=4,nrow=3)

ggsave(filename = 'DEGs/de.gene.picked.umap.grey-red.cluster6.pdf',height=12,width=15,useDingbats=FALSE)
####



##cluster 7
##for(i in de.markers.c7){cat("'",i,"'","=","'","q1","',", sep = "") }
##for(i in de.markers.c7){cat("'",i,"'","=","'","q99","',", sep = "") }

low.q <- c('AC018754.1'='q1','CD96'='q10','GRB14'='q50','GPC5'='q50','AC092920.1'='q50','LINC01588'='q50','LINC01505'='q50','PDE4D'='q60','AC004704.1'='q10','BACE2'='q60')
high.q <- c('AC018754.1'='q99','CD96'='q99','GRB14'='q99','GPC5'='q99','AC092920.1'='q99','LINC01588'='q99','LINC01505'='q99','PDE4D'='q99','AC004704.1'='q99','BACE2'='q99')


for(gene in de.markers.c7){
  cat("gene:",gene,"cutoff:",low.q[[gene]],'-',high.q[[gene]])
  res.p[[gene]] <- featurePlot_single(gene = gene,min.cutoff = low.q[[gene]], max.cutoff = high.q[[gene]])
}

##pick 4 x 3
#for(i in de.markers.c7){cat("res.p[['",i,"']] + ", sep = "") }
options(repr.plot.height=12,repr.plot.width=15)
res.p[['AC018754.1']] + res.p[['CD96']] + res.p[['GRB14']] + res.p[['GPC5']] + res.p[['AC092920.1']] + res.p[['LINC01588']] + res.p[['LINC01505']] + res.p[['PDE4D']] + res.p[['AC004704.1']] + res.p[['BACE2']]  +
plot_layout(ncol=4,nrow=3)

ggsave(filename = 'DEGs/de.gene.picked.umap.grey-red.cluster7.pdf',height=12,width=15,useDingbats=FALSE)
####



##cluster 5
##for(i in de.markers.c5){cat("'",i,"'","=","'","q1","',", sep = "") }
##for(i in de.markers.c5){cat("'",i,"'","=","'","q99","',", sep = "") }

low.q <- c('EPHA1-AS1'='q50','LAMA3'='q50','GPR78'='q50','STAT4'='q50','AJ009632.2'='q50','SLC45A4'='q50','ADAMTSL1'='q50','SLC26A7'='q50','AL691420.1'='q50','ZNF117'='q50')
high.q <- c('EPHA1-AS1'='q99','LAMA3'='q99','GPR78'='q99','STAT4'='q99','AJ009632.2'='q99','SLC45A4'='q99','ADAMTSL1'='q99','SLC26A7'='q99','AL691420.1'='q99','ZNF117'='q99')

for(gene in de.markers.c5){
  cat("gene:",gene,"cutoff:",low.q[[gene]],'-',high.q[[gene]])
  res.p[[gene]] <- featurePlot_single(gene = gene,min.cutoff = low.q[[gene]], max.cutoff = high.q[[gene]])
}

##pick 4 x 3
#for(i in de.markers.c5){cat("res.p[['",i,"']] + ", sep = "") }
options(repr.plot.height=12,repr.plot.width=15)
res.p[['EPHA1-AS1']] + res.p[['LAMA3']] + res.p[['GPR78']] + res.p[['STAT4']] + res.p[['AJ009632.2']] + res.p[['SLC45A4']] + res.p[['ADAMTSL1']] + res.p[['SLC26A7']] + res.p[['AL691420.1']] + res.p[['ZNF117']]  +
plot_layout(ncol=4,nrow=3)

ggsave(filename = 'DEGs/de.gene.picked.umap.grey-red.cluster5.pdf',height=12,width=15,useDingbats=FALSE)
####



##save all the de genes plot 
saveRDS(res.p,'res.p.rds')


##arrange select genes for cluster 7,6,8,2 output 
de.markers.sel <- c('BACE2','PDE4D','PTCHD4','INPP5D','CDK14',"DDX60",'ADGRG6','DTNB' )

##pick 4 x 2
#for(i in de.markers.sel){cat("res.p[['",i,"']] + ", sep = "") }
options(repr.plot.height=7.5,repr.plot.width=15)
res.p[['BACE2']] + res.p[['PDE4D']] + res.p[['PTCHD4']] + res.p[['INPP5D']] + res.p[['CDK14']] + res.p[['DDX60']] + res.p[['ADGRG6']] + res.p[['DTNB']]  +
plot_layout(ncol=4,nrow=2)

ggsave(filename = 'DEGs/de.gene.picked.umap.grey-red.markers.sel.pdf',height=7.5,width=15,useDingbats=FALSE)


#save by file
for (i in de.markers.sel){
    cat ('gene:',i,'\n')
    options(repr.plot.height=5,repr.plot.width=5) 
    print(res.p[[i]])
    ggsave(filename = paste('DEGs/de_genes_sel.split/de.gene.picked.umap.grey-red.markers.',i,  '.pdf',sep=""),height=5,width=5,useDingbats=FALSE)
    
}



####

###########new STB marker?#######
#Karvas, R. M., McInturf, S., Zhou, J., Ezashi, T., Schust, D. J., Roberts, R. M., et al. (2020). Use of a human embryonic stem cell model to discover GABRP, WFDC2, VTCN1 and ACTC1 as markers of early first trimester human trophoblast. Mol. Hum. Reprod

genes <- c('ACTC1', 'GABRP', 'VTCN1','WFDC2')
res.p <- list()
for(gene in genes){
  cat("gene:",gene)#"cutoff:",low.q[[gene]],'-',high.q[[gene]])
  res.p[[gene]] <- featurePlot_single(gene = gene,min.cutoff = 'q10', max.cutoff = 'q80')
}

options(repr.plot.height=10,repr.plot.width=10)
res.p[['ACTC1']] + res.p[['GABRP']] + res.p[['VTCN1']] + res.p[['WFDC2']] +
plot_layout(ncol=2,nrow=2)



######decide cluster 2 type
options(repr.plot.height=8,repr.plot.width=6)
DimPlot(placenta,reduction = 'umap',label = TRUE, label.size = 10)

options(repr.plot.height=5,repr.plot.width=12)
VlnPlot(placenta,features = c('PAPPA','FLT1'),group.by = 'cluster_order' )



########################







