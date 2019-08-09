#MCLclusterrawplot <- function(whichdata, whichprotein, graphtitle): This plots the three points of the three triplicates (not the mean), against the three temperatures with -/+ Fe in two different colors. 
#It also adds a * next to the number that is is statistically significant at the 0.05 level. 
#'protein' is just a place holder,it can be any word

MCLclusterrawplot <- function(whichdata, whichprotein, graphtitle){
   protein <- dplyr::filter(whichdata, grepl (whichprotein, cluster))
  
  T0 <- protein[1, 47]  
  T0 <- round (T0, digits = 3)
  
  T3 <- protein[1, 49]
  T3 <- round (T3, digits = 3)
  
  T6 <- protein[1, 51]
  T6 <- round (T6, digits = 3)
  if (protein [1,46]< 0.05){
    T0 <- paste0(T0, "*")
   } 
  if (protein [1,48]< 0.05){
      T3 <- paste0(T3, "*") 
    }  
  if(protein [1,50]< 0.05){
      T6 <- paste0(T6, "*") } 
    
  protein <- protein [c(3:21)]
  #protein <- protein [-c(2:28, 47:76)]
  names(protein)<- c("name", "A_NoFe_0.5", "B_NoFe_0.5", "C_NoFe_0.5", "A_Fe_0.5", "B_Fe_0.5", "C_Fe_0.5", "A_NoFe_3", "B_NoFe_3", "C_NoFe_3", "A_Fe_3", "B_Fe_3", "C_Fe_3","A_NoFe_6", "B_NoFe_6", "C_NoFe_6", "A_Fe_6", "B_Fe_6", "C_Fe_6")
  protein <-  melt (data = protein,id.vars = "name")
  protein <- separate(protein, col = variable, into = c("Rep", "Fe", "temperature"))
  ggplot(protein, aes(x=temperature,
                      y=value, 
                      colour = Fe ))+
    geom_point(size = 2.5)+
    scale_color_manual(values = c("red", "black"))+
    xlab ("Temperature (°C)")+
    ylab ("Expression-Value")+
    ggtitle(graphtitle)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text( vjust = 1.0, hjust = 1.0))+
    scale_color_manual (name = "", labels = c("+Fe", "-Fe"), values = c("red", "black"))+
   annotate ("text", x = 1:3, y = 0, label= c(T0, T3, T6), color="darkgreen", size = 4, fontface="bold")+
    theme(axis.text.x = element_text (color= "black", #face = "bold", 
                                      size = 14))+
    theme(axis.text.y = element_text (color= "black", #face = "bold", 
                                      size = 14))+
    theme(axis.title.y = element_text (color= "black", face = "bold", size = 16))+
    theme(axis.title.x = element_text (color= "black", face = "bold", size = 16))
}

###`````````````````````````````````````````````````````````````````````````````````````````###

#heatmap function: First it takes a subset from the dataset with whatever name 
#of protein that we want. Then it re-shapes the data by making the protein name
#a row header then it reshapes it into a matrix, then it plots it on a clustered
#heatmap. 

#Arguements:
#group = group of proteins we want to work with
#dataset = which dataset we are using (typically one that has been cleaned)
#proteinname = whatever proteins we want to use

heatmap <- function (dataset, proteinname){
  group <- dplyr::filter (dataset, grepl(proteinname, name))
  group <- group [c(3:21)]
  names(group)<- c("name", "A_noFe_0", "B_noFe_0", "C_noFe_0", "A_Fe_0", "B_Fe_0", "C_Fe_0", "A_noFe_3", "B_noFe_3", "C_noFe_3", "A_Fe_3", "B_Fe_3", "C_Fe_3","A_noFe_6", "B_noFe_6", "C_noFe_6", "A_Fe_6", "B_Fe_6", "C_Fe_6")
  group <- data.frame(group[,-1], row.names = group [,1])
  group <- as.matrix ((group))
  heatmap.2 (group,
             trace = "none", 
             density = "none", 
             col = bluered(100), 
             cexRow = 0.5, 
             cexCol =1,
             margins = c(8,20),
             scale = "row",
             keysize = 1,
             Colv=TRUE)
  

}

###`````````````````````````````````````````````````````````````###

#heatmap_unclustered function:First it takes a subset from the dataset with whatever name 
#of protein that we want. Then it re-shapes the data by making the protein name
#a row header then it reshapes it into a matrix, then it plots it on an unclustered
#heatmap. 

#Arguments:
#group = group of proteins we want to work with
#dataset = which dataset we are using (typically one that has been cleaned)
#proteinname = whatever proteins we want to use

heatmap_unclustered <- function (dataset, proteinname){
  group <- dplyr::filter (dataset, grepl(proteinname, name))
  group <- group [c(3:21)]
  names(group)<- c("name", "A_NoFe_0.5", "B_NoFe_0.5", "C_NoFe_0.5", "A_Fe_0.5", "B_Fe_0.5", "C_Fe_0.5", "A_NoFe_3", "B_NoFe_3", "C_NoFe_3", "A_Fe_3", "B_Fe_3", "C_Fe_3","A_NoFe_6", "B_NoFe_6", "C_NoFe_6", "A_Fe_6", "B_Fe_6", "C_Fe_6")
  group <- data.frame(group[,-1], row.names = group [,1])
  group <- as.matrix ((group))
  heatmap.2 (group, 
             trace = "none",
             density = "none", 
             col = bluered(100),
             cexRow = 0.5, 
             cexCol =1,
             margins = c(8,20),
             scale = "row" ,
             dendrogram = 'none',  
             Colv = FALSE, 
             Rowv = FALSE,
             keysize = 1, 
             colsep=c(3,6,9,12,15), 
             sepcolor = "black", 
             sepwidth = (0.01))
}


# 
# heatmap_unclustered <- function (dataset, proteinname){
#   group <- dplyr::filter (dataset, grepl(proteinname, name))
#   group <- group [c(3:21)]
#   names(group)<- c("name", "A_NoFe_0.5", "B_NoFe_0.5", "C_NoFe_0.5", "A_Fe_0.5", "B_Fe_0.5", "C_Fe_0.5", "A_NoFe_3", "B_NoFe_3", "C_NoFe_3", "A_Fe_3", "B_Fe_3", "C_Fe_3","A_NoFe_6", "B_NoFe_6", "C_NoFe_6", "A_Fe_6", "B_Fe_6", "C_Fe_6")
#   group <- data.frame(group[,-1], row.names = group [,1])
#   group <- as.matrix ((group))
#   heatmap.2 (group, 
#              trace = "none",
#              density = "none", 
#              col = bluered(100),
#              cexRow = 0.5, 
#              cexCol =1,
#              margins = c(8,20),
#              scale = "row" ,
#              dendrogram = 'none',  
#              Colv = FALSE, 
#              Rowv = TRUE,
#              keysize = 1)
# }
###``````````````````````````````````````````````````````````````###



heatmap_mean <- function  (dataset, proteinname){
  group <- dplyr::filter (dataset, grepl(proteinname, name))
  group <- group [c(3:21)]
  colname <- "mean"
  i=0
  while (i < 6){
    colname <- paste(colname, i)
    col1 = 2+(i*3)
    col2 = 4+(i*3)
    group[[colname]] <- rowMeans( group[, c(col1:col2)] ) 
    i <- i+1
  }
  
  group_sub <- group [c(1, 20:25)]
  names(group_sub)<- c("name", "noFe_0C", "Fe_0C", "noFe_3C", "Fe_3C","noFe_6C", "Fe_6C")
  group_sub <- data.frame(group_sub[,-1], row.names = group_sub [,1])
  group_sub <- as.matrix ((group_sub))
  heatmap.2 (group_sub, 
             trace = "none", 
             density = "none",
             col = bluered(100),
             cexRow =0.5,             
             cexCol =1, 
             margins = c(6,20),
             scale = "row", 
             #srtRow = 10, 
             keysize = 1
             
             
  )
  
}







heatmap_unclustered_mean <- function  (dataset, proteinname){
  group <- dplyr::filter (dataset, grepl(proteinname, name))
  group <- group [c(3:21)]
  colname <- "mean"
  i=0
  while (i < 6){
    colname <- paste(colname, i)
    col1 = 2+(i*3)
    col2 = 4+(i*3)
    group[[colname]] <- rowMeans( group[, c(col1:col2)] ) 
    i <- i+1
  }
  
  group_sub <- group [c(1, 20:25)]
  names(group_sub)<- c("name", "noFe_0C", "Fe_0C", "noFe_3C", "Fe_3C","noFe_6C", "Fe_6C")
  group_sub <- data.frame(group_sub[,-1], row.names = group_sub [,1])
  group_sub <- as.matrix ((group_sub))
  heatmap.2 (group_sub, 
             trace = "none", 
             density = "none",
             col = bluered(100),
             cexRow =0.5,             
             cexCol =1, 
             margins = c(6,20),
             scale = "row", 
             dendrogram = 'none',  
             keysize = 1, 
             Colv=FALSE,
             Rowv = FALSE
  )
  
}



heatmap_unclustered <- function (dataset, proteinname){
  group <- dplyr::filter (dataset, grepl(proteinname, name))
  group <- group [c(3:21)]
  names(group)<- c("name", "A_NoFe_0.5", "B_NoFe_0.5", "C_NoFe_0.5", "A_Fe_0.5", "B_Fe_0.5", "C_Fe_0.5", "A_NoFe_3", "B_NoFe_3", "C_NoFe_3", "A_Fe_3", "B_Fe_3", "C_Fe_3","A_NoFe_6", "B_NoFe_6", "C_NoFe_6", "A_Fe_6", "B_Fe_6", "C_Fe_6")
  group <- data.frame(group[,-1], row.names = group [,1])
  group <- as.matrix ((group))
  heatmap.2 (group, 
             trace = "none",
             density = "none", 
             col = bluered(100),
             cexRow = 0.5, 
             cexCol =1,
             margins = c(8,20),
             scale = "row" ,
             dendrogram = 'none',  
             Colv = FALSE, 
             Rowv = FALSE,
             keysize = 1,
             colsep=c(3,6,9,12,15),
             sepcolor = "black",
             sepwidth = c(0.02, 0.02))
}







###```````````````````````````````````````````````````````````````###


volcanoplot <- function(dataset, foldchange, FDR, graphtitle){
  volcano <- dataset[c("name",foldchange, FDR,"Taxa.group" )]
  x <- foldchange
  
  names(volcano) <- c( "gene", "Fold", "FDR", "taxa")
  volcano ["group"] <- "no-fold-change"
  
  if(grepl("Fe_vs_noFe_foldchange", x)) {
    volcano['size'] <- rowMeans( dataset[, c(4:21)])
  } else if (grepl("temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("temp_6C_vs_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:21)])
  } else if (grepl("int_fe_temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("int_fe_temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("noFe_3C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6,10:12 )])
  } else if (grepl("noFe_6C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6, 16:18)])
  } else if (grepl("Fe_3C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(7:9, 13:15)])
  } else if (grepl("Fe_6C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans(dataset[, c(7:9, 19:21)])
  } else if (grepl("noFe_6C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:18,10:12)])
  } else if (grepl("Fe_6C_vs_Fe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(13:15, 19:21)])
  } else if (grepl("Fe_0C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9)])
  } else if (grepl("Fe_3C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:15)])
  } else if (grepl("Fe_6C_vs_noFe_6C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:21)])
  }
  

  volcano [which(volcano['FDR'] < 0.05 & abs(volcano['Fold']) < 2), "group"]<- "nofoldchange_significant" 
  volcano [which(volcano ['FDR'] > 0.05 & abs(volcano['Fold'])> 2), "group"] <- "fold-change_notsignificant"
  volcano [which(volcano ['FDR'] < 0.05 & abs(volcano['Fold']) > 2), "group"] <- "fold-change_significant"
  
  ggplot (data = volcano, aes(x= Fold, y=-log10(FDR), color = taxa, label=gene)) +
    geom_point (alpha=0.5)+
    scale_shape (solid=FALSE)+
    geom_hline (yintercept = 1.301, col = "blue", lty = 2, lwd=0.5 )+
    geom_vline (xintercept = c(-1,1), col = "blue", lty = 2, lwd=0.5 )+
    xlab("Log2 Fold Change")+ ylab("-Log10_P-value")+
    scale_color_manual(values = c("black", "red", "#00FF00"))+
    ggtitle(graphtitle)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
    #geom_label_repel(data = subset(volcano, abs(Fold) >1 & -log10(FDR) >1.3),
     # aes(label=gene), box.padding = 0.2, point.padding = 0.5, color= 'black', size =1)
}

#size = (volcano$size*100))
##-------------------------------------------------------------------------##

volcanoplot_annotated <- function(dataset, foldchange, FDR, graphtitle){
  volcano <- dataset[c("name",foldchange, FDR,"Taxa.group" )]
  x <- foldchange
  
  names(volcano) <- c( "gene", "Fold", "FDR", "taxa")
  volcano ["group"] <- "no-fold-change"
  
  if(grepl("Fe_vs_noFe_foldchange", x)) {
    volcano['size'] <- rowMeans( dataset[, c(4:21)])
  } else if (grepl("temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("temp_6C_vs_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:21)])
  } else if (grepl("int_fe_temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("int_fe_temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("noFe_3C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6,10:12 )])
  } else if (grepl("noFe_6C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6, 16:18)])
  } else if (grepl("Fe_3C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(7:9, 13:15)])
  } else if (grepl("Fe_6C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans(dataset[, c(7:9, 19:21)])
  } else if (grepl("noFe_6C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:18,10:12)])
  } else if (grepl("Fe_6C_vs_Fe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(13:15, 19:21)])
  } else if (grepl("Fe_0C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9)])
  } else if (grepl("Fe_3C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:15)])
  } else if (grepl("Fe_6C_vs_noFe_6C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:21)])
  }
  
  
  volcano [which(volcano['FDR'] < 0.05 & abs(volcano['Fold']) < 2), "group"]<- "nofoldchange_significant" 
  volcano [which(volcano ['FDR'] > 0.05 & abs(volcano['Fold'])> 2), "group"] <- "fold-change_notsignificant"
  volcano [which(volcano ['FDR'] < 0.05 & abs(volcano['Fold']) > 2), "group"] <- "fold-change_significant"
  
  ggplot (data = volcano, aes(x= Fold, y=-log10(FDR), color = taxa, label=gene)) +
    geom_point (alpha=0.5)+ #size = (volcano$size*100))+
    scale_shape (solid=FALSE)+
    geom_hline (yintercept = 1.301, col = "blue", lty = 2, lwd=0.5 )+
    geom_vline (xintercept = c(-1,1), col = "blue", lty = 2, lwd=0.5 )+
    xlab("Log2 Fold Change")+ ylab("-Log10_P-value")+
    scale_color_manual(values = c("black", "red", "#00FF00"))+
    ggtitle(graphtitle)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
  geom_label_repel(data = subset(volcano, abs(Fold) >1 & -log10(FDR) >1.3),
   aes(label=gene), box.padding = 0.2, point.padding = 0.5, color= 'black', size =1)
}


##--------------------------------------------------------------------------------------------------------

volcanoplot_interactive <- function(dataset, Fold, FDR, volcanontitle){
  volcano <- dataset[c("name",Fold, FDR,"Taxa.group")]
  x <- Fold
  
  names(volcano) <- c( "gene", "Fold", "FDR", "taxa")
  volcano ["group"] <- "no-fold-change"
  
  if(grepl("Fe_vs_noFe_foldchange", x)) {
    volcano['size'] <- rowMeans( dataset[, c(4:21)])
  } else if (grepl("temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("temp_6C_vs_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:21)])
  } else if (grepl("int_fe_temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("int_fe_temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("noFe_3C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6,10:12 )])
  } else if (grepl("noFe_6C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6, 16:18)])
  } else if (grepl("Fe_3C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(7:9, 13:15)])
  } else if (grepl("Fe_6C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans(dataset[, c(7:9, 19:21)])
  } else if (grepl("noFe_6C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:18,10:12)])
  } else if (grepl("Fe_6C_vs_Fe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(13:15, 19:21)])
  } else if (grepl("Fe_0C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9)])
  } else if (grepl("Fe_3C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:15)])
  } else if (grepl("Fe_6C_vs_noFe_6C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:21)])
  }
  
  volcano [which(volcano['FDR'] < 0.05 & abs(volcano['Fold']) < 2), "group"]<- "nofoldchange_significant" 
  
  volcano [which(volcano ['FDR'] > 0.05 & abs(volcano['Fold'])> 2), "group"] <- "fold-change_notsignificant"
  
  volcano [which(volcano ['FDR'] < 0.05 & abs(volcano['Fold']) > 2), "group"] <- "fold-change_significant"
  
  vline <- function(x = 0, color = "black") {
    list(
      type = "line", 
      y0 = 0, 
      y1 = 1, 
      yref = "paper",
      x0 = x, 
      x1 = x, 
      line = list(color = color, dash='dot', width=1)
    )
  }
  
  hline <- function(y = 0, color = "black") {
    list(
      type = "line", 
      x0 = 0, 
      x1 = 1, 
      xref = "paper",
      y0 = y, 
      y1 = y, 
      line = list(color = color, dash='dot', width=1)
    )
  }
  
  x <- list (title ="Log2 Fold Change")
  y <- list (title = "-Log10_P-value")
 # ax<- list(showline = TRUE, zeroline=FALSE, linecolor = toRGB("black"),
            #linewidth = 6)
  plot_ly(data = volcano, x = volcano$Fold, y = -log10(volcano$FDR), text = volcano$gene, mode = "markers", symbol = volcano$taxa)%>% #size = volcano$size )%>% 
    layout(title = volcanontitle)%>%
    layout (xaxis = x, yaxis = y)%>%
    layout(shapes = list(vline (-1), vline(1), hline(1.301)))
    #layout (xaxis = ax, yaxis =ax)
}

##..........................................................









##-------------------------------------------------------------------##
##This function plots the log2 fold change against the name of the cluster. 
#thi allows for easy comparison of clusters between groups. 
#dataset:the data you want to enter, foldchange: you select which treatment foldchange you want, same with FDR, and choose a graph title
clustername_vs_foldchange <- function(dataset, foldchange, FDR, graphtitle){
  volcano <- dataset[c("name",foldchange, FDR,"Taxa.group", "cluster" )]
  x <- foldchange
  
  names(volcano) <- c( "gene", "Fold", "FDR", "taxa", "cluster")
  volcano ["group"] <- "no-fold-change"
  
  if(grepl("Fe_vs_noFe_foldchange", x)) {
    volcano['size'] <- rowMeans( dataset[, c(4:21)])
  } else if (grepl("temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("temp_6C_vs_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:21)])
  } else if (grepl("int_fe_temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("int_fe_temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("noFe_3C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6,10:12 )])
  } else if (grepl("noFe_6C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6, 16:18)])
  } else if (grepl("Fe_3C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(7:9, 13:15)])
  } else if (grepl("Fe_6C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans(dataset[, c(7:9, 19:21)])
  } else if (grepl("noFe_6C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:18,10:12)])
  } else if (grepl("Fe_6C_vs_Fe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(13:15, 19:21)])
  } else if (grepl("Fe_0C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9)])
  } else if (grepl("Fe_3C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:15)])
  } else if (grepl("Fe_6C_vs_noFe_6C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:21)])
  }
  
  
  ggplot (data = volcano, aes(x= Fold, y= cluster, color = taxa, label=gene)) +
    geom_point (alpha=0.5 ,size = (volcano$size*50))+
    scale_shape (solid=FALSE)+
    ylab("cluster")+ 
    xlab("log2 fold-change")+
    scale_color_manual(values = c("black", "red", "blue"))+
    ggtitle(graphtitle)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(
            hjust=1,
      size = 9, 
      face = "bold")) +
  geom_text_repel(data = subset(volcano, abs(Fold) >1 & -log10(FDR) >1.3),
  aes(label="*"), box.padding = 0, point.padding = 0.1, color= 'black', size =5)
}

##````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````##
#Does the same thing as culsternamevsfoldchange but for ORFs
ORFname_vs_foldchange <- function(dataset, foldchange, FDR, graphtitle){
  volcano <- dataset[c("name",foldchange, FDR,"Taxa.group", "orf_id" )]
  x <- foldchange
  
  names(volcano) <- c( "gene", "Fold", "FDR", "taxa", "orf_id")
  volcano ["group"] <- "no-fold-change"
  
  if(grepl("Fe_vs_noFe_foldchange", x)) {
    volcano['size'] <- rowMeans( dataset[, c(4:21)])
  } else if (grepl("temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("temp_6C_vs_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:21)])
  } else if (grepl("int_fe_temp_6C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9, 16:21)])
  } else if (grepl("int_fe_temp_3C_vs_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:15)])
  } else if (grepl("noFe_3C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6,10:12 )])
  } else if (grepl("noFe_6C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:6, 16:18)])
  } else if (grepl("Fe_3C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(7:9, 13:15)])
  } else if (grepl("Fe_6C_vs_Fe_0C_foldchange", x)){
    volcano['size'] <- rowMeans(dataset[, c(7:9, 19:21)])
  } else if (grepl("noFe_6C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:18,10:12)])
  } else if (grepl("Fe_6C_vs_Fe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(13:15, 19:21)])
  } else if (grepl("Fe_0C_vs_noFe_0C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(4:9)])
  } else if (grepl("Fe_3C_vs_noFe_3C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(10:15)])
  } else if (grepl("Fe_6C_vs_noFe_6C_foldchange", x)){
    volcano['size'] <- rowMeans( dataset[, c(16:21)])
  }
  
  
  ggplot (data = volcano, aes(x= Fold, y= orf_id, color = taxa, label=gene)) +
    geom_point (alpha=0.5)+# size = (volcano$size * 50 ))+
    scale_shape (solid=FALSE)+
    ylab("ORF_ID")+ 
    xlab("log2 fold-change")+
    scale_color_manual(values = c("black", "red", "blue"))+
    ggtitle(graphtitle)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(
      hjust=1,
      size = 9, 
      face = "bold")) +
    geom_text_repel(data = subset(volcano, abs(Fold) >1 & -log10(FDR) >1.3),
                    aes(label="*"), box.padding = 0, point.padding = 0.1, color= 'black', size =5)
}
###.........................................................................................###


















#This function creates means of the the three treatments and plots the log values 
# xyplotISIP <- function(dataset, title){
#   
#   datasetmelt <-  melt (data = dataset,id.vars = "name", measure.vars = c(2:28))
#   
#   colname <- "mean"
#   i=0
#   while (i < 7){
#     colname <- paste(colname, i)
#     col1 = 2+(i*3)
#     col2 = 4+(i*3)
#     dataset[[colname]] <- rowMeans( dataset[, c(col1:col2)] ) 
#     i <- i+1
#   }
#   
#   datasetmean <- dataset [-c(2:28)]
#   names(datasetmean) <- c("name", "T0","0.5C", "Fe_0.5C", "3C", "Fe_3C", "6C", "Fe_6C" )
#   
#   datasetmeanmelt <-  melt (data = datasetmean,id.vars = "name", measure.vars = c(2:10)) 
#   
#   
#   
#   print (ggplot(datasetmelt, aes(x=variable, y=value, colour = name ))+
#            geom_point(alpha=0.9)+
#            stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.25)+
#            stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
#            xlab ("Treatmet")+
#            ylab ("Expression-Value")+
#            ggtitle(title)+
#            theme_bw()+
#            theme(plot.title = element_text(hjust = 0.5))+
#            theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
#            theme(legend.position = "none"))
#   
#   
#   print(ggplot(datasetmeanmelt, aes(x=variable, y=value, colour = name ))+
#           geom_point(alpha=0.9)+
#           stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.25)+
#           stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
#           xlab ("Treatmet")+
#           ylab ("Expression-Value_Log10")+
#           ggtitle(title)+
#           theme_bw()+
#           theme(plot.title = element_text(hjust = 0.5))+
#           theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
#           scale_y_log10()+
#           theme(legend.position = "none"))
# }



###``````````````````````````````````````````````
#This function creates means of the the three treatments and plots the log values 
xyplotconcise <- function(dataset, title){
  
  datasetmelt <-  melt (data = dataset,id.vars = "name", measure.vars = c(2:28))
  
  colname <- "mean"
  i=0
  while (i < 9){
    colname <- paste(colname, i)
    col1 = 2+(i*3)
    col2 = 4+(i*3)
    dataset[[colname]] <- rowMeans( dataset[, c(col1:col2)] ) 
    i <- i+1
  }
  
  datasetmean <- dataset [-c(2:28)]
  names(datasetmean) <- c("name", "T0", "B12_0.5C", "Fe_B12_0.5C","0.5C", "Fe_0.5C", "3C", "Fe_3C", "6C", "Fe_6C" )
  
  datasetmeanmelt <-  melt (data = datasetmean,id.vars = "name", measure.vars = c(2:10)) 
  
  
  
  print (ggplot(datasetmelt, aes(x=variable, y=value, colour = name ))+
           geom_point(alpha=0.9)+
           stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.25)+
           stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
           xlab ("Treatmet")+
           ylab ("Expression-Value")+
           ggtitle(title)+
           theme_bw()+
           theme(plot.title = element_text(hjust = 0.5))+
           theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
           theme(legend.position = "none"))
  
  
  print(ggplot(datasetmeanmelt, aes(x=variable, y=value, colour = name ))+
          geom_point(alpha=0.9)+
          stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.25)+
          stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
          xlab ("Treatmet")+
          ylab ("Expression-Value_Log10")+
          ggtitle(title)+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5))+
          theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
          scale_y_log10()+
          theme(legend.position = "none"))
  }

###````````````````````````````````````````````````````````````````````###
  
xyplotsmalldata <- function(dataset, title){
  
  datasetmelt <-  melt (data = dataset,id.vars = "name", measure.vars = c(2:28))
  
  colname <- "mean"
  i=0
  while (i < 9){
    colname <- paste(colname, i)
    col1 = 2+(i*3)
    col2 = 4+(i*3)
    dataset[[colname]] <- rowMeans( dataset[, c(col1:col2)] ) 
    i <- i+1
  }
  
  datasetmean <- dataset [-c(2:28)]
  names(datasetmean) <- c("name", "T0", "B12_0.5C", "Fe_B12_0.5C","0.5C", "Fe_0.5C", "3C", "Fe_3C", "6C", "Fe_6C" )
  
  datasetmeanmelt <-  melt (data = datasetmean,id.vars = "name", measure.vars = c(2:10)) 
  
  
  
  print (ggplot(datasetmelt, aes(x=variable, y=value, colour = name ))+
           geom_point()+
           geom_line(aes(group=name))+
           xlab ("Treatmet")+
           ggtitle(title)+
           theme_bw()+
           theme(plot.title = element_text(hjust = 0.5))+
           theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
           theme(legend.position = "none"))

  print(ggplot(datasetmeanmelt, aes(x=variable, y=value, colour = name ))+
          geom_point()+
          geom_line(aes(group=name))+
          xlab ("Treatmet")+
          ylab ("Value(Log10)")+
          ggtitle(title)+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5))+
          theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
          scale_y_log10()+
          theme(legend.position = "none"))
  
  print(ggplot(datasetmeanmelt, aes(x=variable, y=value, colour = name ))+
          geom_point(alpha=0.9)+
          stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.25)+
          stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
          xlab ("Treatmet")+
          ylab ("Expression-Value_Log10")+
          ggtitle(title)+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5))+
          theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
          scale_y_log10()+
          theme(legend.position = "none"))
}
  
###``````````````````````````````````````````````````````````````````````````````###

#simple line and point plot of all the raw flavodoxin data under the different treatments. 
#ggplot(flavodoxinmelt, aes(x=variable, y=value ))+
  #geom_point()+
  #geom_line (aes(group=name))+
  #xlab ("Treatmet")+
  #ggtitle("flavodoxin")+
  #theme_bw()+
  #theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))


#creates a dotplot and plots the means of all the different flavodoxins into a red dot
#ggplot(flavodoxinmelt, aes(x=variable, y=value ))+
  #geom_point(alpha =0.3)+
  #stat_summary(fun.data = mean_cl_normal, geom = "point", colour = "red")+
  #xlab ("Treatmet")+
  #ggtitle("Flavodoxin")+
  #theme_bw()+
  #theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))

#colours the differet flavodoxins and creates a mean line with error bars. 
#ggplot(flavodoxinmelt, aes(x=variable, y=value, colour = name ))+
  #geom_point(alpha=0.9)+
  #stat_summary(fun.y = mean, geom = "line", colour = "red", group =1, size = 1)+
  #stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
  #xlab ("Treatmet")+
  #ggtitle("Flavodoxin")+
  #theme_bw()+
  #theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
  #theme(legend.position = "none")


#this zooms in on the proteins that are expressed a little bit. 
#ggplot(flavodoxinmelt, aes(x=variable, y=value, colour = name ))+
  #geom_line(aes(group=name,alpha=0.9))+
  #stat_summary(fun.y = mean, geom = "line", colour = "red", group =1, size = 1)+
  #stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
  #xlab ("Treatmet")+
  #ggtitle("Flavodoxin")+
  #theme_bw()+
  #theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
  #theme(legend.position = "none")+
  #scale_y_continuous(limits = c(0,0.001))

#log(10) of y-values because of large spread in values
#ggplot(flavodoxinmelt, aes(x=variable, y=value, colour = name ))+
  #geom_line(aes(group=name,alpha=0.9))+
  #stat_summary(fun.y = mean, geom = "line", colour = "red", group =1, size = 1)+
  #stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
  #xlab ("Treatmet")+
  #ylab ("Value(Log10)")+
  #ggtitle("Flavodoxin")+
  #theme_bw()+
  #theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
  #theme(legend.position = "none")+
  #scale_y_log10()















# #x,y plots of the different proteins
# xyplot <- function(dataset,title, yzoom1,yzoom2){
#   title <- title
#   #simple x-yplot 
#   print (ggplot(dataset, aes(x=variable, y=value, colour = name ))+
#            geom_point(alpha=0.9)+
#            stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.5)+
#            stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
#            xlab ("Treatmet")+
#            ggtitle(title)+
#            theme_bw()+
#            theme(plot.title = element_text(hjust = 0.5))+
#            theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
#            theme(legend.position = "none"))
#   
#   #log(10) of y-values because of large spread in values
#   print(ggplot(dataset, aes(x=variable, y=value, colour = name ))+
#           geom_point(alpha=0.9)+
#           stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.5)+
#           stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
#           xlab ("Treatmet")+
#           ylab ("Value(Log10)")+
#           ggtitle(title)+
#           theme_bw()+
#           theme(plot.title = element_text(hjust = 0.5))+
#           theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
#           scale_y_log10()+
#           theme(legend.position = "none"))
#   
#   #zooming in a little 
#   print(ggplot(dataset, aes(x=variable, y=value, colour = name ))+
#           geom_point(aes(group=name,alpha=0.9))+
#           stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.5)+
#           stat_summary(fun.data = mean_se, geom = "errorbar", colour = "black")+
#           xlab ("Treatmet")+
#           ggtitle(title)+
#           # theme_bw()+
#           # theme(plot.title = element_text(hjust = 0.5))+
#           # theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
#           # scale_y_continuous(limits = c(yzoom1,yzoom2))+
#           # theme(legend.position = "none"))
#           
# } 



























# #fixdata:this function creates subsets of the flavodoxin data and plots them as needed. 
# #the loop creates means of the three treatments
# 
# #Arguments:
# #group = group of proteins we want to work with
# #dataset = which dataset we are using (typically one that has been cleaned)
# #proteinname = whatever proteins we want to use
# 
# fixdata <- function(dataset, proteinname, title, title2, yzoom1,yzoom2){
#   
#   datasetmelt <-  melt (data = dataset,id.vars = "name", measure.vars = c(2:28))
#   xyplot(datasetmelt,title, yzoom1, yzoom2)
#   
#   colname <- "mean"
#   i=0
#   while (i < 9){
#     colname <- paste(colname, i)
#     col1 = 2+(i*3)
#     col2 = 4+(i*3)
#     dataset[[colname]] <- rowMeans( dataset[, c(col1:col2)] ) 
#     i <- i+1
#   }
#   
#   datasetmean <- dataset [-c(2:28)]
#   names(datasetmean) <- c("name", "T0", "B12_0.5C", "Fe_B12_0.5C","0.5C", "Fe_0.5C", "3C", "Fe_3C", "6C", "Fe_6C" )
#   
#   heatmap_unclustered(data, datasetmean, proteinname)
#   
#   datasetmeanmelt <-  melt (data = datasetmean,id.vars = "name", measure.vars = c(2:10)) 
#   
#   xyplot(datasetmeanmelt,title2, yzoom1, yzoom2)
# }










# 
# xyplots <- function(dataset, title, title2, yzoom1,yzoom2){
#   
#   datasetmelt <-  melt (data = dataset,id.vars = "name", measure.vars = c(2:28))
#   xyplot(datasetmelt,title, yzoom1, yzoom2)
#   
#   colname <- "mean"
#   i=0
#   while (i < 9){
#     colname <- paste(colname, i)
#     col1 = 2+(i*3)
#     col2 = 4+(i*3)
#     dataset[[colname]] <- rowMeans( dataset[, c(col1:col2)] ) 
#     i <- i+1
#   }
#   
#   datasetmean <- dataset [-c(2:28)]
#   names(datasetmean) <- c("name", "T0", "B12_0.5C", "Fe_B12_0.5C","0.5C", "Fe_0.5C", "3C", "Fe_3C", "6C", "Fe_6C" )
#   
#   datasetmeanmelt <-  melt (data = datasetmean,id.vars = "name", measure.vars = c(2:10)) 
#   
#   xyplot(datasetmeanmelt,title2, yzoom1, yzoom2)
# }
# 


#MCLclusterplot <- function(whichdata, whichprotein, graphtitle): This plots the mean of the three triplicates in a line, against the three temperatures with -/+ Fe in two different colors. 

# MCLclusterplot <- function(whichdata, whichprotein, graphtitle){
#   ISIPMCL <- dplyr::filter(whichdata, grepl (whichprotein, name))
#   ISIPMCL <- ISIPMCL [c(76, 1:75)]
#   ISIPMCL <- ISIPMCL [-c(2:28, 47:76)]
#   names(ISIPMCL)<- c("name", "A_0.5", "B_0.5", "C_0.5", "A_Fe_0.5", "B_Fe_0.5", "C_Fe_0.5", "A_3", "B_3", "C_3", "A_Fe_3", "B_Fe_3", "C_Fe_3","A_6", "B_6", "C_6", "A_Fe_6", "B_Fe_6", "C_Fe_6")
#   
#   colname <- "mean"
#   i=0
#   while (i < 6){
#     colname <- paste(colname, i)
#     col1 = 2+(i*3)
#     col2 = 4+(i*3)
#     ISIPMCL[[colname]] <- rowMeans( ISIPMCL[, c(col1:col2)] )
#     i <- i+1
#   }
#   
#   colname <- "SD"
#   i=0
#   while (i < 6){
#     colname <- paste(colname, i)
#     col1 = 2+(i*3)
#     col2 = 4+(i*3)
#     ISIPMCL[[colname]] <- apply ( ISIPMCL[, c(col1:col2)],1, sd )
#     i <- i+1
#   }
# #   
#   ISIPMCL <- ISIPMCL [-c(2:19)]
#   names(ISIPMCL) <- c("name","mean_NoFe_0.5", "mean_Fe_0.5", "mean_NoFe_3", "mean_Fe_3", "mean_NoFe_6", "mean_Fe_6","SD_NoFe_0.5", "SD_Fe_0.5", "SD_NoFe_3", "SD_Fe_3", "SD_NoFe_6", "SD_Fe_6")
#   
#   ISIPMCL <-  melt (data = ISIPMCL,id.vars = "name")
#   ISIPMCL <- separate(ISIPMCL, col = variable, into = c("Statistic", "Fe", "temperature"))
#   ISIPMCL <- spread(ISIPMCL, Statistic, value)
#   
#   ggplot(ISIPMCL, aes(x=temperature, y=mean, colour = Fe ))+
#     geom_errorbar(aes(ymin = mean-SD, ymax = mean+SD), width = 0.5,   color = "black", size = 0.7)+
#     geom_point(size = 3)+
#     geom_line(size = 1, aes (x=temperature, y=mean, group = interaction(Fe, name)))+
#     scale_color_manual(values = c("red", "black"))+
#     xlab ("Temperature (°C)")+
#     ylab ("Expression-Value")+
#     ggtitle(graphtitle)+
#     theme_bw()+
#     #labs (color = "Treatment") +
#     theme(plot.title = element_text(hjust = 0.5))+
#     theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))
#   #guides(colour = "Treatment")
# }
# 
# ###.....................................................................................................
# 
# 
# ### ORF Plot. Similar to MCL cluster plot, this here plots all the ORFS in a cluster against the treatments
# 
# ORFplots <- function(dataset, whichcluster, graphtitle ){
# 
# cluster <- dplyr::filter(dataset, grepl (whichcluster, cluster))
# 
# #cluster <- cluster [c(107, 1:106)]
# cluster <- cluster [-c(2:64, 83:107)]
# 
# names(cluster)<- c("name", "A_0.5", "B_0.5", "C_0.5", "A_Fe_0.5", "B_Fe_0.5", "C_Fe_0.5", "A_3", "B_3", "C_3", "A_Fe_3", "B_Fe_3", "C_Fe_3","A_6", "B_6", "C_6", "A_Fe_6", "B_Fe_6", "C_Fe_6")
# 
# cluster <-  melt (data = cluster ,id.vars = "name", measure.vars = c(2:19))
# 
# ggplot(cluster, aes(x=variable, y=value, colour = name ))+
#   geom_point()+
#   geom_line(aes(group=name))+
#   # stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.25)+
#   #stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
#   xlab ("Treatmet")+
#   ylab ("Expression-Value")+
#   ggtitle(graphtitle)+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5))+
#   theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))
# }
# 
# 
# 
# 
# ##...................................................................................
# 
# ## this plots specifc contigs
# contigplots <- function(dataset, whichcontig, graphtitle ){
#   
#   contig <- dplyr::filter(dataset, grepl (whichcontig, orf_id))
#   
#   T0 <- contig[1, 83]
#   T3 <- contig[1, 101]
#   T6 <- contig[1, 105]
#   contig <- contig [-c(2:64, 83:107)]
#   
#   names(contig)<- c("name", "A_noFe_0", "B_noFe_0", "C_noFe_0", "A_Fe_0", "B_Fe_0", "C_Fe_0", "A_noFe_3", "B_noFe_3", "C_noFe_3", "A_Fe_3", "B_Fe_3", "C_Fe_3","A_noFe_6", "B_noFe_6", "C_noFe_6", "A_Fe_6", "B_Fe_6", "C_Fe_6")
#   
#   contig <-  melt (data = contig ,id.vars = "name", measure.vars = c(2:19))
#   
#   
#   contig <- separate(data=contig, col = variable, into = c("Replicate", "Fe", "temperature"))
#   
#   ggplot(contig, aes(x=temperature, y=value, colour = Fe, shape = Replicate ))+
#     geom_point(size = 3)+
#     scale_color_manual(values = c("red", "black"))+
#     xlab ("Temperature (°C)")+
#     ylab ("Expression-Value")+
#     ggtitle(graphtitle)+
#     theme_bw()+
#     theme(plot.title = element_text(hjust = 0.5))+
#     theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
#     scale_color_manual (name = "", labels = c("+Fe", "-Fe"), values = c("red", "black"))+
#     annotate ("text", x = 1:3, y = 0, label= c(T0, T3, T6), color="darkgreen", size = 4, fontface="bold")+
#     theme(axis.text.x = element_text (color= "black", face = "bold", size = 11))+
#     theme(axis.text.y = element_text (color= "black", face = "bold"))+
#     theme(axis.title.y = element_text (color= "black", face = "bold", size = 12))+
#     theme(axis.title.x = element_text (color= "black", face = "bold", size = 12))
# }

#``````````````````````````````````````````````````````````````````````````````````````

