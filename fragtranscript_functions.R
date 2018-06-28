#heatmapfunction: group = group of proteins we want to work with
#dataset = which dataset we are using (typically one that has been cleaned)
#proteinname = whatever proteins we want to use

#heatmap function
heatmap <- function (group, dataset, proteinname){
  group <- dplyr::filter (dataset, grepl(proteinname, name))
  group <- data.frame(group[,-1], row.names = group [,1])
  group <- as.matrix ((group))
  heatmap.2 (group, trace = "none", density = "none", col = bluered(100), cexRow = 0.8 , 
             cexCol =0.8, margins = c(5,12), scale = "row")
             }

#heatmap unclustered function
heatmap_unclustered <- function (group, dataset, proteinname){
  group <- dplyr::filter (dataset, grepl(proteinname, name))
  group <- data.frame(group[,-1], row.names = group [,1])
  group <- as.matrix ((group))
  heatmap.2 (group, trace = "none", density = "none", col = bluered(100), cexRow = 0.8, 
             cexCol =1, margins = c(8,9), scale = "row" , dendrogram = 'none',  
             Colv = FALSE, Rowv = FALSE)
             }


#x,y plots of the different proteins
xyplot <- function(dataset,title, yzoom){
  title <- title
  yzoom <- yzoom
  #simple x-yplot 
  print (ggplot(dataset, aes(x=variable, y=value, colour = name ))+
    geom_point(alpha=0.9)+
    stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.5)+
    stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
    xlab ("Treatmet")+
    ggtitle(title)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
    theme(legend.position = "none"))
  
#log(10) of y-values because of large spread in values
  print(ggplot(dataset, aes(x=variable, y=value, colour = name ))+
    geom_point(alpha=0.9)+
    stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.5)+
    stat_summary(fun.data =  mean_se, geom = "errorbar", colour = "black")+
    xlab ("Treatmet")+
    ylab ("Value(Log10)")+
    ggtitle(title)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
    scale_y_log10()+
    theme(legend.position = "none"))

#zooming in a little 
  print(ggplot(dataset, aes(x=variable, y=value, colour = name ))+
    geom_point(aes(group=name,alpha=0.9))+
    stat_summary(fun.y = mean, geom = "line", colour = "black", group =1, size = 1.5)+
    stat_summary(fun.data = mean_se, geom = "errorbar", colour = "black")+
    xlab ("Treatmet")+
    ggtitle(title)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0))+
    scale_y_continuous(limits = c(0,yzoom))+
    theme(legend.position = "none"))
} 



#this function creates subsets of the flavodoxin data and plots them as needed 
#the loop creates means of the three treatments
fixdata <- function(dataset, proteinname, title, title2){
  heatmap_unclustered(data, dataset, proteinname)
  
  datasetmelt <-  melt (data = dataset,id.vars = "name", measure.vars = c(2:28))
  xyplot(datasetmelt,title, 0.001)
  
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
  
  heatmap_unclustered(data, datasetmean, proteinname)
  
  datasetmeanmelt <-  melt (data = datasetmean,id.vars = "name", measure.vars = c(2:10)) 
  
  xyplot(datasetmeanmelt,title2, 0.002)
  
}



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





