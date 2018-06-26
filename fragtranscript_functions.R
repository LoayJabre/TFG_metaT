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

#heatmapfunction: group = group of proteins we want to work with
#dataset = which dataset we are using (typically one that has been cleaned)
#proteinname = whatever proteins we want to use