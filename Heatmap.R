library(maditr)
library(FactoMineR)
library(ggfortify)
library(ggplot2)
library(matlib)
library(readxl)
library(gplots)
library(dendextend)
library(ComplexHeatmap)

# paths of data to use
path <- "C:\\Users\\andre\\OneDrive - Fundación Universitaria San Pablo CEU\\miRNAs\\Resultados_miRNAs\\8_RESULTS_GME.txt"
path_significative <- "C:\\Users\\andre\\OneDrive - Fundación Universitaria San Pablo CEU\\miRNAs\\Resultados_miRNAs\\miRNAs_validation.xlsx"

data <- read.table(path, header = TRUE)
data_validation <- read_excel(path_significative)

miRNAs_sig <- data_validation[,1]
transformed_data <- as.data.frame(dcast(data, miRNA ~ Plate_ID, value.var = "X2..AACt", 
                                        fun.aggregate=sum))
rownames(transformed_data) <- transformed_data[,1]
transformed_data <- transformed_data[,-c(1)]

data_to_heatmap <- transformed_data[c(miRNAs_sig$miRNA),]

# control samples are deleted
data_to_heatmap <- data_to_heatmap[,-c(14:28)]

to_heatmap_Z <- data_to_heatmap

# Z-score normalization
for (row in 1:nrow(to_heatmap_Z)){
  for (col in 1:ncol(to_heatmap_Z)){
    to_heatmap_Z[row, col] <- 
      (as.numeric(to_heatmap_Z[row, col])-mean(as.numeric(data_to_heatmap[row, ])))/sd(as.numeric(data_to_heatmap[row, ]))
  }
  
}

# clustering
col_dend = hclust(dist(t(as.matrix((to_heatmap_Z))), 
                       method = "manhattan"), method = "complete") 

# heatmap
Heatmap(as.matrix((to_heatmap_Z)), 
        row_names_gp = gpar(fontsize = 10),
        name = "Expression values",
        cluster_columns = color_branches(col_dend, k = NULL, col = "black"),
        show_row_dend = FALSE, 
        column_names_gp = gpar(col = c(rep("orange", 11), rep("purple", 10)))
        )
