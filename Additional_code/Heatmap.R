library(maditr)
library(FactoMineR)
library(ggfortify)
library(ggplot2)
library(matlib)
library(readxl)
library(gplots)
library(dendextend)
library(ComplexHeatmap)

# output of miRNAs_normalization.R
data <- read.table("8_RESULTS_GME.txt", header = TRUE)

# output of miRNAs_diff_exp
DE_miRNAs <- read.table("ANOVA_KW_differences_GME-DE.txt", header = T)
DE_miRNAs <- DE_miRNAs[order(abs(DE_miRNAs$UC_vs_ICS_log2.FC.), decreasing = T),]
miRNAs_sig <- unlist(DE_miRNAs[which(DE_miRNAs[,6] < 0.05), 1])

# transform to table
transformed_data <- as.data.frame(dcast(data, miRNA ~ Plate_ID, value.var = "X2..AACt", 
                                        fun.aggregate=sum))

rownames(transformed_data) <- transformed_data[,1]
transformed_data <- transformed_data[,-c(1)]

# only DE miRNA are represented
data_to_heatmap <- transformed_data[c(miRNAs_sig),]

# control samples are deleted
data_to_heatmap <- data_to_heatmap[,-c(14:28)]


# Z-score normalization
to_heatmap_Z <- data_to_heatmap

for (row in 1:nrow(to_heatmap_Z)){
  for (col in 1:ncol(to_heatmap_Z)){
    to_heatmap_Z[row, col] <- 
      (as.numeric(to_heatmap_Z[row, col])-
         mean(as.numeric(data_to_heatmap[row, ])))/
      sd(as.numeric(data_to_heatmap[row, ]))
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
