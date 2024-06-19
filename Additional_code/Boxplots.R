library("ggplot2")
library("gridExtra")
library("ggpubr")
library(readxl)

### FUNCTIONS

data_collection <- function(file_name){
  
  ## data_collection ----------------------------------------------------
  ## Function that gets the data to analyze. 
  ## file_name: name of file that contains the data to be analyzed. 
  ## Returns a list with the data.
  
  output_list <- list()
  
  data <- read.delim(file = file_name, skip = 1)
  
  cts <- c(data[,3])
  miRNAs <- c(data[,2])
  groups <- c(data[,5])
  patients <- c(data[,1])
  
  # control group is not analyzed
  to_remove <- which(groups == "C")
  
  miRNAs <- miRNAs[-to_remove]
  cts <- cts[-to_remove]
  groups <- groups[-to_remove]
  patients <- patients[-to_remove]
  
  # miRNAs chosen for boxplot representation
  
  # output of miRNAs_diff_exp
  DE_miRNAs <- read.table("ANOVA_KW_differences_GME-DE.txt", header = T)
  DE_miRNAs <- DE_miRNAs[order(abs(DE_miRNAs$UC_vs_ICS_log2.FC.), decreasing = T),]
  
  #only DE miRNAs with log2FC >= 1.5
  important_miRNAs <- unlist(DE_miRNAs[which(abs(DE_miRNAs[,7]) >= 1.5 & DE_miRNAs[,6] < 0.05), 1])
  
  # only "important_miRNAs" are used
  to_keep <- c()

  for (i in important_miRNAs){
    to_keep <- c(to_keep, which(miRNAs == i))
  }

  miRNAs <- miRNAs[to_keep]
  cts <- cts[to_keep]
  groups <- groups[to_keep]
  patients <- patients[to_keep]
  
  output_list[["cts"]] <- cts
  output_list[["miRNAs"]] <- miRNAs
  output_list[["groups"]] <- groups
  output_list[["patients"]] <- patients

  return(output_list)
  
}

### CODE

# output of miRNAs_normalization.R
path <- "8_RESULTS_GME.txt"

if (file.exists(path) == FALSE){
  stop("Introduced data file does not exist.")
}

data <- data_collection(path)

par(mfrow = c(3,3))

for (i in unique(data$miRNAs)){

  cts <- data$cts[data$miRNAs == i]
  
  miRNAs <- data$miRNAs[data$miRNAs == i]
  Grupos <- data$groups[data$miRNAs == i]
  
  # group nomenclature changes
  for (j in 1:length(Grupos)){
    
    if (Grupos[j] == "ICS"){
      
      Grupos[j] = "M"
      
    }
    
    else{
      
      Grupos[j] = "SU"
      
    }
  }
  
  df = data.frame(Grupos, cts)
  
  # boxplots
  

  par(cex.lab=1.5)
  boxplot(df$cts ~ factor(Grupos),col = "white", border = c("orange", "purple"),
          cex.main= 2, main = i, cex.ylab = 2, xlab = "", 
          ylab = "Expression values", medlwd = 1, cex.axis = 1.5,
          cex.names = 1.5)
  
  stripchart(df$cts ~ factor(Grupos), vertical = TRUE, pch = 19, add = TRUE,
             col = "black")
  



}




    