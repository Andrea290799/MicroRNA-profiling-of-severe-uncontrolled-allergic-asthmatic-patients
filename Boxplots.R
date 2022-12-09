library("ggplot2")
library("gridExtra")
library("ggpubr")


### FUNCTIONS

data_collection <- function(file_name){
  
  ## data_collection ----------------------------------------------------
  ## Function that gets the data to analyze. 
  ## file_name: name of file that contains the data to be analyzed. 
  ## Returns a list with the data.
  
  output_list <- list()
  
  data <- read.delim(file = file_name, skip = 1)
  
  cts <- output_list[["cts"]]
  miRNAs <- output_list[["miRNAs"]]
  groups <- output_list[["groups"]]
  patients <- output_list[["patients"]]
  
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
  important_miRNAs <- c("hsa-miR-548j-5p", "hsa-miR-554",
                        "hsa-miR-876-3p", "hsa-miR-641",
                        "hsa-miR-877-5p", "hsa-miR-509-3p",
                        "hsa-miR-665", "hsa-miR-32-3p",
                        "hsa-miR-31-5p", "hsa-miR-671-5p",
                        "hsa-miR-33b-5p","hsa-let-7c-5p",
                        "hsa-miR-10b-5p", "hsa-miR-139-3p",
                        "hsa-miR-627-5p", "hsa-miR-30d-5p", 
                        "hsa-miR-185-5p", "hsa-miR-326")
  
  
  to_keep <- c()

  # outliers are deleted for representation
  to_remove <- which(cts > 40)

  miRNAs <- miRNAs[-to_remove]
  cts <- cts[-to_remove]
  groups <- groups[-to_remove]
  patients <- patients[-to_remove]

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

path <- "C:\\Users\\andre\\OneDrive - Fundación Universitaria San Pablo CEU\\miRNAs\\Resultados_miRNAs\\8_RESULTS_GME.txt"

if (file.exists(path) == FALSE){
  stop("Introduced data file does not exist.")
}

data <- data_collection(path)

par(mfrow = c(3,6))

for (i in unique(data$miRNAs)){

  cts <- data$cts[data$miRNAs == i]
  
  miRNAs <- data$miRNAs[data$miRNAs == i]
  Grupos <- data$groups[data$miRNAs == i]
  
  # nomenclature changes
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





    