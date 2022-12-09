library(maditr)
library(FactoMineR)
library(ggfortify)
library(matlib)
library(readxl)
library(ggpubr)
library(ggcorrplot)

# paths of used files
path_miRNAs <- "C:\\Users\\andre\\OneDrive - Fundación Universitaria San Pablo CEU\\miRNAs\\Resultados_miRNAs\\8_RESULTS_GME.txt"
path_metabolites <- "C:\\Users\\andre\\OneDrive - Fundación Universitaria San Pablo CEU\\miRNAs\\Metabolitos.xlsx"
path_significative <- "C:\\Users\\andre\\OneDrive - Fundación Universitaria San Pablo CEU\\miRNAs\\Resultados_miRNAs\\miRNAs_validation.xlsx"

# data obtention: differentially expressed miRNAs and metabolites
miRNA_data <- read.table(path_miRNAs, header = TRUE)
metabolites_data <- as.data.frame(read_excel(path_metabolites, col_names=TRUE))

# miRNA data trasnfromation to a shorter format
miRNA_data <- as.data.frame(dcast(miRNA_data, miRNA ~ Plate_ID, value.var = "X2..AACt", 
                                        fun.aggregate=sum))

# control subjects are not taken into account
miRNA_data <- miRNA_data[,-c(15:29)]

rownames(miRNA_data) <- miRNA_data[,1]
miRNA_data <- miRNA_data[,-c(1)]

data_validation <- read_excel(path_significative)

miRNAs_sig <- data_validation[,1]

miRNA_data <- miRNA_data[c(miRNAs_sig$miRNA),]

rownames(metabolites_data) <- metabolites_data[,1]
metabolites_data <- metabolites_data[,-c(1)]

# vectors for final correlation file
variable_1 <- c()
variable_2 <- c()
pvalue <- c()
r <- c()

for (miRNA in rownames(miRNA_data)){
  for (metabolite in rownames(metabolites_data)){
    
    miRNA_population <- as.numeric(c(miRNA_data[miRNA,]))
    metabolites_population <- as.numeric(c(metabolites_data[metabolite,]))
    
    #normality tests
    shapiro_test_miRNAs <- shapiro.test(miRNA_population)
    shapiro_test_metabolites <- shapiro.test(metabolites_population)
    
    
    if (shapiro_test_miRNAs$p.value > 0.05 &
        shapiro_test_metabolites$p.value > 0.05){
      
      corr_test <- cor.test(miRNA_population, metabolites_population, 
               method = "pearson")
      
      method = "pearson"
      
      
    
    } else {
    
      corr_test <- cor.test(miRNA_population, metabolites_population, 
                                    method = "spearman")
      method = "spearman"
      
      
    }
    
    variable_1 <- c(variable_1,miRNA)
    variable_2 <- c(variable_2,metabolite)
    pvalue <- c(pvalue, corr_test$p.value)
    r <- c(r, corr_test$estimate)
    
    if (corr_test$p.value < 0.05 & abs(corr_test$estimate) > 0){
      
      data_to_plot <- cbind(miRNA_population, metabolites_population)
      
      colnames(data_to_plot) <- c(gsub("-", "_", miRNA), gsub("-", "_", metabolite))

      ggscatter(data_to_plot, x = gsub("-", "_", miRNA), 
                y = gsub("-", "_", metabolite),
                add = "reg.line", conf.int = TRUE,
                cor.coef = TRUE, cor.method = method,
                xlab = miRNA, ylab = metabolite)

      ggsave(paste("C:\\Users\\andre\\OneDrive - Fundación Universitaria San Pablo CEU\\miRNAs\\Correlation_plots\\",
                   variable_1[length(variable_1)], "-", 
                   gsub(":", "-", variable_2[length(variable_2)]),
                   ".pdf", sep=""), device = "pdf")

    }
  }
}

# correlations
final_df <- data.frame(variable_1, variable_2, pvalue, r)

write.table(final_df, "C:\\Users\\andre\\OneDrive - Fundación Universitaria San Pablo CEU\\miRNAs\\Correlaciones_miRNAs-metabolitos.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# short-format correlation information (r)
transformed_data <- as.data.frame(dcast(final_df, variable_1 ~ variable_2, value.var = "r", 
                                        fun.aggregate=sum))

# short-format correlation information (p-val)
transformed_data_pvalues <- as.data.frame(dcast(final_df, variable_1 ~ variable_2, value.var = "pvalue", 
                                        fun.aggregate=sum))

rownames(transformed_data) <- c(transformed_data[,1])
transformed_data <- transformed_data[,-c(1)]

rownames(transformed_data_pvalues) <- c(transformed_data_pvalues[,1])
transformed_data_pvalues <- transformed_data_pvalues[,-c(1)]

# big correlation plot
ggcorrplot(transformed_data, method = "circle", 
           p.mat = as.matrix(transformed_data_pvalues), 
           sig.level = 0.05, insig = "blank")

# small correlation plot
miRNAs_with_targets <- c("hsa-miR-548j-5p", "hsa-miR-554",
                         "hsa-miR-876-3p", "hsa-miR-641",
                         "hsa-miR-877-5p", "hsa-miR-509-3p",
                         "hsa-miR-665", "hsa-miR-32-3p",
                         "hsa-miR-31-5p", "hsa-miR-671-5p",
                         "hsa-miR-33b-5p","hsa-let-7c-5p",
                         "hsa-miR-10b-5p", "hsa-miR-139-3p",
                         "hsa-miR-627-5p", "hsa-miR-30d-5p", 
                         "hsa-miR-185-5p", "hsa-miR-326")

important_metabolites <- c("Arachidonic Acid", "L-Arginine", "L-Leucine", 
                           "Sphingosine Phosphate", "LPC 20:4", "LPC 20:5",
                           "LPC 22:6")

transformed_data <- as.matrix(transformed_data)[miRNAs_with_targets, important_metabolites]
transformed_data_pvalues <- as.matrix(transformed_data_pvalues)[miRNAs_with_targets, important_metabolites]

ggcorrplot(transformed_data, method = "circle", 
           p.mat = as.matrix(transformed_data_pvalues), 
           sig.level = 0.05, insig = "blank")


