library(maditr)
library(FactoMineR)
library(ggfortify)
library(matlib)
library(readxl)
library(ggpubr)
library(ggcorrplot)

# paths of used files

# output of miRNAs_normalization.R
path_miRNAs <- ".\\Resultados_miRNAs\\8_RESULTS_GME.txt"

# table of metabolite abundances with metabolites in rows and patients in columns
path_metabolites <- ".\\Metabolitos.xlsx"

# output of miRNAs_diff_exp
DE_miRNAs <- read.table("ANOVA_KW_differences_GME-DE.txt", header = T)
DE_miRNAs <- DE_miRNAs[order(abs(DE_miRNAs$UC_vs_ICS_log2.FC.), decreasing = T),]
miRNAs_sig <- unlist(DE_miRNAs[which(DE_miRNAs[,6] < 0.05), 1])

# data obtention: miRNAs and metabolites
miRNA_data <- read.table(path_miRNAs, header = TRUE)
metabolites_data <- as.data.frame(read_excel(path_metabolites, col_names=TRUE))

# miRNA data transformation to a shorter format
miRNA_data <- as.data.frame(dcast(miRNA_data, miRNA ~ Plate_ID, value.var = "X2..AACt", 
                                        fun.aggregate=sum))

# control subjects are not taken into account
miRNA_data <- miRNA_data[,-c(15:29)]

rownames(miRNA_data) <- miRNA_data[,1]
miRNA_data <- miRNA_data[,-c(1)]

# only DE miRNAs
miRNA_data <- miRNA_data[miRNAs_sig,]

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
  }
}

# correlations
final_df <- data.frame(variable_1, variable_2, pvalue, r)

write.table(final_df, ".\\Correlaciones_miRNAs-metabolitos.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# short-format correlation information (r)
transformed_data <- as.data.frame(dcast(final_df, variable_1 ~ variable_2, value.var = "r"))

# short-format correlation information (p-val)
transformed_data_pvalues <- as.data.frame(dcast(final_df, variable_1 ~ variable_2, value.var = "pvalue"))

rownames(transformed_data) <- c(transformed_data[,1])
transformed_data <- transformed_data[,-c(1)]

rownames(transformed_data_pvalues) <- c(transformed_data_pvalues[,1])
transformed_data_pvalues <- transformed_data_pvalues[,-c(1)]


# correlation plot

important_metabolites <- rev(c("Sphingosine Phosphate", "Arachidonic Acid", 
                           "LPC 20:4", "LPC 20:5","LPC 22:6", "L-Arginine",
                           "L-Leucine"))

transformed_data <- as.matrix(transformed_data)[, important_metabolites]
transformed_data_pvalues <- as.matrix(transformed_data_pvalues)[, important_metabolites]

# remove miRNAs with no correlations 

miRNAs_no_corr <- c("hsa-let-7f-1-3p", "hsa-miR-139-3p", "hsa-miR-145-5p", 
                    "hsa-miR-18a-5p", "hsa-miR-31-5p", "hsa-miR-323a-3p",
                    "hsa-miR-326", "hsa-miR-33b-5p", "hsa-miR-451a", 
                    "hsa-miR-505-3p", "hsa-miR-521", "hsa-miR-590-3p", 
                    "hsa-miR-631", "hsa-miR-641")

miRNAs_no_corr_indexes <- c()
for (i in miRNAs_no_corr){
  miRNAs_no_corr_indexes <- c(miRNAs_no_corr_indexes,
                              which(rownames(transformed_data)== i))
}


transformed_data <- transformed_data[-miRNAs_no_corr_indexes,]
transformed_data_pvalues <- transformed_data_pvalues[-miRNAs_no_corr_indexes,]

ggcorrplot(transformed_data, method = "circle", 
           p.mat = as.matrix(transformed_data_pvalues), 
           sig.level = 0.05, insig = "blank")

write.table(transformed_data,
            file = ".\\cor_coefs.txt", 
            sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)

write.table(transformed_data_pvalues,
            file = ".\\cor_pvalues.txt", 
            sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)

ggcorrplot(transformed_data, method = "circle", 
           p.mat = as.matrix(transformed_data_pvalues), 
           sig.level = 0.05, insig = "blank")


