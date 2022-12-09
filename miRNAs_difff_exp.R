library("car")
library("ggplot2")
library("rstatix")
library("tidyverse")
library("VennDiagram")
library("WRS2")


### FUNCTIONS

data_collection <- function(file_name){
  
  ## data_collection ----------------------------------------------------
  ## Function that gets the data to analyze. 
  ## file_name: name of file that contains the data to be analyzed. 
  ## Returns a list with the data.
  
  output_list <- list()
  
  path <- getwd()
  data <- read.delim(file = paste(path, file_name, sep = ""), skip = 1)
  
  output_list[["cts"]] <- c(data[,3])
  output_list[["miRNAs"]] <- c(data[,2])
  output_list[["groups"]] <- c(data[,5])
  
  return(output_list)
  
}

statistical_tests_pvalues <- function(cts, groups = NULL, test){

  ## statistical_tests_pvalues -----------------------------------------------
  ## Function that does different statistical tests: 
  ## normality: saphiro.test function
  ## homoscedasticity: leveneTest function
  ## Kruskal-Wallis: kruskal.test function
  ## Pairwise Wilcoxon Rank Sum Tests: pairwise.wilcox.test function
  ## ANOVA: aov function
  ## Pairwise t tests: pairwise.t.test function
  
  ## cts: miRNA cts that will take part of the analysis.
  ## groups: miRNA groups that will take part of the analysis.
  ## It returns a vector with the corresponding p-values.

  # The vector to be returned
  pvalues <- c()

  # NORMALITY
  if (test == "normality"){
  
    # There must be at least 3 values to do this test.
    if (length(cts) > 2){ 
      pvalue_norm <- tryCatch(shapiro.test(unlist(as.numeric(cts)))$p.value, 
                             error = function(e){pvalue_norm = NULL})
        
      if (class(pvalue_norm[1]) != "numeric"){
        pvalue_norm <- NA
      }
    }
        
    else { 
      pvalue_norm <- NA
    }
        
    pvalues <- pvalue_norm

  }
  
  # HOMOCEDASTICITY
  if (test == "homocedasticity"){
    
    # There must be at least 2 groups to do this test.
    if(length(unique(groups))>= 2){ 
      
      pvalue_levene <- leveneTest(y = cts, group = factor(groups), center =
                                   "median")$Pr[1]
    }
    else{
      pvalue_levene <- NA
    }
      
    pvalues <- pvalue_levene
  }
  
  # KRUSCAL-WALLIS
  if (test == "KW"){
    if (length(unique(groups))>= 2){
      pvalue_KW <- tryCatch(kruskal.test(cts, groups)$p.value,
                            error = function(e){pvalue_norm = NA})
    }
    else {pvalue_KW <- NA}
      
     pvalues <-  pvalue_KW
  }
    
  # ANOVA
  if (test == "ANOVA"){
    if (length(unique(groups))>= 2){
      anova_test <- aov(cts ~ groups)
      pvalue_ANOVA <- summary(anova_test)[[1]][["Pr(>F)"]][[1]]
    }
    else {pvalue_ANOVA <- NA}
      
    pvalues <- pvalue_ANOVA
      
  }
  
  # Pairwise Wilcoxon Rank Sum Tests 
  if (test == "post-hoc_KW"){
    if (length(unique(groups)) >= 2){
      pvalue_post_hoc_KW <- pairwise.wilcox.test(cts, groups, 
                                                p.adjust.method = "fdr")$p.value
    }
    pvalues <- pvalue_post_hoc_KW
  }
  
  # Pairwise t tests
  if (test == "post-hoc_ANOVA"){
    if (length(unique(groups)) >= 2){
      pvalue_post_hoc_ANOVA <- pairwise.t.test(cts, groups, 
                                                p.adjust.method = "fdr")$p.value
    }
      
    pvalues <- pvalue_post_hoc_ANOVA
  }
  

  return(pvalues)
  
}


main <- function(file_name, method){

  ## main --------------------------------------------------------------------
  ## Function that executes the main code. 
  ## file_name: name of file that contains the data to be analyzed.
  ## method: "classic" statistical methods. 
  
  if (file.exists(paste(getwd(),file_name, sep="")) == FALSE){
    stop("Introduced data file does not exist.")
  }
  
  data <- data_collection(file_name)
  
  if (method != "classic"){

    stop("Method not recognized. Please, introduce 'classic' method.")

  }
  
  else if (method == "classic"){
    
    # Results file header
    write.table(data.frame("miRNA", "ICS_vs_C_pvalue", "ICS_vs_C_log2(FC)",
                           "UC_vs_C_pvalue", "UC_vs_C_log2(FC)",
                           "UC_vs_ICS_pvalue",
                           "UC_vs_ICS_log2(FC)"), file = 
                  paste("ANOVA_KW_differences", str_sub(file_name, 11, ),
                                                          sep = ""),
                append = FALSE, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    
    # Each miRNA population is studied separately.
    for (i in unique(data$miRNAs)){
      # Spikes are not analyzed
      if (startsWith(i, "hsa") == TRUE | startsWith(i, "cel") == TRUE ){
        pvalues_norm <- c()
        pvalue_levene <- ""
        pvalue_ANOVA <- ""
        results_post_hoc_ANOVA <- c()
        pvalue_KW <- ""
        results_post_hoc_KW <- c()
        
        # Indexes corresponding to miRNA i, regardless the group
        all_miRNA_index <- which(data$miRNAs == i) 
        
        all_miRNA_cts <- data$cts[all_miRNA_index]
        all_miRNA_groups <- data$groups[all_miRNA_index]
        
        for (j in c("ICS", "UC", "C")){
          
          # Indexes corresponding to miRNA i, regarding the group j
          group_miRNA_index <- which(data$miRNAs == i & data$groups == j)
    
          group_miRNA_cts <- data$cts[group_miRNA_index]
          
          # If n >= 7, you can aim to do an ANOVA if the data passes the
          # normality and homoscedasticity tests.
          if (length(group_miRNA_cts) >= 7){ 
            
            # Normality test p-values
            pvalues_norm <- c(pvalues_norm, statistical_tests_pvalues(
              group_miRNA_cts, test = "normality"))
    
          }
        }
        
    
        # If the miRNA data passes the normality test in the 3 groups...
        if (is.null(pvalues_norm) == FALSE){
          if (length(pvalues_norm[which(pvalues_norm > 0.05)]) == 3){ 
      
            # Homocedasticity test p-value
            pvalue_levene <- statistical_tests_pvalues(all_miRNA_cts,
                                                      groups = all_miRNA_groups,
                                                      "homocedasticity")
            
          
            # If the miRNA data passes the homocedasticity test...
            if (pvalue_levene > 0.05){
             
              # ANOVA p-value 
              pvalue_ANOVA <- statistical_tests_pvalues(all_miRNA_cts, groups =
                                                         all_miRNA_groups,
                                                        "ANOVA")
              
              if (pvalue_ANOVA <= 0.05 & is.nan(pvalue_ANOVA) == FALSE){
                
                results_post_hoc_ANOVA = statistical_tests_pvalues(
                  all_miRNA_cts, groups = all_miRNA_groups,"post-hoc_ANOVA")
              }
            }
          }
        }
        
        # If the miRNA data does not pass normality or homocedasticity tests...
        if (pvalue_ANOVA == ""){
          
          # KW test p-value
          pvalue_KW <- statistical_tests_pvalues(all_miRNA_cts, groups =
                                                  all_miRNA_groups, "KW")
          if (is.na(pvalue_KW) == FALSE){
            if (pvalue_KW <= 0.05){
              results_post_hoc_KW = statistical_tests_pvalues(all_miRNA_cts,
                                                              groups = 
                                                                all_miRNA_groups,
                                                              "post-hoc_KW")
            }
          }
        }
        
        # Fold Change
        cts_ICS_mean <- mean(all_miRNA_cts[which(all_miRNA_groups == "ICS")])
        cts_UC_mean <- mean(all_miRNA_cts[which(all_miRNA_groups == "UC")])
        cts_C_mean <- mean(all_miRNA_cts[which(all_miRNA_groups == "C")])
        
        
        # Only if one value is <= 0.05, the results are got
        if (length(which(results_post_hoc_ANOVA <= 0.05)) != 0){ 
          
          write.table(data.frame(i, results_post_hoc_ANOVA[1],
                                 log2(cts_ICS_mean/cts_C_mean), 
                                 results_post_hoc_ANOVA[2],
                                 log2(cts_UC_mean/cts_C_mean),
                                 results_post_hoc_ANOVA[4], 
                                 log2(cts_UC_mean/cts_ICS_mean)),
                      file = paste("ANOVA_KW_differences", 
                                   str_sub(file_name, 11, ),
                                   sep = ""), append = TRUE, sep = "\t",
                      row.names = FALSE, col.names = FALSE, quote = FALSE)
          
        }
        
        else if (length(which(results_post_hoc_KW <= 0.05)) != 0){ 
          
          write.table(data.frame(i, results_post_hoc_KW[1],
                                 log2(cts_ICS_mean/cts_C_mean),
                                 results_post_hoc_KW[2],
                                 log2(cts_UC_mean/cts_C_mean),
                                 results_post_hoc_KW[4],
                                 log2(cts_UC_mean/cts_ICS_mean)),
                      file =  paste("ANOVA_KW_differences",
                                    str_sub(file_name, 11, ),
                                    sep = ""), append = TRUE, sep="\t",
                      row.names = FALSE, col.names = FALSE, quote = FALSE)
          
        }
      }
    }
  }
}

### CODE

main("\\8_RESULTS_GME.txt", "classic")
