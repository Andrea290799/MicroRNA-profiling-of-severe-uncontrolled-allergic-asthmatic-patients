library(readxl)

#data
miEAA_ICS_vs_UC_Analysis_results <- read_excel(".\\Resultados_miRNAs\\miEAA_ICS_vs_UC-Analysis_ results.xlsx", 
                                               sheet = "Tabla top 10")

output_list = list()
all_pathways = c()
all_groups=c()
available_miRNAs = c()

for (i in 1:length(miEAA_ICS_vs_UC_Analysis_results[[1]])){
  miRNAs = c()
  
  pathways=miEAA_ICS_vs_UC_Analysis_results[[1]][i]
  groups = miEAA_ICS_vs_UC_Analysis_results[[2]][i]
  
  no_lista_miRNAs = unlist(strsplit(as.character(miEAA_ICS_vs_UC_Analysis_results[[8]][i]), "; "))
  
  for (miRNA in no_lista_miRNAs){
      
      miRNAs = c(miRNAs, miRNA)
      
  }
  
  all_pathways = c(all_pathways, pathways)
  all_groups = c(all_groups, groups)
  output_list[[pathways]] = miRNAs
  available_miRNAs = c(available_miRNAs, miRNAs)
}

rep_pathways=c()

for (path in all_pathways){
  
  rep_pathways=c(rep_pathways, rep(path, length(output_list[[path]])))
}

require(visNetwork, quietly = TRUE)

nodes = data.frame(id = c(unique(available_miRNAs), unique(all_pathways)), 
                   group =c(rep("miRNAs", length(unique(available_miRNAs))), (all_groups)),
                   label = c(unique(available_miRNAs), unique(all_pathways)), 
                   value = c(rep(1, length(unique(available_miRNAs))), 
                             rep(2, length(unique(all_pathways)))))
# 
# miRNAs_with_targets=c("hsa-miR-548j-5p", "hsa-miR-554",
#                           "hsa-miR-876-3p", "hsa-miR-641",
#                           "hsa-miR-877-5p", "hsa-miR-509-3p",
#                           "hsa-miR-665", "hsa-miR-32-3p",
#                           "hsa-miR-31-5p", "hsa-miR-671-5p",
#                           "hsa-miR-33b-5p","hsa-let-7c-5p",
#                           "hsa-miR-10b-5p", "hsa-miR-139-3p",
#                           "hsa-miR-627-5p", "hsa-miR-30d-5p",
#                           "hsa-miR-185-5p", "hsa-miR-326")
# 
# for (i in 1:length(nodes$id)){
#   
#   if (length(which(miRNAs_with_targets == nodes$id[i])) > 0){
#     nodes$group[i] = "miRNAs_target"
#   }
# }

edges = data.frame(from = available_miRNAs, to = rep_pathways)

colnames(edges) = c("Source", "Target")

write.table(edges,
            file = ".\\edges.csv", 
            sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(nodes,
            file = ".\\nodes.csv", 
            sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
