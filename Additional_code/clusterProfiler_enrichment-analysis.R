library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(dplyr)
library(tidyr)
library(ReactomePA)
library(DOSE) 
library(stats)
library(msigdbr)
library(limma)
library(DAPAR)
library(GOSemSim)
library(enrichplot)
library(ggplot2)

# database of miRNA symbols
edb <- org.Hs.eg.db
keys <- keys(edb, keytype="SYMBOL")

db_sel <- AnnotationDbi::select(edb, keys=keys, columns=c("ENTREZID"),
                                keytype="SYMBOL")

# reading and transforming data 
# output of miRNAs_diff_exp
df = read.table("ANOVA_KW_differences_GME-DE.txt", sep = "\t", header = T) %>%
  filter(UC_vs_ICS_pvalue < 0.05) %>%
  mutate(miRNA = toupper(gsub("-", "", gsub("hsa", "", miRNA)))) %>%
  mutate(miRNA = gsub("*5P", "", miRNA)) %>%
  mutate(miRNA = gsub("*3P", "", miRNA)) %>% 
  mutate(updated_symbols = alias2SymbolTable(miRNA, species="Hs")) %>% 
  left_join(., db_sel, by = c('updated_symbols'='SYMBOL'))


# MIR126 is duplicated because initial names contains 5p and 3p 
df <- df[!duplicated(df$miRNA),]


original_gene_list <- as.numeric(df$UC_vs_ICS_log2.FC.)

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

genes <- sort(gene_list, decreasing = TRUE)

set.seed(333)
univ_h <- univ_AnnotDbPkg("org.Hs.eg.db")

set.seed(333)
egoBP <- enrichGO(gene          = df$ENTREZID,
                  universe      = univ_h,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

# cnetplot
pathways <- data.frame(egoBP) %>%
  filter(Count >= 10) %>% 
  arrange(Count)

sel_pathways <- pathways$Description

pdf("DE_miRNAs-egoBP-cnetplot-count10.pdf", 15, 10)
cnetplot(egoBP, categorySize="pvalue", showCategory = sel_pathways, 
         foldChange=genes, node_label_size = 6)

dev.off()

# treeplot

egoBP_pt <- pairwise_termsim(egoBP)

pdf("DE_miRNAs-egoBP-treeplotp_centroid.pdf", 20, 15)
treeplot(egoBP_pt, showCategory = 50, cluster.params = list(n = 4, method = "centroid"), 
         fontsize = 6, offset.params = list(tiplab = 0.6))
dev.off()

