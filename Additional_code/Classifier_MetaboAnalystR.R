library(MetaboAnalystR)
library(writexl)

mSet<-InitDataObjects("conc", "roc", FALSE)
# table with miRNA expression values, patients in rows and miRNAs in columns
mSet<-Read.TextData(mSet, "40_DE_miRNAs.csv", "rowu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-SanityCheckData(mSet)
mSet<-ImputeMissingVar(mSet, method="min")
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

write_xlsx(as.data.frame(mSet[["dataSet"]]$norm), "imputed_data_by_0.2LODvalue.xlsx")

mSet<-PlotNormSummary(mSet, "norm_0_", "pdf", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "pdf", 72, width=NA)

set.seed(333)

mSet<-SetAnalysisMode(mSet, "explore")
mSet<-PrepareROCData(mSet)
mSet<-PerformCV.explore(mSet, "rf", "rf", 2)
mSet<-PlotProbView(mSet, "cls_prob_0_", "pdf", 72, -1, 0, 0)
mSet<-PlotImpBiomarkers(mSet, "cls_imp_0_", "pdf", 72, -1, "freq", 40);
mSet<-PlotAccuracy(mSet, "cls_accu_0_", "pdf", 72)
mSet<-PlotROC(mSet, "cls_roc_0_", "pdf", 72, 0, "threshold", 0, 0, "fpr", 0.5)


.get.mSet <- function(mSetObj=NA){

    return(mSetObj);

}


GetBestModelIndex <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$analSet$multiROC$best.model.inx;
}

mSet<-GetBestModelIndex(mSet)

