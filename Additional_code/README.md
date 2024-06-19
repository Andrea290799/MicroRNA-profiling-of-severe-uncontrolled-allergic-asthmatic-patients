# Additional code

In this folder there are some scripts used for dowmstream analyses and graphical representation. 

## _Heatmap.R_
This script was used to perform a hierarchical clustering (distance: Manhatan, method: complete), by depicting a heatmap. Used data were expression values of the 40 differentially expressed miRNAs in each group of patients (normalized by Z-score). 

![image](https://user-images.githubusercontent.com/67425702/222125146-7829901e-bcd5-48cc-b6b4-a1a5c0787584.png)

## _Boxplots.R_

This script prepares a 3X3 boxplot figure of the expression values of all those 9 top miRNAs ranked by abs(log2(FC)), comparing mild and severe uncontrolled groups. 

![image](https://github.com/Andrea290799/MicroRNA-profiling-of-severe-uncontrolled-allergic-asthmatic-patients/assets/67425702/a4006a5d-9d2d-4154-8777-e69c73a88139)


## _clusterProfiler_enrichment-analysis.R_

This script performs functional enrichment analysis (over representation analysis, ORA) of the 40 differentially expressed (DE) miRNAs between the severe uncontrolled and mild patients’ groups by using clusterProfiler R package.

![image](https://github.com/Andrea290799/MicroRNA-profiling-of-severe-uncontrolled-allergic-asthmatic-patients/assets/67425702/f0366e67-4fc5-4198-b02a-90c55fb372dc)

![image](https://github.com/Andrea290799/MicroRNA-profiling-of-severe-uncontrolled-allergic-asthmatic-patients/assets/67425702/39bab281-dd80-41fe-8e02-4fe503b03577)


## _Correlations.R_

This script calculates correlations between differentially expresssed miRNAs in mild and severe uncontrolled patients and inflamatory related metabolites. First, normality of both populations is assessed. Pearson correlation test is used when normality is met; otherwise, Spearman correlation test is applied. With r parameter and p-value, a correlation plot is obtained. 

![image](https://user-images.githubusercontent.com/67425702/222125502-d176b2d2-7efd-404e-a974-908a5f8a8656.png)

## _Chordplot.R_

This scripts generates a circos plot of the regulation by the 40 DE miRNAs of their predicted targets, restricting the search to those targets related to correlated metabolites for each DE miRNA. 

![image](https://github.com/Andrea290799/MicroRNA-profiling-of-severe-uncontrolled-allergic-asthmatic-patients/assets/67425702/faed1f9c-4adb-4f34-bf5e-f3a5034dad67)

## _Classifier_MetaboAnalystR.R_

This script uses MetaboAnalystR R package for classifier modelling and performance evaluation trough ROC curves. 

![image](https://github.com/Andrea290799/MicroRNA-profiling-of-severe-uncontrolled-allergic-asthmatic-patients/assets/67425702/12e22db3-54a6-46fb-b09a-ab1ec1da5f13)

![image](https://github.com/Andrea290799/MicroRNA-profiling-of-severe-uncontrolled-allergic-asthmatic-patients/assets/67425702/a0674ff0-6d3f-4fd8-bd23-10279786b899)

![image](https://github.com/Andrea290799/MicroRNA-profiling-of-severe-uncontrolled-allergic-asthmatic-patients/assets/67425702/fbb944aa-70f8-457f-9d42-140cafc90383)

![image](https://github.com/Andrea290799/MicroRNA-profiling-of-severe-uncontrolled-allergic-asthmatic-patients/assets/67425702/7831eef3-0666-422a-81c1-848f8b03b8e1)


