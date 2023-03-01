# Additional code

In this folder there are some scripts used for graphical representation. 

## _Heatmap.R_
This script was used to perform a hierarchical clustering (distance: Manhatan, method: complete), by depicting a heatmap. Used data were expression values of the 40 differentially expressed miRNAs in each group of patients (normalized by Z-score). 

![image](https://user-images.githubusercontent.com/67425702/222125146-7829901e-bcd5-48cc-b6b4-a1a5c0787584.png)

## _Boxplots.R_

This script prepares a 3X4 boxplot figure of the expression values of all those 12 top miRNAs ranked by abs(log2(FC)), comparing mild and severe uncontrolled groups. 

![image](https://user-images.githubusercontent.com/67425702/222125263-ba79c6e1-c5a5-4eaa-9e92-b09bc707099b.png)


## _Data_preparation_to_Gephi.R_
This script prepares the miEAA output file following the format required by Gephi to obtain links between diferentially expressed miRNAs and its associated biological functions. 

![image](https://user-images.githubusercontent.com/67425702/222125436-0c120d3e-3765-4b1e-a4d3-09fb262b2b0e.png)


## _Correlations.R_

This script calculates correlations between differentially expresssed miRNAs in mild and severe uncontrolled patients and inflamatory related metabolites. First, normality of both populations is assessed. Pearson correlation test is used when normality is met; otherwise, Spearman correlation test is applied. With r parameter and p-value, a correlation plot is obtained. 

![image](https://user-images.githubusercontent.com/67425702/222125502-d176b2d2-7efd-404e-a974-908a5f8a8656.png)

