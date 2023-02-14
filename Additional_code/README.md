# Additional code

In this folder there are some scripts used for graphical representation. 

## _Heatmap.R_
This script was used to perform a hierarchical clustering (distance: Manhatan, method: complete), by depicting a heatmap. Used data were expression values of the 40 differentially expressed miRNAs in each group of patients (normalized by Z-score). 

![image](https://user-images.githubusercontent.com/67425702/218722322-6b8604cf-0528-4866-b21f-5869ce78f644.png)

## _Boxplots.R_

This script prepares a 4X4 boxplot figure of the expression values of all those 16 top miRNAs ranked by abs(log2(FC)), comparing mild and severe uncontrolled groups. 

![image](https://user-images.githubusercontent.com/67425702/218722655-e41a49db-1e0c-4927-b776-4f02eaa681ef.png)


## _Data_preparation_to_Gephi.R_
This script prepares the miEAA output file following the format required by Gephi to obtain links between diferentially expressed miRNAs and its associated biological functions. 

![image](https://user-images.githubusercontent.com/67425702/218723022-03de36b9-3f21-4c99-88e5-bfd9bacd7cdd.png)


## _Correlations.R_

This script calculates correlations between differentially expresssed miRNAs in mild and severe uncontrolled patients and inflamatory related metabolites. First, normality of both populations is assessed. Pearson correlation test is used when normality is met; otherwise, Spearman correlation test is applied. With r parameter and p-value, a correlation plot is obtained. 

![image](https://user-images.githubusercontent.com/67425702/218723147-1bb71515-7da5-4227-84f6-cb3696888b8c.png)

