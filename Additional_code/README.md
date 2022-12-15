## _Heatmap.R_
This script was used to perform a hierarchical clustering (distance: Manhatan, method: complete), by depicting a heatmap. Used data were expression values of the 40 differentially expressed miRNAs in each group of patients (normalized by Z-score). 

![image](https://user-images.githubusercontent.com/67425702/207881560-7a6662c1-dce1-48b0-aef8-9496648eac17.png)


## _Boxplots.R_

This script prepares a 6x3 boxplot figure of the expression values of all those miRNAs that have known targets, comparing mild and severe uncontrolled groups. 

![image](https://user-images.githubusercontent.com/67425702/207881935-ba544940-3dcd-40ba-81dc-4d840167ab17.png)


## _Data_preparation_to_Gephi.R_
This script prepares the miEAA output file following the format required by Gephi to obtain links between diferentially expressed miRNAs with known targets and its associated biological functions. 

![image](https://user-images.githubusercontent.com/67425702/207881646-2a522330-f72b-4bb4-9972-1254850da9e9.png)


## _Correlations.R_

This script calculates correlations between differentially expresssed miRNAs in mild and severe uncontrolled patients and inflamatory related metabolites. First, normality of both populations is assessed. Pearson correlation test is used when normality is met; otherwise, Spearman correlation test is applied. With r parameter and p-value, a correlation plot is obtained. 

![image](https://user-images.githubusercontent.com/67425702/207881466-54d8e012-d128-4166-b4fb-ad39533b6a39.png)


