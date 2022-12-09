# MicroRNA profiling of severe uncontrolled allergic asthmatic patients
Repository containing all code used for the development of this work

## _miRNAs_normalization.R_

This script is used to perform RT-qPCR expression values normalization and developes the next functions:
- data_collection: gets the data from all the input files and turn some of its columns into manageable lists.
- undetermined_IPC_delete: Function that deletes the data of those plates that have 2 or 3 of its IPCs as undetermined.
- cts_transformation_40_NA: Function that transforms some Ct values following criteria in image below.

![image](https://user-images.githubusercontent.com/67425702/206719120-02a46280-d95b-46b2-9951-87467cc38e1a.png)

- blank_NA_delete: Function that deletes those miRNAs whose cts are NA or whose miRNA names are Blank.
- interplate_normalization: Function that do interplate normalization on the data.
- GME: Function that applies GME normalization method.

![image](https://user-images.githubusercontent.com/67425702/206719972-6ff2d126-66bf-462b-afc9-a5065cca33df.png)

- ΔΔCt: Function that does ΔΔCt normalization.

![image](https://user-images.githubusercontent.com/67425702/206720501-ad3db4ab-830d-465e-9810-81d28cc9c039.png)

- DOS_ΔΔct: Function that calculates 2^-ΔΔCt.
- main: Function that executes the main code. 



## _miRNAs_diff_exp.R_

This script is used to perform differential expression analysis between groups of patients (control, mild and severe uncontrolled):
- data_collection: gets the data from all the input files and turn some of its columns into manageable lists.
- statistical_tests_pvalues: Function that does different statistical tests. 
- main: Function that executes the main code. It also calculates the Fold Changes of each miRNA in each subjects' comparative. 


## _Boxplots.R_

This script prepares a 6x3 boxplot figure of the expression values of all those miRNAs that have known targets, comparing mild and severe uncontrolled groups. 


## _Correlations.R_

This script calculates correlations between differentially expresssed miRNAs in mild and severe uncontrolled patients and inflamatory related metabolites. First, normality of both populations is assessed. Pearson correlation test is used when normality is met; otherwise, Spearman correlation test is applied. With r parameter and p-value, a correlation plot is obtained. 

## _Heatmap.R_
This script was used to perform a hierarchical clustering (distance: Manhatan, method: complete), by depicting a heatmap. Used data were expression values of the 40 differentially expressed miRNAs in each group of patients (normalized by Z-score). 

## _Data_preparation_to_Gephi.R_
This script prepares the miEAA output file following the format required by Gephi to obtain links between diferentially expressed miRNAs with known targets and its associated biological functions. 

![image](https://user-images.githubusercontent.com/67425702/206741023-e6f04b70-05d2-4b68-bdd8-c0699455836c.png)


~~~
You have to take into account the neccessary changes in paths, file names and variables for the script to properly work
~~~
