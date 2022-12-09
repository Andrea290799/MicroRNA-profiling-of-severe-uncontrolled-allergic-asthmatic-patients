# MicroRNA profiling of severe uncontrolled allergic asthmatic patients
Repository containing all code used for the development of this work

## _miRNAs_normalization.R_

This script is used to perform RT-qPCR expression values normalization and developes the next functions:
- data_collection: gets the data from all the input files and turn some of its columns into manageable lists.
- undetermined_IPC_delete: Function that deletes the data of those plates that have 2 or 3 of its IPCs as undetermined.
- cts_transformation_40_NA: Function that transforms some Ct values following criteria in image below.
![image](https://user-images.githubusercontent.com/67425702/206718923-897b88ec-8af9-494b-b7e6-2b74aa3bb055.png)




** You have to take into account the neccessary changes in paths, file names and variables for the script to properly work
