# sheib
# introduction

This web page provides source code and simulated datasets mentioned in the article "SHEIB: a novel stochastic-based approach for detecting high-order epistatic interactions using bioinformation in genome-wide association studies". The article has been submitted to BMC Bioinformatics but has not been published. Please don't download these files, if you are not reviewers or editors.
The real dataset (WTCCC dataset) in the article will not be uploaded. Because the dataset was downloaded from https://www.wtccc.org.uk/. We have no right to upload it.

# source code

SHEIB is implemented by C++ and compliled by Cygwin64. The source code has been uploaded. We have also uploaded two executable files (sheib was compiled on linux 64 and sheib.exe was compiled on win7 64).

# introduction to the uploaded files

# run on linux 64

compile: g++ -std=c++11 -lpthread *.cpp -o sheib

./sheib -type 0 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in data.txt -out result.txt

./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt

./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt

./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt -AssociatedGenes hrpd_after_3.txt

# run on win7 64

The folder uploaded is a project of CodeBlocks. You can open the sheib.cbp using CodeBlocks.

sheib -type 0 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in data.txt -out result.txt

sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt

sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt

sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt -AssociatedGenes hrpd_after_3.txt

# simulated dataset
The simulated data is too large to upload (3.08GB simulated_data.tar.gz). So we only uploaded the DNME3_100 dataset. If you need all the three simulated datasets, you can send email to 2276632042@qq.com (the first author of the article).
