# sheib
# introduction
This web page provides source code and simulated datasets mentioned in the article "SHEIB: a novel stochastic-based approach for detecting high-order epistatic interactions using bioinformation in genome-wide association studies". The article has been submitted to BMC Bioinformatics but has not been published. Please don't download these files, if you are not reviewers or editors.<br>
The real dataset (WTCCC dataset) in the article will not be uploaded. Because the dataset was downloaded from https://www.wtccc.org.uk/. We have no right to upload it.
# source code
SHEIB is implemented by C++ and compliled by Cygwin64. The source code has been uploaded. We have also uploaded two executable files (sheib was compiled on linux 64 and sheib.exe was compiled on win7 64).
# introduction to the uploaded files
## sheib.zip:
All files like ***.cpp** and ***.h** are source code of SHEIB.<br>
**bd_gwas.tped** and **bd_gwas.tfam** are sample examples of real gwas dataset. There are only 4000 SNPs in the example.<br>
data.txt is a sample example of simulated data. It contains 1000 SNPs and [P1,P2] is the epistatic interaction.<br>
snps_after_2.txt is the gene-mapping data from NCBI.<br>
gene_pairs_after_4.txt is the gene association data from DIP.<br>
hrpd_after_3.txt is the gene association data from HRPD.<br>
mint_after_4.txt is the gene association data from MINT.<br>
reactome_after_4.txt is the gene association data from Reactome.<br>
sheib is an executable file compiled by gcc on Linux 64.<br>
sheib.exe is an executable file compiled by cygwin64 on win7 64.<br>
sheib.cbp, sheib.depend and sheib.layout are files generated by CodeBlocks automatically. Users can use CodeBlocks to open them.<br>
## simulated_datasets
### DME_and_DNME_100
The folder contains 68 2-order epistatic models. For each model, a compressed file is uploaded. In each compressed file, there are 100 simulated data (files). The filename indicates the epistatic interaction in the file.
### DME_and_DNME_1000
This simulated dataset is too large to be uploaded. If you need it, you can send email to 2276632042@qq.com (the first author of the article).<br>
### DNME3_100
The folder contains 40 3-order epistatic models. For each model, a compressed file is uploaded. In each compressed file, there are 100 simulated data (files). The filename indicates the epistatic interaction in the file.
# run on linux 64
compile: g++ -std=c++11 -lpthread *.cpp -o sheib<br>
./sheib -type 0 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in data.txt -out result.txt<br>
./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt<br>
./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt<br>
./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt -AssociatedGenes hrpd_after_3.txt<br>
# run on win7 64
The folder uploaded is a project of CodeBlocks. You can open the sheib.cbp using CodeBlocks.<br>
sheib -type 0 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in data.txt -out result.txt<br>
sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt<br>
sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt<br>
sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt -AssociatedGenes hrpd_after_3.txt<br>
# parameters
-type: default 0, type of the input file. 0 (simulated data), 1 (real gwas data).<br>
-cG: default: 0.05, threshold.<br>
-cGc: default 0.05, threshold.<br>
-o: default -1, maximum order while generating random SNP combinations. -1 means that SHEIB should calculate it based on the number of samples in the gwas data.<br>
-maxGen: default -1, maximum iterations. -1 means that SHEIB will never stop by itself.<br>
-pb: default 0.8, the probability of considering bioinformation while generating random SNP combinations.<br>
-nShow: default 4, it is used to control echo. SHEIB will print epistatic interactions every nShown iteractions. nShown<0 means that SHEIB will not print interactions.<br>
-seed: default 0, random seed.<br>
-rn: default -1, how many epistatic interactions should be collected. -1 means that SHEIB will write all interactions into the result.<br>
-cs: default 0. 1 means that SHEIB will start a control thread. When users type "Enter", the program will be paused and enter an interactive interface. Users can type some commands. type "help" to get a command list and further information. 0 means SHEIB will not start the control thread.<br>
-in: default data.txt, the input gwas file.<br>
-out: default result.txt, the output file.<br>
-SNP2Genes: default null, gene-mapping data.<br>
-AssociatedGenes: default null, gene association data.<br>

