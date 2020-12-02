# sheib
# 1.Introduction
This web page provides source code and simulated datasets mentioned in the article "SHEIB: a novel stochastic-based approach for detecting high-order epistatic interactions using bioinformation in genome-wide association studies". The article has been submitted to BMC Bioinformatics but has not been published. Please don't download these files, if you are not reviewers or editors.<br>
The real dataset (WTCCC dataset) in the article will not be uploaded. Because the dataset was downloaded from https://www.wtccc.org.uk/. We have no right to upload it.
# 2.The uploaded files
**/code/README.md** provides some descriptions of this capsule.<br>
**/code/run.sh** is the script to compile and run sheib on this platform.<br>
**/code/sheib** is an executable file compiled by gcc on Linux 64.<br>
All files like ***.cpp** and ***.h** are source code of SHEIB. It is implemented by C++.<br>
**/data/example.tped** and **/data/example.tfam** are sample examples of real gwas dataset. There are only 4000 SNPs in the example.<br>
**/data/example.txt** is a sample example of simulated data. It contains 1000 SNPs and [P1,P2] is the epistatic interaction.<br>
**/data/bio/snps_after_2.txt** is the gene-mapping data from NCBI.<br>
**/data/bio/gene_pairs_after_4.txt** is the gene association data from DIP.<br>
**/data/bio/hrpd_after_3.txt** is the gene association data from HRPD.<br>
**/data/bio/mint_after_4.txt** is the gene association data from MINT.<br>
**/data/bio/reactome_after_4.txt** is the gene association data from Reactome.<br>
# 3.Execution
compile: g++ -std=c++11 *.cpp -o sheib -lpthread<br>
./sheib -type 0 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in data.txt -out result.txt<br>
./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt<br>
./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt<br>
./sheib -type 1 -cG 0.05 -cGc 0.05 -o -1 -maxGen 4000000 -pb 0.8 -nShow 4000 -seed 0 -rn -1 -cs 0 -in bd_gwas -out result.txt -SNP2Genes snps_after_2.txt -AssociatedGenes hrpd_after_3.txt<br>
# 4.parameters

parameter|default|description
----|----|----
type|0|Type of the input file. 0 (simulated data), 1 (real gwas data).
cG|0.05|Threshold of pvalue of G-test.
cGc|0.05|Threshold of gc.
o|-1|Maximum order while generating random SNP combinations. -1 means that SHEIB should calculate it based on the number of samples in the gwas data.
maxGen|-1|Maximum iterations. -1 means that SHEIB will never stop by itself.
pb|0.8|The probability of considering bioinformation while generating random SNP combinations.
nShow|4|It is used to control echo. SHEIB will print epistatic interactions every nShown iteractions. nShown<0 means that SHEIB will not print interactions.
seed|0|Random seed.
rn|-1|How many epistatic interactions should be collected. -1 means that SHEIB will write all interactions into the result.
cs|0|1 means that SHEIB will start a control thread. When users type "Enter", the program will be paused and enter an interactive interface. Users can type some commands. type "help" to get a command list and further information. 0 means SHEIB will not start the control thread.
in|data.txt|The input gwas file.
out|result.txt|The output file.
SNP2Genes|null|Gene-mapping data.
AssociatedGenes|null|Gene association data.

