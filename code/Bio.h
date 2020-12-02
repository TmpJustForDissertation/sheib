#ifndef BIO_H_INCLUDED
#define BIO_H_INCLUDED
#include "PI.h"
#include <vector>
#include <map>
#include <set>
/*
    header file of class Bio.
    Class Bio is used to save the bioinformation in the memory.
*/
struct Snp{
    std::set<int> * snps=NULL;
    //int * snps=NULL;
    int l=0;
};
class Bio{
public:
    //The number of genes provided in gene-mapping data.
    int nGenes;
    //Initialize the object based on the path provided by users and save the object to pi->bio.
    Bio(PI * pi);
    ~Bio();
    //Gene names provided in gene-mapping data.
    std::vector<std::string> genes;
    //The map converting a gene name to its index[0,nGenes).
    std::map<std::string,int> gene2Index;
    //The map converting SNP index to Gene indexes of genes which the SNP is located in.
    std::vector<int> ** snp2Genes=NULL;
    //The map converting gene index to SNP indexes of SNPs which are located in the gene.
    std::vector<int> ** gene2Snps=NULL;
    //The map converting gene index to gene indexes of genes which are associated with the gene in gene association data.
    std::vector<int> ** associatedGenes=NULL;
    //struct Snp ** snps=NULL;
    //for debug
    void print();
private:
    //max gene name length.
    static const int NAME_SIZE=40;//409600+1;
    PI * pi=NULL;
    //load gene-mapping data into memory.
    bool loadSNP2Genes();
    //load gene association data into memory.
    bool loadAssociatedGenes();
};


#endif // BIO_H_INCLUDED
