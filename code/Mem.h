#ifndef MEM_H_INCLUDED
#define MEM_H_INCLUDED
#define POP_COUNT(BITS,i) BITS[(int)(i&0xFFFF)] + BITS[(int)((i>>16)&0xFFFF)] + BITS[(int)((i>>32)&0xFFFF)] + BITS[(int)((i>>48)&0xFFFF)]
#include <stdio.h>
#include <map>
#include <limits>
#include "PI.h"
/*
    Mem class is used to save the GWAS data in memory by a boolean representation.
    It provides some functions which are calculated base on the boolean representation.
*/
class Mem{
public:
    //initialize the memory based on the path of GWAS file.
    Mem(PI * pi);
    ~Mem();
    //Not be used in sheib now.
    void exchange(int a,int b);
    //For debug.
    void printTable(const int * indexes,const int & l);
    //If the object of Mem is created successfully.
    bool ready=false;
    //Calculate k2 for an SNP combination.
    double computeK2(const int * indexes,const int & l);
    //Calculate p-value of G-test for an SNP combination.
    double computeG(const int * indexes,const int & l);
    //Try to drop an SNP from the combination based on k2.
    void cutK2(int *& indexes,int & l,double & k2);
    //Not used in sheib now. It is a failed version of cutK2.
    void cutK22(int * & indexes, int & l,const double & k2);
    //Randomly generate an SNP combination.
    int ranGen();
    //Names of SNPs in the GWAS data.
    char ** names=NULL;
    //For control system. Print k2 and g for an SNP combination.
    void testAnIndividual(int * & indexes,int & l);
    //Calculate p-value of G-test for each one SNP.
    void fillG1();
    //return a map of (SNP name,SNP index)
    std::map<std::string,int> reversedNames();
private:
    PI * pi=NULL;
    //p-values of each one SNP.
    double * g1=NULL;
    //Buffer size used in reading the file.
    static const int BUFF_SIZE=40000;
    //input file type.
    int type=0;
    //Mean of k2 value of SNP combinations generated.
    double meanK2;
    //The number of SNP combinations generated.
    int nofc=0;
    //Used to calculate the contingency tables.
    int * BITS=NULL;
    //The path of input GWAS file.
    char * filename=NULL;
    //The path of input GWAS tped file.
    char * filenamePed=NULL;
    //The path of input GWAS tfam file.
    char * filenameFam=NULL;
    //Temporary variable used in filling the memory.
    char * chrs=NULL;
    //Temporary variable used in filling the memory.
    char * genos=NULL;
    //Temporary variable used in filling the memory.
    long * poss=NULL;
    //Bit size of the unsigned long long.
    static const int SIZE=sizeof(unsigned long long)*8;
    //Max name length of SNP name.
    static const int SIZE_NAME=40;
    //Temporary variable used in calculation based on the boolean representation.
    static const unsigned long long ONE=0x8000000000000000;
    //The memories used to save the GWAS data.
    unsigned long long * mem1=NULL;
    unsigned long long * mem0=NULL;
    //The number of samples, the number of cases, the number of controls, the number of SNPs, the length of mem1, and the length of mem0.
    int m,m1,m0,n,l1,l0;
    //Input File pointer.
    FILE * fp=NULL;
    //Called in calculating the contingency tables based on the boolean representation.
    int bitCount(unsigned long long i);
    //get the number of SNPs and samples, type=0.
    bool fillMN();
    //Temporary variable used in filling the memory. It is used to save which sample is case or control.
    std::map<int,bool> mSamples;
    //get the number of samples, type=1.
    bool fillM();
    //get the number of SNPs, type=1.
    bool fillN();
    //Fill mem0 and mem1, type=0.
    bool fillMem0();
    //Fill mem0 and mem1, type=1.
    bool fillMem1();
    //Called in calculating the contingency tables based on the boolean representation.
    int popCount(const unsigned long long * a,const int & l);
    //Get the contingency table for an SNP combination based on the boolean representation.
    int * getTable(const int * indexes,const int & l,int & lg);
    //For debug.
    void testMem();
    //For debug.
    void testTable();
    //Calculate max order based on the samples of the GWAS data.
    int recommendOrder();
    //Called by control system. It will print k2 and g for an SNP combination.
    void testK2G(const int * indexes, const int & l);
};



#endif // MEM_H_INCLUDED
