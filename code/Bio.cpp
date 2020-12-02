#include "Bio.h"
#include "PI.h"
#include "Mem.h"
#include <string.h>
#include <string>
#include <list>
/*
    The cpp file of class Bio.
*/
Bio::Bio(PI * pi){
    if(pi->filenameSNP2Genes==NULL){
        if(pi->nShow>0){
            printf("there is no snp2genes map file provided\n");
        }
    }
    else{
        this->pi=pi;
        this->snp2Genes=(std::vector<int> **)malloc(sizeof(std::vector<int *>*)*pi->n);
        for(int i=0;i<pi->n;i++){
            snp2Genes[i]=NULL;
        }
        if(!loadSNP2Genes()){
            printf("Bio::loadSNP2Genes error\n");
            return;
        }
        else{
            if(pi->filenameAssociatedGenes==NULL){
                if(pi->nShow>0){
                    printf("there is no associated genes file provided\n");
                }
            }
            else{
                if(!loadAssociatedGenes()){
                    printf("Bio::loadAssociatedGenes error\n");
                    return;
                }
            }
        }
        pi->bio=this;
    }
}
Bio::~Bio(){
    if(snp2Genes!=NULL){
        for(int i=0;i<pi->n;i++){
            if(snp2Genes[i]!=NULL){
                delete snp2Genes[i];
            }
        }
        free(snp2Genes);
        snp2Genes=NULL;
    }
    if(gene2Snps!=NULL){
        for(int i=0;i<nGenes;i++){
            if(gene2Snps[i]!=NULL){
                delete gene2Snps[i];
            }
        }
        gene2Snps=NULL;
    }
    if(associatedGenes!=NULL){
        for(int i=0;i<nGenes;i++){
            if(associatedGenes[i]!=NULL){
                delete associatedGenes[i];
            }
        }
        associatedGenes=NULL;
    }
}
bool Bio::loadSNP2Genes(){
    std::map<std::string,int> maps=((Mem *)(pi->mem))->reversedNames();
    //¶ÁÈ¡snp2genesÎÄ¼þ
    FILE * fp=fopen(pi->filenameSNP2Genes,"r");
    if(fp==NULL){
        maps.clear();
        return false;
    }
    char name[NAME_SIZE];
    char c;
    int iName=0;
    int snpIndex=-1;
    int iGene=0;
    std::map<std::string,int>::iterator itsi;
    while((c=fgetc(fp))!=EOF){
        if(c=='\t'&&snpIndex==-1){
            name[iName]='\0';
            std::string s(name);
            itsi=maps.find(s);
            if(itsi==maps.end()){
                do{
                    c=fgetc(fp);
                }
                while(c!='\n'&&c!=EOF);
                iName=0;
                continue;
            }
            else{
                snpIndex=itsi->second;
                iName=0;
            }
        }
        else if(c=='\t'&&snpIndex!=-1&&iName>0){
            name[iName]='\0';
            std::string s(name);
            itsi=gene2Index.find(s);
            if(itsi==gene2Index.end()){
                gene2Index.insert(std::pair<std::string,int>(s,iGene));
                if(snp2Genes[snpIndex]==NULL){
                    snp2Genes[snpIndex]=new std::vector<int>();
                }
                snp2Genes[snpIndex]->push_back(iGene);
                genes.push_back(s);
                iGene++;
            }
            else{
                if(snp2Genes[snpIndex]==NULL){
                    snp2Genes[snpIndex]=new std::vector<int>();
                }
                snp2Genes[snpIndex]->push_back(itsi->second);
            }
            iName=0;
        }
        else if(c=='\n'&&snpIndex!=-1){
            if(iName>0){
                name[iName]='\0';
                std::string s(name);
                itsi=gene2Index.find(s);
                if(itsi==gene2Index.end()){
                    gene2Index.insert(std::pair<std::string,int>(s,iGene));
                    if(snp2Genes[snpIndex]==NULL){
                        snp2Genes[snpIndex]=new std::vector<int>();
                    }
                    snp2Genes[snpIndex]->push_back(iGene);
                    genes.push_back(s);
                    iGene++;
                }
                else{
                    if(snp2Genes[snpIndex]==NULL){
                        snp2Genes[snpIndex]=new std::vector<int>();
                    }
                    snp2Genes[snpIndex]->push_back(itsi->second);
                }
                iName=0;
            }
            snpIndex=-1;
        }
        else if(c=='\r'){
            continue;
        }
        else{
            name[iName]=c;
            iName++;
        }
    }
    fclose(fp);
    maps.clear();
    for(int i=0;i<pi->n;i++){
        if(snp2Genes[i]!=NULL){
            snp2Genes[i]->resize(snp2Genes[i]->size());
        }
    }
    nGenes=(int)gene2Index.size();
    gene2Snps=(std::vector<int> **)malloc(sizeof(std::vector<int>*)*nGenes);
    for(int i=0;i<nGenes;i++){
        gene2Snps[i]=NULL;
    }
    for(int i=0;i<pi->n;i++){
        if(snp2Genes[i]!=NULL){
            int tt=snp2Genes[i]->size();
            for(int j=0;j<tt;j++){
                if(gene2Snps[snp2Genes[i]->at(j)]==NULL){
                    gene2Snps[snp2Genes[i]->at(j)]=new std::vector<int>();
                }
                gene2Snps[snp2Genes[i]->at(j)]->push_back(i);
            }
        }
    }
    for(int i=0;i<nGenes;i++){
        if(gene2Snps[i]!=NULL){
            gene2Snps[i]->resize(gene2Snps[i]->size());
        }
    }
    //print();
    return true;
}
bool Bio::loadAssociatedGenes(){
    associatedGenes=(std::vector<int> **)malloc(sizeof(std::vector<int> *)*nGenes);
    for(int i=0;i<nGenes;i++){
        associatedGenes[i]=NULL;
    }
    FILE * fp=fopen(pi->filenameAssociatedGenes,"r");
    int snpIndex=-1;
    int iName=0;
    char c;
    char name[NAME_SIZE];
    std::map<std::string,int>::iterator itsi;
    while((c=fgetc(fp))!=EOF){
        if(c=='\t'&&snpIndex==-1){
            name[iName]='\0';
            std::string s(name);
            itsi=gene2Index.find(s);
            if(itsi==gene2Index.end()){
                do{
                    c=fgetc(fp);
                }
                while(c!='\n'&&c!=EOF);
                iName=0;
                continue;
            }
            else{
                snpIndex=itsi->second;
                iName=0;
            }
        }
        else if(c=='\t'&&snpIndex!=-1&&iName>0){
            name[iName]='\0';
            std::string s(name);
            itsi=gene2Index.find(s);
            if(itsi!=gene2Index.end()){
                if(associatedGenes[snpIndex]==NULL){
                    associatedGenes[snpIndex]=new std::vector<int>();
                }
                associatedGenes[snpIndex]->push_back(itsi->second);
            }
            iName=0;
        }
        else if(c=='\r'){
            continue;
        }
        else if(c=='\n'&&snpIndex!=-1){
            if(iName>0){
                name[iName]='\0';
                std::string s(name);
                itsi=gene2Index.find(s);
                if(itsi!=gene2Index.end()){
                    if(associatedGenes[snpIndex]==NULL){
                        associatedGenes[snpIndex]=new std::vector<int>();
                    }
                    associatedGenes[snpIndex]->push_back(itsi->second);
                }
            }
            iName=0;
            snpIndex=-1;
        }
        else{
            name[iName]=c;
            iName++;
        }
    }
    fclose(fp);
    for(int i=0;i<nGenes;i++){
        if(associatedGenes[i]!=NULL){
            associatedGenes[i]->resize(associatedGenes[i]->size());
        }
    }

    return true;
}
//for debug
void Bio::print(){
    printf("genes\n");
    for(int i=0;i<(int)genes.size();i++){
        printf("%s\n",genes[i].c_str());
    }
    for(int i=0;i<pi->n;i++){
        if(snp2Genes[i]!=NULL){
            printf("%d",i);
            for(int j=0;j<(int)snp2Genes[i]->size();j++){
                printf("\t%s",genes[snp2Genes[i]->at(j)].c_str());
            }
            printf("\n");
        }
    }
    for(int i=0;i<pi->n;i++){
        if(snp2Genes[i]!=NULL){
            printf("%d",i);
            for(int j=0;j<(int)snp2Genes[i]->size();j++){
                printf("\t%d",snp2Genes[i]->at(j));
            }
            printf("\n");
        }
    }
    for(int i=0;i<nGenes;i++){
        if(associatedGenes[i]!=NULL){
            printf("%s",genes[i].c_str());
            for(int j=0;j<(int)associatedGenes[i]->size();j++){
                printf("\t%s",genes[associatedGenes[i]->at(j)].c_str());
            }
            printf("\n");
        }
    }
    for(int i=0;i<nGenes;i++){
        if(gene2Snps[i]!=NULL){
            printf("%s",genes[i].c_str());
            for(int j=0;j<(int)gene2Snps[i]->size();j++){
                printf("\t%d",gene2Snps[i]->at(j));
            }
            printf("\n");
        }
    }
}
