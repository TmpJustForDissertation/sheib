#include "PI.h"
#include "Mem.h"
#include "functions.h"
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
/*
    The implemented source code of PI.h.
*/
bool PI::fillParameters(int argc,char * argv[]){
    /*
    -type 0;
    -cG 0.01;
    -cGc 0.05;
    -o -1;
    -maxGen -1;
    -pe 0.8
    -nShow 4;
    -seed 1
    -rn -1
    -cs 0;
    -SNP2Genes NULL
    -AssociatedGenes NULL
    -in data.txt;
    -out result.txt;
    ./sseh -type 1 -cG 0.01 -cGc 0.01 -o -1 -maxGen 40000 -pe 0.8 -nShow 4 -seed 1 -rn -1 -cs 0 -in  /home/sly/data/WTCCC20181025/backup5/bd_gwas -out result_bd.txt
    */
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i],"-sim")==0){
            i++;
            if(i<argc)
                this->sim=atoi(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-type")==0){
            i++;
            if(i<argc)
                this->type=atoi(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-cG")==0){
            i++;
            if(i<argc)
                this->cG=atof(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-cGc")==0){
            i++;
            if(i<argc)
                this->cGc=atof(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-o")==0){
            i++;
            if(i<argc)
                this->o=atoi(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-maxGen")==0){
            i++;
            if(i<argc)
                this->maxGen=atoi(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-pb")==0){
            i++;
            if(i<argc)
                this->pb=atof(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-nShow")==0){
            i++;
            if(i<argc)
                this->nShow=atoi(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-seed")==0){
            i++;
            if(i<argc)
                this->seed=atoi(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-rn")==0){
            i++;
            if(i<argc)
                this->rn=atoi(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-cs")==0){
            i++;
            if(i<argc)
                this->cs=atoi(argv[i]);
            else
                return false;
        }
        else if(strcmp(argv[i],"-in")==0){
            i++;
            if(i<argc)
                this->filename=argv[i];
            else
                return false;
        }
        else if(strcmp(argv[i],"-out")==0){
            i++;
            if(i<argc)
                this->filenameO=argv[i];
            else
                return false;
        }
        else if(strcmp(argv[i],"-SNP2Genes")==0){
            i++;
            if(i<argc)
                this->filenameSNP2Genes=argv[i];
            else
                return false;
        }
        else if(strcmp(argv[i],"-AssociatedGenes")==0){
            i++;
            if(i<argc)
                this->filenameAssociatedGenes=argv[i];
            else
                return false;
        }
        else{
            return false;
        }
    }
    if(this->filename==NULL){
        this->filename=(char *)"data.txt";
    }
    if(this->filenameO==NULL){
        this->filenameO=(char *)"result.txt";
    }
    return true;
}
void PI::print(){
    printf("parameters:\n");
    printf("%18s  %d\n","type",this->type);
    printf("%18s  %f\n","cG",this->cG);
    printf("%18s  %f\n","cGc",this->cGc);
    printf("%18s  %d\n","o",this->o);
    printf("%18s  %f\n","pb",this->pb);
    printf("%18s  %d\n","maxGen",this->maxGen);
    printf("%18s  %d\n","nShow",this->nShow);
    printf("%18s  %d\n","seed",this->seed);
    printf("%18s  %d\n","rn",this->rn);
    printf("%18s  %d\n","cs",this->cs);
    if(this->filenameSNP2Genes==NULL){
        printf("%18s  NULL\n","SNP2Genes");
    }
    else{
        printf("%18s  %s\n","SNP2Genes",this->filenameSNP2Genes);
    }
    if(this->filenameAssociatedGenes==NULL){
        printf("%18s  NULL\n","AssociatedGenes");
    }
    else{
        printf("%18s  %s\n","AssociatedGenes",this->filenameAssociatedGenes);
    }
    printf("%18s  %s\n","in",this->filename);
    printf("%18s  %s\n","out",this->filenameO);
}
void PI::printDocument(){
    printf("%18s  %14s  %s\n","parameter","default","description");
    printf("%18s  %14s  %s\n","-type","0","type of the input file.");
    printf("%18s  %14s  %s\n","-cG","0.05","the threshold of p-value(G-test) in generating result.");
    printf("%18s  %14s  %s\n","-cGc","0.05","the threshold of p-value(G-test) change in generating result.");
    printf("%18s  %14s  %s\n","-o","-1","the max order of generated individual (SNP combination) in this program.");
    printf("%18s  %14s  %s\n","-pe","0.8","the probability of considering bioinformation while generating a new individual (SNP combination).");
    printf("%18s  %14s  %s\n","-maxGen","-1","max generation which this program runs.");
    printf("%18s  %14s  %s\n","-nShow","4","control the echo of the program.");
    printf("%18s  %14s  %s\n","-seed","0","random seed of the program.");
    printf("%18s  %14s  %s\n","-cs","0","start the control system or not.");
    printf("%18s  %14s  %s\n","-in","data.txt","filename of the gwas data.");
    printf("%18s  %14s  %s\n","-out","result.txt","filename of the results generated by the program.");
    printf("%18s  %14s  %s\n","-SNP2Genes","NULL","filename of gene mapping.");
    printf("%18s  %14s  %s\n","-AssociatedGenes","NULL","filename of associatedGenes.");

}
void PI::addResult(int * x, int l,double g){
    struct Result * r=(struct Result *)malloc(sizeof(struct Result));
    r->o=l;
    r->x=x;
    r->g=g;
    r->g1=(double *)malloc(sizeof(double)*l);
    int * xx=(int *)malloc(sizeof(int));
    for(int i=0;i<l;i++){
        xx[0]=x[i];
        r->g1[i]=((Mem *)mem)->computeG(xx,1);
    }
    free(xx);
    std::pair<std::set<struct Result *>::iterator,bool> ret=this->results.insert(r);
    if(ret.second){
        this->resultList.push_back(r);
        resultList.sort(compResult);
    }
}
void PI::generateResult(){
    FILE * f=fopen(this->filenameO,"w");
    printf("cost time: %fs\n",(double)(clock()-this->startTime)/CLOCKS_PER_SEC);
    int j=0;
    for(this->itrl=this->resultList.begin();this->itrl!=this->resultList.end()&&j!=this->rn;this->itrl++,j++){
        fprintf(f,"[%d",(*(this->itrl))->x[0]);
        for(int i=1;i<(*(this->itrl))->o;i++){
            fprintf(f,",%d",(*(this->itrl))->x[i]);
        }
        fprintf(f,"]\t[%s",((Mem *)(this->mem))->names[(*(this->itrl))->x[0]]);
        for(int i=1;i<(*(this->itrl))->o;i++){
            fprintf(f,",%s",((Mem *)(this->mem))->names[(*(this->itrl))->x[i]]);
        }
        fprintf(f,"]\t[%e",(*(this->itrl))->g1[0]);
        for(int i=1;i<(*(this->itrl))->o;i++){
            fprintf(f,",%e",(*(this->itrl))->g1[i]);
        }
        fprintf(f,"]\t%e\n",(*(this->itrl))->g);
    }
    fclose(f);
}
void PI::getTP(int & tp, int & fp, int & tn, int & fn, int & po){
    po=-1;
    int j=0;
    for(this->itrl=this->resultList.begin();this->itrl!=this->resultList.end();this->itrl++){
        if(solution->o==(*(this->itrl))->o){
            bool b=true;
            for(int i=0;i<solution->o;i++){
                if(solution->x[i]!=(*(this->itrl))->x[i]){
                    b=false;
                    break;
                }
            }
            if(b){
                po=j;
                break;
            }
        }
        j++;
        if(j==this->rn){
            break;
        }
    }
    int p=this->rn;
    if(p<0){
        p=this->resultList.size();
    }
    else if(p>(int)resultList.size()){
        p=resultList.size();
    }
    if(po==-1){
        tp=0;
        fn=1;
        fp=p-tp;
        tn=this->n*(this->n-1)/2-fp-fn-tp;
    }
    else{
        tp=1;
        fn=0;
        fp=p-tp;
        tn=this->n*(this->n-1)/2-fp-fn-tp;
    }
}
PI::~PI(){
    clearMe();
}
void PI::clearMe(){
    for(this->itrl=this->resultList.begin();this->itrl!=this->resultList.end();this->itrl++){
        free(*(this->itrl));
    }
    resultList.clear();
    results.clear();
    if(solution!=NULL){
        free(solution->x);
        free(solution);
        solution=NULL;
    }
    n=-1;
    gen=0;
}
void PI::buildSolution(){
    char * filenameC=(char *)malloc(sizeof(char)*strlen(filename)+1);
    memcpy(filenameC,filename,sizeof(char)*strlen(filename)+1);
    char * start=filenameC;
    char * t=NULL;
    bool c=true;
    while(c){
        c=false;
        t=strchr(start,'/');
        if(t!=NULL){
            start=t+1;
            c=true;
        }
        t=strchr(start,'\\');
        if(t!=NULL){
            start=t+1;
            c=true;
        }
    }
    c=true;
    std::vector<int> tt;
    while(c){
        t=strchr(start,'_');
        if(t!=NULL){
            *t='\0';
            tt.push_back(atoi(start));
            start=t+1;
        }
        else{
            t=strchr(start,'.');
            if(t!=NULL){
                *t='\0';
                tt.push_back(atoi(start));
                c=false;
            }
        }
    }
    this->solution=(struct Result *)malloc(sizeof(struct Result));
    this->solution->o=tt.size();
    this->solution->x=(int *)malloc(sizeof(int)*this->solution->o);
    for(int i=0;i<this->solution->o;i++){
        this->solution->x[i]=tt[i];
    }
}
