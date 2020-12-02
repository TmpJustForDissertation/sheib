#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include "Mem.h"
#include "chi.h"
#include "Bio.h"
/*
    The implemented source of Mem.h.
*/
Mem::~Mem(){
    free(BITS);
    if(names!=NULL){
        for(int i=0;i<n;i++)
            free(names[i]);
        free(names);
        names=NULL;
    }
    if(mem1!=NULL)
        free(mem1);
    if(mem0!=NULL)
        free(mem0);
    if(type==1){
        free(filenamePed);
        free(filenameFam);
        mSamples.clear();
        free(chrs);
        free(poss);
        free(genos);
    }
    if(g1!=NULL){
        free(g1);
        g1=NULL;
    }
    //fclose(fTest);
}
Mem::Mem(PI * pi){
    this->pi=pi;
    type=pi->type;
    ready=false;
    //fTest=fopen("d:/test.txt","w");
    BITS=(int*)malloc(sizeof(int)*0x10000);
    for(int i=0;i<0x10000;i++){
        BITS[i]=bitCount(i);
    }
    if(type==0){
        filename=pi->filename;
        printf("loading %s\n",filename);
        fp=fopen(filename,"r");
        if(fp==NULL){
            printf("file does not exist\n");
        }
        else{
            //确定样本数目和snp数目
            if(!fillMN()){
                printf("file format error\n");
            }
            else{
                if(pi->nShow>0)
                    printf("n:%d case:%d control:%d\n",n,m1,m0);
                pi->n=n;
                if(!fillMem0()){
                    printf("fillMem error\n");
                }
                else{
                    nofc=0;
                    fillG1();
                    pi->mem=this;
                    ready=true;
                    if(pi->nShow>0)
                        printf("mem built\n");
                }
                //testMem();
                //testTable();
                //testMDR();
                //testFitness();
                //testSnpsSet();
                //int ** pop=initPop(2);
                //printTwoArray(pop,numPop,2);

            }
            fclose(fp);
        }
    }
    else if(type==1){
        filename=pi->filename;
        filenamePed=(char *)malloc(sizeof(char)*(strlen(filename)+6));
        filenameFam=(char *)malloc(sizeof(char)*(strlen(filename)+6));
        memcpy(filenamePed,filename,strlen(filename));
        memcpy(filenameFam,filename,strlen(filename));
        memcpy(filenamePed+strlen(filename),".tped",6);
        memcpy(filenameFam+strlen(filename),".tfam",6);
        //确定样本数目以及样本状态。
        fp=fopen(filenameFam,"r");
        if(fp==NULL){
            printf("%s does not exists\n",filenameFam);
        }
        else{
            printf("read tfam file\n");
            if(!fillM()){
                printf("tfam file format error\n");
                fclose(fp);
            }
            else{
                fclose(fp);
                printf("read tped file\n");
                fp=fopen(filenamePed,"r");
                if(!fillN()){
                    printf("tped file format error\n");
                    fclose(fp);
                }
                else{
                    if(pi->nShow>0)
                        printf("#samples:%d\t#snps:%d\nloading data\n",m,n);
                    pi->n=n;
                    if(!fillMem1()){
                        printf("fillMem error\n");
                    }
                    else{
                        printf("loaded\n");
                        nofc=0;
                        fillG1();
                        pi->mem=this;
                        ready=true;
                    }
                }
            }
        }
    }
    if(pi->o==-1){
        pi->o=recommendOrder();
        if(pi->o<=1){
            printf("the number of samples is too small to recommend an order\n");
        }
        else{
            if(pi->nShow>=0)
                printf("the recommended order is %d\n",pi->o);
        }
    }
}
/*
交换两个snp的位置
，主要用于消除模拟数据每次产生的致病snp都是最后几个所带来的偏性。
*/
void Mem::exchange(int a,int b){
    int ll1=3*l1;
    int ll0=3*l0;
    unsigned long long * tmp;
    if(ll1>ll0)
        tmp=(unsigned long long *)malloc(ll1*sizeof(unsigned long long));
    else
        tmp=(unsigned long long *)malloc(ll0*sizeof(unsigned long long));
    memcpy(tmp,mem1+a*ll1,ll1*sizeof(unsigned long long));
    memcpy(mem1+a*ll1,mem1+b*ll1,ll1*sizeof(unsigned long long));
    memcpy(mem1+b*ll1,tmp,ll1*sizeof(unsigned long long));
    memcpy(tmp,mem0+a*ll0,ll0*sizeof(unsigned long long));
    memcpy(mem0+a*ll0,mem0+b*ll0,ll0*sizeof(unsigned long long));
    memcpy(mem0+b*ll0,tmp,ll0*sizeof(unsigned long long));
    char * t=names[a];
    names[a]=names[b];
    names[b]=t;
    free(tmp);
}
double Mem::computeK2(const int * indexes,const int & l){
    int lg;
    int * ca=getTable(indexes,l,lg);
    double k2=0;
    for(int i=0;i<lg;i++){
        double t=0;
        int ca0=ca[i];
        int ca1=ca[i+lg];
        //int ca0=*ADDRESS2(ca,2,lg,0,i);
        //int ca1=*ADDRESS2(ca,2,lg,1,i);
        if(ca0==0&&ca1==0){
            continue;
        }
        else if(ca0>ca1){
            for(int j=2;j<=ca1;j++){
                t+=log(j);
            }
            for(int j=ca0+1;j<=ca0+ca1+1;j++){
                t-=log(j);
            }
        }
        else{
            for(int j=2;j<=ca0;j++){
                t+=log(j);
            }
            for(int j=ca1+1;j<=ca0+ca1+1;j++){
                t-=log(j);
            }
        }
        k2=k2+t;
    }
    k2=-k2;
    free(ca);
    return k2;
}
void Mem::cutK2(int *& indexes,int & l,double & k2){
    int * x=(int *)malloc(sizeof(int)*(l-1));
    for(int i=0;i<l;i++){
        for(int j=0;j<i;j++)
            x[j]=indexes[j];
        for(int j=i+1;j<l;j++)
            x[j-1]=indexes[j];
        double kk2=computeK2(x,l-1);
        if(kk2<k2){
            free(indexes);
            k2=kk2;
            l=l-1;
            indexes=x;
            return;
        }
    }
    free(x);
}
void Mem::cutK22(int * & indexes, int & l, const double & k2){
    int * x=(int *)malloc(sizeof(int)*(l-1));
    bool * flag=(bool *)malloc(sizeof(bool)*l);
    int c=0;
    for(int i=0;i<l;i++){
        for(int j=0;j<i;j++)
            x[j]=indexes[j];
        for(int j=i+1;j<l;j++)
            x[j-1]=indexes[j];
        double kk2=computeK2(x,l-1);
        if(kk2<k2){
            flag[i]=false;
        }
        else{
            flag[i]=true;
            c++;
        }
    }
    free(x);
    if(c==0){
        free(indexes);
        indexes=NULL;
        l=0;
    }
    else if(c<l){
        x=(int *)malloc(sizeof(int)*c);
        for(int i=0,j=0;i<l&&j<c;i++){
            if(flag[i]){
                x[j++]=indexes[i];
            }
        }
        free(indexes);
        indexes=x;
        l=c;
    }
    free(flag);
}
/*
    SNP组合的随机生成过程：
    初始化一个associatedSNPs集合
    随机从n中选出一个SNP作为第一个SNP
        若被选中的SNP有对应的基因，将对应基因的所有SNP，加入到associatedSNPs中，若存在与此基因相关的基因，则将相关基因的SNP亦追加到associatedSNPs中
    当我们为个体生成第二个元素时，有以下两种情况
        若associatedSNPs为空，随机从n中选择......
        若associatedSNPs不为空
            以概率决定执行以下两种操作之一
                随机从n中......
                随机从associatedSNPs中选择
    更新associatedSNPs
*/
int Mem::ranGen(){
    int * x=(int *) malloc(sizeof(int)*pi->o);
    int ti;
    Bio * bio=(Bio *)pi->bio;
    std::set<int> associatedSnps;
    std::set<int>::iterator iti;
    //start to generate the first element
    ti=rand()%n;
    x[0]=ti;
    //try to update associatedSnps
    if(bio!=NULL&&bio->snp2Genes[ti]!=NULL){
        //the SNP(ti) map to genes
        //for each gene, add all its snps to the associatedSnps
        for(std::vector<int>::iterator it=bio->snp2Genes[ti]->begin();it!=bio->snp2Genes[ti]->end();it++){
            //for each gene of this SNP
            if(bio->associatedGenes!=NULL&&bio->associatedGenes[*it]!=NULL){
                for(std::vector<int>::iterator itt=bio->associatedGenes[*it]->begin();itt!=bio->associatedGenes[*it]->end();itt++){
                    //for each associated gene
                    //if there are some SNPs on the gene
                    if(bio->gene2Snps[*itt]!=NULL){
                        for(std::vector<int>::iterator ittt=bio->gene2Snps[*itt]->begin();ittt!=bio->gene2Snps[*itt]->end();ittt++){
                            //add SNPs of the associated gene to the set
                            associatedSnps.insert(*ittt);
                        }
                    }
                }
            }
            for(std::vector<int>::iterator itt=bio->gene2Snps[*it]->begin();itt!=bio->gene2Snps[*it]->end();itt++){
                //add SNPs of the gene to the set
                associatedSnps.insert(*itt);
            }
        }
        //remove SNP(ti) from associatedSnps
        associatedSnps.erase(ti);
    }
    //generate other elements
    for(int i=1;i<pi->o;i++){
        if(associatedSnps.size()==0){
            ti=rand()%(n-i);
            int j;
            for(j=0;j<i;j++){
                if(ti>=x[j]){
                    ti++;
                }
                else{
                    break;
                }
            }
            if(j<i){
                for(int k=i;k>j;k--){
                    x[k]=x[k-1];
                }
                x[j]=ti;
            }
            else{
                x[i]=ti;
            }
        }
        else{
            double ran=(double)rand();
            ran=ran/RAND_MAX;
            if(ran<pi->pb){
                iti=associatedSnps.begin();
                std::advance(iti,rand()%associatedSnps.size());
                ti=*iti;
                int j;
                for(j=0;j<i;j++){
                    if(ti<x[j]){
                        break;
                    }
                }
                if(j<i){
                    for(int k=i;k>j;k--){
                        x[k]=x[k-1];
                    }
                    x[j]=ti;
                }
                else{
                    x[i]=ti;
                }
            }
            else{
                ti=rand()%(n-i);
                int j;
                for(j=0;j<i;j++){
                    if(ti>=x[j]){
                        ti++;
                    }
                    else{
                        break;
                    }
                }
                if(j<i){
                    for(int k=i;k>j;k--){
                        x[k]=x[k-1];
                    }
                    x[j]=ti;
                }
                else{
                    x[i]=ti;
                }
            }
        }
        //try to update associatedSnps
        if(bio!=NULL&&bio->snp2Genes[ti]!=NULL&&i<pi->o){
            //the SNP(ti) map to genes
            //for each gene, add all its snps to the associatedSnps
            for(std::vector<int>::iterator it=bio->snp2Genes[ti]->begin();it!=bio->snp2Genes[ti]->end();it++){
                //for each gene of this SNP
                if(bio->associatedGenes!=NULL&&bio->associatedGenes[*it]!=NULL){
                    for(std::vector<int>::iterator itt=bio->associatedGenes[*it]->begin();itt!=bio->associatedGenes[*it]->end();itt++){
                        //for each associated gene
                        //if there are some SNPs on the gene
                        if(bio->gene2Snps[*itt]!=NULL){
                            for(std::vector<int>::iterator ittt=bio->gene2Snps[*itt]->begin();ittt!=bio->gene2Snps[*itt]->end();ittt++){
                                //add SNPs of the associated gene to the set
                                associatedSnps.insert(*ittt);
                            }
                        }
                    }
                }
                for(std::vector<int>::iterator itt=bio->gene2Snps[*it]->begin();itt!=bio->gene2Snps[*it]->end();itt++){
                    //add SNPs of the gene to the set
                    associatedSnps.insert(*itt);
                }
            }
            //remove elements which has already in x from associatedSnps
            for(int j=0;j<=i;j++)
                associatedSnps.erase(x[j]);
        }
    }
    double k2=computeK2(x,pi->o);
    if(nofc==0){
        meanK2=k2;
        nofc++;
    }
    else{
        meanK2=(meanK2*nofc+k2)/(nofc+1);
        nofc++;
    }
    if(k2<meanK2){
        int l=pi->o;
        int ll;
        while(true){
            if(l==1){
                break;
            }
            ll=l;
            cutK2(x,ll,k2);
            if(ll==l){
                break;
            }
            else{
                l=ll;
            }
        }
        if(l==1){
            free(x);
        }
        else if(l>1){
            double g=computeG(x,l);
            double t=g1[x[0]];
            for(int i=1;i<l;i++){
                t=std::min(g1[x[i]],t);
            }
            if(t==0){
                free(x);
            }
            else{
                double gc=g/t;
                if(g<=pi->cG&&gc<=pi->cGc){
                    pi->addResult(x,l,g);
                }
                else{
                    free(x);
                }
            }
        }
    }
    else{
        free(x);
    }
    /*
    int * x=(int *)malloc(sizeof(int)*pi->o);
    int ti;
    if(pi->bio==NULL){
        ti=rand()%n;
        x[0]=ti;
        for(int i=1;i<pi->o;i++){
            ti=rand()%(n-i);
            int j;
            for(j=0;j<i;j++){
                if(ti>=x[j]){
                    ti++;
                }
                else{
                    break;
                }
            }
            if(j<i){
                for(int k=i;k>j;k--){
                    x[k]=x[k-1];
                }
                x[j]=ti;
            }
            else{
                x[i]=ti;
            }
        }
    }
    else{
        std::set<int> visited;
        bool fail=true;
        while(fail){
            fail=false;
            //随机选择一个gene，再在基因上随机选择一个SNP
            ti=rand()%((Bio *)pi->bio)->nGenes;
            ti=((Bio *)pi->bio)->gene2Snps[ti]->at(rand()%((Bio *)pi->bio)->gene2Snps[ti]->size());
            visited.insert(ti);
            x[0]=ti;
            for(int i=1;i<pi->o;i++){
                //上一个被选中的点是ti，统计ti上有几个备选点
                //首先获得ti对应的所有基因，从中随机选取一个基因
                ti=((Bio *)pi->bio)->snp2Genes[ti]->at(rand()%((Bio *)pi->bio)->snp2Genes[ti]->size());
                //从被选出的gene上随机选择一个SNP
                ti=((Bio *)pi->bio)->gene2Snps[ti]->at(rand()%((Bio *)pi->bio)->gene2Snps[ti]->size());
                if(visited.find(ti)!=visited.end()){
                    fail=true;
                    visited.clear();
                    break;
                }
                for(int j=0;j<i;j++){
                    if(ti<x[j]){
                        for(int k=i;k>j;k--){
                            x[k]=x[k-1];
                        }
                        x[j]=ti;
                    }
                    else if(j==i-1){
                        x[i]=ti;
                    }
                }
                visited.insert(ti);
            }
        }
    }
    double k2=computeK2(x,pi->o);
    if(nofc==0){
        meanK2=k2;
        nofc++;
    }
    else{
        meanK2=(meanK2*nofc+k2)/(nofc+1);
        nofc++;
    }
    if(k2<meanK2){
        int l=pi->o;
        int ll;
        while(true){
            if(l==1){
                break;
            }
            ll=l;
            cutK2(x,ll,k2);
            if(ll==l){
                break;
            }
            else{
                l=ll;
            }
        }
        if(l==1){
            free(x);
        }
        else if(l>1){
            double g=computeG(x,l);
            double t=g1[x[0]];
            for(int i=1;i<l;i++){
                t=std::min(g1[x[i]],t);
            }
            if(t==0){
                free(x);
            }
            else{
                double gc=g/t;
                if(g<=pi->cG&&gc<=pi->cGc){
                    pi->addResult(x,l,g);
                }
                else{
                    free(x);
                }
            }
        }
    }
    else{
        free(x);
    }
    */
    pi->gen++;
    if(pi->gen==pi->maxGen)
        return 1;
    else
        return 0;
}
double Mem::computeG(const int * indexes,const int & l){
    int lg;
    int * ca=getTable(indexes,l,lg);
    double g=0;
    double p0=((double)m0)/m;
    double df=0;
    int c0,c1;
    double e;
    for(int i=0;i<lg;i++){
        c0=ca[i];
        c1=ca[i+lg];
        //c0=*ADDRESS2(ca,2,lg,0,i);
        //c1=*ADDRESS2(ca,2,lg,1,i);
        if(c0!=0){
            e=(c0+c1)*p0;
            g+=c0*log(c0/e);
        }
        if(c1!=0){
            e=(c0+c1)*(1-p0);
            g+=c1*log(c1/e);
        }
        if(c0!=0||c1!=0)
            df++;
    }
    g=1-chi_square_cdf(g,df-1);
    free(ca);
    return g;
}
void Mem::testK2G(const int * indexes,const int & l){
    double g=computeG(indexes,l);
    double k2=computeK2(indexes,l);
    printf("[%d",indexes[0]);
    for(int i=1;i<l;i++){
        printf(",%d",indexes[i]);
    }
    printf("]\tk2:%e\tg:%e\n",k2,g);
    printTable(indexes,l);
}
void Mem::testAnIndividual(int * & x,int & l){
    testK2G(x,l);
    double k2=computeK2(x,l);
    if(k2<meanK2){
        int ll;
        while(true){
            if(l==1){
                break;
            }
            ll=l;
            cutK2(x,ll,k2);
            if(ll==l){
                break;
            }
            else{
                l=ll;
            }
        }
        if(l==1){
            printf("the individual is removed for over cutting.[%d]\n",x[0]);
        }
        else if(l>1){
            double g=computeG(x,l);
            double t=g1[x[0]];
            for(int i=1;i<l;i++){
                t=std::min(g1[x[i]],t);
            }
            if(t==0){
                printf("the individual is removed for zero pvalue of single snp.\n");
            }
            else{
                double gc=g/t;
                if(g<=pi->cG&&gc<=pi->cGc){
                    printf("the individual generates a result.[%e,%e]\n",g,gc);
                }
                else{
                    printf("the individual is removed for larger g or gc.[%e,%e]\n",g,gc);
                }
            }
        }
    }
    else{
        printf("the individual is removed for the larger k2\n");
    }
}
void Mem::fillG1(){
    g1=(double *)malloc(sizeof(double)*n);
    int * x=(int *)malloc(sizeof(int));
    for(int i=0;i<n;i++){
        x[0]=i;
        g1[i]=computeG(x,1);
    }
    free(x);
}
void Mem::printTable(const int * indexes,const int & l){
    int lg;
    int * c=getTable(indexes,l,lg);
    for(int i=0;i<lg;i++){
        printf("%d\t",c[i]);
    }
    printf("\n");
    for(int i=0;i<lg;i++){
        printf("%d\t",c[i+lg]);
    }
    printf("\n");
    free(c);
}
/*
计算i所对应的1的总数。
*/
int Mem::bitCount(unsigned long long i){
    i = i - ((i >> 1) & 0x5555555555555555);
    i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
    i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
    i = i + (i >> 8);
    i = i + (i >> 16);
    i = i + (i >> 32);
    return (int)i & 0x7f;
}
/*
读文件，将样本数和snp数记录下来
*/
bool Mem::fillMN(){
    //printf("fillMN();\n");
    int c;
    m1=0;
    m0=0;
    n=0;
    do{
        c=fgetc(fp);
        if(c=='\t')
            n++;
    }while(c!='\n');
    do{
        fseek(fp,2*n,1);
        c=fgetc(fp);
        if(c=='0')
            m0++;
        else if(c=='1')
            m1++;
        else if(c==EOF)
            break;
        else{
            printf("phenotype should be 0 or 1\n");
            return false;
        }
        do{
            c=fgetc(fp);
        }while(c!='\n');
    }while(c!=EOF);
    m=m0+m1;
    return true;
}
/*
将文件中存储的数据信息存入内存，三维数组，n*3*m
*/
bool Mem::fillMem0(){
    //printf("fillMem();\n");
    //确定需要多少个单元来存储样本
    l1=(int)ceil(((double)m1)/SIZE);
    l0=(int)ceil(((double)m0)/SIZE);
    mem1=(unsigned long long *)calloc(n*3*l1,sizeof(unsigned long long));
    mem0=(unsigned long long *)calloc(n*3*l0,sizeof(unsigned long long));
    rewind(fp);
    names=(char **)malloc(n*sizeof(char *));
    char buff[SIZE_NAME],c;
    for(int i=0;i<n;i++){
        int j=0;
        do{
            c=fgetc(fp);
            if(c=='\t')
                break;
            buff[j]=c;
            j++;
        }while(true);
        buff[j]='\0';
        names[i]=(char *)malloc(strlen(buff)+1);
        memcpy(names[i],buff,(strlen(buff)+1)*sizeof(char));
    }
    while(fgetc(fp)!='\n');
    for(int i=0,i0=0,i1=0;i<m;i++){
        fseek(fp,2*n,1);
        c=fgetc(fp);
        if(c=='1'){
            //this is a positive sample
            int indexOfVector=i1/SIZE;
            int indexInVector=i1-SIZE*indexOfVector;
            i1++;
            fseek(fp,-2*n-1,1);
            for(int j=0;j<n;j++){
                c=fgetc(fp);
                fgetc(fp);
                unsigned long long * a=mem1+(j*3+c-'0')*l1+indexOfVector;
                *a|=ONE>>indexInVector;
            }
        }
        else if(c=='0'){
            //this is a negative sample
            int indexOfVector=i0/SIZE;
            int indexInVector=i0-SIZE*indexOfVector;
            i0++;
            fseek(fp,-2*n-1,1);
            for(int j=0;j<n;j++){
                c=fgetc(fp);
                fgetc(fp);
                unsigned long long * a=mem0+(j*3+c-'0')*l0+indexOfVector;
                *a|=ONE>>indexInVector;
            }
        }
        else{
            printf("phenotype should be 0 or 1\n");
            for(int j=0;j<n;j++){
                free(names[j]);
            }
            free(names);
            names=NULL;
            free(mem1);
            mem1=NULL;
            free(mem0);
            mem0=NULL;
            return false;
        }
        while(fgetc(fp)!='\n');
    }
    if(mem1==NULL||mem0==NULL)
        return false;
    return true;
}
/*
计算a l对应数组中的1的总数。
*/
int Mem::popCount(const unsigned long long * a,const int & l){
    int r=0;
    for(int i=0;i<l;i++){
        r=r+POP_COUNT(BITS,a[i]);
    }
    return r;
}
/*
得到indexes中所存储的snp组合对应的table
*/
int * Mem::getTable(const int * indexes,const int & l,int & lg){
    lg=1;
    for(int i=0;i<l;i++){
        lg=lg*3;
    }
    int * c=(int *)malloc(2*lg*sizeof(int));
    unsigned long long * a=(unsigned long long *)malloc(lg*l0*sizeof(unsigned long long));
    unsigned long long * b=NULL;
    //初始化统计
    //printf("[%d,%d]\n",indexes[0],indexes[1]);
    for(int i=0;i<lg;i+=3){
        memcpy(a+i*l0,mem0+indexes[0]*3*l0,3*l0*sizeof(unsigned long long));
    }
    //逐一进行与操作
    for(int i=1,step=3;i<l;i++,step*=3){
        for(int j=0;j<lg;j++){
            //a+j*l0&mem0+(indexes[i]*3+((j/step)%3))*l0
            b=mem0+(indexes[i]*3+((j/step)%3))*l0;
            for(int k=0;k<l0;k++){
                a[j*l0+k]&=b[k];
            }
        }
    }
    //结果整合
    for(int i=0;i<lg;i++){
        c[i]=popCount(a+i*l0,l0);
    }
    free(a);
    a=(unsigned long long *)malloc(lg*l1*sizeof(unsigned long long));
    b=NULL;
    //初始化统计
    for(int i=0;i<lg;i+=3){
        memcpy(a+i*l1,mem1+indexes[0]*3*l1,3*l1*sizeof(unsigned long long));
    }
    //逐一进行与操作
    for(int i=1,step=3;i<l;i++,step*=3){
        for(int j=0;j<lg;j++){
            //a&mem0+(indexes[i]*3+((j/step)%3))*l0
            b=mem1+(indexes[i]*3+((j/step)%3))*l1;
            for(int k=0;k<l1;k++){
                a[j*l1+k]&=b[k];
            }
        }
    }
    //结果整合
    for(int i=0;i<lg;i++){
        c[i+lg]=popCount(a+i*l1,l1);
    }
    free(a);
    return c;
}

void Mem::testMem(){
    for(int i=0;i<n;i++){
        int cn=0;
        int cp=0;
        for(int j=0;j<3;j++){
            cp+=popCount(mem1+(i*3+j)*l1,l1);
            cn+=popCount(mem0+(i*3+j)*l0,l0);
        }
        printf("positive %d , negative %d\n",cp,cn);
    }
}
void Mem::testTable(){
    int * indexes=(int *)malloc(sizeof(int)*3);
    const int o=3;
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            for(int kk=j+1;kk<n;kk++){
                indexes[0]=i;
                indexes[1]=j;
                indexes[2]=kk;
                int lg;
                int * c=getTable(indexes,o,lg);
                int cp=0,cn=0;
                for(int k=0;k<lg;k++){
                    cp+=c[lg+k];
                    cn+=c[k];
                }
                if(cp!=800||cn!=800)
                    printf("(%d,%d,%d) => positive %d , negative %d\n",i,j,kk,cp,cn);
            }
        }
    }
}
bool Mem::fillM(){
    char c;
    m=m0=m1=0;
    c=fgetc(fp);
    while(c!=EOF){
        for(int i=0;i<5;i++){
            while(fgetc(fp)!='\t');
        }
        c=fgetc(fp);
        if(c=='1'){
            mSamples.insert(std::pair<int,bool>(m,false));
            m0++;
            m++;
        }
        else if(c=='2'){
            mSamples.insert(std::pair<int,bool>(m,true));
            m1++;
            m++;
        }
        else{
            m++;
            printf("error: unsupported phenotype code\n");
            return false;
        }
        c=fgetc(fp);
        while(c!='\n'&&c!=EOF){
            c=fgetc(fp);
        }
        c=fgetc(fp);
    }
    return true;
}
bool Mem::fillN(){
    char buff[BUFF_SIZE];
    n=0;
    char * tt=fgets(buff,BUFF_SIZE,fp);
    while(tt!=NULL&&strcmp(buff,"\n")!=0&&strcmp(buff,"\r\n")!=0){
        n++;
        tt=fgets(buff,BUFF_SIZE,fp);
    }
    /*
    n=0;
    char c=fgetc(fp);
    while(c!=EOF){
        if(c=='\n'){
            n++;
        }
        c=fgetc(fp);
    }
    fseek(fp,-1,1);
    c=fgetc(fp);
    if(c!='\n'){
        n++;
    }
    */
    return true;
}
bool Mem::fillMem1(){
    //printf("fillMem();\n");
    //确定需要多少个单元来存储样本
    rewind(fp);
    l1=(int)ceil(((double)m1)/SIZE);
    l0=(int)ceil(((double)m0)/SIZE);
    mem1=(unsigned long long *)calloc(n*3*l1,sizeof(unsigned long long));
    if(mem1==NULL)
        return false;
    mem0=(unsigned long long *)calloc(n*3*l0,sizeof(unsigned long long));
    if(mem0==NULL)
        return false;
    names=(char **)malloc(n*sizeof(char *));
    if(names==NULL)
        return false;
    chrs=(char *)malloc(n*sizeof(char));
    if(chrs==NULL)
        return false;
    poss=(long *)malloc(sizeof(long)*n);
    if(poss==NULL)
        return false;
    genos=(char *)malloc(n*2*sizeof(char));
    if(genos==NULL)
        return false;
    std::map<int,bool>::iterator iter;
    unsigned long long * a=NULL;
    int indexOfVector,indexInVector;
    //逐行读取文件，填充SNP的memory
    char geno[2];
    int cGeno[2];
    int gg=0;
    int ps=0;
    int p=0;
    char buff[BUFF_SIZE],c;
    for(int i=0;i<n;i++){
        fgets(buff,BUFF_SIZE,fp);
        //chr name
        if(buff[0]=='X'||buff[0]=='x'){
            chrs[i]=23;
            p=1;
            while(buff[p++]!='\t');
            ps=p;
        }
        else if(buff[0]=='Y'||buff[0]=='y'){
            p=1;
            chrs[i]=24;
            while(buff[p++]!='\t');
            ps=p;
        }
        else{
            p=0;
            while(buff[++p]!='\t');
            buff[p]='\0';
            chrs[i]=atoi(buff);
            ps=p+1;
        }
        //rs id
        p=ps;
        while(buff[++p]!='\t');
        buff[p]='\0';
        names[i]=(char *)malloc(p-ps+1);
        memcpy(names[i],buff+ps,(p-ps+1)*sizeof(char));
        ps=p+1;
        //skip
        p=ps;
        while(buff[p++]!='\t');
        ps=p;
        //position
        p=ps;
        while(buff[++p]!='\t');
        buff[p]='\0';
        poss[i]=atol(buff+ps);
        ps=p+1;
        //process the genotypes
        //count genotypes and determine the genotype code
        p=ps;
        geno[0]=geno[1]=0;
        cGeno[0]=cGeno[1]=0;
        for(int j=0;j<m*2;j++){
            c=buff[p++];
            if(geno[0]==0){
                geno[0]=c;
                cGeno[0]=1;
            }
            else{
                if(c==geno[0]){
                    cGeno[0]++;
                }
                else{
                    if(geno[1]==0){
                        geno[1]=c;
                        cGeno[1]=1;
                    }
                    else{
                        if(c==geno[1]){
                            cGeno[1]++;
                        }
                        else{
                            printf("genotype error\n");
                            return false;
                        }
                    }
                }
            }
            p++;
        }
        if(cGeno[0]<cGeno[1]){
            char cc=geno[0];
            geno[0]=geno[1];
            geno[1]=cc;
        }
        genos[i*2]=geno[0];
        genos[i*2+1]=geno[1];
        //save genotype code to mem
        p=ps;
        for(int j=0, j0=0 ,j1=0;j<m;j++){
            gg=0;
            c=buff[p++];
            if(c==geno[1]){
                gg++;
            }
            p++;
            c=buff[p++];
            if(c==geno[1]){
                gg++;
            }
            p++;
            iter=mSamples.find(j);
            if(iter!=mSamples.end()){
                if(iter->second){
                    /*
                    //this is a positive sample
                    int indexOfVector=i1/SIZE;
                    int indexInVector=i1-SIZE*indexOfVector;
                    i1++;
                    fseek(fp,-2*n-1,1);
                    for(int j=0;j<n;j++){
                        c=fgetc(fp);
                        fgetc(fp);
                        unsigned long long * a=mem1+(j*3+c-'0')*l1+indexOfVector;
                        *a|=ONE>>indexInVector;
                    }
                    */
                    indexOfVector=j1/SIZE;
                    indexInVector=j1-SIZE*indexOfVector;
                    a=mem1+i*3*l1+gg*l1+indexOfVector;
                    *a|=ONE>>indexInVector;
                    j1++;
                }
                else{
                    indexOfVector=j0/SIZE;
                    indexInVector=j0-SIZE*indexOfVector;
                    a=mem0+i*3*l0+gg*l0+indexOfVector;
                    *a|=ONE>>indexInVector;
                    j0++;
                }
            }
            else{
                printf("error in mSamples\n");
            }
        }
        //end a line
    }
    if(mem1==NULL||mem0==NULL)
        return false;
    return true;
}
/*
    尝试自动推荐一个探测的order
    m/pow(3,o)>=4
    x<=log(min(m0,m1))/log(3)-1
*/
int Mem::recommendOrder(){
    double t=std::log(std::min(m0,m1))-0.5;
    return std::floor(t);
}
std::map<std::string,int> Mem::reversedNames(){
    std::map<std::string,int> r;
    for(int i=0;i<n;i++){
        std::string s(names[i]);
        r.insert(std::pair<std::string,int>(s,i));
    }
    return r;
}
