#include "PI.h"
#include "Mem.h"
#include <dirent.h>
#include <unistd.h>
#include <time.h>
#include "functions.h"
#include "Bio.h"
#include <cstring>
#include <string>
#define SS '/'
/*
    The implemented source code of functions.h.
    The WORKSPACE is specified to the path of simulated datasets. Users should change it if they want to run the simulated experiments.
*/
#define WORKSPACE "d:/workspace/sseh/simulated_data"
//#define WORKSPACE "/home/sly/data/sseh/simulated_data"
int mainReal(PI * pi){
    srand(pi->seed);
    pi->startTime=clock();
    int status=STATUS_RUNNING;
    pthread_t id;
    struct Paras * paras=NULL;
    if(pi->cs==1){
        paras=(struct Paras *)malloc(sizeof(struct Paras));
        paras->status=&status;
        paras->pi=pi;
    }
    //load data
    Mem mem(pi);
    if(pi->mem==NULL){
        printf("mem object has not been generated, the program terminates.\n");
        return 0;
    }
    Bio bio(pi);
    if(pi->cs==1){
        int r=pthread_create(&id,NULL,control,paras);
        if(r==0){
            printf("control thread started, press enter to enter the control system.\n");
        }
        else{
            printf("control thread error\n");
        }
    }
    int s=0;
    int gen=0;
    while(s==0){
        if(status==STATUS_RUNNING){
            s=mem.ranGen();
            gen++;
            if(pi->nShow>0&&gen%pi->nShow==0){
                printf("\rgen[%d]\n",gen);
                if(pi->results.size()>0){
                    int j=0;
                    for(pi->itrl=pi->resultList.begin();pi->itrl!=pi->resultList.end()&&j!=pi->rn;pi->itrl++,j++){
                        printf("[%d",(*(pi->itrl))->x[0]);
                        for(int i=1;i<(*(pi->itrl))->o;i++){
                            printf(",%d",(*(pi->itrl))->x[i]);
                        }
                        printf("]\t[%s",mem.names[(*(pi->itrl))->x[0]]);
                        for(int i=1;i<(*(pi->itrl))->o;i++){
                            printf(",%s",mem.names[(*(pi->itrl))->x[i]]);
                        }
                        printf("]\t[%e",(*(pi->itrl))->g1[0]);
                        for(int i=1;i<(*(pi->itrl))->o;i++){
                            printf(",%e",(*(pi->itrl))->g1[i]);
                        }
                        printf("]\t%e\n",(*(pi->itrl))->g);
                    }
                }
            }
        }
        else if(status==STATUS_PAUSING){
            printf("command system entered, type 'help' to show the help information\n");
            status=STATUS_PAUSED;
            sleep(4);
        }
        else if(status==STATUS_PAUSED){
            sleep(4);
        }
    }
    pi->generateResult();
    if(pi->cs==1){
        status=STATUS_STOPPING;
        printf("press any key to end\n");
        pthread_join(id,NULL);
        free(paras);
    }
    return 0;
}
void * control(void * arg){
    struct Paras * paras=(struct Paras *)arg;
    char * buff=(char *)malloc(sizeof(char)*400);
    char * tb;
    bool enen;
    while(true){
        //gets(buff);
        enen=true;
        for(int i=0;i<400;i++){
            buff[i]=getchar();
            if(buff[i]=='\n'){
                buff[i]='\0';
                enen=false;
                break;
            }
        }
        if(enen){
            while(getchar()!='\n');
            buff[399]='\0';
            printf("the command is too long.\n");
            continue;
        }
        tb=buff;
        if(*(paras->status)==STATUS_STOPPING){
            break;
        }
        else if(*(paras->status)==STATUS_RUNNING){
            printf("waiting to pause the algorithm and enter the command system......\n");
            *(paras->status)=STATUS_PAUSING;
        }
        else if(*(paras->status)==STATUS_PAUSED){
            printf("processing\n");
            /*
            help
            continue
            w
            test 1,2
            */
            char * sp=strchr(tb,' ');
            if(sp==NULL){
                //无参数命令
                if(strcmp(tb,"help")==0){
                    printf("%-18s  %s\n","COMMAND","DESCRIPTION");
                    printf("%-18s  %s\n","help","show the help information.");
                    printf("%-18s  %s\n","continue","exit the command system and continue the algorithm.");
                    printf("%-18s  %s\n","test","type \"test snp_index1,snp_index2,snp_index3\" to test an SNP combination.");
                    printf("%-18s  %s\n","w","write current results to the file which has been specified by the parameter -out");
                    printf("%-18s  %s\n","paras","show the values of parameters.");
                }
                else if(strcmp(tb,"continue")==0){
                    printf("exit the command system and continue the algorithm\n");
                    *(paras->status)=STATUS_RUNNING;
                }
                else if(strcmp(tb,"w")==0){
                    paras->pi->generateResult();
                }
                else if(strcmp(tb,"paras")==0){
                    paras->pi->print();
                }
                else{
                    printf("%s can't be understood\n",tb);
                }
            }
            else{
                *sp='\0';

                //带参数的命令
                if(strcmp(tb,"test")==0){
                    //test 11,13
                    tb=sp+1;
                    //count the commas
                    int l=1;
                    sp=strchr(tb,',');
                    while(sp!=NULL){
                        l++;
                        sp=strchr(sp+1,',');
                    }
                    int * indexes=(int *)malloc(sizeof(int)*l);
                    for(int i=0;i<l-1;i++){
                        sp=strchr(tb,',');
                        *sp='\0';
                        indexes[i]=atoi(tb);
                        tb=sp+1;
                    }
                    indexes[l-1]=atoi(tb);
                    ((Mem *)(paras->pi->mem))->testAnIndividual(indexes,l);
                    //double tv=((Mem *)paras->mem)->computeFitnessTest(indexes,l);
                    //printf("value of the SNP combination is %f\n",tv);
                    free(indexes);
                }
                else{
                    printf("%s can't be understood\n",tb);
                }
            }
            printf("processed\n");
        }
    }
    free(buff);
    return NULL;
}
void expOnDME100(){
    PI pi;
    char models_dme[8][40]={"DME01_1600_100","DME02_1600_100","DME03_1600_100","DME04_1600_100"\
        ,"DME05_1600_100","DME06_1600_100","DME07_1600_100","DME08_1600_100"};
    pi.maxGen=pi.numPopInDME100*pi.maxGenInDME100;
    pi.nShow=-1;
    pi.rn=1;
    //pi.o=6;
    mainExperiment(&pi,models_dme,sizeof(models_dme)/sizeof(models_dme[0]));
}
void expOnDNME100(){
    PI pi;
    char models_dnme[60][40]={"DNME01_1600_100","DNME02_1600_100","DNME03_1600_100","DNME04_1600_100"\
        ,"DNME05_1600_100","DNME06_1600_100","DNME07_1600_100","DNME08_1600_100"\
        ,"DNME09_1600_100","DNME10_1600_100","DNME11_1600_100","DNME12_1600_100"\
        ,"DNME13_1600_100","DNME14_1600_100","DNME15_1600_100","DNME16_1600_100"\
        ,"DNME17_1600_100","DNME18_1600_100","DNME19_1600_100","DNME20_1600_100"\
        ,"DNME21_1600_100","DNME22_1600_100","DNME23_1600_100","DNME24_1600_100"\
        ,"DNME25_1600_100","DNME26_1600_100","DNME27_1600_100","DNME28_1600_100"\
        ,"DNME29_1600_100","DNME30_1600_100","DNME31_1600_100","DNME32_1600_100"\
        ,"DNME33_1600_100","DNME34_1600_100","DNME35_1600_100","DNME36_1600_100"\
        ,"DNME37_1600_100","DNME38_1600_100","DNME39_1600_100","DNME40_1600_100"\
        ,"DNME41_1600_100","DNME42_1600_100","DNME43_1600_100","DNME44_1600_100"\
        ,"DNME45_1600_100","DNME46_1600_100","DNME47_1600_100","DNME48_1600_100"\
        ,"DNME49_1600_100","DNME50_1600_100","DNME51_1600_100","DNME52_1600_100"\
        ,"DNME53_1600_100","DNME54_1600_100","DNME55_1600_100","DNME56_1600_100"\
        ,"DNME57_1600_100","DNME58_1600_100","DNME59_1600_100","DNME60_1600_100"};
    //char models_dnme[1][40]={"DNME51_1600_100"};
    pi.maxGen=pi.numPopInDNME100*pi.maxGenInDNME100;
    pi.nShow=-1;
    pi.rn=1;
    //pi.o=6;
    mainExperiment(&pi,models_dnme,sizeof(models_dnme)/sizeof(models_dnme[0]));
}
void expOnDNME3100(){
    PI pi;
    char models_dnme[40][40]={"DNME301_1600_100","DNME302_1600_100","DNME303_1600_100","DNME304_1600_100"\
        ,"DNME305_1600_100","DNME306_1600_100","DNME307_1600_100","DNME308_1600_100"\
        ,"DNME309_1600_100","DNME310_1600_100","DNME311_1600_100","DNME312_1600_100"\
        ,"DNME313_1600_100","DNME314_1600_100","DNME315_1600_100","DNME316_1600_100"\
        ,"DNME317_1600_100","DNME318_1600_100","DNME319_1600_100","DNME320_1600_100"\
        ,"DNME321_1600_100","DNME322_1600_100","DNME323_1600_100","DNME324_1600_100"\
        ,"DNME325_1600_100","DNME326_1600_100","DNME327_1600_100","DNME328_1600_100"\
        ,"DNME329_1600_100","DNME330_1600_100","DNME331_1600_100","DNME332_1600_100"\
        ,"DNME333_1600_100","DNME334_1600_100","DNME335_1600_100","DNME336_1600_100"\
        ,"DNME337_1600_100","DNME338_1600_100","DNME339_1600_100","DNME340_1600_100"};
    //char models_dnme[1][40]={"DNME51_1600_100"};
    pi.maxGen=pi.numPopInDNME3100*pi.maxGenInDNME3100;
    pi.nShow=-1;
    pi.rn=1;
    //pi.cGc=0.4;
    //pi.o=6;
    mainExperiment(&pi,models_dnme,sizeof(models_dnme)/sizeof(models_dnme[0]));
}
void expOnDME1000(){
    PI pi;
    char models_dme[8][40]={"DME01_1600_1000","DME02_1600_1000","DME03_1600_1000","DME04_1600_1000"\
        ,"DME05_1600_1000","DME06_1600_1000","DME07_1600_1000","DME08_1600_1000"};
    pi.maxGen=pi.numPopInDME1000*pi.maxGenInDME1000;
    pi.nShow=-1;
    pi.rn=1;
    //pi.o=6;
    mainExperiment(&pi,models_dme,sizeof(models_dme)/sizeof(models_dme[0]));
}
void expOnDNME1000(){
    PI pi;
    char models_dnme[60][40]={"DNME01_1600_1000","DNME02_1600_1000","DNME03_1600_1000","DNME04_1600_1000"\
        ,"DNME05_1600_1000","DNME06_1600_1000","DNME07_1600_1000","DNME08_1600_1000"\
        ,"DNME09_1600_1000","DNME10_1600_1000","DNME11_1600_1000","DNME12_1600_1000"\
        ,"DNME13_1600_1000","DNME14_1600_1000","DNME15_1600_1000","DNME16_1600_1000"\
        ,"DNME17_1600_1000","DNME18_1600_1000","DNME19_1600_1000","DNME20_1600_1000"\
        ,"DNME21_1600_1000","DNME22_1600_1000","DNME23_1600_1000","DNME24_1600_1000"\
        ,"DNME25_1600_1000","DNME26_1600_1000","DNME27_1600_1000","DNME28_1600_1000"\
        ,"DNME29_1600_1000","DNME30_1600_1000","DNME31_1600_1000","DNME32_1600_1000"\
        ,"DNME33_1600_1000","DNME34_1600_1000","DNME35_1600_1000","DNME36_1600_1000"\
        ,"DNME37_1600_1000","DNME38_1600_1000","DNME39_1600_1000","DNME40_1600_1000"\
        ,"DNME41_1600_1000","DNME42_1600_1000","DNME43_1600_1000","DNME44_1600_1000"\
        ,"DNME45_1600_1000","DNME46_1600_1000","DNME47_1600_1000","DNME48_1600_1000"\
        ,"DNME49_1600_1000","DNME50_1600_1000","DNME51_1600_1000","DNME52_1600_1000"\
        ,"DNME53_1600_1000","DNME54_1600_1000","DNME55_1600_1000","DNME56_1600_1000"\
        ,"DNME57_1600_1000","DNME58_1600_1000","DNME59_1600_1000","DNME60_1600_1000"};
    pi.maxGen=pi.numPopInDNME1000*pi.maxGenInDNME1000;
    pi.nShow=-1;
    pi.rn=1;
    //pi.o=6;
    mainExperiment(&pi,models_dnme,sizeof(models_dnme)/sizeof(models_dnme[0]));
}
int mainExperiment(PI * pi,char models[][40],unsigned int lm){
    char workspace[400]=WORKSPACE;
    char * modelName=(char *)malloc(400*sizeof(char));
    for(unsigned int i=0;i<lm;i++){
        memcpy(modelName,workspace,strlen(workspace)*sizeof(char));
        modelName[strlen(workspace)]=SS;
        memcpy(modelName+strlen(workspace)+1,models[i],(strlen(models[i])+1)*sizeof(char));
        forAModel(modelName,pi);
    }
    free(modelName);
    return 0;
}
void forAModel(char * modelName,PI * pi){
    struct dirent * dd=NULL;
    char * filename=(char *)malloc(400*sizeof(char));
    char * filenameO=(char *)malloc(400*sizeof(char));
    memcpy(filenameO,modelName,(strlen(modelName)+1)*sizeof(char));
    //strcat(filenameO+strlen(modelName),".result.working.txt");
    strcat(filenameO+strlen(modelName),pi->tail);
    FILE * fo=fopen(filenameO,"w");
    //printf("%s\n",filenameO);
    //FILE * fo=NULL;
    DIR * d=opendir(modelName);
    while((dd=readdir(d))!=NULL){
        if(strcmp(dd->d_name,".")!=0&&strcmp(dd->d_name,"..")!=0){
            memcpy(filename,modelName,strlen(modelName)*sizeof(char));
            filename[strlen(modelName)]=SS;
            memcpy(filename+strlen(modelName)+1,dd->d_name,(strlen(dd->d_name)+1)*sizeof(char));
            forAFile(filename,fo,pi);
        }
    }
    fclose(fo);
    free(filename);
    free(filenameO);
    closedir(d);
}
void forAFile(char * filename,FILE * fo,PI * pi){
    //printf("%s\n",filename);
    int et=0;
    clock_t start,finish;
    start=clock();
    int tp,tn,fp,fn;
    pi->filename=filename;
    Mem mem(pi);
    pi->buildSolution();
    srand(pi->seed);
    int gen=0;
    while(gen<pi->maxGen){
        mem.ranGen();
        gen++;
    }
    finish=clock();
    int po;
    //et,tp,fp,tn,fn,po,gen
    pi->getTP(tp,fp,tn,fn,po);
    printf("%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\n",pi->maxGen,(double)(finish - start)/CLOCKS_PER_SEC,tp,fp,tn,fn,po,gen);
    if(fo!=NULL)
        fprintf(fo,"%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\n",filename,et,(double)(finish - start)/CLOCKS_PER_SEC,tp,fp,tn,fn,po,gen);
    /*
        在真实情况下不需要的步骤，由于同一个pi被反复使用，pi中在程序中被改变的东西必须重置
    */
    pi->clearMe();
}
bool compResult(struct Result * & r1,struct Result * & r2){
    return r1->g<r2->g;
}

