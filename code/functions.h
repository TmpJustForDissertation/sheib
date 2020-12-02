#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#define STATUS_RUNNING 0
#define STATUS_PAUSING 1
#define STATUS_PAUSED 2
#define STATUS_STOPPING 3
#include <set>
#include "Mem.h"
/*
    The function class is used to implement functions of sheib.
*/
//parameters to start a control thread.
struct Paras{
    PI * pi=NULL;
    int * status;
};
//Compare results to determine which one is more important.
bool compResult(struct Result * & r1,struct Result * & r2);
//The main function.
int mainReal(PI * pi);
//Control thread function. It is used to perform some debug functions.
//Users can type "Enter" to enter the control system if the parameter cs is specified by 1. Typing "help" in control system will list the commands supported in control system.
void * control(void * arg);
//print parameters specified by users.
void print(PI * pi);
//print the document of sheib.
void printDocument();
/*
    experiments on simulated datasets.
*/
void expOnDNME100();
void expOnDNME3100();
void expOnDNME1000();
void expOnDME100();
void expOnDME1000();
//function called in simulated experiments.
void forAFile(char * filename,FILE * fo,PI * pi);
//function called in simulated experiments.
void forAModel(char * modelName,PI * p);
//function called in simulated experiments.
int mainExperiment(PI * pi,char models[][40],unsigned int lm);
#endif // FUNCTIONS_H_INCLUDED
