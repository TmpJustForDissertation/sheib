/*
    由sly1.8修改而来，主要改两处：
    1.将cutk2重新改为递归、循环的形式，原因是，1234，这样的组合，若12完全一样，且和表型关系很强，一次遍历做不到正确的结果。
    2.引入生物信息数据库。
    3.关于生物数据库，我在20181128打算做一次大修
    原来：提供的bio文件，记录SNP与SNP之间的关系
    现在：我想提供两个bio文件，其一记录SNP与gene之间的关系，其二记录gene与gene之间的关系
    我最近又意识到，一点，在利用基因的关系方面，有点难以设计，严谨一些：
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
    此外，
    于20181203，追加编写输出结果到文件的部分。
    此外，
    将bio信息分成两个文件，可以只包含snp2genes文件，而不包含associatedGenes文件。
*/
#include <iostream>
#include "PI.h"
#include "functions.h"
using namespace std;

int main(int argc, char *argv[]){
    printf("rand max : %d this application assumes that this value is more larger than the number of snp.\n",RAND_MAX);
    printf("size of(unsigned long long) : %lu this application assumes that this value is 8.\n",sizeof(unsigned long long));
    PI pi;
    if(pi.fillParameters(argc,argv)){
        //pi.sim=1;
        //pi.cs=1;
        //pi.filename=(char*)"10_44_51.txt";
        //pi.nShow=4;
        //pi.o=6;
        //pi.filenameAssociatedGenes=(char *)"gene_pairs_after_4.txt";
        if(pi.sim==0){
            pi.print();
            mainReal(&pi);
        }
        else{
            expOnDME100();
            expOnDNME100();
            expOnDME1000();
            expOnDNME1000();
            expOnDNME3100();
        }
    }
    else{
        printf("fill parameters error\n");
        pi.printDocument();
    }
    return 0;
}

