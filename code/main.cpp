/*
    ��sly1.8�޸Ķ�������Ҫ��������
    1.��cutk2���¸�Ϊ�ݹ顢ѭ������ʽ��ԭ���ǣ�1234����������ϣ���12��ȫһ�����Һͱ��͹�ϵ��ǿ��һ�α�����������ȷ�Ľ����
    2.����������Ϣ���ݿ⡣
    3.�����������ݿ⣬����20181128������һ�δ���
    ԭ�����ṩ��bio�ļ�����¼SNP��SNP֮��Ĺ�ϵ
    ���ڣ������ṩ����bio�ļ�����һ��¼SNP��gene֮��Ĺ�ϵ�������¼gene��gene֮��Ĺ�ϵ
    ���������ʶ����һ�㣬�����û���Ĺ�ϵ���棬�е�������ƣ��Ͻ�һЩ��
    SNP��ϵ�������ɹ��̣�
    ��ʼ��һ��associatedSNPs����
    �����n��ѡ��һ��SNP��Ϊ��һ��SNP
        ����ѡ�е�SNP�ж�Ӧ�Ļ��򣬽���Ӧ���������SNP�����뵽associatedSNPs�У���������˻�����صĻ�������ػ����SNP��׷�ӵ�associatedSNPs��
    ������Ϊ�������ɵڶ���Ԫ��ʱ���������������
        ��associatedSNPsΪ�գ������n��ѡ��......
        ��associatedSNPs��Ϊ��
            �Ը��ʾ���ִ���������ֲ���֮һ
                �����n��......
                �����associatedSNPs��ѡ��
    ����associatedSNPs
    ���⣬
    ��20181203��׷�ӱ�д���������ļ��Ĳ��֡�
    ���⣬
    ��bio��Ϣ�ֳ������ļ�������ֻ����snp2genes�ļ�����������associatedGenes�ļ���
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

