
#ifndef SCF_H_
#define SCF_H_

#include<vector>
#include<string>
#include<map>
#include<fstream>
#include<algorithm>
#include<cstdlib>


using namespace std;


class asmMerge {

public:

vector<string> r_name;
vector<string> q_name;
map<string,int> ref_len;
map<string,int> q_len;
map<string,vector<int> > ref_st;
map<string,vector<int> > ref_end;
map<string,vector<int> > q_st;
map<string,vector<int> > q_end;
map<string,vector<string> > storeName;
map<string,int> ovlStore;
map<string,string> storeHomolog;
map<string,int> storeHomAln;
map<string,vector<int> > refStart;
map<string,vector<string> > qStoreStart;
map<string,bool> innie;
map<string,char> strandOri;
vector<string> qToRemove; // scaffolded contigs
map<string,int> new_refSt;
map<string,int> storeAlnCov;
};

class fastaSeq{
public:
vector<string> seqName;
map<string,string> seq;
};


string xtractcol(string str,char c, int n);
void fillSeq(fastaSeq & fasta, ifstream& fin);
int ovlCalculator(asmMerge & merge,string & index, vector<int>& ref_st,vector<int> & ref_end, vector<int>& q_st, vector<int>& q_end);
void ovlStoreCalculator(asmMerge & merge);
void findChromPartner(asmMerge & merge);
void storeStart(asmMerge & merge);
void innieChecker(asmMerge & merge);
void joinList(asmMerge & merge, fastaSeq & genome);
unsigned int findElem(vector<int> & v, int & n);
void oriQ(asmMerge & merge);
string revCom(string & str);
int findCoverage(asmMerge & merge, string & tempname,string & tempname2);
bool chkDup(asmMerge & merge,string & r_name, int & new_refSt,string & q_name);
#endif
