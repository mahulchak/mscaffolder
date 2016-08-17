#include<iostream>
#include<fstream>
#include<cstdlib>
#include "scf.h"

using namespace std;


int main(int argc, char * argv[])
{
        if(argc==1)
        {
		cerr<<"Usage: "<<argv[0]<<" -D master_delta.out -f foo.fasta"<<endl;
        	exit(EXIT_FAILURE);
        }

	ifstream Delta,fasta;
	ofstream chromMap;	
	fastaSeq genome;
	asmMerge master; 
	
	string header,ref_name,qu_name,tempname,name;
	
	int qu_st = 0;
	int qu_end = 0;
	int r_st = 0;
	int r_end = 0;
	
	if(string(argv[2]) != "")
	{
		Delta.open(argv[2]);
		while(getline(Delta,header))
       		{
        	        if(header[0] =='>')
                	{

                        	ref_name = xtractcol(header,' ',1);
	                        ref_name = ref_name.substr(1); //blast adds |c| so remove >|c|
        	                master.r_name.push_back(ref_name);
                	        qu_name = xtractcol(header,' ',2);
                	        qu_name = qu_name.substr(1);//remove the leading white sppace
				master.q_name.push_back(qu_name);
				tempname = ref_name.append(qu_name); // tempname is the index for the map. they describe alignment pairs( e.g. BackboneX ctgY). should be unique.
                       		master.ref_len[tempname] = atoi(xtractcol(header,' ',3).c_str());
                	        master.q_len[tempname] = atoi(xtractcol(header,' ',4).c_str());
               		}

       		        if(header[0] != '>' && header.size()>10)
               		{
                       		r_st = atoi(xtractcol(header,' ',1).c_str());
                     		master.ref_st[tempname].push_back(r_st); // storing the coordinates for each alignment hit
               		      	r_end = atoi(xtractcol(header,' ',2).c_str());
                  		master.ref_end[tempname].push_back(r_end);
                     		qu_st = atoi(xtractcol(header,' ',3).c_str());
                  		master.q_st[tempname].push_back(qu_st);
                  		qu_end = atoi(xtractcol(header,' ',4).c_str());
                     		master.q_end[tempname].push_back(qu_end);
              		 }
       		}
	

		Delta.close();
	
		ovlStoreCalculator(master);//calculate the alignment length for each alignment
		findChromPartner(master);	
		innieChecker(master);
		oriQ(master);
		storeStart(master);
		fasta.open(argv[4]);
		fillSeq(genome,fasta);	
		joinList(master,genome);
		chromMap.open("ctgmap.txt");
		//for(unsigned int i =0; i< master.q_name.size();i++)
		int chromCount = 0;
		for(map<string,string>::iterator it = genome.seq.begin();it!= genome.seq.end(); it++)
		{
			//if(find(master.qToRemove.begin(),master.qToRemove.end(),master.q_name[i]) == master.qToRemove.end())
			if(find(master.qToRemove.begin(),master.qToRemove.end(),it->first) == master.qToRemove.end())
			{
				if(!genome.seq[it->first].empty()) //if the sequence is not empty already
				{
					cout<<">U_"<<chromCount++<<endl;
					cout<<it->second<<endl;
					//genome.seq[it->first] = "";
					chromMap<<"Unscaffolded"<<"\t"<<chromCount<<"\t"<<it->first<<endl;
				}
			}
		}
		//chromMap.open("ctgmap.txt");
		for(map<string,string>::iterator it=master.storeHomolog.begin();it!= master.storeHomolog.end();it++)
		{
			chromMap<<"Map"<<"\t"<< it->first << "\t"<< it->second <<endl;
		}
		chromMap.close();
	}

	return 0;
}
