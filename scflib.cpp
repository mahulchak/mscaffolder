#include<iostream>
#include "scf.h"

string xtractcol(string str, char c, int n)
{
int count =0;
int j =0; // j is the index for the string that'll be returned
string elem;
for(unsigned int i=0;i<str.size() && count<n;i++)
        {
        if(str[i] == c)
                {
                count++; //keep a count of the field separator
                }
        if(count == n-1)
                {
                elem.push_back(str[i]);
                j++;
                }
        }
return elem;
}
//////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
int ovlCalculator(vector<int>& q_st, vector<int>& q_end)
{
int ovl = 0;

        for(unsigned int j =0;j<q_st.size();j++)
        {
        ovl = ovl + abs(q_st[j] - q_end[j]);
        }


return ovl;
}
/////////////////////////////////////////////////////////////////////////////////////
void ovlStoreCalculator(asmMerge & merge)
{
	string tempname;
		for(unsigned int i =0;i<merge.r_name.size();i++)
		{
			tempname = merge.r_name[i]+merge.q_name[i];
			merge.ovlStore[tempname] = ovlCalculator(merge.q_st[tempname],merge.q_end[tempname]);
		}
}
////////////////////////////////////////////////////////////////////////////////////
void findChromPartner(asmMerge & merge)
{
	string tempname;
	for(unsigned int i =0;i<merge.q_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		if(merge.storeHomAln[merge.q_name[i]] < merge.ovlStore[tempname]) //if found a longer alignment
		{
			merge.storeHomolog[merge.q_name[i]] = merge.r_name[i];
			merge.storeHomAln[merge.q_name[i]] = merge.ovlStore[tempname];
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////
void storeStart(asmMerge & merge)
{
	string tempname;
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		if(merge.innie[tempname] != true)
		{
if(merge.storeHomolog[merge.q_name[i]] == merge.r_name[i])
{
//cout<<merge.r_name[i]<<"\t"<<merge.q_name[i]<<"\t"<<merge.ref_st[tempname][0]<<endl;
			merge.refStart[merge.r_name[i]].push_back(merge.ref_st[tempname][0]);
			merge.qStoreStart[merge.r_name[i]].push_back(merge.q_name[i]);
}
		}
	}
	
}	

////////////////////////////////////////////////////////////////////////////////
void innieChecker(asmMerge & merge)
{
	string tempname,tempname2;	
	int refEnd1 = 0,refEnd2 =0;
//cout<<merge.r_name.size()<<endl;
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		
		for(unsigned int j =0;j<merge.r_name.size();j++)
		{
                        tempname2 = merge.r_name[j]+merge.q_name[j];
			if((merge.ovlStore[tempname] > merge.ovlStore[tempname2]) && (merge.r_name[i] == merge.r_name[j])) // if the second query length is shorter
			{
if(merge.storeHomolog[merge.q_name[i]] == merge.storeHomolog[merge.q_name[j]])
{
				refEnd1 = merge.ref_end[tempname][merge.ref_st[tempname].size()-1];
				refEnd2 = merge.ref_end[tempname2][merge.ref_st[tempname2].size()-1];
				if(((!(merge.ref_st[tempname][0] > merge.ref_st[tempname2][0])) && (refEnd2<refEnd1)) ||((merge.ref_st[tempname][0] < merge.ref_st[tempname2][0]) && (!(refEnd2>refEnd1))))
				//if((merge.ref_st[tempname][0] < merge.ref_st[tempname2][0]) && (refEnd2<refEnd1))
				{
					merge.innie[tempname2] = true; //a better alignment for this part is tempname
					//cout<<tempname2<<endl;
				}
}
			}
		}
	}
}
					
///////////////////////////////////////////////////////////////////////////////
void fillSeq(fastaSeq & fasta, ifstream& fin) 
{
        string str,str1;
        
        while(getline(fin,str))
        {
                if(str[0] == '>')
                {       
                        fasta.seqName.push_back(str);
//cout<<str<<endl;
                }
                getline(fin,str1);
                fasta.seq[str.substr(1)] = str1; //fasta.seq[str.substr(1)] = str1 needed to remove leading >
		//cout<<fasta.seq[str.substr(1)].size()<<endl;
        }
}
/////////////////////////////////////////////////////////////////////////
void joinList(asmMerge & merge, fastaSeq & genome)
{
	string refName,tempname,scaffold;
	vector<int> multStart;
	vector<int> allStart;
	vector<int> tempStart;
	string filler; 
	unsigned int pos;

	for(int i =0;i<100;i++)
	{
		filler.push_back('N');
	}
//cout<<filler<<endl;
	for (map<string,vector<int> >::iterator it = merge.refStart.begin();it != merge.refStart.end(); it++)
	{
		refName = it->first;
		allStart = merge.refStart[refName];
		sort(allStart.begin(),allStart.end());
//		cout<<refName<<"\t"<<allStart.size()<<endl;
		cout<<">"<<refName<<endl;
		for(unsigned int i =0; i<allStart.size();i++)
		{
			if((i>0) && (allStart[i] != allStart[i-1]))
			{
				pos = findElem(merge.refStart[refName],allStart[i]);
//cout<<pos<<endl;
//cout<<merge.qStoreStart[refName][pos]<<endl;
				tempname = refName + merge.qStoreStart[refName][pos];
				if(merge.strandOri[tempname] == 'F')
				{
					scaffold.append(genome.seq[merge.qStoreStart[refName][pos]]);	
				}
				if(merge.strandOri[tempname] == 'R')
				{
					scaffold.append(revCom(genome.seq[merge.qStoreStart[refName][pos]]));
				}
				if(i != allStart.size()-1) //if it is not the last element
				{
					scaffold.append(filler);
				}
				merge.qToRemove(merge.qStoreStart[refName][pos]);
			}
			if(i == 0)
			{
				pos = findElem(merge.refStart[refName],allStart[i]);
//cout<<pos<<endl;
//cout<<merge.qStoreStart[refName][pos]<<endl;
				tempname = refName + merge.qStoreStart[refName][pos];
				if(merge.strandOri[tempname] == 'F')
                                {
                                        scaffold.append(genome.seq[merge.qStoreStart[refName][pos]]);
                                }
                                if(merge.strandOri[tempname] == 'R')
                                {
                                        scaffold.append(revCom(genome.seq[merge.qStoreStart[refName][pos]]));
//cout<<scaffold.size();
                                }
                 //               scaffold.append(genome.seq[merge.qStoreStart[refName][pos]]); 
				if(i != allStart.size()-1) //if it is not the last element
                                {
                                        scaffold.append(filler);
                                }
				merge.qToRemove(merge.qStoreStart[refName][pos]);

			}
				
		}
		//cout<<scaffold.size()<<endl;
		cout<<scaffold<<endl;
		scaffold.clear();//reset scaffold for the next chromosome
	}	
			 	
}

//////////////returns the index of a vector element//////////////////////////////////////////
unsigned int findElem(vector<int> & v, int & n)
{
unsigned int i =0;

	while(v[i] != n)
	{
		i++;
	}
return i;
}
//////////////////////////check query orientation///////////////////////////////////////////
void oriQ(asmMerge & merge)
{
int oriF =0, oriR =0;
	string tempname;
	for(unsigned int i =0;i<merge.r_name.size();i++)
        {
                tempname = merge.r_name[i] + merge.q_name[i];
		oriR = 0;
		oriF = 0;
		for(unsigned int j = 0;j<merge.q_st[tempname].size();j++)
		{
			if(merge.q_st[tempname][j] > merge.q_end[tempname][j])
			{
				oriR++;
			}
			if(merge.q_st[tempname][j] < merge.q_end[tempname][j])
			{
				oriF++;
			}
		}
		if(oriR > oriF)
		{
			merge.strandOri[tempname] = 'R';
		}
		else
		{
			merge.strandOri[tempname] = 'F';
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////					
string revCom(string & str)
{
string revcom;
unsigned int k = str.size();
        for(int i=k-1; i>-1;--i)
        {
                if(str[i] == 'A')
                {
                revcom.push_back('T');
                }
                if(str[i] == 'T')
                {
                revcom.push_back('A');
                }
                if(str[i] == 'G')
                {
                revcom.push_back('C');
                }
                if(str[i] == 'C')
                {
                revcom.push_back('G');
                }
        }
return revcom;
}
////////////////////////////////////////////////////////////////////////////////////
