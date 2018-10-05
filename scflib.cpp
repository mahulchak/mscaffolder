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
//cout<<merge.r_name[i]<<"\t"<<merge.q_name[i]<<endl;
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////
void storeStart(asmMerge & merge, asmMerge & merge1, char c) //character c tells whether to use user supplied order or not
{
	string tempname;
	
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		if(merge1.innie[tempname] != true) // using the filtering criteria obtained from the mdelta file
		{
			if(merge.storeHomolog[merge.q_name[i]] == merge.r_name[i])
			{
//cout<<merge.r_name[i]<<"\t"<<merge.q_name[i]<<"\t"<<merge.ref_st[tempname][0]<<endl;
				if(merge.new_refSt[tempname] != 0)
				{
					if(c == 'y')
					{
						merge.refStart[merge.r_name[i]].push_back(findElem(merge.qList,merge.q_name[i]));//storing the index from the sorted query name list

					}
					if(c=='n')
					{
						merge.refStart[merge.r_name[i]].push_back(merge.new_refSt[tempname]);
//cout<<merge.r_name[i]<<"\t"<<merge.q_name[i]<<endl;
					}
					merge.qStoreStart[merge.r_name[i]].push_back(merge.q_name[i]);
//cout<<merge.r_name[i]<<"\t"<<merge.q_name[i]<<endl;
//cout<<merge.q_name[i]<<"\t"<<findElem(merge.qList,merge.q_name[i])<<endl;
				}
				if(merge.new_refSt[tempname] == 0)
				{
					if(c=='y')
					{
						merge.refStart[merge.r_name[i]].push_back(findElem(merge.qList,merge.q_name[i]));
					}
					if(c == 'n')
					{
						merge.refStart[merge.r_name[i]].push_back(merge.ref_st[tempname][0]);//take the midpoint instead of the first
					}
					merge.qStoreStart[merge.r_name[i]].push_back(merge.q_name[i]);
				}
			}
		}
	}
	
}	
////////////////////////////////////////////////////////////////////////////////
void innieChecker(asmMerge & merge, asmMerge & merge1)
{
	string tempname,tempname2;	
	int refEnd1 = 0,refEnd2 =0;
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		
		for(unsigned int j =0;j<merge.r_name.size();j++)
		{
                        tempname2 = merge.r_name[j]+merge.q_name[j];
			if((merge.ovlStore[tempname] > merge.ovlStore[tempname2]) && (merge.r_name[i] == merge.r_name[j])) // if the second query length is shorter
			{
				if((merge1.storeHomolog[merge.q_name[i]] == merge1.storeHomolog[merge.q_name[j]])&&(merge1.ref_st[tempname].size()!=0) && (merge1.ref_st[tempname2].size()!=0)) //both of their best reference hits are same
				{	
					refEnd1 = merge1.ref_end[tempname][merge1.ref_st[tempname].size()-1];
					refEnd2 = merge1.ref_end[tempname2][merge1.ref_st[tempname2].size()-1];
					if(((!(merge1.ref_st[tempname][0] > merge1.ref_st[tempname2][0])) && (refEnd2<refEnd1)) ||((merge1.ref_st[tempname][0] < merge1.ref_st[tempname2][0]) && (!(refEnd2>refEnd1))))	
					{

						if(findCoverage(merge,tempname,tempname2) > int(0.8*merge.ovlStore[tempname2]))
			
							{
								merge.innie[tempname2] = true; //a better alignment for this part is tempname
								//cout<<tempname2<<endl;
							}
					}
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////
void innieChecker(asmMerge & merge)
{
        string tempname,tempname2;
        int refEnd1 = 0,refEnd2 =0;
        for(unsigned int i =0;i<merge.r_name.size();i++)
        {
                tempname = merge.r_name[i] + merge.q_name[i];

                for(unsigned int j =0;j<merge.r_name.size();j++)
                {
                        tempname2 = merge.r_name[j]+merge.q_name[j];
                        if((merge.ovlStore[tempname] > merge.ovlStore[tempname2]) && (merge.r_name[i] == merge.r_name[j])) // if the second query length is shorter
                        {
                                if(merge.storeHomolog[merge.q_name[i]] == merge.storeHomolog[merge.q_name[j]]) //both of their best reference hits are same
                                {
                                        refEnd1 = merge.ref_end[tempname][merge.ref_st[tempname].size()-1];
                                        refEnd2 = merge.ref_end[tempname2][merge.ref_st[tempname2].size()-1];
					 if(((!(merge.ref_st[tempname][0] > merge.ref_st[tempname2][0])) && (refEnd2<refEnd1)) ||((merge.ref_st[tempname][0] < merge.ref_st[tempname2][0]) && (!(refEnd2>refEnd1))))
					{
					 	if(findCoverage(merge,tempname,tempname2) > int(0.9*merge.ovlStore[tempname2]))//if more than 90% is covered by another sequence

                                                        {
                                                                merge.innie[tempname2] = true; //a better alignment for this part is tempname
							}
					}
				}
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
void fillSeq(fastaSeq & fasta, ifstream& fin) 
{
        string str,index;
        
        while(getline(fin,str))
        {
                if(str[0] == '>')
                {       
                        fasta.seqName.push_back(str);
			index = str.substr(1);
//cout<<str<<endl;
                }
		if(str[0] != '>')
		{
                	fasta.seq[index].append(str); //fasta.seq[str.substr(1)] = str1 needed to remove leading >
		}
	
        }
}
/////////////////////////////////////////////////////////////////////////
void joinList(asmMerge & merge, fastaSeq & genome, char c)
{
	string refName,tempname,scaffold,revseq;
	vector<int> multStart;
	vector<int> allStart;
	vector<int> tempStart;
	string filler; 
	unsigned int pos;

	for(int i =0;i<100;i++)
	{
		filler.push_back('N');
	}
	for (map<string,vector<int> >::iterator it = merge.refStart.begin();it != merge.refStart.end(); it++)
	{
		refName = it->first;
		allStart = merge.refStart[refName];
		sort(allStart.begin(),allStart.end());
//		cout<<refName<<"\t"<<allStart.size()<<endl;
		cout<<">"<<refName<<endl;
		for(unsigned int i =0; i<allStart.size();i++)
		{
			if((i>0) && (allStart[i] != allStart[i-1])) //if the two contigs don't have the same start site on the reference
			{
				pos = findElem(merge.refStart[refName],allStart[i]);
//cout<<pos<<endl;
//cout<<merge.qStoreStart[refName][pos]<<endl;
				if(c == 'n')
				{
					tempname = refName + merge.qStoreStart[refName][pos]; //find the index corresponding to the query
				}
				if (c == 'y')
				{
					tempname = refName + merge.qList[allStart[i]];
				}
				if(merge.strandOri[tempname] == 'F')
				{	
					if(c == 'n')
					{
						scaffold.append(genome.seq[merge.qStoreStart[refName][pos]]);	
					}
					if(c == 'y')
					{
						scaffold.append(genome.seq[merge.qList[allStart[i]]]);
					}
				}
				if(merge.strandOri[tempname] == 'R')
				{
					if(c == 'n')
					{
						revseq = revCom(genome.seq[merge.qStoreStart[refName][pos]]);
					}
					if(c == 'y')
					{
						revseq = revCom(genome.seq[merge.qList[allStart[i]]]);
					}
					scaffold.append(revseq);
				}
				if(c == 'n')
				{
					merge.qToRemove.push_back(merge.qStoreStart[refName][pos]);
				}
				if(c == 'y')
				{
					merge.qToRemove.push_back(merge.qList[allStart[i]]);
				}
			}
			if((i>0) && (allStart[i] == allStart[i-1]))
			{
//cout<<pos<<endl;
//cout<<merge.qStoreStart[refName][pos]<<endl;
			}
			if(i == 0)
			{
				pos = findElem(merge.refStart[refName],allStart[i]);
//cout<<pos<<endl;
//cout<<merge.qStoreStart[refName][pos]<<endl;
				if(c == 'n')
				{
					tempname = refName + merge.qStoreStart[refName][pos];
				}
				if(c == 'y')
				{
					tempname = refName + merge.qList[allStart[i]];
				}
				if(merge.strandOri[tempname] == 'F')
                                {
					if(c == 'n')
					{
                                       		scaffold.append(genome.seq[merge.qStoreStart[refName][pos]]);
					}
					if(c == 'y')
					{
                                        	scaffold.append(genome.seq[merge.qList[allStart[i]]]);
					}
                                }
                                if(merge.strandOri[tempname] == 'R')
                                {
					if(c == 'n')
					{
                                        	scaffold.append(revCom(genome.seq[merge.qStoreStart[refName][pos]]));
					}
					if(c == 'y')
					{
                                        	scaffold.append(revCom(genome.seq[merge.qList[allStart[i]]]));
					}
                                }
			if(c == 'n')
			{
				merge.qToRemove.push_back(merge.qStoreStart[refName][pos]);
			}
			if(c == 'y')
			{
				merge.qToRemove.push_back(merge.qList[allStart[i]]);
			}

			}
			if(i != allStart.size()-1) //if it is not the last element
                                {
                                        scaffold.append(filler);
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
/////////////////////////////////////////////////////////////////////////////////////////
unsigned int findElem(vector<string> & v, string & n)
{
unsigned int i = 0;
	while((v[i] != n) && (i <v.size()-1))
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
		if(str[i] == 'N')
		{
		revcom.push_back('N');
		}
		if(str[i] == 'a')
                {
                revcom.push_back('t');
                }
                if(str[i] == 't')
                {
                revcom.push_back('a');
                }
                if(str[i] == 'g')
                {
                revcom.push_back('c');
                }
                if(str[i] == 'c')
                {
                revcom.push_back('g');
                }
                if(str[i] == 'n')
                {
                revcom.push_back('n');
                }

        }
return revcom;
}
////////////////////////////////////////////////////////////////////////////////////
int findCoverage(asmMerge & merge, string & tempname,string & tempname2)
{
	int cov = 0;
	int reIndex = 0;
	int qStart = 0, qEnd = 0;

	for(unsigned int j = 0;j<merge.ref_st[tempname2].size();j++)
	{
		qStart = merge.ref_st[tempname2][j];
		qEnd = merge.ref_end[tempname2][j];
		for(unsigned int i = 0;i<merge.ref_st[tempname].size();i++)
		{
			if(!(merge.ref_end[tempname][i] < qStart) && !(merge.ref_st[tempname][i] > qEnd))
			{
				//cov = cov + abs(max(merge.ref_st[tempname][i],qStart) - min(merge.ref_end[tempname][i],merge.ref_end[tempname2][0]));	
				cov = cov + abs(max(merge.ref_st[tempname][i],qStart) - min(merge.ref_end[tempname][i],qEnd));
			}
			if((!(merge.ref_st[tempname][i] > qStart)) && (i<int(0.5*merge.ref_st[tempname].size())))
			{
				reIndex = i;
			}
			
	
		}
	}

	if(cov == 0)
	{
		merge.new_refSt[tempname] = merge.ref_st[tempname][reIndex+1];
//cout<<"innie"<<"\t"<<tempname2<<"\t"<< tempname<<"\t"<<merge.new_refSt[tempname]<<endl;
	}
	return cov;
}
///////////////////////////////////////////////////////////////////////////////////
vector<string> findQlist(ifstream & ctgList)
{
string str,chromName;

vector<string> qChromList;

//int x = 0;

bool ytics=false;

size_t pos1;
        
        while(getline(ctgList,str))
        {
                pos1 = str.find('\t');
                if(str[0] != '*')
                {
                        chromName = str.substr(0,pos1);
                }
                if(str[0] == '*')
                {
                        chromName = str.substr(1,pos1-1);
                }
                if((chromName != "set") && (ytics == true))
                {
//                        x = stoi(str.substr(pos1));
                        qChromList.push_back(chromName);
//cout<<chromName<<chromName.size()<<endl;
                }
                if(chromName == "set")
                {
                        ytics = true;
                }
        }               

return qChromList;
}

