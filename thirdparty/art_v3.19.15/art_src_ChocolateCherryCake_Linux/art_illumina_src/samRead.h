/*
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ART_ILLUMINA -- Artificial Read Transcription 
Copyright(c) 2011 Weichun Huang All Rights Reserved.
_____________________________________________________________________________________________________________
*/
#pragma once

#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <cstdio>
#include <set>
#include <algorithm>
#include "readSeqFile.h"

using namespace std;



class samHeader{
	public:
		samHeader(){ 
			VN="1.4";
		       	SO="unsorted";	
		}

	       	//@HD
		string VN; //=1.4
	       	string SO; //=unsorted	

		//@SQ
		vector<string> SN;
	       	vector<long> LN;
		
	       	//@PG
		string ID;
	       	string PN;
	       	string CL;

		void getRefseqID(char* seqfile);
		void printHeader(ostream& fout);
		void printAlnHeader(ostream& fout);
};

class samRead{

	public:
		samRead(){
			flag=0;
			mapQ=99;
		       	rNext="*";
		       	pNext=0;
			tLen=0;
		};
		string qname;
		int flag;
		string rname;
		long pos;
		short mapQ;
		string cigar;
		string rNext;
		long pNext;
		int tLen;
		string seq;
		string qual;
	       	void reverse_comp();
		void getCigar(string & aln_ref, string& aln_read, bool use_M=false);
		void printRead(ostream& fout);
};
