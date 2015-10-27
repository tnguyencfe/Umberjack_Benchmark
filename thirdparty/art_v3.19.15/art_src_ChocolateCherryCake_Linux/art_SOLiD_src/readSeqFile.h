//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//ART_SOLiD -- Artificial Read Transcription 
//Copyright(c) 2008-2011 Weichun Huang All Rights Reserved.
//___________________________________________________________________________
#pragma once

#include <vector>
#include <fstream>
#include <string>
using namespace std;

class readSeqFile{
  private:
      char* fileName;
      ifstream infile;
  public:
      readSeqFile(char* file_name);
      readSeqFile();
      ~readSeqFile();
      //close opened file and open another file
      bool reSetFile(char* file_name);
  	  //read next num_next seq in fasta format
      int next_seq(vector<string>& nextid,vector<string>& nextseq, int num_next);
      int next_seq(vector<string>& nextid,vector<string>& nextseq);
      int next_seq(string& geneID,string& geneSEQ);
      //set pointer to begining
      bool restart();
 
};
