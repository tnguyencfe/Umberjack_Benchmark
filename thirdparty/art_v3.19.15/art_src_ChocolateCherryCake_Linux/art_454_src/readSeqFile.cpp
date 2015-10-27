//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//ART -- Artificial Read Transcription 
//Copyright(c) 2008-2012 Weichun Huang All Rights Reserved.
//___________________________________________________________________________

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


#include "readSeqFile.h"

readSeqFile::readSeqFile(){
    fileName = (char*) "";
}
readSeqFile::readSeqFile(char* file_name){
    fileName = file_name;
    infile.open(fileName);
    if(! infile.is_open()){
        cerr<<"the file: "<<file_name<<" can not be opened"<<endl;
        exit(0);
    }
}

readSeqFile::~readSeqFile(){
    if(infile.is_open()){
        infile.close();
    }
}

//read next num_next seq in fast format 
int readSeqFile::next_seq(vector<string>& nextid, vector<string>& nextseq, int num_next)
{
    if(!infile.is_open()){
        cerr<<"There is not open file"<<endl;
        return 0;
    }
    string geneID="";
    string geneSEQ="",tmpString="";
    char tempChar;
    int numSeq=0;
    while(!infile.eof() && numSeq<num_next){
        tempChar = 'a'; 
        while(tempChar!='>' && !infile.eof() ){
            infile.get(tempChar);
        }
        if(tempChar!='>'){
            continue;
        }
        getline(infile,geneID);	
        if(geneID.length()<1){
            continue;
        }
        while((!infile.eof()) && (infile.peek()!='>')){ 
            tmpString.clear();
            getline(infile,tmpString);
	    if ((tmpString.size() > 0) && (tmpString[tmpString.size() - 1] == '\r')){ tmpString.resize(tmpString.size() - 1); }
            geneSEQ.append(tmpString);
        }
        nextseq.push_back(geneSEQ);
        nextid.push_back(geneID);
        geneSEQ.clear();
        geneID.clear();
        numSeq++;
    }
    return numSeq;	
}

int readSeqFile::next_seq(vector<string>& nextid,vector<string>& nextseq)
{
    return next_seq(nextid, nextseq,1);
}
//read in one seq
int readSeqFile::next_seq(string& geneID,string& geneSEQ){
    if(!infile.is_open()){
        cerr<<"There is not open file"<<endl;
        return 0;
    }
    geneSEQ.clear();
    string tmpString="";
    char tempChar;
    int numSeq=0;
    //stringstream tmpseq;
    while(!infile.eof() && numSeq==0){
        tempChar = 'a'; //tempChar is not '>'
        while(tempChar!='>' && !infile.eof() ){
            infile.get(tempChar);
        }
        if(tempChar!='>'){
            continue;
        }
        getline(infile,geneID);	
        if(geneID.length()<1){
            continue;
        }
        while((!infile.eof()) && (infile.peek()!='>')){ 
            tmpString.clear();
            getline(infile,tmpString);
	    if ((tmpString.size() > 0) && (tmpString[tmpString.size() - 1] == '\r')){ tmpString.resize(tmpString.size() - 1); }
            geneSEQ.append(tmpString);//v1.1 change from "geneSEQ.append(tmpString)":
        }
        numSeq++;
    }
    return numSeq;
}

//close opened file and open another file
bool readSeqFile::reSetFile(char* file_name){
    if(infile.is_open()){
        infile.close();
    }
    infile.open(file_name);
    if(! infile.is_open()){
        cerr<<"the file: "<<file_name<<" can not be opened";
        return false;
    }
    return true;

}

//set pointer to begining
bool readSeqFile::restart(){
    infile.clear();
    infile.seekg(0, ios::beg);
    return true;
}
