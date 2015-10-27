//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//ART_ILLUMINA -- Artificial Read Transcription 
//Copyright(c) 2008-2014 Weichun Huang All Rights Reserved.
//___________________________________________________________________________
#include "seqRead.h"

int seqRead::get_indel(int read_len){
    //if(ins_rate.size()>=read_len) {cerr<<"fatal error\n";  exit(1)};
    int ins_len=0, del_len=0;
    //deletion
    for(int i=(int)del_rate.size()-1; i>=0; i--){
        if(del_rate[i]>=r_prob()){
            del_len=i+1;
            for(int j=i; j>=0;){
                int pos=(int) floor((read_len-1)*r_prob()); //invalid deletion positions: 0 or read_len-1
                if(pos==0) continue;
                if(indel.count(pos)==0){
                    indel[pos]='-';
                    j--;
                }
            }
            break;
        }
    }

    for(int i=ins_rate.size()-1; i>=0; i--){
        if(ins_rate[i]>=r_prob()){
            ins_len=i+1;
            for(int j=i; j>=0;){
                int pos=(int) floor(r_prob()*read_len);
                if(indel.count(pos)==0){
                    short base=(short)ceil(r_prob()*4);
                    switch(base){
                                case 1:
                                    indel[pos]='A';   break;
                                case 2:
                                    indel[pos]='C';   break;
                                case 3:
                                    indel[pos]='G';   break;
                                case 4:
                                    indel[pos]='T';  
                    }
                    j--;
                }
            }
            break;
        }
    }
    return (ins_len-del_len);
};

void seqRead::ref2read(){
    if(indel.size()==0){
        seq_read=seq_ref;
        return;
    }
    seq_read.clear();
    int k=0;
    for(size_t i=0; i<seq_ref.size();){
        //cout<<i<<"\t"<<k<<endl;
        if(indel.count(k)==0){
            seq_read.push_back(seq_ref[i]); i++; k++; 
        }
        else if(indel[k]=='-'){
            i++;k++;
        }
        else{
            seq_read.push_back(indel[k]); k++;
        }
    }
    while(indel.count(k)>0){
        seq_read.push_back(indel[k]);
        k++;
    }
}

