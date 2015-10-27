//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//ART_454 -- Artificial Read Transcriber 
//Copyright(c) 2008-2014 Weichun Huang, All Rights Reserved.
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#pragma once
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_pow_int.h>

#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <cstdio>
#include <set>


#include "seqRead.h"


using namespace std;

class art{
public:
    //art();
    static const int min_read_len=20;
    string ref_seq;
    string ref_seq_cmp;
    vector<int> homo_plus; 
    vector<int> homo_minus; 

    int read_len;
    long valid_region;
    void ini_set(int readLen){
        read_len=readLen;
        valid_region=ref_seq.size()-readLen;
    };
    void init(){
        ref_seq_cmp.resize(ref_seq.size());
        size_t size=ref_seq.size();
        for(size_t i=0; i<size; i++){
            ref_seq[i]=toupper(ref_seq[i]);
            size_t k=size-i-1;
            switch(ref_seq[i]){
                case 'A':
                    ref_seq_cmp[k]='T'; break;
                case 'C':
                    ref_seq_cmp[k]='G'; break;
                case 'G':
                    ref_seq_cmp[k]='C'; break;
                case 'T':
                    ref_seq_cmp[k]='A'; break;
                default:
                    ref_seq_cmp[k]='N';
            }
        }
    };

    //add initial homopolymer length vector
    void init_fast(){
        size_t size=ref_seq.size();
        ref_seq_cmp.resize(size);
	if (size==0) return;
       	homo_plus.resize(size);
       	homo_minus.resize(size);
        ref_seq[0]=toupper(ref_seq[0]);
        char aBase=ref_seq[0];
	int polymer_len=1;
        for(size_t i=1; i<=size; i++){
		if(i==size){
		       	int idx=0;
		       	while(polymer_len>0){
			       	homo_plus[i-polymer_len]=polymer_len;
			       	homo_minus[size-i+idx]=polymer_len;
			       	idx++;
			       	polymer_len--;
		       	}
		}
		else{
		       	ref_seq[i]=toupper(ref_seq[i]); 
			if(aBase!=ref_seq[i]){
			       	int idx=0;
			       	while(polymer_len>0){
				       	homo_plus[i-polymer_len]=polymer_len;
				       	homo_minus[size-i+idx]=polymer_len;
				       	idx++;
				       	polymer_len--;
			       	}
			       	aBase=ref_seq[i]; 
				polymer_len=1;
		       	}
		       	else{
			       	polymer_len++;
		       	}

		}
	}

        for(size_t i=0; i<size; i++){
            size_t k=size-i-1;
            switch(ref_seq[i]){
                case 'A':
                    ref_seq_cmp[k]='T'; break;
                case 'C':
                    ref_seq_cmp[k]='G'; break;
                case 'G':
                    ref_seq_cmp[k]='C'; break;
                case 'T':
                    ref_seq_cmp[k]='A'; break;
                default:
                    ref_seq_cmp[k]='N';
            }
        }
/*
cerr<<ref_seq<<endl;
for(int kk=0; kk<homo_plus.size(); kk++){
	cerr<<ref_seq[kk]<<"\t"<<homo_plus[kk]<<endl;
}
cerr<<ref_seq_cmp<<endl;
for(int kk=0; kk<homo_minus.size(); kk++){
	cerr<<ref_seq_cmp[kk]<<"\t"<<homo_minus[kk]<<endl;
}
*/

    };

    bool next_read_indel(seqRead& a_read);
    bool next_pair_read_indel(seqRead& read_1, seqRead& read_2);
    bool next_read(seqRead& a_read);
    //bool next_pair_read(seqRead& read_1, seqRead& read_2, int full_len);
    bool next_pair_read(seqRead& read_1, seqRead& read_2);

    //amplicon sequencing
    bool amp_read(seqRead& a_read);
    bool amp_pair_read(seqRead& read_1, seqRead& read_2);

    static gsl_rng* gsl_R;
    static int gaussain_mean;
    static double gaussain_sigma;

    static void ini_read_pair_rand(int mean, double sigma){
        gsl_rng_default_seed=(unsigned int)time(NULL);
        const gsl_rng_type *rndT=gsl_rng_default;
        gsl_R=gsl_rng_alloc(rndT);
        gaussain_mean=mean;
        gaussain_sigma=sigma;
    };

    static void ini_read_pair_rand(int mean, double sigma, unsigned int gseed){
        gsl_rng_default_seed=gseed;
        const gsl_rng_type *rndT=gsl_rng_default;
        gsl_R=gsl_rng_alloc(rndT);
        gaussain_mean=mean;
        gaussain_sigma=sigma;
    }; 
};


