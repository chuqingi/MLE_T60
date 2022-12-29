#pragma once
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef _EXTRACTDATA_H_
#define _EXTRACTDATA_H_
#define pi acos(-1)

/*The WAVE file format*/
typedef struct tagWAVHEADER {
	uint8_t   ChunkID[4];     
	uint32_t  ChunkSize;      
	uint8_t   Format[4];      
	uint8_t   FmtChunkID[4];  
	uint32_t  FmtChunkSize;   
	uint16_t  AudioFormat;    
	uint16_t  NumChannels;    
	uint32_t  SampleRate;     
	uint32_t  ByteRate;       
	uint16_t  BlockAlign;     
	uint16_t  BitsPerSample;
	uint8_t   DataChunkID[4];
	uint32_t  DataChunkSize;
} WAVHEADER;

/*
Parameters for overall processing schemes.
They are not identical to the frame sizesand frame shift used for the
actual RT60 estimation to demonstrate the integration into a given block
processing scheme.
*/
typedef struct tagSIMPAR {
	double remove_from_avg[2];  //first and last values are not taken into account for averaging
	int fs;                     //sampling frequency
	int block_size;             //block size of overall processing scheme
	int overlap;                //frame shift (overlap) of overall processing scheme
}SIMPAR, * PSIMPAR;

/*Initialize RT60 estimation*/
typedef struct tagRTE_HANDLE {
	int fs;                     //sampling frequency
	int down;                   //rate for downsampling before RT60 estimation to reduce computational complexity
	int N_sub;                  //sub-frame length (after downsampling) 273
	int N_shift;                //frame shift (before downsampling) 137
	int La;                     //no. of considered decay rate factors (= no of. RT60s) 11
	int buffer_size;            //buffer size 400
	int no_bins;                //no. of histogram bins 11
	int hist_counter;           //counter increased if histogram is updated
	int nos_min;                //minimal number of subframes to detect a sound decay 3
	int nos_max;                //maximal number of subframes to detect a sound decay 7
	int N;                      //maximal frame length (after downsampling) 1911
	double bin;                 //step-size for RT60 estimation
	double alpha;               //smoothing factor
	double RT_initial;          //initial RT estimate
	double RT_last;             //last RT60 estimate
	double RT_raw;              //raw RT60 estimate obtained by histogram-approach
	double RTml;                //default RT60 estimate (-1 indicates no new RT60 estimate)
	double* Tquant;             //set of qunatized RT60s considered for maximum search
	double* a;                  //corresponding decay rate factors
	double* hist_limits;        //limits of histogram bins
	double* hist_rt;            //histogram with ML estimates 
	int* buffer;	            //buffer with previous indices to update histogram
}RTE_HANDLE, * PRTE_HANDLE;

typedef struct tagVALIND {
	double val;
	int ind;
}VALIND, * PVALIND;

#endif

VALIND minv(double* x, int n);
VALIND maxv(double* x, int n);
double mean(double* x, int n);
double var(double* x, int n);
double* preprocess(short* x_s, int n);
double* linspace(double d1, double d2, double n);
double max_loglf(double* h, int N, double* a, double* Tquant, int len, double* ll);
double MLE_RT(double* x, int n, PSIMPAR psimpar, double** rt_est, int* rt_est_len,
	double** RT_est, int* Rt_est_len, PRTE_HANDLE rte_handle);
double MLE_RT_frame(double* frame, PRTE_HANDLE par);
int MLE_RT_init(int fs, PRTE_HANDLE prte);