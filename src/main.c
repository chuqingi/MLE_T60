/****************************************************************************
Developer: Jiaxin Li
E-mail: 1319376761@qq.com
Github: https://github.com/chuqingi/MLE_T60
Integrated Development Environment(IDE): Microsoft Visual Studio 2019
Title: main.c
Version: 1.0
Description: Maximum-Likelihood (ML) estimation reverberation time (RT60) from reverberant speech
Reference: An improved algorithm for blind reverberation time estimation
****************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include "Common.h"

int main(int argc, char** argv) {
	WAVHEADER FileHeader;
	FILE* fp;
	FILE* re;
	if (argc != 2) {
		fprintf(stderr, "usage: %s revspeech.wav\n", argv[0]);
		return 1;
	}
	fp = fopen(argv[1], "rb");
	if (!fp) {
		fprintf(stderr, "File open is failed!\n");
		return 1;
	}
	re = fopen("Result.txt", "w+");
	fread(&FileHeader, 1, sizeof(WAVHEADER), fp);
	int fs = FileHeader.SampleRate;
	int datalen = FileHeader.DataChunkSize / sizeof(short);
	short* revspeech_s = (short*)calloc(datalen, sizeof(short));
	fread(revspeech_s, sizeof(short), datalen, fp);
	double* revspeech = NULL;
	/*preprocessing*/
	revspeech = preprocess(revspeech_s, datalen);
	SIMPAR simpar;
	simpar.remove_from_avg[0] = datalen / fs / 4.0;
	simpar.remove_from_avg[1] = 0;
	simpar.fs = fs;
	simpar.block_size = lround(20e-3 * fs);
	simpar.overlap = lround(simpar.block_size / 2.0);
	double rt_est_mean = 0;                                       //mean RT60 estimate determined from rt_est
	double* rt_est = NULL;                                        //estimated RT60 over time where the time interval for the estimates is given by rte_handle.N_shift / simpar.fs
	int rt_est_len = 0;
	double* RT_est = NULL;                                        //estimated RT60 over time where the time interval for the estimates is given by simpar.overlap / simpar.fs
	int Rt_est_len = 0;
	RTE_HANDLE rte_handle;                                        //struct containing all buffers, variables and parameters used for frame - wise RT60 estimation
	/*Estimation RT60 based on reverberant speech*/
	rt_est_mean = MLE_RT(revspeech, datalen, &simpar, &rt_est, &rt_est_len, &RT_est, &Rt_est_len, &rte_handle);
	fprintf(re, "RT60:%lfs", rt_est_mean);

	if (rte_handle.Tquant)
		free(rte_handle.Tquant);
	if (rte_handle.a)
		free(rte_handle.a);
	if (rte_handle.buffer)
		free(rte_handle.buffer);
	if (rte_handle.hist_limits)
		free(rte_handle.hist_limits);
	if (rte_handle.hist_rt)
		free(rte_handle.hist_rt);
	if (rt_est)
		free(rt_est);
	if (RT_est)
		free(RT_est);
	if (revspeech_s)
		free(revspeech_s);
	if (revspeech)
		free(revspeech);
	fclose(re);
	fclose(fp);
	return 0;
}
