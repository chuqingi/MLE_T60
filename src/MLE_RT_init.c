/****************************************************************************
Initialize RT60 estimation
provides struct with all parameters and buffers needed for the function
MLE_RT.c to perform a blind RT60 estimation by frame-wise processing
****************************************************************************/
#include "Common.h"

int MLE_RT_init(int fs, PRTE_HANDLE prte) {
	prte->fs = fs;
	double no = fs / 24e3;
	if (fs < 8e3 || fs>24e3) {
		fprintf(stderr, "Algorithm has not been tested for this sampling frequency!\n");
	}
	/* pararmeters for pre-selection of suitable segments*/
	if (fs > 8e3)
		prte->down = 2;
	else
		prte->down = 1;
	prte->N_sub = lround(no * 820 / prte->down);
	prte->N_shift = lround(prte->N_sub * (prte->down / 4.0));
	prte->nos_min = 3;
	prte->nos_max = 7;
	prte->N = prte->nos_max * prte->N_sub;
	/*parameters for ML-estimation*/
	double Tmax = 1.2, Tmin = 0.2;
	prte->bin = 0.1;
	int len = (int)floor((Tmax - Tmin) / prte->bin) + 1;
	prte->Tquant = (double*)calloc(len, sizeof(double));
	prte->a = (double*)calloc(len, sizeof(double));
	if (prte->Tquant == NULL) {
		fprintf(stderr, "Unable to allocate memory!\n");
		return 1; 
	}
	if (prte->a == NULL) {
		fprintf(stderr, "Unable to allocate memory!\n");
		return 1;
	}
	double Ti = Tmin;
	int i = 0;
	double bintmp = prte->bin;
	while (Ti <= Tmax) {
		*(prte->Tquant + i) = Ti;
		*(prte->a + i) = exp(-3 * log(10) / (Ti * fs / prte->down));
		i++;
		Ti = Ti + bintmp;
	}
	if (i == len) {
		prte->La = i;
	}
	else {
		fprintf(stderr, "Tquant array length error!\n");
	}
	/*paramters for histogram-based approach to reduce outliers (order statistics)*/
	prte->buffer_size = lround(no * 1200.0 / prte->down);
	prte->buffer = (int*)calloc(prte->buffer_size, sizeof(int));
	prte->no_bins = prte->La;
	len = (int)floor((Tmax + prte->bin - Tmin) / prte->bin) + 1;
	prte->hist_limits = (double*)calloc(len, sizeof(double));
	i = 0;
	Ti = Tmin - prte->bin / 2.0;
	double tmp = Tmax + prte->bin / 2.0;
	while (Ti <= tmp) {
		prte->hist_limits[i++] = Ti;
		Ti = Ti + bintmp;
	}
	if (i != len) {
		fprintf(stderr, "hist_limits array length error!\n");
	}
	prte->hist_rt = (double*)calloc(prte->no_bins, sizeof(double));
	prte->hist_counter = 0;
	/*paramters for recursive smoothing of final RT60 estimate*/
	prte->alpha = 0.996;
	prte->RT_initial = 0.3;
	prte->RT_last = prte->RT_initial;
	prte->RT_raw = prte->RT_initial;
	return 0;
}