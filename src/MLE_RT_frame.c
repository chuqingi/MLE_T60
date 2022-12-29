/****************************************************************************
RT60 estimation by frame-wise processing
performs blind estimation of the reverberation time (RT60) by frame-wise processing
****************************************************************************/
#include "Common.h"

double MLE_RT_frame(double* frame, PRTE_HANDLE par) {
	int cnt = 0;                                                           //sub-frame counter for pre-selection of possible sound decay
	int i;
	double RTml = -1;
	int len = par->N_sub;
	double* seg = (double*)calloc(len, sizeof(double));
	memcpy(seg, frame, len * sizeof(double));
	/*calculate variance, minimum and maximum of first sub-frame*/
	double var_pre = var(seg, len);
	double min_pre = minv(seg, len).val;
	double max_pre = maxv(seg, len).val;
	double var_cur = 0, max_cur = 0, min_cur = 0;
	int nos_max = par->nos_max;
	int nos_min = par->nos_min;
	int no_bins = par->no_bins;
	int k;
	double* ll = NULL;
	for (k = 1; k < nos_max; ++k) {
		memcpy(seg, &frame[k * len], len * sizeof(double));
		/*calculate variance, minimum and maximum of succeding sub-frame*/
		var_cur = var(seg, len);
		max_cur = maxv(seg, len).val;
		min_cur = minv(seg, len).val;
		/*
		Pre-Selection of suitable speech decays
		if variance, maximum decraease and minimum increase
		=> possible sound decay detected
		*/
		if (var_pre > var_cur && max_pre > max_cur && min_pre < min_cur) {
			cnt++;
			/*current values becomes previous values*/
			var_pre = var_cur;
			max_pre = max_cur;
			min_pre = min_cur;
		}
		else {
			/*minimum length for assumed sound decay achieved?*/
			if (cnt >= nos_min) {
				/*Maximum Likelihood (ML) Estimation of the RT60*/
				RTml = max_loglf(frame, cnt * len, par->a, par->Tquant, par->La, ll);
			}
			break;
		}
		/*maximum frame length achieved?*/
		if (k == nos_max - 1) {
			RTml = max_loglf(frame, cnt * len, par->a, par->Tquant, par->La, ll);
		}
	}
	/*eof sub-frame loop*/

	/*new ML estimate calculated*/
	int index = -1;
	if (RTml >= 0) {
		/*apply order statistics to reduce outliers*/
		par->hist_counter += 1;
		for (i = 0; i < no_bins; ++i) {
			/*find index corresponding to the ML estimate*/
			if (RTml >= par->hist_limits[i] && RTml <= par->hist_limits[i + 1]) {
				index = i;
				break;
			}
		}
		/*update histogram with ML estimates for the RT60*/
		par->hist_rt[index] += 1;
		/*remove old values from histogram*/
		if (par->hist_counter > par->buffer_size + 1)
			par->hist_rt[par->buffer[0]] -= 1;
		memcpy(par->buffer, &par->buffer[1], (par->buffer_size * sizeof(int)) - 1);
		/*update buffer with indices*/
		par->buffer[par->buffer_size - 1] = index;
		/*find index for maximum of the histogram*/
		int idx = maxv(par->hist_rt, par->no_bins).ind;
		/*map index to RT60 value*/
		par->RT_raw = par->Tquant[idx];
	}
	/*final RT60 estimate obtained by recursive smoothing*/
	double RT = par->alpha * par->RT_last + (1 - par->alpha) * par->RT_raw;
	par->RT_last = RT;
	/*intermediate ML estimate for later analysis*/
	par->RTml = RTml;

	if (ll)
		free(ll);
	if (seg)
		free(seg);
	return RT;
}