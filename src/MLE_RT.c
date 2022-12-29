/****************************************************************************
Blind RT60 estimation by means of a maximum-liklihood (ML) estimator
provides RT60 estimates over time for the input sequence x using a frame-wise processing scheme
****************************************************************************/
#include "Common.h"

double MLE_RT(double* x, int n, PSIMPAR psimpar, double** rt_est, int* rt_est_len,
	double** RT_est, int* Rt_est_len, PRTE_HANDLE rte_handle) {
	int num_samples = n;
	int block_size = psimpar->block_size;
	int overlap = psimpar->overlap;
	int fs = psimpar->fs;
	double* remove_from_avg = psimpar->remove_from_avg;
	/*Initialize RT60 estimation*/
	MLE_RT_init(fs, rte_handle);
	/*vectors for all RT60 estimates*/
	int len = (int)ceil(n / rte_handle->N_shift) + 1;
	double* rt_est_tmp = (double*)calloc(len, sizeof(double));
	int i, j;
	double RT = rte_handle->RT_initial;
	for (i = 0; i < len; ++i) {
		rt_est_tmp[i] = RT;
	}
	len = (int)ceil(n / overlap) + 1;
	double* RT_est_tmp = (double*)calloc(len, sizeof(double));
	/*initialize counters*/
	int k = 0;                                                               //frame counter for overall block processing scheme
	int rt_frame_cnt = 0;                                                    //frame counter for RT60 estimation
	int k_rt = rte_handle->N * rte_handle->down;                             //index counter for RT60 estimation
	int cnt;
	int ind0 = 1;
	int down = rte_handle->down;
	len = (int)floor((k_rt - ind0) / down) + 1;
	double* x_seg = (double*)calloc(len, sizeof(double));
	for (cnt = 1; cnt <= num_samples - block_size + 1; cnt += overlap) {
		/*New RT60 Estimation*/
		if (cnt > k_rt) {
			j = 0;
			for (i = ind0; i <= k_rt; i += down) {
				x_seg[j++] = x[i];
			}
			RT = MLE_RT_frame(x_seg, rte_handle);                            //actual RT60 estimation
			k_rt += rte_handle->N_shift;                                     //increase index counter for RT60 estimation
			ind0 = k_rt - rte_handle->N * rte_handle->down + 1;
			rt_est_tmp[rt_frame_cnt++] = RT;                                 //save RT60 estimate over time
		}
		RT_est_tmp[k++] = RT;
	}
	*RT_est = (double*)malloc(k * sizeof(double));
	memcpy(*RT_est, RT_est_tmp, k * sizeof(double));
	*Rt_est_len = k;
	*rt_est = (double*)malloc(rt_frame_cnt * sizeof(double));
	memcpy(*rt_est, rt_est_tmp, rt_frame_cnt * sizeof(double));
	*rt_est_len = rt_frame_cnt;
	/*Mean RT60, averaged over all considered frames*/
	double* fr2sec_idx = NULL;
	fr2sec_idx = linspace(1, num_samples / fs, rt_frame_cnt);
	double rfatmp1 = remove_from_avg[0];
	double rfatmp2 = fr2sec_idx[rt_frame_cnt - 1] - remove_from_avg[1];
	int idx[2] = { 0 };
	for (i = 0; i < rt_frame_cnt; ++i) {
		if (fr2sec_idx[i] > rfatmp1) {
			idx[0] = i;
			break;
		}
	}
	for (i = rt_frame_cnt - 1; i >= 0; --i) {
		if (fr2sec_idx[i] < rfatmp2) {
			idx[1] = i;
			break;
		}
	}
	double rt_est_mean = mean((*rt_est + idx[0]), idx[1] - idx[0] + 1);

	if (x_seg)
		free(x_seg);
	if (RT_est_tmp)
		free(RT_est_tmp);
	if (rt_est_tmp)
		free(rt_est_tmp);
	if (fr2sec_idx)
		free(fr2sec_idx);
	return rt_est_mean;
}