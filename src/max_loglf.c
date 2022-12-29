/****************************************************************************
Maximum-Likelihood (ML) function
****************************************************************************/
#include "Common.h"

double max_loglf(double* h, int N, double* a, double* Tquant, int len, double* ll) {
	ll = (double*)calloc(len, sizeof(double));
	double* h_square = (double*)calloc(N, sizeof(double));
	int i, j;
	for (i = 0; i < N; ++i) {
		if (h_square == NULL) {
			fprintf(stderr, "Unable to allocate memory!\n");
			return 1;
		}
		h_square[i] = h[i] * h[i];
	}
	for (i = 0; i < len; ++i) {
		double Sum = 0;
		for (j = 0; j < N; ++j) {
			Sum += pow(a[i], -2.0 * j) * h_square[j];
		}
		if (ll == NULL) {
			fprintf(stderr, "Unable to allocate memory!\n");
			return 1;
		}
		if (Sum < 1e-12)
			ll[i] = -INFINITY;
		else
			ll[i] = -N / 2.0 * (((double)N - 1) * log(a[i]) + log(2.0 * pi / N * Sum) + 1);
	}
	/*maximum of the log-likelihood function*/
	int idx = maxv(ll, len).ind;
	if (h_square)
		free(h_square);
	/*corresponding ML estimate for the RT60*/
	return Tquant[idx];
}