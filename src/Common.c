/****************************************************************************
Calculate variance, minimum, maximum; removal DC and normalization function
****************************************************************************/
#include "Common.h"

VALIND minv(double* x, int n) {
	VALIND xm = { INFINITY };
	for (int i = 0; i < n; ++i) {
		if (x[i] < xm.val) {
			xm.val = x[i];
			xm.ind = i;
		}
	}
	return xm;
}

VALIND maxv(double* x, int n) {
	VALIND xm = { -INFINITY };
	for (int i = 0; i < n; ++i) {
		if (x[i] > xm.val) {
			xm.val = x[i];
			xm.ind = i;
		}
	}
	return xm;
}

double mean(double* x, int n) {
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		sum += x[i];
	}
	return sum / n;
}

double var(double* x, int n) {
	double xm = mean(x, n);
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		double xtmp = fabs(x[i] - xm);
		sum += xtmp * xtmp;
	}
	return sum / n;
}

double* preprocess(short* x_s, int n) {
	double* x = (double*)calloc(n, sizeof(double));
	double sum = 0, maxv = 0;
	int i;
	for (i = 0; i < n; ++i) {
		if (x == NULL) {
			fprintf(stderr, "Unable to allocate memory!\n");
			return 1;
		}
		x[i] = x_s[i] / 32767.0;
		sum += x[i];
		double tmpabs = fabs(x[i]);
		if (maxv < tmpabs) {
			maxv = tmpabs;
		}
	}
	double meanv = sum / n;
	for (i = 0; i < n; ++i) {
		x[i] -= meanv;
		x[i] /= maxv;
	}
	return x;
}