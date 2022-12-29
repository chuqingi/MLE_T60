/****************************************************************************
Generates n points between d1 and d2
****************************************************************************/
#include "Common.h"

double* linspace(double d1, double d2, double n) {
	double* y = (double*)calloc(n, sizeof(double));
	double n1 = n - 1;
	for (int i = 0; i <= n1; ++i) {
		if (y == NULL) {
			fprintf(stderr, "Unable to allocate memory!\n");
			return 1;
		}
		y[i] = d1 + i * (d2 - d1) / n1;
	}
	return y;
}