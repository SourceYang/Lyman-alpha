#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define M_PI 3.14159265358979323846 // pi

double F(double v){
	double F;
	int imax;
	double temp;
	double f;

	if(fabs(v) <= 8){
		F = v + v * v * v / 3.;
		imax = (int)(v * v) + (int)(fabs(v) * 10);
		temp = v * v;
		f = v * v * v / 3.;
		for(int i = 2; i < imax; i++){
			f *= temp * (2 * i - 1) / i / (2 * i + 1);
			F += f;
		}
		F = exp(-v * v) * F;
	}
	else{
		F = 1 / v;
		temp = 1 / (v * v);
		f = F;
		for(int i = 1; i < 15; i++){
			f *= temp * (2 * i - 1) / 2;
			F += f;
		}

		F = 0.5 * F;
	}

	return F;
}

double H_xa(double x, double a){
	double H_0, H_1, H_2, H_x;

	H_0 = exp(-x * x);
	H_1 = -2 / sqrt(M_PI) * (1 - 2 * x * F(x));
	H_2 = (1 - 2 * x * x) * exp(-x * x);

	H_x = H_0 + a * H_1 + a * a * H_2;

	return H_x;
}
