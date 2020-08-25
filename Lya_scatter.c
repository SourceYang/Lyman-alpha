#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "random.h"

#define C 29979245800 // speed of light, in [cm s^-1]
#define M_HI 1.6735575e-24 // mass of Hydrogen, in [gram]
#define M_p 1.67262193269e-24 // mass of proton, in [gram]
#define k 1.380649e-16 // Boltzmann constant, in [g cm^2 s^-2 K^-1]
#define M_PI 3.14159265358979323846 // pi
#define T_2p 1.6e-9 // lifetime of 1s to 2p, in [s]
#define v_center 2.47e15 // the center frequency of Lyman alpha emission, in [Hz]
#define SCATTER_SCPAR 10.0 // adjustable constant for scattering process

/* calculate the new direction, polarization and frequency of the reemitted Lyman alpha photon */
void scatter(double THI, double *L_mu, double *L_p, double *L_v, long *idum){
	float r, r1, r2; // random numbers

	/* calculate the Voigt distribution */
	double x;
	double lambda_i;
	double lambda_d;
	double a;
	double c;
	double H_x;

	/* calculate the probability of each area */
	double u_1a, u_1b, u_1c, u_1d;
	double h_1, h_2, h_3, h_4, h_5;
	double u_1e, sign_x;
	double w1, w2, w3;
	double Area_a, Area_b, Area_c;

	/* calculate the scattering frequency */
	double u1, u2;
	double R1;
	double p_A, p_B, p_C;

	/* choose scattering angle and calculate the polarization */
	double u, mu; // cos gamma and cos theta
	double gamma;
	double alpha;
	double theta, theta_; // angle absorbed and angle reemitted
	double beta, zeta;
	double x1, y1;
	double E_1;
	double p_NSNS, p_NSEW, p_EWNS, p_EWEW, p_NS, p_EW, p_n;
	double R;

	/* This formula is derived from https://academic.oup.com/mnras/article/369/4/2025/1096494 */
	lambda_i = C / v_center;
	lambda_d = sqrt(2 * k * THI / M_p) * lambda_i / C;
	x = (C / (*L_v + v_center) - lambda_i) / lambda_d;
	a = lambda_i * lambda_i / (4 * M_PI * C * lambda_d * T_2p);
	c = a / sqrt(M_PI);
	if(fabs(x) < 1e-4){
		H_x = (1 - 2 * c) + (-1 + 4 * c) * x * x + (1 / 2. - 5 * c / 3) * pow(x, 4) + (-1 / 6. - 14 * c / 15) * pow(x, 6);
	}
	else{
		H_x = exp(-x * x) - c * exp(-2 * x * x) * (4 * x * x + 3) * (x * x + 1) / (x * x) + c / 2 * (2 * x * x + 3) / pow(x, 4) - c / 2 * (2 * x * x + 3) * exp(-2 * x * x) / pow(x, 4);
	}

	u_1a = -7;
	u_1b = x - SCATTER_SCPAR * a;
	u_1c = x + SCATTER_SCPAR * a;
	u_1d = 7;

	if(u_1b < u_1a){
		u_1b = u_1a;
	}
	if(u_1b > u_1d){
		u_1b = u_1d;
	}
	if(u_1c < u_1a){
		u_1c = u_1a;
	}
	if(u_1c > u_1d){
		u_1c = u_1d;
	}

	h_1 = 1 / (a * a + (x - u_1a) * (x - u_1a));
	h_2 = 1 / (a * a + (x - u_1b) * (x - u_1b));
	h_3 = 1 / (a * a + (x - u_1c) * (x - u_1c));
	h_4 = 1 / (a * a + (x - u_1d) * (x - u_1d));

	w1 = h_1;
	if(h_2 > w1){
		w1 = h_2;
	}
	w3 = h_3;
	if(h_4 > w3){
		w3 = h_4;
	}

	if(u_1b * u_1c < 0){
		w2 = 1 / (a * a);
	}
	else{
		w2 = exp(-u_1b * u_1b);
		if(exp(-u_1c * u_1c) > w2){
			w2 = exp(-u_1c * u_1c);
		}
		w2 = w2 / (a * a);
	}

	if(fabs(x) < 2){
		sign_x = x > 0? 1: -1;
		u_1e = (x - sign_x * sqrt(x * x + 4)) / 2;
		h_5 = 1 / (a * a + (x - u_1a) * (x - u_1e));
		if(u_1e >= u_1a && u_1e < u_1b){
			if(h_5 > w1){
				w1 = h_5;
			}
		}
		else if(u_1e >= u_1b && u_1e < u_1c){
			if(h_5 > w2){
				w2 = h_5;
			}
		}
		else if(u_1e >= u_1c && u_1e <= u_1d){
			if(h_5 > w3){
				w3 = h_5;
			}
		}
	}

	Area_a = (u_1b - u_1a) * w1;
	Area_b = (u_1c - u_1b) * w2;
	Area_c = (u_1d - u_1c) * w3;

	p_A = Area_a / (Area_a + Area_b + Area_c);
	p_B = Area_b / (Area_a + Area_b + Area_c);
	p_C = Area_c / (Area_a + Area_b + Area_c);

	if(fabs(x) < 6.5){
		do{
			r = ran2(idum);
			if(r < p_A){
				u1 = ran2(idum) * (u_1b - u_1a) + u_1a;
				r1 = ran2(idum) * w1;
			}
			else if(r >= p_A && r < p_A + p_B){
				u1 = ran2(idum) * (u_1c - u_1b) + u_1b;
				r1 = ran2(idum) * w2;
			}
			else if(r >= p_A + p_B && r <= 1){
				u1 = ran2(idum) * (u_1d - u_1c) + u_1c;
				r1 = ran2(idum) * w3;
			}
			R1 = exp(-u1 * u1) / (a * a + (x - u1) * (x - u1));
		}while(r1 > R1);
	}
	else{
		r1 = ran2(idum);
		r2 = ran2(idum);
		u1 = 1 / x + sqrt(-log(r1)) * cos(2 * M_PI * r2);
	}

	do{
		mu = *L_mu;
		theta = acos(mu);
		u = 2 * ran2(idum) - 1;
		gamma = acos(u);
		alpha = 2 * M_PI * ran2(idum);
		mu = mu * u + sqrt(1 - mu * mu) * sin(gamma) * cos(alpha);
		theta_ = acos(mu);

		y1 = sin(theta) * sin(alpha);
		x1 = sin(theta) * cos(alpha) * cos(gamma) - cos(theta) * sin(gamma);
		beta = atan2(y1, x1);

		y1 = sin(alpha) * sin(gamma);
		x1 = sin(theta) * cos(gamma) - cos(theta) * cos(alpha) * sin(gamma);
		zeta = atan2(y1, x1);

		if(abs(x - u1) < sqrt(3) * a){
			E_1 = 1;
		}
		else{
			E_1 = 0;
		}

		p_NSNS = sin(theta_) * sin(theta) + cos(theta_) * cos(theta) * cos(zeta);
		p_NSEW = -cos(theta_) * sin(zeta);
		p_EWNS = cos(theta) * sin(zeta);
		p_EWEW = cos(zeta);

		p_NS = (1 - E_1) / (8 * M_PI) + (3 / (16 * M_PI)) * E_1 * ((1 + *L_p) * p_NSNS * p_NSNS) + (1 - *L_p) * p_NSEW * p_NSEW;
		p_EW = (1 - E_1) / (8 * M_PI) + (3 / (16 * M_PI)) * E_1 * ((1 + *L_p) * p_EWNS * p_EWNS) + (1 - *L_p) * p_EWEW * p_EWEW;

		p_n = p_NS + p_EW;

		r = ran2(idum);
		R = (8 * M_PI / 3) * p_n;
	}while(r >= R);

	*L_mu = mu;
	*L_p = ((p_NS - p_EW) / (p_NS + p_EW));

	r1 = ran2(idum);
	r2 = ran2(idum);
	u2 = sqrt(-log(r1)) * cos(2 * M_PI * r2);

	x = x + u1 * cos(gamma) + u2 * sin(gamma) - u1;
	*L_v = C / (lambda_i + x * lambda_d) - v_center;
}

int main(int argc, char **argv){
	int i; // index
	double THI; // setting for the background temperature of hydrogen atom
	double L_mu, L_p, L_v; // settings for direction, polarization and frequency of a certain Lyman alpha photon
	double mu, p, v; // the direction, polarization and frequency of each photon after scattering
	int n_photon; // number of scattering photons with the same setting

	long seed;

	sscanf(argv[1], "%lf", &THI);
	sscanf(argv[2], "%lf", &L_mu);
	sscanf(argv[3], "%lf", &L_p);
	sscanf(argv[4], "%le", &L_v);
	sscanf(argv[5], "%d", &n_photon);
	sscanf(argv[6], "%ld", &seed);

	long *idum = &seed;

	printf("direction polarization frequency\n");
	printf("%lf %lf %le \n", L_mu, L_p, L_v);
	
	printf("direction polarization frequency\n");
	for(int i = 0; i < n_photon; i++){
		mu = L_mu;
		p = L_p;
		v = L_v;
		scatter(THI, &mu, &p, &v, idum);
		printf("%lf %lf %le \n", mu, p, v);
	}

	return 0;
}
