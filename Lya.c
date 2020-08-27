#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "random.h"
#include "functions.h"

#define C 29979245800 // speed of light, in [cm s^-1]
#define M_HI 1.6735575e-24 // mass of Hydrogen, in [gram]
#define M_p 1.67262193269e-24 // mass of proton, in [gram]
#define k 1.380649e-16 // Boltzmann constant, in [g cm^2 s^-2 K^-1]
#define M_PI 3.14159265358979323846 // pi
#define T_2p 1.6e-9 // lifetime of 1s to 2p, in [s]
#define v_center 2.47e15 // the center frequency of Lyman alpha emission, in [Hz]
#define A 6.25e8 // spontaneous decay rate, in [s^-1]
#define DNHI 2.5e16 // column density of Hydrogen atoms, in [cm^-2]
#define H 3.3149362e-17 // Hubble constant in the period of z = 8, in [s^-1]
#define Z 8 // redshift
#define EPS_Z 1e-9 // epsilon z
#define SCATTER_SCPAR 10.0 // adjustable constant for scattering process

#define NGRID 2000 // number of grid
#define N_NU 128 // number of frequency bins
#define N_PHOTON_MAX 1048576 // occupy the max space

/* get a random position of Lyman alpha photon, in grid */
int get_position(double *y1H, double *y1He, double *dy1H, double *Te, double y1H1, long *idum){
	double L_emissivity[NGRID]; // Lyman alpha emssivity
	double T_emissivity = 0; // total Lyman alpha emssivity
	double P_emissivity[NGRID]; // probability of emssivity in each grid
	double Cum_emissivity = 0; // cumulative probability of emssivity
	double q[NGRID]; // reaction coefficient
	int z; // position
	float b;  //random number
	int i; // index

	for(i = 0; i < NGRID; i++){
		q[i] = get_cooling_rate(Te[i], y1H[i], y1He[i]);
		L_emissivity[i] = dy1H[i] * y1H1 * y1H1 * q[i];
		T_emissivity += L_emissivity[i];
	}

	for(i = 0; i < NGRID; i++){
		P_emissivity[i] = L_emissivity[i] / T_emissivity;
	}

	b = ran2(idum);
	i = 0;
	do{
		Cum_emissivity += P_emissivity[i];
		i++;
	} while(b > Cum_emissivity && i <= NGRID);
	z = i - 1;
	i = 0;

	return z;
}

/* get a random frequency of Lyman alpha photon, in delta Hz away from the center */
double get_frequency(double THI, long *idum){
	double sigma, r, theta;  // Voigt Gaussian
	double gamma;  // Voigt Lorentz
	double x;  // variable in Voigt Gaussian
	double y;  // variable in Voigt Lorentz
	double dv;  // frequency away from the center
	float b;  // random number

	sigma = v_center / C * sqrt(k * THI / M_HI);
	b = ran2(idum);
	theta = 2 * M_PI * b;
	b = ran2(idum);
	r = sqrt(-2 * log(1 - b));
	x = r * cos(theta) * sigma;

	gamma = 1 / (4 * M_PI * T_2p);
	b = ran2(idum);
	y = gamma * tan((2 * M_PI * b - M_PI) / 2);

	dv = x + y;

	return dv;
}

/* propagate a Lyman alpha photon until it's absorbed */
double propagate(double dv, double mu, double THI, double y1H1, int U, long *idum){
	/* calculate the Voigt distribution */
	double x;
	double lambda_i;
	double lambda_d;
	double a;
	double c;
	double H_x;

	double sigma_tot; // cross section
	double l_mfp; // mean free path
	float b; // random number
	double s; // distance moved
	double dz; // grid moved relative to the front
	double length; // the length of one grid

	/* This formula is derived from https://academic.oup.com/mnras/article/369/4/2025/1096494 */
	lambda_i = C / v_center;
	lambda_d = sqrt(2 * k * THI / M_p) * lambda_i / C;
	x = (C / (dv + v_center) - lambda_i) / lambda_d;
	a = lambda_i * lambda_i / (4 * M_PI * C * lambda_d * T_2p);
	c = a / sqrt(M_PI);
	if(fabs(x) < 1e-4){
		H_x = (1 - 2 * c) + (-1 + 4 * c) * x * x + (1 / 2. - 5 * c / 3) * pow(x, 4) + (-1 / 6. - 14 * c / 15) * pow(x, 6);
	}
	else{
		H_x = exp(-x * x) - c * exp(-2 * x * x) * (4 * x * x + 3) * (x * x + 1) / (x * x) + c / 2 * (2 * x * x + 3) / pow(x, 4) - c / 2 * (2 * x * x + 3) * exp(-2 * x * x) / pow(x, 4);
	}

	sigma_tot = C / (sqrt(2 * k * THI / M_p) * v_center) / sqrt(M_PI) * 3 * A * C * C * H_x / (8 * M_PI * (dv + v_center) * (dv + v_center));

	l_mfp = 1 / (y1H1 * sigma_tot);
	b = ran2(idum);
	s = l_mfp * log(1 / b);
	dz = s * mu - U * s / C;
	length = DNHI / y1H1;

	return (dz / length);
}

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
	long j, istep; // index
	int n_photon; // number of photons simulated
	double fracflux[N_NU], nu[N_NU], sigH[N_NU], sigHe[N_NU]; // flux and cross section in each bin
	double y1H[NGRID], y1He[NGRID]; // neutral fractions of electron and Helium
	double dy1H[NGRID], dy1He[NGRID]; // delta neutral fractions
	double D_NHI[NGRID]; // delta N_HI
	double dEH[NGRID], EH[NGRID], EHI[NGRID], EHII[NGRID], EHeI[NGRID], EHeII[NGRID]; // delta energy and energy
	double Te[NGRID], THI[NGRID], THII[NGRID], THeI[NGRID], THeII[NGRID]; // temperature in each grid
	double s, dz; // grid moved real distance and delta grid after propagating
	double L_position[N_PHOTON_MAX], L_mu[N_PHOTON_MAX], L_p[N_PHOTON_MAX], L_v[N_PHOTON_MAX]; // the position, direction, polarization and frequency of Lyman alpha photons
	double p_i, p_f; // initial position and final position in one propagating process
	int z; // the position of one Lyman alpha photon
	float b; //random number

	/* using in read data file */
	FILE *fp;
	char* filename = "output.txt";

	int T, U; // blackbody incident temperature and ionization front velocity, will wirte in command line
	double y1H1; // Hydrogen density
	long seed;

	sscanf(argv[1], "%d", &T);
	sscanf(argv[2], "%d", &U);
	sscanf(argv[3], "%lE", &y1H1);
	sscanf(argv[4], "%d", &n_photon);
	sscanf(argv[5], "%ld", &seed);

	long *idum = &seed;

	/* read in the data file */
	fp = fopen(filename, "r");
	if(fp == NULL){
		printf("Cannot read the data file");
		return 1;
	}
	while(fgetc(fp) != 93){
		fgetc(fp);
	}
	while(fgetc(fp) != 100){
			fgetc(fp);
	}
	while(fgetc(fp) != 93){
		fgetc(fp);
	}
	for(j = 0; j < NGRID; j++){
		fscanf(fp, "%d %lE %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &j, &D_NHI[j], &y1H[j], &y1He[j], &EH[j], &Te[j], &EHII[j], &THII[j], &EHI[j], &THI[j], &EHeI[j], &THeI[j], &EHeII[j], &THeII[j], &dEH[j]);
	}
	fclose(fp);

	set_bb(fracflux,T);
	get_ion_rate(y1H, y1He, fracflux, dy1H, dy1He, dEH, Te, istep, U);

	for(int i = 0; i < n_photon; i++){
		L_position[i] = get_position(y1H, y1He, dy1H, Te, y1H1, idum);
		z = (int)L_position[i];
		L_mu[i] = 2 * ran2(idum) - 1;
		L_p[i] = 0;
		L_v[i] = get_frequency(THI[z], idum);

		while((L_position[i] >= 0 && L_position[i] < NGRID) && fabs(L_v[i]) <= 1e13){
			p_i = L_position[i]; // initial position before propagating
			dz = propagate(L_v[i], L_mu[i], THI[z], y1H1, U, idum);
			L_position[i] += dz;
			while((L_position[i] < z || L_position[i] > z + 1) && (L_position[i] >= 0 && L_position[i] < NGRID) && fabs(L_v[i]) <= 1e13){
				if(L_position[i] < z){
					L_position[i] = z - EPS_Z;
					z--;
					s = ((L_position[i] - p_i) * DNHI / y1H1) / (L_mu[i] - U / C);
					L_v[i] -= H * s * (v_center + L_v[i]) / C;
					p_i = L_position[i];
					dz = propagate(L_v[i], L_mu[i], THI[z], y1H1, U, idum);
					L_position[i] += dz;
				}
				else if(L_position[i] > z + 1){
					L_position[i] = z + 1 + EPS_Z;
					z++;
					s = ((L_position[i] - p_i) * DNHI / y1H1) / (L_mu[i] - U / C);
					L_v[i] -= H * s * (v_center + L_v[i]) / C;
					p_i = L_position[i];
					dz = propagate(L_v[i], L_mu[i], THI[z], y1H1, U, idum);
					L_position[i] += dz;
				}
			}
			z = (int)L_position[i];
			p_f = L_position[i]; // final position after propagating
			dz = p_f - p_i;
			s = (dz * DNHI / y1H1) / (L_mu[i] - U / C);
			L_v[i] -= H * s * (v_center + L_v[i]) / C; // Hubble expansion
			fflush(stdout);
			scatter(THI[z], &L_mu[i], &L_p[i], &L_v[i], idum);
			fflush(stdout);
		}
	}

	/* display */
	printf("i, position, direction, polarization, frequency\n");
	for(int i = 0; i < n_photon; i++){
		printf("%d %lf %lf %lf %le \n", i, L_position[i], L_mu[i], L_p[i], L_v[i]);
	}

	return 0;
}