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

#define NGRID 2000 // number of grid
#define N_NU 128 // number of frequency bins

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
	double dz; // grid moved after propagating
	double mu; // the direction of Lyman alpha photons
	int z; // the position of one Lyman alpha photon
	float b; //random number

	/* using in read data file */
	FILE *fp;
	char* filename = "output.txt";

	int T, U; // blackbody incident temperature and ionization front velocity, will wirte in command line
	double y1H1, THI_set; // Hydrogen density and setting for the background temperature of hydrogen atom
	double L_v; // settings for Lyman-alpha frequency
	long seed;

	sscanf(argv[1], "%d", &T);
	sscanf(argv[2], "%d", &U);
	sscanf(argv[3], "%lE", &y1H1);
	sscanf(argv[4], "%lf", &THI_set);
	sscanf(argv[5], "%le", &L_v);
	sscanf(argv[6], "%d", &n_photon);
	sscanf(argv[7], "%ld", &seed);

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
		mu = 2 * ran2(idum) - 1;
		dz = propagate(L_v, mu, THI_set, y1H1, U, idum);
		printf("%lf \n", dz);
	}

	return 0;
}
