#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "random.h"
#include "functions.h"

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
	double L_position[N_PHOTON_MAX], L_mu[N_PHOTON_MAX], L_p[N_PHOTON_MAX], L_v[N_PHOTON_MAX]; // the position, direction, polarization and frequency of Lyman alpha photons
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
		z = get_position(y1H, y1He, dy1H, Te, y1H1, idum);
		printf("%d \n", z);
	}

	return 0;
}