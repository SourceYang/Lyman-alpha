#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "random.h"

#define v_center 2.47e15 // the center frequency of Lyman alpha emission, in [Hz]
#define C 29979245800 // speed of light, in [cm s^-1]
#define k 1.380649e-16 // Boltzmann constant, in [g cm^2 s^-2 K^-1]
#define M_HI 1.6735575e-24 // mass of Hydrogen, in [gram]
#define M_PI 3.14159265358979323846 // pi
#define T_2p 1.6e-9 // lifetime of 1s to 2p, in [s]

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

int main(int argc, char **argv){
	double dv;  // frequency away from the center

	double THI; // background temperature of hydrogen atom
	int n_photon; // number of generating
	long seed; // random seed

	sscanf(argv[1], "%d", &THI);
	sscanf(argv[2], "%d", &n_photon);
	sscanf(argv[3], "%ld", &seed);

	long *idum = &seed;

	for(int i = 0; i < n_photon; i++){
		dv = get_frequency(THI, idum);
		printf("%le \n", dv);
	}
}