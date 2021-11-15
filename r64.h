#include <stdlib.h>
#include <math.h>
#include <stdio.h>

unsigned long long int int64(unsigned long long int j){
	static unsigned long long int v = 4101842887655102017LL;
	static unsigned long long int w = 1;

	v ^= j;

	v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
	w = 4294957665U*(w & 0xffffffff) + (w >> 32);

	return v ^ w;
}

double ran2(unsigned long long int j){
	return 5.42101086242752217E-20 * int64(j);
}