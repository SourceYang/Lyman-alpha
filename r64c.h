#include <stdlib.h>
#include <math.h>
#include <stdio.h>

unsigned long long int int64(unsigned long long int j){
	static unsigned long long int v = 4101842887655102017LL;
	static unsigned long long int w = 1;
	static unsigned long long int u;

	u = j ^ v;
	v = u;
	w = v;

	u = u * 2862933555777941757LL + 7046029254386353087LL;
	v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
	w = 4294957665U*(w & 0xffffffff) + (w >> 32);
	unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
	return (x + v) ^ w;
}

double ran2(unsigned long long int j){
	return 5.42101086242752217E-20 * int64(j);
}