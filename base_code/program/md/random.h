#pragma once

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

class Random {
public:
	long     iy;
	long     iv[NTAB];
	long     idum[1];

	Random(long seed);
	double nextDouble();
    double nextGauss();
};
