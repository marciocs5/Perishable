#include "Master.h"

#ifndef SRC_ADVERSARIAL
#define SRC_ADVERSARIAL

class Adversarial {
protected:

	int n;
	int m;
	int gamma;
	int* dv;
	double* dm;
	double* h;
	double* p;
	double* q;
	Master* mst;

public:
	Adversarial(Master* mst, int n, int m, int gamma, int* dv, double* dm, double* h, double* p, double* q);

	virtual ~Adversarial();

	virtual double solve()=0;

	int getGamma();

};

#endif /* SRC_ADVERSARIAL */
