/*
 * MIP.h
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#ifndef SRC_ADVERSARIAL_MILP_H_
#define SRC_ADVERSARIAL_MILP_H_

#include "Adversarial.h"
#include <ilcplex/cplex.h>
#include <stdlib.h>

class MILP: public Adversarial {
protected:

	double MAX;


	long var;
	int* theta;
	int* beta;
	int* s;
	int* pi;
	int* u;
	int* alpha;
	int* xi;
	int* xi2;
	double* retX;

	void addXiVar(int i);

	void addXi2Var(int i);

	void addTHETAVar(int i);

	void addUVar(int i);

	void addAlphaVar(int i);

	void addBetaVar(int i);

	void addSVar(int i);

	void addPIVar(int i);

	void addThetaRest(int i);

	void addPiRest(int i);

	void addSRest(int i);

	void addURest(int i);

	void addGammaRest();

	CPXENVptr env;
	CPXLPptr lp;


public:
	MILP(Master* mst, int n, int m, int gamma, int* dv, double* dm, double* h, double* p, double* q);

	virtual ~MILP();

	double solve();
};
#endif /* SRC_ADVERSARIAL_MILP_H_ */
