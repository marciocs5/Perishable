/*
 * AffDecRules.h
 *
 *  Created on: 18 de set de 2017
 *      Author: marcio
 */

#ifndef SRC_LBS_AFFDECRULES_H_
#define SRC_LBS_AFFDECRULES_H_

#include <ilcplex/cplex.h>
#include <stdlib.h>
#include "../main/Timer.h"

class AffDecRules {
private:
	CPXENVptr env;
	CPXLPptr lp;

	//problem data
	int n;
	int m;
	int cons_n;//number of constraints
	double MAX;
	double* c;
	double* h;
	double* p;
	double* q;
	double* cap;
	double* dm;
	double* dv;

	Timer* timer;
	double solution;

	long var;
	int theta;
	int* x_main;
	int** x;
	int* y;

	int*** x_affine_pos;
	int*** x_affine_neg;
	int** y_affine_pos;
	int** y_affine_neg;

	int** dualvarspos;//the dual variables
	int** dualvarsneg;//the dual variables
	int* specialdualvars;//the dual variables that are multipleid by gamma

	double** A;//restriction matrix
	double** B;//the idependent terms for each xi
	double** w;

	int gamma;

	void addThetaVar();

	void addXVar();

	void addXMainVar();

	void addYVar();

	void addDualVarPos();

	void addDualVarNeg();

	void addDualVarSpe();

	void addX_affineVarPos();

	void addX_affineVarNeg();

	void addY_affineVarPos();

	void addY_affineVarNeg();

	void addPrimalCons(int i);

	void addDualConsPos(int i, int j);

	void addDualConsNeg(int i, int j);

public:
	AffDecRules(int n, int m, double* c, double* h, double* p, double* q, double* cap,  double* dm, int* dv, int gamma);

	virtual ~AffDecRules();

	void solve();

	double getSolution();

	float getElapsedTime();

	long getVarNumber();

};

#endif /* SRC_LBS_AFFDECRULES_H_ */
