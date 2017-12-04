/*
 * Master.h
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#ifndef SRC_ADVERSARIAL_MASTER_H_
#define SRC_ADVERSARIAL_MASTER_H_

#include <ilcplex/cplex.h>
#include <stdlib.h>

#define MIN(a,b) a > b? b: a
#define POS(a,b) ((a-b) < 0) ? 0 : a-b

class Master {
protected:

	int n;
	int m;
	double MAX;
	double* c;
	double* h;
	double* p;
	double* q;
	double* cap;
	double* dm;
	double* dv;


	long var;
	int theta;
	int* x;
	int* s;
	int* r;
	int* pi;

	int* u;
	int* alpha;

	double* retX;


	double* dem;

	CPXENVptr env;
	CPXLPptr lp;

	void addXVar(int i);

	void addTHETAVar();

	void addRVar(int i);

	void addSVar(int i);

	void addPIVar(int i);

	void addLinVars(int i);

	void addPiRest(int i);

	void addRRest(int i);

	void addSRest(int i);

	void addFULLRest();

	void addPIRest1(int i);

	void addPIRest2(int i);

	void addPIRest3(int i);

	void addURest(int i);

	double bestval;

public:
	Master(int n, int m, double* c, double* h, double* p, double* q, double* cap,  double* dm, int* dv, int gamma);

	double solve(double* solution);//retorna a melhor solu��o para o adversario e solution eh o valor da solucao total

	double getX(int i);

	double getDv(int i);

	double getDm(int i);

	int getM();

	void setDem(int i, double perc);

	void addXiRest();

	virtual ~Master();
};


#endif /* SRC_ADVERSARIAL_MASTER_H_ */
