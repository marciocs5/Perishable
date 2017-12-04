/*
 * HeurMILP1.cpp
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#include "HeurMILP2.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>

HeurMILP2::HeurMILP2(Master* mst, int n, int m, int gamma, int* dv, double* dm, double* h, double* p, double* q) : MILP(mst, n, m, gamma, dv, dm, h, p, q){
}

HeurMILP2::~HeurMILP2(){
}

double HeurMILP2::solve(){
	int status;
	int indices[3*n];
	double values[3*n];
	int index = 0;
	double aux_val;

	for(int i = 0; i < n; i++){

		aux_val = 0.0;
		for(int k = 0; k <= i ; k++){
			aux_val -= dm[k];
			if((k - m) >= 0) aux_val += mst->getX(k-m);
		}
		values[index] = aux_val;
		indices[index] = index;
		index++;
	}

	for(int i = 0; i < n; i++){
		aux_val = 0.0;
		for(int k = 0; k <= i ; k++){
			aux_val += mst->getX(k);
			aux_val -= dm[k];
		}
		values[index] = aux_val;
		indices[index] = index;
		index++;
	}

	status = CPXchgrhs(env, lp, index, indices, values);
	if(status != 0){
		printf("\n PROBLEMS TO CHANGE THE RHS!!!! ERRO CODE = %d\n", status);
		return -1;
	}

	CPXchgprobtype(env, lp, 0);
	CPXchgobjsen(env, lp, CPX_MAX);


	//status = CPXsetdblparam(env, CPX_PARAM_TILIM/*1039*/ , 300.0);
	/*if(status != 0){
		printf("\n PROBLEMS TO ADD TIME LIMIT TO THE MODEL!!!! ERRO CODE = %d\n", status);
		return -1;
	}*/

	//status = CPXwriteprob (env, lp, "MILP.lp", NULL);
	//if(status != 0)printf("\nPROBLEMS TO WRITE THE MODEL!!!! ERRO CODE = %d\n", status);

	status = CPXlpopt (env, lp);
	if(status != 0)printf("\n PROBLEMS TO SOLVE THE MODEL!!!! ERRO CODE = %d\n", status);

	//status = CPXgetstat(env, lp);
	//if(status != 0)printf("\n PROBLEMS TO SOLVE THE MODEL!!!! ERRO CODE = %d\n", status);

	double val;
	//CPXgetobjval(env, lp, &val);
	//if(status != 0)printf("\nPROBLEMS TO RETRIVE OPTIMAL SOLUTION!!!! ERRO CODE = %d\n", status);

	double y[var];
	status = CPXgetx(env, lp, y, 0, CPXgetnumcols(env, lp)-1);

	/****Determining the solution **/
	bool myalpha[n];
	bool mybeta[n];
	for(int i = 0; i < n; i++){
		myalpha[i] = (y[alpha[i]] < 0.5);
		mybeta[i] = (y[beta[i]] < 0.5);
	}

	int cnt = 0;
	int indice[2*n];
	char type[2*n];
	double bound[2*n];

	for(int i = 0; i < n; i++){
		indice[i] = alpha[i];
		indice[i+n] = beta[i];
		type[i] = 'B';
		type[i+n] = 'B';

		if(myalpha[i])bound[i] = 1.0;
		else bound[i] = 0.0;

		if(mybeta[i]) bound[i+n] = 1.0;
		else bound[i+ n] = 0.0;

	}

	CPXchgbds(env, lp, cnt, indice, type, bound);

	status = CPXlpopt (env, lp);
	if(status != 0)printf("\n PROBLEMS TO SOLVE THE MODEL!!!! ERRO CODE = %d\n", status);

	CPXgetobjval(env, lp, &val);
	if(status != 0)printf("\nPROBLEMS TO RETRIVE OPTIMAL SOLUTION!!!! ERRO CODE = %d\n", status);

	status = CPXgetx(env, lp, y, 0, CPXgetnumcols(env, lp)-1);

	/****Determining the solution **/
	double mydv[n];
	for(int i = 0; i < n; i++){
		mst->setDem(i, y[xi[i]] - y[xi2[i]]);
		mydv[i] = (y[xi[i]] - y[xi2[i]])*dv[i];
	}

	double X[n];
	double D[n];
	double PI[n];

	for(int i = 0; i < n; i++){
		X[i] = 0;
		D[i] = 0;
		PI[i] = 0;
		for(int k = 0; k <= i; k++){
			X[i] += mst->getX(k);
			D[i] += dm[k] + mydv[k];
		}

		if((i - m)>0) PI[i] = (X[i-m] - D[i]) > 0? (X[i-m] - D[i]): 0.0;
		else PI[i] = 0.0;
	}

	val = 0;
	val += (X[0] - D[0] ) > 0? h[0]*(X[0] - D[0]) : -p[0]*(X[0] - D[0]);
	for(int i = 1; i < n; i++){
		val += (X[i] - D[i] - PI[i]) > 0? h[i]*(X[i] - D[i] - PI[i]) : -p[i]*(X[i] - D[i] - PI[i]);
		val += q[i]*PI[i];
	}

	//erase the fake bounds and unfix the variables
	for(int i = 0; i < n; i++){
		indice[i] = alpha[i];
		indice[i+n] = alpha[i];
		type[i] = 'L';
		type[i+n] = 'U';
		bound[i] = 1.0;
		bound[i+ n] = 0.0;
	}

	CPXchgbds(env, lp, cnt, indice, type, bound);

	//erase the fake bounds and unfix the variables
	for(int i = 0; i < n; i++){
		indice[i] = beta[i];
		indice[i+n] = beta[i];
		type[i] = 'L';
		type[i+n] = 'U';
		bound[i] = 1.0;
		bound[i+ n] = 0.0;
	}

	CPXchgbds(env, lp, cnt, indice, type, bound);


	return val;
}
