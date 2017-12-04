/*
 * HeurMILP1.cpp
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#include "HeurMILP1.h"

HeurMILP1::HeurMILP1(Master* mst, int n, int m, int gamma, int* dv, double* dm, double* h, double* p, double* q) : MILP(mst, n, m, gamma, dv, dm, h, p, q){
}

HeurMILP1::~HeurMILP1(){
}

double HeurMILP1::solve(){
	int status;
	int indices[2*n];
	double values[2*n];
	int index = 0;
	double aux;

	for(int i = 0; i < n; i++){

		aux = 0.0;
		for(int k = 0; k <= i ; k++){
			aux -= dm[k];
			if((k - m) >= 0) aux += mst->getX(k-m);
		}
		values[index] = aux;
		indices[index] = index;
		index++;
	}

	for(int i = 0; i < n; i++){
		aux = 0.0;
		for(int k = 0; k <= i ; k++){
			aux += mst->getX(k);
			aux -= dm[k];
		}
		values[index] = aux;
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

	double val = 0;
	//CPXgetobjval(env, lp, &val);
	//if(status != 0)printf("\nPROBLEMS TO RETRIVE OPTIMAL SOLUTION!!!! ERRO CODE = %d\n", status);

	double y[var];
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

	return val;
}
