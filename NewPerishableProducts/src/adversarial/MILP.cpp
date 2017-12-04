/*
 * MIP.cpp
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#include "MILP.h"

MILP::MILP(Master* mst, int n, int m, int gamma, int* dv, double* dm, double* h, double* p, double* q) : Adversarial(mst, n, m, gamma, dv, dm, h, p, q){
	int status;
	env = CPXopenCPLEX (&status);
	lp = CPXcreateprob (env, &status, "MILDADVERARIAL");


	MAX = 0;
	for(int i = 0; i < n; i++){
		MAX += (h[i]+ p[i])*(dm[i] + dv[i]);
	}

	var = 0;

	theta = new int[n];
	beta = new int[n];
	s = new int[n];
	pi = new int[n];
	u = new int[n];
	alpha = new int[n];
	xi = new int[n];
	xi2 = new int[n];

	for(int i = 0; i < n; i++){
		addTHETAVar(i);
		theta[i] = var++;
		addBetaVar(i);
		beta[i] = var++;
		addSVar(i);
		s[i] = var++;
		addPIVar(i);
		pi[i] = var++;
		addUVar(i);
		u[i] = var++;
		addAlphaVar(i);
		alpha[i] = var++;
		addXiVar(i);
		xi[i] = var++;
		addXi2Var(i);
		xi2[i] = var++;
	}
	//printf("variaveis  OK!\n");
	//getchar();

	for(int i = 0; i< n; i++){
		addURest(i);
	}

	for(int i = 0; i< n; i++){
		addSRest(i);
	}

	addGammaRest();
	for(int i = 0; i< n; i++){
		addThetaRest(i);
		addPiRest(i);
	}
	//printf("Resti��es OK!\n");
	//getchar();
}

MILP::~MILP(){
	//TODO
}

double MILP::solve(){
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

	CPXchgprobtype(env, lp, 1);
	CPXchgobjsen(env, lp, CPX_MAX);


	//status = CPXsetdblparam(env, CPX_PARAM_TILIM/*1039*/ , 300.0);
	/*if(status != 0){
		printf("\n PROBLEMS TO ADD TIME LIMIT TO THE MODEL!!!! ERRO CODE = %d\n", status);
		return -1;
	}*/

	/*status = CPXwriteprob (env, lp, "MILP.lp", NULL);
	if(status != 0)printf("\nPROBLEMS TO WRITE THE MODEL!!!! ERRO CODE = %d\n", status);*/

	status = CPXmipopt (env, lp);
	if(status != 0)printf("\n PROBLEMS TO SOLVE THE MODEL!!!! ERRO CODE = %d\n", status);

	//printf("Otimização OK \n");
	//getchar();

	//status = CPXgetstat(env, lp);
	//if(status != 0)printf("\n PROBLEMS TO SOLVE THE MODEL!!!! ERRO CODE = %d\n", status);

	double val;
	CPXgetobjval(env, lp, &val);
	if(status != 0)printf("\nPROBLEMS TO RETRIVE OPTIMAL SOLUTION!!!! ERRO CODE = %d\n", status);

	double y[var];
	status = CPXgetx(env, lp, y, 0, CPXgetnumcols(env, lp)-1);

	/*printf("ADVERSARIO!!!\n");
	for(int i = 0; i < n; i++){
		printf("THETA[%d] = %f \t S[%d] = %f \t U[%d] = %f \t Pi[%d] = %f \n", i, y[theta[i]], i , y[s[i]], i, y[u[i]], i, y[pi[i]] );
	}
	printf("%f \n", val);*/


	for(int i = 0; i < n; i++){
		mst->setDem(i, y[xi[i]] - y[xi2[i]]);
	}

	return val;
}

void MILP::addXiVar(int i){
	double obj = 0.0;
	double lb = 0.0;
	double up = 1.0;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "XI+[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
}

void MILP::addXi2Var(int i){
	double obj = 0.0;
	double lb = 0.0;
	double up = 1.0;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "XI-[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
}

void MILP::addTHETAVar(int i){
	double obj = 1.0;
	double lb = 0.0;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "THETA[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
}

void MILP::addUVar(int i){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "U[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
}

void MILP::addAlphaVar(int i){
	double obj = 0.0;
	double lb = 0.0;
	double up = 1.0;
	char type = 'B';
	char* name = new char[30];
	sprintf(name, "ALPHA[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
}


void MILP::addBetaVar(int i){
	double obj = 0.0;
	double lb = 0.0;
	double up = 1.0;
	char type = 'B';
	char* name = new char[30];
	sprintf(name, "BETA[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
}

void MILP::addSVar(int i){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "S[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
}

void MILP::addPIVar(int i){
	double obj = q[i];
	double lb = 0.0;
	double up;
	if(i >= m) up = MAX;
	else up = 0.0;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "Pi[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
}

void MILP::addURest(int i){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'E';

	int rmatind[3*n+1];
	double rmatval[3*n+1];

	rmatind[nzcnt] = u[i];
	rmatval[nzcnt++] = 1.0;

	for(int k = 0; k <= i; k++){
		rmatind[nzcnt] = xi[k];
		rmatval[nzcnt++] = dv[k];
		rmatind[nzcnt] = xi2[k];
		rmatval[nzcnt++] = -dv[k];
		rhs -= dm[k];
		if(k - m >= 0){
			rhs += mst->getX(k-m);
		}
	}

	for(int k = 0; k < i; k++){
		rmatind[nzcnt] = pi[k];
		rmatval[nzcnt++] = 1.0;
	}

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void MILP::addSRest(int i){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'E';

	int rmatind[3*n+1];
	double rmatval[3*n+1];

	for(int k = 0; k <= i; k++){
		rhs += mst->getX(k);
		rmatind[nzcnt] = xi[k];
		rmatval[nzcnt++] = dv[k];
		rmatind[nzcnt] = xi2[k];
		rmatval[nzcnt++] = -dv[k];
		rhs -= dm[k];
	}

	for(int k = 0; k < i; k++){
		rmatind[nzcnt] = pi[k];
		rmatval[nzcnt++] = 1.0;
		rhs -= dm[k];
	}

	rmatind[nzcnt] = s[i];
	rmatval[nzcnt++] = 1.0;

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}


void MILP::addPiRest(int i){
	int nzcnt = 3; //numero de coeficientes nao zero
	double rhs = MAX;// lado direito
	char sense = 'L';
	int rmatind[3] = {alpha[i], u[i], pi[i]};
	double rmatval[3] = {MAX, -1.0, 1.0};
	int rmatbeg[3] = {0, nzcnt-1};
	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);

	int nzcnt2 = 2; //numero de coeficientes nao zero
	double rhs2 = 0.0;// lado direito
	char sense2 = 'L';
	int rmatind2[2] = {alpha[i], pi[i]};
	double rmatval2[2] = {-MAX, 1.0};
	int rmatbeg2[2] = {0, nzcnt-1};
	CPXaddrows(env, lp, 0, 1, nzcnt2, &rhs2, &sense2, rmatbeg2,rmatind2, rmatval2, NULL, NULL);

	int nzcnt3 = 2; //numero de coeficientes nao zero
	double rhs3 = 0.0;// lado direito
	char sense3 = 'G';
	int rmatind3[2] = {u[i], pi[i]};
	double rmatval3[2] = {-1.0, 1.0};
	int rmatbeg3[2] = {0, nzcnt-1};
	CPXaddrows(env, lp, 0, 1, nzcnt3, &rhs3, &sense3, rmatbeg3,rmatind3, rmatval3, NULL, NULL);
}

void MILP::addGammaRest(){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = gamma;// lado direito

	char sense = 'L';

	int rmatind[2*n+1];
	double rmatval[2*n+1];

	for(int k = 0; k < n; k++){
		rmatind[nzcnt] = xi[k];
		rmatval[nzcnt++] = 1.0;
		rmatind[nzcnt] = xi2[k];
		rmatval[nzcnt++] = 1.0;
	}

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void MILP::addThetaRest(int i){
	int nzcnt = 3; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito
	char sense = 'L';
	int rmatind[3] = {beta[i], s[i], theta[i]};
	double rmatval[3] = {-MAX, -h[i], 1.0};
	int rmatbeg[3] = {0, nzcnt-1};
	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);

	int nzcnt2 = 3; //numero de coeficientes nao zero
	double rhs2 = MAX;// lado direito
	char sense2 = 'L';
	int rmatind2[3] = {beta[i], s[i], theta[i]};
	double rmatval2[3] = {MAX, p[i], 1.0};
	int rmatbeg2[3] = {0, nzcnt-1};
	CPXaddrows(env, lp, 0, 1, nzcnt2, &rhs2, &sense2, rmatbeg2,rmatind2, rmatval2, NULL, NULL);
}
