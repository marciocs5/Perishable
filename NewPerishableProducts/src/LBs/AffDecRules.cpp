/*
 * AffDecRules.cpp
 *
 *  Created on: 18 de set de 2017
 *      Author: marcio
 */

#include "AffDecRules.h"

AffDecRules::AffDecRules(int vn, int vm, double* vc, double* vh, double* vp, double* vq, double* vcap,  double* vdm, int* vdv, int vgamma){
	timer = new Timer();

	MAX = 0;
	n = vn;
	m = vm;
	c = new double[n];
	h = new double[n];
	p = new double[n];
	q = new double[n];
	cap = new double[n];
	dm = new double[n];
	dv = new double[n];
	var = 0;
	gamma = vgamma;
	cons_n = 3*n+n*n+1;

	for(int i = 0; i < n; i++){
		c[i] = vc[i];
		h[i] = vh[i];
		p[i] = vp[i];
		q[i] = vq[i];
		cap[i] = vcap[i];
		dm[i] = vdm[i];
		dv[i] = vdv[i];
		MAX += (h[i] + p[i] + c[i])*(dm[i]+dv[i]);
	}

	w = new double*[n];

	for(int i = 0; i < n; i++){
		w[i] = new double[n];

		w[i][i] = 0;
		for(int k = 0; k < n; k++){
			w[i][k] = 0;

			for(int j = i; j < k; j++){
				w[i][k] += h[j];
			}

			for(int j = k; j < i; j++){
				w[i][k] += p[j];
			}
		}
	}

	x_main = new int[n];
	y = new int[n];
	x = new int*[n];

	for(int i = 0; i < n; i++){
		x_main[i] = -1;
		y[i] = -1;
		x[i] = new int[n];
		for(int j = 0; j < n; j++){
			x[i][j] = -1;
		}
	}

	x_affine_pos = new int**[n];
	x_affine_neg = new int**[n];
	for(int i = 0; i < n; i++){
		x_affine_pos[i] = new int*[n];
		x_affine_neg[i] = new int*[n];
		for(int j = 0; j < n; j++){
			x_affine_pos[i][j] = new int[n];
			x_affine_neg[i][j] = new int[n];
			for(int k = 0; k < n; k++){
				x_affine_pos[i][j][k] = -1;
				x_affine_neg[i][j][k] = -1;
			}
		}
	}

	y_affine_pos = new int*[n];
	y_affine_neg = new int*[n];

	for(int i = 0; i < n; i++){
		y_affine_pos[i] = new int[n];
		y_affine_neg[i] = new int[n];

		for(int j = 0; j < n; j++){
			y_affine_pos[i][j] = -1;
			y_affine_neg[i][j] = -1;
		}
	}


	dualvarspos = new int*[cons_n];//the dual variables
	dualvarsneg = new int*[cons_n];//the dual variables
	specialdualvars = new int[cons_n];//the dual variables that are multipleid by gamma

	for(int i = 0; i < cons_n; i++){
		dualvarspos[i] = new int[n];
		dualvarsneg[i] = new int[n];
		specialdualvars[i] = -1;

		for(int j = 0; j < n; j++){
			dualvarspos[i][j] = -1;//the dual variables
			dualvarsneg[i][j] = -1;//the dual variables
		}
	}


	A = new double*[cons_n];//restriction matrix
	B = new double*[cons_n];//the independent terms for each xi

	for(int i = 0; i < cons_n; i++){
		A[i] = new double[n*n + n];
		B[i] = new double[n];
		for(int j = 0; j < n*n + n; j++){
			A[i][j] = 0;
		}
		for(int j = 0; j < n; j++){
			B[i][j] = 0;
		}
	}


	for(int i = 0; i < n; i++){
		//full rest
		for(int j = 0; j < n; j++){
			A[0][i*n + j] = w[i][j];
		}
		if(i > m && i < n - m)A[0][n*n + i] = w[i][i-m] + q[i];

		//sum up to X
		for(int j = 0; j < n; j++){
			A[i+1][i*n + j] = 1.0;
		}
		if(i > m && i < n - m)A[i+1][n*n + i] = 1.0;

		//sum up to the demand
		for(int j = 0; j < n; j++){
			A[n+i+1][j*n + i] = 1.0;

		}
		B[n+i+1][i] = 1.0;

		//Xij positive
		for(int j = 0; j < n;j++){
			A[2*n + n*i+ 1][i*n+j] = 1.0;
		}

		//Xi^{\star} positive
		A[2*n + n*n+ 1+ i][n*n+i] = 1.0;

	}

	int status;
	env = CPXopenCPLEX (&status);
	lp = CPXcreateprob (env, &status, "DUALROBUSTLOTSIZING");



	//printf("REST FULL OK \n");
	//getchar();

	addThetaVar();

	addXVar();

	addXMainVar();

	addYVar();

	addDualVarPos();

	addDualVarNeg();

	addDualVarSpe();

	addX_affineVarPos();

	addX_affineVarNeg();

	addY_affineVarPos();

	addY_affineVarNeg();

	for(int i = 0; i < cons_n; i++)
		addPrimalCons(i);

	for(int i = 0; i < cons_n; i++){
		for(int j = 0; j < n;j++){
			//addDualConsPos(i, j);
			//addDualConsNeg(i, j);
		}
	}

}

AffDecRules::~AffDecRules(){
	//TODO
}


void AffDecRules::addThetaVar(){
	double obj = 1.0;
	double lb = 0.0;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "THETA");
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
	theta = var++;
}

void AffDecRules::addXVar(){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if(j - i < m){
				sprintf(name, "x[%d][%d]", i, j);
				CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
				x[i][j] = var++;
			}else x[i][j] = -1;
		}
	}
}

void AffDecRules::addXMainVar(){
	double obj;
	double lb = 0.0;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = 0; i < n; i++){
		sprintf(name, "X[%d]", i);
		obj = c[i];
		up = cap[i];
		CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
		x_main[i] = var++;
	}
}

void AffDecRules::addYVar(){
	double obj = 0.0;
	double lb = 0.0;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = m; i < n; i++){
		sprintf(name, "Y[%d]", i-m);
		CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
		y[i] = var++;
	}
}

void AffDecRules::addDualVarPos(){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = 0; i < cons_n; i++){
		for(int j = 0; j < n; j++){
			sprintf(name, "z_pos[%d][%d]", i, j);
			CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
			dualvarspos[i][j] = var++;
		}
	}
}

void AffDecRules::addDualVarNeg(){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = 0; i < cons_n; i++){
		for(int j = 0; j < n; j++){
			sprintf(name, "z_neg[%d][%d]", i, j);
			CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
			dualvarsneg[i][j] = var++;
		}
	}

}

void AffDecRules::addDualVarSpe(){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = 0; i < cons_n; i++){
		sprintf(name, "ZESP[%d]", i);
		CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
		specialdualvars[i] = var++;
	}
}

void AffDecRules::addX_affineVarPos(){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = 0; i < n; i++){
		for(int j = 0 ; j < ((i+m > n)?n: i+m); j++){
			for(int k = 0; k <= i; k++){
				sprintf(name, "x_af_pos[%d][%d][%d]", i, j, k);
				CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
				x_affine_pos[i][j][k] = var++;
			}
		}
	}
}

void AffDecRules::addX_affineVarNeg(){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = 0; i < n; i++){
		for(int j = 0 ; j < ((i+m > n)?n: i+m); j++){
			for(int k = 0; k <= i; k++){
				sprintf(name, "x_af_neg[%d][%d][%d]", i, j, k);
				CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
				x_affine_neg[i][j][k] = var++;
			}
		}
	}
}

void AffDecRules::addY_affineVarPos(){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = m; i < n; i++){
		for(int j = 0 ; j <= i; j++){
			sprintf(name, "y_af_pos[%d][%d]", i, j);
			CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
			y_affine_pos[i][j] = var++;
		}
	}
}

void AffDecRules::addY_affineVarNeg(){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];

	for(int i = m; i < n; i++){
		for(int j = 0 ; j <= i; j++){
			sprintf(name, "y_af_pos[%d][%d]", i, j);
			CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
			y_affine_pos[i][j] = var++;
		}
	}
}

void AffDecRules::addPrimalCons(int i){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'G';

	int rmatind[n*n + 4*n +2];
	double rmatval[n*n + 4*n + 2];

	if(i ==0 ){
		rmatind[nzcnt] = theta;
		rmatval[nzcnt++] = -1.0;
	}

	rmatind[nzcnt] = specialdualvars[i];
	rmatval[nzcnt++] = -gamma;


	for(int k = 0; k < n; k++){
		rmatind[nzcnt] = dualvarspos[i][k];
		rmatval[nzcnt++] = -1.0;
		rmatind[nzcnt] = dualvarsneg[i][k];
		rmatval[nzcnt++] = -1.0;

		if(i < n+1 && i > 0){
			rmatind[nzcnt] = x_main[i-1];
			rmatval[nzcnt++] = -1.0;
		}

		if(A[i][n*n + k] > 0 && y[k] >= 0){
			rmatind[nzcnt] = y[k];
			rmatval[nzcnt++] = A[i][n*n + k];
		}

		for(int j = 0; j < n; j++){
			if(A[i][k*n + j] > 0 && x[k][j] >= 0){
				rmatind[nzcnt] = x[k][j];
				rmatval[nzcnt++] = A[i][k*n + j];
			}
		}
		rhs+= B[i][k]*dm[k];
	}

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}



void AffDecRules::addDualConsPos(int i, int j){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'G';

	int rmatind[n*n+n+2];
	double rmatval[n*n+n+2];

	rmatind[nzcnt] = specialdualvars[i];
	rmatval[nzcnt++] = 1.0;

	for(int k = 0; k < n; k++){
		if(A[i][n*n + k+1] > 0 && y_affine_pos[k][j] >= 0){
			rmatind[nzcnt] = y_affine_pos[k][j];
			rmatval[nzcnt++] = A[i][n*n + k+1];
		}
		for(int l = 0; l < n; l++){
			if(A[i][k*n + l +1] > 0 && x_affine_pos[k][l][j] >= 0){
				rmatind[nzcnt] = x_affine_pos[k][l][j];
				rmatval[nzcnt++] = A[i][k*n + l +1];
			}
		}

		rhs+= B[i][k]*dv[k];
	}

	rmatind[nzcnt] = dualvarspos[i][j];
	rmatval[nzcnt++] = 1.0;

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void AffDecRules::addDualConsNeg(int i, int j){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'G';

	int rmatind[n*n+n+2];
	double rmatval[n*n+n+2];

	rmatind[nzcnt] = specialdualvars[i];
	rmatval[nzcnt++] = 1.0;

	for(int k = 0; k < n; k++){
		if(A[i][n*n + k+1] > 0 && y_affine_neg[k][j] >= 0){
			rmatind[nzcnt] = y_affine_neg[k][j];
			rmatval[nzcnt++] = A[i][n*n + k+1];
		}
		for(int l = 0; l < n; l++){
			if(A[i][k*n + l +1] > 0 && x_affine_neg[k][l][j] >= 0){
				rmatind[nzcnt] = x_affine_neg[k][l][j];
				rmatval[nzcnt++] = A[i][k*n + l +1];
			}
		}

		rhs-= B[i][k]*dv[k];
	}

	rmatind[nzcnt] = dualvarsneg[i][j];
	rmatval[nzcnt++] = 1.0;

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void AffDecRules::solve(){

	CPXchgprobtype(env, lp, 0);
	CPXchgobjsen(env, lp, CPX_MIN);
	int status;

	status = CPXwriteprob (env, lp, "AFFINE DECISIONRULES.lp", NULL);
	if(status != 0)printf("\nPROBLEMS TO WRITE THE MODEL!!!! ERRO CODE = %d\n", status);

	//printf("WRITE OK");
	//getchar();

	timer->start();
	status = CPXlpopt(env, lp);
	if(status != 0)printf("\n PROBLEMS TO SOLVE THE MODEL!!!! ERRO CODE = %d\n", status);
	timer->stop();
	//printf("SOL OK");
	//getchar();

	double val;
	CPXgetobjval(env, lp, &val);
	if(status != 0)printf("\nPROBLEMS TO RETRIVE OPTIMAL SOLUTION!!!! ERRO CODE = %d\n", status);
	//printf("GET OK");
	//getchar();

	solution = val;

}

double AffDecRules::getSolution(){
	return solution;
}

float AffDecRules::getElapsedTime(){
	return timer->get();
}

long AffDecRules::getVarNumber(){
	return var;
}
