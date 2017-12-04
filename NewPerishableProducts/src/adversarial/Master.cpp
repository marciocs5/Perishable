/*
 * Master.cpp
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#include "Master.h"

#define MIN(a,b) a > b? b: a
#define POS(a,b) ((a-b) < 0) ? 0 : a-b

Master::Master(int vn, int vm, double* vc, double* vh, double* vp, double* vq, double* vcap,  double* vdm, int* vdv, int vgamma){
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

	x = new int[n];
	s = new int[n];
	r = new int[n];
	pi = new int[n];
	u = new int[n];
	alpha = new int[n];

	retX = new double[n];
	dem = new double[n];

	int status;
	env = CPXopenCPLEX (&status);
	lp = CPXcreateprob (env, &status, "DUALROBUSTLOTSIZING");

	for(int i = 0; i < n; i++){
		c[i] = vc[i];
		h[i] = vh[i];
		p[i] = vp[i];
		q[i] = vq[i];
		cap[i] = vcap[i];
		dm[i] = vdm[i];
		dv[i] = vdv[i];
		dem[i] = vdm[i];
		MAX += (h[i] + p[i] + c[i])*(dm[i]+dv[i]);
	}

	addTHETAVar();

	for(int i = 0; i < n; i++){
		addXVar(i);
	}
}

void Master::addXVar(int i){
	double obj = c[i];
	double lb = 0.0;
	double up = cap[i];
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "X[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
	x[i] = var++;
}

void Master::addTHETAVar(){
	double obj = 1.0;
	double lb = 0.0;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "THETA");
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
	theta=var++;
}

void Master::addSVar(int i){
	double obj = 0.0;
	double lb = 0.0;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "S[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
	s[i] = var++;
}

void Master::addRVar(int i){
	double obj = 0.0;
	double lb = 0.0;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "R[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
	r[i] = var++;
}

void Master::addPIVar(int i){
	double obj = 0.0;
	double lb = 0.0;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "PI[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
	pi[i] = var++;
}

void Master::addLinVars(int i){
	double obj = 0.0;
	double lb = -MAX;
	double up = MAX;
	char type = 'C';
	char* name = new char[30];
	sprintf(name, "U[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
	u[i] = var++;

	obj = 0.0;
	lb = 0.0;
	up = 1.0;
	type = 'B';
	sprintf(name, "alpha[%d]", i);
	CPXnewcols(env, lp, 1, &obj, &lb, &up, &type, &name);
	alpha[i] = var++;
}


double Master::solve(double* solution){//retorna a melhor solu��o para o adversario e solution eh o valor da solucao total
	addXiRest();

	CPXchgprobtype(env, lp, 1);
	CPXchgobjsen(env, lp, CPX_MIN);
	int status;

	/*status = CPXwriteprob (env, lp, "MASTER.lp", NULL);
	if(status != 0)printf("\nPROBLEMS TO WRITE THE MODEL!!!! ERRO CODE = %d\n", status);
	printf("MASTER OK!");
	getchar();*/

	status = CPXmipopt(env, lp);
	if(status != 0)printf("\n PROBLEMS TO SOLVE THE MODEL!!!! ERRO CODE = %d\n", status);

	double val;
	CPXgetobjval(env, lp, &val);
	if(status != 0)printf("\nPROBLEMS TO RETRIVE OPTIMAL SOLUTION!!!! ERRO CODE = %d\n", status);

	*solution = val;

	double y[var];

	status = CPXgetx (env, lp, y, 0, CPXgetnumcols(env, lp)-1);

	//debug!!
	/*printf("MASTER!!!\n");
	for(int i = 0; i < n; i++){
		printf(" S[%d] = %f \t R[%d] = %f \t U[%d] = %f \t Pi[%d] = %f \n", i , y[s[i]],i , y[r[i]], i, y[u[i]], i, y[pi[i]] );
	}
	printf("%f \n", val);*/

	for(int i = 0; i < n; i++) retX[i] = y[x[i]];

	bestval = y[theta];
	return y[theta];
}


double Master::getX(int i){
	return retX[i];
}

void Master::setDem(int i, double perc){
	dem[i] = dm[i] + perc*dv[i];
}

double Master::getDv(int i){
	return dv[i];
}

double Master::getDm(int i){
	return dm[i];
}

int Master::getM(){
	return m;
}

void Master::addXiRest(){

	/*printf("adiconando restricoes para ponto \n");
	for(int i=0; i < n; i++){
		if(dem[i] != dm[i]){
			printf(" %f ||", (dm[i]-dem[i])/ dem[i]);
		}
		else printf("0 ||");
	}
	printf("\n");
	getchar();*/

	for(int i = 0; i < n; i++)
		addSVar(i);

	for(int i = 0; i < n; i++)
		addRVar(i);

	for(int i = 0; i < n; i++)
		addPIVar(i);

	for(int i = 0; i < n; i++)
		addLinVars(i);

	addFULLRest();
	for(int i = 0; i < n; i++)
		addPiRest(i);

	for(int i = 0; i < n; i++)
		addRRest(i);

	for(int i = 0; i < n; i++)
		addSRest(i);
}

Master::~Master(){
	delete x;
	delete s;
	delete r;
	delete pi;
	delete alpha;
	delete u;
	delete retX;
	delete dem;
	delete c;
	delete h;
	delete p;
	delete q;
	delete cap;
	delete dm;
	delete dv;
}


void Master::addRRest(int i){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'G';

	int rmatind[2*n+1];
	double rmatval[2*n+1];

	for(int k = 0; k <= i; k++){
		rmatind[nzcnt] = x[k];
		rmatval[nzcnt++] = 1.0;

		rmatind[nzcnt] = pi[k];
		rmatval[nzcnt++] = -1.0;

		rhs += dem[k];
	}

	rmatind[nzcnt] = r[i];
	rmatval[nzcnt++] = 1.0;

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void Master::addSRest(int i){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'G';

	int rmatind[2*n+1];
	double rmatval[2*n+1];

	for(int k = 0; k <= i; k++){
		rmatind[nzcnt] = x[k];
		rmatval[nzcnt++] = -1.0;

		rmatind[nzcnt] = pi[k];
		rmatval[nzcnt++] = 1.0;
		rhs -= dem[k];
	}

	rmatind[nzcnt] = s[i];
	rmatval[nzcnt++] = 1.0;

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void Master::addFULLRest(){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'G';

	int rmatind[3*n+1];
	double rmatval[3*n+1];

	rmatind[nzcnt] = theta;
	rmatval[nzcnt++] = 1.0;

	for(int k = 0; k < n; k++){
		rmatind[nzcnt] = pi[k];
		rmatval[nzcnt++] = -q[k];
		rmatind[nzcnt] = r[k];
		rmatval[nzcnt++] = -p[k];
		rmatind[nzcnt] = s[k];
		rmatval[nzcnt++] = -h[k];
	}

	int rmatbeg[2] = {0, nzcnt-1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}


void Master::addPiRest(int i){
	addPIRest1(i);
	addPIRest2(i);
	addPIRest3(i);
	addURest(i);
}

void Master::addPIRest1(int i){
	int nzcnt = 2; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'G';

	int rmatind[2] = {pi[i], u[i]};
	double rmatval[2] = {1.0, -1.0};

	int rmatbeg[2] = {0, 1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void Master::addPIRest2(int i){
	int nzcnt = 3; //numero de coeficientes nao zero
	double rhs = MAX;// lado direito

	char sense = 'L';

	int rmatind[3] = {pi[i], u[i], alpha[i]};
	double rmatval[3] = {1.0, -1.0, MAX};

	int rmatbeg[2] = {0, 2};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void Master::addPIRest3(int i){
	int nzcnt = 2; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'L';

	int rmatind[2] = {pi[i], alpha[i]};
	double rmatval[2] = {1.0, -MAX};

	int rmatbeg[2] = {0, 1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}

void Master::addURest(int i){
	int nzcnt = 0; //numero de coeficientes nao zero
	double rhs = 0.0;// lado direito

	char sense = 'E';

	int rmatind[2*n+1];
	double rmatval[2*n+1];

	rmatind[nzcnt] = u[i];
	rmatval[nzcnt++] = 1.0;
	for(int k = 0; k <= i; k++){
		if(k - m >= 0){
			rmatind[nzcnt] = x[k-m];
			rmatval[nzcnt++] = -1.0;
		}
		rhs -= dem[k];
	}

	for(int k = 0; k < i; k++){
		rmatind[nzcnt] = pi[k];
		rmatval[nzcnt++] = 1.0;
	}

	int rmatbeg[2] = {0, nzcnt - 1};

	CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense, rmatbeg,rmatind, rmatval, NULL, NULL);
}
