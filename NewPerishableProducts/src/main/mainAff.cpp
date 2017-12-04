/*
 * mainAff.c
 *
 *  Created on: 20 de set de 2017
 *      Author: marcio
 */

/*
 * main.cpp
 *
 *  Created on: 08/12/2015
 *      Author: santosma
 */

#include "Timer.h"
#include "../LBs/AffDecRules.h"
#include <iomanip>


void readVector(FILE* file, double** vec, long n){
	float fl;
	int ret;
	for(int i = 0; i < n; i++){
		ret = fscanf(file, "%f", &fl);
		if(ret < 0) printf("PROBLEMS TO READ THE FILE\n");
		(*vec)[i] = fl;
	}
}

void readVector(FILE* file, int** vec, long n){
	float fl;
	int ret;
	for(int i = 0; i < n; i++){
		ret = fscanf(file, "%f", &fl);
		if(ret < 0) printf("PROBLEMS TO READ THE FILE\n");
		(*vec)[i] = (int)fl;
	}
}

void InstanceReader(FILE* file, int* n, double** dm, int** dv, double** c, double** cap, double** p, double** h, double** q){
	int ret = fscanf(file, "%d", n);
	if(ret < 0) printf("PROBLEMS TO READ THE FILE \n");

	*c = new double[(*n)];
	readVector(file, c, (*n));

	*cap = new double[(*n)];
	readVector(file, cap, (*n));

	*p = new double[(*n)];
	readVector(file, p, (*n));

	*h = new double[(*n)];
	readVector(file, h, (*n));

	*dm = new double[(*n)];
	readVector(file, dm, (*n));

	*dv = new int[(*n)];
	readVector(file, dv, (*n));

	*q = new double[(*n)];
	readVector(file, q, (*n));
}

int main(int argn, char** args){

	FILE* result = fopen(args[1], "a");
	FILE* file = fopen(args[2], "r");
	int m = atoi(args[3]);
	int gamma = atoi(args[4]);

	int n = 0;
	double* dm = NULL;
	int* dv = NULL;
	double* c = NULL;
	double* cap = NULL;
	double* p = NULL;
	double* h = NULL;
	double* q = NULL;
	InstanceReader(file, &n, &dm, &dv, &c, &cap, &p, &h, &q);

	//printf("READ OK! \n");

	AffDecRules* solver = new AffDecRules(n, m, c,h, p, q, cap, dm, dv, gamma);
	//printf("MODEL OK \n");
	//getchar();

	//printf("MASTER OK \n");
	//getchar();

	solver->solve();
	double solution = solver->getSolution();

	printf("Solution = %f \n", solution);
	//fprintf(result, "%d \t %d \t %d \t %d \t %d \t %f \t %f \t %f \t %f \t %f \t %d \t %d \n", solvernum, heurnum, gamma, n, m, solution, solver->getTotalTime(), solver->getMPTime(), solver->getEXACTAPTime(), solver->getHEURTAPTime(), solver->getEXACTAPnum(), solver->getHEURAPnum());
	//printf("FILE=%s \n GAMMA=%d \n N=%d \n M=%d \n SOL=%f \n TT=%f \n TM=%f \n TADV=%f \n TADVHUR=%f \n NUM_EX=%d \n NUM_HEU = %d \n", args[2], gamma, n, m, solution, solver->getTotalTime(), solver->getMPTime(), solver->getEXACTAPTime(), solver->getHEURTAPTime(), solver->getEXACTAPnum(), solver->getHEURAPnum());

	return 0;
}



