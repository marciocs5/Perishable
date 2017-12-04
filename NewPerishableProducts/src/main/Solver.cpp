/*
 * Solver.cpp
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#include "Solver.h"

#define ZERO 0.05

Solver::Solver(Master* vmp, Adversarial* vap, Adversarial* vheu, int vn){
	n = vn;

	numadv = 0;
	numheu = 0;

	tt = new Timer();
	tmp = new Timer();
	tsp = new Timer();
	thsp = new Timer();

	mp = vmp;
	ap = vap;
	heu = vheu;

	lastsol = 0;
}

double Solver::solve(){
	double* x;
	double* xi = NULL;
	double newsol=0;
	double solution;
	int it = 0;

	numadv = 0;
	numheu = 0;

	tt->start();
	do{
		it++;
		tmp->start();
		//printf("RESOLVENDO O MASTER......................... \n");
		lastsol = mp->solve(&solution);
		//printf("RESOLVI O MASTER! solution = %f sol adv = %f \n", solution, lastsol);
		//getchar();
		tmp->pause();

		if(heu != NULL){
			//printf("WITH HEURISTIC! \n");
			//getchar();
			thsp->start();
			newsol = heu->solve();
			numheu++;
			//printf("RESOLVI A HEUTRISTICA! sol = %f last = %f \n", newsol, lastsol);
			//getchar();
			thsp->pause();

			if(ap != NULL && ((newsol - lastsol)/newsol < ZERO)){
				tsp->start();
				newsol = ap->solve();
				//printf("RESOLVI O ADVERSARIO! sol = %f last = %f \n", newsol, lastsol);
				//getchar();
				numadv++;
				tsp->pause();
			}

		}else{
			//printf("NO HEURISTIC! \n");
			//getchar();

			if(ap != NULL){
				tsp->start();
				//printf("RESOLVENDO O ADVERSARIO....... \n");
				newsol = ap->solve();
				//printf("RESOLVI O ADVERSARIO! sol = %f last = %f \n", newsol, lastsol);
				//getchar();
				numadv++;
				tsp->pause();
			}else{printf("PROBLEMS!!!\n");}

		}

		//printf("Iteration: %d == SOL CUR = %f ||| SOL NEW = %f \n", it, lastsol, newsol);

	}while((newsol - lastsol)/newsol >= ZERO && it <= 100);
	tt->stop();
	tmp->stop();
	tsp->stop();
	thsp->stop();

	return solution;
}

double Solver::getTotalTime(){
	return tt->get();
}

double Solver::getMPTime(){
	return tmp->get();
}

double Solver::getEXACTAPTime(){
	return tsp->get();
}

double Solver::getHEURTAPTime(){
	return thsp->get();
}

int Solver::getEXACTAPnum(){
	return numadv;
}

int Solver::getHEURAPnum(){
	return numheu;
}

double Solver::getProd(){
	double prod = 0;

	for(int i = 0; i < n; i++) prod += mp->getX(i);

	return prod;
}

double Solver::getSpoiled(){
	double myX[n];
	double myD[n];

	double myDev[n];
	int maxdev[n];
	double auxD;
	int auxI;

	for(int i = 0; i < n; i++){
		myDev[i] = mp->getDv(i);
		maxdev[i] = i;
		myX[i] = mp->getX(i);
		myD[i] = mp->getDm(i);
	}

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n-1; j++){
			if(myDev[j] < myDev[j+1]){
				auxD = myDev[j];
				myDev[j] = myDev[j+1];
				myDev[j+1] = auxD;
				auxI = maxdev[j];
				maxdev[j] = maxdev[j+1];
				maxdev[j+1] = auxI;
			}
		}
	}

	int gamma = ap->getGamma();

	for(int i = 0; i < gamma; i++) myD[maxdev[i]] = myD[maxdev[i]] - myDev[i];

	double PIz[n];
	double piz[n];
	double val;
	double Xz[n];
	double Dz[n];

	Xz[0] = myX[0];
	Dz[0] = myD[0];
	PIz[0] = 0;
	piz[0] = 0;
	int m = mp->getM();

	for(int i = 1; i < n; i++){
		Xz[i] = Xz[i-1] + myX[i];
		Dz[i] = Dz[i-1] + myD[i];
		if(i - m < 0) piz[i] = 0;
		else piz[i] = (Xz[i-m] - Dz[i] - PIz[i-1]) > 0?  Xz[i-m] - Dz[i] - PIz[i-1]: 0;
		PIz[i] = PIz[i-1] + piz[i];
		val += piz[i];
	}

	return val;

}

Solver::~Solver(){
	delete tt;
	delete tsp;
	delete tmp;
}
