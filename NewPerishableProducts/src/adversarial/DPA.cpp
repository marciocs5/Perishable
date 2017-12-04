/*
 * DPA.cpp
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#include "DPA.h"
#include <cstdlib>
#include <cstdio>
#include <algorithm>    // std::swap
#include <limits.h>

DPA::DPA(Master* mst, int n, int m, int gamma, int* vdv, double* dm, double* h, double* p, double* q) : Adversarial(mst, n, m, gamma, vdv, dm, h, p, q){

	X = new double[n];
	D = new double[n];

	delta = 0;
	for(int i = 0; i < n; i++){
		delta += dv[i];
	}


	beta = new Dnode**[n];
	if(beta == NULL){
		printf("MEMORY INSUFICIENT!!!\n");
		exit(1);
	}
	for(int i = 0; i< n; i++){
		beta[i] = new Dnode*[gamma+1];
		if(beta[i] == NULL){
			printf("MEMORY INSUFICIENT!!!\n");
			exit(1);
		}
		for(int j = 0; j <= gamma; j++){
			beta[i][j] = new Dnode[delta+1];
			if(beta[i][j] == NULL){
				printf("MEMORY INSUFICIENT!!!\n");
				exit(1);
			}
		}
	}
}

DPA::~DPA(){
	//TODO
}

double DPA::solve(){
	double* ret = new double[n];
	double sol;

	D[0] = dm[0];
	X[0] = mst->getX(0);
	for(int i = 1; i < n; i++){
		D[i] = D[i-1] + dm[i];
		X[i] = X[i-1] + mst->getX(i);
	}

	Dnode* node;
	Dnode* node2;
	Dnode* node3;
	//base cases
	for(int j = 0; j <= gamma; j++){
		for(int i = 0; i <= delta; i++){
			node = &beta[0][j][i];
			node->val = INT_MIN;
			node->xi = false;
			node->next = NULL;
		}
	}

	beta[0][0][0].val = ((X[0] - D[0]) > 0.0 )? h[0]*(X[0] - D[0]) : -p[0]*(X[0] - D[0]);
	beta[0][0][0].xi = false;
	beta[0][0][0].next = NULL;

	if(gamma > 0){
		beta[0][1][dv[0]].val = ((X[0] - D[0] - dv[0]) > 0.0 )? h[0]*(X[0] - D[0] - dv[0]) : -p[0]*(X[0] - D[0] - dv[0]);
		beta[0][1][dv[0]].xi = true;
		beta[0][1][dv[0]].next = NULL;
	}

	double val;
	//printf("k \n");
	//getchar();

	for(int i = 1; i < n; i++){
		for(int j = 0; j <= gamma; j++){
			for(int phi= 0; phi <= delta; phi++){
				val = ((X[i] - D[i] - dv[i]) > 0.0 )? h[i]*(X[i] - D[i] - phi) : -p[i]*(X[i] - D[i] - phi);
				node = &beta[i][j][phi];
				node2 = &beta[i-1][j][phi];

				node->val = node2->val + val;
				node->xi = false;
				node->next = node2;
				if(phi - dv[i] >= 0 && j - 1 >= 0){
					node3 = &beta[i-1][j-1][phi - dv[i]];
					if(node2->val < node3->val){
						node->val = node3->val + val;
						node->xi = true;
						node->next = node3;
					}
				}
			}
		}
	}

	Dnode* solution = &beta[n-1][0][0];
	sol = beta[n-1][0][0].val;
	//looking for solution
	for(int i = 0; i <= gamma; i++){
		for(int j = 0; j <= delta; j++){
			node = &beta[n-1][i][j];
			if(solution->val < node->val) {
				solution = node;
				sol = node->val;
			}
		}
	}

	double Xz[n];
	double Dz[n];
	double PIz[n];
	double piz[n];
	double mydv[n];

	for(int i = 0; i < n; i++){
		mst->setDem(i, 0.0);
		mydv[i] = 0.0;
	}

	for(int i = n-1; i >= 0; i--){
		if(solution->xi == true){
			mst->setDem(i, 1.0);
			mydv[i] = dv[i];
		}
		else{
			mst->setDem(i, 0.0);
		}
		solution = solution->next;
	}

	/*********************/

	val = 0;

	Xz[0] = mst->getX(0);
	Dz[0] = dm[0] + mydv[0];
	PIz[0] = 0;
	piz[0] = 0;
	for(int i = 1; i < n; i++){
		Xz[i] = Xz[i-1] + mst->getX(i);
		Dz[i] = Dz[i-1] + dm[i] + mydv[i];
		if(i - m < 0) piz[i] = 0;
		else piz[i] = (Xz[i-m] - Dz[i] - PIz[i-1]) > 0?  Xz[i-m] - Dz[i] - PIz[i-1]: 0;
		PIz[i] = PIz[i-1] + piz[i];
	}

	/*for(int i = 0; i < n; i++){
		printf("S[%d] = %f || ", i, X[i] - D[i] - PI[i]);
		printf("R[%d] = %f || ", i, -X[i] + D[i] + PI[i]);
		printf("PI[%d] = %f \n", i, pi[i]);
	}*/
	//getchar();

	val = (Xz[0] - Dz[0] ) > 0? h[0]*(Xz[0] - Dz[0]) : -p[0]*(Xz[0] - Dz[0]);
	for(int i = 1; i < n; i++){
		val += (Xz[i] - Dz[i] - PIz[i]) > 0? h[i]*(Xz[i] - Dz[i] - PIz[i]) : -p[i]*(Xz[i] - Dz[i] - PIz[i]);
		if(i > m)val += q[i]*piz[i];
	}

	return val;
}

