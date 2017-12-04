/*
 * Adversarial.cpp
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#include "Adversarial.h"

Adversarial::Adversarial(Master* vmst, int vn, int vm, int vgamma, int* vdv, double* vdm, double* vh, double* vp, double* vq){
	n = vn;
	m = vm;
	gamma = vgamma;
	dv = new int[n];
	dm = new double[n];
	h = new double[n];
	p = new double[n];
	q = new double[n];
	mst = vmst;

	for(int i = 0; i < n; i++){
		dm[i] = vdm[i];
		dv[i] = vdv[i];
		h[i] = vh[i];
		p[i] = vp[i];
		q[i] = vq[i];
	}
}

int Adversarial::getGamma(){
	return gamma;
}

Adversarial::~Adversarial(){
	//TODO
}


