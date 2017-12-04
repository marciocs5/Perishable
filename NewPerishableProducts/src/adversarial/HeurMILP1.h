/*
 * HeurMILP1.h
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#ifndef SRC_ADVERSARIAL_HEURMILP1_H_
#define SRC_ADVERSARIAL_HEURMILP1_H_

#include "Adversarial.h"
#include "MILP.h"

class HeurMILP1: public MILP{
public:
	HeurMILP1(Master* mst, int n, int m, int gamma, int* dv, double* dm, double* h, double* p, double* q);

	virtual ~HeurMILP1();

	double solve();
};

#endif /* SRC_ADVERSARIAL_HEURMILP1_H_ */
