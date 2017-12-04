/*
 * HeurMILP2.h
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#ifndef SRC_ADVERSARIAL_HEURMILP2_H_
#define SRC_ADVERSARIAL_HEURMILP2_H_

#include "Adversarial.h"
#include "MILP.h"

class HeurMILP2: public MILP{
public:
	HeurMILP2(Master* mst, int n, int m, int gamma, int* dv, double* dm, double* h, double* p, double* q);

	virtual ~HeurMILP2();

	double solve();
};

#endif /* SRC_ADVERSARIAL_HEURMILP2_H_ */
