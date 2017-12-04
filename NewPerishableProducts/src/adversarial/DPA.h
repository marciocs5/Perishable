/*
 * DPA.h
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#ifndef SRC_ADVERSARIAL_DPA_H_
#define SRC_ADVERSARIAL_DPA_H_

#include "Adversarial.h"

typedef struct _node{
	double val;
	bool xi;
	struct _node* next;
}Dnode;

class DPA: public Adversarial {
	Dnode*** beta;
	double* X;
	double* D;
	int delta;

	void clearTable();

public:
	DPA(Master* mst, int n, int m, int gamma, int* dv, double* dm, double* h, double* p, double* q);

	virtual ~DPA();

	double solve();
};

#endif /* SRC_ADVERSARIAL_DPA_H_ */
