/*
 * Solver.h
 *
 *  Created on: 02/05/2016
 *      Author: santosma
 */

#ifndef SRC_MAIN_SOLVER_H_
#define SRC_MAIN_SOLVER_H_

#include "../adversarial/Master.h"
#include "../adversarial/Adversarial.h"
#include "Timer.h"

class Solver {
protected:

	int n;

	int numadv;
	int numheu;

	Timer* tt;
	Timer* tmp;
	Timer* tsp;
	Timer* thsp;

	Master* mp;
	Adversarial* ap;
	Adversarial* heu;

	double lastsol;
public:
	Solver(Master* mp, Adversarial* ap, Adversarial* heu, int n);

	double solve();

	double getTotalTime();

	double getMPTime();

	double getEXACTAPTime();

	double getHEURTAPTime();

	int getEXACTAPnum();

	int getHEURAPnum();

	double getProd();

	double getSpoiled();

	virtual ~Solver();
};

#endif /* SRC_MAIN_SOLVER_H_ */
