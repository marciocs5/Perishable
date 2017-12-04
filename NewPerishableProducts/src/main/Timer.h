/*
 * Timer.h
 *
 *  Created on: 26/06/2014
 *      Author: santosma
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <sys/time.h>
#include <unistd.h>

class Timer {
private:
	struct timeval  clockinit;
	struct timeval  clockend;
	struct timezone tz;
	float seconds;
	bool started;
public:
	Timer();

	virtual ~Timer();

	void start();

	float pause();

	float stop();

	float get();

};

#endif /* TIMER_H_ */
