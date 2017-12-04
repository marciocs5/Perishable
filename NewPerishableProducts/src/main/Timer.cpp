/*
 * Timer.cpp
 *
 *  Created on: 26/06/2014
 *      Author: santosma
 */

#include "Timer.h"

Timer::Timer(){
	seconds = 0;
	started = false;
}

Timer::~Timer(){}

void Timer::start(){
	started = true;
	gettimeofday(&clockinit, &tz);
}

float Timer::pause(){
	gettimeofday(&clockend, &tz);
	seconds += ((float)(clockend.tv_sec - clockinit.tv_sec))*1000000 + ((float)(clockend.tv_usec - clockinit.tv_usec));
	gettimeofday(&clockinit, &tz);
	return seconds;
}

float Timer::stop(){
	gettimeofday(&clockend, &tz);
	seconds += ((float)(clockend.tv_sec - clockinit.tv_sec))*1000000 + ((float)(clockend.tv_usec - clockinit.tv_usec));
	return seconds;
}

float Timer::get(){
	if(started)
	return seconds;
	else return 0;
}
