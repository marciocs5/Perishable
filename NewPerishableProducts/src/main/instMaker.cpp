/*
 * instMaker.cpp
 *
 *  Created on: 08/12/2015
 *      Author: santosma
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>

#define PI 3.14159265

double* c = new double[500];
double* cap = new double[500];
double* h = new double[500];
double* p = new double[500];
double* dm = new double[500];
double* dv = new double[500];
double theta = 0.1;
double* q = new double[500];
double capT = 6000;
double qT = 1.0;

void buildInstancePartialRandom(){
	srand(time(NULL));
	for(int i = 0; i <= 300; i++){
		c[i] = 10.0 + 10.0*(rand() % 10 + 1)/10.0;
		cap[i] = capT;
		h[i] = 2.0 + 2.0*(rand() % 10 + 1)/10.0;
		p[i] = 30.0  + 30.0*(rand() % 10 + 1)/10.0;
		dm[i] = 1000.0 + 1000.0*(rand() % 10 + 1)/10.0;
		dv[i] = theta*1000.0 + theta*1000.0*(rand() % 10 + 1)/10.0;
		q[i] = qT + qT*(rand() % 10 + 1)/10.0;
	}
}

void buildInstanceRandom(){
	srand(time(NULL));
	for(int i = 0; i <= 300; i++){
		c[i] = 1.0 + 20.0*(rand() % 10 + 1)/10.0;
		cap[i] = capT;
		h[i] = 1.0 + 4.0*(rand() % 10 + 1)/10.0;
		p[i] = 1.0  + 100.0*(rand() % 10 + 1)/10.0;
		dm[i] = 1.0 + 2000.0*(rand() % 10 + 1)/10.0;
		dv[i] = theta*1.0 + theta*2000.0*(rand() % 10 + 1)/10.0;
		q[i] = 1.0 + qT*(rand() % 10 + 1)/10.0;
	}
}

void buildInstanceStatic(){
	srand(time(NULL));
	for(int i = 0; i <= 300; i++){
		c[i] = 20.0;
		cap[i] = capT;
		h[i] = 4.0 ;
		p[i] = 100.0 ;
		dm[i] = 1000.0 ;
		dv[i] = theta*1000.0;
		q[i] = qT;
	}
}

void buildInstanceDynamic(){
	srand(time(NULL));
	for(int i = 0; i <= 300; i++){
		c[i] = 10.0 + 5.0*sin(i*PI*15/180);
		cap[i] = capT;
		h[i] = 2.0 + 1.0*sin(i*PI*15/180);
		p[i] = 50.0  + 25.0*sin(i*PI*15/180);
		dm[i] = 1000.0 + 500.0*sin(i*PI*15/180);
		dv[i] = theta*1000.0 + theta*500.0*sin(i*PI*15/180);
		q[i] = qT + qT*sin(i*PI*15/180)/2.0;
	}
}

void writeVecInfile(double* vec, FILE* file, long n){
	for(int i = 0; i < n; i++)
		if(vec[i] >= 0)
			fprintf(file,"%.2f\n", vec[i]);
		else fprintf(file,"0.0 \n");
}

void writefile(FILE* file, int i){
	writeVecInfile(c, file, i);
	writeVecInfile(cap, file, i);
	writeVecInfile(p, file, i);
	writeVecInfile(h, file, i);
	writeVecInfile(dm, file, i);
	writeVecInfile(dv, file, i);
	writeVecInfile(q, file, i);
}

int main(){

	double myCapT[2] = {5000.0, 10000000.0};
	double myQT[3] = {2.0, 20.0, 200.0};
	double myTheta[3] = {0.1, 0.2, 0.3};


	for(int i = 10; i <= 100; i = i+10){
		for(int j = 0; j < 2; j++){
			capT = myCapT[j];
			for(int k = 0; k < 3; k++){
				qT = myQT[k];
				for(int l = 0; l < 3; l++){
					theta = myTheta[l];


					FILE* file;
					char name[100];
					sprintf(name, "InstancePARTIALRANDOM_N=%d_CAP=%d_QT=%d_THETA=%d.ins", i,j, k, l);
					file = fopen(name, "w");
					buildInstancePartialRandom();
					fprintf(file, "%d\n", i);
					writefile(file, i);
					fclose(file);

					//sprintf(name, "InstanceRANDOM_N=%d_CAP=%d_QT=%d_THETA=%d.ins", i,j, k, l);
					//file = fopen(name, "w");
					//buildInstanceRandom();
					//fprintf(file, "%d\n", i);
					//writefile(file, i);
					//fclose(file);

					sprintf(name, "InstanceSTATIC_N=%d_CAP=%d_QT=%d_THETA=%d.ins", i,j, k, l);
					file = fopen(name, "w");
					buildInstanceStatic();
					fprintf(file, "%d\n", i);
					writefile(file, i);
					fclose(file);

					sprintf(name, "InstanceDYNAMIC_N=%d_CAP=%d_QT=%d_THETA=%d.ins", i,j, k, l);
					file = fopen(name, "w");
					buildInstanceDynamic();
					fprintf(file, "%d\n", i);
					writefile(file, i);
					fclose(file);
				}
			}
		}
	}
}


