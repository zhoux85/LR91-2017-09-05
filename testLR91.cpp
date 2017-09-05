/* The phase 1 Luo-Rudy(LR1) ventricular model, including four kinds of Euler methods */
/* Xiang Zhou, 2017/08/30 */
/* This code requires a C++ compiler */
/* Detailed list of equations and model description are provided in */
/* Circ Res 1991;68:1501-1526 */

/* IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USER SHOULD PACE THE MODEL UNTIL STEADY-STATE IS REACHED.
THIS REQUIRES APPROXIMATELY 5-20 MINUTES OF PACING DEPENDING ON THE RATE.*/

#include "stdafx.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "LR91.h"
#include "numerical_methods.h"

using std::cout;
using std::endl;

int bcl = 1000; // Basic Cycle Length (ms)
#define beats 1 // Number of Beats

/***for Hybrid_method--step change***/
double dv_min = 0.05;//mv
double dv_max = 0.2;//mv
double dt_max = 0.1;//ms
double dt_min = 0.001;//ms
double old_dt, new_dt;
int cycletmp = 0;

/***for CCL--step change***/
double Dv1_tmp = 0, Dv2 = 0, dt_tmp = 0;
double Voffset = 0.1;//mv, a given membrane potential offset
double dt_univ = 0.1;//ms universal maximum time-step siz

/* Action Potential Duration and Max. Info */
//double vmax[beats];// = { 51.1049, 51.075, 51.0754 }; // Max. Voltage (mV) 
double vmax[beats] = { 42.9639 }; // Max. Voltage (mV) 
double dvdtmax[beats]; // Max. dv/dt (mV/ms) 
double apd[beats]; // Action Potential Duration 
double toneapd[beats]; // Time of dv/dt Max. 
double ttwoapd[beats]; // Time of 90% Repolarization 
double rmbp[beats]; // Resting Membrane Potential 
double nair[beats]; // Intracellular Na At Rest 
double cair[beats]; // Intracellular Ca At Rest 
double kir[beats]; // Intracellular K At Rest 
double caimax[beats]; // Peak Intracellular Ca 
long f_count = 0;
double TNP[beats], CPUtime; // TNP = ttwoapd - toneapd;
double v_onset; //Initial Voltage (mv)
double stimulus_start; //Time to begin stimulus
int i = -1; // Stimulation Counter

double average(double *a);//return the average value of array a
void performance();

/* Creation of Data File */
FILE *ap;
FILE *fmaxs;
FILE *fpara;
FILE *fevaluation;

/* Values are printed to a file called ap */
void prttofile();

int main(int argc, char* argv[]){
	/* Opening of Datafiles */
	ap = fopen("ap", "w");
	fpara = fopen("fpara", "w");
	fmaxs = fopen("fmaxs", "w");
	fevaluation = fopen("fevaluation", "w");

	/* Time Loop Conditions */
	t = 0.0; // Time (ms)
	dt = 0.01; // Time step (ms)
	st = -80.0; // Stimulus (mA)
	tstim = 10.0; // Time to begin stimulus
	stimulus_start = 10.0;
	stimtime = 10.0; // Initial Condition for Stimulus (ms)
	v = -88.654973; // Initial Voltage (mv)
	v_onset = v;

	/* Beginning Ion Concentrations */
	nai = 18; // Initial Intracellular Na (mM)
	nao = 140; // Initial Extracellular Na (mM)
	ki = 145; // Initial Intracellular K (mM)
	ko = 5.4; // Initial Extracellular K (mM)
	cai = 0.0002; // Initial Intracellular Ca (mM)
	cao = 1.8; // Initial Extracellular Ca (mM)

	/* Initial Gate Conditions */
	m = 0.000838;
	h = 0.993336;
	j = 0.995484;
	d = 0.000003;
	f = 0.999745;
	xr = 0.000129;

	steps = (bcl*beats) / dt; // Number of steps
	
	clock_t start, end;
	start = clock();
	/* Beginning of Time Loop */
	for (increment = 0; t <= (bcl*beats); increment++) {
		v_0 = v;
		/***********Forward Euler method***********/
		//vnew = Forward_Euler();

		/***********backward Euler method***********/
		//vnew = Backward_Euler();

		/***********Trapezoid method***********/
		//vnew = Trapezoid();

		/***********Improved Euler method, i.e.£¬one iteration in Trapezoid method**********/
		//vnew = Improved_Euler();

		/***********Fourth-order Runge-Kutta method***********/
		//vnew = RK_method();

		//performance();
		//v = vnew;
		//t = t + dt;
		//prttofile();

		/*********adaptive time step**********************************************/
		/***********Hybrid_method *****************/
		//int cycle = t / 1000;//cycle Counter
		//if (cycle>cycletmp){
		//	cycletmp = cycle;
		//}
		//Hybrid_method(&vnew, &new_dt, cycle);

		/***********CCL method*****************/
		CCL(&vnew, &new_dt);

		performance();
		v = vnew;
		t = t + dt;
		prttofile();
		dt = new_dt;
		/*********adaptive time step**********************************************/
	}
	end = clock();
	double time_used = (double)(end - start) / CLK_TCK;
	printf("time_used__%lf\n", time_used);
	for (i = 0; i<beats; i++){
		apd[i] = ttwoapd[i] - toneapd[i];
		fprintf(fmaxs, "%g\t,", vmax[i]);
	}

	fprintf(fevaluation, "%g\t%g\t%g\t%g\t%g\t%ld\n", average(vmax), average(dvdtmax),
		average(apd), average(TNP), time_used, f_count);
	//	fprintf(fevaluation, "%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%ld\n",
	//		Vmax, dvdt_max, APD90, TNP, time_used, f_count);
	//	fprintf(fpara, "%.5f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
	//		t, v, nai, ki, cai, nao, ko, cao, m, h, j, d, f, xr);

	//cout << "\a\a";
	return (1);
}

double average(double *a)//return the average value of array a
{
	double r = 0;
	int k;
	for (k = 0; k<beats; ++k)
		r += a[k];
	return r / beats;
}

/* Values are printed to a file called ap. The voltage and
currents can be plotted versus time using graphing software. */
void prttofile() {
	if (t>(0) && t<(bcl*beats))	{
		//fprintf(ap, "%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
		//	t, v, nai, ki, cai, ina, isi, ikr, iki, ikp, ib);
		fprintf(ap, "%.5f\t%g\n", t, v);
		//printf("%.5f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
		//	t, v, nai, ki, cai, ina, isi, ikr, iki, ikp, ib);
	}
	//nai, ki, cai are the Intracellular Concentration of nai, ki, cai
}

//performance 
void performance(){
	if (t >= stimulus_start){
		stimulus_start = stimulus_start + bcl;
		i = i + 1; //i+1 when stimulus is applied, t>=10ms
		rmbp[i] = v;// Resting Membrane Potential 
		nair[i] = nai;// Intracellular Na At Rest 
		kir[i] = ki;
		cair[i] = cai;
	}
	dvdt = abs(vnew - v_0) / dt;
	dvdtnew = dvdt;
	if (i >= 0)
	{
		if (vnew>vmax[i])
			vmax[i] = vnew;
		if (dvdtnew>dvdtmax[i]){
			dvdtmax[i] = dvdtnew;
			toneapd[i] = t; //toneapd[beats]; // Time of dv/dt Max. 
		}
		if (vnew >= (vmax[i] - 0.9*(vmax[i] - rmbp[i])))
			ttwoapd[i] = t; // ttwoapd[beats]; // Time of 90% Repolarization 
		if (vnew> vmax[i] * 0.95){
			TNP[i] = TNP[i] + dt;
		}
	}
}
