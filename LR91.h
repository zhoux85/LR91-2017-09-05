/* List of variables and parameters (this code uses all global variables) */
#ifndef _LR91_H
#define _LR91_H

extern int bcl;// Basic Cycle Length (ms)
extern long f_count;//How many times has total current been calculated? 
extern double v_0;//start Voltage

/* Voltage */
extern double v; // Membrane voltage (mV)
extern double vnew; // New Voltage (mV)
extern double dvdt; // Change in Voltage / Change in Time (mV/ms)
extern double dvdtnew; // New dv/dt (mV/ms)

/* Time Step */
extern double dt; // Time step (ms)
extern double t; // Time (ms)
extern int nxstep; // Interval Between Calculating Ion Currents
extern int steps; // Number of Steps
extern int increment; // Loop Control Variable

/* Total Current and Stimulus */
extern double st; // Constant Stimulus (uA/cm^2)
extern double tstim; //Time Stimulus is Applied (ms)//Time to begin stimulus
extern double stimtime; //Time period during which stimulus is applied (ms)
extern double it; // Total current (uA/cm^2)

/* Ion Concentrations */
extern double nai; // Intracellular Na Concentration (mM)
extern double nao; // Extracellular Na Concentration (mM)
extern double cai; // Intracellular Ca Concentration (mM)
extern double cao; // Extracellular Ca Concentration (mM)
extern double ki; // Intracellular K Concentration (mM)
extern double ko; // Extracellular K Concentration (mM)

/* Fast Sodium Current (time dependant) */
extern double m; // Na activation
extern double h; // Na inactivation
extern double j; // Na slow inactivation

/* Current through L-type Ca Channel */
extern double dcai; // Change in myoplasmic Ca concentration (mM)
extern double d; // Voltage dependant activation gate
extern double f; // Voltage dependant inactivation gate

/* Rapidly Activating Potassium Current */
extern double xr; // Rapidly Activating K time-dependant activation  --gate X in LR91

//temporary gate value
extern double m0, h0, j_0;
extern double d0, f0;
extern double xr0;

/* Ion Current Functions */
//// "extern" can be omitted£»Leave it here for easy understanding
extern void comp_ina(); // Calculates Fast Na Current
extern void comp_ical(); // Calculates Currents through L-Type Ca Channel
extern void comp_ikr(); // Calculates Rapidly Activating K Current
extern void comp_iki(); // Calculates Time-Independent K Current
extern void comp_ikp(); // Calculates Plateau K Current
extern void comp_ib(); // Calculates Background Current
extern void comp_it(); // Calculates Total Current

#endif