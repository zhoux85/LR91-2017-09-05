#include "stdafx.h"
#include <math.h>
#include "LR91.h"

/* Voltage */
double v; // Membrane voltage (mV)
double vnew; // New Voltage (mV)
double dvdt; // Change in Voltage / Change in Time (mV/ms)
double dvdtnew; // New dv/dt (mV/ms)

/* Time Step */
double dt; // Time step (ms)
double t; // Time (ms)
double udt; // Universal Time Step
int utsc; // Universal Time Step Counter
int nxstep; // Interval Between Calculating Ion Currents
int steps; // Number of Steps
int increment; // Loop Control Variable

/* Total Current and Stimulus */
double st; // Constant Stimulus (uA/cm^2)
double tstim; //Time Stimulus is Applied (ms)//Time to begin stimulus
double stimtime; //Time period during which stimulus is applied (ms)
double it; // Total current (uA/cm^2)

/* Terms for Solution of Conductance and Reversal Potential */
const double R = 8314; // Universal Gas Constant (J/kmol*K)
const double frdy = 96485; // Faraday's Constant (C/mol)
double temp = 310; // Temperature (K)

/* Ion Concentrations */
double nai; // Intracellular Na Concentration (mM)
double nao; // Extracellular Na Concentration (mM)
double cai; // Intracellular Ca Concentration (mM)
double cao; // Extracellular Ca Concentration (mM)
double ki; // Intracellular K Concentration (mM)
double ko; // Extracellular K Concentration (mM)

/* Fast Sodium Current (time dependant) */
double ina; // Fast Na Current (uA/uF)
double gna; // Max. Conductance of the Na Channel (mS/uF)
double ena; // Reversal Potential of Na (mV)
double am; // Na alpha-m rate constant (ms^-1)
double bm; // Na beta-m rate constant (ms^-1)
double ah; // Na alpha-h rate constant (ms^-1)
double bh; // Na beta-h rate constant (ms^-1)
double aj; // Na alpha-j rate constant (ms^-1)
double bj; // Na beta-j rate constant (ms^-1)
double mtau; // Na activation
double htau; // Na inactivation
double jtau; // Na inactivation
double mss; // Na activation
double hss; // Na inactivation
double jss; // Na slow inactivation
double m; // Na activation
double h; // Na inactivation
double j; // Na slow inactivation

/* Current through L-type Ca Channel */
double dcai; // Change in myoplasmic Ca concentration (mM)
double isi; // Slow inward current (uA/uF)
double esi; // Reversal Potential of si (mV)
double ad; // Ca alpha-d rate constant (ms^-1)
double bd; // Ca beta-d rate constant (ms^-1)
double af; // Ca alpha-f rate constant (ms^-1)
double bf; // Ca beta-f rate constant (ms^-1)

double d; // Voltage dependant activation gate
double dss; // Steady-state value of activation gate d
double taud; // Time constant of gate d (ms^-1)--mistake ???£¿ms£¿
double f; // Voltage dependant inactivation gate
double fss; // Steady-state value of inactivation gate f
double tauf; // Time constant of gate f (ms^-1)
double fca; // Ca dependant inactivation gate -----from LR94

/* Rapidly Activating Potassium Current */
double ikr; // Rapidly Activating K Current (uA/uF)
double gkr; // Channel Conductance of Rapidly Activating K Current (mS/uF)
double ekr; // Reversal Potential of Rapidly Activating K Current (mV)
double ax; // K alpha-x rate constant (ms^-1)
double bx; // K beta-x rate constant (ms^-1)
double xr; // Rapidly Activating K time-dependant activation  --gate X in LR91
double xrss; // Steady-state value of inactivation gate xr  --gate X in LR91
double tauxr; // Time constant of gate xr (ms^-1) --gate X in LR91
double r; // K time-independent inactivation --gate Xi in LR91

/* Potassium Current (time-independent) */
double iki; // Time-independent K current (uA/uF)
double gki; // Channel Conductance of Time Independant K Current (mS/uF)
double eki; // Reversal Potential of Time Independant K Current (mV)
double aki; // K alpha-ki rate constant (ms^-1)
double bki; // K beta-ki rate constant (ms^-1)
double kin; // K inactivation

/* Plateau Potassium Current */
double ikp; // Plateau K current (uA/uF)
double gkp; // Channel Conductance of Plateau K Current (mS/uF)
double ekp; // Reversal Potential of Plateau K Current (mV)
double kp; // K plateau factor

/* Background Current */
double ib; // Background current (uA/uF)

//temporary gate value
double m0, h0, j_0;
double d0, f0;
double xr0;

/********************************************************/
/* Functions that describe the currents begin here */

//Fast sodium current
void comp_ina() {
	gna = 23;
	ena = ((R*temp) / frdy)*log(nao / nai);

	am = 0.32*(v + 47.13) / (1 - exp(-0.1*(v + 47.13)));
	bm = 0.08*exp(-v / 11);
	if (v < -40) {
		ah = 0.135*exp((80 + v) / -6.8);
		bh = 3.56*exp(0.079*v) + 310000 * exp(0.35*v);
		aj = (-127140 * exp(0.2444*v) - 0.00003474*exp(-0.04391*v))*((v + 37.78) / (1 + exp(0.311*(v + 79.23))));
		bj = (0.1212*exp(-0.01052*v)) / (1 + exp(-0.1378*(v + 40.14)));
	}
	else {
		ah = 0;
		bh = 1 / (0.13*(1 + exp((v + 10.66) / -11.1)));
		aj = 0;
		bj = (0.3*exp(-0.0000002535*v)) / (1 + exp(-0.1*(v + 32)));
	}
	mtau = 1 / (am + bm);
	htau = 1 / (ah + bh);
	jtau = 1 / (aj + bj);

	mss = am*mtau;
	hss = ah*htau;
	jss = aj*jtau;

	m0 = mss - (mss - m)*exp(-dt / mtau);
	h0 = hss - (hss - h)*exp(-dt / htau);
	j_0 = jss - (jss - j)*exp(-dt / jtau);

	ina = gna*m0*m0*m0*h0*j_0*(v - ena);
}

//Slow inward current
void comp_ical() {
	esi = 7.7 - 13.0287*log(cai);

	ad = 0.095*exp(-0.01*(v - 5)) / (1 + exp(-0.072*(v - 5)));
	bd = 0.07*exp(-0.017*(v + 44)) / (1 + exp(0.05*(v + 44)));
	af = 0.012*exp(-0.008*(v + 28)) / (1 + exp(0.15*(v + 28)));
	bf = 0.0065*exp(-0.02*(v + 30)) / (1 + exp(-0.2*(v + 30)));

	taud = 1 / (ad + bd);
	tauf = 1 / (af + bf);

	dss = ad*taud;
	fss = af*tauf;

	d0 = dss - (dss - d)*exp(-dt / taud);
	f0 = fss - (fss - f)*exp(-dt / tauf);

	isi = 0.09*d0*f0*(v - esi);

	dcai = -0.0001*isi + 0.07*(0.0001 - cai);

	//if (ca_change == 1){
	//	cai = cai + dcai*dt;
	//}
}

//Time-dependent potassium current, Ik
void comp_ikr() {
	gkr = 0.282*sqrt(ko / 5.4);
	ekr = ((R*temp) / frdy)*log(ko / ki);

	ax = 0.0005*exp(0.083*(v + 50)) / (1 + exp(0.057*(v + 50)));
	bx = 0.0013*exp(-0.06*(v + 20)) / (1 + exp(-0.04*(v + 20)));

	tauxr = 1 / (ax + bx);
	xrss = ax*tauxr;
	xr0 = xrss - (xrss - xr)*exp(-dt / tauxr);

	if (v > -100) {
		r = 2.837*(exp(0.04*(v + 77)) - 1) / ((v + 77)*exp(0.04*(v + 35)));
	}else {
		r = 1;
	}

	ikr = gkr*xr0*r*(v - ekr);
}

//Time-independent potassium current
void comp_iki() {
	gki = 0.6047*(sqrt(ko / 5.4));
	eki = ((R*temp) / frdy)*log(ko / ki);

	aki = 1.02 / (1 + exp(0.2385*(v - eki - 59.215)));
	bki = (0.49124*exp(0.08032*(v - eki + 5.476)) + exp(0.06175*(v - eki - 594.31))) / (1 + exp(-0.5143*(v - eki + 4.753)));
	kin = aki / (aki + bki);

	iki = gki*kin*(v - eki);
}

//Plateau potassium current
void comp_ikp() {
	gkp = 0.0183;
	ekp = eki;

	kp = 1 / (1 + exp((7.488 - v) / 5.98));

	ikp = gkp*kp*(v - ekp);
}

//Background current
void comp_ib() {
	ib = 0.03921*(v + 59.87);
}

/* Total sum of currents is calculated here, if the time is between
stimtime = 0 and stimtime = 0.5 (ms), a stimulus is applied */
void comp_it() {
	if (t >= tstim) { //the initial value of tstim is 10.0ms, time to apply stimulus 
		stimtime = 0;
		tstim = tstim + bcl;//to next stimulation cycle
	}

	//duration of stimulus: 0-0.5, initial value of stimtime = 10.0
	if (stimtime >= 0 && stimtime <= 0.5) {
		it = st + ina + isi + ikr + iki + ikp + ib;
	}
	else {
		it = ina + isi + ikr + iki + ikp + ib;
	}
}
