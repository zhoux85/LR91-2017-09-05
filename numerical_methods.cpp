/* The phase 1 Luo-Rudy(LR1) ventricular model, including four kinds of Euler methods */
/* Xiang Zhou, 2017/08/30 */
/* This code requires a C++ compiler */
/* Detailed list of equations and model description are provided in */
/* Circ Res 1991;68:1501-1526 */

/* IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USER SHOULD PACE THE MODEL UNTIL STEADY-STATE IS REACHED.
THIS REQUIRES APPROXIMATELY 5-20 MINUTES OF PACING DEPENDING ON THE RATE.*/

#include "stdafx.h"

#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include "LR91.h"
#include "numerical_methods.h"

double vb, v_0, v_1, it_0, it_1;

double Forward_Euler(){
	f_count++;
	comp_ina();
	comp_ical();
	comp_ikr();
	comp_iki();
	comp_ikp();
	comp_ib();
	comp_it();
	stimtime = stimtime + dt;
	vnew = v - it*dt;
	cai = cai + dcai*dt;//update Cai
	new_gate();//update Gate values
	return vnew;
}

double Backward_Euler(){
	it_0 = get_init_it();
	stimtime = stimtime + dt;
	v_1 = v_0 - it_0*dt;
	while (1){
		v = v_1;
		it_0 = f_dxy(dt, v);
		v_1 = v_0 - it_0*dt;
		if (fabs(v_1 - v) <= 1e-6){
			vnew = v_1;
			cai = cai + dcai*dt;
			new_gate();
			break;
		}
		//printf("in while__%lf\n", v_1);
	}
	return vnew;
}

double Trapezoid(){
	//Iteration process: y(n+1)^(0),y(n+1)^(1)....y(n+1)^(k),y(n+1)^(k+1)
	it_0 = get_init_it();
	stimtime = stimtime + dt;
	v_1 = v_0 - it_0*dt;  //y(n+1)^(0)
	//it_0 = it;
	while (1){
		v = v_1;
		//ca_change = 1;
		it_1 = f_dxy(dt, v);//get it(Total current) in t+dt, using v_1. 
		vb = v_0 - (dt / 2)*(it_0 + it_1);   //y(n+1)^(1)
		if (fabs(vb - v_1) <= 1e-6){
			vnew = vb;
			cai = cai + dcai*dt;
			new_gate();
			break;
		}
		v_1 = vb;
		//printf("in while__%lf\n", vb);
	}
	return vnew;
}

double Improved_Euler(){
	v_1 = v_0 - it_0*dt;  
	v = v_1;
	it_1 = f_dxy(dt, v);
	vnew = v_0 - (dt / 2)*(it_0 + it_1);//y(n+1)^(1)
	return vnew;
	//printf("new v__%lf\n", vb);
}

//(x,y)is current point, Fourth-order Runge-Kutta method
double RK_method(){
	double K1 = get_init_it();
	double K2 = f_dxy(dt / 2, v + (dt / 2)*K1);
	double K3 = f_dxy(dt / 2, v + (dt / 2)*K2);
	double K4 = f_dxy(dt, v + dt*K3);
	f_dxy(dt, v);//get gate value and dcai when x=dt, y=v
	cai = cai + dcai*dt;//renew Cai
	new_gate();// renew gate value when t+dt
	stimtime = stimtime + dt;
	return v - (dt / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
}

void Hybrid_method(double *vnew, double *new_dt, int cycle){
	*vnew = RK_method();
	double dv = abs(*vnew - v_0);
	dvdt = dv / dt;

	/*********** step change  ***********/
	*new_dt = dt;
	if (dv <= dv_min){ //dv_min=0.05mv
		*new_dt = dv_max / dvdt;
	}else if (dv >= dv_max){//dv_max=0.2mv
		*new_dt = dv_min / dvdt;
	}
	//else{
	//	dt = dv;
	//}

	if (*new_dt>dt_max){
		*new_dt = dt_max;
	}else if (*new_dt<dt_min){
		*new_dt = dt_min;
	}
	/*A 0.5 ms protective zone is enforced on the Hybrid method, and the method chooses a
	minimum time step dt_min within the first 0.5 ms of every beat.*/
	if (t >= 10 + bcl*cycle && t <= (10 + 0.5 + bcl*cycle)){
		*new_dt = 0.001;
	}
}

void CCL(double *vnew, double *new_dt){
	*vnew = RK_method();
	double dv = abs(*vnew - v_0);
	dvdt = dv / dt;

	/*********** step change  ***********/
	*new_dt = dt;

	// then adjust or correct the time step---CCL method
	Dv2 = (dvdt - Dv1_tmp) / dt;
	Dv1_tmp = dvdt;
	double DiscriminantP = 0, DiscriminantN = 0, dtz = 0;
	if (dvdt >= 0){
		DiscriminantP = dvdt*dvdt + 2 * Dv2*Voffset;
		if (Dv2>0){
			*new_dt = (-dvdt + sqrt(DiscriminantP)) / Dv2;
		}else if (Dv2<0){
			dtz = -dvdt / Dv2;
			if (DiscriminantP >= 0){
				*new_dt = (-dvdt + sqrt(DiscriminantP)) / Dv2;
			}else{
				*new_dt = dtz;
			}
		}
	}else{
		DiscriminantN = dvdt*dvdt - 2 * Dv2*Voffset;
		if (Dv2>0){
			dtz = -dvdt / Dv2;
			if (DiscriminantN >= 0){
				*new_dt = (-dvdt - sqrt(DiscriminantP)) / Dv2;
			}else{
				*new_dt = dtz;
			}
		}else if (Dv2<0){
			*new_dt = (-dvdt - sqrt(DiscriminantP)) / Dv2;
		}
	}
	//dt_max = Minimum(2*dt_n, dt_universal)
	if (dt_univ > dt_tmp * 2){
		if (dt_tmp != 0){
			dt_max = dt_tmp * 2;
		}else{
			dt_max = dt_univ;
		}
	}else{
		dt_max = dt_univ;
	}
	if (*new_dt>dt_max){
		*new_dt = dt_max;
	}
	if (*new_dt<dt_min){
		*new_dt = dt_min;
	}
	dt_tmp = *new_dt;
}

void new_gate(){
	m = m0;
	h = h0;
	j = j_0;

	d = d0;
	f = f0;

	xr = xr0;
}

double get_init_it(){
	f_count++;
	double dt_tmp = dt;
	dt = 0;
	comp_ina();
	comp_ical();
	comp_ikr();
	comp_ikr();
	comp_iki();
	comp_ikp();
	comp_ib();
	comp_it();
	dt = dt_tmp;
	return it;
}

// dy/dx=f(x,y)
double f_dxy(double x, double y){
	f_count++;
	double dt_tmp = dt;
	double v_tmp = v;
	dt = x;
	v = y;
	comp_ina();
	comp_ical();
	comp_ikr();
	comp_iki();
	comp_ikp();
	comp_ib();
	comp_it();
	dt = dt_tmp;
	v = v_tmp;
	return it;
}

