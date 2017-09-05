
#ifndef _numerical_methods_H
#define _numerical_methods_H
/* Numerical method  */
extern double dv_min;//mv
extern double dv_max;//mv
extern double dt_max;//ms
extern double dt_min;//ms
extern double dt_min;//ms
extern double old_dt, new_dt;

/***for CCL--step change***/
extern double Dv1_tmp, Dv2, dt_tmp;
extern double Voffset;//mv, a given membrane potential offset
extern double dt_univ;//ms universal maximum time-step siz

extern double f_dxy(double x, double y);  // dy/dx=f(x,y)
extern double get_init_it(); //get_init_Current when dt = 0
extern void new_gate(); // renew gate value when t+dt

extern double Forward_Euler();
extern double Backward_Euler();
extern double Trapezoid();
extern double Improved_Euler();
extern double RK_method();
extern void Hybrid_method(double *vnew, double *new_dt, int cycle);
extern void CCL(double *vnew, double *new_dt);
#endif