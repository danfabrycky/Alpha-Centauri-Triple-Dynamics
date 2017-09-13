/*  This program integrates the N-body problem and outputs the positions of the
    bodies every so often.  It uses the GSL for easy switching between 
    integrators.  This version goes both forward and backward, starting from
    the chosen epoch.  GR is included to precess the innermost orbit. 

    Dan Fabrycky  9 Jan 2008
                  3 Mar 2010
		  23 Mar 2013
    
    on Midway: module load pgi
    compiles with: pgcc -O3 -o tripledyn -lgsl -lgslcblas tripledyn.c

*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>


/* Integration parameters */
#define NPL 2           /* Number of planets to integrate. */
#define DY 1e-13           /* Error allowed in parameter values per timestep. */
#define DELTATPRINT 1e7  /* Print interval */
#define HStart 1e-5      /* Timestep to start.  If you get NaNs, try reducing this. */

/* Some physical constants */
#define G 39.478417604357 /*  Newton's constant, AU^3*yr^-2 */
#define CSQ 3.9992e+09   /* speed of light squared, in (AU/yr)^2 */

int
seteq (double ynew[], const double y[])
{
  int i;
  for(i=0;i<6*NPL;i++) ynew[i]=y[i];
}

int
func (double t, const double y[], double f[],
      void *params)
{
  int i,j;
  int i1,i2;
  double * masses = (double *)params;
  double gmc1, gmc2, gm1, gm2;
  double rc1m3, rc2m3, r12m3, r, eta, scale, rdot, vsq, bracket;

  for(i=0; i<NPL; i++) {
    gmc1 = G*(masses[0]+masses[i+1]);
    rc1m3 = pow(pow(y[i*6+0],2)+pow(y[i*6+1],2)+pow(y[i*6+2],2),-3.0/2);
    if(i==0) {
        /* GR, according to Kidder via Fabrycky review (Non-Kep) */
       r = pow(rc1m3,-1./3);
       eta = masses[0]*masses[1]/pow(masses[0]+masses[1],2);
       scale = - G*(masses[0]+masses[1])*pow(rc1m3,2./3)/CSQ;
       rdot = (y[i*6+0]*y[i*6+3]+y[i*6+1]*y[i*6+4]+y[i*6+2]*y[i*6+5]) / r;
       vsq = pow(y[i*6+3],2)+pow(y[i*6+4],2)+pow(y[i*6+5],2);
       bracket = (1+3*eta)*vsq - 1.5*eta*pow(rdot,2) - 2*(2+eta)*G*(masses[0]+masses[1])/r;
    }
    for(j=0; j<3; j++) {
      f[i*6+j] = y[i*6+3+j];  /* x dot = v */
      f[i*6+3+j] = -gmc1*y[i*6+j]*rc1m3;   /* Pull of the star. */
      if(i==0) { 
	/* GR, according to Kidder via Fabrycky review (Non-Kep) */
	f[i*6+3+j] += scale*( -2*(2-eta)*rdot*y[i*6+3+j] + bracket*y[i*6+j]/r );
      } 

    }   
  }
 
  /* Interaction between each pair of planets. */
  for(i1=0; i1<NPL-1; i1++) {
    gm1 = G*masses[i1+1];
    rc1m3 = pow(pow(y[i1*6+0],2)+pow(y[i1*6+1],2)+pow(y[i1*6+2],2),-3.0/2);
    for(i2=i1+1; i2<NPL; i2++) {
      gm2 = G*masses[i2+1];
      rc2m3 = pow(pow(y[i2*6+0],2)+pow(y[i2*6+1],2)+pow(y[i2*6+2],2),-3.0/2);
      r12m3 = pow(pow(y[i1*6+0]-y[i2*6+0],2)+pow(y[i1*6+1]-y[i2*6+1],2)+pow(y[i1*6+2]-y[i2*6+2],2),-3.0/2);
  
      for(j=0; j<3; j++) f[i1*6+3+j] += -gm2*( (y[i1*6+j]-y[i2*6+j])*r12m3 + y[i2*6+j]*rc2m3 );
      for(j=0; j<3; j++) f[i2*6+3+j] += -gm1*( (y[i2*6+j]-y[i1*6+j])*r12m3 + y[i1*6+j]*rc1m3 );
    }
  }

  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{
  /* This is a part of the code that 
  isn't used for the Alpha Centauri project; 
  it would need to be checked. */
  
  int i,j;
  double * masses = (double *)params;
  double gmc1 = G*(masses[0]+masses[1]);
  double gmc2 = G*(masses[0]+masses[2]);
  double gm1 = G*masses[1];
  double gm2 = G*masses[2];

  double rc1m3 = pow(pow(y[0],2)+pow(y[1],2)+pow(y[2],2),-3.0/2);
  double rc2m3 = pow(pow(y[6],2)+pow(y[7],2)+pow(y[8],2),-3.0/2);
  double r12m3 = pow(pow(y[0]-y[6],2)+pow(y[1]-y[7],2)+pow(y[2]-y[8],2),-3.0/2);
  double rc1m5 = pow(pow(y[0],2)+pow(y[1],2)+pow(y[2],2),-5.0/2);
  double rc2m5 = pow(pow(y[6],2)+pow(y[7],2)+pow(y[8],2),-5.0/2);
  double r12m5 = pow(pow(y[0]-y[6],2)+pow(y[1]-y[7],2)+pow(y[2]-y[8],2),-5.0/2);

  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 6*NPL, 6*NPL);
  gsl_matrix * m = &dfdy_mat.matrix; 

  /* Planet 1. */
  for(i=0;i<=2;i++){
    for(j=0;j<=2;j++) gsl_matrix_set (m, i, j, 0.0);
    for(j=3;j<=5;j++) gsl_matrix_set (m, i, j, 1.0);
    for(j=6;j<=11;j++) gsl_matrix_set (m, i, j, 0.0);
  }
  for(i=3;i<=5;i++){
    for(j=0;j<=2;j++){
      if(i==j+3) {
	gsl_matrix_set (m, i, j, -gmc1*rc1m3 + 3*gmc1*y[j]*y[j]*rc1m5 - gmc2*( r12m3 - 3*y[j]*(y[j]-y[j+6])*r12m5 ) );
      } else {
	gsl_matrix_set (m, i, j, 3*gmc1*y[j]*y[i-3]*rc1m5 + gmc2*3*y[j]*(y[i-3]-y[i+3])*r12m5 );
      }
    }
    for(j=3;j<=5;j++) gsl_matrix_set (m, i, j, 0.0);
    for(j=6;j<=8;j++){
      if(i==j-3) {
	gsl_matrix_set (m, i, j, -gmc2*( -r12m3 + 3*y[j]*(y[j-6]-y[j])*r12m5 + rc2m3 - 3*y[j]*y[j]*rc2m5 ) );
      } else {
	gsl_matrix_set (m, i, j, -gmc2*( 3*y[j]*(y[i-3]-y[i+3])*r12m5 ) - 3*y[j]*y[i+3]*rc2m5 );
      }
    }
    for(j=9;j<=11;j++) gsl_matrix_set (m, i, j, 0.0);
  }


  /* Planet 2. */
  for(i=6;i<=8;i++){
    for(j=0;j<=8;j++) gsl_matrix_set (m, i, j, 0.0);
    for(j=9;j<=11;j++) gsl_matrix_set (m, i, j, 1.0);
  }
  for(i=9;i<=11;i++){
    for(j=0;j<=2;j++){
      if(i==j-9) {
	gsl_matrix_set (m, i, j, -gmc1*( -r12m3 + 3*y[j]*(y[j+6]-y[j])*r12m5 + rc1m3 - 3*y[j]*y[j]*rc1m5 ) );
      } else {
	gsl_matrix_set (m, i, j, -gmc1*( 3*y[j]*(y[i-3]-y[i-9])*r12m5 ) - 3*y[j]*y[i-9]*rc1m5 );
      }
    }
    for(j=3;j<=5;j++) gsl_matrix_set (m, i, j, 0.0);
    for(j=6;j<=8;j++){
      if(i==j+3) {
	gsl_matrix_set (m, i, j, -gmc2*rc2m3 + 3*gmc2*y[j]*y[j]*rc2m5 - gmc1*( r12m3 - 3*y[j]*(y[j]-y[j-6])*r12m5 ) );
      } else {
	gsl_matrix_set (m, i, j, 3*gmc2*y[j]*y[i-3]*rc2m5 + gmc1*3*y[j]*(y[i-3]-y[i-9])*r12m5 );
      }
    }
    for(j=9;j<=11;j++) gsl_matrix_set (m, i, j, 0.0);
  }

  for(i=0;i<=11;i++) dfdt[i] = 0.0;

  return GSL_SUCCESS;
}

int
main (int argc, char *argv[])
{
  const gsl_odeiv_step_type * T 
  /*   = gsl_odeiv_step_bsimp;   14 s */
  = gsl_odeiv_step_rk8pd;   /* 3 s */
/*  = gsl_odeiv_step_rkf45;    14 s */
  /*  = gsl_odeiv_step_rk4;  26 s */
 /*   = gsl_odeiv_step_gear2;  > 60 s */
  /* Timing trials for a transit-timing variation project. Not certain this is 
   the fastest integrator (to reach a given precision) for the current purpose.  */
  
  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc (T, 6*NPL);
  gsl_odeiv_control * c 
    = gsl_odeiv_control_y_new (DY, 0.0);
  gsl_odeiv_evolve * e 
    = gsl_odeiv_evolve_alloc (6*NPL);
  
  long i;
  double mu[NPL+1];
  gsl_odeiv_system sys = {func, jac, 6*NPL, mu};

  char filename[100];
  strcpy(filename, argv[1]);

  double t, tstart; /* this is the given epoch. */
  double t0=0;   /* Just before the first RV or TTV datum. */
  double t1=12e9, tjunk;         /* t1  is time just passed final RV and TTV datum. */

  double junk, hhere;  
  double y[6*NPL], yin[6*NPL];

  FILE *f_in, *f_out, *f_rv;
  f_in = fopen(filename, "r");
  fscanf(f_in, "%ld", &junk);
  fscanf(f_in, "%lf", &tstart);
  fscanf(f_in, "%lf %lf %lf %lf %lf %lf %lf", &mu[0], &yin[0], &yin[1], &yin[2], &yin[3], &yin[4], &yin[5]);
  for(i=0;i<NPL;i++) {
    fscanf(f_in, "%lf %lf %lf %lf %lf %lf %lf", &mu[i+1], &yin[i*6+0], &yin[i*6+1], &yin[i*6+2], &yin[i*6+3], &yin[i*6+4], &yin[i*6+5]);
  }
  fclose(f_in);

  double mtot=mu[0];
  for(i=0;i<NPL;i++) mtot += mu[i+1];


  double yhone[6*NPL], dydt[6*NPL], yerrhone[6*NPL];
  double thone,tstep,dist,vtrans,lasttprint;
  double vstar,astar;


  /* Setup forward integration */
  seteq(y, yin);
  t = tstart;
  double h = HStart;

  while (t < t1 )
    {      
      int status = gsl_odeiv_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, t1,
                                           &h, y);

      if (status != GSL_SUCCESS)
          break;

      if (t>lasttprint+DELTATPRINT) {
	printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", t, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], y[11]);
	lasttprint=t;
      }

    }

   /* Setup backward integration */
  seteq(y, yin);
  t = tstart;
  h = -HStart;
  while (t > t0+1e1 )
    {      
      int status = gsl_odeiv_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, t0,
                                           &h, y);

      if (status != GSL_SUCCESS)
          break;
      if (t<lasttprint-DELTATPRINT) {
	printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", t, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], y[11]);
        lasttprint=t;
      }
    }


  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  return 0;

}
