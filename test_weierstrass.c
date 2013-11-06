/*********************************************************************
 * Filename:      test_weierstrass.c
 * Description:   test Weierstrass implementation
 * Author:        David Dumas <david@dumas.io>
 * Modified at:   Wed Nov  6 13:55:38 2013
 *                
 * Copyright (C) 2013  David Dumas
 *                
 * This program is free software distributed under the GNU General
 * Public License.  See the file COPYING for details.
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "weierstrass.h"

int main(int argc, char *argv[])
{
  int i;
  int j;
  gsl_complex tau_tbl[4] = {{0.0, 1.0},
			    {0.2, 0.3},
			    {0.3, 1.0},
			    {-0.5, 0.1}};
  gsl_complex q,z,tau;
  gsl_complex q14;
  gsl_complex v1,v2,v3,v4;
  gsl_complex v20,v30,v40;
  gsl_complex p,pp,g[2];

  gsl_complex z_tbl[7] = { {0.0, 0.0},
			   {0.5000000000, 0.8660254038},
			   { 0.5000000000, 0.0000000000},
			   { 0.9985675212, -0.0112661793},
			   { 0.4347361266, 0.3543146131},
			   { 3.0, -4.0},
			   { 0.99, 0.99} };

  /* theta tests at z=0 */
  for (i=0;i<4;i++) {
    q = gsl_complex_exp(gsl_complex_mul_imag(tau_tbl[i],M_PI));
    q14 = gsl_complex_exp(gsl_complex_mul_imag(tau_tbl[i],M_PI_4));

    v20 = theta20(q,q14);
    v30 = theta30(q);
    v40 = theta40(q);

    printf("q=(%f,%f) theta20=(%f,%f) theta30=(%f,%f) theta40=(%f,%f)\n",
	   GSL_REAL(q),GSL_IMAG(q),
	   GSL_REAL(v20),GSL_IMAG(v20),
	   GSL_REAL(v30),GSL_IMAG(v30),
	   GSL_REAL(v40),GSL_IMAG(v40));
  }

  /* theta tests at general z */
  for (i=0;i<4;i++) {
    for (j=0;j<6;j++) {
      z = z_tbl[j];
      q = gsl_complex_exp(gsl_complex_mul_imag(tau_tbl[i],M_PI));
      q14 = gsl_complex_exp(gsl_complex_mul_imag(tau_tbl[i],M_PI_4));
      
      v1 = theta1(z,q,q14);
      v2 = theta2(z,q,q14);
      v3 = theta3(z,q);
      v4 = theta4(z,q);
      
      printf("z=(%f,%f) q=(%f,%f) theta1=(%f,%f) theta2=(%f,%f) theta3=(%f,%f) theta4=(%f,%f)\n",
	     GSL_REAL(z),GSL_IMAG(z),
	     GSL_REAL(q),GSL_IMAG(q),
	     GSL_REAL(v1),GSL_IMAG(v1),
	     GSL_REAL(v2),GSL_IMAG(v2),
	     GSL_REAL(v3),GSL_IMAG(v3),
	     GSL_REAL(v4),GSL_IMAG(v4));
    }
  }

  /* P tests */
  for (i=0;i<4;i++) {
    for (j=1;j<7;j++) {
      z = z_tbl[j];
      tau = tau_tbl[i];
      compute_invariants(tau,g);
      p = wP(z,tau,g);
      printf("z=(%f,%f) tau=(%f,%f) p=(%.20f,%.20f)\n",
	     GSL_REAL(z),GSL_IMAG(z),
	     GSL_REAL(tau),GSL_IMAG(tau),
	     GSL_REAL(p),GSL_IMAG(p));
    }
  }

  /* P and P' tests */
  for (i=0;i<4;i++) {
    for (j=1;j<7;j++) {
      z = z_tbl[j];
      tau = tau_tbl[i];
      compute_invariants(tau,g);
      wP_and_prime(z,tau,g,&p,&pp);
      printf("z=(%f,%f) tau=(%f,%f) p=(%.20f,%.20f) pp=(%.20f,%.20f)\n",
	     GSL_REAL(z),GSL_IMAG(z),
	     GSL_REAL(tau),GSL_IMAG(tau),
	     GSL_REAL(p),GSL_IMAG(p),
	     GSL_REAL(pp),GSL_IMAG(pp));
    }
  }

  /* P and P' tests, grid for square torus */
  tau = gsl_complex_rect(0.0,1.0);
  compute_invariants(tau,g);  
  for (i=0;i<19;i++) {
    for (j=0;j<19;j++) {
      z = gsl_complex_rect((j+1)*0.05,(i+1)*0.05);
      wP_and_prime(z,tau,g,&p,&pp);
      printf("%.20f %.20f %.20f %.20f %.20f %.20f\n",
	     GSL_REAL(z),GSL_IMAG(z),
	     GSL_REAL(p),GSL_IMAG(p),
	     GSL_REAL(pp),GSL_IMAG(pp));
    }
  }
}
