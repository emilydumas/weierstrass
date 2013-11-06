/*********************************************************************
 * Filename:      weierstrass.c
 * Description:   the Weierstrass P function
 * Author:        David Dumas <david@dumas.io>
 * Modified at:   Wed Nov  6 11:55:23 2013
 *                
 * Copyright (C) 2013  David Dumas
 *                
 * This program is free software distributed under the GNU General
 * Public License.  See the file COPYING for details.
 ********************************************************************/

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <stdio.h>

#define THETA_ITER_MAX 200

/* (4/3) * pi^4 */
#define _CONST_43PI4 129.8787880453365829819204435849401483329701142

/* 8 pi^6 / 27 */
#define _CONST_827PI6 284.85605735564575912006502034145774781250889720920

/* 4 pi^6 / 9 */
#define _CONST_49PI6 427.28408603346863868009753051218662171876334581380

/* 1/28 */
#define _CONST_1_28 0.035714285714285714285714285714285714285714285714286

/* For printf debugging (!) */
void show(const char *name, gsl_complex value)
{
  printf("%s = %f + I %f\n",name,GSL_REAL(value),GSL_IMAG(value));
}

gsl_complex pow4(gsl_complex x)
{
  return gsl_complex_mul(gsl_complex_mul(x,x),gsl_complex_mul(x,x));
}

gsl_complex theta10(gsl_complex q, gsl_complex q14)
{
  return gsl_complex_rect(0.0,0.0);
}

gsl_complex theta20(gsl_complex q, gsl_complex q14)
{
  int n=0;
  gsl_complex accum = gsl_complex_rect(0.0,0.0);
  gsl_complex q2 = gsl_complex_mul(q,q);
  gsl_complex nextm = q2;
  gsl_complex qpower = gsl_complex_rect(1.0,0.0);

  while ((gsl_complex_abs(qpower) > 2.0*GSL_DBL_EPSILON) && (n < THETA_ITER_MAX)) {
    accum = gsl_complex_add(accum, qpower);
    qpower = gsl_complex_mul(qpower, nextm);
    nextm = gsl_complex_mul(nextm, q2);
    n++;
  }
  if (n >= THETA_ITER_MAX)
    return(gsl_complex_rect(0.0,0.0));
  return gsl_complex_mul_real(gsl_complex_mul(q14,accum),2.0);
}


gsl_complex theta30(gsl_complex q)
{
  int n=0;
  gsl_complex accum = gsl_complex_rect(0.5,0.0);
  gsl_complex q2 = gsl_complex_mul(q,q);
  gsl_complex nextm = gsl_complex_mul(q,q2);
  gsl_complex qpower = q;

  while ((gsl_complex_abs(qpower) > 2.0*GSL_DBL_EPSILON) && (n < THETA_ITER_MAX)) {
    accum = gsl_complex_add(accum, qpower);
    qpower = gsl_complex_mul(qpower, nextm);
    nextm = gsl_complex_mul(nextm, q2);
    n++;
  }
  if (n >= THETA_ITER_MAX)
    return(gsl_complex_rect(0.0,0.0));
  return gsl_complex_mul_real(accum,2.0);
}

gsl_complex theta40(gsl_complex q)
{
  int n=0;
  gsl_complex accum = gsl_complex_rect(0.5,0.0);
  gsl_complex q2 = gsl_complex_mul(q,q);
  gsl_complex nextm = gsl_complex_negative(gsl_complex_mul(q,q2));
  gsl_complex qpower = gsl_complex_negative(q);

  while ((gsl_complex_abs(qpower) > 2.0*GSL_DBL_EPSILON) && (n < THETA_ITER_MAX)) {
    accum = gsl_complex_add(accum, qpower);
    qpower = gsl_complex_mul(qpower, nextm);
    nextm = gsl_complex_mul(nextm, q2);
    n++;
  }
  if (n >= THETA_ITER_MAX)
    return(gsl_complex_rect(0.0,0.0));
  return gsl_complex_mul_real(accum,2.0);
}

gsl_complex theta1(gsl_complex z, gsl_complex q, gsl_complex q14)
{
  int n=0;
  gsl_complex accum = gsl_complex_rect(0.0,0.0);
  gsl_complex q2 = gsl_complex_mul(q,q);
  gsl_complex nextm = gsl_complex_negative(q2);
  gsl_complex qpower = gsl_complex_rect(1.0,0.0);
  gsl_complex term = gsl_complex_rect(1.0,0.0);

  while ((gsl_complex_abs(term) > 2.0*GSL_DBL_EPSILON) && (n < THETA_ITER_MAX)) {
    term = gsl_complex_mul(qpower, gsl_complex_sin(gsl_complex_mul_real(z,2*n+1)));
    accum = gsl_complex_add(accum, term);
    qpower = gsl_complex_mul(qpower, nextm);
    nextm = gsl_complex_mul(nextm, q2);
    n++;
  }
  if (n >= THETA_ITER_MAX)
    return(gsl_complex_rect(0.0,0.0));
  return gsl_complex_mul_real(gsl_complex_mul(q14,accum),2.0);
}

gsl_complex theta2(gsl_complex z, gsl_complex q, gsl_complex q14)
{
  int n=0;
  gsl_complex accum = gsl_complex_rect(0.0,0.0);
  gsl_complex q2 = gsl_complex_mul(q,q);
  gsl_complex nextm = q2;
  gsl_complex qpower = gsl_complex_rect(1.0,0.0);
  gsl_complex term = gsl_complex_rect(1.0,0.0);

  while ((gsl_complex_abs(term) > 2.0*GSL_DBL_EPSILON) && (n < THETA_ITER_MAX)) {
    term = gsl_complex_mul(qpower, gsl_complex_cos(gsl_complex_mul_real(z,2*n+1)));
    accum = gsl_complex_add(accum, term);
    qpower = gsl_complex_mul(qpower, nextm);
    nextm = gsl_complex_mul(nextm, q2);
    n++;
  }
  if (n >= THETA_ITER_MAX)
    return(gsl_complex_rect(0.0,0.0));
  return gsl_complex_mul_real(gsl_complex_mul(q14,accum),2.0);
}

gsl_complex theta3(gsl_complex z, gsl_complex q)
{
  int n=0;
  gsl_complex accum = gsl_complex_rect(0.5,0.0);
  gsl_complex q2 = gsl_complex_mul(q,q);
  gsl_complex nextm = gsl_complex_mul(q,q2);
  gsl_complex qpower = q;
  gsl_complex term = gsl_complex_rect(1.0,0.0);

  while ((gsl_complex_abs(qpower) > 2.0*GSL_DBL_EPSILON) && (n < THETA_ITER_MAX)) {
    term = gsl_complex_mul(qpower, gsl_complex_cos(gsl_complex_mul_real(z,2*(n+1))));
    accum = gsl_complex_add(accum, term);
    qpower = gsl_complex_mul(qpower, nextm);
    nextm = gsl_complex_mul(nextm, q2);
    n++;
  }
  if (n >= THETA_ITER_MAX)
    return(gsl_complex_rect(0.0,0.0));
  return gsl_complex_mul_real(accum,2.0);
}

gsl_complex theta4(gsl_complex z, gsl_complex q)
{
  int n=0;
  gsl_complex accum = gsl_complex_rect(0.5,0.0);
  gsl_complex q2 = gsl_complex_mul(q,q);
  gsl_complex nextm = gsl_complex_negative(gsl_complex_mul(q,q2));
  gsl_complex qpower = gsl_complex_negative(q);
  gsl_complex term = gsl_complex_rect(1.0,0.0);

  while ((gsl_complex_abs(qpower) > 2.0*GSL_DBL_EPSILON) && (n < THETA_ITER_MAX)) {
    term = gsl_complex_mul(qpower, gsl_complex_cos(gsl_complex_mul_real(z,2*(n+1))));
    accum = gsl_complex_add(accum, term);
    qpower = gsl_complex_mul(qpower, nextm);
    nextm = gsl_complex_mul(nextm, q2);
    n++;
  }
  if (n >= THETA_ITER_MAX)
    return(gsl_complex_rect(0.0,0.0));
  return gsl_complex_mul_real(accum,2.0);
}

/*
   In:  tau (lattice parameter)
   Out: g2 -> g[0]
        g3 -> g[1]
*/
void compute_invariants(gsl_complex tau, gsl_complex *g)
{
  gsl_complex q, q14;
  gsl_complex t2,t3,t24,t34;
  gsl_complex g3_term1, g3_term2;
  gsl_complex g2, g3;

  q = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI));
  q14 = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI_4));

  t2=theta20(q,q14);
  t3=theta30(q);
  t24 = pow4(t2);
  t34 = pow4(t3);

  g2 = gsl_complex_mul_real(gsl_complex_sub(gsl_complex_add(gsl_complex_mul(t24,t24),gsl_complex_mul(t34,t34)),gsl_complex_mul(t24,t34)),_CONST_43PI4);

  g3_term1 = gsl_complex_add(gsl_complex_mul(t24,gsl_complex_mul(t24,t24)),gsl_complex_mul(t34,gsl_complex_mul(t34,t34)));
  
  g3_term2 = gsl_complex_mul(gsl_complex_add(t24,t34),gsl_complex_mul(t24,t34));

  g3 = gsl_complex_sub( gsl_complex_mul_real(g3_term1, _CONST_827PI6),
			gsl_complex_mul_real(g3_term2, _CONST_49PI6) );

  g[0] = g2;
  g[1] = g3;
}


/* The Lattes map */
gsl_complex P_doubler(gsl_complex p, const gsl_complex *g)
{
  gsl_complex p2, p3;
  gsl_complex num;
  gsl_complex denom;
  gsl_complex term;

  p2 = gsl_complex_mul(p,p);
  p3 = gsl_complex_mul(p2,p);

  /* denom = 4p^3 - g2p - g3 */
  denom = gsl_complex_sub(gsl_complex_mul_real(p3,4.0),
			  gsl_complex_add(gsl_complex_mul(p,g[0]),g[1]));

  /* num = (p^2 + g2/4)^2 + 2g3p */
  term = gsl_complex_add(p2,gsl_complex_mul_real(g[0],0.25));
  num = gsl_complex_add(gsl_complex_mul(p,gsl_complex_mul_real(g[1],2.0)),
			gsl_complex_mul(term,term));

  return gsl_complex_div(num,denom);
}

/* Compute P using CGL/Lattes iteration */
gsl_complex wP(gsl_complex z, const gsl_complex *g)
{
  int N = 6;
  int i;
  gsl_complex z0;
  gsl_complex z02;
  gsl_complex p;

  z0 = gsl_complex_div_real(z,(double)(1 << N));
  z02 = gsl_complex_mul(z0,z0);

  /* Laurent expansion:  P \approx 1/z^2 + (g2/20)z^2 + (g3/28) z^4 */
  p = gsl_complex_add(gsl_complex_inverse(z02),
		      gsl_complex_add(gsl_complex_mul(z02,gsl_complex_mul_real(g[0],0.05)),
				      gsl_complex_mul(gsl_complex_mul(z02,z02),gsl_complex_mul_real(g[1],_CONST_1_28))));

  for (i=0;i<N;i++) {
    p = P_doubler(p,g);
  }

  return p;
}



/* ---------------------------------------------------------------------- */
/*                   UNUSED FUNCTIONS BELOW                               */
/* ---------------------------------------------------------------------- */

/*
Both P(z,tau) and P'(z,tau) can be written in the form

   a(tau) + b(tau) * M(z,tau)

where a,b are complex constants depending only on the lattice, and
where M is a laurent monomial in theta functions of pi*z.

This function computes the constants a1,b1 (for P) and b2 (for P').
Since a2=0 this constant is not included.
*/

void compute_lattice_coefs(gsl_complex tau, gsl_complex *a1, gsl_complex *b1, gsl_complex *b2)
{
  gsl_complex q, q14;
  gsl_complex t2,t3,t4,t24,t34,tprod;
  gsl_complex g2, g3;
  gsl_complex g3_term1, g3_term2;
  gsl_complex alpha, beta, gamma, Delta, sDelta;
  gsl_complex e1;

  q = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI));
  q14 = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI_4));
 
  t2=theta20(q,q14);
  t3=theta30(q);
  t4=theta40(q);

  t24 = pow4(t2);
  t34 = pow4(t3);
  tprod = gsl_complex_mul(t2,gsl_complex_mul(t3,t4));

  g2 = gsl_complex_mul_real(gsl_complex_sub(gsl_complex_add(gsl_complex_mul(t24,t24),gsl_complex_mul(t34,t34)),gsl_complex_mul(t24,t34)),_CONST_43PI4);

  g3_term1 = gsl_complex_add(gsl_complex_mul(t24,gsl_complex_mul(t24,t24)),gsl_complex_mul(t34,gsl_complex_mul(t34,t34)));
  
  g3_term2 = gsl_complex_mul(gsl_complex_add(t24,t34),gsl_complex_mul(t24,t34));

  g3 = gsl_complex_sub( gsl_complex_mul_real(g3_term1, _CONST_827PI6),
			gsl_complex_mul_real(g3_term2, _CONST_49PI6) );


  Delta = gsl_complex_sub(gsl_complex_mul(g3,g3),gsl_complex_mul_real(gsl_complex_mul(g2,gsl_complex_mul(g2,g2)),1.0/27.0));

  sDelta = gsl_complex_sqrt(Delta);
  /* if (GSL_IMAG(sDelta) < 0.0) */
  /*   sDelta = gsl_complex_negative(sDelta); */

  alpha = gsl_complex_sqrt(gsl_complex_mul_real(g2,1.0/12.0));
  /* if (GSL_IMAG(alpha) < 0.0) */
  /*   alpha = gsl_complex_negative(alpha); */

  beta = gsl_complex_div(gsl_complex_pow_real(gsl_complex_add(g3,sDelta),1.0/3.0),gsl_complex_mul_real(alpha,2.0));

  gamma = gsl_complex_inverse(beta);

  *a1 = gsl_complex_mul(alpha,gsl_complex_add(beta,gamma));
  *b1 = gsl_complex_mul_real(gsl_complex_mul(gsl_complex_mul(t3,t3),gsl_complex_mul(t4,t4)),M_PI*M_PI);
  *b2 = gsl_complex_mul_real(gsl_complex_mul(tprod,tprod),-2.0*M_PI*M_PI*M_PI);

  /* Problems with selecting the wrong root of the Weierstrass polynomial.  Debug... */
  show("tau",tau);
  show("g2",g2);
  show("g3",g3);
  show("Delta",Delta);
  show("sDelta",sDelta);
  show("alpha",alpha);
  show("beta",beta);
  show("gamma",gamma);
  show("a1",*a1);
  show("b1",*b1);
  show("b2",*b2);
}

/* Compute P using theta functions, where lattice coefficients a,b are known */
gsl_complex wP_theta_lc(gsl_complex z, gsl_complex tau, gsl_complex a1, gsl_complex b1)
{
  gsl_complex q,q14;
  gsl_complex t1,t2;

  q = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI));
  q14 = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI_4));
 
  t1 = theta1(gsl_complex_mul_real(z,M_PI),q,q14);
  t2 = theta2(gsl_complex_mul_real(z,M_PI),q,q14);

  return gsl_complex_add(a1,gsl_complex_mul(b1,gsl_complex_div(gsl_complex_mul(t2,t2),gsl_complex_mul(t1,t1))));
}

/* Compute P' using theta functions, where lattice coefficients a,b are known */
gsl_complex wP_prime_theta_lc(gsl_complex z, gsl_complex tau, gsl_complex b2)
{
  gsl_complex q,q14;
  gsl_complex t1,t2,t3,t4;
  gsl_complex piz;

  q = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI));
  q14 = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI_4));
 
  piz = gsl_complex_mul_real(z,M_PI);

  t1 = theta1(piz,q,q14);
  t2 = theta2(piz,q,q14);
  t3 = theta3(piz,q);
  t4 = theta4(piz,q);

  return gsl_complex_mul(b2,gsl_complex_div(gsl_complex_mul(t2,gsl_complex_mul(t3,t4)), gsl_complex_mul(t1,gsl_complex_mul(t1,t1))));
}

