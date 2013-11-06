/*********************************************************************
 * Filename:      weierstrass.c
 * Description:   the Weierstrass P function
 * Author:        David Dumas <david@dumas.io>
 * Modified at:   Wed Nov  6 17:10:38 2013
 *                
 * Copyright (C) 2013  David Dumas
 *                
 * This program is free software distributed under the GNU General
 * Public License.  See the file COPYING for details.
 ********************************************************************/

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#define THETA_ITER_MAX 200

/* (4/3) * pi^4 */
#define _CONST_43PI4 129.8787880453365829819204435849401483329701142

/* 8 pi^6 / 27 */
#define _CONST_827PI6 284.85605735564575912006502034145774781250889720920

/* 4 pi^6 / 9 */
#define _CONST_49PI6 427.28408603346863868009753051218662171876334581380

/* 1/28 */
#define _CONST_1_28 0.035714285714285714285714285714285714285714285714286

/* 1/7 */
#define _CONST_1_7 0.14285714285714285714285714285714285714285714285714

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

/* The extended Lattes map (rational function doubling on the elliptic curve) */
void P_and_Pprime_doubler(gsl_complex *p, gsl_complex *pp, const gsl_complex *g)
{
  gsl_complex pp3;
  gsl_complex ppp, ppp3;


  /* p'' */
  ppp = gsl_complex_sub(gsl_complex_mul_real(gsl_complex_mul(*p,*p),6.0),
			gsl_complex_mul_real(g[0],0.5));
  
  ppp3 = gsl_complex_mul(ppp,gsl_complex_mul(ppp,ppp));
  pp3 = gsl_complex_mul(*pp,gsl_complex_mul(*pp,*pp));

  
  *pp = gsl_complex_sub(gsl_complex_add(gsl_complex_mul_real(gsl_complex_div(gsl_complex_mul(*p,ppp),*pp),3.0),
					gsl_complex_mul_real(gsl_complex_div(ppp3,pp3),-0.25)),
			*pp);
  *p = P_doubler(*p,g);
}

/* Assuming z is in the (1,tau) parallelogram, return the point
   closest to the origin among all translates of z by the lattice. */
gsl_complex near_origin(gsl_complex z, gsl_complex tau)
{
  int i;
  gsl_complex znew;

  znew = gsl_complex_sub_real(z,1.0);
  if (gsl_complex_abs(z) > gsl_complex_abs(znew))
    z = znew;

  znew = gsl_complex_sub(z,tau);
  if (gsl_complex_abs(z) > gsl_complex_abs(znew))
    z = znew;

  znew = gsl_complex_sub(z,gsl_complex_add_real(tau,1.0));
  if (gsl_complex_abs(z) > gsl_complex_abs(znew))
    z = znew;

  return z;
}

/* Compute P using CGL/Lattes iteration */
/* NOTE: Assumes z is in fundamental parallelogram  */
gsl_complex wP(gsl_complex z, gsl_complex tau, const gsl_complex *g)
{
  int N = 6;
  int i;
  gsl_complex z0;
  gsl_complex z02;
  gsl_complex p;

  z = near_origin(z,tau);

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

/* Compute P and P' using CGL/Lattes iteration */
/* NOTE: Assumes z is in fundamental parallelogram  */
void wP_and_prime(gsl_complex z, gsl_complex tau, const gsl_complex *g, gsl_complex *p, gsl_complex *pp)
{
  int N = 6;  /* Enough iterations for good P, not so good P' */
  int i;
  int flip=0;
  gsl_complex z0;
  gsl_complex z02;
  gsl_complex pout, ppout;
  gsl_complex ppsolve;

  z = near_origin(z,tau);

  z0 = gsl_complex_div_real(z,(double)(1 << N));
  z02 = gsl_complex_mul(z0,z0);

  /* Laurent expansion:  P \approx 1/z^2 + (g2/20)z^2 + (g3/28) z^4 */
  pout = gsl_complex_add(gsl_complex_inverse(z02),
			 gsl_complex_add(gsl_complex_mul(z02,gsl_complex_mul_real(g[0],0.05)),
					 gsl_complex_mul(gsl_complex_mul(z02,z02),gsl_complex_mul_real(g[1],_CONST_1_28))));

  /* Laurent expansion:  P' \approx -2/z^3 + g2/10z + g3/7 z^3 */
  ppout = gsl_complex_add(gsl_complex_mul_real(gsl_complex_inverse(gsl_complex_mul(z0,z02)),-2.0),
			  gsl_complex_add(gsl_complex_mul(z0,gsl_complex_mul_real(g[0],0.1)),
					  gsl_complex_mul(gsl_complex_mul(z0,z02),gsl_complex_mul_real(g[1],_CONST_1_7))));

  for (i=0;i<N;i++) {
    P_and_Pprime_doubler(&pout, &ppout, g);
  }

  /* At this point ppout is a decent but not great approximation of P'(z)        */
  /* Instead of using it directly, we use it as a guide for which square root of */
  /* (4P^3 - g2 P - g3) should be selected.                                      */

  ppsolve = gsl_complex_sqrt(
                gsl_complex_sub(
                    gsl_complex_mul_real(gsl_complex_mul(pout,gsl_complex_mul(pout,pout)),4.0),
		    gsl_complex_add(gsl_complex_mul(g[0],pout),g[1])
                )
	    );

  *p = pout;
  if (gsl_complex_abs(gsl_complex_sub(ppsolve,ppout)) < gsl_complex_abs(gsl_complex_add(ppsolve,ppout)))
    *pp = ppsolve;
  else
    *pp = gsl_complex_negative(ppsolve);
}

/* Compute P using CGL/Lattes iteration */
/* NOTE: Assumes z is in fundamental parallelogram  */
gsl_complex wPprime(gsl_complex z, gsl_complex tau, const gsl_complex *g)
{
  gsl_complex p,pp;
  
  wP_and_prime(z,tau,g,&p,&pp);
  return pp;
}

/* Compute P directly from tau */
/* For speed, should compute and store g2,g3 if making many calls with same tau */
gsl_complex wP_tau(gsl_complex z, gsl_complex tau)
{
  gsl_complex g[2];

  compute_invariants(tau,g);
  return wP(z,tau,g);
}

/* Compute P' directly from tau */
/* For speed, should compute and store g2,g3 if making many calls with same tau */
gsl_complex wPprime_tau(gsl_complex z, gsl_complex tau)
{
  gsl_complex g[2];

  compute_invariants(tau,g);
  return wPprime(z,tau,g);
}

/* Compute P' directly from tau */
/* For speed, should compute and store g2,g3 if making many calls with same tau */
void wP_and_prime_tau(gsl_complex z, gsl_complex tau, gsl_complex *p, gsl_complex *pp)
{
  gsl_complex g[2];

  compute_invariants(tau,g);
  wP_and_prime(z,tau,g,p,pp);
}
