/*********************************************************************
 * Filename:      weierstrass.h
 * Description:   the Weierstrass P function
 * Author:        David Dumas <david@dumas.io>
 * Modified at:   Wed Nov  6 17:12:40 2013
 *                
 * Copyright (C) 2013  David Dumas
 *                
 * This program is free software distributed under the GNU General
 * Public License.  See the file COPYING for details.
 ********************************************************************/

#ifndef _WEIERSTRASS_H_
#define _WEIERSTRASS_H_

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

gsl_complex theta10(gsl_complex q, gsl_complex q14);
gsl_complex theta20(gsl_complex q, gsl_complex q14);
gsl_complex theta30(gsl_complex q);
gsl_complex theta40(gsl_complex q);
gsl_complex theta1(gsl_complex z, gsl_complex q, gsl_complex q14);
gsl_complex theta2(gsl_complex z, gsl_complex q, gsl_complex q14);
gsl_complex theta3(gsl_complex z, gsl_complex q);
gsl_complex theta4(gsl_complex z, gsl_complex q);

void compute_invariants(gsl_complex tau, gsl_complex *g);

/* Core functions, with g2,g3 precomputed */
gsl_complex wP(gsl_complex z, gsl_complex tau, const gsl_complex *g);
gsl_complex wPprime(gsl_complex z, gsl_complex tau, const gsl_complex *g);
void wP_and_prime(gsl_complex z, gsl_complex tau, const gsl_complex *g, gsl_complex *p, gsl_complex *pp);

/* Core functions, without g2,g3 precomputed (i.e. _tau = "directly from tau") */
gsl_complex wP_tau(gsl_complex z, gsl_complex tau);
gsl_complex wPprime_tau(gsl_complex z, gsl_complex tau);
void wP_and_prime_tau(gsl_complex z, gsl_complex tau, gsl_complex *p, gsl_complex *pp);

/*  Old method, using theta functions */
/* void compute_lattice_coefs(gsl_complex tau, gsl_complex *a1, gsl_complex *b1, gsl_complex *b2); */
/* gsl_complex wP_theta_lc(gsl_complex z, gsl_complex tau, gsl_complex a1, gsl_complex b1); */
/* gsl_complex wP_prime_theta_lc(gsl_complex z, gsl_complex tau, gsl_complex b2); */


#endif  /* ifndef _WEIERSTRASS_H_ */
