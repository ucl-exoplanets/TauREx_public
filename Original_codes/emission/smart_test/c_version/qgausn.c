/* qgausn.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static logical c_true = TRUE_;

/* Subroutine */ int qgausn_(integer *m, doublereal *gmu, doublereal *gwt)
{
    /* Initialized data */

    static doublereal pi = 0.;
    static integer maxit = 1000;
    static doublereal one = 1.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double asin(doublereal), tan(doublereal), cos(doublereal);

    /* Local variables */
    static integer k;
    static doublereal p, t, x, en;
    static integer nn;
    static doublereal xi, pm1;
    static integer np1;
    static doublereal pm2;
    static integer lim;
    static doublereal tol, tmp, ppr, nnp1, cona;
    static integer iter;
    static doublereal prod, p2pri;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

/*       Compute weights and abscissae for ordinary gaussian quadrature */
/*       (no weight function inside integral) on the interval (0,1) */
/*   INPUT :    M                     order of quadrature rule */
/*   OUTPUT :  GMU(I)  I = 1 TO M,    array of abscissae */
/*             GWT(I)  I = 1 TO M,    array of weights */
/*   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical */
/*                   Integration, Academic Press, New York, pp. 87, 1975. */
/*   METHOD:  Compute the abscissae as roots of the Legendre */
/*            polynomial P-sub-M using a cubically convergent */
/*            refinement of Newton's method.  Compute the */
/*            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note */
/*            that Newton's method can very easily diverge; only a */
/*            very good initial guess can guarantee convergence. */
/*            The initial guess used here has never led to divergence */
/*            even for M up to 1000. */
/*   ACCURACY:  at least 13 significant digits */
/*   INTERNAL VARIABLES: */
/*    ITER      : number of Newton Method iterations */
/*    MAXIT     : maximum allowed iterations of Newton Method */
/*    PM2,PM1,P : 3 successive Legendre polynomials */
/*    PPR       : derivative of Legendre polynomial */
/*    P2PRI     : 2nd derivative of Legendre polynomial */
/*    TOL       : convergence criterion for Legendre poly root iteration */
/*    X,XI      : successive iterates in cubically-convergent version */
/*                of Newtons Method (seeking roots of Legendre poly.) */
/* +---------------------------------------------------------------------+ */
    /* Parameter adjustments */
    --gwt;
    --gmu;

    /* Function Body */
    if (pi == 0.f || tol == 0.f) {
	pi = asin(1.f) * 2.f;
	tol = d1mach_(&c__4) * 10.f;
    }
    if (*m < 1) {
	errmsg_("QGAUSN--Bad value of M", &c_true, (ftnlen)22);
    }
    if (*m == 1) {
	gmu[1] = .5f;
	gwt[1] = 1.f;
	return 0;
    }
    en = (doublereal) (*m);
    np1 = *m + 1;
    nnp1 = (doublereal) (*m * np1);
/* Computing 3rd power */
    i__1 = *m;
    cona = (real) (*m - 1) / (i__1 * (i__1 * i__1) << 3);
    lim = *m / 2;
    i__1 = lim;
    for (k = 1; k <= i__1; ++k) {
/*                                        ** initial guess for k-th root */
/*                                        ** of Legendre polynomial, from */
/*                                        ** Davis/Rabinowitz (2.7.3.3a) */
	t = ((k << 2) - 1) * pi / ((*m << 2) + 2);
	x = cos(t + cona / tan(t));
	iter = 0;
/*                                        ** upward recurrence for */
/*                                        ** Legendre polynomials */
L10:
	++iter;
	pm2 = one;
	pm1 = x;
	i__2 = *m;
	for (nn = 2; nn <= i__2; ++nn) {
	    p = (((nn << 1) - 1) * x * pm1 - (nn - 1) * pm2) / nn;
	    pm2 = pm1;
	    pm1 = p;
/* L20: */
	}
/*                                              ** Newton Method */
/* Computing 2nd power */
	d__1 = x;
	tmp = one / (one - d__1 * d__1);
	ppr = en * (pm2 - x * p) * tmp;
	p2pri = (two * x * ppr - nnp1 * p) * tmp;
	xi = x - p / ppr * (one + p / ppr * p2pri / (two * ppr));
/*                                              ** check for convergence */
	if ((d__1 = xi - x, abs(d__1)) > tol) {
	    if (iter > maxit) {
		errmsg_("QGAUSN--MAX ITERATION COUNT", &c_true, (ftnlen)27);
	    }
	    x = xi;
	    goto L10;
	}
/*                             ** iteration finished--calculate weights, */
/*                             ** abscissae for (-1,1) */
	gmu[k] = -x;
/* Computing 2nd power */
	d__1 = en * pm2;
	gwt[k] = two / (tmp * (d__1 * d__1));
	gmu[np1 - k] = -gmu[k];
	gwt[np1 - k] = gwt[k];
/* L30: */
    }
/*                                    ** set middle abscissa and weight */
/*                                    ** for rules of odd order */
    if (*m % 2 != 0) {
	gmu[lim + 1] = 0.f;
	prod = one;
	i__1 = *m;
	for (k = 3; k <= i__1; k += 2) {
	    prod = prod * k / (k - 1);
/* L40: */
	}
/* Computing 2nd power */
	d__1 = prod;
	gwt[lim + 1] = two / (d__1 * d__1);
    }
/*                                        ** convert from (-1,1) to (0,1) */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	gmu[k] = gmu[k] * .5f + .5f;
	gwt[k] *= .5f;
/* L50: */
    }
    return 0;
} /* qgausn_ */

