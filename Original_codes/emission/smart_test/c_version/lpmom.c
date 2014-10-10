/* lpmom.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int lpmom_(real *func, real *gwt, real *gmu, integer *numu, 
	integer *nmom, real *pmom)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double acos(doublereal);

    /* Local variables */
    static integer n;
    static real pi, pl[1001];
    static integer mom;
    static real pol, plm1[1001];


/* ccccccccccccccccccccccccccccc  l p m o m  cccccccccccccccccccccccccccccc */
/* c                                                                     cc */
/* c    p u r p o s e :                                                  cc */
/* c                                                                     cc */
/* c      compute the legendre polynomial coefficient moments for an     cc */
/* c      arbitrary function, func, that is defined only at the          cc */
/* c      gaussian angles, umu.                                          cc */
/* c                                                                     cc */
/* c     i n p u t :                                                     cc */
/* c                                                                     cc */
/* c   func( n ) :  input angle-dependent quantity at n-th angle         cc */
/* c    gmu( n ) :  cosine of nth gaussian angle.                        cc */
/* c    wgt( n ) :  gaussian weights at each angle                       cc */
/* c        numu :  number of gaussian angles.                           cc */
/* c        nmom :  number of legendre-polynomial momemts (numu - 1)     cc */
/* c                                                                     cc */
/* c                            *** specification of local variables     cc */
/* c                                                                     cc */
/* c   pl( n )   : legendre polynomial p-sub-mom  of argument  gmu(n)    cc */
/* c   plm1( n ) : legendre polynomial p-sub-(mom-1) of argument gmu(n)  cc */
/* c                                                                     cc */
/* ccccccccccccccccccccccccccccc  l p m o m  cccccccccccccccccccccccccccccc */



    /* Parameter adjustments */
    --gmu;
    --gwt;
    --func;

    /* Function Body */
    pi = acos(-1.f);

    i__1 = *nmom;
    for (mom = 0; mom <= i__1; ++mom) {

/* ****        calculate legendre polys. */

	if (mom == 0) {
	    i__2 = *numu;
	    for (n = 1; n <= i__2; ++n) {
		pl[n - 1] = 1.f;
/* L22: */
	    }
	} else if (mom == 1) {
	    i__2 = *numu;
	    for (n = 1; n <= i__2; ++n) {
		plm1[n - 1] = 1.f;
		pl[n - 1] = gmu[n];
/* L24: */
	    }
	} else {
	    i__2 = *numu;
	    for (n = 1; n <= i__2; ++n) {
		pol = (((mom << 1) - 1) * gmu[n] * pl[n - 1] - (mom - 1) * 
			plm1[n - 1]) / mom;
		plm1[n - 1] = pl[n - 1];
		pl[n - 1] = pol;
/* L26: */
	    }
	}

/* ****        get legendre moments by angular quadrature */

	pmom[mom] = 0.f;
	i__2 = *numu;
	for (n = 1; n <= i__2; ++n) {
	    pmom[mom] += gwt[n] * .5f * pl[n - 1] * func[n];
/* L30: */
	}

/* L40: */
    }

    return 0;
} /* lpmom_ */

