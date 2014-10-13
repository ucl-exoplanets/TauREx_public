/* lepoly.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int lepoly_(integer *nmu, integer *m, integer *maxmu, 
	integer *twonm1, real *mu, real *ylm)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    integer ylm_dim1, ylm_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, l, ns;
    static real sqt[1000], tmp1, tmp2;

/*       COMPUTES THE NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL, */
/*       DEFINED IN TERMS OF THE ASSOCIATED LEGENDRE POLYNOMIAL */
/*       PLM = P-SUB-L-SUPER-M AS */
/*             YLM(MU) = SQRT( (L-M)!/(L+M)! ) * PLM(MU) */
/*       FOR FIXED ORDER -M- AND ALL DEGREES FROM L = M TO TWONM1. */
/*       WHEN M.GT.0, ASSUMES THAT Y-SUB(M-1)-SUPER(M-1) IS AVAILABLE */
/*       FROM A PRIOR CALL TO THE ROUTINE. */
/*       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of */
/*                  High-Order Associated Legendre Polynomials, */
/*                  J. Quant. Spectrosc. Radiat. Transfer 10, */
/*                  557-562, 1970.  (hereafter D/A) */
/*       METHOD: Varying degree recurrence relationship. */
/*       NOTE 1: The D/A formulas are transformed by */
/*               setting  M = n-1; L = k-1. */
/*       NOTE 2: Assumes that routine is called first with  M = 0, */
/*               then with  M = 1, etc. up to  M = TWONM1. */
/*       NOTE 3: Loops are written in such a way as to vectorize. */
/*  I N P U T     V A R I A B L E S: */
/*       NMU    :  NUMBER OF ARGUMENTS OF -YLM- */
/*       M      :  ORDER OF -YLM- */
/*       MAXMU  :  FIRST DIMENSION OF -YLM- */
/*       TWONM1 :  MAX DEGREE OF -YLM- */
/*       MU(I)  :  I = 1 TO NMU, ARGUMENTS OF -YLM- */
/*       IF M.GT.0, YLM(M-1,I) FOR I = 1 TO NMU IS REQUIRED */
/*  O U T P U T     V A R I A B L E: */
/*       YLM(L,I) :  L = M TO TWONM1, NORMALIZED ASSOCIATED LEGENDRE */
/*                   POLYNOMIALS EVALUATED AT ARGUMENT -MU(I)- */
/* +---------------------------------------------------------------------+ */
    /* Parameter adjustments */
    ylm_dim1 = *maxmu - 0 + 1;
    ylm_offset = 0 + ylm_dim1;
    ylm -= ylm_offset;
    --mu;

    /* Function Body */
    if (pass1) {
	pass1 = FALSE_;
	for (ns = 1; ns <= 1000; ++ns) {
	    sqt[ns - 1] = sqrt((real) ns);
/* L1: */
	}
    }
/*      IF ( 2*TWONM1 .GT. MAXSQT ) */
/*     $   CALL ERRMSG( 'LEPOLY--NEED TO INCREASE PARAM MAXSQT', .TRUE. ) */
    if (*m == 0) {
/*                             ** UPWARD RECURRENCE FOR ORDINARY */
/*                             ** LEGENDRE POLYNOMIALS */
	i__1 = *nmu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ylm[i__ * ylm_dim1] = 1.f;
	    ylm[i__ * ylm_dim1 + 1] = mu[i__];
/* L10: */
	}
	i__1 = *twonm1;
	for (l = 2; l <= i__1; ++l) {
	    i__2 = *nmu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ylm[l + i__ * ylm_dim1] = (((l << 1) - 1) * mu[i__] * ylm[l - 
			1 + i__ * ylm_dim1] - (l - 1) * ylm[l - 2 + i__ * 
			ylm_dim1]) / l;
/* L20: */
	    }
	}
    } else {
	i__2 = *nmu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*                               ** Y-SUB-M-SUPER-M; DERIVED FROM */
/*                               ** D/A EQS. (11,12) */
/* Computing 2nd power */
	    r__1 = mu[i__];
	    ylm[*m + i__ * ylm_dim1] = -sqt[(*m << 1) - 2] / sqt[(*m << 1) - 
		    1] * sqrt(1.f - r__1 * r__1) * ylm[*m - 1 + i__ * 
		    ylm_dim1];
/*                              ** Y-SUB-(M+1)-SUPER-M; DERIVED FROM */
/*                              ** D/A EQS. (13,14) USING EQS. (11,12) */
	    ylm[*m + 1 + i__ * ylm_dim1] = sqt[*m * 2] * mu[i__] * ylm[*m + 
		    i__ * ylm_dim1];
/* L30: */
	}
/*                                   ** UPWARD RECURRENCE; D/A EQ. (10) */
	i__2 = *twonm1;
	for (l = *m + 2; l <= i__2; ++l) {
	    tmp1 = sqt[l - *m - 1] * sqt[l + *m - 1];
	    tmp2 = sqt[l - *m - 2] * sqt[l + *m - 2];
	    i__1 = *nmu;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ylm[l + i__ * ylm_dim1] = (((l << 1) - 1) * mu[i__] * ylm[l - 
			1 + i__ * ylm_dim1] - tmp2 * ylm[l - 2 + i__ * 
			ylm_dim1]) / tmp1;
/* L40: */
	    }
	}
    }
    return 0;
} /* lepoly_ */

