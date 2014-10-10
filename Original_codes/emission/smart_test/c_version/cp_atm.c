/* cp_atm.f -- translated by f2c (version 20100827).
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

doublereal cp_atm__(integer *ncomp, integer *icomp, real *volmix, real *t)
{
    /* Initialized data */

    static real cp1[2] = { -2.629e-7f,3.47e-7f };
    static real cp2[2] = { 6.059e-4f,-.001269f };
    static real cp3[2] = { -.2232f,1.688f };
    static real cp4[2] = { 1024.f,443.1f };

    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__;
    static real cp, vmix;


/* ccccccccccccccccccccccccccccc  c p _ a t m  ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes the specific heat at constant pressure cc */
/* c    for a specified temperature and atmospheric temperature.        cc */
/* c    specific heats for specific gases are determined by evaluating  cc */
/* c    a cubic polynomial of the form:                                 cc */
/* c                                                                    cc */
/* c                cp(t) = cp1*t**3 + cp2*t**2 + cp3*t + cp4           cc */
/* c                                                                    cc */
/* c    These coefficients are derived from cp/R vs T data (Hilsenrath) cc */
/* c    with the aid of the program: /home/dc/util/fitcp.f              cc */
/* c                                                                    cc */
/* c    The net atmospheric specific heat of a mixture of gases is      cc */
/* c    then determined by taking the volume-mixing-ratio-weighted      cc */
/* c    average of the specific heats of all components.                cc */
/* c                                                                    cc */
/* c    r e f e r e n c e s :                                           cc */
/* c                                                                    cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c    ncomp = number of major atmosphere constituents                 cc */
/* c    icomp = atmosphere constituent index                            cc */
/* c            (1) air  (2) co2  (3) n2  (4) o2  (5) h2  (6) he        cc */
/* c   volmix = volume mixing ratio of each constituent                 cc */
/* c        t = temperature at which cp is needed (K)                   cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c       cp = specific heat at constant pressure at temperature t.    cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccc  c p _ a t m  ccccccccccccccccccccccccccc */






    /* Parameter adjustments */
    --volmix;
    --icomp;

    /* Function Body */

    vmix = 0.f;
    cp = 0.f;

    i__1 = *ncomp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vmix += volmix[i__];
	if (icomp[i__] <= 2) {
	    cp += volmix[i__] * (cp1[icomp[i__] - 1] * *t * *t * *t + cp2[
		    icomp[i__] - 1] * *t * *t + cp3[icomp[i__] - 1] * *t + 
		    cp4[icomp[i__] - 1]);
	}

/* ****          use temperature-independent values */

	if (icomp[i__] == 3) {

/* ****          the following N2 value comes from the CRC */

	    cp += 1059.f;

	}

	if (icomp[i__] == 4) {

/* ****            the following O2 value comes from the CRC pg D165 */

	    cp += 916.29f;

	}

	if (icomp[i__] == 5) {

/* ****        the following value comes from the hilsenrath tables */
/*            for H2 at 80K */

	    cp += 10933.f;

	}

	if (icomp[i__] == 6) {

/* ****               the following He value comes from the CRC pg D165 */

	    cp += 5196.f;

	}
/* L1001: */
    }

    ret_val = cp / vmix;

    return ret_val;
} /* cp_atm__ */

