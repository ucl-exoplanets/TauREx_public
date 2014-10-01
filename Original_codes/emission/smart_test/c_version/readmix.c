/* readmix.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int readmix_(char *mixfile, integer *iumix, integer *ng, 
	integer *iopen, integer *ifrmmix, char *frmmix, integer *ioffmix, 
	integer *nlmax, integer *izmix, integer *imix, integer *icpmix, 
	integer *icmix, real *scpmix, real *scmix, real *wgtatm, real *wgtgas,
	 integer *nlev, real *alt, real *z__, real *rmix, ftnlen mixfile_len, 
	ftnlen frmmix_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double log(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double exp(doublereal);

    /* Local variables */
    extern /* Subroutine */ int xyinterp_(real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer k;
    static real zs[8096], add[20], scz[20];
    static integer icol[20], ncol, ierr, nptj, invz;
    static real zmix[70], rmix0[8096], rmix1[70];
    static integer ncdim, ncmax;
    static real array[20]	/* was [1][20] */;
    static integer iskip, nzdim, nlstd, nxmax, nymax;
    extern /* Subroutine */ int readfil_(integer *, integer *, char *, char *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, real *, real *, real *, integer 
	    *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, "(1a,i5,1a,i5)", 0 };



/* ccccccccccccccccccccccccccc r e a d m i x ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    ********************   p u r p o s e   **********************   cc */
/* c                                                                    cc */
/* c    this subroutine reads reads gas mixing ratios as a function of  cc */
/* c    pressure, converts pressures to pascals, and mixing ratios to   cc */
/* c    mass mixing ratios, if necessary, and interpolates them to a    cc */
/* c    standard pressure grid.                                         cc */
/* c                                                                    cc */
/* c    ********************     i n p u t     *********************    cc */
/* c                                                                    cc */
/* c     iumix - input unit number for atmospheric variables.           cc */
/* c        ng - hitran gas index for this gas.                         cc */
/* c     iopen - open unit flag - 0) don't open unit, 1) open unit.     cc */
/* c   mixfile - name of input file                                     cc */
/* c     izmix - idex of vertical coordinate: 1) pressure, 2) altitude. cc */
/* c      imix - type of mixing ratios: 1) volume, 2) mass (kg/kg)      cc */
/* c    icpmix - column index of vertical coordinate in each record.    cc */
/* c     icmix - column index of mixing ratio in each record.           cc */
/* c     ifrmmix - format of input file:                                cc */
/* c             1) formatted, 2) unformatted, 3) list directed)        cc */
/* c    frmmix - format for formatted data files (ifrmmix=1 only)       cc */
/* c    scpmix - multiplicative factor to convert pressure to pascals   cc */
/* c             or altitudes to km.                                    cc */
/* c     scmix - multiplicative factor to convert mixing ratio to       cc */
/* c             range (0.0 to 1.0)                                     cc */
/* c      ioffmix - number of records to skip at top of file               cc */
/* c      nlev - number of levels in the standard atmosphere p grid.    cc */
/* c     nlmax - maximum number of levels to read from input file.      cc */
/* c     ps(l) - pressure at level l of input rmix atmosphere.          cc */
/* c  rmixs(l) - mixing ratio at level l of input atmosphere.           cc */
/* c    wgtgas - molecular weight of gas (kg/kmole).                    cc */
/* c    wgtatm - molecular weight of background atmosphere (kg/kmole).  cc */
/* c                                                                    cc */
/* c    ********************   o u t p u t     *********************    cc */
/* c                                                                    cc */
/* c      rmix(l) - mass mixing ratio at model level l (kg/kg).         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc r e a d m i x ccccccccccccccccccccccccccccc */




/* ccccccccccccccccccccccccccccc  p a r a m  ccccccccccccccccccccccccccccc */
/* c								      cc */
/* c    p u r p o s e :                                                 cc */
/* c								      cc */
/* c    this include file specifies the dimensions of variables used    cc */
/* c    in the program smt_do.                                          cc */
/* c								      cc */
/* c    q u a n t i t i e s :                                           cc */
/* c								      cc */
/* c      nsp - Maximum number of spectral intervals 		      cc */
/* c       kp - Maximum number of vertical levels.		      cc */
/* c    nmode - Maximum number of aerosol particle modes.               cc */
/* c     ngas - Maximum number of absorbing gases.		      cc */
/* c     ngrp - Maximum number of SMT grouping bins.		      cc */
/* c    mxumu - Maximum number of user output zenith angles in discrete cc */
/* c            ordinate model.					      cc */
/* c    mxphi - Maximum number of output azimuth angles in the discrete cc */
/* c            ordinate model.					      cc */
/* c    mxmom - Maximum number of legendre polynomial moments in the    cc */
/* c            discrete ordinate method.				      cc */
/* c    mxulv - Maximum number of output levels in the discrete 	      cc */
/* c            ordinate model.					      cc */
/* c    mxrad - Maximum number of output radiance levels.		      cc */
/* c     nsol - Maximum number of solar zenith angles.		      cc */
/* c      nex - Maximum number of spectrally-dependent optical          cc */
/* c            (gases, aerosol modes, surface albedos, and solar flux. cc */
/* c     mxpd - maximum number of radiance partial derivatives          cc */
/* c								      cc */
/* ccccccccccccccccccccccccccccc  p a r a m  ccccccccccccccccccccccccccccc */







/* ****   atmospheric structure variables */


/* *****    initialize i/o variables */

    /* Parameter adjustments */
    rmix -= 71;
    --z__;
    --alt;
    --wgtgas;

    /* Function Body */
    ierr = 0;
    ncol = 2;
    ncdim = 20;
    nzdim = 1;
    invz = 2;
    iskip = 0;
    icol[0] = *icpmix;
    icol[1] = *icmix;
    ncmax = icol[1];
    if (icol[0] > icol[1]) {
	ncmax = icol[0];
    }
    scz[0] = *scpmix;
    scz[1] = *scmix;
    add[0] = 0.f;
    add[1] = 0.f;

    if (*iopen == 1) {

/* ****    o p e n   a t m o s p h e r i c   s t r u c t u r e    f i l e */

	readfil_(&c__0, iumix, mixfile, frmmix, ifrmmix, &nzdim, &ncdim, 
		ioffmix, &iskip, &invz, &ncol, icol, &nptj, scz, add, array, &
		ncmax, &ierr, (ftnlen)132, (ftnlen)40);
	if (ierr != 0) {
	    s_wsle(&io___13);
	    do_lio(&c__9, &c__1, "readmix error opening file: ", (ftnlen)28);
	    do_lio(&c__9, &c__1, mixfile, (ftnlen)132);
	    e_wsle();
	    s_stop("", (ftnlen)0);
	}
    }

/* ****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s . */

    if (*ioffmix >= 1) {
	readfil_(&c_n1, iumix, mixfile, frmmix, ifrmmix, &nzdim, &ncdim, 
		ioffmix, &iskip, &invz, &ncol, icol, &nptj, scz, add, array, &
		ncmax, &ierr, (ftnlen)132, (ftnlen)40);

	if (ierr != 0) {
	    s_wsle(&io___14);
	    do_lio(&c__9, &c__1, "End-of-file Encountered while skipping rec"
		    "ords in file ", (ftnlen)55);
	    do_lio(&c__9, &c__1, mixfile, (ftnlen)132);
	    e_wsle();
	    s_stop("", (ftnlen)0);

	}
    }

/* *****   read each p-T record */

    k = 0;

L2401:
    nptj = 1;
    readfil_(&c__1, iumix, mixfile, frmmix, ifrmmix, &nzdim, &ncdim, ioffmix, 
	    &iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &ierr,
	     (ftnlen)132, (ftnlen)40);

    if (ierr != 0) {
	goto L2801;
    }
    ++k;

/* ***      load vectors for vertical coordinate and mixing ratios. */

    if (*izmix == 1) {

/* ****      define the log of pressure */

	if (array[0] > 0.f) {
	    zs[k - 1] = log(array[0]);
	} else {
	    zs[k - 1] = -1e10f;
	}
    } else {

/* ****     define altitude of input */

	zs[k - 1] = array[0];
    }

/* ****     define the log of the mass mixing ratio */

    if (array[1] > 0.f) {
	if (*imix == 1) {
	    rmix0[k - 1] = log(wgtgas[*ng] * array[1] / *wgtatm);
	} else {
	    rmix0[k - 1] = log(array[1]);
	}
    } else {
	rmix0[k - 1] = -100.f;
    }

    if (k < *nlmax || *nlmax == 0) {
	goto L2401;
    }

L2801:
    nlstd = k;

    s_wsfe(&io___19);
    do_fio(&c__1, " Number of input levels for gas ", (ftnlen)32);
    do_fio(&c__1, (char *)&(*ng), (ftnlen)sizeof(integer));
    do_fio(&c__1, " = ", (ftnlen)3);
    do_fio(&c__1, (char *)&nlstd, (ftnlen)sizeof(integer));
    e_wsfe();

/* ****   interpolate the input mixing ratios to the input */
/*       pressure grid.  Assume log of mixing ratio varies linearly */
/*       with the log of pressure. */

    if (*izmix == 1) {
	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    zmix[k - 1] = z__[k];
/* L3201: */
	}
    } else {
	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    zmix[k - 1] = alt[k];
/* L3221: */
	}
    }

    nxmax = 70;
    nymax = 1;

    xyinterp_(zs, rmix0, zmix, rmix1, &nxmax, &nymax, &nlstd, nlev, &c__1);

/* ****     convert output back to mixing ratio */

    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {
	rmix[k + *ng * 70] = exp(rmix1[k - 1]);
/* L4001: */
    }

    return 0;
} /* readmix_ */

