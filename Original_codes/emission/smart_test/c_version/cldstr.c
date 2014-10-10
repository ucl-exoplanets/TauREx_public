/* cldstr.f -- translated by f2c (version 20100827).
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
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__9 = 9;

/* Subroutine */ int cldstr_(integer *nlev, integer *iopen, integer *mode, 
	integer *iuaer, char *aerfile, integer *iofftau, integer *iztau, 
	integer *icptau, integer *ictau, real *scptau, real *sctau, logical *
	lcsh, integer *iccsh, real *sccsh, real *p, real *alt, real *z__, 
	real *dtauaer, ftnlen aerfile_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double log(doublereal), exp(doublereal);

    /* Local variables */
    extern /* Subroutine */ int xyinterp_(real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer i__, k, l, m, n;
    static real z1;
    static integer ic;
    static real dz, add[80];
    static integer icl;
    static real csh[70], dz10, scz[80], altc[70], pbar;
    static integer icol[80];
    static real cshm;
    static integer ncol, ndiv, ierr;
    static real zaer[70];
    static integer nptj, invz;
    static real altc0[7000], tauc0[7000];
    static integer ncdim, ncmax, nclev;
    static real fndiv;
    static integer iform;
    static real array[80]	/* was [1][80] */;
    static integer iskip, nzdim;
    static real dtaus[70], cldab0;
    static integer nxmax, nymax;
    static real dtaus1, cldabs[7000], altaer[70], tauaer[70];
    static integer ncloud;
    static real tautot;
    static char formxy[40];
    extern /* Subroutine */ int readfil_(integer *, integer *, char *, char *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, real *, real *, real *, integer 
	    *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, "(1a,i5,1a,i5)", 0 };
    static cilist io___46 = { 0, 6, 0, "(/,1a,i5,/,/,1a)", 0 };
    static cilist io___49 = { 0, 6, 0, "(4(1pe12.4))", 0 };



/* ccccccccccccccccccccccccccc  c l d s t r  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    ********************   p u r p o s e   *********************    cc */
/* c                                                                    cc */
/* c    this subroutine initializes the aerosol properties of the       cc */
/* c    atmosphere for radiative transfer models.                       cc */
/* c                                                                    cc */
/* c    ********************     i n p u t     *********************    cc */
/* c                                                                    cc */
/* c        iztau - vertical coordinate type 1) pressure, 2) altitude   cc */
/* c      iofftau - number of records to skip above optical depths      cc */
/* c       icptau - column with vertical coordinate                     cc */
/* c        ictau - column with differential optical depth              cc */
/* c       scptau - scaling factor for optical depth vertical           cc */
/* c                coordinate (pressures -> pascals, altitudes -> km)  cc */
/* c        sctau - scaling factor for optical depth                    cc */
/* c         lsch - scale height flag (true -> use scale heights)       cc */
/* c        icsch - column with scale heights.                          cc */
/* c        sccsh - scaling factor to convert scale heights to km       cc */
/* c        nclev - number of discrete cloud layers.                    cc */
/* c       cpl(l) - pressure (bars) at homogeneous aerosol layer boun-  cc */
/* c                daries                                              cc */
/* c     csh(m,l) - cloud particle scale height (km) for particle mode  cc */
/* c                m in homogeneous layer l                            cc */
/* c   dtaus(m,l) - aerosol differential optical depth for mode m in    cc */
/* c                layer l                                             cc */
/* c       ncloud - number of homogeneous aerosol layers                cc */
/* c       nmodes - number of aerosol particle modes                    cc */
/* c         tcld - total fractional cloudiness by all layers           cc */
/* c                                                                    cc */
/* c    ********************    o u t p u t    *********************    cc */
/* c                                                                    cc */
/* c                            - unit 6 -                              cc */
/* c                                                                    cc */
/* c   dtaus(m,l) - aerosol differential optical depth for mode m in    cc */
/* c                layer l                                             cc */
/* c                                                                    cc */
/* c                       - calling arguments -                        cc */
/* c                                                                    cc */
/* c                                                                    cc */
/* c                    - other computed quantities -                   cc */
/* c                                                                    cc */
/* c cpnd(m,l,na) - number density between model half-levels l and l+1  cc */
/* c                for mode m and model atmosphere na (particles/cm**3)cc */
/* c         nstd - spectral interval index where standard atmosphere   cc */
/* c                aerosol optical depths are specified                cc */
/* c dtauaer(m,l) - fractional optical depth for particle mode m, at    cc */
/* c                lth level in the nath model atmosphere              cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  c l d s t r cccccccccccccccccccccccccccccc */




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




/* *****   file structure variables */




/* ****   readfil variables */



    /* Parameter adjustments */
    dtauaer -= 11;
    --z__;
    --alt;
    --p;

    /* Function Body */
    m = *mode;

/* *****    initialize i/o variables */

    ierr = 0;
    nzdim = 1;
    ncdim = 80;
    s_copy(formxy, " ", (ftnlen)40, (ftnlen)1);
    iform = 3;
    invz = 2;
    iskip = 0;
    if (*lcsh) {
	ncol = 3;
	icol[2] = *iccsh;
	scz[2] = *sccsh;
    } else {
	ncol = 2;
	for (k = 1; k <= 70; ++k) {
	    csh[k - 1] = 0.f;
/* L1001: */
	}
    }
    icol[0] = *icptau;
    icol[1] = *ictau;
    scz[0] = *scptau;
    scz[1] = *sctau;
    add[0] = 0.f;
    add[1] = 0.f;
    ncmax = icol[0];
    i__1 = ncol;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (icol[i__ - 1] > ncmax) {
	    ncmax = icol[i__ - 1];
	}
/* L1201: */
    }

    if (*iopen == 1) {

/* ****    o p e n   a e r o s o l   s t r u c t u r e    f i l e */

	readfil_(&c__0, iuaer, aerfile, formxy, &iform, &nzdim, &ncdim, 
		iofftau, &iskip, &invz, &ncol, icol, &nptj, scz, add, array, &
		ncmax, &ierr, (ftnlen)132, (ftnlen)40);
	if (ierr != 0) {
	    s_wsfe(&io___19);
	    do_fio(&c__1, " cldstr error opening file: ", (ftnlen)28);
	    do_fio(&c__1, aerfile, (ftnlen)132);
	    e_wsfe();
	    s_stop("", (ftnlen)0);
	}
    }

    if (*iofftau >= 1) {

/* ****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s . */

	readfil_(&c_n1, iuaer, aerfile, formxy, &iform, &nzdim, &ncdim, 
		iofftau, &iskip, &invz, &ncol, icol, &nptj, scz, add, array, &
		ncmax, &ierr, (ftnlen)132, (ftnlen)40);

	if (ierr != 0) {
	    s_wsle(&io___20);
	    do_lio(&c__9, &c__1, "End-of-file Encountered while skipping rec"
		    "ords in file ", (ftnlen)55);
	    do_lio(&c__9, &c__1, aerfile, (ftnlen)132);
	    e_wsle();
	    s_stop("", (ftnlen)0);

	}
    }

/* *****   read each p-tau record */

    k = 0;
    dtaus[0] = 0.f;

L2402:
    nptj = 1;
    readfil_(&c__1, iuaer, aerfile, formxy, &iform, &nzdim, &ncdim, iofftau, &
	    iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &ierr, 
	    (ftnlen)132, (ftnlen)40);

    if (ierr != 0) {
	goto L2801;
    }
    ++k;

/* ***      load vectors for vertical coordinate and aerosol amounts. */

    if (*iztau == 1) {

/* ****      the vertical coordinate is pressure */
/*          define the log of the cloud-level pressure */

	if (array[0] > 0.f) {
	    zaer[k - 1] = log(array[0]);
	} else {
	    zaer[k - 1] = 1e10f;
	}

    } else {

/* ****     define altitude of input */

	zaer[k - 1] = array[0];

    }

/* ****     define the differential optical depth above level k */

    if (k > 1) {
	dtaus[k - 2] = array[1];
    }

    if (*lcsh) {
	csh[k - 1] = array[2];
    }

    if (k < nclev || nclev == 0) {
	goto L2402;
    }

L2801:
    nclev = k;
    ncloud = nclev - 1;

    s_wsfe(&io___25);
    do_fio(&c__1, " Number of input levels for aerosol mode ", (ftnlen)41);
    do_fio(&c__1, (char *)&(*mode), (ftnlen)sizeof(integer));
    do_fio(&c__1, " = ", (ftnlen)3);
    do_fio(&c__1, (char *)&nclev, (ftnlen)sizeof(integer));
    e_wsfe();

/* ****   find altitudes corresponding to the aerosol layer boundaries */

    if (*iztau == 1) {

/* ****     the input vertical coordinate, zaer, is log pressure. */
/*         interpolate to the altitude grid */

	nxmax = 70;
	nymax = 1;

	xyinterp_(&z__[1], &alt[1], zaer, altc, &nxmax, &nymax, nlev, &nclev, 
		&c__1);

    } else {

	if (zaer[nclev - 1] < zaer[0]) {
	    i__1 = nclev;
	    for (k = 1; k <= i__1; ++k) {
		altc[k - 1] = zaer[k - 1];
/* L3001: */
	    }
	} else {

/* ****       the array must be turned around */

	    i__1 = nclev / 2;
	    for (k = 1; k <= i__1; ++k) {
		z1 = zaer[k - 1];
		dtaus1 = dtaus[k - 1];
		zaer[k - 1] = zaer[nclev - k];
		dtaus[k - 1] = dtaus[nclev - k];
		zaer[nclev - k] = z1;
		dtaus[nclev - k] = dtaus1;
/* L3021: */
	    }

	    i__1 = nclev;
	    for (k = 1; k <= i__1; ++k) {
		altc[k - 1] = zaer[k - 1];
/* L3041: */
	    }
	}
    }

/* ****   interpolate grid to 10 times the vertical resolution.  use the */
/*       cloud particle scale height data to find the optical depth */
/*       profile in each sub-layer. */

    ndiv = 7000 / ncloud - 1;
    fndiv = (real) ndiv;
    icl = 0;
    i__1 = ncloud;
    for (l = 1; l <= i__1; ++l) {
	i__2 = ndiv;
	for (ic = 1; ic <= i__2; ++ic) {
	    ++icl;
	    tauc0[icl - 1] = 0.f;
	    cldabs[icl - 1] = 0.f;
/* L3201: */
	}
/* L3221: */
    }

    icl = 1;
    if (alt[1] > altc[0]) {
	altc0[icl - 1] = alt[1];
	tauc0[icl - 1] = 0.f;
	cldabs[icl - 1] = 0.f;
	++icl;
    }

    altc0[icl - 1] = altc[0];
    tauc0[icl - 1] = 0.f;
    cldabs[icl - 1] = 0.f;

    i__1 = ncloud;
    for (n = 1; n <= i__1; ++n) {

/* ****        determine the thickness of the cloud layer */

	dz = altc[n - 1] - altc[n];

/* ****       determine if the layer has a constant or */
/*           variable number density, and find the absorption */
/*           coefficient per unit length at the base of the layer. */

	dz10 = dz / fndiv;
	cshm = csh[n];
	if (dabs(cshm) > 1e-5f) {

/* ****          use the cloud particle scale heights to distribute */
/*              cloud particles */

	    cldab0 = dtaus[n - 1] / (cshm * (1.f - exp(-dz / cshm)));
	}

	i__2 = ndiv;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++icl;
	    altc0[icl - 1] = altc[n - 1] - (real) i__ * dz10;
	    if (dabs(cshm) > 1e-5f) {

/* ****             use the particle scale height to distribute */
/*                 the particles in the layer */

		cldabs[icl - 1] = cldab0 * exp(-(altc0[icl - 1] - altc[n]) / 
			cshm);
	    } else {

/* ****             the scale height is specified as zero -assume */
/*                 particles have a constant number density. */

		cldabs[icl - 1] = dtaus[n - 1] / dz;

	    }
	    tauc0[icl - 1] = tauc0[icl - 2] + (cldabs[icl - 1] + cldabs[icl - 
		    2]) * .5f * (altc0[icl - 2] - altc0[icl - 1]);
/* L4001: */
	}
/* L4021: */
    }

/* ****   interpolate the cloud absorption coefficients */
/*       to the pressure grid. */

    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {
	altaer[k - 1] = alt[k];
/* L4041: */
    }

    nxmax = 7000;
    nymax = 1;

    xyinterp_(altc0, tauc0, altaer, tauaer, &nxmax, &nymax, &icl, nlev, &
	    nymax);

/* ****   find the differential optical depth in each layer. */

    i__1 = *nlev - 1;
    for (k = 1; k <= i__1; ++k) {
	dtauaer[m + k * 10] = tauaer[k] - tauaer[k - 1];
/* L4061: */
    }

/* ****   print the cumulative optical depth for each particle mode */

    s_wsfe(&io___46);
    do_fio(&c__1, " Optical depth for particle mode:", (ftnlen)33);
    do_fio(&c__1, (char *)&(*mode), (ftnlen)sizeof(integer));
    do_fio(&c__1, "    alt(km)     p (bar)      dtau        tau", (ftnlen)44);
    e_wsfe();
    tautot = 0.f;
    i__1 = *nlev - 1;
    for (k = 1; k <= i__1; ++k) {
	pbar = p[k + 1] * 1e-5f;
	tautot += dtauaer[m + k * 10];
	s_wsfe(&io___49);
	do_fio(&c__1, (char *)&alt[k + 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&pbar, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dtauaer[m + k * 10], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&tautot, (ftnlen)sizeof(real));
	e_wsfe();
/* L4281: */
    }

    return 0;
} /* cldstr_ */

