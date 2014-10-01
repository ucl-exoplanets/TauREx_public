/* mod_atm.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;

/* Subroutine */ int mod_atm__(integer *iuatm, integer *iumix, integer *
	new_pt__, integer *new_mix__, integer *istate, char *atmfile, integer 
	*ifrmatm, char *frmatm, integer *iskatm, integer *icp, integer *ict, 
	real *scp, integer *nlev, integer *nlyr, integer *levout, integer *
	nlout, integer *k_out__, real *pout, real *wgtatm, real *wgtgas, 
	integer *ngases, integer *igas, char *mixfile, integer *ifrmmix, char 
	*frmmix, integer *ioffmix, integer *izmix, integer *imix, integer *
	icpmix, integer *icmix, real *scpmix, real *scmix, real *pd_frac__, 
	real *radius, real *sgrav, real *ratm, real *p, real *t, real *alt, 
	real *z__, real *grav, real *rmix, real *dp_dp__, ftnlen atmfile_len, 
	ftnlen frmatm_len, ftnlen mixfile_len, ftnlen frmmix_len)
{
    /* Initialized data */

    static char gascode[6*50] = "H2O   " "CO2   " "O3    " "N2O   " "CO    " 
	    "CH4   " "O2    " "NO    " "SO2   " "NO2   " "NH3   " "HNO3  " 
	    "OH    " "HF    " "HCl   " "HBr   " "HI    " "ClO   " "OCS   " 
	    "H2CO  " "HOCl  " "N2    " "HCN   " "CH3Cl " "H2O2  " "C2H2  " 
	    "C2H6  " "PH3   " "COF2  " "SF6   " "H2S   " "HCOOH " "HO2   " 
	    "O     " "ClONO2" "NO+   " "HOBr  " "C2H4  " "CH3OH " "H2    " 
	    "He    " "Gas   " "Gas   " "Gas   " "Gas   " "Gas   " "Gas   " 
	    "Gas   " "Gas   " "Gas   ";

    /* System generated locals */
    integer i__1, i__2;
    cllist cl__1;

    /* Builtin functions */
    integer f_clos(cllist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void);
    double log(doublereal);

    /* Local variables */
    static integer k, n, ng, iopen, nlmax, np_mx__;
    static real p_pasc__;
    extern /* Subroutine */ int atmstr_(char *, integer *, integer *, integer 
	    *, char *, integer *, integer *, integer *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    integer *, integer *, ftnlen, ftnlen), readmix_(char *, integer *,
	     integer *, integer *, integer *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    , real *, integer *, real *, real *, real *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 6, 0, "(/,1a,/,1a,10(3x,1a,2x))", 0 };
    static cilist io___10 = { 0, 6, 0, "(13(1pe11.4))", 0 };



/* ccccccccccccccccccccccccccccc  mod_atm  ccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine calls the subroutines atmstr to set up the      cc */
/* c    atmospheric structure, and then call the subroutine readmix     cc */
/* c    to read gas mixing ratios at each model level.                  cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c    atmfile - name of input atmospheric structure file              cc */
/* c      iuatm - unit from which atmospheric structure is read.        cc */
/* c     wgtatm -  mean molecular weight of atmosphere (kg/kmole)       cc */
/* c    ifrmatm - format for reading pressure and temperature:          cc */
/* c               1) formatted 2) unformatted 3) list directed.        cc */
/* c     frmatm - format for reading p and t (ifrmatm = 1 only).        cc */
/* c     iskatm - number of lines to skip at top of atmfile.            cc */
/* c        icp - column containting pressures.                         cc */
/* c        ict - column containting temperattures.                     cc */
/* c        scp - scaling factor to convert pressure to Pascals.        cc */
/* c          p - pressure at each model level (Pascals).               cc */
/* c          t - temperature at each model level (K).                  cc */
/* c       grav - gravitational acceleration at each level (m/s^2).     cc */
/* c     levout - output level flag:                                    cc */
/* c              1) top of the atmosphere only, 2) surface only,       cc */
/* c              3) top of atmosphere and surface,                     cc */
/* c              4) one or more pressure levels',                      cc */
/* c     ngases - number of absorbing gases included in model.          cc */
/* c    mixfile - name of file with gas mixing ratios                   cc */
/* c      iumix - unit from which mixing ratios are read.               cc */
/* c    ifrmmix - format for reading gas mixing ratios:                 cc */
/* c              1) formatted 2) unformatted 3) list directed.         cc */
/* c     frmmix - format for readinggas mixing ratios(ifrmmix=1 only).  cc */
/* c    ioffmix - number of lines to skip at top of mixfile.            cc */
/* c      izmix - type of vertical coordinate:                          cc */
/* c              1) pressure,   2) altitude (km)                       cc */
/* c       imix - type of mixing ratios: 1) volume, 2) mass.            cc */
/* c     icpmix - column containting vertical coordinater (p, z).       cc */
/* c      icmix - column containting gas mixing ratio                   cc */
/* c     scpmix - factor needed to scale vertical coordinate (Pa, km)   cc */
/* c      scmix - factor needed to scale gas mixing ratios.             cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c       ratm - mean atmospheric gas constant (j/kg/k).               cc */
/* c       nlyr - number of layers in the model atmosphere (nlev-1).    cc */
/* c       nlev - number of levels in the model atmosphere (<kp).       cc */
/* c      k_out - index of output level.                                cc */
/* c      p_out - output pressure level(s) (bars).                      cc */
/* c      dp_dp - local pressure gradient near each output level.       cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccc  mod_atm  ccccccccccccccccccccccccccccccc */




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


    /* Parameter adjustments */
    --dp_dp__;
    rmix -= 71;
    --grav;
    --z__;
    --alt;
    --t;
    --p;
    --pd_frac__;
    --scmix;
    --scpmix;
    --icmix;
    --icpmix;
    --imix;
    --izmix;
    --ioffmix;
    frmmix -= 40;
    --ifrmmix;
    mixfile -= 132;
    --igas;
    --wgtgas;
    --pout;
    --k_out__;
    --istate;

    /* Function Body */

/* ****  read the atmospheric structure file: */

    *ratm = 8314.f / *wgtatm;

    if (*new_pt__ == 1) {
	iopen = 1;
	np_mx__ = 70;

	atmstr_(atmfile, iuatm, &iopen, ifrmatm, frmatm, iskatm, &nlmax, icp, 
		ict, scp, ratm, radius, sgrav, &p[1], &t[1], &z__[1], &alt[1],
		 &grav[1], &np_mx__, nlev, (ftnlen)132, (ftnlen)40);

	cl__1.cerr = 0;
	cl__1.cunit = *iuatm;
	cl__1.csta = 0;
	f_clos(&cl__1);

    }

/* ****    define number of layers and number of output levels */

    *nlyr = *nlev - 1;

    if (*levout == 1) {
	k_out__[1] = 1;
	pout[1] = 0.f;
	dp_dp__[1] = 0.f;
    } else {
	if (*levout == 2) {
	    k_out__[1] = *nlev;
	    pout[1] = p[*nlev] * 1e-5f;
	    dp_dp__[1] = 0.f;
	} else {
	    if (*levout == 3) {
		k_out__[1] = 1;
		pout[1] = 0.f;
		dp_dp__[1] = 0.f;
		k_out__[2] = *nlev;
		pout[2] = p[*nlev] * 1e-5f;
		dp_dp__[2] = 0.f;
	    } else {
		i__1 = *nlout;
		for (n = 1; n <= i__1; ++n) {
		    p_pasc__ = pout[n] * 1e5f;
		    i__2 = *nlev;
		    for (k = 1; k <= i__2; ++k) {
			k_out__[n] = k;
			if (p[k] > p_pasc__) {
			    goto L1041;
			}
/* L1021: */
		    }
L1041:
/* L1061: */
		    ;
		}
	    }

	}
    }

    if (*levout > 3) {

/* ****     find the pressure derivative between levels ajacent to each */
/*         output pressure level. */

	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {
	    if (k_out__[k] > 1) {
		dp_dp__[k] = (p[k_out__[k]] - pout[k] * 1e5f) / (p[k_out__[k]]
			 - p[k_out__[k] - 1]);
	    } else {
		dp_dp__[k] = (p[k_out__[k]] - pout[k] * 1e5f) / (p[k_out__[k] 
			+ 1] - p[k_out__[k]]);
	    }
/* L1081: */
	}
    }

    if (*new_mix__ == 1) {

/* ****     read gas mixing ratios */

	i__1 = *ngases;
	for (n = 1; n <= i__1; ++n) {
	    ng = n;
	    nlmax = 0;
	    iopen = 1;

	    readmix_(mixfile + ng * 132, iumix, &ng, &iopen, &ifrmmix[ng], 
		    frmmix + ng * 40, &ioffmix[ng], &nlmax, &izmix[ng], &imix[
		    ng], &icpmix[ng], &icmix[ng], &scpmix[ng], &scmix[ng], 
		    wgtatm, &wgtgas[1], nlev, &alt[1], &z__[1], &rmix[71], (
		    ftnlen)132, (ftnlen)40);

	    cl__1.cerr = 0;
	    cl__1.cunit = *iumix;
	    cl__1.csta = 0;
	    f_clos(&cl__1);

/* L1141: */
	}

    }

/* ****   print the atmospheric structure and gas mixing ratios */

    s_wsfe(&io___9);
    do_fio(&c__1, " Atmospheric Structure and Gas Mixing Ratios", (ftnlen)44);
    do_fio(&c__1, "    alt(km)    p (Pa)      t(K)    g (m/s/s) ", (ftnlen)45)
	    ;
    i__1 = *ngases;
    for (n = 1; n <= i__1; ++n) {
	do_fio(&c__1, gascode + (igas[n] - 1) * 6, (ftnlen)6);
    }
    e_wsfe();
    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&alt[k], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&p[k], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&t[k], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&grav[k], (ftnlen)sizeof(real));
	i__2 = *ngases;
	for (n = 1; n <= i__2; ++n) {
	    do_fio(&c__1, (char *)&rmix[k + n * 70], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L2001: */
    }

    if (istate[1] != 0) {

/* ****     pressure is a varible part of the state vector */
/*         find effective altitude for perturbed surface pressure */

	t[*nlev + 1] = t[*nlev];
	grav[*nlev + 1] = grav[*nlev];
	z__[*nlev + 1] = log((pd_frac__[1] + 1.f) * p[*nlev]);
	alt[*nlev + 1] = alt[*nlev - 1] - *ratm * 5e-4f * (t[*nlev + 1] + t[*
		nlev - 1]) * (z__[*nlev + 1] - z__[*nlev - 1]) / grav[*nlev + 
		1];
	i__1 = *ngases;
	for (n = 1; n <= i__1; ++n) {
	    rmix[*nlev + 1 + n * 70] = rmix[*nlev + n * 70];
/* L3001: */
	}

    }

    return 0;
} /* mod_atm__ */

