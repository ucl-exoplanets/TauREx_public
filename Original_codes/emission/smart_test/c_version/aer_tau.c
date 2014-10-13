/* aer_tau.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int aer_tau__(integer *ne, integer *nlyr, integer *iuaer, 
	integer *nmodes, integer *nzup, integer *nstr, integer *nmom, 
	doublereal *wn, doublereal *wnmin, doublereal *wnmax, real *umu, real 
	*p, real *dtauex, real *g, real *phmom, real *tau_ext__, real *
	tau_sca__, real *g_sca__, doublereal *wnext, real *aerqext, real *
	aerqsca, real *aerg0, real *aerpmom, real *dqextdv, real *dqscadv, 
	real *dpmomdv, real *dg0dv, real *dtausc, real *dtauaer, real *
	p_aer__, real *p_aer_0__, real *tauaer, integer *iflext, integer *
	nsiext, integer *nmomaer, integer *nmom_mx__, integer *ntau_pd0__, 
	integer *istate, integer *ngases, integer *modepd, doublereal *
	wn_eof__, integer *io_end__, integer *io_err__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static real aer_dtau__[70], dtau_tst__;
    static integer k, l, m;
    static real g0[10], dv;
    static integer iu, mom, nze, nst, mode;
    static real qsca[10], qext[10];
    static integer ktau1, nt_pd__;
    static real pmomaer[2010]	/* was [201][10] */, ptauaer;
    extern /* Subroutine */ int aerosol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);


/* cccccccccccccccccccccccccc  a e r _ t a u  cccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine finds the aerosol extinction and scattering     cc */
/* c    cross sections at each model level, and computes the            cc */
/* c    monochromatic aerosol optical depth for each atmospheric layer. cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       nlyr - number of layers in the model atmosphere.             cc */
/* c     nmodes - number of aerosol particle modes.                     cc */
/* c      iuaer - unit number from which aerosol properties are read.   cc */
/* c       nmom - maximum number of legendre polynomial moments         cc */
/* c      nt_pd - number of temperature profiles needed                 cc */
/* c              1 - radiances only, 2 - partial derivaties            cc */
/* c       nzup - index of first upward radiance stream.                cc */
/* c          p - pressure at each model level (Pascals).               cc */
/* c       nstr - number of computational zenith angles (streams).      cc */
/* c        umu - cosine of each upward radiance stream.                cc */
/* c         wn - current wavenumber (cm**-1).                          cc */
/* c     dtauex - differential extintion optical depth in each layer    cc */
/* c     dtausc - differential scattering optical depth in each layer   cc */
/* c      phmom - scattering phase function moments in each layer.      cc */
/* c    aerqext - monochromatic aerosol extinction efficiency at input  cc */
/* c              wavenumber,  wnext.                                   cc */
/* c    aerqsca - monochromatic aerosol scattering efficiency at input  cc */
/* c              wavenumber,  wnext.                                   cc */
/* c      aerg0 - scattering asymmetry parameter at input wavelength    cc */
/* c    aerpmom - monochromatic aerosol phase function moments          cc */
/* c              at input wavenumber,  wnext.                          cc */
/* c    dqextdv - extinction efficiency gradient at wnext.              cc */
/* c    dqextdv - scattering efficiency gradient at wnext.              cc */
/* c      dg0dv - rate of change of g with wavenumber                   cc */
/* c    dpmomdv - phase function moment gradient at wnext.              cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c     io_end: end of file flag for each wavelength constituent:      cc */
/* c             (if io_end ne 0, the end of file has been encountered) cc */
/* c     io_err: i/o error flag  for each wavelength constituent:       cc */
/* c             (if io_err ne 0, an i/o error has been encountered)    cc */
/* c     wn_eof: last wavenumber in each input file.                    cc */
/* c       qext - monochromatic aerosol extinction efficiency at        cc */
/* c              wavenumber, wn.                                       cc */
/* c       qsca - monochromatic aerosol scattering efficiency at        cc */
/* c              wavenumber, wn.                                       cc */
/* c    tau_ext - column integrated extinction optical depth            cc */
/* c              for each variable compoent of the state vector        cc */
/* c    tau_sca - column integrated scattering optical depth            cc */
/* c              for each variable component of the state vector       cc */
/* c    pmomaer - monochromatic phase function moments at wn.           cc */
/* c     wn_eof: last wavenumber in each input file.                    cc */
/* c     tauaer - column-integrated aerleigh optical depth.             cc */
/* c     dtauex - monochromatic extinction optical depth in each layer. cc */
/* c     dtausc - monochromatic scattering optical depth in each layer. cc */
/* c      phmom - scattering phase function moments in each layer.      cc */
/* c      p_aer - pressure level of aerliegh tau=1 for each upward      cc */
/* c              stream (bars).                                        cc */
/* c    p_aer_0 - pressure of Rayleigh tau=1 for a vertical stream.     cc */
/* c       nmom - maximum number of legendre polynomial moments         cc */
/* c      nt_pd - number of temperature profiles needed                 cc */
/* c              1 - radiances only, 2 - partial derivaties            cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  a e r _ t a u  cccccccccccccccccccccccccccc */




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





/* *****   monochromatic optical properties */


/* ****   input aerosol wavelength-dependent variables */


/* ****    enter loop over particle modes */

/* *****   determine whether temperature is a variable part of */
/*        the state vector.  if it is, extinction optical depth must be */
/*        computed for the background and perterbed temperature profile */

    /* Parameter adjustments */
    --io_err__;
    --io_end__;
    wn_eof__ -= 3;
    --modepd;
    --istate;
    --nmom_mx__;
    nmomaer -= 3;
    nsiext -= 3;
    --iflext;
    --tauaer;
    --p_aer__;
    dtauaer -= 11;
    --dtausc;
    --dg0dv;
    dpmomdv -= 201;
    --dqscadv;
    --dqextdv;
    aerpmom -= 403;
    aerg0 -= 3;
    aerqsca -= 3;
    aerqext -= 3;
    wnext -= 3;
    g_sca__ -= 71;
    tau_sca__ -= 71;
    tau_ext__ -= 281;
    phmom -= 201;
    --g;
    dtauex -= 71;
    --p;
    --umu;

    /* Function Body */
    if (istate[2] == 0) {
	nt_pd__ = 1;
    } else {
	nt_pd__ = 2;
    }

    *p_aer_0__ = p[*nlyr + 1];
    i__1 = *nstr;
    for (nze = 1; nze <= i__1; ++nze) {
	p_aer__[nze] = p[*nlyr + 1];
/* L1001: */
    }
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	aer_dtau__[k - 1] = 0.f;
/* L1021: */
    }

    i__1 = *nmodes;
    for (m = 1; m <= i__1; ++m) {

/* ****       update extinction counter */

	++(*ne);

/* *****         find out whether this aerosol is a variable element of */
/*               the state vector */

	if (iflext[*ne] == 1) {
	    mode = m;
	    iu = *iuaer + m;

/* ****        read aerosol optical properties at next wavenumber. */

	    aerosol_(ne, &mode, &iu, wn, wnmin, wnmax, &aerqext[3], &aerqsca[
		    3], &aerg0[3], &aerpmom[403], &dpmomdv[201], &dqextdv[1], 
		    &dqscadv[1], &dg0dv[1], &nmomaer[3], &nmom_mx__[1], &
		    nsiext[3], &wnext[3], &wn_eof__[3], &io_end__[1], &
		    io_err__[1]);

/* ****         turn off the wavenumber grid flag. */

	    iflext[*ne] = 0;

	}

/* ****         update the legendre polynomial mode counter */

	if (nmom_mx__[m] > *nmom) {
	    *nmom = nmom_mx__[m];
	}


/* ****       find the absorption and scattering efficiencies */
/*           and asymmetry factors at wn for this mode */

	dv = (real) (*wn - wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)]);

	qext[m - 1] = aerqext[nsiext[(*ne << 1) + 1] + (m << 1)] + dqextdv[m] 
		* dv;
	qsca[m - 1] = aerqsca[nsiext[(*ne << 1) + 1] + (m << 1)] + dqscadv[m] 
		* dv;
	g0[m - 1] = aerg0[nsiext[(*ne << 1) + 1] + (m << 1)] + dg0dv[m] * dv;
	i__2 = nmom_mx__[m];
	for (mom = 0; mom <= i__2; ++mom) {
	    pmomaer[mom + m * 201 - 201] = aerpmom[nsiext[(*ne << 1) + 1] + (
		    mom + m * 201 << 1)] + dpmomdv[mom + m * 201] * dv;
/* L2001: */
	}

/* ****       find aerosol differential extinction and scattering */
/*           optical depths for each mode. */

	i__2 = *nlyr;
	for (k = 1; k <= i__2; ++k) {
	    aer_dtau__[k - 1] = 0.f;
	    if (dtauaer[m + k * 10] > 0.f) {
		aer_dtau__[k - 1] += qext[m - 1] * dtauaer[m + k * 10];
		dtausc[k] += qsca[m - 1] * dtauaer[m + k * 10];
		g[k] += qsca[m - 1] * dtauaer[m + k * 10] * g0[m - 1];
		i__3 = nt_pd__;
		for (l = 1; l <= i__3; ++l) {
		    dtauex[k + l * 70] += qext[m - 1] * dtauaer[m + k * 10];
/* L2021: */
		}
		i__3 = nmom_mx__[m];
		for (mom = 0; mom <= i__3; ++mom) {
		    phmom[mom + k * 201] += qsca[m - 1] * dtauaer[m + k * 10] 
			    * pmomaer[mom + m * 201 - 201];
/* L2041: */
		}
	    }
/* L2061: */
	}

/* ****      if pressure is a varible part of the state vector */
/*          a perturbed layer is added to the bottom of the atmosphere */
/*          with a surface pressure greater than the nominal value. */
/*          the aerosol optical properties in this layer are identical to */
/*          those in the lowest layer of the atmosphere */

	if (istate[1] != 0) {

	    dtausc[*nlyr + 1] += qsca[m - 1] * dtauaer[m + *nlyr * 10];
	    i__2 = nt_pd__;
	    for (l = 1; l <= i__2; ++l) {
		dtauex[*nlyr + 1 + l * 70] += qext[m - 1] * dtauaer[m + *nlyr 
			* 10];
/* L4001: */
	    }

	}

/* *****         determine whether this gas is a variable element of */
/*              the state vector.  If it is, update extinction counter */
/*              used in calculating optical depth partial derivatives */

	nst = m + *ngases + 2;
	if ((i__2 = istate[nst], abs(i__2)) == 4) {
	    ++(*ntau_pd0__);
	    modepd[*ntau_pd0__] = m;

	    if (*ntau_pd0__ > 0) {
		i__2 = *nlyr;
		for (k = 1; k <= i__2; ++k) {
		    if (dtauaer[m + k * 10] > 0.f) {
			tau_ext__[k + (*ntau_pd0__ * 3 + 1) * 70] = qext[m - 
				1] * dtauaer[m + k * 10];
			tau_sca__[k + *ntau_pd0__ * 70] = qsca[m - 1] * 
				dtauaer[m + k * 10];
			g_sca__[k + *ntau_pd0__ * 70] = g0[m - 1];
		    }
/* L2261: */
		}
	    }
	}
/* L2281: */
    }

/* ****   find the column-integrated aerosol optical depth */

    tauaer[1] = 0.f;
    i__1 = *nlyr + 1;
    for (k = 2; k <= i__1; ++k) {
	tauaer[k] = tauaer[k - 1] + aer_dtau__[k - 2];
/* L3001: */
    }

/* ****   find the layer of aerosol optical absorption optical depth */
/*        unity for each output stream */

    i__1 = *nstr;
    for (nze = *nzup; nze <= i__1; ++nze) {
	ptauaer = 0.f;
	ktau1 = *nlyr + 1;

	i__2 = *nlyr;
	for (k = 1; k <= i__2; ++k) {
	    if (ptauaer == 0.f && tauaer[k + 1] / umu[nze] >= 1.f) {
		ktau1 = k + 1;
		ptauaer = p[k];
	    }
/* L5001: */
	}

/* ****       find the pressure of aerosol extinction optical */
/*           depth unity at each output emission angle */

	dtau_tst__ = (tauaer[ktau1] - tauaer[ktau1 - 1]) / umu[nze];
	if (dtau_tst__ != 0.f) {
	    p_aer__[nze] = (p[ktau1 - 1] + (1.f - tauaer[ktau1 - 1] / umu[nze]
		    ) * (p[ktau1] - p[ktau1 - 1]) / dtau_tst__) * 1e-5f;
	} else {
	    p_aer__[nze] = p[ktau1] * 1e-5f;
	}
/* L5021: */
    }

/* ****    find the layer of gas optical absorption optical depth */
/*        unity for normal incident radiation */

    ptauaer = 0.f;
    ktau1 = *nlyr + 1;

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	if (ptauaer == 0.f && tauaer[k + 1] >= 1.f) {
	    ktau1 = k + 1;
	    ptauaer = p[k];
	}
/* L5041: */
    }

/* ****  find the pressure of aerosol extinction optical */
/*      depth unity at each output emission angle */

    dtau_tst__ = tauaer[ktau1] - tauaer[ktau1 - 1];
    if (dtau_tst__ != 0.f) {
	*p_aer_0__ = (p[ktau1 - 1] + (1.f - tauaer[ktau1 - 1]) * (p[ktau1] - 
		p[ktau1 - 1]) / dtau_tst__) * 1e-5f;
    } else {
	*p_aer_0__ = p[ktau1] * 1e-5f;
    }

    return 0;
} /* aer_tau__ */

