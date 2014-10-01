/* map_spect.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int map_spect__(integer *nlyr, integer *nlay, integer *nstr, 
	integer *nmom, integer *levout, integer *levtau1, integer *nref, 
	integer *nza, doublereal *wn, doublereal *delnu, real *dtauex, real *
	co_pi0__, real *g, real *phmom, real *surf_opt__, real *alb, real *
	taumn, real *tauerr, real *pi0err, real *phferr, real *surferr, 
	doublereal *dnugrp, doublereal *wngrp, real *surfgrp, real *albgrp, 
	real *taugrp, real *pi0grp, real *ggrp, real *pmomgrp, integer *
	nmomgrp, integer *iwngrp, integer *ngroup, integer *ismterr, integer *
	itau1cnt, integer *itaucnt, integer *ipi0cnt, integer *igcnt, integer 
	*isurcnt)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static real smallpi0, smallalb, tauerror;
    static integer i__, j, k, m, nr, nz;
    static real gmn;
    static integer mom;
    static real gmx;
    static doublereal dng0, dng1;
    static integer nbin;
    static real dtmn, dtmx, gmin0, gmax0, pi0mn, pi0mx, g_rng__[35840]	/* 
	    was [70][512] */, ph_mn__, ph_mx__, pi0min0, ph_rng__[7168000]	
	    /* was [200][70][512] */, abserr, smallg, tauabs, pi0max0, dtaumn[
	    35840]	/* was [70][512] */, tautot[140]	/* was [70][2]
	     */, albmin0, errtst, albmax0, ph_min0__, ph_max0__, pi0_rng__[
	    35840]	/* was [70][512] */, taumin0, taumax0, alb_rng__[2048]
	    	/* was [512][4] */, tau_rng__[35840]	/* was [70][512] */;


/* ccccccccccccccccccccccccc  m a p _ s p e c t  ccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine takes monochromatic optical properties and      cc */
/* c    uses a spectral binning technique to identify monochromatic     cc */
/* c    spectral regions with similar absorption properties at all      cc */
/* c    levels of the atmosphere.  These intervals are combined into a  cc */
/* c    smaller number of quasi-monochrmatic bins.                      cc */
/* c                                                                    cc */
/* c    ***  this version of the program uses an exponetial-sum adding  cc */
/* c         method to combine simlilar optical depths.                 cc */
/* c         it also uses an error criterion that depends on the        cc */
/* c         layer absorption optical depth - layers with optical       cc */
/* c         depths less than unity have the most stringent constraint. cc */
/* c         Other values are weighted as:                              cc */
/* c            tauerr =  err*(1 + 0.1dtau + 0.001*dtau**2)             cc */
/* c								      cc */
/* c    note:  this version of the program uses surface optical         cc */
/* c           property gradients.                                      cc */
/* c    note:  this version of the program uses optcial depth           cc */
/* c           gradients and finds corresponding values of pi0grp.      cc */
/* c                                                                    cc */
/* c    i n p u t    p a r a m e t e r s :                              cc */
/* c        wn : wavenumber of input monochromatic interval.            cc */
/* c      nlyr : number of layers in which coefficients are found       cc */
/* c      nlay : number of layers were optical depths are needed:       cc */
/* c              nlyr if pressure is not a varible part of the state   cc */
/* c              vector, nlev if it is (istate(1) = 1)                 cc */
/* c      nmom : number of legendre polynomial coefficients for the     cc */
/* c              scattering phase function expansion.                  cc */
/* c       nref - number of surface optical properties specified at     cc */
/* c              each wavelength.                                      cc */
/* c    tauerr : fractional dtau error allowed in interval bins         cc */
/* c    pi0err : fractional pi0 error allowed in interval bins          cc */
/* c    phferr : maximum fractional asymmetry parameter error allowed   cc */
/* c   surferr : fractional surface optical property error allowed      cc */
/* c             in interval bins                                       cc */
/* c     delnu : width of monochromatic spectral interval               cc */
/* c      dtau : extinction optical depth in spectral interval          cc */
/* c       pi0 : single scattering albedo in spectral interval          cc */
/* c                                                                    cc */
/* c    o u t p u t    v a r i a b l e s :                              cc */
/* c                                                                    cc */
/* c     taugrp - array of total optical depth bins at model            cc */
/* c            levels for each temperature profile                     cc */
/* c    ismterr - error flag: ngroup > ngrp                             cc */
/* c     ngroup - number of spectral bins (0 to ngrp)                   cc */
/* c     iwngrp - index of spectral bin                                 cc */
/* c   itau1cnt - number of bins rejected by poor fit at tau=1          cc */
/* c    itaucnt - number of bins rejected by poor fit away from tau=1   cc */
/* c    ipi0cnt - number of bins rejected by poor pi0 fit               cc */
/* c      igcnt - number of bins rejected by poor g fit                 cc */
/* c    isurcnt - number of bins rejected by poor albedo fit            cc */
/* c    levtau1 - index of tau=1 level                                  cc */
/* c    nmomgrp - maximum number of phase function moments in bin       cc */
/* c     taugrp - mean, min, and max optical depth in each bin          cc */
/* c     piogrp - mean, min, and max single scattering albedo in  bin   cc */
/* c    pmomgrp - mean, min, and max phase function in each bin         cc */
/* c       ggrp - mean, min, and max asymmetry parmeter in each bin     cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  m a p _ s p e c t  ccccccccccccccccccccccccc */




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



/* ****   index of optical depth unity for each bin */


/* ****   spectral binning parameters */


/* ****    internal variables */




/*       monochromatic optical properties */


/* ****   internal variables */


/* ****   variables defining bin limits */



/* ****    define a small number as a minimum tau and pi0 error */

    /* Parameter adjustments */
    --nmomgrp;
    pmomgrp -= 7218111;
    ggrp -= 35911;
    pi0grp -= 35911;
    taugrp -= 35911;
    albgrp -= 5633;
    surfgrp -= 2561;
    wngrp -= 513;
    --dnugrp;
    --alb;
    --surf_opt__;
    phmom -= 201;
    --g;
    co_pi0__ -= 71;
    dtauex -= 71;
    levtau1 -= 513;

    /* Function Body */
    smallpi0 = *pi0err * .005f * *pi0err + 1e-4f;
    smallg = *phferr * .005f * *phferr + 1e-4f;
    smallalb = *surferr * .005f * *surferr + 1e-4f;

/* ****   define the reference level for integrated absorption optical */
/*       depths.  This reference level will be the first level tested */
/*       by the binning scheme.  If radiances are needed only at the top */
/*       of the atmosphere, the optical depth reference is the top */
/*       of the atmosphere. If radiances are needed at the surface, */
/*       the reference level is the surface. */

    j = 1;
    if (*levout == 2) {
	j = 2;
    }

/* *****  compare this monochromatic segment to each smt bin */

    nbin = 0;
    errtst = 1e9f;
    for (m = *ngroup; m >= 1; --m) {
	abserr = 0.f;

/* ****         check errors at tau = 1 */

	dtmn = taugrp[levtau1[m + (j << 9)] + (m + 1536) * 70] - tau_rng__[
		levtau1[m + (j << 9)] + m * 70 - 71];
	if (dtmn < dtaumn[levtau1[m + (j << 9)] + m * 70 - 71]) {
	    dtmn = dtaumn[levtau1[m + (j << 9)] + m * 70 - 71];
	}
	if (dtauex[levtau1[m + (j << 9)] + 70] < dtmn) {
	    ++(*itau1cnt);
	    goto L1712;
	} else {
	    dtmx = taugrp[levtau1[m + (j << 9)] + (m + 1024) * 70] + 
		    tau_rng__[levtau1[m + (j << 9)] + m * 70 - 71];
	    if (dtauex[levtau1[m + (j << 9)] + 70] > dtmx) {
		++(*itau1cnt);
		goto L1712;
	    }
	}
	abserr += (r__1 = dtauex[levtau1[m + (j << 9)] + 70] - taugrp[levtau1[
		m + (j << 9)] + (m + 512) * 70], dabs(r__1)) / tau_rng__[
		levtau1[m + (j << 9)] + m * 70 - 71];

/* ****       check surface albedo match */

	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    if (alb[nz] < albgrp[m + (nz + 30 << 9)] - alb_rng__[m + (nz << 9)
		     - 513]) {
		++(*isurcnt);
		goto L1712;
	    } else {
		if (alb[nz] > albgrp[m + (nz + 20 << 9)] + alb_rng__[m + (nz 
			<< 9) - 513]) {
		    ++(*isurcnt);
		    goto L1712;
		}
	    }
/* L1001: */
	}

/* ****       check optical depth values in all layers below tau = 1 */

	i__1 = *nlyr;
	for (k = levtau1[m + (j << 9)] + 1; k <= i__1; ++k) {

/* ****           define the optical depth limits of this bin */

	    dtmn = taugrp[k + (m + 1536) * 70] - tau_rng__[k + m * 70 - 71];
	    if (dtmn < dtaumn[k + m * 70 - 71]) {
		dtmn = dtaumn[k + m * 70 - 71];
	    }
	    dtmx = taugrp[k + (m + 1024) * 70] + tau_rng__[k + m * 70 - 71];
	    if (dtauex[k + 70] < dtmn || dtauex[k + 70] > dtmx) {

/* ****               this value does not belong in this bin. */

		++(*itaucnt);
		goto L1712;
	    }
	    abserr += (r__1 = dtauex[k + 70] - taugrp[k + (m + 512) * 70], 
		    dabs(r__1)) / tau_rng__[k + m * 70 - 71];
/* L1201: */
	}

/* ****       check optical depth values in all layers above tau = 1 */

	for (k = levtau1[m + (j << 9)] - 1; k >= 1; --k) {

/* ****           define the optical depth limits of this bin */

	    dtmn = taugrp[k + (m + 1536) * 70] - tau_rng__[k + m * 70 - 71];
	    if (dtmn < dtaumn[k + m * 70 - 71]) {
		dtmn = dtaumn[k + m * 70 - 71];
	    }
	    dtmx = taugrp[k + (m + 1024) * 70] + tau_rng__[k + m * 70 - 71];
	    if (dtauex[k + 70] < dtmn || dtauex[k + 70] > dtmx) {

/* ****               this value does not belong in this bin. */

		++(*itaucnt);
		goto L1712;
	    }
	    abserr += (r__1 = dtauex[k + 70] - taugrp[k + (m + 512) * 70], 
		    dabs(r__1)) / tau_rng__[k + m * 70 - 71];
/* L1221: */
	}

/* ****        check single scattering co-albedos in layers below tau = 1 */

	i__1 = *nlyr;
	for (k = levtau1[m + (j << 9)]; k <= i__1; ++k) {

/* ****            if the single scattering co-albedo is larger than */
/*                a minimum value -scattering is important */

	    if ((1.f - co_pi0__[k + 70]) * dtauex[k + 70] > *taumn) {

/* ****             the scattering optical depth is finite. */
/*                 define the single-scattering co-albedo limits */

		pi0mn = pi0grp[k + (m + 1536) * 70] - pi0_rng__[k + m * 70 - 
			71];
		pi0mx = pi0grp[k + (m + 1024) * 70] + pi0_rng__[k + m * 70 - 
			71];
		if (co_pi0__[k + 70] < pi0mn || co_pi0__[k + 70] > pi0mx) {

/* ****             this value does not belong in this bin. */

		    ++(*ipi0cnt);
		    goto L1712;
		}
		abserr += (r__1 = co_pi0__[k + 70] - pi0grp[k + (m + 512) * 
			70], dabs(r__1)) / pi0_rng__[k + m * 70 - 71];
	    }
/* L1401: */
	}

/* ****        check single scattering co-albedos in layers above tau = 1 */

	for (k = levtau1[m + (j << 9)] - 1; k >= 1; --k) {

/* ****            if the single scattering co-albedo is larger than */
/*                a minimum value -scattering is important */

	    if ((1.f - co_pi0__[k + 70]) * dtauex[k + 70] > *taumn) {

/* ****             the scattering optical depth is finite. */
/*                 define the single-scattering co-albedo limits */

		pi0mn = pi0grp[k + (m + 1536) * 70] - pi0_rng__[k + m * 70 - 
			71];
		pi0mx = pi0grp[k + (m + 1024) * 70] + pi0_rng__[k + m * 70 - 
			71];
		if (co_pi0__[k + 70] < pi0mn || co_pi0__[k + 70] > pi0mx) {

/* ****             this value does not belong in this bin. */

		    ++(*ipi0cnt);
		    goto L1712;
		}
		abserr += (r__1 = co_pi0__[k + 70] - pi0grp[k + (m + 512) * 
			70], dabs(r__1)) / pi0_rng__[k + m * 70 - 71];
	    }
/* L1421: */
	}

/* ****        check scattering phase function in layers below tau = 1 */

	i__1 = *nlyr;
	for (k = levtau1[m + (j << 9)]; k <= i__1; ++k) {

/* *****          determine if scattering is important */

	    if ((1.f - co_pi0__[k + 70]) * dtauex[k + 70] > *taumn) {

		gmn = ggrp[k + (m + 1536) * 70] - g_rng__[k + m * 70 - 71];
		gmx = ggrp[k + (m + 1024) * 70] + g_rng__[k + m * 70 - 71];
		if (g[k] < gmn || g[k] > gmx) {

/* ****             this value does not belong in this bin. */

		    ++(*igcnt);
		    goto L1712;
		}

/* ****                check each term in the phase function expansion */

		i__2 = *nstr;
		for (mom = 1; mom <= i__2; ++mom) {
		    ph_mn__ = pmomgrp[mom + (k + (m + 1536) * 70) * 201] - 
			    ph_rng__[mom + (k + m * 70) * 200 - 14201];
		    ph_mx__ = pmomgrp[mom + (k + (m + 1024) * 70) * 201] + 
			    ph_rng__[mom + (k + m * 70) * 200 - 14201];
		    if (phmom[mom + k * 201] < ph_mn__ || phmom[mom + k * 201]
			     > ph_mx__) {

/* ****                       this value does not belong in this bin. */

			++(*igcnt);
			goto L1712;
		    }
/* L1501: */
		}

	    }
/* L1521: */
	}

/* ****        check scattering phase function in layers above tau = 1 */

	for (k = levtau1[m + (j << 9)] - 1; k >= 1; --k) {

/* *****          determine if scattering is important */

	    if ((1.f - co_pi0__[k + 70]) * dtauex[k + 70] > *taumn) {

		gmn = ggrp[k + (m + 1536) * 70] - g_rng__[k + m * 70 - 71];
		gmx = ggrp[k + (m + 1024) * 70] + g_rng__[k + m * 70 - 71];
		if (g[k] < gmn || g[k] > gmx) {

/* ****             this value does not belong in this bin. */

		    ++(*igcnt);
		    goto L1712;
		}

/* ****                check each term in the phase function expansion */

		i__1 = *nstr - 1;
		for (mom = 1; mom <= i__1; ++mom) {
		    ph_mn__ = pmomgrp[mom + (k + (m + 1536) * 70) * 201] - 
			    ph_rng__[mom + (k + m * 70) * 200 - 14201];
		    ph_mx__ = pmomgrp[mom + (k + (m + 1024) * 70) * 201] + 
			    ph_rng__[mom + (k + m * 70) * 200 - 14201];
		    if (phmom[mom + k * 201] < ph_mn__ || phmom[mom + k * 201]
			     > ph_mx__) {

/* ****                       this value does not belong in this bin. */

			++(*igcnt);
			goto L1712;
		    }
/* L1601: */
		}

	    }
/* L1621: */
	}

/* ****       congratulations! you passed all of the tests and you */
/*           qualify to enter this bin */

	if (abserr < errtst) {
	    errtst = abserr;
	    nbin = m;
	}
L1712:
	;
    }

    if (nbin != 0) {
	m = nbin;

/* ****       congratulations! you passed all of the tests and you */
/*           qualify to enter this bin */

	*iwngrp = m;
	dng0 = dnugrp[m] + 1e-7;
	dnugrp[m] += *delnu;
	dng1 = dnugrp[m] + 1e-7;
	wngrp[m + 512] = (dng0 * wngrp[m + 512] + *delnu * *wn) / dng1;
	wngrp[m + 1536] = *wn;

/* ****       modify min and max albedo values */

	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    albgrp[m + (nz + 10 << 9)] = (real) ((dng0 * albgrp[m + (nz + 10 
		    << 9)] + *delnu * alb[nz]) / dng1);
	    if (alb[nz] < albgrp[m + (nz + 20 << 9)]) {
		albgrp[m + (nz + 20 << 9)] = alb[nz];
	    }
	    if (alb[nz] > albgrp[m + (nz + 30 << 9)]) {
		albgrp[m + (nz + 30 << 9)] = alb[nz];
	    }
/* L2001: */
	}
	i__1 = *nref;
	for (nr = 1; nr <= i__1; ++nr) {
	    surfgrp[m + (nr + 4 << 9)] = (real) ((dng0 * surfgrp[m + (nr + 4 
		    << 9)] + *delnu * surf_opt__[nr]) / dng1);
	    if (surf_opt__[nr] < surfgrp[m + (nr + 8 << 9)]) {
		surfgrp[m + (nr + 8 << 9)] = surf_opt__[nr];
	    }
	    if (surf_opt__[nr] > surfgrp[m + (nr + 12 << 9)]) {
		surfgrp[m + (nr + 12 << 9)] = surf_opt__[nr];
	    }
/* L2021: */
	}

	i__1 = *nlay;
	for (k = 1; k <= i__1; ++k) {

	    taugrp[k + (m + 512) * 70] = (real) ((dng0 * taugrp[k + (m + 512) 
		    * 70] + *delnu * dtauex[k + 70]) / dng1);

/* ****           update min and max optical depth values */

	    if (dtauex[k + 70] < taugrp[k + (m + 1024) * 70]) {
		taugrp[k + (m + 1024) * 70] = dtauex[k + 70];
	    }

/* ****           load maximum optical depth and appropriate */
/*               single scattering abledo */

	    if (dtauex[k + 70] > taugrp[k + (m + 1536) * 70]) {
		taugrp[k + (m + 1536) * 70] = dtauex[k + 70];
	    }

/*               use a running average to add co-pi0 contribution. */

	    pi0grp[k + (m + 512) * 70] = (real) ((dng0 * pi0grp[k + (m + 512) 
		    * 70] + *delnu * co_pi0__[k + 70]) / dng1);

/* ****           update min and max co-albedo values */

	    if (co_pi0__[k + 70] < pi0grp[k + (m + 1024) * 70]) {
		pi0grp[k + (m + 1024) * 70] = co_pi0__[k + 70];
	    }
	    if (co_pi0__[k + 70] > pi0grp[k + (m + 1536) * 70]) {
		pi0grp[k + (m + 1536) * 70] = co_pi0__[k + 70];
	    }

/* ****            update the asymmetry factor */

	    ggrp[k + (m + 512) * 70] = (real) ((dng0 * ggrp[k + (m + 512) * 
		    70] + *delnu * g[k]) / dng1);

/* ****            update moments of phase function */

	    i__2 = *nmom;
	    for (mom = 0; mom <= i__2; ++mom) {
		pmomgrp[mom + (k + (m + 512) * 70) * 201] = (real) ((dng0 * 
			pmomgrp[mom + (k + (m + 512) * 70) * 201] + *delnu * 
			phmom[mom + k * 201]) / dng1);
/* L2101: */
	    }

/* ****           update min and max asymmetry factor values */
/*               NOTE: the phase function moments are slaved to the */
/*               variations in the asymmetry parameter */

	    if (g[k] < ggrp[k + (m + 1024) * 70]) {
		ggrp[k + (m + 1024) * 70] = g[k];
		i__2 = *nmom;
		for (mom = 1; mom <= i__2; ++mom) {
		    pmomgrp[mom + (k + (m + 1024) * 70) * 201] = phmom[mom + 
			    k * 201];
/* L2121: */
		}
	    }
	    if (g[k] > ggrp[k + (m + 1536) * 70]) {
		ggrp[k + (m + 1536) * 70] = g[k];
		i__2 = *nmom;
		for (mom = 1; mom <= i__2; ++mom) {
		    pmomgrp[mom + (k + (m + 1536) * 70) * 201] = phmom[mom + 
			    k * 201];
/* L2141: */
		}
	    }

/* L2201: */
	}

/* ****       check to see if the number of phase function moments */
/*           must be updated. */

	if (*nmom > nmomgrp[m]) {
	    nmomgrp[m] = *nmom;
	}

/* ****       you are done. go to the next monochromatic spectral interval */

	return 0;

    }

/* ****   sorry, you're a real misfit. start a bin of your own. */

/* ****   determine if number of bins exceeds dimension bound. */

    if (*ngroup == 512) {

/* ****     no more groups will fit - return */

	*ismterr = 1;
	return 0;

    }

/* ****  increment the group counter */

    ++(*ngroup);

/* ****   initialize bin properties */

    m = *ngroup;
    *iwngrp = *ngroup;
    nmomgrp[m] = *nmom;
    dnugrp[m] = *delnu;
    tautot[0] = 1e-7f;
    tautot[*nlyr + 70] = 1e-7f;
    for (i__ = 1; i__ <= 3; ++i__) {
	wngrp[m + (i__ << 9)] = *wn;
	i__1 = *nref;
	for (nr = 1; nr <= i__1; ++nr) {
	    surfgrp[m + (nr + (i__ << 2) << 9)] = surf_opt__[nr];
/* L3001: */
	}
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    albgrp[m + (nz + i__ * 10 << 9)] = alb[nz];
/* L3021: */
	}
/* L3041: */
    }

    i__1 = *nlay;
    for (k = 1; k <= i__1; ++k) {

/* ****       define the column-integrated absorption optical depth */

	tautot[k] = tautot[k - 1] + co_pi0__[k + 70] * dtauex[k + 70];
	tautot[*nlay - k + 70] = tautot[*nlay - k + 71] + co_pi0__[*nlay - k 
		+ 71] * dtauex[*nlay - k + 71];
/* L3201: */
    }

    i__1 = *nlay;
    for (k = 1; k <= i__1; ++k) {

/* ****        initialize optical depth, single scattering co-albedo, */
/*            assymetry parameter and phase function moments for bin */

	taugrp[k + (m + 512) * 70] = dtauex[k + 70];

	pi0grp[k + (m + 512) * 70] = co_pi0__[k + 70];

	ggrp[k + (m + 512) * 70] = g[k];

	i__2 = *nmom;
	for (mom = 0; mom <= i__2; ++mom) {
	    pmomgrp[mom + (k + (m + 512) * 70) * 201] = phmom[mom + k * 201];
/* L3401: */
	}
/* L3421: */
    }

/* ****   find the index of the tau(abs) =1 level for this bin */

    levtau1[m + 512] = *nlyr;
    levtau1[m + 1024] = 1;
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	if (levtau1[m + 512] == *nlyr && tautot[k] >= 1.f) {
	    levtau1[m + 512] = k;
	}
	if (levtau1[m + 1024] == 1 && tautot[*nlyr - k + 70] >= 1.f) {
	    levtau1[m + 1024] = *nlyr - k + 1;
	}
/* L3601: */
    }

/* *****  define the optical depth range for this bin. */
/*       the criteria is much weaker for large absorption */
/*       optical depths. */

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {

/* ****        define the absorption optical depth */

	tauabs = pi0grp[k + (m + 512) * 70] * taugrp[k + (m + 512) * 70];

/* ****       define the layer-dependent optical depth error for bin. */
/*           the following function varies from ~tauerr at tau < 1, */
/*           to a much larger value at large layer absorption optical */
/*           depths (ie for completely opaque layers). */

	tauerror = taugrp[k + (m + 512) * 70] * *tauerr * (tauabs * .1f + 1.f 
		+ tauabs * .01f * tauabs) + *taumn * 10.f;

/* ****       define the nominal optical depth limits */

	taumax0 = taugrp[k + (m + 512) * 70] + tauerror;
	taumin0 = taugrp[k + (m + 512) * 70] * (1.f - *tauerr) - *taumn;
	tau_rng__[k + m * 70 - 71] = taumax0 - taumin0;
	taugrp[k + (m + 1536) * 70] = taugrp[k + (m + 512) * 70];
	taugrp[k + (m + 1024) * 70] = taugrp[k + (m + 512) * 70];
	dtaumn[k + m * 70 - 71] = taugrp[k + (m + 512) * 70] * (1.f - *tauerr)
		 - *taumn;

/* ****       define the single scattering co-albedo limits */

	pi0max0 = pi0grp[k + (m + 512) * 70] * (*pi0err + 1.f) + smallpi0;
	pi0min0 = pi0grp[k + (m + 512) * 70] * (1.f - *pi0err) - smallpi0;
	pi0_rng__[k + m * 70 - 71] = pi0max0 - pi0min0;
	pi0grp[k + (m + 1536) * 70] = pi0grp[k + (m + 512) * 70];
	pi0grp[k + (m + 1024) * 70] = pi0grp[k + (m + 512) * 70];

/* ****       define the asymmetry parameter limits */

	gmax0 = ggrp[k + (m + 512) * 70] * (*phferr + 1.f) + smallg;
	gmin0 = ggrp[k + (m + 512) * 70] * (1.f - *phferr) - smallg;
	g_rng__[k + m * 70 - 71] = gmax0 - gmin0;
	ggrp[k + (m + 1536) * 70] = ggrp[k + (m + 512) * 70];
	ggrp[k + (m + 1024) * 70] = ggrp[k + (m + 512) * 70];

/* ****       define the phase function limits. */

	pmomgrp[(k + (m + 1024) * 70) * 201] = 1.f;
	pmomgrp[(k + (m + 1536) * 70) * 201] = 1.f;
	i__2 = *nmom;
	for (mom = 1; mom <= i__2; ++mom) {
	    ph_max0__ = pmomgrp[mom + (k + (m + 512) * 70) * 201] * (*phferr 
		    + 1.f) + smallg;
	    ph_min0__ = pmomgrp[mom + (k + (m + 512) * 70) * 201] * (1.f - *
		    phferr) - smallg;
	    ph_rng__[mom + (k + m * 70) * 200 - 14201] = ph_max0__ - 
		    ph_min0__;
	    pmomgrp[mom + (k + (m + 1536) * 70) * 201] = pmomgrp[mom + (k + (
		    m + 512) * 70) * 201];
	    pmomgrp[mom + (k + (m + 1024) * 70) * 201] = pmomgrp[mom + (k + (
		    m + 512) * 70) * 201];
/* L4001: */
	}

/* L4021: */
    }

/* ****     define surface optical property limits */

    i__1 = *nref;
    for (nr = 1; nr <= i__1; ++nr) {
	albmax0 = albgrp[m + (nr + 10 << 9)] * (*surferr + 1.f) + smallalb;
	albmin0 = albgrp[m + (nr + 10 << 9)] * (1.f - *surferr) - smallalb;
	alb_rng__[m + (nr << 9) - 513] = albmax0 - albmin0;
/* L4041: */
    }

/* ****  get the next spectral interval */

    return 0;
} /* map_spect__ */

