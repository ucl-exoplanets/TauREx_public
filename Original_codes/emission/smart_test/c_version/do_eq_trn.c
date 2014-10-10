/* do_eq_trn.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int do_eq_trn__(integer *nlyr, integer *nstr, integer *nmom, 
	integer *numu, integer *nphi, integer *ibcnd, integer *iunits, 
	integer *iref, logical *usrang, logical *lamber, logical *source, 
	real *umu, real *phi, real *umu0, real *phi0, real *accur, real *alb0,
	 real *surf_pr__, real *dtauc, real *ssalb, real *pmom, real *t, real 
	*ttemp, real *btemp, real *temis, real *fbeam, real *fisot, real *
	wng0, real *flup, real *rfldir, real *rfldn, real *uu, real *albmed, 
	real *trnmed, integer *nscat)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, l;
    static real dfdt[70], uavg[70];
    static integer ntau;
    static real utau[70];
    static logical prnt[7];
    static real wvnm;
    static logical plank;
    static integer maxphi;
    static real temper[71];
    static integer maxcly, maxmom;
    static logical onlyfl;
    extern /* Subroutine */ int disort_(integer *, real *, real *, integer *, 
	    real *, real *, real *, logical *, integer *, real *, integer *, 
	    integer *, logical *, integer *, real *, integer *, real *, 
	    integer *, real *, real *, real *, real *, logical *, integer *, 
	    real *, real *, real *, real *, real *, logical *, logical *, 
	    real *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *);
    static integer maxulv, maxumu;
    static logical usrtau;


/* ccccccccccccccccccccccccc  d o _ e q _ t r n  ccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this program initializes discrete ordinate arrays and calls     cc */
/* c    the Stamnes a discrete ordinate model to find radiances         cc */
/* c    and fluxes for smt optical properties.                          cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        ng0 - index  of current spectral bin                        cc */
/* c       nlyr - number of computational model layers                  cc */
/* c       nstr - number of gaussian zenith angles used in D/O code     cc */
/* c       nmom - number of phase function moments used in D/O code     cc */
/* c       numu - number of output zenith angles used in D/O code       cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c      ibcnd - boundy condition flag (0- flx/rad, 1-albedos only)    cc */
/* c     iunits - output units desired for planck function:             cc */
/* c           0: unit flux: b(k) = 1.0                                 cc */
/* c           1: Watts/m**2/cm**-1                                     cc */
/* c           2: Watts/m**2/micron                                     cc */
/* c           3: Watts/m**2/nanometer                                  cc */
/* c           4: Watts/m**2/Angstrom                                   cc */
/* c           5: Watts/m**2/Hz                                         cc */
/* c       iref - bidirectional reflectance options                     cc */
/* c              0 - lambert                                           cc */
/* c              1 - Hapke's BDR model                                 cc */
/* c              2 - Breon's BDR model; combination of Li + Roujean    cc */
/* c              3 - Roujean's BDR model                               cc */
/* c              4 - Cox and Munk glint model                          cc */
/* c     usrang - output radiances at user angles? (logical: T/F)       cc */
/* c     lamber - Include a lambertian surface? (Logical: T/F)          cc */
/* c              note: if lamber = F, use a BRDF is used.              cc */
/* c     source - include thermal fluxes? (logical: T/F)                cc */
/* c        umu - emission zenith angle cosines                         cc */
/* c        phi - emission azimuth angles (degrees)                     cc */
/* c       umu0 - cosine of solar zenith angles                         cc */
/* c       phi0 - solar azimuth angles (degrees)                        cc */
/* c      accur - azimuth convergence accuracy for D/O routine          cc */
/* c       alb0 - surface albedo                                        cc */
/* c    surf_pr - surface proerties for non-lambertian BRDF             cc */
/* c      dtauc - layer optical depth                                   cc */
/* c      ssalb - layer single scattering albedo                        cc */
/* c       pmom - layer particle phase function                         cc */
/* c          t - temperature in each atmospheric                       cc */
/* c      ttemp - temperature of uppermost model layer (space)          cc */
/* c      btemp - surface temperature                                   cc */
/* c      temis - emissivity of uppermost model layer (space)           cc */
/* c      fbeam - intensity of collimated flux at top of atmosphere     cc */
/* c      fisot - thermal flux at top of atmosphere                     cc */
/* c       wng0 - wavenumber of bin                                     cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c       flup - upward (diffuse) flux at each model level             cc */
/* c     rfldir - downward direct flux at each model level              cc */
/* c      rfldn - downward diffuse flux at each model level             cc */
/* c         uu - upward radiance at each level, zenith angle, azimuth  cc */
/* c     albmed - albedo of system                                      cc */
/* c     trnmed - transmissivity of system                              cc */
/* c      nscat - counter for number of scattering calculations         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  d o _ e q _ t r n  ccccccccccccccccccccccccc */




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






/* ***    output variables */


/* ****    set logicals */

    /* Parameter adjustments */
    --trnmed;
    --albmed;
    uu -= 1137;
    --rfldn;
    --rfldir;
    --flup;
    --t;
    pmom -= 201;
    --ssalb;
    --dtauc;
    --surf_pr__;
    --phi;
    --umu;

    /* Function Body */
    onlyfl = FALSE_;
    plank = *source;

/* ****         set dimensions of arrays */

    maxcly = 70;
    maxulv = 70;
    maxumu = 16;
    maxphi = 16;
    maxmom = 200;

/* ****       set wavelength interval (1 cm**-1) */

    wvnm = *wng0;

/* ****   turn off all discr_ord model internal print flags */

    for (l = 1; l <= 7; ++l) {
	prnt[l - 1] = FALSE_;
/* L1001: */
    }

/* ****    set number of arbitrary output levels and angles */

    usrtau = FALSE_;
    ntau = 0;

/* ****   load optical depth, single scattering albedo */
/*       and phase function moment and temperature arrays */

    temper[0] = t[1];
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {

	utau[k - 1] = 0.f;
	temper[k] = t[k + 1];
/* L1221: */
    }

    ++(*nscat);

/* ****   d i s c r e t e   o r d i n a t e    r o u t i n e */

    disort_(nlyr, &dtauc[1], &ssalb[1], nmom, &pmom[201], temper, &wvnm, &
	    usrtau, &ntau, utau, nstr, iunits, usrang, numu, &umu[1], nphi, &
	    phi[1], ibcnd, fbeam, umu0, phi0, fisot, lamber, iref, &surf_pr__[
	    1], alb0, btemp, ttemp, temis, &plank, &onlyfl, accur, prnt, &
	    maxcly, &maxulv, &maxumu, &maxphi, &maxmom, &rfldir[1], &rfldn[1],
	     &flup[1], dfdt, uavg, &uu[1137], &albmed[1], &trnmed[1]);

    return 0;
} /* do_eq_trn__ */

