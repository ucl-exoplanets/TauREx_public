/* source_fcn.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int source_fcn__(integer *ng0, integer *l0, integer *nlyr, 
	integer *nzdn, integer *nzup, integer *nstr, integer *nmom, integer *
	numu, integer *nphi, integer *iref, logical *usrang, logical *lamber, 
	logical *source, real *umu, real *phi, real *umu0, real *phi0, real *
	accur, real *alb0, real *surf_pr__, real *dtauc, real *ssalb, real *
	pmom, real *ts, real *t, real *fbeam, real *wng0, doublereal *
	trnflx_b__, doublereal *refflx_b__, doublereal *trnrad_b__, 
	doublereal *flx_rt__, doublereal *flx_rb__, doublereal *flx_du__, 
	doublereal *flx_dd__, doublereal *flx_dft__, doublereal *flx_uft__, 
	doublereal *flx_dfb__, doublereal *flx_ufb__, doublereal *rad_rt__, 
	doublereal *rad_rb__, doublereal *dnflxsrc, doublereal *upflxsrc, 
	doublereal *radsrc, real *rfldir, real *rfldn, doublereal *flx_ft__, 
	doublereal *flx_fb__, real *uu, integer *nscat)
{
    static integer l;
    extern /* Subroutine */ int do_eq_trn__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, logical *, 
	    logical *, logical *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, integer *);
    static integer ibn;
    static real flup[70];
    static integer ibcnd;
    static real btemp, temis, fisot, ttemp, albmed[16], trnmed[16];
    static integer iunits;
    extern /* Subroutine */ int rad_src__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, doublereal *), flx_src__(
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, real *, 
	    real *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


/* cccccccccccccccccccccccc  s o u r c e _ f c n   ccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes the radiance and flux source functions cc */
/* c    for a given bin for smart                                       cc */
/* c                                                                    cc */
/* c    Note: in the current formulation, if usrang = .true., each      cc */
/* c          downward / upward stream much be accompanied by its       cc */
/* c          conjugate, such that ang(nze) = 180 - ang(numu-nze+1)=    cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        ng0 - index  of current spectral bin                        cc */
/* c         l0 - index of state vector sounding                        cc */
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
/* c     rfldir - downward direct solar flux at each level for this bin cc */
/* c   dnflxsrc - downward flux source at each level for this bin       cc */
/* c   upflxsrc - downward flux source at each level for this bin       cc */
/* c     radsrc - angle-dependent radiance source at each level for     cc */
/* c              this bin                                              cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccc  s o u r c e _ f c n   ccccccccccccccccccccccc */




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



/*      integer nze,naz */





/* ***    output variables */


/* ****   binned layer flux transmittances and absorptances */


/* ****   variables for flux layer adding method */


/* ****  quantities for radiance calculation: downward and upward fluxes */
/*      at layer interfaces */


/* ****   binned layer radiance transmittances and absorptances */


/* ****   variables for radiance layer adding method */


/* ****   internal binned source function variables */

/*      integer k */

    /* Parameter adjustments */
    uu -= 1137;
    flx_fb__ -= 71;
    flx_ft__ -= 71;
    --rfldn;
    --rfldir;
    radsrc -= 18193;
    upflxsrc -= 71;
    dnflxsrc -= 71;
    rad_rb__ -= 18193;
    rad_rt__ -= 18193;
    flx_ufb__ -= 71;
    flx_dfb__ -= 71;
    flx_uft__ -= 71;
    flx_dft__ -= 71;
    flx_dd__ -= 71;
    flx_du__ -= 71;
    flx_rb__ -= 71;
    flx_rt__ -= 71;
    trnrad_b__ -= 6737;
    refflx_b__ -= 421;
    trnflx_b__ -= 421;
    --t;
    pmom -= 201;
    --ssalb;
    --dtauc;
    --surf_pr__;
    --phi;
    --umu;

    /* Function Body */
    ibn = *ng0;
    l = *l0;

/* ****     define thermal source function units */

    iunits = 1;
    ibcnd = 0;

/* ****            set surface temperature and albedo */

    btemp = *ts;

/* ****             set top boundary condition - */
/*                 no emission or reflection */

    fisot = 0.f;
    ttemp = 0.f;
    temis = 1.f;

/* ****     call discrete ordinate driver */

    do_eq_trn__(nlyr, nstr, nmom, numu, nphi, &ibcnd, &iunits, iref, usrang, 
	    lamber, source, &umu[1], &phi[1], umu0, phi0, accur, alb0, &
	    surf_pr__[1], &dtauc[1], &ssalb[1], &pmom[201], &t[1], &ttemp, &
	    btemp, &temis, fbeam, &fisot, wng0, flup, &rfldir[1], &rfldn[1], &
	    uu[1137], albmed, trnmed, nscat);
/*      write(*,'(/,1a,2(1pe14.6),l2))') 'source_fcn: wng0,alb0', */
/*     -                        wng0,alb0 */
/*      write(*,*) nlyr,nstr,nmom,numu,nphi,ibcnd, */
/*     -           iunits,iref,usrang,lamber,umu0,phi0,accur, */
/*     -           ttemp,btemp,temis,fbeam,fisot,(surf_pr(k),k=1,4) */
/*     - */
/*      write(*,'(1a,16(1pe12.4))') 'dtauc ',(dtauc(k),k=1,nlyr) */
/*      write(*,'(1a,16(1pe12.4))') 'ssalb ',(ssalb(k),k=1,nlyr) */
/*      write(*,'(1a,16(1pe12.4))') 'pmom  ',(pmom(1,k),k=1,nlyr) */
/*      write(*,'(1a,16(1pe12.4))') 't     ',(t(k),k=1,nlyr) */
/*      write(*,'(1a,16(1pe12.4))') 'rfldir',(rfldir(k),k=1,nlyr+1) */
/*      write(*,'(1a,16(1pe12.4))') 'rfldn ',(rfldn(k),k=1,nlyr+1) */
/*      write(*,'(1a,16(1pe12.4))') 'flup  ',(flup(k),k=1,nlyr+1) */

/* ****   use the flux adding method to solve for the */
/*       upward flux at the top of each layer, f+, and the */
/*       downward flux at the base of each layer, f-. */

/* ****   find layer source functions for fluxes */

    flx_src__(ng0, l0, nlyr, &trnflx_b__[421], &refflx_b__[421], &flx_rt__[71]
	    , &flx_rb__[71], &flx_du__[71], &flx_dd__[71], &rfldn[1], flup, &
	    flx_fb__[71], &flx_ft__[71], &flx_dft__[71], &flx_uft__[71], &
	    flx_dfb__[71], &flx_ufb__[71], &upflxsrc[71], &dnflxsrc[71]);

/* ****   find radiance source terms */

    rad_src__(ng0, l0, nlyr, nphi, nzdn, numu, &trnrad_b__[6737], &rad_rt__[
	    18193], &rad_rb__[18193], &flx_ft__[71], &flx_fb__[71], &
	    flx_dfb__[71], &flx_uft__[71], &rfldn[1], flup, &uu[1137], &
	    radsrc[18193]);

    return 0;
} /* source_fcn__ */

