/* sm_eq_trn.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sm_eq_trn__(logical *usrang, logical *lamber, logical *
	lplanck, logical *lsolar, integer *ng0, integer *nz0, integer *nlyr, 
	integer *nstr, integer *numu, integer *nzdn, integer *nzup, integer *
	nphi, integer *nmomgrp, integer *nstate, integer *istate, integer *
	n_rad__, integer *iref, integer *nref, integer *nscat, real *tauerr, 
	real *pi0err, real *phferr, real *surferr, real *ws, real *phiw, real 
	*taumn, real *umu, real *umu_f__, real *gwt_f__, real *phi, real *
	umu0nz, real *phi0nz, real *ts, real *t, real *accur, doublereal *
	wngrp, real *surfgrp, real *taugrp, real *pi0grp, real *ggrp, real *
	pmomgrp, real *alb_b__, real *sur_b__, real *dtau_b__, real *
	copi0_b__, real *g_b__, real *pmom_b__, real *dx_b_i__, real *
	dalb_b_i__, doublereal *trnflx_b__, doublereal *refflx_b__, 
	doublereal *absflx_b__, doublereal *dtrnflxdx, doublereal *drefflxdx, 
	doublereal *dabsflxdx, doublereal *trnrad_b__, doublereal *refrad_b__,
	 doublereal *absrad_b__, doublereal *brdf_b__, doublereal *dtrnraddx, 
	doublereal *drefraddx, doublereal *dabsraddx, doublereal *dbrdfdx, 
	doublereal *dnsflxsrc_b__, doublereal *upsflxsrc_b__, doublereal *
	sradsrc_b__, doublereal *ddnsflxdx, doublereal *dupsflxdx, doublereal 
	*dsraddx, doublereal *dntflxsrc_b__, doublereal *uptflxsrc_b__, 
	doublereal *tradsrc_b__, doublereal *ddntflxdx, doublereal *duptflxdx,
	 doublereal *dtraddx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal dnflxsrc[350]	/* was [70][5] */, upflxsrc[350]	
	    /* was [70][5] */;
    static integer k, l;
    extern /* Subroutine */ int surf_brdf__(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, doublereal *, doublereal *);
    static integer l0;
    extern /* Subroutine */ int solar_src__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal dnsflxsrc[350]	/* was [70][5] */, dntflxsrc[350]	
	    /* was [70][5] */, upsflxsrc[350]	/* was [70][5] */, uptflxsrc[
	    350]	/* was [70][5] */;
    static integer nr, nz;
    static real uu[17920]	/* was [16][70][16] */;
    extern /* Subroutine */ int source_fcn__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, logical *, logical *, logical *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, doublereal *, doublereal *, real *, 
	    integer *);
    static integer ibn, mom, naz, nze;
    static real alb0, phi0;
    extern /* Subroutine */ int thermal_src__(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     real *, real *, doublereal *, real *, real *, real *, real *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static real wng0, umu0;
    static integer nlev, nmom;
    static real pmom[14070]	/* was [201][70] */, fbeam, dtauc[70], ssalb[
	    70], rfldn[70];
    static doublereal rad_rb__[89600]	/* was [16][16][70][5] */, flx_dd__[
	    350]	/* was [70][5] */, flx_fb__[350]	/* was [70][5]
	     */, rad_rt__[89600]	/* was [16][16][70][5] */, flx_rb__[
	    350]	/* was [70][5] */, radsrc[89600]	/* was [16][
	    16][70][5] */, flx_du__[350]	/* was [70][5] */;
    static real rfldir[70];
    static doublereal flx_ft__[350]	/* was [70][5] */;
    static logical source;
    static doublereal flx_rt__[350]	/* was [70][5] */;
    extern /* Subroutine */ int rad_ref__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal flx_dfb__[350]	/* was [70][5] */, flx_drb__[350]	
	    /* was [70][5] */, flx_ufb__[350]	/* was [70][5] */, flx_dft__[
	    350]	/* was [70][5] */;
    extern /* Subroutine */ int flx_ref__(integer *, integer *, integer *, 
	    real *, doublereal *, doublereal *, doublereal *, real *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal sradsrc[89600]	/* was [16][16][70][5] */, tradsrc[
	    89600]	/* was [16][16][70][5] */, flx_uft__[350]	/* 
	    was [70][5] */;
    static real surf_pr__[4];
    static doublereal flx_urt__[350]	/* was [70][5] */;


/* cccccccccccccccccccccccc  s m _ e q _ t r n   ccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this routine uses a discrete ordinate method to find radiances  cc */
/* c    and fluxes for spectrally mapped optical properties.            cc */
/* c                                                                    cc */
/* c 3/18/07: This version of sm_eq_trn uses flux layer adding  for     cc */
/* c          both fluxes and radiances.  The same source function are  cc */
/* c          used for both solar and thermal regions.                  cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        ng0 - index  of current spectral bin                        cc */
/* c       nlyr - number of computational model layers                  cc */
/* c       nstr - number of gaussian zenith angles used in D/O code     cc */
/* c       numu - number of output zenith angles used in D/O code       cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c     usrang - output radiances at user angles? (logical: T/F)       cc */
/* c        umu - emission zenith angle cosines                         cc */
/* c        phi - emission azimuth angles (degrees)                     cc */
/* c        nz0 - index of current solar zenith angle (1 to nza)        cc */
/* c       umu0 - cosine of solar zenith angles                         cc */
/* c       phi0 - solar azimuth angles (degrees)                        cc */
/* c      accur - azimuth convergence accuracy for D/O routine          cc */
/* c     lamber - Include a lambertian surface? (Logical: T/F)          cc */
/* c              note: if lamber = F, use a BRDF is used.              cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c      icbnd - disort variable (0- radiances, 1- albedos only)       cc */
/* c      fbeam - intensity of collimated flux at top of atmosphere     cc */
/* c      fisot - thermal flux at top of atmosphere                     cc */
/* c        wng - wavenumber of bin                                     cc */
/* c         ts - surface temperature (K)                               cc */
/* c     nstate - number of variable elements in the state vecror       cc */
/* c     istate - state vector flag indicating which state variables    cc */
/* c              are variable components of the state vector.          cc */
/* c              1 - surface pressure                                  cc */
/* c              2 - surface/atmospheric temperature                   cc */
/* c              3 - gas absorption coeffient                          cc */
/* c              4 - cloud/aerosol optical depth                       cc */
/* c              5 - surface albedo                                    cc */
/* c       iref - bidirectional reflectance options                     cc */
/* c            0 - Lambert                                             cc */
/* c            1 - Hapke's BDR model                                   cc */
/* c            2 - Breon's BDR model; combination of Li + Roujean      cc */
/* c            3 - Roujean's BDR model                                 cc */
/* c            4 - Cox and Munk glint model                            cc */
/* c         ws - wind speed (m/s) for Cox/Munk model                   cc */
/* c       phiw - wind azimuth (deg) for Cox/Munk model                 cc */
/* c      n_rad - radiance calculations flag for each bin               cc */
/* c             1 : radiances for nominal state structure              cc */
/* c             2 : radiances for perturbed surface pressures          cc */
/* c             3 : radiances for perturbed optical depths             cc */
/* c             4 : radiances for perturbed single scattering albedos  cc */
/* c             5 : radiances for perturbed surface reflectance        cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      temis - emissivity of uppermost model layer (space)           cc */
/* c      ttemp - temperature of uppermost model layer (space)          cc */
/* c      sur_b - surface albedo                                        cc */
/* c     dtau_b - differential optical depth in each layer              cc */
/* c    copi0_b - single scattering co-albedo in each layer             cc */
/* c      btemp - surface temperature                                   cc */
/* c     rfldir - downward direct flux at each model level              cc */
/* c      rfldn - downward diffuse flux at each model level             cc */
/* c         uu - upward radiance at each level, zenith angle, azimuth  cc */
/* c     albmed - albedo of system                                      cc */
/* c     trnmed - transmissivity of system                              cc */
/* c      nscat - counter for number of scattering calculations         cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccc  s m _ e q _ t r n   ccccccccccccccccccccccccc */



/* ****   define variables used by discrete ordinate method */


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



/*      integer nlyr0 */
/*      integer nze,naz */

/* ****   spectral binning parameters */


/* ****   spectral binning parameters */



/* ****   emission angle variables */


/* ****   DISORT input variables */


/* ****   DISORT output */


/* *****   state vector variables. */



/* ****   binned layer flux transmittances and absorptances */



/* ****   binned layer radiance transmittances and absorptances */


/* ****   variables for radiance layer adding method */


/* ****   binned solar source function variables and jacobians */



/* ****   binned thermal source function variables and jacobians */



/* ****   internal binned source function variables */


/* ****   internal binned source function variables */


/* ****   source_fcn output variables */


/* ****  quantities for radiance calculation: downward and upward fluxes */
/*      at layer interfaces */


/* ****   specify number of levels, and solar zenith angle */

    /* Parameter adjustments */
    dtraddx -= 89873;
    duptflxdx -= 351;
    ddntflxdx -= 351;
    tradsrc_b__ -= 18193;
    uptflxsrc_b__ -= 71;
    dntflxsrc_b__ -= 71;
    dsraddx -= 89873;
    dupsflxdx -= 351;
    ddnsflxdx -= 351;
    sradsrc_b__ -= 18193;
    upsflxsrc_b__ -= 71;
    dnsflxsrc_b__ -= 71;
    dbrdfdx -= 273;
    dabsraddx -= 4497;
    drefraddx -= 4497;
    dtrnraddx -= 4497;
    brdf_b__ -= 1553;
    absrad_b__ -= 6737;
    refrad_b__ -= 6737;
    trnrad_b__ -= 6737;
    dabsflxdx -= 281;
    drefflxdx -= 281;
    dtrnflxdx -= 281;
    absflx_b__ -= 421;
    refflx_b__ -= 421;
    trnflx_b__ -= 421;
    dalb_b_i__ -= 11;
    dx_b_i__ -= 421;
    pmom_b__ -= 84621;
    g_b__ -= 421;
    copi0_b__ -= 421;
    dtau_b__ -= 421;
    sur_b__ -= 25;
    alb_b__ -= 61;
    pmomgrp -= 7218111;
    ggrp -= 35911;
    pi0grp -= 35911;
    taugrp -= 35911;
    surfgrp -= 2561;
    wngrp -= 513;
    --t;
    --phi;
    --gwt_f__;
    --umu_f__;
    --umu;
    n_rad__ -= 6;
    --istate;
    --nmomgrp;

    /* Function Body */
    nlev = *nlyr + 1;
    nz = *nz0;

    if (*lsolar) {
	umu0 = *umu0nz;
	phi0 = *phi0nz;
    } else {
	umu0 = 0.f;
	phi0 = 0.f;
    }

/* ****    set group counters */

    ibn = *ng0;
    nmom = nmomgrp[ibn];

    wng0 = (real) wngrp[ibn + 512];

/* ****   find the layer reflectance quantities for fluxes and radiances */
/*       for the simplified adding method */

    for (l = 1; l <= 5; ++l) {

	l0 = l;

	if (l == 1 || n_rad__[l + ibn * 5] != 0) {

/* ****       define the surface optical properties */

	    alb0 = alb_b__[nz + (l + ibn * 5) * 10];
	    trnflx_b__[nlev + (l + ibn * 5) * 70] = 0.f;
	    refflx_b__[nlev + (l + ibn * 5) * 70] = alb0;

	    if (*lamber) {

		i__1 = *nphi;
		for (naz = 1; naz <= i__1; ++naz) {
		    i__2 = *numu;
		    for (nze = 1; nze <= i__2; ++nze) {
			brdf_b__[nze + (naz + (l + ibn * 5 << 4) << 4)] = 
				alb0;
			if (l == 5) {
			    dbrdfdx[nze + (naz + (ibn << 4) << 4)] = (alb0 - 
				    brdf_b__[nze + (naz + (ibn * 5 + 1 << 4) 
				    << 4)]) * dalb_b_i__[nz + ibn * 10];
			}
/* L1221: */
		    }
/* L1241: */
		}

	    } else {

/* ****          use disort to find brdf of surface */

		surf_brdf__(usrang, lamber, &l0, ng0, nz0, nstr, numu, nphi, &
			nmom, &n_rad__[6], nref, &umu[1], &phi[1], &
			dalb_b_i__[11], &alb_b__[61], &sur_b__[25], ws, phiw, 
			&umu0, &phi0, &brdf_b__[1553], &dbrdfdx[273]);

	    }


/* ****         find the layer transmittances, reflectances and */
/*             absorptances for combined layers */


	    flx_ref__(&l0, nlyr, ng0, &alb0, &trnflx_b__[421], &refflx_b__[
		    421], &absflx_b__[421], &dx_b_i__[421], flx_dd__, 
		    flx_du__, flx_urt__, flx_drb__, flx_rb__, flx_rt__);

/* ****          find the radiance sources */

	    rad_ref__(&l0, nlyr, ng0, numu, nzdn, nzup, nphi, &trnrad_b__[
		    6737], &refrad_b__[6737], &brdf_b__[1553], flx_urt__, 
		    flx_drb__, rad_rb__, rad_rt__);

	}

/* L1401: */
    }

/* ******  o p t i c a l    d e p t h   p r o f i l e   l o o p   ******* */

/* ****   perform radiance calculations for the nomical case and */
/*       for each variable component of the optical depth structure */
/*       and set the radiance calculation flag, n_rad, where */
/*       n_rad = 1 : radiances for nominal state structure */
/*             = 2 : radiances for perturbed optical depths */
/*             = 3 : radiances for perturbed single scattering albedos */
/*             = 4 : radiances for perturbed phase functions */
/*             = 5 : radiances for perturbed surface albedos */

    for (l = 1; l <= 5; ++l) {
	l0 = l;

	if (l == 1 || n_rad__[l + ibn * 5] != 0) {

/* ****        define the surface optical properties */

	    alb0 = alb_b__[nz + (l + ibn * 5) * 10];

/* ****        set the surface optical properties for disort */

	    i__1 = *nref;
	    for (nr = 1; nr <= i__1; ++nr) {
		surf_pr__[nr - 1] = sur_b__[nr + (l + ibn * 5 << 2)];
/* L2021: */
	    }

	    if (! (*lamber)) {

		if (*iref == 4) {

/* ****            set wind speed and direction for Cox/Munk model */

		    surf_pr__[2] = *ws;
		    surf_pr__[3] = *phiw;
		}

	    }

/* ****            set number of l-p phase function moments */

	    nmom = nmomgrp[ibn];

	    i__1 = *nlyr;
	    for (k = 1; k <= i__1; ++k) {

/* ****            load the optical depth and single scattering albedos */
/*                for use in the discrete ordinate algorithm */

		dtauc[k - 1] = dtau_b__[k + (l + ibn * 5) * 70];
		ssalb[k - 1] = 1.f - copi0_b__[k + (l + ibn * 5) * 70];
		if (ssalb[k - 1] < 0.f) {
		    ssalb[k - 1] = 0.f;
		}

/* ****            load the phase function moments at each level */

		pmom[k * 201 - 201] = 1.f;
		i__2 = nmomgrp[ibn];
		for (mom = 1; mom <= i__2; ++mom) {
		    pmom[mom + k * 201 - 201] = pmom_b__[mom + (k + (l + ibn *
			     5) * 70) * 201];
/* L2121: */
		}
/* L2141: */
	    }

	}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

	if (*lsolar) {

/* ****         find solar fluxes, radiances and source functions */
/*             in each layer */

	    if (l == 1 || n_rad__[l + ibn * 5] != 0) {

/* ****             turn off thermal sources and set a unit direct flux */
/*                 at the top of the atmosphere */

		source = FALSE_;
		fbeam = 1.f;

/* ****             find the solar flux and radince sources in each layer */
/*                 note: source_fcn output variables dnsflxsrc, */
/*                 upsflxsrc, and sradsrc are specific to this call */

		source_fcn__(ng0, &l0, nlyr, nzdn, nzup, nstr, &nmom, numu, 
			nphi, iref, usrang, lamber, &source, &umu[1], &phi[1],
			 &umu0, &phi0, accur, &alb0, surf_pr__, dtauc, ssalb, 
			pmom, ts, &t[1], &fbeam, &wng0, &trnflx_b__[421], &
			refflx_b__[421], &trnrad_b__[6737], flx_rt__, 
			flx_rb__, flx_du__, flx_dd__, flx_dft__, flx_uft__, 
			flx_dfb__, flx_ufb__, rad_rt__, rad_rb__, dnflxsrc, 
			upflxsrc, radsrc, rfldir, rfldn, flx_ft__, flx_fb__, 
			uu, nscat);

	    }

	    solar_src__(ng0, &l0, nz0, nlyr, nzdn, nzup, numu, nphi, &n_rad__[
		    6], &alb0, &umu0, rfldir, rfldn, &dx_b_i__[421], &
		    dalb_b_i__[11], &trnrad_b__[6737], rad_rt__, rad_rb__, 
		    flx_ft__, flx_fb__, &trnflx_b__[421], flx_rt__, flx_rb__, 
		    flx_du__, flx_dd__, flx_dft__, flx_uft__, flx_dfb__, 
		    flx_ufb__, dnflxsrc, upflxsrc, radsrc, dnsflxsrc, 
		    upsflxsrc, sradsrc, &dnsflxsrc_b__[71], &upsflxsrc_b__[71]
		    , &sradsrc_b__[18193], &ddnsflxdx[351], &dupsflxdx[351], &
		    dsraddx[89873]);

	}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

	if (*lplanck) {

/* ****         find thermal fluxes, radiances and source functions */
/*             in each layer */

	    if (l == 1 || n_rad__[l + ibn * 5] != 0) {

/* ****             turn on thermal sources and set the direct flux */
/*                 at the top of the atmosphere to zero */

		source = TRUE_;
		fbeam = 0.f;

/*                 note: source_fcn output variables dntflxsrc, */
/*                 uptflxsrc, and tradsrc are specific to this call */

		source_fcn__(ng0, &l0, nlyr, nzdn, nzup, nstr, &nmom, numu, 
			nphi, iref, usrang, lamber, &source, &umu[1], &phi[1],
			 &umu0, &phi0, accur, &alb0, surf_pr__, dtauc, ssalb, 
			pmom, ts, &t[1], &fbeam, &wng0, &trnflx_b__[421], &
			refflx_b__[421], &trnrad_b__[6737], flx_rt__, 
			flx_rb__, flx_du__, flx_dd__, flx_dft__, flx_uft__, 
			flx_dfb__, flx_ufb__, rad_rt__, rad_rb__, dnflxsrc, 
			upflxsrc, radsrc, rfldir, rfldn, flx_ft__, flx_fb__, 
			uu, nscat);

	    }

/* ****             find thermal flux and radince sources in each layer */

	    thermal_src__(ng0, &l0, nz0, nlyr, nzdn, nzup, numu, nphi, &
		    n_rad__[6], &wng0, &alb0, &brdf_b__[1553], ts, &t[1], &
		    dx_b_i__[421], &dalb_b_i__[11], &trnrad_b__[6737], 
		    rad_rt__, rad_rb__, flx_uft__, flx_dfb__, dnflxsrc, 
		    upflxsrc, radsrc, dntflxsrc, uptflxsrc, tradsrc, &
		    copi0_b__[421], &absflx_b__[421], &absrad_b__[6737], &
		    dntflxsrc_b__[71], &uptflxsrc_b__[71], &tradsrc_b__[18193]
		    , &ddntflxdx[351], &duptflxdx[351], &dtraddx[89873]);

	}

/* L2601: */
    }

    return 0;
} /* sm_eq_trn__ */

