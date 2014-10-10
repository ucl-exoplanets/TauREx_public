/* map_back.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int map_back__(logical *lsolar, logical *lplanck, logical *
	lamber, logical *usrang, integer *iugrp, integer *iuthrm, integer *
	iusol0, integer *iusol2, integer *iutrn, integer *iu_flux__, integer *
	iuout, integer *npd, integer *ipd, integer *iunits, integer *irad, 
	integer *iu_pd__, integer *iutpd, integer *iuspd, integer *iupdrad, 
	integer *ifrmout, integer *nza, integer *nz0, integer *nstr, integer *
	numu, integer *nzup, integer *nzdn, integer *nphi, integer *nlyr, 
	integer *nlay, integer *levout, integer *nlout, integer *k_out__, 
	integer *igs, integer *iref, integer *modepd, integer *nt_pd__, 
	integer *nstate, integer *istate, integer *ntau_pd__, integer *
	isptype, integer *islit, doublereal *width, doublereal *dwn, real *
	points, doublereal *wnmin, doublereal *wnmax, doublereal *wnsmin, 
	doublereal *wnsmax, doublereal *wn_tol__, real *units, real *umu0nz, 
	real *phi0nz, real *umu_f__, real *gwt_f__, real *umu, real *phi, 
	real *pd_pert__, real *p, real *t, real *ts, real *dp_dp__, real *
	rmix, real *dtauaer, real *alb_b__, real *dtau_b__, real *copi0_b__, 
	real *g_b__, doublereal *trnflx_b__, doublereal *dtrnflxdx, 
	doublereal *refflx_b__, doublereal *drefflxdx, doublereal *absflx_b__,
	 doublereal *dabsflxdx, doublereal *refrad_b__, doublereal *drefraddx,
	 doublereal *absrad_b__, doublereal *dabsraddx, doublereal *brdf_b__, 
	doublereal *dbrdfdx, doublereal *dnsflxsrc_b__, doublereal *ddnsflxdx,
	 doublereal *upsflxsrc_b__, doublereal *dupsflxdx, doublereal *
	dntflxsrc_b__, doublereal *ddntflxdx, doublereal *uptflxsrc_b__, 
	doublereal *duptflxdx, doublereal *sradsrc_b__, doublereal *dsraddx, 
	doublereal *tradsrc_b__, doublereal *dtraddx, doublereal *dirsoflx, 
	doublereal *dnsoflx, doublereal *upsoflx, doublereal *dnthflx, 
	doublereal *upthflx)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer f_rew(alist *), f_clos(cllist *), f_open(olist *), s_wsue(cilist *
	    ), do_uio(integer *, char *, ftnlen), e_wsue(void), s_rsue(cilist 
	    *), e_rsue(void), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real dn_s_flx__[1400]	/* was [70][10][2] */, dn_t_flx__[140]
	    	/* was [70][2] */, up_s_flx__[1400]	/* was [70][10][2] */,
	     up_t_flx__[140]	/* was [70][2] */;
    extern /* Subroutine */ int flux_int__(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, real *, real *, real *, real *
	    , real *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static real tau_aer_0__, tau_gas_0__, tau_ray_0__;
    static integer k, n;
    static real pd_absflx__[700]	/* was [70][10] */, dir_s_flx__[1400]	
	    /* was [70][10][2] */, pd_refflx__[700]	/* was [70][10] */;
    static integer n1;
    static real pd_trndir__[700]	/* was [70][10] */, pd_trnflx__[700]	
	    /* was [70][10] */;
    static integer ne;
    static doublereal wn;
    static integer nz;
    static real pd_dns_src__[700]	/* was [70][10] */, pd_dnt_src__[700]	
	    /* was [70][10] */;
    static integer ne0;
    static real pd_ups_src__[700]	/* was [70][10] */, pd_upt_src__[700]	
	    /* was [70][10] */;
    extern /* Subroutine */ int find_units__(integer *, integer *, doublereal 
	    *, real *);
    static doublereal wn0, wn1;
    static real alb[10], rad[768]	/* was [16][16][3] */;
    static integer ifl;
    static doublereal dvi;
    static integer npd1;
    static real sol0[2];
    static doublereal wns1;
    static integer iend[113], ibin, iuin;
    static doublereal dist[113];
    static integer nout;
    static doublereal dist0, delnu, wn_io__;
    extern /* Subroutine */ int solar_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, real *, real 
	    *);
    static doublereal wnvar[2];
    static real p_aer0__[16], p_gas0__[16];
    static doublereal wnout[2260]	/* was [2][113][10] */, wntst;
    static real p_ray0__[16], pd_rad__[537600]	/* was [16][16][3][70][10] */,
	     dn_flx__[3];
    extern /* Subroutine */ int pd_out__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, logical *, 
	    logical *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *);
    static real dsoldv, up_flx__[3];
    static integer ninter, iflout[1130]	/* was [113][10] */, nsivar[2];
    static real p_aer_0__, solflx, p_gas_0__;
    static integer nsiout[2260]	/* was [2][113][10] */;
    static real p_ray_0__;
    static integer iend_ne__;
    extern /* Subroutine */ int map_rad__(logical *, logical *, logical *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, real *, real *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *);
    static real abs_flx__[70], ref_flx__[70], dir_flx__[3], dns_src__[70];
    static doublereal dwn_max__;
    static real dnt_src__[70];
    extern /* Subroutine */ int rad_out__(logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *);
    static real trn_dir__[70];
    static doublereal distmin;
    static real tau_ext__[2100]	/* was [70][3][10] */, trn_flx__[70], 
	    ups_src__[70], upt_src__[70];
    static doublereal wn_next__;

    /* Fortran I/O blocks */
    static cilist io___29 = { 0, 0, 0, 0, 0 };
    static cilist io___30 = { 0, 0, 1, 0, 0 };
    static cilist io___32 = { 0, 6, 0, "(1a,i5,8(1pe14.7))", 0 };
    static cilist io___33 = { 0, 6, 0, "(/,1a,/,2a,2(1pe14.6))", 0 };
    static cilist io___34 = { 0, 6, 0, "(/,7a)", 0 };
    static cilist io___35 = { 0, 6, 0, "(/,6i5,6(1pe14.6))", 0 };
    static cilist io___74 = { 0, 6, 0, "(1a,8(1pe14.7))", 0 };



/* ccccccccccccccccccccccccc   m a p _ b a c k   ccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine maps binned radiances back to a high-           cc */
/* c    resolution spectral grid.                                       cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c      iugrp - unit number of file with binned radiances             cc */
/* c      iutrn - unit number for output transmission file              cc */
/* c    iu_flux - unit number for output level-dependent fluxes         cc */
/* c     iusol0 - unit number for solar radiance scratch file           cc */
/* c     iusol2 - unit number with solar radiances                      cc */
/* c     iuthrm - unit number for thermal radiances                     cc */
/* c      iuout - unit number of output radiance file                   cc */
/* c      iutpd - unit numbers for output level-dependent               cc */
/* c              thermal fluxes and thier partial derivatives          cc */
/* c      iuspd - unit numbers for output level-dependent               cc */
/* c              solar fluxes and thier partial derivatives            cc */
/* c    iupdrad - unit numbers for output level-dependent radiances     cc */
/* c              and thier partial derivatives                         cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c     iunits - index of output radiance units:                       cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) ergs/s/cm**2/sr/cm-1                               cc */
/* c              5) photons/s/m**2/sr/micron                           cc */
/* c        nz0 - index of this solar zenith angle (1 to nsol)          cc */
/* c        nza - number of solar zenith angles                         cc */
/* c    ifrmout - output file format (1) ascii or (2) binary            cc */
/* c     levout - output level index (1) top of atmosphere,             cc */
/* c              (2) surface, (3) arbitrary level                      cc */
/* c      dp_dp - vertical pressure gradient (p-pout)/dp                cc */
/* c      k_out - index of output level (1 - nlev)                      cc */
/* c       nlyr - number of computational model layers                  cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       numu - number of zenith angles in input file                 cc */
/* c    isptype - output spectrum type:                                 cc */
/* c              1) Full-resolution spectrum                           cc */
/* c              2) Sampled spectrum smoothed with a slit function     cc */
/* c       irad - index of output file type:                            cc */
/* c              1) fluxes, radiances, and heating rates,              cc */
/* c                 at computational azimuths and zenith angles,       cc */
/* c              2) fluxes, radiances, heating rates, and transmission cc */
/* c                functions at computational zenith angles,           cc */
/* c              3) fluxes, radiances, heating rates, and contribution cc */
/* c                 functions at computational zenith angles,          cc */
/* c              4) fluxes, radiances, heating rates, transmission     cc */
/* c                 functions and and contribution functions           cc */
/* c                 at computational zenith angles,                    cc */
/* c              5) fluxes, radiances, and heating rates,              cc */
/* c                 at computational azimuths and zenith angles,       cc */
/* c              6) fluxes, radiances and transmission functions       cc */
/* c                 at arbitrary zenith angles,                        cc */
/* c              7) fluxes, radiances, and contribution functions      cc */
/* c                 at arbitrary zenith angles,                        cc */
/* c              8) fluxes, radiances, transmission functions,and      cc */
/* c                 contribution functions at arbitrary zenith angles. cc */
/* c     nstate - number of variable elements in the state vecror       cc */
/* c     istate - state vector flag indicating which state variables    cc */
/* c              are variable components of the state vector.          cc */
/* c              1 - surface pressure                                  cc */
/* c              2 - surface/atmospheric temperature                   cc */
/* c              3 - gas absorption coeffient                          cc */
/* c              4 - cloud/aerosol optical depth                       cc */
/* c              5 - surface albedo                                    cc */
/* c      islit - index of slit function type:                          cc */
/* c              1) boxcar                                             cc */
/* c              2) triangular (approximate slit spectrometer)         cc */
/* c      width - half-width at half max of slit function  (cm**-1)     cc */
/* c        dwn - output sampling resolution (cm**-1)                   cc */
/* c         wn - current wavenumber (cm**-1)                           cc */
/* c     solflx - solar flux at this wavenumber                         cc */
/* c      wnmin - minimum wavenumber in output spectral grid (cm**-1)   cc */
/* c      wnmax - maximum wavenumber in output spectral grid (cm**-1)   cc */
/* c     wnsmin - minimum wavenumber in this binned segment (cm**-1)    cc */
/* c     wnsmax - maximum wavenumber in this binned segment (cm**-1)    cc */
/* c         ts - surface temperature (K)                               cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c         wl - current wavenumber (cm**-1)                           cc */
/* c         wn - current wavenumber (cm**-1)                           cc */
/* c    dir_flx - downward direct flux (irradiance) at wavenumber wn    cc */
/* c     dn_flx - total downward flux (irradiance) at wavenumber wn     cc */
/* c     up_flx - upward flux (irradiance) at wavenumber wn             cc */
/* c        rad - radiance at each output zenith and azimuth angle      cc */
/* c  trn_ray_0 - normal-incidence rayleigh-scattering transmission     cc */
/* c  trn_gas_0 - normal-incidence gas transmission                     cc */
/* c  trn_ray_0 - normal-incidence aerosol transmission                 cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc   m a p _ b a c k   ccccccccccccccccccccccccc */




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




/* ****   spectrally-dependent output flux and radiance */


/* ****    number and type of partial derivatives */


/* ****   counters used in map_back */


/* ****    units for flux and radiance jacobians */


/* ****   spectral binning parameters */


/* ****   spectral binning parameters */


/* ****   binned layer radiance transmittances and absorptances */


/* ****   binned layer flux transmittances and absorptances */


/* ****    aerosol optical depth for each particle mode */


/* ****   binned solar source function variables and jacobians */



/* ****   binned thermal source function variables and jacobians */






/* ****     monochormatic radiance interpolation variables */


/* ****   solar variables */



/* ****    output fluxes and radiances */





/* ****    flux transmission and absorption functions at wavenuber, wn */


/* ****    flux transmission and absorption partial */
/*        derivatives for simplified adding method at wavenumber wn */


/* ****   perturbations used for partial derivatives */


/* ****   layer-dependent source terms for solar and thermal fluxes */


/* ****   output flux and radiance partial derivatives */


/* ****   slit function variables */


/* ****   double precision variables for flux integration */


/* ****   specify the solar zenith angle index */

    /* Parameter adjustments */
    --upthflx;
    --dnthflx;
    upsoflx -= 71;
    dnsoflx -= 71;
    dirsoflx -= 71;
    dtraddx -= 89873;
    tradsrc_b__ -= 18193;
    dsraddx -= 89873;
    sradsrc_b__ -= 18193;
    duptflxdx -= 351;
    uptflxsrc_b__ -= 71;
    ddntflxdx -= 351;
    dntflxsrc_b__ -= 71;
    dupsflxdx -= 351;
    upsflxsrc_b__ -= 71;
    ddnsflxdx -= 351;
    dnsflxsrc_b__ -= 71;
    dbrdfdx -= 273;
    brdf_b__ -= 1553;
    dabsraddx -= 4497;
    absrad_b__ -= 6737;
    drefraddx -= 4497;
    refrad_b__ -= 6737;
    dabsflxdx -= 281;
    absflx_b__ -= 421;
    drefflxdx -= 281;
    refflx_b__ -= 421;
    dtrnflxdx -= 281;
    trnflx_b__ -= 421;
    g_b__ -= 421;
    copi0_b__ -= 421;
    dtau_b__ -= 421;
    alb_b__ -= 61;
    dtauaer -= 11;
    rmix -= 71;
    --dp_dp__;
    --t;
    --p;
    --pd_pert__;
    --phi;
    --umu;
    --gwt_f__;
    --umu_f__;
    --points;
    --istate;
    --modepd;
    --igs;
    --k_out__;
    iupdrad -= 10513;
    iuspd -= 11;
    iutpd -= 11;
    --ipd;

    /* Function Body */
    nz = *nz0;

/* **** rewind thermal flux/radiance file if necessary */

    if (*lplanck) {
	al__1.aerr = 0;
	al__1.aunit = *iuthrm;
	f_rew(&al__1);
    }

/* ****    set a maximum spectral interval step size */

    dwn_max__ = (*wnmax - *wnmin) * .01;
    if (dwn_max__ < *wnmin * .01) {
	dwn_max__ = *wnmin * .01;
    }
/* *****   initialize the spectral grid flags for the binned */
/*        radiances and the solar fluxes */

    if (*lsolar) {

/*         Two spectral grids (solar fluxes and spectral radiance) */
/*         must be combined, so the number of outoput grirds, nout = 2. */

	nout = 2;
	if (*wnsmin == *wnmin) {

/* ****      set the value of the variable n1 = 1 so that the */
/*          index counters, nsiout, for both the solar fluxes */
/*          and the binned radiances are reset */

	    n1 = 1;

/* ****      since this is the first spectral interval, set the */
/*          "last solar flux" quantity, wns1 = wnsmin */

	    wns1 = *wnsmin;
	} else {

/* ****      Set the the variable n1 = 2 such that the index */
/*          counters for the binned radiances are reset, but */
/*          the index counters for the solar flux are not reset */

	    n1 = 2;

/* ****       test to see if the last solar flux value is beyond */
/*           the current wavenumber, wn.  If it is not, set the */
/*           mininum wavenumber value for this interval to a */
/*           value that is just beyond the last solar flux value. */

	    if (nz == 1) {
		wn_next__ = (*wn_tol__ + 1.) * wns1;
		if (wn_next__ < *wnsmin) {
		    *wnsmin = wn_next__;
		}
	    }
	}

/* ****     wavnumber dependent solar fluxes are written to */
/*         unit iusol0 as the program processes the first solar */
/*         zenith angle.  These values are used for all */
/*         other solar zenith angles. If this is the first */
/*         solar zenith angle, close and reopen the solar */
/*         flux unit. if multiple solar zenith angles are */
/*         used rewind this unit */

	if (*nza > 1) {
	    if (nz == 1) {

		cl__1.cerr = 0;
		cl__1.cunit = *iusol0;
		cl__1.csta = 0;
		f_clos(&cl__1);
		o__1.oerr = 0;
		o__1.ounit = *iusol0;
		o__1.ofnm = 0;
		o__1.orl = 0;
		o__1.osta = "scratch";
		o__1.oacc = 0;
		o__1.ofm = "unformatted";
		o__1.oblnk = 0;
		f_open(&o__1);

	    } else {
		al__1.aerr = 0;
		al__1.aunit = *iusol0;
		f_rew(&al__1);
	    }
	}

    } else {

/* ****     solar fluxes are not included.  Thermal fluxes are the only */
/*         output quantity, so nout = 1, and n1 = 1 */

	nout = 1;
	n1 = 1;

    }

/* ****  initialize wavelength pointers */

    wn = *wnsmin;
    wn0 = *wnsmin;
    wn1 = *wnsmin;
    wn_io__ = -999.;

    if (*wnsmin == *wnmin) {

/* ****     intialize output flags and spectral grid pointers. */

	i__1 = nout;
	for (ne = n1; ne <= i__1; ++ne) {

/* ****      set output flag */

	    iflout[ne + nz * 113 - 114] = 1;

/* ****         set end of bin file flag */

	    iend[ne - 1] = 0;

/* ****      initialize spectral counters for radiances and solar flux. */

	    nsiout[(ne + nz * 113 << 1) - 228] = 0;
	    nsiout[(ne + nz * 113 << 1) - 227] = 0;
/* L1001: */
	}

/* ****     set partial derivative pointers */

	i__1 = *nstate;
	for (n = 1; n <= i__1; ++n) {
	    if (istate[n] != 0) {

		for (k = 1; k <= 70; ++k) {
		    pd_ups_src__[k + n * 70 - 71] = 0.f;
		    pd_dns_src__[k + n * 70 - 71] = 0.f;
		    pd_upt_src__[k + n * 70 - 71] = 0.f;
		    pd_dnt_src__[k + n * 70 - 71] = 0.f;
/* L1201: */
		}

	    }
/* L1221: */
	}

/* ****     set the index variable for optical depths: if ipd = 1, */
/*         the first varible component of the state vector is pressure. */
/*         if ipd = 2, the second variable component of the state vector */
/*         is temperature. */

	npd1 = 0;
	if (abs(ipd[1]) > 0 && abs(ipd[1]) < 3) {
	    npd1 = 1;
	}
	if (abs(ipd[2]) > 0 && abs(ipd[2]) < 3) {
	    ++npd1;
	}

    }

    ninter = 0;

/* ****        s p e c t r a l    i n t e r v a l    l o o p */

/* ****   initialize spectral interval counters */

L2001:
    ne = 0;

    ++ninter;

/* ****       f i n d    t o a   s o l a r   f l u x e s */

    if (*lsolar) {

/* ****          increment spectral counter and check wn grid flag */

	++ne;

/* ****          if this is the first solar zenith angle, read the */
/*              toa solar flux values. */

	if (nz == 1) {

	    if (iflout[ne + nz * 113 - 114] == 1) {
		nsivar[0] = nsiout[(ne + nz * 113 << 1) - 228];
		nsivar[1] = nsiout[(ne + nz * 113 << 1) - 227];
		wnvar[0] = wnout[(ne + nz * 113 << 1) - 228];
		wnvar[1] = wnout[(ne + nz * 113 << 1) - 227];

/* ****            read solar fluxes at next wavenumber */

		solar_(iusol2, &ne, nsivar, iend, wnmax, &wn, wnvar, sol0, &
			dsoldv);

		nsiout[(ne + nz * 113 << 1) - 228] = nsivar[0];
		nsiout[(ne + nz * 113 << 1) - 227] = nsivar[1];
		wnout[(ne + nz * 113 << 1) - 228] = wnvar[0];
		wnout[(ne + nz * 113 << 1) - 227] = wnvar[1];

/* ****             turn off the wavenumber grid flag. */

		iflout[ne + nz * 113 - 114] = 0;

	    }

/* ****           interpolate toa solar fluxes to this wn */

	    solflx = sol0[0] + dsoldv * (real) (wn - wnout[(ne + nz * 113 << 
		    1) - 228]);
	    if (solflx < 0.f) {
		solflx = 0.f;
	    }

/* ****           if there is more than one solar zenith angle,save */
/*              toa solar fluxes and spectral derivative test values */

	    if (*nza > 1) {
		io___29.ciunit = *iusol0;
		s_wsue(&io___29);
		do_uio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) - 228], (
			ftnlen)sizeof(integer));
		do_uio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) - 227], (
			ftnlen)sizeof(integer));
		do_uio(&c__1, (char *)&wnout[nsiout[(ne + nz * 113 << 1) - 
			228] + (ne + nz * 113 << 1) - 229], (ftnlen)sizeof(
			doublereal));
		do_uio(&c__1, (char *)&wnout[nsiout[(ne + nz * 113 << 1) - 
			227] + (ne + nz * 113 << 1) - 229], (ftnlen)sizeof(
			doublereal));
		do_uio(&c__1, (char *)&wn, (ftnlen)sizeof(doublereal));
		do_uio(&c__1, (char *)&solflx, (ftnlen)sizeof(real));
		e_wsue();
	    }

	} else {

/* ****           read previously mapped solar fluxes */

	    io___30.ciunit = *iusol0;
	    i__1 = s_rsue(&io___30);
	    if (i__1 != 0) {
		goto L2201;
	    }
	    i__1 = do_uio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) - 228], 
		    (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L2201;
	    }
	    i__1 = do_uio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) - 227], 
		    (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L2201;
	    }
	    i__1 = do_uio(&c__1, (char *)&wnout[nsiout[(ne + nz * 113 << 1) - 
		    228] + (ne + nz * 113 << 1) - 229], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L2201;
	    }
	    i__1 = do_uio(&c__1, (char *)&wnout[nsiout[(ne + nz * 113 << 1) - 
		    227] + (ne + nz * 113 << 1) - 229], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L2201;
	    }
	    i__1 = do_uio(&c__1, (char *)&wntst, (ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L2201;
	    }
	    i__1 = do_uio(&c__1, (char *)&solflx, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L2201;
	    }
	    i__1 = e_rsue();
	    if (i__1 != 0) {
		goto L2201;
	    }
	    goto L2221;

/* ****             print end of file error if necessary and stop */

L2201:
	    s_wsfe(&io___32);
	    do_fio(&c__1, "map_back: read past EOF for solar flux unit", (
		    ftnlen)43);
	    do_fio(&c__1, (char *)&(*iusol0), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&wn, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&wntst, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&wn0, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*wnsmax), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&wns1, (ftnlen)sizeof(doublereal));
	    e_wsfe();

	    s_stop("", (ftnlen)0);

/* ****            check wavelength grid alignment */

L2221:
	    if ((d__1 = wntst - wn, abs(d__1)) > *wn_tol__ * 2. * wn) {
		s_wsfe(&io___33);
		do_fio(&c__1, "Error in map_back:", (ftnlen)18);
		do_fio(&c__1, "wavelength grids not aligned for different ", (
			ftnlen)43);
		do_fio(&c__1, "solar zenith angles: wn, wntst: ", (ftnlen)32);
		do_fio(&c__1, (char *)&wn, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&wntst, (ftnlen)sizeof(doublereal));
		e_wsfe();
		s_wsfe(&io___34);
		do_fio(&c__1, "nz,ninter", (ftnlen)9);
		do_fio(&c__1, "nsiout(1,1,nz),nsiout(2,1,nz)", (ftnlen)29);
		do_fio(&c__1, "nsiout(1,2,nz),nsiout(2,2,nz)", (ftnlen)29);
		do_fio(&c__1, "wnout(nsiout(1,1,nz),1,nz)", (ftnlen)26);
		do_fio(&c__1, "wnout(nsiout(2,1,nz),1,nz)", (ftnlen)26);
		do_fio(&c__1, "wnout(nsiout(1,2,nz),2,nz)", (ftnlen)26);
		do_fio(&c__1, "wnout(nsiout(2,2,nz),2,nz)", (ftnlen)26);
		do_fio(&c__1, "wntst,solflx", (ftnlen)12);
		e_wsfe();
		s_wsfe(&io___35);
		do_fio(&c__1, (char *)&nz, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ninter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nsiout[(nz * 113 + 1 << 1) - 228], (
			ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nsiout[(nz * 113 + 1 << 1) - 227], (
			ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nsiout[(nz * 113 + 2 << 1) - 228], (
			ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nsiout[(nz * 113 + 2 << 1) - 227], (
			ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&wnout[nsiout[(nz * 113 + 1 << 1) - 228]
			 + (nz * 113 + 1 << 1) - 229], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&wnout[nsiout[(nz * 113 + 1 << 1) - 227]
			 + (nz * 113 + 1 << 1) - 229], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&wnout[nsiout[(nz * 113 + 2 << 1) - 228]
			 + (nz * 113 + 2 << 1) - 229], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&wnout[nsiout[(nz * 113 + 2 << 1) - 227]
			 + (nz * 113 + 2 << 1) - 229], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&wntst, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&solflx, (ftnlen)sizeof(real));
		e_wsfe();
		s_stop("", (ftnlen)0);

	    }
	}
    }

/* ****        m a p    s p e c t r a l    b i n s   t o   t h i s   wn */

    ++ne;
    ifl = iflout[ne + nz * 113 - 114];
    iend_ne__ = iend[ne - 1];

    if (wn_io__ < *wnsmax) {
	ne0 = ne;

/* ***            interpolate binned properties to wavenumber wn */

	map_rad__(lsolar, lplanck, lamber, usrang, iugrp, iuthrm, &ifl, &
		iend_ne__, &ne0, nz0, nphi, nstr, numu, nzup, nzdn, nza, nlyr,
		 nlay, nlout, levout, &k_out__[1], iref, &ibin, nsiout, &igs[
		1], npd, &npd1, &ipd[1], ntau_pd__, nt_pd__, &modepd[1], irad,
		 wnmin, wnmax, &wn, &wn_io__, wnout, wn_tol__, &dvi, umu0nz, 
		phi0nz, &solflx, &umu_f__[1], &gwt_f__[1], &umu[1], &phi[1], 
		ts, &t[1], &p[1], &dp_dp__[1], &rmix[71], alb, &dtauaer[11], 
		tau_ext__, &pd_pert__[1], &alb_b__[61], &dtau_b__[421], &
		copi0_b__[421], &g_b__[421], &p_ray_0__, &p_gas_0__, &
		p_aer_0__, p_ray0__, p_gas0__, p_aer0__, &tau_ray_0__, &
		tau_gas_0__, &tau_aer_0__, &trnflx_b__[421], &dtrnflxdx[281], 
		&refflx_b__[421], &drefflxdx[281], &absflx_b__[421], &
		dabsflxdx[281], &refrad_b__[6737], &drefraddx[4497], &
		absrad_b__[6737], &dabsraddx[4497], &brdf_b__[1553], &dbrdfdx[
		273], &dnsflxsrc_b__[71], &ddnsflxdx[351], &upsflxsrc_b__[71],
		 &dupsflxdx[351], &dntflxsrc_b__[71], &ddntflxdx[351], &
		uptflxsrc_b__[71], &duptflxdx[351], &tradsrc_b__[18193], &
		dtraddx[89873], &sradsrc_b__[18193], &dsraddx[89873], 
		up_s_flx__, dn_s_flx__, dir_s_flx__, up_t_flx__, dn_t_flx__, 
		dns_src__, ups_src__, pd_dns_src__, pd_ups_src__, dnt_src__, 
		upt_src__, pd_dnt_src__, pd_upt_src__, pd_rad__, pd_trndir__, 
		pd_trnflx__, pd_refflx__, pd_absflx__, up_flx__, dn_flx__, 
		dir_flx__, rad, trn_dir__, trn_flx__, ref_flx__, abs_flx__);

/* ****           turn off the wavenumber grid flag. */

	if (iend[ne - 1] == 0) {
	    iflout[ne + nz * 113 - 114] = 0;
	} else {

/* *****           end if iugrp file encountered */

	    s_wsfe(&io___74);
	    do_fio(&c__1, "End of bin file: wn,wn_io", (ftnlen)25);
	    do_fio(&c__1, (char *)&wn, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&wn_io__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*wnsmax), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    *wnsmax = wn_io__;
	    if (nz == *nza) {
		wns1 = wn;
	    }
	    return 0;

	}

    }

/* ****           e f f e c t i v e    i n t e r v a l    w i d t h */

/* *****       find the next wavenumber where monochromatic properties */
/*            are specified.  Check all input data sets and select the */
/*            wavenumber that is closest to the current wavenumber. */

/*            note: the minimum distance must exceed the round-off */
/*                tolerance of the computer, wn_tol */

    distmin = 1e30;
    dist0 = *wn_tol__ * wn;

    i__1 = nout;
    for (ne0 = 1; ne0 <= i__1; ++ne0) {
	dist[ne0 - 1] = wnout[nsiout[(ne0 + nz * 113 << 1) - 227] + (ne0 + nz 
		* 113 << 1) - 229] - wn;
	if (dist[ne0 - 1] < distmin) {
	    if (dist[ne0 - 1] > dist0) {
		distmin = dist[ne0 - 1];
	    } else {
		distmin = dist0;
	    }
	}
/* L6001: */
    }

    if (distmin > dwn_max__) {
	distmin = dwn_max__;
    }

/* ****        set the next wavenumber and the effective spectral */
/*            interval width, dnu.  The spectral interval is assumed */
/*            to extend between the current wavenumber, wn, and half-way */
/*            to the previous value, wn0, and the next value, wn1. */

    wn1 = wn + distmin;
    delnu = wn - wn0;

/* ****         add this contribution to spectrally-integrated fluxes */

    flux_int__(lsolar, lplanck, nlyr, &nz, &ninter, &delnu, dir_s_flx__, 
	    dn_s_flx__, up_s_flx__, dn_t_flx__, up_t_flx__, &dirsoflx[71], &
	    dnsoflx[71], &upsoflx[71], &dnthflx[1], &upthflx[1]);

/* ****          find the factor needed to convert fluxes and radiances */
/*              from w/m**2/cm-1 to the desired output units */

    iuin = 1;

    find_units__(&iuin, iunits, &wn, units);

/* ****           print these radiances */

    rad_out__(lsolar, iuout, iutrn, iu_flux__, ifrmout, irad, levout, nlout, 
	    nlyr, nz0, nza, nphi, numu, nzup, isptype, islit, &wn, wnmin, 
	    wnmax, width, dwn, units, &p[1], &solflx, up_s_flx__, dn_s_flx__, 
	    dir_s_flx__, up_t_flx__, dn_t_flx__, up_flx__, dn_flx__, 
	    dir_flx__, rad, p_ray0__, p_gas0__, p_aer0__, &tau_ray_0__, &
	    tau_gas_0__, &tau_aer_0__, &p_ray_0__, &p_gas_0__, &p_aer_0__, &
	    points[1]);

    if (*nstate > 0) {

	pd_out__(nz0, npd, &npd1, &ipd[1], &igs[1], ifrmout, levout, nlout, 
		nlyr, nphi, numu, nzup, iu_pd__, &iutpd[11], &iuspd[11], &
		modepd[1], &iupdrad[10513], iunits, isptype, islit, lsolar, 
		lplanck, &wn, wnmin, wnmax, width, dwn, &points[1], units, 
		umu0nz, ts, &t[1], &p[1], &rmix[71], alb, &dtauaer[11], &
		solflx, up_s_flx__, dn_s_flx__, dir_s_flx__, up_t_flx__, 
		dn_t_flx__, &pd_pert__[1], dns_src__, pd_dns_src__, ups_src__,
		 pd_ups_src__, dnt_src__, pd_dnt_src__, upt_src__, 
		pd_upt_src__, trn_dir__, trn_flx__, ref_flx__, abs_flx__, 
		pd_trndir__, pd_trnflx__, pd_refflx__, pd_absflx__, rad, 
		pd_rad__);

    }

/* *****      specify the next wavenumber grid point and reset */
/*           wavenumber grid flags. */

    wn0 = wn;
    wn = wn1;

    i__1 = nout;
    for (ne = 1; ne <= i__1; ++ne) {
	dist[ne - 1] = wnout[nsiout[(ne + nz * 113 << 1) - 227] + (ne + nz * 
		113 << 1) - 229] - wn;
	if (dist[ne - 1] <= 0.f && iend[ne - 1] == 0) {

/* ****              wn is beyond latest input point for this constituent. */
/*                  set flag to read next point. */

	    iflout[ne + nz * 113 - 114] = 1;
	} else {
	    iflout[ne + nz * 113 - 114] = 0;
	}
/* L7001: */
    }

/* ****         get monochromatic optical properties for next wavenumber. */

    if (wn0 < *wnsmax) {
	goto L2001;
    }

/* ****   if this is the last zenith angle, set the end-of interval wn. */

    if (nz == *nza) {
	wns1 = wn1;
    }

    return 0;
} /* map_back__ */

