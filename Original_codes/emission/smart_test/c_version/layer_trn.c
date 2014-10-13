/* layer_trn.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int layer_trn__(logical *usrang0, integer *l0, integer *ng0, 
	integer *nlyr0, integer *nstr, integer *numu0, integer *nmom, integer 
	*nphi, integer *n_rad__, real *dtau_b__, real *copi0_b__, real *
	pmom_b__, real *dx_b_i__, real *umu, real *phi, real *umu_f__, real *
	gwt_f__, doublereal *trnflx_b__, doublereal *dtrnflxdx, doublereal *
	refflx_b__, doublereal *drefflxdx, doublereal *absflx_b__, doublereal 
	*dabsflxdx, doublereal *trnrad_b__, doublereal *dtrnraddx, doublereal 
	*refrad_b__, doublereal *drefraddx, doublereal *absrad_b__, 
	doublereal *dabsraddx)
{
    /* Initialized data */

    static integer ifirst = 0;

    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer k, l, m;
    static real uu[17920]	/* was [16][70][16] */, alb;
    static integer ibn, mom, naz, nze;
    static real phi0, phi1[16], umu0, umu1[16], dfdt[70];
    static integer iref;
    static doublereal dtau;
    static real uavg[70];
    static integer ntau;
    static real pmom[14070]	/* was [201][70] */, flup[70], utau[70];
    static logical prnt[7];
    static integer nlyr, numu;
    static real wvnm, fbeam;
    static integer ibcnd;
    static real accur, dtauc[70], ssalb[70];
    static logical plank;
    static real btemp, rfldn[70], temis, fisot, ttemp;
    static integer numuf;
    static real albmed[16];
    static logical lamber;
    static real rfldir[70];
    static integer maxphi;
    static real trnmed[16], temper[71];
    static integer maxcly;
    static logical usrang;
    static integer maxmom;
    static logical onlyfl;
    extern /* Subroutine */ int disort_(integer *, real *, real *, integer *, 
	    real *, real *, real *, logical *, integer *, real *, integer *, 
	    integer *, logical *, integer *, real *, integer *, real *, 
	    integer *, real *, real *, real *, real *, logical *, integer *, 
	    real *, real *, real *, real *, real *, logical *, logical *, 
	    real *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *);
    static integer maxulv, maxumu, iunits;
    static logical usrtau;
    static real surf_pr__[4];


/* ccccccccccccccccccccccccc  l a y e r _ t r n   cccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes the effective layer radiance and flux  cc */
/* c    transmittances, reflectances and and absorptances for each      cc */
/* c    spectral bin.  These quantities are used in the spectral        cc */
/* c    mapping and jacobian calculations.                              cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       nlyr0: number of layers on the atmosphere.                   cc */
/* c        nstr: number of computational radiance streams.             cc */
/* c         ng0: spectral bin number                                   cc */
/* c      dtau_b: layer optical depth for bin, ibn0.                    cc */
/* c     copi0_b: layer single scattering albedo for bin, ibn0.         cc */
/* c      pmom_b: layer particle phase function for bin, ibn0.          cc */
/* c         umu: zenith angle of each stream.                          cc */
/* c         gwt: gaussian weight for each stream.                      cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c   trnrad_b: layer transmittance for each radiance stream.          cc */
/* c   refrad_b: layer reflectance for each radiance stream.            cc */
/* c   absrad_b: layer absorptance for each rediance stream             cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  l a y e r _ t r n   cccccccccccccccccccccccc */




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









/* ****    local variables used in disort */




/* *****   state vector variables. */


/* ***    disort output variables */


/* ****   binned layer flux transmittances and absorptances */


/* ****   binned layer radiance transmittances and absorptances */



    /* Parameter adjustments */
    dabsraddx -= 4497;
    absrad_b__ -= 6737;
    drefraddx -= 4497;
    refrad_b__ -= 6737;
    dtrnraddx -= 4497;
    trnrad_b__ -= 6737;
    dabsflxdx -= 281;
    absflx_b__ -= 421;
    drefflxdx -= 281;
    refflx_b__ -= 421;
    dtrnflxdx -= 281;
    trnflx_b__ -= 421;
    --gwt_f__;
    --umu_f__;
    --phi;
    --umu;
    dx_b_i__ -= 421;
    pmom_b__ -= 84621;
    copi0_b__ -= 421;
    dtau_b__ -= 421;
    n_rad__ -= 6;

    /* Function Body */

    if (ifirst == 0) {

/* ****      initialzie disort varibles */

	ifirst = 1;
	numu = *numu0;
	iunits = 1;

/* ****      set logicals */

	usrang = FALSE_;
	onlyfl = FALSE_;
	plank = FALSE_;
	lamber = TRUE_;
	usrtau = FALSE_;
	iref = 0;
	wvnm = 1.f;

/* ****    set variable for isotropic illumniation from top */

	ibcnd = 1;

/* ****         set dimensions of arrays */

	maxcly = 2;
	maxulv = 2;
	maxumu = 16;
	maxphi = 16;
	maxmom = 200;

/* ****     turn off all discr_ord model internal print flags */

	for (m = 1; m <= 7; ++m) {
	    prnt[m - 1] = FALSE_;
/* L1001: */
	}

/* ****     set surface albedos to zero */

	alb = 0.f;
	for (m = 1; m <= 4; ++m) {
	    surf_pr__[m - 1] = 0.f;
/* L1021: */
	}

	ntau = 0;
	accur = 1e-4f;

/* ****     set number of layers to 1 */

	nlyr = 1;

    }

    ibn = *ng0;
    l = *l0;

    if (l < 5 && n_rad__[l + ibn * 5] != 0) {

/* ****    find the layer transittances and reflectances for each layer */

/*      write(*,'(/,1a)') 'above disort layer loop: ' */
/*      write(*,'(1a,16(1pe12.4))') 'umu: ',(umu(nze),nze=1,nstr) */
/*      write(*,'(1a,16(1pe12.4))') 'umu_f:',(umu_f(nze),nze=1,nstr) */
	i__1 = *nlyr0;
	for (k = 1; k <= i__1; ++k) {

/* ****       find layer transittances and reflectances for each sounding */

	    dtauc[0] = dtau_b__[k + (l + ibn * 5) * 70];
	    ssalb[0] = 1.f - copi0_b__[k + (l + ibn * 5) * 70];
	    if (ssalb[0] < 0.f) {
		ssalb[0] = 0.f;
	    }
	    pmom[0] = 1.f;
	    i__2 = *nmom;
	    for (mom = 1; mom <= i__2; ++mom) {
		pmom[mom] = pmom_b__[mom + (k + (l + ibn * 5) * 70) * 201];
/* L2001: */
	    }

	    numuf = *nstr;
	    usrang = FALSE_;

/* ****         call disort to find layer transmissivities and albedos. */
/*             note: in this mode, it changes the values of umu, so a */
/*             dummy variable should be used. */

	    disort_(&nlyr, dtauc, ssalb, nmom, pmom, temper, &wvnm, &usrtau, &
		    ntau, utau, nstr, &iunits, &usrang, &numuf, umu1, nphi, 
		    phi1, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    iref, surf_pr__, &alb, &btemp, &ttemp, &temis, &plank, &
		    onlyfl, &accur, prnt, &maxcly, &maxulv, &maxumu, &maxphi, 
		    &maxmom, rfldir, rfldn, flup, dfdt, uavg, uu, albmed, 
		    trnmed);

/* ****           repack radiance variables */

	    i__2 = *nstr / 2;
	    for (nze = 1; nze <= i__2; ++nze) {
		trnmed[nze + *nstr / 2 - 1] = trnmed[nze - 1];
		albmed[nze + *nstr / 2 - 1] = albmed[nze - 1];
/* L2201: */
	    }
	    i__2 = *nstr / 2;
	    for (nze = 1; nze <= i__2; ++nze) {
		trnmed[nze - 1] = trnmed[*nstr - nze];
		albmed[nze - 1] = albmed[*nstr - nze];
/* L2221: */
	    }

/* ****         find the flux transmittance and absorptance by */
/*             integrating transmittance over zenith angle using */
/*             gaussian quadrature */

	    trnflx_b__[k + (l + ibn * 5) * 70] = 0.f;
	    refflx_b__[k + (l + ibn * 5) * 70] = 0.f;
	    absflx_b__[k + (l + ibn * 5) * 70] = 0.f;
	    i__2 = *nstr / 2;
	    for (nze = 1; nze <= i__2; ++nze) {
		trnflx_b__[k + (l + ibn * 5) * 70] += gwt_f__[nze] * 2.f * 
			umu_f__[nze] * trnmed[nze - 1];
		refflx_b__[k + (l + ibn * 5) * 70] += gwt_f__[nze] * 2.f * 
			umu_f__[nze] * albmed[nze - 1];
		absflx_b__[k + (l + ibn * 5) * 70] += gwt_f__[nze] * 2.f * 
			umu_f__[nze] * (1.f - (trnmed[nze - 1] + albmed[nze - 
			1]));
/* L2241: */
	    }

	    if (*usrang0) {

/* ****         find layer transmissivities and reflectivities for */
/*             the user-specified zenith angles */

		numu = *numu0;
		usrang = TRUE_;
		i__2 = numu;
		for (nze = 1; nze <= i__2; ++nze) {
		    umu1[nze - 1] = umu[nze];
/* L2301: */
		}
		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {
		    phi1[naz - 1] = phi[nze];
/* L2321: */
		}

/* ****         call disort to find layer transmissivities and albedos. */
/*             note: in this mode, it changes the values of umu, so a */
/*             dummy variable should be used. */

		disort_(&nlyr, dtauc, ssalb, nmom, pmom, temper, &wvnm, &
			usrtau, &ntau, utau, nstr, &iunits, &usrang, &numu, 
			umu1, nphi, phi1, &ibcnd, &fbeam, &umu0, &phi0, &
			fisot, &lamber, &iref, surf_pr__, &alb, &btemp, &
			ttemp, &temis, &plank, &onlyfl, &accur, prnt, &maxcly,
			 &maxulv, &maxumu, &maxphi, &maxmom, rfldir, rfldn, 
			flup, dfdt, uavg, uu, albmed, trnmed);

	    }

/* ****         load layer radiance tranmittances */

	    i__2 = *numu0;
	    for (nze = 1; nze <= i__2; ++nze) {
		dtau = dtau_b__[k + (l + ibn * 5) * 70] / (r__1 = umu[nze], 
			dabs(r__1));
		trnrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] = exp(-dtau);
		refrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] = albmed[nze 
			- 1];
		absrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] = 1. - (
			trnmed[nze - 1] + albmed[nze - 1]);
/* L2401: */
	    }
/* L2461: */
	}
/*      write(*,'(/,1a)') 'below disort layer loop: ' */
/*      write(*,'(1a,16(1pe12.4))') 'umu: ',(umu(nze),nze=1,nstr) */
/*      write(*,'(1a,16(1pe12.4))') 'umu_f:',(umu_f(nze),nze=1,nstr) */
/*      do 1 nze=1,numu */
/*       write(*,'(/,1a,2i5,1pe12.4)') */
/*     -      'layer_trn: l,nze,umu(nze)',l,nze,umu(nze) */
/* 1       write(*,'(1a,16(1pe12.4))') 'trnrad_b', */
/*     - (trnrad_b(nze,k,l,ibn),k=1,nlyr0) */

    } else {

/* ****      either l=5 (perturbed surface albedos) or nrad=0. */
/*          in either case, the layer transmittance, reflectance, */
/*          and absorptance are set to values for nominal case */

	i__1 = *nlyr0;
	for (k = 1; k <= i__1; ++k) {
	    trnflx_b__[k + (ibn * 5 + 5) * 70] = trnflx_b__[k + (ibn * 5 + 1) 
		    * 70];
	    refflx_b__[k + (ibn * 5 + 5) * 70] = refflx_b__[k + (ibn * 5 + 1) 
		    * 70];
	    absflx_b__[k + (ibn * 5 + 5) * 70] = absflx_b__[k + (ibn * 5 + 1) 
		    * 70];

	    i__2 = numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		trnrad_b__[nze + (k + (ibn * 5 + 5) * 70 << 4)] = trnrad_b__[
			nze + (k + (ibn * 5 + 1) * 70 << 4)];
		refrad_b__[nze + (k + (ibn * 5 + 5) * 70 << 4)] = refrad_b__[
			nze + (k + (ibn * 5 + 1) * 70 << 4)];
		absrad_b__[nze + (k + (ibn * 5 + 5) * 70 << 4)] = absrad_b__[
			nze + (k + (ibn * 5 + 1) * 70 << 4)];
/* L2601: */
	    }
/* L2621: */
	}

    }

    if (l > 1 && l < 5) {

/* ****        find the transmittance and absorptance jacobians. */
/*            these are only needed for optical depth (l=2), */
/*            single scattering albedo (l=3), and scattering */
/*            phase function (l=4) */

	i__1 = *nlyr0;
	for (k = 1; k <= i__1; ++k) {
	    dtrnflxdx[k + (l - 1 + ibn * 3) * 70] = (trnflx_b__[k + (l + ibn *
		     5) * 70] - trnflx_b__[k + (ibn * 5 + 1) * 70]) * 
		    dx_b_i__[k + (l + ibn * 5) * 70];
	    drefflxdx[k + (l - 1 + ibn * 3) * 70] = (refflx_b__[k + (l + ibn *
		     5) * 70] - refflx_b__[k + (ibn * 5 + 1) * 70]) * 
		    dx_b_i__[k + (l + ibn * 5) * 70];
	    dabsflxdx[k + (l - 1 + ibn * 3) * 70] = (absflx_b__[k + (l + ibn *
		     5) * 70] - absflx_b__[k + (ibn * 5 + 1) * 70]) * 
		    dx_b_i__[k + (l + ibn * 5) * 70];
	    i__2 = numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		dtrnraddx[nze + (k + (l - 1 + ibn * 3) * 70 << 4)] = (
			trnrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] - 
			trnrad_b__[nze + (k + (ibn * 5 + 1) * 70 << 4)]) * 
			dx_b_i__[k + (l + ibn * 5) * 70];
		drefraddx[nze + (k + (l - 1 + ibn * 3) * 70 << 4)] = (
			refrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] - 
			refrad_b__[nze + (k + (ibn * 5 + 1) * 70 << 4)]) * 
			dx_b_i__[k + (l + ibn * 5) * 70];
		dabsraddx[nze + (k + (l - 1 + ibn * 3) * 70 << 4)] = (
			absrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] - 
			absrad_b__[nze + (k + (ibn * 5 + 1) * 70 << 4)]) * 
			dx_b_i__[k + (l + ibn * 5) * 70];
/* L2801: */
	    }

/* L2821: */
	}

    }

    return 0;
} /* layer_trn__ */

