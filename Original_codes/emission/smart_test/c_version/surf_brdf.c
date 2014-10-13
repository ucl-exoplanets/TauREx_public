/* surf_brdf.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int surf_brdf__(logical *usrang, logical *lamber, integer *
	l0, integer *ng0, integer *nz0, integer *nstr, integer *numu, integer 
	*nphi, integer *nmom, integer *n_rad__, integer *nref, real *umu, 
	real *phi, real *dalb_b_i__, real *alb_b__, real *sur_b__, real *ws, 
	real *phiw, real *umu0, real *phi0, doublereal *brdf_b__, doublereal *
	dbrdfdx)
{
    /* Initialized data */

    static integer ifirst = 0;

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double acos(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer l, m;
    static real pi;
    static integer nr, nz;
    static real uu[17920]	/* was [16][70][16] */, alb;
    static integer ibn, mom, naz, nze;
    static real dfdt[70];
    static integer iref;
    static real uavg[70];
    static integer ntau;
    static real pmom[14070]	/* was [201][70] */, flup[70], utau[70];
    static logical prnt[7];
    static integer nlyr;
    static real wvnm, fbeam;
    static integer ibcnd;
    static real accur, dtauc[70], ssalb[70];
    static logical plank;
    static real rfldn[70], btemp, temis, fisot, ttemp, albmed[16], rfldir[70];
    static integer maxphi;
    static real trnmed[16], temper[71];
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
    static integer iunits, maxulv, maxumu;
    static logical usrtau;
    static real surf_pr__[4];

    /* Fortran I/O blocks */
    static cilist io___47 = { 0, 6, 0, "(/,1a,3i5)", 0 };
    static cilist io___48 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___49 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___50 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___51 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___52 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };



/* ccccccccccccccccccccccccc  s u r f _ b r d f   cccccccccccccccccccccccc */
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
/* c    brdf_b: surface bi-directional reflection distribution function cc */
/* c   dbrdfdx: partial derivative of brdf_b with respect to albedo     cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  s u r f _ b r d f   cccccccccccccccccccccccc */




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


/* ****   binned layer radiance transmittances and absorptances */


    /* Parameter adjustments */
    dbrdfdx -= 273;
    brdf_b__ -= 1553;
    sur_b__ -= 25;
    alb_b__ -= 61;
    dalb_b_i__ -= 11;
    --phi;
    --umu;
    n_rad__ -= 6;

    /* Function Body */

    if (ifirst == 0) {

/* ****      initialzie disort varibles */

	ifirst = 1;
	iunits = 1;
	ibcnd = 0;

/* ****      set logicals */

	*usrang = FALSE_;
	onlyfl = FALSE_;
	plank = FALSE_;
	usrtau = FALSE_;
	iref = 0;
	wvnm = 1.f;
	pi = acos(-1.f);

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

	ntau = 0;
	accur = 1e-4f;

/* ****     set number of layers to 1 */

	nlyr = 1;

    }

    ibn = *ng0;
    l = *l0;
    nz = *nz0;

    if (l == 1 || l == 5 && n_rad__[l + ibn * 5] != 0) {

/* ****     set surface temperature and atmospheric temperature */

	btemp = 296.f;
	temper[0] = 296.f;
	temper[1] = 296.f;
	wvnm = 1e3f;

/* ****             set top boundary condition - no emission or reflection */

	fisot = 0.f;
	ttemp = 0.f;
	temis = 1.f;
	fbeam = 1.f;

	dtauc[0] = 0.f;
	ssalb[0] = 0.f;
	pmom[0] = 1.f;
	i__1 = *nmom;
	for (mom = 1; mom <= i__1; ++mom) {
	    pmom[mom] = 0.f;
/* L2001: */
	}

/* ****    define the surface optical properties */

	alb = alb_b__[nz + (l + ibn * 5) * 10];

/* ****    set the surface optical properties for disort */

	i__1 = *nref;
	for (nr = 1; nr <= i__1; ++nr) {
	    surf_pr__[nr - 1] = sur_b__[nr + (l + ibn * 5 << 2)];
/* L2021: */
	}

	if (iref == 4) {

/* ****      set wind speed and direction for Cox/Munk model */

	    surf_pr__[2] = *ws;
	    surf_pr__[3] = *phiw;
	}

	disort_(&nlyr, dtauc, ssalb, nmom, pmom, temper, &wvnm, &usrtau, &
		ntau, utau, nstr, &iunits, usrang, numu, &umu[1], nphi, &phi[
		1], &ibcnd, &fbeam, umu0, phi0, &fisot, lamber, &iref, 
		surf_pr__, &alb, &btemp, &ttemp, &temis, &plank, &onlyfl, &
		accur, prnt, &maxcly, &maxulv, &maxumu, &maxphi, &maxmom, 
		rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed);

	if (l == 1) {

/* ****         load brdf's: note-values are independent of atmospheric */
/*             optical depth, single scattering albedo, and particle */
/*             phase function, so values of l=1,4 are the same. */

	    for (m = 1; m <= 4; ++m) {
		i__1 = *nphi;
		for (naz = 1; naz <= i__1; ++naz) {
		    i__2 = *numu;
		    for (nze = 1; nze <= i__2; ++nze) {
			brdf_b__[nze + (naz + (m + ibn * 5 << 4) << 4)] = pi *
				 uu[nze + (naz * 70 + 2 << 4) - 1137] / 
				rfldir[1];
/* L2421: */
		    }
/* L2441: */
		}
/* L2461: */
	    }

	} else {

/* ****        find the brdf jacobians. */

	    i__1 = *nphi;
	    for (naz = 1; naz <= i__1; ++naz) {
		i__2 = *numu;
		for (nze = 1; nze <= i__2; ++nze) {
		    brdf_b__[nze + (naz + (l + ibn * 5 << 4) << 4)] = pi * uu[
			    nze + (naz * 70 + 2 << 4) - 1137] / rfldir[1];
		    dbrdfdx[nze + (naz + (ibn << 4) << 4)] = (brdf_b__[nze + (
			    naz + (l + ibn * 5 << 4) << 4)] - brdf_b__[nze + (
			    naz + (ibn * 5 + 1 << 4) << 4)]) * dalb_b_i__[nz 
			    + ibn * 10];
/* L2801: */
		}
/* L2821: */
	    }

	}

    }
    s_wsfe(&io___47);
    do_fio(&c__1, "surf_brdf: nz,l,ibn: ", (ftnlen)21);
    do_fio(&c__1, (char *)&nz, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&l, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ibn, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___48);
    do_fio(&c__1, "rfldir    ", (ftnlen)10);
    do_fio(&c__1, (char *)&rfldir[nlyr], (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___49);
    do_fio(&c__1, "rfldn     ", (ftnlen)10);
    do_fio(&c__1, (char *)&rfldn[nlyr], (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___50);
    do_fio(&c__1, "uu(1)     ", (ftnlen)10);
    i__1 = *numu;
    for (nze = 1; nze <= i__1; ++nze) {
	do_fio(&c__1, (char *)&uu[nze - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___51);
    do_fio(&c__1, "uu(2)     ", (ftnlen)10);
    i__1 = *numu;
    for (nze = 1; nze <= i__1; ++nze) {
	do_fio(&c__1, (char *)&uu[nze + 15], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___52);
    do_fio(&c__1, "brdf_b: ", (ftnlen)8);
    i__1 = *numu;
    for (nze = 1; nze <= i__1; ++nze) {
	do_fio(&c__1, (char *)&brdf_b__[nze + ((l + ibn * 5 << 4) + 1 << 4)], 
		(ftnlen)sizeof(doublereal));
    }
    e_wsfe();

    return 0;
} /* surf_brdf__ */

