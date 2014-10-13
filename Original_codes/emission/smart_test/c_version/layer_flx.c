/* layer_flx.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int layer_flx__(logical *usrang, logical *lamber, integer *
	l0, integer *nlyr0, integer *ng0, integer *nstr, integer *numu, 
	integer *nzdn, integer *nzup, integer *nmom, integer *nphi, integer *
	iref, real *umu, real *phi, real *dtau_b__, real *copi0_b__, real *
	pmom_b__, real *umu0, real *phi0, real *accur)
{
    /* Initialized data */

    static integer ifirst = 0;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer k, l, m;
    static real uu[17920]	/* was [16][70][16] */, alb, rad[17920]	/* 
	    was [16][16][70] */;
    static integer ibn;
    static real dir[70];
    static integer mom, naz, nze;
    static real dfdt[70], flxd[350]	/* was [70][5] */, uavg[70];
    static integer nlev, ntau;
    static real pmom[14070]	/* was [201][70] */, flup[70], utau[70], flxu[
	    350]	/* was [70][5] */;
    static logical prnt[7];
    static integer nlyr;
    static real wvnm, fbeam;
    static integer ibcnd;
    static real dtauc[70], ssalb[70];
    static logical plank;
    static real rfldn[70], btemp, temis, fisot, ttemp, albmed[16], deltau, 
	    rfldir[70];
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
    static integer maxulv, maxumu, iunits;
    static logical usrtau;
    static real deltaui[70];
    static doublereal dflxddx[70], dflxudx[70];
    static real surf_pr__[4];

    /* Fortran I/O blocks */
    static cilist io___48 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___49 = { 0, 6, 0, "(1a,16(1pe12.4),4(/16(1pe12.4)))", 0 }
	    ;
    static cilist io___50 = { 0, 6, 0, "(1a,16(1pe12.4),4(/16(1pe12.4)))", 0 }
	    ;
    static cilist io___51 = { 0, 6, 0, "(1a,16(1pe12.4),4(/16(1pe12.4)))", 0 }
	    ;
    static cilist io___56 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___57 = { 0, 6, 0, "(1a,16(1pe12.4),4(/16(1pe12.4)))", 0 }
	    ;
    static cilist io___58 = { 0, 6, 0, "(1a,16(1pe12.4),4(/16(1pe12.4)))", 0 }
	    ;
    static cilist io___59 = { 0, 6, 0, "(1a,16(1pe12.4),4(/16(1pe12.4)))", 0 }
	    ;



/* ccccccccccccccccccccccccc  l a y e r _ f l x   cccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes the radiances and fluxes layer by      cc */
/* c    layer to test the adding method.                                cc */
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
/* c    trn_rad: layer transmittance for each radiance stream.          cc */
/* c    abs_rad: layer absorptance for each rediance stream             cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  l a y e r _ f l x   cccccccccccccccccccccccc */




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











/* ***   define a dummy cosine of zenith angle variable for call */
/*      to DISORT for the flux transmission calculation. */

/* *****   state vector variables. */


/* ***    disort output variables */



    /* Parameter adjustments */
    pmom_b__ -= 84621;
    copi0_b__ -= 421;
    dtau_b__ -= 421;
    --phi;
    --umu;

    /* Function Body */

/*      write(*,*) 'layer_flx',l0,nlyr0,ng0, */
/*     -                     nstr,numu,nzdn,nzup,nmom,nphi,iref */
    if (ifirst == 0) {

/* ****      initialzie disort varibles */

	ifirst = 1;
	nlev = *nlyr0 + 1;

/* ****      set logicals */

	onlyfl = FALSE_;
	plank = FALSE_;
	*lamber = TRUE_;
	usrtau = FALSE_;
	wvnm = 1.f;
	fbeam = 1.f;
	fisot = 0.f;
	iunits = 0;

/* ****    set variable for isotropic illumniation from top */

	ibcnd = 0;

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

/* ****     set number of layers to 1 */

	nlyr = 1;

    }

    ibn = *ng0;
    l = *l0;

/* ****     set surface albedos to zero */

    alb = 0.f;
    surf_pr__[0] = alb;
    for (m = 2; m <= 4; ++m) {
	surf_pr__[m - 1] = 0.f;
/* L1021: */
    }

/* ****    find the layer transittances and reflectances for each layer */

    dir[0] = *umu0 * fbeam;
    flxd[l * 70 - 70] = 0.f;
    i__1 = *nlyr0;
    for (k = 1; k <= i__1; ++k) {

/* ****    find the layer transittances and reflectances for each sounding */

	dtauc[0] = dtau_b__[k + (l + ibn * 5) * 70];
	ssalb[0] = 1.f - copi0_b__[k + (l + ibn * 5) * 70];
	pmom[0] = 1.f;
	i__2 = *nmom;
	for (mom = 1; mom <= i__2; ++mom) {
	    pmom[mom] = pmom_b__[mom + (k + (l + ibn * 5) * 70) * 201];
/* L2001: */
	}

	fbeam = 1.f;

	disort_(&nlyr, dtauc, ssalb, nmom, pmom, temper, &wvnm, &usrtau, &
		ntau, utau, nstr, &iunits, usrang, numu, &umu[1], nphi, &phi[
		1], &ibcnd, &fbeam, umu0, phi0, &fisot, lamber, iref, 
		surf_pr__, &alb, &btemp, &ttemp, &temis, &plank, &onlyfl, 
		accur, prnt, &maxcly, &maxulv, &maxumu, &maxphi, &maxmom, 
		rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed);

	flxu[k + l * 70 - 71] = flup[0];
	flxd[k + 1 + l * 70 - 71] = rfldn[1];
	dir[k] = rfldir[1];
	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {
	    i__3 = *nzdn;
	    for (nze = 1; nze <= i__3; ++nze) {
		rad[nze + (naz + (k + 1 << 4) << 4) - 273] = uu[nze + (naz * 
			70 + 2 << 4) - 1137];
/* L2201: */
	    }
	    i__3 = *numu;
	    for (nze = *nzup; nze <= i__3; ++nze) {
		rad[nze + (naz + (k << 4) << 4) - 273] = uu[nze + (naz * 70 + 
			1 << 4) - 1137];
/* L2221: */
	    }
/* L2241: */
	}

/* L2441: */
    }
    k = nlev;
    i__1 = *nphi;
    for (naz = 1; naz <= i__1; ++naz) {
	i__2 = *numu;
	for (nze = *nzup; nze <= i__2; ++nze) {
	    rad[nze + (naz + (k << 4) << 4) - 273] = uu[nze + (naz * 70 + 1 <<
		     4) - 1137];
/* L2261: */
	}
/* L2281: */
    }
    flxu[nlev + l * 70 - 71] = flup[1];
    if (ibn == 1) {
	s_wsfe(&io___48);
	do_fio(&c__1, "from layer_flx: l = ", (ftnlen)20);
	do_fio(&c__1, (char *)&l, (ftnlen)sizeof(integer));
	e_wsfe();
/*      write(*,'(1a,16(1pe12.4))') 'dir:        ',(dir(k),k=1,nlev) */
	s_wsfe(&io___49);
	do_fio(&c__1, "dtau_b:     ", (ftnlen)12);
	i__1 = nlev - 1;
	for (k = 1; k <= i__1; ++k) {
	    do_fio(&c__1, (char *)&dtau_b__[k + (l + ibn * 5) * 70], (ftnlen)
		    sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___50);
	do_fio(&c__1, "flxd:       ", (ftnlen)12);
	i__1 = nlev;
	for (k = 1; k <= i__1; ++k) {
	    do_fio(&c__1, (char *)&flxd[k + l * 70 - 71], (ftnlen)sizeof(real)
		    );
	}
	e_wsfe();
	s_wsfe(&io___51);
	do_fio(&c__1, "flxu:       ", (ftnlen)12);
	i__1 = nlev;
	for (k = 1; k <= i__1; ++k) {
	    do_fio(&c__1, (char *)&flxu[k + l * 70 - 71], (ftnlen)sizeof(real)
		    );
	}
	e_wsfe();
/*      write(*,*) */
/*      do 1 nze=1,numu */
/*      naz = 1 */
/* 1     write(*,'(1a,i2,16(1pe12.4))') 'layer rad ', */
/*     -  nze,((rad(nze,naz,k)),k=1,nlev) */
	if (l > 1) {
	    i__1 = nlev - 1;
	    for (k = 1; k <= i__1; ++k) {
		deltau = dtau_b__[k + (l + ibn * 5) * 70] - dtau_b__[k + (ibn 
			* 5 + 1) * 70];
		if (deltau != 0.f) {
		    deltaui[k - 1] = 1.f / deltau;
		    dflxddx[k] = (flxd[k + 1 + l * 70 - 71] - flxd[k]) * 
			    deltaui[k - 1];
		    dflxudx[k - 1] = (flxu[k + l * 70 - 71] - flxu[k - 1]) * 
			    deltaui[k - 1];
		} else {
		    deltaui[k - 1] = 0.f;
		    dflxddx[k - 1] = 0.f;
		    dflxudx[k - 1] = 0.f;
		}
/* L3001: */
	    }
	    s_wsfe(&io___56);
	    do_fio(&c__1, "from layer_flx partials: l = ", (ftnlen)29);
	    do_fio(&c__1, (char *)&l, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___57);
	    do_fio(&c__1, "dflxddx:   ", (ftnlen)11);
	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&dflxddx[k - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    s_wsfe(&io___58);
	    do_fio(&c__1, "dflxudx:   ", (ftnlen)11);
	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&dflxudx[k - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    s_wsfe(&io___59);
	    do_fio(&c__1, "deltaui:   ", (ftnlen)11);
	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&deltaui[k - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
    }

    return 0;
} /* layer_flx__ */

