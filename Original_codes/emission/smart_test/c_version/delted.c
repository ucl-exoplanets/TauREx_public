/* delted.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int delted_(integer *nlay, integer *nza, doublereal *alb, 
	doublereal *u0, doublereal *solflx, doublereal *dtau, doublereal *om, 
	doublereal *g, doublereal *flxu, doublereal *flxd, doublereal *dir)
{
    /* Initialized data */

    static doublereal prec = 1e-10;
    static doublereal emax1 = 8.;
    static doublereal emax2 = 70.;

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer l, m;
    static doublereal a1[75], a2[75], g1[75], g2[75], g3[75], g4[75], da[75], 
	    db[75], dc[75], dd[75], gp[75], du[75], sk[75], rl[75], tl[75], 
	    rs[75], ru[75];
    static integer np2;
    static doublereal dfl[75], ekt[75], ufl[75];
    static integer nlm;
    static doublereal tau[75], gsq[75], omp[75], etu[75], skt;
    static integer nlev;
    static doublereal tdmu[75], e2ktm, denom, dtaup[75], dflux[75], uflux[75],
	     vlarge;


/* cccccccccccccccccccccccccc  d e l t e d  cccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes the upward, downward and net solar     cc */
/* c    fluxes in an inhomogeneous absorbing, scattering atmosphere     cc */
/* c    the delta-eddington approximation is used to find the diffuse   cc */
/* c    reflectivity and transmissivity and the total upward and        cc */
/* c    downward fluxes for each of the nlay homogeneous layers.        cc */
/* c    the adding method is then used to combine these layers.  if     cc */
/* c    any layer is thicker than dtau = *emax1*, it is assumed to be   cc */
/* c    semi-infinite.  layers thicker than dtau = *emax2*, are treated cc */
/* c    as infinite layers.                                             cc */
/* c                                                                    cc */
/* c    note: to account for a diffuse flux at the top of the atmos-    cc */
/* c          phere, the user must initialize flxd(1,m) to the appro-   cc */
/* c          priate value at all m zenith angles.  if there is no      cc */
/* c          downward diffuse flux at the top of the atmosphere, the   cc */
/* c          user must initialize flxd(1,m) to zero in the calling     cc */
/* c          program.                                                  cc */
/* c                                                                    cc */
/* c      nlay - number of homogeneous model layers.                    cc */
/* c       nza - number of zenith angles                                cc */
/* c       alb - surface albedo                                         cc */
/* c        u0 - array of cosines solar zenith angles                   cc */
/* c    solflx - normal incidence solar flux at top of atmosphere.      cc */
/* c      dtau - array of normal incidence optical depths in each       cc */
/* c             homogeneous model layer.                               cc */
/* c        om - array of single scattering albedos for each homo-      cc */
/* c             geneous model layer.                                   cc */
/* c         g - array of assymetry parameters for each homogeneous     cc */
/* c             model layer.                                           cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      flxu - upward flux at nlay+1 layer boundries.                 cc */
/* c             (flxu(l) refers to the upward flux at the top          cc */
/* c              of layer l)                                           cc */
/* c      flxd - downward flux at nlay+1 layer boundries.               cc */
/* c             (flxd(l) refers to the downward flux at the bottom     cc */
/* c              of layer l-1)                                         cc */
/* c       dir - direct solar flux at the top of each layer and at      cc */
/* c             surface.                                               cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  d e l t e d  cccccccccccccccccccccccccccccc */


/* ********** common blocks used in delta-eddington routines. */





    /* Parameter adjustments */
    dir -= 76;
    flxd -= 76;
    flxu -= 76;
    --g;
    --om;
    --dtau;
    --u0;

    /* Function Body */


    nlev = *nlay + 1;
    np2 = *nlay + 2;
    nlm = *nlay - 1;
    vlarge = exp(emax2);
    tau[0] = 0.f;

/* ****   scale the optical depths, single scattering albedos and the */
/*       scattering assymetry factors for use in the delta-eddington */
/*       approximation.  initialize other quantities.  use wiscombe's */
/*       trick to subtract a small value from the single scattering */
/*       for the case of a conservative atmosphere */

    i__1 = *nlay;
    for (l = 1; l <= i__1; ++l) {
	gsq[l - 1] = g[l] * g[l];
	gp[l - 1] = g[l] / (g[l] + 1.f);
	if (1.f - om[l] < prec) {
	    om[l] = 1.f - prec;
	}
	dtaup[l - 1] = (1.f - om[l] * gsq[l - 1]) * dtau[l];
	tau[l] = tau[l - 1] + dtaup[l - 1];
	omp[l - 1] = (1.f - gsq[l - 1]) * om[l] / (1.f - om[l] * gsq[l - 1]);
	g1[l - 1] = (7.f - omp[l - 1] * (gp[l - 1] * 3.f + 4.f)) * .25f;
	g2[l - 1] = (1.f - omp[l - 1] * (4.f - gp[l - 1] * 3.f)) * -.25f;
	sk[l - 1] = sqrt(g1[l - 1] * g1[l - 1] - g2[l - 1] * g2[l - 1]);
/* L1081: */
    }

/* **** diffuse transmittance and reflectance of each homogeneous layer */

/*     the equations for rl and tl were derived by solving the homogeneous */
/*     part of the equation of transfer (ie. no source term) */

    i__1 = *nlay;
    for (l = 1; l <= i__1; ++l) {
	skt = sk[l - 1] * dtaup[l - 1];
	if (skt > emax2) {
	    goto L1141;
	}
	ekt[l - 1] = exp(skt);
	if (skt > emax1) {
	    goto L1121;
	}
	e2ktm = ekt[l - 1] * ekt[l - 1] - 1.f;
	denom = g1[l - 1] * e2ktm + sk[l - 1] * (e2ktm + 2.f);
	rl[l - 1] = g2[l - 1] * e2ktm / denom;
	tl[l - 1] = sk[l - 1] * 2.f * ekt[l - 1] / denom;
	goto L1163;

/*     semi-infinite layers */

L1121:
	rl[l - 1] = g2[l - 1] / (g1[l - 1] + sk[l - 1]);
	tl[l - 1] = sk[l - 1] * 2.f / (ekt[l - 1] * (g1[l - 1] + sk[l - 1]));
	goto L1163;

/*     infintie layers */

L1141:
	rl[l - 1] = g2[l - 1] / (g1[l - 1] + sk[l - 1]);
	tl[l - 1] = 0.f;
	ekt[l - 1] = vlarge;
L1163:
	;
    }

/* ****   set the "reflectivity" and "transmissivity" of the surface. */

    rl[nlev - 1] = *alb;
    tl[nlev - 1] = 0.f;

/* ****   use adding method to find the reflectance and transmittance */
/*       of combined layers.  add downward from the top and upward */
/*       from the bottom at the same time. */

    rs[0] = rl[0];
    ru[nlev - 1] = *alb;
    i__1 = *nlay;
    for (l = 1; l <= i__1; ++l) {
	dd[l - 1] = 1.f / (1.f - rs[l - 1] * rl[l]);
	rs[l] = rl[l] + tl[l] * tl[l] * rs[l - 1] * dd[l - 1];
	du[nlev - l - 1] = 1.f / (1.f - rl[nlev - l - 1] * ru[np2 - l - 1]);
	ru[nlev - l - 1] = rl[nlev - l - 1] + tl[nlev - l - 1] * tl[nlev - l 
		- 1] * ru[np2 - l - 1] * du[nlev - l - 1];
/* L1221: */
    }

/* ****    z e n i t h    a n g l e    l o o p */

    i__1 = *nza;
    for (m = 1; m <= i__1; ++m) {

/* ****       initialize the direct solar flux at the top of each layer */
/*           and other quantities needed to find the diffuse fluxes. */

	dir[m * 75 + 1] = u0[m] * *solflx;
	i__2 = *nlay;
	for (l = 1; l <= i__2; ++l) {
	    g3[l - 1] = (2.f - gp[l - 1] * 3.f * u0[m]) * .25f;
	    g4[l - 1] = 1.f - g3[l - 1];
	    db[l - 1] = sk[l - 1] * u0[m] + 1.f;
	    a1[l - 1] = g1[l - 1] * g4[l - 1] + g2[l - 1] * g3[l - 1];
	    a2[l - 1] = g1[l - 1] * g3[l - 1] + g2[l - 1] * g4[l - 1];
	    dc[l - 1] = 2.f - db[l - 1];
	    ufl[l - 1] = 0.f;
	    dfl[l - 1] = 0.f;
	    dir[l + 1 + m * 75] = 0.f;
	    tdmu[l - 1] = dtaup[l - 1] / u0[m];
	    etu[l - 1] = 1.f / vlarge;
	    if (tdmu[l - 1] <= emax2) {
		etu[l - 1] = exp(-tdmu[l - 1]);
	    }
	    da[l - 1] = db[l - 1] * dc[l - 1] * ((sk[l - 1] + g1[l - 1]) * 
		    ekt[l - 1] + (sk[l - 1] - g1[l - 1]) / ekt[l - 1]);
/* L1502: */
	}
	i__2 = nlev;
	for (l = 1; l <= i__2; ++l) {
	    tdmu[l - 1] = tau[l - 1] / u0[m];
/* L1521: */
	}
/* L1551: */
	i__2 = nlev;
	for (l = 2; l <= i__2; ++l) {
	    if (tdmu[l - 1] <= emax2) {
		dir[l + m * 75] = *solflx * u0[m] * exp(-tdmu[l - 1]);
	    }
/* L1561: */
	}

/* ****       find the upward flux at the top and the downward */
/*           flux at the bottom of each homogeneous model layer. */

	ufl[nlev - 1] = *alb * dir[nlev + m * 75];
	dfl[nlev - 1] = 0.f;
	i__2 = *nlay;
	for (l = 1; l <= i__2; ++l) {
	    ufl[l - 1] = omp[l - 1] * dir[l + m * 75] / da[l - 1] * (dc[l - 1]
		     * (a2[l - 1] + sk[l - 1] * g3[l - 1]) * ekt[l - 1] - db[
		    l - 1] * (a2[l - 1] - sk[l - 1] * g3[l - 1]) / ekt[l - 1] 
		    - sk[l - 1] * 2.f * (g3[l - 1] - a2[l - 1] * u0[m]) * etu[
		    l - 1]);
	    dfl[l - 1] = -dir[l + 1 + m * 75] * omp[l - 1] / da[l - 1] * (db[
		    l - 1] * (a1[l - 1] + sk[l - 1] * g4[l - 1]) * ekt[l - 1] 
		    - dc[l - 1] * (a1[l - 1] - sk[l - 1] * g4[l - 1]) / ekt[l 
		    - 1] - sk[l - 1] * 2.f * (g4[l - 1] + a1[l - 1] * u0[m]) /
		     etu[l - 1]);
/* L1621: */
	}

/* ****       use adding method to find upward and downward fluxes */
/*           for combined layers.  start at top */

	dflux[0] = dfl[0] + tl[0] * flxd[m * 75 + 1];
	i__2 = nlm;
	for (l = 1; l <= i__2; ++l) {
/*              dflux(l+1) = tl(l+1)*(rs(l)*(rl(l+1)*dflux(l) + */
/*     *                     ufl(l+1))*dd(l) + dflux(l)) + dfl(l+1) */
	    dflux[l] = dfl[l] + tl[l] * (dflux[(integer) (l + rs[l - 1] * ufl[
		    l]) - 1] * dd[l - 1]);
/* L1701: */
	}
	flxd[nlev + m * 75] = (dflux[*nlay - 1] + rs[*nlay - 1] * dir[nlev + 
		m * 75] * *alb) / (1.f - rs[*nlay - 1] * *alb);

/* ****       use adding method to find upward and downward fluxes */
/*           for combined layers.  start at bottom. */

	uflux[nlev - 1] = ufl[nlev - 1];
	i__2 = *nlay;
	for (l = 1; l <= i__2; ++l) {
	    uflux[nlev - l - 1] = tl[nlev - l - 1] * (uflux[np2 - l - 1] + ru[
		    np2 - l - 1] * dfl[nlev - l - 1]) * du[nlev - l - 1] + 
		    ufl[nlev - l - 1];
/* L1721: */
	}

/* ****       find the total upward and downward fluxes at interfaces */
/*           between inhomogeneous layers. */

	flxu[m * 75 + 1] = uflux[0] + ru[0] * flxd[m * 75 + 1];
	flxu[nlev + m * 75] = *alb * (flxd[nlev + m * 75] + dir[nlev + m * 75]
		);
	i__2 = nlm;
	for (l = 1; l <= i__2; ++l) {
	    flxu[nlev - l + m * 75] = (uflux[nlev - l - 1] + ru[nlev - l - 1] 
		    * dflux[*nlay - l - 1]) / (1.f - ru[nlev - l - 1] * rs[*
		    nlay - l - 1]);
/* L1741: */
	}
	i__2 = nlm;
	for (l = 1; l <= i__2; ++l) {
	    flxd[l + 1 + m * 75] = dflux[l - 1] + rs[l - 1] * flxu[l + 1 + m *
		     75];
/* L1782: */
	}
    }

    return 0;
} /* delted_ */

