/* bdrf.f -- translated by f2c (version 20100827).
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

static doublereal c_b2 = .5;
static doublereal c_b9 = 2.;

doublereal bdref_(real *mu, real *mup, real *dphi, real *surf_pr__, integer *
	iref)
{
    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double cos(doublereal), pow_dd(doublereal *, doublereal *), acos(
	    doublereal), tan(doublereal), sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int sunglint_(real *, real *, real *, real *, 
	    real *, real *, real *, real *);
    static real b, h__, p, w, b0;
    extern doublereal f2_(real *, real *, real *);
    static real h0, k0, k1, k2, fi, hh, ni, nr, ts, tv;
    extern doublereal f1_2__(real *, real *, real *), f1_3__(real *, real *, 
	    real *);
    static real rog, azw, mup1, th_s__, th_v__, wspd, gamma, theta, ctheta;

/*      Supplies surface bi-directional reflectivity. */

/*    REFMODE  : bidirectional reflectance options */
/*             0 - Lambert */
/*             1 - Hapke's BDR model */
/*             2 - Breon's BDR model; combination of Li + Roujean */
/*             3 - Roujean's BDR model */
/*             4 - Cox and Munk glint model */

/*      NOTE 1: Bidirectional reflectivity in DISORT is defined */
/*              by Eq. 39 in STWL. */
/*      NOTE 2: Both MU and MUP (cosines of reflection and incidence */
/*              angles) are positive. */

/*  INPUT: */

/*    MU     : Cosine of angle of reflection (positive) */

/*    MUP    : Cosine of angle of incidence (positive) */

/*    DPHI   : Difference of azimuth angles of incidence and reflection */
/*                (radians) */
/*   SURF_PR : Wavelength dependent surface properties array */
/*             IREF= 0 - Lambert albedo */
/*             IREF= 1 - Hapke : single scatter albedo, W, and */
/*                               angular width factr HH */
/*             IREF= 2 - Breon's BDR model: k0, k1, k2 */
/*             IREF= 3 - Roujean's BDR model: k0, k1, k2 */
/*             IREF= 4 - Cox and Munk glint model: n, k, ws, phiw */

/*  LOCAL VARIABLES: */

/*    IREF   : bidirectional reflectance options */

/*    B0     : empirical factor to account for the finite size of */
/*             particles in Hapke's BDR model */

/*    B      : term that accounts for the opposition effect */
/*             (retroreflectance, hot spot) in Hapke's BDR model */

/*    CTHETA : cosine of phase angle in Hapke's BDR model */

/*    GAMMA  : albedo factor in Hapke's BDR model */

/*    H0     : H( mu0 ) in Hapke's BDR model */

/*    H      : H( mu ) in Hapke's BDR model */

/*    HH     : angular width parameter of opposition effect in Hapke's */
/*             BDR model */

/*    P      : scattering phase function in Hapke's BDR model */

/*    THETA  : phase angle (radians); the angle between incidence and */
/*             reflection directions in Hapke's BDR model */

/*    W      : single scattering albedo in Hapke's BDR model */

/*    ws     : wind speed (m/s) in the Cox Munk model */

/*    phiw   : wind direction (radians) in the Cox Munk model */

/*   Called by- DREF, SURFAC */
/* +-------------------------------------------------------------------+ */



/*     .. Scalar Arguments .. */


/*     .. Local Scalars .. */


/*     .. Intrinsic Functions .. */


    /* Parameter adjustments */
    --surf_pr__;

    /* Function Body */
    if (*iref == 1) {
/*        ** Hapke's BRDF model (times Pi/Mu0) */
/*        ** (Hapke, B., Theory of reflectance */
/*        ** and emittance spectroscopy, Cambridge */
/*        ** University Press, 1993, Eq. 8.89 on */
/*        ** page 233. Parameters are from */
/*        ** Fig. 8.15 on page 231, expect for w.) */
/* Computing 2nd power */
	r__1 = *mu;
	d__1 = (doublereal) (1.f - r__1 * r__1);
/* Computing 2nd power */
	r__2 = *mup;
	d__2 = (doublereal) (1.f - r__2 * r__2);
	ctheta = *mu * *mup + pow_dd(&d__1, &c_b2) * pow_dd(&d__2, &c_b2) * 
		cos(*dphi);
	theta = acos(ctheta);
	p = ctheta * .5f + 1.f;
	hh = surf_pr__[2];
/* 0.06 */
	b0 = 1.f;
	b = b0 * hh / (hh + tan(theta / 2.f));
	w = surf_pr__[1];
/* 0.6 */
	gamma = sqrt(1.f - w);
	h0 = (*mup * 2.f + 1.f) / (*mup * 2.f * gamma + 1.f);
	h__ = (*mu * 2.f + 1.f) / (*mu * 2.f * gamma + 1.f);
	ret_val = w / 4.f / (*mu + *mup) * ((b + 1.f) * p + h0 * h__ - 1.f);
    }

    if (*iref == 2) {
	k0 = surf_pr__[1];
/* 0.3530000 */
	k1 = surf_pr__[2];
/* 5.7000000E-02 */
	k2 = surf_pr__[3];
/* 0.3590000 */
	th_v__ = acos(*mu);
	th_s__ = acos(*mup);

	ret_val = k0 + k1 * f1_2__(&th_s__, &th_v__, dphi) + k2 * f2_(&th_s__,
		 &th_v__, dphi);

    }

    if (*iref == 3) {

/*             Roujean's BDR model */

	k0 = surf_pr__[1];
/* 0.3530000 */
	k1 = surf_pr__[2];
/* 5.7000000E-02 */
	k2 = surf_pr__[3];
/* 0.3590000 */
	th_v__ = acos(*mu);
	th_s__ = acos(*mup);
	ret_val = k0 + k1 * f1_3__(&th_s__, &th_v__, dphi) + k2 * f2_(&th_s__,
		 &th_v__, dphi);

    }

    if (*iref == 4) {

/*             Cox and Munk glint model from 6s */
/* input parameters:   wspd=speed of the wind (in m/s) */
/*                     nr=index of refraction of the sea water */
/*                     ni=extinction coefficient of the sea water */
/*                     azw=azim. of the sun - azim. of the wind (in deg.) */
/*                     ts=solar zenith angle (in deg.) */
/*                     tv=view zenith angle (in deg.) */
/*                     fi=relative azimuth (sun-viewing angle) */
/* output parameters:  rog=reflectance of the sun glint */

/* *****  note, for near glancing angles (mup < 0.03) this routine */
/*       produces flux albedo values > 1.0.  To fix this, set all */
/*       near glancing agles to mup > 0.3 */

	if (*mup >= .08f) {
	    mup1 = *mup;
	} else {
	    mup1 = .08f;
	}
	nr = surf_pr__[1];
	ni = surf_pr__[2];
	wspd = surf_pr__[3];
	azw = surf_pr__[4];
	ts = acos(mup1) * 180.f / 3.14159265358979323f;
	tv = acos(*mu) * 180.f / 3.14159265358979323f;
	fi = *dphi * 180.f / 3.14159265358979323f;

	sunglint_(&wspd, &nr, &ni, &azw, &ts, &tv, &fi, &rog);
/*       write(*,'(10(1pe12.4))') wspd,nr,ni,azw,ts,tv,fi,rog */
	if (rog < 0.f) {
	    rog = 0.f;
	}
	ret_val = rog;

    }
    return ret_val;
} /* bdref_ */


/* FUNCTION BDREF */
doublereal f1_2__(real *theta_s__, real *theta_v__, real *phi)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double tan(doublereal), cos(doublereal), sin(doublereal), acos(doublereal)
	    , sqrt(doublereal);

    /* Local variables */
    static real q, t, x, kk, dumy, invpi, costt, cos_ksi__;


/*     Li-Sparse F1 */

/*  INPUT: */

/*    THETA_S : Cosine of angle of incidence (positive) */

/*    THETA_V : Cosine of angle of reflection (positive) */

/*    PHI     : azimuth angles of incidence and reflection (radians) */

/* +-------------------------------------------------------------------+ */

/*     .. Scalar Arguments .. */

    dumy = 1.2217310000000001f;
    if (*theta_s__ > dumy) {
	*theta_s__ = dumy;
    }
    if (*theta_v__ > dumy) {
	*theta_v__ = dumy;
    }
    invpi = .31830988618379069f;
/* Computing 2nd power */
    r__1 = tan(*theta_s__);
/* Computing 2nd power */
    r__2 = tan(*theta_v__);
    x = r__1 * r__1 + r__2 * r__2 - tan(*theta_s__) * 2 * tan(*theta_v__) * 
	    cos(*phi);
    cos_ksi__ = cos(*theta_s__) * cos(*theta_v__) + sin(*theta_s__) * sin(*
	    theta_v__) * cos(*phi);
    kk = acos(cos(*theta_s__) * cos(*theta_v__) + sin(*theta_s__) * sin(*
	    theta_v__) * cos(*phi));
    kk = acos(cos_ksi__);
/* Computing 2nd power */
    r__1 = tan(*theta_s__) * tan(*theta_v__) * sin(*phi);
    costt = sqrt(x + r__1 * r__1) * 2.f / (1.f / cos(*theta_s__) + 1.f / cos(*
	    theta_v__));
    if (costt > 1.f) {
	costt = 1.f;
    }
    t = acos(costt);
    q = invpi * (t - sin(t) * cos(t)) * (1.f / cos(*theta_s__) + 1.f / cos(*
	    theta_v__));
    ret_val = q - 1.f / cos(*theta_v__) - 1.f / cos(*theta_s__) + (cos(kk) + 
	    1) * .5f * (1.f / cos(*theta_s__)) * (1.f / cos(*theta_v__));
    return ret_val;
} /* f1_2__ */

/* FUNCTION F1_2 */
doublereal f1_3__(real *theta_s__, real *theta_v__, real *phi)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double tan(doublereal), cos(doublereal), sin(doublereal), sqrt(doublereal)
	    ;

    /* Local variables */
    static real x, dumy, invpi, inv2pi;

/*     Roujean F1 */

/*  INPUT: */

/*    THETA_S : Cosine of angle of incidence (positive) */

/*    THETA_V : Cosine of angle of reflection (positive) */

/*    PHI     : azimuth angles of incidence and reflection (radians) */

/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
    dumy = 1.2217310000000001f;
    if (*theta_s__ > dumy) {
	*theta_s__ = dumy;
    }
    if (*theta_v__ > dumy) {
	*theta_v__ = dumy;
    }
    invpi = .31830988618379069f;
    inv2pi = .15915494309189535f;
/* Computing 2nd power */
    r__1 = tan(*theta_s__);
/* Computing 2nd power */
    r__2 = tan(*theta_v__);
    x = r__1 * r__1 + r__2 * r__2 - tan(*theta_s__) * 2 * tan(*theta_v__) * 
	    cos(*phi);
    ret_val = inv2pi * ((3.14159265358979323f - *phi) * cos(*phi) + sin(*phi))
	     * tan(*theta_s__) * tan(*theta_v__) - invpi * (tan(*theta_s__) + 
	    tan(*theta_v__) + sqrt(x));
    return ret_val;
} /* f1_3__ */

/* FUNCTION F1_3 */
doublereal f2_(real *theta_s__, real *theta_v__, real *phi)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), acos(doublereal);

    /* Local variables */
    static real x, dumy, cos_ksi__;

/*     Roujean F2 */

/*  INPUT: */

/*    THETA_S : Cosine of angle of incidence (positive) */

/*    THETA_V : Cosine of angle of reflection (positive) */

/*    PHI     : azimuth angles of incidence and reflection (radians) */

/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */

    dumy = 1.2217310000000001f;
    if (*theta_s__ > dumy) {
	*theta_s__ = dumy;
    }
    if (*theta_v__ > dumy) {
	*theta_v__ = dumy;
    }
    cos_ksi__ = cos(*theta_s__) * cos(*theta_v__) + sin(*theta_s__) * sin(*
	    theta_v__) * cos(*phi);
    x = acos(cos(*theta_s__) * cos(*theta_v__) + sin(*theta_s__) * sin(*
	    theta_v__) * cos(*phi));
    x = acos(cos_ksi__);
    ret_val = 4 / ((cos(*theta_s__) + cos(*theta_v__)) * 9.4247779607693793f) 
	    * ((1.5707963267948966f - x) * cos(x) + sin(x)) - 
	    .33333333333333331f;
    return ret_val;
} /* f2_ */

/* FUNCTION F2 */
/* Subroutine */ int sunglint_(real *wspd, real *nr, real *ni, real *azw, 
	real *ts, real *tv, real *fi, real *rog)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double acos(doublereal), cos(doublereal), sin(doublereal), sqrt(
	    doublereal), atan(doublereal), exp(doublereal);

    /* Local variables */
    static real r1;
    static doublereal c21, c03, c40, c04, c22, cs, pi, cv, xe, ss, xn, sv, zx,
	     zy, xe2, xn2, fac, phi, phw, coef, tilt, proba, sigmac;
    static real coschi, sinchi;
    static doublereal sigmau, cos2chi;
    extern /* Subroutine */ int fresnel_(real *, real *, real *, real *, real 
	    *);
    static doublereal tantilt;


/* input parameters:   wspd=speed of the wind (in m/s) */
/*                     nr=index of refraction of the sea water */
/*                     ni=extinction coefficient of the sea water */
/*                     azw=azim. of the sun - azim. of the wind (in deg.) */
/*                     ts=solar zenith angle (in deg.) */
/*                     tv=view zenith angle (in deg.) */
/*                     fi=relative azimuth (sun-satellite) */
/* output parameters:  rog=reflectance of the sun glint */



    pi = acos(-1.);
    fac = pi / 180.;
    phw = *azw * fac;
    cs = cos(*ts * fac);
    cv = cos(*tv * fac);
    ss = sin(*ts * fac);
    sv = sin(*tv * fac);
    phi = *fi * fac;
    zx = -sv * sin(phi) / (cs + cv);
    zy = (ss + sv * cos(phi)) / (cs + cv);
    tantilt = sqrt(zx * zx + zy * zy);
    tilt = atan(tantilt);

/*  Anisotropic Gaussian distribution */
/*    phw=phi_sun-phi_wind */

    sigmac = *wspd * .00192f + .003;
    sigmau = *wspd * .00316;
    c21 = .01 - *wspd * .0086;
    c03 = .04 - *wspd * .033;
    c40 = .4;
    c22 = .12;
    c04 = .23;
    xe = (cos(phw) * zx + sin(phw) * zy) / sqrt(sigmac);
    xn = (-sin(phw) * zx + cos(phw) * zy) / sqrt(sigmau);
    xe2 = xe * xe;
    xn2 = xn * xn;
    coef = 1 - c21 / 2.f * (xe2 - 1) * xn - c03 / 6.f * (xn2 - 3) * xn;
    coef += c40 / 24. * (xe2 * xe2 - xe2 * 6 + 3);
    coef += c04 / 24. * (xn2 * xn2 - xn2 * 6 + 3);
    coef += c22 / 4. * (xe2 - 1) * (xn2 - 1);
    proba = coef / 2. / pi / sqrt(sigmau) / sqrt(sigmac) * exp(-(xe2 + xn2) / 
	    2.);

/* Compute Fresnel's coefficient R1 */

    cos2chi = cv * cs + sv * ss * cos(phi);
    if (cos2chi > 1.f) {
	cos2chi = .99999999999f;
    }
    if (cos2chi < -1.f) {
	cos2chi = -.99999999999f;
    }
    coschi = sqrt((cos2chi + 1.) * .5);
    sinchi = sqrt((1. - cos2chi) * .5);

    fresnel_(nr, ni, &coschi, &sinchi, &r1);

/* Compute Reflectance of the sun glint */

/* Computing 4th power */
    d__1 = cos(tilt), d__1 *= d__1;
    *rog = pi * .25f * r1 * proba / cs / cv / (d__1 * d__1);
/*      write(*,'(10(1pe12.4))') pi,r1,proba,cs,cv,tilt,cos(tilt),Rog */

    return 0;
} /* sunglint_ */



/* Subroutine */ int fresnel_(real *nr, real *ni, real *coschi, real *sinchi, 
	real *r1)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static real u, v, a1, a2, b1, b2, rl2, rr2;


/* to compute the Fresnel's coefficient of reflection (see for */
/* example M. Born and E. Wolf, Principles of Optics, Pergamon Press, fifth */
/* edition, 1975, pp 628 */
/* input parameters: nr=index of refraction of the sea water */
/*                   ni=extinction coefficient of the sea water */
/*                   coschi & sinchi=cosine and sine of the incident radiation */
/*                                   with respect of the wave facet normal. */
/* output parameter: R1=Fresnel's coefficient for reflection */

/* absolute value for a1 to get v=0 when ni=0 */
    a1 = (r__1 = *nr * *nr - *ni * *ni - *sinchi * *sinchi, dabs(r__1));
    d__1 = (doublereal) (*nr * *nr - *ni * *ni - *sinchi * *sinchi);
    a2 = sqrt(pow_dd(&d__1, &c_b9) + *nr * 4 * *nr * *ni * *ni);
    u = sqrt((a1 + a2) * .5f);
    v = sqrt((-a1 + a2) * .5f);
/* Computing 2nd power */
    r__1 = *coschi - u;
/* Computing 2nd power */
    r__2 = *coschi + u;
    rr2 = (r__1 * r__1 + v * v) / (r__2 * r__2 + v * v);
    b1 = (*nr * *nr - *ni * *ni) * *coschi;
    b2 = *nr * 2 * *ni * *coschi;
/* Computing 2nd power */
    r__1 = b1 - u;
/* Computing 2nd power */
    r__2 = b2 + v;
/* Computing 2nd power */
    r__3 = b1 + u;
/* Computing 2nd power */
    r__4 = b2 - v;
    rl2 = (r__1 * r__1 + r__2 * r__2) / (r__3 * r__3 + r__4 * r__4);
    *r1 = (rr2 + rl2) / 2.f;
    return 0;
} /* fresnel_ */

