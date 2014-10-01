/* fresnel.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int fresnel_(real *rn1, real *rk1, real *ang1, real *rn2, 
	real *rk2, real *ang2, real *rpar, real *rper, real *tpar, real *tper)
{
    /* System generated locals */
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *);
    double cos(doublereal);
    void c_sqrt(complex *, complex *), r_cnjg(complex *, complex *);
    double sin(doublereal), asin(doublereal);

    /* Local variables */
    static complex n1, n2, rp, tp, rs, ts;
    static real arg, cos1;
    static complex cos2, refrac;


/* cccccccccccccccccccccccccc  f r e s n e l  ccccccccccccccccccccccccccccc */
/* c                                                                     cc */
/* c    p u r p o s e :                                                  cc */
/* c                                                                     cc */
/* c    this subroutine evaluates fresnel's reflection formula for       cc */
/* c    arbitrary angle and refractive index, and returns the            cc */
/* c    parallel and perpendictular polarized components of the          cc */
/* c    reflectance.                                                     cc */
/* c                                                                     cc */
/* c    i n p u t :                                                      cc */
/* c                                                                     cc */
/* c     rn1 - real part of the complex refractive index of 1st medium   cc */
/* c     rk1 - complex part of the complex refractive index of 1st       cc */
/* c           medium                                                    cc */
/* c     rn2 - real part of the complex refractive index of 2nd medium   cc */
/* c     rk2 - complex part of the complex refractive index of 2nd       cc */
/* c           medium                                                    cc */
/* c    ang1 - angle of incidence of the wave to the surface normal      cc */
/* c           (in radians) in  first media                              cc */
/* c                                                                     cc */
/* c    o u t p u t :                                                    cc */
/* c                                                                     cc */
/* c    ang2 - angle of reflection of the wave to the surface normal     cc */
/* c           (in radians) in second media                              cc */
/* c    rpar - parallel-polarized component of the reflectance           cc */
/* c    rper - perpendictular-polarized component of the reflectance     cc */
/* c    tpar - parallel-polarized component of the transmittance         cc */
/* c    tper - perpendictular-polarized component of the transmittance   cc */
/* c                                                                     cc */
/* c    a u t h o r :                                                    cc */
/* c                                                                     cc */
/* c    dave crisp  4/29/89                                              cc */
/* c      update 8/4/90 to fix tpar,tper                                 cc */
/* c                                                                     cc */
/* c    r e f e r e n c e s :                                            cc */
/* c  Wendlandt and Hecht, Reflectance Spectroscopy, p.17                cc */
/* c                                                                     cc */
/* c  Bennett J.M. and H.E. Bennet, "Polarization," Handbook of optics,  cc */
/* c      W. Driscoll ed., Mcgraw-Hill Book Company, NY, pg 10-7, 1978.  cc */
/* c                                                                     cc */
/* c  Schanda, E., Physical Fundamentals of Remote Sensing, Springer     cc */
/* c      Verlag, Berlin, pg 28, 1986, 187 pp.                           cc */
/* c                                                                     cc */
/* cccccccccccccccccccccccccc  f r e s n e l  ccccccccccccccccccccccccccccc */



    if (*ang1 > 1.55f) {

/* *****   if ang1 is close to 90 degrees, assume reflectance is unity. */
/*        set up the correct values for rpar,rper, tpar, tper, */
/*        solve for ang2 and return */

	*rpar = 1.f;
	*rper = 1.f;
	*tpar = 0.f;
	*tper = 0.f;

    } else {

/* *****   define the complex index of refraction for the two media */

	q__2.r = *rk1 * 0.f, q__2.i = *rk1 * 1.f;
	q__1.r = *rn1 - q__2.r, q__1.i = -q__2.i;
	n1.r = q__1.r, n1.i = q__1.i;
	q__2.r = *rk2 * 0.f, q__2.i = *rk2 * 1.f;
	q__1.r = *rn2 - q__2.r, q__1.i = -q__2.i;
	n2.r = q__1.r, n2.i = q__1.i;

	c_div(&q__1, &n1, &n2);
	refrac.r = q__1.r, refrac.i = q__1.i;
	cos1 = cos(*ang1);
	q__4.r = refrac.r * refrac.r - refrac.i * refrac.i, q__4.i = refrac.r 
		* refrac.i + refrac.i * refrac.r;
/* Computing 2nd power */
	r__2 = cos1;
	r__1 = r__2 * r__2;
	q__5.r = 1.f - r__1, q__5.i = 0.f;
	q__3.r = q__4.r * q__5.r - q__4.i * q__5.i, q__3.i = q__4.r * q__5.i 
		+ q__4.i * q__5.r;
	q__2.r = 1.f - q__3.r, q__2.i = 0.f - q__3.i;
	c_sqrt(&q__1, &q__2);
	cos2.r = q__1.r, cos2.i = q__1.i;

/* ****     define the parallel-polarized component of the reflectance */

	q__3.r = cos1 * n2.r, q__3.i = cos1 * n2.i;
	q__4.r = n1.r * cos2.r - n1.i * cos2.i, q__4.i = n1.r * cos2.i + n1.i 
		* cos2.r;
	q__2.r = q__3.r - q__4.r, q__2.i = q__3.i - q__4.i;
	q__6.r = cos1 * n2.r, q__6.i = cos1 * n2.i;
	q__7.r = n1.r * cos2.r - n1.i * cos2.i, q__7.i = n1.r * cos2.i + n1.i 
		* cos2.r;
	q__5.r = q__6.r + q__7.r, q__5.i = q__6.i + q__7.i;
	c_div(&q__1, &q__2, &q__5);
	rp.r = q__1.r, rp.i = q__1.i;
	r_cnjg(&q__2, &rp);
	q__1.r = rp.r * q__2.r - rp.i * q__2.i, q__1.i = rp.r * q__2.i + rp.i 
		* q__2.r;
	*rpar = q__1.r;

/* ****     define the perpendictular-polarized component of  reflectance */

	q__3.r = cos1 * n1.r, q__3.i = cos1 * n1.i;
	q__4.r = n2.r * cos2.r - n2.i * cos2.i, q__4.i = n2.r * cos2.i + n2.i 
		* cos2.r;
	q__2.r = q__3.r - q__4.r, q__2.i = q__3.i - q__4.i;
	q__6.r = cos1 * n1.r, q__6.i = cos1 * n1.i;
	q__7.r = n2.r * cos2.r - n2.i * cos2.i, q__7.i = n2.r * cos2.i + n2.i 
		* cos2.r;
	q__5.r = q__6.r + q__7.r, q__5.i = q__6.i + q__7.i;
	c_div(&q__1, &q__2, &q__5);
	rs.r = q__1.r, rs.i = q__1.i;
	r_cnjg(&q__2, &rs);
	q__1.r = rs.r * q__2.r - rs.i * q__2.i, q__1.i = rs.r * q__2.i + rs.i 
		* q__2.r;
	*rper = q__1.r;

/* ****   define the parallel-polarized component of the transmittance */

	q__3.r = n1.r * 2.f, q__3.i = n1.i * 2.f;
	q__2.r = cos1 * q__3.r, q__2.i = cos1 * q__3.i;
	q__5.r = cos1 * n2.r, q__5.i = cos1 * n2.i;
	q__6.r = n1.r * cos2.r - n1.i * cos2.i, q__6.i = n1.r * cos2.i + n1.i 
		* cos2.r;
	q__4.r = q__5.r + q__6.r, q__4.i = q__5.i + q__6.i;
	c_div(&q__1, &q__2, &q__4);
	tp.r = q__1.r, tp.i = q__1.i;
	q__4.r = n2.r * cos2.r - n2.i * cos2.i, q__4.i = n2.r * cos2.i + n2.i 
		* cos2.r;
	q__5.r = cos1 * n1.r, q__5.i = cos1 * n1.i;
	c_div(&q__3, &q__4, &q__5);
	q__2.r = q__3.r * tp.r - q__3.i * tp.i, q__2.i = q__3.r * tp.i + 
		q__3.i * tp.r;
	r_cnjg(&q__6, &tp);
	q__1.r = q__2.r * q__6.r - q__2.i * q__6.i, q__1.i = q__2.r * q__6.i 
		+ q__2.i * q__6.r;
	*tpar = q__1.r;

/* ****     define perpendictular-polarized component of  transmittance */

	q__3.r = n1.r * 2.f, q__3.i = n1.i * 2.f;
	q__2.r = cos1 * q__3.r, q__2.i = cos1 * q__3.i;
	q__5.r = cos1 * n1.r, q__5.i = cos1 * n1.i;
	q__6.r = n2.r * cos2.r - n2.i * cos2.i, q__6.i = n2.r * cos2.i + n2.i 
		* cos2.r;
	q__4.r = q__5.r + q__6.r, q__4.i = q__5.i + q__6.i;
	c_div(&q__1, &q__2, &q__4);
	ts.r = q__1.r, ts.i = q__1.i;
	q__4.r = n2.r * cos2.r - n2.i * cos2.i, q__4.i = n2.r * cos2.i + n2.i 
		* cos2.r;
	q__5.r = cos1 * n1.r, q__5.i = cos1 * n1.i;
	c_div(&q__3, &q__4, &q__5);
	q__2.r = q__3.r * ts.r - q__3.i * ts.i, q__2.i = q__3.r * ts.i + 
		q__3.i * ts.r;
	r_cnjg(&q__6, &ts);
	q__1.r = q__2.r * q__6.r - q__2.i * q__6.i, q__1.i = q__2.r * q__6.i 
		+ q__2.i * q__6.r;
	*tper = q__1.r;

    }

    arg = refrac.r * sin(*ang1);
    if (arg > -1.f && arg < 1.f) {
	*ang2 = asin(refrac.r * sin(*ang1));
    } else {
	*ang2 = 1.5707963f;
    }

    return 0;
} /* fresnel_ */

