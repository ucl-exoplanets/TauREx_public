/* plkavg.f -- translated by f2c (version 20100827).
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

static integer c__2 = 2;
static integer c__4 = 4;
static logical c_true = TRUE_;
static logical c_false = FALSE_;

doublereal plkavg_(real *wnumlo, real *wnumhi, real *t)
{
    /* Initialized data */

    static real c2 = 1.438786f;
    static real sigma = 5.67032e-8f;
    static real vcut = 1.5f;
    static real vcp[7] = { 10.25f,5.7f,3.9f,2.9f,2.3f,1.9f,0.f };
    static real pi = 0.f;

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double asin(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static real d__[2];
    static integer i__, k, m, n;
    static real p[2], v[2], hh, ex, mv, del, val, exm, vsq, val0, conc;
    static integer mmax;
    static real vmax, epsil;
    extern doublereal r1mach_(integer *);
    static real sigdpi, oldval;
    static integer smallv;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

/*        computes planck function integrated between two wavenumbers */

/*  input :  wnumlo : lower wavenumber (inv cm) of spectral interval */

/*           wnumhi : upper wavenumber */

/*           t      : temperature (k) */

/*  output : plkavg : integrated planck function ( watts/sq m ) */
/*                      = integral (wnumlo to wnumhi) of */
/*                        2h c**2  nu**3 / ( exp(hc nu/kt) - 1) */
/*                        (where h=plancks constant, c=speed of */
/*                         light, nu=wavenumber, t=temperature, */
/*                         and k = boltzmann constant) */

/*  reference : specifications of the physical world: new value */
/*                 of the fundamental constants, dimensions/n.b.s., */
/*                 jan. 1974 */

/*  method :  for wnumlo close to wnumhi, a simpson-rule quadrature */
/*            is done to avoid ill-conditioning; otherwise */

/*            (1)  for wnumlo or wnumhi small, */
/*                 integral(0 to wnumlo/hi) is calculated by expanding */
/*                 the integrand in a power series and integrating */
/*                 term by term; */

/*            (2)  otherwise, integral(wnumlo/hi to infinity) is */
/*                 calculated by expanding the denominator of the */
/*                 integrand in powers of the exponential and */
/*                 integrating term by term. */

/*  accuracy :  at least 6 significant digits, assuming the */
/*              physical constants are infinitely accurate */

/*  errors which are not trapped: */

/*      * power or exponential series may underflow, giving no */
/*        significant digits.  this may or may not be of concern, */
/*        depending on the application. */

/*      * simpson-rule special case is skipped when denominator of */
/*        integrand will cause overflow.  in that case the normal */
/*        procedure is used, which may be inaccurate if the */
/*        wavenumber limits (wnumlo, wnumhi) are close together. */

/*  local variables */

/*        a1,2,... :  power series coefficients */
/*        c2       :  h * c / k, in units cm*k (h = plancks constant, */
/*                      c = speed of light, k = boltzmann constant) */
/*        d(i)     :  exponential series expansion of integral of */
/*                       planck function from wnumlo (i=1) or wnumhi */
/*                       (i=2) to infinity */
/*        epsil    :  smallest number such that 1+epsil .gt. 1 on */
/*                       computer */
/*        ex       :  exp( - v(i) ) */
/*        exm      :  ex**m */
/*        mmax     :  no. of terms to take in exponential series */
/*        mv       :  multiples of v(i) */
/*        p(i)     :  power series expansion of integral of */
/*                       planck function from zero to wnumlo (i=1) or */
/*                       wnumhi (i=2) */
/*        pi       :  3.14159... */
/*        sigma    :  stefan-boltzmann constant (w/m**2/k**4) */
/*        sigdpi   :  sigma / pi */
/*        smallv   :  number of times the power series is used (0,1,2) */
/*        v(i)     :  c2 * (wnumlo(i=1) or wnumhi(i=2)) / temperature */
/*        vcut     :  power-series cutoff point */
/*        vcp      :  exponential series cutoff points */
/*        vmax     :  largest allowable argument of exp function */

/*   called by- disort */
/*   calls- r1mach, errmsg */
/* ---------------------------------------------------------------------- */
/*     .. parameters .. */
/*     .. */
/*     .. scalar arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. local arrays .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
/*     .. statement functions .. */
/*     .. */
/*     .. statement function definitions .. */
/*     .. */
    if (pi == 0.f) {
	pi = asin(1.f) * 2.f;
	vmax = log(r1mach_(&c__2));
	epsil = r1mach_(&c__4);
	sigdpi = sigma / pi;
/* Computing 4th power */
	r__1 = pi, r__1 *= r__1;
	conc = 15.f / (r__1 * r__1);
    }
    if (*t < 0.f || *wnumhi <= *wnumlo || *wnumlo < 0.f) {
	errmsg_("plkavg--temperature or wavenums. wrong", &c_true, (ftnlen)38)
		;
    }
    if (*t < 1e-4f) {
	ret_val = 0.f;
	return ret_val;
    }
    v[0] = c2 * *wnumlo / *t;
    v[1] = c2 * *wnumhi / *t;
    if (v[0] > epsil && v[1] < vmax && (*wnumhi - *wnumlo) / *wnumhi < .01f) {
/*                          ** wavenumbers are very close.  get integral */
/*                          ** by iterating simpson rule to convergence. */
	hh = v[1] - v[0];
	oldval = 0.f;
/* Computing 3rd power */
	r__1 = v[0];
/* Computing 3rd power */
	r__2 = v[1];
	val0 = r__1 * (r__1 * r__1) / (exp(v[0]) - 1.f) + r__2 * (r__2 * r__2)
		 / (exp(v[1]) - 1.f);
	for (n = 1; n <= 10; ++n) {
	    del = hh / (n << 1);
	    val = val0;
	    i__1 = (n << 1) - 1;
	    for (k = 1; k <= i__1; ++k) {
		r__1 = v[0] + k * del;
/* Computing 3rd power */
		r__2 = r__1;
		val += (k % 2 + 1 << 1) * (r__2 * (r__2 * r__2) / (exp(r__1) 
			- 1.f));
/* L10: */
	    }
	    val = del / 3.f * val;
	    if ((r__1 = (val - oldval) / val, dabs(r__1)) <= 1e-6f) {
		goto L30;
	    }
	    oldval = val;
/* L20: */
	}
	errmsg_("plkavg--simpson rule didnt converge", &c_false, (ftnlen)35);
L30:
/* Computing 4th power */
	r__1 = *t, r__1 *= r__1;
	ret_val = sigdpi * (r__1 * r__1) * conc * val;
	return ret_val;
    }
/*                          *** general case *** */
    smallv = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	if (v[i__ - 1] < vcut) {
/*                                   ** use power series */
	    ++smallv;
/* Computing 2nd power */
	    r__1 = v[i__ - 1];
	    vsq = r__1 * r__1;
	    p[i__ - 1] = conc * vsq * v[i__ - 1] * (v[i__ - 1] * (v[i__ - 1] *
		     (vsq * (vsq * (vsq * -7.5156325156325161e-8f + 
		    3.6743092298647855e-6f) - 1.9841269841269841e-4f) + 
		    .016666666666666666f) - .125f) + .33333333333333331f);
	} else {
/*                      ** use exponential series */
	    mmax = 0;
/*                                ** find upper limit of series */
L40:
	    ++mmax;
	    if (v[i__ - 1] < vcp[mmax - 1]) {
		goto L40;
	    }
	    ex = exp(-v[i__ - 1]);
	    exm = 1.f;
	    d__[i__ - 1] = 0.f;
	    i__1 = mmax;
	    for (m = 1; m <= i__1; ++m) {
		mv = m * v[i__ - 1];
		exm = ex * exm;
/* Computing 4th power */
		i__2 = m, i__2 *= i__2;
		d__[i__ - 1] += exm * (mv * (mv * (mv + 3.f) + 6.f) + 6.f) / (
			i__2 * i__2);
/* L50: */
	    }
	    d__[i__ - 1] = conc * d__[i__ - 1];
	}
/* L60: */
    }
/*                              ** handle ill-conditioning */
    if (smallv == 2) {
/*                                    ** wnumlo and wnumhi both small */
	ret_val = p[1] - p[0];
    } else if (smallv == 1) {
/*                                    ** wnumlo small, wnumhi large */
	ret_val = 1.f - p[0] - d__[1];
    } else {
/*                                    ** wnumlo and wnumhi both large */
	ret_val = d__[0] - d__[1];
    }
/* Computing 4th power */
    r__1 = *t, r__1 *= r__1;
    ret_val = sigdpi * (r__1 * r__1) * ret_val;
/*      if( plkavg.eq.0.0 ) */
/*     &    call errmsg('plkavg--returns zero; possible underflow', */
/*     &    .false.) */
    return ret_val;
} /* plkavg_ */

