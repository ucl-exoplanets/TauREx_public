/* R1MACH.f -- translated by f2c (version 20100827).
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

static logical c_true = TRUE_;
static integer c__1 = 1;
static real c_b55 = 0.f;

doublereal r1mach_(integer *i__)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    real ret_val, r__1;
    static real equiv_4[6];

    /* Builtin functions */
    double r_lg10(real *);

    /* Local variables */
    static real s, eps;
#define log10 ((integer *)equiv_4 + 4)
#define large ((integer *)equiv_4 + 1)
#define rmach (equiv_4)
#define small ((integer *)equiv_4)
#define diver ((integer *)equiv_4 + 3)
#define right ((integer *)equiv_4 + 2)
    static real ratio;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static real epsnew;

/*  SINGLE-PRECISION MACHINE CONSTANTS */
/*  ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT, */
/*  BASE-B FORM */
/*     (SIGN) (B**E) * ( (X(1)/B) + ... + (X(T)/B**T) ) */
/*  WHERE 0.LE.X(I).LT.B, I=1,...,T, X(1).GT.0; EMIN.LE.E.LE.EMAX. */
/*  FOR EXAMPLE, IN BASE 10, THE X(I) WOULD JUST BE THE DIGITS */
/*  FOLLOWING THE DECIMAL POINT, AND CLEARLY THE FIRST DIGIT (X(1)) */
/*  WOULD HAVE TO BE NONZERO OR THE DECIMAL POINT COULD BE MOVED */
/*  OVER WITH A CORRESPONDING CHANGE IN 'E'.  THEN */
/*  R1MACH(1) = B**(EMIN-1), SMALLEST POSITIVE MAGNITUDE */
/*                 (E = EMIN, X(1) = 1, ALL OTHER X(I) = 0). */
/*  R1MACH(2) = B**EMAX*(1 - B**(-T)), LARGEST MAGNITUDE */
/*                 (E = EMAX, X(I) = B-1). */
/*  R1MACH(3) = 1/B**T, SMALLEST RELATIVE SPACING. */
/*  R1MACH(4) = 1/B**(T-1), MACHINE PRECISION (X(1)=X(T)=1,ALL OTHER */
/*                 X(I) = 0); SMALLEST POSITIVE EPS SUCH THAT 1+EPS.NE.1 */
/*  R1MACH(5) = LOG10(B) */
/*  REFERENCE: FOX P.A., HALL A.D., SCHRYER N.L.,'FRAMEWORK FOR A */
/*               PORTABLE LIBRARY', ACM TRANSACTIONS ON MATHEMATICAL */
/*               SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188. */
/*  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT, */
/*  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY */
/*  DELETING THE C FROM COLUMN 1. */
/*  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED. */
/*  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.) */
/*  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST */
/*  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE. */
/*  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED */
/*  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING */
/*  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD */
/*  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO */
/*  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER */
/*  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS. */
/* IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T */
/* 3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T */
/* PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300). */
/*      DATA SMALL(1)/8388608/, LARGE(1)/2139095039/, */
/*     $     RIGHT(1)/864026624/, DIVER(1)/872415232/, */
/*     $     LOG10(1)/ 1050288283/, SC/987/ */
/* AMDAHL MACHINES. */
/*      DATA SMALL(1)/1048576/, LARGE(1)/2147483647/, */
/*     $     RIGHT(1)/990904320/, DIVER(1)/1007681536/, */
/*     $     LOG10(1)/1091781651/, SC/987/ */
/* BURROUGHS 1700 SYSTEM. */
/*      DATA RMACH/Z400800000,Z5FFFFFFFF,Z4E9800000,Z4EA800000, */
/*     $             Z500E730E8/, SC/987/ */
/* BURROUGHS 5700/6700/7700 SYSTEMS. */
/*      DATA RMACH/O1771000000000000,O0777777777777777,O1311000000000000, */
/*     $             O1301000000000000,O1157163034761675/, SC/987/ */
/* FTN4 ON CDC 6000/7000 SERIES. */
/*      DATA RMACH/00564000000000000000B,37767777777777777776B, */
/*     $ 16414000000000000000B,16424000000000000000B, */
/*     $ 17164642023241175720B/, SC/987/ */
/* FTN5 ON CDC 6000/7000 SERIES. */
/*      DATA RMACH/O"00564000000000000000",O"37767777777777777776", */
/*     $ O"16414000000000000000",O"16424000000000000000", */
/*     $ O"17164642023241175720"/, SC/987/ */
/* CONVEX C-1. */
/*      DATA RMACH/'00800000'X,'7FFFFFFF'X,'34800000'X, */
/*     $ '35000000'X,'3F9A209B'X/, SC/987/ */
/* CRAY 1, XMP, 2, AND 3. */
/*      DATA RMACH/200034000000000000000B,577767777777777777776B, */
/*     $ 377224000000000000000B,377234000000000000000B, */
/*     $ 377774642023241175720B/, SC/987/ */
/* DATA GENERAL ECLIPSE S/200. */
/* NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE - */
/* STATIC RMACH(5) */
/*      DATA SMALL/20K,0/, LARGE/77777K,177777K/, RIGHT/35420K,0/, */
/*     $  DIVER/36020K,0/, LOG10/40423K,42023K/, SC/987/ */
/* HARRIS SLASH 6 AND SLASH 7. */
/*      DATA SMALL/'20000000,'00000201/, LARGE/'37777777,'00000177/, */
/*     $  RIGHT/'20000000,'00000352/, DIVER/'20000000,'00000353/, */
/*     $  LOG10/'23210115,'00000377/, SC/987/ */
/* HONEYWELL DPS 8/70 SERIES. */
/*      DATA RMACH/O402400000000,O376777777777,O714400000000, */
/*     $ O716400000000,O776464202324/, SC/987/ */
/* IBM 360/370 SERIES, */
/* XEROX SIGMA 5/7/9 AND SEL SYSTEMS 85/86. */
/*      DATA RMACH/Z00100000,Z7FFFFFFF,Z3B100000,Z3C100000, */
/*     $ Z41134413/, SC/987/ */
/* INTERDATA 8/32 WITH UNIX SYSTEM FORTRAN 77 COMPILER. */
/* FOR INTERDATA FORTRAN VII COMPILER REPLACE */
/* Z'S SPECIFYING HEX CONSTANTS WITH Y'S. */
/*      DATA RMACH/Z'00100000',Z'7EFFFFFF',Z'3B100000',Z'3C100000', */
/*     $ Z'41134413'/, SC/987/ */
/* PDP-10 (KA OR KI PROCESSOR). */
/*      DATA RMACH/"000400000000,"377777777777,"146400000000, */
/*     $ "147400000000,"177464202324/, SC/987/ */
/* PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS */
/* (EXPRESSED IN INTEGER AND OCTAL). */
/*      DATA SMALL(1)/8388608/, LARGE(1)/2147483647/, */
/*     $  RIGHT(1)/880803840/, DIVER(1)/889192448/, */
/*     $  LOG10(1)/1067065499/, SC/987/ */
/*      DATA RMACH/O00040000000,O17777777777,O06440000000, */
/*     $ O06500000000,O07746420233/, SC/987/ */
/* PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS */
/* (EXPRESSED IN INTEGER AND OCTAL). */
/*      DATA SMALL/128,0/, LARGE/32767,-1/, RIGHT/13440,0/, */
/*     $  DIVER/13568,0/, LOG10/16282,8347/, SC/987/ */
/*      DATA SMALL/O000200,O000000/, LARGE/O077777,O177777/, */
/*     $  RIGHT/O032200,O000000/, DIVER/O032400,O000000/, */
/*     $  LOG10/O037632,O020233/, SC/987/ */
/* SEQUENT BALANCE 8000. */
/*      DATA SMALL(1)/$00800000/, LARGE(1)/$7F7FFFFF/, */
/*     $  RIGHT(1)/$33800000/, DIVER(1)/$34000000/, */
/*     $  LOG10(1)/$3E9A209B/, SC/987/ */
/* UNIVAC 1100 SERIES. */
/*      DATA RMACH/O000400000000,O377777777777,O146400000000, */
/*     $ O147400000000,O177464202324/, SC/987/ */
/* VAX UNIX F77 COMPILER. */
/*      DATA SMALL(1)/128/, LARGE(1)/-32769/, RIGHT(1)/13440/, */
/*     $  DIVER(1)/13568/, LOG10(1)/547045274/, SC/987/ */
/* VAX-11 WITH FORTRAN IV-PLUS COMPILER. */
/*      DATA RMACH/Z00000080,ZFFFF7FFF,Z00003480,Z00003500, */
/*     $ Z209B3F9A/, SC/987/ */
/* VAX/VMS VERSION 2.2. */
/*      DATA RMACH/'80'X,'FFFF7FFF'X,'3480'X,'3500'X, */
/*     $ '209B3F9A'X/, SC/987/ */
    if (pass1) {
	pass1 = FALSE_;
/*         IF (SC.NE.987) */
/*     $       CALL ERRMSG( 'R1MACH--NO DATA STATEMENTS ACTIVE',.TRUE.) */
/*         sc = 987 */
	rmach[0] = slamch_("S", (ftnlen)1);
	rmach[1] = slamch_("O", (ftnlen)1);
	rmach[2] = slamch_("E", (ftnlen)1);
	rmach[3] = slamch_("P", (ftnlen)1);
	r__1 = slamch_("B", (ftnlen)1);
	rmach[4] = r_lg10(&r__1);
/*                      ** CALCULATE MACHINE PRECISION */
	epsnew = .01f;
L10:
	eps = epsnew;
	epsnew /= 1.1f;
/*                               ** IMPORTANT TO STORE 'S' SINCE MAY BE */
/*                               ** KEPT IN HIGHER PRECISION IN REGISTERS */
	s = epsnew + 1.f;
	if (s > 1.f) {
	    goto L10;
	}
	ratio = eps / rmach[3];
/* 	write(*,'(1a,1p5e15.7)') 's, epsnew,eps,rmach(4),eps/rmach(4)', */
/*     - s, epsnew,eps,rmach(4),eps/rmach(4) */
	if (ratio < .5f || ratio > 2.f) {
	    errmsg_("R1MACH--TABULATED PRECISION WRONG", &c_true, (ftnlen)33);
	}
    }
    if (*i__ < 1 || *i__ > 5) {
	errmsg_("R1MACH--ARGUMENT OUT OF BOUNDS", &c_true, (ftnlen)30);
    }
    ret_val = rmach[*i__ - 1];
    return ret_val;
} /* r1mach_ */

#undef right
#undef diver
#undef small
#undef rmach
#undef large
#undef log10


/* ..................... */
doublereal slamch_(char *cmach, ftnlen cmach_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Builtin functions */
    double pow_ri(real *, integer *);

    /* Local variables */
    static real t;
    static integer it;
    static real rnd, eps, base;
    static integer beta;
    static real emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static real rmin, rmax, rmach, small, sfmin;
    extern /* Subroutine */ int slamc2_(integer *, integer *, logical *, real 
	    *, integer *, real *, integer *, real *), errmsg_(char *, logical 
	    *, ftnlen);

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  SLAMCH determines single precision machine parameters. */
/*  Arguments */
/*  ========= */
/*  CMACH   (input) CHARACTER*1 */
/*          Specifies the value to be returned by SLAMCH: */
/*          = 'E' or 'e',   SLAMCH := eps */
/*          = 'S' or 's ,   SLAMCH := sfmin */
/*          = 'B' or 'b',   SLAMCH := base */
/*          = 'P' or 'p',   SLAMCH := eps*base */
/*          = 'N' or 'n',   SLAMCH := t */
/*          = 'R' or 'r',   SLAMCH := rnd */
/*          = 'M' or 'm',   SLAMCH := emin */
/*          = 'U' or 'u',   SLAMCH := rmin */
/*          = 'L' or 'l',   SLAMCH := emax */
/*          = 'O' or 'o',   SLAMCH := rmax */
/*          where */
/*          eps   = relative machine precision */
/*          sfmin = safe minimum, such that 1/sfmin does not overflow */
/*          base  = base of the machine */
/*          prec  = eps*base */
/*          t     = number of (base) digits in the mantissa */
/*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
/*          emin  = minimum exponent before (gradual) underflow */
/*          rmin  = underflow threshold - base**(emin-1) */
/*          emax  = largest exponent before overflow */
/*          rmax  = overflow threshold  - (base**emax)*(1-eps) */
/* ===================================================================== */
/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */
    if (first) {
	first = FALSE_;
	slamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (real) beta;
	t = (real) it;
	if (lrnd) {
	    rnd = 1.f;
	    i__1 = 1 - it;
	    eps = pow_ri(&base, &i__1) / 2;
	} else {
	    rnd = 0.f;
	    i__1 = 1 - it;
	    eps = pow_ri(&base, &i__1);
	}
	prec = eps * base;
	emin = (real) imin;
	emax = (real) imax;
	sfmin = rmin;
	small = 1.f / rmax;
	if (small >= sfmin) {
/*           Use SMALL plus a bit, to avoid the possibility of rounding */
/*           causing overflow when computing  1/sfmin. */
	    sfmin = small * (eps + 1.f);
	}
    }
    if (*(unsigned char *)cmach == 'E' || *(unsigned char *)cmach == 'e') {
	rmach = eps;
    } else if (*(unsigned char *)cmach == 'S' || *(unsigned char *)cmach == 
	    's') {
	rmach = sfmin;
    } else if (*(unsigned char *)cmach == 'B' || *(unsigned char *)cmach == 
	    'b') {
	rmach = base;
    } else if (*(unsigned char *)cmach == 'P' || *(unsigned char *)cmach == 
	    'p') {
	rmach = prec;
    } else if (*(unsigned char *)cmach == 'U' || *(unsigned char *)cmach == 
	    'u') {
	rmach = rmin;
    } else if (*(unsigned char *)cmach == 'O' || *(unsigned char *)cmach == 
	    'o') {
	rmach = rmax;
    } else if (*(unsigned char *)cmach == 'N' || *(unsigned char *)cmach == 
	    'n') {
	rmach = t;
    } else if (*(unsigned char *)cmach == 'R' || *(unsigned char *)cmach == 
	    'r') {
	rmach = rnd;
    } else if (*(unsigned char *)cmach == 'M' || *(unsigned char *)cmach == 
	    'm') {
	rmach = emin;
    } else if (*(unsigned char *)cmach == 'L' || *(unsigned char *)cmach == 
	    'l') {
	rmach = emax;
    } else {
	errmsg_("SLAMCH--invalid argument", &c_true, (ftnlen)24);
    }
    ret_val = rmach;
    return ret_val;
/*     End of SLAMCH */
} /* slamch_ */

/* *********************************************************************** */
/* Subroutine */ int slamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static real a, b, c__, f, t1, t2;
    static integer lt;
    static real one, qtr;
    static logical lrnd;
    static integer lbeta;
    static real savec;
    static logical lieee1;
    extern doublereal slamc3_(real *, real *);

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  SLAMC1 determines the machine parameters given by BETA, T, RND, and */
/*  IEEE1. */
/*  Arguments */
/*  ========= */
/*  BETA    (output) INTEGER */
/*          The base of the machine. */
/*  T       (output) INTEGER */
/*          The number of ( BETA ) digits in the mantissa. */
/*  RND     (output) LOGICAL */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/*          be a reliable guide to the way in which the machine performs */
/*          its arithmetic. */
/*  IEEE1   (output) LOGICAL */
/*          Specifies whether rounding appears to be done in the IEEE */
/*          'round to nearest' style. */
/*  Further Details */
/*  =============== */
/*  The routine is based on the routine  ENVRON  by Malcolm and */
/*  incorporates suggestions by Gentleman and Marovich. See */
/*     Malcolm M. A. (1972) Algorithms to reveal properties of */
/*        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */
/*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
/*        that reveal properties of floating point arithmetic units. */
/*        Comms. of the ACM, 17, 276-277. */
/* ===================================================================== */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */
    if (first) {
	first = FALSE_;
	one = 1.f;
/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, */
/*        IEEE1, T and RND. */
/*        Throughout this routine  we use the function  SLAMC3  to ensure */
/*        that relevant values are  stored and not held in registers,  or */
/*        are not affected by optimizers. */
/*        Compute  a = 2.0**m  with the  smallest positive integer m such */
/*        that */
/*           fl( a + 1.0 ) = a. */
	a = 1.f;
	c__ = 1.f;
L10:
	if (c__ == one) {
	    a *= 2;
	    c__ = slamc3_(&a, &one);
	    r__1 = -a;
	    c__ = slamc3_(&c__, &r__1);
	    goto L10;
	}
/*        Now compute  b = 2.0**m  with the smallest positive integer m */
/*        such that */
/*           fl( a + b ) .gt. a. */
	b = 1.f;
	c__ = slamc3_(&a, &b);
L20:
	if (c__ == a) {
	    b *= 2;
	    c__ = slamc3_(&a, &b);
	    goto L20;
	}
/*        Now compute the base.  a and c  are neighbouring floating point */
/*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so */
/*        their difference is beta. Adding 0.25 to c is to ensure that it */
/*        is truncated to beta and not ( beta - 1 ). */
	qtr = one / 4;
	savec = c__;
	r__1 = -a;
	c__ = slamc3_(&c__, &r__1);
	lbeta = c__ + qtr;
/*        Now determine whether rounding or chopping occurs,  by adding a */
/*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */
	b = (real) lbeta;
	r__1 = b / 2;
	r__2 = -b / 100;
	f = slamc3_(&r__1, &r__2);
	c__ = slamc3_(&f, &a);
	if (c__ == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	r__1 = b / 2;
	r__2 = b / 100;
	f = slamc3_(&r__1, &r__2);
	c__ = slamc3_(&f, &a);
	if (lrnd && c__ == a) {
	    lrnd = FALSE_;
	}
/*        Try and decide whether rounding is done in the  IEEE  'round to */
/*        nearest' style. B/2 is half a unit in the last place of the two */
/*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit */
/*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change */
/*        A, but adding B/2 to SAVEC should change SAVEC. */
	r__1 = b / 2;
	t1 = slamc3_(&r__1, &a);
	r__1 = b / 2;
	t2 = slamc3_(&r__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;
/*        Now find  the  mantissa, t.  It should  be the  integer part of */
/*        log to the base beta of a,  however it is safer to determine  t */
/*        by powering.  So we find t as the smallest positive integer for */
/*        which */
/*           fl( beta**t + 1.0 ) = 1.0. */
	lt = 0;
	a = 1.f;
	c__ = 1.f;
L30:
	if (c__ == one) {
	    ++lt;
	    a *= lbeta;
	    c__ = slamc3_(&a, &one);
	    r__1 = -a;
	    c__ = slamc3_(&c__, &r__1);
	    goto L30;
	}
    }
    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;
/*     End of SLAMC1 */
} /* slamc1_ */

/* *********************************************************************** */
/* Subroutine */ int slamc2_(integer *beta, integer *t, logical *rnd, real *
	eps, integer *emin, real *rmin, integer *emax, real *rmax)
{
    /* Initialized data */

    static logical first = TRUE_;
    static logical iwarn = FALSE_;

    /* Format strings */
    static char fmt_9999[] = "(//\002 WARNING. The value EMIN may be incorre"
	    "ct:-\002,\002  EMIN = \002,i8,/\002 If, after inspection, the va"
	    "lue EMIN looks\002,\002 acceptable please comment out \002,/\002"
	    " the IF block as marked within the code of routine\002,\002 SLAM"
	    "C2,\002,/\002 otherwise supply EMIN explicitly.\002,/)";

    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double pow_ri(real *, integer *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static real a, b, c__;
    static integer i__, lt;
    static real one, two;
    static logical ieee;
    static real half;
    static logical lrnd;
    static real leps, zero;
    static integer lbeta;
    static real rbase;
    static integer lemin, lemax, gnmin;
    static real small;
    static integer gpmin;
    static real third, lrmin, lrmax, sixth;
    static logical lieee1;
    extern /* Subroutine */ int slamc1_(integer *, integer *, logical *, 
	    logical *);
    extern doublereal slamc3_(real *, real *);
    extern /* Subroutine */ int slamc4_(integer *, real *, integer *), 
	    slamc5_(integer *, integer *, integer *, logical *, integer *, 
	    real *);
    static integer ngnmin, ngpmin;

    /* Fortran I/O blocks */
    static cilist io___69 = { 0, 6, 0, fmt_9999, 0 };


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  SLAMC2 determines the machine parameters specified in its argument */
/*  list. */
/*  Arguments */
/*  ========= */
/*  BETA    (output) INTEGER */
/*          The base of the machine. */
/*  T       (output) INTEGER */
/*          The number of ( BETA ) digits in the mantissa. */
/*  RND     (output) LOGICAL */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/*          be a reliable guide to the way in which the machine performs */
/*          its arithmetic. */
/*  EPS     (output) REAL */
/*          The smallest positive number such that */
/*             fl( 1.0 - EPS ) .LT. 1.0, */
/*          where fl denotes the computed value. */
/*  EMIN    (output) INTEGER */
/*          The minimum exponent before (gradual) underflow occurs. */
/*  RMIN    (output) REAL */
/*          The smallest normalized number for the machine, given by */
/*          BASE**( EMIN - 1 ), where  BASE  is the floating point value */
/*          of BETA. */
/*  EMAX    (output) INTEGER */
/*          The maximum exponent before overflow occurs. */
/*  RMAX    (output) REAL */
/*          The largest positive number for the machine, given by */
/*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point */
/*          value of BETA. */
/*  Further Details */
/*  =============== */
/*  The computation of  EPS  is based on a routine PARANOIA by */
/*  W. Kahan of the University of California at Berkeley. */
/* ===================================================================== */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */
    if (first) {
	first = FALSE_;
	zero = 0.f;
	one = 1.f;
	two = 2.f;
/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of */
/*        BETA, T, RND, EPS, EMIN and RMIN. */
/*        Throughout this routine  we use the function  SLAMC3  to ensure */
/*        that relevant values are stored  and not held in registers,  or */
/*        are not affected by optimizers. */
/*        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */
	slamc1_(&lbeta, &lt, &lrnd, &lieee1);
/*        Start to find EPS. */
	b = (real) lbeta;
	i__1 = -lt;
	a = pow_ri(&b, &i__1);
	leps = a;
/*        Try some tricks to see whether or not this is the correct  EPS. */
	b = two / 3;
	half = one / 2;
	r__1 = -half;
	sixth = slamc3_(&b, &r__1);
	third = slamc3_(&sixth, &sixth);
	r__1 = -half;
	b = slamc3_(&third, &r__1);
	b = slamc3_(&b, &sixth);
	b = dabs(b);
	if (b < leps) {
	    b = leps;
	}
	leps = 1.f;
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    r__1 = half * leps;
/* Computing 5th power */
	    r__3 = two, r__4 = r__3, r__3 *= r__3;
/* Computing 2nd power */
	    r__5 = leps;
	    r__2 = r__4 * (r__3 * r__3) * (r__5 * r__5);
	    c__ = slamc3_(&r__1, &r__2);
	    r__1 = -c__;
	    c__ = slamc3_(&half, &r__1);
	    b = slamc3_(&half, &c__);
	    r__1 = -b;
	    c__ = slamc3_(&half, &r__1);
	    b = slamc3_(&half, &c__);
	    goto L10;
	}
	if (a < leps) {
	    leps = a;
	}
/*      write(*,*) ' EPS from PARANOIA:', LEPS */
/*        Computation of EPS complete. */
/*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)). */
/*        Keep dividing  A by BETA until (gradual) underflow occurs. This */
/*        is detected when we cannot recover the previous A. */
	rbase = one / lbeta;
	small = one;
	for (i__ = 1; i__ <= 3; ++i__) {
	    r__1 = small * rbase;
	    small = slamc3_(&r__1, &zero);
/* L20: */
	}
	a = slamc3_(&one, &small);
	slamc4_(&ngpmin, &one, &lbeta);
	r__1 = -one;
	slamc4_(&ngnmin, &r__1, &lbeta);
	slamc4_(&gpmin, &a, &lbeta);
	r__1 = -a;
	slamc4_(&gnmin, &r__1, &lbeta);
	ieee = FALSE_;
	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual underflow; */
/*              e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual underflow; */
/*              e.g., IEEE standard followers ) */
	    } else {
		lemin = min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }
	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
		lemin = max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow; */
/*              e.g., CYBER 205 ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }
	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - min(ngpmin,ngnmin) == 3) {
		lemin = max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflow; */
/*              no known machine ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }
	} else {
/* Computing MIN */
	    i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
	    lemin = min(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = TRUE_;
	}
/* ** */
/* Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = TRUE_;
	    s_wsfe(&io___69);
	    do_fio(&c__1, (char *)&lemin, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
/* ** */
/*        Assume IEEE arithmetic if we found denormalised  numbers above, */
/*        or if arithmetic seems to round in the  IEEE style,  determined */
/*        in routine SLAMC1. A true IEEE machine should have both  things */
/*        true; however, faulty machines may have one or the other. */
	ieee = ieee || lieee1;
/*        Compute  RMIN by successive division by  BETA. We could compute */
/*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during */
/*        this computation. */
	lrmin = 1.f;
	i__1 = 1 - lemin;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__1 = lrmin * rbase;
	    lrmin = slamc3_(&r__1, &zero);
/* L30: */
	}
/*        Finally, call SLAMC5 to compute EMAX and RMAX. */
	slamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }
    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;
    return 0;
/*     End of SLAMC2 */
} /* slamc2_ */

/* *********************************************************************** */
doublereal slamc3_(real *a, real *b)
{
    /* System generated locals */
    real ret_val;

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  SLAMC3  is intended to force  A  and  B  to be stored prior to doing */
/*  the addition of  A  and  B ,  for use in situations where optimizers */
/*  might hold one of these in a register. */
/*  Arguments */
/*  ========= */
/*  A, B    (input) REAL */
/*          The values A and B. */
/* ===================================================================== */
/*     .. Executable Statements .. */
    ret_val = *a + *b;
    return ret_val;
/*     End of SLAMC3 */
} /* slamc3_ */

/* *********************************************************************** */
/* Subroutine */ int slamc4_(integer *emin, real *start, integer *base)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static real a;
    static integer i__;
    static real b1, b2, c1, c2, d1, d2, one, zero, rbase;
    extern doublereal slamc3_(real *, real *);

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  SLAMC4 is a service routine for SLAMC2. */
/*  Arguments */
/*  ========= */
/*  EMIN    (output) EMIN */
/*          The minimum exponent before (gradual) underflow, computed by */
/*          setting A = START and dividing by BASE until the previous A */
/*          can not be recovered. */
/*  START   (input) REAL */
/*          The starting point for determining EMIN. */
/*  BASE    (input) INTEGER */
/*          The base of the machine. */
/* ===================================================================== */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */
    a = *start;
    one = 1.f;
    rbase = one / *base;
    zero = 0.f;
    *emin = 1;
    r__1 = a * rbase;
    b1 = slamc3_(&r__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	r__1 = a / *base;
	b1 = slamc3_(&r__1, &zero);
	r__1 = b1 * *base;
	c1 = slamc3_(&r__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d1 += b1;
/* L20: */
	}
	r__1 = a * rbase;
	b2 = slamc3_(&r__1, &zero);
	r__1 = b2 / rbase;
	c2 = slamc3_(&r__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
    return 0;
/*     End of SLAMC4 */
} /* slamc4_ */

/* *********************************************************************** */
/* Subroutine */ int slamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, real *rmax)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__;
    static real y, z__;
    static integer try__, lexp;
    static real oldy;
    static integer uexp, nbits;
    extern doublereal slamc3_(real *, real *);
    static real recbas;
    static integer exbits, expsum;

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  SLAMC5 attempts to compute RMAX, the largest machine floating-point */
/*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum */
/*  approximately to a power of 2.  It will fail on machines where this */
/*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625, */
/*  EMAX = 28718).  It will also fail if the value supplied for EMIN is */
/*  too large (i.e. too close to zero), probably with overflow. */
/*  Arguments */
/*  ========= */
/*  BETA    (input) INTEGER */
/*          The base of floating-point arithmetic. */
/*  P       (input) INTEGER */
/*          The number of base BETA digits in the mantissa of a */
/*          floating-point value. */
/*  EMIN    (input) INTEGER */
/*          The minimum exponent before (gradual) underflow. */
/*  IEEE    (input) LOGICAL */
/*          A logical flag specifying whether or not the arithmetic */
/*          system is thought to comply with the IEEE standard. */
/*  EMAX    (output) INTEGER */
/*          The largest exponent before overflow */
/*  RMAX    (output) REAL */
/*          The largest machine floating-point number. */
/* ===================================================================== */
/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
/*     First compute LEXP and UEXP, two powers of 2 that bound */
/*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum */
/*     approximately to the bound that is closest to abs(EMIN). */
/*     (EMAX is the exponent of the required number RMAX). */
    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }
/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
/*     than or equal to EMIN. EXBITS is the number of bits needed to */
/*     store the exponent. */
    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }
/*     EXPSUM is the exponent range, approximately equal to */
/*     EMAX - EMIN + 1 . */
    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;
/*     NBITS is the total number of bits needed to store a */
/*     floating-point number. */
    if (nbits % 2 == 1 && *beta == 2) {
/*        Either there are an odd number of bits used to store a */
/*        floating-point number, which is unlikely, or some bits are */
/*        not used in the representation of numbers, which is possible, */
/*        (e.g. Cray machines) or the mantissa has an implicit bit, */
/*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the */
/*        most likely. We have to assume the last alternative. */
/*        If this is true, then we need to reduce EMAX by one because */
/*        there must be some way of representing zero in an implicit-bit */
/*        system. On machines like Cray, we are reducing EMAX by one */
/*        unnecessarily. */
	--(*emax);
    }
    if (*ieee) {
/*        Assume we are on an IEEE machine which reserves one exponent */
/*        for infinity and NaN. */
	--(*emax);
    }
/*     Now create RMAX, the largest machine number, which should */
/*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */
/*     First compute 1.0 - BETA**(-P), being careful that the */
/*     result is less than 1.0 . */
    recbas = 1.f / *beta;
    z__ = *beta - 1.f;
    y = 0.f;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ *= recbas;
	if (y < 1.f) {
	    oldy = y;
	}
	y = slamc3_(&y, &z__);
/* L20: */
    }
    if (y >= 1.f) {
	y = oldy;
    }
/*     Now multiply by BETA**EMAX to get RMAX. */
    i__1 = *emax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__1 = y * *beta;
	y = slamc3_(&r__1, &c_b55);
/* L30: */
    }
    *rmax = y;
    return 0;
/*     End of SLAMC5 */
} /* slamc5_ */

