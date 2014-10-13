/* D1MACH.f -- translated by f2c (version 20100827).
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
static doublereal c_b55 = 0.;

doublereal d1mach_(integer *i__)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    doublereal ret_val, d__1;
    static doublereal equiv_4[6];

    /* Builtin functions */
    double d_lg10(doublereal *);

    /* Local variables */
    static doublereal s, eps;
#define log10 ((integer *)equiv_4 + 8)
#define dmach (equiv_4)
#define large ((integer *)equiv_4 + 2)
#define small ((integer *)equiv_4)
#define diver ((integer *)equiv_4 + 6)
#define right ((integer *)equiv_4 + 4)
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static doublereal epsnew;

/*  DOUBLE-PRECISION MACHINE CONSTANTS (SEE R1MACH FOR DOCUMENTATION) */
/*  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST */
/*  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE. */
/* IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T 3B SERIES AND */
/* MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T PC 7300), */
/* IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST. */
/*      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2146435071,-1/, */
/*     $  (RIGHT(N),N=1,2)/1017118720,0/, (DIVER(N),N=1,2)/1018167296,0/, */
/*     $  (LOG10(N),N=1,2)/1070810131,1352628735/, SC/987/ */
/* IEEE ARITHMETIC MACHINES AND 8087-BASED MICROS, SUCH AS THE IBM PC */
/* AND AT&T 6300, IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST. */
/*      DATA (SMALL(N),N=1,2)/0,1048576/, (LARGE(N),N=1,2)/-1,2146435071/, */
/*     $  (RIGHT(N),N=1,2)/0,1017118720/, (DIVER(N),N=1,2)/0,1018167296/, */
/*     $  (LOG10(N),N=1,2)/1352628735,1070810131/, SC/987/ */
/* AMDAHL MACHINES. */
/*      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2147483647,-1/, */
/*     $ (RIGHT(N),N=1,2)/856686592,0/, (DIVER(N),N=1,2)/ 873463808,0/, */
/*     $ (LOG10(N),N=1,2)/1091781651,1352628735/, SC/987/ */
/* BURROUGHS 1700 SYSTEM. */
/*      DATA (SMALL(N),N=1,2)/ZC00800000,Z000000000/, */
/*     $ (LARGE(N),N=1,2)/ZDFFFFFFFF,ZFFFFFFFFF/, */
/*     $ (RIGHT(N),N=1,2)/ZCC5800000,Z000000000/, */
/*     $ (DIVER(N),N=1,2)/ZCC6800000,Z000000000/, */
/*     $ (LOG10(N),N=1,2)/ZD00E730E7,ZC77800DC0/, SC/987/ */
/* BURROUGHS 5700 SYSTEM. */
/*      DATA (SMALL(N),N=1,2)/O1771000000000000,O0000000000000000/, */
/*     $  (LARGE(N),N=1,2)/O0777777777777777,O0007777777777777/, */
/*     $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/, */
/*     $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/, */
/*     $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/ */
/* BURROUGHS 6700/7700 SYSTEMS. */
/*      DATA (SMALL(N),N=1,2)/O1771000000000000,O7770000000000000/, */
/*     $  (LARGE(N),N=1,2)/O0777777777777777,O7777777777777777/, */
/*     $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/, */
/*     $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/, */
/*     $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/ */
/* FTN4 ON THE CDC 6000/7000 SERIES. */
/*      DATA */
/*     $  (SMALL(N),N=1,2)/00564000000000000000B,00000000000000000000B/, */
/*     $  (LARGE(N),N=1,2)/37757777777777777777B,37157777777777777774B/, */
/*     $  (RIGHT(N),N=1,2)/15624000000000000000B,00000000000000000000B/, */
/*     $  (DIVER(N),N=1,2)/15634000000000000000B,00000000000000000000B/, */
/*     $  (LOG10(N),N=1,2)/17164642023241175717B,16367571421742254654B/, */
/*     $  SC/987/ */
/* FTN5 ON THE CDC 6000/7000 SERIES. */
/*      DATA */
/*     $(SMALL(N),N=1,2)/O"00564000000000000000",O"00000000000000000000"/, */
/*     $(LARGE(N),N=1,2)/O"37757777777777777777",O"37157777777777777774"/, */
/*     $(RIGHT(N),N=1,2)/O"15624000000000000000",O"00000000000000000000"/, */
/*     $(DIVER(N),N=1,2)/O"15634000000000000000",O"00000000000000000000"/, */
/*     $(LOG10(N),N=1,2)/O"17164642023241175717",O"16367571421742254654"/, */
/*     $ SC/987/ */
/* CONVEX C-1 */
/*      DATA (SMALL(N),N=1,2)/'00100000'X,'00000000'X/, */
/*     $  (LARGE(N),N=1,2)/'7FFFFFFF'X,'FFFFFFFF'X/, */
/*     $  (RIGHT(N),N=1,2)/'3CC00000'X,'00000000'X/, */
/*     $  (DIVER(N),N=1,2)/'3CD00000'X,'00000000'X/, */
/*     $  (LOG10(N),N=1,2)/'3FF34413'X,'509F79FF'X/, SC/987/ */
/* CRAY 1, XMP, 2, AND 3. */
/*      DATA */
/*     $ (SMALL(N),N=1,2)/201354000000000000000B,000000000000000000000B/, */
/*     $ (LARGE(N),N=1,2)/577767777777777777777B,000007777777777777776B/, */
/*     $ (RIGHT(N),N=1,2)/376434000000000000000B,000000000000000000000B/, */
/*     $ (DIVER(N),N=1,2)/376444000000000000000B,000000000000000000000B/, */
/*     $ (LOG10(N),N=1,2)/377774642023241175717B,000007571421742254654B/, */
/*     $ SC/987/ */
/* DATA GENERAL ECLIPSE S/200 */
/* NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE - */
/* STATIC DMACH(5) */
/*      DATA SMALL/20K,3*0/, LARGE/77777K,3*177777K/, */
/*     $  RIGHT/31420K,3*0/, DIVER/32020K,3*0/, */
/*     $  LOG10/40423K,42023K,50237K,74776K/, SC/987/ */
/* HARRIS SLASH 6 AND SLASH 7 */
/*      DATA (SMALL(N),N=1,2)/'20000000,'00000201/, */
/*     $  (LARGE(N),N=1,2)/'37777777,'37777577/, */
/*     $  (RIGHT(N),N=1,2)/'20000000,'00000333/, */
/*     $  (DIVER(N),N=1,2)/'20000000,'00000334/, */
/*     $  (LOG10(N),N=1,2)/'23210115,'10237777/, SC/987/ */
/* HONEYWELL DPS 8/70 SERIES. */
/*      DATA (SMALL(N),N=1,2)/O402400000000,O000000000000/, */
/*     $  (LARGE(N),N=1,2)/O376777777777,O777777777777/, */
/*     $  (RIGHT(N),N=1,2)/O604400000000,O000000000000/, */
/*     $  (DIVER(N),N=1,2)/O606400000000,O000000000000/, */
/*     $  (LOG10(N),N=1,2)/O776464202324,O117571775714/, SC/987/ */
/* IBM 360/370 SERIES, XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86. */
/*      DATA (SMALL(N),N=1,2)/Z00100000,Z00000000/, */
/*     $  (LARGE(N),N=1,2)/Z7FFFFFFF,ZFFFFFFFF/, */
/*     $  (RIGHT(N),N=1,2)/Z33100000,Z00000000/, */
/*     $  (DIVER(N),N=1,2)/Z34100000,Z00000000/, */
/*     $  (LOG10(N),N=1,2)/Z41134413,Z509F79FF/, SC/987/ */
/* INTERDATA 8/32 WITH THE UNIX SYSTEM FORTRAN 77 COMPILER. */
/* FOR THE INTERDATA FORTRAN VII COMPILER REPLACE */
/* THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S. */
/*      DATA (SMALL(N),N=1,2)/Z'00100000',Z'00000000'/, */
/*     $  (LARGE(N),N=1,2)/Z'7EFFFFFF',Z'FFFFFFFF'/, */
/*     $  (RIGHT(N),N=1,2)/Z'33100000',Z'00000000'/, */
/*     $  (DIVER(N),N=1,2)/Z'34100000',Z'00000000'/, */
/*     $  (LOG10(N),N=1,2)/Z'41134413',Z'509F79FF'/, SC/987/ */
/* PDP-10 (KA PROCESSOR). */
/*      DATA (SMALL(N),N=1,2)/"033400000000,"000000000000/, */
/*     $  (LARGE(N),N=1,2)/"377777777777,"344777777777/, */
/*     $  (RIGHT(N),N=1,2)/"113400000000,"000000000000/, */
/*     $  (DIVER(N),N=1,2)/"114400000000,"000000000000/, */
/*     $  (LOG10(N),N=1,2)/"177464202324,"144117571776/, SC/987/ */
/* PDP-10 (KI PROCESSOR). */
/*      DATA (SMALL(N),N=1,2)/"000400000000,"000000000000/, */
/*     $  (LARGE(N),N=1,2)/"377777777777,"377777777777/, */
/*     $  (RIGHT(N),N=1,2)/"103400000000,"000000000000/, */
/*     $  (DIVER(N),N=1,2)/"104400000000,"000000000000/, */
/*     $  (LOG10(N),N=1,2)/"177464202324,"047674776746/, SC/987/ */
/* PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS */
/* (EXPRESSED IN INTEGER AND OCTAL). */
/*      DATA (SMALL(N),N=1,2)/8388608,0/, (LARGE(N),N=1,2)/2147483647,-1/, */
/*     $  (RIGHT(N),N=1,2)/612368384,0/, (DIVER(N),N=1,2)/620756992,0/, */
/*     $  (LOG10(N),N=1,2)/1067065498,-2063872008/, SC/987/ */
/*      DATA (SMALL(N),N=1,2)/O00040000000,O00000000000/, */
/*     $  (LARGE(N),N=1,2)/O17777777777,O37777777777/, */
/*     $  (RIGHT(N),N=1,2)/O04440000000,O00000000000/, */
/*     $  (DIVER(N),N=1,2)/O04500000000,O00000000000/, */
/*     $  (LOG10(N),N=1,2)/O07746420232,O20476747770/, SC/987/ */
/* PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS */
/* (EXPRESSED IN INTEGER AND OCTAL). */
/*      DATA SMALL/128,3*0/, LARGE/32767,3*-1/, RIGHT/9344,3*0/, */
/*     $  DIVER/9472,3*0/, LOG10/16282,8346,-31493,-12296/, SC/987/ */
/*      DATA SMALL/O000200,3*O000000/, LARGE/O077777,3*O177777/, */
/*     $  RIGHT/O022200,3*O000000/, DIVER/O022400,3*O000000/, */
/*     $  LOG10/O037632,O020232,O102373,O147770/, SC/987/ */
/* PRIME 50 SERIES SYSTEMS WITH 32-BIT INTEGERS AND 64V MODE */
/* INSTRUCTIONS, SUPPLIED BY IGOR BRAY. */
/*      DATA (SMALL(N),N=1,2)/:10000000000,:00000100001/, */
/*     $  (LARGE(N),N=1,2)/:17777777777,:37777677775/, */
/*     $  (RIGHT(N),N=1,2)/:10000000000,:00000000122/, */
/*     $  (DIVER(N),N=1,2)/:10000000000,:00000000123/, */
/*     $  (LOG10(N),N=1,2)/:11504046501,:07674600177/, SC/987/ */
/* SEQUENT BALANCE 8000 */
/*      DATA (SMALL(N),N=1,2)/$00000000, $00100000/, */
/*     $  (LARGE(N),N=1,2)/$FFFFFFFF, $7FEFFFFF/, */
/*     $  (RIGHT(N),N=1,2)/$00000000, $3CA00000/, */
/*     $  (DIVER(N),N=1,2)/$00000000, $3CB00000/, */
/*     $  (LOG10(N),N=1,2)/$509F79FF, $3FD34413/, SC/987/ */
/* UNIVAC 1100 SERIES. */
/*      DATA (SMALL(N),N=1,2)/O000040000000,O000000000000/, */
/*     $  (LARGE(N),N=1,2)/O377777777777,O777777777777/, */
/*     $  (RIGHT(N),N=1,2)/O170540000000,O000000000000/, */
/*     $  (DIVER(N),N=1,2)/O170640000000,O000000000000/, */
/*     $  (LOG10(N),N=1,2)/O177746420232,O411757177572/, SC/987/ */
/* VAX UNIX F77 COMPILER */
/*      DATA (SMALL(N),N=1,2)/128,0/, (LARGE(N),N=1,2)/-32769,-1/, */
/*     $  (RIGHT(N),N=1,2)/9344,0/, (DIVER(N),N=1,2)/9472,0/, */
/*     $  (LOG10(N),N=1,2)/546979738,-805796613/, SC/987/ */
/* VAX-11 WITH FORTRAN IV-PLUS COMPILER */
/*      DATA (SMALL(N),N=1,2)/Z00000080,Z00000000/, */
/*     $  (LARGE(N),N=1,2)/ZFFFF7FFF,ZFFFFFFFF/, */
/*     $  (RIGHT(N),N=1,2)/Z00002480,Z00000000/, */
/*     $  (DIVER(N),N=1,2)/Z00002500,Z00000000/, */
/*     $  (LOG10(N),N=1,2)/Z209A3F9A,ZCFF884FB/, SC/987/ */
/* VAX/VMS VERSION 2.2 */
/*      DATA (SMALL(N),N=1,2)/'80'X,'0'X/, */
/*     $  (LARGE(N),N=1,2)/'FFFF7FFF'X,'FFFFFFFF'X/, */
/*     $  (RIGHT(N),N=1,2)/'2480'X,'0'X/, (DIVER(N),N=1,2)/'2500'X,'0'X/, */
/*     $  (LOG10(N),N=1,2)/'209A3F9A'X,'CFF884FB'X/, SC/987/ */
    if (pass1) {
	pass1 = FALSE_;
/*         IF (SC.NE.987) */
/*     $       CALL ERRMSG( 'D1MACH--NO DATA STATEMENTS ACTIVE',.TRUE.) */
/*         sc = 987 */
	dmach[0] = dlamch_("S", (ftnlen)1);
	dmach[1] = dlamch_("O", (ftnlen)1);
	dmach[2] = dlamch_("E", (ftnlen)1);
	dmach[3] = dlamch_("P", (ftnlen)1);
	d__1 = dlamch_("B", (ftnlen)1);
	dmach[4] = d_lg10(&d__1);
/*                        ** CALCULATE MACHINE PRECISION */
	epsnew = .01;
L10:
	eps = epsnew;
	epsnew /= 1.1;
/*                               ** IMPORTANT TO STORE 'S' SINCE MAY BE */
/*                               ** KEPT IN HIGHER PRECISION IN REGISTERS */
	s = epsnew + 1.;
/* 	write(*,'(1a,1p5e15.7)') 's, epsnew,eps,dmach(4),eps/dmach(4)', */
/*     - s, epsnew,eps,dmach(4),eps/dmach(4) */
	if (s > 1.) {
	    goto L10;
	}
	if (eps / dmach[3] < .5 || eps / dmach[3] > 2.) {
	    errmsg_("D1MACH--TABULATED PRECISION WRONG", &c_true, (ftnlen)33);
	}
    }
    if (*i__ < 1 || *i__ > 5) {
	errmsg_("D1MACH--ARGUMENT OUT OF BOUNDS", &c_true, (ftnlen)30);
    }
    ret_val = dmach[*i__ - 1];
    return ret_val;
} /* d1mach_ */

#undef right
#undef diver
#undef small
#undef large
#undef dmach
#undef log10


/* .................................... */
doublereal dlamch_(char *cmach, ftnlen cmach_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal t;
    static integer it;
    static doublereal rnd, eps, base;
    static integer beta;
    static doublereal emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static doublereal rmin, rmax, rmach, small, sfmin;
    extern /* Subroutine */ int dlamc2_(integer *, integer *, logical *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *), 
	    errmsg_(char *, logical *, ftnlen);

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  DLAMCH determines double precision machine parameters. */
/*  Arguments */
/*  ========= */
/*  CMACH   (input) CHARACTER*1 */
/*          Specifies the value to be returned by DLAMCH: */
/*          = 'E' or 'e',   DLAMCH := eps */
/*          = 'S' or 's ,   DLAMCH := sfmin */
/*          = 'B' or 'b',   DLAMCH := base */
/*          = 'P' or 'p',   DLAMCH := eps*base */
/*          = 'N' or 'n',   DLAMCH := t */
/*          = 'R' or 'r',   DLAMCH := rnd */
/*          = 'M' or 'm',   DLAMCH := emin */
/*          = 'U' or 'u',   DLAMCH := rmin */
/*          = 'L' or 'l',   DLAMCH := emax */
/*          = 'O' or 'o',   DLAMCH := rmax */
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
	dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (doublereal) beta;
	t = (doublereal) it;
	if (lrnd) {
	    rnd = 1.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1) / 2;
	} else {
	    rnd = 0.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1);
	}
	prec = eps * base;
	emin = (doublereal) imin;
	emax = (doublereal) imax;
	sfmin = rmin;
	small = 1. / rmax;
	if (small >= sfmin) {
/*           Use SMALL plus a bit, to avoid the possibility of rounding */
/*           causing overflow when computing  1/sfmin. */
	    sfmin = small * (eps + 1.);
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
/*     End of DLAMCH */
} /* dlamch_ */

/* *********************************************************************** */
/* Subroutine */ int dlamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a, b, c__, f, t1, t2;
    static integer lt;
    static doublereal one, qtr;
    static logical lrnd;
    static integer lbeta;
    static doublereal savec;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static logical lieee1;

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  DLAMC1 determines the machine parameters given by BETA, T, RND, and */
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
	one = 1.;
/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, */
/*        IEEE1, T and RND. */
/*        Throughout this routine  we use the function  DLAMC3  to ensure */
/*        that relevant values are  stored and not held in registers,  or */
/*        are not affected by optimizers. */
/*        Compute  a = 2.0**m  with the  smallest positive integer m such */
/*        that */
/*           fl( a + 1.0 ) = a. */
	a = 1.;
	c__ = 1.;
L10:
	if (c__ == one) {
	    a *= 2;
	    c__ = dlamc3_(&a, &one);
	    d__1 = -a;
	    c__ = dlamc3_(&c__, &d__1);
	    goto L10;
	}
/*        Now compute  b = 2.0**m  with the smallest positive integer m */
/*        such that */
/*           fl( a + b ) .gt. a. */
	b = 1.;
	c__ = dlamc3_(&a, &b);
L20:
	if (c__ == a) {
	    b *= 2;
	    c__ = dlamc3_(&a, &b);
	    goto L20;
	}
/*        Now compute the base.  a and c  are neighbouring floating point */
/*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so */
/*        their difference is beta. Adding 0.25 to c is to ensure that it */
/*        is truncated to beta and not ( beta - 1 ). */
	qtr = one / 4;
	savec = c__;
	d__1 = -a;
	c__ = dlamc3_(&c__, &d__1);
	lbeta = (integer) (c__ + qtr);
/*        Now determine whether rounding or chopping occurs,  by adding a */
/*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */
	b = (doublereal) lbeta;
	d__1 = b / 2;
	d__2 = -b / 100;
	f = dlamc3_(&d__1, &d__2);
	c__ = dlamc3_(&f, &a);
	if (c__ == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	d__1 = b / 2;
	d__2 = b / 100;
	f = dlamc3_(&d__1, &d__2);
	c__ = dlamc3_(&f, &a);
	if (lrnd && c__ == a) {
	    lrnd = FALSE_;
	}
/*        Try and decide whether rounding is done in the  IEEE  'round to */
/*        nearest' style. B/2 is half a unit in the last place of the two */
/*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit */
/*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change */
/*        A, but adding B/2 to SAVEC should change SAVEC. */
	d__1 = b / 2;
	t1 = dlamc3_(&d__1, &a);
	d__1 = b / 2;
	t2 = dlamc3_(&d__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;
/*        Now find  the  mantissa, t.  It should  be the  integer part of */
/*        log to the base beta of a,  however it is safer to determine  t */
/*        by powering.  So we find t as the smallest positive integer for */
/*        which */
/*           fl( beta**t + 1.0 ) = 1.0. */
	lt = 0;
	a = 1.;
	c__ = 1.;
L30:
	if (c__ == one) {
	    ++lt;
	    a *= lbeta;
	    c__ = dlamc3_(&a, &one);
	    d__1 = -a;
	    c__ = dlamc3_(&c__, &d__1);
	    goto L30;
	}
    }
    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;
/*     End of DLAMC1 */
} /* dlamc1_ */

/* *********************************************************************** */
/* Subroutine */ int dlamc2_(integer *beta, integer *t, logical *rnd, 
	doublereal *eps, integer *emin, doublereal *rmin, integer *emax, 
	doublereal *rmax)
{
    /* Initialized data */

    static logical first = TRUE_;
    static logical iwarn = FALSE_;

    /* Format strings */
    static char fmt_9999[] = "(//\002 WARNING. The value EMIN may be incorre"
	    "ct:-\002,\002  EMIN = \002,i8,/\002 If, after inspection, the va"
	    "lue EMIN looks\002,\002 acceptable please comment out \002,/\002"
	    " the IF block as marked within the code of routine\002,\002 DLAM"
	    "C2,\002,/\002 otherwise supply EMIN explicitly.\002,/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal a, b, c__;
    static integer i__, lt;
    static doublereal one, two;
    static logical ieee;
    static doublereal half;
    static logical lrnd;
    static doublereal leps, zero;
    static integer lbeta;
    static doublereal rbase;
    static integer lemin, lemax, gnmin;
    static doublereal small;
    static integer gpmin;
    static doublereal third, lrmin, lrmax, sixth;
    extern /* Subroutine */ int dlamc1_(integer *, integer *, logical *, 
	    logical *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static logical lieee1;
    extern /* Subroutine */ int dlamc4_(integer *, doublereal *, integer *), 
	    dlamc5_(integer *, integer *, integer *, logical *, integer *, 
	    doublereal *);
    static integer ngnmin, ngpmin;

    /* Fortran I/O blocks */
    static cilist io___68 = { 0, 6, 0, fmt_9999, 0 };


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  DLAMC2 determines the machine parameters specified in its argument */
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
/*  EPS     (output) DOUBLE PRECISION */
/*          The smallest positive number such that */
/*             fl( 1.0 - EPS ) .LT. 1.0, */
/*          where fl denotes the computed value. */
/*  EMIN    (output) INTEGER */
/*          The minimum exponent before (gradual) underflow occurs. */
/*  RMIN    (output) DOUBLE PRECISION */
/*          The smallest normalized number for the machine, given by */
/*          BASE**( EMIN - 1 ), where  BASE  is the floating point value */
/*          of BETA. */
/*  EMAX    (output) INTEGER */
/*          The maximum exponent before overflow occurs. */
/*  RMAX    (output) DOUBLE PRECISION */
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
	zero = 0.;
	one = 1.;
	two = 2.;
/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of */
/*        BETA, T, RND, EPS, EMIN and RMIN. */
/*        Throughout this routine  we use the function  DLAMC3  to ensure */
/*        that relevant values are stored  and not held in registers,  or */
/*        are not affected by optimizers. */
/*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */
	dlamc1_(&lbeta, &lt, &lrnd, &lieee1);
/*        Start to find EPS. */
	b = (doublereal) lbeta;
	i__1 = -lt;
	a = pow_di(&b, &i__1);
	leps = a;
/*        Try some tricks to see whether or not this is the correct  EPS. */
	b = two / 3;
	half = one / 2;
	d__1 = -half;
	sixth = dlamc3_(&b, &d__1);
	third = dlamc3_(&sixth, &sixth);
	d__1 = -half;
	b = dlamc3_(&third, &d__1);
	b = dlamc3_(&b, &sixth);
	b = abs(b);
	if (b < leps) {
	    b = leps;
	}
	leps = 1.;
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    d__1 = half * leps;
/* Computing 5th power */
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
	    d__5 = leps;
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
	    c__ = dlamc3_(&d__1, &d__2);
	    d__1 = -c__;
	    c__ = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c__);
	    d__1 = -b;
	    c__ = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c__);
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
	    d__1 = small * rbase;
	    small = dlamc3_(&d__1, &zero);
/* L20: */
	}
	a = dlamc3_(&one, &small);
	dlamc4_(&ngpmin, &one, &lbeta);
	d__1 = -one;
	dlamc4_(&ngnmin, &d__1, &lbeta);
	dlamc4_(&gpmin, &a, &lbeta);
	d__1 = -a;
	dlamc4_(&gnmin, &d__1, &lbeta);
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
	    s_wsfe(&io___68);
	    do_fio(&c__1, (char *)&lemin, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
/* ** */
/*        Assume IEEE arithmetic if we found denormalised  numbers above, */
/*        or if arithmetic seems to round in the  IEEE style,  determined */
/*        in routine DLAMC1. A true IEEE machine should have both  things */
/*        true; however, faulty machines may have one or the other. */
	ieee = ieee || lieee1;
/*        Compute  RMIN by successive division by  BETA. We could compute */
/*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during */
/*        this computation. */
	lrmin = 1.;
	i__1 = 1 - lemin;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__1 = lrmin * rbase;
	    lrmin = dlamc3_(&d__1, &zero);
/* L30: */
	}
/*        Finally, call DLAMC5 to compute EMAX and RMAX. */
	dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
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
/*     End of DLAMC2 */
} /* dlamc2_ */

/* *********************************************************************** */
doublereal dlamc3_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing */
/*  the addition of  A  and  B ,  for use in situations where optimizers */
/*  might hold one of these in a register. */
/*  Arguments */
/*  ========= */
/*  A, B    (input) DOUBLE PRECISION */
/*          The values A and B. */
/* ===================================================================== */
/*     .. Executable Statements .. */
    ret_val = *a + *b;
    return ret_val;
/*     End of DLAMC3 */
} /* dlamc3_ */

/* *********************************************************************** */
/* Subroutine */ int dlamc4_(integer *emin, doublereal *start, integer *base)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal a;
    static integer i__;
    static doublereal b1, b2, c1, c2, d1, d2, one, zero, rbase;
    extern doublereal dlamc3_(doublereal *, doublereal *);

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  DLAMC4 is a service routine for DLAMC2. */
/*  Arguments */
/*  ========= */
/*  EMIN    (output) EMIN */
/*          The minimum exponent before (gradual) underflow, computed by */
/*          setting A = START and dividing by BASE until the previous A */
/*          can not be recovered. */
/*  START   (input) DOUBLE PRECISION */
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
    one = 1.;
    rbase = one / *base;
    zero = 0.;
    *emin = 1;
    d__1 = a * rbase;
    b1 = dlamc3_(&d__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	d__1 = a / *base;
	b1 = dlamc3_(&d__1, &zero);
	d__1 = b1 * *base;
	c1 = dlamc3_(&d__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d1 += b1;
/* L20: */
	}
	d__1 = a * rbase;
	b2 = dlamc3_(&d__1, &zero);
	d__1 = b2 / rbase;
	c2 = dlamc3_(&d__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
    return 0;
/*     End of DLAMC4 */
} /* dlamc4_ */

/* *********************************************************************** */
/* Subroutine */ int dlamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, doublereal *rmax)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal y, z__;
    static integer try__, lexp;
    static doublereal oldy;
    static integer uexp, nbits;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static doublereal recbas;
    static integer exbits, expsum;

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */
/*     .. Scalar Arguments .. */
/*     .. */
/*  Purpose */
/*  ======= */
/*  DLAMC5 attempts to compute RMAX, the largest machine floating-point */
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
/*  RMAX    (output) DOUBLE PRECISION */
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
    recbas = 1. / *beta;
    z__ = *beta - 1.;
    y = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ *= recbas;
	if (y < 1.) {
	    oldy = y;
	}
	y = dlamc3_(&y, &z__);
/* L20: */
    }
    if (y >= 1.) {
	y = oldy;
    }
/*     Now multiply by BETA**EMAX to get RMAX. */
    i__1 = *emax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = y * *beta;
	y = dlamc3_(&d__1, &c_b55);
/* L30: */
    }
    *rmax = y;
    return 0;
/*     End of DLAMC5 */
} /* dlamc5_ */

