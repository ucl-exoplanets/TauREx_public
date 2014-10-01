      DOUBLE PRECISION FUNCTION D1MACH(I)

C  DOUBLE-PRECISION MACHINE CONSTANTS (SEE R1MACH FOR DOCUMENTATION)

C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.

      INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4), SC
      DOUBLE PRECISION DMACH(5), EPS, EPSNEW, S
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
      EQUIVALENCE (DMACH(1),SMALL(1)), (DMACH(2),LARGE(1)),
     $            (DMACH(3),RIGHT(1)), (DMACH(4),DIVER(1)),
     $            (DMACH(5),LOG10(1))

      LOGICAL  PASS1
      SAVE     PASS1
      DATA     PASS1/.TRUE./

C IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T 3B SERIES AND
C MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T PC 7300),
C IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.

c      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2146435071,-1/,
c     $  (RIGHT(N),N=1,2)/1017118720,0/, (DIVER(N),N=1,2)/1018167296,0/,
c     $  (LOG10(N),N=1,2)/1070810131,1352628735/, SC/987/

C IEEE ARITHMETIC MACHINES AND 8087-BASED MICROS, SUCH AS THE IBM PC
C AND AT&T 6300, IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST.

C      DATA (SMALL(N),N=1,2)/0,1048576/, (LARGE(N),N=1,2)/-1,2146435071/,
C     $  (RIGHT(N),N=1,2)/0,1017118720/, (DIVER(N),N=1,2)/0,1018167296/,
C     $  (LOG10(N),N=1,2)/1352628735,1070810131/, SC/987/

C AMDAHL MACHINES.

C      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2147483647,-1/,
C     $ (RIGHT(N),N=1,2)/856686592,0/, (DIVER(N),N=1,2)/ 873463808,0/,
C     $ (LOG10(N),N=1,2)/1091781651,1352628735/, SC/987/

C BURROUGHS 1700 SYSTEM.

C      DATA (SMALL(N),N=1,2)/ZC00800000,Z000000000/,
C     $ (LARGE(N),N=1,2)/ZDFFFFFFFF,ZFFFFFFFFF/,
C     $ (RIGHT(N),N=1,2)/ZCC5800000,Z000000000/,
C     $ (DIVER(N),N=1,2)/ZCC6800000,Z000000000/,
C     $ (LOG10(N),N=1,2)/ZD00E730E7,ZC77800DC0/, SC/987/

C BURROUGHS 5700 SYSTEM.

C      DATA (SMALL(N),N=1,2)/O1771000000000000,O0000000000000000/,
C     $  (LARGE(N),N=1,2)/O0777777777777777,O0007777777777777/,
C     $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/,
C     $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/,
C     $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/

C BURROUGHS 6700/7700 SYSTEMS.

C      DATA (SMALL(N),N=1,2)/O1771000000000000,O7770000000000000/,
C     $  (LARGE(N),N=1,2)/O0777777777777777,O7777777777777777/,
C     $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/,
C     $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/,
C     $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/

C FTN4 ON THE CDC 6000/7000 SERIES.

C      DATA
C     $  (SMALL(N),N=1,2)/00564000000000000000B,00000000000000000000B/,
C     $  (LARGE(N),N=1,2)/37757777777777777777B,37157777777777777774B/,
C     $  (RIGHT(N),N=1,2)/15624000000000000000B,00000000000000000000B/,
C     $  (DIVER(N),N=1,2)/15634000000000000000B,00000000000000000000B/,
C     $  (LOG10(N),N=1,2)/17164642023241175717B,16367571421742254654B/,
C     $  SC/987/

C FTN5 ON THE CDC 6000/7000 SERIES.

C      DATA
C     $(SMALL(N),N=1,2)/O"00564000000000000000",O"00000000000000000000"/,
C     $(LARGE(N),N=1,2)/O"37757777777777777777",O"37157777777777777774"/,
C     $(RIGHT(N),N=1,2)/O"15624000000000000000",O"00000000000000000000"/,
C     $(DIVER(N),N=1,2)/O"15634000000000000000",O"00000000000000000000"/,
C     $(LOG10(N),N=1,2)/O"17164642023241175717",O"16367571421742254654"/,
C     $ SC/987/

C CONVEX C-1

C      DATA (SMALL(N),N=1,2)/'00100000'X,'00000000'X/,
C     $  (LARGE(N),N=1,2)/'7FFFFFFF'X,'FFFFFFFF'X/,
C     $  (RIGHT(N),N=1,2)/'3CC00000'X,'00000000'X/,
C     $  (DIVER(N),N=1,2)/'3CD00000'X,'00000000'X/,
C     $  (LOG10(N),N=1,2)/'3FF34413'X,'509F79FF'X/, SC/987/

C CRAY 1, XMP, 2, AND 3.

C      DATA
C     $ (SMALL(N),N=1,2)/201354000000000000000B,000000000000000000000B/,
C     $ (LARGE(N),N=1,2)/577767777777777777777B,000007777777777777776B/,
C     $ (RIGHT(N),N=1,2)/376434000000000000000B,000000000000000000000B/,
C     $ (DIVER(N),N=1,2)/376444000000000000000B,000000000000000000000B/,
C     $ (LOG10(N),N=1,2)/377774642023241175717B,000007571421742254654B/,
C     $ SC/987/

C DATA GENERAL ECLIPSE S/200
C NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C STATIC DMACH(5)

C      DATA SMALL/20K,3*0/, LARGE/77777K,3*177777K/,
C     $  RIGHT/31420K,3*0/, DIVER/32020K,3*0/,
C     $  LOG10/40423K,42023K,50237K,74776K/, SC/987/

C HARRIS SLASH 6 AND SLASH 7

C      DATA (SMALL(N),N=1,2)/'20000000,'00000201/,
C     $  (LARGE(N),N=1,2)/'37777777,'37777577/,
C     $  (RIGHT(N),N=1,2)/'20000000,'00000333/,
C     $  (DIVER(N),N=1,2)/'20000000,'00000334/,
C     $  (LOG10(N),N=1,2)/'23210115,'10237777/, SC/987/

C HONEYWELL DPS 8/70 SERIES.

C      DATA (SMALL(N),N=1,2)/O402400000000,O000000000000/,
C     $  (LARGE(N),N=1,2)/O376777777777,O777777777777/,
C     $  (RIGHT(N),N=1,2)/O604400000000,O000000000000/,
C     $  (DIVER(N),N=1,2)/O606400000000,O000000000000/,
C     $  (LOG10(N),N=1,2)/O776464202324,O117571775714/, SC/987/

C IBM 360/370 SERIES, XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.

C      DATA (SMALL(N),N=1,2)/Z00100000,Z00000000/,
C     $  (LARGE(N),N=1,2)/Z7FFFFFFF,ZFFFFFFFF/,
C     $  (RIGHT(N),N=1,2)/Z33100000,Z00000000/,
C     $  (DIVER(N),N=1,2)/Z34100000,Z00000000/,
C     $  (LOG10(N),N=1,2)/Z41134413,Z509F79FF/, SC/987/

C INTERDATA 8/32 WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.

C      DATA (SMALL(N),N=1,2)/Z'00100000',Z'00000000'/,
C     $  (LARGE(N),N=1,2)/Z'7EFFFFFF',Z'FFFFFFFF'/,
C     $  (RIGHT(N),N=1,2)/Z'33100000',Z'00000000'/,
C     $  (DIVER(N),N=1,2)/Z'34100000',Z'00000000'/,
C     $  (LOG10(N),N=1,2)/Z'41134413',Z'509F79FF'/, SC/987/

C PDP-10 (KA PROCESSOR).

C      DATA (SMALL(N),N=1,2)/"033400000000,"000000000000/,
C     $  (LARGE(N),N=1,2)/"377777777777,"344777777777/,
C     $  (RIGHT(N),N=1,2)/"113400000000,"000000000000/,
C     $  (DIVER(N),N=1,2)/"114400000000,"000000000000/,
C     $  (LOG10(N),N=1,2)/"177464202324,"144117571776/, SC/987/

C PDP-10 (KI PROCESSOR).

C      DATA (SMALL(N),N=1,2)/"000400000000,"000000000000/,
C     $  (LARGE(N),N=1,2)/"377777777777,"377777777777/,
C     $  (RIGHT(N),N=1,2)/"103400000000,"000000000000/,
C     $  (DIVER(N),N=1,2)/"104400000000,"000000000000/,
C     $  (LOG10(N),N=1,2)/"177464202324,"047674776746/, SC/987/

C PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS
C (EXPRESSED IN INTEGER AND OCTAL).

C      DATA (SMALL(N),N=1,2)/8388608,0/, (LARGE(N),N=1,2)/2147483647,-1/,
C     $  (RIGHT(N),N=1,2)/612368384,0/, (DIVER(N),N=1,2)/620756992,0/,
C     $  (LOG10(N),N=1,2)/1067065498,-2063872008/, SC/987/

C      DATA (SMALL(N),N=1,2)/O00040000000,O00000000000/,
C     $  (LARGE(N),N=1,2)/O17777777777,O37777777777/,
C     $  (RIGHT(N),N=1,2)/O04440000000,O00000000000/,
C     $  (DIVER(N),N=1,2)/O04500000000,O00000000000/,
C     $  (LOG10(N),N=1,2)/O07746420232,O20476747770/, SC/987/

C PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS
C (EXPRESSED IN INTEGER AND OCTAL).

C      DATA SMALL/128,3*0/, LARGE/32767,3*-1/, RIGHT/9344,3*0/,
C     $  DIVER/9472,3*0/, LOG10/16282,8346,-31493,-12296/, SC/987/

C      DATA SMALL/O000200,3*O000000/, LARGE/O077777,3*O177777/,
C     $  RIGHT/O022200,3*O000000/, DIVER/O022400,3*O000000/,
C     $  LOG10/O037632,O020232,O102373,O147770/, SC/987/

C PRIME 50 SERIES SYSTEMS WITH 32-BIT INTEGERS AND 64V MODE
C INSTRUCTIONS, SUPPLIED BY IGOR BRAY.

C      DATA (SMALL(N),N=1,2)/:10000000000,:00000100001/,
C     $  (LARGE(N),N=1,2)/:17777777777,:37777677775/,
C     $  (RIGHT(N),N=1,2)/:10000000000,:00000000122/,
C     $  (DIVER(N),N=1,2)/:10000000000,:00000000123/,
C     $  (LOG10(N),N=1,2)/:11504046501,:07674600177/, SC/987/

C SEQUENT BALANCE 8000

C      DATA (SMALL(N),N=1,2)/$00000000, $00100000/,
C     $  (LARGE(N),N=1,2)/$FFFFFFFF, $7FEFFFFF/,
C     $  (RIGHT(N),N=1,2)/$00000000, $3CA00000/,
C     $  (DIVER(N),N=1,2)/$00000000, $3CB00000/,
C     $  (LOG10(N),N=1,2)/$509F79FF, $3FD34413/, SC/987/

C UNIVAC 1100 SERIES.

C      DATA (SMALL(N),N=1,2)/O000040000000,O000000000000/,
C     $  (LARGE(N),N=1,2)/O377777777777,O777777777777/,
C     $  (RIGHT(N),N=1,2)/O170540000000,O000000000000/,
C     $  (DIVER(N),N=1,2)/O170640000000,O000000000000/,
C     $  (LOG10(N),N=1,2)/O177746420232,O411757177572/, SC/987/

C VAX UNIX F77 COMPILER

C      DATA (SMALL(N),N=1,2)/128,0/, (LARGE(N),N=1,2)/-32769,-1/,
C     $  (RIGHT(N),N=1,2)/9344,0/, (DIVER(N),N=1,2)/9472,0/,
C     $  (LOG10(N),N=1,2)/546979738,-805796613/, SC/987/

C VAX-11 WITH FORTRAN IV-PLUS COMPILER

C      DATA (SMALL(N),N=1,2)/Z00000080,Z00000000/,
C     $  (LARGE(N),N=1,2)/ZFFFF7FFF,ZFFFFFFFF/,
C     $  (RIGHT(N),N=1,2)/Z00002480,Z00000000/,
C     $  (DIVER(N),N=1,2)/Z00002500,Z00000000/,
C     $  (LOG10(N),N=1,2)/Z209A3F9A,ZCFF884FB/, SC/987/

C VAX/VMS VERSION 2.2

C      DATA (SMALL(N),N=1,2)/'80'X,'0'X/,
C     $  (LARGE(N),N=1,2)/'FFFF7FFF'X,'FFFFFFFF'X/,
C     $  (RIGHT(N),N=1,2)/'2480'X,'0'X/, (DIVER(N),N=1,2)/'2500'X,'0'X/,
C     $  (LOG10(N),N=1,2)/'209A3F9A'X,'CFF884FB'X/, SC/987/

      IF( PASS1 )  THEN

         PASS1 = .FALSE.
c         IF (SC.NE.987)
c     $       CALL ERRMSG( 'D1MACH--NO DATA STATEMENTS ACTIVE',.TRUE.)
c         sc = 987
         dmach(1) = DLAMCH( 'S' )
         dmach(2) = DLAMCH( 'O' )
         dmach(3) = DLAMCH( 'E' )
         dmach(4) = DLAMCH( 'P' )
         dmach(5) = dlog10(DLAMCH( 'B' ))

C                        ** CALCULATE MACHINE PRECISION
         EPSNEW = 0.01D0
   10    EPS = EPSNEW
            EPSNEW = EPSNEW / 1.1D0
C                               ** IMPORTANT TO STORE 'S' SINCE MAY BE
C                               ** KEPT IN HIGHER PRECISION IN REGISTERS
            S = 1.D0 + EPSNEW
c	write(*,'(1a,1p5e15.7)') 's, epsnew,eps,dmach(4),eps/dmach(4)',
c     - s, epsnew,eps,dmach(4),eps/dmach(4)
            IF( S.GT.1.D0 ) GO TO 10
         IF( EPS/DMACH(4).LT.0.5D0 .OR. EPS/DMACH(4).GT.2.D0 )
     $       CALL ERRMSG( 'D1MACH--TABULATED PRECISION WRONG',.TRUE.)

      END IF

      IF (I.LT.1.OR.I.GT.5)
     $    CALL ERRMSG( 'D1MACH--ARGUMENT OUT OF BOUNDS',.TRUE.)
      D1MACH = DMACH(I)
      RETURN
      END
c....................................
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..

*  Purpose
*  =======

*  DLAMCH determines double precision machine parameters.

*  Arguments
*  =========

*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax

*          where

*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)

* =====================================================================

*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..

      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN

*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.

            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF

      IF( CMACH.EQ.'E' .OR. CMACH.EQ.'e' ) THEN
         RMACH = EPS
      ELSE IF( CMACH.EQ.'S' .OR. CMACH.EQ.'s' ) THEN
         RMACH = SFMIN
      ELSE IF( CMACH.EQ.'B' .OR. CMACH.EQ.'b' ) THEN
         RMACH = BASE
      ELSE IF( CMACH.EQ.'P' .OR. CMACH.EQ.'p' ) THEN
         RMACH = PREC
      ELSE IF( CMACH.EQ.'U' .OR. CMACH.EQ.'u' ) THEN
         RMACH = RMIN
      ELSE IF( CMACH.EQ.'O' .OR. CMACH.EQ.'o' ) THEN
         RMACH = RMAX
      ELSE IF( CMACH.EQ.'N' .OR. CMACH.EQ.'n' ) THEN
         RMACH = T
      ELSE IF( CMACH.EQ.'R' .OR. CMACH.EQ.'r' ) THEN
         RMACH = RND
      ELSE IF( CMACH.EQ.'M' .OR. CMACH.EQ.'m' ) THEN
         RMACH = EMIN
      ELSE IF( CMACH.EQ.'L' .OR. CMACH.EQ.'l' ) THEN
         RMACH = EMAX
      ELSE
         CALL ErrMsg('SLAMCH--invalid argument',.True.)
      END IF

      DLAMCH = RMACH
      RETURN

*     End of DLAMCH

      END

************************************************************************

      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..

*  Purpose
*  =======

*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.

*  Arguments
*  =========

*  BETA    (output) INTEGER
*          The base of the machine.

*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.

*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.

*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.

*  Further Details
*  ===============

*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See

*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.

*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.

* =====================================================================

*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..

      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1

*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.

*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.

*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that

*           fl( a + 1.0 ) = a.

         A = 1
         C = 1

   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF

*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that

*           fl( a + b ) .gt. a.

         B = 1
         C = DLAMC3( A, B )

   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF

*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).

         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR

*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.

         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( LRND .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.

*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.

         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND

*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which

*           fl( beta**t + 1.0 ) = 1.0.

         LT = 0
         A = 1
         C = 1

   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF

      END IF

      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN

*     End of DLAMC1

      END

************************************************************************

      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..

*  Purpose
*  =======

*  DLAMC2 determines the machine parameters specified in its argument
*  list.

*  Arguments
*  =========

*  BETA    (output) INTEGER
*          The base of the machine.

*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.

*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.

*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that

*             fl( 1.0 - EPS ) .LT. 1.0,

*          where fl denotes the computed value.

*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.

*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.

*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.

*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.

*  Further Details
*  ===============

*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.

* =====================================================================

*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..

      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2

*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.

*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.

*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.

         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )

*        Start to find EPS.

         B = LBETA
         A = B**( -LT )
         LEPS = A

*        Try some tricks to see whether or not this is the correct  EPS.

         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS

         LEPS = 1

   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF

         IF( A.LT.LEPS )
     $      LEPS = A
c      write(*,*) ' EPS from PARANOIA:', LEPS

*        Computation of EPS complete.

*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.

         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.

         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF

         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF

         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF

         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***

*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.

         IEEE = IEEE .OR. LIEEE1

*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.

         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE

*        Finally, call DLAMC5 to compute EMAX and RMAX.

         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF

      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX

      RETURN

 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )

*     End of DLAMC2

      END

************************************************************************

      DOUBLE PRECISION FUNCTION DLAMC3( A, B )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..

*  Purpose
*  =======

*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.

*  Arguments
*  =========

*  A, B    (input) DOUBLE PRECISION
*          The values A and B.

* =====================================================================

*     .. Executable Statements ..

      DLAMC3 = A + B

      RETURN

*     End of DLAMC3

      END

************************************************************************

      SUBROUTINE DLAMC4( EMIN, START, BASE )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..

*  Purpose
*  =======

*  DLAMC4 is a service routine for DLAMC2.

*  Arguments
*  =========

*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.

*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.

*  BASE    (input) INTEGER
*          The base of the machine.

* =====================================================================

*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..

      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF

      RETURN

*     End of DLAMC4

      END

************************************************************************

      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..

*  Purpose
*  =======

*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.

*  Arguments
*  =========

*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.

*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.

*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.

*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.

*  EMAX    (output) INTEGER
*          The largest exponent before overflow

*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.

* =====================================================================

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..

*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).

      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF

*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.

      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF

*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .

      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P

*     NBITS is the total number of bits needed to store a
*     floating-point number.

      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN

*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.

         EMAX = EMAX - 1
      END IF

      IF( IEEE ) THEN

*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.

         EMAX = EMAX - 1
      END IF

*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .

*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .

      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY

*     Now multiply by BETA**EMAX to get RMAX.

      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE

      RMAX = Y
      RETURN

*     End of DLAMC5

      END
