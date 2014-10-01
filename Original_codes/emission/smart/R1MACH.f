      REAL FUNCTION R1MACH(I)

C  SINGLE-PRECISION MACHINE CONSTANTS

C  ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C  BASE-B FORM

C     (SIGN) (B**E) * ( (X(1)/B) + ... + (X(T)/B**T) )

C  WHERE 0.LE.X(I).LT.B, I=1,...,T, X(1).GT.0; EMIN.LE.E.LE.EMAX.
C  FOR EXAMPLE, IN BASE 10, THE X(I) WOULD JUST BE THE DIGITS
C  FOLLOWING THE DECIMAL POINT, AND CLEARLY THE FIRST DIGIT (X(1))
C  WOULD HAVE TO BE NONZERO OR THE DECIMAL POINT COULD BE MOVED
C  OVER WITH A CORRESPONDING CHANGE IN 'E'.  THEN

C  R1MACH(1) = B**(EMIN-1), SMALLEST POSITIVE MAGNITUDE
C                 (E = EMIN, X(1) = 1, ALL OTHER X(I) = 0).
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), LARGEST MAGNITUDE
C                 (E = EMAX, X(I) = B-1).
C  R1MACH(3) = 1/B**T, SMALLEST RELATIVE SPACING.
C  R1MACH(4) = 1/B**(T-1), MACHINE PRECISION (X(1)=X(T)=1,ALL OTHER
C                 X(I) = 0); SMALLEST POSITIVE EPS SUCH THAT 1+EPS.NE.1
C  R1MACH(5) = LOG10(B)

C  REFERENCE: FOX P.A., HALL A.D., SCHRYER N.L.,'FRAMEWORK FOR A
C               PORTABLE LIBRARY', ACM TRANSACTIONS ON MATHEMATICAL
C               SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.

C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  DELETING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)

C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.

C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.

      INTEGER SMALL(2), LARGE(2), RIGHT(2), DIVER(2), LOG10(2), SC
      REAL RMACH(5)
      REAL              SLAMCH
      EXTERNAL          SLAMCH

      EQUIVALENCE (RMACH(1),SMALL(1)), (RMACH(2),LARGE(1)),
     $            (RMACH(3),RIGHT(1)), (RMACH(4),DIVER(1)),
     $            (RMACH(5),LOG10(1))

      LOGICAL  PASS1
      SAVE     PASS1
      DATA     PASS1/.TRUE./

C IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C 3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).

c      DATA SMALL(1)/8388608/, LARGE(1)/2139095039/,
c     $     RIGHT(1)/864026624/, DIVER(1)/872415232/,
c     $     LOG10(1)/ 1050288283/, SC/987/

C AMDAHL MACHINES.

C      DATA SMALL(1)/1048576/, LARGE(1)/2147483647/,
C     $     RIGHT(1)/990904320/, DIVER(1)/1007681536/,
C     $     LOG10(1)/1091781651/, SC/987/

C BURROUGHS 1700 SYSTEM.

C      DATA RMACH/Z400800000,Z5FFFFFFFF,Z4E9800000,Z4EA800000,
C     $             Z500E730E8/, SC/987/

C BURROUGHS 5700/6700/7700 SYSTEMS.

C      DATA RMACH/O1771000000000000,O0777777777777777,O1311000000000000,
C     $             O1301000000000000,O1157163034761675/, SC/987/

C FTN4 ON CDC 6000/7000 SERIES.

C      DATA RMACH/00564000000000000000B,37767777777777777776B,
C     $ 16414000000000000000B,16424000000000000000B,
C     $ 17164642023241175720B/, SC/987/

C FTN5 ON CDC 6000/7000 SERIES.

C      DATA RMACH/O"00564000000000000000",O"37767777777777777776",
C     $ O"16414000000000000000",O"16424000000000000000",
C     $ O"17164642023241175720"/, SC/987/

C CONVEX C-1.

C      DATA RMACH/'00800000'X,'7FFFFFFF'X,'34800000'X,
C     $ '35000000'X,'3F9A209B'X/, SC/987/

C CRAY 1, XMP, 2, AND 3.

C      DATA RMACH/200034000000000000000B,577767777777777777776B,
C     $ 377224000000000000000B,377234000000000000000B,
C     $ 377774642023241175720B/, SC/987/

C DATA GENERAL ECLIPSE S/200.
C NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C STATIC RMACH(5)

C      DATA SMALL/20K,0/, LARGE/77777K,177777K/, RIGHT/35420K,0/,
C     $  DIVER/36020K,0/, LOG10/40423K,42023K/, SC/987/

C HARRIS SLASH 6 AND SLASH 7.

C      DATA SMALL/'20000000,'00000201/, LARGE/'37777777,'00000177/,
C     $  RIGHT/'20000000,'00000352/, DIVER/'20000000,'00000353/,
C     $  LOG10/'23210115,'00000377/, SC/987/

C HONEYWELL DPS 8/70 SERIES.

C      DATA RMACH/O402400000000,O376777777777,O714400000000,
C     $ O716400000000,O776464202324/, SC/987/

C IBM 360/370 SERIES,
C XEROX SIGMA 5/7/9 AND SEL SYSTEMS 85/86.

C      DATA RMACH/Z00100000,Z7FFFFFFF,Z3B100000,Z3C100000,
C     $ Z41134413/, SC/987/

C INTERDATA 8/32 WITH UNIX SYSTEM FORTRAN 77 COMPILER.
C FOR INTERDATA FORTRAN VII COMPILER REPLACE
C Z'S SPECIFYING HEX CONSTANTS WITH Y'S.

C      DATA RMACH/Z'00100000',Z'7EFFFFFF',Z'3B100000',Z'3C100000',
C     $ Z'41134413'/, SC/987/

C PDP-10 (KA OR KI PROCESSOR).

C      DATA RMACH/"000400000000,"377777777777,"146400000000,
C     $ "147400000000,"177464202324/, SC/987/

C PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS
C (EXPRESSED IN INTEGER AND OCTAL).

C      DATA SMALL(1)/8388608/, LARGE(1)/2147483647/,
C     $  RIGHT(1)/880803840/, DIVER(1)/889192448/,
C     $  LOG10(1)/1067065499/, SC/987/

C      DATA RMACH/O00040000000,O17777777777,O06440000000,
C     $ O06500000000,O07746420233/, SC/987/

C PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS
C (EXPRESSED IN INTEGER AND OCTAL).

C      DATA SMALL/128,0/, LARGE/32767,-1/, RIGHT/13440,0/,
C     $  DIVER/13568,0/, LOG10/16282,8347/, SC/987/

C      DATA SMALL/O000200,O000000/, LARGE/O077777,O177777/,
C     $  RIGHT/O032200,O000000/, DIVER/O032400,O000000/,
C     $  LOG10/O037632,O020233/, SC/987/

C SEQUENT BALANCE 8000.

C      DATA SMALL(1)/$00800000/, LARGE(1)/$7F7FFFFF/,
C     $  RIGHT(1)/$33800000/, DIVER(1)/$34000000/,
C     $  LOG10(1)/$3E9A209B/, SC/987/

C UNIVAC 1100 SERIES.

C      DATA RMACH/O000400000000,O377777777777,O146400000000,
C     $ O147400000000,O177464202324/, SC/987/

C VAX UNIX F77 COMPILER.

C      DATA SMALL(1)/128/, LARGE(1)/-32769/, RIGHT(1)/13440/,
C     $  DIVER(1)/13568/, LOG10(1)/547045274/, SC/987/

C VAX-11 WITH FORTRAN IV-PLUS COMPILER.

C      DATA RMACH/Z00000080,ZFFFF7FFF,Z00003480,Z00003500,
C     $ Z209B3F9A/, SC/987/

C VAX/VMS VERSION 2.2.

C      DATA RMACH/'80'X,'FFFF7FFF'X,'3480'X,'3500'X,
C     $ '209B3F9A'X/, SC/987/

      IF( PASS1 )  THEN

         PASS1 = .FALSE.
c         IF (SC.NE.987)
c     $       CALL ERRMSG( 'R1MACH--NO DATA STATEMENTS ACTIVE',.TRUE.)

c         sc = 987
         rmach(1) = SLAMCH( 'S' )
         rmach(2) = SLAMCH( 'O' )
         rmach(3) = SLAMCH( 'E' )
         rmach(4) = SLAMCH( 'P' )
         rmach(5) = alog10(SLAMCH( 'B' ))

C                      ** CALCULATE MACHINE PRECISION
         EPSNEW = 0.01
   10    EPS = EPSNEW
            EPSNEW = EPSNEW / 1.1
C                               ** IMPORTANT TO STORE 'S' SINCE MAY BE
C                               ** KEPT IN HIGHER PRECISION IN REGISTERS
            S = 1.0 + EPSNEW
            IF( S.GT.1.0 ) GO TO 10
         RATIO = EPS / RMACH(4)
c	write(*,'(1a,1p5e15.7)') 's, epsnew,eps,rmach(4),eps/rmach(4)',
c     - s, epsnew,eps,rmach(4),eps/rmach(4)

         IF( RATIO.LT.0.5 .OR. RATIO.GT.2.0 )
     $       CALL ERRMSG( 'R1MACH--TABULATED PRECISION WRONG',.TRUE.)

      END IF

      IF (I.LT.1.OR.I.GT.5)
     $    CALL ERRMSG( 'R1MACH--ARGUMENT OUT OF BOUNDS',.TRUE.)

      R1MACH = RMACH(I)

      RETURN
      END
c.....................
      REAL             FUNCTION SLAMCH( CMACH )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..

*  Purpose
*  =======

*  SLAMCH determines single precision machine parameters.

*  Arguments
*  =========

*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by SLAMCH:
*          = 'E' or 'e',   SLAMCH := eps
*          = 'S' or 's ,   SLAMCH := sfmin
*          = 'B' or 'b',   SLAMCH := base
*          = 'P' or 'p',   SLAMCH := eps*base
*          = 'N' or 'n',   SLAMCH := t
*          = 'R' or 'r',   SLAMCH := rnd
*          = 'M' or 'm',   SLAMCH := emin
*          = 'U' or 'u',   SLAMCH := rmin
*          = 'L' or 'l',   SLAMCH := emax
*          = 'O' or 'o',   SLAMCH := rmax

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
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAMC2
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
         CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
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

      SLAMCH = RMACH
      RETURN

*     End of SLAMCH

      END

************************************************************************

      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )

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

*  SLAMC1 determines the machine parameters given by BETA, T, RND, and
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
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
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

*        Throughout this routine  we use the function  SLAMC3  to ensure
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
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF

*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that

*           fl( a + b ) .gt. a.

         B = 1
         C = SLAMC3( A, B )

   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF

*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).

         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR

*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.

         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( LRND .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.

*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.

         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
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
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF

      END IF

      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN

*     End of SLAMC1

      END

************************************************************************

      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL               EPS, RMAX, RMIN
*     ..

*  Purpose
*  =======

*  SLAMC2 determines the machine parameters specified in its argument
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

*  EPS     (output) REAL
*          The smallest positive number such that

*             fl( 1.0 - EPS ) .LT. 1.0,

*          where fl denotes the computed value.

*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.

*  RMIN    (output) REAL
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.

*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.

*  RMAX    (output) REAL
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
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAMC1, SLAMC4, SLAMC5
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

*        Throughout this routine  we use the function  SLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.

*        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.

         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )

*        Start to find EPS.

         B = LBETA
         A = B**( -LT )
         LEPS = A

*        Try some tricks to see whether or not this is the correct  EPS.

         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS

         LEPS = 1

   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
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
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
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
*        in routine SLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.

         IEEE = IEEE .OR. LIEEE1

*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.

         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE

*        Finally, call SLAMC5 to compute EMAX and RMAX.

         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
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
     $      ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )

*     End of SLAMC2

      END

************************************************************************

      REAL             FUNCTION SLAMC3( A, B )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      REAL               A, B
*     ..

*  Purpose
*  =======

*  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.

*  Arguments
*  =========

*  A, B    (input) REAL
*          The values A and B.

* =====================================================================

*     .. Executable Statements ..

      SLAMC3 = A + B

      RETURN

*     End of SLAMC3

      END

************************************************************************

      SUBROUTINE SLAMC4( EMIN, START, BASE )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL               START
*     ..

*  Purpose
*  =======

*  SLAMC4 is a service routine for SLAMC2.

*  Arguments
*  =========

*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.

*  START   (input) REAL
*          The starting point for determining EMIN.

*  BASE    (input) INTEGER
*          The base of the machine.

* =====================================================================

*     .. Local Scalars ..
      INTEGER            I
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. Executable Statements ..

      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF

      RETURN

*     End of SLAMC4

      END

************************************************************************

      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )

*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992

*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL               RMAX
*     ..

*  Purpose
*  =======

*  SLAMC5 attempts to compute RMAX, the largest machine floating-point
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

*  RMAX    (output) REAL
*          The largest machine floating-point number.

* =====================================================================

*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL               OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
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
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY

*     Now multiply by BETA**EMAX to get RMAX.

      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE

      RMAX = Y
      RETURN

*     End of SLAMC5

      END
