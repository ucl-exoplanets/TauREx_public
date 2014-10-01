      SUBROUTINE  QGAUSN( M, GMU, GWT )

C       Compute weights and abscissae for ordinary gaussian quadrature
C       (no weight function inside integral) on the interval (0,1)

C   INPUT :    M                     order of quadrature rule

C   OUTPUT :  GMU(I)  I = 1 TO M,    array of abscissae
C             GWT(I)  I = 1 TO M,    array of weights

C   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
C                   Integration, Academic Press, New York, pp. 87, 1975.

C   METHOD:  Compute the abscissae as roots of the Legendre
C            polynomial P-sub-M using a cubically convergent
C            refinement of Newton's method.  Compute the
C            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
C            that Newton's method can very easily diverge; only a
C            very good initial guess can guarantee convergence.
C            The initial guess used here has never led to divergence
C            even for M up to 1000.

C   ACCURACY:  at least 13 significant digits

C   INTERNAL VARIABLES:

C    ITER      : number of Newton Method iterations
C    MAXIT     : maximum allowed iterations of Newton Method
C    PM2,PM1,P : 3 successive Legendre polynomials
C    PPR       : derivative of Legendre polynomial
C    P2PRI     : 2nd derivative of Legendre polynomial
C    TOL       : convergence criterion for Legendre poly root iteration
C    X,XI      : successive iterates in cubically-convergent version
C                of Newtons Method (seeking roots of Legendre poly.)
C+---------------------------------------------------------------------+
      DOUBLE PRECISION     CONA, GMU(*), GWT(*), PI, T
      INTEGER  ITER, LIM, M, MAXIT, NP1
      DOUBLE   PRECISION  D1MACH
      DOUBLE   PRECISION  EN, NNP1, ONE, P, PM1, PM2, PPR, P2PRI, PROD,
     $                    TMP, TOL, TWO, X, XI
      SAVE     PI, TOL
      DATA     PI / 0.0 /,  MAXIT / 1000 /,  ONE / 1.D0 /,  TWO / 2.D0 /


      IF ( PI.EQ.0.0 .or. tol .eq. 0.0 )  THEN
         PI = 2. * ASIN(1.0)
         TOL = 10. * D1MACH(4)
      END IF

      IF ( M.LT.1 )  CALL ErrMsg( 'QGAUSN--Bad value of M', .TRUE. )
      IF ( M.EQ.1 )  THEN
         GMU( 1 ) = 0.5
         GWT( 1 ) = 1.0
         RETURN
      END IF

      EN   = M
      NP1  = M + 1
      NNP1 = M * NP1
      CONA = FLOAT( M-1 ) / ( 8 * M**3 )

      LIM  = M / 2
      DO 30  K = 1, LIM
C                                        ** initial guess for k-th root
C                                        ** of Legendre polynomial, from
C                                        ** Davis/Rabinowitz (2.7.3.3a)
         T = ( 4*K - 1 ) * PI / ( 4*M + 2 )
         X = COS ( T + CONA / TAN( T ) )
         ITER = 0
C                                        ** upward recurrence for
C                                        ** Legendre polynomials
   10    ITER = ITER + 1
         PM2 = ONE
         PM1 = X
         DO 20 NN = 2, M
            P   = ( ( 2*NN - 1 ) * X * PM1 - ( NN-1 ) * PM2 ) / NN
            PM2 = PM1
            PM1 = P
   20    CONTINUE
C                                              ** Newton Method
         TMP   = ONE / ( ONE - X**2 )
         PPR   = EN * ( PM2 - X * P ) * TMP
         P2PRI = ( TWO * X * PPR - NNP1 * P ) * TMP
         XI    = X - ( P / PPR ) * ( ONE +
     $               ( P / PPR ) * P2PRI / ( TWO * PPR ) )

C                                              ** check for convergence
         IF ( ABS(XI-X) .GT. TOL ) THEN
            IF( ITER.GT.MAXIT )
     $          CALL ERRMSG( 'QGAUSN--MAX ITERATION COUNT', .TRUE. )
            X = XI
            GO TO 10
         END IF
C                             ** iteration finished--calculate weights,
C                             ** abscissae for (-1,1)
         GMU( K ) = - X
         GWT( K ) = TWO / ( TMP * ( EN * PM2 )**2 )
         GMU( NP1 - K ) = - GMU( K )
         GWT( NP1 - K ) =   GWT( K )
  30  CONTINUE
C                                    ** set middle abscissa and weight
C                                    ** for rules of odd order
      IF ( MOD( M,2 ) .NE. 0 )  THEN
         GMU( LIM + 1 ) = 0.0
         PROD = ONE
         DO 40 K = 3, M, 2
            PROD = PROD * K / ( K-1 )
  40     CONTINUE
         GWT( LIM + 1 ) = TWO / PROD**2
      END IF
C                                        ** convert from (-1,1) to (0,1)
      DO 50  K = 1, M
         GMU( K ) = 0.5 * GMU( K ) + 0.5
         GWT( K ) = 0.5 * GWT( K )
  50  CONTINUE

      RETURN
      END
