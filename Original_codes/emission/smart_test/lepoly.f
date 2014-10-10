      SUBROUTINE  LEPOLY( NMU, M, MAXMU, TWONM1, MU, YLM )

C       COMPUTES THE NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL,
C       DEFINED IN TERMS OF THE ASSOCIATED LEGENDRE POLYNOMIAL
C       PLM = P-SUB-L-SUPER-M AS

C             YLM(MU) = SQRT( (L-M)!/(L+M)! ) * PLM(MU)

C       FOR FIXED ORDER -M- AND ALL DEGREES FROM L = M TO TWONM1.
C       WHEN M.GT.0, ASSUMES THAT Y-SUB(M-1)-SUPER(M-1) IS AVAILABLE
C       FROM A PRIOR CALL TO THE ROUTINE.

C       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
C                  High-Order Associated Legendre Polynomials,
C                  J. Quant. Spectrosc. Radiat. Transfer 10,
C                  557-562, 1970.  (hereafter D/A)

C       METHOD: Varying degree recurrence relationship.

C       NOTE 1: The D/A formulas are transformed by
C               setting  M = n-1; L = k-1.
C       NOTE 2: Assumes that routine is called first with  M = 0,
C               then with  M = 1, etc. up to  M = TWONM1.
C       NOTE 3: Loops are written in such a way as to vectorize.

C  I N P U T     V A R I A B L E S:

C       NMU    :  NUMBER OF ARGUMENTS OF -YLM-
C       M      :  ORDER OF -YLM-
C       MAXMU  :  FIRST DIMENSION OF -YLM-
C       TWONM1 :  MAX DEGREE OF -YLM-
C       MU(I)  :  I = 1 TO NMU, ARGUMENTS OF -YLM-
C       IF M.GT.0, YLM(M-1,I) FOR I = 1 TO NMU IS REQUIRED

C  O U T P U T     V A R I A B L E:

C       YLM(L,I) :  L = M TO TWONM1, NORMALIZED ASSOCIATED LEGENDRE
C                   POLYNOMIALS EVALUATED AT ARGUMENT -MU(I)-
C+---------------------------------------------------------------------+
      REAL     MU(*), YLM( 0:MAXMU,* )
      INTEGER  M, NMU, TWONM1
      PARAMETER  ( MAXSQT = 1000 )
      REAL     SQT( MAXSQT )
      LOGICAL  PASS1
      SAVE  SQT, PASS1
      DATA  PASS1 / .TRUE. /


      IF ( PASS1 )  THEN
         PASS1 = .FALSE.
         DO 1  NS = 1, MAXSQT
            SQT( NS ) = SQRT( FLOAT(NS) )
    1    CONTINUE
      ENDIF

c      IF ( 2*TWONM1 .GT. MAXSQT )
c     $   CALL ERRMSG( 'LEPOLY--NEED TO INCREASE PARAM MAXSQT', .TRUE. )

      IF ( M .EQ. 0 )  THEN
C                             ** UPWARD RECURRENCE FOR ORDINARY
C                             ** LEGENDRE POLYNOMIALS
         DO  10  I = 1, NMU
            YLM( 0,I ) = 1.
            YLM( 1,I ) = MU( I )
10       CONTINUE
         DO  20  L = 2, TWONM1
            DO  20  I = 1, NMU
               YLM( L,I ) = ( ( 2*L-1 ) * MU(I) * YLM( L-1,I )
     $                      - ( L-1 ) * YLM( L-2,I ) ) / L
20       CONTINUE

      ELSE

         DO  30  I = 1, NMU
C                               ** Y-SUB-M-SUPER-M; DERIVED FROM
C                               ** D/A EQS. (11,12)

            YLM( M,I) = - SQT( 2*M-1 ) / SQT( 2*M )
     $                  * SQRT( 1. - MU(I)**2 ) * YLM( M-1,I )

C                              ** Y-SUB-(M+1)-SUPER-M; DERIVED FROM
C                              ** D/A EQS. (13,14) USING EQS. (11,12)

            YLM( M+1,I ) = SQT( 2*M+1 ) * MU(I) * YLM( M,I )
30       CONTINUE
C                                   ** UPWARD RECURRENCE; D/A EQ. (10)
         DO  40  L = M+2, TWONM1
            TMP1 = SQT( L-M ) * SQT( L+M )
            TMP2 = SQT( L-M-1 ) * SQT( L+M-1 )
            DO  40  I = 1, NMU
               YLM( L,I ) = ( ( 2*L-1 ) * MU(I) * YLM( L-1,I )
     $                        - TMP2 * YLM( L-2,I ) ) / TMP1
40       CONTINUE

      END IF

      RETURN
      END
