c
c----------------------------------------------------------------------
c
      SUBROUTINE  SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )
C
C         FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION 
C         AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.
C
C     INPUT:
C
C        ABD     REAL(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     EXAMPLE:  IF THE ORIGINAL MATRIX IS
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN
C
C            *  *  *  +  +  +  , * = NOT USED
C            *  * 13 24 35 46  , + = USED FOR PIVOTING
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *
C
C
C     ROUTINES CALLED:  FROM LINPACK: SGBFA
C                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
C                       FROM FORTRAN: ABS, DMAX1, MAX0, MIN0, SIGN
C
      implicit none
      INTEGER  LDA, N, ML, MU, IPVT(*)
      REAL*8     ABD(LDA,*), Z(*)
      REAL*8     RCOND
C
      REAL*8     SDOT, EK, T, WK, WKM
      REAL*8     ANORM, S, SASUM, SM, YNORM
      INTEGER  IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
C
C
C                       ** COMPUTE 1-NORM OF A
      ANORM = 0.0E0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM, SASUM(L,ABD(IS,J), 1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C                                               ** FACTOR
      CALL SGBFA(ABD, LDA, N, ML, MU, IPVT, INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C                     ** SOLVE TRANS(U)*W = E
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
C
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
         IF (ABS(EK-Z(K)) .GT. ABS(ABD(M,K))) THEN
            S = ABS(ABD(M,K))/ABS(EK-Z(K))
            CALL SSCAL(N, S, Z, 1)
            EK = S*EK
         ENDIF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (ABD(M,K) .NE. 0.0E0) THEN
            WK  = WK /ABD(M,K)
            WKM = WKM/ABD(M,K)
         ELSE
            WK  = 1.0E0
            WKM = 1.0E0
         ENDIF
         KP1 = K + 1
         JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
         MM = M
         IF (KP1 .LE. JU) THEN
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + ABS(Z(J)+WKM*ABD(MM,J))
               Z(J) = Z(J) + WK*ABD(MM,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
            ENDIF
         ENDIF
         Z(K) = WK
  100 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C
C                         ** SOLVE TRANS(L)*Y = W
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(ML, N-K)
         IF (K .LT. N) Z(K) = Z(K) + SDOT(LM, ABD(M+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0 / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
         ENDIF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C
      YNORM = 1.0E0
C                         ** SOLVE L*V = Y
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN0(ML, N-K)
         IF (K .LT. N) CALL SAXPY(LM, T, ABD(M+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0 / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
  140 CONTINUE
C
      S = 1.0E0/SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C                           ** SOLVE  U*Z = W
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .GT. ABS(ABD(M,K))) THEN
            S = ABS(ABD(M,K)) / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
         IF (ABD(M,K) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
         IF (ABD(M,K) .EQ. 0.0E0) Z(K) = 1.0E0
         LM = MIN0(K, M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL SAXPY(LM, T, ABD(LA,K), 1, Z(LZ), 1)
  160 CONTINUE
C                              ** MAKE ZNORM = 1.0
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE  SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )
C
C         FACTORS A REAL BAND MATRIX BY ELIMINATION.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     INPUT:  SAME AS 'SGBCO'
C
C     ON RETURN:
C
C        ABD,IPVT    SAME AS 'SGBCO'
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
C                       FROM FORTRAN: MAX0, MIN0
C
      implicit none
      INTEGER  LDA, N, ML, MU, IPVT(*), INFO
      REAL*8     ABD(LDA,*)
C
      REAL*8     T
      integer tmpone
      INTEGER  I,ISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C
      M = ML + MU + 1
      INFO = 0
C                        ** ZERO INITIAL FILL-IN COLUMNS
      J0 = MU + 2
      J1 = MIN0(N, M) - 1
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
      JZ = J1
      JU = 0
C
C                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      NM1 = N - 1
      DO 120 K = 1, NM1
         KP1 = K + 1
C                                  ** ZERO NEXT FILL-IN COLUMN
         JZ = JZ + 1
         IF (JZ .LE. N) THEN
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0E0
   40       CONTINUE
         ENDIF
C                                  ** FIND L = PIVOT INDEX
         LM = MIN0(ML, N-K)
         tmpone=1
         L = ISAMAX(LM+1, ABD(M,K),tmpone) + M - 1
         IPVT(K) = L + K - M
C
         IF (ABD(L,K) .EQ. 0.0E0) THEN
C                                      ** ZERO PIVOT IMPLIES THIS COLUMN 
C                                      ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
C                                ** INTERCHANGE IF NECESSARY
            IF (L .NE. M) THEN
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
            ENDIF
C                                   ** COMPUTE MULTIPLIERS
            T = -1.0E0 / ABD(M,K)
            CALL SSCAL(LM, T, ABD(M+1,K), 1)
C
C                               ** ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
            MM = M
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .NE. MM) THEN
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
               ENDIF
               CALL SAXPY(LM, T, ABD(M+1,K), 1, ABD(MM+1,J), 1)
   80       CONTINUE
C
         ENDIF
C
  120 CONTINUE
C
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE  SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )
C
C         SOLVES THE REAL BAND SYSTEM
C            A * X = B  OR  TRANSPOSE(A) * X = B
C         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     INPUT:
C
C        ABD     REAL(LDA, N)
C                THE OUTPUT FROM SBGCO OR SGBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SBGCO OR SGBFA.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
C        OR SGBFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
C                       FROM FORTRAN: MIN0
C
      implicit none
      INTEGER  LDA, N, ML, MU, IPVT(*), JOB
      REAL*8     ABD(LDA,*), B(*)
C
      REAL*8     SDOT,T
      INTEGER  K,KB,L,LA,LB,LM,M,NM1
C
C
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .EQ. 0) THEN
C                               ** JOB = 0 , SOLVE  A * X = B
C                               ** FIRST SOLVE L*Y = B
         IF (ML .NE. 0) THEN
            DO 20 K = 1, NM1
               LM = MIN0(ML, N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .NE. K) THEN
                  B(L) = B(K)
                  B(K) = T
               ENDIF
               CALL SAXPY( LM, T, ABD(M+1,K), 1, B(K+1), 1 )
   20       CONTINUE
         ENDIF
C                           ** NOW SOLVE  U*X = Y
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / ABD(M,K)
            LM = MIN0(K, M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL SAXPY(LM, T, ABD(LA,K), 1, B(LB), 1)
   40    CONTINUE
C
      ELSE
C                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                  ** FIRST SOLVE  TRANS(U)*Y = B
         DO 60 K = 1, N
            LM = MIN0(K, M) - 1
            LA = M - LM
            LB = K - LM
            T = SDOT(LM, ABD(LA,K), 1, B(LB), 1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C                                  ** NOW SOLVE TRANS(L)*X = Y
         IF (ML .NE. 0) THEN
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML, N-K)
               B(K) = B(K) + SDOT(LM, ABD(M+1,K), 1, B(K+1), 1)
               L = IPVT(K)
               IF (L .NE. K) THEN
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
               ENDIF
   80       CONTINUE
         ENDIF
C
      ENDIF
C
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE  SGECO( A, LDA, N,IPVT, RCOND, Z )
C
C         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C         AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
C         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     ROUTINES CALLED:  FROM LINPACK: SGEFA
C                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
C                       FROM FORTRAN: ABS, DMAX1, SIGN
C
      implicit none
      INTEGER  LDA, N, IPVT(*)
      REAL*8     A(LDA,*), Z(*)
      REAL*8     RCOND
C
      REAL*8     SDOT,EK,T,WK,WKM
      REAL*8     ANORM,S,SASUM,SM,YNORM
      INTEGER  INFO,J,K,KB,KP1,L
C
C
C                        ** COMPUTE 1-NORM OF A
      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = DMAX1( ANORM, SASUM(N,A(1,J),1) )
   10 CONTINUE
C                                      ** FACTOR
      CALL SGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C                        ** SOLVE TRANS(U)*W = E
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
C
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
         IF (ABS(EK-Z(K)) .GT. ABS(A(K,K))) THEN
            S = ABS(A(K,K)) / ABS(EK-Z(K))
            CALL SSCAL(N, S, Z, 1)
            EK = S*EK
         ENDIF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .NE. 0.0E0) THEN
            WK  = WK  / A(K,K)
            WKM = WKM / A(K,K)
         ELSE
            WK  = 1.0E0
            WKM = 1.0E0
         ENDIF
         KP1 = K + 1
         IF (KP1 .LE. N) THEN
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
            ENDIF
         ENDIF
         Z(K) = WK
  100 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C                                ** SOLVE TRANS(L)*Y = W
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + SDOT(N-K, A(K+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
         ENDIF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C                                 ** SOLVE L*V = Y
      YNORM = 1.0E0
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL SAXPY(N-K, T, A(K+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
  140 CONTINUE
C
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
C                                  ** SOLVE  U*Z = V
      YNORM = S*YNORM
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
            S = ABS(A(K,K))/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
         IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         T = -Z(K)
         CALL SAXPY(K-1, T, A(1,K), 1, Z(1), 1)
  160 CONTINUE
C                                   ** MAKE ZNORM = 1.0
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE  SGEFA( A, LDA, N, IPVT, INFO )
C
C         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
C
C     INPUT:  SAME AS 'SGECO'
C
C     ON RETURN:
C
C        A,IPVT  SAME AS 'SGECO'
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
C
      implicit none
      INTEGER  LDA, N, IPVT(*), INFO
      REAL*8     A(LDA,*)
C
      REAL*8     T
      integer tmpone
      INTEGER  ISAMAX,J,K,KP1,L,NM1
C
C
C                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      INFO = 0
      NM1 = N - 1
      tmpone=1
      DO 60 K = 1, NM1
         KP1 = K + 1
C                                            ** FIND L = PIVOT INDEX
         L = ISAMAX( N-K+1, A(K,K), tmpone) + K-1
         IPVT(K) = L
C
         IF (A(L,K) .EQ. 0.0E0) THEN
C                                     ** ZERO PIVOT IMPLIES THIS COLUMN 
C                                     ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
C                                     ** INTERCHANGE IF NECESSARY
            IF (L .NE. K) THEN
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
            ENDIF
C                                     ** COMPUTE MULTIPLIERS
            T = -1.0E0 / A(K,K)
            CALL SSCAL( N-K, T, A(K+1,K), 1 )
C
C                              ** ROW ELIMINATION WITH COLUMN INDEXING
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .NE. K) THEN
                  A(L,J) = A(K,J)
                  A(K,J) = T
               ENDIF
               CALL SAXPY( N-K, T, A(K+1,K), 1, A(K+1,J), 1 )
   30       CONTINUE
C
         ENDIF
C
   60 CONTINUE
C
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE  SGESL( A, LDA, N,IPVT, B, JOB )
C
C         SOLVES THE REAL SYSTEM
C            A * X = B  OR  TRANS(A) * X = B
C         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
C
      implicit none
      INTEGER  LDA, N, IPVT(*), JOB
      REAL*8     A(LDA,*), B(*)
C
      REAL*8     SDOT,T
      INTEGER  K,KB,L,NM1
C
C
      NM1 = N - 1
      IF (JOB .EQ. 0) THEN
C                                 ** JOB = 0 , SOLVE  A * X = B
C                                     ** FIRST SOLVE  L*Y = B
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .NE. K) THEN
               B(L) = B(K)
               B(K) = T
            ENDIF
            CALL SAXPY( N-K, T, A(K+1,K), 1, B(K+1), 1 )
   20    CONTINUE
C                                    ** NOW SOLVE  U*X = Y
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / A(K,K)
            T = -B(K)
            CALL SAXPY( K-1, T, A(1,K), 1, B(1), 1 )
   40    CONTINUE
C
      ELSE
C                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                    ** FIRST SOLVE  TRANS(U)*Y = B
         DO 60 K = 1, N
            T = SDOT( K-1, A(1,K), 1, B(1), 1 )
            B(K) = (B(K) - T) / A(K,K)
   60    CONTINUE
C                                    ** NOW SOLVE  TRANS(L)*X = Y
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + SDOT( N-K, A(K+1,K), 1, B(K+1), 1 )
            L = IPVT(K)
            IF (L .NE. K) THEN
               T = B(L)
               B(L) = B(K)
               B(K) = T
            ENDIF
   80    CONTINUE
C
      ENDIF
C
      RETURN
      END
c
c----------------------------------------------------------------------
c
      REAL*8 FUNCTION  SASUM( N, SX, INCX )
C
C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
C
C --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
C
      implicit none
      integer n,incx,i,m
      REAL*8 SX(*)
c,SASUM
C
C
      SASUM = 0.0
      IF( N.LE.0 )  RETURN
      IF( INCX.NE.1 ) THEN
C                                          ** NON-UNIT INCREMENTS
          DO 10 I = 1, 1+(N-1)*INCX, INCX
             SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      ELSE
C                                          ** UNIT INCREMENTS
         M = MOD(N,6)
         IF( M.NE.0 ) THEN
C                             ** CLEAN-UP LOOP SO REMAINING VECTOR 
C                             ** LENGTH IS A MULTIPLE OF 6.
            DO 30  I = 1, M
              SASUM = SASUM + ABS(SX(I))
   30       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 6
           SASUM = SASUM + ABS(SX(I))   + ABS(SX(I+1)) + ABS(SX(I+2))
     $                   + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50    CONTINUE
      ENDIF
C
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE     SAXPY( N, SA, SX, INCX, SY, INCY )
C
C          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)
C
C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
C
C --OUTPUT--
C       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH 
C                 SA*SX(LX+I*INCX) + SY(LY+I*INCY), 
C            WHERE LX = 1          IF INCX .GE. 0,
C                     = (-INCX)*N  IF INCX .LT. 0
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      implicit none
      integer N,incx,incy,i,ix,iy,m
      REAL*8 SX(*), SY(*), SA
C
C
      IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN
C
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SY(I) = SY(I) + SA * SX(I)
   10     CONTINUE
C
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
C
C                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,4)
         IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
              SY(I) = SY(I) + SA * SX(I)
   20       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 4
            SY(I)   = SY(I)   + SA * SX(I)
            SY(I+1) = SY(I+1) + SA * SX(I+1)
            SY(I+2) = SY(I+2) + SA * SX(I+2)
            SY(I+3) = SY(I+3) + SA * SX(I+3)
   30    CONTINUE
C
      ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            SY(IY) = SY(IY) + SA*SX(IX)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
C
      ENDIF
C
      RETURN
      END
c
c----------------------------------------------------------------------
c
      REAL*8 FUNCTION  SDOT( N, SX, INCX, SY, INCY )
C
C          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'
C
C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
C
C --OUTPUT--
C     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
C            WHERE  LX = 1          IF INCX .GE. 0, 
C                      = (-INCX)*N  IF INCX .LT. 0,
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      implicit none
c      real*8 sdot
      integer n,incy,incx,i,ix,iy,m
      REAL*8 SX(*), SY(*)
C
C
      SDOT = 0.0
      IF( N.LE.0 )  RETURN
C
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SDOT = SDOT + SX(I) * SY(I)
   10     CONTINUE
C
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
C
C                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,5)
         IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
               SDOT = SDOT + SX(I) * SY(I)
   20       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 5
            SDOT = SDOT + SX(I)*SY(I)     + SX(I+1)*SY(I+1)
     $                  + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3)
     $                  + SX(I+4)*SY(I+4)
   30    CONTINUE
C
      ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            SDOT = SDOT + SX(IX) * SY(IY)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
C
      ENDIF
C
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE     SSCAL( N, SA, SX, INCX )
C
C         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)
C
C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
C            SA  SINGLE PRECISION SCALE FACTOR
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
C
C --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX) 
C                FOR I = 0 TO N-1
C
      implicit none
      integer incx,n,i,m
      REAL*8 SA, SX(*)
C
C
      IF( N.LE.0 ) RETURN
C
      IF( INCX.NE.1 ) THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SX(I) = SA * SX(I)
   10     CONTINUE
C
      ELSE
C
         M = MOD(N,5)
         IF( M.NE.0 ) THEN
C                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                           ** IS A MULTIPLE OF 5.
            DO 30  I = 1, M
               SX(I) = SA * SX(I)
   30       CONTINUE
         ENDIF
C                             ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 5
            SX(I)   = SA * SX(I)
            SX(I+1) = SA * SX(I+1)
            SX(I+2) = SA * SX(I+2)
            SX(I+3) = SA * SX(I+3)
            SX(I+4) = SA * SX(I+4)
   50    CONTINUE
C
      ENDIF
C
      RETURN
      END
c
c----------------------------------------------------------------------
c
      INTEGER FUNCTION  ISAMAX( N, SX, INCX )
C
C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
C
C --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
C                         ABS(SX(1+(I-1)*INCX))
C
      implicit none
      integer n,incx,i,ii
      REAL*8 SX(*), SMAX, XMAG
C
C gg init by default
      isamax=0
      IF( N.LE.0 ) THEN
         ISAMAX = 0
      ELSE IF( N.EQ.1 ) THEN
         ISAMAX = 1
      ELSE
         SMAX = 0.0
         II = 1
         DO 20  I = 1, 1+(N-1)*INCX, INCX
            XMAG = ABS(SX(I))
            IF( SMAX.LT.XMAG ) THEN
               SMAX = XMAG
               ISAMAX = II
            ENDIF
            II = II + 1
   20    CONTINUE
      ENDIF
C
      RETURN
      END
c
c----------------------------------------------------------------------
c
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      implicit none
        integer lda,n,ipvt(*)
      real*8 a(lda,*),z(*)
      real*8 rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine dgesl(a,lda,n,ipvt,b,job)
      implicit none
        integer lda,n,ipvt(*),job
       real*8 a(lda,*),b(*)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine dgefa(a,lda,n,ipvt,info)
      implicit none
        integer lda,n,ipvt(*),info
      real*8 a(lda,*)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
      real*8 dx(*),dy(*),da
      integer i,incx,incy,m,mp1,ix,iy,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
      real*8 da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
c
c----------------------------------------------------------------------
c
      real*8 function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
c      real*8 ddot
      real*8 dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
c
c----------------------------------------------------------------------
c
      real*8 function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
      real*8 dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
c
c----------------------------------------------------------------------
c
      function idamax(n,dx,incx)
      implicit none
        integer idamax
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      real*8 dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
      implicit none
        integer lda,n,ml,mu,ipvt(*)
      real*8 abd(lda,*),z(*),rcond
c
c     dgbco factors a double precision band matrix by gaussian
c     elimination and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgbfa is slightly faster.
c     to solve  a*x = b , follow dgbco by dgbsl.
c     to compute  inverse(a)*c , follow dgbco by dgbsl.
c     to compute  determinant(a) , follow dgbco by dgbdi.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     example..  if the original matrix is
c
c           11 12 13  0  0  0
c           21 22 23 24  0  0
c            0 32 33 34 35  0
c            0  0 43 44 45 46
c            0  0  0 54 55 56
c            0  0  0  0 65 66
c
c      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abd should contain
c
c            *  *  *  +  +  +  , * = not used
c            *  * 13 24 35 46  , + = used for pivoting
c            * 12 23 34 45 56
c           11 22 33 44 55 66
c           21 32 43 54 65  *
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgbfa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,max0,min0,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer is,info,j,ju,k,kb,kp1,l,la,lm,lz,m,mm
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      l = ml + 1
      is = l + mu
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(l,abd(is,j),1))
         if (is .gt. ml + 1) is = is - 1
         if (j .le. mu) l = l + 1
         if (j .ge. n - ml) l = l - 1
   10 continue
c
c     factor
c
      call dgbfa(abd,lda,n,ml,mu,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      m = ml + mu + 1
      ju = 0
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(abd(m,k))) go to 30
            s = dabs(abd(m,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (abd(m,k) .eq. 0.0d0) go to 40
            wk = wk/abd(m,k)
            wkm = wkm/abd(m,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         ju = min0(max0(ju,mu+ipvt(k)),n)
         mm = m
         if (kp1 .gt. ju) go to 90
            do 60 j = kp1, ju
               mm = mm - 1
               sm = sm + dabs(z(j)+wkm*abd(mm,j))
               z(j) = z(j) + wk*abd(mm,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               mm = m
               do 70 j = kp1, ju
                  mm = mm - 1
                  z(j) = z(j) + t*abd(mm,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         lm = min0(ml,n-k)
         if (k .lt. n) z(k) = z(k) + ddot(lm,abd(m+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         lm = min0(ml,n-k)
         if (k .lt. n) call daxpy(lm,t,abd(m+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = w
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(abd(m,k))) go to 150
            s = dabs(abd(m,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (abd(m,k) .ne. 0.0d0) z(k) = z(k)/abd(m,k)
         if (abd(m,k) .eq. 0.0d0) z(k) = 1.0d0
         lm = min0(k,m) - 1
         la = m - lm
         lz = k - lm
         t = -z(k)
         call daxpy(lm,t,abd(la,k),1,z(lz),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      implicit none
        integer lda,n,ml,mu,ipvt(*),job
      real*8 abd(lda,*),b(*)
c
c     dgbsl solves the double precision band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran min0
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
      implicit none
        integer lda,n,ml,mu,ipvt(*),info
      real*8 abd(lda,*)
c
c     dgbfa factors a double precision band matrix by elimination.
c
c     dgbfa is usually called by dgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c     fortran max0,min0
c
c     internal variables
c
      double precision t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0d0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end
c
c----------------------------------------------------------------------
c
