      SUBROUTINE  DISORT(ien, NLYR, DTAUC, SSALB, PMOM,  USRTAU, NTAU, 
     $                    UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI,
     $                    FBEAM, UMU0, SRC, SRCU, PHI0, FISOT, LAMBER,
     $                    ALBEDO, HL, BSRC, TSRC, DELTAM, INTSRC,
     $                    ONLYFL, ACCUR, PRNT, MAXCLY, MAXULV,
     $                    MAXUMU, MAXCMU, MAXPHI, RFLDIR, RFLDN,
     $                    FLUP, DFDT, UAVG, UU, U0U, ipk, linear,
     $ ischeck )
C
C **********************************************************************
C       PLANE-PARALLEL DISCRETE ORDINATES RADIATIVE TRANSFER PROGRAM
C         MODIFIED TO  SOLVE FOR AN ANISOTROPIC INTERNAL SOURCE,
C                      ARVE KYLLING JUNE 1990
C             ( SEE DISORT.DOC FOR COMPLETE DOCUMENTATION )
*
*       Modified to fit the needs of the electron transport problem
*       by making FISOT an array, distinguish between "convex" and
*       "concave" source terms (Oystein's suggestion), and including
*       some more checks on the source term to catch overflows.  The
*       logical variable linear causes the use of a linear source
*       for the case linear=.true., and the exp-lin source for 
*       linear=.false.  Modifications are in lower case.  Some write
*       statements for debugging purposes are left as comment lines.
*       The relation between index and angle for the array SRC is
*       changed from Arve's original convention, to make it consistent
*       with the convention used for the intensities.  Note that the
*       arrays XR0 and XR1 are left as they were.
*               Dirk Lummerzheim, October 1990
*       Some arrays are in double precision in order to avoid under-
*       flows (jl, 1992).
*       The underflows problems have been cured through the xra 
*       computation which is now simplified. The d.p. process is kept, 
*       since it increases the precision and do not make any harm.
*
C **********************************************************************
C
C+---------------------------------------------------------------------+
C------------------    I/O VARIABLE SPECIFICATIONS     -----------------
C+---------------------------------------------------------------------+
C
c      IMPLICIT INTEGER (M)
      IMPLICIT NONE
      INTEGER MXCLY,MXCMU,NNLYRI,MI9M2,MXUMU,MXULV,MXPHI,MI,ipk
c       parameter(MXCLY=1001,MXULV=1001*2-1,MXCMU=32,MXUMU=32,
c     .   MXPHI=1,MI=16,MI9M2=9*16-2,NNLYRI=32*1001)
        PARAMETER ( MXCLY =1001, MXULV =mxcly*2-1, MXCMU = 32, 
     .            MXUMU = mxcmu,
     $            MXPHI = 1, MI = MXCMU/2, MI9M2 = 9*MI-2,
     $            NNLYRI = MXCMU*MXCLY )


      LOGICAL  DELTAM, LAMBER, INTSRC, ONLYFL, PRNT(7), USRANG, USRTAU
      LOGICAL ischeck
      INTEGER  MAXCLY, MAXUMU, MAXULV, MAXCMU, MAXPHI, NLYR,
     $         NUMU, NSTR, NPHI, NTAU
      REAL*8  ACCUR, ALBEDO, BSRC, DTAUC(NTAU), FBEAM, 
     .         fisot(-ipk:-1),
     $         HL( 0:MAXCMU ), PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ),
     $         PHI0, SSALB( MAXCLY ), SRC( 3*MAXCLY, MAXCMU ), 
     $         SRCU( 3*MAXCLY, MAXUMU), TSRC,
     $         UMU( MAXUMU ), UMU0, UTAU( MAXULV )
C
      REAL*8  RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     $         UAVG( MAXULV ), DFDT( MAXULV ), U0U( MAXUMU, MAXULV ),
     $         UU( MAXUMU, MAXULV, MAXPHI )
C
C+---------------------------------------------------------------------+
C      ROUTINES CALLED (IN ORDER):  ZEROAL, CHEKIN, SETDIS, PRTINP,
C                                   LEPOLY, SURFAC, SOLEIG, UPBEAM,
C                                   UPISOT, TERPEV, TERPSO, SETMTX,
C                                   SOLVE0, FLUXES, USRINT, PRAVIN,
C                                   PRTINT
C+---------------------------------------------------------------------+
C
C  INDEX CONVENTIONS (FOR ALL DO-LOOPS AND ALL VARIABLE DESCRIPTIONS):
C
C     IU     :  FOR USER POLAR ANGLES
C
C  IQ,JQ,KQ  :  FOR COMPUTATIONAL POLAR ANGLES ('QUADRATURE ANGLES')
C
C   IQ/2     :  FOR HALF THE COMPUTATIONAL POLAR ANGLES (JUST THE ONES
C               IN EITHER 0-90 DEGREES, OR 90-180 DEGREES)
C
C     J      :  FOR USER AZIMUTHAL ANGLES
C
C     K,L    :  FOR LEGENDRE EXPANSION COEFFICIENTS OR, ALTERNATIVELY,
C               SUBSCRIPTS OF ASSOCIATED LEGENDRE POLYNOMIALS
C
C     LU     :  FOR USER LEVELS
C
C     LC     :  FOR COMPUTATIONAL LAYERS (EACH HAVING A DIFFERENT
C               SINGLE-SCATTER ALBEDO AND/OR PHASE FUNCTION)
C
C    LEV     :  FOR COMPUTATIONAL LEVELS
C
C    MAZ     :  FOR AZIMUTHAL COMPONENTS IN FOURIER COSINE EXPANSION
C               OF INTENSITY AND PHASE FUNCTION
C
C+---------------------------------------------------------------------+
C               I N T E R N A L    V A R I A B L E S
C
C   AMB(IQ/2,IQ/2)    FIRST MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEM
C                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN 'SOLEIG')
C
C   APB(IQ/2,IQ/2)    SECOND MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEM
C                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN 'SOLEIG')
C
C   ARRAY(IQ,IQ)      SCRATCH MATRIX FOR 'SOLEIG', 'UPBEAM' AND 'UPISOT'
C                     (SEE EACH SUBROUTINE FOR DEFINITION)
C
C   B()               RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO
C                     *SOLVE0,1*;  RETURNS AS SOLUTION VECTOR
C                     VECTOR CAPITAL-L, THE CONSTANTS OF INTEGRATION
C
C   BDR(IQ/2,0:IQ/2)  BOTTOM-BOUNDARY BIDIRECTIONAL REFLECTIVITY FOR A
C                     GIVEN AZIMUTHAL COMPONENT.  FIRST INDEX ALWAYS
C                     REFERS TO A COMPUTATIONAL ANGLE.  SECOND INDEX:
C                     IF ZERO, REFERS TO INCIDENT BEAM ANGLE -UMU0-;
C                     IF NON-ZERO, REFERS TO A COMPUTATIONAL ANGLE.
C
C   BEM(IQ/2)         BOTTOM-BOUNDARY DIRECTIONAL EMISSIVITY AT COMPU-
C                     TATIONAL ANGLES.
C
C   CBAND()           MATRIX OF LEFT-HAND SIDE OF THE LINEAR SYSTEM
C                     EQ. SC(5), SCALED BY EQ. SC(12);  IN BANDED
C                     FORM REQUIRED BY LINPACK SOLUTION ROUTINES
C
C   CC(IQ,IQ)         CAPITAL-C-SUB-IJ IN EQ. SS(5)
C
C   CMU(IQ)           COMPUTATIONAL POLAR ANGLES (GAUSSIAN)
C
C   CWT(IQ)           QUADRATURE WEIGHTS CORRESP. TO -CMU-
C
C   DELM0             KRONECKER DELTA, DELTA-SUB-M0, WHERE 'M' = MAZ
C                     IS THE NUMBER OF THE FOURIER COMPONENT IN THE
C                     AZIMUTH COSINE EXPANSION
C
C   EMU(IU)           BOTTOM-BOUNDARY DIRECTIONAL EMISSIVITY AT USER
C                     ANGLES.
C
C   EVAL(IQ)          TEMPORARY STORAGE FOR EIGENVALUES OF EQ. SS(12)
C
C   EVECC(IQ,IQ)      COMPLETE EIGENVECTORS OF SS(7) ON RETURN FROM
C                     *SOLEIG* ; STORED PERMANENTLY IN -GC-
C
C   EXPBEA(LC)        TRANSMISSION OF DIRECT BEAM IN DELTA-M OPTICAL
C                     DEPTH COORDINATES
C
C   FLYR(LC)          TRUNCATED FRACTION IN DELTA-M METHOD
C
C   GL(K,LC)          PHASE FUNCTION LEGENDRE POLY. EXPANSION
C                     COEFFICIENTS, CALCULATED FROM 'PMOM' BY
C                     INCLUDING SINGLE-SCATTERING ALBEDO, FACTOR
C                     2K+1, AND (IF DELTAM=TRUE) THE DELTA-M
C                     SCALING
C
C   GC(IQ,IQ,LC)      EIGENVECTORS AT POLAR QUADRATURE ANGLES,
C                     LITTLE-G  IN EQ. SC(1)
C
C   GU(IU,IQ,LC)      EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
C                     ( LITTLE-G  IN EQS. SC(3) AND S1(8-9), i.e.
C                     CAPITAL-G WITHOUT THE CAPITAL-L FACTOR )
C
C   HLPR()            LEGENDRE COEFFICIENTS OF BOTTOM BIDIRECTIONAL
C                     REFLECTIVITY (AFTER INCLUSION OF 2K+1 FACTOR)
C
C   IPVT(LC*IQ)       INTEGER VECTOR OF PIVOT INDICES FOR LINPACK
C                     ROUTINES
C
C   KK(IQ,LC)         EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C
C   KCONV             COUNTER IN AZIMUTH CONVERGENCE TEST
C
C   LAYRU(LU)         COMPUTATIONAL LAYER IN WHICH USER OUTPUT LEVEL
C                     -UTAU(LU)- IS LOCATED
C
C   LL(IQ,LC)         CONSTANTS OF INTEGRATION CAPITAL-L IN EQ. SC(1),
C                     OBTAINED BY SOLVING SCALED VERSION OF EQ. SC(5)
C
C   LYRCUT            TRUE, RADIATION IS ASSUMED ZERO BELOW LAYER
C                     -NCUT- BECAUSE OF ALMOST COMPLETE ABSORPTION
C
C   NAZ               NUMBER OF AZIMUTHAL COMPONENTS CONSIDERED
C
C   NCUT              COMPUTATIONAL LAYER NUMBER IN WHICH ABSORPTION
C                     OPTICAL DEPTH FIRST EXCEEDS -ABSCUT-
C
C   OPRIM(LC)         SINGLE SCATTERING ALBEDO AFTER DELTA-M SCALING
C
C   PASS1             TRUE ON FIRST ENTRY, FALSE THEREAFTER
C
C   PSIO(IQ),         SUM JUST AFTER SQUARE BRACKET IN  EQ. SD(9)
C    PSI1(IQ)         FOR Z0 AND Z1 RESPECTIVELY
C
C   RMU(IU,0:IQ)      BOTTOM-BOUNDARY BIDIRECTIONAL REFLECTIVITY FOR A
C                     GIVEN AZIMUTHAL COMPONENT.  FIRST INDEX ALWAYS
C                     REFERS TO A USER ANGLE.  SECOND INDEX:
C                     IF ZERO, REFERS TO INCIDENT BEAM ANGLE -UMU0-;
C                     IF NON-ZERO, REFERS TO A COMPUTATIONAL ANGLE.
C
C   TAUC(0:LC)        CUMULATIVE OPTICAL DEPTH (UN-DELTA-M-SCALED)
C
C   TAUCPR(0:LC)      CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED IF
C                     DELTAM = TRUE, OTHERWISE EQUAL TO -TAUC-)
C
C   UUM(IU,LU,MAZ)    COMPONENTS OF THE INTENSITY (U-SUPER-M) WHEN
C                     EXPANDED IN FOURIER COSINE SERIES IN AZIMUTH ANGLE
C
C   U0C(IQ,LU)        AZIMUTHALLY-AVERAGED INTENSITY
C
C   UTAUPR(LU)        OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
C                     COORDINATES;  EQUAL TO  -UTAU(LU)- IF NO DELTA-M
C
C   WK()              SCRATCH ARRAY
C
C   XR0(LC)           X-SUB-ZERO IN EXPANSION OF THERMAL SOURCE FUNC-
C                     TION PRECEDING EQ. SS(14) (HAS NO MU-DEPENDENCE)
C
C   XR1(LC)           X-SUB-ONE IN EXPANSION OF THERMAL SOURCE FUNC-
C                     TION;  SEE  EQS. SS(14-16)
C
C   XRA(LC)           X-SUB-TWO IN EXPANSION OF THERMAL SOURCE FUNC-
C                     TION;  SEE  EQS. SS(14-16)
C
C   YLM0(L)           NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                     OF SUBSCRIPT 'L' AT THE BEAM ANGLE (NOT SAVED
C                     AS FUNCTION OF SUPERSCIPT 'M')
C
C   YLMC(L,IQ)        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                     OF SUBSCRIPT 'L' AT THE COMPUTATIONAL ANGLES
C                     (NOT SAVED AS FUNCTION OF SUPERSCIPT 'M')
C
C   YLMU(L,IU)        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                     OF SUBSCRIPT 'L' AT THE USER ANGLES
C                     (NOT SAVED AS FUNCTION OF SUPERSCIPT 'M')
C
C   Z()               SCRATCH ARRAY USED IN * SOLVE0,1* TO SOLVE A
C                     LINEAR SYSTEM FOR THE CONSTANTS OF INTEGRATION
C
C   Z0(IQ)            SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)
C
C   Z0U(IU,LC)        Z-SUB-ZERO IN EQ. SS(16) INTERPOLATED TO USER
C                     ANGLES FROM AN EQUATION DERIVED FROM SS(16)
C
C   Z1(IQ)            SOLUTION VECTORS Z-SUB-ONE  OF EQ. SS(16)
C
C   ZA(LC)            THE ALFA COEFFICIENT AT EACH LAYER, EQ.
C
C   Z1U(IU,LC)        Z-SUB-ONE IN EQ. SS(16) INTERPOLATED TO USER
C                     ANGLES FROM AN EQUATION DERIVED FROM SS(16)
C
C   ZBEAM(IU,LC)      PARTICULAR SOLUTION FOR BEAM SOURCE
C
C   ZJ(IQ)            RIGHT-HAND SIDE VECTOR CAPITAL-X-SUB-ZERO IN
C                     EQ. SS(19), ALSO THE SOLUTION VECTOR CAPITAL
C                     -Z-SUB-ZERO AFTER SOLVING THAT SYSTEM
C
C   ZZ(IQ,LC)         PERMANENT STORAGE FOR THE BEAM SOURCE VECTORS -ZJ-
C
C   ZPLK0(IQ,LC)      PERMANENT STORAGE FOR THE INTERNAL SOURCE
C                     VECTORS  -Z0-  OBTAINED BY SOLVING  EQ. SS(16)
C
C   ZPLK1(IQ,LC)      PERMANENT STORAGE FOR THE INTERNAL SOURCE
C                     VECTORS  -Z1-  OBTAINED BY SOLVING  EQ. SS(16)
C
C   ZPLKA(LC)         PERMANENT STORAGE FOR THE INTERNAL SOURCE
C                     VECTORS  -ZA-  IN EQ. SS(16)
C
C+---------------------------------------------------------------------+
C   LOCAL SYMBOLIC DIMENSIONS:
C
C       MXCLY  = MAX NO. OF COMPUTATIONAL LAYERS
C       MXULV  = MAX NO. OF OUTPUT LEVELS
C       MXCMU  = MAX NO. OF COMPUTATION POLAR ANGLES
C       MXUMU  = MAX NO. OF OUTPUT POLAR ANGLES
C       MXPHI  = MAX NO. OF OUTPUT AZIMUTHAL ANGLES
C+---------------------------------------------------------------------+
        INTEGER NN
c      INTEGER MXCLY,ien,ipk,iounit,I,NNLYRI,MXULV,MXCMU,MXUMU,NXPHI
c      INTEGER MI,MI9M2
C
      REAL*8 tempsqrt,templog,pi,SGN,rpd,epsil,exptest,AZTERM,AZERR
      REAL*8 COSPHI,DELTAT,exxp,DELM0,RATIO,ABS

      INTEGER L,KCONV,NCUT
      INTEGER I,J,LU,LS,IQ,LC,MAZ,NAZ,NCOL,IU,ien
      LOGICAL LYRCUT, PASS1,pass2,linear
      INTEGER IPVT( NNLYRI ), LAYRU( MXULV )
CAKYB
      REAL*8    AMB( MI,MI ), APB( MI,MI ), ARRAY( MXCMU,MXCMU ),  
     $        B( NNLYRI ), BDR( MI,0:MI ), BEM( MI ),  
     $        CBAND( MI9M2,NNLYRI ), CC( MXCMU,MXCMU ), CMU( MXCMU ),  
     $        CWT( MXCMU ), EMU( MXUMU ), EVAL( MI ),
     $        EVECC( MXCMU, MXCMU ), EXPBEA( 0:MXCLY ), FLYR( MXCLY ),
     $        FLDN( MXULV ), FLDIR( MXULV ), GL( 0:MXCMU,MXCLY ),
     $        GC( MXCMU,MXCMU,MXCLY ), GU( MXUMU,MXCMU,MXCLY ),
     $        HLPR( 0:MXCMU ), KK( MXCMU,MXCLY ), LL( MXCMU,MXCLY ),
     $        OPRIM( MXCLY ), PHIRAD( MXPHI ), PSI0( MXCMU ),  
     $        PSI1(MXCMU), RMU( MXUMU,0:MI ), TAUC( 0:MXCLY ), 
     $        TAUCPR( 0:MXCLY ), U0C( MXCMU,MXULV ), UTAUPR( MXULV ),
     $        UUM( MXUMU,MXULV,0:MXCMU ), WK( MXCMU ), 
     $        XR0( MXCMU, MXCLY ),XR1( MXCMU,MXCLY ), XRA( MXCLY), 
     $        XR0U( MXUMU, MXCLY ),XR1U( MXUMU,MXCLY ), XRAU( MXCLY), 
     $        YLM0( 0:MXCMU ), YLMC( 0:MXCMU,MXCMU ),
     $        YLMU( 0:MXCMU,MXUMU ), Z( NNLYRI ), Z0( MXCMU ), 
     $        Z0U( MXUMU,MXCLY ), Z1( MXCMU ),  ZA(MXCLY), 
     $        Z1U( MXUMU,MXCLY ), ZJ( MXCMU ), 
     $        ZZ( MXCMU,MXCLY ), ZPLK0( MXCMU,MXCLY ),
     $        ZPLK1( MXCMU,MXCLY ), ZPLKA(MXCLY), ZBEAM( MXUMU,MXCLY ) 
C
      real*8  dARRAY( MXCMU,mxcmu ),dWK(mxcmu),
     .          dzj(mxcmu),dz0(mxcmu),dz1(mxcmu)
      real*8 dCBAND( MI9M2,NNLYRI ),dZ( NNLYRI ),dB( NNLYRI )
c
      real*8   AAD( MI,MI ), EVALD( MI ) , EVECCD( MI,MI ),
     $                   WKD( MXCMU )
C
      SAVE  PASS1, pi, epsil, rpd, exptest
      DATA  PASS1 / .TRUE. /pass2/.true./
C
C

c      write(*,*)'testgg2 SRC(8,7)',SRC(8,7),SRC(89,7)
c      print*,"pmom=",PMOM
      IF ( PASS1 )  THEN
         PI = 2. * ASIN(1.0)
         EPSIL = 1e-5           ! 10.*R1MACH(3)
         RPD = PI / 180.0
*         exptest=log(2.)*2**13         ! for Cray Y-MP
         exptest=log(2.)*2**6           ! for SUN (single precision)
         PASS1 = .FALSE.
      END IF

c      print*,'fisot',FISOT
cc      print*,"FISOT"
c     do tmpgg=-ipk,-1
c      print*,FISOT(tmpgg)
c      enddo
C
c   1  CONTINUE
c      IF ( PRNT(1) )  WRITE( 6,1010 )  HEADER
C
C                         ** ZERO SOME ARRAYS (NOT STRICTLY NECESSARY,
C                         ** BUT OTHERWISE UNUSED PARTS OF ARRAYS
C                         ** COLLECT GARBAGE)
      DO 10 I = 1, NNLYRI
         IPVT(I) = 0
10    CONTINUE
      CALL  ZEROAL( AMB, APB, ARRAY, CC, CMU, CWT, EVAL, EVECC,
     $              GC, GU, HLPR, KK, LL, PSI0, PSI1, WK, XR0, XR1, XRA,
     $              XR0U, XR1U, XRAU, YLM0, YLMC, YLMU, Z, Z0, Z1, ZA, 
     $              ZJ, ZZ, ZPLK0, ZPLK1, ZPLKA, Z0U, Z1U, ZBEAM, 
     $              MI, MXCMU, MXCLY, NNLYRI, MXUMU )
      call dzeroit (dWK,mxcmu)
      call dzeroit (dzj,mxcmu)
      call dzeroit (dz0,mxcmu)
      call dzeroit (dz1,mxcmu)
      call dzeroit (dZ,NNLYRI)
      CALL  dZEROIT(darray,MXCMU*mxcmu)
      call dzeroit (dCBAND,MI9M2*NNLYRI )
C
C                                  ** CALCULATE CUMULATIVE OPTICAL DEPTH
C                                  ** AND DITHER SINGLE-SCATTER ALBEDO
C                                  ** TO IMPROVE NUMERICAL BEHAVIOR OF
C                                  ** EIGENVALUE/VECTOR COMPUTATION
      TAUC( 0 ) = 0.
      CALL  ZEROIT( TAUC(0), MXCLY+1 )
      DO 20  LC = 1, NLYR
*         IF( SSALB(LC).EQ.1.0 )  SSALB(LC) = 1.0 - EPSIL
         ssalb(lc) = min(ssalb(lc),1.0-epsil)
         TAUC(LC) = TAUC(LC-1) + DTAUC(LC)
20    CONTINUE
C                                ** CHECK INPUT DIMENSIONS AND VARIABLES
C
         if(ischeck) then
      CALL  CHEKIN( NLYR, DTAUC, SSALB, PMOM, 
     $              SRC, USRTAU, NTAU, UTAU, NSTR, USRANG,
     $              NUMU, UMU, NPHI, PHI, FBEAM, UMU0,
     $              PHI0, FISOT, LAMBER, ALBEDO, HL, BSRC,
     $              TSRC, INTSRC, ONLYFL, ACCUR, MAXCLY,
     $              MAXULV, MAXUMU, MAXCMU, MAXPHI, MXCLY,
     $              MXULV,  MXUMU,  MXCMU,  MXPHI, TAUC, ipk )
        endif
C
C                                 ** PERFORM VARIOUS SETUP OPERATIONS
C
      CALL  SETDIS( CMU, CWT, DELTAM, DTAUC, EXPBEA, FBEAM, FLYR,
     $              GL, HL, HLPR, LAMBER, LAYRU, LYRCUT, MAXCLY,
     $              MAXUMU, MAXCMU, mxcly, MXCMU, NCUT, NLYR, NTAU, NN,
     $              NSTR, INTSRC, NUMU, ONLYFL, OPRIM, PMOM, SRC, 
     $              SRCU, SSALB, TAUC, TAUCPR, UTAU, UTAUPR, UMU, 
     $              UMU0, USRTAU, USRANG)
C
C                                             ** PRINT INPUT INFORMATION
      IF ( PRNT(1) )
     $     CALL PRTINP( NLYR, DTAUC, SSALB, PMOM, NTAU, UTAU, NSTR, 
     $                  NUMU, UMU, NPHI, PHI, FBEAM, UMU0, PHI0,
     $                  FISOT, LAMBER, ALBEDO, HL, BSRC, TSRC,
     $                  DELTAM, INTSRC, ONLYFL, ACCUR, FLYR, LYRCUT, 
     $                  OPRIM, TAUC, TAUCPR, MAXCMU, PRNT(7), ipk )
C
      IF ( .NOT. INTSRC )  THEN
         BSRC = 0.0
         TSRC = 0.0
         CALL  ZEROIT( SRC, 3*MAXCLY*MAXCMU )
      END IF
C
C
C ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
C ========  (EQ STWJ 5)
C
      KCONV = 0
      NAZ = NSTR-1
C                                            ** AZIMUTH-INDEPENDENT CASE
C
      IF ( FBEAM.EQ.0.0 .OR. (1.-UMU0).LT.1.E-5 .OR. ONLYFL .OR.
     $     (NUMU.EQ.1.AND.(1.-UMU(1)).LT.1.E-5 ) )
     $   NAZ = 0
C
      CALL  ZEROIT( UU, MAXUMU*MAXULV*MAXPHI )
      DO  200  MAZ = 0, NAZ
C
      IF ( MAZ.EQ.0 )  DELM0 = 1.0
      IF ( MAZ.GT.0 )  DELM0 = 0.0
C                                  ** GET NORMALIZED ASSOCIATED LEGENDRE
C                          ** POLYNOMIALS FOR INCIDENT BEAM ANGLE COSINE
      IF ( FBEAM.GT.0.0 )
     $     CALL  LEPOLY( 1, MAZ, MXCMU, NSTR-1, -UMU0, YLM0)
C
C                                  ** GET NORMALIZED ASSOCIATED LEGENDRE
C                                      ** POLYNOMIALS FOR COMPUTATIONAL
C                                      ** AND USER POLAR ANGLE COSINES
      IF ( .NOT.ONLYFL .AND. USRANG )
     $ CALL  LEPOLY( NUMU, MAZ, MXCMU, NSTR-1, UMU, YLMU )
       CALL  LEPOLY( NN,   MAZ, MXCMU, NSTR-1, CMU, YLMC )
C
C                       ** EVALUATE NORMALIZED ASSOCIATED LEGENDRE
C                       ** POLYNOMIALS WITH NEGATIVE -CMU- FROM THOSE
C                       ** WITH POSITIVE -CMU-; DAVE/ARMSTRONG EQ. (15)
      SGN  = - 1.0
      DO  50  L = MAZ, NSTR-1
         SGN = - SGN
         DO  50  IQ = NN+1, NSTR
            YLMC( L,IQ ) = SGN * YLMC( L,IQ-NN )
 50   CONTINUE
C                                 ** SPECIFY USER'S BOTTOM REFLECTIVITY
C                                 ** AND EMISSIVITY PROPERTIES
      IF ( .NOT.LYRCUT )
     $   CALL  SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER,
     $                 MI, MAZ, MXCMU, MXUMU, NN, NUMU, NSTR, ONLYFL,
     $                 UMU, USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM,
     $                 RMU)
C
C ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============
C
c       if (ien.le.18 .and. ien.ge.16)then
c       if (ien.eq.17)then
c       write(47,*)
c       write(47,*)'ien ',ien
c       write(47,*) 'rcond in UPISOT, (xr1(i),z1(i),i=1,nn)'
c       endif
      DO 100  LC = 1, NCUT
C
C                        ** SOLVE EIGENFUNCTION PROBLEM IN EQ. STWJ(8B);
C                        ** RETURN EIGENVALUES AND EIGENVECTORS
C
         CALL  SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL(0,LC), MI, MAZ,
     $                 MXCMU, NN, NSTR,  YLMC, CC, EVECC, EVAL,
     $                 KK(1,LC), GC(1,1,LC),  AAD, WKD, EVECCD,
     .                 EVALD )
C
C                                  ** CALCULATE PARTICULAR SOLUTIONS OF
C                                  ** EQ.SS(18) FOR INCIDENT BEAM SOURCE
         IF ( FBEAM.GT.0.0 )
     $        CALL  UPBEAM( CC, CMU, DELM0, FBEAM, GL(0,LC),
     $                   IPVT, MAZ, MXCMU, NN, NSTR, PI, UMU0, 
     $                   YLM0, YLMC, ZZ(1,LC),darray,dwk,dzj)
C
         IF ( INTSRC .AND. MAZ.EQ.0 ) THEN
C
C                              ** CALCULATE PARTICULAR SOLUTIONS OF
C                              ** EQ. SS(15) FOR INTERNAL EMISSION SOURCE
C
            DELTAT = TAUCPR(LC) - TAUCPR(LC-1)
            XRA( LC ) = 0.0
C                              ** XRA = ALFA NOT ANGLE  DEPENDENT
C
            IQ = 1
            LS = 3*(LC-1)+1
            IF ( DELTAT.GT.0.0 ) THEN
               if( src(ls+2,iq) .gt. 0.) then
                   tempsqrt=(src(ls+1,iq)/src(ls+2,iq))**2 -
     .                  src(ls,iq)/src(ls+2,iq)
                   templog=src(ls+1,iq)/src(ls+2,iq)
                   if( tempsqrt  .gt. 0.0         .and.
     .                 src(lc,iq).lt.src(ls+2,iq)) then
                                xra(lc) = (2./deltat) *
     .                                  log( templog + sqrt(tempsqrt))
            elseif(tempsqrt.gt.0.0.and.src(lc,iq).gt.src(ls+2,iq)) then
                         if(templog.gt.sqrt(tempsqrt)) then
                    xra(lc) = (2./deltat)*log( templog - sqrt(tempsqrt))
                          else
                              xra(lc)=0.
                         endif
                   else
                            xra(lc)=0.
                   endif
               ELSEIF ( SRC(LS+2,IQ).EQ.0.0 ) THEN
                   xra(lc)=0.
               ENDIF
            ENDIF
            if(linear) xra(lc)=0.               ! test with alpha=0
            if(xra(lc)*taucpr(lc).ge.exptest.or.
     .         xra(lc)*taucpr(lc-1).ge.exptest) then
*                 write(6,*) 'XRA is too large ',iq,lc,xra(lc),
*     .                 taucpr(lc),taucpr(lc-1)
                  xra(lc)=0.
            end if
            DO 110 IQ = 1, nn
               XR1( IQ,LC ) = 0.0
               xr1( iq+nn,lc ) = 0.0
               IF ( DELTAT.GT.0.0 ) THEN
                  XR1(IQ,LC) = (1./DELTAT)*
     $                   ( SRC(LS+2,IQ+nn)*EXxP( XRA(LC)*TAUCPR(LC) ) -
     $                     SRC(LS,IQ+nn)*EXxP( XRA(LC)*TAUCPR(LC-1) ) )
c47      format(7(1pe10.2))
                  xr1(iq+nn,lc) = (1./deltat)*
     $                   ( src(ls+2,nn+1-iq)*exp( xra(lc)*taucpr(lc) ) -
     $                     src(ls,nn+1-iq)*exp( xra(lc)*taucpr(lc-1) ) )
               else
                write(6,*) 'deltat = ',deltat,' !'
                write(6,*) 'Programme stoppe a #E =', ien
                stop
               ENDIF
              XR0( IQ, LC ) = SRC(LS,IQ+nn)*EXxP( XRA(LC)*TAUCPR(LC-1) )
     $                         - XR1(IQ,LC)*TAUCPR(LC-1)
               xr0(iq+nn,lc) = src(ls,nn+1-iq)*exp(xra(lc)*taucpr(lc-1))
     $                         - xr1(iq+nn,lc)*taucpr(lc-1)
*           halftau=abs(taucpr(lc)+taucpr(lc-1))/2.
*           qdumi(iq+nn)=exp(-xra(lc)*halftau)*
*     .         (xr1(iq,lc)*halftau+xr0(iq,lc))
*           qdumi(1-iq+nn)=exp(-xra(lc)*halftau)*
*     .         (xr1(iq+nn,lc)*halftau+xr0(iq+nn,lc))
110         CONTINUE
*       write(6,10003) 'qcenter ',lc,halftau,
*     .         (qdumi(iq),iq=1,nstr)
*       write(6,10003) 'top     ',ls,taucpr(lc-1),
*     .         (src(ls,iq),iq=1,nstr)
*       write(6,10003) 'center  ',ls,halftau,
*     .         (src(ls+1,iq),iq=1,nstr)
*       write(6,10003) 'bottom  ',ls,taucpr(lc),
*     .         (src(ls+2,iq),iq=1,nstr)
c 10003   format(1x,a,i5,f9.4,(t25,8(1pe12.3)))
C
         IF ( USRANG ) THEN
            XRAU(LC) = XRA(LC)
            DO 111 IU = 1, NUMU
               XR1U( IU, LC ) = 0.0
               IF ( DELTAT.GT.0.0 ) THEN
                  XR1U(IU,LC) = (1./DELTAT)*
     $                  ( SRCU(LS+2,IU)*EXxP( XRAU(LC)*TAUCPR(LC) )-
     $                    SRCU(LS,IU)*EXxP( XRAU(LC)*TAUCPR(LC-1) )  )
               ENDIF
               XR0U( IU, LC ) =SRCU(LS,IU)*EXxP( XRAU(LC)*TAUCPR(LC-1) )
     $                          - XR1U(IU,LC)*TAUCPR(LC-1)
111         CONTINUE
         ENDIF
C
            CALL UPISOT( CC, CMU, IPVT, MXCMU, NN, NSTR, 
     $                    XR0(1,LC), XR1(1,LC), XRA(LC), 
     $                    ZA(lc), ZPLK0(1,LC), ZPLK1(1,LC), 
     $                   ZPLKA(LC),dARRAY,dWK,dz0,dz1)
         END IF
C
C
         IF ( .NOT.ONLYFL .AND. USRANG ) THEN
C                                            ** INTERPOLATE EIGENVECTORS
C                                            ** TO USER ANGLES
C
            CALL  TERPEV( CWT, EVECC, GL(0,LC), GU(1,1,LC), MAZ, MXCMU,
     $                    MXUMU, NN, NSTR, NUMU, WK, YLMC, YLMU)
C
C                                            ** INTERPOLATE SOURCE TERMS
C                                            ** TO USER ANGLES
C
            CALL  TERPSO( CWT, DELM0, FBEAM, GL(0,LC), MAZ,
     $                    MXCMU, INTSRC, NUMU, NSTR, 
     $                    PI, YLM0, YLMC, YLMU, PSI0, PSI1, 
     $                    XR0U(1,LC), XR1U(1,LC), Z0, Z1, ZJ, 
     $                    ZBEAM(1,LC), Z0U(1,LC), Z1U(1,LC) )   
         END IF
C
100   CONTINUE
C
C ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============
C
C                      ** SET COEFFICIENT MATRIX OF EQUATIONS COMBINING
C                      ** BOUNDARY AND LAYER INTERFACE CONDITIONS
C
      CALL  SETMTX( BDR, CBAND, CMU, CWT, DELM0, GC, KK, LAMBER,
     $              LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI,
     $              NN, NSTR, TAUCPR, WK )
C
C                      ** SOLVE FOR CONSTANTS OF INTEGRATION IN HOMO-
C                      ** GENEOUS SOLUTION (GENERAL BOUNDARY CONDITIONS)
C
      CALL  SOLVE0( B, BDR, BEM, BSRC, CBAND, CMU, CWT, EXPBEA,
     $              FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT,
     $              MAZ, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR,
     $              NNLYRI, PI, TSRC, TAUCPR, UMU0, ZZ,
     $              ZPLK0, ZPLK1, ZPLKA, ipk, dcband,dz,db)
C
C                                  ** COMPUTE UPWARD AND DOWNWARD FLUXES
      IF ( MAZ.EQ.0 )
     $     CALL FLUXES( CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
     $                  MXCMU, MXULV, NCUT, NN, NSTR, NTAU, PI,
     $                  PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,
     $                  XR0, XR1, XRA, ZZ, ZPLK0, ZPLK1, ZPLKA, 
     $                  DFDT, FLUP, FLDN, FLDIR, RFLDIR, RFLDN, 
     $                  UAVG, U0C, MAXULV)
C
      IF ( ONLYFL )  THEN
         IF( MAXUMU.GE.NSTR )  THEN
C                                         ** SAVE AZIM-AVGD INTENSITIES
C                                         ** AT QUADRATURE ANGLES
            DO 120 LU = 1, NTAU
               DO 120 IQ = 1, NSTR
                  U0U( IQ,LU ) = U0C( IQ,LU )
120         CONTINUE
         ELSE
               CALL  ZEROIT( U0U, MAXUMU*MAXULV )
         ENDIF
         GO TO 210
      ENDIF
C
      IF ( USRANG ) THEN
C                                     ** COMPUTE AZIMUTHAL INTENSITY
C                                     ** COMPONENTS AT USER ANGLES
C
*         CALL  USRINT( BSRC, CMU, CWT, DELM0, EMU, EXPBEA,
*     $                 FBEAM, FISOT, GC, GU, KK, LAMBER, LAYRU, LL,
*     $                 LYRCUT, MAZ, MXCMU, MXULV, MXUMU, NCUT,
*     $                 NLYR, NN, NSTR, INTSRC, NUMU, NTAU, PI, RMU,
*     $                 TAUCPR, TSRC, UMU, UMU0, UTAUPR, WK, ZBEAM,  
*     $                 Z0U, Z1U, ZA, ZZ, ZPLK0, ZPLK1, ZPLKA, UUM)
        print*,'>>> DISORT:  USRINT = TRUE is not allowed'
C
      ELSE
C                                     ** COMPUTE AZIMUTHAL INTENSITY
C                                     ** COMPONENTS AT QUADRATURE ANGLES
C
         CALL  CMPINT( FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZ,
     $                 MXCMU, MXULV, MXUMU, NCUT, NN, NSTR,
     $                 INTSRC, NTAU, TAUCPR, UMU0, UTAUPR,
     $                 ZZ, ZPLK0, ZPLK1, ZPLKA, UUM )
      END IF
C
      IF( MAZ.EQ.0 ) THEN
C
         DO  140  J = 1, NPHI
            PHIRAD( J ) = RPD * ( PHI(J) - PHI0 )
 140     CONTINUE
C                               ** SAVE AZIMUTHALLY AVERAGED INTENSITIES
         DO 160  LU = 1, NTAU
            DO 160  IU = 1, NUMU
               U0U( IU,LU ) = UUM( IU,LU,0 )
 160     CONTINUE
C                              ** PRINT AZIMUTHALLY AVERAGED INTENSITIES
C                              ** AT USER ANGLES
         IF ( PRNT(4) )
     $        CALL PRAVIN( UMU, NUMU, MAXUMU, UTAU, NTAU, U0U)
C
      END IF
C                                ** INCREMENT INTENSITY BY CURRENT
C                                ** AZIMUTHAL COMPONENT (FOURIER
C                                ** COSINE SERIES);  EQ SD(2)
      AZERR = 0.0
      DO 180  J = 1, NPHI
         COSPHI = COS( MAZ * PHIRAD(J) )
         DO 180  LU = 1, NTAU
            DO 180  IU = 1, NUMU
               AZTERM = UUM( IU,LU,MAZ ) * COSPHI
               UU( IU,LU,J ) = UU( IU,LU,J ) + AZTERM
               AZERR = DMAX1( RATIO( ABS(AZTERM), ABS(UU(IU,LU,J)) ),
     .                   AZERR )
180   CONTINUE
      IF ( AZERR.LE.ACCUR )  KCONV = KCONV + 1
      IF ( KCONV.GE.2 )      GOTO 210
C
200   CONTINUE
C
C ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============
C
C
C                                                 ** PRINT INTENSITIES
C
 210  IF ( PRNT(5) .AND. .NOT.ONLYFL )
     $     CALL  PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI,
     $                   MAXULV, MAXUMU)
C
      RETURN
C
c1010  FORMAT ( ////, 1X, 120('*'), /, 25X,
c     $  'DISCRETE ORDINATES RADIATIVE TRANSFER PROGRAM,',/, 20X,
c     $  'SPECIAL VERSION WITH GENERAL ANISOTROPIC INTERNAL SOURCE',
c     $  /, 1X, A, /, 1X, 120('*') )
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  ASYMTX( A, EVEC, EVAL, M, IA, IEVEC, IER,
     $                    AAD, EVECD, EVALD, WKD )
C
C    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======
C
C       SOLVES EIGENFUNCTION PROBLEM FOR REAL ASYMMETRIC MATRIX
C       FOR WHICH IT IS KNOWN A PRIORI THAT THE EIGENVALUES ARE REAL.
C
C       THIS IS AN ADAPTATION OF A SUBROUTINE 'EIGRF' IN THE IMSL
C       LIBRARY TO USE REAL INSTEAD OF COMPLEX ARITHMETIC, ACCOUNTING
C       FOR THE KNOWN FACT THAT THE EIGENVALUES AND EIGENVECTORS IN
C       THE DISCRETE ORDINATE SOLUTION ARE REAL.  OTHER CHANGES INCLUDE
C       PUTTING ALL THE CALLED SUBROUTINES IN-LINE, IN DELETING
C       THE PERFORMANCE INDEX CALCULATION, IN UPDATING MANY DO-LOOPS
C       TO FORTRAN77, AND IN CALCULATING THE MACHINE PRECISION
C       'TOL' INSTEAD OF SPECIFYING IT IN A DATA STATEMENT.
C
C       'EIGRF' IS BASED PRIMARILY ON 'EISPACK' ROUTINES.  THE MATRIX
C       IS FIRST BALANCED USING THE PARLETT-REINSCH ALGORITHM.  THEN
C       THE MARTIN-WILKINSON ALGORITHM IS APPLIED.
C
C       REFERENCES:
C          DONGARRA, J. AND C. MOLER, EISPACK -- A PACKAGE FOR SOLVING
C             MATRIX EIGENVALUE PROBLEMS, IN COWELL, ED., 1984:
C             SOURCES AND DEVELOPMENT OF MATHEMATICAL SOFTWARE,
C             PRENTICE-HALL, ENGLEWOOD CLIFFS, NJ
C         PARLETT AND REINSCH, 1969: BALANCING A MATRIX FOR CALCULATION
C             OF EIGENVALUES AND EIGENVECTORS, NUM. MATH. 13, 293-304
C         WILKINSON, J., 1965: THE ALGEBRAIC EIGENVALUE PROBLEM,
C             CLARENDON PRESS, OXFORD
C
C   I N P U T    V A R I A B L E S:
C
C        A    :  INPUT ASYMMETRIC MATRIX, DESTROYED AFTER SOLVED
C        M    :  ORDER OF -A-
C       IA    :  FIRST DIMENSION OF -A-
C    IEVEC    :  FIRST DIMENSION OF -EVEC-
C
C   O U T P U T    V A R I A B L E S:
C
C       EVEC  :  (UNNORMALIZED) EIGENVECTORS OF -A-
C                   ( COLUMN J CORRESPONDS TO EVAL(J) )
C
C       EVAL  :  (UNORDERED) EIGENVALUES OF -A- ( DIMENSION AT LEAST M )
C
C       IER   :  IF .NE. 0, SIGNALS THAT EVAL(IER) FAILED TO CONVERGE;
C                   IN THAT CASE EIGENVALUES IER+1,IER+2,...,M  ARE
C                   CORRECT BUT EIGENVALUES 1,...,IER ARE SET TO ZERO.
C
C   S C R A T C H   V A R I A B L E S:
C
C       WK    :  WORK AREA ( DIMENSION AT LEAST 2*M )
C       AAD    :  DOUBLE PRECISION STAND-IN FOR -A-
C       EVECD :  DOUBLE PRECISION STAND-IN FOR -EVEC-
C       EVALD :  DOUBLE PRECISION STAND-IN FOR -EVAL-
C       WKD   :  DOUBLE PRECISION STAND-IN FOR -WK-
C+---------------------------------------------------------------------+
C
      implicit none
      integer ia,ievec
      REAL*8  A( IA,* ),   EVAL(*), EVEC( IEVEC,* )
      REAL*8  AAD( IA,* ), WKD(*), EVALD(*), EVECD( IA,* )
      REAL*8 Z,ZERO,W,UU,TOL,VV,X,T,Y,SGN,S,ROW,Q,P,R,SCALE,RNORM,ONE
      REAL*8 G,H,C3,C6,C4,C5,REPL,C1,C2,discri,COL,F
      INTEGER L,LLL,N2,KA,J,II,I,LB,N,K,IN,KKK,IER,n1,M
*      DOUBLE PRECISION  D1MACH
      LOGICAL           NOCONV, NOTLAS
      DATA     C1 / 0.4375D0 /, C2/ 0.5D0 /, C3/ 0.75D0 /, C4/ 0.95D0 /,
     $         C5/ 16.D0 /, C6/ 256.D0 /, ZERO / 0.D0 /, ONE / 1.D0 /
C
C       gg initialization
      r=0
      p=0
      q=0
      TOL = 1d-30               ! D1MACH(3)
      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M )
     $     CALL ERRMSG( 'ASYMTX--BAD INPUT VARIABLE(S)', 
     .          .TRUE. )
C
C                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES
      IF ( M.EQ.1 )  THEN
         EVAL(1) = A(1,1)
         EVEC(1,1) = 1.0
         RETURN
C
      ELSE IF ( M.EQ.2 )  THEN
         DISCRI = ( A(1,1) - A(2,2) )**2 + 4. * A(1,2) * A(2,1)
         IF ( DISCRI.LT.0.0 )
     $        CALL ERRMSG( 'ASYMTX--COMPLEX EVALS IN 2X2 CASE',
     .          .TRUE. )
         SGN = 1.0
         IF ( A(1,1).LT.A(2,2) )  SGN = - 1.0
         EVAL(1) = 0.5 * ( A(1,1) + A(2,2) + SGN*SQRT(DISCRI) )
         EVAL(2) = 0.5 * ( A(1,1) + A(2,2) - SGN*SQRT(DISCRI) )
         EVEC(1,1) = 1.0
         EVEC(2,2) = 1.0
         IF ( A(1,1).EQ.A(2,2) .AND. (A(2,1).EQ.0.0.OR.A(1,2).EQ.0.0) )
     $        THEN
            RNORM = ABS(A(1,1))+ ABS(A(1,2)) + ABS(A(2,1)) + ABS(A(2,2))
            W = TOL * RNORM
            EVEC(2,1) = A(2,1) / W
            EVEC(1,2) = - A(1,2) / W
         ELSE
            EVEC(2,1) = A(2,1) / ( EVAL(1) - A(2,2) )
            EVEC(1,2) = A(1,2) / ( EVAL(2) - A(1,1) )
         ENDIF
         RETURN
      END IF
C                               ** PUT S.P. MATRIX INTO D.P. MATRIX
      DO 1  J = 1, M
         DO 1  K = 1, M
            AAD( J,K ) = ( A(J,K) )
    1 CONTINUE
C                                        ** INITIALIZE OUTPUT VARIABLES
      IER = 0
      DO 20 I = 1, M
         EVALD(I) = ZERO
         DO 10 J = 1, M
            EVECD(I,J) = ZERO
10       CONTINUE
         EVECD(I,I) = ONE
20    CONTINUE
C                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
C                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
C                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                  ** AND PUSH THEM DOWN
      RNORM = ZERO
      L  = 1
      K  = M
C
30    KKK = K
         DO 70  J = KKK, 1, -1
            ROW = ZERO
            DO 40 I = 1, K
               IF ( I.NE.J ) ROW = ROW + ABS( AAD(J,I) )
40          CONTINUE
            IF ( ROW.EQ.ZERO ) THEN
               WKD(K) = J
               IF ( J.NE.K ) THEN
                  DO 50 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,K)
                     AAD(I,K) = REPL
50                CONTINUE
                  DO 60 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(K,I)
                     AAD(K,I) = REPL
60                CONTINUE
               END IF
               K = K - 1
               GO TO 30
            END IF
70       CONTINUE
C                                     ** SEARCH FOR COLUMNS ISOLATING AN
C                                       ** EIGENVALUE AND PUSH THEM LEFT
80    LLL = L
         DO 120 J = LLL, K
            COL = ZERO
            DO 90 I = L, K
               IF ( I.NE.J ) COL = COL + ABS( AAD(I,J) )
90          CONTINUE
            IF ( COL.EQ.ZERO ) THEN
               WKD(L) = J
               IF ( J.NE.L ) THEN
                  DO 100 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,L)
                     AAD(I,L) = REPL
100               CONTINUE
                  DO 110 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(L,I)
                     AAD(L,I) = REPL
110               CONTINUE
               END IF
               L = L + 1
               GO TO 80
            END IF
120      CONTINUE
C                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K
      DO 130 I = L, K
         WKD(I) = ONE
130   CONTINUE
C
140   NOCONV = .FALSE.
         DO 200 I = L, K
            COL = ZERO
            ROW = ZERO
            DO 150 J = L, K
               IF ( J.NE.I ) THEN
                  COL = COL + ABS( AAD(J,I) )
                  ROW = ROW + ABS( AAD(I,J) )
               END IF
150         CONTINUE
            F = ONE
            G = ROW / C5
            H = COL + ROW
160         IF ( COL.LT.G ) THEN
               F   = F * C5
               COL = COL * C6
               GO TO 160
            END IF
            G = ROW * C5
170         IF ( COL.GE.G ) THEN
               F   = F / C5
               COL = COL / C6
               GO TO 170
            END IF
C                                                         ** NOW BALANCE
            IF ( (COL+ROW)/F .LT. C4*H ) THEN
               WKD(I)  = WKD(I) * F
               NOCONV = .TRUE.
               DO 180 J = L, M
                  AAD(I,J) = AAD(I,J) / F
180            CONTINUE
               DO 190 J = 1, K
                  AAD(J,I) = AAD(J,I) * F
190            CONTINUE
            END IF
200      CONTINUE
C
      IF ( NOCONV ) GO TO 140
C                                  ** IS -A- ALREADY IN HESSENBERG FORM?
      IF ( K-1 .LT. L+1 ) GO TO 350
C                                   ** TRANSFER -A- TO A HESSENBERG FORM
      DO 290 N = L+1, K-1
         H        = ZERO
         WKD(N+M) = ZERO
         SCALE    = ZERO
C                                                        ** SCALE COLUMN
         DO 210 I = N, K
            SCALE = SCALE + ABS(AAD(I,N-1))
210      CONTINUE
         IF ( SCALE.NE.ZERO ) THEN
            DO 220 I = K, N, -1
               WKD(I+M) = AAD(I,N-1) / SCALE
               H = H + WKD(I+M)**2
220         CONTINUE
            G = - SIGN( SQRT(H), WKD(N+M) )
            H = H - WKD(N+M) * G
            WKD(N+M) = WKD(N+M) - G
C                                                 ** FORM (I-(U*UT)/H)*A
            DO 250 J = N, M
               F = ZERO
               DO 230  I = K, N, -1
                  F = F + WKD(I+M) * AAD(I,J)
230            CONTINUE
               DO 240 I = N, K
                  AAD(I,J) = AAD(I,J) - WKD(I+M) * F / H
240            CONTINUE
250         CONTINUE
C                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 280 I = 1, K
               F = ZERO
               DO 260  J = K, N, -1
                  F = F + WKD(J+M) * AAD(I,J)
260            CONTINUE
               DO 270 J = N, K
                  AAD(I,J) = AAD(I,J) - WKD(J+M) * F / H
270            CONTINUE
280         CONTINUE
            WKD(N+M)  = SCALE * WKD(N+M)
            AAD(N,N-1) = SCALE * G
         END IF
290   CONTINUE
C
      DO 340  N = K-2, L, -1
         N1 = N + 1
         N2 = N + 2
         F  = AAD(N1,N)
         IF ( F.NE.ZERO ) THEN
            F  = F * WKD(N1+M)
            DO 300 I = N2, K
               WKD(I+M) = AAD(I,N)
300         CONTINUE
            IF ( N1.LE.K ) THEN
               DO 330 J = 1, M
                  G = ZERO
                  DO 310 I = N1, K
                     G = G + WKD(I+M) * EVECD(I,J)
310               CONTINUE
                  G = G / F
                  DO 320 I = N1, K
                     EVECD(I,J) = EVECD(I,J) + G * WKD(I+M)
320               CONTINUE
330            CONTINUE
            END IF
         END IF
340   CONTINUE
C
350   CONTINUE
      N = 1
      DO 370 I = 1, M
         DO 360 J = N, M
            RNORM = RNORM + ABS(AAD(I,J))
360      CONTINUE
         N = I
         IF ( I.LT.L .OR. I.GT.K ) EVALD(I) = AAD(I,I)
370   CONTINUE
      N = K
      T = ZERO
C                                         ** SEARCH FOR NEXT EIGENVALUES
380   IF ( N.LT.L ) GO TO 530
      IN = 0
      N1 = N - 1
      N2 = N - 2
C                          ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
390   CONTINUE
      DO 400 I = L, N
         LB = N+L - I
         IF ( LB.EQ.L ) GO TO 410
         S = ABS( AAD(LB-1,LB-1) ) + ABS( AAD(LB,LB) )
         IF ( S.EQ.ZERO ) S = RNORM
         IF ( ABS(AAD(LB,LB-1)) .LE. TOL*S ) GO TO 410
400   CONTINUE
C
410   X = AAD(N,N)
      IF ( LB.EQ.N ) THEN
C                                        ** ONE EIGENVALUE FOUND
         AAD(N,N)  = X + T
         EVALD(N) = AAD(N,N)
         N = N1
         GO TO 380
      END IF
C
      Y = AAD(N1,N1)
      W = AAD(N,N1) * AAD(N1,N)
      IF ( LB.EQ.N1 ) THEN
C                                        ** TWO EIGENVALUES FOUND
         P = (Y-X) * C2
         Q = P**2 + W
         Z = SQRT( ABS(Q) )
         AAD(N,N) = X + T
         X = AAD(N,N)
         AAD(N1,N1) = Y + T
C                                        ** REAL PAIR
         Z = P + SIGN(Z,P)
         EVALD(N1) = X + Z
         EVALD(N)  = EVALD(N1)
         IF ( Z.NE.ZERO ) EVALD(N) = X - W / Z
         X = AAD(N,N1)
C                                  ** EMPLOY SCALE FACTOR IN CASE
C                                  ** X AND Z ARE VERY SMALL
         R = SQRT( X*X + Z*Z )
         P = X / R
         Q = Z / R
C                                             ** ROW MODIFICATION
         DO 420 J = N1, M
            Z = AAD(N1,J)
            AAD(N1,J) = Q * Z + P * AAD(N,J)
            AAD(N,J)  = Q * AAD(N,J) - P * Z
420      CONTINUE
C                                             ** COLUMN MODIFICATION
         DO 430 I = 1, N
            Z = AAD(I,N1)
            AAD(I,N1) = Q * Z + P * AAD(I,N)
            AAD(I,N)  = Q * AAD(I,N) - P * Z
430      CONTINUE
C                                          ** ACCUMULATE TRANSFORMATIONS
         DO 440 I = L, K
            Z = EVECD(I,N1)
            EVECD(I,N1) = Q * Z + P * EVECD(I,N)
            EVECD(I,N)  = Q * EVECD(I,N) - P * Z
440      CONTINUE
C
         N = N2
         GO TO 380
      END IF
C
      IF ( IN.EQ.30 ) THEN
C                    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR
C                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE
         IER = N
         GO TO 670
      END IF
C                                                          ** FORM SHIFT
      IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
         T = T + X
         DO 450 I = L, N
            AAD(I,I) = AAD(I,I) - X
450      CONTINUE
         S = ABS(AAD(N,N1)) + ABS(AAD(N1,N2))
         X = C3 * S
         Y = X
         W = - C1 * S**2
      END IF
C
      IN = IN + 1
C                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
C
      DO 460 J = LB, N2
         I = N2+LB - J
         Z = AAD(I,I)
         R = X - Z
         S = Y - Z
         P = ( R * S - W ) / AAD(I+1,I) + AAD(I,I+1)
         Q = AAD(I+1,I+1) - Z - R - S
         R = AAD(I+2,I+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF ( I.EQ.LB ) GO TO 470
         UU = ABS( AAD(I,I-1) ) * ( ABS(Q) + ABS(R) )
         VV = ABS(P)*( ABS(AAD(I-1,I-1)) + ABS(Z) + ABS(AAD(I+1,I+1)) )
         IF ( UU .LE. TOL*VV ) GO TO 470
460   CONTINUE
C
470   CONTINUE
      AAD(I+2,I) = ZERO
      DO 480 J = I+3, N
         AAD(J,J-2) = ZERO
         AAD(J,J-3) = ZERO
480   CONTINUE
C
C             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N
C
      DO 520 KA = I, N1
         NOTLAS = KA.NE.N1
         IF ( KA.EQ.I ) THEN
            S = SIGN( SQRT( P*P + Q*Q + R*R ), P )
            IF ( LB.NE.I ) AAD(KA,KA-1) = - AAD(KA,KA-1)
         ELSE
            P = AAD(KA,KA-1)
            Q = AAD(KA+1,KA-1)
            R = ZERO
            IF ( NOTLAS ) R = AAD(KA+2,KA-1)
            X = ABS(P) + ABS(Q) + ABS(R)
            IF ( X.EQ.ZERO ) GO TO 520
            P = P / X
            Q = Q / X
            R = R / X
            S = SIGN( SQRT( P*P + Q*Q + R*R ), P )
            AAD(KA,KA-1) = - S * X
         END IF
         P = P + S
         X = P / S
         Y = Q / S
         Z = R / S
         Q = Q / P
         R = R / P
C                                                    ** ROW MODIFICATION
         DO 490 J = KA, M
            P = AAD(KA,J) + Q * AAD(KA+1,J)
            IF ( NOTLAS ) THEN
               P = P + R * AAD(KA+2,J)
               AAD(KA+2,J) = AAD(KA+2,J) - P * Z
            END IF
            AAD(KA+1,J) = AAD(KA+1,J) - P * Y
            AAD(KA,J)   = AAD(KA,J)   - P * X
490      CONTINUE
C                                                 ** COLUMN MODIFICATION
         DO 500 II = 1, MIN0(N,KA+3)
            P = X * AAD(II,KA) + Y * AAD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * AAD(II,KA+2)
               AAD(II,KA+2) = AAD(II,KA+2) - P * R
            END IF
            AAD(II,KA+1) = AAD(II,KA+1) - P * Q
            AAD(II,KA)   = AAD(II,KA) - P
500      CONTINUE
C                                          ** ACCUMULATE TRANSFORMATIONS
         DO 510 II = L, K
            P = X * EVECD(II,KA) + Y * EVECD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * EVECD(II,KA+2)
               EVECD(II,KA+2) = EVECD(II,KA+2) - P * R
            END IF
            EVECD(II,KA+1) = EVECD(II,KA+1) - P * Q
            EVECD(II,KA)   = EVECD(II,KA) - P
510      CONTINUE
C
520   CONTINUE
      GO TO 390
C                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR
530   CONTINUE
      IF ( RNORM.NE.ZERO ) THEN
         DO 560  N = M, 1, -1
            N2 = N
            AAD(N,N) = ONE
            DO 550  I = N-1, 1, -1
               W = AAD(I,I) - EVALD(N)
               IF ( W.EQ.ZERO ) W = TOL * RNORM
               R = AAD(I,N)
               DO 540 J = N2, N-1
                  R = R + AAD(I,J) * AAD(J,N)
540            CONTINUE
               AAD(I,N) = - R / W
               N2 = I
550         CONTINUE
560      CONTINUE
C                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS
C
         DO 580 I = 1, M
            IF ( I.LT.L .OR. I.GT.K ) THEN
               DO 570 J = I, M
                  EVECD(I,J) = AAD(I,J)
570            CONTINUE
            END IF
580      CONTINUE
C                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF ( K.NE.0 ) THEN
            DO 600  J = M, L, -1
               DO 600 I = L, K
                  Z = ZERO
                  DO 590 N = L, MIN0(J,K)
                     Z = Z + EVECD(I,N) * AAD(N,J)
590               CONTINUE
                  EVECD(I,J) = Z
600         CONTINUE
         END IF
C
      END IF
C
      DO 620 I = L, K
         DO 620 J = 1, M
            EVECD(I,J) = EVECD(I,J) * WKD(I)
620   CONTINUE
C                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO 640  I = L-1, 1, -1
         J = WKD(I)
         IF ( I.NE.J ) THEN
            DO 630 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
630         CONTINUE
         END IF
640   CONTINUE
C
      DO 660 I = K+1, M
         J = WKD(I)
         IF ( I.NE.J ) THEN
            DO 650 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
650         CONTINUE
         END IF
660   CONTINUE
C                         ** PUT RESULTS INTO OUTPUT ARRAYS
  670 CONTINUE
      DO 680 J = 1, M
         EVAL( J ) = EVALD(J)
         DO 680 K = 1, M
            EVEC( J,K ) = EVECD(J,K)
680   CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  CHEKIN( NLYR, DTAUC, SSALB, PMOM, SRC,
     $                    USRTAU, NTAU, UTAU, NSTR, USRANG,
     $                    NUMU, UMU, NPHI, PHI, FBEAM, UMU0,
     $                    PHI0, FISOT, LAMBER, ALBEDO, HL, BSRC,
     $                    TSRC, INTSRC, ONLYFL, ACCUR, MAXCLY,
     $                    MAXULV, MAXUMU, MAXCMU, MAXPHI, MXCLY,
     $                    MXULV,  MXUMU,  MXCMU,  MXPHI, TAUC, 
     $                    ipk )
C
C           CHECKS THE INPUT DIMENSIONS AND VARIABLES
C
      implicit none
      integer ipk
      LOGICAL  WRTBAD, WRTDIM
      LOGICAL  LAMBER, INTSRC, ONLYFL, USRANG, USRTAU, INPERR
      INTEGER  MAXCLY, MAXUMU, MAXULV, MAXCMU, MAXPHI, NLYR,
     $         NUMU, NSTR, NPHI, NTAU, MXCMU, MXUMU, MXPHI, MXCLY,
     $         MXULV
      REAL*8     ACCUR, ALBEDO, BSRC, DTAUC( MAXCLY ), FBEAM, 
     $          fisot(-ipk:-1),
     $         HL( 0:MAXCMU ), PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ),
     $         PHI0, SRC(3*MAXCLY,*), SSALB( MAXCLY ), TSRC,
     $         UMU( MAXUMU ), UMU0, UTAU( MAXULV ),
     $         TAUC( 0:* )
      REAL*8 RMU
      INTEGER BOUCLE,LU,IQ,LS,IU,J,K,LC
      REAL*8 fmin,flxalb,dref
C
C
      INPERR = .FALSE.
      IF ( NLYR.LT.1 ) INPERR = WRTBAD( 'NLYR' )
      IF ( NLYR.GT.MAXCLY ) INPERR = WRTBAD( 'MAXCLY' )
C
      DO 10  LC = 1, NLYR
         IF ( DTAUC(LC).LT.0.0 ) INPERR = WRTBAD( 'DTAUC' )
         IF ( SSALB(LC).LT.0.0 .OR. SSALB(LC).GT.1.0 )
     $        INPERR = WRTBAD( 'SSALB' )
         IF ( INTSRC )  THEN
            LS = 3*(LC-1) + 1
            DO 11 IQ = 1, NSTR 
              IF( SRC(LS,IQ).LT.0.0 ) THEN 
                      INPERR = WRTBAD( 'SRC' )
              
                        write(*,*)"SRC,LS,IQ"
                        write(*,*)SRC(LS,IQ),LS,IQ
      
                ENDIF
              IF( SRC(LS+1,IQ).LT.0.0 ) THEN
                      INPERR = WRTBAD( 'SRC')

                        write(*,*)"SRC,LS+1,IQ"
                        write(*,*)SRC(LS+1,IQ),LS+1,IQ

              ENDIF
             IF( SRC(LS+2,IQ).LT.0.0 ) THEN
                     INPERR = WRTBAD( 'SRC')

                        write(*,*)"SRC,LS+2,IQ"
                        write(*,*)SRC(LS+2,IQ),LS+2,IQ


             ENDIF
11          CONTINUE
            DO 12 IU = 1, NUMU
              IF( SRC(LS,IU).LT.0.0 ) INPERR = WRTBAD( 'SRCU')
              IF( SRC(LS+1,IU).LT.0.0 ) INPERR = WRTBAD( 'SRCU' )
              IF( SRC(LS+2,IU).LT.0.0 ) INPERR = WRTBAD( 'SRCU' )
12          CONTINUE
         ENDIF
         DO 5  K = 0, NSTR
            IF( PMOM(K,LC).LT.-1.0 .OR. PMOM(K,LC).GT.1.0 ) then
               INPERR = WRTBAD( 'PMOM' )
               print*,'k,lc,pmom:',k,lc,pmom(k,lc)
            end if
 5       CONTINUE
10    CONTINUE
C
      IF ( USRTAU )  THEN
c              print*, "NTAU : ",NTAU,"NSTR",NSTR
         IF ( NTAU.LT.1 ) INPERR = WRTBAD( 'NTAU' )
         IF ( MAXULV.LT.NTAU ) INPERR = WRTBAD( 'MAXULV')
         DO 20  LU = 1, NTAU
            IF( ABS(UTAU(LU)-TAUC(NLYR)).LE.1.E-4) UTAU(LU) = TAUC(NLYR)
            
            IF( ABS((UTAU(LU)-TAUC(NLYR))/TAUC(NLYR)).LE.1.E-4)
     *  UTAU(LU) = TAUC(NLYR)
            IF( UTAU(LU).LT.0.0 .OR. UTAU(LU).GT.TAUC(NLYR) ) then
            print*,"Probleme utau"
            if (UTAU(LU).GT.TAUC(NLYR)) print*,"Merde sur utau lu"
            print*,tauc(nlyr),nlyr
            if (UTAU(LU).LT.0.) print*,"Merde sur utau<0."
            print*,"valeur utau lu",utau(LU),LU

               INPERR = WRTBAD( 'UTAU' )

       endif
20       CONTINUE
      ELSE
         IF ( MAXULV.LT.NLYR+1 ) INPERR = WRTBAD( 'MAXULV' )
      END IF
C
      IF ( NSTR.LT.2 .OR. MOD(NSTR,2).NE.0 )
     $                  INPERR = WRTBAD( 'NSTR' )
      IF ( NSTR.GT.MAXCMU ) INPERR = WRTBAD( 'MAXCMU' )
C
      IF ( USRANG )  THEN
         IF ( NUMU.LT.0 ) INPERR = WRTBAD( 'NUMU' )
         IF ( .NOT.ONLYFL .AND. NUMU.EQ.0 )
     $                  INPERR = WRTBAD( 'NUMU'  )
         IF ( NUMU.GT.MAXUMU ) INPERR = WRTBAD( 'MAXUMU' )
         DO 30  IU = 1, NUMU
            IF( UMU(IU).LT.-1.0 .OR. UMU(IU).GT.1.0 .OR. UMU(IU).EQ.0.0)
     $           INPERR = WRTBAD( 'UMU' )
            IF( IU.GT.1 .AND. UMU(IU).LT.UMU(IU-1) )
     $           INPERR = WRTBAD( 'UMU' )
30       CONTINUE
      ELSE
         IF( MAXUMU.LT.NSTR ) INPERR = WRTBAD( 'MAXUMU' )
      END IF
C
      IF ( .NOT.ONLYFL )  THEN
         IF ( NPHI.LE.0 ) INPERR = WRTBAD( 'NPHI' )
         IF ( NPHI.GT.MAXPHI ) INPERR = WRTBAD( 'MAXPHI' )
         DO 40  J = 1, NPHI
            IF ( PHI(J).LT.0.0 .OR. PHI(J).GT.360.0 )
     $           INPERR = WRTBAD( 'PHI' )
40       CONTINUE
      END IF
C
      IF ( FBEAM.LT.0.0 ) INPERR = WRTBAD( 'FBEAM' )
      IF ( FBEAM.GT.0.0 .AND. ( UMU0.LE.0.0 .OR. UMU0.GT.1.0 ) )
     $     INPERR = WRTBAD( 'UMU0' )
      IF ( FBEAM.GT.0.0 .AND. ( PHI0.LT.0.0 .OR. PHI0.GT.360.0 ) )
     $     INPERR = WRTBAD( 'PHI0')
      IF ( fmin(fisot,ipk).LT.0.0 ) INPERR = WRTBAD( 'FISOT' )
      IF ( LAMBER )  THEN
         IF ( ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0 )
     $        INPERR = WRTBAD( 'ALBEDO' )
      ELSE
C                    ** MAKE SURE FLUX ALBEDO AT DENSE MESH OF INCIDENT
C                    ** ANGLES DOES NOT ASSUME UNPHYSICAL VALUES
C
C         DO 50  RMU = 0.0, 1.0, 0.01
C        Modification pour g95
        RMU=0.00
         DO 50  BOUCLE = 0, 100, 1
          
            FLXALB = DREF( RMU, HL, NSTR)
            IF ( FLXALB.LT.0.0 .OR. FLXALB.GT.1.0 )
     $           INPERR = WRTBAD( 'HL' )
        RMU=RMU+0.01
50      CONTINUE
      ENDIF
      IF ( INTSRC )  THEN
         IF ( BSRC.LT.0.0 ) INPERR = WRTBAD( 'BSRC' )
         IF ( TSRC.LT.0.0 ) INPERR = WRTBAD( 'TSRC' )
      END IF
C
      IF ( ACCUR.LT.0.0 .OR. ACCUR.GT.1.E-2 )
     $     INPERR = WRTBAD( 'ACCUR' )
C
      IF ( MXCLY.LT.NLYR ) INPERR = WRTDIM( 'MXCLY', NLYR)
      IF ( USRTAU .AND. MXULV.LT.NTAU )
     $     INPERR = WRTDIM( 'MXULV', NTAU )
      IF ( .NOT.USRTAU .AND. MXULV.LT.NLYR+1 )
     $     INPERR = WRTDIM( 'MXULV', NLYR+1 )
      IF ( MXCMU.LT.NSTR ) INPERR = WRTDIM( 'MXCMU', NSTR)
      IF ( USRANG .AND. MXUMU.LT.NUMU )
     $     INPERR = WRTDIM( 'MXUMU', NUMU )
      IF ( .NOT.USRANG .AND. MXUMU.LT.NSTR )
     $      INPERR = WRTDIM( 'MXUMU', NSTR )
      IF ( .NOT.ONLYFL .AND. MXPHI.LT.NPHI )
     $      INPERR = WRTDIM( 'MXPHI', NPHI )
C
      IF ( INPERR )
     $   CALL ERRMSG( 'DISORT--INPUT AND/OR DIMENSION ERRORS', .TRUE.
     $                   )
C
      RETURN
      END
      real*8 function fmin(array,n)
      implicit none
      integer n,i
*       return minimum value of an array
      real*8 array(n)
      fmin=array(1)
      do 100 i=2,n,1
      fmin=min(fmin,array(i))
100   continue
      return
      end
*
*----------------------------------------------------------------------
*
      SUBROUTINE  CMPINT( FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZ,
     $                    MXCMU, MXULV, MXUMU, NCUT, NN, NSTR,
     $                    INTSRC, NTAU, TAUCPR, UMU0, UTAUPR,
     $                    ZZ, ZPLK0, ZPLK1, ZPLKA, UUM )
C
C       CALCULATES THE FOURIER INTENSITY COMPONENTS AT THE QUADRATURE
C       ANGLES FOR AZIMUTHAL EXPANSION TERMS (MAZ) IN EQ. SD(2)
C
C    I N P U T    V A R I A B L E S:
C
C       KK      :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       GC      :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       LL      :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
C                  BY SOLVING SCALED VERSION OF EQ. SC(5);
C                  EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
C       LYRCUT  :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       MAZ     :  ORDER OF AZIMUTHAL COMPONENT
C       NCUT    :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION
C                  OPTICAL DEPTH EXCEEDS -ABSCUT-
C       NN      :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
C       TAUCPR  :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       UTAUPR  :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
C                  COORDINATES;  EQUAL TO -UTAU- IF NO DELTA-M
C       ZZ      :  BEAM SOURCE VECTORS IN EQ. SS(19)
C       ZPLK0   :  INTERNAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
C       ZPLK1   :  INTERNAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)
C       ZPLKA   :  INTERNAL SOURCE VECTORS -ZA-, ALFA IN EQ. 
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C    O U T P U T   V A R I A B L E S:
C
C       UUM     :  FOURIER COMPONENTS OF THE INTENSITY IN EQ.  SD(12)
C                   ( AT POLAR QUADRATURE ANGLES )
C
C    I N T E R N A L   V A R I A B L E S:
C
C       FACT    :  EXxP( - UTAUPR / UMU0 )
C       ZINT    :  INTENSITY OF M=0 CASE, IN EQ. SC(1)
C+----------------------------------------------------------------------
C
      implicit none
      integer lu,iq,jq,ntau,nstr,nn,ncut,mxumu,mxulv,maz,mxcmu
      integer lyu
        real*8 fbeam
      real*8 zint,umu0,exxp
       LOGICAL  LYRCUT, INTSRC
       INTEGER  LAYRU(*)
       REAL*8     UUM( MXUMU, MXULV, 0:* )
       REAL*8     GC( MXCMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ),
     $          TAUCPR( 0:* ), UTAUPR(*), ZZ( MXCMU, *),
     $          ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* ), ZPLKA( * ) 
C
C
C                                                  ** ZERO OUTPUT ARRAY
       CALL ZEROIT( UUM, MXUMU*MXULV*(MXCMU + 1) )
C
C                                       ** LOOP OVER USER LEVELS
       DO 100  LU = 1, NTAU
C
          LYU = LAYRU(LU)
          IF ( LYRCUT .AND. LYU.GT.NCUT )  GO TO 100
C
          DO 20  IQ = 1, NSTR
             ZINT = 0.0
             DO 10  JQ = 1, NN
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                   EXxP( - KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU)) )
10           CONTINUE
             DO 11  JQ = NN+1, NSTR
                ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                 EXxP( - KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU-1)) )
11           CONTINUE
C
             UUM(IQ,LU,MAZ) = ZINT
             IF ( FBEAM.GT.0.0 )
     $            UUM(IQ,LU,MAZ) = ZINT + ZZ(IQ,LYU)
     $                                    * EXxP( - UTAUPR(LU) / UMU0 )
             IF ( INTSRC .AND. MAZ.EQ.0 )
     $            UUM(IQ,LU,MAZ) = UUM(IQ,LU,MAZ) + 
     $                             EXxP(-ZPLKA(LYU)*UTAUPR(LU))*
     $                         (ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU))
20        CONTINUE
C
100   CONTINUE
C
      RETURN
      END
      REAL*8 FUNCTION  DREF( MU, HL, NSTR )
C
C        EXACT FLUX ALBEDO FOR GIVEN ANGLE OF INCIDENCE, GIVEN
C        A BIDIRECTIONAL REFLECTIVITY CHARACTERIZED BY ITS
C        LEGENDRE COEFFICIENTS ( NOTE** THESE WILL ONLY AGREE
C        WITH BOTTOM-BOUNDARY ALBEDOS CALCULATED BY 'DISORT' IN
C        THE LIMIT AS NUMBER OF STREAMS GO TO INFINITY, BECAUSE
C        'DISORT' EVALUATES THE INTEGRAL 'CL' ONLY APPROXIMATELY,
C        BY QUADRATURE, WHILE THIS ROUTINE CALCULATES IT EXACTLY. )
C
C      INPUT :   MU     COSINE OF INCIDENCE ANGLE
C                HL     LEGENDRE COEFFICIENTS OF BIDIRECTIONAL REF'Y
C              NSTR     NUMBER OF ELEMENTS OF 'HL' TO CONSIDER
C
C      INTERNAL VARIABLES (P-SUB-L IS THE L-TH LEGENDRE POLYNOMIAL) :
C
C              CL    INTEGRAL FROM 0 TO 1 OF  MU * P-SUB-L(MU)
C                       (VANISHES FOR  L = 3, 5, 7, ... )
C              PL    P-SUB-L
C            PLM1    P-SUB-(L-1)
C            PLM2    P-SUB-(L-2)
C
              implicit none
                  integer maxtrm
      PARAMETER  ( MAXTRM = 100 )
      LOGICAL      PASS1
      REAL*8         MU, HL( 0:* )
      real*8 CL,PL,PLM1,PLM2
      integer L,NSTR
        REAL*8         C( MAXTRM )
      DATA  PASS1 / .TRUE. /
      save c
C
C
      IF ( PASS1 )  THEN
         PASS1 = .FALSE.
         CL = 0.125
         C(2) = 10. * CL
         DO 1  L = 4, MAXTRM, 2
            CL = - CL * (L-3) / (L+2)
            C(L) = 2. * (2*L+1) * CL
    1    CONTINUE
      END IF
C
      IF ( NSTR.GT.MAXTRM )  CALL
     $     ERRMSG( 'DREF--PARAMETER MAXTRM TOO SMALL', .TRUE. )
C
      DREF = HL(0) - 2.*HL(1) * MU
      PLM2 = 1.0
      PLM1 = - MU
      DO 10  L = 2, NSTR-1
C                                ** LEGENDRE POLYNOMIAL RECURRENCE
C
         PL = ( (2*L-1) * (-MU) * PLM1 - (L-1) * PLM2 ) / L
         IF( MOD(L,2).EQ.0 )  DREF = DREF + C(L) * HL(L) * PL
         PLM2 = PLM1
         PLM1 = PL
   10 CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  FLUXES( CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
     $                    MXCMU, MXULV, NCUT, NN, NSTR, NTAU, PI,
     $                    PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,
     $                    XR0, XR1, XRA, ZZ, ZPLK0, ZPLK1, ZPLKA, 
     $                    DFDT, FLUP, FLDN, FLDIR, RFLDIR, RFLDN, 
     $                    UAVG, U0C, MAXULV)
C
C       CALCULATES THE RADIATIVE FLUXES, MEAN INTENSITY, AND FLUX
C       DERIVATIVE WITH RESPECT TO OPTICAL DEPTH FROM THE M=0 INTENSITY
C       COMPONENTS (THE AZIMUTHALLY-AVERAGED INTENSITY)
C
C    I N P U T     V A R I A B L E S:
C
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       LAYRU    :  LAYER NUMBER OF USER LEVEL -UTAU-
C       LL       :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
C                   BY SOLVING SCALED VERSION OF EQ. SC(5);
C                   EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
C       NCUT     :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION
C                     OPTICAL DEPTH EXCEEDS -ABSCUT-
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       UTAUPR   :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
C                     COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M
C       XR0      :  EXPANSION OF INERNAL SOURCE FUNCTION IN EQ. SS(14)
C       XR1      :  EXPANSION OF INTERNAL SOURCE FUNCTION EQS. SS(14)
C       XRA      :  EXPANSION OF INTERNAL SOURCE FUNCTION EQS. SS(14)
C       ZZ       :  BEAM SOURCE VECTORS IN EQ. SS(19)
C       ZPLK0    :  INTERNAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
C       ZPLK1    :  INTERNAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)
C       ZPLKA    :  INTERNAL SOURCE VECTORS -ZA-, ALFA IN EQ.
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T     V A R I A B L E S:
C
C       U0C      :  AZIMUTHALLY AVERAGED INTENSITIES
C                   ( AT POLAR QUADRATURE ANGLES )
C       (RFLDIR, RFLDN, FLUP, DFDT, UAVG ARE 'DISORT' OUTPUT VARIABLES)
C
C   I N T E R N A L       V A R I A B L E S:
C
C       DIRINT   :  DIRECT INTENSITY ATTENUATED
C       FDNTOT   :  TOTAL DOWNWARD FLUX (DIRECT + DIFFUSE)
C       FLDIR    :  DIRECT-BEAM FLUX (DELTA-M SCALED)
C       FLDN     :  DIFFUSE DOWN-FLUX (DELTA-M SCALED)
C       FNET     :  NET FLUX (TOTAL-DOWN - DIFFUSE-UP)
C       FACT     :  EXxP( - UTAUPR / UMU0 )
C       SORC     :  INTERNAL SOURCE FUNCTION INTEGRATED OVER ALL
C                   ANGLES
C       ZINT     :  INTENSITY OF m = 0 CASE, IN EQ. SC(1)
C+---------------------------------------------------------------------+
C
      implicit none
      real*8 sorc,zint,fnet,fact,exxp,fdntot,dirint,ang1,ang2,pi,umu0
      integer lu,lyu,iq,jq,maxulv,mxcmu,ncut,mxulv,nn,ntau,nstr
      real*8 fbeam
      LOGICAL LYRCUT, PRNT(*)
      REAL*8  DFDT(*), FLUP(*), FLDIR(*), FLDN(*), RFLDIR(*), RFLDN(* ),
     $        U0C( MXCMU,MXULV ), UAVG(*)
      INTEGER LAYRU(*)
      REAL*8    CMU(*), CWT(*), GC( MXCMU,MXCMU,* ), KK( MXCMU,* ),
     $        LL( MXCMU,* ), SSALB(*), TAUCPR( 0:* ),
     $        UTAU(*), UTAUPR(*), XR0(MXCMU,*), 
     $        XR1(MXCMU,*), XRA(*), ZZ( MXCMU,* ),
     $        ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* ), ZPLKA( * ) 
C
C   GG : init fact at 0
      fact=0

      IF ( PRNT(2) )  WRITE( 6,1010 )
C                                          ** ZERO DISORT OUTPUT ARRAYS
      CALL  ZEROIT( U0C, MXULV*MXCMU )
      CALL  ZEROIT( RFLDIR, MAXULV )
      CALL  ZEROIT( FLDIR,  MXULV )
      CALL  ZEROIT( RFLDN,  MAXULV )
      CALL  ZEROIT( FLDN,   MXULV )
      CALL  ZEROIT( FLUP,   MAXULV )
      CALL  ZEROIT( UAVG,   MAXULV )
      CALL  ZEROIT( DFDT,   MAXULV )
C                                        ** LOOP OVER USER LEVELS
      DO 100  LU = 1, NTAU
C
         LYU = LAYRU(LU)
C
         IF ( LYRCUT .AND. LYU.GT.NCUT ) THEN
C                                                ** NO RADIATION REACHES
C                                                ** THIS LEVEL
            FDNTOT = 0.0
            FNET   = 0.0
            SORC = 0.0
            GO TO 90
         END IF
C
         IF ( FBEAM.GT.0.0 )  THEN
            FACT  = EXxP( - UTAUPR(LU) / UMU0 )
            DIRINT = FBEAM * FACT
            FLDIR(  LU ) = UMU0 * ( FBEAM * FACT )
            RFLDIR( LU ) = UMU0 * FBEAM * EXxP( - UTAU( LU ) / UMU0 )
         ELSE
            DIRINT = 0.0
            FLDIR(  LU ) = 0.0
            RFLDIR( LU ) = 0.0
         END IF
C
         DO 20  IQ = 1, NN
C
            ZINT = 0.0
            DO 10  JQ = 1, NN
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                EXxP(- KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) )
10          CONTINUE
            DO 11  JQ = NN+1, NSTR
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                EXxP(- KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)) )
11          CONTINUE
C
            U0C( IQ,LU ) = ZINT
            IF ( FBEAM.GT.0.0 )  U0C( IQ,LU ) = ZINT + ZZ(IQ,LYU) * FACT
            U0C( IQ,LU ) = U0C( IQ,LU ) +
     $                             EXxP(-ZPLKA(LYU)*UTAUPR(LU))*
     $                         (ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU))
*       if(lu.le.10.and.iq.eq.1)
*     . write(6,1000) lu,lyu,utaupr(lu),taucpr(lyu-1),taucpr(lyu),
*     .         iq,zplka(lyu),zplk0(iq,lyu),zplk1(iq,lyu),
*     .          EXxP(-ZPLKA(LYU)*UTAUPR(LU))*
*     .                    (ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU)),
*     .         zint,u0c(iq,lu)
*1000   format(1x,2i4,3(1pe11.2),i4,3e12.3,4x,3(1pe12.3))
            UAVG(LU) = UAVG(LU) + CWT(NN+1-IQ) * U0C( IQ,LU )
            FLDN(LU) = FLDN(LU) + CWT(NN+1-IQ)*CMU(NN+1-IQ) * U0C(IQ,LU)
20       CONTINUE
C
         DO 40  IQ = NN+1, NSTR
C
            ZINT = 0.0
            DO 30  JQ = 1, NN
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                EXxP(- KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) )
30          CONTINUE
            DO 31  JQ = NN+1, NSTR
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                EXxP(- KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)) )
31          CONTINUE
C
            U0C( IQ,LU ) = ZINT
            IF ( FBEAM.GT.0.0 )  U0C( IQ,LU ) = ZINT + ZZ(IQ,LYU) * FACT
            U0C( IQ,LU ) = U0C( IQ,LU ) + 
     $                             EXxP(-ZPLKA(LYU)*UTAUPR(LU))*
     $                         (ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU))
            UAVG(LU) = UAVG(LU) + CWT(IQ-NN) * U0C( IQ,LU )
            FLUP(LU) = FLUP(LU) + CWT(IQ-NN) * CMU(IQ-NN) * U0C( IQ,LU )
40       CONTINUE
C
         FLUP( LU )  = 2.0 * PI * FLUP( LU )
         FLDN( LU )  = 2.0 * PI * FLDN( LU )
         FDNTOT = FLDN( LU ) + FLDIR( LU )
         FNET   = FDNTOT - FLUP( LU )
         RFLDN( LU ) = FDNTOT - RFLDIR( LU )
         UAVG( LU ) = ( 2.0 * PI * UAVG(LU) + DIRINT ) / ( 4.*PI )
         SORC = 0.0
         DO 50 IQ = 1, NN
                SORC = SORC + CWT(IQ)*
     $                           EXxP(-XRA(LYU)*UTAUPR(LU))*
     $                          (XR0(IQ,LYU) + XR1(IQ,LYU)* UTAUPR(LU))
50       CONTINUE
         DO 51 IQ = NN+1, NSTR
                SORC = SORC + CWT(IQ-NN)*
     $                           EXxP(-XRA(LYU)*UTAUPR(LU))*
     $                          (XR0(IQ,LYU) + XR1(IQ,LYU)* UTAUPR(LU))
51       CONTINUE
C                             **  FACTOR OF 0.5 BECAUSE OF THE
C                             **  USE OF DOUBLE GAUSSIAN QUADRATURE
         SORC = 0.5 * SORC    
         DFDT( LU ) = ( 1.0-SSALB(LYU) ) * 4.*PI* UAVG(LU) 
     $                 - 4. * PI * SORC 
 90      IF( PRNT(2) )  WRITE( 6,1020 ) UTAU(LU), LYU, RFLDIR(LU),
     $                                 RFLDN(LU), FDNTOT, FLUP(LU),
     $                                 FNET, UAVG(LU), SORC, DFDT(LU)
100   CONTINUE
C
      IF ( PRNT(3) )  THEN
         WRITE ( 6,1100 )
         DO 200  LU = 1, NTAU
            WRITE( 6,1110 )  UTAU( LU )
            DO  200  IQ = 1, NN
               ANG1 = 180./PI * ACOS( CMU(2*NN-IQ+1) )
               ANG2 = 180./PI * ACOS( CMU(IQ) )
               WRITE( 6,1120 ) ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),
     $                         ANG2, CMU(IQ),        U0C(IQ+NN,LU)
200      CONTINUE
      END IF
C
1010  FORMAT( //, 21X,
     $ '<----------------------- FLUXES ----------------------->', /,
     $ '   OPTICAL  COMPU    DOWNWARD    DOWNWARD    DOWNWARD     ',
     $ ' UPWARD                    MEAN     INTERNAL  D(NET FLUX)', /,
     $ '     DEPTH  LAYER      DIRECT     DIFFUSE       TOTAL     ',
     $ 'DIFFUSE         NET   INTENSITY      SOURCE   / D(OP DEP)', / )
1020  FORMAT( F10.4, I7, 1P,7E12.3, E14.3 )
1100  FORMAT( //, ' ******** AZIMUTHALLY AVERAGED INTENSITIES',
     $      ' ( AT POLAR QUADRATURE ANGLES ) *******' )
1110  FORMAT( /, ' OPTICAL DEPTH =', F10.4, //,
     $  '     ANGLE (DEG)   COS(ANGLE)     INTENSITY',
     $  '     ANGLE (DEG)   COS(ANGLE)     INTENSITY' )
1120  FORMAT( 2( 0P,F16.4, F13.5, 1P,E14.3 ) )
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  LEPOLY( NMU, M, MAXMU, TWONM1, MU, YLM)
C
C       COMPUTES THE NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL,
C       DEFINED IN TERMS OF THE ASSOCIATED LEGENDRE POLYNOMIAL
C       PLM = P-SUB-L-SUPER-M AS
C
C             YLM(MU) = SQRT( (L-M)!/(L+M)! ) * PLM(MU)
C
C       FOR FIXED ORDER -M- AND ALL DEGREES FROM L = M TO TWONM1.
C       WHEN M.GT.0, ASSUMES THAT Y-SUB(M-1)-SUPER(M-1) IS AVAILABLE
C       FROM A PRIOR CALL TO THE ROUTINE.
C
C       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
C                  High-Order Associated Legendre Polynomials,
C                  J. Quant. Spectrosc. Radiat. Transfer 10,
C                  557-562, 1970.  (hereafter D/A)
C
C       METHOD: Varying degree recurrence relationship.
C
C       NOTE 1: The D/A formulas are transformed by
C               setting  M = n-1; L = k-1.
C       NOTE 2: Assumes that routine is called first with  M = 0,
C               then with  M = 1, etc. up to  M = TWONM1.
C       NOTE 3: Loops are written in such a way as to vectorize.
C
C  I N P U T     V A R I A B L E S:
C
C       NMU    :  NUMBER OF ARGUMENTS OF -YLM-
C       M      :  ORDER OF -YLM-
C       MAXMU  :  FIRST DIMENSION OF -YLM-
C       TWONM1 :  MAX DEGREE OF -YLM-
C       MU(I)  :  I = 1 TO NMU, ARGUMENTS OF -YLM-
C       IF M.GT.0, YLM(M-1,I) FOR I = 1 TO NMU IS REQUIRED
C
C  O U T P U T     V A R I A B L E:
C
C       YLM(L,I) :  L = M TO TWONM1, NORMALIZED ASSOCIATED LEGENDRE
C                   POLYNOMIALS EVALUATED AT ARGUMENT -MU(I)-
C+---------------------------------------------------------------------+
      implicit none
        integer maxmu,maxsqt,i,l,ns
        real*8 tmp1,tmp2
        REAL*8     MU(*), YLM( 0:MAXMU,* )
      INTEGER  M, NMU, TWONM1
      PARAMETER  ( MAXSQT = 1000 )
      REAL*8     SQT( MAXSQT )
      LOGICAL  PASS1
      SAVE  SQT, PASS1
      DATA  PASS1 / .TRUE. /
C
C
      IF ( PASS1 )  THEN
         PASS1 = .FALSE.
         DO 1  NS = 1, MAXSQT
            SQT( NS ) = SQRT( FLOAT(NS) )
    1    CONTINUE
      ENDIF
C
      IF ( 2*TWONM1 .GT. MAXSQT )
     $   CALL ERRMSG( 'LEPOLY--NEED TO INCREASE PARAM MAXSQT', .TRUE.
     $            )
C
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
C
      ELSE
C
         DO  30  I = 1, NMU
C                               ** Y-SUB-M-SUPER-M; DERIVED FROM
C                               ** D/A EQS. (11,12)
C
            YLM( M,I) = - SQT( 2*M-1 ) / SQT( 2*M )
     $                  * SQRT( 1. - MU(I)**2 ) * YLM( M-1,I )
C
C                              ** Y-SUB-(M+1)-SUPER-M; DERIVED FROM
C                              ** D/A EQS. (13,14) USING EQS. (11,12)
C
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
C
      END IF
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  PRAVIN( UMU, NUMU, MAXUMU, UTAU, NTAU, U0U)
C
C        PRINT AZIMUTHALLY AVERAGED INTENSITIES AT USER ANGLES
C
      implicit none
      integer iounit,np,iu,lu,npass,lenfmt,iumin,ntau,numu,maxumu
        integer iumax      
        REAL*8     UMU(*), UTAU(*), U0U( MAXUMU,* )
C
C
        iounit=6
      WRITE ( iounit, '(//,A)' )
     $         ' *********  AZIMUTHALLY AVERAGED INTENSITIES '
     $       // '(USER POLAR ANGLES)  *********'
      LENFMT = 8
      NPASS = 1 + NUMU / LENFMT
      IF ( MOD(NUMU,LENFMT) .EQ. 0 )  NPASS = NPASS - 1
      DO 10  NP = 1, NPASS
         IUMIN = 1 + LENFMT * (NP-1)
         IUMAX = MIN0( LENFMT*NP, NUMU )
         WRITE ( 6,101 )  ( UMU(IU), IU = IUMIN, IUMAX )
         DO 10  LU = 1, NTAU
            WRITE( 6,102 ) UTAU(LU), ( U0U(IU,LU), IU=IUMIN,IUMAX)
 10   CONTINUE
C
      RETURN
C
101   FORMAT( /, 3X,'OPTICAL   POLAR ANGLE COSINES',
     $        /, 3X,'  DEPTH', 8F14.5 )
102   FORMAT( 0P,F10.4, 1P,8E14.4 )
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  PRTINP( NLYR, DTAUC, SSALB, PMOM, NTAU, UTAU, 
     $                    NSTR, NUMU, UMU, NPHI, PHI, FBEAM, 
     $                    UMU0, PHI0, FISOT, LAMBER, ALBEDO, HL, 
     $                    BSRC, TSRC, DELTAM, INTSRC, ONLYFL, 
     $                    ACCUR, FLYR, LYRCUT, OPRIM, TAUC, TAUCPR,
     $                    MAXCMU, PRTMOM, ipk)
C
C        PRINT VALUES OF INPUT VARIABLES
C
      implicit none
      integer lc,lu,iu,j,ipk,MAXCMU,NLYR,NSTR,NTAU,NUMU
      integer NPHI,k
      real*8 fbeam
      real*8 yessct,bsrc,tsrc,accur,phi0,albedo,umu0
      LOGICAL  DELTAM, LAMBER, LYRCUT, INTSRC, ONLYFL, PRTMOM
      REAL*8     UMU(*), FLYR(*), DTAUC(*), OPRIM(*), PHI(*),
     $         PMOM( 0:MAXCMU,* ), SSALB(*), UTAU(*), TAUC( 0:* ),
     $         TAUCPR( 0:* ), HL( 0:MAXCMU ), fisot(-ipk:-1)
C
C
      WRITE( 6,1010 )  NSTR, NLYR
      WRITE( 6,1030 )  NTAU, (UTAU(LU), LU = 1, NTAU)
      IF ( .NOT.ONLYFL )
     $      WRITE( 6,1040 )  NUMU, ( UMU(IU), IU = 1, NUMU )
      IF ( .NOT.ONLYFL )
     $      WRITE( 6,1050 )  NPHI, ( PHI(J), J = 1, NPHI )
      IF ( .NOT. INTSRC )  WRITE( 6,1100 )
      WRITE( 6,1060 ) FBEAM, UMU0, PHI0, (FISOT(k),k=-nstr/2,-1,1)
      IF ( LAMBER )   WRITE( 6,1080 ) ALBEDO
      IF ( .NOT.LAMBER )  WRITE( 6,1090 ) ( HL(K), K = 0, NSTR )
      IF ( INTSRC )  WRITE( 6,1110 ) BSRC, TSRC
      IF ( DELTAM )      WRITE( 6,1120 )
      IF ( .NOT.DELTAM ) WRITE( 6,1130 )
      IF ( ONLYFL )  THEN
         WRITE( 6,1140 )
      ELSE
         WRITE( 6,1150 )
      ENDIF
      WRITE( 6,1160 )  ACCUR
      IF ( LYRCUT )  WRITE( 6,1170 )
      IF( INTSRC       )  WRITE ( 6,1190 )
      IF( .NOT. INTSRC )  WRITE ( 6,1191 )
      YESSCT = 0.0
      DO 10 LC = 1, NLYR
         YESSCT = YESSCT + SSALB(LC)
         IF( INTSRC )
     $       WRITE( 6,1200 )  LC, DTAUC(LC), TAUC(LC), SSALB(LC),
     $                    FLYR(LC), TAUCPR(LC)-TAUCPR(LC-1), TAUCPR(LC),
     $                    OPRIM(LC), PMOM(1,LC)
         IF( .NOT. INTSRC )
     $       WRITE( 6,1200 )  LC, DTAUC(LC), TAUC(LC), SSALB(LC),
     $                    FLYR(LC), TAUCPR(LC)-TAUCPR(LC-1), TAUCPR(LC),
     $                    OPRIM(LC), PMOM(1,LC)
 10   CONTINUE
C
      IF( PRTMOM .AND. YESSCT.GT.0.0 )  THEN
         WRITE( 6, '(/,A)' )  ' LAYER   PHASE FUNCTION MOMENTS'
         DO 20 LC = 1, NLYR
            IF( SSALB(LC).GT.0.0 )
     $          WRITE( 6,1300 )  LC, ( PMOM(K,LC), K = 0, NSTR )
 20      CONTINUE
      ENDIF
C
      RETURN
C
1010  FORMAT ( /, ' NO. STREAMS =', I4,
     $  '     NO. COMPUTATIONAL LAYERS =', I4 )
1030  FORMAT( I4,' USER OPTICAL DEPTHS :',10F10.4, /, (26X,10F10.4) )
1040  FORMAT( I4,' USER POLAR ANGLE COSINES :',10F9.5,/,(31X,10F9.5) )
1050  FORMAT( I4,' USER AZIMUTHAL ANGLES :', 10F9.2, /, (28X,10F9.2) )
1060  FORMAT( '    INCIDENT BEAM WITH INTENSITY =', 1P,E11.3, ' AND',
     $ ' POLAR ANGLE COSINE = ', 0P,F8.5,'  AND AZIMUTH ANGLE =', F7.2,
     $ /,'    PLUS INCIDENT INTENSITY =', 1P,(E11.3) )
1080  FORMAT( '    BOTTOM ALBEDO (LAMBERTIAN) =', 0P,F8.4 )
1090  FORMAT( '    LEGENDRE COEFFS OF BOTTOM BIDIRECTIONAL',
     $ ' REFLECTIVITY :', /, (10X,10F9.5) )
1100  FORMAT( ' NO INTERNAL SOURCES' )
1110  FORMAT('    BOTTOM SOURCE =', F10.2, '     TOP SOURCE =', F10.2 )
1120  FORMAT( ' USES DELTA-M METHOD' )
1130  FORMAT( ' DOES NOT USE DELTA-M METHOD' )
1140  FORMAT( ' CALCULATE FLUXES AND AZIM-AVERAGED INTENSITIES ONLY' )
1150  FORMAT( ' CALCULATE FLUXES AND INTENSITIES' )
1160  FORMAT( ' RELATIVE CONVERGENCE CRITERION FOR AZIMUTH SERIES =',
     $   1P,E11.2 )
1170  FORMAT( ' SETS RADIATION = 0 BELOW ABSORPTION OPTICAL DEPTH 10' )
1190  FORMAT( /, 37X, '<------------- DELTA-M --------------->', /,
     $'                   TOTAL    SINGLE                           ',
     $               'TOTAL    SINGLE', /,
     $'       OPTICAL   OPTICAL   SCATTER   TRUNCATED   ',
     $   'OPTICAL   OPTICAL   SCATTER    ASYMM', /,
     $'         DEPTH     DEPTH    ALBEDO    FRACTION     ',
     $     'DEPTH     DEPTH    ALBEDO   FACTOR' )
1191  FORMAT( /, 37X, '<------------- DELTA-M --------------->', /,
     $'                   TOTAL    SINGLE                           ',
     $               'TOTAL    SINGLE', /,
     $'       OPTICAL   OPTICAL   SCATTER   TRUNCATED   ',
     $   'OPTICAL   OPTICAL   SCATTER    ASYMM', /,
     $'         DEPTH     DEPTH    ALBEDO    FRACTION     ',
     $     'DEPTH     DEPTH    ALBEDO   FACTOR' )
1200  FORMAT( I4, 2F10.4, F10.5, F12.5, 2F10.4, F10.5, F9.4 )
1300  FORMAT( I6, 10F11.6, /, (6X,10F11.6) )
C
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI,
     $                    MAXULV, MAXUMU)
C
C         PRINTS THE INTENSITY AT USER POLAR AND AZIMUTHAL ANGLES
C
C     ALL ARGUMENTS ARE DISORT INPUT OR OUTPUT VARIABLES
C
C+---------------------------------------------------------------------+
      implicit none
       integer MAXULV,MAXUMU,NTAU,NUMU,NPHI,jmax,j,lu
       integer jmin,lenfmt,np,npass,iu
        REAL*8   PHI(*), UMU(*), UTAU(*), UU( MAXUMU, MAXULV, * )
C
C
      WRITE ( 6, '(//,A)' )
     $         ' *********  I N T E N S I T I E S  *********'
      LENFMT = 10
      NPASS = 1 + NPHI / LENFMT
      IF ( MOD(NPHI,LENFMT) .EQ. 0 )  NPASS = NPASS - 1
      DO 10  LU = 1, NTAU
         DO 10  NP = 1, NPASS
            JMIN = 1 + LENFMT * (NP-1)
            JMAX = MIN0( LENFMT*NP, NPHI )
            WRITE( 6,101 )  ( PHI(J), J = JMIN, JMAX )
            DO 10  IU = 1, NUMU
               IF( IU.EQ.1 )  WRITE( 6,102 )  UTAU(LU), UMU(IU),
     $           ( UU( IU,LU,J ), J = JMIN, JMAX )
               IF( IU.GT.1 )  WRITE( 6,103 )  UMU(IU),
     $           ( UU( IU,LU,J ), J = JMIN, JMAX )
10    CONTINUE
C
      RETURN
C
101   FORMAT( /, 3X,'          POLAR   AZIMUTH ANGLES (DEGREES)',
     $        /, 3X,'OPTICAL   ANGLE',
     $        /, 3X,' DEPTH   COSINE', 10F11.2 )
102   FORMAT( F10.4, F8.4, 1P,10E11.3 )
103   FORMAT( 10X,   F8.4, 1P,10E11.3 )
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  QGAUSN( M, GMU, GWT )
C
C       COMPUTE WEIGHTS AND ABSCISSAE FOR ORDINARY GAUSSIAN QUADRATURE
C       (NO WEIGHT FUNCTION INSIDE INTEGRAL) ON THE INTERVAL (0,1)
C
C       REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
C                   Integration, Academic Press, New York, pp. 87, 1975.
C
C          METHOD:  Compute the abscissae as roots of the Legendre
C                   Polynomial P-SUB-N using a cubically convergent
C                   refinement of Newton's method.  Compute the
C                   weights from EQ. 2.7.3.8 of Davis/Rabinowitz.
C
C        ACCURACY:  at least 13 significant digits
C
C
C  I N P U T :    M       ORDER OF QUADRATURE RULE
C
C  O U T P U T :  GMU(I)  I = 1 TO M,    ARRAY OF ABSCISSAE
C                 GWT(I)  I = 1 TO M,    ARRAY OF WEIGHTS
C
C  I N T E R N A L    V A R I A B L E S:
C
C    PM2,PM1,P : 3 SUCCESSIVE LEGENDRE POLYNOMIALS
C    PPR       : DERIVATIVE OF LEGENDRE POLYNOMIAL
C    P2PRI     : 2ND DERIVATIVE OF LEGENDRE POLYNOMIAL
C    TOL       : CONVERGENCE CRITERION FOR LEGENDRE POLY ROOT ITERATION
C    X,XI      : SUCCESSIVE ITERATES IN CUBICALLY-
C                CONVERGENT VERSION OF NEWTON'S METHOD
C                ( SEEKING ROOTS OF LEGENDRE POLYNOMIAL )
C+---------------------------------------------------------------------+
      implicit none
        REAL*8     CONA, GMU(*), GWT(*), PI, T
      INTEGER  LIM, M, NP1,K,NN
*      DOUBLE   PRECISION  D1MACH
       real*8 EN, NNP1, P, PM1, PM2, PPR, P2PRI, PROD,
     $                    TMP, TOL, X, XI, one, two
*      real               EN, NNP1, P, PM1, PM2, PPR, P2PRI, PROD,
*     $                    TMP, TOL, X, XI, one, two
      DATA     PI / 0.0 /
      save one,two,tol,pi
C
C
      IF ( PI.EQ.0.0 )  THEN
         PI = 2. * ASIN(1.0)
         TOL = 1D-13            ! 10. * D1MACH(3)
         one=1.
         two=2.
      END IF
C
      IF ( M.LE.1 )  THEN
         M = 1
         GMU( 1 ) = 0.5
         GWT( 1 ) = 1.0
         RETURN
      END IF
C
      EN   = M
      NP1  = M + 1
      NNP1 = M * NP1
      CONA = FLOAT( M-1 ) / ( 8 * M**3 )
C                                        ** INITIAL GUESS FOR K-TH ROOT
C                                        ** OF LEGENDRE POLYNOMIAL, FROM
C                                        ** DAVIS/RABINOWITZ (2.7.3.3A)
      LIM  = M / 2
      DO 30  K = 1, LIM
         T = ( 4*K - 1 ) * PI / ( 4*M + 2 )
         X = COS ( T + CONA / TAN( T ) )
C                                        ** RECURSION RELATION FOR
C                                        ** LEGENDRE POLYNOMIALS
10       PM2 = one
         PM1 = X
         DO 20 NN = 2, M
            P   = ( ( 2*NN - 1 ) * X * PM1 - ( NN-1 ) * PM2 ) / NN
            PM2 = PM1
            PM1 = P
20       CONTINUE
C
         TMP   = one / ( one - X**2 )
         PPR   = EN * ( PM2 - X * P ) * TMP
         P2PRI = ( two * X * PPR - NNP1 * P ) * TMP
         XI    = X - ( P / PPR ) * ( one +
     $               ( P / PPR ) * P2PRI / ( two * PPR ) )
C
C                                              ** CHECK FOR CONVERGENCE
         IF ( ABS(XI-X) .GT. TOL ) THEN
            X = XI
            GO TO 10
         END IF
C                          ** ITERATION FINISHED--CALC. WEIGHTS,
C                          ** ABSCISSAE FOR (-1,1)
         GMU( K ) = - X
         GWT( K ) = two / ( TMP * ( EN * PM2 )**2 )
         GMU( NP1 - K ) = - GMU( K )
         GWT( NP1 - K ) =   GWT( K )
30    CONTINUE
C                                    ** SET MIDDLE ABSCISSA AND WEIGHT
C                                    ** FOR RULES OF ODD ORDER
      IF ( MOD( M,2 ) .NE. 0 )  THEN
         GMU( LIM + 1 ) = 0.0
         PROD = one
         DO 40 K = 3, M, 2
            PROD = PROD * K / ( K-1 )
40       CONTINUE
         GWT( LIM + 1 ) = two / PROD**2
      END IF
C                                        ** CONVERT FROM (-1,1) TO (0,1)
      DO 50  K = 1, M
         GMU( K ) = 0.5 * GMU( K ) + 0.5
         GWT( K ) = 0.5 * GWT( K )
50    CONTINUE
C
      RETURN
      END
      REAL*8 FUNCTION  RATIO( A, B )
      implicit none
      real*8 A,B
C
C        CALCULATE RATIO  A/B  WITH OVER- AND UNDER-FLOW PROTECTION
C
         IF ( ABS(A).LT.1.0E-8 .AND. ABS(B).LT.1.0E-8 )  THEN
            RATIO = 1.0
         ELSE IF ( B.EQ.0.0 )  THEN
            RATIO = 1.E+20
         ELSE
            RATIO = A / B
         END IF
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  SETDIS( CMU, CWT, DELTAM, DTAUC, EXPBEA, FBEAM, FLYR,
     $                    GL, HL, HLPR, LAMBER, LAYRU, LYRCUT, MAXCLY,
     $                    MAXUMU, MAXCMU, mxcly, MXCMU, NCUT, NLYR,
     $                    NTAU, NN,
     $                    NSTR, INTSRC, NUMU, ONLYFL, OPRIM, PMOM, SRC,
     $                    SRCU, SSALB, TAUC, TAUCPR, UTAU, UTAUPR, UMU, 
     $                    UMU0, USRTAU, USRANG)
C
C          PERFORM MISCELLANEOUS SETTING-UP OPERATIONS
C
C       ROUTINES CALLED:  ERRMSG, QGAUSN, ZEROIT
C
C       INPUT :  ALL ARE 'DISORT' INPUT VARIABLES (SEE DOC FILE)
C
C       OUTPUT:  NTAU,UTAU   IF USRTAU = FALSE
C                NUMU,UMU    IF USRANG = FALSE
C                CMU,CWT     COMPUTATIONAL POLAR ANGLES AND
C                               CORRESPONDING QUADRATURE WEIGHTS
C                EXPBEA      TRANSMISSION OF DIRECT BEAM
C                FLYR        TRUNCATED FRACTION IN DELTA-M METHOD
C                GL          PHASE FUNCTION LEGENDRE COEFFICIENTS MULTI-
C                              PLIED BY (2L+1) AND SINGLE-SCATTER ALBEDO
C                HLPR        LEGENDRE MOMENTS OF SURFACE BIDIRECTIONAL
C                              REFLECTIVITY, TIMES 2K+1
C                LAYRU       COMPUTATIONAL LAYER IN WHICH -UTAU- FALLS
C                LYRCUT      FLAG AS TO WHETHER RADIATION WILL BE ZEROED
C                              BELOW LAYER -NCUT-
C                NCUT        COMPUTATIONAL LAYER WHERE ABSORPTION
C                              OPTICAL DEPTH FIRST EXCEEDS -ABSCUT-
C                NN          NSTR / 2
C                OPRIM       DELTA-M-SCALED SINGLE-SCATTER ALBEDO
C                TAUCPR      DELTA-M-SCALED OPTICAL DEPTH
C                UTAUPR      DELTA-M-SCALED VERSION OF -UTAU-
C
      implicit none
      integer NUMU,NN,NTAU,MAXCLY,MAXUMU,maxcmu,mxcly
      integer ncut,nlyr,lu,lc,k,ls,iu,mxcmu,nstr,iq
      real*8 fbeam
      real*8 f,abscut,abstau,umu0,exxp
      LOGICAL  DELTAM, LAMBER, LYRCUT, INTSRC, ONLYFL, USRTAU, USRANG
      INTEGER  LAYRU(*)
      REAL*8     CMU(*), CWT(*), DTAUC(*), EXPBEA(0:*), FLYR(*),
     $         GL(0:MXCMU,*), HL(0:*), HLPR(0:*), OPRIM(*),
     $         PMOM(0:MAXCMU,*), SRC(3*MAXCLY,*), SRCU(3*MAXCLY,*),
     $         SSALB(*), TAUC(0:*), TAUCPR(0:*),
     $         UTAU(*), UTAUPR(*), UMU(*)
      DATA  ABSCUT / 10. /
C
C
      IF ( .NOT.USRTAU ) THEN
C                              ** SET OUTPUT LEVELS AT COMPUTATIONAL
C                              ** LAYER BOUNDARIES
         NTAU = NLYR + 1
         DO 30  LC = 0, NTAU-1
            UTAU(LC+1) = TAUC(LC)
30       CONTINUE
      END IF
C                        ** APPLY DELTA-M SCALING AND MOVE DESCRIPTION
C                        ** OF COMPUTATIONAL LAYERS TO LOCAL VARIABLES
      EXPBEA( 0 ) = 1.0
      CALL  ZEROIT( TAUCPR(0), MXCLY+1 )
      CALL  ZEROIT( EXPBEA(1), MXCLY )
      CALL  ZEROIT( FLYR, MXCLY )
      CALL  ZEROIT( GL, (MXCMU+1)*MXCLY )
      CALL  ZEROIT( OPRIM, MXCLY )
      ABSTAU = 0.0
      DO  60  LC = 1, NLYR
         PMOM(0,LC) = 1.0
         IF ( ABSTAU.LT.ABSCUT )  NCUT = LC
         ABSTAU = ABSTAU + ( 1. - SSALB(LC) ) * DTAUC(LC)
C
         IF ( .NOT.DELTAM )  THEN
            OPRIM(LC) = SSALB(LC)
            TAUCPR(LC) = TAUC(LC)
            DO 40  K = 0, NSTR-1
               GL(K,LC) = (2*K+1) * OPRIM(LC) * PMOM(K,LC)
 40         CONTINUE
            F = 0.0
         ELSE
C                                    ** DO DELTA-M TRANSFORMATION
            F = PMOM( NSTR,LC )
            OPRIM(LC) = SSALB(LC) * ( 1. - F ) / ( 1. - F * SSALB(LC) )
            TAUCPR(LC) = TAUCPR(LC-1) + ( 1. - F*SSALB(LC) ) * DTAUC(LC)
            DO 50  K = 0, NSTR-1
               GL(K,LC) = (2*K+1) * OPRIM(LC) * (PMOM(K,LC)-F) / (1.-F)
 50         CONTINUE
            LS = 3*(LC-1) + 1
            DO 51 IQ = 1, NSTR
               SRC(LS,IQ)   = SRC(LS,IQ)   / ( 1. - F * SSALB(LC) )
               SRC(LS+1,IQ) = SRC(LS+1,IQ) / ( 1. - F * SSALB(LC) )
               SRC(LS+2,IQ) = SRC(LS+2,IQ) / ( 1. - F * SSALB(LC) )
51          CONTINUE
            DO 52 IU = 1, NUMU
               SRCU(LS,IU)   = SRCU(LS,IU)   / ( 1. - F * SSALB(LC) )
               SRCU(LS+1,IU) = SRCU(LS+1,IU) / ( 1. - F * SSALB(LC) )
               SRCU(LS+2,IU) = SRCU(LS+2,IU) / ( 1. - F * SSALB(LC) )
52          CONTINUE
         ENDIF
C
         FLYR(LC) = F
         EXPBEA(LC) = 0.0
         IF ( FBEAM.GT.0.0 )  EXPBEA(LC) = EXxP( - TAUCPR(LC) / UMU0 )
60    CONTINUE
C                      ** IF NO INTERNAL SOURCES, CUT OFF MEDIUM BELOW
C                      ** ABSORPTION OPTICAL DEPTH = ABSCUT ( NOTE THAT
C                      ** DELTA-M TRANSFORMATION LEAVES ABSORPTION
C                      ** OPTICAL DEPTH INVARIANT ).  NOT WORTH THE
C                      ** TROUBLE FOR ONE-LAYER PROBLEMS, THOUGH.
      LYRCUT = .FALSE.
      IF ( ABSTAU.GE.ABSCUT .AND. .NOT. INTSRC
     $     .AND. NLYR.GT.1 )  LYRCUT =.TRUE.
      IF ( .NOT.LYRCUT )  NCUT = NLYR
C
C                             ** SET ARRAYS DEFINING LOCATION OF USER
C                             ** OUTPUT LEVELS WITHIN DELTA-M-SCALED
C                             ** COMPUTATIONAL MESH
      DO 90  LU = 1, NTAU
         DO 70 LC = 1, NLYR
            IF ( UTAU(LU).GE.TAUC(LC-1) .AND. UTAU(LU).LE.TAUC(LC) )
     $           GO TO 80
70       CONTINUE
         LC = NLYR
C
80       UTAUPR(LU) = UTAU(LU)
         IF(DELTAM) UTAUPR(LU) = TAUCPR(LC-1) + (1.-SSALB(LC)*FLYR(LC))
     $                                        * (UTAU(LU) - TAUC(LC-1))
         LAYRU(LU) = LC
90    CONTINUE
C                      ** CALCULATE COMPUTATIONAL POLAR ANGLE COSINES
C                      ** AND ASSOCIATED QUADRATURE WEIGHTS FOR GAUSSIAN
C                      ** QUADRATURE ON THE INTERVAL (0,1) (UPWARD)
      NN = NSTR / 2
      CALL  QGAUSN( NN, CMU, CWT )
C                                  ** DOWNWARD (NEG) ANGLES AND WEIGHTS
      DO 100  IQ = 1, NN
         CMU(IQ+NN) = - CMU(IQ)
         CWT(IQ+NN) =   CWT(IQ)
100   CONTINUE
C
      IF ( FBEAM.GT.0.0 )  THEN
C                               ** COMPARE BEAM ANGLE TO COMPU'L ANGLES
         DO 110  IQ = 1, NN
            IF ( ABS(UMU0-CMU(IQ))/UMU0 .LT. 1.E-4 )  CALL ERRMSG
     $         ( 'SETDIS--BEAM ANGLE=COMPUTATIONAL ANGLE; CHANGE NSTR',
     $            .TRUE. )
  110    CONTINUE
      END IF
C
      IF ( .NOT.USRANG .OR. (ONLYFL .AND. MAXUMU.GE.NSTR) )  THEN
C
C                                   ** SET OUTPUT POLAR ANGLES TO
C                                   ** COMPUTATIONAL POLAR ANGLES
            NUMU = NSTR
            DO 120  IU = 1, NN
               UMU(IU) = - CMU(NN+1-IU)
120         CONTINUE
            DO 121  IU = NN+1, NSTR
               UMU(IU) = CMU(IU-NN)
121         CONTINUE
      END IF
C
      IF ( .NOT.LYRCUT .AND. .NOT.LAMBER )  THEN
         DO 160  K = 0, NSTR
            HLPR(K) = (2*K+1) * HL(K)
160      CONTINUE
      END IF
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  SETMTX( BDR, CBAND, CMU, CWT, DELM0, GC, KK, LAMBER,
     $                    LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI,
     $                    NN, NSTR, TAUCPR, WK )
C
C        CALCULATE COEFFICIENT MATRIX FOR THE SET OF EQUATIONS
C        OBTAINED FROM THE BOUNDARY CONDITIONS AND THE CONTINUITY-
C        OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS;  STORE IN THE
C        SPECIAL BANDED-MATRIX FORMAT REQUIRED BY LINPACK ROUTINES
C
C     ROUTINES CALLED:  ZEROIT
C
C     I N P U T      V A R I A B L E S:
C
C       BDR      :  SURFACE BIDIRECTIONAL REFLECTIVITY
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       DELM0    :  KRONECKER DELTA, DELTA-SUB-M0
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       NN       :  NUMBER OF STREAMS IN A HEMISPHERE (NSTR/2)
C       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T     V A R I A B L E S:
C
C       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
C                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
C                   BY LINPACK SOLUTION ROUTINES
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
C
C   I N T E R N A L    V A R I A B L E S:
C
C       IROW     :  POINTS TO ROW IN  -CBAND-
C       JCOL     :  POINTS TO POSITION IN LAYER BLOCK
C       LDA      :  ROW DIMENSION OF -CBAND-
C       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
C       NSHIFT   :  FOR POSITIONING NUMBER OF ROWS IN BAND STORAGE
C       WK       :  TEMPORARY STORAGE FOR 'EXP' EVALUATIONS
C ---------------------------------------------------------------------+
      implicit none
      integer nstr,nn,nnlyri,ncol,ncut,nshift,nncol,mi,mi9m2
      integer ncd,lda,lc,k,jcol,jq,iq,irow,mxcmu
      real*8 sum,exxp,expa,delm0
        LOGICAL LAMBER, LYRCUT
      REAL*8    BDR( MI,0:* ), CBAND( MI9M2,NNLYRI ), CMU(*), CWT(*),
     $        GC( MXCMU,MXCMU,* ), KK( MXCMU,* ), TAUCPR( 0:* ), WK(*)
C
C
      CALL  ZEROIT( CBAND, MI9M2*NNLYRI )
      NCD    = 3*NN - 1
      LDA    = 3*NCD + 1
      NSHIFT = LDA - 2*NSTR + 1
      NCOL   = 0
C                         ** USE CONTINUITY CONDITIONS OF EQ. STWJ(17)
C                         ** TO FORM COEFFICIENT MATRIX IN STWJ(20);
C                         ** EMPLOY SCALING TRANSFORMATION STWJ(22)
      DO 30  LC = 1, NCUT
C
         DO 4  IQ = 1, NN
            WK(IQ) = EXxP(KK(IQ,LC) * (TAUCPR(LC) - TAUCPR(LC-1)) )
 4       CONTINUE
C
         JCOL = 0
         DO 10  IQ = 1, NN
            NCOL = NCOL + 1
            IROW = NSHIFT - JCOL
            DO 5  JQ = 1, NSTR
               CBAND(IROW+NSTR,NCOL) =   GC(JQ,IQ,LC)
               CBAND(IROW,     NCOL) = - GC(JQ,IQ,LC) * WK(IQ)
               IROW = IROW + 1
 5          CONTINUE
            JCOL = JCOL + 1
10       CONTINUE
C
         DO 20  IQ = NN+1, NSTR
            NCOL = NCOL + 1
            IROW = NSHIFT - JCOL
            DO 15  JQ = 1, NSTR
               CBAND(IROW+NSTR,NCOL) =   GC(JQ,IQ,LC) * WK(NSTR+1-IQ)
               CBAND(IROW,     NCOL) = - GC(JQ,IQ,LC)
               IROW = IROW + 1
15          CONTINUE
            JCOL = JCOL + 1
20       CONTINUE
C
30    CONTINUE
C                  ** USE TOP BOUNDARY CONDITION OF STWJ(20A) FOR
C                  ** FIRST LAYER
      JCOL = 0
      DO 40  IQ = 1, NN
         EXPA = EXxP( KK(IQ,1) * TAUCPR(1) )
         IROW = NSHIFT - JCOL + NN
         DO 35  JQ = NN, 1, -1
            CBAND(IROW,JCOL+1) = GC(JQ,IQ,1) * EXPA
            IROW = IROW+1
35       CONTINUE
         JCOL = JCOL+1
40    CONTINUE
C
      DO 50  IQ = NN+1, NSTR
         IROW = NSHIFT - JCOL + NN
         DO 45  JQ = NN, 1, -1
            CBAND(IROW,JCOL+1) = GC(JQ,IQ,1)
            IROW = IROW+1
45       CONTINUE
         JCOL = JCOL+1
50    CONTINUE
C                           ** USE BOTTOM BOUNDARY CONDITION OF
C                           ** STWJ(20C) FOR LAST LAYER
      NNCOL = NCOL - NSTR
      JCOL  = 0
      DO 70  IQ = 1, NN
         NNCOL = NNCOL + 1
         IROW  = NSHIFT - JCOL + NSTR
C
         DO 60  JQ = NN+1, NSTR
            IF ( LYRCUT .OR. (LAMBER .AND. DELM0.EQ.0) ) THEN
C
C                          ** NO AZIMUTHAL-DEPENDENT INTENSITY IF LAM-
C                          ** BERT SURFACE; NO INTENSITY COMPONENT IF
C                          ** TRUNCATED BOTTOM LAYER
C
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT)
            ELSE
               SUM = 0.0
               DO 55  K = 1, NN
                  SUM = SUM + CWT(K) * CMU(K) * BDR(JQ-NN,K)
     $                        * GC(NN+1-K,IQ,NCUT)
55             CONTINUE
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT) - (1.+DELM0) * SUM
            END IF
C
            IROW = IROW + 1
60       CONTINUE
         JCOL = JCOL + 1
70    CONTINUE
C
      DO 90  IQ = NN+1, NSTR
         NNCOL = NNCOL + 1
         IROW  = NSHIFT - JCOL + NSTR
         EXPA = WK(NSTR+1-IQ)
C
         DO 80  JQ = NN+1, NSTR
C
            IF ( LYRCUT .OR. (LAMBER .AND. DELM0.EQ.0) ) THEN
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT) * EXPA
            ELSE
               SUM = 0.0
               DO 75  K = 1, NN
                  SUM = SUM + CWT(K) * CMU(K) * BDR(JQ-NN,K)
     $                        * GC(NN+1-K,IQ,NCUT)
75             CONTINUE
               CBAND(IROW,NNCOL) = ( GC(JQ,IQ,NCUT)
     $                               - (1.+DELM0) * SUM ) * EXPA
            END IF
C
            IROW = IROW + 1
80       CONTINUE
         JCOL = JCOL + 1
90    CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MI, MAZ,
     $                    MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL,
     $                    KK, GC, AAD, WKD, EVECCD, EVALD )
C
C         SOLVES EIGENVALUE/VECTOR PROBLEM NECESSARY TO CONSTRUCT
C         HOMOGENEOUS PART OF DISCRETE ORDINATE SOLUTION; STWJ(8B)
C         ** NOTE ** EIGENVALUE PROBLEM IS DEGENERATE WHEN SINGLE
C                    SCATTERING ALBEDO = 1;  PRESENT WAY OF DOING IT
C                    SEEMS NUMERICALLY MORE STABLE THAN ALTERNATIVE
C                    METHODS THAT WE TRIED
C
C     ROUTINES CALLED:  ASYMTX
C
C   I N P U T     V A R I A B L E S:
C
C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
C                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO)
C       CMU    :  COMPUTATIONAL POLAR ANGLE COSINES
C       CWT    :  WEIGHTS FOR QUADRATURE OVER POLAR ANGLE COSINE
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       NN     :  HALF THE TOTAL NUMBER OF STREAMS
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE QUADRATURE ANGLES -CMU-
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T    V A R I A B L E S:
C
C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5); NEEDED IN SS(15&18)
C       EVAL   :  -NN- EIGENVALUES OF EQ. SS(12) ON RETURN FROM 'ASYMTX'
C                    BUT THEN SQUARE ROOTS TAKEN
C       EVECC  :  -NN- EIGENVECTORS  (G+) - (G-)  ON RETURN
C                    FROM 'ASYMTX' ( COLUMN J CORRESPONDS TO -EVAL(J)- )
C                    BUT THEN  (G+) + (G-)  IS CALCULATED FROM SS(10),
C                    G+  AND  G-  ARE SEPARATED, AND  G+  IS STACKED ON
C                    TOP OF  G-  TO FORM -NSTR- EIGENVECTORS OF SS(7)
C       GC     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVECTORS, BUT
C                    IN AN ORDER CORRESPONDING TO -KK-
C       KK     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVALUES OF SS(7),
C                    BUT RE-ORDERED WITH NEGATIVE VALUES FIRST ( SQUARE
C                    ROOTS OF -EVAL- TAKEN AND NEGATIVES ADDED )
C
C   I N T E R N A L   V A R I A B L E S:
C
C       AMB,APB :  MATRICES (ALPHA-BETA), (ALPHA+BETA) IN REDUCED
C                    EIGENVALUE PROBLEM
C       ARRAY   :  COMPLETE COEFFICIENT MATRIX OF REDUCED EIGENVALUE
C                    PROBLEM: (ALFA+BETA)*(ALFA-BETA)
C       GPPLGM  :  (G+) + (G-) (CF. EQS. SS(10-11))
C       GPMIGM  :  (G+) - (G-) (CF. EQS. SS(10-11))
C       WK      :  SCRATCH ARRAY REQUIRED BY 'ASYMTX'
C+---------------------------------------------------------------------+
      implicit none
      integer kq,jq,ier,iq,l,nstr,nn,mxcmu,maz,mi
      real*8 sum,gpplgm,gpmigm,alpha,beta
      REAL*8    AMB( MI,* ), APB( MI,* ), ARRAY( MI,* ), CC( MXCMU,* ),
     $        CMU(*), CWT(*), EVAL(*), EVECC( MXCMU,* ), GC( MXCMU,* ),
     $        GL(0:*), KK(*),  YLMC( 0:MXCMU,* )
      real*8   EVECCD( MI,* ), EVALD(*), WKD(*), AAD( MI,* )
C
C
C                             ** CALCULATE QUANTITIES IN EQS. SS(5-6)
      DO 40 IQ  = 1, NN
C
         DO 20  JQ = 1, NSTR
            SUM = 0.0
            DO 10  L = MAZ, NSTR-1
               SUM = SUM + GL(L) * YLMC(L,IQ) * YLMC(L,JQ)
10          CONTINUE
            CC(IQ,JQ) = 0.5 * SUM * CWT(JQ)
20       CONTINUE
C
         DO 30  JQ = 1, NN
C                             ** FILL REMAINDER OF ARRAY USING SYMMETRY
C                             ** RELATIONS  C(-MUI,MUJ) = C(MUI,-MUJ)
C                             ** AND        C(-MUI,-MUJ) = C(MUI,MUJ)
C
            CC(IQ+NN,JQ) = CC(IQ,JQ+NN)
            CC(IQ+NN,JQ+NN) = CC(IQ,JQ)
C                                      ** GET FACTORS OF COEFF. MATRIX
C                                      ** OF REDUCED EIGENVALUE PROBLEM
            ALPHA =   CC(IQ,JQ) / CMU(IQ)
            BETA = CC(IQ,JQ+NN) / CMU(IQ)
            AMB(IQ,JQ) = ALPHA - BETA
            APB(IQ,JQ) = ALPHA + BETA
30       CONTINUE
         AMB(IQ,IQ) = AMB(IQ,IQ) - 1.0 / CMU(IQ)
         APB(IQ,IQ) = APB(IQ,IQ) - 1.0 / CMU(IQ)
C
40    CONTINUE
C                      ** FINISH CALCULATION OF COEFFICIENT MATRIX OF
C                      ** REDUCED EIGENVALUE PROBLEM:  GET MATRIX
C                      ** PRODUCT (ALFA+BETA)*(ALFA-BETA); SS(12)
      DO 70  IQ = 1, NN
         DO 70  JQ = 1, NN
            SUM = 0.
            DO 60  KQ = 1, NN
               SUM = SUM + APB(IQ,KQ) * AMB(KQ,JQ)
60          CONTINUE
            ARRAY(IQ,JQ) = SUM
70    CONTINUE
C                      ** FIND (REAL) EIGENVALUES AND EIGENVECTORS
C
      CALL  ASYMTX( ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, 
     .          AAD, EVECCD, EVALD, WKD)
C
      IF ( IER.GT.0 )  THEN
         WRITE( 6, '(//,A,I4,A)' ) ' ASYMTX--EIGENVALUE NO. ', IER,
     $     '  DIDNT CONVERGE.  LOWER-NUMBERED EIGENVALUES WRONG.'
         CALL  ERRMSG( 'ASYMTX--CONVERGENCE PROBLEMS', .TRUE. )
      END IF
C
      DO 75  IQ = 1, NN
         EVAL(IQ) = SQRT( ABS( EVAL(IQ) ) )
         KK( IQ+NN ) = EVAL(IQ)
C                                             ** ADD NEGATIVE EIGENVALUE
         KK( NN+1-IQ ) = - EVAL(IQ)
75    CONTINUE
C                          ** FIND EIGENVECTORS (G+) + (G-) FROM SS(10)
C                          ** AND STORE TEMPORARILY IN -APB- ARRAY
      DO 90  JQ = 1, NN
         DO 90  IQ = 1, NN
            SUM = 0.
            DO 80  KQ = 1,NN
               SUM = SUM + AMB(IQ,KQ) * EVECC(KQ,JQ)
80          CONTINUE
            APB(IQ,JQ) = SUM / EVAL(JQ)
90    CONTINUE
C
      DO 100  JQ = 1, NN
         DO 100  IQ = 1, NN
            GPPLGM = APB(IQ,JQ)
            GPMIGM = EVECC(IQ,JQ)
C                                ** RECOVER EIGENVECTORS G+,G- FROM
C                                ** THEIR SUM AND DIFFERENCE; STACK THEM
C                                ** TO GET EIGENVECTORS OF FULL SYSTEM
C                                ** SS(7) (JQ = EIGENVECTOR NUMBER)
C
            EVECC(IQ,      JQ) = 0.5 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,   JQ) = 0.5 * ( GPPLGM - GPMIGM )
C
C                                ** EIGENVECTORS CORRESPONDING TO
C                                ** NEGATIVE EIGENVALUES (CORRESP. TO
C                                ** REVERSING SIGN OF 'K' IN SS(10) )
            GPPLGM = - GPPLGM
            EVECC(IQ,   JQ+NN) = 0.5 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,JQ+NN) = 0.5 * ( GPPLGM - GPMIGM )
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )
100   CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  SOLVE0(B, BDR, BEM, BSRC,CBAND, CMU, CWT, EXPBEA,
     $                    FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT,
     $                    MAZ, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR,
     $                    NNLYRI, PI, TSRC, TAUCPR, UMU0,  ZZ,
     $                    ZPLK0, ZPLK1, ZPLKA, ipk,dcband,dz,db)
C
C        CONSTRUCT RIGHT-HAND SIDE VECTOR -B- FOR GENERAL BOUNDARY
C        CONDITIONS STWJ(17) AND SOLVE SYSTEM OF EQUATIONS OBTAINED
C        FROM THE BOUNDARY CONDITIONS AND THE
C        CONTINUITY-OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS.
C        THERMAL EMISSION CONTRIBUTES ONLY IN AZIMUTHAL INDEPENDENCE.
C
C     ROUTINES CALLED:  SGBCO, SGBSL, ZEROIT
C
C     I N P U T      V A R I A B L E S:
C
C       BDR      :  SURFACE BIDIRECTIONAL REFLECTIVITY
C       BEM      :  SURFACE BIDIRECTIONAL EMISSIVITY
C       BSRC     :  BOTTOM BOUNDARY EMISSION
C       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
C                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
C                   BY LINPACK SOLUTION ROUTINES
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       EXPBEA   :  TRANSMISSION OF INCIDENT BEAM, EXxP(-TAUCPR/UMU0)
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       MAZ      :  ORDER OF AZIMUTHAL COMPONENT
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
C       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
C       TSRC     :  TOP BOUNDARY EMISSION
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       ZZ       :  BEAM SOURCE VECTORS IN EQ. SS(19)
C       ZPLK0    :  INTERNAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
C       ZPLK1    :  INTERNAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)
C       ZPLKA    :  INTERNAL SOURCE VECTORS -ZA-, ALFA IN EQ.
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T     V A R I A B L E S:
C
C       B        :  RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO
C                   *SGBSL*; RETURNS AS SOLUTION VECTOR OF EQ.
C                   SC(12), CONSTANTS OF INTEGRATION WITHOUT
C                   EXPONENTIAL TERM
C      LL        :  PERMANENT STORAGE FOR -B-, BUT RE-ORDERED
C
C   I N T E R N A L    V A R I A B L E S:
C
C       IPVT     :  INTEGER VECTOR OF PIVOT INDICES
C       IT       :  POINTER FOR POSITION IN  -B-
C       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL
C       RCOND    :  INDICATOR OF SINGULARITY FOR -CBAND-
C       Z        :  SCRATCH ARRAY REQUIRED BY *SGBCO*
C+---------------------------------------------------------------------+
C
      implicit none
      integer ipk,mi
      real*8 sum,umu0,tsrc,pi,bsrc,exxp
      integer ncd,jq,it,iq,j,lc,i,ipnt,nnlyri,nstr,nn,ncut
      integer ncol,maz,mxcmu,mi9m2
      real*8 fbeam
      LOGICAL  LAMBER, LYRCUT
      INTEGER  IPVT(*)
      REAL*8     B(*), BDR( MI,0:* ), BEM(*), CBAND( MI9M2,NNLYRI ),
     $         CMU(*), CWT(*), EXPBEA(0:*), LL( MXCMU,* ),
     $         TAUCPR( 0:* ), ZZ( MXCMU,* ), ZPLK0( MXCMU,* ),
     $         ZPLK1( MXCMU,* ), ZPLKA( * ) , fisot(-ipk:-1)
        real*8 dcband(mi9m2,nnlyri),dz(*),drcond,db(nnlyri)
C
C
      CALL  ZEROIT( B, NNLYRI )
      CALL  dZEROIT( dB, NNLYRI )
C                             ** CONSTRUCT -B-,  STWJ(20A,C) FOR
C                             ** PARALLEL BEAM + BOTTOM REFLECTION +
C                             ** EMISSION AT TOP AND/OR BOTTOM
C
      IF ( MAZ.GT.0 .AND. FBEAM.GT.0.0 )  THEN
C
C                                         ** AZIMUTH-DEPENDENT CASE
C                                         ** (NEVER CALLED IF FBEAM = 0)
         IF ( LYRCUT .OR. LAMBER ) THEN
C
C               ** NO AZIMUTHAL-DEPENDENT INTENSITY FOR LAMBERT SURFACE;
C               ** NO INTENSITY COMPONENT FOR TRUNCATED BOTTOM LAYER
C
            DO 10  IQ = 1, NN
C                                                     ** TOP BOUNDARY
               B(IQ) = - ZZ(NN+1-IQ,1)
C                                                  ** BOTTOM BOUNDARY
               B(NCOL-NN+IQ) = - ZZ(IQ+NN,NCUT) * EXPBEA(NCUT)
10          CONTINUE
C
         ELSE
C
            DO 20  IQ = 1, NN
               B(IQ) = - ZZ(NN+1-IQ,1)
C
               SUM = 0.
               DO 15  JQ = 1, NN
                  SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)
     $                        * ZZ(NN+1-JQ,NCUT) * EXPBEA(NCUT)
15             CONTINUE
               B(NCOL-NN+IQ) = SUM
               IF ( FBEAM.GT.0.0 )
     $              B(NCOL-NN+IQ) = SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI
     $                                 - ZZ(IQ+NN,NCUT) ) * EXPBEA(NCUT)
20          CONTINUE
         END IF
C                             ** CONTINUITY CONDITION FOR LAYER
C                             ** INTERFACES OF EQ. STWJ(20B)
         IT = NN
         DO 40  LC = 1, NCUT-1
            DO 30  IQ = 1, NSTR
               IT    = IT + 1
               B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)
30          CONTINUE
40       CONTINUE
C
      ELSE
C                                   ** AZIMUTH-INDEPENDENT CASE
         IF ( FBEAM.EQ.0.0 )  THEN
C
            DO 50 IQ = 1, NN
C                                      ** TOP BOUNDARY
C
               B(IQ) = - ZPLK0(NN+1-IQ,1) + FISOT(-iq) + TSRC
50          CONTINUE
C
            IF ( LYRCUT ) THEN
C                               ** NO INTENSITY COMPONENT FOR TRUNCATED
C                               ** BOTTOM LAYER
               DO 60 IQ = 1, NN
C                                      ** BOTTOM BOUNDARY
C
                  B(NCOL-NN+IQ) = 
     $                       -  EXxP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $               (ZPLK0(IQ+NN,NCUT)+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
60             CONTINUE
C
            ELSE
C
               DO 80 IQ = 1, NN
C
                  SUM = 0.
                  DO 70 JQ = 1, NN
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)
     $                          * ( 
     $                EXxP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $       (ZPLK0(NN+1-JQ,NCUT)+ZPLK1(NN+1-JQ,NCUT)*TAUCPR(NCUT)))
70                CONTINUE
                  B(NCOL-NN+IQ) = 2.*SUM + BEM(IQ) * BSRC
     $                        - EXxP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $               (ZPLK0(IQ+NN,NCUT)+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
80             CONTINUE
            END IF
c47          format(i5,8(1pe9.1))
C                             ** CONTINUITY CONDITION FOR LAYER
C                             ** INTERFACES, STWJ(20B)
            IT = NN
            DO 100  LC = 1, NCUT-1
               DO 90  IQ = 1, NSTR
                  IT    = IT + 1
                  B(IT) =
     $              +  EXxP(-ZPLKA(LC+1)*TAUCPR(LC))*
     $       (ZPLK0(IQ,LC+1)+ZPLK1(IQ,LC+1)*TAUCPR(LC))
     $              -  EXxP(-ZPLKA(LC)*TAUCPR(LC))*
     $       (ZPLK0(IQ,LC)+ZPLK1(IQ,LC)*TAUCPR(LC))
90             CONTINUE
100         CONTINUE
C
         ELSE
C
            DO 150 IQ = 1, NN
               B(IQ) = - ZZ(NN+1-IQ,1) - ZPLK0(NN+1-IQ,1) +
     .                  FISOT(-iq) + TSRC
150         CONTINUE
C
            IF ( LYRCUT ) THEN
C
               DO 160 IQ = 1, NN
                  B(NCOL-NN+IQ) = - ZZ(IQ+NN,NCUT) * EXPBEA(NCUT)
     $                        - EXxP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $               (ZPLK0(IQ+NN,NCUT)+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
160            CONTINUE
C
            ELSE
C
               DO 180 IQ = 1, NN
C
                  SUM = 0.
                  DO 170 JQ = 1, NN
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)
     $                          * ( ZZ(NN+1-JQ,NCUT) * EXPBEA(NCUT)
     $              +  EXxP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $       (ZPLK0(NN+1-JQ,NCUT)+ZPLK1(NN+1-JQ,NCUT)*TAUCPR(NCUT)))
170               CONTINUE
                  B(NCOL-NN+IQ) = 2.*SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI
     $                                 - ZZ(IQ+NN,NCUT) ) * EXPBEA(NCUT)
     $                            + BEM(IQ) * BSRC
     $                       -  EXxP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $               (ZPLK0(IQ+NN,NCUT)+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
180            CONTINUE
            END IF
C
            IT = NN
            DO 200  LC = 1, NCUT-1
               DO 190  IQ = 1, NSTR
                  IT    = IT + 1
                  B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)
     $              +  EXxP(-ZPLKA(LC+1)*TAUCPR(LC))*
     $       (ZPLK0(IQ,LC+1)+ZPLK1(IQ,LC+1)*TAUCPR(LC))
     $              -  EXxP(-ZPLKA(LC)*TAUCPR(LC))*
     $       (ZPLK0(IQ,LC)+ZPLK1(IQ,LC)*TAUCPR(LC))
190            CONTINUE
200         CONTINUE
C
         END IF
C
      END IF
C                     ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITION
C                     ** OF BAND MATRIX -CBAND- AND TEST IF IT IS NEARLY
C                     ** SINGULAR (NOTE: -CBAND- IS DESTROYED)
C                     ** (-CBAND- IS IN LINPACK PACKED FORMAT)
        do i=1,nnlyri
          db(i) = (b(i))
          do j=1,mi9m2
            dcband(j,i) = (cband(j,i))
          enddo
        enddo
      dRCOND = 0.0
      NCD = 3*NN - 1
      CALL  dGBCO( dCBAND, MI9M2, NCOL, NCD, NCD, IPVT, dRCOND, dZ )
*      write(6,'(a,1pe20.5)') ' rcond in SOLVE0=',rcond
      IF ( 1.0+dRCOND .EQ. 1.0 )  CALL  ERRMSG
     $   ( 'SOLVE0--SGBCO SAYS MATRIX NEAR SINGULAR',.FALSE.)
C
C                   ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -CBAND-
C                   ** AND R.H. SIDE(S) -B- AFTER -CBAND- HAS BEEN L-U
C                   ** DECOMPOSED.  SOLUTION IS RETURNED IN -B-.
C
      CALL  dGBSL( dCBAND, MI9M2, NCOL, NCD, NCD, IPVT, dB, 0 )
C
C                   ** ZERO -CBAND- (IT MAY CONTAIN 'FOREIGN'
C                   ** ELEMENTS UPON RETURNING FROM LINPACK);
C                   ** NECESSARY TO PREVENT ERRORS
C
      CALL  dZEROIT( dCBAND, MI9M2*NNLYRI )
C
      DO 220  LC = 1, NCUT
         IPNT = LC*NSTR - NN
         DO 220  IQ = 1, NN
            LL(NN+1-IQ,LC) = sngl(dB(IPNT+1-IQ))
            LL(IQ+NN,  LC) = sngl(dB(IQ+IPNT))
220   CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER,
     $                    MI, MAZ, MXCMU, MXUMU, NN, NUMU, NSTR, ONLYFL,
     $                    UMU, USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM,
     $                    RMU )
C
C       SPECIFIES USER'S SURFACE BIDIRECTIONAL PROPERTIES, STWJ(21)
C
C   I N P U T     V A R I A B L E S:
C
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0
C       HLPR   :  LEGENDRE MOMENTS OF SURFACE BIDIRECTIONAL REFLECTIVITY
C                    (WITH 2K+1 FACTOR INCLUDED)
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       NN     :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
C       YLM0   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE BEAM ANGLE
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
C                 AT THE QUADRATURE ANGLES
C       YLMU   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
C                 AT THE USER ANGLES
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C    O U T P U T     V A R I A B L E S:
C
C       BDR :  SURFACE BIDIRECTIONAL REFLECTIVITY (COMPUTATIONAL ANGLES)
C       RMU :  SURFACE BIDIRECTIONAL REFLECTIVITY (USER ANGLES)
C       BEM :  SURFACE DIRECTIONAL EMISSIVITY (COMPUTATIONAL ANGLES)
C       EMU :  SURFACE DIRECTIONAL EMISSIVITY (USER ANGLES)
C
C    I N T E R N A L     V A R I A B L E S:
C
C       DREF      DIRECTIONAL REFLECTIVITY
C       NMUG   :  NUMBER OF ANGLE COSINE QUADRATURE POINTS
C                 ON (0,1) FOR INTEGRATING BIDIRECTIONAL REFLECTIVITY
C                 TO GET DIRECTIONAL EMISSIVITY (IT IS NECESSARY TO USE
C                 A QUADRATURE SET DISTINCT FROM THE COMPUTATIONAL
C                 ANGLES, BECAUSE THE COMPUTATIONAL ANGLES MAY NOT BE
C                 DENSE ENOUGH -- I.E. 'NSTR' MAY BE TOO SMALL-- TO GIVE
C                 AN ACCURATE APPROXIMATION FOR THE INTEGRATION).
C       GMU    :  THE 'NMUG' ANGLE COSINE QUADRATURE POINTS ON (0,1)
C       GWT    :  THE 'NMUG' ANGLE COSINE QUADRATURE WEIGHTS ON (0,1)
C       YLMG   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
C                 AT THE 'NMUG' QUADRATURE ANGLES
C+---------------------------------------------------------------------+
      implicit none
      integer nmug,mi,maz,mxcmu,nn,numu,nstr,iq
      integer k,iu,jg,jq,maxstr,mxumu
      real*8 albedo,delm0,dref,sum,sgn
      real*8 fbeam

        LOGICAL  LAMBER, ONLYFL, USRANG
      REAL*8     BDR( MI,0:* ), BEM(*), EMU(*),
     $         HLPR(0:*), RMU( MXUMU,0:* ), UMU(*),
     $         YLM0(0:*), YLMC( 0:MXCMU,* ), YLMU( 0:MXCMU,* )
      PARAMETER  ( NMUG = 10, MAXSTR = 100 )
      LOGICAL  PASS1
      REAL*8     GMU( NMUG ), GWT( NMUG ), YLMG( 0:MAXSTR, NMUG )
      DATA  PASS1 / .TRUE. /
      save ylmg
C
C
      IF ( PASS1 )  THEN
         PASS1 = .FALSE.
         CALL QGAUSN( NMUG, GMU, GWT )
C
         CALL LEPOLY( NMUG, 0, MAXSTR, MAXSTR, GMU, YLMG)
C                       ** CONVERT LEGENDRE POLYS. TO NEGATIVE -GMU-
         SGN  = - 1.0
         DO 1  K = 0, MAXSTR
            SGN = - SGN
            DO 1  JG = 1, NMUG
               YLMG( K,JG ) = SGN * YLMG( K,JG )
 1       CONTINUE
C
      END IF
C
      CALL  ZEROIT( BDR, MI*(MI+1) )
      CALL  ZEROIT( BEM, MI )
C
      IF ( LAMBER .AND. MAZ.EQ.0 ) THEN
C
         DO 20 IQ = 1, NN
            BEM(IQ) = 1.0 - ALBEDO
            DO 20 JQ = 0, NN
               BDR(IQ,JQ) = ALBEDO
20       CONTINUE
C
      ELSE IF ( .NOT.LAMBER ) THEN
C                                  ** COMPUTE SURFACE BIDIRECTIONAL
C                                  ** PROPERTIES AT COMPUTATIONAL ANGLES
         DO 60 IQ = 1, NN
C
            DO 40 JQ = 1, NN
              SUM = 0.0
              DO 30 K = MAZ, NSTR-1
                 SUM = SUM + HLPR(K) * YLMC(K,IQ) * YLMC(K,JQ+NN)
30            CONTINUE
              BDR(IQ,JQ) = (2.-DELM0) * SUM
40          CONTINUE
C
            IF ( FBEAM.GT.0.0 )  THEN
               SUM = 0.0
               DO 50 K = MAZ, NSTR-1
                  SUM = SUM + HLPR(K) * YLMC(K,IQ) * YLM0(K)
50             CONTINUE
               BDR(IQ,0) = (2.-DELM0) * SUM
            ENDIF
C
60       CONTINUE
C
         IF ( MAZ.EQ.0 ) THEN
C
            IF ( NSTR.GT.MAXSTR )  CALL
     $           ERRMSG( 'SURFAC--PARAMETER MAXSTR TOO SMALL', .TRUE.
     $            )
C
C                              ** INTEGRATE BIDIRECTIONAL REFLECTIVITY
C                              ** AT REFLECTION POLAR ANGLES -CMU- AND
C                              ** INCIDENT ANGLES -GMU- TO GET
C                              ** DIRECTIONAL EMISSIVITY AT
C                              ** COMPUTATIONAL ANGLES -CMU-.
            DO 100  IQ = 1, NN
               DREF = 0.0
               DO 90  JG = 1, NMUG
                  SUM = 0.0
                  DO 80  K = 0, NSTR-1
                     SUM = SUM + HLPR(K) * YLMC(K,IQ) * YLMG(K,JG)
80                CONTINUE
                  DREF = DREF + 2.* GWT(JG) * GMU(JG) * SUM
90             CONTINUE
               BEM(IQ) = 1.0 - DREF
100         CONTINUE
C
         END IF
C
      END IF
C                                       ** COMPUTE SURFACE BIDIRECTIONAL
C                                       ** PROPERTIES AT USER ANGLES
C
      IF ( .NOT.ONLYFL .AND. USRANG )  THEN
C
         CALL  ZEROIT( EMU, MXUMU )
         CALL  ZEROIT( RMU, MXUMU*(MI+1) )
C
         DO 170 IU = 1, NUMU
            IF ( UMU(IU).GT.0.0 )  THEN
C
               IF ( LAMBER .AND. MAZ.EQ.0 )  THEN
                  DO 110 IQ = 0, NN
                     RMU(IU,IQ) = ALBEDO
110               CONTINUE
                  EMU(IU) = 1.0 - ALBEDO
C
               ELSE IF ( .NOT.LAMBER ) THEN
                  DO 130 IQ = 1, NN
                     SUM = 0.0
                     DO 120 K = MAZ, NSTR-1
                        SUM = SUM + HLPR(K) * YLMU(K,IU) * YLMC(K,IQ+NN)
120                  CONTINUE
                     RMU(IU,IQ) = (2.-DELM0) * SUM
130               CONTINUE
C
                  IF ( FBEAM.GT.0.0 )  THEN
                     SUM = 0.0
                     DO 140 K = MAZ, NSTR-1
                        SUM = SUM + HLPR(K) * YLMU(K,IU) * YLM0(K)
140                  CONTINUE
                     RMU(IU,0) = (2.-DELM0) * SUM
                  END IF
C
                  IF ( MAZ.EQ.0 ) THEN
C
C                               ** INTEGRATE BIDIRECTIONAL REFLECTIVITY
C                               ** AT REFLECTION ANGLES -UMU- AND
C                               ** INCIDENT ANGLES -GMU- TO GET
C                               ** DIRECTIONAL EMISSIVITY AT
C                               ** USER ANGLES -UMU-.
                     DREF = 0.0
                     DO 160 JG = 1, NMUG
                        SUM = 0.0
                        DO 150 K = 0, NSTR-1
                           SUM = SUM + HLPR(K) * YLMU(K,IU) * YLMG(K,JG)
150                     CONTINUE
                        DREF = DREF + 2.* GWT(JG) * GMU(JG) * SUM
160                  CONTINUE
C
                     EMU(IU) = 1.0 - DREF
                  END IF
C
               END IF
            END IF
170      CONTINUE
C
      END IF
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  TERPEV( CWT, EVECC, GL, GU, MAZ, MXCMU, MXUMU,
     $                    NN, NSTR, NUMU, WK, YLMC, YLMU )
C
C         INTERPOLATE EIGENVECTORS TO USER ANGLES; EQ SD(8)
C
      implicit none
      integer mxumu,nstr,numu,mxcmu,nn
      integer l,jq,maz,iu,iq
      real*8 sum
      REAL*8 CWT(*), EVECC( MXCMU,* ), GL(0:*), GU(  MXUMU,* ), WK(*),
     $      YLMC(  0:MXCMU,* ), YLMU(  0:MXCMU,* )
C
C
      DO 50  IQ = 1, NSTR
C
         DO 20  L = MAZ, NSTR-1
C                                       ** INNER SUM IN SD(8) TIMES ALL
C                                   ** FACTORS IN OUTER SUM BUT PLM(MU)
            SUM = 0.0
            DO 10  JQ = 1, NSTR
               SUM = SUM + CWT(JQ) * YLMC(L,JQ) * EVECC(JQ,IQ)
10          CONTINUE
            WK(L+1) = 0.5 * GL(L) * SUM
20       CONTINUE
C                                    ** FINISH OUTER SUM IN SD(8)
C                                    ** AND STORE EIGENVECTORS
         DO 40  IU = 1, NUMU
            SUM = 0.
            DO 30  L = MAZ, NSTR-1
               SUM = SUM + WK(L+1) * YLMU(L,IU)
30          CONTINUE
            IF ( IQ.LE.NN )  GU( IU, IQ+NN     ) = SUM
            IF ( IQ.GT.NN )  GU( IU, NSTR+1-IQ ) = SUM
40       CONTINUE
C
50    CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  TERPSO( CWT, DELM0, FBEAM, GL, MAZ, MXCMU, 
     $                    INTSRC, NUMU, NSTR,  PI, YLM0, YLMC,
     $                    YLMU, PSI0, PSI1, XR0U, XR1U, Z0, Z1, 
     $                    ZJ, ZBEAM, Z0U, Z1U)
C
C         INTERPOLATES SOURCE FUNCTIONS TO USER ANGLES
C
C    I N P U T      V A R I A B L E S:
C
C       CWT    :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0
C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
C                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO)
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       OPRIM  :  SINGLE SCATTERING ALBEDO
C       XR0    :  EXPANSION OF INTERNAL SOURCE FUNCTION
C       XR1    :  EXPANSION OF INTERNAL SOURCE FUNCTION EQS.SS(14-16)
C       XRA    :  EXPANSION OF INTERNAL SOURCE FUNCTION EQS.SS(14-16)
C       YLM0   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE BEAM ANGLE
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE QUADRATURE ANGLES
C       YLMU   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE USER ANGLES
C       Z0     :  SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)
C       Z1     :  SOLUTION VECTORS Z-SUB-ONE  OF EQ. SS(16)
C       ZJ     :  SOLUTION VECTOR CAPITAL -Z-SUB-ZERO AFTER SOLVING
C                 EQ. SS(19)
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C    O U T P U T     V A R I A B L E S:
C
C       ZBEAM  :  INCIDENT-BEAM SOURCE FUNCTION AT USER ANGLES
C       Z0U,Z1U:  COMPONENTS OF A EXPONENTIAL-LINEARN-IN-OPTICAL-
C         ZA      DEPTH-DEPENDENT SOURCE (APPROXIMATING THE
C                             INTERNAL EMISSION SOURCE)
C
C   I N T E R N A L       V A R I A B L E S:
C
C       PSI0,  :  SUM JUST AFTER SQUARE BRACKET IN  EQ. SD(9)
C        PSI1     WITH Z0 AND Z1 RESPECTIVELY
C+---------------------------------------------------------------------+
      IMPLICIT NONE
      INTEGER MXCMU,NUMU,NSTR,iu,jq,iq,maz
      real*8 sum1,sum,psum1,sum0,psum,psum0,fact,pi,delm0
      real*8 fbeam
        LOGICAL  INTSRC
      REAL*8     CWT(*), GL(0:*), PSI0(*), PSI1(*), XR0U(*), XR1U(*),
     $         YLM0(0:*), YLMC( 0:MXCMU,* ), YLMU( 0:MXCMU,*), 
     $         Z0(*), Z1(*), ZJ(*), ZBEAM(*), Z0U(*), Z1U(*)
C
C
      IF ( FBEAM.GT.0.0 )  THEN
C                                  ** BEAM SOURCE TERMS; EQ. SD(9)
         DO 20  IQ = MAZ, NSTR-1
            PSUM = 0.
            DO 10  JQ = 1, NSTR
               PSUM = PSUM + CWT(JQ) * YLMC(IQ,JQ) * ZJ(JQ)
10          CONTINUE
            PSI0(IQ+1) = 0.5 * GL(IQ) * PSUM
20       CONTINUE
C
         FACT = ( 2. - DELM0 ) * FBEAM / (4.0*PI)
         DO 40  IU = 1, NUMU
            SUM = 0.
            DO 30 IQ = MAZ, NSTR-1
               SUM = SUM + YLMU(IQ,IU) *
     $                    ( PSI0(IQ+1) + FACT * GL(IQ) * YLM0(IQ) )
30          CONTINUE
            ZBEAM(IU) = SUM
40       CONTINUE
      END IF
C
      IF ( INTSRC .AND. MAZ.EQ.0 )  THEN
C
C                                ** INTERNAL SOURCE TERMS, STWJ(27C)
         DO 80  IQ = MAZ, NSTR-1
            PSUM0 = 0.0
            PSUM1 = 0.0
            DO 70  JQ = 1, NSTR
               PSUM0 = PSUM0 + CWT(JQ) * YLMC(IQ,JQ) * Z0(JQ)
               PSUM1 = PSUM1 + CWT(JQ) * YLMC(IQ,JQ) * Z1(JQ)
 70         CONTINUE
            PSI0(IQ+1) = 0.5 * GL(IQ) * PSUM0
            PSI1(IQ+1) = 0.5 * GL(IQ) * PSUM1
 80       CONTINUE
C
          DO 100  IU = 1, NUMU
            SUM0 = 0.0
            SUM1 = 0.0
            DO 90   IQ = MAZ, NSTR-1
               SUM0 = SUM0 + YLMU(IQ,IU) *  PSI0(IQ+1)
               SUM1 = SUM1 + YLMU(IQ,IU) *  PSI1(IQ+1)
90          CONTINUE
            Z0U(IU) = SUM0 + XR0U(IU) 
            Z1U(IU) = SUM1 + XR1U(IU)
100      CONTINUE
C
      END IF
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  UPBEAM(  CC, CMU, DELM0, FBEAM, GL, IPVT, MAZ,
     $                    MXCMU, NN, NSTR, PI, UMU0, YLM0, YLMC,
     $                    ZZ,darray,dwk,dzj)
C
C         FINDS THE INCIDENT-BEAM PARTICULAR SOLUTION  OF SS(18)
C
C     ROUTINES CALLED:  SGECO, SGESL
C
C   I N P U T    V A R I A B L E S:
C
C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5)
C       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0
C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
C                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO)
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       YLM0   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE BEAM ANGLE
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE QUADRATURE ANGLES
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T    V A R I A B L E S:
C
C       ZJ     :  RIGHT-HAND SIDE VECTOR CAPITAL-X-SUB-ZERO IN SS(19);
C                 ALSO THE SOLUTION VECTOR CAPITAL-Z-SUB-ZERO
C                 AFTER SOLVING THAT SYSTEM
C       ZZ     :  PERMANENT STORAGE FOR -ZJ-, BUT RE-ORDERED
C
C   I N T E R N A L    V A R I A B L E S:
C
C       ARRAY  :  COEFFICIENT MATRIX IN LEFT-HAND SIDE OF EQ. SS(19)
C       IPVT   :  INTEGER VECTOR OF PIVOT INDICES REQUIRED BY *LINPACK*
C       WK     :  SCRATCH ARRAY REQUIRED BY *LINPACK*
C+---------------------------------------------------------------------+
C
      implicit none
      integer mxcmu,nn,nstr,maz
      real*8 sum,umu0,pi,delm0
      real*8 fbeam
      integer jq,job,k,iq
      INTEGER  IPVT(*)
      REAL*8     CC( MXCMU,* ), CMU(*), GL(0:*),
     $          YLM0(0:*), YLMC( 0:MXCMU,* ),  ZZ(*)
c
       real*8 dARRAY( MXCMU,mxcmu ),dWK(mxcmu),drcond,
     .          dzj(mxcmu)
C
C
      DO 40  IQ = 1, NSTR
C
         DO 10  JQ = 1, NSTR
            dARRAY(IQ,JQ) = (- CC(IQ,JQ))
10       CONTINUE
         dARRAY(IQ,IQ) = 1. + (CMU(IQ) / UMU0) + dARRAY(IQ,IQ)
C
         SUM = 0.
         DO 20  K = MAZ, NSTR-1
            SUM = SUM + GL(K) * YLMC(K,IQ) * YLM0(K)
20       CONTINUE
         dZJ(IQ) = (( 2. - DELM0 ) * FBEAM * SUM / (4.0*PI))
40    CONTINUE
C                  ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITION
C                  ** OF -ARRAY- AND SEE IF IT IS NEARLY SINGULAR
C                  ** (NOTE:  -ARRAY- IS DESTROYED)
      dRCOND = 0.0
      CALL  dGECO( dARRAY, MXCMU, NSTR, IPVT, dRCOND, dWK )
      IF ( 1.0+dRCOND .EQ. 1.0 )  CALL  ERRMSG
     $   ( 'UPBEAM--dGECO SAYS MATRIX NEAR SINGULAR',.FALSE.)
C
C                ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -ARRAY-
C                ** (ASSUMED ALREADY L-U DECOMPOSED) AND R.H. SIDE(S)
C                ** -ZJ-;  RETURN SOLUTION(S) IN -ZJ-
      JOB = 0
      CALL  dGESL( dARRAY, MXCMU, NSTR, IPVT, dZJ, JOB )
C
      DO 50  IQ = 1, NN
         ZZ( IQ+NN )   = sngl(dZJ( IQ ))
         ZZ( NN+1-IQ ) = sngl(dZJ( IQ+NN ))
50    CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  UPISOT(  CC, CMU, IPVT, MXCMU, NN, NSTR, 
     $                    XR0, XR1, XRA,  ZA, ZPLK0, ZPLK1, 
     $                    ZPLKA, dARRAY,dWK,dz0,dz1)
C
C       FINDS THE PARTICULAR SOLUTION OF THERMAL RADIATION OF SS(15)
C
C     ROUTINES CALLED:  SGECO, SGESL
C
C   I N P U T     V A R I A B L E S:
C
C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5)
C       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       OPRIM  :  DELTA-M SCALED SINGLE SCATTERING ALBEDO
C       XR0    :  EXPANSION OF INTERNAL SOURCE FUNCTION
C       XR1    :  EXPANSION OF INTERNAL SOURCE FUNCTION EQS. SS(14-16)
C       XRA    :  EXPANSION OF INTERNAL SOURCE FUNCTION EQS. SS(14-16)
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C    O U T P U T    V A R I A B L E S:
C
C       Z0     :  SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)
C       Z1     :  SOLUTION VECTORS Z-SUB-ONE  OF EQ. SS(16)
C       ZA     :  ALFA COEFFICIENT IN EQ.
C       ZPLK0, :  PERMANENT STORAGE FOR -Z0,Z1,ZA-, BUT RE-ORDERED
C        ZPLK1,
C        ZPLKA
C
C   I N T E R N A L    V A R I A B L E S:
C
C       ARRAY  :  COEFFICIENT MATRIX IN LEFT-HAND SIDE OF EQ. SS(16)
C       IPVT   :  INTEGER VECTOR OF PIVOT INDICES REQUIRED BY *LINPACK*
C       WK     :  SCRATCH ARRAY REQUIRED BY *LINPACK*
C+---------------------------------------------------------------------+
C
      implicit none
      integer nn,nstr,mxcmu,jq,iq
      INTEGER IPVT(*)
      REAL*8     CC( MXCMU,* ), CMU(*), 
     $        XR0(*), XR1(*), XRA, ZA, 
     $        ZPLK0(*), ZPLK1(*), ZPLKA
      real*8 dARRAY( MXCMU,mxcmu ),dWK(mxcmu),drcond,
     .          dz0(mxcmu),dz1(mxcmu)
C
C
c       write(6,'(a8,8(1pe9.1))') ' xr0/1  ',
c    .          (xr0(i),i=1,nstr),(xr1(i),i=1,nstr)
      DO 20 IQ = 1, NSTR
C
         DO 10 JQ = 1, NSTR
            dARRAY(IQ,JQ) = (- CC(IQ,JQ))
10       CONTINUE
         dARRAY(IQ,IQ) = 1.0 + (XRA*CMU(IQ)) + dARRAY(IQ,IQ)
C
         ZA     = XRA
         dZ1(IQ) = (XR1(IQ))
20    CONTINUE
c       if (ien.eq.17)then
c       write(47,*) 'array avant DGECO'
c       write(47,47)darray(1,1),darray(1,2),darray(1,3),darray(1,4)
c       write(47,47)darray(2,1),darray(2,2),darray(2,3),darray(2,4)
c       write(47,47)darray(3,1),darray(3,2),darray(3,3),darray(3,4)
c       write(47,47)darray(4,1),darray(4,2),darray(4,3),darray(4,4)
c       write(47,*)
c       endif
C                       ** SOLVE LINEAR EQUATIONS: SAME AS IN *UPBEAM*,
C                       ** EXCEPT -ZJ- REPLACED BY -Z0-
      dRCOND = 0.0
      CALL  dGECO( dARRAY, MXCMU, NSTR, IPVT, dRCOND, dWK )
      IF ( 1.0+dRCOND .EQ. 1.0 )  CALL  ERRMSG
     $   ( 'UPISOT--dGECO SAYS MATRIX NEAR SINGULAR',.FALSE.)
C
      CALL  dGESL( dARRAY, MXCMU, NSTR, IPVT, dZ1, 0 )
C
      DO 30 IQ = 1, NSTR
         dZ0(IQ) = (XR0(IQ)) + (CMU(IQ)) * dZ1(IQ)
30    CONTINUE
        
      CALL  dGESL( dARRAY, MXCMU, NSTR, IPVT, dZ0, 0 )
C
      DO 40  IQ = 1, NN
         ZPLKA            = ZA
         ZPLK0( IQ+NN )   = sngl(dZ0( IQ ))
         ZPLK1( IQ+NN )   = sngl(dZ1( IQ ))
         ZPLK0( NN+1-IQ ) = sngl(dZ0( IQ+NN ))
         ZPLK1( NN+1-IQ ) = sngl(dZ1( IQ+NN ))
40    CONTINUE
c       if (ien.le.18 .and. ien.ge.16)then
c       if (ien.eq.17)then
c       write(47,'(1pd10.2,4(2(1pd9.1),2x))')drcond,
c    .           (xr1(i),dz1(i),i=1,nn)
c       write(47,'(10x,8(2(1pd9.1),2x))') (xr1(i+nn),dz1(i+nn),i=1,nn)
c       write(47,*) 'array apres DGECO'
c       write(47,47)darray(1,1),darray(1,2),darray(1,3),darray(1,4)
c       write(47,47)darray(2,1),darray(2,2),darray(2,3),darray(2,4)
c       write(47,47)darray(3,1),darray(3,2),darray(3,3),darray(3,4)
c       write(47,47)darray(4,1),darray(4,2),darray(4,3),darray(4,4)
c       write(47,48)
c       endif
c47     format (t10,4(1pd12.4))
c48     format (70('-'))
c
*       write(6,'(1x,a7,(t12,1p8e12.3))')'zplk0/1',
*     .         (zplk0(i),i=1,nstr),(zplk1(i),i=1,nstr)
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  USRINT( BSRC, CMU, CWT, DELM0, EMU, EXPBEA,
     $                    FBEAM, FISOT, GC, GU, KK, LAMBER, LAYRU, LL,
     $                    LYRCUT, MAZ, MXCMU, MXULV, MXUMU, NCUT,
     $                    NLYR, NN, NSTR, INTSRC, NUMU, NTAU, PI, RMU,
     $                    TAUCPR, TSRC, UMU, UMU0, UTAUPR, WK, ZBEAM,  
     $                    Z0U, Z1U, ZA, ZZ, ZPLK0, ZPLK1, ZPLKA, UUM)
C
C       COMPUTES INTENSITY COMPONENTS AT USER OUTPUT ANGLES
C       FOR AZIMUTHAL EXPANSION TERMS IN EQ. SD(2)
C
C   I N P U T    V A R I A B L E S:
C
C       BSRC   :  EMISSION FROM BOTTOM BOUNDARY
C       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT    :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0
C       EMU    :  SURFACE DIRECTIONAL EMISSIVITY (USER ANGLES)
C       EXPBEA :  TRANSMISSION OF INCIDENT BEAM, EXxP(-TAUCPR/UMU0)
C       GC     :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       GU     :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
C                 (i.e., g IN EQ. SC(1) )
C       KK     :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       LAYRU  :  LAYER NUMBER OF USER LEVEL -UTAU-
C       LL     :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
C                 BY SOLVING SCALED VERSION OF EQ. SC(5);
C                 EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
C       LYRCUT :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       NCUT   :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
C       NN     :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
C       RMU    :  SURFACE BIDIRECTIONAL REFLECTIVITY (USER ANGLES)
C       TAUCPR :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       TSRC   :  EMISSION FROM TOP BOUNDARY
C       UTAUPR :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
C                    COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M
C       Z0U    :  Z-SUB-ZERO IN EQ. SS(16) INTERPOLATED TO USER
C                 ANGLES FROM AN EQUATION DERIVED FROM SS(16)
C       Z1U    :  Z-SUB-ONE IN EQ. SS(16) INTERPOLATED TO USER
C                 ANGLES FROM AN EQUATION DERIVED FROM SS(16)
C       ZA     :  ALFA IN EQ.
C       ZZ     :  BEAM SOURCE VECTORS IN EQ. SS(19)
C       ZPLK0  :  INTERNAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
C       ZPLK1  :  INTERNAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)
C       ZPLKA  :  INTERNAL SOURCE VECTORS -ZA , ALFA IN EQ.
C       ZBEAM  :  INCIDENT-BEAM SOURCE VECTORS
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T    V A R I A B L E S:
C
C       UUM  :  AZIMUTHAL COMPONENTS OF THE INTENSITY IN EQ. STWJ(5)
C
C   I N T E R N A L    V A R I A B L E S:
C
C       BNDDIR :  DIRECT INTENSITY DOWN AT THE BOTTOM BOUNDARY
C       BNDDFU :  DIFFUSE INTENSITY DOWN AT THE BOTTOM BOUNDARY
C       BNDINT :  INTENSITY ATTENUATED AT BOTH BOUNDARIES, STWJ(25-6)
C       DTAU   :  OPTICAL DEPTH OF A COMPUTATIONAL LAYER
C       LYREND :  END LAYER OF INTEGRATION
C       LYRSTR :  START LAYER OF INTEGRATION
C       PALINT :  INTENSITY COMPONENT FROM PARALLEL BEAM
C       PLKINT :  INTENSITY COMPONENT FROM INTERNAL SOURCE
C       WK     :  SCRATCH VECTOR FOR SAVING 'EXP' EVALUATIONS
C       ALL THE EXPONENTIAL FACTORS ( EXP1, EXPN,... etc.)
C       COME FROM THE SUBSTITUTION OF CONSTANTS OF INTEGRATION IN
C       EQ. SC(12) INTO EQS. S1(8-9).  THEY ALL HAVE NEGATIVE
C       ARGUMENTS SO THERE SHOULD NEVER BE OVERFLOW PROBLEMS.
C+---------------------------------------------------------------------+
C
      implicit none
      integer mxumu,mxcmu,nlyr,nn,nstr
      real*8 bsrc,delm0,fisot,exxp,pi,tsrc,umu0
      real*8 fbeam
      integer maz,mxulv,ncut,numu,ntau
      LOGICAL  LAMBER, LYRCUT, INTSRC, NEGUMU
      real*8 plkint,sgn,palint,exp1,exp2,expn
      real*8 dtau2,dtau,dtau1,den,denom,bndint,bnddir
      real*8 bnddfu,dfuint,exp0,alfa
      integer lyu,lu,lc,lyrend,lyrstr,iu,iq,jq
      INTEGER  LAYRU(*)
      REAL*8  CMU(*), CWT(*), EMU(*), EXPBEA(0:*), GC( MXCMU,MXCMU,* ),
     $         GU( MXUMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ),
     $         RMU( MXUMU,0:* ), TAUCPR( 0:* ), UUM( MXUMU,MXULV,0:* ),
     $         UMU(*), UTAUPR(*), WK(*), Z0U( MXUMU,* ), Z1U( MXUMU,* ),
     $         ZA( * ), ZBEAM( MXUMU,* ), ZZ( MXCMU,* ),  
     $         ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* ), ZPLKA( * ) 
C
C GG : initialization of exp1 and exp2
        exp1=0
        exp2=0
      CALL  ZEROIT( UUM, MXUMU*MXULV*(MXCMU+1) )
C
C                          ** INCORPORATE CONSTANTS OF INTEGRATION INTO
C                          ** INTERPOLATED EIGENVECTORS
      DO 10  LC = 1, NCUT
         DO  10  IQ = 1, NSTR
            DO 10  IU = 1, NUMU
               GU(IU,IQ,LC) = GU(IU,IQ,LC) * LL(IQ,LC)
10    CONTINUE
C                           ** LOOP OVER LEVELS AT WHICH INTENSITIES
C                           ** ARE DESIRED ('USER OUTPUT LEVELS')
      DO 200  LU = 1, NTAU
C
         EXP0 = EXxP( - UTAUPR(LU) / UMU0 )
         LYU = LAYRU(LU)
C                              ** LOOP OVER POLAR ANGLES AT WHICH
C                              ** INTENSITIES ARE DESIRED
         DO 100  IU = 1, NUMU
            IF ( LYRCUT .AND. LYU.GT.NCUT )  GO TO 100
            NEGUMU = UMU(IU).LT.0.0
            IF( NEGUMU )  THEN
               LYRSTR = 1
               LYREND = LYU - 1
               SGN = - 1.0
            ELSE
               LYRSTR = LYU + 1
               LYREND = NCUT
               SGN = 1.0
            END IF
C                          ** FOR DOWNWARD INTENSITY, INTEGRATE FROM TOP
C                          ** TO 'LYU-1' IN EQ. S1(8); FOR UPWARD,
C                          ** INTEGRATE FROM BOTTOM TO 'LYU+1' IN S1(9)
            PALINT = 0.0
            PLKINT = 0.0
            DO 30  LC = LYRSTR, LYREND
C
               DTAU = TAUCPR(LC) - TAUCPR(LC-1)
               EXP1 =  EXxP( (UTAUPR(LU) - TAUCPR(LC-1)) / UMU(IU) )
               EXP2 =  EXxP( (UTAUPR(LU) - TAUCPR( LC )) / UMU(IU) )
               DEN  =  SGN*1./(ZA(LC)*UMU(IU)+1.)
C
               IF ( INTSRC .AND. MAZ.EQ.0 )
     $           PLKINT = PLKINT +
     $    (Z0U(IU,LC)*DEN*(EXxP(-ZA(LC)*TAUCPR(LC-1))*EXP1
     $         -EXxP(-ZA(LC)*TAUCPR(LC))*EXP2)
     $         +Z1U(IU,LC)*DEN*(
     $         (TAUCPR(LC-1)+SGN*DEN*UMU(IU))*EXxP(-ZA(LC)*
     $          TAUCPR(LC-1))*EXP1
     $       -(TAUCPR(LC)+SGN*DEN*UMU(IU))*EXxP(-ZA(LC)*TAUCPR(LC))*
     $          EXP2))
C
               IF ( FBEAM.GT.0.0 )  THEN
                  DENOM = 1.0 + UMU(IU) / UMU0
                  IF ( ABS(DENOM).LT.0.0001 ) THEN
C                                                   ** L'HOSPITAL LIMIT
                     EXPN = ( DTAU / UMU0 ) * EXP0
                  ELSE
                     EXPN = ( EXP1 * EXPBEA(LC-1) - EXP2 * EXPBEA(LC) )
     $                      * SGN / DENOM
                  END IF
                  PALINT = PALINT + ZBEAM(IU,LC) * EXPN
               ENDIF
C                                                   ** -KK- IS NEGATIVE
               DO 20  IQ = 1, NN
                  WK(IQ) = EXxP( KK(IQ,LC) * DTAU )
                  DENOM = 1.0 + UMU(IU) * KK(IQ,LC)
                  IF ( ABS(DENOM).LT.0.0001 ) THEN
C                                                   ** L'HOSPITAL LIMIT
                     EXPN = DTAU / UMU(IU) * EXP2
                  ELSE
                     EXPN = SGN * ( EXP1 * WK(IQ) - EXP2 ) / DENOM
                  END IF
                  PALINT = PALINT + GU(IU,IQ,LC) * EXPN
20             CONTINUE
C                                                   ** -KK- IS POSITIVE
               DO 21  IQ = NN+1, NSTR
                  DENOM = 1.0 + UMU(IU) * KK(IQ,LC)
                  IF ( ABS(DENOM).LT.0.0001 ) THEN
C                                                   ** L'HOSPITAL LIMIT
                     EXPN = - DTAU / UMU(IU) * EXP1
                  ELSE
                     EXPN = SGN *( EXP1 - EXP2 * WK(NSTR+1-IQ) ) / DENOM
                  END IF
                  PALINT = PALINT + GU(IU,IQ,LC) * EXPN
21             CONTINUE
C
30          CONTINUE
C                           ** CALCULATE CONTRIBUTION FROM USER
C                           ** OUTPUT LEVEL TO NEXT COMPUTATIONAL LEVEL
C
            DTAU1 = UTAUPR(LU) - TAUCPR(LYU-1)
            DTAU2 = UTAUPR(LU) - TAUCPR(LYU)
            IF( ABS(DTAU1).LT.1.E-6 .AND. NEGUMU )  GO TO 50
            IF( ABS(DTAU2).LT.1.E-6 .AND. (.NOT.NEGUMU) )  GO TO 50
            IF( NEGUMU ) EXP1 = EXxP( DTAU1 / UMU(IU) )
            IF( .NOT.NEGUMU ) EXP2 = EXxP( DTAU2 / UMU(IU) )
C
            IF ( FBEAM.GT.0.0 )  THEN
               DENOM = 1.0 + UMU(IU) / UMU0
               IF ( ABS(DENOM).LT.0.0001 ) THEN
                  EXPN =  ( DTAU1 / UMU0 ) * EXP0
               ELSE IF ( NEGUMU ) THEN
                  EXPN = ( EXP0 - EXPBEA(LYU-1) * EXP1 ) / DENOM
               ELSE
                  EXPN = ( EXP0 - EXPBEA(LYU) * EXP2 ) / DENOM
               END IF
               PALINT = PALINT + ZBEAM(IU,LYU) * EXPN
            ENDIF
C                                                   ** -KK- IS NEGATIVE
            DTAU = TAUCPR(LYU) - TAUCPR(LYU-1)
            DO 40  IQ = 1, NN
               DENOM = 1.0 + UMU(IU) * KK(IQ,LYU)
               IF ( ABS(DENOM).LT.0.0001 ) THEN
                  EXPN = - DTAU2 / UMU(IU) * EXP2
               ELSE IF ( NEGUMU ) THEN
                  EXPN = ( EXxP( - KK(IQ,LYU) * DTAU2 ) -
     $                     EXxP( KK(IQ,LYU) * DTAU ) * EXP1 ) / DENOM
               ELSE
                  EXPN = ( EXxP( - KK(IQ,LYU) * DTAU2 ) - EXP2 ) / DENOM
               END IF
               PALINT = PALINT + GU(IU,IQ,LYU) * EXPN
40          CONTINUE
C                                                   ** -KK- IS POSITIVE
            DO 41  IQ = NN+1, NSTR
               DENOM = 1.0 + UMU(IU) * KK(IQ,LYU)
               IF ( ABS(DENOM).LT.0.0001 ) THEN
                  EXPN = - DTAU1 / UMU(IU) * EXP1
               ELSE IF ( NEGUMU ) THEN
                  EXPN = ( EXxP(- KK(IQ,LYU) * DTAU1 ) - EXP1 ) / DENOM
               ELSE
                  EXPN = ( EXxP( - KK(IQ,LYU) * DTAU1 ) -
     $                     EXxP( - KK(IQ,LYU) * DTAU ) * EXP2 ) / DENOM
               END IF
               PALINT = PALINT + GU(IU,IQ,LYU) * EXPN
41          CONTINUE
C
            IF ( INTSRC .AND. MAZ.EQ.0 )  THEN
              IF ( NEGUMU ) THEN
                 EXPN = EXP1
                 ALFA = ZA(LYU)
                 DEN = (-1./(ALFA*UMU(IU)+1.))
                 PLKINT = PLKINT + 
     $                    Z0U(IU,LYU) * DEN * (
     $                   -EXxP(-ALFA*UTAUPR(LU))+ 
     $                   exp(-alfa*taucpr(lyu-1))*EXPN )
     $                   +  Z1U(IU,LYU)*DEN*(
     $                  -(UTAUPR(LU)-UMU(IU)*DEN)*EXxP(-ALFA*UTAUPR(LU))
     $                   +(TAUCPR(LYU-1)-UMU(IU)*DEN) *
     $                   exp(-alfa*taucpr(lyu-1))*EXPN)
              ELSE
                 EXPN = EXP2
                 ALFA = ZA(LYU)
                 DEN = (1./(ALFA*UMU(IU)+1.))
                 PLKINT = PLKINT + 
     $                    Z0U(IU,LYU) * DEN *
     $                   (EXxP(-ALFA*UTAUPR(LU))-EXxP(-ALFA*TAUCPR(LYU))
     $                                           *EXPN)
     $                 +  Z1U(IU,LYU)*DEN*(
     $                  (UTAUPR(LU) +UMU(IU)*DEN)*EXxP(-ALFA*UTAUPR(LU))
     $                -(TAUCPR(LYU)+UMU(IU)*DEN)*EXxP(-ALFA*TAUCPR(LYU))
     $                                             *EXPN)
              END IF
            END IF
C                            ** CALCULATE INTENSITY COMPONENTS
C                            ** ATTENUATED AT BOTH BOUNDARIES.
C                            ** NOTE:: NO AZIMUTHAL INTENSITY
C                            ** COMPONENT FOR ISOTROPIC SURFACE
50          BNDINT = 0.0
            IF ( NEGUMU .AND. MAZ.EQ.0 ) THEN
              BNDINT = ( FISOT + TSRC ) * EXxP( UTAUPR(LU) / UMU(IU) )
            ELSE IF ( .NOT.NEGUMU ) THEN
              IF ( LYRCUT .OR. (LAMBER .AND. MAZ.GT.0) )  GO TO 90
              DO 60  JQ = NN+1, NSTR
               WK(JQ) = EXxP(-KK(JQ,NLYR)*(TAUCPR(NLYR)-TAUCPR(NLYR-1)))
60            CONTINUE
              BNDDFU = 0.0
              DO 80  IQ = NN, 1, -1
                 DFUINT = 0.0
                 DO 70  JQ = 1, NN
                    DFUINT = DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)
70               CONTINUE
                 DO 71  JQ = NN+1, NSTR
                    DFUINT = DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)
     $                                * WK(JQ)
71               CONTINUE
                 IF ( FBEAM.GT.0.0 )
     $                DFUINT = DFUINT + ZZ(IQ,NLYR) * EXPBEA(NLYR)
                 DFUINT = DFUINT + DELM0 * ( 
     $              EXxP(-ZPLKA(NLYR)*TAUCPR(NLYR))*
     $              (ZPLK0(IQ,NLYR)+ZPLK1(IQ,NLYR)*TAUCPR(NLYR)))
                 BNDDFU = BNDDFU + ( 1. + DELM0 ) * RMU(IU,NN+1-IQ)
     $                           * CMU(NN+1-IQ) * CWT(NN+1-IQ) * DFUINT
80            CONTINUE
C
              BNDDIR = 0.0
              IF (FBEAM.GT.0.0) BNDDIR = UMU0*FBEAM/PI * RMU(IU,0)
     $                                   * EXPBEA(NLYR)
              BNDINT = ( BNDDFU + BNDDIR + DELM0 * EMU(IU) * BSRC )
     $                 * EXxP( (UTAUPR(LU)-TAUCPR(NLYR)) / UMU(IU) )
            END IF
C
90          UUM( IU, LU, MAZ ) = PALINT + PLKINT + BNDINT
C
100      CONTINUE
200   CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  ZEROAL( AMB, APB, ARRAY, CC, CMU, CWT, EVAL, EVECC,
     $                    GC, GU, HLPR, KK, LL, PSI0, PSI1, WK, XR0, 
     $                    XR1, XRA,XR0U, XR1U, XRAU, YLM0, YLMC, YLMU, 
     $                    Z, Z0, Z1, ZA, ZJ, ZZ, ZPLK0, ZPLK1, ZPLKA, 
     $                    Z0U, Z1U, ZBEAM, MI, MXCMU, MXCLY, 
     $                    NNLYRI, MXUMU )
C
C            ZERO ARRAYS
C
      implicit none
      integer mxcmu,mxcly,mi,mxumu,nnlyri
      REAL*8  AMB(MI,*), APB(MI,*), ARRAY(MXCMU,*), CC(MXCMU,*),
     $      CMU(*), CWT(*), EVAL(*), EVECC(MXCMU,*), GC(MXCMU,MXCMU,*),
     $      GU(MXUMU,MXCMU,*), HLPR(0:*), KK(MXCMU,*), LL(MXCMU,*),
     $      PSI0(*), PSI1(*), WK(*), XR0(MXCMU,*), XR1(MXCMU,*), 
     $      XRA(*), XR0U(MXUMU,*),XR1U(MXUMU,*),XRAU(*), YLM0(0:*),
     $      YLMC(0:MXCMU,*), YLMU(0:MXCMU,*), Z(*), Z0(*), Z1(*), ZA(*),
     $      Z0U(MXUMU,*), Z1U(MXUMU,*),ZJ(*), ZZ(MXCMU,*),
     $      ZPLK0(MXCMU,*), ZPLK1(MXCMU,*), ZPLKA(*),
     $      ZBEAM(MXUMU,*)  
C
C
      CALL  ZEROIT( XR0, MXCMU*MXCLY )
      CALL  ZEROIT( XR1, MXCMU*MXCLY )
      CALL  ZEROIT( XRA, MXCLY )
      CALL  ZEROIT( XR0U, MXUMU*MXCLY )
      CALL  ZEROIT( XR1U, MXUMU*MXCLY )
      CALL  ZEROIT( XRAU, MXCLY )
      CALL  ZEROIT( CMU, MXCMU )
      CALL  ZEROIT( CWT, MXCMU )
      CALL  ZEROIT( PSI0, MXCMU )
      CALL  ZEROIT( PSI1, MXCMU )
      CALL  ZEROIT( EVAL,MI )
      CALL  ZEROIT( WK,  MXCMU )
      CALL  ZEROIT( Z0,  MXCMU )
      CALL  ZEROIT( Z1,  MXCMU )
      CALL  ZEROIT( ZA,  MXCLY )
      CALL  ZEROIT( ZJ,  MXCMU )
      CALL  ZEROIT( HLPR, MXCMU+1 )
      CALL  ZEROIT( YLM0, MXCMU+1 )
      CALL  ZEROIT( ARRAY, MXCMU**2 )
      CALL  ZEROIT( CC,    MXCMU**2 )
      CALL  ZEROIT( EVECC, MXCMU**2 )
      CALL  ZEROIT( YLMC, (MXCMU+1)*MXCMU )
      CALL  ZEROIT( YLMU, (MXCMU+1)*MXUMU )
      CALL  ZEROIT( AMB, MI**2 )
      CALL  ZEROIT( APB, MI**2 )
      CALL  ZEROIT( KK,     MXCMU*MXCLY )
      CALL  ZEROIT( LL,     MXCMU*MXCLY )
      CALL  ZEROIT( ZZ,     MXCMU*MXCLY )
      CALL  ZEROIT( ZPLK0,  MXCMU*MXCLY )
      CALL  ZEROIT( ZPLK1,  MXCMU*MXCLY )
      CALL  ZEROIT( ZPLKA,  MXCLY )
      CALL  ZEROIT( Z0U,   MXUMU*MXCLY )
      CALL  ZEROIT( Z1U,   MXUMU*MXCLY )
      CALL  ZEROIT( ZBEAM, MXUMU*MXCLY )
      CALL  ZEROIT( GC, MXCMU**2*MXCLY )
      CALL  ZEROIT( GU, MXUMU*MXCMU*MXCLY )
      CALL  ZEROIT( Z, NNLYRI )
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  dZEROIT( A, LENGTH )
C
C     ZEROS A double precision ARRAY -A- HAVING -LENGTH- ELEMENTS
C
      implicit none
      real*8 A(*)
      integer length,l
C
      DO 10  L = 1, LENGTH
         A( L ) = 0.0
10    CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  ZEROIT( A, LENGTH )
C
C         ZEROS A REAL ARRAY -A- HAVING -LENGTH- ELEMENTS
C
      implicit none
      integer length,l
      REAL*8  A(*)
C
      DO 10  L = 1, LENGTH
         A( L ) = 0.0
10    CONTINUE
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE  ERRMSG( MESSAG, FATAL)
C
C        PRINT OUT A WARNING OR ERROR MESSAGE;  ABORT IF ERROR
C
      implicit none
      LOGICAL       FATAL, ONCE
      CHARACTER*(*) MESSAG
      INTEGER       MAXMSG, NUMMSG
      SAVE          MAXMSG, NUMMSG, ONCE
      DATA NUMMSG / 0 /,  MAXMSG / 100 /,  ONCE / .FALSE. /
C
C
      IF ( FATAL )  THEN
c         WRITE ( 6, '(/,2A)' )  ' ******* ERROR >>>>>>  ', MESSAG
         WRITE (*,*)  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF
C
      NUMMSG = NUMMSG + 1
      IF ( NUMMSG.GT.MAXMSG )  THEN
         IF ( .NOT.ONCE )  WRITE (*,99 )
c         IF ( .NOT.ONCE )  WRITE ( 6,99 )
         ONCE = .TRUE.
      ELSE
c   WRITE ( 6, '(/,2A)' )  ' ******* WARNING >>>>>>  ', MESSAG
         WRITE (*,*)  ' ******* WARNING >>>>>>  ', MESSAG
      ENDIF
C
      RETURN
C
   99 FORMAT( ///,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     $   'THEY WILL NO LONGER BE PRINTED  <<<<<<<', /// )
      END
      LOGICAL FUNCTION  WRTBAD ( VARNAM )
      implicit none
C
C          WRITE NAMES OF ERRONEOUS VARIABLES AND RETURN 'TRUE'
C
C      INPUT :   VARNAM = NAME OF ERRONEOUS VARIABLE TO BE WRITTEN
C                         ( CHARACTER, ANY LENGTH )
C ----------------------------------------------------------------------
      CHARACTER*(*)  VARNAM
      INTEGER        MAXMSG, NUMMSG
      SAVE  NUMMSG, MAXMSG
      DATA  NUMMSG / 0 /,  MAXMSG / 50 /
C
C
      WRTBAD = .TRUE.
      NUMMSG = NUMMSG + 1
c      WRITE ( 6, '(3A)' )  ' ****  INPUT VARIABLE  ', VARNAM,
      WRITE (*, '(3A)' )  ' ****  INPUT VARIABLE  ', VARNAM,
     $                     '  IN ERROR  ****'
      IF ( NUMMSG.EQ.MAXMSG )
     $   CALL  ERRMSG ( 'TOO MANY INPUT ERRORS.  ABORTING...$', .TRUE.
     $                  )
      RETURN
      END
      LOGICAL FUNCTION  WRTDIM ( DIMNAM, MINVAL )
C
C          WRITE NAME OF TOO-SMALL SYMBOLIC DIMENSION AND
C          THE VALUE IT SHOULD BE INCREASED TO;  RETURN 'TRUE'
C
C      INPUT :  DIMNAM = NAME OF SYMBOLIC DIMENSION WHICH IS TOO SMALL
C                        ( CHARACTER, ANY LENGTH )
C               MINVAL = VALUE TO WHICH THAT DIMENSION SHOULD BE
C                        INCREASED (AT LEAST)
C ----------------------------------------------------------------------
      CHARACTER*(*)  DIMNAM
      INTEGER        MINVAL
C
C
c      WRITE ( 6, '(3A,I7)' ) ' ****  SYMBOLIC DIMENSION  ', DIMNAM,
      WRITE ( *, '(3A,I7)' ) ' ****  SYMBOLIC DIMENSION  ', DIMNAM,
     $                     '  SHOULD BE INCREASED TO AT LEAST ', MINVAL
      WRTDIM = .TRUE.
      RETURN
      END
c
c----------------------------------------------------------------------
c
        function exxp(a)
        implicit none
        real*8 exxp,a
        if(a.lt.-69.) then
                exxp=0.
        else
                exxp=exp(a)
        end if
        return
        end
