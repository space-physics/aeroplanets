C From:	PACF2::HEDIN        15-FEB-1990 14:21:14.82
C To:	NSSDCA::BILITZA
C CC:
C Subj:	VENUS MODEL

C Dieter,   Here is the Venus subroutine as well as a test program
C similar to that for MSIS.  I haven't run it in some time.  I am sure it
C has not been worked over as much as MSIS to make it compatible with
C other computers. But, Give it a try.
C Al'
C -----------------------------------------------------------------------
C VTS3   main Venus subroutine
C -----------------------------------------------------------------------
      SUBROUTINE VTS3(ALT,YRD,XLAT,TLOC,F107A,F107,MAS ,D,T)
C       VENUS ATMOSPHERE MODEL (03/05/81)
C         USES LOCAL TIME AND LATITUDE TO ESTIMATE SOLAR ZENITH
C         ANGLE FOR SYMMETRICAL TERMS

C    INPUT VARIABLES
C      ALT - ALTITUDE (KM)
C      XLAT - LATITUDE (DEG)
C      XLOC - LOCAL HOUR ANGLE (HRS)
C       MAS - MASS REQUIRED, 48 FOR ALL
C      YRD            NOT CURRENTLY USED

C    OUTPUT VARIABLES
C      D(1) - TOTAL MASS DENSITY (GM/CM3)
C        2  - CO2 NUMBER DENSITY (CM-3)
C        3  - O
C        4  - CO
C        5  - HE
C        6  - N
C        7  - N2
C      T(1) - EXOSPHERIC TEMPERATURE
C        2  - TEMPERATURE AT ALTITUDE
      INTEGER, intent(in) :: YRD
      real, intent(in) :: alt, xlat, tloc, f107a, f107
      integer, intent(in) :: mas
      real, intent(out) :: D(7),T(2)

      real :: MT(7)

      COMMON/PARM/PT(50),PD(50,6)
      COMMON/LOWER/PTM(8),PDM(8,6)
      COMMON/VTS1C/TLB,S,DB44,DB16,DB29,DB04,DB14,DB28
      COMMON/PARMB/GSURF,RE
      COMMON/VFIT/ZAA,FRR
      DATA IFL/0/,MT/48,44,16,29,4,14,28/
      REAL truc

      MASS = ABS(MAS)
      TB=0
      IF(MAS.EQ.-48) TB=1
      IF(IFL.EQ.1) GOTO 10
      CALL PARMSV
      IFL=1
10    CONTINUE
C       TEMPERATURE STRUCTURE
      TINFV=GLOBV(YRD,XLAT,TLOC,F107A,F107,PT(1))
      TINF=PTM(1)*(1.+TINFV)*PT(1)
      T(1)=TINF
      ZA=PTM(5)*PT(31)
      ZAA=ZA
      T0=PTM(3)*PT(47)
      S=PTM(4)*PT(46)
      SP=PT(49)
C      SET TA BETWEEN T0 AND TINF
      FR=PT(16)*(1.+GLOBV(YRD,XLAT,TLOC,F107A,F107,PT(16)))
      FRR=FR
      FR=1./(1.+EXP(-4.*(FR-.5)))
      TA=TINF+FR*(T0-TINF)
      ZM=PTM(7)*PT(48)
      AM=PTM(8)*PD(31,1)*(1.+GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(31,1)))
      IF(MASS.EQ.0) GO TO 50
      ENTRY VTSB(D,T)
      XMR=PDM(2,2)*PD(16,2)*(1.+GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(16,2))
     &)
      IF(XMR.LT.1.E-3) XMR=1.E-3
      H1=PDM(7,2)*PD(47,2)
      ZHO=PDM(6,2)*PD(46,2)
      YMR=XMR*CCORVENUS(PDM(6,1),PDM(5,2),H1,ZHO)
C       CALCULATE MEAN MASS
      XMM=(44.+44.*YMR+PDM(2,6)*28.)/(1.+2.*YMR+PDM(2,6))
      PDM(5,1)=XMM
      ZHM=PDM(4,1)/PD(31,2)
13    CONTINUE
      DO 15 J=1,7
      IF(MASS.EQ.MT(J)) GO TO 20
15    CONTINUE
      GO TO 90
C       CO2 DENSITY
20    CONTINUE
      DB44=PDM(1,1)*EXP(GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(1,1)))*PD(1,1)
      D(2)=DENSS(ALT,DB44,TINF,TLB,44.,0.,T(2),PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      ZH44=PDM(3,1)*PD(16,1)
C        GET MIXING DENSITY AT ZLB
      XMD=44.-XMM
      B44=DENSS(ZH44,DB44,TINF,TLB,XMD,-1.,TZ,PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      DM44=DENSS(ALT,B44,TINF,TLB,XMM,0.,TZ,PTM(6),S,T0,ZA,TA,ZM,AM,SP)
      D(2)=DNETVENUS(D(2),DM44,ZHM,XMM,44.)
22    CONTINUE
C       O DENSITY
25    DB16=PDM(1,2)*EXP(GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(1,2)))*PD(1,2)
      D(3)=DENSS(ALT,DB16,TINF,TLB,16.,0.,T(2),PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      DM16=DM44*XMR
      D(3)=DNETVENUS(D(3),DM16,ZHM,XMM,16.)
      D(3)=D(3)*CCORVENUS(ALT,PDM(5,2),H1,ZHO)
      IF(MASS.EQ.44.OR .MASS.EQ.16.OR .MASS.EQ.04) GO TO 27
C        GET O TURBOPAUSE ESTIMATE
      DD16=DENSS(ZH44,DB16,TINF,TLB,16.,0.,TZ,PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      DMZ44=DENSS(ZH44,B44,TINF,TLB,XMM,0.,TZ,PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      CALL TURBO(DD16,DMZ44*XMR,ZH44,ZH16,16.,XMM,TZ)
      PDM(3,2)=ZH16
27    CONTINUE
      GO TO (30,90,90,30,35,40,45),J
C       CO DENSITY
30    DB29=PDM(1,3)*EXP(GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(1,3)))*PD(1,3)
      D(4)=DENSS(ALT,DB29,TINF,TLB,28.,0.,T(2),PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
c      write(*,*),'retour de globv'
      truc=EXP(GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(1,3)))
c      write(*,*),'GLOBV  ',truc,' DB29', DB29
c      write(*,*),' e ',YRD,' ',XLAT,' ',TLOC,' ',F107A,' ',F107,' '

      DM29=DM44*XMR
      D(4)=DNETVENUS(D(4),DM29,ZHM,XMM,28.)
      D(4)=D(4)*CCORVENUS(ALT,PDM(5,2),H1,ZHO)
      IF(TB.EQ.0.) GO TO 32
      DD29=DENSS(ZH44,DB29,TINF,TLB,28.,0.,TZ,PTM(6),S,T0,ZA,TA,
     & ZM,AM,SP)
      CALL TURBO(DD29,DMZ44*XMR,ZH44,PDM(3,3),28.,XMM,TZ)
32    CONTINUE
      IF(MASS.EQ.29) GO TO 40
      IF(MASS.NE.48) GO TO 90
C      HE DENSITY
35    DB04=PDM(1,4)*EXP(GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(1,4)))*PD(1,4)
      D(5)=DENSS(ALT,DB04,TINF,TLB,4.,-.6,T(2),PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      DM04=DM44*PDM(2,4)*PD(16,4)
      D(5)=DNETVENUS(D(5),DM04,ZHM,XMM,4.)
      IF(TB.EQ.0.) GO TO 39
      DD04=DENSS(ZH44,DB04,TINF,TLB,4.,-.6,TZ,PTM(6),S,T0,ZA,TA,
     & ZM,AM,SP)
      CALL TURBO(DD04,DMZ44*PDM(2,4)*PD(16,4),ZH44,PDM(3,4),4.,XMM,TZ)
39    CONTINUE
      IF(MASS.NE.48) GO TO 90
C      N DENSITY
40    DB14=PDM(1,5)*EXP(GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(1,5)))*PD(1,5)
      D(6)=DENSS(ALT,DB14,TINF,TLB,14.,0.0,T(2),PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      ZH14=ZH16
      PDM(3,5)=ZH14
      XMD=14.-XMM
      B14=DENSS(ZH14,DB14,TINF,TLB,XMD,-1.,TZ,PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      DM14=DENSS(ALT,B14,TINF,TLB,XMM,0.,TZ,PTM(6),S,T0,ZA,TA,ZM,AM,SP)
      D(6)=DNETVENUS(D(6),DM14,ZHM,XMM,14.)
      D(6)=D(6)*CCORVENUS(ALT,PDM(5,5),PDM(7,5),PDM(6,5))
42    CONTINUE
      IF(MASS.NE.48) GO TO 90
C       N2 DENSITY
45    DB28=PDM(1,6)*EXP(GLOBV(YRD,XLAT,TLOC,F107A,F107,PD(1,6)))*PD(1,6)
      D(7)=DENSS(ALT,DB28,TINF,TLB,28.,0.,T(2),PTM(6),S,T0,ZA,TA,ZM,AM,
     &SP)
      DM28=DM44*PDM(2,6)
      D(7)=DNETVENUS(D(7),DM28,ZHM,XMM,28.)
      IF(TB.EQ.0.) GO TO 47
      DD28=DENSS(ZH44,DB28,TINF,TLB,28.,0.,TZ,PTM(6),S,T0,ZA,TA,
     & ZM,AM,SP)
      CALL TURBO(DD28,DMZ44*PDM(2,6),ZH44,PDM(3,6),28.,XMM,TZ)
47    CONTINUE
      IF(MASS.EQ.28) GO TO 40
      IF(MASS.NE.48) GO TO 90
C      TOTAL MASS DENSITY
      D(1) = 1.66E-24*(44.*D(2)+16.*D(3)+28.*D(4)+
     & 4.*D(5)+14.*D(6)+28.*D(7))
      GO TO 90
50    DD=DENSS(ALT,1.,TINF,TLB,0.,0.,T(2),PTM(6),S,T0,ZA,TA,ZM,AM,SP)
90    CONTINUE
      DO 91 I=1,7
C      MULTIPLICATION FACTOR FOR ONMS DENSITIES
91    D(I)=D(I)*PDM(8,1)

      END subroutine vts3

C------------------------------------------------------------------
      SUBROUTINE PARMSV
C        SET INITIAL PARAMETERS (04/24/81)
      COMMON/PARM/PT(50),PD(50,6)
      COMMON/LOWER/PTM(8),PDM(8,6)
      COMMON/PARMB/GSURF,RE
      DIMENSION PT1(50),PA1(50),PB1(50),PC1(50),PD1(50),PE1(50),
     & PF1(50),PP(50,6),PM(8,7)
      EQUIVALENCE (PA1(1),PP(1,1)),(PB1(1),PP(1,2)),(PC1(1),PP(1,3)),
     & (PD1(1),PP(1,4)),(PE1(1),PP(1,5)),(PF1(1),PP(1,6))
C         TEMPERATURE
      DATA PT1/
     *  1.14575E00, 5.73938E-01,-1.05623E-01,-1.53896E-01,-7.11596E-03,
     * -1.82894E-01, 4.41052E-03, 1.22197E-01, 3.20351E-04,-9.28368E-03,
     * -2.32161E-05, 0.0        , 0.0        , 1.00000E-03, 6.00000E-04,
     *  7.93416E-01, 1.30399E-01, 8.82217E-02,-4.98274E-01,-2.05990E-02,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  1.12000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  7.23933E-01, 9.74046E-01, 1.00000E00, 7.24953E-02, 0.0        /
C         CO2 DENSITY
      DATA PA1/
     *  7.76049E-01, 2.93750E00,-1.47210E-01,-1.57306E-01,-6.69599E-02,
     * -8.18055E-01,-1.06697E-02, 3.00201E-01, 7.96075E-04, 3.24607E-01,
     *  9.23450E-05, 0.0        , 0.0        ,-2.37584E-03,-1.34238E-04,
     *  1.00000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  1.13277E00,-1.87582E00,-3.38029E-01,-9.31956E-01, 9.04382E-02,
     *  1.67238E00, 7.32745E-03, 8.28310E-01, 1.69341E-03,-6.84008E-01,
     * -1.00458E-04, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 2.01296E-02, 0.0        /
C         O DENSITY
      DATA PB1/
     *  1.07134E00, 7.13595E-01,-3.61877E-02, 2.82843E-01, 4.85210E-03,
     * -1.95412E-01,-1.76002E-03,-3.34167E-01,-9.68110E-04, 3.87223E-01,
     *  3.88084E-05, 0.0        , 0.0        , 2.84044E-03, 1.20219E-03,
     *  3.00713E-02, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  3.90127E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  1.12140E00, 1.30508E00, 0.0        , 0.0        , 9.80888E-04/
C         CO DENSITY & TLB
      DATA PC1/
     *  9.84586E-01, 1.92394E00,-8.16346E-02, 1.85261E-01,-4.62701E-03,
     * -4.08268E-01,-1.68582E-03,-2.05573E-01,-1.15921E-03, 4.93592E-01,
     * -2.59753E-05, 0.0        , 0.0        , 1.39529E-03, 5.53068E-04,
     *  1.00000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  1.00000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        /
C         HE DENSITY
      DATA PD1/
     *  9.61871E-01,-1.42734E00, 5.93447E-01, 9.36320E-02, 1.39517E-01,
     *  8.39837E-01,-3.60608E-03,-3.57368E-01,-1.38972E-03,-1.96184E-02,
     *  8.86656E-05, 0.0        , 0.0        , 2.15513E-03,-7.68169E-04,
     *  3.03416E-02, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  1.00000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        /
C          N DENSITY
      DATA PE1/
     *  1.07218E00, 1.07061E00,-1.92920E-01, 1.72877E-01,-4.19517E-02,
     * -2.37737E-01,-3.55792E-04,-2.46647E-01,-8.06493E-04, 4.72518E-01,
     *  8.04218E-06, 0.0        , 0.0        , 4.85444E-03, 1.24067E-03,
     *  1.00000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  1.00000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        /
C         N2 DENSITY
      DATA PF1/
     *  8.11373E-01, 1.15335E00, 1.26869E-01, 2.25488E-01, 4.20755E-02,
     * -2.21562E-01, 6.93173E-03,-4.49676E-01,-2.56643E-04, 5.91909E-01,
     *  1.22099E-05, 0.0        , 0.0        ,-1.42414E-03,-6.52404E-04,
     *  1.00000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  1.00000E00, 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     *  0.0        , 0.0        , 0.0        , 1.88000E00, 0.0        /
C         LOWER BOUNDARY
      DATA PM /      1.99171E02, 1.87600E02, 1.77000E02, 9.21400E-02,
     *  1.25000E02, 1.50000E02, 1.40000E02, 1.00000E-01,
     *               5.00530E08, 0.0        , 1.30000E02, 4.00000E01,
     *  4.34577E01, 1.00000E02, 0.0        , 1.63000E00,
     *               1.04912E09, 1.00000E-01, 1.12862E02, 0.0        ,
     *  0.0        , 1.20000E02, 5.80000E00, 1.00000E00,
     *               2.95223E08, 0.0        , 1.26638E02, 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        ,
     *               5.64985E06, 2.00000E-05, 1.23157E02, 1.00000E01,
     *  0.0        , 0.0        , 0.0        , 0.0        ,
     *               9.30114E06, 0.0        , 1.12862E02, 0.0        ,
     *  0.0        , 1.24000E02, 1.45000E00, 0.0        ,
     *               1.39749E08, 3.50000E-02, 1.29748E02, 0.0        ,
     *  0.0        , 0.0        , 0.0        , 0.0        /
      GSURF = 887.4
      RE = 6053.
      DO 10 J=1,50
10    PT(J) = PT1(J)
      DO 20 K=1,6
      DO 20 J=1,50
20    PD(J,K) = PP(J,K)
      DO 30 J=1,8
30    PTM(J) = PM(J,1)
      DO 40 K=1,6
      DO 40 J=1,8
40    PDM(J,K) = PM(J,K+1)

      END SUBROUTINE PARMSV

C--------------------------------------------------------------------
      FUNCTION DENSS(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,T0,ZA,TA,ZM,AM,
     & SP)
C        CALCULATE DENSITY (12/12/80)
C         TLB CALCULATED FROM TA
      COMMON/PARMB/GSURF,RE
      COMMON/FIT/TAF
      DATA RGAS/831.4/
      ZETA(ZZ)=(ZZ-ZLB)*(RE+ZLB)/(RE+ZZ)
      DENSS=DLB
        TAF=(TA-T0)/(TINF-T0)
        IF(TAF.GT.0.AND.TAF.LT.1.)  GO TO 7
C         WRITE(6,110) TA,TINF,TLB,S2,ALT,ZA,T0
          IF(TAF.LE.0.)TA=T0+0.0001*(TINF-T0)
          IF(TAF.GE.1.)TA=TINF-.0001*(TINF-T0)
7       CONTINUE
110     FORMAT(' DENSS ',1P7E11.2)
5       Z=AMAX1(ALT,ZA)
        ZG2=ZETA(Z)
        ZGA=ZETA(ZA)
      ZG=ZETA(ALT)
        TLB=TINF-(TINF-TA)*EXP(S2*ZGA)
        T2=TINF-(TINF-TLB)*EXP(-S2*ZG2)
        TT=T2
        IF(ALT.GE.ZA) GO TO 10
          S1=-S2*(TINF-TA)/(TA-T0)
          ZG1=ZG-ZGA
C           CALCULATE TEMPERATURE BELOW ZA
          IF(S1*ZG1.LT. 10.) GO TO 8
            T1=T0
            GO TO 9
8         T1=T0-(T0-TA)*EXP(-S1*ZG1)
9         TT=T1
10      CONTINUE
        TZ=TT
        TZL=TLB
        IF(AM.EQ.0.) GO TO 15
          ZGM=ZETA(ZM)
          EXZ=EXP(-(SP*(ZG -ZGM)))
          DIT=4.*AM*EXZ/(1.+EXZ)**2/TINF
          TZ=TT/(1.+TT*DIT)
          EXM=EXP(SP*ZGM)
          DITL=4.*AM*EXM/(1.+EXM)**2/TINF
          TZL=TLB/(1.+TLB*DITL)
15        IF(XM.EQ.0.) GO TO 350
C          CALCULATE DENSITY ABOVE ZA
          GLB=GSURF/(1.+ZLB/RE)**2
          GAMMA2=XM*GLB/(S2*RGAS*TINF)
          DENS2=DLB*(TLB/T2)**(GAMMA2)*EXP(-S2*GAMMA2*ZG2)
          DENSS=DENS2
          IF(ALT.GE.ZA) GO TO 300
C            CALCULATE DENSITY BELOW ZA
            GAMMA1=XM*GLB/(S1*RGAS*T0)
            DENS1=DENS2*(TA/TT)**(GAMMA1)*EXP(-S1*GAMMA1*ZG1)
            DENSS=DENS1
300      CONTINUE
         IF(AM.EQ.0.) GO TO 340
          GAMMAP=XM*GLB/(SP*RGAS*TINF)
           DENSM=EXP(GAMMAP*4.*AM*(EXZ/(1.+EXZ)-EXM/(1.+EXM)))
           DENSS=DENSM*DENSS
340      CONTINUE
         DENSS=DENSS*(TZL/TZ)**(1.+ALPHA)
350    CONTINUE

      END FUNCTION DENSS


      FUNCTION DNETVENUS(DD,DM,ZHM,XMM,XM)
C       8/20/80;7/25/83
C       TURBOPAUSE CORRECTION
      A=ZHM/(XMM-XM)
      ADM=0.
      IF(DM.GT.0.) ADM=ALOG(DM)
      ADD=0.
      IF(DD.GT.0.) ADD=ALOG(DD)
      YLOG=A*(ADM-ADD)
      IF(YLOG.LT.-10.) GO TO 10
      IF(YLOG.GT.10.)  GO TO 20
        DNETVENUS=DD*(1.+EXP(YLOG))**(1/A)
        GO TO 50
10    CONTINUE
        DNETVENUS=DD
        GO TO 50
20    CONTINUE
        DNETVENUS=DM
        GO TO 50
50    CONTINUE

      END FUNCTION DNETVENUS


      pure real FUNCTION CCORVENUS(ALT, R,H1,ZH)
C        CHEMISTRY/DISSOCIATION CORRECTION
      real, intent(in) :: alt, R,H1,ZH

      E=(ZH-ALT)/H1
      IF(E.GT.170.) GO TO 10
      IF(E.LT.-170.) GO TO 20
        EX=EXP(E)
        CCORVENUS=(1.+R*EX)/(1.+EX)
        return
10      CCORVENUS=R
        return
20      CCORVENUS=1.

      END FUNCTION CCORVENUS


      SUBROUTINE TURBO(DD,DM,ZB,ZH,XM,XMM,TZ)
C        ESTIMATE TURBOPAUSE HEIGHT
      COMMON/PARMB/GSURF,RE
      DATA RGAS/831.4/
      GZB=GSURF/(1.+ZB/RE)**2
      ZH=ZB+RGAS*TZ/GZB/(XM-XMM)*ALOG(DD/DM)

      END SUBROUTINE TURBO

C--------------------------------------------------------------------
      FUNCTION GLOBV(YRD,LAT,TLOC,F107A,F107,P)
C        CALCULATE GLOBAL VARIATIONS (12/12/80)
C         USES SOLAR ZENITH ANGLE FOR SYMMETRICAL TERMS
      INTEGER YRD
      REAL LAT
      DIMENSION PLG(6,6),P(15),SV(15),T(15)
      COMMON/CSW/SW(15),ISW
      DATA DGTR/1.74533E-2/, XL/-1000./,SV/15*1./,TLL/1000./,HR/.2618/
      IF(ISW.NE.64999) CALL TSELECVENUS(SV)
c       comment par gg,0 par GG
c      IF( XL.EQ.LAT.AND.TLL.EQ.TLOC) then
c      PLG(2,1)=0.
c      PLG(3,1)=0.
c      PLG(4,1)=0.
c      PLG(5,1)=0.
c      PLG(6,1)=0.
c      PLG(2,2)=0.
c      PLG(3,3)=0.
c      PLG(4,4)=0.
c      PLG(5,5)=0.
c      PLG(6,6)=0.
c      write(*,*),"Passage par les boucles"
c      GO TO 15
c      endif
C      CALCULATE LEGENDRE POLYNOMIALS
      S=COS(LAT*DGTR)
C        COSINE OF SOLAR ZENITH ANGLE
      C=COS(LAT*DGTR)*COS(HR*(TLOC-12.))
      S2=S*S
      C2=C*C
      C4=C2*C2
      PLG(2,1)=C
      PLG(3,1)=0.5*(3.*C2-1.)
      PLG(4,1)=0.5*(5.*C2-3.)*C
      PLG(5,1)=(35.*C4-30.*C2+3.)/8.
      PLG(6,1)=(63.*C4-70.*C2+15.)*C/8.
      PLG(2,2)=S
      PLG(3,3)=3.*S2
      PLG(4,4)=15.*S2*S
      PLG(5,5)=105.*S2*S2
      PLG(6,6)=945.*S2*S2*S
       XL=LAT
c15    CONTINUE
c      IF(TLL.EQ.TLOC) GO TO 16
      STLOC=SIN(HR*TLOC)
      S2TLOC=SIN(2.*HR*TLOC)
      S3TLOC=SIN(3.*HR*TLOC)
      S4TLOC=SIN(4.*HR*TLOC)
      S5TLOC=SIN(5*HR*TLOC)
c16    CONTINUE
      TLL=TLOC
      T(1)=P(14)*(F107A-200.)+P(15)*(F107-F107A)
C
C  THE FOLLOWING TEST INADVERTENTLY TREATS HELIUM DIFFERENTLY THAN THE
C  OTHER SPECIES WITH REGARD TO THE F107 EFFECT.
C  THE .GT. SHOULD BE .NE. BUT IS NOT CHANGED HERE TO BE CONSISTANT WITH
C  THE COEFFICIENTS AND TABLES AS ORIGINALLY DERIVED AND PUBLISHED.  7/25/83
C
      IF(P(2).GT.0.) GO TO 18
        F=1.
        GO TO 20
18    F=1.+T(1)/(P(2)-P(4)+P(6)-P(8)+P(10))
20    T(2)=(P(2)*PLG(2,1)
     &    +P(3)*PLG(2,2)*STLOC)*F
      T(3)=(P(4)*PLG(3,1)
     &    +P(5)*PLG(3,3)*S2TLOC)*F
      T(4)=(P(6)*PLG(4,1)
     &    +P(7)*PLG(4,4)*S3TLOC)*F
      T(5)=(P(8)*PLG(5,1)+P(9)*PLG(5,5)*S4TLOC)*F
      T(6)=(P(10)*PLG(6,1)+P(11)*PLG(6,6)*S5TLOC)*F
      G=0.
      DO 50 I=1,6
50    G=G+SW(I)*T(I)
C      write(*,*),'GLOB G ',G,' SW1 ',SW(1),' T1 ',T(I)
      IF(G.gt.10) then
c      G=0
c      DO 42 I=1,6
c        G=G+SW(I)*T(I)
c42      write(*,*),'GLOB G ',G,' SW ',SW(I),' T ',T(I)
      T(2)=(P(2)*PLG(2,1)+P(3)*PLG(2,2)*STLOC)*F
      write(*,*)'T2',T(2),' P2 ',P(2),' P3', P(3)
      write(*,*)'F ',F,' STLOC ',STLOC,' PLG22 ', PLG(2,2)
      write(*,*)' PLG21 ', PLG(2,1) ,' S ', S, 'LAT',LAT
      write(*,*)'DGTR',DGTR

      endif
      GLOBV=G

      END FUNCTION GLOBV


      SUBROUTINE TSELECVENUS(SV)
      DIMENSION SV(15)
      COMMON/CSW/SW(15),ISW
      DO 100 I=1,15
100   SW(I)=SV(I)
      ISW=64999

      END SUBROUTINE TSELECVENUS
C --------------------------------------------------------------------------
C TVTS3.FOR  test driver for VTS3
C --------------------------------------------------------------------------

!      PROGRAM vts
!      DIMENSION D(7),T(2),RES(62),RES1(62),RES2(62),RES3(62),RES4(62),
!     . RES5(62),RES6(62),RES7(62),RES8(62)
!        INTEGER I,J
!        integer resu
!        DATA resu/35/
!        open(unit=resu,file='res.txt')
c        resu=
c        REAL altgg
c      WRITE(*,*),'Bonjour les gens'
c      CALL VTS3(250.,0,0.,0.,200.,200.,48,D,T)
c      WRITE(6,133) D,T
c      WRITE(*,*),'Bonjour les gens'
c      CALL VTS3(250.,0,0.,12.,200.,200.,48,D,T)
c      WRITE(6,133) D,T
c      WRITE(*,*),'Bonjour les gens'
c      CALL VTS3(110.,78343,-31.3,6.8,166.,190.,48,D,T)
c      WRITE(6,133) D,T
c      WRITE(*,*),'Bonjour les gens'
c      CALL VTS3(130.,78343,-37.9,8.5,166.,190.,48,D,T)
c      WRITE(6,133) D,T
c      WRITE(*,*),'Bonjour les gens'
c      CALL VTS3(150.3,79097,16.3,5.3,179.2,200.3,48,D,T)
c      WRITE(6,133) D,T

c       On va travailler sur les altitudes croissantes

!        altgg=405.
!        DO 4242 I=1,62,1
c         write(*,*),I
!         altgg=altgg-5.
!         CALL VTS3(altgg,0,45.,12.,200.,200.,48,D,T)
!         write(*,*),altgg,D,T
!         RES(I)=altgg
!         RES1(I)=T(2)
!         RES2(I)=D(2)
!         RES3(I)=D(3)
!         RES4(I)=D(4)
!         RES5(I)=D(5)
!         RES6(I)=D(6)
!         RES7(I)=D(7)

!4242    continue
!        write(*,*),'Fin gg'
!        rewind(resu)
C        DO 4243 J=0,11
C         write(resu,*),RES(5*J+1),RES(5*J+2),RES(5*J+3),RES(5*J+4),
C     .  RES(5*J+5)
C4243    continue
c         write(resu,*),RES(61),RES(62)
!        write(resu,*),'Altitudes'
!        write(resu,1042),RES
!        write(resu,*),'Temperature neutre'
!        write(resu,1020),RES1
!        write(resu,*),'CO2'
!        write(resu,1020),RES2
!        write(resu,*),'O'
!        write(resu,1020),RES3
!        write(resu,*),'CO'
!        write(resu,1020),RES4
!        write(resu,*),'HE'
!        write(resu,1020),RES5
!        write(resu,*),'N'
!        write(resu,1020),RES5
!        write(resu,*),'N2'
!        write(resu,1020),RES5

c        write(resu,1020),RES
c        write(resu,*),RES


!1042    format(5f10.2)

!1020    format(5(1pe14.6))


!        STOP
c133   FORMAT(1X,1P7E10.2/2E10.3)
! !     END
C -------------------------------------------------------------------------
C output of TVTS3
C ----------------------------------------------------------------------------
C   3.70E-18  1.09E-07  2.78E+03  2.39E-02  5.46E+05  7.53E+01  1.23E-02
C 1.338E+02 1.338E+02
C   5.52E-16  3.58E+03  1.95E+07  2.93E+05  9.63E+05  4.88E+05  5.90E+04
C 3.081E+02 3.080E+02
C   4.96E-09  6.60E+13  3.32E+10  3.83E+10  3.92E+08  3.73E+05  3.07E+12
C 2.542E+02 1.709E+02
C   3.27E-11  4.12E+11  1.47E+10  1.27E+10  2.29E+07  3.44E+08  3.46E+10
C 2.933E+02 1.844E+02
C   6.83E-14  2.79E+08  1.20E+09  1.98E+08  3.22E+07  9.32E+06  1.39E+08
C 1.754E+02 1.441E+02
