c
C IRIS12.FOR ---------------------------------------- OCTOBER 1991
C
C*****************************************************************
C CHANGES FROM  IRIS11.FOR  TO   IRIS12.FOR:
C    - CIRA-1986 INSTEAD OF CIRA-1972 FOR NEUTRAL TEMPERATURE
C    - 10/30/91 VNER FOR NIGHTTIME LAY-VERSION:  ABS(..)
C    - 10/30/91 XNE(..) IN CASE OF LAY-VERSION
C    - 10/30/91 CHANGE SSIN=F/T TO IIQU=0,1,2
C    - 10/30/91 Te > Ti > Tn ENFORCED IN FINAL PROFILE
C    - 10/30/91 SUB ALL NAMES WITH 6 OR MORE CHARACTERS
C    - 10/31/91 CORRECTED HF1 IN HST SEARCH:  NE(HF1)>NME
C    - 11/14/91 C1=0 IF NO F1-REGION
C    - 11/14/91 CORRECTED HHMIN AND HZ FOR LIN. APP.
C
C*****************************************************************
C********* INTERNATIONAL REFERENCE IONOSPHERE (IRI). *************
C*****************************************************************
C****************    OCTOBER 1991     ****************************
C****************     SUBROUTINE      ****************************
C*****************************************************************
C
C
       SUBROUTINE IRIS12(JF,JMAG,ALATI,ALONG,RZ12,MMDD,DHOUR,
     &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
c
cf2py intent(out) outf,oarr
C-----------------------------------------------------------------
C INTERNATIONAL REFERENCE IONOSPHERE 1991
C
C INPUT:  JMAG=0/1      GEODETIC/GEOMAGNETIC LATITUDE AND LONGITUDE
C         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C         RZ12 (-COV)   12-MONTHS-RUNNING MEAN OF SOLAR SUNSPOT NUMBER
C                          (OR EQUIVALENT F10.7 SOLAR RADIO FLUX AS
C                          NEGATIVE NUMBER)
C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL 
C                          HOURS
C         HEIBEG,       BEGIN, END, AND STEPWIDTH OF HEIGHT RANGE  
C          HEIEND,HEISTP   IN KM (MAXIMUM NUMBER OF STEPS IS 50 !!)
C         JF(1:12)      TRUE/FALSE FLAGS FOR SEVERAL OPTIONS
C          JF(1)=.TRUE.[.FALSE.]   ELECTRON DENSITY IS [NOT] CALCULATED
C          JF(2)=T[F]    TEMPERATURES ARE [NOT] CALCULATED
C          JF(3)=T[F]    ION COMPOSITION IS [NOT] CALCULATED
C          JF(4)=T[F]    B0 FROM TABLE [FROM GULYEAVA 1987]
C          JF(5)=T[F]    F2 PEAK FROM CCIR [FROM URSI]
C          JF(6)=T[F]    ION COMP. STANDARD [DANILOV-YAICHNIKOV-1985]
C          JF(7)=T[F]    STAND. IRI TOPSIDE [IRI-79]
C          JF(8)=T[F]    NMF2 PEAK MODEL [INPUT VALUES]
C          JF(9)=T[F]    HMF2 PEAK MODEL [INPUT VALUES]
C          JF(10)=T[F]   TE MODEL [TE-NE MODEL WITH NE INPUT]
C          JF(11)=T[F]   NE STANDARD [LAY-FUNCTIONS VERSION]
C          JF(12)=T[F]   MESSAGE ARE WRITTEN TO UNIT=6 [=12]
C
C  JF(1:11)=.TRUE. GENERATES THE STANDARD IRI-90 PARAMETERS.
C  IF YOU SET JF(8)=.FALSE., THAN YOU HAVE TO PROVIDE THE F2 PEAK
C  NMF2/M-3 OR FOF2/MHZ IN OARR(1). SIMILARLY, IF YOU SET JF(9)=
C  .FALSE., THAN YOU HAVE TO PROVIDE THE F2 PEAK HEIGHT HMF2/KM IN
C  OARR(2). IF YOU SET JF(10)=.FALSE., THAN YOU HAVE TO PROVIDE THE
C  ELECTRON DENSITY IN M-3 AT 300KM AND/OR 400KM AND/OR 600KM IN
C  OARR(3), OARR(4), AND OARR(5). IF YOU WANT TO USE THIS OPTION AT
C  ONLY ONE OF THE THREE ALTITUDES, THAN SET THE DENSITIES AT THE
C  OTHER TWO TO ZERO.
C
C  OUTPUT:  OUTF(1:11,1:50)   IRI PROFILES
C              OUTF(1,*)  ELECTRON DENSITY/M-3
C              OUTF(2,*)  NEUTRAL TEMPERATURE/K
C              OUTF(3,*)  ION TEMPERATURE/K
C              OUTF(4,*)  ELECTRON TEMPERATURE/K
C              OUTF(5,*)  PERCENTAGE OF O+ IONS IN %
C              OUTF(6,*)  PERCENTAGE OF H+ IONS IN %
C              OUTF(7,*)  PERCENTAGE OF HE+ IONS IN %
C              OUTF(8,*)  PERCENTAGE OF O2+ IONS IN %
C              OUTF(9,*)  PERCENTAGE OF NO+ IONS IN %
C                 AND, IF JF(6)=.FALSE.:
c       inversion d'indice dans les commentaires... Corrige (jl, 2006)
C              OUTF(10,*)  PERCENTAGE OF CLUSTER IONS IN %
C              OUTF(11,*)  PERCENTAGE OF N+ IONS IN %
C
C            OARR(1:30)   ADDITIONAL OUTPUT PARAMETERS         
C              OARR(1) = NMF2/M-3        OARR(2) = HMF2/KM
C              OARR(3) = NMF1/M-3        OARR(4) = HMF1/KM
C              OARR(5) = NME/M-3         OARR(6) = HME/KM
C              OARR(7) = NMD/M-3         OARR(8) = HMD/KM
C              OARR(9) = HHALF/KM        OARR(10) = B0/KM
C              OARR(11) =VALLEY-BASE/M-3 OARR(12) = VALLEY-TOP/KM
C              OARR(13) = TE-PEAK/K      OARR(14) = TE-PEAK HEIGHT/KM
C              OARR(15) = TE-MOD(300KM)  OARR(16) = TE-MOD(400KM)/K
C              OARR(17) = TE-MOD(600KM)  OARR(17) = TE-MOD(1400KM)/K
C              OARR(18) = TE-MOD(3000KM) OARR(19) = TE(120KM)=TN=TI/K
C              OARR(20) = TI-MOD(430KM)  OARR(21) = X/KM, WHERE TE=TI
C              OARR(22) = SOLAR ZENITH ANGLE/DEG
C              OARR(23) = SUN DECLINATION/DEG
C              OARR(24) = DIP            OARR(25) = DIP LATITUDE
C              OARR(26) = MODIFIED DIP LATITUDE
C              OARR(27:30) FREE
C-------------------------------------------------------------------
C*** THIS PROGRAM PRODUCES PROFILES OF                         ***
C***      ELECTRON DENSITY                                     ***
C***      NEUTRAL TEMPERATURE (CIRA 86)                        ***
C***      ELECTRON TEMPERATURE                                 ***
C***      ION TEMPERATURE                                      ***
C***      RELATIVE PERCENTAGE DENSITIES OF THE IONS            ***
C***           ATOMIC OXYGEN, HYDROGEN, HELIUM,                ***
C***           MOLECULAR OXYGEN AND NITROGEN OXYD (NO+)        ***
C*****************************************************************
C*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
C***     ELECTRON DENSITY         60/80 KM       1000 KM       ***
C***     TEMPERATURES              120 KM        3000 KM       ***
C***     ION DENSITIES             100 KM        1000 KM       ***
C*****************************************************************
C*     --------------------ADDRESSES------------------------     *
C*     I  PROF. K. RAWER              DR. D. BILITZA       I     *
C*     I  HERRENSTR. 43               GSFC/NSSDC CODE 633  I     *
C*     I  D-7801 MARCH                GREENBELT MD 20771   I     *
C*     I  F.R.G.                      USA                  I     *
C*     -----------------------------------------------------     *
C*****************************************************************
C*****************************************************************
C*****************************************************************
C*********       ALL ANGLES ARE IN DEGREE           **************
C*********       ALL DENSITIES ARE IN M-3           **************
C*********       ALL ALTITUDES ARE IN KM            **************
C*********     ALL TEMPERATURES ARE IN KELVIN       **************
C*********     ALL TIMES ARE IN DECIMAL HOURS       **************
C*****************************************************************
C********************  OPTIONS  **********************************
C*****************************************************************
C* FOR HMF2=0 OR FOF2=0 THE F2 PEAK VALUES ARE CALCULATED WITH   *
C* THE CCIR OR URSI MODELS. THE CCIR COEFFICIENT SET FOR THE     *
C* MONTH "mm" IS EXPECTED IN THE BINARY FILE "CCIRmm.BIN" AND    *
C* THE URSI SET IN "URSImm.BIN". IF YOU USE THE ASCII CODED      *
C* FILES "CCIRmm.ASC", YOU HAVE TO INCORPORATE THE CHANGES       *
C* INDICTED IN PROGRAM SECTION ENTITLED "READ CCIR COEFFICIENT   *
C* SET FOR CHOSEN MONTH."                                        * 
C*****************************************************************
C*****************************************************************
C*****************************************************************
      INTEGER           EGNR,AGNR,DAYNR,DDO,DO2,SEASON,SEADAY
      REAL              LATI,LONGI,MO2,MO,MODIP,NMF2,MAGBR
      REAL              NMF1,NME,NMD,NEI,MM,MLAT,MLONG,NOBO2
      CHARACTER*80      FILNAM
      DIMENSION  F(3),RIF(4),E(4),XDELS(4),DNDS(4)
      DIMENSION  FF0(988),XM0(441),F2(13,76,2),FM3(9,49,2)
      DIMENSION  AMP(4),HXL(4),SCL(4),B0B1(5)
      DIMENSION  CTN(3),CTNN(3),XSM(4),MM(5),DTI(4)
      DIMENSION  AHH(7),STTE(6),DTE(5),ATE(7),TEA(6),HOA(3),XNAR(3)
      DIMENSION  PG1O(80),PG2O(32),PG3O(80),PF1O(12),PF2O(4),PF3O(12)
      DIMENSION  HO(4),MO(5),DDO(4),HO2(2),MO2(3),DO2(2),DION(7)
      DIMENSION  OUTF(11,50),OARR(30)
      LOGICAL           EXT,SCHALT,NIGHT,TCON(3)
      LOGICAL           F1REG,FOF2IN,HMF2IN,URSIF2,LAYVER,DY,GULB0
      LOGICAL           NODEN,NOTEM,NOION,TENEOP
      LOGICAL           OLD79,TOPSI,BOTTO,BELOWE,JF(12),URSIFO
      COMMON    /BLOCK1/HMF2,NMF2,HMF1           /CONSTiri/UMR
     &          /BLOCK2/B0,B1,C1      /BLOCK3/HZ,T,HST,STR
     &          /BLOCK4/HME,NME,HEF   /BLOCK5/NIGHT,E
     &          /BLOCK6/HMD,NMD,HDX   /BLOCK7/D1,XKK,FP30,FP3U,FP1,FP2
     &          /BLOCK8/HS,TNHS,XSM,MM,DTI,MXSM         
     &          /BLOTN/XSM1,TEXOS,TLBDH,SIGMA /BLOTE/AHH,ATE1,STTE,DTE
     &          /BLO10/BETA,ETA,DELTA,ZETA       /ARGEXP/ARGMAX
      EXTERNAL          XE1,XE2,XE3,XE4,XE5,XE6,TEDER
      DATA  HOA  /300.,400.,600./,   XNAR       /3*0.0/,   
     &      XDELS   /3*5.,10./,      DNDS   /.016,.01,2*.016/,
     &      DDO   /9,5,5,25/,        DO2        /5,5/,
     &      B0B1  /.755566,.778596,.797332,.812928,.826146/
C
C PROGAM CONSTANTS
C
        icalls=icalls+1
        ARGMAX=88.0
        UMR=ATAN(1.0)*4./180.
        ALOG2=ALOG(2.)
        ALG100=ALOG(100.)
        ISTART=1
        NUMHEI=INT((HEIEND-HEIBEG)/HEISTP)+1
        IF(NUMHEI.GT.50) NUMHEI=50
C
C Code inserted to aleviate block data problem for PC version.
C Thus avoiding DATA statement with parameters from COMMON block.
C
        AHH(1)=120.
        AHH(2)=0.
        AHH(3)=300.
        AHH(4)=400.
        AHH(5)=600.
        AHH(6)=1400.
        AHH(7)=3000.
        DTE(1)=5.
        DTE(2)=5.
        DTE(3)=10.
        DTE(4)=20.
        DTE(5)=20.
        DTI(1)=10.
        DTI(2)=10.
        DTI(3)=20.
        DTI(4)=20.
C
C FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
C AGNR=OUTPUT (OUTPUT IS DISPLAYED OR STORED IN FILE OUTPUT.IRI)...
C IUCCIR=UNIT NUMBER FOR CCIR COEFFICIENTS ........................
C
      MONITO=6
      IUCCIR=10
      KONSOL=6
      IF(.not.JF(12)) KONSOL=12

c
c selection of density and ion composition options ..................
c

      NODEN=(.NOT.JF(1))
      NOTEM=(.NOT.JF(2))
      NOION=(.NOT.JF(3))
      DY=(.NOT.JF(6))
      LAYVER=(.NOT.JF(11))
      OLD79=(.NOT.JF(7))
      GULB0=(.NOT.JF(4))
c
c f peak density ....................................................
c
      FOF2IN=(.NOT.JF(8))
       IF(FOF2IN) THEN
          AFOF2=OARR(1)
          IF(AFOF2.GT.100.) AFOF2=SQRT(AFOF2/1.24E10)
          ENDIF
      URSIF2=(.NOT.JF(5))
c
c f peak altitude ..................................................
c
      HMF2IN=(.NOT.JF(9))
       IF(HMF2IN) AHMF2=OARR(2)
c
C TE-NE MODEL OPTION ..............................................
C
      TENEOP=(.NOT.JF(10))
        IF(TENEOP) THEN
           DO 8154 JXNAR=1,3
              XNAR(JXNAR)=OARR(JXNAR+2)
              TCON(JXNAR)=.FALSE. 
8154          IF(XNAR(JXNAR).GT.0.) TCON(JXNAR)=.TRUE. 
           ENDIF

      if(icalls.gt.1) goto 8201
c       write(*,*) '*** IRI parameters are being calculated ***'
      if(NODEN) goto 2889
        if(LAYVER) write(*,*) 'Ne, E-F: The LAY-Version is ',
     &    'prelimenary. Erroneous profile features can occur.'
        if(GULB0) write(*,*) 'Ne, B0: Bottomside thickness is ',
     &    'obtained with Gulyaeva-1987 model.'
        if(OLD79) write(*,*) 'Ne: Using IRI-79. Correction',
     &    ' of equatorial topside is not included.'
c       if(HMF2IN) write(*,*) 'Ne, hmF2: Input values are used.'
        if(FOF2IN) then
c         write(*,*) 'Ne, foF2: Input values are used.'
          goto 2889
          endif
c       if(URSIF2) then
c         write(*,*) 'Ne, foF2: URSI model is used.'
c       else
c         write(*,*) 'Ne, foF2: CCIR model is used.'
c       endif
2889  if((.not.NOION).and.(DY)) 
     &     write(*,*) 'Ion Com.: Using Danilov-Yaichnikov-1985.'
      if((.not.NOTEM).and.(TENEOP))
     &     write(*,*) 'Te: Temperature-density correlation is used.'
8201    continue
C
C CALCULATION OF MEAN F10.7CM SOLAR RADIO FLUX (COV)................
C CALCULATION OF RESTRICTED SOLAR ACTIVITIES (RG,COVG)..............
C
      IF(RZ12.GT.0.0) THEN
        R=RZ12
        COV=63.75+R*(0.728+R*0.00089)
      ELSE
        COV=-RZ12
        R=33.52*(SQRT(COV+85.12)-12.2)
      ENDIF
      RG=R
      COVG=COV
      IF(R.GT.150.) RG=150.
      IF(COV.GT.193.) COVG=193.
C
C CALCULATION OF GEOG. OR GEOM. COORDINATES IN DEG....................
C CALCULATION OF MAGNETIC INCLINATION (DIP), DECLINATION (DEC)........
C   DIP LATITUDE (MAGBR) AND MODIFIED DIP (MODIP). ALL IN DEGREE......
C
        IF(JMAG.GT.0) THEN
           MLAT=ALATI
           MLONG=ALONG
        ELSE
           LATI=ALATI
           LONGI=ALONG
        ENDIF
        CALL GGM(JMAG,LONGI,LATI,MLONG,MLAT)
        ABSLAT=ABS(LATI)
        CALL FIELDG(LATI,LONGI,300.0,XMA,YMA,ZMA,BET,DIP,DEC,MODIP)
        MAGBR=ATAN(0.5*TAN(DIP*UMR))/UMR
        ABSMLT=ABS(MLAT)
        ABSMDP=ABS(MODIP)
        ABSMBR=ABS(MAGBR)
C
C CALCULATION OF SEASON (SUMMER=2, WINTER=4)..........................
C CALCULATION OF DAY OF YEAR AND SUN DECLINATION......................
C
        if(MMDD.lt.0) then
                DAYNR=-MMDD
                call MODA(1,MONTH,IDAY,DAYNR)
        else
                MONTH=MMDD/100
                IDAY=MMDD-MONTH*100
                call MODA(0,MONTH,IDAY,DAYNR)
        endif
      SEASON=INT((DAYNR+45.0)/92.0)
      IF(SEASON.LT.1) SEASON=4
      NSESON=SEASON
      seaday=daynr
      IF(LATI.GT.0.0) GOTO 5592
        SEASON=SEASON-2
        IF(SEASON.LT.1) SEASON=SEASON+4
        seaday=daynr+183
        if(seaday.gt.366) seaday=seaday-366
C
C CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG).........................
C NOON VALUE (XHINON).................................................
C
5592  IF(DHOUR.GT.24.1) THEN
        UT=DHOUR-25.
        HOUR=UT+LONGI/15.
        IF(HOUR.GT.24.) HOUR=HOUR-24.
      ELSE
        HOUR=DHOUR
        UT=HOUR-LONGI/15.
        IF(UT.LT.0.) UT=UT+24.
      ENDIF

        CALL SOCO(DAYNR,HOUR,LATI,LONGI,SUNDEC,XHI,SAX,SUX)
        CALL SOCO(DAYNR,12.0,LATI,LONGI,SUNDE1,XHINON,SAXNON,SUXNON)

                NIGHT=.FALSE.
        if(abs(sax).gt.25.0) then
                if(sax.lt.0.0) NIGHT=.TRUE.
                goto 1334
                endif
        if(SAX.le.SUX) goto 1386
        if((hour.gt.sux).and.(hour.lt.sax)) night=.true.
        goto 1334
1386    IF((HOUR.GT.SUX).OR.(HOUR.LT.SAX)) NIGHT=.TRUE.
C
C CALCULATION OF ELECTRON DENSITY PARAMETERS................
C
1334  HNEA=65.
      IF(NIGHT) HNEA=80.
      HNEE=2000.
      IF(NODEN) GOTO 4933
      DELA=4.32
      IF(ABSMDP.GE.18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)
      DELL=1+EXP(-(ABSLAT-20.)/10.)
C!!!!!!! F-REGION PARAMETERS AND E-PEAK !!!!!!!!!!!!!!!!!!!!!!!!!!
      FOE=FOEEDI(COV,XHI,XHINON,ABSLAT)
      NME=1.24E10*FOE*FOE
      HME=105.0
      IF((FOF2IN).AND.(HMF2IN)) GOTO 501
      IF(URSIF2.NEqv.URSIFO) GOTO 7797 
      IF((MONTH.EQ.MONTHO).AND.(RG.EQ.RGO)) GOTO 4292
      IF(MONTH.EQ.MONTHO) GOTO 4291
C
C READ CCIR COEFFICIENT SET FOR CHOSEN MONTH....................
C
7797    WRITE(FILNAM,104)MONTH+10
c104     FORMAT('dir.iri90/CCIR',I2,'.BIN')
c1344    OPEN(IUCCIR,FILE=FILNAM(1:lenc(filnam)),STATUS='OLD',
c    &          ERR=8448,FORM='UNFORMATTED')
        lread=1
c       READ(IUCCIR) F2,FM3
C !!!!!!!!!!!!!!! FOR ASCII CODED CCIR COEFFIECENTS FILES
C !!!!!!!!!!!!!!! SUBSTITUTE THE LAST 3 STATEMENTS BY:
 104     FORMAT('datairi/ccir',I2,'.asc')
 1344    OPEN(IUCCIR,FILE=FILNAM(1:lenc(filnam)),STATUS='OLD',
     &          ERR=8448,FORM='FORMATTED')
         READ(IUCCIR,4689) F2,FM3
 4689    FORMAT(1X,4E15.8)
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CLOSE(IUCCIR)
C
C READ URSI COEFFICIENT SET FOR CHOSEN MONTH....................
C
        if(URSIF2) then
c         WRITE(FILNAM,1144) MONTH+10
c1144    FORMAT('dir.iri90/URSI',I2,'.BIN')
c1244      OPEN(IUCCIR,FILE=FILNAM(1:lenc(filnam)),
c    &          STATUS='OLD',ERR=8448,FORM='UNFORMATTED')
         WRITE(FILNAM,1144) MONTH+10
1144    FORMAT('datairi/ursi',I2,'.asc')
1244      OPEN(IUCCIR,FILE=FILNAM(1:lenc(filnam)),
     &          STATUS='OLD',ERR=8448,FORM='FORMATTED')
          lread=2
          READ(IUCCIR,*) F2
          CLOSE(IUCCIR)
          endif
        MONTHO=MONTH
        URSIFO=URSIF2
        GOTO 4291

8448    WRITE(MONITO,8449) FILNAM
c8449   FORMAT(1X////,
c    &    ' !!!!   The file ',A45,/,
c    &    ' !!!!   is not in your directory,'/
c    &    ' !!!!   try a different diskette (enter: 1),'/
c    &    ' !!!!   or exit (enter: 0)')
c       read(*,*) idisk
c       if(idisk.eq.1) goto (1344,1244) lread
8449    FORMAT(1X////,
     &    ' !!!!   The file ',A45,/,
     &    ' !!!!   is not in your directory,'/
     &    ' !!!!   IRI is NOT excecuted,')
        GOTO 3330
C
C LINEAR INTERPOLATION IN SOLAR ACTIVITY
C
4291    RR2=RG/100.
        RR1=1.-RR2
        DO 20 I=1,76
        DO 20 J=1,13
        K=J+13*(I-1)
20      FF0(K)=F2(J,I,1)*RR1+F2(J,I,2)*RR2
        DO 30 I=1,49
        DO 30 J=1,9
        K=J+9*(I-1)
30      XM0(K)=FM3(J,I,1)*RR1+FM3(J,I,2)*RR2
        RGO=RG

4292  CALL F2OUT(MODIP,LATI,LONGI,FF0,XM0,UT,YFOF2,XM3000)

501     IF(FOF2IN) THEN
          FOF2=AFOF2
        ELSE
          FOF2=YFOF2
        ENDIF
        NMF2=1.24E10*FOF2*FOF2

        IF(HMF2IN) THEN
          HMF2=AHMF2
        ELSE
          HMF2=HMF2ED(MAGBR,RG,FOF2/FOE,XM3000)
        ENDIF

        TOPSI=(HEIEND.GT.HMF2)
        BOTTO=((HEIEND.GE.HME).AND.(HEIBEG.LE.HMF2))
        BELOWE=(HEIBEG.LT.HME)
c
c topside profile parameters .............................
c
        IF(.NOT.TOPSI) GOTO 1501
      COS2=COS(MLAT*UMR)
      COS2=COS2*COS2
      FLU=(COVG-40.0)/30.0
      IF(OLD79) then
        ETA1=-0.0070305*COS2
      else
        EX=EXP(-MLAT/15.)
        EX1=EX+1
        EPIN=4.*EX/(EX1*EX1)
        ETA1=-0.02*EPIN
      endif
      ETA=0.058798+ETA1+FLU*(-0.014065+0.0069724*COS2)+
     &(0.0024287+0.0042810*COS2-0.00015280*FOF2)*FOF2
      ZETA=0.078922-0.0046702*COS2+FLU*(-0.019132+0.0076545*COS2)+
     &(0.0032513+0.0060290*COS2-0.00020872*FOF2)*FOF2
      BETA=-128.03+20.253*COS2+FLU*(-8.0755-0.65896*COS2)+(0.44041
     &+0.71458*COS2-0.042966*FOF2)*FOF2
      Z=EXP(94.45/BETA)
      Z1=Z+1
      Z2=Z/(BETA*Z1*Z1)
      DELTA=(ETA/Z1-ZETA/2.0)/(ETA*Z2+ZETA/400.0)
c
c bottomside profile parameters .............................
C
1501    HMF1=HMF2
        HZ=HMF2
        HEF=HME
        IF(.not.BOTTO) GOTO 2727  
        B1=3.0
C!!!!!!! INTERPOLATION FOR B0 OUT OF ARRAY B0F !!!!!!!!!!!!!!!!!!!!!
        if(GULB0) then
          call ROGUL(SEADAY,XHI,SEAX,GRAT)
          if(NIGHT) GRAT=0.91-HMF2/4000.
          B0CNEW=HMF2*(1.-GRAT)
          B0=B0CNEW/B0B1(1)
        else
          B0 = B0POL(HOUR,SAX,SUX,SEASON,RG,DELA)
        endif
C!!!!!!! F1-REGION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      F1REG=.FALSE.
      HMF1=0.
      PNMF1=0.
      C1=0.  
      IF(NIGHT.OR.(SEASON.EQ.4)) GOTO 150
        FOF1=FOF1ED(ABSMBR,R,XHI)
        IF(FOF1.LT.1.E-3) GOTO 150
          F1REG=.TRUE.
          C1=.09+.11/DELA
          PNMF1=1.24E10*FOF1*FOF1
150   NMF1=PNMF1
C!!!!!!! PARAMETER FOR E AND VALLEY-REGION !!!!!!!!!!!!!!!!!!!!!
      XDEL=XDELS(SEASON)/DELA
      DNDHBR=DNDS(SEASON)/DELA
      HDEEP=HPOL(HOUR,10.5/DELA,28.,SAX,SUX,1.,1.)
      WIDTH=HPOL(HOUR,17.8/DELA,45.+22./DELA,SAX,SUX,1.,1.)
      DEPTH=HPOL(HOUR,XDEL,81.,SAX,SUX,1.,1.)
      DLNDH=HPOL(HOUR,DNDHBR,.06,SAX,SUX,1.,1.)
      IF(DEPTH.LT.1.0) GOTO 600
        IF(NIGHT) DEPTH=-DEPTH
        CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
        IF(.NOT.EXT) GOTO 667
          WRITE(KONSOL,650)
650   FORMAT(1X,'*NE* E-REGION VALLEY CAN NOT BE MODELLED')
600       WIDTH=.0
667   HEF=HME+WIDTH
      VNER = (1. - ABS(DEPTH) / 100.) * NME
c
c Parameters below E  .............................
c
2727    IF(.not.BELOWE) GOTO 2726
C!!!!!!!D-REGION PARAMETER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NMD=XMDED(XHI,R,4.0E8)
      HMD=HPOL(HOUR,81.0,88.0,SAX,SUX,1.,1.)
      F(1)=HPOL(HOUR,0.02+0.03/DELA,0.05,SAX,SUX,1.,1.)
      F(2)=HPOL(HOUR,4.6,4.5,SAX,SUX,1.,1.)
      F(3)=HPOL(HOUR,-11.5,-4.0,SAX,SUX,1.,1.)
      FP1=F(1)
      FP2=-FP1*FP1/2.0
      FP30=(-F(2)*FP2-FP1+1.0/F(2))/(F(2)*F(2))
      FP3U=(-F(3)*FP2-FP1-1.0/F(3))/(F(3)*F(3))
      HDX=HMD+F(2)
      X=HDX-HMD
      XDX=NMD*EXP(X*(FP1+X*(FP2+X*FP30)))
      DXDX=XDX*(FP1+X*(2.0*FP2+X*3.0*FP30))
      X=HME-HDX
      XKK=-DXDX*X/(XDX*ALOG(XDX/NME))
      D1=DXDX/(XDX*XKK*X**(XKK-1.0))
C
C SEARCH FOR HMF1 ..................................................
C
2726    IF(.not.BOTTO) GOTO 4933
        if(LAYVER) goto 6153
924     IF(.not.F1REG) GOTO 380
        XE2H=XE2(HEF)
      CALL REGFA1(HEF,HMF2,XE2H,NMF2,0.001,NMF1,XE2,SCHALT,HMF1)
        IF(.not.SCHALT) GOTO 380
          WRITE(KONSOL,11)
11    FORMAT(1X,'*NE* HMF1 IS NOT EVALUATED BY THE FUNCTION XE2')
        IREGFA=1
c
c change B1 and try again ..........................................
c
9244    IF(B1.GT.4.5) GOTO (7398,8922) IREGFA
                B1=B1+0.5
                WRITE(KONSOL,902) B1-0.5,B1
902   FORMAT(6X,'CORR.: B1(OLD)=',F4.1,' B1(NEW)=',F4.1)
                IF(GULB0) then
                        ib1=int(b1*2.-5.)
                        B0=B0CNEW/b0b1(ib1)
                        endif
                GOTO 924
c
c omit F1 feature ....................................................
c
7398  WRITE(KONSOL,9269)
9269  FORMAT(1X,'CORR.: NO F1 REGION, B1=3, C1=0.0')
        HMF1=0.
        NMF1=0.
        C1=0.0
        B1=3.
        F1REG=.FALSE.

C
C SEARCH FOR HST [NE3(HST)=NME] ..........................................
C
380     RRRR=0.5
        IF(F1REG) then
                hf1=hmf1
                xf1=nmf1
                GOTO 3972
                ENDIF
        RATHH=0.5
3973    hf1=hef+(hmf2-hef)*RATHH
        xf1=xe3(hf1)
        IF(XF1.LT.NME) THEN
                RATHH=RATHH+.1
                GOTO 3973
                ENDIF
3972    h=hf1
        deh=10.
        XXMIN=XF1
        HHMIN=HF1
3895    h=h-deh
        if(h.lt.HEF) then
          h=h+2*deh
          deh=deh/10.
          if(deh.lt.1.) goto 3885
          endif
        XE3H=XE3(h)
        IF(XE3H.LT.XXMIN) then
          XXMIN=XE3H
          HHMIN=h
          endif
        if(XE3H.gt.NME) goto 3895
      CALL REGFA1(h,HF1,XE3H,XF1,0.001,NME,XE3,SCHALT,HST)
        STR=HST
        IF(.not.SCHALT) GOTO 360
3885    WRITE(KONSOL,100)
100   FORMAT(1X,'*NE* HST IS NOT EVALUATED BY THE FUNCTION XE3')
        IREGFA=2
        IF(XXMIN/NME.LT.1.3) GOTO 9244
c
c assume linear interpolation between HZ and HEF ..................
c
8922    HZ=HHMIN+(HF1-HHMIN)*RRRR
        XNEHZ=XE3(HZ)
        if(xnehz-nme.lt.0.001) then
          RRRR=RRRR+.1
          GOTO 8922
          endif
        WRITE(KONSOL,901) HZ,HEF
901   FORMAT(6X,'CORR.: LIN. APP. BETWEEN HZ=',F5.1,
     &          ' AND HEF=',F5.1)
        T=(XNEHZ-NME)/(HZ-HEF)
        HST=-333.
        GOTO 4933
c
c calculate HZ, D and T ............................................
c
360     HZ=(HST+HF1)/2.0
        D=HZ-HST
        T=D*D/(HZ-HEF-D)
        GOTO 4933
C
C LAY-functions for middle ionosphere
C
6153    HMF1M=165.+0.6428*XHI
        HHALF = GRAT * HMF2
        HV1R = HME + WIDTH
        HV2R = HME + HDEEP
        HHMF2 = HMF2
        CALL INILAY(NIGHT,NMF2,NMF1,NME,VNER,HHMF2,HMF1M,HME,
     &                  HV1R,HV2R,HHALF,HXL,SCL,AMP,IIQU)
        IF(IIQU.EQ.1) WRITE(KONSOL,7733)
7733    FORMAT('*NE* LAY amplitudes found with 2nd choice of HXL(1).')
        IF(IIQU.EQ.2) WRITE(KONSOL,7722)
7722    FORMAT('*NE* LAY amplitudes could not be found.')

C---------- CALCULATION OF NEUTRAL TEMPERATURE PARAMETER-------

4933  HTA=120.0
      HTE=3000.0
        IF(NOTEM) GOTO 240
      SEC=UT*3600.
      CALL CIRA86(DAYNR,SEC,LATI,LONGI,HOUR,COV,TEXOS,TN120,SIGMA)
        IF(HOUR.NE.0.0) THEN
      SECNI=(24.-LONGI/15.)*3600.
      CALL CIRA86(DAYNR,SECNI,LATI,LONGI,0.,COV,TEXNI,TN1NI,SIGNI)
        ELSE
      TEXNI=TEXOS
      TN1NI=TN120
      SIGNI=SIGMA
        ENDIF
      TLBDH=TEXOS-TN120
      TLBDN=TEXNI-TN1NI
C
C--------- CALCULATION OF ELECTRON TEMPERATURE PARAMETER--------
C
881   CONTINUE

C !!!!!!!!!! TE(120KM)=TN(120KM) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ATE(1)=TN120

C !!!!!!!!!! TE-MAXIMUM (JICAMARCA,ARECIBO) !!!!!!!!!!!!!!!!!!!!
      HMAXD=60.*EXP(-(MLAT/22.41)**2)+210.
      HMAXN=150.
      AHH(2)=HPOL(HOUR,HMAXD,HMAXN,SAX,SUX,1.,1.)
      TMAXD=800.*EXP(-(MLAT/33.)**2)+1500.
      TMAXN=TN(HMAXN,TEXNI,TLBDN,SIGNI)+20
      ATE(2)=HPOL(HOUR,TMAXD,TMAXN,SAX,SUX,1.,1.)

C !!!!!!!!!! TE(300,400KM)=TE-AE-C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! TE(1400,3000KM)=TE-ISIS !!!!!!!!!!!!!!!!!!!!!!!!!!!
        DIPLAT=MAGBR
      CALL TEBA(DIPLAT,HOUR,NSESON,TEA)
      ATE(3)=TEA(1)
      ATE(4)=TEA(2)
      ATE(6)=TEA(3)
      ATE(7)=TEA(4)

C !!!!!!!!!! TE(600KM)=TE-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ETT=EXP(-MLAT/11.35)
      TET=2900.-5600.*ETT/((ETT+1)**2.)
      TEN=839.+1161./(1.+EXP(-(ABSMLT-45.)/5.))
      ATE(5)=HPOL(HOUR,TET,TEN,SAX,SUX,1.5,1.5)

C !!!!!!!!!! OPTION TO USE TE-NE-RELATION !!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! AT 300, 400 OR 600 KM  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(TENEOP) THEN
        DO 3395 I=1,3
3395      IF(TCON(I)) ATE(I+2)=TEDE(HOA(I),XNAR(I),-COV)
        ENDIF
C !!!!!!!!!! TE'S ARE CORRECTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! ALSO TE > TN ENFORCED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TNAHH2=TN(AHH(2),TEXOS,TLBDH,SIGMA)
      IF(ATE(2).LT.TNAHH2) ATE(2)=TNAHH2
      STTE1=(ATE(2)-ATE(1))/(AHH(2)-AHH(1))
      DO 1901 I=2,6
       TNAHHI=TN(AHH(I+1),TEXOS,TLBDH,SIGMA)
       IF(ATE(I+1).LT.TNAHHI) ATE(I+1)=TNAHHI
       STTE2=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
       ATE(I)=ATE(I)-(STTE2-STTE1)*DTE(I-1)*ALOG2
1901  STTE1=STTE2
C !!!!!!!!!! GRADIENTS ARE CALCULATED WITH !!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! CORRECTED REGION BOUNDARIES !!!!!!!!!!!!!!!!!!!!!!
      DO 1902 I=1,6
1902  STTE(I)=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
      ATE1=ATE(1)
887   CONTINUE
C
C------------ CALCULATION OF ION TEMPERATURE PARAMETERS--------
C
C !!!!!!!!!! TI(430KM,DAY)=TI-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!!!
      XSM1=430.0
      XSM(1)=XSM1
      Z1=EXP(-0.09*MLAT)
      Z2=Z1+1.
      TID1 = 1240.0 - 1400.0 * Z1 / ( Z2 * Z2 )
      MM(2)=HPOL(HOUR,3.0,0.0,SAX,SUX,1.,1.)
C !!!!!!!!!!  TI < TE   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TED1=TEA(6)+30.
        IF(TID1.GT.TED1) TID1=TED1

C !!!!!!!!!! TI(430KM,NIGHT)=TI-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!
      Z1=ABSMLT
      Z2=Z1*(0.47+Z1*0.024)*UMR
      Z3=COS(Z2)
      TIN1=1200.0-300.0*SIGN(1.0,Z3)*SQRT(ABS(Z3))
C !!!!!!!!!! TN < TI < TE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TEN1=TEA(5)
        TNN1=TN(XSM1,TEXNI,TLBDN,SIGNI)
        IF(TEN1.LT.TNN1) TEN1=TNN1
        IF(TIN1.GT.TEN1) TIN1=TEN1
        IF(TIN1.LT.TNN1) TIN1=TNN1

C !!!!!!!!!! TI(430KM,LT) FROM STEP FUNCTION !!!!!!!!!!!!!!!!!!
        TI1=TIN1  
        IF(TID1.GT.TIN1) TI1=HPOL(HOUR,TID1,TIN1,SAX,SUX,1.,1.)

C !!!!!!!!!! TANGENT ON TN DETERMINES HS !!!!!!!!!!!!!!!!!!!!!!
        TI13=TEDER(130.)
        TI50=TEDER(500.)
      CALL REGFA1(130.0,500.0,TI13,TI50,0.01,TI1,TEDER,SCHALT,HS)
      IF(SCHALT) HS=200.
      TNHS=TN(HS,TEXOS,TLBDH,SIGMA)
      MM(1)=DTNDH(HS,TEXOS,TLBDH,SIGMA)
      IF(SCHALT) MM(1)=(TI1-TNHS)/(XSM1-HS)
      MXSM=2

C !!!!!!!!!! XTETI ALTITTUDE WHERE TE=TI !!!!!!!!!!!!!!!!!!!!!!
2391    XTTS=500.
        X=500.
2390    X=X+XTTS
        IF(X.GE.AHH(7)) GOTO 240
        TEX=ELTE(X)
        TIX=TI(X)
        IF(TIX.LT.TEX) GOTO 2390
        X=X-XTTS
        XTTS=XTTS/10.
        IF(XTTS.GT.0.1) GOTO 2390
        XTETI=X+XTTS*5.

C !!!!!!!!!! TI=TE ABOVE XTETI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MXSM=3
        MM(3)=STTE(6)
        XSM(2)=XTETI
        IF(XTETI.GT.AHH(6)) GOTO 240
        MXSM=4
        MM(3)=STTE(5)
        MM(4)=STTE(6)
        XSM(3)=AHH(6)
        IF(XTETI.GT.AHH(5)) GOTO 240
        MXSM=5
        DTI(1)=5.
        DTI(2)=5.
        MM(3)=STTE(4)
        MM(4)=STTE(5)
        MM(5)=STTE(6)
        XSM(3)=AHH(5)
        XSM(4)=AHH(6)
C
C CALCULATION OF ION DENSITY PARAMETER..................
C
240   IF(NOION) GOTO 141
      HNIA=100.
      HNIE=2000.
        if(DY) goto 141
C
C INPUT OF THE ION DENSITY PARAMETER ARRAYS PF1O,PF2O AND PF3O......
C
      RIF(1)=2.
      IF(ABSLAT.LT.30.0) RIF(1)=1.
      RIF(2)=2.
      IF(COV.LT.100.0) RIF(2)=1.
      RIF(3)=SEASON
      IF(SEASON.EQ.1) RIF(3)=3.
      RIF(4)=1.
      IF(NIGHT) RIF(4)=2.
      CALL KOEFP1(PG1O)
      CALL KOEFP2(PG2O)
      CALL KOEFP3(PG3O)
      CALL SUFE(PG1O,RIF,12,PF1O)
      CALL SUFE(PG2O,RIF, 4,PF2O)
      CALL SUFE(PG3O,RIF,12,PF3O)
c
c calculate O+ profile parameters
c
      IF(ABS(XHI).LE.90.0) THEN
        ZZZ1=COS(XHI*UMR)
      ELSE
        ZZZ1=0.0
      ENDIF
      msumo=4
      RDOMAX=100.0
      MO(1)=EPSTEP(PF1O(1),PF1O(2),PF1O(3),PF1O(4),ZZZ1)
      MO(2)=EPSTEP(PF1O(5),PF1O(6),PF1O(7),PF1O(8),ZZZ1)
      MO(3)=0.0
      HO(1)=EPSTEP(PF1O(9),PF1O(10),PF1O(11),PF1O(12),ZZZ1)
      HO(2)=290.0
      IF((RIF(2).EQ.2.).AND.(RIF(3).EQ.2.)) HO(2)=237.0
      HO(4)=PF2O(1)
        ho05=pf2o(4)
      MO(4)=PF2O(2)
      MO(5)=PF2O(3)
c
c adjust gradient MO(4) of O+ profile segment above F peak
c
7100    HO(3)=(ALG100-MO(5)*(HO(4)-ho05))/MO(4)+HO(4)
        IF(HO(3).LE.HO(2)+20.) THEN
                MO(4)=MO(4)-0.001
                GOTO 7100
                endif
        hfixo=(ho(2)+ho(3))/2.
c
c find height H0O of maximum O+ relative density
c
      DELX=5.0
      X=HO(2)
      YMAXX=0.0
7102  X=X+DELX
      Y=RPID(X,HFIXO,RDOMAX,msumo,MO,DDO,HO)
      IF(Y.LE.YMAXX) then
        if(delx.le.0.1) GOTO 7104
        x=x-delx
        delx=delx/5.
      ELSE
        YMAXX=Y
      ENDIF
      GOTO 7102 
7104  H0O=X-DELX/2.
7101    if(y.lt.100.0) goto 7103
          rdomax=rdomax-0.01
        y=rpid(h0o,hfixo,rdomax,msumo,mo,ddo,ho)
        goto 7101
7103    yo2h0o=100.-y
        yoh0o=y
c
c calculate parameters for O2+ profile
c
        hfixo2  = pf3o(1)
        rdo2mx = pf3o(2) 
      DO 7105 L=1,2
                I = L * 2
                HO2(L)=PF3O(1+I)+PF3O(2+I)*ZZZ1
7105            MO2(L+1)=PF3O(7+I)+PF3O(8+I)*ZZZ1
      MO2(1)=PF3O(7)+PF3O(8)*ZZZ1
        if(hfixo2.gt.ho2(1)) then
           ymo2z=mo2(2)
        else
           ymo2z=mo2(1)
        endif
        aldo21=alog(rdo2mx)+ymo2z*(ho2(1)-hfixo2)
        hfixo2=(ho2(2)+ho2(1))/2.
        rdo2mx=exp(aldo21+mo2(2)*(hfixo2-ho2(1)))
c
c make sure that rd(O2+) is less or equal 100-rd(O+) at O+ maximum
c
7106  Y=RPID(H0O,hfixo2,rdo2mx,2,MO2,DO2,HO2)
      IF(Y.GT.yo2h0o) then
        MO2(3)=MO2(3)-0.02
        GOTO 7106
        endif
c
C use ratio of NO+ to O2+ density at O+ maximum to calculate
c NO+ density above the O+ maximum (H0O)
c
      IF(y.LT.1.) then
        NOBO2=0.0
      ELSE
        NOBO2= (yo2h0o-y)/y
      ENDIF
C
C CALCULATION FOR THE REQUIRED HEIGHT RANGE.......................
C
141   IF(.NOT.F1REG) HMF1=HZ

      DO 7397 KI=1,11
      DO 7397 KK=1,50
7397    OUTF(KI,KK)=-1.

      HEIGHT=HEIBEG
      KK=1

300   IF(NODEN) GOTO 330
      IF((HEIGHT.GT.HNEE).OR.(HEIGHT.LT.HNEA)) GOTO 330
        IF(LAYVER) THEN
          ELEDE=-9.
          IF(IIQU.LT.2) ELEDE=XEN(HEIGHT,HMF2,NMF2,HME,4,HXL,SCL,AMP)
        ELSE
          ELEDE=XE(HEIGHT)
        ENDIF
        OUTF(1,KK)=ELEDE
330   IF(NOTEM) GOTO 7108
      IF((HEIGHT.GT.HTE).OR.(HEIGHT.LT.HTA)) GOTO 7108
        TNH=TN(HEIGHT,TEXOS,TLBDH,SIGMA)
        TIH=TNH
        IF(HEIGHT.GE.HS) TIH=TI(HEIGHT)
        TEH=ELTE(HEIGHT)
        IF(TIH.LT.TNH) TIH=TNH
        IF(TEH.LT.TIH) TEH=TIH
        OUTF(2,KK)=TNH
        OUTF(3,KK)=TIH
        OUTF(4,KK)=TEH
7108  IF(NOION) GOTO 7118
      IF((HEIGHT.GT.HNIE).OR.(HEIGHT.LT.HNIA)) GOTO 7118
        if(DY) then
      call IONCOM(HEIGHT,XHI*UMR,LATI*UMR,COV,MONTH,DION)
      ROX=DION(1)
      RHX=DION(2)
      RNX=DION(3)
      RHEX=DION(4)
      RNOX=DION(5)
      RO2X=DION(6)
      RCLUST=DION(7)
        else
      ROX=RPID(HEIGHT,HFIXO,RDOMAX,msumo,MO,DDO,HO)
      RO2X=RPID(HEIGHT,HFIXO2,rdo2mx,2,MO2,DO2,HO2)
      CALL RDHHE(HEIGHT,H0O,ROX,RO2X,NOBO2,10.,RHX,RHEX)
      RNOX=RDNO(HEIGHT,H0O,RO2X,ROX,NOBO2)
      RNX=-1.
      RCLUST=-1.
        endif
      OUTF(5,KK)=ROX
      OUTF(6,KK)=RHX
      OUTF(7,KK)=RHEX
      OUTF(8,KK)=RO2X
      OUTF(9,KK)=RNOX
      OUTF(10,KK)=RCLUST
      OUTF(11,KK)=RNX

7118    HEIGHT=HEIGHT+HEISTP
        KK=KK+1
        IF(KK.LE.NUMHEI) GOTO 300
C
C ADDITIONAL PARAMETER FIELD OARR
C
        IF(NODEN) GOTO 6192
      OARR(1)=NMF2
      OARR(2)=HMF2
      OARR(3)=NMF1
      OARR(4)=HMF1
      OARR(5)=NME
      OARR(6)=HME
      OARR(7)=NMD
      OARR(8)=HMD
      OARR(9)=HHALF
      OARR(10)=B0
      OARR(11)=VNER
      OARR(12)=HEF
6192    IF(NOTEM) GOTO 6092
      OARR(13)=ATE(2)
      OARR(14)=AHH(2)
      OARR(15)=ATE(3)
      OARR(16)=ATE(4)
      OARR(17)=ATE(5)
      OARR(18)=ATE(6)
      OARR(19)=ATE(7)
      OARR(20)=ATE(1)
      OARR(21)=TI1
      OARR(22)=XTETI
6092  OARR(23)=XHI
      OARR(24)=SUNDEC
      OARR(25)=DIP
      OARR(26)=MAGBR
      OARR(27)=MODIP

3330  CONTINUE
      RETURN
      END
