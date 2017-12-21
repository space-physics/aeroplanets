      Subroutine CaltoJul_V05(iY,iM,iD,ihour,imin,sec,xJD)
!f2py intent(out) xJD
C     Compute Julian day (Real*8) by method of Meeus, Astronomical
C       Algorithms, 2nd Edition, 1998, page 61. Inputs are year iY,
C       month iM, day of month iD, and time of day in hours, minutes,
C       and seconds (all integer except seconds).  Output is Real*8
C       Julian day, xJD.
C
      Implicit None
      Double Precision sec,xJD,D
      Integer iY,iM,iD,ihour,imin,Y,M,A,B
      Y = iY
      M = iM
C...  Consider Jan or Feb as if months 13 and 14 of previous year
      If (iM.le.2)Then
        Y = iY - 1
        M = iM + 12
      Endif
C...  Compute day of month plus fractional part
      D = iD + ihour/2.4d1 + imin/1.440d3 + sec/8.64D4
      A = IDint(Y/100.0d0)
      B = 2 - A + IDint(A/4.0d0)
C...  Compute Julian day with fractional part
      xJD = IDint(365.25d0*(Y+4716)) + IDint(30.6001d0*(M+1))
     &  + D + B - 1524.5d0
      Return
      End
C----------------------------------------------------------------------
      Subroutine Datastep_V05(I,CHGT,CLAT,CLON,CSEC,DAY0,RHOd,RHOu,
     & RHOv,EOF,DELHGT,DELLAT,DELLON,DELTIME,TEMP,PRES,DENSLO,DENS,
     & DENSHI,DENSP,EWWIND,EWpert,NSWIND,NSpert,Hrho,HSCALE,dsunlat,
     & dsunlon,dsunLs,dradau,dowlt,LonEast,corlim,DENSTOT,IERT,IUTC,
     & fmolCO2,fmolN2,fmolO,fmolCO,fmolHe,fmolN,fmolH,AMz,pertstep,
     & corlmin,iupdate,ALS,SZA,owlt,sunlat,sunlon,VenusAU,TLOCAL,
     & profnear,proffar,nprof, DATADIR)
!f2py intent(inout) iupdate
!f2py intent(out) TEMP,PRES,DENSLO,DENS,DENSHI,DENSP,EWWIND,EWpert,NSWIND,NSpert,Hrho,HSCALE
!f2py intent(out) corlim,DENSTOT,fmolCO2,fmolN2,fmolO,fmolCO,fmolHe,fmolN,fmolH,AMz,ALS,SZA,owlt,sunlat,sunlon,VenusAU,TLOCAL
      ! modif GG include DATADIR
C...  Subroutine to compute atmosphere at next trajectory or profile
C     position, including perturbations, with appropriate correlation
C     maintained over step size
      Implicit None
      Double Precision CHGT,CLAT,CLON,DELHGT,DELLAT,DELLON,DTR,Rref,
     &  NSWIND,NSpert,NStot,TEMP,PRES,DENSLO,DENS,DENSHI,DENSP,EWWIND,
     &  EWpert,Hrho,HSCALE,corlim,DENSTOT,pi,gz,yo,yco,yn,amo,amco,amn,
     &  ysum,yh,amh,yhe,amhe,yco2,amco2,devtot,ytot,amz,presmb,
     &  patsurf,sigmalevel,devhi,varx,vary,pvra,alogdens,var,sza,
     &  yhp,yhep,yco2p,devav,preshgt,olon,ppnd_V05,random_V05,ogz,
     &  dvra,devlo,tvra,rhou,rhov,ewtot,yop,ycop,ynp,
     &  ssp2,fc,fa,fb,fact,wlon,z1,rhod,densrand,rpscale,sigd,sigu,
     &  sigv,trajlat,trajlon,correl,facthi,dplus,factlo,dminus,z2,
     &  rsc,vls,hls,delns,delew,delz,trajhgt,fmolh,fmolhe,
     &  pertstep,corlmin,corbeta,yn2,amn2,yn2p,profnear,proffar,
     &  CSEC,DAY0,DELTIME,dsunlat,dsunlon,dsunLs,dradau,dowlt,profwgt,
     &  TrajSEC,DAY,VenusAU,owlt,sunlat,sunlon,ALS,EOT,ttsec,dt,
     &  TLOCAL,DayVenus,zeta,fmolCO2,fmolN2,fmolO,fmolCO,fmolN,
     &  zlo,plo,dlo,tlo,zetalo,Rlo,zmd,pmd,dmd,tmd,yco2md,yn2md,yomd,
     &  ycomd,yhemd,ynmd,zhi,phi,dhi,thi,yco2hi,yn2hi,yohi,
     &  ycohi,yhehi,ynhi,yhhi
      INTEGER EOF,LonEast,logscale,nvarx,nvary,iu0,iup,maxfiles,L,
     &  ifault,i,npos,iupdate,IERT,IUTC,nprof
c      Character*60 lstfl,outfl
      Character*300 DATADIR !+modif GG +300!!!!
      Character*8 densunits
C...  Venus model VIRA data at low ,middle, and high altitudes
      Common /VIRAdata_V05/zlo(81),plo(81,5),dlo(81,5),tlo(81,5),
     &  zetalo(81,5),Rlo(81,5),zmd(11),pmd(11,2),dmd(11,2),tmd(11,2),
     &  yco2md(11,2),yn2md(11,2),yomd(11,2),ycomd(11,2),yhemd(11,2),
     &  ynmd(11,2),zhi(21),phi(21,7),dhi(21,7),
     &  thi(21,7),yco2hi(21,7),yn2hi(21,7),yohi(21,7),ycohi(21,7),
     &  yhehi(21,7),ynhi(21,7),yhhi(21,7)
      COMMON /DATACOM_V05/rpscale,NPOS,NVARX,NVARY,logscale,iu0,iup,
     & maxfiles
c      COMMON /FILENAME_V05/lstfl,outfl
C...  Molecular weights for atmospheric constituents
      Data amh,amco2,amhe,amn2,amo,amco,amn/1.01d0,44.0d0,4.0d0,28.0d0,
     &  16.0d0,28.0d0,14.0d0/
      z1 =0. ! gg avoids warning about uninitialization
      Call ReadVIRA_V05(DATADIR) ! Modif GG
      DTR = DAtan(1.0d0)/45.0d0
      pi = 4.0d0*dAtan(1.0d0)
      EOF = 0
C...  Venus day length, seconds
      DayVenus = 1.00872d7
C...  Get radius RSC to current position
      Call rgplanet_V05(CLAT,Rref,CHGT,gz,2)
      RSC = Rref + CHGT
C...  Set vertical and horizontal scale parameters
C     Vertical scale approximates VLS/Hrho = 2
      VLS = 32.0d0  - 22.0d0*CHGT/65.0d0
      If (VLS.lt.10.0d0)VLS = 10.0d0
      If (CHGT.gt.150.0d0)VLS = 10.0d0 + 0.6d0*(CHGT-150.0d0)
C...  Horizontal scale selected to be fairly consistent with wavelength
C     estimates from Fig. 6 of Bougher and Borucki, J. Geophys. Res.
C     99(E2), 3759-3776 (1994), Fig. 2 of Mayr et al., J. Geophys. Res.
C     93(A10), 11247-262 (1988), Kasprzak et al. J. Geophys. Res.
C     93(A10), 11237-246 (1988), and Kasprzak et al. Geophys. Res.
C     Lett. 20 2755-2758 (1993)
      HLS = 60.0d0*VLS
C...  Relative displacements between previous and current position
      If (I.ne.0)Then
        DELNS = DTR*RSC*(DELLAT)/HLS
        DELEW = -DTR*RSC*dCOS(DTR*CLAT)*DELLON/HLS
        DELZ = DELHGT/VLS
      Else
        DELNS = 0.0
        DELEW = 0.0
        DELZ = 0.0
      Endif
      NPOS = 1
      IF(NPOS.le.0)then
C...    Read new position if trajectory data file is being used
        READ (7,*,ERR=9998,END=999)TrajSEC,TrajHGT,TrajLAT,TrajLON
C...    Convert negative longitudes
        If (TrajLON.lt.0.0d0)TrajLON = TrajLON + 360.0d0
C...    Convert to East Longitude if LonEast = 0
        If (LonEast.eq.0)TrajLON = 360.0d0 - TrajLON
        If (I.gt.0)Then
C...      Compute displacement magnitude of new from previous position
          DELHGT = TrajHGT - CHGT
          DELLAT = dAbs(TrajLAT - CLAT)
          DELLON = dAbs(TrajLON - CLON)
          DELTIME = TrajSEC - CSEC
        Endif
C...    Correct DELLON for cases near 0/360 longitude discontinuity
        If (DELLON.gt.180.0d0)DELLON = 360.0d0 - DELLON
        If (DELLON.lt.0.0d0)DELLON = DELLON + 360.0d0
C...    Correct DELLON and DELLAT near polar discontinuities
        If (DELLON.gt.90.0d0.and.(dAbs(TrajLAT).ge.70.0d0.or.
     &   dAbs(CLAT).ge.70.0d0))Then
          DELLON = dAbs(180.0d0 - DELLON)
          DELLAT = 180.0d0 - dAbs(TrajLAT) - dAbs(CLAT)
        Endif
C...    Relative displacements between previous and current position
        DELNS = DTR*RSC*(DELLAT)/HLS
        DELEW = -DTR*RSC*dCOS(DTR*CLAT)*DELLON/HLS
        DELZ = DELHGT/VLS
C...    Set current position to new position
        CHGT = TrajHGT
        CLAT = TrajLAT
        CLON = TrajLON
        CSEC = TrajSEC
      Else If (I.gt.0) then
        CHGT = CHGT + DELHGT
        CLAT = CLAT + DELLAT
        CLON = CLON + DELLON
        CSEC = CSEC + DELTIME
      Endif
C...  Update Julian day
      DAY = DAY0 + CSEC/8.64d4
C...  Get radius RSC to new position
      Call rgplanet_V05(CLAT,Rref,CHGT,ogz,2)
      RSC = Rref + CHGT
C...  Correct latitude and longitude if position crosses either pole
      IF(dABS(CLAT).GT.90.0d0)then
        CLAT = DSIGN(180.0d0,CLAT) - CLAT
        CLON = CLON + 180.0d0
        DELLAT = -DELLAT
      Endif
      IF(CLON.LT.0.0d0)CLON = CLON + 360.0d0
      IF(CLON.GE.360.0d0) CLON = CLON - 360.0d0
C...  Sun and Venus positions at new time
C...  Use input ephemeris values, if supplied
      If (dradau.gt.0.0)Then
        SUNLON = dsunlon
        SUNLAT = dsunlat
        ALS = dsunLs
        VenusAU = dradau
        owlt = dowlt
      Else
C...    Use built-in Venus ephemeris routine
C...    Convert to Terrestrial (Dynamical) Time, if necessary
        ttsec = 0.0
        If (iutc.eq.1)Then
C...      Get terrestrial dynamical time offset (seconds)
          dt = (DAY - 2451545.0d0)/36525.0d0
C...      Terrestrial time offset (in seconds) TT = UTC + ttsec
          ttsec = (64.184d0 + 95.0d0*dt + 35.0d0*dt**2)/86400.0d0
        Endif
        Call venephem_V05(DAY+ttsec,sunlat,sunlon,ALS,VenusAU,owlt,EOT)
C...    Convert to Venus-Event Time, if necessary
        If(iert.eq.1)Call venephem_V05(DAY+ttsec-owlt/1440.0d0,sunlat,
     &    sunlon,ALS,VenusAU,owlt,EOT)
      Endif
C...  Local true solar time in Venus hours (1/24th Venus day)
      TLOCAL = 12.0d0 + (SUNLON - CLON)/15.0d0
      IF (TLOCAL .LT. 0.)TLOCAL = TLOCAL + 24.
      IF (TLOCAL .GT. 24.)TLOCAL = TLOCAL - 24.
C...  Get solar zenith angle
      Call  SolZenAng_V05(sunlat,sunlon,clat,clon,sza)
C...  Evaluate Venus engineering model atmospheric parameters
      Call Venusatm_V05(CHGT,CLAT,CLON,TLOCAL,sza,HSCALE,TEMP,DENS,
     &  FACTHI,FACTLO,PRES,Hrho,AMz,EWWIND,NSWIND,ytot,yco2,yn2,yo,yco,
     &  yhe,yn,yh,profnear,proffar,nprof,profwgt)
C...  Evaluate correlation and step size relative to accuracy limit
      If (iupdate.ge.0)pertstep=pertstep+dAbs(DELNS) + dAbs(DELEW)
     &  + dAbs(DELZ)+dSqrt(EWWIND**2+NSWIND**2)*DELTIME/(1000.0d0*HLS)
      corlim = -pertstep/dlog(0.995d0)
      If (corlim.le.corlmin.or.iupdate.lt.0)Then
        CORREL = 1.0d0
        corbeta = 0.0
        If (iupdate.lt.0)Then
          iupdate = -1
        Else
          iupdate = 0
        Endif
      Else
C...    Get uniform and Gaussian random numbers from PPND_V05 and
C       RANDOM_V05
  480   Z2 = RANDOM_V05(L)
        IF (L .EQ. 1)GOTO 480
        Z1 = PPND_V05(Z2,IFAULT)
        IF (IFAULT .EQ. 1)STOP ' PPND ERROR'
        CORREL = dExp(-pertstep)
        corbeta = DSqrt(1.0d0 - CORREL**2)
        pertstep = 0.0
        iupdate = 1
      Endif
C...  Compute Density plus ~ 1 sigma
      DENSHI = DENS*FACTHI
      DPLUS = DENSHI - DENS
C...  Compute Density minus ~ 1 sigma
      DENSLO = DENS*FACTLO
      DMINUS = DENS - DENSLO
C...  Current random density perturbation value, correlated with
C...  previous random density perturbation
      RHOd = CORREL*RHOd + corbeta*Z1
      IF(RHOd.LT.0.0d0)DensRand = RHOd*DMINUS*rpscale
      IF(RHOd.GE.0.0d0)DensRand = RHOd*DPLUS*rpscale
C...  Add random density perturbation
      DENSP = DENS + DensRand
C...  Check upper and lower bounds on density perturbations
      If (DENSP .lt. 0.1d0*DENS)DENSP = 0.1d0*DENS
      If (DENSP .gt. 10.0d0*DENS)DENSP = 10.0d0*DENS
C...  Save as total density, for output
      DENSTOT = DENSP
C...  Standard deviation in random density perturbation (% of
C...  unperturbed mean) for output
      SIGD = rpscale*50.0d0*dAbs(DENSHI-DENSLO)/DENS
C...  Adjust random DENSHI, DENSLO for rpscale
      DENSHI = DENS + rpscale*(DENSHI - DENS)
      DENSLO = DENS + rpscale*(DENSLO - DENS)
      If (DENSLO .lt. 0.1d0*DENS)DENSLO = 0.1d0*DENS
C...  Convert random density perturbation to % of (unperturbed) mean
      DensRand = 100.0d0*(DENSP - DENS)/DENS
C...  Convert DENSP to density perturbation % of (unpertubed) mean
      DENSP = DensRand
C...  Standard deviation for wind perturbations, from approximation to
C     VIRA model
      If (CHGT.le.30.0d0)Then
        SIGU = 1.0d0 + 7.0d0*CHGT/30.0d0
      Else If (CHGT.le.45.0d0)Then
        SIGU = 8.0d0
      Else If (CHGT.le.60.0d0)Then
        SIGU = 8.0d0 + 7.0d0*(CHGT-45.0d0)/15.0d0
      Else If (CHGT.le.160.0d0)Then
        SIGU = 9.0d0 + CHGT/10.0d0
      Else
        SIGU = 25.0d0
      Endif
      SIGU = SIGU*rpscale
      SIGV = SIGU/2.0d0
C...  Compute EW and NS wind perturbations and total wind
 586  If (corbeta.ne.0.0)Then
        Z2 = RANDOM_V05(L)
        If (L.eq.1)Goto 586
        Z1 = PPND_V05(Z2,ifault)
      Endif
      RHOu = CORREL*RHOu + corbeta*Z1
C...  EW component of perturbation in wind and total wind
      EWpert = RHOu*SIGU
      EWtot = EWWIND + EWpert
 587  If (corbeta.ne.0.0)Then
        Z2 = RANDOM_V05(L)
        If (L.eq.1)Goto 587
        Z1 = PPND_V05(Z2,ifault)
      Endif
      RHOv = CORREL*RHOv + corbeta*Z1
C...  NS component of perturbation in wind and total wind
      NSpert = RHOv*SIGV
      NStot = NSWIND + NSpert
C...  Limit winds to sound speed
C     Assume specific heat ratio = 1.286 (CO2 - N2 mixture)
      ssp2 = 1.286d0*PRES/DENS
      fc = ssp2 - EWtot**2 - NStot**2
      If (fc.lt.0.0d0)Then
C...    Find multiplier factor for wind perturbations that keep total
C       wind speed <= speed of sound
        fc = ssp2 - EWWIND**2 - NSWIND**2
        fa = EWpert**2 + NSpert**2
        If (fa.le.0.0d0)fa = 1.0d0
        fb = EWpert*EWWIND + NSpert*NSWIND
        fact = (-fb + dSqrt(dAbs(fb**2 + fa*fc)))/fa
        If (fact.lt.0.0d0)fact = 0.0d0
        If (fact.gt.1.0d0)fact = 1.0d0
C...    Recompute perturbations and total winds with required factor
        EWpert = EWpert*fact
        NSpert = NSpert*fact
        EWtot = EWWIND + EWpert
        NStot = NSWIND + NSpert
      Endif
C...  Write descriptively formatted data on LIST.txt file
      wlon = 360.0d0 - CLON
      zeta = PRES*AMz/(DENS*8314.472*TEMP)
      If (CHGT.gt.500.)zeta = 1.0
      If(iup.gt.0)Write(iup,590)CSEC,CSEC/DayVenus,ALS,CHGT,SZA,
     & RSC,Rref,owlt,HSCALE,Hrho,zeta,CLAT,CLON,wlon,SUNLAT,
     & VenusAU,SUNLON,TLOCAL
  590 FORMAT(' Time (rel. to T0) =',F11.1,' sec. (',F7.3,' Venus',
     & ' Days)    Ls =',F6.1/' Height Above Reference Ellipsoid =',F9.3
     & ,' km ',9x,'SZA =',F7.2,' deg'/' Total Radius (Ref Radius) = ',
     & F9.3,' (',F9.3,') km   OWLT =',F7.2,' Min'/' Scale Heights:',
     & ' H(p) =',F6.2,'       H(rho) =',F6.2,' km    zeta =',F7.4   /
     & ' Latitude = ',F7.2,'  degrees       Longitude = ',F7.2,' E (',
     & F7.2,' W) deg.'/ ' Sun Latitude = ',F8.2,' deg.      Venus',
     & ' Orbital Radius  =',F7.4,' AU'/' Sun Longitude = ',F7.2,
     & ' deg.E     Local True Solar Time =',F6.2,' Venus hr')
C...  Compute percent deviations from Venus average values
      Call Venusref_V05(CHGT,tvra,pvra,dvra)
      If (dvra.le.0.0d0)Then
        devlo = -99.9d0
        devav = -99.9d0
        devhi = -99.9d0
        devtot = -99.9d0
      Else
        devlo = 100.0d0*(DENSLO-dvra)/dvra
        devav = 100.0d0*(DENS-dvra)/dvra
        devhi = 100.0d0*(DENSHI-dvra)/dvra
        devtot = 100.0d0*(DENSTOT-dvra)/dvra
      Endif
      densunits = 'kg/m**3 '
C...  Convert density units to kg/km**3 if logscale = 3
      If (logscale .eq.3)Then
        DENS = 1.0D9*DENS
        DENSLO = 1.0D9*DENSLO
        DENSHI = 1.0D9*DENSHI
        DENSTOT = 1.0D9*DENSTOT
        densunits = 'kg/km**3'
      Endif
C...  Get mole fractions (volume fractions) from number densities
      fmolCO2 = 100.0*yco2/ytot
      fmolN2 = 100.0*yn2/ytot
      fmolH = 100.0d0*yh/ytot
      fmolHe = 100.0d0*yhe/ytot
      fmolO = 100.0d0*yO/ytot
      fmolCO = 100.0d0*yCO/ytot
      fmolN = 100.0d0*yN/ytot
C...  Get mass fractions from number densities
      ysum = yh*amh + yhe*amhe + yco2*amco2 + yn2*amn2 + yO*amO +
     &  yCO*amCO + yN*amN
      yhp = 100.0d0*yh*amh/ysum
      yhep = 100.0d0*yhe*amhe/ysum
      yco2p = 100.0d0*yco2*amco2/ysum
      yn2p = 100.0*yn2*amn2/ysum
      yOp = 100.0*yO*amO/ysum
      yCOp = 100.0*yCO*amCO/ysum
      yNp = 100.0*yn*amn/ysum
C...  Write formatted output to list file
      If(iup.gt.0)Then
        Write(iup,600)TEMP,PRES,profwgt,DENSLO,DENS,DENSHI,densunits,
     &  devlo,devav,devhi,DENSTOT,densunits,DensRand,iupdate,EWWIND,
     &  EWpert,EWtot,NSWIND,NSpert,NStot
        Write(iup,640)yco2,yn2,yo,yco,yco2p,yn2p,yop,ycop,fmolco2,
     &    fmoln2,fmolo,fmolco
  640   Format(' CO2 =',1p,E10.3,' N2 =',E10.3,'  O =',E10.3,' CO =',
     &    E10.3,0p,' #/m**3'/F16.3,3F15.3,' % by mass'/F16.3,
     &    3F15.3,' % by volume')
        Write(iup,620)yhe,yn,yh,ytot,yhep,ynp,yhp,AMz,fmolhe,fmoln,
     &    fmolh
        If (I.gt.0.and.corlim.lt.1.0d0.and.dAbs(CLAT).lt.89.99d0)
     &    Write(iup,610)corlim
        Write(iup,650)
      Endif
  600 FORMAT(' Temperature = ',F7.1,' K',5X,' Pressure =',1p,E10.3,
     & ' N/m**2   profwgt =',F6.3,/ ' Density (Low, Avg., High)  =',
     & 3E12.3,1X,A8/ ' Departure from Venus Avg = ',0p,F10.1,' %',
     & 2(F10.1,' %')/' Tot.Dens. =',1p,E10.3,0p,1X,A8,
     & '    Dens.Pert. =',F7.2,' % of mean  iupdate=',I2,/
     & ' Eastward Wind  (Mean, Perturbed, Total) = ', 3F7.1,' m/s'/
     & ' Northward Wind (Mean, Perturbed, Total) = ',3F7.1,' m/s')
  610 Format(' Warning: Step size smaller than accuracy limit by a ',
     &  'factor of',F6.3)
  620 Format(1p,'  He =',E10.3,'  N =',E10.3,'  H =',E10.3,'   Total=',
     &   E10.3,' #/m**3',0p/F16.3,2F15.3,'   % by mass   ',
     &   '  MolWgt=',F6.3/F16.3,2F15.3,'   % volume (or mole)',
     &   ' fraction')
  650 FORMAT(' -----------------------------------------------------',
     & '----------------------')
C...  Compute pressure in millibars
      PRESmb = PRES/100.0d0
C...  Get pressure at surface, sigma = p/psurf, and pressure altitude
C     = -Ln(sigma)
      Patsurf = plo(1,1)
      sigmalevel = PRES/Patsurf
      preshgt = -dlog(sigmalevel)
C...  Select plot output variables from codes NVARX and NVARY
      If(NVARX.EQ.9)VARX = PRESmb
      If(NVARY.EQ.9)VARY = PRESmb
      If(NVARX.EQ.10)VARX = preshgt
      If(NVARY.EQ.10)VARY = preshgt
      If(NVARX.EQ.11)VARX = sigmalevel
      If(NVARY.EQ.11)VARY = sigmalevel
C...  Output deviations from Venus average if logscale = 2
      If (logscale.eq.2)Then
        DENSLO = devlo
        DENS = devav
        DENSHI = devhi
        DENSTOT = devtot
        If (pvra.le.0.0d0)Then
          PRES = -99.9d0
        Else
          PRES = 100.0d0*(PRES-pvra)/pvra
        Endif
      Endif
C...  Write parameters on plot format files
      IF(NVARX.EQ.1)VARX = CHGT
      IF(NVARX.EQ.2)VARX = RSC
      IF(NVARX.EQ.3)VARX = CLAT
      IF(NVARX.EQ.4)Then
        VARX = CLON
        If (LonEast.eq.0)VARX = 360.0d0 - CLON
      Endif
      If(NVARX.eq.13)Then
        VARX = CLON
        If (VARX.gt.180.0d0)VARX = VARX - 360.0d0
        If (LonEast.eq.0)VARX = -VARX
      Endif
      If (NVARX.eq.5)VARX = CSEC
      If (NVARX.eq.6)VARX = CSEC/DayVenus
      If (NVARX.eq.7)VARX = ALS
      If (NVARX.eq.8)VARX = TLOCAL
      If (NVARX.eq.12)VARX = SZA
      Alogdens = 0.0d0
C...  Change density and pressure units as determined by LOGSCALE
      If (logscale .ne. 2)Alogdens = dlog10(DENS)
      If (logscale .eq. 1)then
        DENS = Alogdens
        PRES = dlog10(PRES)
        DENSLO = dlog10(DENSLO)
        DENSHI = dLOG10(DENSHI)
        DENSTOT = dlog10(DENSTOT)
      Endif
C...  Write plotable output files
      If (NVARY .eq. 0.and.iup.gt.0)then
        Write(21,796)VARX,DENSLO,DENS,DENSHI,DENSTOT,RSC,Rref,
     &    ogz,LOGSCALE,profwgt
        Write(22,790)VARX,SIGD,DensRand,corlim,SIGU,SIGV,iupdate,
     &    TEMP*(1.0-SIGD/100.),TEMP,TEMP*(1.0+SIGD/100.),
     &    TEMP*(1.0-DensRand/100.)
        Write(23,791)VARX,EWWIND,EWpert,EWtot,NSWIND,NSpert,NStot,
     &    iupdate
        Write(24,798)VARX,TEMP,PRES,TEMP-273.15d0,PRESmb,Hrho,HSCALE,
     &    AMz,fmolco2,fmoln2,fmolo,fmolco,fmolhe,fmoln,fmolh,LOGSCALE
      Else If (iup.gt.0) Then
        IF(NVARY.EQ.1)VARY = CHGT
        IF(NVARY.EQ.2)VARY = RSC
        IF(NVARY.EQ.3)VARY = CLAT
        IF(NVARY.EQ.4)Then
          VARY = CLON
          If (LonEast.eq.0)VARY = 360.0d0 - CLON
        Endif
        If(NVARY.eq.13)Then
          VARY = CLON
          If (VARY.gt.180.0d0)VARY = VARY - 360.0d0
          If (LonEast.eq.0)VARY = -VARY
        Endif
        If (NVARY.eq.5)VARY = CSEC
        If (NVARY.eq.6)VARY = CSEC/DayVenus
        If (NVARY.eq.7)VARY = ALS
        If (NVARY.eq.8)VARY = TLOCAL
        If (NVARY.eq.12)VARY = SZA
        Write(21,797)VARX,VARY,DENSLO,DENS,DENSHI,DENSTOT,
     &    RSC,Rref,ogz,LOGSCALE,profwgt
        Write(22,795)VARX,VARY,SIGD,DensRand,corlim,SIGU,SIGV,iupdate,
     &    TEMP*(1.0-SIGD/100.),TEMP,TEMP*(1.0+SIGD/100.),
     &    TEMP*(1.0-DensRand/100.)
        Write(23,792)VARX,VARY,EWWIND,EWpert,EWtot,NSWIND,NSpert,NStot,
     &    iupdate
        Write(24,799)VARX,VARY,TEMP,PRES,TEMP-273.15d0,PRESmb,Hrho,
     &    HSCALE,AMz,fmolco2,fmoln2,fmolo,fmolco,fmolhe,fmoln,fmolh,
     &    LOGSCALE
      Endif
  790 FORMAT(G13.5,F6.2,F10.3,1p,E10.3,0p,2F7.2,I5,3x,4F7.2)
  791 Format(G13.5,6F8.2,I5)
  792 Format(2G13.5,6F8.2,I5)
  795 Format(2G13.5,F6.2,F10.3,1p,E10.3,0p,2F7.2,I5,3x,4F7.2)
  796 FORMAT(G13.5,1p,4E11.3,0p,2F9.2,F7.3,I5,F9.3)
  797 Format(2G13.5,1p,4E11.3,0p,2F9.2,F7.3,I5,F9.3)
  798 FORMAT(G13.5,2(F7.1,1p,E11.3,0p),2F7.2,F6.2,F6.1,5F5.1,F6.1,I5)
  799 Format(2G13.5,2(F7.1,1p,E11.3,0p),2F7.2,F6.2,F6.1,5F5.1,F6.1,I5)
C...  Write non-descriptively formatted data on OUTPUT file
      VAR = CHGT
      IF(NVARX.EQ.2.or.NVARY.eq.2)VAR = RSC
      If(iup.gt.0)Then
        OLON = CLON
        If (LonEast.eq.0)OLON = 360.0d0 - CLON
        If (NVARX.eq.13.or.NVARY.eq.13)Then
          If(OLON.gt.180.0d0)OLON = OLON - 360.0d0
        Endif
        If(logscale.eq.0.or.logscale.eq.3)Then
          WRITE(25,800)CSEC,VAR,CLAT,OLON,dens,TEMP,EWWIND,NSWIND,SIGD,
     &      ALS,SZA,yco2p,yn2p,yop,ycop,yhep,ynp,yhp
        Else If (logscale.eq.1)Then
          WRITE(25,810)CSEC,VAR,CLAT,OLON,dens,TEMP,EWWIND,NSWIND,SIGD,
     &      ALS,SZA,yco2p,yn2p,yop,ycop,yhep,ynp,yhp
        Else
          WRITE(25,805)CSEC,VAR,CLAT,OLON,dens,TEMP,EWWIND,NSWIND,SIGD,
     &      ALS,SZA,yco2p,yn2p,yop,ycop,yhep,ynp,yhp
        Endif
      Endif
  800 FORMAT(F10.1,F8.2,F7.2,F8.2,1P,E9.2,0P,3F7.1,5F6.1,4F5.1,F6.1)
  805 FORMAT(F10.1,F8.2,F7.2,F8.2,F9.2,3F7.1,5F6.1,4F5.1,F6.1)
  810 FORMAT(F10.1,F8.2,F7.2,F8.2,F9.3,3F7.1,5F6.1,4F5.1,F6.1)
      Return
  999 EOF = 1
      If (NPOS.LE.0)Rewind(7)
      Return
 9998 Stop ' Error termination reading trajectory data file!'
      END
C-----------------------------------------------------------------------
      Double Precision Function PPND_V05(p, ifault)
C
C     Algorithm AS 111 Appl. Statist. (1977) Vol. 26, p. 118
C
C     Produces normal deviate corresponding to lower tail area of p.
C     Returns ifault = 1 in input p >= 1 or <= 0, ifault = 0
C     otherwise.  If ifault = 1, PPND_V05 value is set to 0.
C     Double precision version with error epsilon = 2 ** (-31).
C     For Single precision version, change DOUBLE PRECISION to REAL
C     in the FUNCTION statement and the declaration of variables;
C     change E0 to D0 in the DATA statements and change dABS, dLOG
C     and dSQRT to ABS, ALOG and SQRT in the assignment statements.
C     The hash sums are the sums of the moduli of the coefficients.
C     They have no inherent meanings, but are included for use in
C     checking transpositions.
C
      Implicit None
      Double Precision zero, split, half, one, a0, a1, a2, a3,
     & b1, b2, b3, b4, c0, c1, c2, c3, d1, d2, p, q, r
      Integer ifault
C
      Data zero, half, one, split /0.0D0, 0.5D0, 1.0D0, 0.42D0/
C
      Data a0 /       2.50662823884D0/,
     & a1 /     -18.61500062529D0/,
     & a2 /      41.39119773534D0/,
     & a3 /     -25.44106049637D0/,
     & b1 /      -8.47351093090D0/,
     & b2 /      23.08336743743D0/,
     & b3 /     -21.06224101826D0/,
     & b4 /       3.13082909833D0/
C
C     Hash sum for a & b = 143.70383558076
C
      Data c0 /      -2.78718931138D0/,
     & c1 /      -2.29796479134D0/,
     & c2 /       4.85014127135D0/,
     & c3 /       2.32121276858D0/,
     & d1 /       3.54388924762D0/,
     & d2 /       1.63706781897D0/
C
C     Hash sum for c & d = 17.43746520924
C
C
      ifault = 0
      q = p - half
      If (dABS(q) .gt. split) goto 1
      r = q * q
      PPND_V05 = q * (((a3 * r + a2) * r + a1) * r + a0) /
     & ((((b4 * r + b3) * r + b2) * r + b1) * r + one)
      Return
    1 r = p
      If (q .gt.zero) r = one - p
      If (r .lt.zero) goto 2
      r = dSQRT(-dLOG(r))
      PPND_V05 = (((c3 * r + c2) * r + c1) * r + c0) /
     & ((d2 * r + d1) * r + one)
      If (q .lt. zero) PPND_V05 = -PPND_V05
      Return
    2 ifault = 1
      PPND_V05 = zero
      Return
      End
C-----------------------------------------------------------------------
      Double Precision Function RANDOM_V05(L)
C
C     Algorithm AS 183 Appl. Statist. (1982) Vol. 31, p.188
C
C     Returns a pseudo-random number rectangularly distributed
C     between 0 and 1.
C
C     IX, IY and IZ should be set to integer values between
C     1 and 30,000 before first entry.
C
C     Integer arithmetic up to 30323 is required.
C
C     Returns L = 0 unless random = 0 or random = 1, in which
C     case L = 1
C
      Implicit None
      Double Precision one,zero
      Integer IX,IY,IZ,L
      Common /RANDCOM_V05/ IX, IY, IZ
      Data one,zero/1.0D0,0.0D0/
C     IX = 171 * Mod(IX, 177) -  2 * (IX / 177)
C     IY = 172 * Mod(IY, 176) - 35 * (IY / 176)
C     IZ = 170 * Mod(IZ, 178) - 63 * (IZ / 178)
C
C     If (IX .lt. 0) IX = IX + 30269
C     If (IY .lt. 0) IY = IY + 30307
C     If (IZ .lt. 0) IZ = IZ + 30323
C
C     If integer arithmetic up to 5,212,632 is not available,
C     the preceding 6 statements may be used instead of the following 3
C
      IX = Mod(171 * IX, 30269)
      IY = Mod(172 * IY, 30307)
      IZ = Mod(170 * IZ, 30323)
C
C     On some machines, this may slightly decrease the speed.
C     The results should be identical.
C
      Random_V05 = dmod(dble(IX) / 30269.0d0 + dble(IY) / 30307.0d0 +
     & dble(IZ) / 30323.0d0, one)
      L = 0
      If (Random_V05 .le. zero .or. Random_V05 .ge. one)L = 1
      Return
      End
C-----------------------------------------------------------------------
        Subroutine rgplanet_V05(dlat,Rref,hgt,gz,iplanet)
C...    Computes planetary radius (Rref, km), and gravity (gz,
C       m/s**2) at latitude dlat (deg) and height (hgt, km) for
C       planet with code iplanet (see list below)
        Implicit None
        Double Precision dlat,Rref,hgt,gz,GM,requ,rpol,J2,omega,per,
     &    pi,flat,s2lat,P2
        Integer iplanet
        pi = 4.0D0*Datan(1.0D0)
C
C...    Planet number code
C       Venus = 2
C       Earth = 3
C       Mars = 4
C       Jupiter = 5
C       Saturn = 6
C       Uranus = 7
C       Neptune = 8
C       Titan = 606 (sixth satellite of Saturn)
C
C...    Planetary constants
C       GM = mass constant (m**3/s**2)
C       requ, rpol = equatorial and polar radii (km)
C       J2 = first gravity field non-spherical term (unitless)
C       per = rotation period (sec) (negative for retrograde)
        If (iplanet.eq.2)Then
C...      Venus
          GM = 3.2485863D14
          requ = 6051.893d0
          rpol = 6051.893d0
          J2 = 0.0d0
          per = -20996798.0d0
        Else If (iplanet.eq.3)Then
C...      Earth
          GM = 3.986005D14
          requ = 6378.1370D0
          rpol = 6356.7523D0
          J2 = 0.00108263D0
          per = 86164.100D0
        Else If (iplanet.eq.4)Then
C...      Mars
          GM = 4.2828314258D13
          requ = 3396.0D0
          rpol = 3378.32D0
          J2 = 0.001958616128D0
          per = 88642.66D0
        Else If (iplanet.eq.5)Then
C...      Jupiter
          GM = 1.26686537D17
          requ = 71492.0d0
          rpol = 66854.0d0
          J2 = 0.014687D0
          per = 35729.685d0
        Else If (iplanet.eq.6)Then
C...      Saturn
          GM = 3.79312845D16
          requ = 60268.0d0
          rpol = 54364.0d0
          J2 = 0.016358D0
          per = 38362.4d0
        Else If (iplanet.eq.7)Then
C...      Uranus
          GM = 5.793947D15
          requ = 25559.0d0
          rpol = 24973.0d0
          J2 = 0.003515D0
          per = -62064.0d0
        Else If (iplanet.eq.8)Then
C...      Neptune
          GM = 6.835107D15
          requ = 24766.0D0
          rpol = 24342.0D0
          J2 = 0.003540D0
          per = 57996.0D0
        Else If (iplanet.eq.606)Then
C...      Titan
          GM = 8.97803D12
          requ = 2575.5D0
          rpol = 2575.5D0
          J2 = 0.0D0
          per = 1377684.0d0
        Else
          Stop ' Invalid planet code in call to rgplanet_V05'
        Endif
C       Flattening term
        flat = (requ - rpol)/rpol
C       Rotation rate (rad/sec)
        omega = 2.0D0*pi/DAbs(per)
C...    Planetary radius at latitude
        s2lat = (Dsin(pi*dlat/1.8D2))**2
        Rref = requ/Sqrt(1.0D0 + flat*(flat+2.0D0)*s2lat)
C...    Gravity at height hgt
        P2 = 1.5D0*s2lat - 0.5D0
        gz = 1.0D-6*GM/(Rref+hgt)**2
C...    J2 correction term
        gz = gz - gz*3.0D0*J2*P2*((requ/(Rref+hgt))**4)
C...    Rotation correction term
        gz = gz - 1.0D3*(Rref+hgt)*(omega*Dcos(pi*dlat/1.8D2))**2
        Return
        End
C-----------------------------------------------------------------------
      Subroutine Winds_V05(z,clat,time,u,v)
C...  Compute zonal and meridional winds (u and v) at height z,
C     latitude clat, and local solar time.  Approximations to data in
C     Fig. 3, page 466 and Fig. 5, page 469 of "Venus II", and Fig. 8,
C     page 696 of "Venus".
C
C     Note: Many references adopt the convention for Venus that
C     super-rotating (westward, or retrograde) zonal winds are
C     positive.  We retain the traditional right-handed coordinate
C     convention, whereby zonal winds are positive eastward and
C     meridional winds are positive northward.
C
      Implicit None
      Double Precision xlat,z,u,v,clat,ue,pi180,ve,ut,time,usas,vsas
      pi180 = dAtan(1.0d0)/45.0d0
      xlat = clat
C...  Retrograde zonal wind magnitude versus height (- for westward
C     superrotation)
      ue = -1.0d0 - 99.0d0*z/80.0d0
C...  Decrease in ue above 80 km parameterized from page 333 of
C     Lellouch et al., Icarus 110, 315-319 (1994) and Fig 2. of Hou and
C     Farrell, J. Atmos. Sci. 44, 1049-1061 (1987).
      If (z.ge.80.0d0) ue = -100.0d0*(1.0d0 - (z-80.0d0)/30.0d0)
      If (ue.gt.0.0)ue = 0.0
C...  Meridional wind magnitude versus height for retrograde wind
      ve = 0.1d0*ue*dCos(2.25d0*pi180*z)
C...  Sub-solar to anti-solar diurnal wind (u peaks at terminator),
C     parameterized from Fig. 2 of Zhang et al. J. Geophys. Res.
C     101(E10), 23,195-205 (1996) and Fig 4. of Bougher et al. Icarus,
C     73, 545-573 (1988)
      ut = -3.0d0*(z - 81.0d0)
      If (ut.gt.0.0)ut = 0.0
      If (ut.lt.-237.0d0)ut = -237.0d0
C...  Time variation of sub-solar to anti-solar wind components
      usas = ut*dSin(pi180*15.0d0*(time - 12.0d0))
      vsas = -ut*dCos(pi180*15.0d0*(time - 12.0d0))
C...  Latitude variation of zonal and meridional wind
      u = (ue + usas)*(1.0d0 - (xlat/90.0)**4)
      v = ve*dSin(2.0d0*pi180*xlat) + vsas*dSin(pi180*xlat)
      Return
      End
C-----------------------------------------------------------------------
      Subroutine Venusatm_V05(CHGT,CLAT,CLON,TLOCAL,sza,HSCALE,TEMP,
     &  DENS,FACTHI,FACTLO,PRES,Hrho,AMz,EWWIND,NSWIND,ytot,yco2,yn2,
     &  yo,yco,yhe,yn,yh,profnear,proffar,nprof,profwgt)
C...  Evaluates Venus atmospheric parameters by interpolating VIRA data
C     to given height, latitude, and local solar time or solar zenith
C     angle, as necessary
      Implicit None
      Double Precision CHGT,CLAT,NSWIND,HSCALE,TEMP,DENS,FACTHI,sza,
     &  FACTLO,PRES,Hrho,AMz,EWWIND,ytot,yh,yco2,yhe,yn2,yo,yco,yn,
     &  TLOCAL,yco2md,yn2md,yomd,R0,zfact,ynb,yn1,z,T1,T2,p1,p2,R1,R2,
     &  d1,d2,yh1,yh2,yhe1,yhe2,Rgas,zlo,plo,dlo,tlo,zetalo,Rlo,ycoa,
     &  zmd,pmd,dmd,tmd,AVn,ycomd,yhemd,ynmd,zhi,phi,dhi,thi,yco2hi,
     &  yn2hi,yohi,yco21,yco22,ycohi,yhehi,ynhi,yhhi,yn21,yn22,yo1,yo2,
     &  yco1,ycob,zeta,zeta1,zeta2,pa,da,ta,zetaa,Ra,yco2a,yn2a,yoa,
     &  yna,yhea,tf,yco2f,yn2f,yof,ycof,yhef,ynf,yhf,df,pf,AMzp,
     &  profnear,proffar,profwgt,tout,pout,dout,uout,vout,yratio,CLON
      Integer j,jm,nprof
C...  Venus model VIRA data at low ,middle, and high altitudes
      Common /VIRAdata_V05/zlo(81),plo(81,5),dlo(81,5),tlo(81,5),
     &  zetalo(81,5),Rlo(81,5),zmd(11),pmd(11,2),dmd(11,2),tmd(11,2),
     &  yco2md(11,2),yn2md(11,2),yomd(11,2),ycomd(11,2),yhemd(11,2),
     &  ynmd(11,2),zhi(21),phi(21,7),dhi(21,7),
     &  thi(21,7),yco2hi(21,7),yn2hi(21,7),yohi(21,7),ycohi(21,7),
     &  yhehi(21,7),ynhi(21,7),yhhi(21,7)
C.... Universal gas constant
      R0 = 8314.472d0
C...  Avogadro's number
      AVn =  6.02214199d26
C...  Convert to height z (km) and insure 0 < z < 250 km
      z = CHGT
C...  Limit heights to > 0 km
      If (z.lt.0.0d0)z = 0.0d0
C...  For heights 0-98 km use interpolation of low altitude VIRA data
      If (z.le.zlo(80))Then
        Do 100 j = 2,80
C...      Find height index for interpolation
          If (z.le.zlo(j))Then
            jm = j - 1
C...        Get low altitude VIRA data at interpolation heights
            Call Lowterp_V05(jm,clat,p1,d1,t1,zeta1,R1,yco21,yn21,yo1,
     &        yco1)
            Call Lowterp_V05(j,clat,p2,d2,t2,zeta2,R2,yco22,yn22,yo2,
     &        ycob)
C...        Height interpolation factor
            zfact = (z-zlo(jm))/(zlo(j)-zlo(jm))
C...        Linear interpolation of temperature
            TEMP = t1 + (t2 - t1)*zfact
C...        Linear interpolation of compressibility and gas constant
            zeta = zeta1 + (zeta2 - zeta1)*zfact
            Rgas = R1 + (R2 - R1)*zfact
C...        Logarithmic interpolation of pressure
            PRES = dExp(dLog(p1) + dLog(p2/p1)*zfact)
C...        Density from compressible gas law
            DENS = PRES/(zeta*Rgas*TEMP)
C...        Mean molecular weight
            AMz = R0/Rgas
C...        Logarithmic interpolation of CO2 and N2 number densities
            yco2 = dExp(dLog(yco21) + dLog(yco22/yco21)*zfact)
            yn2 = dExp(dLog(yn21) + dLog(yn22/yn21)*zfact)
C...        Logarithmic interpolation of O and CO, avoiding problems
C           with transition from 0 values below 86 km
            If (yo2.le.0.0)Then
              yo = 0.0
            Else
              If (yo1.le.0.0)yo1 = 0.1d0*yo2
              yo = dExp(dLog(yo1) + dLog(yo2/yo1)*zfact)
            Endif
            If (ycob.le.0.0)Then
              yco = 0.0
            Else
              If (yco1.le.0.0)yco1 = 0.1d0*ycob
              yco = dExp(dLog(yco1) + dLog(ycob/yco1)*zfact)
            Endif
C...        Helium, and atomic H and N not present in low atmosphere
            yhe = 0.0
            yn = 0.0
            yh = 0.0
C...        Total number density
            ytot = yco2 + yn2 + yo + yco
C...        Scale heights for pressure and density
            HSCALE = (zlo(j)-zlo(jm))/dLog(p1/p2)
            Hrho = (zlo(j)-zlo(jm))/dLog(d1/d2)
            Goto 900
          Endif
 100     Enddo
      Else If (z.le.zlo(81))Then
C...    For height 98-100 km, interpolate using combination of low and
C       middle altitude VIRA data
C...    Get low altitude VIRA data at 98 km
        Call Lowterp_V05(80,clat,p1,d1,t1,zeta1,R1,yco21,yn21,yo1,yco1)
C...    Get low altitude VIRA data at 100 km
        Call Lowterp_V05(81,clat,pa,da,ta,zetaa,Ra,yco2a,yn2a,yoa,ycoa)
C...    Get middle altitude VIRA data at 100 km
        Call Midterp_V05(1,TLOCAL,clat,p2,d2,t2,yco22,yn22,yo2,ycob,
     &    yhe2,ynb)
C...    Assume He and H values for logarithmic interpolation
        yhe1 = 0.1d0*yhe2
        yn1 = 0.1d0*ynb
C...    Logarithmic average pressure and density at 100 km from low
C       and middle altitude VIRA data
        p2 = dSqrt(pa*p2)
        d2 = dSqrt(da*d2)
C...    Linear average temperature at 100 km from low and middle
C       altitude VIRA data
        t2 = (ta + t2)/2.0d0
C...    Logarithmic average number densities at 100 km from low
C       and middle altitude VIRA data
        yco22 = dSqrt(yco2a*yco22)
        yn22 = dSqrt(yn2a*yn22)
        yo2 = dSqrt(yoa*yo2)
        ycob = dSqrt(ycoa*ycob)
C...    Height interpolation factor
        zfact = (z - zlo(80))/(zlo(81) - zlo(80))
C...    Linear interpolation of temperature
        TEMP = t1 + (t2 - t1)*zfact
C...    Gas constant at 98 and 100 km
        R1 = p1/(d1*t1)
        R2 = p2/(d2*t2)
C...    Linear interpolation of gas constant
        Rgas = R1 + (R2 - R1)*zfact
C...    Logarithmic interpolation of pressure
        PRES = dExp(dLog(p1) + dLog(p2/p1)*zfact)
C...    Density from perfect gas law (zeta = 1)
        DENS = PRES/(Rgas*TEMP)
C...    Mean molecular weight
        AMz = R0/Rgas
C...    Logarithmic interpolation of number densities
        yco2 = dExp(dLog(yco21) + dLog(yco22/yco21)*zfact)
        yn2 = dExp(dLog(yn21) + dLog(yn22/yn21)*zfact)
        yo = dExp(dLog(yo1) + dLog(yo2/yo1)*zfact)
        yco = dExp(dLog(yco1) + dLog(ycob/yco1)*zfact)
        yhe = dExp(dLog(yhe1) + dLog(yhe2/yhe1)*zfact)
        yn = dExp(dLog(yn1) + dLog(ynb/yn1)*zfact)
C.......Atomic hydrogen not present at this height
        yh = 0.0
C...    Total number density
        ytot = yco2 + yn2 + yo + yco + yhe + yn
C...    Scale heights for pressure and density
        HSCALE = (zlo(81)-zlo(80))/dLog(p1/p2)
        Hrho = (zlo(81)-zlo(80))/dLog(d1/d2)
        Goto 900
      Else If (z.le.zmd(2))Then
C...    For height 100-105 km, interpolate using combination of low and
C       middle altitude VIRA data
C...    Get low altitude VIRA data at 100 km
        Call Lowterp_V05(81,clat,p1,d1,t1,zeta1,R1,yco21,yn21,yo1,yco1)
C...    Get middle altitude VIRA data at 100 km
        Call Midterp_V05(1,TLOCAL,clat,pa,da,ta,yco2a,yn2a,yoa,ycoa,
     &    yhea,yna)
C...    Get middle altitude VIRA data at 105 km
        Call Midterp_V05(2,TLOCAL,clat,p2,d2,t2,yco22,yn22,yo2,ycob,
     &    yhe2,ynb)
C...    Logarithmic average pressure and density at 100 km from low
C       and middle altitude VIRA data
        p1 = dSqrt(pa*p1)
        d1 = dSqrt(da*d1)
C...    Linear average temperature at 100 km from low and middle
C       altitude VIRA data
        t1 = (ta + t1)/2.0d0
C...    Logarithmic average number densities at 100 km from low
C       and middle altitude VIRA data
        yco21 = dSqrt(yco2a*yco21)
        yn21 = dSqrt(yn2a*yn21)
        yo1 = dSqrt(yoa*yo1)
        yco1 = dSqrt(ycoa*yco1)
C...    Take 100 km He and N number densities from middle atmosphere
C       VIRA data
        yhe1 = yhea
        yn1 = yna
C...    Height interpolation factor
        zfact = (z - zmd(1))/(zmd(2) - zmd(1))
C...    Linear interpolation of temperature
        TEMP = t1 + (t2 - t1)*zfact
C...    Gas constant at 100 and 105 km
        R1 = p1/(d1*t1)
        R2 = p2/(d2*t2)
C...    Linear interpolation of gas constant
        Rgas = R1 + (R2 - R1)*zfact
C...    Logarithmic interpolation of pressure
        PRES = dExp(dLog(p1) + dLog(p2/p1)*zfact)
C...    Density from perfect gas law (zeta = 1)
        DENS = PRES/(Rgas*TEMP)
C...    Mean molecular weight
        AMz = R0/Rgas
C...    Logarithmic interpolation of number densities
        yco2 = dExp(dLog(yco21) + dLog(yco22/yco21)*zfact)
        yn2 = dExp(dLog(yn21) + dLog(yn22/yn21)*zfact)
        yo = dExp(dLog(yo1) + dLog(yo2/yo1)*zfact)
        yco = dExp(dLog(yco1) + dLog(ycob/yco1)*zfact)
        yhe = dExp(dLog(yhe1) + dLog(yhe2/yhe1)*zfact)
        yn = dExp(dLog(yn1) + dLog(ynb/yn1)*zfact)
C...    Atomic hydrogen not present at this height
        yh = 0.0
C...    Total number density
        ytot = yco2 + yn2 + yo + yco + yhe + yn
C...    Scale heights for pressure and density
        HSCALE = (zmd(2)-zmd(1))/dLog(p1/p2)
        Hrho = (zmd(2)-zmd(1))/dLog(d1/d2)
        Goto 900
      Else If (z.le.zmd(10))Then
C...    For height 105-145 km, interpolate using middle altitude
C       VIRA data
        Do 200 j = 2,10
C...      Find height index for interpolation
          If (z.le.zmd(j))Then
            jm = j - 1
C...        Get middle altitude VIRA data at interpolation heights
            Call Midterp_V05(jm,TLOCAL,clat,p1,d1,t1,yco21,yn21,yo1,
     &        yco1,yhe1,yn1)
            Call Midterp_V05(j,TLOCAL,clat,p2,d2,t2,yco22,yn22,yo2,
     &        ycob,yhe2,ynb)
C...        Height interpolation factor
            zfact = (z-zmd(jm))/(zmd(j)-zmd(jm))
C...        Linear interpolation of temperature
            TEMP = t1 + (t2 - t1)*zfact
C...        Gas constant at two interpolation heights
            R1 = p1/(d1*t1)
            R2 = p2/(d2*t2)
C...        Linear interpolation of gas constant
            Rgas = R1 + (R2 - R1)*zfact
C...        Logarithmic interpolation of pressure
            PRES = dExp(dLog(p1) + dLog(p2/p1)*zfact)
C...        Density from perfect gas law
            DENS = PRES/(Rgas*TEMP)
C...        Mean molecular weight
            AMz = R0/Rgas
C...        Logarithmic interpolation of number densities
            yco2 = dExp(dLog(yco21) + dLog(yco22/yco21)*zfact)
            yn2 = dExp(dLog(yn21) + dLog(yn22/yn21)*zfact)
            yo = dExp(dLog(yo1) + dLog(yo2/yo1)*zfact)
            yco = dExp(dLog(yco1) + dLog(ycob/yco1)*zfact)
            yhe = dExp(dLog(yhe1) + dLog(yhe2/yhe1)*zfact)
            yn = dExp(dLog(yn1) + dLog(ynb/yn1)*zfact)
C...        Atomic hydrogen not present in this height range
            yh = 0.0
C...        Total number density
            ytot = yco2 + yn2 + yo + yco + yhe + yn
C...        Scale heights for pressure and density
            HSCALE = (zmd(j)-zmd(jm))/dLog(p1/p2)
            Hrho = (zmd(j)-zmd(jm))/dLog(d1/d2)
            Goto 900
           Endif
 200    Enddo
      Else If (z.le.zmd(11))Then
C...    For height 145-150 km, interpolate using combination of middle
C       and high altitude VIRA data
C...    Get middle altitude VIRA data at 145 km
        Call Midterp_V05(10,TLOCAL,clat,p1,d1,t1,yco21,yn21,yo1,yco1,
     &    yhe1,yn1)
C...    Get middle altitude VIRA data at 150 km
        Call Midterp_V05(11,TLOCAL,clat,pa,da,ta,yco2a,yn2a,yoa,ycoa,
     &    yhea,yna)
C...    Get high altitude VIRA data at 150 km
        Call Highterp_V05(1,sza,p2,d2,t2,yco22,yn22,yo2,ycob,yhe2,ynb,
     &    yh2)
C...    Logarithmic average pressure and density at 150 km from middle
C       and high altitude VIRA data
        p2 = dSqrt(pa*p2)
        d2 = dSqrt(da*d2)
C...    Linear average temperature at 150 km from middle and high
C       altitude VIRA data
        t2 = (ta + t2)/2.0d0
C...    Logarithmic average number densities at 150 km from middle
C       and high altitude VIRA data
        yco22 = dSqrt(yco2a*yco22)
        yn22 = dSqrt(yn2a*yn22)
        yo2 = dSqrt(yoa*yo2)
        ycob = dSqrt(ycoa*ycob)
        yhe2 = dSqrt(yhea*yhe2)
        ynb = dSqrt(yna*ynb)
C...    Assume H number density at 145 km for logarithmic interpolation
        yh1 = 0.1d0*yh2
C...    Height interpolation factor
        zfact = (z - zmd(10))/(zmd(11) - zmd(10))
C...    Linear interpolation of temperature
        TEMP = t1 + (t2 - t1)*zfact
C...    Gas constant at 145 and 150 km
        R1 = p1/(d1*t1)
        R2 = p2/(d2*t2)
C...    Linear interpolation of gas constant
        Rgas = R1 + (R2 - R1)*zfact
C...    Logarithmic interpolation of pressure
        PRES = dExp(dLog(p1) + dLog(p2/p1)*zfact)
C...    Density from perfect gas law
        DENS = PRES/(Rgas*TEMP)
C...    Mean molecular weight
        AMz = R0/Rgas
C...    Logarithmic interpolation of number densities
        yco2 = dExp(dLog(yco21) + dLog(yco22/yco21)*zfact)
        yn2 = dExp(dLog(yn21) + dLog(yn22/yn21)*zfact)
        yo = dExp(dLog(yo1) + dLog(yo2/yo1)*zfact)
        yco = dExp(dLog(yco1) + dLog(ycob/yco1)*zfact)
        yhe = dExp(dLog(yhe1) + dLog(yhe2/yhe1)*zfact)
        yn = dExp(dLog(yn1) + dLog(ynb/yn1)*zfact)
        yh = dExp(dLog(yh1) + dLog(yh2/yh1)*zfact)
C...    Total number density
        ytot = yco2 + yn2 + yo + yco + yhe + yn + yh
C...    Scale heights for pressure and density
        HSCALE = (zmd(11)-zmd(10))/dLog(p1/p2)
        Hrho = (zmd(11)-zmd(10))/dLog(d1/d2)
        Goto 900
      Else If (z.le.zhi(2))Then
C...    For height 150-155 km, interpolate using combination of middle
C       and high altitude VIRA data
C...    Get middle altitude VIRA data at 150 km
        Call Midterp_V05(11,TLOCAL,clat,pa,da,ta,yco2a,yn2a,yoa,ycoa,
     &    yhea,yna)
C...    Get high altitude VIRA data at 150 km
        Call Highterp_V05(1,sza,p1,d1,t1,yco21,yn21,yo1,yco1,yhe1,yn1,
     &    yh1)
C...    Get high altitude VIRA data at 155 km
        Call Highterp_V05(2,sza,p2,d2,t2,yco22,yn22,yo2,ycob,yhe2,ynb,
     &    yh2)
C...    Logarithmic average pressure and density at 150 km from middle
C       and high altitude VIRA data
        p1 = dSqrt(pa*p1)
        d1 = dSqrt(da*d1)
C...    Linear average temperature at 150 km from middle and high
C       altitude VIRA data
        t1 = (ta + t1)/2.0d0
C...    Logarithmic average number densities at 150 km from middle
C       and high altitude VIRA data
        yco21 = dSqrt(yco2a*yco21)
        yn21 = dSqrt(yn2a*yn21)
        yo1 = dSqrt(yoa*yo1)
        yco1 = dSqrt(ycoa*yco1)
        yhe1 = dSqrt(yhea*yhe1)
        yn1 = dSqrt(yna*yn1)
C...    Height interpolation factor
        zfact = (z - zhi(1))/(zhi(2) - zhi(1))
C...    Linear interpolation of temperature
        TEMP = t1 + (t2 - t1)*zfact
C...    Gas constant at 150 and 155 km
        R1 = p1/(d1*t1)
        R2 = p2/(d2*t2)
C...    Linear interpolation of gas constant
        Rgas = R1 + (R2 - R1)*zfact
C...    Logarithmic interpolation of pressure
        PRES = dExp(dLog(p1) + dLog(p2/p1)*zfact)
C...    Density from perfect gas law
        DENS = PRES/(Rgas*TEMP)
C...    Mean molecular weight
        AMz = R0/Rgas
C...    Logarithmic interpolation of number densities
        yco2 = dExp(dLog(yco21) + dLog(yco22/yco21)*zfact)
        yn2 = dExp(dLog(yn21) + dLog(yn22/yn21)*zfact)
        yo = dExp(dLog(yo1) + dLog(yo2/yo1)*zfact)
        yco = dExp(dLog(yco1) + dLog(ycob/yco1)*zfact)
        yhe = dExp(dLog(yhe1) + dLog(yhe2/yhe1)*zfact)
        yn = dExp(dLog(yn1) + dLog(ynb/yn1)*zfact)
        yh = dExp(dLog(yh1) + dLog(yh2/yh1)*zfact)
C...    Total number density
        ytot = yco2 + yn2 + yo + yco + yhe + yn + yh
C...    Scale heights for pressure and density
        HSCALE = (zhi(2)-zhi(1))/dLog(p1/p2)
        Hrho = (zhi(2)-zhi(1))/dLog(d1/d2)
        Goto 900
      Else If (z.le.zhi(21))Then
C...    For height 155-250 km, interpolate using high altitude
C       VIRA data
        Do 300 j = 2,21
C...      Find height index for interpolation
          If (z.le.zhi(j))Then
            jm = j - 1
C...        Get high altitude VIRA data at interpolation heights
            Call Highterp_V05(jm,sza,p1,d1,t1,yco21,yn21,yo1,yco1,yhe1,
     &        yn1,yh1)
            Call Highterp_V05(j,sza,p2,d2,t2,yco22,yn22,yo2,ycob,yhe2,
     &        ynb,yh2)
C...        Height interpolation factor
            zfact = (z-zhi(jm))/(zhi(j)-zhi(jm))
C...        Linear interpolation of temperature
            TEMP = t1 + (t2 - t1)*zfact
C...        Gas constant at two interpolation heights
            R1 = p1/(d1*t1)
            R2 = p2/(d2*t2)
C...        Linear interpolation of gas constant
            Rgas = R1 + (R2 - R1)*zfact
            PRES = dExp(dLog(p1) + dLog(p2/p1)*zfact)
C...        Density from perfect gas law
            DENS = PRES/(Rgas*TEMP)
C...        Mean molecular weight
            AMz = R0/Rgas
C...        Logarithmic interpolation of number densities
            yco2 = dExp(dLog(yco21) + dLog(yco22/yco21)*zfact)
            yn2 = dExp(dLog(yn21) + dLog(yn22/yn21)*zfact)
            yo = dExp(dLog(yo1) + dLog(yo2/yo1)*zfact)
            yco = dExp(dLog(yco1) + dLog(ycob/yco1)*zfact)
            yhe = dExp(dLog(yhe1) + dLog(yhe2/yhe1)*zfact)
            yn = dExp(dLog(yn1) + dLog(ynb/yn1)*zfact)
            yh = dExp(dLog(yh1) + dLog(yh2/yh1)*zfact)
C...        Total number density
            ytot = yco2 + yn2 + yo + yco + yhe + yn + yh
C...        Scale heights for pressure and density
            HSCALE = (zhi(j)-zhi(jm))/dLog(p1/p2)
            Hrho = (zhi(j)-zhi(jm))/dLog(d1/d2)
            Goto 900
          Endif
 300    Enddo
      Else
C...  Evaluate atmosphere with constant temperature thermosphere
C...    Get high altitude VIRA data at 250 km
        Call Highterp_V05(21,sza,pf,df,tf,yco2f,yn2f,yof,ycof,yhef,ynf,
     &    yhf)
C...    Get thermospheric model values at heights z+1 km and z
        TEMP = tf
        Call Thermos_V05(z+1.0d0,clat,250.0d0,yco2f,yn2f,yof,ycof,yhef,
     &    ynf,yhf,tf,PRES,DENS,AMzp,yco2,yn2,yo,yco,yhe,yn,yh,ytot,
     &    HSCALE)
        Call Thermos_V05(z,clat,250.0d0,yco2f,yn2f,yof,ycof,yhef,ynf,
     &    yhf,tf,PRES,DENS,AMz,yco2,yn2,yo,yco,yhe,yn,yh,ytot,HSCALE)
C...    Get density scale height from pressure scale height and dM/dz
        Hrho = HSCALE/(1.0d0 - HSCALE*(AMzp-AMz)/AMz)
      Endif
C...  Get wind components from approximation to VIRA data
 900  Call Winds_V05(z,CLAT,TLOCAL,EWWIND,NSWIND)
C...  Compute high and low factors for perturbations: parameterized
C     from data in Figs. 1-8(d) and 1-12(a) of Seiff et al. VIRA data,
C     pp. 247, 259 & 278 of "Venus", pp. 200 & 201 of "The Planet
C     Venus", p. 283. of "Venus II", and Fig. 6 of Hinson and Jenkins
C     Icarus 114, 310-327 (1995).
      If (z.lt.50.0d0)Then
        FACTHI = 1.02d0
      Else
        FACTHI = 1.02d0 + 0.0013d0*(z-50.0d0)
      Endif
      If (FACTHI.gt.1.15d0)FACTHI = 1.15d0
      FACTLO = 1.0d0/FACTHI
C---  Use weighted average profile data if profnear > 0. Weight=1 if
C     lat-lon radius < profnear. Weight=0 if lat-lon radius > proffar.
      If (profnear.gt.0.0d0)Then
        Call ProfTerp_V05(z,CLAT,CLON,TEMP,PRES,DENS,EWWIND,NSWIND,
     &    tout,pout,dout,uout,vout,nprof,profnear,proffar,profwgt)
        yratio = dout/DENS
        yco2 = yco2*yratio
        yn2 = yn2*yratio
        yo = yo*yratio
        yco = yco*yratio
        yhe = yhe*yratio
        yn = yn*yratio
        yh = yh*yratio
        ytot = ytot*yratio
        TEMP = tout
        PRES = pout
        DENS = dout
        EWWIND = uout
        NSWIND = vout
      Endif
      Return
      End
C-----------------------------------------------------------------------
      Subroutine Venusref_V05(CHGT,tvra,pvra,dvra)
C...  Evaluates Venus reference atmospheric parameters from approximate
C     average VIRA data
      Implicit None
      Double Precision CHGT,tvra,pvra,dvra,z,T1,T2,p1,p2,R1,R2,delz,
     &  dztot,R,HSCALE,zref,pref,dref,tref,AMz,yco2,yn2,yo,yco,yhe,yn,
     &  yh,totnd
      Integer iz,i
C...  Common block for VIRA reference data arrays
      Common /VIRAref_V05/zref(111),pref(111),dref(111),tref(111)
      i = 0 ! gg avoids warning about uninitialization
C...  Convert to height z (km)
      z = CHGT
C...  Check that height is in bounds
      If (z.lt.zref(1))z = zref(1)
      If (z.gt.zref(111))Then
C...    Evaluate thermosphere at height z for SZA = 90 deg and 250 km
C       boundary values
        Call Thermos_V05(z,0.0d0,250.0d0,3.07D6,8.61D8,1.05D12,1.13D9,
     &    3.07D12,2.18D10,1.05D12,230.0d0,pvra,dvra,AMz,yco2,yn2,yo,
     &    yco,yhe,yn,yh,totnd,hscale)
        tvra = 230.0d0
        Return
      Endif
C...  Find height index for vertical interpolation
      Do 10 iz = 1,110
        If (zref(iz).le.z.and.zref(iz+1).ge.z)Then
          i = iz
          Goto 20
        Endif
  10  Enddo
C...  Get Venus average temperature, pressure, and gas law constant
C     at upper and lower height indexes
  20  T1 = tref(i)
      T2 = tref(i+1)
      p1 = pref(i)
      p2 = pref(i+1)
      R1 = pref(i)/(dref(i)*tref(i))
      R2 = pref(i+1)/(dref(i+1)*tref(i+1))
C...  Linear height interpolation on temperature
      delz = z - zref(i)
      dztot = zref(i+1) - zref(i)
      tvra = T1 + (T2 - T1)*delz/dztot
C...  Pressure scale height and vertical pressure interpolation
      HSCALE = dztot/dLog(p1/p2)
      pvra = P1*dExp(-delz/HSCALE)
C...  Linear height interpolation for gas constant
      R = R1 + (R2 - R1)*delz/dztot
C...  Density from perfect gas law
      dvra = pvra/(R*tvra)
      Return
      End
C-----------------------------------------------------------------------
      Subroutine Rescale_V05(x)
C...  Puts x into range 0 - 360
      Double precision x
      x = x/360.0d0 - Dint(x/360.0d0) + 1.0d0
      x = (x - Dint(x))*360.0d0
      Return
      End
C----------------------------------------------------------------------
      Subroutine Shiftdif_V05(x)
C...  Shifts difference x to be +/- and close to 0.0
      Double Precision x
      If (x.gt.180.0d0)Then
        x = x - 360.0d0
      Else If (x.lt.-180.0d0)Then
        x = x + 360.0d0
      Endif
      Return
      End
C----------------------------------------------------------------------
      Subroutine venephem_V05(xday,sunlat,sunlon,sunLsubs,radius,owlt,
     &  EOT)
C...  Computes sunlat, sunlon= latitude and longitude of sub-solar
C     point on the surface, sunLsubs= planetocentric longitude of Sun
C     (Ls), radius= current orbital radius from Sun to Venus, heliolon=
C     Venus heliocentric longitude, owlt= Venus-Earth one-way light
C     time (minutes), and EOT= equation of time (deg), calculated from
C     Julian day and time, xday.  Notes: input xday is NOT UTC, but
C     Terrestrial (Dynamical) Venus-Event Time (NOT Earth-Receive Time).
C     Venus Local Mean Solar Time (hrs ) = Local True Solar Time (hrs)
C     minus EOT (in hrs). Output is for Terrestrial (Dynamical)
C     Venus Event Time (corresponding to input xday).
C
C     Equations for "moderately accurate" Venus solar time, seasonal
C     parameters, and one-way Venus-Earth light time, adapted from
C     Allison and McEwen, Planet. Space Sci., 48, 215-235 (2000),
C     and Allison, Geophys Res. Lett., 24(16), 1967-1970 (1997).
C
      Implicit None
      Double Precision xday,sunlat,sunlon,sunLsubs,radius,pi180,dt,
     &  anomM,alphFMS,alsrad,helilon,alphs,EOT,pmr,gE,rE,dlat1,
     &  helonE,Vm,owlt,yranom,anom0,yrtrop,veqlon0,perlon0,ecc,obl,
     &  dlat2,dlat3,veqlon1,inc,anlon0,siday,rad0,ecc2,ecc3,ecc4,
     &  ecc5,ecc6,Vm0,Vmday0,xE,yE,xpl,ypl,zpl,eqcenter,argper,
     &  trueanom,coslat
      pi180 = Datan(1.0d0)/45.0d0
C...  Days since 2000 January 1.5
      dt = xday - 2451545.0d0
C
C.....................................................................
C
C...  Planetary orbit parameters
C
C     Semi-major axis (AU) = mean distance from Sun
      rad0 = 0.723330d0
C     Anomalistic year (days, perihelion-to-perihelion)
      yranom = 224.7d0
C     Tropical year (days, for rate of fictitious mean sun)
      yrtrop = 224.6994d0
C     Mean anomaly for J2000 (degrees)
      anom0 = 50.4084d0
C     Heliocentric longitude of perihelion at J2000 (deg)
      perlon0 = 131.5709d0
C     Terms for heliocentric longitude at Ls=0 (deg)
      veqlon0 = 237.841d0
      veqlon1 = 9.77d-6
C     Eccentricity and powers
      ecc = 0.006773d0 - 1.302D-9*dt
      ecc2 = ecc**2
      ecc3 = ecc2*ecc
      ecc4 = ecc3*ecc
      ecc5 = ecc4*ecc
      ecc6 = ecc5*ecc
C     Obliquity angle (radians)
      obl = (177.36d0 + 0.0d0*dt)*pi180
C     Inclination (radians)
      inc = (3.3946d0 +2.75D-8*dt)*pi180
C     Longitude of ascending node at J2000 (deg)
      anlon0 = 76.6799d0
C     Sidereal period of rotation (Earth days)
      siday = -243.02011d0
C     Heliocentric lon of prime meridian (deg) at Julian day Vmday0
      Vm0 = 75.1289d0
      Vmday0 = 2451545.0d0
C     Difference terms, planetocentric to planetographic lat (deg)
      dlat1 = -0.002d0
      dlat2 = 0.002d0
      dlat3 = 0.0d0
C
C.....................................................................
C
C...  Mean anomaly (radians)
C...  Allison & McEwen (2000) equation (16)
      anomM = (anom0 + (360.0d0/yranom)*dt)*pi180
C...  Right ascension of fictitious mean sun (deg)
C...  Allison & McEwen (2000) equation (17)
      alphFMS = perlon0 - veqlon0 + anom0 + (360.0d0/yrtrop)*dt
C...  Venus equation of center, A&M eqn. (4) (degrees)
      eqcenter = ((2.0d0*ecc - 0.25d0*ecc3 +
     &  (5.0d0/96.0d0)*ecc5)*Dsin(anomM) +
     &  (1.25d0*ecc2 - (11.0d0/24.0d0)*ecc4 +
     &  (17.0d0/192.0d0)*ecc6)*Dsin(2.0d0*anomM) +
     &  ((13.0d0/12.0d0)*ecc3 - (43.0d0/63.0d0)*ecc5)*
     &  Dsin(3.0d0*anomM) + ((103.0d0/96.0d0)*ecc4 -
     &  (451.0d0/480.0d0)*ecc6)*Dsin(4.0d0*anomM) +
     &  ((1097.0d0/960.0d0)*ecc5)*Dsin(5.0d0*anomM) +
     &  ((12323.0d0/960.0d0)*ecc6)*Dsin(6.0d0*anomM))/pi180
C...  True planetocentric solar longitude (Ls), A&M eqns. (2) and (4)
      sunLsubs = alphFMS + eqcenter
      Call Rescale_V05(sunLsubs)
C...  Ls angle in radians
      alsrad = sunLsubs*pi180
C...  Sub-solar latitude of sun (planetographic solar declination),
C     Allison (1997) eqn. (5) with empirical Ls and 3*Ls terms
      sunlat = DAsin(Sin(obl)*Dsin(alsrad))/pi180 + dlat1*Dsin(alsrad)
     &  + dlat2*Dcos(alsrad) + dlat3*Dsin(3.0d0*alsrad)
C...  Solar right ascension, un-numbered equation, A&M page 217
      alphs = Datan2(Dcos(obl)*dSin(alsrad),dCos(alsrad))/pi180
C...  Venus orbital radius, Astronomical Almanac page E4
      radius = rad0*(1.0d0 - ecc2)/(1.0d0 + ecc*DCos(anomM + alsrad -
     &  alphFMS*pi180))
C.... Approximate Venus heliocentric longitude, A&M eqn, (11)
      helilon = sunLsubs + veqlon0 - veqlon1*dt-(Dtan(0.5d0*inc)**2)*
     &  Dsin(2.0d0*(alsrad + (veqlon0 - anlon0)*pi180))/pi180
      Call Rescale_V05(helilon)
C...  Equation of time (deg)
      EOT = alphFMS - alphs
      Call Rescale_V05(EOT)
      Call Shiftdif_V05(EOT)
      Call Rescale_V05(alphs)
C...  Earth heliocentric distance and longitude, Allison eqns (20)-
C     (22)
      gE = (357.528d0 + 0.9856003d0*dt)*pi180
      rE = 1.00014d0 - 0.01671d0*Dcos(gE) - 0.00014d0*Dcos(2.0d0*gE)
      helonE = 100.472d0 + 0.9856474d0*dt + 1.915d0*Dsin(gE)
     &  + 0.020d0*Dsin(2.0d0*gE)
C...  Earth Cartesian coordinates
      xE = rE*dCos(helonE*pi180)
      yE = rE*dSin(helonE*pi180)
C...  Venus true anolmaly (radians)
      trueanom = eqcenter*pi180 + anomM
C...  Venus argument of perihelion (radians)
      argper = (54.8910d0 + 1.38374d-5*dt)*pi180
C...  Venus Cartesian coordinates
      zpl = radius*dSin(trueanom + argper)*dSin(inc)
      coslat = dSqrt(1.0d0 - (zpl/radius)**2)
      xpl = radius*dCos((helilon+3.82394d-5*dt)*pi180)*coslat
      ypl = radius*dSin((helilon+3.82394d-5*dt)*pi180)*coslat
C...  One-way light time (minutes), Allison eqn.(19)
      owlt = dSqrt((xpl-xE)**2+(ypl-yE)**2+zpl**2)*499.005d0/60.0d0
C...  Venus (Heliocentric) prime meridian, Allison eqn (11)
      Vm = Vm0 + (360.0d0/dAbs(siday))*(xday - Vmday0)
C...  Sub-solar longitude from true solar time at prime meridian,
C     A&M page 217
      pmr = (Vm - alphs)/360.0d0
      sunlon = (pmr - Dint(pmr))*360.0d0 + 180.0d0
      Call Rescale_V05(sunlon)
      Return
      End
C----------------------------------------------------------------------
      Subroutine SolZenAng_V05(sunlat,sunlon,sitelat,sitelon,sza)
C...  Solar zenith angle (sza, degrees) from latitude and longitude
C     of sun and site (degrees)
      Implicit None
      Double Precision sunlat,sunlon,sitelat,sitelon,sza,pi180,csza
      pi180 = dAtan(1.0d0)/45.0d0
C...  Cosine of solar zenith angle
      csza = dSin(sunlat*pi180)*dSin(sitelat*pi180) +
     &  dCos(sunlat*pi180)*dCos(sitelat*pi180)*
     &  dCos(pi180*(sunlon-sitelon))
C...  Solar zenith angle
      sza = dAcos(csza)/pi180
      Return
      End
C----------------------------------------------------------------------
      Subroutine Lowterp_V05(i,clat,p,d,t,zeta,Rgas,yco2,yn2,yo,yco)
C...  Interpolation routine for low altitude VIRA data (0-100 km)
C     Computes pressure, density, temperature (p,d,t), compressibility
C     (zeta), gas constant (Rgas), and number densities (yco2,yn2,yo,
C     and yco) from VIRA low-altitude index (i) and current latitude
C     (clat)
      Implicit None
      Double Precision clat,p,d,t,zeta,Rgas,zlo,plo,dlo,tlo,zetalo,
     &  Rlo,zmd,pmd,dmd,tmd,yco2md,yn2md,yomd,ycomd,yhemd,ynmd,R0,
     &  zhi,phi,dhi,thi,yco2hi,yn2hi,yohi,ycohi,yhehi,ynhi,yhhi,AM,
     &  vlat(5),xlat,flat,yco2,yn2,yo,yco,fco2,fn2,fo,fco,AVn,AMx
      Integer i,j,jm
C...  Venus model VIRA data at low ,middle, and high altitudes
      Common /VIRAdata_V05/zlo(81),plo(81,5),dlo(81,5),tlo(81,5),
     &  zetalo(81,5),Rlo(81,5),zmd(11),pmd(11,2),dmd(11,2),tmd(11,2),
     &  yco2md(11,2),yn2md(11,2),yomd(11,2),ycomd(11,2),yhemd(11,2),
     &  ynmd(11,2),zhi(21),phi(21,7),dhi(21,7),
     &  thi(21,7),yco2hi(21,7),yn2hi(21,7),yohi(21,7),ycohi(21,7),
     &  yhehi(21,7),ynhi(21,7),yhhi(21,7)
C...  Latitudes for low-altitude VIRA data
      Data vlat/30.0d0,45.0d0,60.0d0,75.0d0,85.0d0/
C...  Universal gas constant
      R0 = 8314.472d0
C...  Avogadro's number
      AVn = 6.02214199d26
      xlat = dAbs(clat)
C...  Use VIRA data for lat=0-30 if abs(lat)<30
      If (xlat.le.vlat(1))Then
        p = plo(i,1)
        d = dlo(i,1)
        t = tlo(i,1)
        zeta = zetalo(i,1)
        Rgas = Rlo(i,1)
C...  Use VIRA data for lat=85-90 if abs(lat)>85
      Else If (xlat.ge.vlat(5))Then
        p = plo(i,5)
        d = dlo(i,5)
        t = tlo(i,5)
        zeta = zetalo(i,5)
        Rgas = Rlo(i,5)
      Else
C...    Find latitude interpolation index j
        Do 10 j = 2,5
          If (xlat.le.vlat(j))Then
            jm = j - 1
C...        Get latitude interpolation factor
            flat = (xlat - vlat(jm))/(vlat(j) - vlat(jm))
C...        Logarithmic interpolation on pressure
            p = dExp(dLog(plo(i,jm)) + dLog(plo(i,j)/plo(i,jm))*flat)
C...        Linear interpolation on temperature, zeta, and gas constant
            t = tlo(i,jm) + (tlo(i,j)-tlo(i,jm))*flat
            zeta = zetalo(i,jm) + (zetalo(i,j)-zetalo(i,jm))*flat
            Rgas = Rlo(i,jm) + (Rlo(i,j)-Rlo(i,jm))*flat
C,,,        Density from perfect gas law (with compressibility zeta)
            d = p/(zeta*Rgas*t)
            Goto 20
          Endif
  10    Enddo
      Endif
  20  AM = R0/Rgas
C...  Use uniform mixing ratios of CO2&N2 for heights up to 82 km
      If (i.lt.73)Then
        fco2 = 0.965d0
        fn2 = 0.035d0
        fco = 0.0
        fo = 0.0
      Else
C...    Use height-variable CO2/N2/O/CO mixing ratios from 82 to 100 km
        AMx = AM
        If (AMx.gt.43.44d0)AMx = 43.44d0
C...    Get mixing ratios from mean molecular weight
        fco2 = -1.723936d0 + 6.19d-2*AMx
        fn2 = 2.515424d0 - 5.71d-2*AMx
        fo = 3.4752d-2 - 8.0d-4*AMx
        fco = 0.17376d0 - 4.0d-3*AMx
        If (fo.le.0.0)fo = 0.0
        If (fco.le.0.0)fco = 0.0
      Endif
C...  Get number densities from mixing ratios and mass density
      yco2 = fco2*d*AVn/AM
      yn2 = fn2*d*AVn/AM
      yco = fco*d*AVn/AM
      yo = fo*d*AVn/AM
      Return
      End
C----------------------------------------------------------------------
      Subroutine Midterp_V05(i,time,clat,p,d,t,yco2,yn2,yo,yco,yhe,yn)
C...  Interpolation routine for middle altitude VIRA data (100-150 km)
C     Computes pressure, density, temperature (p,d,t), and number
C     densities (yco2,yn2,yo, yco, yhe and yn) from VIRA midddle-
C     altitude index (i), current local solar time (time) and current
C     latitude (clat)
      Implicit None
      Double Precision time,p,d,t,zlo,plo,dlo,tlo,zetalo,
     &  Rlo,zmd,pmd,dmd,tmd,yco2md,yn2md,yomd,ycomd,yhemd,ynmd,
     &  zhi,phi,dhi,thi,yco2hi,yn2hi,yohi,ycohi,yhehi,ynhi,yhhi,
     &  yco2,yn2,yo,yco,yhe,yn,ftime,A,B,R1,R2,Rgas,clat,alat,flat
      Integer i
C...  Venus model VIRA data at low ,middle, and high altitudes
      Common /VIRAdata_V05/zlo(81),plo(81,5),dlo(81,5),tlo(81,5),
     &  zetalo(81,5),Rlo(81,5),zmd(11),pmd(11,2),dmd(11,2),tmd(11,2),
     &  yco2md(11,2),yn2md(11,2),yomd(11,2),ycomd(11,2),yhemd(11,2),
     &  ynmd(11,2),zhi(21),phi(21,7),dhi(21,7),
     &  thi(21,7),yco2hi(21,7),yn2hi(21,7),yohi(21,7),ycohi(21,7),
     &  yhehi(21,7),ynhi(21,7),yhhi(21,7)
C...  Diurnal time variation parameter
      ftime = dSin(dAtan(1.0d0)*(time - 6.0)/3.0d0)
C...  Latitude factor (suppresses diurnal variation near poles)
      flat = 1.0d0
      alat = dAbs(clat)
      If (alat.ge.70.0d0)flat = 1.0d0 - ((alat-70.0d0)/20.0d0)**2
C...  Mean (A) and diurnal amplitude (B) for temperature
      A = 0.5d0*(tmd(i,2) + tmd(i,1))
      B = 0.5d0*(tmd(i,2) - tmd(i,1))*flat
C...  Time-dependent temperature
      t = A + B*ftime
C...  Mean (A) and diurnal amplitude (B) for pressure
      A = 0.5d0*(dLog(pmd(i,2)*pmd(i,1)))
      B = 0.5d0*(dLog(pmd(i,2)/pmd(i,1)))*flat
C...  Time-dependent pressure
      p = dExp(A + B*ftime)
C...  Gas constant at LST = 0 and 12 hr
      R1 = pmd(i,1)/(dmd(i,1)*tmd(i,1))
      R2 = pmd(i,2)/(dmd(i,2)*tmd(i,2))
C...  Mean (A) and diurnal amplitude (B) for gas constant
      A = 0.5d0*(R2 + R1)
      B = 0.5d0*(R2 - R1)*flat
C...  Time-dependent gas constant
      Rgas = A + B*ftime
C...  Density from perfect gas law
      d = p/(Rgas*t)
C...  Mean (A) and diurnal amplitude (B) for CO2 number density
      A = 0.5d0*(dLog(yco2md(i,2)*yco2md(i,1)))
      B = 0.5d0*(dLog(yco2md(i,2)/yco2md(i,1)))*flat
C...  Time-dependent CO2 number density
      yco2 = dExp(A + B*ftime)
C...  Mean (A) and diurnal amplitude (B) for N2 number density
      A = 0.5d0*(dLog(yn2md(i,2)*yn2md(i,1)))
      B = 0.5d0*(dLog(yn2md(i,2)/yn2md(i,1)))*flat
C...  Time-dependent N2 number density
      yn2 = dExp(A + B*ftime)
C...  Mean (A) and diurnal amplitude (B) for O number density
      A = 0.5d0*(dLog(yomd(i,2)*yomd(i,1)))
      B = 0.5d0*(dLog(yomd(i,2)/yomd(i,1)))*flat
C...  Time-dependent O number density
      yo = dExp(A + B*ftime)
C...  Mean (A) and diurnal amplitude (B) for CO number density
      A = 0.5d0*(dLog(ycomd(i,2)*ycomd(i,1)))
      B = 0.5d0*(dLog(ycomd(i,2)/ycomd(i,1)))*flat
C...  Time-dependent CO number density
      yco = dExp(A + B*ftime)
C...  Mean (A) and diurnal amplitude (B) for He number density
      A = 0.5d0*(dLog(yhemd(i,2)*yhemd(i,1)))
      B = 0.5d0*(dLog(yhemd(i,2)/yhemd(i,1)))*flat
C...  Time-dependent He number density
      yhe = dExp(A + B*ftime)
C...  Mean (A) and diurnal amplitude (B) for N number density
      A = 0.5d0*(dLog(ynmd(i,2)*ynmd(i,1)))
      B = 0.5d0*(dLog(ynmd(i,2)/ynmd(i,1)))*flat
C...  Time-dependent N number density
      yn = dExp(A + B*ftime)
      Return
      End
C----------------------------------------------------------------------
      Subroutine Highterp_V05(i,sza,p,d,t,yco2,yn2,yo,yco,yhe,yn,yh)
C...  Interpolation routine for high altitude VIRA data (150-250 km)
C     Computes pressure, density, temperature (p,d,t), and number
C     densities (yco2,yn2,yo, yco, yhe, yn, and yh) from VIRA high-
C     altitude index (i), and current solar zenith angle (sza)
      Implicit None
      Double Precision sza,p,d,t,zlo,plo,dlo,tlo,zetalo,R1,R2,Rgas,
     &  Rlo,zmd,pmd,dmd,tmd,yco2md,yn2md,yomd,ycomd,yhemd,ynmd,
     &  zhi,phi,dhi,thi,yco2hi,yn2hi,yohi,ycohi,yhehi,ynhi,yhhi,
     &  vsza(7),xsza,yco2,yn2,yo,yco,yhe,yn,yh,fsza
      Integer i,j,jm
C...  Venus model VIRA data at low ,middle, and high altitudes
      Common /VIRAdata_V05/zlo(81),plo(81,5),dlo(81,5),tlo(81,5),
     &  zetalo(81,5),Rlo(81,5),zmd(11),pmd(11,2),dmd(11,2),tmd(11,2),
     &  yco2md(11,2),yn2md(11,2),yomd(11,2),ycomd(11,2),yhemd(11,2),
     &  ynmd(11,2),zhi(21),phi(21,7),dhi(21,7),
     &  thi(21,7),yco2hi(21,7),yn2hi(21,7),yohi(21,7),ycohi(21,7),
     &  yhehi(21,7),ynhi(21,7),yhhi(21,7)
C,,,  Solar zenith angles for high-altitude VIRA data
      Data vsza/15.0d0,34.0d0,61.0d0,90.0d0,119.0d0,146.0d0,165.0d0/
      xsza = dAbs(sza)
C...  Use VIRA data for 0-15 deg if abs(sza) < 15 deg
      If (xsza.le.vsza(1))Then
        p = phi(i,1)
        d = dhi(i,1)
        t = thi(i,1)
        yco2 = yco2hi(i,1)
        yn2 = yn2hi(i,1)
        yo = yohi(i,1)
        yco = ycohi(i,1)
        yhe = yhehi(i,1)
        yn = ynhi(i,1)
        yh = yhhi(i,1)
C...  Use VIRA data for 165 deg if abs(sza) > 165 deg
      Else If (xsza.ge.vsza(7))Then
        p = phi(i,7)
        d = dhi(i,7)
        t = thi(i,7)
        yco2 = yco2hi(i,7)
        yn2 = yn2hi(i,7)
        yo = yohi(i,7)
        yco = ycohi(i,7)
        yhe = yhehi(i,7)
        yn = ynhi(i,7)
        yh = yhhi(i,7)
      Else
C...    Find sza interpolation index j
        Do 10 j = 2,7
          If (xsza.le.vsza(j))Then
            jm = j - 1
C...        Solar zenith angle interpolation factor
            fsza = (xsza - vsza(jm))/(vsza(j) - vsza(jm))
C...        Logarithmic interpolation on pressure
            p = dExp(dLog(phi(i,jm)) + dLog(phi(i,j)/phi(i,jm))*fsza)
C...        Linear interpolation on temperature and gas constant
            t = thi(i,jm) + (thi(i,j)-thi(i,jm))*fsza
            R1 = phi(i,jm)/(dhi(i,jm)*thi(i,jm))
            R2 = phi(i,j)/(dhi(i,j)*thi(i,j))
            Rgas = R1 + (R2-R1)*fsza
C...        Density from perfect gas law
            d = p/(Rgas*t)
C...        Logarithmic interpolation for number densities
            yco2 = dExp(dLog(yco2hi(i,jm)) +
     &        dLog(yco2hi(i,j)/yco2hi(i,jm))*fsza)
            yn2 = dExp(dLog(yn2hi(i,jm)) +
     &        dLog(yn2hi(i,j)/yn2hi(i,jm))*fsza)
            yo=dExp(dLog(yohi(i,jm))+dLog(yohi(i,j)/yohi(i,jm))*fsza)
            yco = dExp(dLog(ycohi(i,jm)) +
     &        dLog(ycohi(i,j)/ycohi(i,jm))*fsza)
            yhe = dExp(dLog(yhehi(i,jm)) +
     &        dLog(yhehi(i,j)/yhehi(i,jm))*fsza)
            yn=dExp(dLog(ynhi(i,jm))+dLog(ynhi(i,j)/ynhi(i,jm))*fsza)
            yh=dExp(dLog(yhhi(i,jm))+dLog(yhhi(i,j)/yhhi(i,jm))*fsza)
            Return
          Endif
  10    Enddo
      Endif
      Return
      End
C----------------------------------------------------------------------
      Subroutine Thermos_V05(z,clat,zf,yco2f,yn2f,yof,ycof,yhef,ynf,
     &  yhf,Tf,presz,densz,AMz,yco2,yn2,yo,yco,yhe,yn,yh,totnd,hscale)
C...  Special Venus thermosphere model asssuming constant exospheric
C     temperature above height zf (250 km)
C...  Evaluates pressure (presz), density (densz), molecular weight
C     (AMz), scale height (hscale), and number densities yxxx for 7
C     species, given conditions at height zf (= 250 km), for altitude
C     (z) and latitude (clat)
      Implicit None
      Double Precision z,zf,yco2f,yn2f,yof,ycof,yhef,ynf,yhf,Tf,presz,
     &  densz,AMz,yco2,yn2,yo,yco,yhe,yn,yh,totnd,amh,amco2,amhe,amn2,
     &  amo,amco,amn,zzf,bk,R0,Rf,Rsfc,gf,clat,y,ppco2f,ppn2f,ppof,
     &  ppcof,pphef,ppnf,pphf,Tinf,Tz,Hh,Hco2,Hhe,Hn2,Ho,Hco,Hn,ymin,
     &  pph,ppco2,pphe,ppn2,ppo,ppco,ppn,hscale,gz
C...  Molecular weights for atmospheric constituents
      Data amh,amco2,amhe,amn2,amo,amco,amn/1.01d0,44.0d0,4.0d0,28.0d0,
     &  16.0d0,28.0d0,14.0d0/
C...  Height above 250 km
      zzf = z - zf
C...  Boltzmann constant
      bk = 1.38065D-23
C...  Universal gas constant
      R0 = 8314.472D0
C...  Surface radius (Rsfc) and gravity (gf) at height zf
      Call rgplanet_V05(clat,Rsfc,zf,gf,2)
C...  Radius at height zf
      Rf = Rsfc + zf
C...  Partial pressures of constituents at height zf
      ppco2f = bk*yco2f*Tf
      ppn2f = bk*yn2f*Tf
      ppof = bk*yof*Tf
      ppcof = bk*ycof*Tf
      pphef = bk*yhef*Tf
      ppnf = bk*ynf*Tf
      pphf = bk*yhf*Tf
C...  Height parameter
      y = zzf*Rf/(Rf + zzf)
C...  Assume constant temperature versus height
      Tinf = tf
      Tz = tf
C...  Get scale height, partial pressure, and number density of each
C     constituent
      Hco2 = R0*Tinf/(1000.0d0*amco2*gf)
      ppco2 = ppco2f*Exp(-y/Hco2)
      yco2 = ppco2/(bk*Tz)
      Hn2 = R0*Tinf/(1000.0d0*amn2*gf)
      ppn2 = ppn2f*Exp(-y/Hn2)
      yn2 = ppn2/(bk*Tz)
      Ho = R0*Tinf/(1000.0d0*amo*gf)
      ppo = ppof*Exp(-y/Ho)
      yo = ppo/(bk*Tz)
      Hco = R0*Tinf/(1000.0d0*amco*gf)
      ppco = ppcof*Exp(-y/Hco)
      yco = ppco/(bk*Tz)
      Hhe = R0*Tinf/(1000.0d0*amhe*gf)
      pphe = pphef*Exp(-y/Hhe)
      yhe = pphe/(bk*Tz)
      Hn = R0*Tinf/(1000.0d0*amn*gf)
      ppn = ppnf*Exp(-y/Hn)
      yn = ppn/(bk*Tz)
      Hh = R0*Tinf/(1000.0d0*amh*gf)
      pph = pphf*Exp(-y/Hh)
      yh = pph/(bk*Tz)
C...  Limit number densities to > 1 per km**3
      ymin = 1.0d-9
      If (yco2.lt.ymin)yco2=ymin
      If (yn2.lt.ymin)yn2=ymin
      If (yco.lt.ymin)yco=ymin
C...  Get total pressure from partial pressures
      presz = ppco2 + ppn2 + ppo + ppco + pphe + ppn + pph
C...  Get total number density from constituent number densities
      totnd = yco2 + yn2 + yo + yco + yhe + yn + yh
C...  Get mean molecular weight
      Amz = (yco2*amco2+yn2*amn2+yo*amo+yco*amco+yhe*amhe+yn*amn+
     &  yh*amh)/totnd
C...  Get density from perfect gas law
      densz = presz*Amz/(R0*Tz)
C...  Get gravity at height z
      Call rgplanet_V05(clat,Rsfc,z,gz,2)
C...  Get pressure scale height
      hscale = R0*Tz/(1000.0d0*AMz*gz)
      Return
      End
C----------------------------------------------------------------------
        Subroutine ProfTerp_V05(chgt,clat,clon,tin,pin,din,uin,vin,
     &    ptemp,ppres,pdens,puwin,pvwin,nprof,profnear,proffar,profwgt)
C---    Interpolates profile data to current position (chgt,clat,clon)
C       and weights results (with factor profwgt) with input values
C       (tin,pin,din,uin,vin), yielding weighted average (ptemp,ppres,
C       pdens,puwin,pvwin).  Input profnear is lat-lon radius over
C       which profile is weighted with 1.0; proffar is lat-lon radius
C       beyond which profile is given zero weight.
        Implicit Double Precision (A-H,O-Z)
        Parameter (npmax = 100000)
        Common /pterp_V05/phgt(npmax),plat(npmax),plon(npmax),
     &    ptmp(npmax),pprs(npmax),pden(npmax),puwn(npmax),pvwn(npmax)

        i1 = 0 ! gg avoids warnings about i1 and i2
        i2 = 0
C---    Calculate pi/2
        pi2 = 2.0d0*dAtan(1.0d0)
        If (chgt.lt.phgt(1))Then
C---    Profile weighting zero for height below 1st profile data point
          profwgt = 0.0d0
          pilat = plat(1)
          pilon = plon(1)
          ptemp = ptmp(1)
          ppres = pprs(1)
          pdens = pden(1)
          puwin = puwn(1)
          pvwin = pvwn(1)
        Else If (chgt.gt.phgt(nprof))Then
C---    Profile weighting zero for height above last profile data point
          profwgt = 0.0d0
          pilat = plat(nprof)
          pilon = plon(nprof)
          ptemp = ptmp(nprof)
          ppres = pprs(nprof)
          pdens = pden(nprof)
          puwin = puwn(nprof)
          pvwin = pvwn(nprof)
        Else
C---    Find index values i1 and i2=i1+1 bracketing current height
          Do 20 i = 1,nprof-1
            If (chgt.ge.phgt(i).and.chgt.le.phgt(i+1))Then
              i1 = i
              i2 = i + 1
              Goto 25
            Endif
   20     Enddo
C---      Compute factor for linear height interpolation
   25     factor = (chgt - phgt(i1))/(phgt(i2)-phgt(i1))
C---      Linear height interpolation for lat,lon,temperature,winds
          pilat = plat(i1) + factor*(plat(i2)-plat(i1))
          pilon = plon(i1) + factor*(plon(i2)-plon(i1))
          ptemp = ptmp(i1) + factor*(ptmp(i2)-ptmp(i1))
          puwin = puwn(i1) + factor*(puwn(i2)-puwn(i1))
          pvwin = pvwn(i1) + factor*(pvwn(i2)-pvwn(i1))
C---      Power-law interpolation for pressure (unless profile pressure
C         is zero, for which zero weight will be used)
          pdens = 0.0d0
          If (pden(i1).gt.0.0d0)pdens = pden(i1)*
     &       (pden(i2)/pden(i1))**factor
C---      Power-law interpolation for density (unless profile density
C         is zero, for which zero weight will be used)
          ppres = 0.0d0
          If (pprs(i1).gt.0.0d0)ppres = pprs(i1)*
     &      (pprs(i2)/pprs(i1))**factor
C---      Initialize weighting factor components for height and lat-lon
          facthgt = 1.0d0
          factll = 1.0d0
          If (chgt.le.phgt(2))Then
C---      Sin-squared variation of height weighting from 0 at 1st point
C         to 1 at 2nd point
            facthgt=(chgt-phgt(1))/(phgt(2)-phgt(1))
            facthgt = (dSin(pi2*facthgt))**2
          Else If (chgt.ge.phgt(nprof-1))Then
C---      Sin-squared variation of height weighting from 0 at next-to-
C         last point to 1 at last point
            facthgt=(chgt-phgt(nprof))/(phgt(nprof-1)-phgt(nprof))
            facthgt = (dSin(pi2*facthgt))**2
          Endif
C---      Compute absolute lat-lon difference of current position from
C         profile lat-lon
          dlat = dAbs(clat - pilat)
          dlon = dAbs(clon - pilon)
C---      Adjust lon difference for wrap at lon 360
          If (dlon.gt.180.0d0)dlon = 360.0d0 - dlon
C---      Lat-lon radius of current position from profile lat-lon
          radius = dSqrt(dlat**2 + dlon**2)
C---      Use weight=0 if radius>proffar, weight=1 if radius<profnear,
C         with sin-squared variation between proffar and profnear
          If (radius.ge.proffar)Then
            factll = 0.0d0
          Else If (radius.le.profnear)Then
            factll = 1.0d0
          Else
            factll = (proffar-radius)/(proffar - profnear)
            factll = (dSin(pi2*factll))**2
          Endif
C---      Total weight = product of weights for lat-lon and height
          profwgt = factll*facthgt
        Endif
        tpdwgt = profwgt
        uvwgt = profwgt
C---    Set profile weight to zero for p,d, & t if profile values are 0
        If (ptemp*ppres*pdens.eq.0.0d0)tpdwgt=0.0d0
C---    Set profile weight to zero for u & v if profile values are 0
        If (dAbs(puwin)+dAbs(pvwin).eq.0.0d0)uvwgt = 0.0d0
C---    Apply weighted averaging of profile values with input values
        ptemp = tpdwgt*ptemp + (1.0d0 - tpdwgt)*tin
        ppres = tpdwgt*ppres + (1.0d0 - tpdwgt)*pin
        pdens = tpdwgt*pdens + (1.0d0 - tpdwgt)*din
        puwin = uvwgt*puwin + (1.0d0 - uvwgt)*uin
        pvwin = uvwgt*pvwin + (1.0d0 - uvwgt)*vin
        Return
        End
C----------------------------------------------------------------------
        Subroutine RdProf_V05(profile,nprof,LonEast)
C---    Reads alternate profile data file profile. Returns number of
C       lines of data (nprof).  Converts input longitudes from East to
C       West if LonEast = 1
        Implicit Double Precision (A-H,O-Z)
        Character*60 profile
        Character*1 dummy
        Parameter (npmax = 100000)
        Common /pterp_V05/phgt(npmax),plat(npmax),plon(npmax),
     &    ptmp(npmax),pprs(npmax),pden(npmax),puwn(npmax),pvwn(npmax)
C---    Compute string length for profile file name
        lenprof = index(profile,' ')-1
        If (lenprof.lt.1.or.lenprof.gt.60)lendir = 60
C---    Open profile data file
        Open(33,file=profile(1:lenprof),status='old')
C---    Read and ignore header line
        Read(33,5)dummy
   5    Format(A1)
        n = 0
C---    Start of loop to read profile data
  10    Read(33,*,End=99)zhgt,xlat,xlon,t,p,d,u,v
C---    Convert negative longitudes
        If (xlon.lt.0.0d0)xlon = xlon + 360.0d0
C---    Convert to West Longitude if LonEast = 1
        If (LonEast.eq.1)xlon = 360.0d0 - xlon
C---    Count number of lines read
        n = n + 1
C---    Store profile data in arrays, for common pterp_V05
        phgt(n) = zhgt
C---    Stop if two successive heights are the same
        If (n.gt.1.and.phgt(n).eq.phgt(n-1))
     &    Stop ' Consecutive profile heights cannot be same'
        plat(n) = xlat
        plon(n) = xlon
        ptmp(n) = t
        pprs(n) = p
        pden(n) = d
        puwn(n) = u
        pvwn(n) = v
        nprof = n
C---    Cycle back to read another line of profile data
        Goto 10
C---    Close profile input file when end-of-file encountered
  99    Close(33)
        Return
        End
C----------------------------------------------------------------------
      Subroutine ReadVIRA_V05(DATADIR)
C...  Reads VIRA data tables and stores in common block VIRAdata_V05.
C     Computes reference VIRA data values and stores in common block
C     VIRAref_V05.  DATADIR is first part of path to VIRA data
C     directory.
      Implicit None
      Integer lendir,i,j
      Character*300 DATADIR
      Character*1 dummy
      Double Precision zlo,plo,dlo,tlo,zetalo,Rlo,zmd,pmd,dmd,tmd,
     &  yco2md,yn2md,yomd,ycomd,yhemd,ynmd,zhi,phi,dhi,xlat,sza,xlst,
     &  thi,yco2hi,yn2hi,yohi,ycohi,yhehi,ynhi,yhhi,z,p,d,t,Rgas,
     &  zeta,yco2,yn2,yo,yco,yhe,yn,yh,avmw,R0,AvN,zref,pref,
     &  dref,tref
C...  Common block for VIRA data arrays
C     Low altitude data (0-100 km) indexed 81 by 5, with 81 heights
C     and 5 latitudes (30,45,60,75,85).
C     Middle altitude data (100-150 km) indexed 11 by 2, with 11
C     heights and 2 times (LST = 0, 12).
C     High altitude data (150-250 km) indexed 21 by 7, with 21
C     heights and 7 solar zenith angles (15,34,61,90,119,146,165).
C     Low, middle, high altitude data designated lo, md, hi, with -
C     z = height (km)
C     p = pressure (N/m**2)
C     d = density (kg/m**3)
C     t = temperature (K)
C     R = gas constant (SI units)
C     zeta = compressibility factor [ = p/(d*R*t) ]
C     yco2 = CO2 number density (#/m**3)
C     yn2 = N2 number density (#/m**3)
C     yo = O number density (#/m**3)
C     yco = CO number density (#/m**3)
C     yhe = He number density (#/m**3)
C     yn = N number density (#/m**3)
C     yh = H number density (#/m**3)
      Common /VIRAdata_V05/zlo(81),plo(81,5),dlo(81,5),tlo(81,5),
     &  zetalo(81,5),Rlo(81,5),zmd(11),pmd(11,2),dmd(11,2),tmd(11,2),
     &  yco2md(11,2),yn2md(11,2),yomd(11,2),ycomd(11,2),yhemd(11,2),
     &  ynmd(11,2),zhi(21),phi(21,7),dhi(21,7),
     &  thi(21,7),yco2hi(21,7),yn2hi(21,7),yohi(21,7),ycohi(21,7),
     &  yhehi(21,7),ynhi(21,7),yhhi(21,7)
C...  Common block for VIRA reference data arrays
      Common /VIRAref_V05/zref(111),pref(111),dref(111),tref(111)
C...  Compute character string length of DATADIR path name
      lendir = Index(DATADIR,' ')-1
      If (lendir.lt.1.or.lendir.gt.300)lendir = 300
      Open(26,file=DATADIR(1:lendir)//'VIRALow.txt',status='old')
      Open(27,file=DATADIR(1:lendir)//'VIRAMid.txt',status='old')
      Open(28,file=DATADIR(1:lendir)//'VIRAHi.txt',status='old')
C...  Read and ignore header information
      Read(26,10)dummy
      Read(27,10)dummy
      Read(28,10)dummy
  10  Format(A1)
      R0 = 8314.472d0
      AvN = 6.02214199d26
C.... Read low altitude (0-100 km) VIRA data and store in Common arrays
      Do 30 j = 1,5
      Do 20 i = 1,81
        Read(26,*)z,xlat,d,p,t,zeta,Rgas
        avmw = R0/Rgas
C...    Store low altitude VIRA data in arrays
        zlo(i) = z
        plo(i,j) = zeta*d*Rgas*t
        dlo(i,j) = d
        tlo(i,j) = t
        zetalo(i,j) = zeta
        Rlo(i,j) = Rgas
  20  Enddo
  30  Enddo
C.... Read middle altitude (100-150 km) VIRA data and store in Common
C     arrays
      Do 50 j = 1,2
      Do 40 i = 1,11
        Read(27,*)z,xlst,yco2,yn2,yo,yco,yhe,yn,d,p,t,avmw
C...    Store middle altitude VIRA data in arrays
        zmd(i) = z
        pmd(i,j) = p
        dmd(i,j) = d
        tmd(i,j) = t
        yco2md(i,j) = yco2
        yn2md(i,j) = yn2
        yomd(i,j) = yo
        ycomd(i,j) = yco
        yhemd(i,j) = yhe
        ynmd(i,j) = yn
  40  Enddo
  50  Enddo
C.... Read high altitude (150-250 km) VIRA data; store in Common arrays
      Do 70 j = 1,7
      Do 60 i = 1,21
        Read(28,*)z,sza,yco2,yn2,yo,yco,yhe,yn,yh,d,p,t,avmw
C...    Store high altitude VIRA data in arrays
        zhi(i) = z
        phi(i,j) = p
        dhi(i,j) = d
        thi(i,j) = t
        yco2hi(i,j) = yco2
        yn2hi(i,j) = yn2
        yohi(i,j) = yo
        ycohi(i,j) = yco
        yhehi(i,j) = yhe
        ynhi(i,j) = yn
        yhhi(i,j) = yh
   60 Enddo
   70 Enddo
C...  Close VIRA data files
      Close(26)
      Close(27)
      Close(28)
C...  Save reference values of pressure, density, temperature
C...  Reference values for low altitudes (0-98 km) (0-30 lat data)
      Do 80 i = 1,80
        zref(i) = zlo(i)
        pref(i) = plo(i,1)
        dref(i) = dlo(i,1)
        tref(i) = tlo(i,1)
  80  Enddo
C...  Reference values for 100 km altitude (average from low and
C     middle altitude data sets at 100 km)
      zref(81) = zmd(1)
      pref(81) = dSqrt(dSqrt(pmd(1,1)*pmd(1,2))*plo(81,1))
      dref(81) = dSqrt(dSqrt(dmd(1,1)*dmd(1,2))*dlo(81,1))
      tref(81) = 0.5d0*(0.5d0*(tmd(1,1)+tmd(1,2)) + tlo(81,1))
C...  Reference values for middle altitudes (105-145 km) (average
C     of LST = 0 and LST = 12)
      Do 90 i = 2,10
        zref(80+i) = zmd(i)
        pref(80+i) = dSqrt(pmd(i,1)*pmd(i,2))
        dref(80+i) = dSqrt(dmd(i,1)*dmd(i,2))
        tref(80+i) = 0.5d0*(tmd(i,1)+tmd(i,2))
  90  Enddo
C...  Reference values for 150 km (average of middle and high altitude
C     data sets at 150 km)
      zref(91) = zhi(1)
      pref(91) = dSqrt(phi(1,4)*dSqrt(pmd(11,1)*pmd(11,2)))
      dref(91) = dSqrt(dhi(1,4)*dSqrt(dmd(11,1)*dmd(11,2)))
      tref(91) = 0.5d0*(thi(1,4) + 0.5d0*(tmd(11,1)+tmd(11,2)))
C...  Reference values for high altitudes (155-250 km) (solar zenith
C     angle 90 data from high altitude range)
      Do 100 i = 2,21
        zref(90+i) = zhi(i)
        pref(90+i) = phi(i,4)
        dref(90+i) = dhi(i,4)
        tref(90+i) = thi(i,4)
 100  Enddo
      Return
      End
C---------------------------------------------------------------------- 
