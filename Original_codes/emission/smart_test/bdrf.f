      REAL FUNCTION BDREF(MU, MUP, DPHI, SURF_PR, IREF)

c      Supplies surface bi-directional reflectivity.
c
c    REFMODE  : bidirectional reflectance options
c             0 - Lambert
c             1 - Hapke's BDR model
c             2 - Breon's BDR model; combination of Li + Roujean
c             3 - Roujean's BDR model
c             4 - Cox and Munk glint model
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MUP (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
c
c    MU     : Cosine of angle of reflection (positive)
c
c    MUP    : Cosine of angle of incidence (positive)
c
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c   SURF_PR : Wavelength dependent surface properties array 
c             IREF= 0 - Lambert albedo
c             IREF= 1 - Hapke : single scatter albedo, W, and 
c                               angular width factr HH 
c             IREF= 2 - Breon's BDR model: k0, k1, k2
c             IREF= 3 - Roujean's BDR model: k0, k1, k2
c             IREF= 4 - Cox and Munk glint model: n, k, ws, phiw
c
c  LOCAL VARIABLES:
c
c    IREF   : bidirectional reflectance options
c
c    B0     : empirical factor to account for the finite size of
c             particles in Hapke's BDR model
c
c    B      : term that accounts for the opposition effect
c             (retroreflectance, hot spot) in Hapke's BDR model
c
c    CTHETA : cosine of phase angle in Hapke's BDR model
c
c    GAMMA  : albedo factor in Hapke's BDR model
c
c    H0     : H( mu0 ) in Hapke's BDR model
c
c    H      : H( mu ) in Hapke's BDR model
c
c    HH     : angular width parameter of opposition effect in Hapke's
c             BDR model
c
c    P      : scattering phase function in Hapke's BDR model
c
c    THETA  : phase angle (radians); the angle between incidence and
c             reflection directions in Hapke's BDR model
c
c    W      : single scattering albedo in Hapke's BDR model
c
c    ws     : wind speed (m/s) in the Cox Munk model
c
c    phiw   : wind direction (radians) in the Cox Munk model
c
c   Called by- DREF, SURFAC
c +-------------------------------------------------------------------+
c
      implicit none
c
      REAL PI
      PARAMETER (PI=3.14159265358979323)
c
c     .. Scalar Arguments ..
c
      INTEGER   IREF
      REAL      DPHI, MU, MUP, SURF_PR(4)
      REAL      th_s,th_v,K0,K1,K2,f1_2,f1_3,f2
      REAL      nr,ni,wspd,azw,ts,tv,fi,rog,mup1
c
c     .. Local Scalars ..
c
      REAL      B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W
c
c     .. Intrinsic Functions ..
c
      INTRINSIC COS, SQRT
c
      IF ( IREF.EQ.1 ) THEN

c        ** Hapke's BRDF model (times Pi/Mu0)
c        ** (Hapke, B., Theory of reflectance
c        ** and emittance spectroscopy, Cambridge
c        ** University Press, 1993, Eq. 8.89 on
c        ** page 233. Parameters are from
c        ** Fig. 8.15 on page 231, expect for w.)

         CTHETA = MU * MUP + (1.-MU**2)**.5 * (1.-MUP**2)**.5
     &            * COS( DPHI )
         THETA = ACOS( CTHETA )
         P    = 1. + 0.5 * CTHETA
         HH   = SURF_PR(2) ! 0.06
         B0   = 1.0
         B    = B0 * HH / ( HH + TAN( THETA/2.) )
         W = SURF_PR(1) ! 0.6
         GAMMA = SQRT( 1. - W )
         H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
         H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )
         BDREF = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      END IF
c
      IF ( IREF.EQ.2 ) THEN
       k0 =  SURF_PR(1) ! 0.3530000
       k1 =  SURF_PR(2) ! 5.7000000E-02
       k2 =  SURF_PR(3) ! 0.3590000
       th_v = acos(mu)
       th_s = acos(mup)
c
       BDREF = K0 + K1*F1_2(th_s,th_v,dphi) + K2*F2(th_s,th_v,dphi)
c
      END IF
c
      IF ( IREF.EQ.3 ) THEN
c
c             Roujean's BDR model
c
       k0 =  SURF_PR(1) ! 0.3530000
       k1 =  SURF_PR(2) ! 5.7000000E-02
       k2 =  SURF_PR(3) ! 0.3590000
       th_v = acos(mu)
       th_s = acos(mup)

       BDREF = K0 + K1*F1_3(th_s,th_v,dphi) + K2*F2(th_s,th_v,dphi)
c
      END IF
c
      IF ( IREF.EQ.4 ) THEN
c
c             Cox and Munk glint model from 6s
C input parameters:   wspd=speed of the wind (in m/s)
C                     nr=index of refraction of the sea water
C                     ni=extinction coefficient of the sea water
c                     azw=azim. of the sun - azim. of the wind (in deg.)
C                     ts=solar zenith angle (in deg.)
C                     tv=view zenith angle (in deg.)
C                     fi=relative azimuth (sun-viewing angle)
C output parameters:  rog=reflectance of the sun glint
c
c*****  note, for near glancing angles (mup < 0.03) this routine
c       produces flux albedo values > 1.0.  To fix this, set all
c       near glancing agles to mup > 0.3
c
       if(mup .ge. 0.08) then
         mup1 = mup
       else
         mup1 = 0.08
       endif
       nr = SURF_PR(1)
       ni = SURF_PR(2)
       wspd = SURF_PR(3)
       azw = SURF_PR(4)
       ts = 180.0*acos(mup1)/PI
       tv = 180.0*acos(mu)/PI
       fi = 180.0*dphi/PI
c
       call sunglint(wspd,nr,ni,azw,ts,tv,fi,rog)
c       write(*,'(10(1pe12.4))') wspd,nr,ni,azw,ts,tv,fi,rog
       if(rog .lt. 0.0) rog = 0.0
       BDREF = rog
c

      END IF

      END ! FUNCTION BDREF

c
      REAL FUNCTION F1_2(THETA_S, THETA_V, PHI )
c
c     Li-Sparse F1
c
c  INPUT:
c
c    THETA_S : Cosine of angle of incidence (positive)
c
c    THETA_V : Cosine of angle of reflection (positive)
c
c    PHI     : azimuth angles of incidence and reflection (radians)
c
c +-------------------------------------------------------------------+
c
c     .. Scalar Arguments ..
c
      IMPLICIT NONE
      REAL PI
      PARAMETER (PI=3.14159265358979323)

      REAL PHI, THETA_S, THETA_V, INVPI, X, DUMY
      REAL Q, t, KK, cos_ksi, costt

      DUMY=70 * 1.74533E-2
      IF(THETA_S .GT. DUMY)  THETA_S=DUMY
      IF(THETA_V .GT. DUMY)  THETA_V=DUMY
      INVPI=1/PI

      X = TAN(THETA_S)**2+TAN(THETA_V)**2
     >  - 2*TAN(THETA_S)*TAN(THETA_V)*COS(PHI)
       cos_ksi=COS(THETA_S)*COS(THETA_V)
     >+ SIN(THETA_S)*SIN(THETA_V)*COS(PHI)
       KK = ACOS(COS(THETA_S)*COS(THETA_V)
     >  + SIN(THETA_S)*SIN(THETA_V)*COS(PHI))
       KK=ACOS(cos_ksi)

      costt = 2.*sqrt(X+(TAN(THETA_S)*TAN(THETA_V)*sin(PHI))**2)/
     >  (1./cos(THETA_S)+1./cos(THETA_V))

      If(costt .gt. 1) then
       costt =1
      endif
      t = ACOS(costt)

      Q = INVPI*(t-sin(t)*cos(t))*(1./cos(THETA_S)+1./cos(THETA_V))

      F1_2 = Q - 1./cos(THETA_V) - 1./cos(THETA_S) +0.5*(1+cos(KK))
     >    *(1./cos(THETA_S))*(1./cos(THETA_V))

      END ! FUNCTION F1_2

      REAL FUNCTION F1_3(THETA_S, THETA_V, PHI )

c     Roujean F1
c
c  INPUT:
c
c    THETA_S : Cosine of angle of incidence (positive)
c
c    THETA_V : Cosine of angle of reflection (positive)
c
c    PHI     : azimuth angles of incidence and reflection (radians)
c
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..
      IMPLICIT NONE
      REAL PI
      PARAMETER (PI=3.14159265358979323)

      REAL PHI, THETA_S, THETA_V, INV2PI, INVPI, X, DUMY

      DUMY=70 * 1.74533E-2
      IF(THETA_S .GT. DUMY)  THETA_S=DUMY
      IF(THETA_V .GT. DUMY)  THETA_V=DUMY
      INVPI=1/PI
      INV2PI=1/(2*PI)

      X = TAN(THETA_S)**2+TAN(THETA_V)**2
     >  - 2*TAN(THETA_S)*TAN(THETA_V)*COS(PHI)
      F1_3=INV2PI*((PI-PHI)*COS(PHI)+SIN(PHI))*TAN(THETA_S)
     >  *TAN(THETA_V) - INVPI*(TAN(THETA_S)+TAN(THETA_V)+SQRT(X))

       END ! FUNCTION F1_3

      REAL FUNCTION F2(THETA_S, THETA_V, PHI )

c     Roujean F2
c
c  INPUT:
c
c    THETA_S : Cosine of angle of incidence (positive)
c
c    THETA_V : Cosine of angle of reflection (positive)
c
c    PHI     : azimuth angles of incidence and reflection (radians)
c
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      IMPLICIT NONE
      REAL PI
      PARAMETER (PI=3.14159265358979323)
c
      REAL PHI, THETA_S, THETA_V, X, DUMY
      REAL COS_KSI 

      DUMY=70 * 1.74533E-2
      IF(THETA_S .GT. DUMY) THETA_S=DUMY
      IF(THETA_V .GT. DUMY) THETA_V=DUMY

      cos_ksi=COS(THETA_S)*COS(THETA_V)
     >+ SIN(THETA_S)*SIN(THETA_V)*COS(PHI)
      X = ACOS(COS(THETA_S)*COS(THETA_V)
     >  + SIN(THETA_S)*SIN(THETA_V)*COS(PHI))
      X=ACOS(cos_ksi)
       F2 = 4/(3*PI*(COS(THETA_S)+COS(THETA_V)))
     >   * ((PI/2-X)*COS(X)+SIN(X))-1./3.

      END ! FUNCTION F2

      subroutine sunglint(wspd,nr,ni,azw,ts,tv,fi,rog)
C
C input parameters:   wspd=speed of the wind (in m/s)
C                     nr=index of refraction of the sea water
C                     ni=extinction coefficient of the sea water
c                     azw=azim. of the sun - azim. of the wind (in deg.)
C                     ts=solar zenith angle (in deg.)
C                     tv=view zenith angle (in deg.)
C                     fi=relative azimuth (sun-satellite)
C output parameters:  rog=reflectance of the sun glint
C
      implicit none
c
      real r1
      real coschi,sinchi
      real wspd,nr,ni,azw,ts,tv,fi,rog
      double precision pi,fac,cs,cv,ss,sv,phi,zx,zy,tantilt,tilt,proba,
     -                 xe,xn,xe2,xn2,phw
      double precision coef,cos2chi
      double precision sigmaC,sigmaU,C21,C03,C40,C04,C22
c
      pi=dacos(-1.0d0)
      fac=pi/180.0d0
      phw=azw*fac
      cs=dcos(ts*fac)
      cv=dcos(tv*fac)
      ss=dsin(ts*fac)
      sv=dsin(tv*fac)
      phi=fi*fac
      Zx=-sv*sin(phi)/(cs+cv)
      Zy=(ss+sv*cos(phi))/(cs+cv)
      tantilt=dsqrt(zx*zx+zy*zy)
      tilt=datan(tantilt)
c
c  Anisotropic Gaussian distribution
c    phw=phi_sun-phi_wind
c
      sigmaC=0.003d0+0.00192*wspd
      sigmaU=0.00316d0*wspd
      C21=0.01d0-0.0086d0*wspd
      C03=0.04d0-0.033d0*wspd
      C40=0.40d0
      C22=0.12d0
      C04=0.23d0
      xe=(dcos(phw)*Zx+dsin(phw)*Zy)/dsqrt(SigmaC)
      xn=(-dsin(phw)*Zx+dcos(phw)*Zy)/dsqrt(SigmaU)
      xe2=xe*xe
      xn2=xn*xn
      coef=1-C21/2.*(xe2-1)*xn-C03/6.*(xn2-3)*xn
      coef=coef+c40/24.0d0*(xe2*xe2-6*xe2+3)
      coef=coef+C04/24.0d0*(xn2*xn2-6*xn2+3)
      coef=coef+C22/4.0d0*(xe2-1)*(xn2-1)
      proba=coef/2.0d0/pi/sqrt(sigmaU)/dsqrt(sigmaC)*
     -     dexp(-(xe2+xn2)/2.0d0)
c
c Compute Fresnel's coefficient R1
c
      cos2chi=cv*cs+sv*ss*dcos(phi)
      if (cos2chi.gt.1.0)cos2chi=0.99999999999
      if (cos2chi.lt.-1.0)cos2chi=-0.99999999999
      coschi=dsqrt(0.50d0*(1.0d0+cos2chi))
      sinchi=dsqrt(0.50d0*(1.0d0-cos2chi))
c
      Call Fresnel(nr,ni,coschi,sinchi,r1)
c
C Compute Reflectance of the sun glint
c
      Rog=0.25*pi*r1*proba/cs/cv/(dcos(tilt)**4)
c      write(*,'(10(1pe12.4))') pi,r1,proba,cs,cv,tilt,cos(tilt),Rog
c
      return
      end
C
C
      Subroutine Fresnel(nr,ni,coschi,sinchi,R1)
C
C to compute the Fresnel's coefficient of reflection (see for
C example M. Born and E. Wolf, Principles of Optics, Pergamon Press, fifth
C edition, 1975, pp 628
C input parameters: nr=index of refraction of the sea water
C                   ni=extinction coefficient of the sea water
C                   coschi & sinchi=cosine and sine of the incident radiation
C                                   with respect of the wave facet normal.
C output parameter: R1=Fresnel's coefficient for reflection
C
      real nr,ni,a1,a2,u,v,Rr2,Rl2,b1,b2,R1,coschi,sinchi
c absolute value for a1 to get v=0 when ni=0
      a1=abs(nr*nr-ni*ni-sinchi*sinchi)
      a2=sqrt((nr*nr-ni*ni-sinchi*sinchi)**2.+4*nr*nr*ni*ni)
      u=sqrt(0.5*(a1+a2))
      v=sqrt(0.5*(-a1+a2))
      Rr2=((coschi-u)**2+v*v)/((coschi+u)**2+v*v)
      b1=(nr*nr-ni*ni)*coschi
      b2=2*nr*ni*coschi
      Rl2=((b1-u)**2+(b2+v)**2)/((b1+u)**2+(b2-v)**2)
      R1=(Rr2+Rl2)/2.
      return
      end
C
