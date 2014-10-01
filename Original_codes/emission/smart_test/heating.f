      subroutine heating(iuheat,nlyr,nza,umu0_1,wnmin,wnmax,p,t,alt,
     -                    aid_lr,dirsoflx,dnsoflx,upsoflx,
     -                    dnthflx,upthflx,soheat,thheat)
c
cccccccccccccccccccccccccccc  h e a t i n g  ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes flux divergences adn heating rates     cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc      nlyr - number of of layers in the model atmosphere            cc
cc       nza - number of solar zenith angles for which heating rates  cc
cc             are needed.                                            cc
cc     flxdn - spectrally-integrated downward flux (W/m**2).          cc
cc     flxup - spectrally-integrated upward flux (W/m**2).            cc
cc    aid_lr - aidiabatic lapse rate (g/cp) in mks units              cc
cc         p - pressure at each level (pascals)                       cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      soheat - solar radiative heating rate (K/sec)                 cc
cc      thheat - thermal radiative heating rate (K/sec)               cc
cc                                                                    cc
cccccccccccccccccccccccccccc  h e a t i n g  ccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer iuheat,nlyr,nza
      integer nz,k
c
      real umu0_1(nsol)
      double precision soflxn(mxulv),thflxn(mxulv)
c
c****   atmospheric structure variables
c
      real p(kp),t(kp),alt(kp)
c
c****   double precision variables for heating rate integration
c
      double precision wnmin,wnmax
      double precision aid_lr(mxulv),dirsoflx(mxulv,nsol),
     -       dnsoflx(mxulv,nsol),upsoflx(mxulv,nsol),dnthflx(mxulv),
     -       upthflx(mxulv),soheat(mxulv,nsol),thheat(mxulv)
c
      real pi,pbar,tbar,altkm,soheat_d,thheat_d,solang
c
      data pi /3.141593 /
c
c****     define solar heating rates at each solar zenith angle
c
      do 1201 nz=1,nza
c
c****      define the net flux at each level:
c
          do 1001 k=1,nlyr+1
              soflxn(k) = dirsoflx(k,nz) + 
     -                    dnsoflx(k,nz) - upsoflx(k,nz)
1001      continue 
c
c****       define the heating rate in each layer
c
          do 1021 k=1,nlyr
c
              soheat(k,nz) = aid_lr(k)*
     -                       (soflxn(k) - soflxn(k+1))/(p(k+1) - p(k))
1021      continue
1201  continue
c
c****    define the net flux at each level:
c
      do 2001 k=1,nlyr+1
              thflxn(k) = dnthflx(k) - upthflx(k)
2001  continue 
c
c***** define thermal cooling rates (independent of zenith angle)
c
      do 2021 k=1,nlyr
          thheat(k) = aid_lr(k)*
     -                (thflxn(k+1) - thflxn(k))/(p(k+1) - p(k))
2021  continue
      rewind(iuheat)
      write(iuheat,'(/,/,1a)')
     - ' R a d i a t i v e     H e a t i n g    R a t e s '
c
      do 6221 nz=1,nza
          solang = 180.*acos(umu0_1(nz))/pi
          write(iuheat,
     -     '(/,/,1a,2(1pe14.6,1a),/,1a,1pe13.5,1a,/,/,3a,/,3a)') 
     -     'Spectral Range =', wnmin,' -',wnmax,' wavenumbers',
     -     'Solar Zenith Angle =', solang,' degrees',
     -     '   pressure  temperature   altitude     solar Q ',  
     -     '   thermal Q  dir sol flx  dn sol flx   up sol flx',
     -     '   dn th flx   up th flx     p(flux)     t(flux)',
     -     '     (bar)       (K)         (km)      (K/day) ', 
     -     '    (K/day)   (watts/m*m)  (watts/m*m)  (watts/m*m)',
     -     '   (watts/m*m)  (watts/m*m)     (bar)       (K)'
c
          pbar = 1.e-5*p(1)
          write(iuheat,'(3(1pe12.4),24x,5(1pe13.5),2(1pe13.4))')
     -              pbar,t(1),alt(1),dirsoflx(1,nz),dnsoflx(1,nz),
     -              upsoflx(1,nz),dnthflx(1),upthflx(1),
     -              1.e-5*p(1),t(1)
          do 6201 k=1,nlyr 
              pbar = 5.e-6*(p(k) + p(k+1))
              tbar = 0.5*(t(k) + t(k+1))
              altkm = 0.5*(alt(k) + alt(k+1))
              soheat_d = 86400.*real(soheat(k,nz))
              thheat_d = 86400.*real(thheat(k))
              write(iuheat,'(5(1pe12.4),5(1pe13.5),2(1pe13.4))')
     -              pbar,tbar,altkm,soheat_d,thheat_d,
     -              dirsoflx(k+1,nz),dnsoflx(k+1,nz),upsoflx(k+1,nz),
     -              dnthflx(k+1),upthflx(k+1),1.e-5*p(k+1),t(k+1)
6201      continue
6221  continue
c
      return
      end
