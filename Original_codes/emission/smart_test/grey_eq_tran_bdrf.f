      subroutine grey_eq_trn(nlyr,nstr,numu,nphi,nzup,nzdn,iref,
     -                       wnio,dtau,tatm,umu,phi,umu0nz,phi0nz,
     -                       fbeam,surf_opt,btemp,ttemp,temis,
     -                       lsolar,lplanck,usrang,lamber,
     -                       upsflx,dnsflx,dirsflx,sol_rad,
     -                       uptflx,dntflx,th_rad)
c
ccccccccccccccccccccccc  g r e y _ e q _ t r n  cccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine solves the equation of transfer for radiances   cc
cc    and fluxes in non-scattering atmospheres.  It assumes that the  cc
cc    planck function varies linearly with optical depth throughout   cc
cc    each level.                                                     cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc       nlyr - number of computational model layers                  cc
cc       dtau - differential extintion optical depth in each layer    cc
cc       tatm - temperature of each atmospheric level                 cc
cc         wn - wavelnumber                                           cc
cc       nstr - number of gaussian zenith angles used in D/O code     cc
cc     usrang - output radiances at user angles? (logical: T/F)       cc
cc       numu - number of output zenith angles used in D/O code       cc
cc        umu - emission zenith angle cosines                         cc
cc       nphi - number of output azimuth angles                       cc
cc        phi - emission azimuths read from input header              cc
cc       umu0 - cosine of solar zenith angles                         cc
cc       phi0 - solar azimuth angle cosines read from header          cc
cc      fbeam - intensity of collimated flux at top of atmosphere     cc
cc     lamber - Include a lambertian surface? (Logical: T/F)          cc
cc                 note: this version assumes this                    cc
cc   surf_opt - surface albedo and other surface optical properties   cc
cc      btemp - surface temperature                                   cc
cc      ttemp - temperature of uppermost model layer (space)          cc
cc      temis - emissivity of uppermost model layer (space)           cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     dirsflx - downward direct solar flux at each model level       cc
cc      dnsflx - downward diffuse solar flux at each model level      cc
cc      upsflx - upward (diffuse) solar flux at each model level      cc
cc     sol_rad - solar radiance at each model level                   cc
cc      dntflx - downward diffuse thermal flux at each model level    cc
cc      uptflx - upward (diffuse) thermal flux at each model level    cc
cc      th_rad - thermal radiance at each model level                 cc
cc                                                                    cc
ccccccccccccccccccccccc  g r e y _ e q _ t r n  cccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical usrang,lamber,lplanck,lsolar
c
      logical pass1
c
      integer nlyr, numu, nstr, nzdn, nzup, nphi, iref
c
      integer k,naz,nze
c
      real phi(mxphi),umu(mxumu),phi0nz,umu0nz, 
     -     tatm(kp),dtau(kp),surf_opt(4),
     -     btemp0(1),bs0(1),ttemp0(1),bt0(1),
     -     wnio,fbeam,btemp,ttemp,temis
c
c       define local variables
c
      real b(kp),dbdt(kp)
      real pi,twopi,dphi,fa,bdr,bdr_flx
      real c1nu3,bs,bt
      real bdref,dref
c
      double precision dtauk,trsol(kp), tr(kp,mxmom), 
     -       trans(kp,mxmom),tr1(kp,mxmom), 
     -       uth( mxmom,kp), rad(kp,mxmom )
c
      double precision tr_flx(mxulv,mxmom),tr1_flx(mxulv,mxmom),
     -       trans_flx(kp,mxmom),uu_flx(mxumu,kp,mxphi), 
     -       u0_s(mxumu,kp),rad_flx(kp,mxmom ),uth_flx(mxmom,kp)
c
      real g_phi( mxphi ), gwt_phi( mxphi ),
     -     umu_ze( mxmom ), gwt_ze( mxmom )
c
c****   spectrally-dependent output flux and radiance
c
      real upsflx(mxrad),dnsflx(mxrad),dirsflx(mxrad),
     -     sol_rad(mxumu,mxphi,mxrad),
     -     uptflx(mxrad),dntflx(mxrad),
     -     th_rad(mxumu,mxphi,mxrad)
c     ..
      save pass1, g_phi, gwt_phi, umu_ze, gwt_ze, pi, twopi
      data pass1 / .true. /
c
      if( pass1 ) then
        pass1 = .False.
        pi = 2.*asin( 1.0 )
        twopi = 2.0*pi
c
c****      find the gaussian angles needed to perform the integration
c          over azimuth
c
        call qgausn( mxphi/2, g_phi, gwt_phi )
c
        do 1001 k = 1, mxphi / 2
            g_phi( k + mxphi/2 ) = -g_phi( k )
            gwt_phi( k + mxphi/2 ) = gwt_phi( k )
1001    continue
c
c****   find the gaussian angles needed to perform the flux integration
c
        call qgausn(nstr/2,umu_ze,gwt_ze)
c
c****      restack the values into order used by disort
c
        do 1101 nze=1,nstr/2
           gwt_ze(nze+nstr/2) = twopi*gwt_ze(nze)*umu_ze(nze)
           umu_ze(nze+nstr/2) = umu_ze(nze)
1101    continue
        do 1121 nze=1,nstr/2
            gwt_ze(nze) = gwt_ze(nstr-nze+1)
            umu_ze(nze) = -umu_ze(nstr-nze+1)
1121    continue
c
      endif
c
c****   find the layer transmission for each downward zenith angle
c
      do 1402 nze=1,nzdn
          do 1402 k=1,nlyr
              dtauk = dtau(k)/umu(nze)
              tr(k,nze) = dexp(dtauk)
              tr1(k,nze) = 1.0d0 - tr(k,nze)
1402  continue
c
c****   find the layer transmission for each upward zenith angle
c
      do 1422 nze=nzup,numu
          do 1422 k=1,nlyr
              dtauk = dtau(k)/umu(nze)
              tr(k,nze) = dexp(-dtauk)
              tr1(k,nze) = 1.0d0 - tr(k,nze)
1422  continue 
c
      if(usrang) then
c
c****    calculate the transmittances for use in finding solar fluxes
c 
c****    find the layer transmission for each downward zenith angle
c
        do 1502 nze=1,nstr/2
            do 1502 k=1,nlyr
                dtauk = dtau(k)/umu_ze(nze)
                tr_flx(k,nze) = dexp(dtauk)
                tr1_flx(k,nze) = 1.0d0 - tr_flx(k,nze)
1502    continue
c
c****     find the layer transmission for each upward zenith angle
c
        do 1522 nze=nstr/2+1,nstr
            do 1522 k=1,nlyr
                dtauk = dtau(k)/umu_ze(nze)
                tr_flx(k,nze) = dexp(-dtauk)
                tr1_flx(k,nze) = 1.0d0 - tr_flx(k,nze)
1522    continue
c
      endif
c
c****  initialize radiances and fluxes
c
      do 1641 k=1,nlyr+1
          dnsflx(k) = 0.0
          dirsflx(k) = 0.0
          upsflx(k) = 0.0
          dntflx(k) = 0.0
          uptflx(k) = 0.0
c
          do 1602 nze=1,mxumu
              u0_s(nze,k) = 0.0
              do 1602 naz=1,mxphi
                  sol_rad(nze,naz,k) = 0.0
                  th_rad(nze,naz,k) = 0.0
1602      continue
c
          if(usrang) then
            do 1622 nze=1,nstr
                do 1622 naz=1,mxphi
                    uu_flx(nze,k,naz) = 0.0
1622       continue
c
          endif
c
1641  continue
c
c****             s o l a r    r a d i a n c e s
c
      if(lsolar) then
c
c****      find the transmission along the solar zenith angle
c
        trsol(1) = 1.0
        do 2001 k=2,nlyr+1
            dtauk = dtau(k-1)/umu0nz
            trsol(k) = trsol(k-1)*dexp(-dtauk)
2001    continue
c
c****     find the transmission between the top of the atmosphere and 
c         each model level
c
        do 2022 nze=1,nzdn
              trans(1,nze) = 1.0d0
              do 2022 k=2,nlyr+1
                  trans(k,nze) = trans(k-1,nze)*tr(k-1,nze) 
2022    continue
c
c****      find the transmission along each stream from the top of 
c          the atmosphere, to the surface, and back up.  
c
        do 2042 nze=nzup,numu
            trans(nlyr+1,nze) = trsol(nlyr+1)
            do 2042 k=nlyr,1,-1
                trans(k,nze) = trans(k+1,nze)*tr(k,nze)
2042    continue
c
        if(usrang) then
c
c****       find the transmission between the top of the atmosphere and 
c           each model level
c
          do 2122 nze=1,nzdn
              trans_flx(1,nze) = 1.0d0
              do 2122 k=2,nlyr+1
                  trans_flx(k,nze) = trans_flx(k-1,nze)*tr_flx(k-1,nze) 
2122      continue
c
c****      find the transmission along each stream from the top of 
c          the atmosphere, to the surface, and back up.  
c
          do 2142 nze=nzup,nstr
              trans_flx(nlyr+1,nze) = trsol(nlyr+1)
              do 2142 k=nlyr,1,-1
                  trans_flx(k,nze) = trans_flx(k+1,nze)*tr_flx(k,nze)
2142      continue
        else
c
c****        use the transmission values along the usual gaussian angles
c
          do 2202 nze=nzup,nstr
              do 2202 k=1,nlyr+1
                  trans_flx(k,nze) = trans(k,nze)
2202      continue
        endif
c
c****     find the direct solar irradiances at each level
c         Note: the downward diffuse radiance should be zero 
c         at each zenith angle and azimuth for the non-scattering case
c
        do 2501 k=1,nlyr+1
            dirsflx(k) = umu0nz*fbeam*trsol(k)
2501    continue
c
c****     find the upward scattered solar radiances at each level.
c
        if(lamber) then
c
c****        u p w a r d   r a d i a n c e s - l a m b e r t
c
          fa = umu0nz*fbeam*surf_opt(1)/pi
          do 2641 nze=nzup,numu
              do 2621 k=1,nlyr+1
                  sol_rad(nze,1,k) = fa*trans(k,nze)
                  do 2601 naz=2,nphi
                      sol_rad(nze,naz,k) = sol_rad(nze,1,k)
2601              continue
2621          continue
2641      continue
c
c****        u p w a r d   s o l a r    f l u x - l a m b e r t
c
c****      compute the azimuthally-averaged radiances
c          along the models computational streams
c
c****      This version assumes the surface is lambertian
c
          do 2681 nze=nzup,nstr
              do 2661 k=1,nlyr+1
                  u0_s(nze,k) = umu0nz*fbeam*surf_opt(1)*
     -                          trans_flx(k,nze)/pi
2661          continue
2681      continue
c
        else
c
c****       u p w a r d    r a d i a n c e s - b r d f
c
          do 2741 nze=nzup,numu
              do 2721 naz=1,nphi
                  dphi = pi*abs(phi(naz) - phi0nz)/180.
                  bdr = umu0nz*fbeam*
     -                  bdref(umu(nze),umu0nz,dphi,surf_opt,iref)/pi
                  do 2701 k=1,nlyr+1
                      sol_rad(nze,naz,k) = bdr*trans(k,nze)
2701              continue
2721          continue
2741      continue
c
c****        u p w a r d   s o l a r    f l u x  - b r d f
c
c****      compute the azimuthally-averaged radiances
c          along the models computational streams
c
c****      This version assumes the surface is lambertian
c
          bdr_flx = umu0nz*fbeam*dref(umu0nz,surf_opt,iref)/pi
c
          do 2781 nze=nzup,nstr
              do 2761 k=1,nlyr+1
                  u0_s(nze,k) = bdr*trans_flx(k,nze)
2761          continue
2781      continue
c
        endif
c
c****      use gaussian quadrature to compute the upward fluxes
c
        do 2902 nze=nstr/2+1,nstr
            do 2902 k=1,nlyr+1
                upsflx(k) = upsflx(k) + gwt_ze(nze)*u0_s(nze,k)
2902    continue
c
      endif
c
c****    t h e r m a l   r a d i a n c e s   a n d    f l u x e s
c
      if(lplanck) then
c
c****   evaluate planck function at surface and at each model level
c
c****           s u r f a c e    t h e r m a l    e m i s s i o n
c
	btemp0(1) = btemp
c
        call planck(1,1,1,1,wnio,btemp0,bs0,c1nu3)
c
        bs = bs0(1)
c
c****           t o p    o f   a t m o s p h e r e    e m i s s i o n
c
        if(ttemp .gt. 0.0 .and. temis .gt. 0.0) then
          ttemp0(1) = ttemp
c
          call planck(1,1,1,1,wnio,ttemp0,bt0,c1nu3)
c
          bt = bt0(1)
        else
          bt = 0.0
        endif
c
c****        a t m o s p h e r i c     t h e r m a l    e m i s s i o n
c
        call planck(1,1,1,1,wnio,tatm(1),b(1),c1nu3)
c
        do 3021 k=1,nlyr
c
            call planck(1,1,1,1,wnio,tatm(k),b(k+1),c1nu3)
c
            if(dtau(k) .ne. 0.0) then
              dbdt(k) = (b(k+1) - b(k))/dtau(k)
            else
              dbdt(k) = 0.0
            endif
c
3021    continue
c
c****     evaluate angle-dependent radiance contribution by each layer
c
c
c****     downward radiance at base of this layer
c
        do 3102 nze=1,nzdn
            do 3102 k=1,nlyr
                rad(k,nze) = b(k)*tr1(k,nze) + b(k+1) - b(k) +
     -                       umu(nze)*dbdt(k)*tr1(k,nze)
3102    continue
c
c****      upward radiance at top of this layer
c
        do 3122 nze=nzdn+1,nzup
            do 3122 k=1,nlyr
                rad(k,nze) = b(k)*tr1(k,nze) - 
     -                       (b(k+1) - b(k))*tr(k,nze) + 
     -                        umu(nze)*dbdt(k)*tr1(k,nze)
3122    continue
c
        if(usrang) then
c
c****     find downward and upward radiances on the gaussian
c         computational grid for use in the flx calculations
c
c****     downward radiance at base of this layer
c
          do 3142 nze=1,nstr/2
              do 3142 k=1,nlyr
                  rad_flx(k,nze) = b(k)*tr1_flx(k,nze) + 
     -                             b(k+1) - b(k) +
     -                             umu_ze(nze)*dbdt(k)*tr1_flx(k,nze)
3142      continue
c
c****        upward radiance at top of this layer
c
          do 3162 nze=nstr/2+1,nstr
              do 3162 k=1,nlyr
                  rad_flx(k,nze) = b(k)*tr1_flx(k,nze) - 
     -                   (b(k+1) - b(k))*tr_flx(k,nze) + 
     -                    umu_ze(nze)*dbdt(k)*
     -                    tr1_flx(k,nze)
3162      continue
        else
          do 3182 nze=1,nstr
              do 3182 k=1,nlyr
                  rad_flx(k,nze) = rad(k,nze)
3182      continue       
        endif
c
c****     d o w n w a r d     t h e r m a l     r a d i a n c e
c
c****     find the downward diffuse radiances at each level and angle 
c         by integrating downward from the top of the atmosphere.
c
        do 3221 nze=1,nzdn
c
c****         set radiance at top of atmosphere to top BC
c
            uth(nze,1) = temis*bt
c
c****             specify downward radiances at all other levels
c
            do 3201 k=1,nlyr
                uth(nze,k+1) = rad(k,nze) + 
     -                                  tr(k,nze)*uth(nze,k)
3201        continue
c
3221    continue
c
c****      find the downward irradiance (flux) at each level
c 
        if(usrang) then
c
c****       find downward radiances on the gaussian computational grid
c
          do 3321 nze=1,nstr/2
c
c****           set radiance at top of atmosphere to top BC
c
              uth_flx(nze,1) = temis*bt
c
c****               specify downward radiances at all other levels
c
              do 3301 k=1,nlyr
                      uth_flx(nze,k+1) = rad_flx(k,nze) + 
     -                                   tr_flx(k,nze)*uth_flx(nze,k)
3301          continue
c
3321      continue
        else
c
          do 3342 nze=1,nstr/2
              do 3342 k=1,nlyr
                  uth_flx(nze,k+1) = uth(nze,k+1)
3342      continue
c
        endif
c
c****      use gaussian quadrature to integrate over the upper 
c          hemisphere (note, this explicitly assumes that the
c          downward thermal radiances are azimuthally uniform
c
        do 3422 nze=1,nstr/2
            do 3422 k=1,nlyr+1
                dntflx(k) = dntflx(k) + gwt_ze(nze)*uth_flx(nze,k)
3422    continue
c
c****     u p w a r d     t h e r m a l    r a d i a n c e
c
c****     find the upward diffuse radiances at each level and angle
c         by integrating upward from the surface.
c
        do 3621 nze=nzup,numu
c
c****         specify the surface boundary condition
c             note: This version still assumes a lamber surface 
c             at thermal wavelengths
c
c****           assume that the the surface is lambertian 
c
            uth(nze,nlyr+1) = (1.0 - surf_opt(1))*bs + 
     -                          surf_opt(1)*dntflx(nlyr+1)/pi
c
c****         find radiance radiances at all other levels
c
            do 3601 k=nlyr,1,-1
                uth(nze,k) = rad(k,nze) + tr(k,nze)*uth(nze,k+1)
3601        continue
c
3621    continue
c
c****      find the upward irradiance (flux) at each level
c
c****      use gaussian quadrature to integrate over the lower hemisphere
c
        do 3822 nze=nzup,nstr
            do 3822 k=1,nlyr+1
                uptflx(k) = uptflx(k) + gwt_ze(nze)*uth_flx(nze,k)
3822    continue
c
      endif
c
      return
      end
