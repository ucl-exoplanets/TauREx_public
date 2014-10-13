      subroutine grey_eq_trn(lsolar,lplanck,usrang,lamber,
     -                       nlev,nz0,nstr,numu,nphi,
     -                       nzup,nzdn,iref,wn_io,dtau,copi0,tatm,
     -                       umu,umu_f,gwt_f,phi,umu0nz,phi0nz,
     -                       fbeam,alb0,btemp,ttemp,temis,
     -                       trndir,trnflx,absflx,trn_cmu,abs_cmu,
     -                       trnrad,absrad,
     -                       upsflx,dnsflx,dirsflx,sol_rad,
     -                       uptflx,dntflx,th_rad,
     -                       dn_s_src,up_s_src,dn_t_src,up_t_src)
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
cc       alb0 - surface albedo and other surface optical properties   cc
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
      logical lplanck,lsolar,lamber,usrang
c
c****   input integer variables
c
      integer nlev,nz0,nstr,numu,nzdn,nzup,nphi,iref
c
      integer k,kk,naz,nze,nlyr
      integer ipass1
      save ipass1
c
c****   input real variables
c
      double precision wn_io
c
      real phi(mxphi),umu(mxumu),phi0nz,umu0nz, 
     -     tatm(kp),dtau(kp),copi0(kp),alb0,
     -     btemp0(1),ttemp0(1),fbeam,btemp,ttemp,temis
c
c       define local real variables
c
      real pi,twopi,wnio
      save pi,twopi
c
      real b(kp),dbdt(kp)
      real bs,bt,fa
      double precision bb(kp)
c
c****    layer transmission values for simplified adding method
c
      double precision trnrad(mxumu,kp),absrad(mxumu,kp),
     -                 trn_cmu(mxumu,kp),abs_cmu(mxumu,kp)
      double precision trnflx(kp),absflx(kp),trndir(kp)

      double precision trsol(kp),trans(kp,mxumu),thrad(mxumu,kp)
c
      double precision trans_flx(kp,mxumu),
     -                 uth_flx(mxumu,kp),dth_flx(mxumu,kp)
c
      real rad_flx(kp,mxumu)
c
      real umu_f(mxumu),gwt_f(mxumu)
c
c****   spectrally-dependent output flux and radiance
c
      real upsflx(mxrad),dnsflx(mxrad),dirsflx(mxrad),
     -     sol_rad(mxumu,mxphi,mxrad),
     -     uptflx(mxrad),dntflx(mxrad),
     -     th_rad(mxumu,mxphi,mxrad)
c
c****   interpolated solar and thermal source functions
c
      double precision up_s_src(kp),dn_s_src(kp)
      double precision up_t_src(kp),dn_t_src(kp)
c
      double precision bb_flx_up(kp),bb_flx_dn(kp)
c
c****   define the speed of light (mks)
c
      data ipass1 / 0 /
c
      nlyr = nlev - 1
      if( ipass1 .eq. 0 ) then
c
        ipass1 = 1
        pi = acos(-1.0 )
        twopi = 2.0*pi
c
      endif
c
c****         find the flux transmittance and absorptance by 
c             integrating transmittance over zenith angle using 
c             gaussian quadrature
c
      do 1221 k=1,nlyr
          trnflx(k) = 0.0
          absflx(k) = 0.0
          do 1201 nze=1,nstr/2
              trnflx(k) = trnflx(k) + 
     -                    2.0*gwt_f(nze)*umu_f(nze)*trn_cmu(nze,k)
              absflx(k) = absflx(k) + 
     -                    2.0*gwt_f(nze)*umu_f(nze)*abs_cmu(nze,k)
1201      continue
1221  continue
c
c***********************************************************************
c
c****             s o l a r    r a d i a n c e s
c
c***********************************************************************
c
      if(lsolar) then
c
c****  initialize radiances and fluxes
c
      do 1841 k=1,nlev
          dn_s_src(k) = 0.0
          up_s_src(k) = 0.0
          dnsflx(k) = 0.0
          dirsflx(k) = 0.0
          upsflx(k) = 0.0
c
          do 1821 nze=1,numu
              rad_flx(k,nze) = 0.0
              do 1801 naz=1,nphi
                  sol_rad(nze,naz,k) = 0.0
1801          continue
1821      continue
1841  continue
c
c****      find the transmission along the solar zenith angle
c 
        trsol(1) = 1.0d0
        do 2001 k=2,nlev
            trsol(k) = trsol(k-1)*trndir(k)
2001    continue
c
c****      find the transmission along each stream from the top of 
c          the atmosphere, to the surface, and back up.  
c
        do 2041 nze=nzup,numu
            trans(nlev,nze) = trsol(nlev)
            do 2021 k=nlyr,1,-1
                trans(k,nze) = trans(k+1,nze)*trnrad(nze,k)
2021        continue
2041    continue
c
c****      find the transmission along each upward stream from 
c          surface to the top of the atmosphere  
c
        do 2081 nze=nstr/2+1,nstr
            trans_flx(nlev,nze) = trsol(nlev)
            do 2061 k=nlyr,1,-1
                trans_flx(k,nze) = trans_flx(k+1,nze)*trn_cmu(nze,k)
2061        continue
2081    continue
c
c****     find the direct downward solar irradiances at each level
c         Note: the downward diffuse radiance should be zero 
c         at each zenith angle and azimuth for the non-scattering case
c
        do 2101 k=1,nlev
            dirsflx(k) = umu0nz*fbeam*real(trsol(k))
2101    continue
c
c****     find the upward scattered solar radiances at each level.
c
c****        u p w a r d   r a d i a n c e s - l a m b e r t
c
        fa = umu0nz*fbeam*alb0/pi
        do 2241 nze=nzup,numu
            do 2221 k=1,nlev
                sol_rad(nze,1,k) = fa*real(trans(k,nze))
                do 2201 naz=2,nphi
                    sol_rad(nze,naz,k) = sol_rad(nze,1,k)
2201            continue
2221        continue
2241    continue
c
c****        u p w a r d   s o l a r    f l u x - l a m b e r t
c
c****      compute the azimuthally-averaged radiances
c          along the models computational streams
c
        do 2281 nze=nstr/2+1,nstr
            do 2261 k=1,nlev
                rad_flx(k,nze) = fa*real(trans_flx(k,nze))
2261        continue
2281    continue
c
c****      use gaussian quadrature to compute the upward fluxes
c
        do 2321 k=1,nlev
            do 2301 nze=nstr/2+1,nstr
                upsflx(k) = upsflx(k) + twopi*umu_f(nze)*gwt_f(nze)*
     -                      rad_flx(k,nze)
2301        continue
2321    continue
c
c****      initialize the effective 'source' for each layer
c
c         note: the downward diffuse solar fluxes are zero by definition
c
c****     initialize the upward fluxes at the surface
c
        up_s_src(nlev) = upsflx(nlev)
c
        do 2401 k=2,nlev
c
c****         find upward flux source at top of each layer
c
            kk = nlev - k + 1
            up_s_src(kk) = upsflx(kk) - upsflx(kk+1)*trnflx(kk)
c
2401    continue
c
      endif
c***********************************************************************
c
c****    t h e r m a l   r a d i a n c e s   a n d    f l u x e s
c
c***********************************************************************
      if(lplanck .and. nz0 .eq. 1) then
c
c****  initialize radiances and fluxes
c
        do 3041 k=1,nlev
            dn_t_src(k) = 0.0
            up_t_src(k) = 0.0
            dntflx(k) = 0.0
            uptflx(k) = 0.0
c
            do 3021 nze=1,numu
                do 3001 naz=1,nphi
                    th_rad(nze,naz,k) = 0.0
3001            continue
3021        continue
3041    continue
c
c****   evaluate planck function at surface and at each model level
c
c****           s u r f a c e    t h e r m a l    e m i s s i o n
c
        wnio = real(wn_io)
        btemp0(1) = btemp
c
        call planck(1,1,1,1,wnio,btemp0,bb)
c
        bs = real(bb(1))
c
c****           t o p    o f   a t m o s p h e r e    e m i s s i o n
c
        if(ttemp .gt. 0.0 .and. temis .gt. 0.0) then
          ttemp0(1) = ttemp
c
          call planck(1,1,1,1,wnio,ttemp0,bb)
c
          bt = real(bb(1))
        else
          bt = 0.0
        endif
c
c****        a t m o s p h e r i c     t h e r m a l    e m i s s i o n
c
        call planck(1,nlev,1,1,wnio,tatm,bb)
c
        do 3101 k=1,nlev
            b(k) = real(bb(k))
3101    continue
c
        do 3121 k=1,nlyr
c
            if(copi0(k)*dtau(k) .ne. 0.0) then
              dbdt(k) = (b(k+1) - b(k))/(copi0(k)*dtau(k))
            else
              dbdt(k) = 0.0
            endif
c
3121    continue
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
            thrad(nze,1) = temis*bt
c
c****             specify downward radiances at all other levels
c
            do 3201 k=1,nlyr
                thrad(nze,k+1) = real(b(k+1) - b(k) + (b(k) + 
     -                           umu(nze)*dbdt(k))*absrad(nze,k) + 
     -                           trnrad(nze,k)*thrad(nze,k))
c Danie Liang
                if(dtau(k).eq.0) thrad(nze,k+1) = thrad(nze,k)
3201        continue
c
3221    continue
c
c****       find downward radiances on the gaussian computational grid
c
        do 3321 nze=1,nstr/2
c
c****         set radiance at top of atmosphere to top BC
c
            dth_flx(nze,1) = temis*bt
c
c****             specify downward radiances at all other levels
c
            do 3301 k=1,nlyr
                dth_flx(nze,k+1) =b(k+1) - b(k) + (b(k) + 
     -                             umu_f(nze)*dbdt(k))*
     -                             abs_cmu(nze,k) + 
     -                             trn_cmu(nze,k)*dth_flx(nze,k)
c Danie Liang
                if(dtau(k).eq.0) dth_flx(nze,k+1) = dth_flx(nze,k)
                if(dth_flx(nze,k+1).lt.0) dth_flx(nze,k+1) = 0.
3301        continue
c
3321    continue
c
c****      use gaussian quadrature to integrate over the upper 
c          hemisphere (note, this explicitly assumes that the
c          downward thermal radiances are azimuthally uniform
c
        do 3351 nze=1,nstr/2
            do 3341 k=1,nlev
c Danie Liang, change size of -twopi to +twopi
                dntflx(k) = dntflx(k) + twopi*umu_f(nze)*
     -                        gwt_f(nze)*real(dth_flx(nze,k))
3341        continue
3351    continue
c
c****     u p w a r d     t h e r m a l    r a d i a n c e
c
c****     find the upward diffuse radiances at each level and angle
c         by integrating upward from the surface.
c
c****        use a lambertian surface
c
        do 3421 nze=nzup,numu
c
c****           specify the surface boundary condition
c               note: This version still assumes a lamber surface 
c               at thermal wavelengths
c
c****             assume that the the surface is lambertian 
c
            thrad(nze,nlev) = (1.0 - alb0)*bs + 
     -                            alb0*dntflx(nlev)/pi
c
c****           find radiance radiances at all other levels
c
            do 3401 k=nlyr,1,-1
c Danie Liang
                thrad(nze,k) = b(k)*absrad(nze,k) - 
     -                           (b(k+1) - b(k))*trnrad(nze,k) + 
     -                           umu(nze)*dbdt(k)*absrad(nze,k) + 
     -                           trnrad(nze,k)*thrad(nze,k+1)
c                thrad(nze,k) = b(k+1) - b(k) + (b(k) - 
c     -                           umu(nze)*dbdt(k))*absrad(nze,k) + 
c     -                           trnrad(nze,k)*thrad(nze,k+1)
3401        continue
c
3421    continue
c
c*****    load output arrays
c
        do 3641 k=1,nlev
            do 3621 naz=1,nphi
                do 3601 nze=1,numu
                    th_rad(nze,naz,k) = real(thrad(nze,k))
3601            continue
3621        continue 
3641    continue
c 
c****     u p w a r d     t h e r m a l    f l u x
c
c****       find downward radiances on gaussian computational grid
c
        do 3721 nze=nstr/2+1,nstr
c
c****         set radiance at top of atmosphere to top BC
c
            uth_flx(nze,nlev) = (1.0 - alb0)*bs + 
     -                             alb0*dntflx(nlev)/pi
c
c****             specify downward radiances at all other levels
c
            do 3701 k=nlyr,1,-1
c                uth_flx(nze,k) = b(k+1) - b(k) + (b(k) - 
c     -                           umu_f(nze)*dbdt(k))*abs_cmu(nze,k) + 
c     -                           trn_cmu(nze,k)*uth_flx(nze,k+1)
c Danie Liang
                uth_flx(nze,k) = b(k)*absrad(nze,k) -
     -                           (b(k+1) - b(k))*trn_cmu(nze,k) +
     -                           umu_f(nze)*dbdt(k)*abs_cmu(nze,k) +
     -                           trn_cmu(nze,k)*uth_flx(nze,k+1)
3701        continue
c
3721    continue
c
c****      use gaussian quadrature to integrate over the upper 
c          hemisphere (note, this explicitly assumes that the
c          downward thermal radiances are azimuthally uniform
c
        do 3751 nze=nstr/2+1,nstr
            do 3741 k=1,nlev
                uptflx(k) = uptflx(k) + twopi*umu_f(nze)*
     -                        gwt_f(nze)*real(uth_flx(nze,k))
3741        continue
3751    continue
c
c****    t h e r m a l   s o u r c e    f u n c t i o n s
c
c****           define the black body source functions for flux and 
c               radiance, assuming a linear-in-tau formulation.
c
        bb_flx_dn(1) = 0.0
        bb_flx_up(nlev) = pi*bs*(1.0d0 - alb0)
        do 3801 k=1,nlyr
            bb_flx_up(k) = 0.5*pi*(b(k) + b(k+1))*absflx(k)
            bb_flx_dn(k+1) = 0.5*pi*(b(k) + b(k+1))*absflx(k)
3801    continue
c
c****     compute the effective 'source function' for flux in each 
c         layer for use in climate modeling
c
        dn_t_src(1) = 0.0
        up_t_src(nlev) = uptflx(nlev) - bb_flx_up(nlev)
c
        do 3821 k=2,nlev
c
c****         find downward source at bottom of each layer
c
            dn_t_src(k) = (dntflx(k) - dntflx(k-1)*trnflx(k-1)) -
     -                                  bb_flx_dn(k)
c
c****         find upward source at top of each layer
c
            kk = nlev - k + 1
            up_t_src(kk) = uptflx(kk) - bb_flx_up(kk) -
     -                                 (uptflx(kk+1)*trnflx(kk)) 
c
3821    continue
c
      endif
c
      return
      end
