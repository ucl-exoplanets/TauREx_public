      subroutine solar_src(ng0,l0,nz0,nlyr,nzdn,nzup,numu,nphi,n_rad,
     -                     alb0,umu0,rfldir,rfldn,dx_b_i,dalb_b_i,
     -                     trnrad_b,rad_rt,rad_rb,flx_ft,flx_fb,
     -                     trnflx_b,flx_rt,flx_rb,flx_du,flx_dd,
     -                     flx_dft,flx_uft,flx_dfb,flx_ufb,
     -                     dnflxsrc,upflxsrc,radsrc,
     -                     dnsflxsrc,upsflxsrc,sradsrc,
     -                     dnsflxsrc_b,upsflxsrc_b,sradsrc_b,
     -                     ddnsflxdx,dupsflxdx,dsraddx)
c
cccccccccccccccccccccccccc  s o l a r _ s r c  ccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine scales binned solar radiance and flux sources   cc
cc    and then loads them and and their jacobians into the arrays     cc
cc    used by interp_rad                                              cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc      rfldir - downward direct solar flux at each level for bin     cc
cc dnsflxsrc - downward flux source at each level for this bin        cc
cc upsflxsrc - downward flux source at each level for this bin        cc
cc   sradsrc - angle-dependent radiance source at each level for      cc
cc               this bin                                             cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc dnsflxsrc_b - downward solar flux source at each level for bin     cc
cc upsflxsrc_b - downward solar flux source at each level for bin     cc
cc   sradsrc_b - angle-dependent solar radiance source at each level  cc
cc              for this bin                                          cc
cc  ddnsflxdx - downward solar flux jacobian at each level for  bin   cc
cc  dupsflxdx - upward solar flux jacobian at each level for bin      cc
cc    dsraddx - solar radiance jacobian at each level for bin         cc
cc                                                                    cc
cccccccccccccccccccccccccc  s o l a r _ s r c  ccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlyr,nphi,numu,l0,ng0,nz0,nzdn,nzup,n_rad(5,ngrp)
      integer ibn,l,k,nze,naz,nlev,nz
c
c****   internal binned source function variables
c
      double precision upsflxsrc(kp,5),dnsflxsrc(kp,5),
     -                 sradsrc(mxumu,mxphi,kp,5)
c
c****   binned solar source function variables and jacobians
c
      real rfldir( mxulv ), rfldn( mxulv )
      real alb0, umu0
      real dx_b_i(kp,5,ngrp),dalb_b_i(nsol,ngrp)
c
c****   source_fcn output variables
c
      double precision upflxsrc(kp,5),dnflxsrc(kp,5),
     -                 radsrc(mxumu,mxphi,kp,5)
c
c****   scaled solar source function variables and jacobians
c
      double precision upsflxsrc_b(kp,ngrp),dnsflxsrc_b(kp,ngrp),
     -                 sradsrc_b(mxumu,mxphi,kp,ngrp)
c
      double precision ddnsflxdx(kp,4,ngrp),dupsflxdx(kp,4,ngrp),
     -                 dsraddx(mxumu,mxphi,kp,4,ngrp)
c
c****  quantities for radiance calculation: downward and upward fluxes
c      at layer interfaces
c
      double precision flx_ft(kp,5),flx_fb(kp,5)
      double precision flx_dfb(kp,5),flx_uft(kp,5)
      double precision flx_dft(kp,5),flx_ufb(kp,5)
c
c****   binned layer radiance transmittances and absorptances
c
      double precision trnrad_b(mxumu,kp,5,ngrp)
c
c****   variables for radiance layer adding method
c
      double precision rad_rb(mxumu,mxphi,kp,5),rad_rt(mxumu,mxphi,kp,5)
c
      double precision pi

      double precision trnflx_b(kp,5,ngrp)
      double precision flx_rb(kp,5),flx_rt(kp,5),
     -                 flx_dd(kp,5),flx_du(kp,5)
      double precision test1(kp),test2(kp)
c
      real dflx,uflx
c
c****   define pi 
c
      data pi /3.141592654d0 /
c
      l = l0
      ibn = ng0
      nz = nz0
      nlev = nlyr + 1
c
      if(l .eq. 1 .or. n_rad(l,ibn) .ne. 0) then
c
c****     set values of fluxes at top and bottom of atmosphere.
c
        dnsflxsrc(1,l) = 0.0
        upsflxsrc(nlev,l) = alb0*rfldir(nlev)
c
        do 1201 k=1,nlyr
c
c****         define the flux scaling factor
c
            dflx = rfldir(k)/umu0 + flx_dft(k,l)
c
            if(dflx .gt. 0.0) then
c
              test1(k+1) = dnflxsrc(k+1,l)/(rfldir(k)/umu0 +
     -            (flx_ft(k,l) + flx_rt(k,l)*upflxsrc(k,l))*flx_dd(k,l))
              uflx = rfldir(k+1) + flx_ufb(k+1,l)
              test2(k+1) = dnflxsrc(k+1,l)/(rfldir(k)/umu0 +
     -       (flx_ft(k,l) + flx_rt(k,l)*upflxsrc(k,l)/uflx)*flx_dd(k,l))

              dnsflxsrc(k+1,l) = dnflxsrc(k+1,l)/dflx
c
            else
c
                dnsflxsrc(k+1,l) = 0.0
c
            endif
c
            uflx = rfldir(k+1) + flx_ufb(k+1,l)
c
            if(uflx .ne. 0.0) then
c
              upsflxsrc(k,l) = upflxsrc(k,l)/uflx
c
            else
c
                upsflxsrc(k,l) = 0.0
c
            endif
c
1201    continue
c
c****      scale radiances 
c             note: the upward radiances at the surface (k=nlev) are
c                   not scaled since trnrad_b(nlev) = 0.0
c
c****     set values of radiances at top and bottom of atmosphere
c         where they are not scaled
c
        do 2041 naz=1,nphi
            do 2001 nze=1,nzdn
                sradsrc(nze,naz,1,l) = radsrc(nze,naz,1,l)
2001        continue
            do 2021 nze=nzup,numu
                sradsrc(nze,naz,nlev,l) = radsrc(nze,naz,nlev,l)
2021        continue
2041    continue
c
c****     scale values at other layers
c
        do 2261 k=1,nlyr
            do 2241 naz=1,nphi
                do 2201 nze=1,nzdn
                    sradsrc(nze,naz,k+1,l) = radsrc(nze,naz,k+1,l) + 
     -                                       trnrad_b(nze,k,l,ibn)*
     -                                       rad_rt(nze,naz,k,l)*
     -                                       flx_uft(k,l)/pi
2201            continue
c
                do 2221 nze=nzup,numu
                    sradsrc(nze,naz,k,l) = radsrc(nze,naz,k,l) + 
     -                                     trnrad_b(nze,k,l,ibn)*
     -                                     rad_rb(nze,naz,k+1,l)*
     -                                     flx_dfb(k+1,l)/pi
2221            continue
2241        continue
2261    continue
c
      endif
c
c****  store flux and radiance source variables in bin arrays
c
      if(l0 .eq. 1) then
c
c****           load nominal flux and radiance source terms
c
        do 2641 k=1,nlev
c                    
            dnsflxsrc_b(k,ibn) = dnsflxsrc(k,1)
            upsflxsrc_b(k,ibn) = upsflxsrc(k,1)
c
            do 2621 naz=1,nphi
                do 2601 nze=1,numu
                     sradsrc_b(nze,naz,k,ibn) = sradsrc(nze,naz,k,1)
2601            continue
2621        continue
2641    continue
c
      else
c
        if(n_rad(l,ibn) .ne. 0) then
c
c****         flux derivatives at top of atmosphere and surface
c
          ddnsflxdx(1,l-1,ibn) = 0.0d0
c
          do 2781 naz=1,nphi
c
c****         downward radiance derivatives at top of atmosphere
c
              do 2721 nze=1,nzdn
                  dsraddx(nze,naz,1,l-1,ibn) = 0.0
2721          continue
c
c****           upward radiance streams derivatives at the surface
c
              if(l .lt. 5) then
                do 2741 nze=nzup,numu
                    dsraddx(nze,naz,nlev,l-1,ibn) = 0.0
2741            continue
c
              else
c
                do 2761 nze=nzup,numu
                    dsraddx(nze,naz,nlev,l-1,ibn) = dalb_b_i(nz,ibn)*
     -                        (sradsrc(nze,naz,nlev,l) - 
     -                         sradsrc(nze,naz,nlev,1))
2761            continue
              endif
2781      continue
c
c***        find flux and radiance derivatives at other levels
c
          do 2861 k=1,nlyr
c
c****          flux derivatives each level
c
              ddnsflxdx(k+1,l-1,ibn) = dx_b_i(k,l,ibn)*
     -                   (dnsflxsrc(k+1,l) - dnsflxsrc(k+1,1))
              dupsflxdx(k,l-1,ibn) = dx_b_i(k,l,ibn)*
     -                   (upsflxsrc(k,l) - upsflxsrc(k,1))
c
              do 2841 naz=1,nphi
c
c****             downward radiance derivatives at each level
c
                  do 2801 nze=1,nzdn
                      dsraddx(nze,naz,k+1,l-1,ibn) = dx_b_i(k,l,ibn)*
     -                             (sradsrc(nze,naz,k+1,l) - 
     -                              sradsrc(nze,naz,k+1,1))
2801              continue
c
c****                 upward radiance streams derivatives at  surface
c
                  do 2821 nze=nzup,numu
                      dsraddx(nze,naz,k,l-1,ibn) = dx_b_i(k,l,ibn)*
     -                             (sradsrc(nze,naz,k,l) - 
     -                              sradsrc(nze,naz,k,1))
2821              continue
2841          continue
2861      continue
c
        else
c
c****       n_rad = 0, not a variable part of the state vector
c
          do 2941 k=1,nlev
              ddnsflxdx(k,l-1,ibn) = 0.0d0
              dupsflxdx(k,l-1,ibn) = 0.0d0
              do 2921 naz=1,nphi
                  do 2901 nze=1,numu
                      dsraddx(nze,naz,k,l-1,ibn) = 0.0d0
2901              continue
2921          continue
2941      continue
c
        endif
c
      endif
c
      return
      end
