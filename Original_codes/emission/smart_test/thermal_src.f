      subroutine thermal_src(ng0,l0,nz0,nlyr,nzdn,nzup,numu,nphi,n_rad,
     -                       wng0,alb0,brdf_b,ts,t,dx_b_i,dalb_b_i,
     -                       trnrad_b,rad_rt,rad_rb,flx_uft,flx_dfb,
     -                       dnflxsrc,upflxsrc,radsrc,
     -                       dntflxsrc,uptflxsrc,tradsrc,
     -                       copi0_b,absflx_b,absrad_b,
     -                       dntflxsrc_b,uptflxsrc_b,tradsrc_b,
     -                       ddntflxdx,duptflxdx,dtraddx,bb_rad)
c
ccccccccccccccccccccccccc  t h e r m a l _ s r c  cccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine scales binned thermal radiance and flux sources cc
cc    and then loads them and and their jacobians into the arrays     cc
cc    used by interp_rad                                              cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc         ng0 - spectral bin index                                   cc
cc   dntflxsrc - downward flux source at each level for this bin      cc
cc   uptflxsrc - downward flux source at each level for this bin      cc
cc     tradsrc - angle-dependent radiance source at each level for    cc
cc               this bin                                             cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc dntflxsrc_b - downward thermal flux source at each level for bin   cc
cc uptflxsrc_b - downward thermal flux source at each level for bin   cc
cc   tradsrc_b - angle-dependent thermal radiance source at each      cc
cc              level for this bin                                    cc
cc  ddntflxdx - downward thermal flux jacobian at each level for bin  cc
cc  duptflxdx - upward thermal flux jacobian at each level for bin    cc
cc    dtraddx - thermal radiance jacobian at each level for bin       cc
cc                                                                    cc
ccccccccccccccccccccccccc  t h e r m a l _ s r c  cccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlyr,nphi,numu,l0,ng0,nz0,nzdn,nzup,n_rad(5,ngrp)
      integer ibn,l,k,nze,naz,nlev,iunits,nz
c
      double precision absflx_b(kp,5,ngrp)
c
c****   source_fcn output variables
c
      double precision upflxsrc(kp,5),dnflxsrc(kp,5),
     -                 radsrc(mxumu,mxphi,kp,5)
c
c****   scaled thermal source function variables
c
      double precision uptflxsrc(kp,5),dntflxsrc(kp,5),
     -                 tradsrc(mxumu,mxphi,kp,5)
c
c****   internal black body sources
c
      double precision bb_sur,bb(kp),bb_flx_up(kp),bb_flx_dn(kp),
     -                 bb_rad(mxumu,mxphi,kp)
c
c****   binned thermal source function variables and jacobians
c
      double precision uptflxsrc_b(kp,ngrp),dntflxsrc_b(kp,ngrp),
     -                 tradsrc_b(mxumu,mxphi,kp,ngrp)
c
      double precision ddntflxdx(kp,4,ngrp),duptflxdx(kp,4,ngrp),
     -                 dtraddx(mxumu,mxphi,kp,4,ngrp)
c
c****  quantities for radiance calculation: downward and upward fluxes
c      at layer interfaces
c
      double precision flx_dfb(kp,5),flx_uft(kp,5)
c
c****   binned layer radiance transmittances and absorptances
c
      double precision trnrad_b(mxumu,kp,5,ngrp),
     -                 absrad_b(mxumu,kp,5,ngrp),
     -                 brdf_b(mxumu,mxphi,5,ngrp)
c
c****   variables for radiance layer adding method
c
      double precision rad_rb(mxumu,mxphi,kp,5),rad_rt(mxumu,mxphi,kp,5)
c
      real dx_b_i(kp,5,ngrp),dalb_b_i(nsol,ngrp),wng0,ts,t(kp),alb0
      real copi0_b(kp,5,ngrp)
c
      real tsurf0(1)
c
      double precision pi
c
      data pi / 3.141592654d0 /
c
      l = l0
      ibn = ng0
      nz = nz0
      nlev = nlyr + 1
      iunits = 1
c
      if(l .eq. 1 .or. n_rad(l,ibn) .ne. 0) then
c
c****   evaluate planck function at the surface
c
        tsurf0(1) = ts
c
        call planck(1,1,1,iunits,wng0,tsurf0,bb)
c
        bb_sur = bb(1)
c
c****    find the atmospheric planck functions at this wavenumber
c
        call planck(1,nlev,1,iunits,wng0,t,bb)
c
c****         define the black body source functions for flux and
c             radiance, assuming a linear-in-tau formulation.
c
        bb_flx_dn(1) = pi*bb(1)*copi0_b(1,l,ibn)*absflx_b(1,l,ibn)
        bb_flx_up(nlev) = pi*bb_sur*(1.0d0 - alb0)
c
        do 1001 k=1,nlyr
          bb_flx_dn(k+1) = 0.5d0*pi*(bb(k) + bb(k+1))*
     -                     copi0_b(k,l,ibn)*absflx_b(k,l,ibn)
          bb_flx_up(k) = 0.5d0*pi*(bb(k) + bb(k+1))*
     -                   copi0_b(k,l,ibn)*absflx_b(k,l,ibn)
1001    continue
c
c****    thermal sources for downward radiance streams at top of
c        atmosphere and upward radiance streams at surface   
c
        do 1241 naz=1,nphi
            do 1201 nze=1,nzdn
c
c****            define a non-zero thermal source for scaling top level
c            
                bb_rad(nze,naz,1) = bb(1)*copi0_b(1,l,ibn)*
     -                                 absrad_b(nze,1,l,ibn)
1201        continue
c
c****         define the surface blackbody emission
c 
            do 1221 nze=nzup,numu
                bb_rad(nze,naz,nlev) = bb_sur*
     -                                 (1.0d0 - brdf_b(nze,naz,l,ibn))
1221        continue
1241    continue            
c
c****             thermal source for downward streams
c
        do 1461 k=1,nlyr
            naz = 1
            do 1401 nze=1,nzdn
                 bb_rad(nze,naz,k+1) = 0.5d0*(bb(k) + bb(k+1))*
     -                                 copi0_b(k,l,ibn)*
     -                                 absrad_b(nze,k,l,ibn)
1401        continue
c
c****         set values at all other azimuth angles 
c             (no azimuth dependence above surface)
c
            do 1441 naz=2,nphi
                do 1421 nze=1,nzdn
                    bb_rad(nze,naz,k+1) = bb_rad(nze,1,k+1)
1421            continue
1441        continue
1461    continue
c
c****       thermal source for upward streams
c
        do 1661 k=1,nlyr
            naz = 1
            do 1601 nze=nzup,numu
                bb_rad(nze,naz,k) = 0.5d0*(bb(k) + bb(k+1))*
     -                              copi0_b(k,l,ibn)*
     -                              absrad_b(nze,k,l,ibn)
1601        continue
            do 1641 naz=2,nphi
                do 1621 nze=nzup,numu
                    bb_rad(nze,naz,k) = bb_rad(nze,1,k)
1621            continue
1641        continue
1661    continue
c
c****     normalize the fluxes and radiances by the black body 
c         emission in the layer.
c
        do 2001 k=1,nlev
c
c****          scale downward fluxes
c
            if(bb_flx_dn(k) .gt. 0.0d0) then
              dntflxsrc(k,l) = (dnflxsrc(k,l) - bb_flx_dn(k))/
     -                             bb_flx_dn(k)
            else
              dntflxsrc(k,l) = 0.0d0
            endif
c              dntflxsrc(k,l) = dnflxsrc(k,l)
c
c****         scale upward fluxes
c
            if(bb_flx_up(k) .gt. 0.0d0) then
              uptflxsrc(k,l) = (upflxsrc(k,l) - bb_flx_up(k))/
     -                           bb_flx_up(k)
            else
              uptflxsrc(k,l) = 0.0d0
            endif
c            uptflxsrc(k,l) = upflxsrc(k,l)
2001    continue
c
c****     scale radiance source terms at top of atmosphere and surface
c
        do 2061 naz=1,nphi
c
c****          scale downward radiance sources at top of atmosphere
c
            k = 1
            do 2021 nze=1,nzdn
                   tradsrc(nze,naz,k,l) = radsrc(nze,naz,k,l) - 
     -                                    bb_rad(nze,naz,k)
2021        continue
c
c****              scale upward radiance sources at the surface
c
            k = nlev
            do 2041 nze=nzup,numu
                  tradsrc(nze,naz,k,l) = radsrc(nze,naz,k,l) - 
     -                                     bb_rad(nze,naz,k)
2041        continue
2061    continue
c
c****     scale values at other levels
c
        do 2161 k=1,nlyr
c
            do 2141 naz=1,nphi
c
c****              scale downward radiance sources
c
                do 2101 nze=1,nzdn
                      tradsrc(nze,naz,k+1,l) = radsrc(nze,naz,k+1,l) + 
     -                                          trnrad_b(nze,k,l,ibn)*
     -                                          rad_rt(nze,naz,k,l)*
     -                                          flx_uft(k,l)/pi - 
     -                                          bb_rad(nze,naz,k)
2101            continue
c
c****              scale upward radiance sources
c
                do 2121 nze=nzup,numu
                      tradsrc(nze,naz,k,l) = radsrc(nze,naz,k,l) + 
     -                                        trnrad_b(nze,k,l,ibn)*
     -                                        rad_rb(nze,naz,k+1,l)*
     -                                        flx_dfb(k+1,l)/pi - 
     -                                        bb_rad(nze,naz,k)
2121            continue
2141          continue
c
2161    continue
c
      endif
c
      if(l0 .eq. 1) then
c
c****           load nominal flux and radiance source terms
c
        do 2241 k=1,nlev
c                    
            dntflxsrc_b(k,ibn) = dntflxsrc(k,1)
            uptflxsrc_b(k,ibn) = uptflxsrc(k,1)
c
            do 2221 naz=1,nphi
                do 2201 nze=1,numu
                     tradsrc_b(nze,naz,k,ibn) = tradsrc(nze,naz,k,1)
2201            continue
2221        continue
2241    continue
c
      else
c
        if(n_rad(l,ibn) .ne. 0) then
c
c****         flux derivatives at top of atmosphere and surface
c
          ddntflxdx(1,l-1,ibn) = 0.0d0
          duptflxdx(nlev,l-1,ibn) = dx_b_i(nlev,l,ibn)*
     -          (uptflxsrc(nlev,l) - uptflxsrc(nlev,1))
c
          do 2461 naz=1,nphi
c
c****         downward radiance derivatives at top of atmosphere
c
              do 2401 nze=1,nzdn
                  dtraddx(nze,naz,1,l-1,ibn) = 0.0
2401          continue
c
c****           upward radiance streams derivatives at the surface
c
              if(l .lt. 5) then
                do 2421 nze=nzup,numu
                    dtraddx(nze,naz,nlev,l-1,ibn) = 0.0
2421            continue
c
              else
c
                do 2441 nze=nzup,numu
                    dtraddx(nze,naz,nlev,l-1,ibn) = dalb_b_i(nz,ibn)*
     -                      (tradsrc(nze,naz,nlev,l) - 
     -                       tradsrc(nze,naz,nlev,1))
2441            continue
c
              endif
c
2461      continue
c
c***        find flux and radiance derivatives at other levels
c
          do 2661 k=1,nlyr
c
c****          flux derivatives each level
c
              ddntflxdx(k+1,l-1,ibn) = dx_b_i(k,l,ibn)*
     -                   (dntflxsrc(k+1,l) - dntflxsrc(k+1,1))
              duptflxdx(k,l-1,ibn) = dx_b_i(k,l,ibn)*
     -                   (uptflxsrc(k,l) - uptflxsrc(k,1))
c
              do 2641 naz=1,nphi
c
c****             downward radiance derivatives at each level
c
                  do 2601 nze=1,nzdn
                      dtraddx(nze,naz,k+1,l-1,ibn) = dx_b_i(k,l,ibn)*
     -                             (tradsrc(nze,naz,k+1,l) - 
     -                              tradsrc(nze,naz,k+1,1))
2601              continue
c
c****                 upward radiance streams derivatives at  surface
c
                  do 2621 nze=nzup,numu
                      dtraddx(nze,naz,k,l-1,ibn) = dx_b_i(k,l,ibn)*
     -                             (tradsrc(nze,naz,k,l) - 
     -                              tradsrc(nze,naz,k,1))
2621              continue
2641          continue
2661      continue
c
        else
c
c****       n_rad = 0, not a variable part of the state vector
c
          do 2841 k=1,nlev
              ddntflxdx(k,l-1,ibn) = 0.0d0
              duptflxdx(k,l-1,ibn) = 0.0d0
              do 2821 naz=1,nphi
                  do 2801 nze=1,numu
                      dtraddx(nze,naz,k,l-1,ibn) = 0.0d0
2801              continue
2821          continue
2841      continue
c
        endif
c
      endif
c
      return
      end
