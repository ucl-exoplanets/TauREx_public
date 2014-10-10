      subroutine interp_rad(lsolar,lplanck,
     -           ng0,l_1,l_2,nlev,nz0,nphi,numu,nzup,nzdn,
     -           wn_io,umu0nz,umu,alb_b,dtau_b,copi0_b,g_b,
     -           trnflx_b,dtrnflxdx,refflx_b,drefflxdx,
     -           absflx_b,dabsflxdx,refrad_b,drefraddx,
     -           absrad_b,dabsraddx,brdf_b,dbrdfdx,
     -           tsurf,tatm,alb0,dtau,copi0,g0,dalb,deltau,delpi0,delg,
     -           dnsflxsrc_b,ddnsflxdx,upsflxsrc_b,dupsflxdx,
     -           dntflxsrc_b,ddntflxdx,uptflxsrc_b,duptflxdx,
     -           sradsrc_b,dsraddx,tradsrc_b,dtraddx,
     -           bb,bb_sur,bb_flx_up,bb_flx_dn,bb_rad,
     -           dn_s_src,up_s_src,rad_s_src,
     -           dn_t_src,up_t_src,rad_t_src,
     -           trndir,trnflx,refflx,absflx,trnrad,refrad,absrad,
     -           upsflx,dnsflx,dirsflx,sol_rad,uptflx,dntflx,th_rad)
c
ccccccccccccccccccccccc   i n t e r p _ r a d   cccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine linearly interpolates the binned radiance and   cc
cc    flux source functions for each layer to the optical depth,      cc
cc    single scatterng albedo, asymmetry parameter, surface pressure, cc
cc    and surface albedo at the wavelength of interest.               cc
cc                                                                    cc
cc    notes:                                                          cc
cc           This version uses a full surface bdrf, but only the      cc
cc           first element in the surface optical property vector     cc
cc           is used in the interpolation.  The assumption is that    cc
cc           variation in the other elements are correlated with      cc
cc           first element.                                           cc
cc           This version uses optical depth interpolation.           cc
cc                                                                    cc
cc    NOTE: this version has been modified to produce layer dependent cc
cc    flux source functions and their partial derivatives, as well    cc
cc    as the corresponding flux values.                               cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc         l1 - starting vertical index for source calculation        cc
cc         l2 - ending vertical index for source calculation          cc
cc        ng0 - spectral mapping bin index for this spectral point    cc
cc       nlev - number of model vertical levels                       cc
cc        nz0 - index of this solar zenith angle (1 to nsol)          cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc      wn_io - current wavenumber (cm**-1)                           cc
cc       nzup - number of upward steams for radiances                 cc
cc       nzdn - number of downward streams for radiances              cc
cc       numu - total number of streams for radiances                 cc
cc       nphi - total number of azimuths for radiances                cc
cc dnsflxsrc_b - downward flux source at each level for this bin      cc
cc upsflxsrc_b - downward flux source at each level for this bin      cc
cc   sradsrc_b - angle-dependent radiance source at each level for    cc
cc              this bin                                              cc
cc  ddnsflxdx - downward flux jacobian at each level for this bin     cc
cc  dupsflxdx - upward flux jacobian at each level for this bin       cc
cc    dsraddx - radiance jacobian at each level for this bin          cc
cc dntflxsrc_b - downward flux source at each level for this bin      cc
cc uptflxsrc_b - downward flux source at each level for this bin      cc
cc   tradsrc_b - angle-dependent radiance source at each level for    cc
cc              this bin                                              cc
cc  ddntflxdx - downward flux jacobian at each level for this bin     cc
cc  duptflxdx - upward flux jacobian at each level for this bin       cc
cc    dtraddx - radiance jacobian at each level for this bin          cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     uptflx - wn-dependent upward thermal flux                      cc
cc     dntflx - wn-dependent downward thermal flux                    cc
cc     th_rad - wn-depndent, angle-dependent thermal radiances        cc
cc     upsflx - wn-dependent upward solar flux                        cc
cc     dnsflx - wn-dependent downward direct+diffuse solar flux       cc
cc    dirsflx - wn-dependent downward direct solar flux               cc
cc    sol_rad - wn-depndent, angle-dependent solar radiances          cc
cc                                                                    cc
ccccccccccccccccccccccc   i n t e r p _ r a d   cccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical lplanck,lsolar
c
      integer l_1,l_2,nlev,nz0,nphi,numu,nzup,nzdn
c
      integer k,naz,nze,ibn,lev1,l2tau
c
c****   binned optical properties
c
      integer ng0,nlyr
c
c*****   state vector variables.
c
      real alb_b(nsol,5,ngrp),dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     g_b(kp,5,ngrp)
c
c****   binned layer radiance transmittances and absorptances
c
      double precision refrad_b(mxumu,kp,5,ngrp),
     -                 drefraddx(mxumu,kp,3,ngrp),
     -                 absrad_b(mxumu,kp,5,ngrp),
     -                 dabsraddx(mxumu,kp,3,ngrp),
     -                 brdf_b(mxumu,mxphi,5,ngrp),
     -                 dbrdfdx(mxumu,mxphi,ngrp)
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),dtrnflxdx(kp,3,ngrp),
     -                 refflx_b(kp,5,ngrp),drefflxdx(kp,3,ngrp),
     -                 absflx_b(kp,5,ngrp),dabsflxdx(kp,3,ngrp)
c
c****   layer flux transmittances, absorptances, and fluxes
c
      double precision tl(kp),rl(kp),fup(kp),fdn(kp)
c
c****   surface bi-directional reflection function
c
      double precision brdf(mxumu,mxphi),dir(kp)
c
      real solflx
c
c****   output upward and downward radiances
c
      double precision rad0(mxumu,mxphi,kp)
c
c****   binned solar source function variables and jacobians
c
      double precision upsflxsrc_b(kp,ngrp),dnsflxsrc_b(kp,ngrp),
     -                 sradsrc_b(mxumu,mxphi,kp,ngrp)
c
      double precision ddnsflxdx(kp,4,ngrp),dupsflxdx(kp,4,ngrp),
     -                 dsraddx(mxumu,mxphi,kp,4,ngrp)
c
c****   binned thermal source function variables and jacobians
c
      double precision uptflxsrc_b(kp,ngrp),dntflxsrc_b(kp,ngrp),
     -                 tradsrc_b(mxumu,mxphi,kp,ngrp)
c
      double precision ddntflxdx(kp,4,ngrp),duptflxdx(kp,4,ngrp),
     -                 dtraddx(mxumu,mxphi,kp,4,ngrp)
c
      double precision wn_io
c
      real umu0nz,umu(mxumu)
c
c****   vertical structure variables at this wavelength
c
      real alb0,tsurf,tatm(kp),dtau(kp),copi0(kp),g0(kp)
c
c****    layer transmission values for simplified adding method
c
      double precision trnrad(mxumu,kp),refrad(mxumu,kp),
     -                 absrad(mxumu,kp)
c
      double precision trnflx(kp),refflx(kp),absflx(kp),trndir(kp)
c
c****   interpolated solar and thermal source functions
c
      double precision up_s_src(kp),dn_s_src(kp),
     -                 rad_s_src(mxumu,mxphi,kp)
      double precision up_t_src(kp),dn_t_src(kp),
     -                 rad_t_src(mxumu,mxphi,kp)
c
      double precision rad_src0(mxumu,mxphi,kp)
c
c****    quantitities used on radiance adding
c
      double precision urt(kp),drb(kp),uft(kp),dfb(kp)
      double precision ft(kp),fb(kp),flxu(kp),flxd(kp)
c
      real dalb,deltau(kp),delpi0(kp),delg(kp)
c
c****   planck function thermal source variables
c
      double precision bb_sur,bb(kp),bb_flx_up(kp),bb_flx_dn(kp),
     -                 bb_rad(mxumu,mxphi,kp)
c
c****   spectrally-dependent output flux and radiance
c
      real upsflx(mxrad),dnsflx(mxrad),dirsflx(mxrad),
     -     sol_rad(mxumu,mxphi,mxrad),
     -     uptflx(mxrad),dntflx(mxrad),
     -     th_rad(mxumu,mxphi,mxrad)
      
      double precision pi,dtau0
      real tsurf0(1),wnio
c
c****   define pi 
c
      data pi /3.141592654d0 /
c
c****   specify the bin index
c
      ibn = ng0
      nlyr = nlev - 1
c
c****  set the level range counters
c
      if(l_1 .eq. 1) then
        lev1 = 2
      else
        lev1 = l_1
      endif
      if(l_2 .lt. nlev) then
        l2tau = l_2
      else
        l2tau = nlev-1
      endif
c
c****    find the surface optical properties differences.
c        for surface optical properties, use only the first element 
c        of the surface optical property vector for interpolation.
c
c      dalb = 0.0
      dalb = alb0 - alb_b(nz0,1,ibn)
c
c****     find the optical depth and single scattering albedo 
c         and asymmetry parameter differences.
c
      do 1001 k=l_1,l2tau
          deltau(k) = dtau(k) - dtau_b(k,1,ibn)
          delpi0(k) = copi0(k) - copi0_b(k,1,ibn)
          delg(k) = g0(k) - g_b(k,1,ibn)
c          deltau(k) = 0.0
c          delpi0(k) = 0.0
c          delg(k) = 0.0
1001  continue
c
c****     Layer transmittances are used in mapping fluxes and
c         radiances back to the high resolution spectrum.
c         find the transmission through each layer for each stream.
c
      do 1221 k=l_1,l2tau
c
c****            find the flux transmittances and absorptances
c
          trnflx(k) = trnflx_b(k,1,ibn) +
     -                dtrnflxdx(k,1,ibn)*deltau(k) +
     -                dtrnflxdx(k,2,ibn)*delpi0(k) +
     -                dtrnflxdx(k,3,ibn)*delg(k)
          refflx(k) = refflx_b(k,1,ibn) +
     -                drefflxdx(k,1,ibn)*deltau(k) +
     -                drefflxdx(k,2,ibn)*delpi0(k) +
     -                drefflxdx(k,3,ibn)*delg(k)
c
c****            find the radiance transmittances and absorptances
c
          do 1201 nze=1,numu
              dtau0 = dtau(k)/abs(umu(nze))
              trnrad(nze,k) = dexp(-dtau0)
              refrad(nze,k) = refrad_b(nze,k,1,ibn) +
     -                        drefraddx(nze,k,1,ibn)*deltau(k) +
     -                        drefraddx(nze,k,2,ibn)*delpi0(k) +
     -                        drefraddx(nze,k,3,ibn)*delg(k)
1201      continue
1221  continue
c
      refflx(nlev) = alb0
      trnflx(nlev) = 0.0
      do 1261 naz=1,nphi
          do 1241 nze=1,numu
              brdf(nze,naz) = brdf_b(nze,naz,1,ibn) + 
     -                        dbrdfdx(nze,naz,ibn)*dalb
1241      continue
1261  continue
c
      if(lplanck) then
        do 1421 k=l_1,l2tau
c
c****            find the flux transmittances and absorptances
c
            absflx(k) = absflx_b(k,1,ibn) +
     -                  dabsflxdx(k,1,ibn)*deltau(k) +
     -                  dabsflxdx(k,2,ibn)*delpi0(k) +
     -                  dabsflxdx(k,3,ibn)*delg(k)
            do 1401 nze=1,numu
c
c****            find the radiance transmittances and absorptances
c
                absrad(nze,k) = absrad_b(nze,k,1,ibn) +
     -                          dabsraddx(nze,k,1,ibn)*deltau(k) +
     -                          dabsraddx(nze,k,2,ibn)*delpi0(k) + 
     -                          dabsraddx(nze,k,3,ibn)*delg(k)
1401        continue
1421    continue
c
      endif
c
      if(lsolar) then
c
c****    find the direct fluxes, dir
c
        solflx = 1.0
        dir(1) = umu0nz*solflx
        do 1601 k=1,nlyr
            dir(k+1) = dir(k)*trndir(k)
1601    continue
c
      endif
c
c****    i n t e r p o l a t e   s o u r c e   f u n c t i o n s
c
      call interp_src(lsolar,lplanck,ng0,nlev,l_1,l_2,
     -                      lev1,l2tau,nzup,nzdn,numu,nphi,
     -                      dnsflxsrc_b,upsflxsrc_b,sradsrc_b,
     -                      ddnsflxdx,dupsflxdx,dsraddx,
     -                      dntflxsrc_b,uptflxsrc_b,tradsrc_b,
     -                      ddntflxdx,duptflxdx,dtraddx,
     -                      dalb,deltau,delpi0,delg,alb0,wn_io,
     -                      up_s_src,dn_s_src,rad_s_src,
     -                      up_t_src,dn_t_src,rad_t_src)
c
      if(lsolar) then
c
c**********************************************************************
c
c*****         f i n d   s o l a r    r a d i a n c e s
c
c**********************************************************************
c
c****      initialize solar radiances and fluxes at all levels
c
        do 2041 k=1,nlev
            upsflx(k) = 0.0
            dnsflx(k) = 0.0
            dirsflx(k) = 0.0
c
c****        initialize output radiances
c
            do 2021 naz=1,nphi
                do 2001 nze=1,numu
                    sol_rad(nze,naz,k) = 0.0
2001            continue
2021        continue
2041    continue
c
c****     use the flux adding method to find the 
c         upward and downward fluxes
c
        fdn(1) = dn_s_src(1)
        do 2201 k=1,nlyr
            rl(k) = refflx(k)
            tl(k) = trnflx(k)
            fdn(k+1) = dn_s_src(k+1)
            fup(k) = up_s_src(k)
2201    continue
c
        rl(nlev) = alb0
        tl(nlev) = 0.0d0
c
c*****    define the surface flux source function - for fluxes
c         the reflected solar gives and exact solution (everything
c         else is accouted for in the other terms).
c
        fup(nlev) = alb0*dir(nlev)

c
        call add_sflx(nlyr,tl,rl,umu0nz,dir,fup,fdn,urt,drb,
     -                    uft,dfb,ft,fb,flxu,flxd)
c
c****         set direct and diffuse fluxes at all levels
c
        do 2441 k=1,nlev
            upsflx(k) = real(flxu(k))
            dnsflx(k) = real(flxd(k))
            dirsflx(k) = real(dir(k))
2441    continue
c
c****     find upward and downward radiances at all levels
c
c****     load the radiance source terms
c
        do 2641 naz=1,nphi
            do 2621 nze=1,numu
                do 2601 k=1,nlev
                    rad_src0(nze,naz,k) = rad_s_src(nze,naz,k)
2601            continue
2621        continue
2641    continue
c
        call add_rad(nlyr,nzdn,nzup,numu,nphi,trnrad,refrad,brdf,
     -                    rad_src0,ft,fb,urt,drb,uft,dfb,rad0)
c
        do 2841 naz=1,nphi
            do 2821 nze=1,numu
                do 2801 k=1,nlev
                    sol_rad(nze,naz,k) = real(rad0(nze,naz,k))
2801            continue

2821        continue
2841    continue
c
c****    exit solar block
c
      endif
c
c***********************************************************************
c
c****           f i n d    t h e r m a l    r a d i a n c e s
c
c***********************************************************************
c
      if(lplanck .and. nz0 .eq. 1) then
c
c****     find the surface planck source at this wavenumber
c
        wnio = real(wn_io)
        if(l_2 .eq. nlev) then
          tsurf0(1) = tsurf
          call planck(1,1,1,1,wnio,tsurf0,bb)
          bb_sur = bb(1)
          bb_flx_up(nlev) = pi*bb_sur*(1.0d0 - alb0)
        endif
c
c****    find the atmospheric planck functions at this wavenumber
c
        call planck(1,l_2,1,1,wnio,tatm,bb)
c
c****         define the black body source functions for flux and 
c             radiance, assuming a linear-in-tau formulation.
c
        if(l_1 .eq. 1) bb_flx_dn(1) = pi*bb(1)*copi0(1)*absflx(1)
        do 3001 k=l_1,l2tau
            bb_flx_dn(k+1) = 0.5d0*pi*(bb(k) + bb(k+1))*
     -                       copi0(k)*absflx(k)
            bb_flx_up(k) = 0.5d0*pi*(bb(k) + bb(k+1))*
     -                     copi0(k)*absflx(k)
3001    continue
c
c****    thermal sources for downward radiance streams at top of
c        atmosphere and upward radiance streams at surface   
c
        do 3241 naz=1,nphi
            if(l_1 .eq. 1) then
c
c****          the following creates a non-zero scaling factor for 
c              the top layer
c
              do 3201 nze=1,nzdn
                   bb_rad(nze,naz,1) = bb(1)*copi0(1)*absrad(nze,1)
3201          continue
c
            endif
c
            if(l_2 .eq. nlev) then
c
c****           set upward thermal source at surface
c
              do 3221 nze=nzup,numu
                  bb_rad(nze,naz,nlev) = bb_sur*
     -                                   (1.0d0 - brdf(nze,naz))
3221          continue
c
            endif
3241    continue            
c
c****             thermal source for downward streams
c
        do 3461 k=l_1,l2tau
            naz = 1
            do 3401 nze=1,nzdn
                 bb_rad(nze,naz,k+1) = 0.5d0*(bb(k) + bb(k+1))*
     -                                 copi0(k)*absrad(nze,k)
3401        continue
c
c****         set values at all other azimuth angles 
c             (no azimuth dependence above surface)
c
            do 3441 naz=2,nphi
                do 3421 nze=1,nzdn
                    bb_rad(nze,naz,k+1) = bb_rad(nze,1,k+1)
3421            continue
3441        continue
3461    continue
c
c****       thermal source for upward streams
c
        do 3661 k=l_1,l2tau
            naz = 1
            do 3601 nze=nzup,numu
                bb_rad(nze,naz,k) = 0.5d0*(bb(k) + bb(k+1))*
     -                              copi0(k)*absrad(nze,k)
3601        continue
            do 3641 naz=2,nphi
                do 3621 nze=nzup,numu
                    bb_rad(nze,naz,k) = bb_rad(nze,1,k)
3621            continue
3641        continue
3661    continue
c
c****     initialize thermal fluxes and radiances at all levels
c
        do 4241 k=1,nlev
            uptflx(k) = 0.0
            dntflx(k) = 0.0
c
c****     initialize thermal radiances at all levels
c
            do 4221 naz=1,nphi
                do 4201 nze=1,numu
                    th_rad(nze,naz,k) = 0.0
4201            continue
4221        continue
4241    continue
c
c****  use the flux adding method to find the upward and downward fluxes
c
        do 4401 k=1,nlyr
            rl(k) = refflx(k)
            tl(k) = trnflx(k)
            fdn(k) = bb_flx_dn(k)*dn_t_src(k) + bb_flx_dn(k)
            fup(k) = bb_flx_up(k)*up_t_src(k) + bb_flx_up(k)
4401    continue
        k = nlev
        rl(k) = alb0
        tl(k) = 0.0d0
        fdn(k) = bb_flx_dn(k)*dn_t_src(k) + bb_flx_dn(k)
        fup(k) = bb_flx_up(k)*up_t_src(k) + bb_flx_up(k)
c
        call add_tflx(nlyr,tl,rl,fup,fdn,urt,drb,
     -                    uft,dfb,ft,fb,flxu,flxd)
c
c****         set downward direct and diffuse fluxes at all levels
c
        do 4421 k=1,nlev
            uptflx(k) = real(flxu(k))
            dntflx(k) = real(flxd(k))
4421    continue

c      if(nz0 .eq. 1 .and.  (l_1 .eq. 1 .or. l_1 .eq. 10)) then
c       write(*,'(1a,16(1pe12.4))') 'dn_t_src:      ',
c     - (dn_t_src(k),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'up_t_src:      ',
c     - (up_t_src(k),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'dntflx:      ',
c     - (dntflx(k),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'uptflx:      ',
c     - (uptflx(k),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'bb_flx_dn    ',
c     - (bb_flx_dn(k),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'bb_flx_up    ',
c     - (bb_flx_up(k),k=1,nlev)
c      do 7 nze=1,numu
c7      write(*,'(1a,i3,16(1pe12.4))') 'bb_rad ',nze,
c     - (bb_rad(nze,1,k),k=1,nlev)
c      endif
c
c****     rescale the radiance source terms
c
        do 4641 naz=1,nphi
c
c****        scale downward thermal sources at top of atmosphere
c
           k = 1
           do 4601 nze=1,nzdn
c                rad_src0(nze,naz,k) = bb_rad(nze,naz,k)*
c     -                                  rad_t_src(nze,naz,k) +
c     -                                  bb_rad(nze,naz,k)
                rad_src0(nze,naz,k) = rad_t_src(nze,naz,k) +
     -                                  bb_rad(nze,naz,k)
4601        continue
c
c****        scale upward thermal sources at surface
c
            k = nlev
            do 4621 nze=nzup,numu
c                rad_src0(nze,naz,k) = bb_rad(nze,naz,k)*
c     -                                rad_t_src(nze,naz,k) +
c     -                                bb_rad(nze,naz,k)
                rad_src0(nze,naz,k) = rad_t_src(nze,naz,k) +
     -                                bb_rad(nze,naz,k)
4621        continue
4641    continue
c
c*****    scale upward and downward thermal sources at other levels
c
        do 4861 k=1,nlyr
            do 4841 naz=1,nphi
                do 4801 nze=1,nzdn
c                        rad_src0(nze,naz,k+1) = bb_rad(nze,naz,k)*
c     -                                         rad_t_src(nze,naz,k+1) +
c     -                                          bb_rad(nze,naz,k)
                        rad_src0(nze,naz,k+1) = rad_t_src(nze,naz,k+1) +
     -                                          bb_rad(nze,naz,k)
4801            continue
                do 4821 nze=nzup,numu
c                        rad_src0(nze,naz,k) = bb_rad(nze,naz,k)*
c     -                                        rad_t_src(nze,naz,k) +
c     -                                        bb_rad(nze,naz,k)
                        rad_src0(nze,naz,k) = rad_t_src(nze,naz,k)+
     -                                        bb_rad(nze,naz,k)
4821            continue
4841        continue
4861    continue
c
        call add_rad(nlyr,nzdn,nzup,numu,nphi,trnrad,refrad,brdf,
     -                    rad_src0,ft,fb,urt,drb,uft,dfb,rad0)
c
        do 4941 naz=1,nphi
            do 4921 nze=1,numu
                do 4901 k=1,nlev
                    th_rad(nze,naz,k) = real(rad0(nze,naz,k))
4901            continue
4921        continue
4941    continue
c
      endif
c
      return
      end
