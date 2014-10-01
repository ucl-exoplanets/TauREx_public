      subroutine interp_src(lsolar,lplanck,ng0,nlev,l_1,l_2,
     -                      lev1,l2tau,nzup,nzdn,numu,nphi,
     -                      dnsflxsrc_b,upsflxsrc_b,sradsrc_b,
     -                      ddnsflxdx,dupsflxdx,dsraddx,
     -                      dntflxsrc_b,uptflxsrc_b,tradsrc_b,
     -                      ddntflxdx,duptflxdx,dtraddx,
     -                      dalb,deltau,delpi0,delg,alb0,wn_io,
     -                      up_s_src,dn_s_src,rad_s_src,
     -                      up_t_src,dn_t_src,rad_t_src)
c
ccccccccccccccccccccccc   i n t e r p _ s r c   cccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine uses linear interpolation to scale the solar    cc
cc    and thermal radiative sources to this spetral point.            cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc        ng0 - spectral mapping bin index for this spectral point    cc
cc       nlev - number of model vertical levels                       cc
cc         l1 - starting vertical index for source calculation        cc
cc         l2 - ending vertical index for source calculation          cc
cc       lev1 - first vertical level for source interpolation         cc
cc      l2tau - last vertical level for source function interpolation cc
cc       nzup - number of upward steams for radiances                 cc
cc       nzdn - number of downward streams for radiances              cc
cc       numu - total number of streams for radiances                 cc
cc       nphi - total number of azimuths for radiances                cc
cc dnsflxsrc_b - downward solar flux source at each level for bin     cc
cc upsflxsrc_b - downward solar flux source at each level for bin     cc
cc   sradsrc_b - angle-dependent solar radiance source at each level  cc
cc              for this bin                                          cc
cc  ddnsflxdx - downward solar flux jacobian at each level for  bin   cc
cc  dupsflxdx - upward solar flux jacobian at each level for bin      cc
cc    dsraddx - solar radiance jacobian at each level for bin         cc
cc dntflxsrc_b - downward thermal flux source at each level for bin   cc
cc uptflxsrc_b - downward thermal flux source at each level for bin   cc
cc   tradsrc_b - angle-dependent thermal radiance source at each      cc
cc              level for this bin                                    cc
cc  ddntflxdx - downward thermal flux jacobian at each level for bin  cc
cc  duptflxdx - upward thermal flux jacobian at each level for bin    cc
cc    dtraddx - thermal radiance jacobian at each level for bin       cc
cc       wnio - current wavenumber (cm**-1)                           cc
cc       dalb - difference between albedo at wnio and binned albedo   cc
cc     deltau - difference between optical depth in each layer and    cc
cc              binned value                                          cc
cc     delpi0 - difference between single scattering albedo in each   cc
cc              layer and binned value                                cc
cc       delg - difference between scattering asymmetry parameter in  cc
cc              each layer and binned value                           cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc   up_s_src - interpolated/scaled upward solar flux source for      cc
cc              each level at this spectral point                     cc
cc   dn_s_src - interpolated/scaled downward solar flux source for    cc
cc              each level at this spectral point                     cc
cc  rad_s_src - interpolated/scaled solar radinace source for each    cc
cc              level at this spectral point                          cc
cc   up_s_src - interpolated/scaled upward thermal flux source for    cc
cc              each level at this spectral point                     cc
cc   dn_s_src - interpolated/scaled downward thermal flux source for  cc
cc              each level at this spectral point                     cc
cc  rad_s_src - interpolated/scaled thermal radinace source for each  cc
cc              level at this spectral point                          cc
cc                                                                    cc
ccccccccccccccccccccccc   i n t e r p _ s r c   cccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical lplanck,lsolar
c
      integer ng0,nlev,l_1,l_2,lev1,l2tau,nzup,nzdn,numu,nphi
      integer nze,naz,k,ibn
c      integer kk
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
c****   interpolated solar and thermal source functions
c
      double precision up_s_src(kp),dn_s_src(kp),
     -                 rad_s_src(mxumu,mxphi,kp)
      double precision up_t_src(kp),dn_t_src(kp),
     -                 rad_t_src(mxumu,mxphi,kp)
c
      real dalb,deltau(kp),delpi0(kp),delg(kp)
c
      real alb0
      double precision wn_io
c
c****   specify the bin index
c
      ibn = ng0
c
      if(lsolar) then
c
c**********************************************************************
c
c*****         f i n d   s o l a r    s o u r c e s
c
c**********************************************************************
c
        if(l_1 .eq. 1) then
c
c****        set downward flux and radiance source functions
c            at the top boundary
c
          k = 1
          dn_s_src(k) = 0.0
          do 1021 naz=1,nphi
              do 1001 nze=1,nzdn
                  rad_s_src(nze,naz,k)= 0.0
1001          continue
1021      continue
c
        endif
c
        if(l_2 .eq. nlev) then
c
c****      set upward radiances at surface
c
          do 1061 naz=1,nphi
              do 1041 nze=nzup,numu
                  rad_s_src(nze,naz,nlev) = 
     -                       sradsrc_b(nze,naz,nlev,ibn) + 
     -                       dsraddx(nze,naz,nlev,4,ibn)*dalb
1041          continue
1061      continue
c
        endif
c
c****    find the downward radiance sources for all other levels
c
        do 1281 k=lev1,l2tau+1
c
c****              interpolate wavenumber-dependent downward fluxes 
c 
            dn_s_src(k) = dnsflxsrc_b(k,ibn) + 
     -                              ddnsflxdx(k,1,ibn)*deltau(k-1) + 
     -                              ddnsflxdx(k,2,ibn)*delpi0(k-1) + 
     -                              ddnsflxdx(k,3,ibn)*delg(k-1) + 
     -                              ddnsflxdx(k,4,ibn)*dalb
c
c****              interpolate wavenumber-dependent downward radiances 
c 
            do 1221 naz=1,nphi
                do 1201 nze=1,nzdn
                    rad_s_src(nze,naz,k) = sradsrc_b(nze,naz,k,ibn) + 
     -                       dsraddx(nze,naz,k,1,ibn)*deltau(k-1) +
     -                       dsraddx(nze,naz,k,2,ibn)*delpi0(k-1) +
     -                       dsraddx(nze,naz,k,3,ibn)*delg(k-1) + 
     -                       dsraddx(nze,naz,k,4,ibn)*dalb
1201            continue
1221        continue
c
c****        interpolate upward radiances and fluxes
c
            up_s_src(k-1) = upsflxsrc_b(k-1,ibn) + 
     -                      dupsflxdx(k-1,1,ibn)*deltau(k-1) + 
     -                      dupsflxdx(k-1,2,ibn)*delpi0(k-1) + 
     -                      dupsflxdx(k-1,3,ibn)*delg(k-1) + 
     -                      dupsflxdx(k-1,4,ibn)*dalb
c
c****           interpolate wavenumber-dependent upward radiances 
c 
            do 1261 naz=1,nphi
                do 1241 nze=nzup,numu
                    rad_s_src(nze,naz,k-1) = 
     -                        sradsrc_b(nze,naz,k-1,ibn) +
     -                         dsraddx(nze,naz,k-1,1,ibn)*deltau(k-1) +
     -                         dsraddx(nze,naz,k-1,2,ibn)*delpi0(k-1) +
     -                         dsraddx(nze,naz,k-1,3,ibn)*delg(k-1) +
     -                         dsraddx(nze,naz,k-1,4,ibn)*dalb
1241            continue      
1261        continue      
c
1281    continue
c
      endif
c
      if(lplanck) then
c
c**********************************************************************
c
c*****         f i n d   t h e r m a l    s o u r c e s
c
c**********************************************************************
c
        if(l_1 .eq. 1) then
c
c****        set downward flux and radiance source functions
c            at the top boundary
c
          k = 1
          dn_t_src(k) = 0.0
          do 2221 naz=1,nphi
              do 2201 nze=1,nzdn
                  rad_t_src(nze,naz,k)= 0.0
2201          continue
2221      continue
c
        endif
c
        if(l_2 .eq. nlev) then
c
c****      set upward fluxes and radiances at surface
c
          k = nlev
c
          up_t_src(k) = uptflxsrc_b(k,ibn) + duptflxdx(k,4,ibn)*dalb
c     
          do 2461 naz=1,nphi
              do 2441 nze=nzup,numu
                  rad_t_src(nze,naz,k) = tradsrc_b(nze,naz,k,ibn) + 
     -                                   dtraddx(nze,naz,k,4,ibn)*dalb
2441          continue
2461      continue
c
        endif
c
c****    find the downward radiance sources for all other levels
c
        do 2681 k=lev1,l2tau+1
c
c****              interpolate wavenumber-dependent downward fluxes 
c 
            dn_t_src(k) = dntflxsrc_b(k,ibn) + 
     -                    ddntflxdx(k,1,ibn)*deltau(k-1) + 
     -                    ddntflxdx(k,2,ibn)*delpi0(k-1) + 
     -                    ddntflxdx(k,3,ibn)*delg(k-1) + 
     -                    ddntflxdx(k,4,ibn)*dalb
c
c****              interpolate wavenumber-dependent downward radiances 
c 
            do 2621 naz=1,nphi
                do 2601 nze=1,nzdn
                    rad_t_src(nze,naz,k) = 
     -                           tradsrc_b(nze,naz,k,ibn) + 
     -                           dtraddx(nze,naz,k,1,ibn)*deltau(k-1) + 
     -                           dtraddx(nze,naz,k,2,ibn)*delpi0(k-1) + 
     -                           dtraddx(nze,naz,k,3,ibn)*delg(k-1) + 
     -                           dtraddx(nze,naz,k,4,ibn)*dalb
2601            continue
2621        continue
c
c****        interpolate upward radiances and flux sourceses
c
            up_t_src(k-1) = uptflxsrc_b(k-1,ibn) + 
     -                      duptflxdx(k-1,1,ibn)*deltau(k-1) + 
     -                      duptflxdx(k-1,2,ibn)*delpi0(k-1) + 
     -                      duptflxdx(k-1,3,ibn)*delg(k-1) + 
     -                      duptflxdx(k-1,4,ibn)*dalb
c
c****           interpolate wavenumber-dependent upward radiances 
c 
            do 2661 naz=1,nphi
                do 2641 nze=nzup,numu
                    rad_t_src(nze,naz,k-1) = 
     -                     tradsrc_b(nze,naz,k-1,ibn) +
     -                     dtraddx(nze,naz,k-1,1,ibn)*deltau(k-1) +
     -                     dtraddx(nze,naz,k-1,2,ibn)*delpi0(k-1) +
     -                     dtraddx(nze,naz,k-1,3,ibn)*delg(k-1) + 
     -                     dtraddx(nze,naz,k-1,4,ibn)*dalb
2641            continue      
2661        continue      
c
2681    continue

c       write(*,'(/,1a,2i5,1pe14.6)') 'interp_src: ibn,l_1, wn_io',
c     -  ibn,l_1,wn_io
c       write(*,'(1a,16(1pe12.4))') 'dn_t_src  ',
c     - (dn_t_src(k),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'dntflxsrc_b: ',
c     - (dntflxsrc_b(k,ibn),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'ddntf*dtau',
c     - ((ddntflxdx(k-1,1,ibn)*deltau(k-1)),k=2,nlev)
c       write(*,'(1a,16(1pe12.4))') 'ddntf*dpi0',
c     - ((ddntflxdx(k-1,2,ibn)*delpi0(k-1)),k=2,nlev)
c       write(*,'(1a,16(1pe12.4))') 'ddntf*dg  ',
c     - ((ddntflxdx(k-1,3,ibn)*delg(k-1)),k=2,nlev)
c       write(*,'(1a,16(1pe12.4))') 'ddntf*dalb',
c     - ((ddntflxdx(k-1,4,ibn)*dalb),k=2,nlev)
c  
c       write(*,'(/,1a,16(1pe12.4))') 'up_t_src: ',
c     - (up_t_src(k),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'uptflxsrc_b: ',
c     - (uptflxsrc_b(k,ibn),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'duptflx*dtau',
c     - ((duptflxdx(k-1,1,ibn)*deltau(k-1)),k=2,nlev)
c       write(*,'(1a,16(1pe12.4))') 'duptflx*dpi0',
c     - ((duptflxdx(k-1,2,ibn)*delpi0(k-1)),k=2,nlev)
c       write(*,'(1a,16(1pe12.4))') 'duptflx*dg  ',
c     - ((duptflxdx(k-1,3,ibn)*delg(k-1)),k=2,nlev)
c       write(*,'(1a,16(1pe12.4))') 'duptflx*dalb',
c     - ((duptflxdx(k-1,4,ibn)*dalb),k=2,nlev)
c
      endif
c
      return
      end
