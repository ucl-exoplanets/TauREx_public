      subroutine pd_out(nz0,npd,npd1,ipd,igs,ifrmout,
     -                    levout,nlout,nlyr,nphi,numu,nzup,
     -                    iu_pd,iutpd,iuspd,modepd,iupdrad,
     -                    iunits,isptype,islit,lsolar,lplanck,
     -                    wn,wnmin,wnmax,width,dwn,points,units,
     -                    umu0nz,ts,t,p,rmix,alb,dtauaer,
     -                    solflx,up_s_flx,dn_s_flx,dir_s_flx,
     -                    up_t_flx,dn_t_flx,pd_pert,
     -                    dns_src,pd_dns_src,ups_src,pd_ups_src,
     -                    dnt_src,pd_dnt_src,upt_src,pd_upt_src,
     -                    trn_dir,trn_flx,ref_flx,abs_flx,
     -                    pd_trndir,pd_trnflx,pd_refflx,pd_absflx,
     -                    rad,pd_rad)
c
cccccccccccccccccccccccccccc   p d _ o u t   ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine prints monochromatic partial derivatives,       cc
cc    or values smoothed with a slit function.                        cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc     iunits - index of output radiance units:                       cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) ergs/s/cm**2/sr/cm-1                               cc
cc              5) photons/s/m**2/sr/micron                           cc
cc      iutpd - unit numbers for output level-dependent               cc
cc              thermal fluxes and thier partial derivatives          cc
cc      iuspd - unit numbers for output level-dependent               cc
cc              solar fluxes and thier partial derivatives            cc
cc    iupdrad - unit numbers for output level-dependent radiances     cc
cc              and thier partial derivatives                         cc
cc         nz - index of this solar zenith angle (1 to nsol)          cc
cc    ifrmout - output file format (1) ascii or (2) binary            cc
cc              (2) surface, (3) arbitrary level                      cc
cc       nphi - number of output azimuth angles                       cc
cc       numu - number of zenith angles in input file                 cc
cc    isptype - output spectrum type:                                 cc
cc              1) Full-resolution spectrum                           cc
cc              2) Sampled spectrum smoothed with a slit function     cc
cc      islit - index of slit function type:                          cc
cc              1) boxcar                                             cc
cc              2) triangular (approximate slit spectrometer)         cc
cc      width - half-width at half max of slit function  (cm**-1)     cc
cc        dwn - output sampling resolution (cm**-1)                   cc
cc         wn - current wavenumber (cm**-1)                           cc
cc      wnmin - minimum wavenumber in output spectral grid (cm**-1)   cc
cc      wnmax - maximum wavenumber in output spectral grid (cm**-1)   cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc         wl - current wavenumber (cm**-1)                           cc
cc       wnio - current wavenumber (cm**-1)                           cc
cc     pd_rad - radiance partial derivatives at wn                    cc
cc                                                                    cc
cccccccccccccccccccccccccccc   p d _ o u t   ccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ncd,nfmx
      parameter (ncd=mxrad*mxrad+mxumu*mxphi+6, 
     -           nfmx = (mxlout+mxpd+1)*nsol+1)
c
      logical lsolar,lplanck
c
      integer npd1,igs(mxpd),modepd(mxpd)
      integer ifrmout,nlout,nlyr,nphi,numu,
     -         iunits,isptype,iscwt,islit,levout,nzup
      integer nz0,nlev,nfpd
      integer iquit_pd(nfmx)
      save iquit_pd
      integer nn,k,l,n,nl,naz,nze,npd0,l1(mxpd),l2(mxpd)
      integer iuo,iord
      integer nza_1(mxlout)
c
c****    units for flux and radiance jacobians
c 
      integer npd,ipd(mxpd),
     -        iu_pd,iutpd(mxpd,nsol),iuspd(mxpd,nsol),
     -        iupdrad(mxumu,mxphi,mxpd,mxlout,nsol)
c
      double precision wn,wnmin,wnmax,width,dwn
c
      real solflx,umu0nz,units
      real ts,t(kp),p(kp),rmix(kp,ngas),alb(nsol)
c
      real pd_pert(mxpd),dtauaer(nmode,kp)
c
c****   slit function variables
c
      real spect(ncd),points(nfmx)
c
c****   layer-dependent source terms for solar and thermal fluxes
c
      real dns_src(mxrad),ups_src(mxrad),dnt_src(mxrad),upt_src(mxrad)
c
      real rad(mxumu,mxphi,mxlout)
c
c****    flux transmission and absorption functions at wavenuber, wn
c
      real trn_dir(mxrad),trn_flx(mxrad),ref_flx(mxrad),
     -     abs_flx(mxrad)
c
c****   output flux and radiance partial derivatives
c
      real pd_ups_src(mxrad,mxpd),pd_dns_src(mxrad,mxpd),
     -     pd_upt_src(mxrad,mxpd),pd_dnt_src(mxrad,mxpd),
     -     pd_rad(mxumu,mxphi,mxlout,mxrad,mxpd)
c
c****    output fluxes and radiances
c
      real dir_s_flx(mxrad,nsol,2),up_s_flx(mxrad,nsol,2),
     -     dn_s_flx(mxrad,nsol,2),up_t_flx(mxrad,2),dn_t_flx(mxrad,2)
c
      real dn_flx_trn(mxrad),up_flx_trn(mxrad),
     -     dn_flx(mxrad),up_flx(mxrad)
c
c****    flux transmission and absorption partial 
c        derivatives for simplified adding method at wavenumber wn
c
      real pd_trndir(mxrad,mxpd),
     -     pd_trnflx(mxrad,mxpd),
     -     pd_refflx(mxrad,mxpd),
     -     pd_absflx(mxrad,mxpd)
c
      real wl,var,x0(kp,mxpd),rad0,pd_rad0(mxrad)
      real wnio,wnmins,wnmaxs
c
c****   specify the file index for the first output file 
c       (required for slit function programs)
c
      nlev = nlyr + 1
c
c****   set io variables
c
      if(isptype .eq. 2) then
        do 1001 n=1,mxpd
            if(wn .le. wnmin) then
              iquit_pd(n)= 0
            else
              if(wn .ge. wnmax .and. iquit_pd(n) .eq. 0) 
     -           iquit_pd(n) = 1
            endif
1001    continue
      endif
c
c*****  define the single precision wavenumber and wavelength
c
      wnio = real(wn)
      if(wn .ne. 0.0) then
        wl = real(1.e4/wn)
      else
        wl = 1.e6
      endif
c
c****   initialize state vector varibles
c
      do 1081 n=1,npd
c
          if(abs(ipd(n)) .eq. 1) then
            l1(n) = nlev
            l2(n) = nlev
            x0(nlev,n) = p(nlev)
          else
            if(abs(ipd(n)) .eq. 2) then
              l1(n) = 1
              l2(n) = nlev
              do 1021 l=l1(n),l2(n)
                  x0(l,n) = t(l)
1021          continue
              x0(nlev+1,n) = ts
            else
              if(abs(ipd(n)) .eq. 3) then
                l1(n) = 1
                l2(n) = nlev
                do 1041 l=l1(n),l2(n)
                    x0(l,n) = rmix(l,igs(n-npd1))
1041            continue
              else
                if(abs(ipd(n)) .eq. 4) then
                  l1(n) = 1
                  l2(n) = nlev-1
                  do 1061 l=l1(n),l2(n)
                      x0(l,n) = dtauaer(modepd(n-npd1),l)
1061              continue
                else
                  l1(n) = nlev
                  l2(n) = nlev
                  x0(nlev,n) = alb(nz0)
                endif
              endif
            endif
          endif
c
1081  continue
c
c****   F l u x    J a c o b i a n s
c
      do 2201 n=1,npd
c
          npd0 = n
c
          if(ipd(n) .lt. 0) then
c
            if(lplanck .and. nz0 .eq. 1) then
c
c****          write out thermal flux partials at each level
c
              nfpd = iutpd(n,nz0) - iu_pd + 5
              var = alb(nz0)
c
              do 2001 k=1,nlev
                  dnt_src(k) = units*dnt_src(k)
                  upt_src(k) = units*upt_src(k)
2001          continue
c
              do 2021 k=l1(n),l2(n)
                    pd_dnt_src(k,n) = units*pd_dnt_src(k,n)
                    pd_upt_src(k,n) = units*pd_upt_src(k,n)
2021          continue     
c
c***            corrections for transmission flux products
c
              do 2041 k=1,nlev-1
                  dn_flx(k) = dn_t_flx(k,2)
                  up_flx(k) = up_t_flx(k+1,2)
                  dn_flx_trn(k) = dn_t_flx(k,2)*trn_flx(k)
                  up_flx_trn(k) = up_t_flx(k+1,2)*trn_flx(k)
2041          continue   
c
c****            write thermal fluxes and their jacobians
c              
              call write_pd_flx(
     -                        iutpd(n,nz0),npd0,ipd,nlev,l1,l2,ifrmout,
     -                        iunits,isptype,islit,nfpd,iquit_pd,
     -                        wn,wnmin,wnmax,width,dwn,points,
     -                        umu0nz,var,p,x0,pd_pert,
     -                        up_flx,dn_flx,up_flx_trn,dn_flx_trn,
     -                        trn_flx,ref_flx,abs_flx,
     -                        pd_trnflx,pd_refflx,pd_absflx,
     -                        dnt_src,upt_src,pd_dnt_src,pd_upt_src)
c
            endif
c
            if(lsolar) then
c
c****           write downward diffuse solar flux partials at each level
c
              nfpd = iuspd(n,nz0) - iu_pd + 5 + nsol
              var = units*solflx
c
              do 2101 k=1,nlev
                  dns_src(k) = units*dns_src(k)
                  ups_src(k) = units*ups_src(k)
2101          continue
c
              do 2121 k=l1(n),l2(n)
                    pd_dns_src(k,n) = units*pd_dns_src(k,n)
                    pd_ups_src(k,n) = units*pd_ups_src(k,n)
2121          continue
c
              do 2141 k=1,nlev-1
                  dn_flx(k) = dn_s_flx(k,nz0,2) + dir_s_flx(k,nz0,2)
                  up_flx(k) = up_s_flx(k+1,nz0,2)
                  dn_flx_trn(k) = dn_s_flx(k,nz0,2)*trn_flx(k)
                  up_flx_trn(k) = up_s_flx(k+1,nz0,2)*trn_flx(k)
2141          continue   
c
c****            write solar fluxes and their jacobians
c              
              call write_pd_flx(
     -                        iuspd(n,nz0),npd0,ipd,nlev,l1,l2,ifrmout,
     -                        iunits,isptype,islit,nfpd,iquit_pd,
     -                        wn,wnmin,wnmax,width,dwn,points,
     -                        umu0nz,var,p,x0,pd_pert,
     -                        up_flx,dn_flx,up_flx_trn,dn_flx_trn,
     -                        trn_flx,ref_flx,trn_dir,
     -                        pd_trnflx,pd_refflx,pd_trndir,
     -                        dns_src,ups_src,pd_dns_src,pd_ups_src)
            endif
c
          endif
c
2201  continue
c
c****    R a d i a n c e    J a c o b i a n s
c
c****    define the first output zenith angle for radiances
c        (only upward radiances are saved at the top fo the atmosphere)
c
      do 3001 nl=1,nlout
          if(levout .eq. 1 .or. (levout .eq. 3 .and. nl .eq. 1)) then
             nza_1(nl) = nzup
          else
             nza_1(nl) = 1
          endif
3001  continue
c
      do 3261 n=1,npd
c
          if(ipd(n) .gt. 0) then
c
            do 3241 nl=1,nlout
c
                do 3221 naz=1,nphi
c
                    do 3201 nze=nza_1(nl),numu
c  
c****                     set the unit number
c
                        iuo = iupdrad(nze,naz,n,nl,nz0)
c
                        if(wn .eq. wnmin) then
c
c*****                      print print out the pressure profile if 
c                           this is the first point in the file  
c
                          if(ifrmout .eq. 1) then
                            write(iuo,'(5i10)') 
     -                           nlev,l1(n),l2(n),ipd(n),iunits
                            write(iuo,'(2(1pe14.6),10(1pe12.4),/,
     -                                        10(12(1pe12.4),/))') 
     -                          pd_pert(n),(x0(l,n),l=l1(n),l2(n))
                            write(iuo,'(2(1pe14.6),10(1pe12.4),/,
     -                                        10(12(1pe12.4),/))')
     -                                  wnmin,wnmax,(p(k),k=1,nlev)
                          else
                            wnmins = real(wnmin)
                            wnmaxs = real(wnmax)
                            write(iuo) nlev,l1(n),l2(n),ipd(n),iunits
                            write(iuo) pd_pert(n),
     -                                 (x0(l,n),l=l1(n),l2(n))
                            write(iuo) wnmins,wnmaxs,(p(k),k=1,nlev)
                          endif
c
                        endif
c
c****                     print level-dependent partials at wavelength
c
                        if(isptype .eq. 1) then
c
c****                       write monochromatic data at full resolution
c
                          if(ifrmout .eq. 1) then
c
                            if(levout .eq. 1 .or. (levout .eq. 3 .and. 
     -                                             nl .eq. 1)) then
c
                              rad0 = units*rad(nze,naz,1)
                              do 3101 l=l1(n),l2(n)
                                  pd_rad0(l) = units*
     -                                         pd_rad(nze,naz,nl,l,n)
 3101                         continue
                              write(iuo,'(2(1pe14.6),10(1pe12.4),/,
     -                                        10(12(1pe12.4),/))') 
     -                           wl,wn,rad0,(pd_rad0(l),l=l1(n),l2(n))
                            else
c
                              if(levout .eq. 2 .or. (levout .eq. 3 .and.
     -                                           nl .eq. 2)) then
                                rad0 = units*rad(nze,naz,nl)
                                do 3121 l=l1(n),l2(n)
                                    pd_rad0(l) = units*
     -                                           pd_rad(nze,naz,nl,l,n)
3121                            continue
c
                                write(iuo,'(2(1pe14.6),10(1pe12.4),/,
     -                                          10(12(1pe12.4),/))') 
     -                           wl,wn,rad0,(pd_rad0(l),l=l1(n),l2(n))
                              endif
                            endif
                          else
c
                            if(levout .eq. 1 .or. (levout .eq. 3 .and. 
     -                                             nl .eq. 1)) then
c
                              rad0 = units*rad(nze,naz,1)
                              do 3141 l=l1(n),l2(n)
                                  pd_rad0(l) = units*
     -                                         pd_rad(nze,naz,nl,l,n)
 3141                         continue
                              write(iuo) wl,wnio,rad0,
     -                            (pd_rad0(l),l=l1(n),l2(n))
                            else
                              if(levout .eq. 2 .or. (levout .eq. 3 .and.
     -                                           nl .eq. 2)) then
                                rad0 = units*rad(nze,naz,nl)
                                do 3161 l=l1(n),l2(n)
                                    pd_rad0(l) = units*
     -                                           pd_rad(nze,naz,nl,l,n)
3161                            continue
c
                                write(iuo) wl,wnio,rad0,
     -                            (pd_rad0(l),l=l1(n),l2(n))
                              endif
                            endif
c
                          endif
c
                        else
c
c****                       convolve the data with a slit function.
c                           set slit function index
c
                          nfpd = iuo - iu_pd + 5 + nsol
c
c****                      increment counting variables
c
                          points(nfpd) =  points(nfpd) + 1.0
c
                          nn = 1
                          if(levout .eq. 1 .or. (levout .eq. 3 .and. 
     -                                           nl .eq. 1)) then
                            spect(nn) = units*rad(nze,naz,1)
                          else
                            if(levout .eq. 2 .or. (levout .eq. 3 .and. 
     -                                           nl .eq. 2))
     -                      spect(nn) = units*rad(nze,naz,nl)
                          endif
c
                          do 3181 l=l1(n),l2(n)
                              nn = nn + 1
                              spect(nn) = units*pd_rad(nze,naz,nl,l,n)
3181                      continue
c
                          iord = 2
                          iscwt = 0
c
c****                        call the slit function program
c
                          call rad_slit(iuo,nfpd,ifrmout,
     -                        islit,iscwt,iord,nn,
     -                        wnmin,wnmax,dwn,width,
     -                        points,wn,spect,iquit_pd(nfpd))
c

                        endif
c
3201                continue
3221            continue
3241        continue
c
          endif
c
3261  continue
c
      return
      end
