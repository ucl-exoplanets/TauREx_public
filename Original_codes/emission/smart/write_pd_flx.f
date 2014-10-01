      subroutine write_pd_flx(iuo,npd0,ipd,nlev,l1,l2,ifrmout,
     -                        iunits,isptype,islit,nfpd,iquit_pd,
     -                        wn,wnmin,wnmax,width,dwn,points,
     -                        umu0nz,var,p,x0,pd_pert,
     -                        up_flx,dn_flx,up_flx_trn,dn_flx_trn,
     -                        trnflx1,trnflx2,trnflx3,
     -                        pd_trnflx1,pd_trnflx2,pd_trnflx3,
     -                        dn_src,up_src,pd_dn_src,pd_up_src)
c
ccccccccccccccccccccc  w r i t e _ p d _ f l x   ccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine prints partial derivatives for fluxes           cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        iuo - unit number for output level-dependent partials       cc
cc     iunits - index of output radiance units:                       cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) ergs/s/cm**2/sr/cm-1                               cc
cc              5) photons/s/m**2/sr/micron                           cc
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
ccccccccccccccccccccc  w r i t e _ p d _ f l x   ccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ncd,nfmx
      parameter (ncd=mxrad*mxrad+mxumu*mxphi+6, 
     -           nfmx = (mxlout+mxpd+1)*nsol+1)
c
      integer iuo,ifrmout,npd0,ipd(mxpd),nlev,l1(mxpd),l2(mxpd)
      integer iunits,isptype,islit,nfpd,iquit_pd(nfmx)
      integer k,l,n,nn,iord,iscwt,l2t
c
      double precision wn,wnmin,wnmax,width,dwn
c
      real var,umu0nz
      real p(kp),x0(kp,mxpd)
      real pd_pert(mxpd)
c
      real spect(ncd),points(nfmx)
c
      real wl,wnio,wnmins,wnmaxs
c
c****    flux source functions and jacobians
c
      real dn_src(mxrad),up_src(mxrad),
     -     pd_dn_src(mxrad,mxpd),pd_up_src(mxrad,mxpd)
c
c****    flux transmission and absorption functions and jacobians
c
      real trnflx1(mxrad),trnflx2(mxrad),trnflx3(mxrad),
     -     pd_trnflx1(mxrad,mxpd),pd_trnflx2(mxrad,mxpd),
     -     pd_trnflx3(mxrad,mxpd)
c
      real dn_flx(mxrad),up_flx(mxrad)
      real dn_flx_trn(mxrad),up_flx_trn(mxrad)
c
      n = npd0
c
      if(wn .eq. wnmin) then
c
c****    write a header, including number of model levels (nlev), 
c        the range of levels used for this jacobian (l1,l2), 
c        the jacobian type index (ipd), the units index (iunits)
c
        if(ifrmout .eq. 1) then
          write(iuo,'(5i10)') nlev,l1(n),l2(n),ipd(n),iunits
c
c****       write the amplitude of the perturbation used for this
c           jacobian, pd_pert, and the nominal value of this state 
c           vector element at each level where the jacobian is defined
c
          write(iuo,'(9(1pe14.6))')
     -                        umu0nz,pd_pert(n),(x0(k,n),k=l1(n),l2(n))
c
c****       write the mininum and maximum wavenumber incldued in the
c           jacobian file, and the pressure at each model level
c
          write(iuo,'(9(1pe14.6))') wnmin,wnmax,(p(k),k=1,nlev)
        else
          wnmins = real(wnmin)
          wnmaxs = real(wnmax)
          write(iuo) nlev,l1(n),l2(n),ipd(n),iunits
          write(iuo) umu0nz,pd_pert(n),(x0(k,n),k=l1(n),l2(n))
          write(iuo) wnmins,wnmaxs,(p(k),k=1,nlev)
        endif
c
        if(l2(n) .ge. nlev) then
          l2t = nlev - 1
        else
          l2t = l2(n)
        endif
c
      endif
c
c*****  define the wavelength
c
      wnio = real(wn)
      if(wn .ne. 0.0d0) then
        wl = real(1.0d4/wn)
      else
        wl = 1.e6
      endif
c
      if(isptype .eq. 1) then    
c
c****     create output at full monochromatic resolution
c
        if(ifrmout .eq. 1) then
c
c****        create ascii formatted output
c
          write(iuo,'(9(1pe14.6))') wl,wn,var,
     -                 (dn_src(k),k=1,nlev),
     -                 (up_src(k),k=1,nlev),
     -                 (pd_dn_src(l,n),l=l1(n),l2(n)),
     -                 (pd_up_src(l,n),l=l1(n),l2(n)),
     -                 (dn_flx(k),k=1,nlev-1),
     -                 (up_flx(k),k=1,nlev-1),
     -                 (dn_flx_trn(k),k=1,nlev-1),
     -                 (up_flx_trn(k),k=1,nlev-1),
     -                 (trnflx1(k),k=1,nlev-1),
     -                 (trnflx2(k),k=1,nlev-1),
     -                 (trnflx3(k),k=1,nlev-1),
     -                 (pd_trnflx1(l,n),l=l1(n),l2t),
     -                 (pd_trnflx2(l,n),l=l1(n),l2t),
     -                 (pd_trnflx3(l,n),l=l1(n),l2t)
c
        else
c
c****        create binary output
c
          write(iuo) wl,wnio,var,
     -                 (dn_src(k),k=1,nlev),
     -                 (up_src(k),k=1,nlev),
     -                 (pd_dn_src(l,n),l=l1(n),l2(n)),
     -                 (pd_up_src(l,n),l=l1(n),l2(n)),
     -                 (dn_flx(k),k=1,nlev-1),
     -                 (up_flx(k),k=1,nlev-1),
     -                 (dn_flx_trn(k),k=1,nlev-1),
     -                 (up_flx_trn(k),k=1,nlev-1),
     -                 (trnflx1(k),k=1,nlev-1),
     -                 (trnflx2(k),k=1,nlev-1),
     -                 (trnflx3(k),k=1,nlev-1),
     -                 (pd_trnflx1(l,n),l=l1(n),l2t),
     -                 (pd_trnflx2(l,n),l=l1(n),l2t),
     -                 (pd_trnflx3(l,n),l=l1(n),l2t)
c
        endif
c
      else
c
c****     increment counting variables
c
        points(nfpd) =  points(nfpd) + 1.0
c
        nn = 1
        spect(nn) = var
c
c****      downward and upward sources and their jacobians
c
        do 2001 k=1,nlev
            nn = nn + 1
            spect(nn) = dn_src(k)
2001    continue
        do 2021 k=1,nlev
            nn = nn + 1
            spect(nn) = up_src(k)
2021    continue
        do 2041 l=l1(n),l2(n)
            nn = nn + 1
            spect(nn) = pd_dn_src(l,n)
2041    continue
        do 2061 l=l1(n),l2(n)
            nn = nn + 1
            spect(nn) = pd_up_src(l,n)
2061    continue
c
c***      product of flux and transmission
c
        do 2101 k=1,nlev-1
            nn = nn + 1
            spect(nn) = dn_flx(k)
2101    continue
        do 2121 k=1,nlev-1
            nn = nn + 1
            spect(nn) = up_flx(k)
2121    continue
        do 2141 k=1,nlev-1
            nn = nn + 1
            spect(nn) = dn_flx_trn(k)
2141    continue
        do 2161 k=1,nlev-1
            nn = nn + 1
            spect(nn) = up_flx_trn(k)
2161    continue
c
c****     transmission variables and their jacobians
c
        do 2221 k=1,nlev-1
            nn = nn + 1
            spect(nn) = trnflx1(k)
2221    continue
        do 2241 k=1,nlev-1
            nn = nn + 1
            spect(nn) = trnflx2(k)
2241    continue
        do 2261 k=1,nlev-1
            nn = nn + 1
            spect(nn) = trnflx3(k)
2261    continue
        do 2341 l=l1(n),l2t
            nn = nn + 1
            spect(nn) = pd_trnflx1(l,n)
2341    continue
        do 2361 l=l1(n),l2t
            nn = nn + 1
            spect(nn) = pd_trnflx2(l,n)
2361    continue
        do 2381 l=l1(n),l2t
            nn = nn + 1
            spect(nn) = pd_trnflx3(l,n)
2381    continue
c
        iord = 2
        iscwt = 0
c
        call rad_slit(iuo,nfpd,ifrmout,
     -                islit,iscwt,iord,nn,
     -                wnmin,wnmax,dwn,width,
     -                points,wn,spect,iquit_pd(nfpd))
c
      endif
c
      return
      end
