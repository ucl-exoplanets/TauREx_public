      subroutine skip_wn(iwnmode,ij0,nlay,iwnflg,ntau_pd,numu,irad,nza,
     -                   dwnmx,small_tau,small_pi0,small_g,small_alb,
     -                   tauerr,pi0err,phferr,surferr,
     -                   wn,dtauex,dtausc,co_pi0,g,phmom,surf_opt,alb,
     -                   wni,dtauexi,dtausci,co_pi0i,gi,phmomi,
     -                   surf_opti,albi,tau_ext,tau_sca,g_sca,
     -                   tau_exti,tau_scai,g_scai,
     -                   p_ray_0,p_ray,tau_ray,p_ray_0i,p_rayi,tau_rayi,
     -                   p_gas_0,p_gas,tau_gas,p_gas_0i,p_gasi,tau_gasi,
     -                   p_aer_0,p_aer,tau_aer,p_aer_0i,p_aeri,tau_aeri)
c
ccccccccccccccccccccccccc    s k i p _ w n    cccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine determines whether the optical properties at    cc
cc    a spectral grid point will make a contribution to the spectrum, cc
cc    when compared to its neighboring points.  If not, it is skipped.cc 
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                      nza                              cc
cc    iwnmode - operating mode:                                       cc
cc              (0) setup                                             cc
cc              (1) check to see if point x(2) is needed              cc
cc              (2) all 3 points are good                             cc
cc              (3) overwrite point x(2) with point x(3)              cc
cc        ij0 - wn index counter (1,2,3)                              cc
cc       nlay - number of levels (nlev or nlev+1 if pressure is a     cc
cc              variable part of the state vector                     cc
cc     iwnflg - flag stating if this point is needed:                 cc
cc              (0) no, (1) yes                                       cc
cc      dwnmx -  maximum spectral interval between output points      cc
cc  small_tau - small optical depth                                   cc
cc  small_pi0 - small single scattering albedo                        cc
cc    small_g - small asymmetry parameter                             cc
cc     tauerr - optical depth relative binning error (0. to ~0.8)     cc
cc     pi0err - co-single scattering albedo absolute binning error    cc
cc     phferr - asymmetry factor absolute binning error               cc
cc    surferr - surface optical property binning error                cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    input values, reset to skip values at specified wavenumber      cc
cc                                                                    cc
ccccccccccccccccccccccccc    s k i p _ w n    cccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ij0,iwnflg,iwnmode,nlay,ntau_pd,numu,nza
      integer k,l,mom,n,nze,irad,m,nz
c
c****   atmospheric optical properties
c
      real dtauex(kp,2),co_pi0(kp,2),g(kp),phmom(0:mxmom,kp),
     -     surf_opt(4),alb(nsol)
c
c****   atmospheric scattering optical depth and surface properties
c
      real dtausc(kp)
c
      double precision wn,wni(3),dwni,dwnmx
c
c****    pressure of tau = 1 and column optical depths
c
      real p_ray_0,p_gas_0,p_aer_0
      real tau_ray,tau_gas,tau_aer
      real p_ray(mxumu),p_gas(mxumu),p_aer(mxumu)
c
c****    optical depths of each variable component of state vector
c
      real tau_ext(kp,3,mxpd),tau_sca(kp,mxpd),g_sca(kp,mxpd)
c
      real dtauexi(kp,2,3),dtausci(kp,3),co_pi0i(kp,2,3),gi(kp,3),
     -     phmomi(0:mxmom,kp,3),surf_opti(4,3),albi(nsol,3),
     -     small_tau,small_pi0,small_g,small_alb
c
      real tau_exti(kp,3,mxpd,3),tau_scai(kp,mxpd,3),g_scai(kp,mxpd,3)
c
      real p_ray_0i(3),p_rayi(mxumu,3),tau_rayi(3),
     -     p_gas_0i(3),p_gasi(mxumu,3),tau_gasi(3),
     -     p_aer_0i(3),p_aeri(mxumu,3),tau_aeri(3)
c
      real tauerr,pi0err,phferr,surferr
c
      double precision diff
c
c****    check to see if this is a setup step (iwnmode = 0),
c        a wavelength skip test is to be performed (iwnmode = 1),
c        or a wavelength skip step is performed (iwnmode = 2)
c
        if(iwnmode .eq. 1) then
c
c****       determine if the optical properties are changing
c           rapidly enough to justify an RT calculation for this point
c
          wni(ij0) = wn
c
          if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -          .or. irad .eq. 8) then
            p_ray_0i(ij0) = p_ray_0
            p_gas_0i(ij0) = p_gas_0
            p_aer_0i(ij0) = p_aer_0
            tau_rayi(ij0) = tau_ray
            tau_gasi(ij0) = tau_gas
            tau_aeri(ij0) = tau_aer
            do 2001 nze=1,numu
                p_rayi(nze,ij0) = p_ray(nze)
                p_gasi(nze,ij0) = p_gas(nze)
                p_aeri(nze,ij0) = p_aer(nze)
2001        continue
c
          endif
c
          do 2061 k=1,nlay
              dtauexi(k,1,ij0) = dtauex(k,1)
              dtauexi(k,2,ij0) = dtauex(k,2)
              dtausci(k,ij0) = dtausc(k)
              co_pi0i(k,1,ij0) = co_pi0(k,1)
              co_pi0i(k,2,ij0) = co_pi0(k,2)
              gi(k,ij0) = g(k)
              do 2021 mom=0,mxmom
                  phmomi(mom,k,ij0) = phmom(mom,k)
2021          continue
              do 2041 n=1,ntau_pd
                  tau_exti(k,1,n,ij0) = tau_ext(k,1,n)
                  tau_exti(k,2,n,ij0) = tau_ext(k,2,n)
                  tau_exti(k,3,n,ij0) = tau_ext(k,3,n)
                  tau_scai(k,n,ij0) = tau_sca(k,n)
                  g_scai(k,n,ij0) = g_sca(k,n)
2041          continue
2061      continue
          do 2081 l=1,4
              surf_opti(l,ij0) = surf_opt(l)
2081      continue 
          do 2091 nz=1,nza
              albi(nz,ij0) = alb(nz)
2091      continue 
c
          if(ij0 .lt. 3) return
c
          iwnflg = 0
c
          if(wni(3) - wni(1) .gt. dwnmx) iwnflg = 1
          dwni = (wni(2) - wni(1))/(wni(3) - wni(1))
          do 3021 k=1,nlay
              diff = dtauexi(k,1,2) - (dtauexi(k,1,1) + 
     -                 (dtauexi(k,1,3) - dtauexi(k,1,1))*dwni)
              if(abs(diff) .gt. 
     -             0.01*tauerr*(dtauexi(k,1,2) + small_tau)) iwnflg = 1
c
              diff = dtauexi(k,2,2) - (dtauexi(k,2,1) + 
     -                 (dtauexi(k,2,3) - dtauexi(k,2,1))*dwni)
              if(abs(diff) .gt. 
     -             0.01*tauerr*(dtauexi(k,2,2) + small_tau)) iwnflg = 1
c
              diff = dtausci(k,2) - (dtausci(k,1) + 
     -                 (dtausci(k,3) - dtausci(k,1))*dwni)
              if(abs(diff) .gt. 
     -             0.01*tauerr*(dtausci(k,2) + small_tau)) iwnflg = 1
c
              diff = co_pi0i(k,1,2) - (co_pi0i(k,1,1) + 
     -                 (co_pi0i(k,1,3) - co_pi0i(k,1,1))*dwni)
              if(abs(diff) .gt. 
     -             0.01*pi0err*(co_pi0i(k,1,2) + small_pi0)) iwnflg = 1
c
              diff = co_pi0i(k,2,2) - (co_pi0i(k,2,1) + 
     -                 (co_pi0i(k,2,3) - co_pi0i(k,2,1))*dwni)
              if(abs(diff) .gt. 
     -             0.01*pi0err*(co_pi0i(k,2,2) + small_pi0)) iwnflg = 1
c
              diff = gi(k,2) - (gi(k,1) + 
     -                 (gi(k,3) - gi(k,1))*dwni)
              if(abs(diff) .gt. 
     -             0.01*phferr*(gi(k,2) + small_g)) iwnflg = 1
c
              do 3001 mom=0,mxmom
                  diff = phmomi(mom,k,2) - (phmomi(mom,k,1) + 
     -                     (phmomi(mom,k,3) - phmomi(mom,k,1))*dwni)
                  if(abs(diff) .gt. 
     -                 0.01*phferr*(phmomi(mom,k,2) + small_g)) 
     -                 iwnflg = 1
3001          continue
3021      continue
          do 3041 l=1,4
              surf_opti(l,ij0) = surf_opt(l)
              diff = surf_opti(l,2) - (surf_opti(l,1) + 
     -                 (surf_opti(l,3) - surf_opti(l,1))*dwni)
              if(abs(diff) .gt. 
     -             0.01*surferr*(surf_opti(l,2) + small_alb)) iwnflg = 1
3041      continue
          do 3061 nz=1,nza
              albi(nz,ij0) = alb(nz)
              diff = albi(nz,2) - (albi(nz,1) + 
     -                 (albi(nz,3) - albi(nz,1))*dwni)
              if(abs(diff) .gt. 
     -             0.01*surferr*(albi(l,2) + small_alb)) iwnflg = 1
3061      continue
c
        else
c
          if(iwnmode .eq. 2) then
c
c*****       All 3 points are good. put x(1) into working arrays
c            and process that point, and then pack x(2) into x(1)
c            and x(3) into x(2)
c
            wn = wni(1)
            wni(1) = wni(2)
            wni(2) = wni(3)
c
            if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -          .or. irad .eq. 8) then
              p_ray_0 = p_ray_0i(1)
              p_ray_0i(1) = p_ray_0i(2)
              p_ray_0i(2) = p_ray_0i(3)
              p_gas_0 = p_gas_0i(1)
              p_gas_0i(1) = p_gas_0i(2)
              p_gas_0i(2) = p_gas_0i(3)
              p_aer_0 = p_aer_0i(1)
              p_aer_0i(1) = p_aer_0i(2)
              p_aer_0i(2) = p_aer_0i(3)
              tau_ray = tau_rayi(1)
              tau_rayi(1) = tau_rayi(2)
              tau_rayi(2) = tau_rayi(3)
              tau_gas = tau_gasi(1)
              tau_gasi(1) = tau_gasi(2)
              tau_gasi(2) = tau_gasi(3)
              tau_aer = tau_aeri(1)
              tau_aeri(1) = tau_aeri(2)
              tau_aeri(2) = tau_aeri(3)
              do 4001 nze=1,numu
                  p_ray(nze) = p_rayi(nze,1)
                  p_rayi(nze,1) = p_rayi(nze,2)
                  p_rayi(nze,2) = p_rayi(nze,3)
                  p_gas(nze) = p_gasi(nze,1)
                  p_gasi(nze,1) = p_gasi(nze,2)
                  p_gasi(nze,2) = p_gasi(nze,3)
                  p_aer(nze) = p_aeri(nze,1)
                  p_aeri(nze,1) = p_aeri(nze,2)
                  p_aeri(nze,2) = p_aeri(nze,3)
4001          continue
c
            endif
c
            do 4061 k=1,nlay
                dtauex(k,1) = dtauexi(k,1,1)
                dtauexi(k,1,1) = dtauexi(k,1,2)
                dtauexi(k,1,2) = dtauexi(k,1,3)
                dtauex(k,2) = dtauexi(k,2,1)
                dtauexi(k,2,1) = dtauexi(k,2,2)
                dtauexi(k,2,2) = dtauexi(k,2,3)
                dtausc(k) = dtausci(k,1)
                dtausci(k,1) = dtausci(k,2)
                dtausci(k,2) = dtausci(k,3)
                co_pi0(k,1) = co_pi0i(k,1,1)
                co_pi0i(k,1,1) = co_pi0i(k,1,2)
                co_pi0i(k,1,2) = co_pi0i(k,1,3)
                co_pi0(k,2) = co_pi0i(k,2,1)
                co_pi0i(k,2,1) = co_pi0i(k,2,2)
                co_pi0i(k,2,2) = co_pi0i(k,2,3)
                g(k) = gi(k,1)
                gi(k,1) = gi(k,2)
                gi(k,2) = gi(k,3)
                do 4021 mom=0,mxmom
                    phmom(mom,k) = phmomi(mom,k,1)
                    phmomi(mom,k,1) = phmomi(mom,k,2)
                    phmomi(mom,k,2) = phmomi(mom,k,3)
4021            continue
                do 4041 n=1,ntau_pd
                    do 4031 m=1,3
                        tau_ext(k,m,n) = tau_exti(k,m,n,1)
                        tau_exti(k,m,n,1) = tau_exti(k,m,n,2)
                        tau_exti(k,m,n,2) = tau_exti(k,m,n,3)
4031                continue
                    tau_sca(k,n) = tau_scai(k,n,1)
                    tau_scai(k,n,1) = tau_scai(k,n,2)
                    tau_scai(k,n,2) = tau_scai(k,n,3)
                    g_sca(k,n) = g_scai(k,n,1)
                    g_scai(k,n,1) = g_scai(k,n,2)
                    g_scai(k,n,2) = g_scai(k,n,3)
4041            continue
4061        continue
            do 4081 l=1,4
                surf_opt(l) = surf_opti(l,1)
                surf_opti(l,1) = surf_opti(l,2)
                surf_opti(l,2) = surf_opti(l,3)
4081        continue
            do 4091 nz=1,nza
                alb(nz) = albi(nz,1)
                albi(nz,1) = albi(nz,2)
                albi(nz,2) = albi(nz,3)
4091        continue
c
          else
c
c****        iwnmode = 3: skip point x(2).  Load x(3) into x(2)
c
            wni(2) = wni(3)
c
            if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -          .or. irad .eq. 8) then
              p_ray_0i(2) = p_ray_0i(3)
              p_gas_0i(2) = p_gas_0i(3)
              p_aer_0i(2) = p_aer_0i(3)
              tau_rayi(2) = tau_rayi(3)
              tau_gasi(2) = tau_gasi(3)
              tau_aeri(2) = tau_aeri(3)
              do 5001 nze=1,numu
                  p_rayi(nze,2) = p_rayi(nze,3)
                  p_gasi(nze,2) = p_gasi(nze,3)
                  p_aeri(nze,2) = p_aeri(nze,3)
5001          continue
            endif
c
            do 5061 k=1,nlay
                dtauexi(k,1,2) = dtauexi(k,1,3)
                dtauexi(k,2,2) = dtauexi(k,2,3)
                dtausci(k,2) = dtausci(k,3)
                co_pi0i(k,1,2) = co_pi0i(k,1,3)
                co_pi0i(k,2,2) = co_pi0i(k,2,3)
                gi(k,2) = gi(k,3)
                do 5021 mom=0,mxmom
                    phmomi(mom,k,2) = phmomi(mom,k,3)
5021            continue
                do 5041 n=1,ntau_pd
                    tau_exti(k,1,n,2) = tau_exti(k,1,n,3)
                    tau_exti(k,2,n,2) = tau_exti(k,2,n,3)
                    tau_exti(k,3,n,2) = tau_exti(k,3,n,3)
                    tau_scai(k,n,2) = tau_scai(k,n,3)
                    g_scai(k,n,2) = g_scai(k,n,3)
5041            continue
5061        continue
            do 5081 l=1,4
                surf_opti(l,2) =  surf_opti(l,3)
5081        continue
            do 5091 nz=1,nza
                albi(nz,2) =  albi(nz,3)
5091        continue
c
          endif
c
        endif
c
      return
      end
          
