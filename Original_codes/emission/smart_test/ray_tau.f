      subroutine ray_tau(ne,nlyr,ncomp,icomp,volmix,p,grav,
     -                   wgtatm,nzup,nstr,iflext,nsiext,istate,
     -                   pd_frac,wnext,wn,umu,dtauex,dtausc,phmom,
     -                   tauray,p_ray,p_ray_0)
c
ccccccccccccccccccccccccccc  r a y _ t a u  cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine finds the Rayleigh scattering cross sections    cc
cc    at each model level, and computes the effective monochromatic   cc
cc    Rayleigh scattering optical depth for each atmospheric layer.   cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc       nlyr - number of layers in the model atmosphere         .    cc
cc      ncomp - number of rayleigh-scattering constituentes           cc
cc      icomp - index of each rayleigh scattering constituent         cc
cc         a0 - avogodro's number.                                    cc
cc     volmix - volume mixing ratio of each rayleigh scatterer        cc
cc       nzup - index of first upward radiance stream.                cc
cc          p - pressure at each model level (Pascals).               cc
cc       grav - gravitational acceleration at each level (m/s^2).     cc
cc       nstr - number of computational zenith angles (streams).      cc
cc        umu - cosine of each upward radiance stream.                cc
cc         wn - current wavenumber (cm**-1).                          cc
cc     iflext - index for extinction sources                          cc
cc     nsiext - index counter for extinction sources                  cc
cc      wnext - wavenumber of the extincition source                  cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     tauray - column-integrated rayleigh optical depth.             cc
cc     dtauex - monochromatic extinction optical depth in each layer. cc
cc     dtausc - monochromatic scattering optical depth in each layer. cc
cc      phmom - scattering phase function moments in each layer.      cc
cc      p_ray - pressure level of rayliegh tau=1 for each upward      cc
cc             stream (bars).                                         cc
cc    p_ray_0 - pressure of Rayleigh tau=1 for a vertical stream.     cc
cc     iflext - index for extinction sources                          cc
cc     nsiext - index counter for extinction sources                  cc
cc      wnext - wavenumber of the extincition source                  cc
cc                                                                    cc
ccccccccccccccccccccccccccc  r a y _ t a u  cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer icomp(6),ne,nlyr,ncomp,nstr,nzup,istate(nex)
c
c*****   wavenumber grid variables
c
      integer iflext(nex),nsiext(2,nex)
      integer nt_pd,k,l,ktau1,nze,nlay
c
      double precision wn,wnext(2,nex)
c
      real wgtatm
      real a0,wl,sigray,dtray,dtau_tst,ptauray
      real raylei
c
c****   partial derivative fractional change
c
      real pd_frac(nex)
c
c****  rayleigh scattering optical properties
c
      real tauray(kp),dtauex(kp,2),dtausc(kp),umu(mxumu),
     -     p_ray(mxumu),p_ray_0,volmix(ncomp),p(kp),grav(kp),
     -     phmom(0:mxmom,kp)
c
c****    define avagadro's number (mks units: molecules/kmole)
c
      data a0/6.02297e26/
c
c*****   determine whether temperature is a variable part of
c        the state vector.  if it is, extinction optical depth must be 
c        computed for the background and perterbed temperature profile
c
      if(istate(2) .eq. 0) then
        nt_pd = 1
      else
        nt_pd = 2
      endif
c
c****   define the local wavelength (microns)
c
      wl = real(1.0d4/wn)
c
c****       find the rayleigh scattering cross section at wl
c
      sigray = raylei(wl,ncomp,icomp,volmix)
      tauray(1) = 0.0
c
      do 4021 k=1,nlyr
c
c****         find rayleigh scattering optical depth
c
          dtray = a0*sigray*(p(k+1) - p(k))/
     -            (wgtatm*0.5*(grav(k)+grav(k+1)))
          tauray(k+1) = tauray(k) + dtray
          dtausc(k) = dtausc(k) + dtray
          do 4001 l=1,nt_pd
              dtauex(k,l) = dtauex(k,l) + dtray
4001      continue
c
c****        define rayleigh scattering phase function moments
c
          phmom(0,k) = phmom(0,k) + dtray
          phmom(1,k) = phmom(1,k) + 0.1*dtray 
4021  continue
c
c****    if surface pressure is a variable part of the state vector,
c        find the optical depth for a layer below the base of the 
c        atmosphere that is 1% thicker than the nominal layer
c
      if(istate(1) .ne. 0) then
        nlay = nlyr + 1
        k = nlay
        dtray = a0*sigray*((1.0 + pd_frac(1))*p(k) - p(k-1))/
     -            (wgtatm*0.5*(grav(k) + grav(k-1)))
        dtausc(k) = dtausc(k) + dtray
        dtauex(k,1) = dtauex(k,1) + dtray
c      write(*,'(1a,i5,3(1pe14.6))') 'ray_tau: nlay,pd_frac,dtauex:',
c     -      nlay,pd_frac(1),dtauex(nlyr,1),dtauex(nlay,1)      
c
c****        define rayleigh scattering phase function moments
c
          phmom(0,k) = phmom(0,k) + dtray
          phmom(1,k) = phmom(1,k) + 0.1*dtray 
      endif
c
c****       find pressure of rayleigh scattering optical depth unity:
c
      ktau1 = nlyr + 1
      ptauray = 0.0
      do 4121 nze=nzup,nstr
          do 4101 k=1,nlyr
              if(ptauray .eq. 0.0 .and. 
     -           tauray(k+1)/umu(nze) .ge. 1.0) then
             ktau1 = k+1
             ptauray = p(k)
              endif 
4101      continue
          dtau_tst = (tauray(ktau1) - tauray(ktau1-1))/umu(nze)
          if(dtau_tst .ne. 0.0) then
            p_ray(nze) = 1.e-5*(p(ktau1-1) + 
     -               (1.0 - tauray(ktau1-1)/umu(nze))*
     -               (p(ktau1) - p(ktau1-1))/dtau_tst)
          else
            p_ray(nze) = 1.e-5*p(ktau1)
          endif
4121  continue
c
c****       find pressure of rayleigh scattering optical depth unity 
c           for a vertical path:
c
      ktau1 = nlyr+1
      ptauray = 0.0
      do 4141 k=1,nlyr
          if(ptauray .eq. 0.0 .and. 
     -       tauray(k+1) .ge. 1.0) then
             ktau1 = k+1
             ptauray = p(k)
          endif 
4141      continue
      dtau_tst = tauray(ktau1) - tauray(ktau1-1)
      if(dtau_tst .ne. 0.0) then
        p_ray_0 = 1.e-5*(p(ktau1-1) + 
     -           (1.0 - tauray(ktau1-1))*
     -           (p(ktau1) - p(ktau1-1))/dtau_tst)
      else
        p_ray_0 = 1.e-5*p(ktau1)
      endif
c
c****       update extinction counter and rayleigh flag
c
      ne = ne + 1
      iflext(ne) = 0
c
c*****       find the next required rayleigh scattering wavenumber
c
      nsiext(2,ne) = 2
      if(tauray(ktau1) .lt. 0.01) then
        wnext(2,ne) = 1.05d0*wn
      else
        if(tauray(ktau1) .lt. 0.1) then
          wnext(2,ne) = 1.010d0*wn
        else
          wnext(2,ne) = 1.0010d0*wn
        endif
      endif
c
      return
      end
