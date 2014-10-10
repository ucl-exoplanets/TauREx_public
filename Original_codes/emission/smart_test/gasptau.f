      subroutine gasptau(ne,iuabc,ng0,ngt0,igpab,nlev,
     -                           nsiext,nlay,istate,nt_pd,ntaug,
     -                           wn,wnmin,wnmax,wnext,p,t,alt,rmix,
     -                           d_pres,d_temp,d_rmix,
     -                           dtaug,dtaugdv,wneof,
     -                           io_end,io_err)
c
cccccccccccccccccccccccccc  g a s p t a u  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine reads gas p-induced absorption coefficients &   cc
cc    interpolates them to the appropriate pressure, and temperature  cc
cc    grids.                                                          cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc         ne - index of extinction source                            cc
cc        ng0 - gas index (1 - ngases)                                cc
cc      nt_pd - number of temperature values where absorption         cc
cc              coefficients are needed (2 are needed if temperature  cc
cc      wnmin - minimum wavenumber (cm**-1)                           cc
cc      wnmax - maximum wavenumber (cm**-1)                           cc
cc     d_pres - perturbed surface pressure used when pressure is a    cc
cc              variable part of the state vector                     cc
cc     d_temp - fractional temperature perturbation used when         cc
cc              temperature is a variable componen of the state       cc
cc              vector.                                               cc
cc     d_rmix - perturbed gas mixing ratio used when rmix is a        cc
cc              variable part of the state vector                     cc
cc      ntaug - number of optical depths needed for this gas:         cc
cc              =1 if this gas is not a variable part of state vector cc
cc              =3 if this gas is a variable part of state vector     cc
cc              (nominal, rmix perturbed at top of layer, rmix        cc
cc               perturbed at bottom of layer).                       cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     io_end - end-of file flag                                      cc
cc     io_err - flag specifying error in file                         cc
cc      wneof - wavenumber at end of file                             cc
cc      wnpab - wavenumber of each pressure-induced optical depth     cc
cc        pab - pressure-induced gas optical depth at each input wn   cc
cc       npab - number of wavenumbers were pab's are specified        cc
cc      dtaug - gas optical depth in each layer and input wavelength  cc
cc    dtaugdv - wavenumber gradient of each gas optical depth         cc
cc                                                                    cc
cccccccccccccccccccccccccc  g a s p t a u  ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ngp,nsp
      parameter (ngp=5,nsp=10000)
c
      integer ne,ng0,ngt0,nlev,nlay,nt_pd,ntaug
      integer igpab(ngtmax,ngas),istate(nex)
      integer nsipr(ngtmax,ngp),npab(ngtmax,ngp)
      integer nsiext(2,nex),iuabc(ngtmax,ngas)
      integer io_err(nex),io_end(nex)
      integer nlyr,ng,ngt,ngpab,k,l,m,n,ni,ns,ntpab,nxmax,nymax,mm
      integer iuabc0,npab0
c
      double precision wn,wnmin,wnmax,wnext(2,nex),wneof(2)
c
      real pab0(kp,nsp),tpab0(kp)
      real pab1(kp),rho(kp,3)
      real dtaupab(kp,nsp,ngtmax,ngp,2,3)
      save dtaupab
      real pab(kp),rgas,convrt,tt,tt0(1)
      double precision wnpab(nsp,ngtmax,ngp)
      save wnpab
      double precision wnpab0(nsp)
      real dvi
c
c****   atmospheric structure variables
c
      real p(kp),t(kp),alt(kp),rmix(kp,ngas)
c
c****   gas absorption coefficients
c
      real dtaug(kp,2,nex,2,3),dtaugdv(kp,nex,2,3)
c
c****   partial derivative fractional change
c
      real d_pres(2),d_temp(2),d_rmix(2,ngas)
c
      data rgas /8314.0/
c
      nlyr = nlev - 1
      ng = ng0
      ngt = ngt0
      ngpab = igpab(ngt,ng)
c
      if(nsiext(2,ne) .eq. 0) then
c
c****     no pressure induced absorption coefficents have been
c         read in yet.  define conversion from cgs to cm-amagat
c         (1 amagat = density at stp = 0.044616 kmoles/m**3)
c
        convrt = 1./0.044616
c
c****     find the density (amagats) and absorber amount (cm-amagat**2)
c
c            write(*,*) '08/18, output rmix to see if it is mmr'
c
        do 1021 l=1,nt_pd
            do 1001 k=1,nlev
                tt = d_temp(l)*t(k)
                rho(k,l) = convrt*rmix(k,ng)*p(k)/(rgas*tt)
ccccccccccccc
c            write(*,*) rmix(k,ng)
ccccccccccccc
1001        continue
c
1021    continue
c
        if(istate(1) .ne. 0) then
c
c****        pressure is a variable part of the state vector:
c            create a perturbed layer at  base of the atmosphere with
c            perturbed pressure, but the same temperature. 

          rho(nlev+1,1) = convrt*rmix(nlev,ng)*
     -                    d_pres(2)*p(nlev)/(rgas*t(nlev))
        endif
c
c****     read wavelength-dependent pressure-induced
c         absorption coefficient for each gas.
c
        iuabc0 = iuabc(ngt,ng)
c
        call readpab(iuabc0,wnmin,wnmax,npab0,ntpab,
     -               tpab0,wnpab0,pab0,io_end(ne),io_err(ne),wneof)
c
        npab(ngt,ngpab) = npab0
        if(npab(ngt,ngpab) .eq. 0) then
c
c****      there is no continuum absorption in this spectral range
c
          io_end(ne) = 1
          nsiext(1,ne) = 1
          nsiext(2,ne) = 2
          wneof(1) = -9999.0d0       
          wneof(2) = -9999.0d0       
          wnext(1,ne) = wnmin
          wnext(2,ne) = wnmax
c
c****      inialize gas optical depth - use nlyr+1 layers just in case
c          istate(1) = 1.
c
          do 1061 m=1,ntaug
              do 1041 k=1,nlay
                  dtaug(k,1,ne,1,m) = 0.0
                  dtaug(k,2,ne,1,m) = 0.0
                  dtaugdv(k,ne,1,m) =  0.0
                  dtaug(k,1,ne,2,m) = 0.0
                  dtaug(k,2,ne,2,m) = 0.0
                  dtaugdv(k,ne,2,m) =  0.0
1041          continue
1061      continue
c
          return
c
        endif
c
c****     interpolate absorption coefficients to the model atmosphere
c         temperature grid
c
        if(ntpab .gt. 1) then
c
c****      pressured induced absorption coefficients are defined at  
c          ntpab temperatures
c
          nxmax = kp
          nymax = 1
c
          do 1281 n=1,npab(ngt,ngpab)
              do 1261 l=1,nt_pd
                  do 1221 k=2,nlev  
                      do 1201 mm=1,ntpab
                          pab1(mm) = pab0(mm,n)
1201                  continue
                      tt0(1) = d_temp(l)*t(k)
c
                      call xyinterp(tpab0,pab1,tt0,pab(k),
     -                      nxmax,nymax,ntpab,nymax,nymax)
c
1221              continue
c
                  do 1241 k=1,nlyr
                      dtaupab(k,n,ngt,ngpab,l,1) = 50000.*
     -                   ((rho(k+1,l)**2)*pab(k+1) + 
     -                    (rho(k,l)**2)*pab(k))*(alt(k) - alt(k+1))
c
                      if(ntaug .gt. 1) then
c
c****                    find the (rho**2*pab)dz for mixing ratio 
c                        perturbation at top of layer k (level k)
c
                        dtaupab(k,n,ngt,ngpab,l,2) = 50000.*
     -                         ((rho(k+1,l)**2)*pab(k+1) +
     -                         ((d_rmix(2,ng)*rho(k,l))**2)*pab(k))*
     -                         (alt(k) - alt(k+1))
c
c****                     find the (rho**2)dz for mixing ratio 
c                         perturbation at bottom of layer k (level k+1)
c
                        dtaupab(k,n,ngt,ngpab,l,3) = 50000.*
     -                     (((d_rmix(2,ng)*rho(k+1,l))**2)*pab(k+1) + 
     -                       (rho(k,l)**2)*pab(k))*
     -                       (alt(k) - alt(k+1))
                      endif
c
1241              continue
c
1261          continue
c
              if(istate(1) .ne. 0) then
c
c****             pressure is a variable part of the state vector.
c                 find optical depth for perturbed surface pressure.
c                 (use pab values for the nomial surface temperature)
c
                dtaupab(nlay,n,ngt,ngpab,1,1) = 50000.*
     -                   ((rho(nlev+1,l)**2)*pab(nlev) + 
     -                    (rho(nlev-1,l)**2)*pab(nlev-1))*
     -                    (alt(nlev-1) - alt(nlev+1))
              endif
c
1281      continue
c
        else
c
          do 1641 n=1,npab(ngt,ngpab)
              do 1631 m=1,ntaug
                  do 1621 l=1,nt_pd
                      do 1601 k=1,nlyr
                          dtaupab(k,n,ngt,ngpab,l,m) = 50000.*
     -                       pab0(1,n)*(rho(k+1,l)**2 + rho(k,l)**2)*
     -                       (alt(k) - alt(k+1))
1601                  continue
c
1621              continue
c
1631          continue
c
              if(istate(1) .ne. 0) then
c
c****             pressure is a variable part of the state vector.
c                 find optical depth for perturbed surface pressure.
c                 (use pab values for the nomial surface temperature)
c
                dtaupab(nlay,n,ngt,ngpab,1,1) = 50000.*
     -                 pab0(1,n)*(rho(nlev+1,l)**2 + rho(nlev-1,l)**2)*
     -                 (alt(nlev-1) - alt(nlev+1))
                
              endif
c
1641      continue
c
        endif
c
        do 1801 n=1,npab(ngt,ngpab)
            wnpab(n,ngt,ngpab) = wnpab0(n)
1801    continue
c
        wneof(1) = wnpab(1,ngt,ngpab)
        if(npab(ngt,ngpab) .gt. 0) then
            wneof(2) = wnpab(npab(ngt,ngpab),ngt,ngpab)
        else
            wneof(2) = -99999.0d0
        endif
c
c****     initialize spectral counters
c
        nsiext(1,ne) = 1
        nsiext(2,ne) = 2
        nsipr(ngt,ngpab) = 0
c
c****   find the spectral index for the first wn
c
        if(wnpab(1,ngt,ngpab) .gt. wnmin) then
          wnext(1,ne) = 0.0d0
          wnext(2,ne) = wnmin
c
c****       initialize optical depths. Define nlyr+1 just in case
c           pressure is a variable part of the state vector.
c
          do 2011 m=1,ntaug
              do 2001 k=1,nlay
                  dtaug(k,1,ne,1,m) = 0.0
                  dtaug(k,2,ne,1,m) = 0.0
                  dtaugdv(k,ne,1,m) = 0.0
                  dtaug(k,1,ne,2,m) = 0.0
                  dtaug(k,2,ne,2,m) = 0.0
                  dtaugdv(k,ne,2,m) = 0.0
2001          continue
2011      continue
c
        else
c
c****           find the first coefficient within wnmin - wnmax
c
          do 2021 ns=1,npab(ngt,ngpab)
              if(wnpab(ns,ngt,ngpab) .gt. wnmin) go to 2041 
              nsipr(ngt,ngpab) = ns
2021      continue
c
2041      wnext(2,ne) = wnpab(nsipr(ngt,ngpab),ngt,ngpab)
c
          do 2091 m=1,ntaug
              do 2081 l=1,nt_pd
                  do 2061 k=1,nlay
                      dtaug(k,2,ne,l,m) = 
     -                      dtaupab(k,nsipr(ngt,ngpab),ngt,ngpab,l,m)
2061              continue
c
2081          continue
c
2091      continue
c
        endif
c
      endif
c
c****    increment pressure-induced absorption counter.
c
3001  nsipr(ngt,ngpab) = nsipr(ngt,ngpab) + 1
c
c****    swap spectral interval counters
c
      ni = nsiext(1,ne)
      nsiext(1,ne) = nsiext(2,ne)
      nsiext(2,ne) = ni
c
c****  find optical depth at next input wavenumber. 
c
      if(nsipr(ngt,ngpab) .le. npab(ngt,ngpab)) then
        wnext(ni,ne) = wnpab(nsipr(ngt,ngpab),ngt,ngpab)
        dvi = real(1.0d0/(wnext(2,ne) - wnext(1,ne)))
        do 3051 m=1,ntaug
           do 3041 l=1,nt_pd
               do 3021 k=1,nlay
                   dtaug(k,ni,ne,l,m) = 
     -                   dtaupab(k,nsipr(ngt,ngpab),ngt,ngpab,l,m)
                   dtaugdv(k,ne,l,m) =  (dtaug(k,2,ne,l,m) - 
     -                                   dtaug(k,1,ne,l,m))*dvi
3021           continue
3041       continue
3051    continue
        wneof(2) = wnext(nsiext(2,ne),ne)
        if(wnext(ni,ne) .lt. wn) go to 3001
c
      else
c
c****     you've run out of coefficients.  set further values to zero.
c
        io_end(ne) = 1
        io_err(ne) = 0
        wneof(2) = wnext(nsiext(1,ne),ne)        
c
        wnext(ni,ne) = wnmax
c
c****    set all optical depths to zero.  
c
        do 3071 m=1,ntaug
            do 3061 k=1,nlay
                dtaug(k,1,ne,1,m) = 0.0
                dtaug(k,2,ne,1,m) = 0.0
                dtaugdv(k,ne,1,m) =  0.0
                dtaug(k,1,ne,2,m) = 0.0
                dtaug(k,2,ne,2,m) = 0.0
                dtaugdv(k,ne,2,m) =  0.0
3061        continue
3071    continue
c
      endif
c
      return
      end
