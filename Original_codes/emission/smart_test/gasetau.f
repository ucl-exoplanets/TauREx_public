      subroutine gasetau(ne,iuabc,ng0,ngt0,nlev,nlay,
     -                   nsiext,nt_pd,istate,ntaug,
     -                   wn,wnmin,wnmax,wnext,
     -                   p,t,alt,rmix,ratm,
     -                   d_pres,d_temp,d_rmix,
     -                   dtaug,dtaugdv,wneof,io_end,io_err)
c
cccccccccccccccccccccccccc  g a s e t a u  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine reads wavelength-dependent gas absorption       cc
cc    cross-sections for electronic transitions, computes the gas     cc
cc    optical depth, and interpolates these values to the             cc
cc    appropriate pressure, and temperature grids.                    cc
cc                                                                    cc
cc    note: this version does not include temperature dependence.     cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc         ne - index of extinction source                            cc
cc        ng0 - gas index (1 - ngases)                                cc
cc      nt_pd - number of temperature values where absorption         cc
cc              coefficients are needed (2 are needed if temperature  cc
cc     istate - state vector flag indicating which state variables    cc
cc              are variable components of the state vector.          cc
cc       nlay - number of layers were optical depths are needed:      cc
cc              nlev-1 if pressure is not a varible part of the state cc
cc              vector, nlev if it is (istate(1) = 1)                 cc
cc      nt_pd - number of temperature values where absorption         cc
cc              coefficients are needed (2 are needed if temperature  cc
c               is a variable component of the state vector.          cc 
cc      wnmin - minimum wavenumber (cm**-1)                           cc
cc      wnmax - maximum wavenumber (cm**-1)                           cc
cc     d_pres - perturbed surface pressure used when pressure is a    cc
cc              variable part of the state vector                     cc
cc     d_temp - fractional temperature perturbation used when         cc
cc              temperature is a variable componen of the state       cc
cc              vector.                                               cc
cc     d_rmix - perturbed gas mixing ratio used when rmix is a        cc
cc              variable part of the state vector                     cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     io_end - end-of file flag                                      cc
cc     io_err - flag specifying error in file                         cc
cc      wneof - wavenumber at end of file                             cc
cc     wnxsec - wavenumber of each absorption cross section           cc
cc       xsec - cross-section at each input wn                        cc
cc      nxsec - number of wavenumbers were xsec's are specified       cc
cc      dtaug - gas optical depth in each layer and input wavelength  cc
cc    dtaugdv - wavenumber gradient of each gas optical depth         cc
cc                                                                    cc
cccccccccccccccccccccccccc  g a s e t a u  ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nsp
      parameter (nsp=14000)
c
      integer nsiext(2,nex),iuabc(ngtmax,ngas),istate(nex)
      integer nsix(ngtmax,ngas),nxsec(ngtmax,ngas),nt_pd
      save nsix,nxsec
      integer io_err(nex),io_end(nex)
      integer ne,ng0,ngt0,nlev,nlay,ntaug
      integer ng,ngt,nxmx,iquit,n,l,k,ns,ni,m
c
      double precision wn,wnmin,wnmax,wnext(2,nex),wneof(2)
c
      double precision wlxsec(nsp)
      real xsec0(nsp),atoms(kp,2)
c
c****   atmospheric structure variables
c
      real p(kp),t(kp),alt(kp),rmix(kp,ngas),ratm
      real dtaug(kp,2,nex,2,3),dtaugdv(kp,nex,2,3)
      real a0,wgtgs,convert,dvi
c
c****    uv gas absorption cross-sections
c
      double precision wnxsec(nsp,ngtmax,ngas)
      save wnxsec
      real xsec(nsp,ngtmax,ngas),rhodz(kp,ngas,2,3)
      save xsec,rhodz
c
c****   partial derivative fractional change
c
      real d_pres(2),d_temp(2),d_rmix(2,ngas)
c
c****    specify avagodro's number (kg/kmole). To convert from
c        density to number density - multiply by a0/wgt
c
      data a0 /6.02e26/
c
      ng = ng0
      ngt = ngt0
      nxmx = nsp
      iquit = 0
c
      if(nsiext(2,ne) .eq. 0) then
c
c****     read wavelength-dependent gas absorption cross-sections
c
        call readxsec(iuabc(ngt,ng),nxsec(ngt,ng),wgtgs,
     -                wnmin,wnmax,wlxsec,xsec0,nxmx,iquit)
c
c****     check ordering of absorption coefficients.  Reverse order
c         if necessary
c
        if(wlxsec(1) .lt. wlxsec(nxsec(ngt,ng))) then
          do 1001 n=1,nxsec(ngt,ng)
              if(wlxsec(nxsec(ngt,ng)-n+1) .gt. 0.0d0) then
                wnxsec(n,ngt,ng) = 1.0d4/wlxsec(nxsec(ngt,ng)-n+1)
              else
                wnxsec(n,ngt,ng) = 1.0d8
              endif
              xsec(n,ngt,ng) = xsec0(nxsec(ngt,ng)-n+1)
1001      continue
c
        else
c
          do 1021 n=1,nxsec(ngt,ng)
              if(wlxsec(n) .gt. 0.0d0) then
                wnxsec(n,ngt,ng) = 1.0d4/wlxsec(n)
              else
                wnxsec(n,ngt,ng) = 1.0d8
              endif
              xsec(n,ngt,ng) = xsec0(n)
1021      continue
        endif
c
c****     set start-of-buffer wavenumber
c
        if(wneof(1) .le. 0.0d0 .or. wnxsec(n,ngt,ng) .le. wneof(1))
     -  wneof(1) = wnxsec(n,ngt,ng)
c
c****      define factor to convert  rho(kg/m**3)*dz(km) to
c          N(atoms/cm**3)*dz(cm) = (a0/wgt)*(1e-6m**3/cm**3)*1e5cm/km
c          or 0.1*a0/wgt - note wgt divides out of rho and convert.
c
        convert = 0.1*a0/(ratm*wgtgs)
        do 2221 l=1,nt_pd
            do 2001 k=1,nlev
                atoms(k,1) = convert*rmix(k,ng)*
     -                       p(k)/(d_temp(l)*t(k))
                if(ntaug .gt. 1) then
                  atoms(k,2) = convert*d_rmix(2,ng)*rmix(k,ng)*
     -                         p(k)/(d_temp(l)*t(k))
                endif
2001        continue
c
            if(istate(1) .ne. 0) then
c
c****         if pressure is a variable part of the state vector 
c             create a perturbed layer at the base of the atmosphere 
c             with 1% more pressure than the surface layer, but the  
c             same temperature. (needed only for l=1)
c
              k = nlev
              atoms(k+1,1) = convert*d_pres(2)*
     -                        rmix(k,ng)*p(k)/(d_temp(l)*t(k))
            endif
c
c****         find rho*dz in each layer
c
            do 2201 k=1,nlev
                rhodz(k,ng,l,1) = 0.5*(atoms(k+1,1) + atoms(k,1))*
     -                              (alt(k) - alt(k+1))
                if(ntaug .gt. 1) then
c
c****                  find the sigma*rho*dz for mixing ratio 
c                      perturbation at top of layer k (level k)
c
                  rhodz(k,ng,l,2) = 0.5*(atoms(k+1,1) + atoms(k,2))*
     -                                (alt(k) - alt(k+1))
c
c****                  find the sigma*rho*dz for mixing ratio 
c                      perturbation at bottom of layer k (level k+1)
c
                  rhodz(k,ng,l,3) = 0.5*(atoms(k+1,2) + atoms(k,1))*
     -                                (alt(k) - alt(k+1))
                endif
2201        continue
c
            if(istate(1) .ne. 0) then
c
c****           find sigma*rho*dz for lower layer with peturbed pressure
c
              rhodz(nlev+1,ng,l,1) = 0.5*(atoms(nlev+1,1) + 
     -                                   atoms(nlev-1,1))*
     -                                   (alt(nlev-1) - alt(nlev+1))
              if(ntaug .gt. 1) then
                rhodz(nlev+1,ng,l,2) = rhodz(nlev+1,ng,l,1)
                rhodz(nlev+1,ng,l,3) = rhodz(nlev+1,ng,l,1)
              endif
c
            endif
c
2221    continue
c
c****     initialize spectral counters
c
        nsiext(1,ne) = 1
        nsiext(2,ne) = 2
        nsix(ngt,ng) = 0
c  
c****      find the spectral index for the first wn
c
        if(wnxsec(1,ngt,ng) .gt. wnmin) then
          wnext(1,ne) = 0.0d0
          wnext(2,ne) = wnmin
c
c****      inialize gas optical depth - use nlev layers just in case
c          istate(1) = +/-1.
c
          do 2421 m=1,ntaug
              do 2401 k=1,nlay
                  dtaug(k,1,ne,1,m) = 0.0
                  dtaug(k,2,ne,1,m) = 0.0
                  dtaugdv(k,ne,1,m) = 0.0
                  dtaug(k,1,ne,2,m) = 0.0
                  dtaug(k,2,ne,2,m) = 0.0
                  dtaugdv(k,ne,2,m) = 0.0
2401          continue
2421      continue
c
        else
c
c****           find the first coefficient within wnmin - wnmax
c
          do 2601 ns=1,nxsec(ngt,ng)
              if(wnxsec(ns,ngt,ng) .gt. wnmin) go to 2801 
              nsix(ngt,ng) = ns
2601      continue
c
2801      wnext(2,ne) = wnxsec(nsix(ngt,ng),ngt,ng)
          do 2861 l=1,nt_pd
              do 2841 m=1,ntaug
                  do 2821 k=1,nlay
	write(96,*)wn,k,ng,l,m,k,rhodz(k,ng,l,m)
                      dtaug(k,2,ne,l,m) = xsec(nsix(ngt,ng),ngt,ng)*
     -                                    rhodz(k,ng,l,m)
2821              continue
2841          continue
2861      continue
        endif
c
      endif
c
c****    increment continuum absorption counter.
c
3001  nsix(ngt,ng) = nsix(ngt,ng) + 1
c
c****    swap spectral interval counters
c
      ni = nsiext(1,ne)
      nsiext(1,ne) = nsiext(2,ne)
      nsiext(2,ne) = ni
c
c****  find optical depth at next input wavenumber
c
      if(nsix(ngt,ng) .le. nxsec(ngt,ng)) then
        wnext(ni,ne) = wnxsec(nsix(ngt,ng),ngt,ng)
        dvi = real(1.0d0/(wnext(2,ne) - wnext(1,ne)))
        do 3041 l=1,nt_pd
c
c****         find an additional value if pressure is a variable
c             part of the state vector
c
            do 3031 m=1,ntaug
                do 3021 k=1,nlay
                    dtaug(k,ni,ne,l,m) = xsec(nsix(ngt,ng),ngt,ng)*
     -                                   rhodz(k,ng,l,m)
                    dtaugdv(k,ne,l,m) =  (dtaug(k,2,ne,l,m) - 
     -                                    dtaug(k,1,ne,l,m))*dvi
3021            continue
3031        continue
3041    continue
        wneof(2) = wnext(nsiext(2,ne),ne)
        if(wnext(ni,ne) .lt. wn) go to 3001
c
      else
c
c****     you've run out of coefficients.  punt!
c
        io_end(ne) = 1
        io_err(ne) = 0
        wneof(2) = wnext(nsiext(1,ne),ne)
c
        wnext(ni,ne) = wnmax
c
        do 3061 k=1,nlay
            dtaug(k,ni,ne,1,m) = 0.0
            dtaugdv(k,ne,1,m) =  0.0
            dtaug(k,ni,ne,2,m) = 0.0
            dtaugdv(k,ne,2,m) =  0.0
3061    continue
      endif
c
      return
      end
