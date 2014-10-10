      subroutine aer_tau(ne,nlyr,iuaer,nmodes,nzup,nstr,nmom,
     -                   wn,wnmin,wnmax,umu,p,dtauex,g,phmom,
     -                   tau_ext,tau_sca,g_sca,wnext,
     -                   aerqext,aerqsca,aerg0,aerpmom,
     -                   dqextdv,dqscadv,dpmomdv,dg0dv,
     -                   dtausc,dtauaer,p_aer,p_aer_0,tauaer,
     -                   iflext,nsiext,nmomaer,nmom_mx,ntau_pd0,
     -                   istate,ngases,modepd,wn_eof,io_end,io_err)
c
ccccccccccccccccccccccccccc  a e r _ t a u  cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine finds the aerosol extinction and scattering     cc
cc    cross sections at each model level, and computes the            cc
cc    monochromatic aerosol optical depth for each atmospheric layer. cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc       nlyr - number of layers in the model atmosphere.             cc
cc     nmodes - number of aerosol particle modes.                     cc
cc      iuaer - unit number from which aerosol properties are read.   cc
cc       nmom - maximum number of legendre polynomial moments         cc
cc      nt_pd - number of temperature profiles needed                 cc
cc              1 - radiances only, 2 - partial derivaties            cc
cc       nzup - index of first upward radiance stream.                cc
cc          p - pressure at each model level (Pascals).               cc
cc       nstr - number of computational zenith angles (streams).      cc
cc        umu - cosine of each upward radiance stream.                cc
cc         wn - current wavenumber (cm**-1).                          cc
cc     dtauex - differential extintion optical depth in each layer    cc
cc     dtausc - differential scattering optical depth in each layer   cc
cc      phmom - scattering phase function moments in each layer.      cc
cc    aerqext - monochromatic aerosol extinction efficiency at input  cc
cc              wavenumber,  wnext.                                   cc
cc    aerqsca - monochromatic aerosol scattering efficiency at input  cc
cc              wavenumber,  wnext.                                   cc
cc      aerg0 - scattering asymmetry parameter at input wavelength    cc
cc    aerpmom - monochromatic aerosol phase function moments          cc
cc              at input wavenumber,  wnext.                          cc
cc    dqextdv - extinction efficiency gradient at wnext.              cc
cc    dqextdv - scattering efficiency gradient at wnext.              cc
cc      dg0dv - rate of change of g with wavenumber                   cc
cc    dpmomdv - phase function moment gradient at wnext.              cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     io_end: end of file flag for each wavelength constituent:      cc
cc             (if io_end ne 0, the end of file has been encountered) cc
cc     io_err: i/o error flag  for each wavelength constituent:       cc
cc             (if io_err ne 0, an i/o error has been encountered)    cc
cc     wn_eof: last wavenumber in each input file.                    cc
cc       qext - monochromatic aerosol extinction efficiency at        cc
cc              wavenumber, wn.                                       cc
cc       qsca - monochromatic aerosol scattering efficiency at        cc
cc              wavenumber, wn.                                       cc
cc    tau_ext - column integrated extinction optical depth            cc
cc              for each variable compoent of the state vector        cc
cc    tau_sca - column integrated scattering optical depth            cc
cc              for each variable component of the state vector       cc
cc    pmomaer - monochromatic phase function moments at wn.           cc
cc     wn_eof: last wavenumber in each input file.                    cc
cc     tauaer - column-integrated aerleigh optical depth.             cc
cc     dtauex - monochromatic extinction optical depth in each layer. cc
cc     dtausc - monochromatic scattering optical depth in each layer. cc
cc      phmom - scattering phase function moments in each layer.      cc
cc      p_aer - pressure level of aerliegh tau=1 for each upward      cc
cc              stream (bars).                                        cc
cc    p_aer_0 - pressure of Rayleigh tau=1 for a vertical stream.     cc
cc       nmom - maximum number of legendre polynomial moments         cc
cc      nt_pd - number of temperature profiles needed                 cc
cc              1 - radiances only, 2 - partial derivaties            cc
cc                                                                    cc
ccccccccccccccccccccccccccc  a e r _ t a u  cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ne,nlyr,iuaer,nmodes,nzup,nstr,nmom
      integer nmomaer(2,nmode),nmom_mx(nmode)
      integer iflext(nex),nsiext(2,nex),istate(nex)
      integer ntau_pd0,ngases,modepd(mxpd)
      integer nt_pd,m,mode, iu,mom,nst,k,l,nze,ktau1
      integer io_err(nex),io_end(nex)
c
      double precision wn,wnmin,wnmax,wnext(2,nex),wn_eof(2,nex)
c
      real tauaer(kp),dtausc(kp),dtauaer(nmode,kp)
      real p(kp),p_aer(mxumu),p_aer_0
      real qext(nmode),qsca(nmode),g0(nmode),aer_dtau(kp)
      real umu(mxumu),pmomaer(0:mxmom,nmode)
      real ptauaer,dtau_tst,dv
c
c*****   monochromatic optical properties
c
      real dtauex(kp,2),g(kp),phmom(0:mxmom,kp)
      real tau_ext(kp,3,mxpd),tau_sca(kp,mxpd),g_sca(kp,mxpd)
c
c****   input aerosol wavelength-dependent variables
c
      real aerqext(2,nmode),aerqsca(2,nmode),aerg0(2,nmode),
     -     aerpmom(2,0:mxmom,nmode),dpmomdv(0:mxmom,nmode),
     -     dqextdv(nmode),dqscadv(nmode),dg0dv(nmode)
c
c****    enter loop over particle modes
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
      p_aer_0 = p(nlyr+1)
      do 1001 nze=1,nstr
          p_aer(nze) = p(nlyr+1)
1001  continue
      do 1021 k=1,nlyr
          aer_dtau(k) = 0.0
1021  continue
c
      do 2281 m=1,nmodes
c
c****       update extinction counter
c
          ne = ne + 1
c
c*****         find out whether this aerosol is a variable element of 
c               the state vector
c
          if(iflext(ne) .eq. 1) then
            mode = m
            iu = iuaer + m
c
c****        read aerosol optical properties at next wavenumber.
c
            call aerosol(ne,mode,iu,wn,wnmin,wnmax,
     -                   aerqext,aerqsca,aerg0,aerpmom,
     -                   dpmomdv,dqextdv,dqscadv,dg0dv,
     -                   nmomaer,nmom_mx,nsiext,wnext,
     -                   wn_eof,io_end,io_err)
c
c****         turn off the wavenumber grid flag.
c
            iflext(ne) = 0
c
          endif
c
c****         update the legendre polynomial mode counter
c
          if(nmom_mx(m) .gt. nmom) nmom = nmom_mx(m)
c
c
c****       find the absorption and scattering efficiencies 
c           and asymmetry factors at wn for this mode
c
          dv = real(wn - wnext(nsiext(1,ne),ne))
c
          qext(m) = aerqext(nsiext(1,ne),m) + dqextdv(m)*dv
          qsca(m) = aerqsca(nsiext(1,ne),m) + dqscadv(m)*dv
          g0(m) = aerg0(nsiext(1,ne),m) + dg0dv(m)*dv
          do 2001 mom=0,nmom_mx(m)
              pmomaer(mom,m) = aerpmom(nsiext(1,ne),mom,m) + 
     -                        dpmomdv(mom,m)*dv
2001      continue
c
c****       find aerosol differential extinction and scattering 
c           optical depths for each mode.
c
          do 2061 k=1,nlyr
              aer_dtau(k) = 0.0
              if(dtauaer(m,k) .gt. 0.) then
                aer_dtau(k) = aer_dtau(k) + qext(m)*dtauaer(m,k)
                dtausc(k) = dtausc(k) + qsca(m)*dtauaer(m,k)
                g(k) = g(k) + qsca(m)*dtauaer(m,k)*g0(m)
                do 2021 l=1,nt_pd
                    dtauex(k,l) = dtauex(k,l) + qext(m)*dtauaer(m,k)
2021            continue
                do 2041 mom=0,nmom_mx(m)
                    phmom(mom,k) = phmom(mom,k) + 
     -                             qsca(m)*dtauaer(m,k)*
     -                             pmomaer(mom,m)
2041            continue
              endif
2061      continue
c
c****      if pressure is a varible part of the state vector 
c          a perturbed layer is added to the bottom of the atmosphere 
c          with a surface pressure greater than the nominal value.
c          the aerosol optical properties in this layer are identical to
c          those in the lowest layer of the atmosphere 
c
          if(istate(1) .ne. 0) then
c
            dtausc(nlyr+1) = dtausc(nlyr+1) + qsca(m)*dtauaer(m,nlyr)
            do 4001 l=1,nt_pd
                dtauex(nlyr+1,l) = dtauex(nlyr+1,l) + 
     -                             qext(m)*dtauaer(m,nlyr)
4001        continue
c
          endif               
c
c*****         determine whether this gas is a variable element of 
c              the state vector.  If it is, update extinction counter
c              used in calculating optical depth partial derivatives
c
          nst = m+ngases+2
          if(iabs(istate(nst)) .eq. 4) then
            ntau_pd0 = ntau_pd0 + 1
            modepd(ntau_pd0) = m
c
            if(ntau_pd0 .gt. 0) then
              do 2261 k=1,nlyr
                  if(dtauaer(m,k) .gt. 0.) then
                    tau_ext(k,1,ntau_pd0) = qext(m)*dtauaer(m,k)
                    tau_sca(k,ntau_pd0) = qsca(m)*dtauaer(m,k)
                    g_sca(k,ntau_pd0) =  g0(m)
                  endif
2261          continue
            endif
          endif
2281  continue
c
c****   find the column-integrated aerosol optical depth
c
      tauaer(1) = 0.0
      do 3001 k=2,nlyr+1
          tauaer(k) = tauaer(k-1) + aer_dtau(k-1)
3001  continue
c
c****   find the layer of aerosol optical absorption optical depth
c        unity for each output stream
c
      do 5021 nze=nzup,nstr
          ptauaer = 0.0
          ktau1 = nlyr+1
c
          do 5001 k=1,nlyr
              if(ptauaer .eq. 0.0 .and. 
     -           tauaer(k+1)/umu(nze) .ge. 1.0) then
                 ktau1 = k+1
                 ptauaer = p(k)
              endif 
5001      continue
c
c****       find the pressure of aerosol extinction optical 
c           depth unity at each output emission angle
c
          dtau_tst = (tauaer(ktau1) - tauaer(ktau1-1))/umu(nze)
          if(dtau_tst .ne. 0.0) then
            p_aer(nze) = 1.e-5*(p(ktau1-1) + 
     -                   (1.0 - tauaer(ktau1-1)/umu(nze))*
     -                   (p(ktau1) - p(ktau1-1))/dtau_tst)
          else
            p_aer(nze) = 1.e-5*p(ktau1)
          endif
5021  continue
c
c****    find the layer of gas optical absorption optical depth
c        unity for normal incident radiation
c
      ptauaer = 0.0
      ktau1 = nlyr+1
c
      do 5041 k=1,nlyr
          if(ptauaer .eq. 0.0 .and. tauaer(k+1) .ge. 1.0) then
             ktau1 = k+1
             ptauaer = p(k)
           endif 
5041  continue
c
c****  find the pressure of aerosol extinction optical 
c      depth unity at each output emission angle
c
      dtau_tst = tauaer(ktau1) - tauaer(ktau1-1)
      if(dtau_tst .ne. 0.0) then
        p_aer_0 = 1.e-5*(p(ktau1-1) + 
     -            (1.0 - tauaer(ktau1-1))*
     -            (p(ktau1) - p(ktau1-1))/dtau_tst)
      else
        p_aer_0 = 1.e-5*p(ktau1)
      endif
c
      return
      end
