      subroutine gas_tau(iuabc,igtrn,igpab,ngases,ngastyp,
     -                  npg,ntg,nbgas,ne,nlev,nzup,nt_pd,ntaug,
     -                  nstr,iflext,nsiext,nlay,ntau_pd0,istate,
     -                  wn,wnmin,wnmax,ratm,umu,p,t,alt,grav,
     -                  z,rmix,pd_frac,wnext,dtauex,tau_ext,    
     -                  taugas,p_gas,p_gas_0,wn_eof,io_end,io_err)
c
ccccccccccccccccccccccccccc  g a s _  t a u cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine finds the gas absorption coefficients at each   cc
cc    model level, and computes the effective monochromatic gas       cc
cc    absorption optical depth for each atmospheric layer.            cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc      iuabc - unit numbers for gas line parameters                  cc
cc       nlyr - number of layers in the model atmosphere (nlev-1).    cc
cc       nlev - number of levels in the model atmosphere (<kp).       cc
cc     ngases - number of absorbing gases included in model.          cc
cc         ne - extinction/source index.                              cc
cc     iflext - extinction flag: 0) do not read next spectral value,  cc
cc             1) read next spectral value.                           cc
cc      igtrn - type of energy transition: (1) line,                  cc
cc             (2) pressure-induced, (3) electronic                   cc
cc       ratm - mean atmospheric gas constant (j/kg/k).               cc
cc       nzup - index of first upward radiance stream.                cc
cc       nstr - number of computational zenith angles (streams).      cc
cc        umu - cosine of each upward radiance stream.                cc
cc         wn - current wavenumber (cm**-1).                          cc
cc      ntaug - number of optical depths needed for this gas:         cc
cc              =1 if this gas is not a variable part of state vector cc
cc              =3 if this gas is a variable part of state vector     cc
cc              (nominal, rmix perturbed at top of layer, rmix        cc
cc               perturbed at bottom of layer).
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     io_end: end of file flag for each wavelength constituent:      cc
cc             (if io_end ne 0, the end of file has been encountered) cc
cc     io_err: i/o error flag  for each wavelength constituent:       cc
cc             (if io_err ne 0, an i/o error has been encountered)    cc
cc     wn_eof: last wavenumber in each input file.                    cc
cc     taugas - column-integrated gas optical depth.                  cc
cc      dtaug - monochromatic gas optical depth in each layer.        cc
cc     dtauex - monochromatic extinction optical depth in each layer. cc
cc    dtaugdv - derivative of dtaug with respect to wavenumber.       cc
cc      p_gas - pressure level of gas tau=1 for each upward stream    cc
cc              (bars).                                               cc
cc    p_gas_0 - pressure of gas tau=1 for a vertical stream.          cc
cc     d_pres - perturbed surface pressure used when pressure is a    cc
cc              variable part of the state vector                     cc
cc     d_temp - fractional temperature perturbation used when         cc
cc              temperature is a variable componen of the state       cc
cc              vector.                                               cc
cc     d_rmix - perturbed gas mixing ratio used when rmix is a        cc
cc              variable part of the state vector                     cc
cc                                                                    cc
ccccccccccccccccccccccccccc  g a s _  t a u cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer iuabc(ngtmax,ngas),igtrn(ngtmax,ngas),igpab(ngtmax,ngas)
      integer ngases,ngastyp(ngas)
      integer npg(ngtmax,ngas),ntg(ngtmax,ngas),
     -          nbgas(ngtmax,ngas)
      integer ne,nlev,nzup,nstr
      integer iflext(nex),nsiext(2,nex)
      integer io_err(nex),io_end(nex)
      integer ntau_pd0,istate(nex)
      integer ntaug
c
      integer nlyr,nlay,nt_pd,k,ng,ng0,ngt,ngt0,l,nze,ktau1,ni,nee,m
c
      double precision wn,wnmin,wnmax
c
      real p_gas(mxumu),p_gas_0,ratm,umu(mxumu)
c
      double precision wnext(2,nex),wn_eof(2,nex)
c
      real p(kp),t(kp),alt(kp),grav(kp),z(kp),rmix(kp,ngas)
      real pd_frac(nex)
      real dtauex(kp,2),tau_ext(kp,3,mxpd),taugas(kp)
c
      real dtaug(kp,2,nex,2,3),dtaugdv(kp,nex,2,3)
      save dtaug,dtaugdv
c
      double precision wneof(2)
c
      real d_pres(2),d_temp(2),d_rmix(2,ngas)
      real dv,ptaugas,dtau_tst,dtaugas(kp,2,3),dtau0
c
      nlyr = nlev - 1
c
c*****  set the wavenumber start counter
c
      wneof(1) = wn_eof(1,ne)
c
c*****   determine whether pressure is a variable part of
c        the state vector.  if it is, gas absorption must be 
c        computed for the background and perterbed pressure profile
c
      d_pres(1) = 1.0
      if(istate(1) .ne. 0) then
         d_pres(2) = 1.0 + pd_frac(1)
      endif
c
c*****   determine whether temperature is a variable part of
c        the state vector.  if it is, gas absorption must be 
c        computed for the background and perterbed temperature profile
c
      d_temp(1) = 1.0
      if(istate(2) .ne. 0) then
         d_temp(2) = 1.0 + pd_frac(2)
      endif
c
      do 1001 k=1,nlev
          taugas(k) = 0.0
1001  continue
      if(istate(1) .ne. 0) taugas(nlev+1) = 0.0
c
c Danie Liang, 10/05/2007
c      if(wn .eq. wnmin) then
      if(wn .lt. wnmin) then
        do 1191 m=1,3
            do 1181 ni=1,2
                do 1161 nee=1,nex
                    do 1101 k=1,nlev
                         dtaugdv(k,nee,ni,m) = 0.0
1101                continue
                    do 1141 l=1,2
                        do 1121 k=1,nlev
                            dtaug(k,l,nee,ni,m) = 0.0 
1121                    continue
1141                continue
1161            continue
1181        continue
1191    continue
c
      endif
c
c****   enter the loop over absorbing gases
c
      do 1801 ng=1,ngases
          ng0 = ng
c
c*****         determine whether this gas is a variable element of 
c              the state vector.  If it is, update extinction counter
c              used in calculating partial derivatives
c
          ntaug = 1
          d_rmix(1,ng) = 1.0
          if(iabs(istate(ng+2)) .eq. 3) then
            ntau_pd0 = ntau_pd0 + 1
            ntaug = 3
            d_rmix(2,ng) = (1.0 + pd_frac(ng+2))
c
            do 1221 m=1,3
                do 1201 k=1,nlay
                    tau_ext(k,m,ntau_pd0) = 0.0
1201            continue
1221        continue
c
          endif
c
c****        update extinction index and check extinction flag
c
          do 1701 ngt=1,ngastyp(ng)
              ngt0 = ngt
              ne = ne + 1
              if(iflext(ne) .eq. 1) then
c
c****              determine the type of transition: (1) line, 
c                  (2) pressure-induced, (3) electronic
c
                if(igtrn(ngt,ng) .eq. 1) then
c
c****              read the next line absorption coefficient
c
                  call gaslbltau(ne,iuabc,ng0,ngt0,nlev,npg,ntg,nbgas,
     -                     nsiext,nlay,istate,nt_pd,ntaug,
     -                     wn,wnmin,wnmax,wnext,p,t,z,grav,rmix,
     -                     d_pres,d_temp,d_rmix,
     -                     dtaug,dtaugdv,wneof,
     -                     io_end,io_err)
c
                else
c
                  if(igtrn(ngt,ng) .eq. 2) then
c
c****                 read next p-induced transition
c
                    call gasptau(ne,iuabc,ng0,ngt0,igpab,nlev,
     -                           nsiext,nlay,istate,nt_pd,ntaug,
     -                           wn,wnmin,wnmax,wnext,p,t,alt,rmix,
     -                           d_pres,d_temp,d_rmix,
     -                           dtaug,dtaugdv,wneof,
     -                           io_end,io_err)
c
                  else
c
c****                 read next electronic transition cross-section
c
                    call gasetau(ne,iuabc,ng0,ngt0,nlev,nlay,
     -                           nsiext,nt_pd,istate,ntaug,
     -                           wn,wnmin,wnmax,wnext,
     -                           p,t,alt,rmix,ratm,
     -                           d_pres,d_temp,d_rmix,
     -                           dtaug,dtaugdv,wneof,io_end,io_err)
c
                  endif
                endif
c
c****             turn off the extinction flag for this gas
c
              iflext(ne) = 0
c
              endif
              wn_eof(1,ne) = wneof(1)
              wn_eof(2,ne) = wneof(2)
c
c****          interpolate gas optical depths to wavenumber wn
c
              dv = real(wn - wnext(nsiext(1,ne),ne))
c
              do 1441 m=1,ntaug
                  do 1421 l=1,nt_pd
                       do 1401 k=1,nlay
                           dtau0 = dtaug(k,nsiext(1,ne),ne,l,m) + 
     -                               dtaugdv(k,ne,l,m)*dv
                           if(dtau0 .gt. 0.0) then
                             dtaugas(k,l,m) = dtau0
                           else
                             dtaugas(k,l,m) = 0.0
                           endif
1401                   continue
1421              continue
1441          continue
              do 1481 l=1,nt_pd
                  do 1461 k=1,nlay
                      dtauex(k,l) = dtauex(k,l) + dtaugas(k,l,1)
1461              continue
1481          continue
c
              do 1501 k=1,nlay
                  taugas(k+1) = taugas(k+1) + dtaugas(k,1,1)
1501          continue
c
c****                if this is a variable component of the state 
c                    vector, save differential extinction optical depth 
c
              if(iabs(istate(ng+2)) .eq. 3) then
c
                do 1621 m=1,ntaug
                    do 1601 k=1,nlay
                        tau_ext(k,m,ntau_pd0) = 
     -                           tau_ext(k,m,ntau_pd0) + dtaugas(k,1,m)
1601                continue
c
1621            continue
c
              endif
c
1701      continue
c
1801  continue
c
c****      find the layer of gas optical absorption optical depth 
c          unity at each of the output emission angles
c
      do 2041 nze=nzup,nstr
          ptaugas = 0.0
          ktau1 = 0
c
          do 2001 k=1,nlyr
              if(ptaugas .eq. 0.0 .and. 
     -           taugas(k+1)/umu(nze) .ge. 1.0) then
                ktau1 = k+1
                ptaugas = p(k+1)
              endif 
2001      continue
c
c****       find the pressure of gas absorption optical depth unity:
c
          if(ktau1 .gt. 1) then
c
            dtau_tst = (taugas(ktau1) - taugas(ktau1-1))/umu(nze)
            if(dtau_tst .ne. 0.0) then
              p_gas(nze) = 1.e-5*(p(ktau1) - 
     -               (taugas(ktau1)/umu(nze) - 1.0)*
     -               (p(ktau1) - p(ktau1-1))/dtau_tst)

            else
              p_gas(nze) =  1.e-5*p(ktau1)
            endif
c
          else
c
c****        tau < 1 in column.  set to surface value
c
            ktau1 = nlev
            p_gas(nze) = 1.e-5*p(ktau1)
c
          endif
2041  continue
c
c****       find the pressure of tau = 1 for a vertical path
c
      ptaugas = 0.0
      ktau1 = 0
c
      do 3001 k=1,nlyr
          if(ptaugas .eq. 0.0 .and. 
     -       taugas(k+1) .ge. 1.0) then
             ktau1 = k+1
             ptaugas = p(k+1)
          endif 
3001  continue
c
c****      find the pressure of gas absorption optical depth unity:
c          note: ktau can't equal to 1.0 because taugas(1) = 0.0
c
      if(ktau1 .gt. 1) then
c
        dtau_tst = taugas(ktau1) - taugas(ktau1-1)
        if(dtau_tst .ne. 0.0) then
          p_gas_0 = 1.e-5*(p(ktau1) - (taugas(ktau1) - 1.0)*
     -              (p(ktau1) - p(ktau1-1))/dtau_tst)
        else
          ktau1 = nlev
          p_gas_0 = 1.e-5*p(ktau1)
        endif
c
      else
c
c****      tau < 1 in column.  set to surface value
c
        p_gas_0 = 1.e-5*p(nlev)
c
      endif
c
c Danie Liang
      do k=1,nlay
      write(97,*)wn,p(k),taugas(k+1)
      enddo       

      return
      end
