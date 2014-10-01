      subroutine gaslbltau(ne,iuabc,ng0,ngt0,nlev,npg,ntg,nbgas,
     -                     nsiext,nlay,istate,nt_pd,ntaug,
     -                     wn,wnmin,wnmax,wnext,p,t,z,grav,rmix,
     -                     d_pres,d_temp,d_rmix,
     -                     dtaug,dtaugdv,wneof,
     -                     io_end,io_err)
c
cccccccccccccccccccccccc  g a s l b l t a u  ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine reads gas absorption coefficients, inter-       cc
cc    polates them to the appropriate pressure and temperature        cc
cc    grids, and computes the gas optical depth, and it derivative    cc
cc    with respect to wavenumber                                      cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc    ntemp - Maximum number of temperature profiles for gas          cc
cc            absorption coefficients.				      cc
cc         ne - index of extinction source                            cc
cc         ng - gas index (1 - ngases)                                cc
cc        ngt - number of different absorption coefficeint types      cc
cc              for this gas.                                         cc
cc       nlev - number of levels in the model atmosphere (<kp).       cc
cc       nlay - number of layers were optical depths are needed:      cc
cc              nlev-1 if pressure is not a varible part of the state cc
cc              vector, nlev if it is (istate(1) = 1)                 cc
cc      nt_pd - number of temperature values where absorption         cc
cc              coefficients are needed (2 are needed if temperature  cc
c               is a variable component of the state vector.          cc
cc      ntaug - number of optical depth profile of each gas:          cc
cc                - if istate not equal to 3, ntaug = 1               cc
cc                - if istate equal to 3, ntaug = 3                   cc
cc      wnmin - minimum wavenumber (cm**-1)                           cc
cc      wnmax - maximum wavenumber (cm**-1)                           cc
cc         p0 - pressures of input absorption coefficients            cc
cc         t0 - temperatures of input absorption coefficients         cc
cc      rmix0 - mixing ratios of input absorption coefficients        cc
cc         wn - wavenumber of each gas absorption coefficient         cc
cc        abs - input absorption coefficient                          cc
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
cc      dtaug - gas optical depth in each layer and input wavelength  cc
cc    dtaugdv - wavenumber gradient of each gas optical depth         cc
cc                                                                    cc
cccccccccccccccccccccccc  g a s l b l t a u  ccccccccccccccccccccccccccc
c
      implicit none
c
      integer ntemp
      parameter (ntemp=3)
c
      include 'param.inc'
c
      integer npg(ngtmax,ngas),ntg(ngtmax,ngas),
     -           nbgas(ngtmax,ngas),nt_pd,ntaug,istate(nex)
      integer nsiext(2,nex),iuabc(ngtmax,ngas)
      integer ne,ng0,ngt0
      integer io_err(nex),io_end(nex)
      integer nlev,nlay
      integer nlyr,k,i,l,nxmax,nymax,ni,n1,n2,itl,m,ng,ngt
c
      integer int_p(ngtmax,ngas),int_t(ngtmax,ngas),
     -                 itlev(kp,ngtmax,ngas)
      save int_p,int_t,itlev
c
      real ppas(kp,ngtmax,ngas),tatm1(kp,ntemp,ngtmax,ngas)
      save ppas,tatm1
c
      double precision wn,wnmin,wnmax,wnext(2,nex),wneof(2)
c
      double precision dv
      real dvi
      real tol,p5,dt0,wns,dnu0,dt02,stp,stm,tt,delt
      real abc0(kp,3,2),abc1(kp,3),abc2(kp,2,3),z1(kp),
     -        p0(kp),tatm(kp,ntemp),rg(kp,ngas),at(kp),bt(kp)
c
c****    atmospheric structure and gas mixing ratios 
c
      real p(kp),t(kp),z(kp),grav(kp),rmix(kp,ngas)
      real dtaug(kp,2,nex,2,3),dtaugdv(kp,nex,2,3)
      real sigma
      real pstar(1)
c
c****   partial derivative fractional change
c
      real d_pres(2),d_temp(2),d_rmix(2,ngas)
c
      tol = 1.e-7
      p5 = 0.5
      nlyr = nlev - 1
      ng = ng0
      ngt = ngt0
c
      if(nsiext(2,ne) .eq. 0) then
c
c****    initialize gas spectral interval counters
c
        nsiext(1,ne) = 1
        nsiext(2,ne) = 2
        wnext(1,ne) = wnmin
        wnext(2,ne) = wnmin
c
c****   open read the atmospheric structure information
c
        read(iuabc(ngt,ng)) (p0(k),k=1,npg(ngt,ng))
        read(iuabc(ngt,ng)) 
     -      ((rg(k,i),i=1,nbgas(ngt,ng)),k=1,npg(ngt,ng))
        read(iuabc(ngt,ng)) 
     -      ((tatm(k,l),k=1,npg(ngt,ng)),l=1,ntg(ngt,ng))
c
c*****    determine if absorption coeficients must be interpolated
c         to a new pressure grid:
c
        if(npg(ngt,ng) .ne. nlev) then
c
c****      the input atmosphere has a different number of levels.
c          a pressure interpolation is required.
c
          int_p(ngt,ng) = 1
        else
c
c****     the input atmosphere has the same number of levels as the
c         absorption coefficient file.  Check each level
c
          int_p(ngt,ng) = 0
          do 1201 k=1,nlev
              if(abs(p0(k) - p(k)) .gt. 0.01*p(k)) int_p(ngt,ng) = 1
1201      continue      
        endif
c
        if(int_p(ngt,ng) .ne. 0) then
c
c****      interpolate temperatures to the output p grid.
c
          do 2421 k=1,npg(ngt,ng)
              ppas(k,ngt,ng) = p0(k)
              if(ppas(k,ngt,ng) .ne. 0.) then
                z1(k) = alog(ppas(k,ngt,ng))
              else
                z1(k) = - 80.
              endif
2421      continue
c
c****       interpolate absorption coefficients to output p grid.
c
          nxmax = kp
          nymax = ntemp
c
          call xyinterp(z1,tatm(1,1),z,tatm1(1,1,ngt,ng),
     -                  nxmax,nymax,npg(ngt,ng),nlev,ntg(ngt,ng))
        else
          do 2461 l=1,ntg(ngt,ng)
              do 2441 k=1,nlev
                  tatm1(k,l,ngt,ng) = tatm(k,l)
2441          continue
2461      continue
        endif
c
c****    determine if a temperature interpolation is needed.
c
        if(ntg(ngt,ng) .gt. 1) then
          int_t(ngt,ng) = 0
          do 2501 k=1,nlev
c
c****            find the index of the temperature that is closest
c
              dt0 = tatm1(k,2,ngt,ng) - tatm1(k,1,ngt,ng)
              itlev(k,ngt,ng) = nint((t(k) - tatm1(k,1,ngt,ng))/dt0) + 1
c
c****           check to see if itlev is with range of tatm
c
              if(itlev(k,ngt,ng) .lt. 1) then
                itlev(k,ngt,ng) = 1
              else
                if(itlev(k,ngt,ng) .gt. ntg(ngt,ng)) 
     -             itlev(k,ngt,ng) = ntg(ngt,ng)
              endif
c
c****           determine if a temperature interpolation is needed
c
              if(abs(tatm1(k,itlev(k,ngt,ng),ngt,ng) - t(k)) .gt. 
     -           0.00005*t(k) .or. istate(2) .ne. 0) int_t(ngt,ng) = 1
2501      continue
        else
          int_t(ngt,ng) = 0
          do 2521 k=1,nlev
              itlev(k,ngt,ng) = 1
2521      continue
        endif
c
c****   initialize absorption coefficients for this gas
c
        ni = nsiext(2,ne)
        wnext(2,ne) = wnmin
        do 2621 m=1,ntaug
            do 2601 k=1,nlay
                dtaug(k,1,ne,1,m) = 0.0
                dtaug(k,2,ne,1,m) = 0.0
                dtaug(k,1,ne,2,m) = 0.0
                dtaug(k,2,ne,2,m) = 0.0
2601        continue
2621    continue
c
      endif
c
c****    swap spectral interval counters
c
2801  ni = nsiext(1,ne)
      nsiext(1,ne) = nsiext(2,ne)
      nsiext(2,ne) = ni
c
c****     read the next absorption coefficient for this gas
c
      read(iuabc(ngt,ng),err=4201,end=4001) wns,dnu0,
     -    ((abc0(k,l,ni),k=1,npg(ngt,ng)),l=1,ntg(ngt,ng))
c
c****   determine if this segment is in desired spectral window
c
      wnext(ni,ne) = wns
      if(wnext(ni,ne) .le. wnmin) wneof(1) = wnext(ni,ne)
      if(wnext(ni,ne) .lt. wn) go to 2801
c
      wneof(2) = wnext(ni,ne)
c
c****    interpolate absorption coefficients to appropriate output grid.
c
      if(wnext(nsiext(1,ne),ne) .le. wnmin) then
        n1 = 1
      else
        n1 = 2
      endif
      n2 = 2
c
      nxmax = kp
      nymax = 3
c
      do 3481 ni=n1,n2
          if(int_p(ngt,ng) .eq. 1) then
c
c****    interpolate absorption coefficients to p grid
c
            call xyinterp(ppas(1,ngt,ng),abc0(1,1,nsiext(ni,ne)),
     -                    p,abc1,
     -                    nxmax,nymax,npg(ngt,ng),nlev,ntg(ngt,ng))
c
          else
c
c****        no p interpolation needed.  Just load abc1 array
c
            do 3121 l=1,ntg(ngt,ng)
                do 3101 k=1,nlay+1
                    abc1(k,l) = abc0(k,l,nsiext(ni,ne))
3101            continue
3121        continue
c
          endif
c
          if(istate(1) .ne. 0) then
c
c****         interpolate abc to layer (1.0+p_pert(1))*p(nlev)
c
            pstar(1) = d_pres(2)*p(nlev)
            call xyinterp(ppas(1,ngt,ng),abc0(1,1,nsiext(ni,ne)),
     -                    pstar,abc1(nlev+1,1),
     -                    nxmax,nymax,npg(ngt,ng),1,ntg(ngt,ng))  
c
c****         define tatm1 at nlev+1
c
            do 3141 l=1,ntg(ngt,ng)
                tatm1(nlev+1,l,ngt,ng) = tatm1(nlev,l,ngt,ng)
3141        continue
c
          endif
c
c****       interpolate absorption coefficients to output T grid.
c           (assume that log(abc) is quadratic in t.)
c           Weight results by mass-mixing-ratio/gravity (cgs units)
c
          if(int_t(ngt,ng) .ge. 1) then
            do 3241 k=1,nlay+1
                dt0 = tatm1(k,2,ngt,ng) - tatm1(k,1,ngt,ng)
                dt02 = dt0*dt0
                do 3201 l=2,ntg(ngt,ng)-1
                    stp = alog(abc1(k,l+1)/abc1(k,l))
                    stm= alog(abc1(k,l-1)/abc1(k,l))
                    at(l) = p5*(stp - stm)/dt0
                    bt(l) = p5*(stp + stm)/dt02
3201            continue
c
                do 3221 l=1,nt_pd
                    tt = d_temp(l)*t(k)
                    itl = nint((tt - tatm1(k,1,ngt,ng))/dt0) + 1
                    if(itl .lt. 2) itl = 2
                    if(itl .gt. ntg(ngt,ng)-1) itl = ntg(ngt,ng) - 1
                    delt = tt - tatm1(k,itl,ngt,ng)
                    sigma = 0.1*abc1(k,itl)*exp(at(itl)*delt +
     -                            bt(itl)*delt*delt)
                    abc2(k,l,1) = sigma*rmix(k,ng)/grav(k)
                    if(ntaug .gt. 1) then
                      abc2(k,l,2) = sigma*
     -                              d_rmix(2,ng)*rmix(k,ng)/grav(k)
                    endif                      
3221            continue
c
3241        continue
c
          else
c
c****       values are only specified at one temperature
c
            do 3401 k=1,nlay+1
                abc2(k,1,1) = 0.1*abc1(k,itlev(k,ngt,ng))*
     -                          rmix(k,ng)/grav(k)
c
                if(ntaug .gt. 1) then
c
c****               find the absorption for perturbed profile
c
                  abc2(k,1,2) = 0.1*abc1(k,itlev(k,ngt,ng))*
     -                            d_rmix(2,ng)*rmix(k,ng)/grav(k)
                
                endif
3401        continue
c
          endif
c
c****         find the gas optical depth in each layer.
c
          do 3441 l=1,nt_pd
              do 3421 k=1,nlyr
                  dtaug(k,nsiext(ni,ne),ne,l,1) = p5*
     -                       (abc2(k+1,l,1) + abc2(k,l,1))*
     -                       (p(k+1) - p(k))
c
                  if(ntaug .gt. 1) then
c
c****                  find the optical depth for mixing ratio 
c                      perturbation at top of layer k (level k)
c
                    dtaug(k,nsiext(ni,ne),ne,l,2) = p5*
     -                       (abc2(k+1,l,1) + abc2(k,l,2))*
     -                       (p(k+1) - p(k))
c
c****                  find the optical depth for mixing ratio 
c                      perturbation at bottom of layer k (level k+1)
c
                    dtaug(k,nsiext(ni,ne),ne,l,3) = p5*
     -                       (abc2(k+1,l,2) + abc2(k,l,1))*
     -                       (p(k+1) - p(k))
                  endif
3421          continue
c
3441      continue
c
          if(istate(1) .ne. 0) then
c
c****         surface pressure is a variable component of the 
c             state vector - find perturbation optical depth
c             for the lowest pressure level. 
c
            k = nlev
            dtaug(k,nsiext(ni,ne),ne,1,1) = p5*
     -                   (abc2(k+1,1,1) + abc2(k-1,1,1))*
     -                   (d_pres(2)*p(k) - p(k-1))
          endif
c
3481  continue
c
c****     compute the derivative of the optical depth w.r.t. wavenumber
c
      dv = wnext(2,ne) - wnext(1,ne)
      if(abs(dv) .gt. tol*wnext(1,ne)) then
        dvi = real(1.0d0/dv)
        do 3841 m=1,ntaug
            do 3821 l=1,nt_pd
                do 3801 k=1,nlay
                    dtaugdv(k,ne,l,m) = (dtaug(k,2,ne,l,m) - 
     -                                   dtaug(k,1,ne,l,m))*dvi
3801            continue
3821        continue
3841    continue
c
      else
c
        do 3881 m=1,ntaug
            do 3861 k=1,nlay
                dtaugdv(k,ne,1,m) = 0.0
                dtaugdv(k,ne,2,m) = 0.0
3861        continue  
3881    continue
      endif
c
      return
c
c***   end of file encountered
c
4001  io_end(ne) = 1
c
      wneof(2) = wnext(nsiext(1,ne),ne)
c
c****     set all gas optical depths and derivatives to zero
c
      wnext(nsiext(2,ne),ne) = 1.01d0*wnmax
      do 4041 m=1,ntaug
          do 4021 k=1,nlay
              dtaug(k,nsiext(1,ne),ne,1,m) = 0.0
              dtaug(k,nsiext(2,ne),ne,1,m) = 0.0
              dtaugdv(k,ne,1,m) = 0.0
              dtaug(k,nsiext(1,ne),ne,2,m) = 0.0
              dtaug(k,nsiext(2,ne),ne,2,m) = 0.0
              dtaugdv(k,ne,2,m) = 0.0
4021      continue
4041  continue
c
      return
c 
c****   file read error encountered
c
4201  io_err(ne) = 1
      wneof(2) = wnext(nsiext(1,ne),ne)
      write(*,*) 'gaslbltau: error reading unit iuabc=',iuabc(ngt,ng),
     -            ' at wavenumber',wneof(2)
c
      wnext(nsiext(2,ne),ne) = 1.01d0*wnmax
c
c****     set all gas optical depths and derivatives to zero
c
      do 4241 m=1,ntaug
          do 4221 k=1,nlay
              dtaug(k,nsiext(1,ne),ne,1,m) = 0.0
              dtaug(k,nsiext(2,ne),ne,1,m) = 0.0
              dtaugdv(k,ne,1,m) = 0.0
              dtaug(k,nsiext(1,ne),ne,2,m) = 0.0
              dtaug(k,nsiext(2,ne),ne,2,m) = 0.0
              dtaugdv(k,ne,2,m) = 0.0
4221      continue
4241  continue
c 
      return
      end
