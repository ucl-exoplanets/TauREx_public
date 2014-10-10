      subroutine albedo(ne,iusur,iflext,nsiext,iref,nref,nza,
     -                lamber,umu0,wn,wnmin,wnmax,wnext,
     -                surf_0,dsurfdv,ws,phiw,surf_opt,alb,
     -                wn_eof,io_end,io_err)
c
cccccccccccccccccccccccccc   a l b e d o   cccccccccccccccccccccccccccc
cc                                                                   cc
cc    p u r p o s e :                                                cc
cc                                                                   cc
cc    this subroutine interpolates surface properties to the         cc
cc    common wavenumber grid used by line-by-line programs.          cc
cc                                                                   cc
cc     note: this version calculates the lamber albedo when a        cc
cc           non-lambertian albedo is selected.                      cc
cc                                                                   cc
cc    i n p u t :                                                    cc
cc                                                                   cc
cc       ne: optical property index counter (prior to albedo)        cc
cc   iflext: spectral flag: 1 => find next albedo point              cc
cc   nsiext: spectral index counter                                  cc
cc    iusur: input unit number for surface albedo file.              cc
cc     nref: number of wavelength dependent surface properties       cc
cc       wn: current wavenumber (cm-1)                               cc
cc    wnext: wavenumbers where albedos are specified (cm-1)          cc
cc    wnmax: maximum wavenumber (cm-1)                               cc
cc                                                                   cc
cc    o u t p u t :                                                  cc
cc                                                                   cc
cc   surf_0: surface properties at spectral point, wnext(ni,ne)      cc
cc  surf_opt: interpolated surface albedo at spectral point wn       cc
cc                                                                   cc
cccccccccccccccccccccccccc   a l b e d o   cccccccccccccccccccccccccccc
c
      implicit none
c      
      logical lamber
c
      include 'param.inc'
c
      integer ne,iusur
c
c*****   albedo i/o flags and counters
c
      integer iflext(nex),nsiext(2,nex),nref,iref,nza
      integer ni,n,nz
      integer io_err(nex),io_end(nex)
c
c****   output albedos and spectral gradients
c
      double precision wn,wnmin,wnmax,wnext(2,nex),wn_eof(2,nex)
c
      real surf_0(2,4),dsurfdv(4),ws,phiw,surf_opt(4),alb(nsol)
      real umu0(nsol)
c
      double precision dv,dvi
c
c****    variables for the output wavenumber grid
c
      real dref
c
c****   increment index counter and check new-value flag
c
      ne = ne + 1
c
c****   initialize values if this is the first interval
c
      if(nsiext(2,ne) .eq. 0) then
c
c****      initialize spectral interval counter for this aerosols mode
c
        nsiext(1,ne) = 1
        nsiext(2,ne) = 2
c
c****      initialize the aerosol optical properties
c
        wnext(2,ne) = 0.0d0
        do 1001 n=1,nref
            surf_0(2,n) = 0.0d0
1001    continue
c
      endif
c
      if(iflext(ne) .eq. 1) then
c
c****    swap spectral interval counters
c
2001    ni = nsiext(1,ne)
        nsiext(1,ne) = nsiext(2,ne)
        nsiext(2,ne) = ni
c
c****     read the next surface optical property value for this gas
c
        read(iusur,err=4201,end=4001) wnext(ni,ne),
     -      (surf_0(ni,n),n=1,nref)
c
c****       determine if this segment is in desired spectral window
c
        if(wnext(ni,ne) .le. wnmin) wn_eof(1,ne) = wnext(ni,ne)
c
        if(wnext(ni,ne) .le. wn) go to 2001
c
        wn_eof(2,ne) = wnext(ni,ne)
c
c****    find the derivative of each quantity wrt wavenumber
c
        dvi = 1.0d0/(wnext(2,ne) - wnext(1,ne))
        do 2021 n=1,nref
            dsurfdv(n) = real((surf_0(2,n) - surf_0(1,n))*dvi)
2021    continue
c
        iflext(ne) = 0
c
      endif
c
c****       interpolate  surface properties to this wn
c
      dv = wn - wnext(nsiext(1,ne),ne)
      do 2041 n=1,nref
          surf_opt(n) = real(surf_0(nsiext(1,ne),n) + dsurfdv(n)*dv)
2041  continue
c
      if(lamber) then
c
c****     a lambert albedo is used
c
        do 2061 nz=1,nza
            alb(nz) = surf_opt(1)
2061    continue
c
      else
c
        if(iref .eq. 4) then
c
c****      set wind speed and direction for Cox/Munk model
c
          surf_opt(3) = ws
          surf_opt(4) = phiw
        endif
c
c****    use the disort dref routine to find the flux albedo 
c        for each solar zenith angle
c
        do 2081 nz=1,nza
            alb(nz) = dref(umu0(nz),surf_opt,iref)
2081    continue
c
      endif
c
      return
c
c****   end of surface optical property file:
c
4001  wnext(ni,ne) = wnmax
      do 4021 n=1,nref
          surf_0(ni,n) = 0.0
          dsurfdv(n) = 0.0
4021  continue
      iflext(ne) = 0
      write(*,'(/,1a,1pe14.6)') 
     - 'End of surface albedo file at wavenumber ',wn
      io_end(ne) = 1
c
      return
c 
c****   error in surface optical property file
c
4201  write(*,'(/,1a,1pe12.4)') 
     - 'Error in surface albedo file after wavenumber',
     -  wnext(nsiext(1,ne),ne)
      io_err(ne) = 1
c
      return
      end
