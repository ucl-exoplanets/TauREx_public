      subroutine solar(iusol,ne,nsivar,iend,wnmax,wn,wnvar,
     -                 sol0,dsoldv)
c
ccccccccccccccccccccccccccccc  s o l a r   ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine interpolates solar fluxes  to the standard      cc
cc    standard wavenumber grid used by line-by-line programs.         cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc      iusol - unit number of file with solar fluxes                 cc
cc      wnmax - maximum wavenumber                                    cc
cc     nsivar - spectral interval initialization flag                 cc
cc      wnvar - wavenumber at which solar fluxes are specified        cc
cc         wn - current wavenumber                                    cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      wnvar - wavenumber at which solar fluxes are specified        cc
cc       sol0 - solar flux at wavenumber wnvar                        cc
cc     dsoldv - derivative of solar flux with respect to wavenumber   cc
cc       iend - end of file flag                                      cc
cc                                                                    cc
ccccccccccccccccccccccccccccc  s o l a r   ccccccccccccccccccccccccccccc
c
      implicit none
c 
      include 'param.inc'
c
      integer nsivar(2),iend(nex),iusol,ne
      integer ni
c
      double precision wn,wnmax,wnvar(2)
c
      real sol0(2),dsoldv
c
      if(nsivar(2) .eq. 0) then
c
c****    initialize spectral interval counter for this aerosols mode
c
        nsivar(1) = 1
        nsivar(2) = 2
c
c****     initialize the aerosol optical properties
c
        wnvar(2) = 0.0d0
        sol0(2) = 0.0
c
      endif
c
c****    swap spectral interval counters
c
2001  ni = nsivar(1)
      nsivar(1) = nsivar(2)
      nsivar(2) = ni
c
c****     read the next solar flux value
c
      read(iusol,err=4201,end=4001) wnvar(ni),sol0(ni)
c
c****     determine if this segment is in desired spectral window
c
      if(wnvar(ni) .le. wn) go to 2001
c
c****    find the derivative of each quantity wrt wavenumber
c
      dsoldv = real((sol0(2) - sol0(1))/(wnvar(2) - wnvar(1)))
c
      return
c
4001  write(*,'(/,1a,1pe12.4)') 'End of solar flux file at wavenumber',
     -            wnvar(nsivar(1))
c
      wnvar(ni) = 1.0001d0*wnmax
      sol0(ni) = 0.0
      iend(ne) = 1
c
      return
c 
4201  write(*,'(/,2a,1pe12.4)') 'Error in solar flux file ',
     -           'after wavenumber',wnvar(nsivar(1))
c
      wnvar(ni) = 1.01d0*wnmax
      sol0(ni) = 0.0
      iend(ne) = 1
c
      return
      end
