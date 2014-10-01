      subroutine readxsec(iu,nxsec,wgtgs,wnmin,wnmax,wl,xsec,nxmx,iquit)
c
cccccccccccccccccccccccccc  r e a d x s e c  ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine reads wavelength-dependent gas absorption       cc
cc    cross-sections.  This version of the code                       cc
cc    assumes no pressure or temperature dependence for these         cc
cc    coefficients.                                                   cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        iu - unit number for absorption coefficeints.               cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     nxsec - number of wavenumbers where absorption coefficients    cc
cc             are specified                                          cc
cc        wl - wavenumbers where absorption coefficients are given    cc
cc      xsec - absorption coefficient at each wavenumber              cc
cc                                                                    cc
cccccccccccccccccccccccccc  r e a d x s e c  ccccccccccccccccccccccccccc
c
      implicit none
c
      integer nskip,n,iu,nxsec,iquit,nxmx
c
      double precision wnmin,wnmax,wl(nxmx)
      real wgtgs,xsec(nxmx)
      real xsec0
      double precision wlmin,wlmax,wl0
c
c****    skip over the file header
c
      wlmin = 1.0d4/wnmax
      wlmax = 1.0d4/wnmin
c
      nskip = 2
      do 1001 n=1,nskip
          read(iu,*)
1001  continue
c
c****    read the gas molecular weight (kg/kmole)
c     
      read(iu,*) wgtgs
c
      nskip = 4
      do 1021 n=1,nskip
          read(iu,*)
1021  continue
c
c****   read the gas absorption cross-section (cm**2)
c
      n = 0
2001  if(n .lt. nxmx .and. iquit .eq. 0) then
        read(iu,*,end=2022,err=2042) wl0,xsec0
c
        if(wl0 .ge. wlmin .and. wl0 .le. wlmax) then
          n = n + 1
          wl(n) = wl0
          xsec(n) = xsec0
        else
          if(n .lt. 2) then 
            n = 1
            wl(n) = wl0
            xsec(n) = xsec0
          else 
            n = n + 1
            wl(n) = wl0
            xsec(n) = xsec0
            go to 2022
          endif
        endif
c
        go to 2001
      else
        if(n .gt. nxmx) then
          write(*,'(/,1a,1pe14.6,/2a,i5)') 
     -    'Error in readxsec after wl = ',wl(n),
     -    'Number of elements in array exceeds dimension bound:',
     -    ' nxmx = ',nxmx
        endif
      endif
c
2022  nxsec = n
c
      iquit = 1
      return
c
2042  iquit = -1
      write(*,'(/,1a,1pe14.6)') 'Error in readxsec after wl = ',wl(n)  
c
      return
      end
