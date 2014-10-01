      subroutine rad_slit(iuout,nz0,ifrm,islit,iscwt,iord,nspt,
     -                    wnmin,wnmax,dwn,width,
     -                    points,wnz,z,iquit)
c
cccccccccccccccccccccccccc  r a d _  s l i t  cccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine uses a weighted-average scheme to map randomly  cc
cc    spaced data onto a regular, but arbitrarily-spaced output grid. cc
cc                                                                    cc
cc    the output of this routine is actually the weighted sum of all  cc
cc    contributions at each point on the output grid, and the sum of  cc
cc    all weights that contribute at that point.  the weighted sum at cc
cc    each point must be divided by the sum of all weights at that    cc
cc    point by the calling program after all input points are         cc
cc    included.                                                       cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc     ifrm - output file format: 1) ascii, 2) binary.                cc
cc    iord - ordering of input grid: 1) wavelength, 2) wavenumber.    cc
cc    iscwt - spectral resolution flag: 0) constant,                  cc
cc                                      1) linear with wn             cc
cc    width - half-width at half max for weighting function.          cc
cc            these quantities are also used to normalize the         cc
cc            spatial/temporal domain of the grid.                    cc
cc      wn1 - wavenumber of each input data point.                    cc
cc            (only the first 2 are used for the output map)          cc
cc        z - vector containing values of input quantity at wn.       cc
cc     nspt - number of z values at each wavenumber.                  cc
cc    wn_io - output wavenumber grid.                                 cc
cc    wnmin - minimum wavenumber in input and output grid             cc
cc    wnmax - maximum wavenumber in input and output grid             cc
cc      dwn - spacing of output spectral grid (cm**-1)                cc
cc                                                                    cc
cc     islit - index specifying type of weighting function.           cc
cc                 1) boxcar (constant) weighting within r0, ry, r0   cc
cc                 2) triangular (approximate slit spectrometer)      cc
cc     nfcn - number of values in fcn table.                          cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc        wn - wavenumber at each output sample point.                cc
cc     spect - value of each ouput spectral quantity at point wn.     cc
cc                                                                    cc
cc    note: this version of slit is designed to accomodate multiple   cc
cc          streams and solar zenith angles.                          cc
cc                                                                    cc
cccccccccccccccccccccccccc  r a d _  s l i t  cccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nfmx,ncd,nv0,nptmx
c
      parameter (ncd=mxrad*mxrad+mxumu*mxphi+6)
      parameter (nv0 = 21)
      parameter (nptmx = 2)
      parameter (nfmx = (mxlout+mxpd+1)*nsol+1)
c
      integer iuout,nz0,ifrm,islit,iscwt,iord,nspt,iquit
      integer nz,nwn,n,j,nc,nwnwt,ndwn,nw,i0,in01,nn,i1,in
c
      integer nfcn
      save nfcn
c
      integer in0(nfmx),in1(nfmx),in2(nfmx),
     -        i0new(nfmx),i1new(nfmx),i2new(nfmx),nout(nfmx)
      save in0,in1,in2,i0new,i1new,i2new,nout
c
      double precision wnmin,wnmax,dwn,width,wnz
c
      real points(nfmx),z(ncd)
c
      real spect(ncd)
c
      double precision dzdwn(ncd),dnuwt,dwn_dist
c
      real dwnz,wni_o(nfmx),wl
      save wni_o
c
      real fcn(nptmx),dfdnu(nptmx)
      save fcn,dfdnu
c
      double precision dnu,dnui,delnu(nptmx),hwhm,dist,wt
      save dnu,dnui,delnu
c
      double precision wni0(nfmx),delwn(nfmx),a(ncd,nv0,nfmx),
     -                 w(nv0,nfmx),wnmin0(nfmx),wnmax0(nfmx),
     -                 wnm1(nfmx),wn0(nfmx),wn1(nfmx),dnu0(nfmx)
      save wni0,delwn,a,w,wnmin0,wnmax0,wnm1,wn0,wn1,dnu0
c
      real zm1(ncd,nfmx),z0(ncd,nfmx)
      save zm1,z0
c
      real wn_io(nv0,nfmx),spectc(ncd)
      save wn_io,spectc
c
      nz = nz0
c
      if (iquit .lt. 0) return
c
      if(points(nz) .le. 1.) then
c
c****    initialize index pointers for wrap-around output buffer.
c
        in0(nz) = 1
        in2(nz) = 0
        wni0(nz) = wnmin
        nout(nz) = 0
c
c****     find index of last contributing element.  Assume the ouput 
c         grid is equally-spaced with spacing, dwn.
c
        nwn = int((wnmax - wnmin)/dwn) + 1
        in1(nz) = nwn
        if(in1(nz) .gt. nv0) in1(nz) = nv0
        delwn(nz) = 0.0
        wn0(nz) = -9999.0
c
c****      initialize output arrays for range extending from i0 to i1:
c
        do 1021 n=in0(nz),in1(nz)
            wn_io(n,nz) = real(wnmin + delwn(nz))
            delwn(nz) = delwn(nz) + dwn
            w(n,nz) = 0.0
            do 1001 j=1,nspt
                a(j,n,nz) = 0.0
1001        continue
1021    continue
c
c****     initialize internal variables
c
        wn0(nz) = wnz
        wnm1(nz) = wnz
        do 1041 nc=1,nspt
            z0(nc,nz) = z(nc)
            zm1(nc,nz) = z(nc)
1041    continue
c
c****          you can't do any more, so return
c
c          write(*,'(/,1a,i5)') 
c     -      'Returning from rad_slit after initialing unit: ',iuout
c
        return
c
      else
c
c****      This value's contribution to the output spectrum depends on 
c          the spectral interval width that it represents.  Because this
c          width is not uniform, it must be estimated from the 
c          distance between adjacent intervals.  
c
        if(points(nz) .eq. 2) then
c
c****         now you know the width of interval 1.  express all other 
c             interval widths relative to this width.
c
          dnu0(nz) = abs(wnz - wn0(nz))
          dnuwt = 1.0
        else
c
          if(dnu0(nz) .ne. 0.0) then
            dnuwt = 0.5*(wnz - wnm1(nz))/dnu0(nz)
          else
            dnuwt = 0.0
          endif
c
        endif
c
      endif
c
c****     find the distance between point wn0 and the 
c         previous point, wnm1.
c
      dwn_dist = wn0(nz) - wnm1(nz)
c
c****    determine if the slit function must be re-defined
c
      if(points(nz) .le. 2 .or. 
     -   (iscwt .eq. 1 .and. dwn_dist .gt. 0.005*wn0(nz))) then
c
c****     define the half-width of the weighting function
c
        if(iscwt .eq. 1) then
           if(wn0(nz) .gt. 0.) then
             hwhm = 0.5*width*wn0(nz)
           else
             hwhm = 0.5*width*wnmin
           endif
        else
           hwhm = width
        endif
c
c****           d e f i n e    s l i t    f u n c t i o n 
c
c****     define weighting functions and their derivative.
c
        if(islit .eq. 1) then
c
c****          b o x c a r    f u n c t i o n
c
          nfcn = 2
          dnu = hwhm
          dnui = 1.0d0/dnu
          delnu(1) = 0.0
          fcn(1) = 1.0
          dfdnu(1) = 0.0
          delnu(2) = hwhm
          fcn(2) = 1.0
        endif
c
        if(islit .eq. 2) then
c
c****          t r i a n g u l a r   f u n c t i o n
c
          nfcn = 2
          dnu = 2.0*hwhm
          dnui = 1.0d0/dnu
          delnu(1) = 0.0
          fcn(1) = 1.0
          dfdnu(1) = real(-1.0/dnu)
          delnu(2) = 2.0*hwhm
          fcn(2) = 0.0
        endif
c
c****    determine if output grid is wide enough to 
c        include weighting function.
c
        nwnwt = int(2.0*delnu(nfcn)/dwn)
c
        if(nwnwt .gt. nv0) then
          write(*,'(1x,1a,/,1x,1a,/,1x,1a/,1x,1a)') 
     -      'Number of output grid points needed to resolve ',
     -      'weighting function exceeds dimension bound (nv0).',
     -      'Use a smaller number of points to resolve ',
     -      'weighting function.'
          write(*,'(1a,2i10,3(1pe13.5))') 
     -       'nwnwt, nfcn, wnz, delnu, dwn',
     -        nwnwt,nfcn,wnz,delnu(nfcn),dwn 
          stop
        endif
      endif
c
c*****     s p e c t r a l    c o n v o l u t i o n   l o o p
c
c****    find the number of output intervals needed between 
c        wn0 and wnm1
c
      ndwn = int(abs(dwn_dist)/dwn) + 1
c
c****    if more than one spectral interval is needed, find 
c        the spectral gradient
c
      if(ndwn .gt. 1) then
        do 2001 nc=1,nspt
            dzdwn(nc) = (z0(nc,nz) - zm1(nc,nz))/dwn_dist
2001    continue
c
      else
        do 2021 nc=1,nspt
            dzdwn(nc) = 0.
2021    continue
      endif
c
c****    enter loop over output spectral intervals
c
      do 2801 nw=1,ndwn
c
          dwnz = real(float(ndwn - nw)*dwn)
          wn1(nz) = wn0(nz) - dwnz
          do 2101 nc=1,nspt
              spect(nc) = real(z0(nc,nz) - dzdwn(nc)*dwnz)
2101      continue
c
          if(wn1(nz) .ge. wnmax .or. iquit .eq. 1) go to 3002
c
          if(points(nz) .le. 1) then
c
c****        if the input spectral grid does not extend beyond end of 
c            output grid, set wn1 to wnmin - delnu.  This is a fudge 
c              to fix a computational problem.
c
            if(wn1(nz) .gt. wnmin - delnu(nfcn)) 
     -         wn1(nz) = wnmin - delnu(nfcn)
          endif
c
c****       find min and max wavenumber of the output grid where this 
c           input value contributes.
c
          wnmin0(nz) = wn1(nz) - delnu(nfcn)
          wnmax0(nz) = wn1(nz) + delnu(nfcn)
c
c****      find index of the output grid that corresponds to wnmin0. 
c          (note: This approach assumes an equally-spaced output grid.)
c 
          i0 = int((wnmin0(nz) - wni0(nz))/dwn)
c
c****        determine if this interval is beyond shortwave limit of
c            output wrap-around buffer.
c
          if(i0 .lt. 0) then
            i0 = 0
          else
c
c****      move the initial point in the wrap-around output buffer by 
c          i0 points.  write any completed output intervals to disk.
c
            if(nout(nz) .eq. 0) then
              in01 = in0(nz)
              nout(nz) = 1
            else
              in01 = in0(nz) + 1
            endif
c
            i0new(nz) = in0(nz) + i0
            if(i0new(nz) .le. nv0) then
              if(i0 .gt. 1) then
                do 2241 n=in01,i0new(nz)
                    if(w(n,nz) .ne. 0.0) then
                      do 2201 j=1,nspt
                          spectc(j) = real(a(j,n,nz)/w(n,nz))
                          a(j,n,nz) = 0.0
2201                  continue
                      w(n,nz) = 0.0
                    else
                      do 2221 j=1,nspt
                          spectc(j) = 0.0
                          a(j,n,nz) = 0.0
2221                  continue
                    endif
                    wl = 1.e4/wn_io(n,nz)
                    wni_o(nz) = wn_io(n,nz)
                    if(iord .eq. 1) then
                      if(ifrm .eq. 1) then
                        write(iuout,'(9(1pe14.6))')
     -                        wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                      else
                        write(iuout) wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                      endif
                    else
                      if(ifrm .eq. 1) then
                        write(iuout,'(9(1pe14.6))')
     -                        wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                      else
                        write(iuout) wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                      endif
                    endif
2241            continue
c
c****            reset the in0 pointer to in0 + i0
c
                in0(nz) = i0new(nz)
                wni0(nz) = wn_io(in0(nz),nz)
              endif
c
            else
c
c****          write all values between in0 and nv0 to disk.  wrap the
c              starting interval around the end of the output buffer
c              and write values between 1 and i0new - nv0 to disk.
c
              do 2341 n=in01,nv0
                  if(w(n,nz) .ne. 0.0) then
                    do 2301 j=1,nspt
                        spectc(j) = real(a(j,n,nz)/w(n,nz))
                        a(j,n,nz) = 0.0
2301                continue
                    w(n,nz) = 0.0
                  else
                    do 2321 j=1,nspt
                        spectc(j) = 0.0
                        a(j,n,nz) = 0.0
2321                continue
                  endif
                  wl = 1.0d4/wn_io(n,nz)
                  wni_o(nz) = wn_io(n,nz)
                  if(iord .eq. 1) then
                    if(ifrm .eq. 1) then
                      write(iuout,'(9(1pe14.6))')
     -                      wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                    else
                      write(iuout) wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                    endif
                  else
                    if(ifrm .eq. 1) then
                      write(iuout,'(9(1pe14.6))')
     -                      wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                   else
                      write(iuout) wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                    endif
                  endif
c
2341          continue
c
c****            write values between beginning of the buffer and in0
c
              in0(nz) = i0new(nz) - nv0
              do 2481 nn=1,in0(nz)
                  if(nn .le. nv0) then
                    n = nn
                    if(w(n,nz) .ne. 0.0) then
                      do 2401 j=1,nspt
                          spectc(j) = real(a(j,n,nz)/w(n,nz))
                          a(j,n,nz) = 0.0
2401                  continue
                      w(n,nz) = 0.0
                    else
                      do 2421 j=1,nspt
                          spectc(j) = 0.0
                          a(j,n,nz) = 0.0
2421                  continue                  
                    endif
c
                  else
c
                    n = 1
                    wn_io(n,nz) = real(wnmin + delwn(nz))
                    delwn(nz) = delwn(nz) + dwn
                    do 2461 j=1,nspt
                        spectc(j) = 0.0
                        a(j,n,nz) = 0.0
2461                continue
                  endif
c
                  wl = 1.0e4/wn_io(n,nz)
                  wni_o(nz) = wn_io(n,nz)
                  if(iord .eq. 1) then
                    if(ifrm .eq. 1) then
                      write(iuout,'(9(1pe14.6))')
     -                      wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                    else
                      write(iuout) wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                    endif
                  else
                    if(ifrm .eq. 1) then
                      write(iuout,'(9(1pe14.6))')
     -                      wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                    else
                      write(iuout) wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                    endif
                  endif
c
2481          continue
c
c****           reset end-of-buffer counters
c
              if(in0(nz) .le. nv0) then
                in1(nz) = in2(nz)
              else
                in0(nz) = 1
                in1(nz) = 0
              endif
              in2(nz) = 0
c
              wni0(nz) = wn_io(in0(nz),nz)
c
            endif
          endif
c
c****       find the index of the maximum output wavenumber bin affected
c           by this input quantity.
c
          i1 = int((wnmax0(nz) - wni0(nz))/dwn)
c
c****       determine if interval is beyond end of wrap-around buffer.
c 
          i1new(nz) = in0(nz) + i1
          if(i1new(nz) .le. nv0) then
c
c****      create new output wavenumbers and initialize output arrays.
c
            do 2521 n=in1(nz)+1,i1new(nz)
                 wn_io(n,nz) = real(wnmin + delwn(nz))
                 delwn(nz) = delwn(nz) + dwn
                 w(n,nz) = 0.0
                 do 2501 j=1,nspt
                     a(j,n,nz) = 0.0
2501            continue
2521        continue
c
            if(i1new(nz) .gt. in1(nz)) in1(nz) = i1new(nz)
            in2(nz) = 0
c          
          else
c
c****         initialize wavenumbers and spectral quantities between 
c             current in1 and end of wrap-around buffer.
c
            do 2621 n=in1(nz)+1,nv0
                wn_io(n,nz) = real(wnmin + delwn(nz))
                delwn(nz) = delwn(nz) + dwn
                w(n,nz) = 0.0
                do 2601 j=1,nspt
                    a(j,n,nz) = 0.0
2601            continue
2621        continue
c
c****         define new termination index (in2) and initialize values 
c             between the beginning of the buffer and this point.
c
            in1(nz) = nv0
            i2new(nz) = i1new(nz) - nv0 
            if(i2new(nz) .gt. nv0) i2new(nz) = nv0
            do 2661 n=in2(nz)+1,i2new(nz)
                wn_io(n,nz) = real(wnmin + delwn(nz))
                delwn(nz) = delwn(nz) + dwn
                w(n,nz) = 0.0
                do 2641 j=1,nspt
                    a(j,n,nz) = 0.0
2641            continue    
2661        continue    
            in2(nz) = i2new(nz)  
            if(in2(nz) .ge. in0(nz)) in2(nz) = in0(nz)      
c
          endif
c
c****      add distance-weighted average contribution to each output bin
c
          if(nout(nz) .eq. 0) then
            in01 = in0(nz)
            nout(nz) = 1
          else
            in01 = in0(nz) + 1
          endif
c
          do 2721 n=in01,in1(nz)
              dist = abs(wn1(nz) - wn_io(n,nz))
              in = int(dnui*dist) + 1
              if(in .lt. nptmx) then
                wt = (fcn(in) + dfdnu(in)*(dist - delnu(in)))*dnuwt
                if(wt .ne. 0.0) then
                  w(n,nz) = w(n,nz) + wt
                  do 2701 j=1,nspt
                      a(j,n,nz) = a(j,n,nz) + wt*spect(j)
2701              continue
                endif
              endif
2721      continue
c
          do 2761 n=1,in2(nz)
              dist = abs(wn1(nz) - wn_io(n,nz))
              in = int(dnui*dist) + 1
              if(in .lt. nptmx) then
                wt = (fcn(in) + dfdnu(in)*(dist - delnu(in)))*dnuwt
                if(wt .ne. 0.0) then
                  w(n,nz) = w(n,nz) + wt
                  do 2741 j=1,nspt
                      a(j,n,nz) = a(j,n,nz) + wt*spect(j)
2741              continue
                endif
              endif
2761      continue
c
          wnm1(nz) = wn0(nz)
c
2801  continue
c
c****   reset old wavenumber and spectral variables
c
      wnm1(nz) = wn0(nz)
      wn0(nz) = wnz
      do 2901 nc=1,nspt
          zm1(nc,nz) = z0(nc,nz)
          z0(nc,nz) = z(nc)
2901  continue
c
c****      g e t    t h e    n e x t    i n p u t    r e c o r d :
c
      if(wn1(nz) .le. wnmax .and. iquit .eq. 0) return
c
c****      flush remaining buffers and quit
c
3002  continue
c
c      write(*,'(1x,1a,i5,1a,i5,1a,1pe14.6)') 
c     - 'rad_slit: flushing buffers and closing unit: ',
c     - iuout,' for variable: ',nz0,' at wavenumber: ',wnz
c
      do 3221 n=in0(nz)+1,in1(nz)
          if(wn_io(n,nz) .le. wnmax .and. 
     -       wn_io(n,nz) .gt. wni_o(nz)) then
            if(w(n,nz) .ne. 0.0) then
              do 3201 j=1,nspt
                  spectc(j) = real(a(j,n,nz)/w(n,nz))
3201          continue
              wni_o(nz) = wn_io(n,nz)
              wl = 1.0d4/wn_io(n,nz)
              if(iord .eq. 1) then
                if(ifrm .eq. 1) then
                  write(iuout,'(9(1pe14.6))')
     -                wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                else
                  write(iuout) wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                endif
              else
                if(ifrm .eq. 1) then
                  write(iuout,'(9(1pe14.6))')
     -                wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                else
                  write(iuout) wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                endif
              endif
            endif
          endif
c
3221  continue
c
      do 3261 n=1,in2(nz)
          if(wn_io(n,nz) .le. wnmax .and. 
     -       wn_io(n,nz) .gt. wni_o(nz)) then
            if(w(n,nz) .ne. 0.0) then
              do 3241 j=1,nspt
                  spectc(j) = real(a(j,n,nz)/w(n,nz))
3241          continue
              wl = 1.0d4/wn_io(n,nz)
              wni_o(nz) = wn_io(n,nz)
              if(iord .eq. 1) then
                if(ifrm .eq. 1) then
                  write(iuout,'(9(1pe14.6))')
     -                 wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                else
                  write(iuout) wn_io(n,nz),wl,(spectc(j),j=1,nspt)
                endif
              else
                if(ifrm .eq. 1) then
                  write(iuout,'(9(1pe14.6))')
     -                   wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                else
                  write(iuout) wl,wn_io(n,nz),(spectc(j),j=1,nspt)
                endif
              endif
            endif
          endif
c
3261  continue
c
      iquit = -1
      close (iuout)
c
      return
      end
