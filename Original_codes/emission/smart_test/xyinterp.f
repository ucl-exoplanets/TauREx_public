      subroutine xyinterp(xd,zd,xi,zi,nxmax,nymax,nxd,nxi,ny)
c
cccccccccccccccccccccccccc  x y i n t e r p  ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e:                                                  cc
cc                                                                    cc
cc    this subroutine uses linear interpolation for an unequally-     cc
cc    spaced input grid to interpolate a 2-d file to a new grid.      cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc    xd    : vector of length nxd. x points where input array is     cc
cc            known.                                                  cc
cc    zd    : nxd array of values defined at points xd.               cc
cc    nxd   : number of points in input x-array                       cc
cc    xi    : vector of length nxi. x points where array is needed.   cc
cc    zi    : approximate value of array zd at points xi, yi.         cc
cc    nxi   : number of points in x-vector xi                         cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    zi    : value of input function, zd, at points xi.              cc
cc                                                                    cc
cccccccccccccccccccccccccc  x y i n t e r p  ccccccccccccccccccccccccccc
c
      implicit none
c
      integer nxmax,nymax,nxd,nxi,ny
      integer i,j,iori,iord,lstd,l,ll
c
      real xd(nxd),zd(nxmax,nymax),xi(nxi),zi(nxmax,nymax)
      real xp1,xp2
c
      if(nxd .eq. 1) then
        do 1021 j=1,ny
            do 1001 i=1,nxi
                zi(i,j) = zd(1,j)
1001        continue
1021    continue
        return
      endif
c
c****   determine the ordering of the array xi
c
      iori = 0
      if(nxi .gt. 1) then
        if(xi(nxi) - xi(1) .lt. 0.) iori = 1
      endif
c
c****   determine the ordering of the array xd
c
      iord = 0
      if(xd(nxd) - xd(1) .lt. 0.) iord = 1
c
c   use linear interpolation to estimate the value of
c   zi at the desired x values, xi.
c
      if(iord .eq. 0) then
        lstd = 2
        xp1 = xd(1)
        xp2 = xd(2)
c
c****      enter loop over interpolated positions.
c
        do 2062 ll=1,nxi
            l = ll
            if(iori .eq. 1) l = nxi-ll+1
c
c****          determine if xi is within the bounds of xd
c
            if(xi(l) .ge. xd(1) .and. xi(l) .le. xd(nxd)) then
c
c****              determine if xi is between points xp1 and xp2
c
2001          if(xi(l) .ge. xp1 .and. xi(l) .le. xp2) then
                do 2021 j=1,ny
                    zi(l,j) = zd(lstd-1,j) + 
     -                      (zd(lstd,j) - zd(lstd-1,j))*
     -                      (xi(l) - xd(lstd-1))/
     -                      (xd(lstd) - xd(lstd-1))
2021            continue
              else
c
c****             update counters and xp1, xp2
c
                lstd = lstd + 1
                xp1 = xp2
                xp2 = xd(lstd)
                go to 2001
              endif
            else
c
c****           xi(l) is outside the range xd(1)-xd(nxd).  extrapolate
c
              if(xi(l) .lt. xd(1)) then
                xp1 = xd(1)
                xp2 = xd(2)
                lstd = 2
              else
                xp1 = xd(nxd-1)
                xp2 = xd(nxd)
                lstd = nxd
              endif
c
              do 2041 j=1,ny
                  zi(l,j) = zd(lstd-1,j) + 
     -                      (zd(lstd,j) - zd(lstd-1,j))*
     -                      (xi(l) - xd(lstd-1))/
     -                      (xd(lstd) - xd(lstd-1))
2041          continue
            endif
2062    continue
      else
c
c****    the input xd array is monotonically decreasing.
c
        lstd = 2
        xp1 = xd(1)
        xp2 = xd(2)
c
c****      enter loop over interpolated positions.
c
        do 3062 ll=1,nxi
            l = ll
            if(iori .eq. 0) l = nxi-ll+1
c
c****          determine if xi is within the bounds of xd
c
            if(xi(l) .le. xd(1) .and. xi(l) .ge. xd(nxd)) then
c
c****              determine if xi is between points xp1 and xp2
c
3001          if(xi(l) .le. xp1 .and. xi(l) .ge. xp2) then
                do 3021 j=1,ny
                    zi(l,j) = zd(lstd-1,j) + 
     -                      (zd(lstd,j) - zd(lstd-1,j))*
     -                      (xi(l) - xd(lstd-1))/
     -                      (xd(lstd) - xd(lstd-1))
3021            continue
              else
c
c****             update counters and xp1, xp2
c
                lstd = lstd + 1
                xp1 = xp2
                xp2 = xd(lstd)
                go to 3001
              endif
            else
c
c****           xi(l) is outside the range xd(1)-xd(nxd).  extrapolate
c
              if(xi(l) .gt. xd(1)) then
                xp1 = xd(1)
                xp2 = xd(2)
                lstd = 2
              else
                xp1 = xd(nxd-1)
                xp2 = xd(nxd)
                lstd = nxd
              endif
c
              do 3041 j=1,ny
                  zi(l,j) = zd(lstd-1,j) + 
     -                      (zd(lstd,j) - zd(lstd-1,j))*
     -                      (xi(l) - xd(lstd-1))/
     -                      (xd(lstd) - xd(lstd-1))
3041          continue
            endif
3062    continue
      endif
c
      return
      end

