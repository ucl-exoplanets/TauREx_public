      subroutine planck(itype,nlev,iuin,iunits,w,t,b)
c
cccccccccccccccccccccccccccc  p l a n c k  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine evaluates the planck function at 'nlev'         cc
cc    temperatures at a specific wavelenght or wavenumber.            cc
cc    note: if w = 0.0, this version returns a unit value for b       cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc   itype - calculation type: 1) radiance, 2) flux                   cc
cc    nlev - number of levels.                                        cc
cc    iuin - input frequency/wavelength unit:                         cc
cc           1: wavenumber (cm**-1)                                   cc
cc           2: wavelength (microns)                                  cc
cc           3: wavelength (nanometers)                               cc
cc           4: wavelength (Angstroms)                                cc
cc   iunits - output units desired for planck function:               cc
cc           0: unit flux: b(k) = 1.0                                 cc
cc           1: Watts/m**2/cm**-1                                     cc
cc           2: Watts/m**2/micron                                     cc
cc           3: Watts/m**2/nanometer                                  cc
cc           4: Watts/m**2/Angstrom                                   cc
cc           5: Watts/m**2/Hz                                         cc
cc       w - wavenumber (cm**-1) or wavelength                        cc
cc       t - temperature at each level                                cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc       b - planck function (units specified by iunits)               cc
cc    c1fv - numerator of planck function in desired units            cc
cc                                                                    cc
cccccccccccccccccccccccccccc  p l a n c k  ccccccccccccccccccccccccccccc
c
      implicit none
c
      integer itype,nlev,iuin,iunits
      integer k
c
      double precision b(nlev)
      real w,t(nlev)
c
      double precision h,bk,c,twopi,c2,fv,wlm,ang,c2fv,arg,c1fv
c
c****   define Planks constant, h (J s), Boltzmann's constant, k (J/K)
c       and the speed of light, c (m/sec).
c       
      data h/6.6262d-34/, bk/1.3806d-23/, c/2.998d+8/
c
      data twopi/6.283185d0/
c
c****  evaluate the 2nd radiation constant and convert 
c      input wavenumber or wavelength to Hz
c
      if(iunits .eq. 0) then
c
c****   return a unit flux
c
        do 1001 k=1,nlev
            b(k) = 1.0
1001    continue
c
        return
c
      endif
c
      c2 = h/bk
      if(iuin .eq. 1) then
        fv = 100.0d0*c*w
        wlm = 1.0d-2/w
      else
        if(iuin .eq. 2) then
          fv = 1.0d6*c/w
          wlm = 1.0d-2*w
        else
          if(iuin .eq. 3) then
            fv = 1.0d9*c/w
            wlm = 1.0d-9*w
          else
            fv = 1.e10*c/w
            wlm = 1.e-10*w
          endif
        endif      
      endif
      c2fv = c2*fv
c
c****   specify the output units of the Planck function.
c
      if(itype .eq. 1) then
        ang = 2.0d0
      else
        ang = twopi
      endif
      if(iunits .eq. 1) then
        c1fv = 100.d0*ang*h*fv*fv*fv/c
      else
        if(iunits .eq. 2) then
          c1fv = 1.0d-6*ang*h*c*c/(wlm**5)
        else
          if(iunits .eq. 3) then
            c1fv = 1.0d-9*ang*h*c*c/(wlm**5)
          else
            if(iunits .eq. 4) then
              c1fv = 1.0d-10*ang*h*c*c/(wlm**5)
            else
              c1fv = ang*h*fv*fv*fv/(c*c)
            endif
          endif
        endif
      endif
c
c compute the Planck function based on the temperature profile.
c
      do 4201 k = 1,nlev
         arg = c2fv/t(k)
         if(arg .lt. 200.d0) then
            b(k) = c1fv/(dexp(arg) - 1.d0)
         else
c
c****        use the Rayleigh-Jeans limit
c
           b(k) = c1fv*exp(-arg)
c
         endif
4201  continue
c
      return
      end
