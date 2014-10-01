      subroutine planck(itype,nlev,iuin,iuout,w,t,b)
c
cccccccccccccccccccccccccccc  p l a n c k  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine evaluates the planck function at 'nlev'         cc
cc    temperatures at a specific wavelenght or wavenumber.            cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc   itype - calculation type: 1) radinace, 2) flux                   cc
cc    nlev - number of levels.                                        cc
cc    iuin - input frequency/wavelength unit:                         cc
cc           1: wavenumber (cm**-1)                                   cc
cc           2: wavelength (microns)                                  cc
cc           3: wavelength (nanometers)                               cc
cc           4: wavelength (Angstroms)                                cc
cc   iuout - output units desired for planck function:                cc
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
cc       b - planck function (units specified by iuout)               cc
cc    c1fv - numerator of planck function in desired units            cc
cc                                                                    cc
cccccccccccccccccccccccccccc  p l a n c k cccccccccccccccccccccccccccccc
c
      implicit none
c
      integer itype,nlev,iuin,iuout
      integer j
c
      double precision b(nlev)
      real w,t(nlev)
c
      double precision h,bk,c,twopi,c2,fv,wlm,ang,c2fv,arg,c1fv
c
c****   define Planks constant, h (J s), Boltzmann's constant, k (J/K)
c       and the speed of light, c (m/sec).
c       
      data h/6.6262e-34/, bk/1.3806e-23/, c/2.998e+8/
c
      data twopi/6.283185/
c
c****  evaluate the 2nd radiation constant and convert 
c      input wavenumber or wavelength to Hz
c
      write(*,*) 'in planck: iuin,h, bk',iuin,h,bk
      c2 = h/bk
      fv = 100.0d0*c*w
      c2fv = c2*fv
c
c****   specify the output units of the Planck function.
c
      if(itype .eq. 1) then
        ang = 2.0d0
      else
        ang = twopi
      endif
        c1fv = 100.d0*ang*h*fv*fv*fv/c
      write(*,*) 'c2fv, c1fv',c2fv,c1fv
c
c compute the Planck function based on the temperature profile.
c
      do 4201 j = 1,nlev
         arg = c2fv/t(j)
         write(*,*) 'j,c2fv,t(j),arg=',j,c2fv,t(j),arg
         if(arg .gt. 0.0 .and. arg .lt. 200.d0) then
            b(j) = c1fv/(dexp(arg) - 1.d0)
         else
c
c****        use the Rayleigh-Jeans limit
c
           b(j) = c1fv*exp(-arg)
c
         endif
4201  continue
c
      return
      end
