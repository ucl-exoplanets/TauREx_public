      function raylei(wleff,ncomp,icomp,volmix)
c
cccccccccccccccccccccccccccccc  r a y l e i  ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes the rayleigh scattering cross section  cc
cc    per molecule (meters**-2)  for atmospheres composed of (1) air, cc
cc    (2) co2, (3) n2, (4) o2, (5) h2, (6) he, or any combination of  cc
cc    these gases.                                                    cc
cc    to find the rayleigh scattering optical depth at any level of   cc
cc    the atmosphere, these cross sections must be multiplied by      cc
cc    the pathlength-integrated number density, n(z)*dz.  in a        cc
cc    hydrostatic atmosphere, n(z)*dz = a0*dp/(u0*grav), where u0     cc
cc    is the mean molecular mass (kg/mole), and a0 is avagadro's      cc
cc    number(kg/kmole)                                                cc
cc                                                                    cc
cc    r e f e r e n c e s :                                           cc
cc                                                                    cc
cc    e.j. mccartney, optics of the atmosphere, wiley, p. 187-215,    cc
cc          1976.                                                     cc
cc    a.t. young, revised depolarization corrections for atmospheric  cc
cc          extinction, appl. opt. 19, 3427-3428, 1980.               cc
cc    c.w. allen, astrophysical quantities, athlone press, p. 87,     cc
cc         1964.                                                      cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc    ncomp = number of major atmosphere constituents                 cc
cc    icomp = atmosphere constituent index                            cc
cc            (1) air  (2) co2  (3) n2  (4) o2 (5) h2 (6) he          cc
cc   volmix = volume mixing ratio of each constituent                 cc
cc    wleff = effective wavelength (microns)                          cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    raylei - raleigh scattering cross section per molecule (m**2)   cc
cc                                                                    cc
cccccccccccccccccccccccccccccc  r a y l e i  ccccccccccccccccccccccccccc
c
      implicit none
c
      double precision delta(6),a(6),b(6),wl2i,r,aniso,r2,sum
c
      integer ncomp,icomp(ncomp)
      integer i
      real volmix(ncomp)
      real raylei,wleff
      double precision cnst
c
c       depolarization factors for air, co2, n2, o2 (young, 1980)
c
      data delta/0.02790d0,0.0780d0,0.0210d0,0.0580d0,0.0000d0,0.0000d0/
c
c***    wavelength dependence coefficients for the refractive index
c       (allen, 1964) (note: wavelengths must be in microns)
c
      data a/2.871d-04,4.39d-04,2.906d-04,2.663d-04,1.358d-4,3.48d-5/
      data b/5.67d-03,6.4d-03,7.7d-03,5.07d-03,7.52d-3,2.3d-3/
c
c****   Define the constant multiplier, cnst.  this 
c	constant is equal to 24*pi**3/(1.e-24*L**2), 
c	where L = loschmidt's number (mks units) ,
c	(L = 2.687e25 molecules/m3) and the factor 1.e-24 is 
c	needed to convert wl**4 from microns to meters.
c
      data cnst /1.031d-24/
c
      wl2i=1.0d0/(wleff*wleff)
      sum = 0.
c
      do 1201 i=1,ncomp
            r = 1.0d0 + a(icomp(i))*(1.0d0 + b(icomp(i))*wl2i)
            aniso = (6.0d0 + 3.0d0*delta(icomp(i)))/
     -              (6.0d0 - 7.0d0*delta(icomp(i)))
            r2 = r*r + 2.0d0
            sum = sum + volmix(i)*aniso*
     -            ((r*r - 1.)/r2)**2
1201  continue
c
      raylei = real(cnst*wl2i*wl2i*sum)
c Danie Liang
      raylei = 0.
c
      return
      end
