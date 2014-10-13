      subroutine fresnel(rn1,rk1,ang1,rn2,rk2,ang2,rpar,rper,tpar,tper)
c
ccccccccccccccccccccccccccc  f r e s n e l  ccccccccccccccccccccccccccccc
cc                                                                     cc
cc    p u r p o s e :                                                  cc
cc                                                                     cc
cc    this subroutine evaluates fresnel's reflection formula for       cc
cc    arbitrary angle and refractive index, and returns the            cc
cc    parallel and perpendictular polarized components of the          cc
cc    reflectance.                                                     cc
cc                                                                     cc
cc    i n p u t :                                                      cc
cc                                                                     cc
cc     rn1 - real part of the complex refractive index of 1st medium   cc
cc     rk1 - complex part of the complex refractive index of 1st       cc
cc           medium                                                    cc
cc     rn2 - real part of the complex refractive index of 2nd medium   cc
cc     rk2 - complex part of the complex refractive index of 2nd       cc
cc           medium                                                    cc
cc    ang1 - angle of incidence of the wave to the surface normal      cc
cc           (in radians) in  first media                              cc
cc                                                                     cc
cc    o u t p u t :                                                    cc
cc                                                                     cc
cc    ang2 - angle of reflection of the wave to the surface normal     cc
cc           (in radians) in second media                              cc
cc    rpar - parallel-polarized component of the reflectance           cc
cc    rper - perpendictular-polarized component of the reflectance     cc
cc    tpar - parallel-polarized component of the transmittance         cc
cc    tper - perpendictular-polarized component of the transmittance   cc
cc                                                                     cc
cc    a u t h o r :                                                    cc
cc                                                                     cc
cc    dave crisp  4/29/89                                              cc
cc      update 8/4/90 to fix tpar,tper                                 cc
cc                                                                     cc
cc    r e f e r e n c e s :                                            cc
cc  Wendlandt and Hecht, Reflectance Spectroscopy, p.17                cc
cc                                                                     cc
cc  Bennett J.M. and H.E. Bennet, "Polarization," Handbook of optics,  cc
cc      W. Driscoll ed., Mcgraw-Hill Book Company, NY, pg 10-7, 1978.  cc
cc                                                                     cc
cc  Schanda, E., Physical Fundamentals of Remote Sensing, Springer     cc
cc      Verlag, Berlin, pg 28, 1986, 187 pp.                           cc
cc                                                                     cc
ccccccccccccccccccccccccccc  f r e s n e l  ccccccccccccccccccccccccccccc
c
      implicit none
c
      complex refrac,rs,rp,ts,tp,cos2,n1,n2
      real rn1,rk1,rn2,rk2
      real ang1,ang2,rpar,rper,tpar,tper,cos1,arg
c
      if (ang1 .gt. 1.55) then
c
c*****   if ang1 is close to 90 degrees, assume reflectance is unity.
c        set up the correct values for rpar,rper, tpar, tper, 
c        solve for ang2 and return
c
        rpar = 1.0
        rper = 1.0
        tpar = 0.0
        tper = 0.0
c
      else
c
c*****   define the complex index of refraction for the two media
c
        n1 = rn1 - (0.,1.)*rk1
        n2 = rn2 - (0.,1.)*rk2
c
        refrac = n1/n2
        cos1 = cos(ang1)
        cos2 = csqrt((1.,0.) - refrac*refrac*((1.,0) - cos1**2))
c
c****     define the parallel-polarized component of the reflectance
c
        rp = (n2*cos1 - n1*cos2)/(n2*cos1 + n1*cos2)
        rpar = rp*conjg(rp)
c
c****     define the perpendictular-polarized component of  reflectance
c
        rs = (n1*cos1 - n2*cos2)/(n1*cos1 + n2*cos2)
        rper = rs*conjg(rs)
c
c****   define the parallel-polarized component of the transmittance
c
        tp = 2.0*n1*cos1/(n2*cos1 + n1*cos2)
        tpar = (n2*cos2/(n1*cos1))*tp*conjg(tp)
c
c****     define perpendictular-polarized component of  transmittance
c
        ts = 2.*n1*cos1/(n1*cos1 + n2*cos2)
        tper = (n2*cos2/(n1*cos1))*ts*conjg(ts)
c
      endif
c
      arg = real(refrac)*sin(ang1)
      if(arg .gt. -1.0 .and. arg .lt. 1.0) then
        ang2 = asin(real(refrac)*sin(ang1))
      else
        ang2 = 0.5*3.1415926
      endif
c
      return
      end
