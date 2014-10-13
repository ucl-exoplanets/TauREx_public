c
cccccccccccccccccccccccccccccc  p a r a m  ccccccccccccccccccccccccccccc
cc								      cc
cc    p u r p o s e :                                                 cc
cc								      cc
cc    this include file specifies the dimensions of variables used    cc
cc    in the program smt_do.                                          cc
cc								      cc
cc    q u a n t i t i e s :                                           cc
cc								      cc
cc      nsp - Maximum number of spectral intervals 		      cc
cc       kp - Maximum number of vertical levels.		      cc
cc    nmode - Maximum number of aerosol particle modes.               cc
cc     ngas - Maximum number of absorbing gases.		      cc
cc     ngrp - Maximum number of SMT grouping bins.		      cc
cc    mxumu - Maximum number of user output zenith angles in discrete cc
cc            ordinate model.					      cc
cc    mxphi - Maximum number of output azimuth angles in the discrete cc
cc            ordinate model.					      cc
cc    mxmom - Maximum number of legendre polynomial moments in the    cc
cc            discrete ordinate method.				      cc
cc    mxulv - Maximum number of output levels in the discrete 	      cc
cc            ordinate model.					      cc
cc    mxrad - Maximum number of output radiance levels.		      cc
cc     nsol - Maximum number of solar zenith angles.		      cc
cc      nex - Maximum number of spectrally-dependent optical          cc
cc            (gases, aerosol modes, surface albedos, and solar flux. cc
cc     mxpd - maximum number of radiance partial derivatives          cc  
cc								      cc
cccccccccccccccccccccccccccccc  p a r a m  ccccccccccccccccccccccccccccc
c
      integer kp, nmode, ngas, ngrp, mxrad, nsol, nex,
     -        mxlout, mxulv, mxumu, mxphi, ngtmax, mxmom, mxpd
c
      parameter (kp=70, nmode=10, ngas=50, ngtmax=3,  
     -           ngrp=512, mxumu=16, mxphi=16, mxmom=200, 
     -           mxulv=kp, mxrad=kp, nsol=10, mxlout=3,
     -           nex=2*ngas + nmode + 3, mxpd = 10)
