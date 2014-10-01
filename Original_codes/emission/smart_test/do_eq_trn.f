      subroutine do_eq_trn(nlyr,nstr,nmom,numu,nphi,ibcnd,
     -                     iunits,iref,usrang,lamber,source,
     -                     umu,phi,umu0,phi0,accur,
     -                     alb0,surf_pr,dtauc,ssalb,pmom,t,
     -                     ttemp,btemp,temis,fbeam,fisot,wng0,
     -                     flup,rfldir,rfldn,uu,albmed,trnmed,nscat)
c
cccccccccccccccccccccccccc  d o _ e q _ t r n  ccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this program initializes discrete ordinate arrays and calls     cc
cc    the Stamnes a discrete ordinate model to find radiances         cc
cc    and fluxes for smt optical properties.                          cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        ng0 - index  of current spectral bin                        cc
cc       nlyr - number of computational model layers                  cc
cc       nstr - number of gaussian zenith angles used in D/O code     cc
cc       nmom - number of phase function moments used in D/O code     cc
cc       numu - number of output zenith angles used in D/O code       cc
cc       nphi - number of output azimuth angles                       cc
cc      ibcnd - boundy condition flag (0- flx/rad, 1-albedos only)    cc
cc     iunits - output units desired for planck function:             cc
cc           0: unit flux: b(k) = 1.0                                 cc
cc           1: Watts/m**2/cm**-1                                     cc
cc           2: Watts/m**2/micron                                     cc
cc           3: Watts/m**2/nanometer                                  cc
cc           4: Watts/m**2/Angstrom                                   cc
cc           5: Watts/m**2/Hz                                         cc
cc       iref - bidirectional reflectance options                     cc
cc              0 - lambert                                           cc
cc              1 - Hapke's BDR model                                 cc
cc              2 - Breon's BDR model; combination of Li + Roujean    cc
cc              3 - Roujean's BDR model                               cc
cc              4 - Cox and Munk glint model                          cc
cc     usrang - output radiances at user angles? (logical: T/F)       cc
cc     lamber - Include a lambertian surface? (Logical: T/F)          cc
cc              note: if lamber = F, use a BRDF is used.              cc
cc     source - include thermal fluxes? (logical: T/F)                cc
cc        umu - emission zenith angle cosines                         cc
cc        phi - emission azimuth angles (degrees)                     cc
cc       umu0 - cosine of solar zenith angles                         cc
cc       phi0 - solar azimuth angles (degrees)                        cc
cc      accur - azimuth convergence accuracy for D/O routine          cc
cc       alb0 - surface albedo                                        cc
cc    surf_pr - surface proerties for non-lambertian BRDF             cc
cc      dtauc - layer optical depth                                   cc
cc      ssalb - layer single scattering albedo                        cc
cc       pmom - layer particle phase function                         cc
cc          t - temperature in each atmospheric                       cc
cc      ttemp - temperature of uppermost model layer (space)          cc
cc      btemp - surface temperature                                   cc
cc      temis - emissivity of uppermost model layer (space)           cc
cc      fbeam - intensity of collimated flux at top of atmosphere     cc
cc      fisot - thermal flux at top of atmosphere                     cc
cc       wng0 - wavenumber of bin                                     cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc       flup - upward (diffuse) flux at each model level             cc
cc     rfldir - downward direct flux at each model level              cc
cc      rfldn - downward diffuse flux at each model level             cc
cc         uu - upward radiance at each level, zenith angle, azimuth  cc
cc     albmed - albedo of system                                      cc
cc     trnmed - transmissivity of system                              cc
cc      nscat - counter for number of scattering calculations         cc
cc                                                                    cc
cccccccccccccccccccccccccc  d o _ e q _ t r n  ccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical usrang,lamber,source
      logical prnt(7), onlyfl, plank, usrtau
c
      integer ibcnd, maxcly, maxmom, maxphi, maxulv, maxumu
      integer nlyr,nmom, nphi, nstr, ntau, numu, iunits
      integer iref,nscat
      integer l,k
c
      real accur, alb0, btemp, fbeam, fisot, phi0, temis, ttemp,
     &          umu0, wvnm,t(kp)
c
      real phi(mxphi),umu(mxumu),surf_pr(4),dtauc(kp),ssalb(kp),
     -       pmom(0:mxmom,kp),temper(0:kp),
     -       utau(mxulv),wng0
c
c***    output variables
c
      real rfldir( mxulv ), rfldn( mxulv ), flup( mxulv ),
     -       dfdt( mxulv ), uavg( mxulv ),uu( mxumu, mxulv, mxphi ), 
     -       albmed( mxumu ),trnmed( mxumu )
c
c****    set logicals
c
      onlyfl = .false.
      plank = source
c
c****         set dimensions of arrays
c
      maxcly = kp
      maxulv = mxulv
      maxumu = mxumu
      maxphi = mxphi
      maxmom = mxmom
c
c****       set wavelength interval (1 cm**-1)
c
      wvnm = wng0
c
c****   turn off all discr_ord model internal print flags
c
      do 1001 l=1,7
          prnt(l) = .false.
1001  continue
c
c****    set number of arbitrary output levels and angles
c 
      usrtau = .false.
      ntau = 0
c
c****   load optical depth, single scattering albedo
c       and phase function moment and temperature arrays
c
      temper(0) = t(1)
      do 1221 k=1,nlyr
c
          utau(k) = 0.0
          temper(k) = t(k+1)
1221  continue
c
        nscat = nscat + 1
c
c****   d i s c r e t e   o r d i n a t e    r o u t i n e
c
      call disort( nlyr, dtauc, ssalb, nmom, pmom, temper, 
     -             wvnm, usrtau, ntau, utau, nstr,iunits,
     -             usrang, numu, umu, nphi, phi, ibcnd, fbeam,  
     -             umu0, phi0, fisot, lamber, iref, 
     -             surf_pr, alb0, btemp, ttemp, temis, 
     -             plank, onlyfl, accur, prnt, 
     -             maxcly, maxulv, maxumu, maxphi, maxmom, 
     -             rfldir, rfldn, flup, dfdt, uavg,uu, 
     -             albmed, trnmed )
c
      return
      end
