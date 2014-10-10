      subroutine init_iu(iuheat,iustat,iuatm,iumie,iumix,iusur,
     -                    iugrp,iuthrm,iutrn,iuaer,iuflx,iuout,
     -                    iugas,iusol0,iusol1,iusol2,iu_pd)
c
cccccccccccccccccccccccccc   i n i t _ i u    cccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    This subroutine initializes the i/o unit numbers for SMART      cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc    nsol - maximum number of solar zenith angles                    cc
cc    nmode - maximum number of particle modes                        cc
cc    mxlout - maximum number of output radiance/flux levels          cc
cc    ngas - maximum number of absorbing gases                        cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    iuheat - output unit for heating and cooling rates.             cc
cc    iustat - output unit for binning statistics                     cc
cc    iuatm - input unit for input atmospheric thermal structure      cc
cc    iumie - first input unit for aerosol optical properties         cc
cc    iumix - first input unit for gas mixing ratios                  cc
cc    iusur - input unit for surface optical properties               cc
cc    iugrp - scratch i/o unit for binned optcal properties           cc
cc    iuthrm - scratch i/o unit for binned thermal radiances          cc
cc    iutrn - output unit with transmission/pressure file             cc
cc    iuaer - i/o unit with combined aerosol optical properties       cc
cc    iuflx - output unit with flx vs pressure                        cc
cc    iuout - output unit for wavenumber-dependent radiances/fluxes   cc
cc    iugas - first input unit with gas gas optical properties        cc
cc    iusol0 - solar flux scratch file used by map_back               cc
cc    iusol1 - first of two solar flux input units                    cc
cc    iusol2 - second of two solar flux input units                   cc
cc    iu_pd - output unit for flux and radiance jacobians             cc
cc                                                                    cc
cccccccccccccccccccccccccc   i n i t _ i u    cccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
c****   define the output unit numbers for heating rates
c        (iuheat), and binning statistics (iustat)
c
      integer iuatm, iumix,iumie, iusur
      integer iustat,iuaer,iugas,iu_pd,iuheat
      integer iugrp,iutrn,iuflx,iusol0,iusol1,iusol2,iuthrm,iuout
c
      iuheat = 9
      iustat = 10
c
c****  define the unit numbers for the input atmospheric thermal 
c      structure, aerosol optical properties, gas mixing ratios, 
c      and surface albedos.
c
      iuatm = 11
      iumie = 12
      iumix = 13
      iusur = 14
c
c****  define a scratch unit for bin numbers at each wavenumber
c
      iugrp = 15
c
c****  define a scratch unit for thermal fluxes, iuthrm
c
      iuthrm = 18
c
c****   define an output unit for optical depths, iutrn
c
      iutrn = 19
c
c****   define input unit for aerosol optical properties, iuaer
c       note: nmode units are needed.
c
      iuaer = 20
c
c****   define the first output unit for the level-dependent fluxes.
c       note: there are up to nsol of these units
c
      iuflx = iuaer + nmode
c
c****   define units for output flux/radiance files
c       note: up to nsol units are needed
c
      iuout = iuflx + nsol
c
c****   define input units for gas absorption coefficeints
c       note: ngas units are needed (this can be a large number).
c
      iugas = iuout + nsol*mxlout
c
c****  note: for solar flux files ordered in increasing wavenumber,
c      only 2 unit numbers are needed.  for large files ordered in
c      increasing wavelength, a large number of files are needed
c      to reverse the order of the solar flux files.  The unit iusol1 
c      is released after the subroutine 'readsol' completes, and is
c      then reused in 'backmap'
c
      iusol0 = iugas + ngas + 1
      iusol1 = iusol0 + 1
      iusol2 = iusol1 + 1
c
c****    define unit number for partial derivitives
c
      iu_pd = iusol2 + 1
c
      return
      end
