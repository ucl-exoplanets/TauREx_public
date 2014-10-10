      subroutine flux_int(lsolar,lplanck,nlyr,nz0,ninter,delnu,
     -                    dir_s_flx,dn_s_flx,up_s_flx,
     -                    dn_t_flx,up_t_flx,
     -                    dirsoflx,dnsoflx,upsoflx,
     -                    dnthflx,upthflx)
c 
cccccccccccccccccccccccccc   f l u x _ i n t   ccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine integrates level-dependent solar and thermal    cc
cc    fluxes over wavenumber.  These spectrally-integrated values     cc
cc    are used to find heating rates.                                 cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc       nlyr - number of computational model layers                  cc
cc        nzo - index of current solar zenith angle                   cc
cc      delnu - current spectral resolution (cm**-1)                  cc
cc     ninter - counter that records first time routine is called     cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    dir_s_flx - direct downward solar flux at wavenumber wn         cc
cc     dn_s_flx - total downward solar flux at wavenumber wn          cc
cc     up_s_flx - upward solar flux (irradiance) at wavenumber wn     cc
cc     dn_t_flx - downward thermal flux (irradiance) at wavenumber wn cc
cc     up_t_flx - upward thermal flux (irradiance) at wavenumber wn   cc
cc                                                                    cc
cccccccccccccccccccccccccc   f l u x _ i n t   ccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical lplanck,lsolar
c
      integer nlyr,nz0,ninter
      integer nz,k
c
c****   spectrally-dependent output flux and radiance
c
      double precision delnu
c
c****   spectrally dependent fluxes
c
      real dir_s_flx(mxrad,nsol,2),up_s_flx(mxrad,nsol,2),
     -     dn_s_flx(mxrad,nsol,2),up_t_flx(mxrad,2),dn_t_flx(mxrad,2)
c
c****   double precision integration variables for heating rates
c
      double precision dirsoflx(mxulv,nsol),dnsoflx(mxulv,nsol),
     -       upsoflx(mxulv,nsol),dnthflx(mxulv),upthflx(mxulv)
c
c****   specify the solar zenith angle index
c
      nz = nz0
c
      if(ninter .gt. 1) then  
c
c****    add contributions to solar heating rates:
c
        if(lsolar) then
c
c****         o u t p u t    l e v e l   l o o p
c
          do 7001 k=1,nlyr+1
c
c****          define spectrally-integrated solar fluxes for 
c              heating rate calculation.
c
              dirsoflx(k,nz) = dirsoflx(k,nz) + delnu*
     -                         0.5*(dir_s_flx(k,nz,2) + 
     -                              dir_s_flx(k,nz,1))
              dir_s_flx(k,nz,1) = dir_s_flx(k,nz,2)
              dnsoflx(k,nz) = dnsoflx(k,nz) + delnu*
     -                        0.5*(dn_s_flx(k,nz,2) + 
     -                             dn_s_flx(k,nz,1))
              dn_s_flx(k,nz,1) = dn_s_flx(k,nz,2)
              upsoflx(k,nz) = upsoflx(k,nz) + delnu*
     -                        0.5*(up_s_flx(k,nz,2) + 
     -                             up_s_flx(k,nz,1))
              up_s_flx(k,nz,1) = up_s_flx(k,nz,2)
c
7001      continue
c
        endif
c
c****    add contributions to thermal cooling rates:
c
        if(lplanck .and. nz .eq. 1) then
c
c****      define spectrally-integrated thermal fluxes for 
c          heating rate calculation.
c
          do 4541 k=1,nlyr+1
              dnthflx(k) = dnthflx(k) + delnu*0.5*
     -                     (dn_t_flx(k,2) + 
     -                      dn_t_flx(k,1))
              dn_t_flx(k,1) = dn_t_flx(k,2)
              upthflx(k) = upthflx(k) + delnu*0.5*
     -                     (up_t_flx(k,2) +
     -                      up_t_flx(k,1))
              up_t_flx(k,1) = up_t_flx(k,2)
4541      continue  
        endif
      else
c
c*****     initialize flux integration arrays
c
        do 4601 k=1,nlyr+1
            up_s_flx(k,nz,1) = up_s_flx(k,nz,2)
            dn_s_flx(k,nz,1) = dn_s_flx(k,nz,2)
            dir_s_flx(k,nz,1) = dir_s_flx(k,nz,2)
4601    continue
        do 4621 k=1,nlyr+1
            dn_t_flx(k,1) = dn_t_flx(k,2)
            up_t_flx(k,1) = up_t_flx(k,2)
4621    continue  
      endif
c
      return
      end
