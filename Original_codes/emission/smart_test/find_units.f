      subroutine find_units(iuin,iunits,wn,units)
c
cccccccccccccccccccccc  f i n d _ u n i t s  ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine sets the scaling factor for converting radiance cc
cc    units between several choices, specified by the index, iunits.  cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc       iuin - index of input radiance units:                        cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) ergs/s/cm**2/sr/cm-1                               cc
cc              5) photons/s/m**2/sr/micron                           cc
cc     iunits - index of output radiance units:                       cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) ergs/s/cm**2/sr/cm-1                               cc
cc              5) photons/s/m**2/sr/micron                           cc
cc         wn - wavenumber (cm-1)                                     cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      units - conversion factor needed to convert units from        cc
cc              iuin to iunits.                                       cc
cc                                                                    cc
cccccccccccccccccccccc  f i n d _ u n i t s  ccccccccccccccccccccccccccc
c
      implicit none
c
      integer iuin,iunits
c
      double precision wn
c
      real units
c
      double precision hc
c
c****   define the product of the speed of light, c (2.99792458e+8 m/s)
c       and Plancks Constant, h (6.626068e-34 m**2 kg/s)
c
      data hc /1.9864421d-25/
c
c****   scale output units to W/m**2/sr/cm**-1
c
      if(iuin .ne. iunits .and. wn .ne. 0.0d0) then
c
        if(iuin .eq. 1) then
c
c****       input units are W/m**2/cm**-1
c
          if(iunits .eq. 2) then
c
c****         convert from W/m**2/cm**-1 to W/m**2/micron
c
            units = real(1.0d-4*wn*wn)
c
          else
c
            if(iunits .eq. 3) then
c
c****           convert from W/m**2/cm**-1 to W/m**2/nm
c
              units = real(1.0d-7*wn*wn)
c
            else
c
              if(iunits .eq. 4) then
c
c****           convert from W/m**2/cm**-1 to erg/s/cm**2/cm**-1
c
                units = 1000.0
c
              else
c
c****           convert from W/m**2/cm**-1 to photons/s/m**2/micron
c
                units = real(1.0d-6*wn/hc)
c
              endif
c
            endif
          endif
c
        else
c
          if(iuin .eq. 2) then
c
c****          input units are W/m**2/micron
c
           if(iunits .eq. 1) then
c
c****           convert from W/m**2/micron to W/m**2/cm**-1
c
              units = real(1.0d4/(wn*wn))
c
            else
c
              if(iunits .eq. 3) then
c
c****             convert from W/m**2/micron to W/m**2/nm
c
                units = 0.001
c
              else
c
                if(iunits .eq. 4) then
c
c****               convert from W/m**2/micron to erg/s/cm**2/cm**-1
c
                  units = real(1.0d7/(wn*wn))
c
                else
c
c****               convert from W/m**2/micron to photons/s/m**2/micron
c
                  units = real(1.0d-2/(wn*hc))
c
                endif
c                  
              endif
c
            endif
c
          else
c
            if(iuin .eq. 3) then
c
c****           the input units are W/m**2/nm
c
              if(iunits .eq. 1) then
c
c****             convert from W/m**2/nm to W/m**2/cm**-1
c
                 units = real(1.0d7/(wn*wn))
c
              else
c
                if(iunits .eq. 2) then
c
c****               convert from W/m**2/nm to W/m**2/micron
c
                  units = 1000.0
c
                else
c
                  if(iunits .eq. 4) then
c
c****                 convert from W/m**2/nm to erg/s/cm**2/cm**-1
c
                    units = real(1.0d10/(wn*wn))
c                    
                  else
c
c****                 convert from W/m**2/nm to photons/s/m**2/micron
c
                    units = real(10.0d0/(wn*hc))
c
                  endif
c
                endif
c
              endif
c
            else
c
              if(iuin .eq. 4) then
c
c****            the input units are erg/s/cm**2/cm**-1
c
                if(iunits .eq. 1) then
c
c****               convert from erg/s/cm**2/cm**-1 to W/m**2/cm**-1
c
                  units = 0.001
c
                else
c
                  if(iunits .eq. 2) then
c
c****                convert from erg/s/cm**2/cm**-1 to W/m**2/micron
c
                    units = real(1.0d-7*wn*wn)
c
                  else
c
                    if(iunits .eq. 3) then
c
c****                   convert from erg/s/cm**2/cm**-1 to W/m**2/nm
c
                      units = real(1.0d-10*wn*wn)
c
                    else                
c
c****                   convert from erg/s/cm**2/cm**-1 to 
c                       photons/sec/m**2/micron
c
                      units = real(1.0d-9*wn/hc)
c
                    endif
c
                  endif
c
                endif
c         
              else
c
c****             the input units are photons/sec/m**2/micron
c
                if(iunits .eq. 1) then
c
c****              convert from photons/sec/m**2/micron to W/m**2/cm**-1
c
                  units = real(1.0d6*hc/wn)
c
                else
c
                  if(iunits .eq. 2) then
c
c****                 convert from photons/sec/m**2/micron to 
c                     W/m**2/micron
c
                    units = real(100.0d0*wn*hc)
c
                  else
c
                    if(iunits .eq. 3) then
c
c****                   convert from photons/sec/m**2/micron to 
c                       W/m**2/nm
c
                      units = real(0.1d0*wn*hc)
c
                    else
c
c****                   convert from photons/sec/m**2/micron to 
c                       erg/sec/cm**2/sr/cm-1
c
                      units = real(1.0d9*hc/wn)
c
                    endif  
c
                  endif
c
                endif
c
              endif
c
            endif
c
          endif
c
        endif
c
      else
c
        units = 1.0
c
      endif
c
      return
      end
