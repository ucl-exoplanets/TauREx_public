      subroutine rad_trn(usrang,lplanck,l_1,l_2,ng0,nstr,numu,
     -                     umu,umu_f,gwt_f,dtau,copi0,g0,
     -                     trnrad,absrad,trn_cmu,abs_cmu)
c
ccccccccccccccccccccccccccc  r a d _ t r n   ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes the effective layer radiance and flux  cc
cc    transmittances and absorptances for each spectral bin.  These   cc
cc    quantities are used in the spectral mapping and jacobian        cc
cc    calculations.                                                   cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc     usrang: logical variable.  If false, use computational angles, cc
cc             if true, find radiances for user specified angles.     cc
cc        l_1: index of first layer for which values are needed       cc
cc        l_2: index of last layer for which values are needed        cc
cc        ng0: spectral bin number: 0=> single scat, >0=> mult scat   cc
cc       nstr: number of zenith angle streams.                        cc
cc       numu: number of output zenith angles.                        cc
cc       dtau: layer optical depth                                    cc
cc      copi0: layer single scattering albedo                         cc
cc         g0: layer asymmetry parameter                              cc
cc        umu: zenith angle of each stream.                           cc
cc      umu_f: gaussian point for each stream.                        cc
cc      gwt_f: gaussian weight for each stream.                       cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    trn_cmu: layer transmittance for each radiance stream.          cc
cc    abs_cmu: layer absorptance for each rediance stream             cc
cc     trnrad: layer transmittance for each user radiance stream.     cc
cc     absrad: layer absorptance for each user radiance stream.       cc
cc                                                                    cc
ccccccccccccccccccccccccccc  r a d _ t r n   ccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical usrang,lplanck
c
      integer l_1,l_2,nstr,numu,ng0
c
      integer k,nze
c
c***    set up direction cosines for flux integration
c
      real umu_f(mxumu),gwt_f(mxumu),umu(mxumu)
c
c*****   dummy state vector optical properties.
c
      real dtau(kp),copi0(kp),g0(kp)
c
c****    layer transmission values for simplified adding method
c
      double precision trnrad(mxumu,kp),absrad(mxumu,kp),
     -                 trn_cmu(mxumu,kp),abs_cmu(mxumu,kp)
      double precision dtauexe,dtauabs
c
c****    find the layer transittances and reflectances for each layer
c
      do 2621 k=l_1,l_2
c
c****       define the delta-m scaled extinction optical depth
c
          dtauexe = dtau(k)
c          dtauexe = dtau(k)*(1.0 - (1.0 - copi0(k))*g0(k)*g0(k))
c
c****       define the absorption optical depth
c
          dtauabs = dtau(k)*copi0(k)
c
c           define the direct transmission and absorption 
c           along each computational stream
c 
          do 2001 nze=1,nstr/2
              trn_cmu(nze+nstr/2,k) = dexp(-dtauexe/umu_f(nze+nstr/2))
              trn_cmu(nstr/2-nze+1,k) = trn_cmu(nze+nstr/2,k)
2001      continue
          if(ng0 .eq. 0 .or. lplanck) then
            do 2021 nze=1,nstr/2
                abs_cmu(nze+nstr/2,k) = copi0(k)*(1.0d0 -
     -                              dexp(-dtauabs/umu_f(nze+nstr/2)))
                abs_cmu(nstr/2-nze+1,k) = abs_cmu(nze+nstr/2,k)
2021        continue
          endif
c
          if(.not. usrang) then
c
c             define the direct transmission along each stream
c 
            do 2201 nze=1,numu
                trnrad(nze,k) = trn_cmu(nze,k)
                absrad(nze,k) = abs_cmu(nze,k)
2201        continue
c
          else
c
c****         usrang = .true.  Recompute radiances on the 
c             computational grid and integrate these over 
c             angle to get fluxes
c
c             define the direct transmission along each stream
c 
            do 2401 nze=1,numu
                trnrad(nze,k) = exp(-dtauexe/abs(umu(nze)))
                absrad(nze,k) = copi0(k)*(1.0d0 - 
     -                          exp(-dtauabs/abs(umu(nze))))
2401        continue
c
          endif
c
2621  continue
c
      return
      end
