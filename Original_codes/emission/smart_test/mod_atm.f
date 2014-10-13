      subroutine mod_atm(iuatm,iumix,new_pt,new_mix,istate,
     -           atmfile,ifrmatm,frmatm,iskatm,icp,ict,scp,
     -           nlev,nlyr,levout,nlout,k_out,pout,wgtatm,wgtgas,
     -           ngases,igas,mixfile,ifrmmix,frmmix,ioffmix,
     -           izmix,imix,icpmix,icmix,scpmix,scmix,pd_frac,
     -           radius,sgrav,ratm,p,t,alt,z,grav,rmix,dp_dp)
c
cccccccccccccccccccccccccccccc  mod_atm  ccccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine calls the subroutines atmstr to set up the      cc
cc    atmospheric structure, and then call the subroutine readmix     cc
cc    to read gas mixing ratios at each model level.                  cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc    atmfile - name of input atmospheric structure file              cc
cc      iuatm - unit from which atmospheric structure is read.        cc
cc     wgtatm -  mean molecular weight of atmosphere (kg/kmole)       cc
cc    ifrmatm - format for reading pressure and temperature:          cc
cc               1) formatted 2) unformatted 3) list directed.        cc
cc     frmatm - format for reading p and t (ifrmatm = 1 only).        cc
cc     iskatm - number of lines to skip at top of atmfile.            cc
cc        icp - column containting pressures.                         cc
cc        ict - column containting temperattures.                     cc
cc        scp - scaling factor to convert pressure to Pascals.        cc
cc          p - pressure at each model level (Pascals).               cc
cc          t - temperature at each model level (K).                  cc
cc       grav - gravitational acceleration at each level (m/s^2).     cc
cc     levout - output level flag:                                    cc
cc              1) top of the atmosphere only, 2) surface only,       cc
cc              3) top of atmosphere and surface,                     cc
cc              4) one or more pressure levels',                      cc
cc     ngases - number of absorbing gases included in model.          cc
cc    mixfile - name of file with gas mixing ratios                   cc
cc      iumix - unit from which mixing ratios are read.               cc
cc    ifrmmix - format for reading gas mixing ratios:                 cc
cc              1) formatted 2) unformatted 3) list directed.         cc
cc     frmmix - format for readinggas mixing ratios(ifrmmix=1 only).  cc
cc    ioffmix - number of lines to skip at top of mixfile.            cc
cc      izmix - type of vertical coordinate:                          cc
cc              1) pressure,   2) altitude (km)                       cc
cc       imix - type of mixing ratios: 1) volume, 2) mass.            cc
cc     icpmix - column containting vertical coordinater (p, z).       cc
cc      icmix - column containting gas mixing ratio                   cc
cc     scpmix - factor needed to scale vertical coordinate (Pa, km)   cc
cc      scmix - factor needed to scale gas mixing ratios.             cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc       ratm - mean atmospheric gas constant (j/kg/k).               cc
cc       nlyr - number of layers in the model atmosphere (nlev-1).    cc
cc       nlev - number of levels in the model atmosphere (<kp).       cc
cc      k_out - index of output level.                                cc
cc      p_out - output pressure level(s) (bars).                      cc
cc      dp_dp - local pressure gradient near each output level.       cc
cc                                                                    cc
cccccccccccccccccccccccccccccc  mod_atm  ccccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      character*132 atmfile,mixfile(ngas)
      character*6 gascode(ngas)
      character*40 frmmix(ngas),frmatm
c
      integer np_mx,istate(nex)
c
      integer iuatm,ifrmatm,iskatm,icp,ict,nlyr,levout,nlout, ngases,
     -        iumix,new_pt,new_mix
      integer iopen,k,n,ng,nlmax
      integer k_out(mxlout),ifrmmix(ngas),ioffmix(ngas), izmix(ngas),
     -        igas(ngas),imix(ngas),icpmix(ngas),icmix(ngas)
c
      real ratm,radius,scp,p_pasc
      real pout(mxlout),dp_dp(mxlout)
      real scmix(ngas),scpmix(ngas)
      real pd_frac(nex)
c
c****   atmospheric structure variables
c
      real p(kp),t(kp),alt(kp),grav(kp),z(kp),rmix(kp,ngas),
     -     wgtgas(ngas),wgtatm,sgrav
      integer nlev
c
      data gascode /'H2O','CO2','O3','N2O','CO',
     -              'CH4','O2','NO','SO2','NO2',
     -              'NH3','HNO3','OH','HF','HCl',
     -              'HBr','HI','ClO','OCS','H2CO',
     -              'HOCl','N2','HCN','CH3Cl','H2O2',
     -              'C2H2','C2H6','PH3','COF2','SF6',
     -              'H2S','HCOOH','HO2','O','ClONO2',
     -              'NO+','HOBr','C2H4','CH3OH','H2','He',
     -               9*'Gas'/
c
c****  read the atmospheric structure file:
c
      ratm = 8314./wgtatm
c
      if(new_pt .eq. 1) then
        iopen = 1
        np_mx = kp
c
        call atmstr(atmfile,iuatm,iopen,ifrmatm,frmatm,
     -              iskatm,nlmax,icp,ict,scp,ratm,
     -              radius,sgrav,p,t,z,alt,grav,np_mx,nlev)
c
        close(iuatm)
c
      endif
c
c****    define number of layers and number of output levels
c
      nlyr = nlev - 1
c
      if(levout .eq. 1) then
        k_out(1) = 1
        pout(1) = 0.0
        dp_dp(1) = 0.0
      else
        if(levout .eq. 2) then
          k_out(1) = nlev
          pout(1) = 1.e-5*p(nlev)
          dp_dp(1) = 0.0
        else
          if(levout .eq. 3) then
            k_out(1) = 1
            pout(1) = 0.0
            dp_dp(1) = 0.0
            k_out(2) = nlev
            pout(2) = 1.e-5*p(nlev)
            dp_dp(2) = 0.0
          else 
            do 1061 n=1,nlout
                p_pasc = 1.e5*pout(n)
                do 1021 k=1,nlev
                    k_out(n) = k
                    if(p(k) .gt. p_pasc) go to 1041
1021                continue            
1041            continue
1061         continue
           endif
c
        endif
      endif
c
      if(levout .gt. 3) then
c
c****     find the pressure derivative between levels ajacent to each
c         output pressure level.
c
        do 1081 k=1,nlout
            if(k_out(k) .gt. 1) then
              dp_dp(k) = (p(k_out(k)) - 1.e5*pout(k))/
     -                   (p(k_out(k)) - p(k_out(k)-1))
            else
              dp_dp(k) = (p(k_out(k)) - 1.e5*pout(k))/
     -                   (p(k_out(k)+1) - p(k_out(k)))
            endif
1081    continue
      endif
c
      if(new_mix .eq. 1) then
c
c****     read gas mixing ratios
c
        do 1141 n=1,ngases
            ng = n
            nlmax = 0
            iopen = 1
c
            call readmix(mixfile(ng),iumix,ng,iopen,ifrmmix(ng),
     -                   frmmix(ng),ioffmix(ng),nlmax,izmix(ng),
     -                   imix(ng),icpmix(ng),icmix(ng),
     -                   scpmix(ng),scmix(ng),wgtatm,wgtgas,
     -                   nlev,alt,z,rmix)
c
            close(iumix)
c
1141    continue
c
      endif
c
c****   print the atmospheric structure and gas mixing ratios
c
      write(*,'(/,1a,/,1a,10(3x,1a,2x))') 
     - ' Atmospheric Structure and Gas Mixing Ratios',
     - '    alt(km)    p (Pa)      t(K)    g (m/s/s) ',
     - (gascode(igas(n)),n=1,ngases)
      do 2001 k=1,nlev
          write(*,'(13(1pe11.4))') 
     -    alt(k),p(k),t(k),grav(k),(rmix(k,n),n=1,ngases)
2001  continue
c
      if(istate(1) .ne. 0) then
c
c****     pressure is a varible part of the state vector
c         find effective altitude for perturbed surface pressure
c
        t(nlev+1) = t(nlev)
        grav(nlev+1) = grav(nlev)
        z(nlev+1) = alog((1.0 + pd_frac(1))*p(nlev))
        alt(nlev+1) = alt(nlev-1) - 0.0005*ratm*
     -                  (t(nlev+1) + t(nlev-1))*
     -                  (z(nlev+1) - z(nlev-1))/grav(nlev+1)
        do 3001 n=1,ngases
            rmix(nlev+1,n) = rmix(nlev,n)
3001    continue
c
      endif
c
      return
      end
