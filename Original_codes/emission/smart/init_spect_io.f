      subroutine init_spect_io(nza,nlout,irad,ifrmout,sza0,clev,
     -                         iuout,iuflx,iuheat,iustat,iutrn,
     -                         name,len,radfile,heatfile,statfile,
     -                         trnfile,flxfile)
c
cccccccccccccccccccccc  i n i t _ s p e c t _ i o  ccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    This subroutine initiazies output file names for program smart  cc
cc                                                                    cc
cc    i n p u t:                                                      cc
cc                                                                    cc
cc        nza - number of solar zenith angles                         cc
cc      nlout - number of levels where radiance spectra are output    cc
cc       irad - index of output file type:                            cc
cc              1) wavelength-dependent fluxes and radiances at the   cc
cc                 computational azimuths and zenith angles at the    cc
cc                 specified output levels, and spectrally-integrated cc
cc                 fluxes and heating rates at each computational     cc
cc                 level.                                             cc
cc              2) wavelength-dependent fluxes, radiances, and        cc
cc                 transmission values at computational zenith angles cc
cc                 and specified output levels, and spectrally-       cc
cc                 integrated fluxes and heating rates at each        cc
cc                 computational level.                               cc
cc              3) wavelength-dependent fluxes and radiances at the   cc
cc                 computational azimuths and zenith angles at the    cc
cc                 specified output levels, wavelength-dependent,     cc
cc                 level-dependent pressure-weighted, flux            cc
cc                 divergences, and spectrally integrated fluxes      cc
cc                 and heating rates at each computational level.     cc
cc              4) wavelength-dependent fluxes, radiances,            cc
cc                 transmission values, wavelength-dependent,         cc
cc                 level-dependent pressure-weighted, flux            cc
cc                 divergences, and spectrally-integrated fluxes      cc
cc                 and heating rates at each computational level.     cc
cc              5) fluxes, radiances, and heating rates,              cc
cc                 at arbitrary azimuths and zenith angles,           cc
cc              6) fluxes, radiances and transmission functions       cc
cc                 at arbitrary zenith angles,                        cc
cc              7) fluxes, radiances, and contribution functions      cc
cc                 at arbitrary zenith angles,                        cc
cc              8) fluxes, radiances, transmission functions,and      cc
cc                 contribution functions at arbitrary zenith angles. cc
cc     ifmout - index of output file format (1) ascii, (2) binary,    cc
cc                                          (3) binary, no header     cc
cc       sza0 - solar zenith angles (degrees)                         cc
cc       clev - string variable used to indicate output level         cc
cc      iuout - unit number for output radiance file                  cc
cc      iutrn - unit number for output transmission/pressure file     cc
cc      iuflx - unit number for output level-dependent fluxes         cc
cc     iuheat - unit number for output solar heating rates            cc
cc     iustat - unit number for output binning statistics             cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      name - character string with parsed output file name          cc
cc       len - length of output file name (stripped of blanks)        cc
cc    radfile - name output file for flux/radiance spectra            cc
cc   statfile - name of output file with binning statistics           cc
cc   heatfile - name of output file with heating/cooling rates        cc
cc    trnfile - name of output file with tau=1 trans/pressure         cc
cc    flxfile - name of output file with level-dependent flux spectra cc
cc                                                                    cc
cccccccccccccccccccccc  i n i t _ s p e c t _ i o  ccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      character*132 radfile(mxlout,nsol),heatfile,statfile,trnfile,
     -             flxfile(nsol)
c
      character*132 filein
      character*1 name(132)
      character*4 outext,flxext,clev(mxlout)
c
      integer nlout,irad,nza,len,ifrmout
      integer iuflx,iuout,iuheat,iustat,iutrn
      integer nl,nz,j,iuf,iuo,iover,iza,nlb,ntb,ic,icmx
c
      real sza0(nsol)
c
      write(*,'(/,1a,/,1a)') ' Output Files: ',
     -      '  Unit    File Type  Filename'
c
      icmx = 10
      ic = 0
7021  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') 
     - ' Enter prefix of output radiance and heating rate file names: '
      read(*,'(1a)') filein
c
c****   parse the input file name
c
      call charsp(filein, name, len,132,nlb,ntb) 
c
      do 7221 nl=1,nlout
          if(nza .gt. 1) then
            do 7201 nz=1,nza
                radfile(nl,nz) = ' '
                iza = nint(sza0(nz))
                write(outext,'(1a2,i2.2)') '.r',iza
                write(radfile(nl,nz),'(132a)') 
     -               (name(j),j=1,len),clev(nl),outext
c                write(*,'(1x,2a)') 'radfile =  ',radfile(nl,nz) 
7201        continue          
          else
            radfile(nl,1) = ' '
            write(radfile(nl,1),'(132a)') 
     -           (name(j),j=1,len),clev(nl),'.rad'
c            write(*,'(1x,2a)') 'radfile = ',radfile(nl,1)
          endif
7221  continue
c            
      if(irad .eq. 3 .or. irad .eq. 4 .or. irad .eq. 7 
     -         .or. irad .eq. 8) then
        if(nza .gt. 1) then
c
c****       create an output flux file name for each zenith angle
c
          do 7241 nz=1,nza
              iza = nint(sza0(nz))
              write(flxext,'(1a2,i2.2)') '.f',iza
              write(flxfile(nz),'(132a)') (name(j),j=1,len),flxext
7241      continue          
        else
c
          flxfile(1) = ' '
          write(flxfile(1),'(132a)') (name(j),j=1,len),'.flx'
        endif
      endif
c
      heatfile = ' '
      write(heatfile,'(132a)') (name(j),j=1,len),'.hrt'
c
      statfile = ' '
      write(statfile,'(132a)') (name(j),j=1,len),'.stat'
      if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -   .or. irad .eq. 8) then
        trnfile = ' '
        write(trnfile,'(132a)') (name(j),j=1,len),'.trn'
      endif
c
c****    open units for output fluxes and radiances.  For 
c        level-dependent fluxes, a separate unit is needed for each
c        output solar zenith angle.  For radiances, a separate unit
c        is needed for each output level and solar zenith angle.
c
      open(iuheat,file=heatfile,form='formatted',
     -             status='new',err=7408)
c
      call charsp(heatfile,name,len,132,nlb,ntb) 
c
      write(*,'(i5,a16,132a)') iuheat,' Heating Rates  ',
     -     (name(j),j=1,len)
c
      open(iustat,file=statfile,form='formatted',
     -             status='new',err=7408)
c
      call charsp(statfile,name,len,132,nlb,ntb) 
c
      write(*,'(i5,a16,132a)') iustat,' Statistics  ',
     -      (name(j),j=1,len)
c
      if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -   .or. irad .eq. 8) then
        if(ifrmout .eq. 1) then
          open(iutrn,file=trnfile,form='formatted',
     -             status='new',err=7408)
        else
          open(iutrn,file=trnfile,form='unformatted',
     -             status='new',err=7408)
        endif
c
        call charsp(trnfile,name,len,132,nlb,ntb) 
c
        write(*,'(i5,a16,132a)') iutrn,' Transmission  ',
     -       (name(j),j=1,len)
      endif
c
      do 7321 nz=1,nza
c
          if(irad .eq. 3 .or. irad .eq. 4 .or. irad .eq. 7 
     -         .or. irad .eq. 8) then
c
c****          open unit for output fluxes 
c
            iuf = iuflx + nz - 1
            if(ifrmout .eq. 1) then
              open(iuf,file=flxfile(nz),form='formatted',
     -             status='new',err=7408)
            else
              open(iuf,file=flxfile(nz),form='unformatted',
     -             status='new',err=7408)
            endif
c
            call charsp(flxfile(nz),name,len,132,nlb,ntb) 
c
            write(*,'(i5,a16,132a)') iuf,' Flux  ',(name(j),j=1,len)
          endif
c
c****       open output radiance file:
c
          do 7301 nl=1,nlout
              iuo = iuout + (nz-1)*nlout + nl - 1
              if(ifrmout .eq. 1) then
                open(iuo,file=radfile(nl,nz),form='formatted',
     -               status='new',err=7408)
              else
                open(iuo,file=radfile(nl,nz),form='unformatted',
     -               status='new',err=7408)
              endif
c
              call charsp(radfile(nl,nz),name,len,132,nlb,ntb) 
c
              write(*,'(i5,a16,132a)') iuo,' Radiances  ',
     -             (name(j),j=1,len)
7301      continue
7321  continue
      go to 7601
c
7408  write(*,'(3(/,1a))') 
     - ' Output file already exists.  Choose option:',
     - ' 1) enter a new file name',' 2) overwrite existing file' 
      read(*,*) iover
      write(*,'(1x,1a,i5)') 'iover =',iover
      if(iover .eq. 1) go to 7021
c
c****   overwrite existing file
c
7601  close(iuheat)
      open(iuheat,file=heatfile,form='formatted',status='unknown')
c
      call charsp(heatfile,name,len,132,nlb,ntb) 
c
      write(*,'(i5,a16,132a)') iuheat,' Heating Rates  ',
     -     (name(j),j=1,len)
c
      close(iustat)
      open(iustat,file=statfile,form='formatted',status='unknown')
c
      call charsp(statfile,name,len,132,nlb,ntb) 
c
      write(*,'(i5,a16,132a)') iustat,' Statistics  ',
     -      (name(j),j=1,len)
c
      if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -   .or. irad .eq. 8) then
        close(iutrn)
        if(ifrmout .eq. 1) then
          open(iutrn,file=trnfile,form='formatted',status='unknown')
        else
          open(iutrn,file=trnfile,form='unformatted',status='unknown')
        endif
c
        call charsp(trnfile,name,len,132,nlb,ntb) 
c
        write(*,'(i5,a16,132a)') iutrn,' Transmission  ',
     -       (name(j),j=1,len)
      endif
c
      do 7661 nz=1,nza
          iuf = iuflx + nz - 1
          close(iuf)
          if(irad .eq. 3 .or. irad .eq. 4 .or. irad .eq. 7 
     -       .or. irad .eq. 8) then
            if(ifrmout .eq. 1) then
              open(iuf,file=flxfile(nz),form='formatted',
     -             status='unknown')
            else
              open(iuf,file=flxfile(nz),form='unformatted',
     -             status='unknown')
            endif
c
            call charsp(flxfile(nz),name,len,132,nlb,ntb) 
c
            write(*,'(i5,a16,132a)') iuf,' Flux  ',(name(j),j=1,len)
          endif
          do 7621 nl=1,nlout
              iuo = iuout + (nz-1)*nlout + nl - 1
              close(iuo)
              if(ifrmout .eq. 1) then
                open(iuo,file=radfile(nl,nz),form='formatted',
     -               status='unknown')
              else
                open(iuo,file=radfile(nl,nz),form='unformatted',
     -           status='unknown')
              endif
c
              call charsp(radfile(nl,nz),name,len,132,nlb,ntb) 
c
              write(*,'(i5,a16,132a)') iuo,' Radiances  ',
     -             (name(j),j=1,len)
7621      continue
7661  continue
c
      return
      end
