      subroutine aerosol(ne,mode,iu,wn,wnmin,wnmax,
     -                   aerqext,aerqsca,aerg0,aerpmom,
     -                   dpmomdv,dqextdv,dqscadv,dg0dv,
     -                   nmomaer,nmom_mx,nsiext,wnext,
     -                   wn_eof,io_end,io_err)
c
cccccccccccccccccccccccccc  a e r o s o l  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine interpolates aerosol optical properties to the  cc
cc    standard wavenumber grid used by line-by-line programs.         cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc         ne : index of radiatively active component                 cc
cc       mode : aerosol particle mode index                           cc
cc       nmom : number of legendre polynomial coefficients for the    cc
cc              scattering phase function expansion.                  cc
cc          iu: unit for aerosol optical property i/o                 cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    aerqext(ni,mode): aerosol extincition efficiency                cc
cc    aerqsca(ni,mode): aerosol scattering efficiency                 cc
cc      aerg0(ni,mode): scattering asymmetry parameter                cc
cc    dqextdv(mode): rate of change of qext with wavenumber           cc
cc    dqscadv(mode): rate of change of qsca with wavenumber           cc
cc      dg0dv(mode): rate of change of g with wavenumber              cc
cc    aerpmom(ni,mom,mode) : legendre polynomial momemt of aerosol    cc
cc                         scattering phase function                  cc
cc    dpmomdv(0,mode) : rate of change of pmom with wavenumber        cc
cc                                                                    cc
cccccccccccccccccccccccccc  a e r o s o l  ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
c****   aerosol moment variables
c
      integer nmomaer(2,nmode),nmom_mx(nmode)
c
c****    variables for the output wavenumber grid
c
      integer ne,mode,iu
      integer nsiext(2,nex)
      integer io_err(nex),io_end(nex)
      integer mom,ni
c
      double precision wn,wnmin,wnmax,wnext(2,nex),wn_eof(2,nex)
c
      real dvi
c
c****   input aerosol wavelength-dependent variables
c
      real aerqext(2,nmode),aerqsca(2,nmode),aerg0(2,nmode),
     -     aerpmom(2,0:mxmom,nmode),dpmomdv(0:mxmom,nmode),
     -     dqextdv(nmode),dqscadv(nmode),dg0dv(nmode)
c
      if(nsiext(2,ne) .eq. 0) then
c        write(*,*) 'intializing aerosol properties for mode: ',mode
c
c****    initialize spectral interval counter for this aerosols mode
c
        nsiext(1,ne) = 1
        nsiext(2,ne) = 2
c
c****     initialize the aerosol optical properties
c
        wnext(2,ne) = 0.0d0
        aerqext(2,mode) = 0.0
        aerqsca(2,mode) = 0.0
        aerg0(2,mode) = 0.0
        aerpmom(2,0,mode) = 1.0
        nmomaer(2,mode) = 1
        do 1201 mom=1,mxmom
            aerpmom(2,mom,mode) = 0.0
1201    continue
c
      endif
c
c****    swap spectral interval counters
c
2001  ni = nsiext(1,ne)
      nsiext(1,ne) = nsiext(2,ne)
      nsiext(2,ne) = ni
c
c****      read aerosol extinction efficiencies and phase functions
c
      read(iu,err=4201,end=4001) nmomaer(ni,mode),wnext(ni,ne),
     -       aerqext(ni,mode),aerqsca(ni,mode),aerg0(ni,mode),
     -       (aerpmom(ni,mom,mode),mom=0,nmomaer(ni,mode))
c
c****     determine if this segment is in desired spectral window
c
      if(wnext(ni,ne) .le. wnmin) wn_eof(1,ne) = wnext(ni,ne)
      if(wnext(ni,ne) .le. wn) go to 2001
      wn_eof(2,ne) = wnext(ni,ne)
c
c****    find the derivative of each quantity wrt wavenumber
c
      dvi = real(1.0d0/(wnext(2,ne) - wnext(1,ne)))
      dqextdv(mode) = (aerqext(2,mode) - aerqext(1,mode))*dvi
      dqscadv(mode) = (aerqsca(2,mode) - aerqsca(1,mode))*dvi
      dg0dv(mode) = (aerg0(2,mode) - aerg0(1,mode))*dvi
c
c****   find the mean values of the phase function moments
c        - the number of phase function moments is the maximum 
c          of the two at ni = 1,2
c
      nmom_mx(mode) = nmomaer(1,mode)
      if(nmomaer(2,mode) .gt. nmom_mx(mode)) 
     -   nmom_mx(mode) = nmomaer(2,mode)
c
      do 3101 mom=1,nmom_mx(mode)
          dpmomdv(mom,mode) = (aerpmom(2,mom,mode) - 
     -                         aerpmom(1,mom,mode))*dvi
3101  continue
c
      return
c
4001  write(*,*) 'End of aerosol file for mode',mode,
     -           'after wavenumber',wnext(nsiext(1,ne),ne)
c
      io_end(ne) = 1
      wnext(ni,ne) = 1.010d0*wnmax
      aerqext(ni,mode) = 0.0
      aerqsca(ni,mode) = 0.0
      aerg0(ni,mode) = 0.0
      dqextdv(mode) = 0.0
      dqscadv(mode) = 0.0
      dg0dv(mode) = 0.0
      aerpmom(ni,0,mode) = 1.0
      nmomaer(ni,mode) = 1
      dpmomdv(0,mode) = 0.
      do 4101 mom=1,mxmom
          aerpmom(ni,mom,mode) = 0.0
          dpmomdv(mom,mode) = 0.0
4101  continue
c
      return
c 
4201  write(*,*) 'Error in aerosol file for mode',mode,
     -           'after wavenumber',wnext(nsiext(1,ne),ne)
c
      io_err(ne) = 1
      wnext(ni,ne) = 1.010d0*wnmax
      aerqext(ni,mode) = 0.0
      aerqsca(ni,mode) = 0.0
      aerg0(ni,mode) = 0.0
      dqextdv(mode) = 0.0
      dqscadv(mode) = 0.0
      dg0dv(mode) = 0.0
      aerpmom(ni,0,mode) = 1.0
      nmomaer(ni,mode) = 1
      dpmomdv(0,mode) = 0.0
      do 4221 mom=1,mxmom
          aerpmom(ni,mom,mode) = 0.0
          dpmomdv(mom,mode) = 0.0
4221  continue
c
      return
      end
