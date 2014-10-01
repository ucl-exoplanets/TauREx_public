      subroutine charsp(charin,chout,length,nl,nlb,ntb)
c
ccccccccccccccccccccccccccc  c h a r s p  ccccccccccccccccccccccccccccc
cc                                                                   cc
cc    p u r p o s e :                                                cc
cc                                                                   cc
cc    this subroutine takes a character of arbitrary length and      cc
cc    splits it into an array of characters with unit length.  the   cc
cc    number of leading and trailing blanks and the length of the    cc
cc    input character, less these spaces is then found.              cc
cc                                                                   cc
cc    i n p u t :                                                    cc
cc                                                                   cc
cc    charin : input character of length nl                          cc
cc        nl : length of the input character and dimension of chout  cc
cc                                                                   cc
cc    o u t p u t :                                                  cc
cc                                                                   cc
cc     chout : output character array containing charin characters   cc
cc             stripped of leading and trailing blanks               cc
cc    length : number of elements in charin less leading and         cc
cc             trailing blanks                                       cc
cc       nlb : number of leading blanks in charin                    cc
cc       ntb : number of trailing blanks in charin                   cc
cc                                                                   cc
ccccccccccccccccccccccccccc  c h a r s p  ccccccccccccccccccccccccccccc
c
      character*(*) charin
      character*1 chout(nl)
c
c****  find the number of leading blanks
c
      nlb = 0
      do 1001 i=1,nl
          if(charin(i:i) .ne. ' ' .and. 
     -       ichar(charin(i:i)) .ne. 0) go to 1201
          nlb = nlb + 1
1001  continue
c
c****  find the number of trailing blanks
c
1201  ntb = 0
      do 1221 i=1,nl
          ii = nl-i+1
          if(charin(ii:ii) .ne. ' ' .and. 
     -       ichar(charin(ii:ii)) .ne. 0) go to 1401
          ntb = ntb + 1
1221  continue
c
c****  find the length of charin less leading and trailing blanks
c
1401  length = nl - nlb - ntb
      l1 = length
c      do 1421 l=1,nl
c          if(charin(l:l) .eq. '|' .or. charin(l:l) .eq. '~')
c     -    length = length - 1
c1421  continue
c
c****  compose output character array chout
c
      do 1601 i=1,l1
          ii = i+nlb
          chout(i) = charin(ii:ii)
1601  continue
c
      return
      end
