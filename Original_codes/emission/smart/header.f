      subroutine header(iu,ifrm,var,value,comment,lenvar,lenval,lenrec)
c
ccccccccccccccccccccccccccccc  h e a d e r  ccccccccccccccccccccccccccccc
cc                                                                     cc
cc    p u r p o s e :                                                  cc
cc                                                                     cc
cc    this subroutine writes an 80-columm character header that        cc
cc    includes input values for the smt routine                        cc
cc                                                                     cc
cc    i n p u t :                                                      cc
cc                                                                     cc
cc        iu - i/o unit number					       cc
cc      ifrm - i/o format: 1) ascii, 2) binary			       cc
cc       var - character specifying output variable name	       cc
cc     value - value of output variable				       cc
cc   comment - character string describing value		       cc
cc    lenvar - length of output variable name			       cc
cc    lenrec - output record length				       cc
cc                                                                     cc
cc    o u t p u t :                                                    cc
cc                                                                     cc
cc     var, value, and comment					       cc
cc                                                                     cc
ccccccccccccccccccccccccccccc  h e a d  e r  cccccccccccccccccccccccccccc
c
      character*(*) value
      character*(*) var
      character*(*) comment
      character*(127) record
      character*1 name0(127),name1(127),name2(127),
     -            name3(127),space(127)
c
      do 1001 i=1,127
          space(i) = ' '
          name0(i) = ' '
          name1(i) = ' '
          name2(i) = ' '
          name3(i) = ' '
1001  continue
c
      lenvr = len(var)
      lenvl = len(value)
      lencm = len(comment)
c
      call charsp(var,name0,lvar,lenvr,nlb,ntb)
      call charsp(value,name1,length,lenvl,nlb,ntb)
      call charsp(comment,name2,lcom,lencm,nlb,ntb)
c
c****    find the number of spaces between the value and the comment
c 
      l0 = lenval - length
      if(l0 .lt. 1) then
        l0 = 1
        l1 = lenrec - (length + lenvar + 4)
        if(lcom .gt. l1) lcom = l1
      endif
c
      record = ' '
c
      write(record,'(1x,126a)') (name0(i),i=1,lenvar),
     - '= ',(name1(i),i=1,length),(space(i),i=1,l0),
     - '/ ',(name2(i),i=1,lcom)
c
      call charsp(record,name3,len3,127,nlb,ntb)
c
      if(ifrm .eq. 1) then
        write(iu,'(126a)') (name3(i),i=1,lenrec)
      else
        write(iu) (name3(i),i=1,lenrec)
      endif
c
      return
      end
