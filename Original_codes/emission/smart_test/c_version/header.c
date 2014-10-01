/* header.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__127 = 127;

/* Subroutine */ int header_(integer *iu, integer *ifrm, char *var, char *
	value, char *comment, integer *lenvar, integer *lenval, integer *
	lenrec, ftnlen var_len, ftnlen value_len, ftnlen comment_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_wsfe(cilist *), e_wsfe(void), s_wsue(cilist *), do_uio(
	    integer *, char *, ftnlen), e_wsue(void);

    /* Local variables */
    static integer i__, l0, l1, nlb, ntb, len3, lcom, lvar;
    static char name0[1*127], name1[1*127], name2[1*127], name3[1*127], space[
	    1*127];
    static integer lencm, lenvl, lenvr;
    static char record[127];
    extern /* Subroutine */ int charsp_(char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer length;

    /* Fortran I/O blocks */
    static icilist io___18 = { 0, record, 0, "(1x,126a)", 127, 1 };
    static cilist io___20 = { 0, 0, 0, "(126a)", 0 };
    static cilist io___21 = { 0, 0, 0, 0, 0 };



/* cccccccccccccccccccccccccccc  h e a d e r  ccccccccccccccccccccccccccccc */
/* c                                                                     cc */
/* c    p u r p o s e :                                                  cc */
/* c                                                                     cc */
/* c    this subroutine writes an 80-columm character header that        cc */
/* c    includes input values for the smt routine                        cc */
/* c                                                                     cc */
/* c    i n p u t :                                                      cc */
/* c                                                                     cc */
/* c        iu - i/o unit number					       cc */
/* c      ifrm - i/o format: 1) ascii, 2) binary			       cc */
/* c       var - character specifying output variable name	       cc */
/* c     value - value of output variable				       cc */
/* c   comment - character string describing value		       cc */
/* c    lenvar - length of output variable name			       cc */
/* c    lenrec - output record length				       cc */
/* c                                                                     cc */
/* c    o u t p u t :                                                    cc */
/* c                                                                     cc */
/* c     var, value, and comment					       cc */
/* c                                                                     cc */
/* cccccccccccccccccccccccccccc  h e a d  e r  cccccccccccccccccccccccccccc */


    for (i__ = 1; i__ <= 127; ++i__) {
	*(unsigned char *)&space[i__ - 1] = ' ';
	*(unsigned char *)&name0[i__ - 1] = ' ';
	*(unsigned char *)&name1[i__ - 1] = ' ';
	*(unsigned char *)&name2[i__ - 1] = ' ';
	*(unsigned char *)&name3[i__ - 1] = ' ';
/* L1001: */
    }

    lenvr = i_len(var, var_len);
    lenvl = i_len(value, value_len);
    lencm = i_len(comment, comment_len);

    charsp_(var, name0, &lvar, &lenvr, &nlb, &ntb, var_len, (ftnlen)1);
    charsp_(value, name1, &length, &lenvl, &nlb, &ntb, value_len, (ftnlen)1);
    charsp_(comment, name2, &lcom, &lencm, &nlb, &ntb, comment_len, (ftnlen)1)
	    ;

/* ****    find the number of spaces between the value and the comment */

    l0 = *lenval - length;
    if (l0 < 1) {
	l0 = 1;
	l1 = *lenrec - (length + *lenvar + 4);
	if (lcom > l1) {
	    lcom = l1;
	}
    }

    s_copy(record, " ", (ftnlen)127, (ftnlen)1);

    s_wsfi(&io___18);
    i__1 = *lenvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, name0 + (i__ - 1), (ftnlen)1);
    }
    do_fio(&c__1, "= ", (ftnlen)2);
    i__2 = length;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_fio(&c__1, name1 + (i__ - 1), (ftnlen)1);
    }
    i__3 = l0;
    for (i__ = 1; i__ <= i__3; ++i__) {
	do_fio(&c__1, space + (i__ - 1), (ftnlen)1);
    }
    do_fio(&c__1, "/ ", (ftnlen)2);
    i__4 = lcom;
    for (i__ = 1; i__ <= i__4; ++i__) {
	do_fio(&c__1, name2 + (i__ - 1), (ftnlen)1);
    }
    e_wsfi();

    charsp_(record, name3, &len3, &c__127, &nlb, &ntb, (ftnlen)127, (ftnlen)1)
	    ;

    if (*ifrm == 1) {
	io___20.ciunit = *iu;
	s_wsfe(&io___20);
	i__1 = *lenrec;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, name3 + (i__ - 1), (ftnlen)1);
	}
	e_wsfe();
    } else {
	io___21.ciunit = *iu;
	s_wsue(&io___21);
	i__1 = *lenrec;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_uio(&c__1, name3 + (i__ - 1), (ftnlen)1);
	}
	e_wsue();
    }

    return 0;
} /* header_ */

