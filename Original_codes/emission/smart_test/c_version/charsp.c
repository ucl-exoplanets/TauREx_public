/* charsp.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int charsp_(char *charin, char *chout, integer *length, 
	integer *nl, integer *nlb, integer *ntb, ftnlen charin_len, ftnlen 
	chout_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, l1, ii;


/* cccccccccccccccccccccccccc  c h a r s p  ccccccccccccccccccccccccccccc */
/* c                                                                   cc */
/* c    p u r p o s e :                                                cc */
/* c                                                                   cc */
/* c    this subroutine takes a character of arbitrary length and      cc */
/* c    splits it into an array of characters with unit length.  the   cc */
/* c    number of leading and trailing blanks and the length of the    cc */
/* c    input character, less these spaces is then found.              cc */
/* c                                                                   cc */
/* c    i n p u t :                                                    cc */
/* c                                                                   cc */
/* c    charin : input character of length nl                          cc */
/* c        nl : length of the input character and dimension of chout  cc */
/* c                                                                   cc */
/* c    o u t p u t :                                                  cc */
/* c                                                                   cc */
/* c     chout : output character array containing charin characters   cc */
/* c             stripped of leading and trailing blanks               cc */
/* c    length : number of elements in charin less leading and         cc */
/* c             trailing blanks                                       cc */
/* c       nlb : number of leading blanks in charin                    cc */
/* c       ntb : number of trailing blanks in charin                   cc */
/* c                                                                   cc */
/* cccccccccccccccccccccccccc  c h a r s p  ccccccccccccccccccccccccccccc */


/* ****  find the number of leading blanks */

    /* Parameter adjustments */
    --chout;

    /* Function Body */
    *nlb = 0;
    i__1 = *nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&charin[i__ - 1] != ' ' && *(unsigned char *)&
		charin[i__ - 1] != 0) {
	    goto L1201;
	}
	++(*nlb);
/* L1001: */
    }

/* ****  find the number of trailing blanks */

L1201:
    *ntb = 0;
    i__1 = *nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *nl - i__ + 1;
	if (*(unsigned char *)&charin[ii - 1] != ' ' && *(unsigned char *)&
		charin[ii - 1] != 0) {
	    goto L1401;
	}
	++(*ntb);
/* L1221: */
    }

/* ****  find the length of charin less leading and trailing blanks */

L1401:
    *length = *nl - *nlb - *ntb;
    l1 = *length;
/*      do 1421 l=1,nl */
/*          if(charin(l:l) .eq. '|' .or. charin(l:l) .eq. '~') */
/*     -    length = length - 1 */
/* 1421  continue */

/* ****  compose output character array chout */

    i__1 = l1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = i__ + *nlb;
	*(unsigned char *)&chout[i__] = *(unsigned char *)&charin[ii - 1];
/* L1601: */
    }

    return 0;
} /* charsp_ */

