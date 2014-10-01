/* readfil.f -- translated by f2c (version 20100827).
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

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__132 = 132;
static integer c__4 = 4;

/* Subroutine */ int readfil_(integer *i__, integer *iunit, char *filez, char 
	*formz, integer *iform, integer *nzd, integer *ncd, integer *ioff, 
	integer *iskip, integer *invz, integer *nval, integer *icol, integer *
	npt, real *scz, real *addz, real *z__, integer *ncmax, integer *ierr, 
	ftnlen filez_len, ftnlen formz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    cilist ci__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfi(icilist *), do_fio(integer *, char *, ftnlen)
	    , e_wsfi(void), s_wsfe(cilist *), e_wsfe(void), f_clos(cllist *), 
	    f_open(olist *), s_rsle(cilist *), e_rsle(void), s_rsfe(cilist *),
	     e_rsfe(void), s_cmp(char *, char *, ftnlen, ftnlen), s_rsue(
	    cilist *), e_rsue(void), do_uio(integer *, char *, ftnlen);

    /* Local variables */
    static integer j, k, l, n;
    static char fl[1];
    static integer ii, nlb, len, ntb;
    static char file[132], name__[1*132];
    static integer nptj, iflag;
    static real array[80];
    extern /* Subroutine */ int charsp_(char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal darray[80];
    static integer iarray[80];

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 6, 0, 0, 0 };
    static icilist io___8 = { 0, file, 0, "(132a)", 132, 1 };
    static cilist io___11 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___15 = { 0, 0, 1, 0, 0 };
    static cilist io___20 = { 0, 6, 0, "(2a)", 0 };
    static cilist io___21 = { 0, 6, 0, "(/,132a)", 0 };
    static cilist io___22 = { 0, 0, 1, 0, 0 };
    static cilist io___23 = { 1, 0, 1, 0, 0 };
    static cilist io___24 = { 1, 0, 1, 0, 0 };
    static cilist io___25 = { 1, 0, 1, 0, 0 };
    static cilist io___27 = { 1, 0, 1, 0, 0 };
    static cilist io___28 = { 0, 6, 0, "(/,132a)", 0 };
    static cilist io___29 = { 0, 0, 1, 0, 0 };
    static cilist io___30 = { 1, 0, 1, 0, 0 };
    static cilist io___31 = { 1, 0, 1, 0, 0 };
    static cilist io___32 = { 0, 6, 0, "(/,1a,i5,132a)", 0 };



/* ccccccccccccccccccccccccc  r e a d f i l  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine reads data records from a formatted or          cc */
/* c    unformatted file.                                               cc */
/* c                                                                    cc */
/* c         i - i/o mode -1) skip records, 0) open file, +1) read file cc */
/* c     iunit - unit number for input                                  cc */
/* c     formz - format for formatted data files                        cc */
/* c     filez - name of input file                                     cc */
/* c       nzd - maximum number of input records                        cc */
/* c       ncd - maximum number of input columns                        cc */
/* c         i - data point index                                       cc */
/* c     iunit - unit number of coordinate grid input                   cc */
/* c       npt - number of points in z array                            cc */
/* c     ncmax - number of variables to read in each record             cc */
/* c      nval - number of variables to keep in each record             cc */
/* c      icol - column index of each value to keep in each record      cc */
/* c      invz - data order inversion flag (1) yes (2) no               cc */
/* c     iform - format mode (formatted, unfor., list directed)         cc */
/* c       scz - multiplicative scaling factor for each input quantity  cc */
/* c      addz - additive offset factor for each input quantity         cc */
/* c      ioff - number of records to skip at top of file               cc */
/* c     iskip - number of records to skip between records              cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c       z(npt,icol) - value of quantity at each point                cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  r e a d f i l  ccccccccccccccccccccccccccccc */






    /* Parameter adjustments */
    z_dim1 = *nzd;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --addz;
    --scz;
    --icol;

    /* Function Body */
    *ierr = 0;
    *ncmax = icol[1];
    i__1 = *nval;
    for (n = 1; n <= i__1; ++n) {
	if (icol[n] > *ncmax) {
	    *ncmax = icol[n];
	}
/* L1001: */
    }
    if (*ncmax > 80) {
	s_wsle(&io___2);
	do_lio(&c__9, &c__1, " Number of input values on each record exceeds "
		, (ftnlen)47);
	do_lio(&c__9, &c__1, "dimension bound in readfil: ", (ftnlen)28);
	do_lio(&c__9, &c__1, "ncmax =", (ftnlen)7);
	do_lio(&c__3, &c__1, (char *)&(*ncmax), (ftnlen)sizeof(integer));
	e_wsle();
	*ierr = 1;
	return 0;
    }

/* ****   truncate file name if necessary */

    charsp_(filez, name__, &len, &c__132, &nlb, &ntb, (ftnlen)132, (ftnlen)1);
    s_wsfi(&io___8);
    i__1 = len;
    for (ii = 1; ii <= i__1; ++ii) {
	do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
    }
    e_wsfi();

    if (*iform == 1) {

/* ****     f o r m a t t e d    d a t a    f i l e */

	if (*i__ == 0) {
	    nptj = 0;
	    s_wsfe(&io___11);
	    do_fio(&c__1, " readfil reading formatted file: ", (ftnlen)33);
	    do_fio(&c__1, file, (ftnlen)132);
	    e_wsfe();
	    cl__1.cerr = 0;
	    cl__1.cunit = *iunit;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    o__1.oerr = 1;
	    o__1.ounit = *iunit;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = file;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = "formatted";
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L5203;
	    }
	    iflag = 0;

/* *****        determine if the quantites are real or integer */

	    for (l = 1; l <= 40; ++l) {
		*(unsigned char *)fl = *(unsigned char *)&formz[l - 1];
		if (*(unsigned char *)fl == 'i' || *(unsigned char *)fl == 
			'i') {
		    iflag = 1;
		}
/* L1201: */
	    }
	} else {
	    if (*i__ < 0) {

/* ****         skip over records at beginning of file */

		i__1 = *ioff;
		for (l = 1; l <= i__1; ++l) {
		    io___15.ciunit = *iunit;
		    i__2 = s_rsle(&io___15);
		    if (i__2 != 0) {
			goto L5010;
		    }
		    i__2 = e_rsle();
		    if (i__2 != 0) {
			goto L5010;
		    }
/* L2001: */
		}
	    } else {

/* ****         read in data values and scale. */

		i__1 = *npt;
		for (l = 1; l <= i__1; ++l) {
		    j = l;
		    if (*invz != 2) {
			j = *npt - l + 1;
		    }
		    if (iflag == 1) {
			ci__1.cierr = 1;
			ci__1.ciend = 1;
			ci__1.ciunit = *iunit;
			ci__1.cifmt = formz;
			i__2 = s_rsfe(&ci__1);
			if (i__2 != 0) {
			    goto L100001;
			}
			i__3 = *ncmax;
			for (k = 1; k <= i__3; ++k) {
			    i__2 = do_fio(&c__1, (char *)&iarray[k - 1], (
				    ftnlen)sizeof(integer));
			    if (i__2 != 0) {
				goto L100001;
			    }
			}
			i__2 = e_rsfe();
L100001:
			if (i__2 < 0) {
			    goto L5010;
			}
			if (i__2 > 0) {
			    goto L5110;
			}
			i__2 = *nval;
			for (k = 1; k <= i__2; ++k) {
			    z__[j + k * z_dim1] = scz[k] * (real) iarray[icol[
				    k] - 1] + addz[k];
/* L2201: */
			}
		    } else {
			ci__1.cierr = 1;
			ci__1.ciend = 1;
			ci__1.ciunit = *iunit;
			ci__1.cifmt = formz;
			i__2 = s_rsfe(&ci__1);
			if (i__2 != 0) {
			    goto L100002;
			}
			i__3 = *ncmax;
			for (k = 1; k <= i__3; ++k) {
			    i__2 = do_fio(&c__1, (char *)&array[k - 1], (
				    ftnlen)sizeof(real));
			    if (i__2 != 0) {
				goto L100002;
			    }
			}
			i__2 = e_rsfe();
L100002:
			if (i__2 < 0) {
			    goto L5010;
			}
			if (i__2 > 0) {
			    goto L5110;
			}
			i__2 = *nval;
			for (k = 1; k <= i__2; ++k) {
			    z__[j + k * z_dim1] = scz[k] * array[icol[k] - 1] 
				    + addz[k];
/* L2221: */
			}
		    }
		    i__2 = *iskip;
		    for (k = 1; k <= i__2; ++k) {
			ci__1.cierr = 1;
			ci__1.ciend = 1;
			ci__1.ciunit = *iunit;
			ci__1.cifmt = formz;
			i__3 = s_rsfe(&ci__1);
			if (i__3 != 0) {
			    goto L100003;
			}
			i__3 = e_rsfe();
L100003:
			if (i__3 < 0) {
			    goto L5010;
			}
			if (i__3 > 0) {
			    goto L5110;
			}
/* L2421: */
		    }
		    ++nptj;
/* L2441: */
		}
	    }
	}

    } else {
	if (*iform == 2) {

/* ****      read an unformatted file */

	    if (*i__ == 0) {
		nptj = 0;
		iflag = 2;
		s_wsfe(&io___20);
		do_fio(&c__1, "formz = ", (ftnlen)8);
		do_fio(&c__1, formz, (ftnlen)6);
		e_wsfe();
		if (s_cmp(formz, "integer", (ftnlen)7, (ftnlen)7) == 0) {
		    iflag = 1;
		}
		if (s_cmp(formz, "real*8", (ftnlen)6, (ftnlen)6) == 0 || 
			s_cmp(formz, "double", (ftnlen)6, (ftnlen)6) == 0) {
		    iflag = 3;
		}

		s_wsfe(&io___21);
		do_fio(&c__1, " readfil reading from unformatted file: ", (
			ftnlen)40);
		i__1 = len;
		for (ii = 1; ii <= i__1; ++ii) {
		    do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
		}
		e_wsfe();
		cl__1.cerr = 0;
		cl__1.cunit = *iunit;
		cl__1.csta = 0;
		f_clos(&cl__1);
		o__1.oerr = 1;
		o__1.ounit = *iunit;
		o__1.ofnmlen = 132;
		o__1.ofnm = file;
		o__1.orl = 0;
		o__1.osta = "old";
		o__1.oacc = 0;
		o__1.ofm = "unformatted";
		o__1.oblnk = 0;
		i__1 = f_open(&o__1);
		if (i__1 != 0) {
		    goto L5203;
		}

	    } else {
		if (*i__ < 0) {

/* ****         skip records at top of file. */

		    i__1 = *ioff;
		    for (l = 1; l <= i__1; ++l) {
			io___22.ciunit = *iunit;
			i__2 = s_rsue(&io___22);
			if (i__2 != 0) {
			    goto L5010;
			}
			i__2 = e_rsue();
			if (i__2 != 0) {
			    goto L5010;
			}
/* L3001: */
		    }
		} else {
		    i__1 = *npt;
		    for (l = 1; l <= i__1; ++l) {
			j = l;
			if (*invz != 2) {
			    j = *npt - l + 1;
			}

/* ****               read the next data point */

			if (iflag == 1) {
			    io___23.ciunit = *iunit;
			    i__2 = s_rsue(&io___23);
			    if (i__2 != 0) {
				goto L100004;
			    }
			    i__3 = *ncmax;
			    for (k = 1; k <= i__3; ++k) {
				i__2 = do_uio(&c__1, (char *)&iarray[k - 1], (
					ftnlen)sizeof(integer));
				if (i__2 != 0) {
				    goto L100004;
				}
			    }
			    i__2 = e_rsue();
L100004:
			    if (i__2 < 0) {
				goto L5010;
			    }
			    if (i__2 > 0) {
				goto L5110;
			    }
			    i__2 = *nval;
			    for (k = 1; k <= i__2; ++k) {
				z__[j + k * z_dim1] = scz[k] * (real) iarray[
					icol[k] - 1] + addz[k];
/* L3101: */
			    }
			} else {
			    if (iflag == 2) {
				io___24.ciunit = *iunit;
				i__2 = s_rsue(&io___24);
				if (i__2 != 0) {
				    goto L100005;
				}
				i__3 = *ncmax;
				for (k = 1; k <= i__3; ++k) {
				    i__2 = do_uio(&c__1, (char *)&array[k - 1]
					    , (ftnlen)sizeof(real));
				    if (i__2 != 0) {
					goto L100005;
				    }
				}
				i__2 = e_rsue();
L100005:
				if (i__2 < 0) {
				    goto L5010;
				}
				if (i__2 > 0) {
				    goto L5110;
				}
				i__2 = *nval;
				for (k = 1; k <= i__2; ++k) {
				    z__[j + k * z_dim1] = scz[k] * array[icol[
					    k] - 1] + addz[k];
/* L3201: */
				}
			    } else {
				io___25.ciunit = *iunit;
				i__2 = s_rsue(&io___25);
				if (i__2 != 0) {
				    goto L100006;
				}
				i__3 = *ncmax;
				for (k = 1; k <= i__3; ++k) {
				    i__2 = do_uio(&c__1, (char *)&darray[k - 
					    1], (ftnlen)sizeof(doublereal));
				    if (i__2 != 0) {
					goto L100006;
				    }
				}
				i__2 = e_rsue();
L100006:
				if (i__2 < 0) {
				    goto L5010;
				}
				if (i__2 > 0) {
				    goto L5110;
				}
				i__2 = *nval;
				for (k = 1; k <= i__2; ++k) {
				    z__[j + k * z_dim1] = scz[k] * darray[
					    icol[k] - 1] + addz[k];
/* L3221: */
				}
			    }
			}
			i__2 = *iskip;
			for (k = 1; k <= i__2; ++k) {
			    io___27.ciunit = *iunit;
			    i__3 = s_rsue(&io___27);
			    if (i__3 != 0) {
				goto L100007;
			    }
			    i__3 = e_rsue();
L100007:
			    if (i__3 < 0) {
				goto L5010;
			    }
			    if (i__3 > 0) {
				goto L5110;
			    }
/* L3421: */
			}
			++nptj;
/* L3441: */
		    }
		}
	    }

	} else {

/* ****     read a list-directed file */

	    if (*i__ == 0) {

/* ****         open the file */

		nptj = 0;
		s_wsfe(&io___28);
		do_fio(&c__1, " readfil reading from list directed file: ", (
			ftnlen)42);
		i__1 = len;
		for (ii = 1; ii <= i__1; ++ii) {
		    do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
		}
		e_wsfe();
		cl__1.cerr = 0;
		cl__1.cunit = *iunit;
		cl__1.csta = 0;
		f_clos(&cl__1);
		o__1.oerr = 1;
		o__1.ounit = *iunit;
		o__1.ofnmlen = 132;
		o__1.ofnm = file;
		o__1.orl = 0;
		o__1.osta = "old";
		o__1.oacc = 0;
		o__1.ofm = "formatted";
		o__1.oblnk = 0;
		i__1 = f_open(&o__1);
		if (i__1 != 0) {
		    goto L5203;
		}

	    } else {
		if (*i__ < 0) {

/* ****         skip records at top of file. */

		    i__1 = *ioff;
		    for (l = 1; l <= i__1; ++l) {
			io___29.ciunit = *iunit;
			i__2 = s_rsle(&io___29);
			if (i__2 != 0) {
			    goto L5010;
			}
			i__2 = e_rsle();
			if (i__2 != 0) {
			    goto L5010;
			}
/* L4001: */
		    }
		} else {
		    i__1 = *npt;
		    for (l = 1; l <= i__1; ++l) {
			j = l;
			if (*invz != 2) {
			    j = *npt - l + 1;
			}
			io___30.ciunit = *iunit;
			i__2 = s_rsle(&io___30);
			if (i__2 != 0) {
			    goto L100008;
			}
			i__3 = *ncmax;
			for (k = 1; k <= i__3; ++k) {
			    i__2 = do_lio(&c__4, &c__1, (char *)&array[k - 1],
				     (ftnlen)sizeof(real));
			    if (i__2 != 0) {
				goto L100008;
			    }
			}
			i__2 = e_rsle();
L100008:
			if (i__2 < 0) {
			    goto L5010;
			}
			if (i__2 > 0) {
			    goto L5110;
			}
			i__2 = *nval;
			for (k = 1; k <= i__2; ++k) {
			    z__[j + k * z_dim1] = scz[k] * array[icol[k] - 1] 
				    + addz[k];
/* L4201: */
			}
			i__2 = *iskip;
			for (k = 1; k <= i__2; ++k) {
			    io___31.ciunit = *iunit;
			    i__3 = s_rsle(&io___31);
			    if (i__3 != 0) {
				goto L100009;
			    }
			    i__3 = e_rsle();
L100009:
			    if (i__3 < 0) {
				goto L5010;
			    }
			    if (i__3 > 0) {
				goto L5110;
			    }
/* L4421: */
			}
			++nptj;
/* L4441: */
		    }
		}
	    }
	}
    }

    return 0;

/* ****   end of file encountered */

L5010:
    *ierr = 1;
    *npt = nptj;
    cl__1.cerr = 0;
    cl__1.cunit = *iunit;
    cl__1.csta = 0;
    f_clos(&cl__1);

    return 0;

/* ****    error encountered in file */

L5110:
    *ierr = -1;
    cl__1.cerr = 0;
    cl__1.cunit = *iunit;
    cl__1.csta = 0;
    f_clos(&cl__1);

    return 0;

/* ****  could not open file */

L5203:
    *ierr = -2;
    s_wsfe(&io___32);
    do_fio(&c__1, "I/0 Error - could not open unit ", (ftnlen)32);
    do_fio(&c__1, (char *)&(*iunit), (ftnlen)sizeof(integer));
    do_fio(&c__1, " for file: ", (ftnlen)11);
    i__1 = len;
    for (ii = 1; ii <= i__1; ++ii) {
	do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = *iunit;
    cl__1.csta = 0;
    f_clos(&cl__1);

    return 0;
} /* readfil_ */

