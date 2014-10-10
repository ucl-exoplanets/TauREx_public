/* init_pd_iu.f -- translated by f2c (version 20100827).
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

static integer c__132 = 132;
static integer c__9 = 9;
static integer c__1 = 1;

/* Subroutine */ int init_pd_iu__(integer *nza, integer *numu, integer *nphi, 
	integer *nlout, integer *ifrmout, integer *nza_1__, char *radfile, 
	logical *lplanck, logical *lsolar, char *clev, integer *nstate, 
	integer *istate, integer *iu_pd__, char *j_ext__, integer *iutpd, 
	integer *iuspd, integer *iupdrad, ftnlen radfile_len, ftnlen clev_len,
	 ftnlen j_ext_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    olist o__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , f_open(olist *), s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer j, n, nl, ir, nz;
    static char nm1[1*132];
    static integer nlb, len, npd, ntb, naz, nze;
    static char name__[1*132];
    static integer lenf;
    static char sflx[11];
    static integer ipdct, lenpd;
    static char pdfile[132];
    extern /* Subroutine */ int charsp_(char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char pd_ext__[20], pd_name__[1*9];

    /* Fortran I/O blocks */
    static icilist io___12 = { 0, pdfile, 0, "(132a)", 132, 1 };
    static cilist io___16 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static icilist io___18 = { 0, sflx, 0, "(1a9,i2.2)", 11, 1 };
    static icilist io___19 = { 0, pdfile, 0, "(132a)", 132, 1 };
    static cilist io___20 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static icilist io___26 = { 0, pd_ext__, 0, "(1a4,i3.3,1a4,9a)", 20, 1 };
    static icilist io___27 = { 0, pdfile, 0, "(132a)", 132, 1 };
    static cilist io___28 = { 0, 6, 0, "(i5,a16,132a)", 0 };



/* ccccccccccccccccccccccc  i n i t _ p d _ i u   cccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    This subroutine initializes the output units for radiance and   cc */
/* c    and flux jacobians.                                             cc */
/* c                                                                    cc */
/* c    i n p u t:                                                      cc */
/* c                                                                    cc */
/* c        nza - number of solar zenith angles                         cc */
/* c       numu - number of output zenith angles                        cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c        npd - number of variable elements of the state vector       cc */
/* c      nlout - number of levels where radiance spectra are output    cc */
/* c     ifmout - index of output file format (1) ascii, (2) binary,    cc */
/* c                                          (3) binary, no header     cc */
/* c      nza_1 - index of first stream to be printed out               cc */
/* c       name - string vector with parsed name of output file         cc */
/* c        len - length of name (no leading or trailing blanks)        cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c     nstate - number of elements in the state vector                cc */
/* c     istate - state vector flag indicating which state variables    cc */
/* c              are variable components of the state vector.          cc */
/* c              0 - not a variable component of the state vector      cc */
/* c              1 - surface pressure                                  cc */
/* c              2 - surface/atmospheric temperature                   cc */
/* c              3 - gas absorption coeffient                          cc */
/* c              4 - cloud/aerosol optical depth                       cc */
/* c              5 - surface albedo                                    cc */
/* c      iu_pd - starting index of partial derivative output units     cc */
/* c       clev - character vector with out level index                 cc */
/* c      j_ext - character varible with state vector type (istate)     cc */
/* c                                                                    cc */
/* c    o u t p u t:                                                    cc */
/* c                                                                    cc */
/* c       iutpd - unit number for each thermal flux */
/* c                                                                    cc */
/* ccccccccccccccccccccccc  i n i t _ p d _ i u   cccccccccccccccccccccccc */




/* ccccccccccccccccccccccccccccc  p a r a m  ccccccccccccccccccccccccccccc */
/* c								      cc */
/* c    p u r p o s e :                                                 cc */
/* c								      cc */
/* c    this include file specifies the dimensions of variables used    cc */
/* c    in the program smt_do.                                          cc */
/* c								      cc */
/* c    q u a n t i t i e s :                                           cc */
/* c								      cc */
/* c      nsp - Maximum number of spectral intervals 		      cc */
/* c       kp - Maximum number of vertical levels.		      cc */
/* c    nmode - Maximum number of aerosol particle modes.               cc */
/* c     ngas - Maximum number of absorbing gases.		      cc */
/* c     ngrp - Maximum number of SMT grouping bins.		      cc */
/* c    mxumu - Maximum number of user output zenith angles in discrete cc */
/* c            ordinate model.					      cc */
/* c    mxphi - Maximum number of output azimuth angles in the discrete cc */
/* c            ordinate model.					      cc */
/* c    mxmom - Maximum number of legendre polynomial moments in the    cc */
/* c            discrete ordinate method.				      cc */
/* c    mxulv - Maximum number of output levels in the discrete 	      cc */
/* c            ordinate model.					      cc */
/* c    mxrad - Maximum number of output radiance levels.		      cc */
/* c     nsol - Maximum number of solar zenith angles.		      cc */
/* c      nex - Maximum number of spectrally-dependent optical          cc */
/* c            (gases, aerosol modes, surface albedos, and solar flux. cc */
/* c     mxpd - maximum number of radiance partial derivatives          cc */
/* c								      cc */
/* ccccccccccccccccccccccccccccc  p a r a m  ccccccccccccccccccccccccccccc */




/* ****    units for flux and radiance jacobians */


/* ****   local integer variables */


/* ****       create an output files for each partial derivative */

    /* Parameter adjustments */
    iupdrad -= 10513;
    iuspd -= 11;
    iutpd -= 11;
    j_ext__ -= 9;
    --istate;
    clev -= 4;
    radfile -= 528;
    --nza_1__;

    /* Function Body */
    ipdct = -1;

/* ****   enter loop over solar zenith angle */

    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {

	npd = 0;

/* ****       enter loop over partial derivative */

	i__2 = *nstate;
	for (n = 1; n <= i__2; ++n) {

/* ****            flux partial derivtives are only calculated */
/*                if istate is less than zero */

	    if (istate[n] < 0) {

		charsp_(radfile + (nz * 3 + 1) * 132, name__, &len, &c__132, &
			nlb, &ntb, (ftnlen)132, (ftnlen)1);

		++npd;

		charsp_(j_ext__ + n * 9, pd_name__, &lenpd, &c__9, &nlb, &ntb,
			 (ftnlen)9, (ftnlen)1);

/* ****             increment the unit counter.  This counter has to */
/*                 be incremented even if thermal fluxes aren't found */
/*                 in this run to support the separate solar and thermal */
/*                 flux calculations in vpl_climate */

		if (nz == 1) {
		    ++ipdct;
		    if (*lplanck) {

/* ****                 open units for partial derivatives for the */
/*                     diffuse thermal flux */

			s_wsfi(&io___12);
			i__3 = len;
			for (j = 1; j <= i__3; ++j) {
			    do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
			}
			do_fio(&c__1, "_tflx", (ftnlen)5);
			i__4 = lenpd;
			for (j = 1; j <= i__4; ++j) {
			    do_fio(&c__1, pd_name__ + (j - 1), (ftnlen)1);
			}
			e_wsfi();

			iutpd[npd + nz * 10] = *iu_pd__ + ipdct;
			if (*ifrmout == 1) {
			    o__1.oerr = 0;
			    o__1.ounit = iutpd[npd + nz * 10];
			    o__1.ofnmlen = 132;
			    o__1.ofnm = pdfile;
			    o__1.orl = 0;
			    o__1.osta = "unknown";
			    o__1.oacc = 0;
			    o__1.ofm = "formatted";
			    o__1.oblnk = 0;
			    f_open(&o__1);
			} else {
			    o__1.oerr = 0;
			    o__1.ounit = iutpd[npd + nz * 10];
			    o__1.ofnmlen = 132;
			    o__1.ofnm = pdfile;
			    o__1.orl = 0;
			    o__1.osta = "unknown";
			    o__1.oacc = 0;
			    o__1.ofm = "unformatted";
			    o__1.oblnk = 0;
			    f_open(&o__1);
			}

			charsp_(pdfile, nm1, &lenf, &c__132, &nlb, &ntb, (
				ftnlen)132, (ftnlen)1);

			s_wsfe(&io___16);
			do_fio(&c__1, (char *)&iutpd[npd + nz * 10], (ftnlen)
				sizeof(integer));
			do_fio(&c__1, " Jacobian  ", (ftnlen)11);
			i__3 = lenf;
			for (j = 1; j <= i__3; ++j) {
			    do_fio(&c__1, nm1 + (j - 1), (ftnlen)1);
			}
			e_wsfe();

		    }
		}

/* ****             increment the unit counter.  This counter has to */
/*                 be incremented even if solar fluxes aren't found */
/*                 in this run to support the separate solar and thermal */
/*                 flux calculations in vpl_climate */

		++ipdct;
		if (*lsolar) {

/* ****               open unit for partial derivatives for downward */
/*                   solar flux */

		    s_wsfi(&io___18);
		    do_fio(&c__1, "_sflx_sza", (ftnlen)9);
		    do_fio(&c__1, (char *)&nz, (ftnlen)sizeof(integer));
		    e_wsfi();
		    s_wsfi(&io___19);
		    i__3 = len;
		    for (j = 1; j <= i__3; ++j) {
			do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
		    }
		    do_fio(&c__1, sflx, (ftnlen)11);
		    i__4 = lenpd;
		    for (j = 1; j <= i__4; ++j) {
			do_fio(&c__1, pd_name__ + (j - 1), (ftnlen)1);
		    }
		    e_wsfi();

		    iuspd[npd + nz * 10] = *iu_pd__ + ipdct;
		    if (*ifrmout == 1) {
			o__1.oerr = 0;
			o__1.ounit = iuspd[npd + nz * 10];
			o__1.ofnmlen = 132;
			o__1.ofnm = pdfile;
			o__1.orl = 0;
			o__1.osta = "unknown";
			o__1.oacc = 0;
			o__1.ofm = "formatted";
			o__1.oblnk = 0;
			f_open(&o__1);
		    } else {
			o__1.oerr = 0;
			o__1.ounit = iuspd[npd + nz * 10];
			o__1.ofnmlen = 132;
			o__1.ofnm = pdfile;
			o__1.orl = 0;
			o__1.osta = "unknown";
			o__1.oacc = 0;
			o__1.ofm = "unformatted";
			o__1.oblnk = 0;
			f_open(&o__1);
		    }


		    charsp_(pdfile, nm1, &lenf, &c__132, &nlb, &ntb, (ftnlen)
			    132, (ftnlen)1);

		    s_wsfe(&io___20);
		    do_fio(&c__1, (char *)&iuspd[npd + nz * 10], (ftnlen)
			    sizeof(integer));
		    do_fio(&c__1, " Jacobian  ", (ftnlen)11);
		    i__3 = lenf;
		    for (j = 1; j <= i__3; ++j) {
			do_fio(&c__1, nm1 + (j - 1), (ftnlen)1);
		    }
		    e_wsfe();

		}

	    }

	    if (istate[n] > 0) {

		++npd;

		charsp_(j_ext__ + n * 9, pd_name__, &lenpd, &c__9, &nlb, &ntb,
			 (ftnlen)9, (ftnlen)1);

/* ****            r a d i a n c e   p a r t i a l   d e r i v a t i v e s */

/* ***             enter loop over output level */

		i__3 = *nlout;
		for (nl = 1; nl <= i__3; ++nl) {

		    charsp_(radfile + (nl + nz * 3) * 132, name__, &len, &
			    c__132, &nlb, &ntb, (ftnlen)132, (ftnlen)1);

/* ****                 enter loops over emission azimuth and zenith angle */

		    ir = 0;
		    i__4 = *nphi;
		    for (naz = 1; naz <= i__4; ++naz) {

			i__5 = *numu;
			for (nze = nza_1__[nl]; nze <= i__5; ++nze) {

			    ++ir;
			    ++ipdct;
			    s_wsfi(&io___26);
			    do_fio(&c__1, "_rad", (ftnlen)4);
			    do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer)
				    );
			    do_fio(&c__1, clev + (nl << 2), (ftnlen)4);
			    i__6 = lenpd;
			    for (j = 1; j <= i__6; ++j) {
				do_fio(&c__1, pd_name__ + (j - 1), (ftnlen)1);
			    }
			    e_wsfi();

/* ****                          compute the unit number */

			    iupdrad[nze + (naz + (npd + (nl + nz * 3) * 10 << 
				    4) << 4)] = *iu_pd__ + ipdct;

			    s_wsfi(&io___27);
			    i__6 = len;
			    for (j = 1; j <= i__6; ++j) {
				do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
			    }
			    do_fio(&c__1, pd_ext__, (ftnlen)20);
			    e_wsfi();
			    if (*ifrmout == 1) {
				o__1.oerr = 0;
				o__1.ounit = iupdrad[nze + (naz + (npd + (nl 
					+ nz * 3) * 10 << 4) << 4)];
				o__1.ofnmlen = 132;
				o__1.ofnm = pdfile;
				o__1.orl = 0;
				o__1.osta = "unknown";
				o__1.oacc = 0;
				o__1.ofm = "formatted";
				o__1.oblnk = 0;
				f_open(&o__1);
			    } else {
				o__1.oerr = 0;
				o__1.ounit = iupdrad[nze + (naz + (npd + (nl 
					+ nz * 3) * 10 << 4) << 4)];
				o__1.ofnmlen = 132;
				o__1.ofnm = pdfile;
				o__1.orl = 0;
				o__1.osta = "unknown";
				o__1.oacc = 0;
				o__1.ofm = "unformatted";
				o__1.oblnk = 0;
				f_open(&o__1);
			    }

			    charsp_(pdfile, nm1, &lenf, &c__132, &nlb, &ntb, (
				    ftnlen)132, (ftnlen)1);

			    s_wsfe(&io___28);
			    do_fio(&c__1, (char *)&iupdrad[nze + (naz + (npd 
				    + (nl + nz * 3) * 10 << 4) << 4)], (
				    ftnlen)sizeof(integer));
			    do_fio(&c__1, " Jacobian  ", (ftnlen)11);
			    i__6 = lenf;
			    for (j = 1; j <= i__6; ++j) {
				do_fio(&c__1, nm1 + (j - 1), (ftnlen)1);
			    }
			    e_wsfe();
/* L1201: */
			}
/* L1221: */
		    }
/* L1241: */
		}

	    }

/* L1261: */
	}
/* L1281: */
    }

    return 0;
} /* init_pd_iu__ */

