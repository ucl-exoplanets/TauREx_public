/* init_spect_io.f -- translated by f2c (version 20100827).
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
static integer c__132 = 132;
static integer c__3 = 3;

/* Subroutine */ int init_spect_io__(integer *nza, integer *nlout, integer *
	irad, integer *ifrmout, real *sza0, char *clev, integer *iuout, 
	integer *iuflx, integer *iuheat, integer *iustat, integer *iutrn, 
	char *name__, integer *len, char *radfile, char *heatfile, char *
	statfile, char *trnfile, char *flxfile, ftnlen clev_len, ftnlen 
	name_len, ftnlen radfile_len, ftnlen heatfile_len, ftnlen 
	statfile_len, ftnlen trnfile_len, ftnlen flxfile_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    icilist ici__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_rsfe(cilist *), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_nint(real *), s_wsfi(icilist *), e_wsfi(void), f_open(olist *), 
	    s_rsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer j, ic, nl, nz, nlb, iuf, iza, ntb, iuo, icmx, iover;
    static char filein[132];
    extern /* Subroutine */ int charsp_(char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char flxext[4], outext[4];

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___4 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___5 = { 0, 5, 0, "(1a)", 0 };
    static icilist io___13 = { 0, outext, 0, "(1a2,i2.2)", 4, 1 };
    static icilist io___16 = { 0, flxext, 0, "(1a2,i2.2)", 4, 1 };
    static cilist io___17 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___18 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___19 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___21 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___23 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___24 = { 0, 6, 0, "(3(/,1a))", 0 };
    static cilist io___25 = { 0, 5, 0, 0, 0 };
    static cilist io___27 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___28 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___29 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___30 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___31 = { 0, 6, 0, "(i5,a16,132a)", 0 };
    static cilist io___32 = { 0, 6, 0, "(i5,a16,132a)", 0 };



/* ccccccccccccccccccccc  i n i t _ s p e c t _ i o  ccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    This subroutine initiazies output file names for program smart  cc */
/* c                                                                    cc */
/* c    i n p u t:                                                      cc */
/* c                                                                    cc */
/* c        nza - number of solar zenith angles                         cc */
/* c      nlout - number of levels where radiance spectra are output    cc */
/* c       irad - index of output file type:                            cc */
/* c              1) wavelength-dependent fluxes and radiances at the   cc */
/* c                 computational azimuths and zenith angles at the    cc */
/* c                 specified output levels, and spectrally-integrated cc */
/* c                 fluxes and heating rates at each computational     cc */
/* c                 level.                                             cc */
/* c              2) wavelength-dependent fluxes, radiances, and        cc */
/* c                 transmission values at computational zenith angles cc */
/* c                 and specified output levels, and spectrally-       cc */
/* c                 integrated fluxes and heating rates at each        cc */
/* c                 computational level.                               cc */
/* c              3) wavelength-dependent fluxes and radiances at the   cc */
/* c                 computational azimuths and zenith angles at the    cc */
/* c                 specified output levels, wavelength-dependent,     cc */
/* c                 level-dependent pressure-weighted, flux            cc */
/* c                 divergences, and spectrally integrated fluxes      cc */
/* c                 and heating rates at each computational level.     cc */
/* c              4) wavelength-dependent fluxes, radiances,            cc */
/* c                 transmission values, wavelength-dependent,         cc */
/* c                 level-dependent pressure-weighted, flux            cc */
/* c                 divergences, and spectrally-integrated fluxes      cc */
/* c                 and heating rates at each computational level.     cc */
/* c              5) fluxes, radiances, and heating rates,              cc */
/* c                 at arbitrary azimuths and zenith angles,           cc */
/* c              6) fluxes, radiances and transmission functions       cc */
/* c                 at arbitrary zenith angles,                        cc */
/* c              7) fluxes, radiances, and contribution functions      cc */
/* c                 at arbitrary zenith angles,                        cc */
/* c              8) fluxes, radiances, transmission functions,and      cc */
/* c                 contribution functions at arbitrary zenith angles. cc */
/* c     ifmout - index of output file format (1) ascii, (2) binary,    cc */
/* c                                          (3) binary, no header     cc */
/* c       sza0 - solar zenith angles (degrees)                         cc */
/* c       clev - string variable used to indicate output level         cc */
/* c      iuout - unit number for output radiance file                  cc */
/* c      iutrn - unit number for output transmission/pressure file     cc */
/* c      iuflx - unit number for output level-dependent fluxes         cc */
/* c     iuheat - unit number for output solar heating rates            cc */
/* c     iustat - unit number for output binning statistics             cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      name - character string with parsed output file name          cc */
/* c       len - length of output file name (stripped of blanks)        cc */
/* c    radfile - name output file for flux/radiance spectra            cc */
/* c   statfile - name of output file with binning statistics           cc */
/* c   heatfile - name of output file with heating/cooling rates        cc */
/* c    trnfile - name of output file with tau=1 trans/pressure         cc */
/* c    flxfile - name of output file with level-dependent flux spectra cc */
/* c                                                                    cc */
/* ccccccccccccccccccccc  i n i t _ s p e c t _ i o  ccccccccccccccccccccc */




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






    /* Parameter adjustments */
    flxfile -= 132;
    radfile -= 528;
    --name__;
    clev -= 4;
    --sza0;

    /* Function Body */
    s_wsfe(&io___1);
    do_fio(&c__1, " Output Files: ", (ftnlen)15);
    do_fio(&c__1, "  Unit    File Type  Filename", (ftnlen)29);
    e_wsfe();

    icmx = 10;
    ic = 0;
L7021:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___4);
    do_fio(&c__1, " Enter prefix of output radiance and heating rate file na"
	    "mes: ", (ftnlen)62);
    e_wsfe();
    s_rsfe(&io___5);
    do_fio(&c__1, filein, (ftnlen)132);
    e_rsfe();

/* ****   parse the input file name */

    charsp_(filein, name__ + 1, len, &c__132, &nlb, &ntb, (ftnlen)132, (
	    ftnlen)1);

    i__1 = *nlout;
    for (nl = 1; nl <= i__1; ++nl) {
	if (*nza > 1) {
	    i__2 = *nza;
	    for (nz = 1; nz <= i__2; ++nz) {
		s_copy(radfile + (nl + nz * 3) * 132, " ", (ftnlen)132, (
			ftnlen)1);
		iza = i_nint(&sza0[nz]);
		s_wsfi(&io___13);
		do_fio(&c__1, ".r", (ftnlen)2);
		do_fio(&c__1, (char *)&iza, (ftnlen)sizeof(integer));
		e_wsfi();
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 132;
		ici__1.iciunit = radfile + (nl + nz * 3) * 132;
		ici__1.icifmt = "(132a)";
		s_wsfi(&ici__1);
		i__3 = *len;
		for (j = 1; j <= i__3; ++j) {
		    do_fio(&c__1, name__ + j, (ftnlen)1);
		}
		do_fio(&c__1, clev + (nl << 2), (ftnlen)4);
		do_fio(&c__1, outext, (ftnlen)4);
		e_wsfi();
/*                write(*,'(1x,2a)') 'radfile =  ',radfile(nl,nz) */
/* L7201: */
	    }
	} else {
	    s_copy(radfile + (nl + 3) * 132, " ", (ftnlen)132, (ftnlen)1);
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 132;
	    ici__1.iciunit = radfile + (nl + 3) * 132;
	    ici__1.icifmt = "(132a)";
	    s_wsfi(&ici__1);
	    i__2 = *len;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, name__ + j, (ftnlen)1);
	    }
	    do_fio(&c__1, clev + (nl << 2), (ftnlen)4);
	    do_fio(&c__1, ".rad", (ftnlen)4);
	    e_wsfi();
/*            write(*,'(1x,2a)') 'radfile = ',radfile(nl,1) */
	}
/* L7221: */
    }

    if (*irad == 3 || *irad == 4 || *irad == 7 || *irad == 8) {
	if (*nza > 1) {

/* ****       create an output flux file name for each zenith angle */

	    i__1 = *nza;
	    for (nz = 1; nz <= i__1; ++nz) {
		iza = i_nint(&sza0[nz]);
		s_wsfi(&io___16);
		do_fio(&c__1, ".f", (ftnlen)2);
		do_fio(&c__1, (char *)&iza, (ftnlen)sizeof(integer));
		e_wsfi();
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 132;
		ici__1.iciunit = flxfile + nz * 132;
		ici__1.icifmt = "(132a)";
		s_wsfi(&ici__1);
		i__2 = *len;
		for (j = 1; j <= i__2; ++j) {
		    do_fio(&c__1, name__ + j, (ftnlen)1);
		}
		do_fio(&c__1, flxext, (ftnlen)4);
		e_wsfi();
/* L7241: */
	    }
	} else {

	    s_copy(flxfile + 132, " ", (ftnlen)132, (ftnlen)1);
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 132;
	    ici__1.iciunit = flxfile + 132;
	    ici__1.icifmt = "(132a)";
	    s_wsfi(&ici__1);
	    i__1 = *len;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, name__ + j, (ftnlen)1);
	    }
	    do_fio(&c__1, ".flx", (ftnlen)4);
	    e_wsfi();
	}
    }

    s_copy(heatfile, " ", (ftnlen)132, (ftnlen)1);
    ici__1.icierr = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = 132;
    ici__1.iciunit = heatfile;
    ici__1.icifmt = "(132a)";
    s_wsfi(&ici__1);
    i__1 = *len;
    for (j = 1; j <= i__1; ++j) {
	do_fio(&c__1, name__ + j, (ftnlen)1);
    }
    do_fio(&c__1, ".hrt", (ftnlen)4);
    e_wsfi();

    s_copy(statfile, " ", (ftnlen)132, (ftnlen)1);
    ici__1.icierr = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = 132;
    ici__1.iciunit = statfile;
    ici__1.icifmt = "(132a)";
    s_wsfi(&ici__1);
    i__1 = *len;
    for (j = 1; j <= i__1; ++j) {
	do_fio(&c__1, name__ + j, (ftnlen)1);
    }
    do_fio(&c__1, ".stat", (ftnlen)5);
    e_wsfi();
    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
	s_copy(trnfile, " ", (ftnlen)132, (ftnlen)1);
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 132;
	ici__1.iciunit = trnfile;
	ici__1.icifmt = "(132a)";
	s_wsfi(&ici__1);
	i__1 = *len;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, name__ + j, (ftnlen)1);
	}
	do_fio(&c__1, ".trn", (ftnlen)4);
	e_wsfi();
    }

/* ****    open units for output fluxes and radiances.  For */
/*        level-dependent fluxes, a separate unit is needed for each */
/*        output solar zenith angle.  For radiances, a separate unit */
/*        is needed for each output level and solar zenith angle. */

    o__1.oerr = 1;
    o__1.ounit = *iuheat;
    o__1.ofnmlen = 132;
    o__1.ofnm = heatfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = "formatted";
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L7408;
    }

    charsp_(heatfile, name__ + 1, len, &c__132, &nlb, &ntb, (ftnlen)132, (
	    ftnlen)1);

    s_wsfe(&io___17);
    do_fio(&c__1, (char *)&(*iuheat), (ftnlen)sizeof(integer));
    do_fio(&c__1, " Heating Rates  ", (ftnlen)16);
    i__1 = *len;
    for (j = 1; j <= i__1; ++j) {
	do_fio(&c__1, name__ + j, (ftnlen)1);
    }
    e_wsfe();

    o__1.oerr = 1;
    o__1.ounit = *iustat;
    o__1.ofnmlen = 132;
    o__1.ofnm = statfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = "formatted";
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L7408;
    }

    charsp_(statfile, name__ + 1, len, &c__132, &nlb, &ntb, (ftnlen)132, (
	    ftnlen)1);

    s_wsfe(&io___18);
    do_fio(&c__1, (char *)&(*iustat), (ftnlen)sizeof(integer));
    do_fio(&c__1, " Statistics  ", (ftnlen)13);
    i__1 = *len;
    for (j = 1; j <= i__1; ++j) {
	do_fio(&c__1, name__ + j, (ftnlen)1);
    }
    e_wsfe();

    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
	if (*ifrmout == 1) {
	    o__1.oerr = 1;
	    o__1.ounit = *iutrn;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = trnfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = "formatted";
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L7408;
	    }
	} else {
	    o__1.oerr = 1;
	    o__1.ounit = *iutrn;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = trnfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = "unformatted";
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L7408;
	    }
	}

	charsp_(trnfile, name__ + 1, len, &c__132, &nlb, &ntb, (ftnlen)132, (
		ftnlen)1);

	s_wsfe(&io___19);
	do_fio(&c__1, (char *)&(*iutrn), (ftnlen)sizeof(integer));
	do_fio(&c__1, " Transmission  ", (ftnlen)15);
	i__1 = *len;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, name__ + j, (ftnlen)1);
	}
	e_wsfe();
    }

    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {

	if (*irad == 3 || *irad == 4 || *irad == 7 || *irad == 8) {

/* ****          open unit for output fluxes */

	    iuf = *iuflx + nz - 1;
	    if (*ifrmout == 1) {
		o__1.oerr = 1;
		o__1.ounit = iuf;
		o__1.ofnmlen = 132;
		o__1.ofnm = flxfile + nz * 132;
		o__1.orl = 0;
		o__1.osta = "new";
		o__1.oacc = 0;
		o__1.ofm = "formatted";
		o__1.oblnk = 0;
		i__2 = f_open(&o__1);
		if (i__2 != 0) {
		    goto L7408;
		}
	    } else {
		o__1.oerr = 1;
		o__1.ounit = iuf;
		o__1.ofnmlen = 132;
		o__1.ofnm = flxfile + nz * 132;
		o__1.orl = 0;
		o__1.osta = "new";
		o__1.oacc = 0;
		o__1.ofm = "unformatted";
		o__1.oblnk = 0;
		i__2 = f_open(&o__1);
		if (i__2 != 0) {
		    goto L7408;
		}
	    }

	    charsp_(flxfile + nz * 132, name__ + 1, len, &c__132, &nlb, &ntb, 
		    (ftnlen)132, (ftnlen)1);

	    s_wsfe(&io___21);
	    do_fio(&c__1, (char *)&iuf, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " Flux  ", (ftnlen)7);
	    i__2 = *len;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, name__ + j, (ftnlen)1);
	    }
	    e_wsfe();
	}

/* ****       open output radiance file: */

	i__2 = *nlout;
	for (nl = 1; nl <= i__2; ++nl) {
	    iuo = *iuout + (nz - 1) * *nlout + nl - 1;
	    if (*ifrmout == 1) {
		o__1.oerr = 1;
		o__1.ounit = iuo;
		o__1.ofnmlen = 132;
		o__1.ofnm = radfile + (nl + nz * 3) * 132;
		o__1.orl = 0;
		o__1.osta = "new";
		o__1.oacc = 0;
		o__1.ofm = "formatted";
		o__1.oblnk = 0;
		i__3 = f_open(&o__1);
		if (i__3 != 0) {
		    goto L7408;
		}
	    } else {
		o__1.oerr = 1;
		o__1.ounit = iuo;
		o__1.ofnmlen = 132;
		o__1.ofnm = radfile + (nl + nz * 3) * 132;
		o__1.orl = 0;
		o__1.osta = "new";
		o__1.oacc = 0;
		o__1.ofm = "unformatted";
		o__1.oblnk = 0;
		i__3 = f_open(&o__1);
		if (i__3 != 0) {
		    goto L7408;
		}
	    }

	    charsp_(radfile + (nl + nz * 3) * 132, name__ + 1, len, &c__132, &
		    nlb, &ntb, (ftnlen)132, (ftnlen)1);

	    s_wsfe(&io___23);
	    do_fio(&c__1, (char *)&iuo, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " Radiances  ", (ftnlen)12);
	    i__3 = *len;
	    for (j = 1; j <= i__3; ++j) {
		do_fio(&c__1, name__ + j, (ftnlen)1);
	    }
	    e_wsfe();
/* L7301: */
	}
/* L7321: */
    }
    goto L7601;

L7408:
    s_wsfe(&io___24);
    do_fio(&c__1, " Output file already exists.  Choose option:", (ftnlen)44);
    do_fio(&c__1, " 1) enter a new file name", (ftnlen)25);
    do_fio(&c__1, " 2) overwrite existing file", (ftnlen)27);
    e_wsfe();
    s_rsle(&io___25);
    do_lio(&c__3, &c__1, (char *)&iover, (ftnlen)sizeof(integer));
    e_rsle();
    s_wsfe(&io___27);
    do_fio(&c__1, "iover =", (ftnlen)7);
    do_fio(&c__1, (char *)&iover, (ftnlen)sizeof(integer));
    e_wsfe();
    if (iover == 1) {
	goto L7021;
    }

/* ****   overwrite existing file */

L7601:
    cl__1.cerr = 0;
    cl__1.cunit = *iuheat;
    cl__1.csta = 0;
    f_clos(&cl__1);
    o__1.oerr = 0;
    o__1.ounit = *iuheat;
    o__1.ofnmlen = 132;
    o__1.ofnm = heatfile;
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = "formatted";
    o__1.oblnk = 0;
    f_open(&o__1);

    charsp_(heatfile, name__ + 1, len, &c__132, &nlb, &ntb, (ftnlen)132, (
	    ftnlen)1);

    s_wsfe(&io___28);
    do_fio(&c__1, (char *)&(*iuheat), (ftnlen)sizeof(integer));
    do_fio(&c__1, " Heating Rates  ", (ftnlen)16);
    i__1 = *len;
    for (j = 1; j <= i__1; ++j) {
	do_fio(&c__1, name__ + j, (ftnlen)1);
    }
    e_wsfe();

    cl__1.cerr = 0;
    cl__1.cunit = *iustat;
    cl__1.csta = 0;
    f_clos(&cl__1);
    o__1.oerr = 0;
    o__1.ounit = *iustat;
    o__1.ofnmlen = 132;
    o__1.ofnm = statfile;
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = "formatted";
    o__1.oblnk = 0;
    f_open(&o__1);

    charsp_(statfile, name__ + 1, len, &c__132, &nlb, &ntb, (ftnlen)132, (
	    ftnlen)1);

    s_wsfe(&io___29);
    do_fio(&c__1, (char *)&(*iustat), (ftnlen)sizeof(integer));
    do_fio(&c__1, " Statistics  ", (ftnlen)13);
    i__1 = *len;
    for (j = 1; j <= i__1; ++j) {
	do_fio(&c__1, name__ + j, (ftnlen)1);
    }
    e_wsfe();

    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
	cl__1.cerr = 0;
	cl__1.cunit = *iutrn;
	cl__1.csta = 0;
	f_clos(&cl__1);
	if (*ifrmout == 1) {
	    o__1.oerr = 0;
	    o__1.ounit = *iutrn;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = trnfile;
	    o__1.orl = 0;
	    o__1.osta = "unknown";
	    o__1.oacc = 0;
	    o__1.ofm = "formatted";
	    o__1.oblnk = 0;
	    f_open(&o__1);
	} else {
	    o__1.oerr = 0;
	    o__1.ounit = *iutrn;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = trnfile;
	    o__1.orl = 0;
	    o__1.osta = "unknown";
	    o__1.oacc = 0;
	    o__1.ofm = "unformatted";
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}

	charsp_(trnfile, name__ + 1, len, &c__132, &nlb, &ntb, (ftnlen)132, (
		ftnlen)1);

	s_wsfe(&io___30);
	do_fio(&c__1, (char *)&(*iutrn), (ftnlen)sizeof(integer));
	do_fio(&c__1, " Transmission  ", (ftnlen)15);
	i__1 = *len;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, name__ + j, (ftnlen)1);
	}
	e_wsfe();
    }

    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {
	iuf = *iuflx + nz - 1;
	cl__1.cerr = 0;
	cl__1.cunit = iuf;
	cl__1.csta = 0;
	f_clos(&cl__1);
	if (*irad == 3 || *irad == 4 || *irad == 7 || *irad == 8) {
	    if (*ifrmout == 1) {
		o__1.oerr = 0;
		o__1.ounit = iuf;
		o__1.ofnmlen = 132;
		o__1.ofnm = flxfile + nz * 132;
		o__1.orl = 0;
		o__1.osta = "unknown";
		o__1.oacc = 0;
		o__1.ofm = "formatted";
		o__1.oblnk = 0;
		f_open(&o__1);
	    } else {
		o__1.oerr = 0;
		o__1.ounit = iuf;
		o__1.ofnmlen = 132;
		o__1.ofnm = flxfile + nz * 132;
		o__1.orl = 0;
		o__1.osta = "unknown";
		o__1.oacc = 0;
		o__1.ofm = "unformatted";
		o__1.oblnk = 0;
		f_open(&o__1);
	    }

	    charsp_(flxfile + nz * 132, name__ + 1, len, &c__132, &nlb, &ntb, 
		    (ftnlen)132, (ftnlen)1);

	    s_wsfe(&io___31);
	    do_fio(&c__1, (char *)&iuf, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " Flux  ", (ftnlen)7);
	    i__2 = *len;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, name__ + j, (ftnlen)1);
	    }
	    e_wsfe();
	}
	i__2 = *nlout;
	for (nl = 1; nl <= i__2; ++nl) {
	    iuo = *iuout + (nz - 1) * *nlout + nl - 1;
	    cl__1.cerr = 0;
	    cl__1.cunit = iuo;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    if (*ifrmout == 1) {
		o__1.oerr = 0;
		o__1.ounit = iuo;
		o__1.ofnmlen = 132;
		o__1.ofnm = radfile + (nl + nz * 3) * 132;
		o__1.orl = 0;
		o__1.osta = "unknown";
		o__1.oacc = 0;
		o__1.ofm = "formatted";
		o__1.oblnk = 0;
		f_open(&o__1);
	    } else {
		o__1.oerr = 0;
		o__1.ounit = iuo;
		o__1.ofnmlen = 132;
		o__1.ofnm = radfile + (nl + nz * 3) * 132;
		o__1.orl = 0;
		o__1.osta = "unknown";
		o__1.oacc = 0;
		o__1.ofm = "unformatted";
		o__1.oblnk = 0;
		f_open(&o__1);
	    }

	    charsp_(radfile + (nl + nz * 3) * 132, name__ + 1, len, &c__132, &
		    nlb, &ntb, (ftnlen)132, (ftnlen)1);

	    s_wsfe(&io___32);
	    do_fio(&c__1, (char *)&iuo, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " Radiances  ", (ftnlen)12);
	    i__3 = *len;
	    for (j = 1; j <= i__3; ++j) {
		do_fio(&c__1, name__ + j, (ftnlen)1);
	    }
	    e_wsfe();
/* L7621: */
	}
/* L7661: */
    }

    return 0;
} /* init_spect_io__ */

