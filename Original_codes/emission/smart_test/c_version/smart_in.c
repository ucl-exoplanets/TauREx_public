/* smart_in.f -- translated by f2c (version 20100827).
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
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__132 = 132;
static integer c__9 = 9;
static integer c__50 = 50;
static integer c__6 = 6;
static integer c__10 = 10;
static integer c__5 = 5;
static integer c__8 = 8;
static integer c__16 = 16;

/* Subroutine */ int smart_in__(integer *iuatm, char *atmfile, integer *
	ifrmatm, char *frmatm, integer *iskatm, integer *icp, integer *ict, 
	real *scp, integer *ngases, integer *ngastyp, integer *igas, integer *
	iugas, integer *iuabc, integer *igtrn, integer *igpab, real *wgtgas, 
	char *gasfile, integer *iumix, char *mixfile, integer *ifrmmix, char *
	frmmix, integer *ioffmix, integer *izmix, integer *imix, integer *
	icpmix, integer *icmix, real *scpmix, real *scmix, integer *nmodes, 
	integer *iuaer, char *aerfile, integer *iumie, char *miefile, integer 
	*ioffmie, integer *ioffmom, integer *iang, integer *icwlq, integer *
	icqext, integer *icqsca, integer *icg1, integer *icg2, integer *icf, 
	doublereal *wnaer, integer *iofftau, integer *iztau, integer *icptau, 
	integer *ictau, real *scptau, real *sctau, logical *lcsh, integer *
	iccsh, real *sccsh, integer *iusur, char *surfile, integer *ifrmsur, 
	char *frmsur, integer *ioffsur, real *ws, real *phiw, integer *iwnsur,
	 integer *icwnsur, integer *icalb, real *scwalb, real *scalb, integer 
	*iusol1, char *solfile, integer *ifrms0, char *frms0, integer *ioffs0,
	 integer *ixs0, integer *icwns0, integer *icsol, real *scwns0, real *
	scsol, integer *iuins0, integer *ncomp, integer *icomp, real *volmix, 
	real *ts, real *au, real *radius, real *sgrav, real *wgtatm, 
	doublereal *wnmin, doublereal *wnmax, integer *isptype, integer *
	islit, doublereal *width, doublereal *dwn, integer *nstr, logical *
	usrang, integer *numu, real *umu, integer *nphi, real *phi, integer *
	nza, real *umu0, real *phi0, integer *levout, integer *nlout, real *
	pout, real *accur, logical *lamber, logical *lplanck, logical *lsolar,
	 integer *isource, integer *irad, integer *iref, integer *nref, 
	integer *iuout, integer *iunits, integer *ifrmout, integer *iuheat, 
	integer *iustat, integer *iutrn, integer *iuflx, real *tauerr, real *
	pi0err, real *phferr, real *surferr, real *pd_frac__, integer *nstate,
	 integer *istate, integer *iu_pd__, integer *iutpd, integer *iuspd, 
	integer *iupdrad, ftnlen atmfile_len, ftnlen frmatm_len, ftnlen 
	gasfile_len, ftnlen mixfile_len, ftnlen frmmix_len, ftnlen 
	aerfile_len, ftnlen miefile_len, ftnlen surfile_len, ftnlen 
	frmsur_len, ftnlen solfile_len, ftnlen frms0_len)
{
    /* Initialized data */

    static char cgas[6*50] = "   H2O" "   CO2" "    O3" "   N2O" "    CO" 
	    "   CH4" "    O2" "    NO" "   SO2" "   NO2" "   NH3" "  HNO3" 
	    "    OH" "    HF" "   HCl" "   HBr" "    HI" "   ClO" "   OCS" 
	    "  H2CO" "  HOCl" "    N2" "   HCN" " CH3Cl" "  H2O2" "  C2H2" 
	    "  C2H6" "   PH3" "  COF2" "   SF6" "   H2S" " HCOOH" "   HO2" 
	    "     O" "ClONO2" "   NO+" "  HOBr" "  C2H4" "CH3OH " "H2    " 
	    "He    " "Other " "Other " "Other " "Other " "Other " "Other " 
	    "Other " "Other " "Other ";

    /* System generated locals */
    integer i__1, i__2, i__3;
    icilist ici__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    double acos(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsle(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsfe(cilist *), e_rsfe(void), s_wsfi(icilist *), e_wsfi(void), 
	    f_open(olist *), f_clos(cllist *), s_wsle(cilist *), e_wsle(void);
    double cos(doublereal);

    /* Local variables */
    static char heatfile[132], statfile[132];
    static integer i__, j, k, m, n, ic;
    static real pi;
    static integer nl, ir;
    extern /* Subroutine */ int init_pd_iu__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, char *, logical *, logical *, 
	    char *, integer *, integer *, integer *, char *, integer *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static real ang;
    static integer nlb, len, ntb, ngt, nze, ist;
    static real sza0[10];
    static char name__[1*132], clev[4*3];
    static integer icmx, nzdn, nzup, nza_1__[3], npabs;
    static char j_ext__[9*113];
    extern /* Subroutine */ int init_spect_io__(integer *, integer *, integer 
	    *, integer *, real *, char *, integer *, integer *, integer *, 
	    integer *, integer *, char *, integer *, char *, char *, char *, 
	    char *, char *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static real sctau0;
    static char filein[132];
    extern /* Subroutine */ int charsp_(char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char radfile[132*3*10], flxfile[132*10], trnfile[132];
    static integer isurtyp;

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___7 = { 0, 6, 0, "(1a,3(/,1a))", 0 };
    static cilist io___8 = { 1, 5, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, "(3(1x,1a,i5),2a)", 0 };
    static cilist io___12 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___13 = { 0, 5, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___15 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___16 = { 0, 6, 0, "(/,2a,/,1a)", 0 };
    static cilist io___17 = { 0, 6, 0, "(1a,3(/,1a))", 0 };
    static cilist io___18 = { 1, 5, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, "(3(1x,1a,i5),2a)", 0 };
    static cilist io___20 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___21 = { 0, 5, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___23 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___24 = { 0, 6, 0, "(2(/,1a))", 0 };
    static cilist io___25 = { 1, 5, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___27 = { 0, 6, 0, "(2(/,1a))", 0 };
    static cilist io___28 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___29 = { 0, 6, 0, "(1x,2a)", 0 };
    static cilist io___30 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___31 = { 0, 5, 0, "(1a)", 0 };
    static icilist io___37 = { 0, filein, 0, "(132a)", 132, 1 };
    static cilist io___39 = { 0, 6, 0, "(/,1a,i5,2a)", 0 };
    static cilist io___40 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___41 = { 1, 5, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___43 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___44 = { 1, 5, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, "(1x,2(1a,i5))", 0 };
    static cilist io___46 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___47 = { 1, 5, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___49 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___50 = { 1, 5, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___52 = { 0, 6, 0, "(4(/,1a))", 0 };
    static cilist io___53 = { 1, 5, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___57 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___58 = { 1, 5, 0, 0, 0 };
    static cilist io___59 = { 0, 6, 0, "(1x,1a,i5,10a)", 0 };
    static cilist io___60 = { 0, 6, 0, "(/,10a)", 0 };
    static cilist io___61 = { 0, 6, 0, "(1a,3(/,1a))", 0 };
    static cilist io___62 = { 1, 5, 0, 0, 0 };
    static cilist io___63 = { 0, 6, 0, "(3(1x,1a,i5),2a)", 0 };
    static cilist io___64 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___65 = { 0, 5, 0, 0, 0 };
    static cilist io___66 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___67 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___68 = { 0, 6, 0, "(/,10a)", 0 };
    static cilist io___69 = { 1, 5, 0, 0, 0 };
    static cilist io___70 = { 0, 6, 0, 0, 0 };
    static cilist io___71 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___73 = { 0, 6, 0, "(/,2a,/,3(/,1a))", 0 };
    static cilist io___74 = { 1, 5, 0, 0, 0 };
    static cilist io___75 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___76 = { 0, 6, 0, "(/,10a)", 0 };
    static cilist io___77 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___78 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___79 = { 0, 6, 0, "(/,2a,/,2x,1a)", 0 };
    static cilist io___80 = { 0, 6, 0, "(2(/,1a))", 0 };
    static cilist io___81 = { 1, 5, 0, 0, 0 };
    static cilist io___82 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___83 = { 0, 6, 0, "(2(/,1a))", 0 };
    static cilist io___84 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___85 = { 0, 6, 0, "(1x,2a)", 0 };
    static cilist io___86 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___87 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___88 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___89 = { 0, 6, 0, "(1x,10a)", 0 };
    static cilist io___90 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___91 = { 1, 5, 0, 0, 0 };
    static cilist io___92 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___93 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___94 = { 1, 5, 0, 0, 0 };
    static cilist io___95 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___96 = { 0, 6, 0, "(/,10a)", 0 };
    static cilist io___97 = { 1, 5, 0, 0, 0 };
    static cilist io___98 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___99 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___100 = { 0, 6, 0, "(2a,/,2(/,1a))", 0 };
    static cilist io___101 = { 1, 5, 0, 0, 0 };
    static cilist io___102 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___103 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___104 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___105 = { 1, 5, 0, 0, 0 };
    static cilist io___106 = { 0, 6, 0, "(1x,1a,1pe13.6)", 0 };
    static cilist io___107 = { 0, 6, 0, "(1x,1a,1pe13.6)", 0 };
    static cilist io___108 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___109 = { 1, 5, 0, 0, 0 };
    static cilist io___110 = { 0, 6, 0, "(1x,1a,1pe13.6)", 0 };
    static cilist io___111 = { 0, 6, 0, "(/,/,1a)", 0 };
    static cilist io___112 = { 1, 5, 0, 0, 0 };
    static cilist io___113 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___114 = { 0, 6, 0, "(1a,i5)", 0 };
    static cilist io___116 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___117 = { 0, 6, 0, "(1a,3(/,1a))", 0 };
    static cilist io___118 = { 1, 5, 0, 0, 0 };
    static cilist io___119 = { 0, 6, 0, "(3(1x,1a,i5),2a)", 0 };
    static cilist io___120 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___121 = { 0, 5, 0, 0, 0 };
    static cilist io___122 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___123 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___124 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___125 = { 0, 5, 0, "(1a)", 0 };
    static icilist io___126 = { 0, filein, 0, "(132a)", 132, 1 };
    static cilist io___127 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___128 = { 0, 6, 0, "(1a,i5)", 0 };
    static cilist io___129 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___130 = { 1, 5, 0, 0, 0 };
    static cilist io___131 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___132 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___133 = { 1, 5, 0, 0, 0 };
    static cilist io___134 = { 0, 6, 0, "(1x,1a,4i5)", 0 };
    static cilist io___135 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___136 = { 0, 6, 0, "(3(/1a))", 0 };
    static cilist io___137 = { 1, 5, 0, 0, 0 };
    static cilist io___138 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___139 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___140 = { 1, 5, 0, 0, 0 };
    static cilist io___141 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___142 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___143 = { 1, 5, 0, 0, 0 };
    static cilist io___144 = { 0, 6, 0, "(1x,1a,2i5)", 0 };
    static cilist io___145 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___146 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___147 = { 0, 6, 0, "(1x,2a)", 0 };
    static cilist io___148 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___149 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___150 = { 1, 5, 0, 0, 0 };
    static cilist io___151 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___152 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___153 = { 0, 5, 0, 0, 0 };
    static cilist io___154 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___155 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___156 = { 1, 5, 0, 0, 0 };
    static cilist io___157 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___158 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___159 = { 1, 5, 0, 0, 0 };
    static cilist io___160 = { 0, 6, 0, "(1x,1a,2i5)", 0 };
    static cilist io___161 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___162 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___163 = { 1, 5, 0, 0, 0 };
    static cilist io___165 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___166 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___167 = { 1, 5, 0, 0, 0 };
    static cilist io___168 = { 0, 6, 0, "(1x,1a,l6)", 0 };
    static cilist io___169 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___170 = { 1, 5, 0, 0, 0 };
    static cilist io___171 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___172 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___173 = { 1, 5, 0, 0, 0 };
    static cilist io___174 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___175 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___176 = { 1, 5, 0, 0, 0 };
    static cilist io___178 = { 0, 6, 0, "(/,1a,4(/,1a))", 0 };
    static cilist io___179 = { 1, 5, 0, 0, 0 };
    static cilist io___180 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___181 = { 0, 6, 0, "(1x,1a,l6)", 0 };
    static cilist io___182 = { 0, 6, 0, "(/,2a,/,1a)", 0 };
    static cilist io___183 = { 0, 6, 0, "(1a,3(/,1a))", 0 };
    static cilist io___184 = { 1, 5, 0, 0, 0 };
    static cilist io___185 = { 0, 6, 0, "(3(1x,1a,i5),2a)", 0 };
    static cilist io___186 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___187 = { 1, 5, 0, 0, 0 };
    static cilist io___188 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___189 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___190 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___191 = { 1, 5, 0, 0, 0 };
    static cilist io___192 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___193 = { 0, 6, 0, "(2(/,1a))", 0 };
    static cilist io___194 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___195 = { 0, 6, 0, "(1x,2a)", 0 };
    static cilist io___196 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___197 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___198 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___199 = { 0, 6, 0, "(2a)", 0 };
    static cilist io___200 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___201 = { 1, 5, 0, 0, 0 };
    static cilist io___202 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___203 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___204 = { 1, 5, 0, 0, 0 };
    static cilist io___206 = { 0, 6, 0, "(1x,1a,5i5)", 0 };
    static cilist io___207 = { 0, 6, 0, "(2(/,1a))", 0 };
    static cilist io___208 = { 1, 5, 0, 0, 0 };
    static cilist io___209 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___210 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___211 = { 1, 5, 0, 0, 0 };
    static cilist io___212 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___213 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___214 = { 1, 5, 0, 0, 0 };
    static cilist io___215 = { 0, 6, 0, "(1x,1a,4(1pe13.5))", 0 };
    static cilist io___216 = { 0, 6, 0, "(/,1x,1a)", 0 };
    static cilist io___217 = { 1, 5, 0, 0, 0 };
    static cilist io___218 = { 0, 6, 0, "(1x,1a,1pe12.4,1a,1pe12.4,1a)", 0 };
    static cilist io___219 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___220 = { 1, 5, 0, 0, 0 };
    static cilist io___221 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___222 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___223 = { 1, 5, 0, 0, 0 };
    static cilist io___224 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___225 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___226 = { 1, 5, 0, 0, 0 };
    static cilist io___227 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___228 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___229 = { 1, 5, 0, 0, 0 };
    static cilist io___230 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___231 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___232 = { 1, 5, 0, 0, 0 };
    static cilist io___233 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___234 = { 0, 6, 0, 0, 0 };
    static cilist io___236 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___237 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___238 = { 1, 5, 0, 0, 0 };
    static cilist io___239 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___240 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___241 = { 1, 5, 0, 0, 0 };
    static cilist io___242 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___243 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___244 = { 1, 5, 0, 0, 0 };
    static cilist io___245 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___246 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___247 = { 0, 6, 0, "(5(/,1a))", 0 };
    static cilist io___248 = { 1, 5, 0, 0, 0 };
    static cilist io___249 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___250 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___251 = { 1, 5, 0, 0, 0 };
    static cilist io___252 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___253 = { 0, 6, 0, "(2(/,1a))", 0 };
    static cilist io___254 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___255 = { 0, 6, 0, "(1x,2a)", 0 };
    static cilist io___256 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___257 = { 0, 5, 0, "(1a)", 0 };
    static cilist io___258 = { 0, 6, 0, "(1x,2a)", 0 };
    static cilist io___259 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___260 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___261 = { 1, 5, 0, 0, 0 };
    static cilist io___262 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___263 = { 0, 6, 0, "(5(/,1a))", 0 };
    static cilist io___264 = { 1, 5, 0, 0, 0 };
    static cilist io___265 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___266 = { 0, 6, 0, "(4(/,1a))", 0 };
    static cilist io___267 = { 1, 5, 0, 0, 0 };
    static cilist io___268 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___269 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___270 = { 1, 5, 0, 0, 0 };
    static cilist io___271 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___272 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___273 = { 1, 5, 0, 0, 0 };
    static cilist io___274 = { 0, 6, 0, "(1x,1a,2i5)", 0 };
    static cilist io___275 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___276 = { 1, 5, 0, 0, 0 };
    static cilist io___277 = { 0, 6, 0, "(1x,1a,2i5)", 0 };
    static cilist io___278 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___279 = { 1, 5, 0, 0, 0 };
    static cilist io___280 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___281 = { 0, 6, 0, "(/,2a,i5)", 0 };
    static cilist io___282 = { 0, 6, 0, "(/,2a,/,1a,/,1a,/1a)", 0 };
    static cilist io___283 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___284 = { 1, 5, 0, 0, 0 };
    static cilist io___286 = { 0, 6, 0, "(1x,1a,2(1pe13.5))", 0 };
    static cilist io___287 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___288 = { 0, 5, 0, 0, 0 };
    static cilist io___289 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___290 = { 0, 6, 0, "(16(/,1a))", 0 };
    static cilist io___291 = { 0, 6, 0, "(1a,19(/,1a))", 0 };
    static cilist io___292 = { 1, 5, 0, 0, 0 };
    static cilist io___293 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___294 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___295 = { 1, 5, 0, 0, 0 };
    static cilist io___296 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___297 = { 0, 6, 0, "(/,2a,i5)", 0 };
    static cilist io___298 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___299 = { 0, 6, 0, "(1a,i5)", 0 };
    static cilist io___300 = { 1, 5, 0, 0, 0 };
    static cilist io___302 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___303 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___304 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___305 = { 1, 5, 0, 0, 0 };
    static cilist io___306 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___307 = { 0, 6, 0, "(/,2a,i5)", 0 };
    static cilist io___308 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___309 = { 0, 6, 0, "(1a,i5)", 0 };
    static cilist io___310 = { 1, 5, 0, 0, 0 };
    static cilist io___311 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___312 = { 0, 6, 0, "(6(/1a))", 0 };
    static cilist io___313 = { 1, 5, 0, 0, 0 };
    static cilist io___314 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___316 = { 0, 6, 0, "(/,1a,/,1a,i5)", 0 };
    static cilist io___317 = { 1, 5, 0, 0, 0 };
    static cilist io___318 = { 0, 6, 0, 0, 0 };
    static cilist io___319 = { 0, 6, 0, 0, 0 };
    static cilist io___321 = { 0, 6, 0, "(/,1a,i5)", 0 };
    static cilist io___322 = { 1, 5, 0, 0, 0 };
    static cilist io___323 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___324 = { 0, 6, 0, "(6(/,1a))", 0 };
    static cilist io___325 = { 1, 5, 0, 0, 0 };
    static cilist io___326 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___327 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___328 = { 1, 5, 0, 0, 0 };
    static cilist io___329 = { 0, 6, 0, "(1x,1a,2(1pe14.6))", 0 };
    static cilist io___330 = { 0, 6, 0, "(4(/,1a))", 0 };
    static cilist io___331 = { 1, 5, 0, 0, 0 };
    static cilist io___332 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___333 = { 0, 6, 0, "(9(/,1a))", 0 };
    static cilist io___334 = { 1, 5, 0, 0, 0 };
    static cilist io___335 = { 0, 6, 0, "(1x,1a,i5)", 0 };
    static cilist io___336 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___337 = { 1, 5, 0, 0, 0 };
    static cilist io___338 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___339 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___340 = { 1, 5, 0, 0, 0 };
    static cilist io___341 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___342 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___343 = { 1, 5, 0, 0, 0 };
    static cilist io___344 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___345 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___346 = { 1, 5, 0, 0, 0 };
    static cilist io___347 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___348 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___349 = { 1, 5, 0, 0, 0 };
    static cilist io___350 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___351 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___352 = { 1, 5, 0, 0, 0 };
    static cilist io___353 = { 0, 6, 0, "(1x,1a,1pe13.5)", 0 };
    static cilist io___354 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___355 = { 1, 5, 0, 0, 0 };
    static cilist io___356 = { 0, 6, 0, "(1a,i5)", 0 };



/* cccccccccccccccccccccccccc  s m a r t _ i n  cccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    This subroutine reads list-directed input for the program smart cc */
/* c                                                                    cc */
/* c    i n p u t:                                                      cc */
/* c                                                                    cc */
/* c      iusol1 - unit number with solar radiances                     cc */
/* c     iuthrm - unit number for thermal radiances                     cc */
/* c      iuout - unit number for output radiance file                  cc */
/* c      iutrn - unit number for output transmission/pressure file     cc */
/* c      iuflx - unit number for output level-dependent fluxes         cc */
/* c     iuheat - unit number for output solar heating rates            cc */
/* c     iustat - unit number for output binning statistics             cc */
/* c      iutpd - unit numbers for output level-dependent               cc */
/* c              thermal fluxes and thier partial derivatives          cc */
/* c      iuspd - unit numbers for output level-dependent               cc */
/* c              solar fluxes and thier partial derivatives            cc */
/* c    iupdrad - unit numbers for output level-dependent radiances     cc */
/* c              and thier partial derivatives                         cc */
/* c    atmfile - name of input atmospheric structure file              cc */
/* c    aerfile - name of input aerosol vertical structure file         cc */
/* c    miefile - name of file with aerosol optical properties vs. wn   cc */
/* c    solfile - name of file with wn-dependent solar fluxes           cc */
/* c    surfile - name of file with wn-dependent surface optical prop.  cc */
/* c    mixfile - name of file with gas mixing ratios                   cc */
/* c    gasfile - name of file with gas absorption coeffiecients vs. wn cc */
/* c    radfile - name output file for flux/radiance spectra            cc */
/* c   heatfile - name of output file with heating/cooling rates        cc */
/* c   statfile - name of output file with binning statistics           cc */
/* c    trnfile - name of output file with tau=1 trans/pressure         cc */
/* c    flxfile - name of output file with level-dependent flux spectra cc */
/* c     nstate - number of variable elements in the state vecror       cc */
/* c     istate - state vector flag indicating which state variables    cc */
/* c              are variable components of the state vector.          cc */
/* c              0 - not a variable component of the state vector      cc */
/* c              1 - surface pressure                                  cc */
/* c              2 - surface/atmospheric temperature                   cc */
/* c              3 - gas absorption coeffient                          cc */
/* c              4 - cloud/aerosol optical depth                       cc */
/* c              5 - surface albedo                                    cc */
/* c     nmodes - number of discrete aerosol partical modes             cc */
/* c      ncomp - number of rayleigh-scattering constituentes           cc */
/* c      icomp - index of each rayleigh scattering constituent         cc */
/* c     volmix - volume mixing ratio of each rayleigh scatterer        cc */
/* c         ts - sufrace temperature (K)                               cc */
/* c         au - distance to the sun (in AU's)                         cc */
/* c     tauerr - optical depth relative binning error (0. to ~0.8)     cc */
/* c     pi0err - co-single scattering albedo absolute binning error    cc */
/* c     phferr - asymmetry factor absolute binning error               cc */
/* c    surferr - surface optical property binning error                cc */
/* c       pout - output pressure level (bars)                          cc */
/* c    isource - index of source function type: (1) solar only         cc */
/* c              (2) thermal only, (3) both                            cc */
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
/* c       iref - bidirectional reflectance options                     cc */
/* c              0 - Lambert                                           cc */
/* c              1 - Hapke's BDR model                                 cc */
/* c              2 - Breon's BDR model; combination of Li + Roujean    cc */
/* c              3 - Roujean's BDR model                               cc */
/* c              4 - Cox and Munk glint model                          cc */
/* c       sza0 - solar zenith angles (degrees)                         cc */
/* c       umu0 - cosine of solar zenith angles                         cc */
/* c       phi0 - solar azimuth angles (degrees)                        cc */
/* c       nstr - number of gaussian zenith angles used in D/O code     cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       nlev - number of output levels                               cc */
/* c       nlyr - number of computational model layers                  cc */
/* c      nlout - number of levels where radiance spectra are output    cc */
/* c     levout - output level index (1) top of atmosphere,             cc */
/* c              (2) surface, (3) arbitrary level                      cc */
/* c        nza - number of solar zenith angles                         cc */
/* c        phi - emission azimuth angles (degrees)                     cc */
/* c        umu - emission zenith angle cosines                         cc */
/* c      accur - azimuth convergence accuracy for D/O routine          cc */
/* c     lamber - Include a lambertian surface? (Logical: T/F)          cc */
/* c              note: if lamber = F, a bi-directional reflection      cc */
/* c              function is required.                                 cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c     usrang - output radiances at user angles? (logical: T/F)       cc */
/* c     iunits - index of output radiance units:                       cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) ergs/s/cm**2/sr/cm-1                               cc */
/* c              5) photons/s/m**2/sr/micron                           cc */
/* c      wnmin - minimum wavenumber of desired spectral window         cc */
/* c      wnmax - maximum wavenumber of desired spectral window         cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    same as input                                                   cc */
/* c                                                                    cc */
/* c    m o d i f i c a t i o n    h i s t o r y                        cc */
/* c                                                                    cc */
/* c    8/19/00: an aerosol optical depth scaling factor, sctau(nmode), cc */
/* c             was added to the input prompt for the optical depth    cc */
/* c             pressure scaling factor. In its current form, the code cc */
/* c             should work whether sctau is specified or not, as long cc */
/* c             as the input record includes a non-numeric character.  cc */
/* c    4/15/06  modified to support radiance and flux jacobian files   cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  s m a r t _ i n  cccccccccccccccccccccccccc */




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




    /* Parameter adjustments */
    iupdrad -= 10513;
    iuspd -= 11;
    iutpd -= 11;
    --istate;
    --pd_frac__;
    --pout;
    --phi0;
    --umu0;
    --phi;
    --umu;
    --volmix;
    --icomp;
    --scalb;
    --icalb;
    --sccsh;
    --iccsh;
    --lcsh;
    --sctau;
    --scptau;
    --ictau;
    --icptau;
    --iztau;
    --iofftau;
    --wnaer;
    --icf;
    --icg2;
    --icg1;
    --icqsca;
    --icqext;
    --icwlq;
    --iang;
    --ioffmom;
    --ioffmie;
    miefile -= 132;
    aerfile -= 132;
    --scmix;
    --scpmix;
    --icmix;
    --icpmix;
    --imix;
    --izmix;
    --ioffmix;
    frmmix -= 40;
    --ifrmmix;
    mixfile -= 132;
    gasfile -= 528;
    --wgtgas;
    igpab -= 4;
    igtrn -= 4;
    iuabc -= 4;
    --igas;
    --ngastyp;

    /* Function Body */

    pi = acos(-1.f);

/* ****  set maximimum number of retries for input errors */

    icmx = 10;

    for (n = 1; n <= 113; ++n) {
	istate[n] = 0;
/* L1001: */
    }

/* ****   initialize the number of values in the state structure: */
/*       The first 2 values are pressure and temperature */

    ic = 0;
L1101:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    *nstate = 1;
    istate[1] = 0;
    s_wsfe(&io___6);
    do_fio(&c__1, " Compute Jacobians for Pressure?", (ftnlen)32);
    e_wsfe();
    s_wsfe(&io___7);
    do_fio(&c__1, " Enter index of Jacobian Type: ", (ftnlen)31);
    do_fio(&c__1, " 0) No Jacobians", (ftnlen)16);
    do_fio(&c__1, " 1) Radiance Jacobians", (ftnlen)22);
    do_fio(&c__1, " 2) Flux Jacobians", (ftnlen)18);
    e_wsfe();
    i__1 = s_rsle(&io___8);
    if (i__1 != 0) {
	goto L1101;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&ist, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L1101;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L1101;
    }

    if (ist != 0) {
	s_copy(j_ext__, ".j_pres", (ftnlen)9, (ftnlen)7);
	if (ist == 1) {
	    istate[*nstate] = 1;
	} else {
	    istate[*nstate] = -1;
	}

	s_wsfe(&io___11);
	do_fio(&c__1, "ist =", (ftnlen)5);
	do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	do_fio(&c__1, " nstate=", (ftnlen)8);
	do_fio(&c__1, (char *)&(*nstate), (ftnlen)sizeof(integer));
	do_fio(&c__1, " istate =", (ftnlen)9);
	do_fio(&c__1, (char *)&istate[*nstate], (ftnlen)sizeof(integer));
	do_fio(&c__1, " j_ext =", (ftnlen)8);
	do_fio(&c__1, j_ext__ + (*nstate - 1) * 9, (ftnlen)9);
	e_wsfe();

	s_wsfe(&io___12);
	do_fio(&c__1, " Enter the fractional change in pressure (0.0-1.0): ", 
		(ftnlen)52);
	e_wsfe();
	s_rsle(&io___13);
	do_lio(&c__4, &c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(real)
		);
	e_rsle();
	s_wsfe(&io___14);
	do_fio(&c__1, "pd_frac =", (ftnlen)9);
	do_fio(&c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(real));
	e_wsfe();
    } else {
	s_wsfe(&io___15);
	do_fio(&c__1, "ist =", (ftnlen)5);
	do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	e_wsfe();
    }

    *nstate = 2;
    istate[*nstate] = 0;
    pd_frac__[*nstate] = 0.f;
L1111:
    s_wsfe(&io___16);
    do_fio(&c__1, " Compute Jacobians for Temperature?", (ftnlen)35);
    e_wsfe();
    s_wsfe(&io___17);
    do_fio(&c__1, " Enter index of Jacobian Type: ", (ftnlen)31);
    do_fio(&c__1, " 0) No Jacobians", (ftnlen)16);
    do_fio(&c__1, " 1) Radiance Jacobians", (ftnlen)22);
    do_fio(&c__1, " 2) Flux Jacobians", (ftnlen)18);
    e_wsfe();
    i__1 = s_rsle(&io___18);
    if (i__1 != 0) {
	goto L1111;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&ist, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L1111;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L1111;
    }

    if (ist != 0) {
	s_copy(j_ext__ + (*nstate - 1) * 9, ".j_temp", (ftnlen)9, (ftnlen)7);
	if (ist == 1) {
	    istate[*nstate] = 2;
	} else {
	    istate[*nstate] = -2;
	}
	s_wsfe(&io___19);
	do_fio(&c__1, "ist =", (ftnlen)5);
	do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	do_fio(&c__1, " nstate=", (ftnlen)8);
	do_fio(&c__1, (char *)&(*nstate), (ftnlen)sizeof(integer));
	do_fio(&c__1, " istate =", (ftnlen)9);
	do_fio(&c__1, (char *)&istate[*nstate], (ftnlen)sizeof(integer));
	do_fio(&c__1, " j_ext =", (ftnlen)8);
	do_fio(&c__1, j_ext__ + (*nstate - 1) * 9, (ftnlen)9);
	e_wsfe();
	s_wsfe(&io___20);
	do_fio(&c__1, " Enter the fractional change in Temperature (0.0-1.0)"
		": ", (ftnlen)55);
	e_wsfe();
	s_rsle(&io___21);
	do_lio(&c__4, &c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(real)
		);
	e_rsle();
	s_wsfe(&io___22);
	do_fio(&c__1, "pd_frac =", (ftnlen)9);
	do_fio(&c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(real));
	e_wsfe();
    } else {
	s_wsfe(&io___23);
	do_fio(&c__1, "ist =", (ftnlen)5);
	do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/* ****  enter input file names and other data for this run. */

    icmx = 5;
    ic = 0;
L1121:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___24);
    do_fio(&c__1, " Enter atmospheric structure file format index (1 - 3):  ",
	     (ftnlen)57);
    do_fio(&c__1, "  1) formatted 2) unformatted 3) list directed ", (ftnlen)
	    47);
    e_wsfe();
    i__1 = s_rsle(&io___25);
    if (i__1 != 0) {
	goto L1121;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ifrmatm), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L1121;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L1121;
    }
    s_wsfe(&io___26);
    do_fio(&c__1, "ifrmatm =", (ftnlen)9);
    do_fio(&c__1, (char *)&(*ifrmatm), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*ifrmatm == 1) {
	s_wsfe(&io___27);
	do_fio(&c__1, " Enter format for reading pressure and temperature", (
		ftnlen)50);
	do_fio(&c__1, " [enclosed in parenthesis, ie; (2f3.5) ]: ", (ftnlen)
		42);
	e_wsfe();
	s_rsfe(&io___28);
	do_fio(&c__1, frmatm, (ftnlen)40);
	e_rsfe();
	s_wsfe(&io___29);
	do_fio(&c__1, "frmatm = ", (ftnlen)9);
	do_fio(&c__1, frmatm, (ftnlen)40);
	e_wsfe();
    } else {
	s_copy(frmatm, " ", (ftnlen)40, (ftnlen)1);
    }

    ic = 0;
L1132:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_copy(atmfile, " ", (ftnlen)132, (ftnlen)1);
    s_wsfe(&io___30);
    do_fio(&c__1, " Enter name of the Atmospheric Structure file: ", (ftnlen)
	    47);
    e_wsfe();
    s_rsfe(&io___31);
    do_fio(&c__1, atmfile, (ftnlen)132);
    e_rsfe();

    charsp_(atmfile, name__, &len, &c__132, &nlb, &ntb, (ftnlen)132, (ftnlen)
	    1);

    s_copy(filein, " ", (ftnlen)132, (ftnlen)1);
    s_wsfi(&io___37);
    i__1 = len;
    for (j = 1; j <= i__1; ++j) {
	do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
    }
    e_wsfi();
    s_wsfe(&io___39);
    do_fio(&c__1, " unit: ", (ftnlen)7);
    do_fio(&c__1, (char *)&(*iuatm), (ftnlen)sizeof(integer));
    do_fio(&c__1, " atmfile:", (ftnlen)9);
    do_fio(&c__1, filein, len);
    e_wsfe();

    if (*ifrmatm != 2) {
	o__1.oerr = 1;
	o__1.ounit = *iuatm;
	o__1.ofnmlen = len;
	o__1.ofnm = filein;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "formatted";
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L1132;
	}
    } else {
	o__1.oerr = 1;
	o__1.ounit = *iuatm;
	o__1.ofnmlen = len;
	o__1.ofnm = filein;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L1132;
	}
    }

    cl__1.cerr = 0;
    cl__1.cunit = *iuatm;
    cl__1.csta = 0;
    f_clos(&cl__1);

    ic = 0;
L1141:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___40);
    do_fio(&c__1, " Enter number of lines to skip at top of this file:", (
	    ftnlen)51);
    e_wsfe();
    i__1 = s_rsle(&io___41);
    if (i__1 != 0) {
	goto L1141;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*iskatm), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L1141;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L1141;
    }
    s_wsfe(&io___42);
    do_fio(&c__1, "iskatm =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*iskatm), (ftnlen)sizeof(integer));
    e_wsfe();
    ic = 0;
L1161:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___43);
    do_fio(&c__1, " Enter columns containing p and T:", (ftnlen)34);
    e_wsfe();
    i__1 = s_rsle(&io___44);
    if (i__1 != 0) {
	goto L1161;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*icp), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L1161;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ict), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L1161;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L1161;
    }
    s_wsfe(&io___45);
    do_fio(&c__1, " icp =", (ftnlen)6);
    do_fio(&c__1, (char *)&(*icp), (ftnlen)sizeof(integer));
    do_fio(&c__1, "  ict =", (ftnlen)7);
    do_fio(&c__1, (char *)&(*ict), (ftnlen)sizeof(integer));
    e_wsfe();
    ic = 0;
L1181:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___46);
    do_fio(&c__1, " Enter scaling factor to convert pressure to Pascals:", (
	    ftnlen)53);
    e_wsfe();
    i__1 = s_rsle(&io___47);
    if (i__1 != 0) {
	goto L1181;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*scp), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L1181;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L1181;
    }
    s_wsfe(&io___48);
    do_fio(&c__1, "scp =", (ftnlen)5);
    do_fio(&c__1, (char *)&(*scp), (ftnlen)sizeof(real));
    e_wsfe();

    ic = 0;
L1201:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___49);
    do_fio(&c__1, " Enter the surface temperature, ts (K): ", (ftnlen)40);
    e_wsfe();
    i__1 = s_rsle(&io___50);
    if (i__1 != 0) {
	goto L1201;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*ts), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L1201;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L1201;
    }
    s_wsfe(&io___51);
    do_fio(&c__1, "ts = ", (ftnlen)5);
    do_fio(&c__1, (char *)&(*ts), (ftnlen)sizeof(real));
    e_wsfe();

/* ****                  a b s o r b i n g    g a s e s */

    ic = 0;
L1301:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___52);
    do_fio(&c__1, " Enter the number of different absorbing gases: ", (ftnlen)
	    48);
    e_wsfe();
    i__1 = s_rsle(&io___53);
    if (i__1 != 0) {
	goto L1301;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ngases), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L1301;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L1301;
    }
    s_wsfe(&io___54);
    do_fio(&c__1, "ngases =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*ngases), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ngases > 50) {
	s_wsle(&io___55);
	do_lio(&c__9, &c__1, "Number of gases exceeds dimension bound: ngas=",
		 (ftnlen)46);
	do_lio(&c__3, &c__1, (char *)&c__50, (ftnlen)sizeof(integer));
	e_wsle();
	s_stop("", (ftnlen)0);
    }

    npabs = 0;

    i__1 = *ngases;
    for (n = 1; n <= i__1; ++n) {

/* ****       increment state structure counter */

	++(*nstate);

	ic = 0;
L1401:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___57);
	do_fio(&c__1, " Enter the HITRAN gas code for gas # ", (ftnlen)37);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	i__2 = s_rsle(&io___58);
	if (i__2 != 0) {
	    goto L1401;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&igas[n], (ftnlen)sizeof(integer))
		;
	if (i__2 != 0) {
	    goto L1401;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1401;
	}

	charsp_(cgas + (igas[n] - 1) * 6, name__, &len, &c__6, &nlb, &ntb, (
		ftnlen)6, (ftnlen)1);

	s_wsfe(&io___59);
	do_fio(&c__1, "igas =", (ftnlen)6);
	do_fio(&c__1, (char *)&igas[n], (ftnlen)sizeof(integer));
	do_fio(&c__1, " name = ", (ftnlen)8);
	i__2 = len;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
	}
	e_wsfe();
	ic = 0;
L1421:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	istate[*nstate] = 0;
	s_wsfe(&io___60);
	do_fio(&c__1, " Compute Jacobians for ", (ftnlen)23);
	i__2 = len;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
	}
	e_wsfe();
	s_wsfe(&io___61);
	do_fio(&c__1, " Enter index of Jacobian Type: ", (ftnlen)31);
	do_fio(&c__1, " 0) No Jacobians", (ftnlen)16);
	do_fio(&c__1, " 1) Radiance Jacobians", (ftnlen)22);
	do_fio(&c__1, " 2) Flux Jacobians", (ftnlen)18);
	e_wsfe();
	i__2 = s_rsle(&io___62);
	if (i__2 != 0) {
	    goto L1421;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&ist, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L1421;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1421;
	}

	if (ist != 0) {
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 9;
	    ici__1.iciunit = j_ext__ + (*nstate - 1) * 9;
	    ici__1.icifmt = "(1a3,6a)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, ".j_", (ftnlen)3);
	    i__2 = len;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
	    }
	    e_wsfi();
	    if (ist == 1) {
		istate[*nstate] = 3;
	    } else {
		istate[*nstate] = -3;
	    }
	    s_wsfe(&io___63);
	    do_fio(&c__1, "ist =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " nstate=", (ftnlen)8);
	    do_fio(&c__1, (char *)&(*nstate), (ftnlen)sizeof(integer));
	    do_fio(&c__1, " istate =", (ftnlen)9);
	    do_fio(&c__1, (char *)&istate[*nstate], (ftnlen)sizeof(integer));
	    do_fio(&c__1, " j_ext =", (ftnlen)8);
	    do_fio(&c__1, j_ext__ + (*nstate - 1) * 9, (ftnlen)9);
	    e_wsfe();

	    s_wsfe(&io___64);
	    do_fio(&c__1, " Enter fractional change in the mixing ratio (0.0"
		    "-1.0): ", (ftnlen)56);
	    e_wsfe();
	    s_rsle(&io___65);
	    do_lio(&c__4, &c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(
		    real));
	    e_rsle();
	    s_wsfe(&io___66);
	    do_fio(&c__1, "pd_frac =", (ftnlen)9);
	    do_fio(&c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(real));
	    e_wsfe();
	} else {
	    s_wsfe(&io___67);
	    do_fio(&c__1, "ist =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	    e_wsfe();
	}

	ic = 0;
L1432:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___68);
	do_fio(&c__1, " Enter  number of absorption coefficient types for ", (
		ftnlen)51);
	i__2 = len;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
	}
	e_wsfe();
	i__2 = s_rsle(&io___69);
	if (i__2 != 0) {
	    goto L1432;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&ngastyp[n], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L1432;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1432;
	}
	if (ngastyp[n] > 3) {
	    s_wsle(&io___70);
	    do_lio(&c__9, &c__1, " Error: number of absorption coeffienct ty"
		    "pes ", (ftnlen)46);
	    do_lio(&c__9, &c__1, "exceeds maximum dimension bound: ngtmax=", (
		    ftnlen)40);
	    do_lio(&c__3, &c__1, (char *)&c__3, (ftnlen)sizeof(integer));
	    e_wsle();
	    goto L1432;
	}
	s_wsfe(&io___71);
	do_fio(&c__1, "igas =", (ftnlen)6);
	do_fio(&c__1, (char *)&ngastyp[n], (ftnlen)sizeof(integer));
	e_wsfe();

	i__2 = ngastyp[n];
	for (ngt = 1; ngt <= i__2; ++ngt) {

	    ic = 0;
L1441:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___73);
	    do_fio(&c__1, " Enter index of transition type for ", (ftnlen)36);
	    do_fio(&c__1, cgas + (igas[n] - 1) * 6, (ftnlen)6);
	    do_fio(&c__1, "  (1) vibration-rotation transitions (line absorb"
		    "ers), ", (ftnlen)55);
	    do_fio(&c__1, "  (2) pressure-induced transitions (density-squar"
		    "ed), ", (ftnlen)54);
	    do_fio(&c__1, "  (3) cross-sections per molecule (UV and VIS con"
		    "tinuum)", (ftnlen)56);
	    e_wsfe();
	    i__3 = s_rsle(&io___74);
	    if (i__3 != 0) {
		goto L1441;
	    }
	    i__3 = do_lio(&c__3, &c__1, (char *)&igtrn[ngt + n * 3], (ftnlen)
		    sizeof(integer));
	    if (i__3 != 0) {
		goto L1441;
	    }
	    i__3 = e_rsle();
	    if (i__3 != 0) {
		goto L1441;
	    }
	    s_wsfe(&io___75);
	    do_fio(&c__1, "igtrn =", (ftnlen)7);
	    do_fio(&c__1, (char *)&igtrn[ngt + n * 3], (ftnlen)sizeof(integer)
		    );
	    e_wsfe();

/* ****           create a gas absorption coefficient input unit number. */

	    iuabc[ngt + n * 3] = *iugas + (n - 1) * 3 + ngt - 1;

	    ic = 0;
L1462:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___76);
	    do_fio(&c__1, " Enter name of absorption coefficient file for ", (
		    ftnlen)47);
	    i__3 = len;
	    for (j = 1; j <= i__3; ++j) {
		do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
	    }
	    e_wsfe();
	    s_rsfe(&io___77);
	    do_fio(&c__1, gasfile + (ngt + n * 3) * 132, (ftnlen)132);
	    e_rsfe();

	    s_wsfe(&io___78);
	    do_fio(&c__1, " opening gas unit:iuabc(n)=", (ftnlen)27);
	    do_fio(&c__1, (char *)&iuabc[ngt + n * 3], (ftnlen)sizeof(integer)
		    );
	    e_wsfe();
	    if (igtrn[ngt + n * 3] == 1) {
		o__1.oerr = 1;
		o__1.ounit = iuabc[ngt + n * 3];
		o__1.ofnmlen = 132;
		o__1.ofnm = gasfile + (ngt + n * 3) * 132;
		o__1.orl = 0;
		o__1.osta = "old";
		o__1.oacc = 0;
		o__1.ofm = "unformatted";
		o__1.oblnk = 0;
		i__3 = f_open(&o__1);
		if (i__3 != 0) {
		    goto L1462;
		}
	    } else {
		o__1.oerr = 1;
		o__1.ounit = iuabc[ngt + n * 3];
		o__1.ofnmlen = 132;
		o__1.ofnm = gasfile + (ngt + n * 3) * 132;
		o__1.orl = 0;
		o__1.osta = "old";
		o__1.oacc = 0;
		o__1.ofm = "formatted";
		o__1.oblnk = 0;
		i__3 = f_open(&o__1);
		if (i__3 != 0) {
		    goto L1462;
		}

		if (igtrn[ngt + n * 3] == 2) {
		    ++npabs;
		    igpab[ngt + n * 3] = npabs;
		}
	    }

	    s_wsfe(&io___79);
	    do_fio(&c__1, " Absorption coefficient file for ", (ftnlen)33);
	    do_fio(&c__1, cgas + (igas[n] - 1) * 6, (ftnlen)6);
	    do_fio(&c__1, gasfile + (ngt + n * 3) * 132, (ftnlen)132);
	    e_wsfe();

/* L1471: */
	}

	ic = 0;
L1481:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___80);
	do_fio(&c__1, " Enter gas mixing ratio file format index (1 - 3):  ", 
		(ftnlen)52);
	do_fio(&c__1, "  1) formatted 2) unformatted 3) list directed ", (
		ftnlen)47);
	e_wsfe();
	i__2 = s_rsle(&io___81);
	if (i__2 != 0) {
	    goto L1481;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&ifrmmix[n], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L1481;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1481;
	}
	s_wsfe(&io___82);
	do_fio(&c__1, "ifrmmix =", (ftnlen)9);
	do_fio(&c__1, (char *)&ifrmmix[n], (ftnlen)sizeof(integer));
	e_wsfe();

	if (ifrmmix[n] == 1) {
	    s_wsfe(&io___83);
	    do_fio(&c__1, " Enter format for reading pressure and rmix", (
		    ftnlen)43);
	    do_fio(&c__1, " [enclosed in parenthesis, ie; (2f3.5) ]: ", (
		    ftnlen)42);
	    e_wsfe();
	    s_rsfe(&io___84);
	    do_fio(&c__1, frmmix + n * 40, (ftnlen)40);
	    e_rsfe();
	    s_wsfe(&io___85);
	    do_fio(&c__1, "frmmix = ", (ftnlen)9);
	    do_fio(&c__1, frmmix + n * 40, (ftnlen)40);
	    e_wsfe();
	} else {
	    s_copy(frmmix + n * 40, " ", (ftnlen)40, (ftnlen)1);
	}

	ic = 0;
L1602:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___86);
	do_fio(&c__1, " Enter name of gas mixing ratio file for ", (ftnlen)41)
		;
	do_fio(&c__1, cgas + (igas[n] - 1) * 6, (ftnlen)6);
	e_wsfe();
	s_rsfe(&io___87);
	do_fio(&c__1, mixfile + n * 132, (ftnlen)132);
	e_rsfe();

	s_wsfe(&io___88);
	do_fio(&c__1, " Opening gas mixing ratio unit: iumix=", (ftnlen)38);
	do_fio(&c__1, (char *)&(*iumix), (ftnlen)sizeof(integer));
	e_wsfe();
	if (ifrmmix[n] != 2) {
	    o__1.oerr = 1;
	    o__1.ounit = *iumix;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = mixfile + n * 132;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = "formatted";
	    o__1.oblnk = 0;
	    i__2 = f_open(&o__1);
	    if (i__2 != 0) {
		goto L1602;
	    }
	} else {
	    o__1.oerr = 1;
	    o__1.ounit = *iumix;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = mixfile + n * 132;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = "unformatted";
	    o__1.oblnk = 0;
	    i__2 = f_open(&o__1);
	    if (i__2 != 0) {
		goto L1602;
	    }
	}

	s_wsfe(&io___89);
	i__2 = len;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
	}
	do_fio(&c__1, " mixfile(n) = ", (ftnlen)14);
	do_fio(&c__1, mixfile + n * 132, (ftnlen)132);
	e_wsfe();

	cl__1.cerr = 0;
	cl__1.cunit = *iumix;
	cl__1.csta = 0;
	f_clos(&cl__1);

	ic = 0;
L1621:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___90);
	do_fio(&c__1, " Enter number of lines to skip at top of file:", (
		ftnlen)46);
	e_wsfe();
	i__2 = s_rsle(&io___91);
	if (i__2 != 0) {
	    goto L1621;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&ioffmix[n], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L1621;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1621;
	}
	s_wsfe(&io___92);
	do_fio(&c__1, "ioffmix =", (ftnlen)9);
	do_fio(&c__1, (char *)&ioffmix[n], (ftnlen)sizeof(integer));
	e_wsfe();

	ic = 0;
L1641:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___93);
	do_fio(&c__1, " Enter type of vertical coordinate: ", (ftnlen)36);
	do_fio(&c__1, "  1) pressure,   2) altitude (km)", (ftnlen)33);
	e_wsfe();
	i__2 = s_rsle(&io___94);
	if (i__2 != 0) {
	    goto L1641;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&izmix[n], (ftnlen)sizeof(integer)
		);
	if (i__2 != 0) {
	    goto L1641;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1641;
	}
	s_wsfe(&io___95);
	do_fio(&c__1, "izmix =", (ftnlen)7);
	do_fio(&c__1, (char *)&izmix[n], (ftnlen)sizeof(integer));
	e_wsfe();
	ic = 0;
L1661:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___96);
	do_fio(&c__1, " Enter number of column containing z and rmix for ", (
		ftnlen)50);
	i__2 = len;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
	}
	e_wsfe();
	i__2 = s_rsle(&io___97);
	if (i__2 != 0) {
	    goto L1661;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&icpmix[n], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L1661;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&icmix[n], (ftnlen)sizeof(integer)
		);
	if (i__2 != 0) {
	    goto L1661;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1661;
	}
	s_wsfe(&io___98);
	do_fio(&c__1, "icpmix =", (ftnlen)8);
	do_fio(&c__1, (char *)&icpmix[n], (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___99);
	do_fio(&c__1, "icmix =", (ftnlen)7);
	do_fio(&c__1, (char *)&icmix[n], (ftnlen)sizeof(integer));
	e_wsfe();
	ic = 0;
L1681:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___100);
	do_fio(&c__1, " Enter type of mixing ratios for ", (ftnlen)33);
	do_fio(&c__1, cgas + (igas[n] - 1) * 6, (ftnlen)6);
	do_fio(&c__1, "  1) volume mixing ratio,", (ftnlen)25);
	do_fio(&c__1, "  2) mass mixing ratio", (ftnlen)22);
	e_wsfe();
	i__2 = s_rsle(&io___101);
	if (i__2 != 0) {
	    goto L1681;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&imix[n], (ftnlen)sizeof(integer))
		;
	if (i__2 != 0) {
	    goto L1681;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1681;
	}
	s_wsfe(&io___102);
	do_fio(&c__1, "imix =", (ftnlen)6);
	do_fio(&c__1, (char *)&imix[n], (ftnlen)sizeof(integer));
	e_wsfe();
	ic = 0;
L1801:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	if (izmix[n] == 1) {
	    s_wsfe(&io___103);
	    do_fio(&c__1, " Enter factors to convert pressure to Pascals ", (
		    ftnlen)46);
	    do_fio(&c__1, " and mixing ratios to range 0-1: ", (ftnlen)33);
	    e_wsfe();
	} else {
	    s_wsfe(&io___104);
	    do_fio(&c__1, " Enter multiplicative factors to convert altitude"
		    "s to km ", (ftnlen)57);
	    do_fio(&c__1, " and mixing ratios to range 0-1: ", (ftnlen)33);
	    e_wsfe();
	}
	i__2 = s_rsle(&io___105);
	if (i__2 != 0) {
	    goto L1801;
	}
	i__2 = do_lio(&c__4, &c__1, (char *)&scpmix[n], (ftnlen)sizeof(real));
	if (i__2 != 0) {
	    goto L1801;
	}
	i__2 = do_lio(&c__4, &c__1, (char *)&scmix[n], (ftnlen)sizeof(real));
	if (i__2 != 0) {
	    goto L1801;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1801;
	}
	s_wsfe(&io___106);
	do_fio(&c__1, "scpmix =", (ftnlen)8);
	do_fio(&c__1, (char *)&scpmix[n], (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___107);
	do_fio(&c__1, "scmix  =", (ftnlen)8);
	do_fio(&c__1, (char *)&scmix[n], (ftnlen)sizeof(real));
	e_wsfe();
	if (imix[n] == 1) {
	    ic = 0;
L1821:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___108);
	    do_fio(&c__1, " Enter molecular weight of gas (kg/Kmole): ", (
		    ftnlen)43);
	    e_wsfe();
	    i__2 = s_rsle(&io___109);
	    if (i__2 != 0) {
		goto L1821;
	    }
	    i__2 = do_lio(&c__4, &c__1, (char *)&wgtgas[n], (ftnlen)sizeof(
		    real));
	    if (i__2 != 0) {
		goto L1821;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L1821;
	    }
	    s_wsfe(&io___110);
	    do_fio(&c__1, "wgtgas =", (ftnlen)8);
	    do_fio(&c__1, (char *)&wgtgas[n], (ftnlen)sizeof(real));
	    e_wsfe();
	}

/* L1841: */
    }

/* ****       a e r o s o l    o p t i c a l    p r o p e r t i e s */

/* ****   enter the number of particle modes. */

    ic = 0;
L2001:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___111);
    do_fio(&c__1, " Enter number of aerosol particle modes: ", (ftnlen)41);
    e_wsfe();
    i__1 = s_rsle(&io___112);
    if (i__1 != 0) {
	goto L2001;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*nmodes), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L2001;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L2001;
    }
    s_wsfe(&io___113);
    do_fio(&c__1, "nmodes =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*nmodes), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*nmodes > 10) {
	s_wsfe(&io___114);
	do_fio(&c__1, " Error: Number of modes exceeds dimension bound: nmod"
		"e=", (ftnlen)55);
	do_fio(&c__1, (char *)&c__10, (ftnlen)sizeof(integer));
	e_wsfe();
	s_stop("", (ftnlen)0);
    }

/* ****   initialize i/o variables for each mode */

    for (m = 1; m <= 10; ++m) {
	s_copy(miefile + m * 132, " ", (ftnlen)132, (ftnlen)1);
	ioffmie[m] = 0;
	ioffmom[m] = 0;
	icwlq[m] = 0;
	icqext[m] = 0;
	icqsca[m] = 0;
	icg1[m] = 0;
	icg2[m] = 0;
	icf[m] = 0;
/* L2101: */
    }

/* ****   read in file names for aerosol optical properties */

    i__1 = *nmodes;
    for (m = 1; m <= i__1; ++m) {
	++(*nstate);
	ic = 0;
L2201:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	istate[*nstate] = 0;
	s_wsfe(&io___116);
	do_fio(&c__1, " Compute Jacobians for aerosol mode: ", (ftnlen)37);
	do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___117);
	do_fio(&c__1, " Enter index of Jacobian Type: ", (ftnlen)31);
	do_fio(&c__1, " 0) No Jacobians", (ftnlen)16);
	do_fio(&c__1, " 1) Radiance Jacobians", (ftnlen)22);
	do_fio(&c__1, " 2) Flux Jacobians", (ftnlen)18);
	e_wsfe();
	i__2 = s_rsle(&io___118);
	if (i__2 != 0) {
	    goto L2201;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&ist, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L2201;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2201;
	}

	if (ist != 0) {

	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 9;
	    ici__1.iciunit = j_ext__ + (*nstate - 1) * 9;
	    ici__1.icifmt = "(1a7,i2.2)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, ".j_aer", (ftnlen)6);
	    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	    e_wsfi();
	    if (ist == 1) {
		istate[*nstate] = 4;
	    } else {
		istate[*nstate] = -4;
	    }
	    s_wsfe(&io___119);
	    do_fio(&c__1, "ist =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " nstate=", (ftnlen)8);
	    do_fio(&c__1, (char *)&(*nstate), (ftnlen)sizeof(integer));
	    do_fio(&c__1, " istate =", (ftnlen)9);
	    do_fio(&c__1, (char *)&istate[*nstate], (ftnlen)sizeof(integer));
	    do_fio(&c__1, " j_ext =", (ftnlen)8);
	    do_fio(&c__1, j_ext__ + (*nstate - 1) * 9, (ftnlen)9);
	    e_wsfe();
	    s_wsfe(&io___120);
	    do_fio(&c__1, " Enter fractional change in optical depth (0.0-1."
		    "0): ", (ftnlen)53);
	    e_wsfe();
	    s_rsle(&io___121);
	    do_lio(&c__4, &c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(
		    real));
	    e_rsle();
	    s_wsfe(&io___122);
	    do_fio(&c__1, "pd_frac =", (ftnlen)9);
	    do_fio(&c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(real));
	    e_wsfe();
	} else {
	    s_wsfe(&io___123);
	    do_fio(&c__1, "ist =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	    e_wsfe();
	}

	ic = 0;
L2221:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___124);
	do_fio(&c__1, " Enter prefix of the mie scattering files for mode # ",
		 (ftnlen)53);
	do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	e_wsfe();
	s_rsfe(&io___125);
	do_fio(&c__1, miefile + m * 132, (ftnlen)132);
	e_rsfe();

	charsp_(miefile + m * 132, name__, &len, &c__132, &nlb, &ntb, (ftnlen)
		132, (ftnlen)1);

	s_copy(filein, " ", (ftnlen)132, (ftnlen)1);
	s_wsfi(&io___126);
	i__2 = len;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, name__ + (j - 1), (ftnlen)1);
	}
	do_fio(&c__1, ".mie", (ftnlen)4);
	e_wsfi();
	s_wsfe(&io___127);
	do_fio(&c__1, " miefile = ", (ftnlen)11);
	do_fio(&c__1, filein, (ftnlen)132);
	e_wsfe();

	s_wsfe(&io___128);
	do_fio(&c__1, " opening unit iumie=", (ftnlen)20);
	do_fio(&c__1, (char *)&(*iumie), (ftnlen)sizeof(integer));
	e_wsfe();
	o__1.oerr = 1;
	o__1.ounit = *iumie;
	o__1.ofnmlen = 132;
	o__1.ofnm = filein;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "formatted";
	o__1.oblnk = 0;
	i__2 = f_open(&o__1);
	if (i__2 != 0) {
	    goto L2221;
	}
	cl__1.cerr = 0;
	cl__1.cunit = *iumie;
	cl__1.csta = 0;
	f_clos(&cl__1);

	ic = 0;
L2241:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___129);
	do_fio(&c__1, " Enter number of records to skip above mie file: ", (
		ftnlen)49);
	e_wsfe();
	i__2 = s_rsle(&io___130);
	if (i__2 != 0) {
	    goto L2241;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&ioffmie[m], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L2241;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2241;
	}
	s_wsfe(&io___131);
	do_fio(&c__1, "ioffmie =", (ftnlen)9);
	do_fio(&c__1, (char *)&ioffmie[m], (ftnlen)sizeof(integer));
	e_wsfe();
	ic = 0;
L2261:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___132);
	do_fio(&c__1, " Enter column numbers for wl, qext, qsca, and g1 : ", (
		ftnlen)51);
	e_wsfe();
	i__2 = s_rsle(&io___133);
	if (i__2 != 0) {
	    goto L2261;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&icwlq[m], (ftnlen)sizeof(integer)
		);
	if (i__2 != 0) {
	    goto L2261;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&icqext[m], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L2261;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&icqsca[m], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L2261;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&icg1[m], (ftnlen)sizeof(integer))
		;
	if (i__2 != 0) {
	    goto L2261;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2261;
	}
	s_wsfe(&io___134);
	do_fio(&c__1, "icwlq, icqext, icqsca, icg1 =", (ftnlen)29);
	do_fio(&c__1, (char *)&icwlq[m], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&icqext[m], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&icqsca[m], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&icg1[m], (ftnlen)sizeof(integer));
	e_wsfe();

	ic = 0;
L2281:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___135);
	do_fio(&c__1, " Choose phase function type ", (ftnlen)28);
	e_wsfe();
	s_wsfe(&io___136);
	do_fio(&c__1, " 1) full phase function", (ftnlen)23);
	do_fio(&c__1, " 2) Henyey Greenstein,  ", (ftnlen)24);
	do_fio(&c__1, " 3) Double Henyey-Greenstein,  ", (ftnlen)31);
	e_wsfe();
	i__2 = s_rsle(&io___137);
	if (i__2 != 0) {
	    goto L2281;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&iang[m], (ftnlen)sizeof(integer))
		;
	if (i__2 != 0) {
	    goto L2281;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2281;
	}
	s_wsfe(&io___138);
	do_fio(&c__1, "iang =", (ftnlen)6);
	do_fio(&c__1, (char *)&iang[m], (ftnlen)sizeof(integer));
	e_wsfe();
	ioffmom[m] = 0;
	if (iang[m] == 1) {
	    ic = 0;
L2401:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___139);
	    do_fio(&c__1, " Enter number of records to skip above phase fcn "
		    "file: ", (ftnlen)55);
	    e_wsfe();
	    i__2 = s_rsle(&io___140);
	    if (i__2 != 0) {
		goto L2401;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ioffmom[m], (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L2401;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L2401;
	    }
	    s_wsfe(&io___141);
	    do_fio(&c__1, "ioffmom =", (ftnlen)9);
	    do_fio(&c__1, (char *)&ioffmom[m], (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    if (iang[m] == 3) {
		ic = 0;
L2421:
		++ic;
		if (ic > icmx) {
		    s_stop("", (ftnlen)0);
		}
		s_wsfe(&io___142);
		do_fio(&c__1, " Enter column numbers for g2 and f (0-1): ", (
			ftnlen)42);
		e_wsfe();
		i__2 = s_rsle(&io___143);
		if (i__2 != 0) {
		    goto L2421;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&icg2[m], (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L2421;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&icf[m], (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L2421;
		}
		i__2 = e_rsle();
		if (i__2 != 0) {
		    goto L2421;
		}
		s_wsfe(&io___144);
		do_fio(&c__1, "icg2, icf =", (ftnlen)11);
		do_fio(&c__1, (char *)&icg2[m], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icf[m], (ftnlen)sizeof(integer));
		e_wsfe();
	    } else {
		icg2[m] = 0;
		icf[m] = 0;
	    }
	}

/* ****       a e r o s o l    v e r t i c a l    s t r u c t u r e */

	ic = 0;
L2461:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___145);
	do_fio(&c__1, " Enter name of file with optical depths for mode: ", (
		ftnlen)50);
	do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	e_wsfe();
	s_rsfe(&io___146);
	do_fio(&c__1, aerfile + m * 132, (ftnlen)132);
	e_rsfe();
	s_wsfe(&io___147);
	do_fio(&c__1, "aerfile =", (ftnlen)9);
	do_fio(&c__1, aerfile + m * 132, (ftnlen)132);
	e_wsfe();

	s_wsfe(&io___148);
	do_fio(&c__1, " opening unit iuaer=", (ftnlen)20);
	do_fio(&c__1, (char *)&(*iuaer), (ftnlen)sizeof(integer));
	e_wsfe();
	o__1.oerr = 1;
	o__1.ounit = *iuaer;
	o__1.ofnmlen = 132;
	o__1.ofnm = aerfile + m * 132;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "formatted";
	o__1.oblnk = 0;
	i__2 = f_open(&o__1);
	if (i__2 != 0) {
	    goto L2461;
	}
	cl__1.cerr = 0;
	cl__1.cunit = *iuaer;
	cl__1.csta = 0;
	f_clos(&cl__1);

/* ****       a e r o s o l   o p t i c a l    d e p t h s */

	ic = 0;
L2601:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___149);
	do_fio(&c__1, " Enter standard wavenumber where optical", (ftnlen)40);
	do_fio(&c__1, " depths are given: ", (ftnlen)19);
	e_wsfe();
	i__2 = s_rsle(&io___150);
	if (i__2 != 0) {
	    goto L2601;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&wnaer[m], (ftnlen)sizeof(
		doublereal));
	if (i__2 != 0) {
	    goto L2601;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2601;
	}
	s_wsfe(&io___151);
	do_fio(&c__1, "wnaer =", (ftnlen)7);
	do_fio(&c__1, (char *)&wnaer[m], (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___152);
	do_fio(&c__1, " Enter number of records to skip above optical depths"
		": ", (ftnlen)55);
	e_wsfe();
	s_rsle(&io___153);
	do_lio(&c__3, &c__1, (char *)&iofftau[m], (ftnlen)sizeof(integer));
	e_rsle();
	s_wsfe(&io___154);
	do_fio(&c__1, "iofftau =", (ftnlen)9);
	do_fio(&c__1, (char *)&iofftau[m], (ftnlen)sizeof(integer));
	e_wsfe();
	ic = 0;
L2621:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___155);
	do_fio(&c__1, " Enter type of vertical coordinate: ", (ftnlen)36);
	do_fio(&c__1, "  1) pressure,   2) altitude (km)", (ftnlen)33);
	e_wsfe();
	i__2 = s_rsle(&io___156);
	if (i__2 != 0) {
	    goto L2621;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&iztau[m], (ftnlen)sizeof(integer)
		);
	if (i__2 != 0) {
	    goto L2621;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2621;
	}
	s_wsfe(&io___157);
	do_fio(&c__1, "iztau =", (ftnlen)7);
	do_fio(&c__1, (char *)&iztau[m], (ftnlen)sizeof(integer));
	e_wsfe();
	ic = 0;
L2641:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___158);
	do_fio(&c__1, " Enter column numbers with vertical coordinate ", (
		ftnlen)47);
	do_fio(&c__1, " and differential optical depth:", (ftnlen)32);
	e_wsfe();
	i__2 = s_rsle(&io___159);
	if (i__2 != 0) {
	    goto L2641;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&icptau[m], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L2641;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&ictau[m], (ftnlen)sizeof(integer)
		);
	if (i__2 != 0) {
	    goto L2641;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2641;
	}
	s_wsfe(&io___160);
	do_fio(&c__1, "icptau, ictau =", (ftnlen)15);
	do_fio(&c__1, (char *)&icptau[m], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ictau[m], (ftnlen)sizeof(integer));
	e_wsfe();
	ic = 0;
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	if (iztau[m] == 1) {
	    s_wsfe(&io___161);
	    do_fio(&c__1, " Enter multiplicative factors to convert", (ftnlen)
		    40);
	    do_fio(&c__1, " pressure to Pascals and to scale optical depth: ",
		     (ftnlen)49);
	    e_wsfe();
	} else {
	    s_wsfe(&io___162);
	    do_fio(&c__1, " Enter multiplicative factors to convert", (ftnlen)
		    40);
	    do_fio(&c__1, " altitudes to km and to scale optical depth: ", (
		    ftnlen)45);
	    e_wsfe();
	}
	scptau[m] = 1.f;
	sctau[m] = 1.f;
	i__2 = s_rsle(&io___163);
	if (i__2 != 0) {
	    goto L2661;
	}
	i__2 = do_lio(&c__4, &c__1, (char *)&scptau[m], (ftnlen)sizeof(real));
	if (i__2 != 0) {
	    goto L2661;
	}
	i__2 = do_lio(&c__4, &c__1, (char *)&sctau0, (ftnlen)sizeof(real));
	if (i__2 != 0) {
	    goto L2661;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2661;
	}
	sctau[m] = sctau0;
L2661:
	s_wsfe(&io___165);
	do_fio(&c__1, "scptau =", (ftnlen)8);
	do_fio(&c__1, (char *)&scptau[m], (ftnlen)sizeof(real));
	do_fio(&c__1, "sctau =", (ftnlen)7);
	do_fio(&c__1, (char *)&sctau[m], (ftnlen)sizeof(real));
	e_wsfe();

/* ****       a e r o s o l    s c a l e    h e i g h t s */

	ic = 0;
L2801:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___166);
	do_fio(&c__1, " Are aerosol scale heights specified for", (ftnlen)40);
	do_fio(&c__1, " each layer (true/false)? ", (ftnlen)26);
	e_wsfe();
	i__2 = s_rsle(&io___167);
	if (i__2 != 0) {
	    goto L2801;
	}
	i__2 = do_lio(&c__8, &c__1, (char *)&lcsh[m], (ftnlen)sizeof(logical))
		;
	if (i__2 != 0) {
	    goto L2801;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L2801;
	}
	s_wsfe(&io___168);
	do_fio(&c__1, " lcsh =", (ftnlen)7);
	do_fio(&c__1, (char *)&lcsh[m], (ftnlen)sizeof(logical));
	e_wsfe();
	if (lcsh[m]) {
	    ic = 0;
L2821:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___169);
	    do_fio(&c__1, " Enter column numbers with aerosol scale height:", 
		    (ftnlen)48);
	    e_wsfe();
	    i__2 = s_rsle(&io___170);
	    if (i__2 != 0) {
		goto L2821;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&iccsh[m], (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L2821;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L2821;
	    }
	    s_wsfe(&io___171);
	    do_fio(&c__1, "iccsh =", (ftnlen)7);
	    do_fio(&c__1, (char *)&iccsh[m], (ftnlen)sizeof(integer));
	    e_wsfe();
	    ic = 0;
L2841:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___172);
	    do_fio(&c__1, " Enter factor to convert scale heights to km: ", (
		    ftnlen)46);
	    e_wsfe();
	    i__2 = s_rsle(&io___173);
	    if (i__2 != 0) {
		goto L2841;
	    }
	    i__2 = do_lio(&c__4, &c__1, (char *)&sccsh[m], (ftnlen)sizeof(
		    real));
	    if (i__2 != 0) {
		goto L2841;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L2841;
	    }
	    s_wsfe(&io___174);
	    do_fio(&c__1, "sccsh =", (ftnlen)7);
	    do_fio(&c__1, (char *)&sccsh[m], (ftnlen)sizeof(real));
	    e_wsfe();
	}

/* L2861: */
    }

/* ****       s u r f a c e    o p t i c a l    p r o p e r t i e s */

    ic = 0;
L3001:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___175);
    do_fio(&c__1, " Enter the index of the surface BRDF type: ", (ftnlen)43);
    do_fio(&c__1, " (0) Lambertian, (1) Non-Lambertian: ", (ftnlen)37);
    e_wsfe();
    i__1 = s_rsle(&io___176);
    if (i__1 != 0) {
	goto L3001;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&isurtyp, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L3001;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3001;
    }
    if (isurtyp == 0) {
	*lamber = TRUE_;
    } else {
	*lamber = FALSE_;
	s_wsfe(&io___178);
	do_fio(&c__1, " Enter the index of the surface BRDF mode: ", (ftnlen)
		43);
	do_fio(&c__1, "  1) Hapke BDR model", (ftnlen)20);
	do_fio(&c__1, "  2) Breon BDR model; combination of Li + Roujean", (
		ftnlen)49);
	do_fio(&c__1, "  3) Roujean BDR model", (ftnlen)22);
	do_fio(&c__1, "  4) Cox and Munk glint model", (ftnlen)29);
	e_wsfe();
	i__1 = s_rsle(&io___179);
	if (i__1 != 0) {
	    goto L3001;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*iref), (ftnlen)sizeof(integer))
		;
	if (i__1 != 0) {
	    goto L3001;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L3001;
	}
	s_wsfe(&io___180);
	do_fio(&c__1, "iref =", (ftnlen)6);
	do_fio(&c__1, (char *)&(*iref), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    s_wsfe(&io___181);
    do_fio(&c__1, "lamber =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*lamber), (ftnlen)sizeof(logical));
    e_wsfe();

    ++(*nstate);
    ic = 0;
L3011:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    istate[*nstate] = 0;
    s_wsfe(&io___182);
    do_fio(&c__1, " Compute Jacobians for Surface Reflectance?", (ftnlen)43);
    e_wsfe();
    s_wsfe(&io___183);
    do_fio(&c__1, " Enter index of Jacobian Type: ", (ftnlen)31);
    do_fio(&c__1, " 0) No Jacobians", (ftnlen)16);
    do_fio(&c__1, " 1) Radiance Jacobians", (ftnlen)22);
    do_fio(&c__1, " 2) Flux Jacobians", (ftnlen)18);
    e_wsfe();
    i__1 = s_rsle(&io___184);
    if (i__1 != 0) {
	goto L3011;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&ist, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L3011;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3011;
    }

    if (ist != 0) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 9;
	ici__1.iciunit = j_ext__ + (*nstate - 1) * 9;
	ici__1.icifmt = "(1a7)";
	s_wsfi(&ici__1);
	do_fio(&c__1, ".j_surf", (ftnlen)7);
	e_wsfi();
	if (ist == 1) {
	    istate[*nstate] = 5;
	} else {
	    istate[*nstate] = -5;
	}
	s_wsfe(&io___185);
	do_fio(&c__1, "ist =", (ftnlen)5);
	do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	do_fio(&c__1, " nstate=", (ftnlen)8);
	do_fio(&c__1, (char *)&(*nstate), (ftnlen)sizeof(integer));
	do_fio(&c__1, " istate =", (ftnlen)9);
	do_fio(&c__1, (char *)&istate[*nstate], (ftnlen)sizeof(integer));
	do_fio(&c__1, " j_ext =", (ftnlen)8);
	do_fio(&c__1, j_ext__ + (*nstate - 1) * 9, (ftnlen)9);
	e_wsfe();

L3021:
	s_wsfe(&io___186);
	do_fio(&c__1, " Enter the fractional change in reflectance (0.0-1.0)"
		": ", (ftnlen)55);
	e_wsfe();
	i__1 = s_rsle(&io___187);
	if (i__1 != 0) {
	    goto L3021;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&pd_frac__[*nstate], (ftnlen)
		sizeof(real));
	if (i__1 != 0) {
	    goto L3021;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L3021;
	}
	s_wsfe(&io___188);
	do_fio(&c__1, "pd_frac =", (ftnlen)9);
	do_fio(&c__1, (char *)&pd_frac__[*nstate], (ftnlen)sizeof(real));
	e_wsfe();
    } else {
	s_wsfe(&io___189);
	do_fio(&c__1, "ist =", (ftnlen)5);
	do_fio(&c__1, (char *)&ist, (ftnlen)sizeof(integer));
	e_wsfe();
    }

    ic = 0;
L3031:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___190);
    do_fio(&c__1, " Enter the surface reflectance file format index (1 - 3): "
	    , (ftnlen)58);
    do_fio(&c__1, "  1) formatted 2) unformatted 3) list directed ", (ftnlen)
	    47);
    e_wsfe();
    i__1 = s_rsle(&io___191);
    if (i__1 != 0) {
	goto L3031;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ifrmsur), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L3031;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3031;
    }
    s_wsfe(&io___192);
    do_fio(&c__1, "ifrmsur =", (ftnlen)9);
    do_fio(&c__1, (char *)&(*ifrmsur), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*ifrmsur == 1) {
	s_wsfe(&io___193);
	do_fio(&c__1, " Enter format for reading surface properties", (ftnlen)
		44);
	do_fio(&c__1, " [enclosed in parenthesis, ie; (2f3.5) ]: ", (ftnlen)
		42);
	e_wsfe();
	s_rsfe(&io___194);
	do_fio(&c__1, frmsur, (ftnlen)40);
	e_rsfe();
	s_wsfe(&io___195);
	do_fio(&c__1, "frmsur =", (ftnlen)8);
	do_fio(&c__1, frmsur, (ftnlen)40);
	e_wsfe();
    } else {
	s_copy(frmsur, " ", (ftnlen)40, (ftnlen)1);
    }

    ic = 0;
L3042:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___196);
    do_fio(&c__1, " Enter name surface optical property file: ", (ftnlen)43);
    e_wsfe();
    s_rsfe(&io___197);
    do_fio(&c__1, surfile, (ftnlen)132);
    e_rsfe();

    s_wsfe(&io___198);
    do_fio(&c__1, " opening unit iusur=", (ftnlen)20);
    do_fio(&c__1, (char *)&(*iusur), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ifrmsur != 2) {
	o__1.oerr = 1;
	o__1.ounit = *iusur;
	o__1.ofnmlen = 132;
	o__1.ofnm = surfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "formatted";
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L3042;
	}
    } else {
	o__1.oerr = 1;
	o__1.ounit = *iusur;
	o__1.ofnmlen = 132;
	o__1.ofnm = surfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L3042;
	}
    }

    cl__1.cerr = 0;
    cl__1.cunit = *iusur;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_wsfe(&io___199);
    do_fio(&c__1, " surfile =", (ftnlen)10);
    do_fio(&c__1, surfile, (ftnlen)132);
    e_wsfe();

    ic = 0;
L3061:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___200);
    do_fio(&c__1, " Enter number of lines to skip at top of this file:", (
	    ftnlen)51);
    e_wsfe();
    i__1 = s_rsle(&io___201);
    if (i__1 != 0) {
	goto L3061;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ioffsur), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L3061;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3061;
    }
    s_wsfe(&io___202);
    do_fio(&c__1, "ioffsur =", (ftnlen)9);
    do_fio(&c__1, (char *)&(*ioffsur), (ftnlen)sizeof(integer));
    e_wsfe();

/* ****    initialize the number of surface reflectance properties */
/*        at each wavelength used in the subroutine bdrf */
/*             iref= 0, nref = 1: Lambert albedo */
/*             iref= 1, nref = 1: Hapke : single scatter albedo, W, and */
/*                                angular width factor HH */
/*             iref= 2, nref = 3: - Breon's BDR model: k0, k1, k2 */
/*             iref= 3, nref = 3: - Roujean's BDR model: k0, k1, k2 */
/*             iref= 4, nref = 2: - Cox and Munk glint model: n, k */

    if (*iref == 0) {
	*nref = 1;
    } else {
	if (*iref == 1) {
	    *nref = 2;
	} else {
	    if (*iref == 2 || *iref == 3) {
		*nref = 3;
	    } else {
		*nref = 2;
	    }
	}
    }
    ic = 0;
L3081:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___203);
    do_fio(&c__1, " Enter columns with spectral quantity and each", (ftnlen)
	    46);
    do_fio(&c__1, " of the surface optical properties:", (ftnlen)35);
    e_wsfe();
    i__1 = s_rsle(&io___204);
    if (i__1 != 0) {
	goto L3081;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*icwnsur), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L3081;
    }
    i__2 = *nref;
    for (ir = 1; ir <= i__2; ++ir) {
	i__1 = do_lio(&c__3, &c__1, (char *)&icalb[ir], (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L3081;
	}
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3081;
    }
    s_wsfe(&io___206);
    do_fio(&c__1, "icwnsur, icalb =", (ftnlen)16);
    do_fio(&c__1, (char *)&(*icwnsur), (ftnlen)sizeof(integer));
    i__1 = *nref;
    for (ir = 1; ir <= i__1; ++ir) {
	do_fio(&c__1, (char *)&icalb[ir], (ftnlen)sizeof(integer));
    }
    e_wsfe();

    ic = 0;
L3221:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___207);
    do_fio(&c__1, " Enter type of spectral quantity:", (ftnlen)33);
    do_fio(&c__1, "  1) wavelength,  2) wavenumber", (ftnlen)31);
    e_wsfe();
    i__1 = s_rsle(&io___208);
    if (i__1 != 0) {
	goto L3221;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*iwnsur), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L3221;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3221;
    }
    s_wsfe(&io___209);
    do_fio(&c__1, "iwnsur =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*iwnsur), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*iwnsur == 1) {
	ic = 0;
L3241:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___210);
	do_fio(&c__1, " Enter factor needed to convert wavelengths to micron"
		"s: ", (ftnlen)56);
	e_wsfe();
	i__1 = s_rsle(&io___211);
	if (i__1 != 0) {
	    goto L3241;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*scwalb), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L3241;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L3241;
	}
	s_wsfe(&io___212);
	do_fio(&c__1, "scwalb =", (ftnlen)8);
	do_fio(&c__1, (char *)&(*scwalb), (ftnlen)sizeof(real));
	e_wsfe();
    } else {
	*scwalb = 1.f;
    }

    ic = 0;
L3261:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___213);
    do_fio(&c__1, " Enter factor needed to scale surface optical properties: "
	    , (ftnlen)58);
    e_wsfe();
    i__1 = s_rsle(&io___214);
    if (i__1 != 0) {
	goto L3261;
    }
    i__2 = *nref;
    for (ir = 1; ir <= i__2; ++ir) {
	i__1 = do_lio(&c__4, &c__1, (char *)&scalb[ir], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L3261;
	}
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3261;
    }
    s_wsfe(&io___215);
    do_fio(&c__1, "scalb =", (ftnlen)7);
    i__1 = *nref;
    for (ir = 1; ir <= i__1; ++ir) {
	do_fio(&c__1, (char *)&scalb[ir], (ftnlen)sizeof(real));
    }
    e_wsfe();

    if (*iref == 4) {

/* ****    set wind speed and direction for Cox/Munk model */

L3301:
	s_wsfe(&io___216);
	do_fio(&c__1, "Enter the wind speed (m/s) and azimuth (deg): ", (
		ftnlen)46);
	e_wsfe();
	i__1 = s_rsle(&io___217);
	if (i__1 != 0) {
	    goto L3301;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*ws), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L3301;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*phiw), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L3301;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L3301;
	}
	s_wsfe(&io___218);
	do_fio(&c__1, "Wind speed =", (ftnlen)12);
	do_fio(&c__1, (char *)&(*ws), (ftnlen)sizeof(real));
	do_fio(&c__1, " m/s   Wind Azimuth =", (ftnlen)21);
	do_fio(&c__1, (char *)&(*phiw), (ftnlen)sizeof(real));
	do_fio(&c__1, " deg", (ftnlen)4);
	e_wsfe();
    } else {
	*ws = 0.f;
	*phiw = 0.f;
    }

/* ****  read physical properties of planet and atmosphere. */

    ic = 0;
L3401:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___219);
    do_fio(&c__1, " Enter the planets distance from sun (au):", (ftnlen)42);
    e_wsfe();
    i__1 = s_rsle(&io___220);
    if (i__1 != 0) {
	goto L3401;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*au), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L3401;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3401;
    }
    s_wsfe(&io___221);
    do_fio(&c__1, "au =", (ftnlen)4);
    do_fio(&c__1, (char *)&(*au), (ftnlen)sizeof(real));
    e_wsfe();

    ic = 0;
L3441:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___222);
    do_fio(&c__1, " Enter gravitational acceleration at surface (m/s**2): ", (
	    ftnlen)55);
    e_wsfe();
    i__1 = s_rsle(&io___223);
    if (i__1 != 0) {
	goto L3441;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*sgrav), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L3441;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3441;
    }
    s_wsfe(&io___224);
    do_fio(&c__1, "sgrav =", (ftnlen)7);
    do_fio(&c__1, (char *)&(*sgrav), (ftnlen)sizeof(real));
    e_wsfe();

    ic = 0;
L3461:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___225);
    do_fio(&c__1, " Enter the radius of the planet (km): ", (ftnlen)38);
    e_wsfe();
    i__1 = s_rsle(&io___226);
    if (i__1 != 0) {
	goto L3461;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*radius), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L3461;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3461;
    }
    s_wsfe(&io___227);
    do_fio(&c__1, "radius =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*radius), (ftnlen)sizeof(real));
    e_wsfe();

    ic = 0;
L3481:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___228);
    do_fio(&c__1, " Enter the atmospheric mean molecular weight (kg/kmole): ",
	     (ftnlen)57);
    e_wsfe();
    i__1 = s_rsle(&io___229);
    if (i__1 != 0) {
	goto L3481;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*wgtatm), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L3481;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3481;
    }
    s_wsfe(&io___230);
    do_fio(&c__1, "wgtatm =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*wgtatm), (ftnlen)sizeof(real));
    e_wsfe();

    ic = 0;
L3601:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___231);
    do_fio(&c__1, " Enter the number of major atmospheric components:", (
	    ftnlen)50);
    e_wsfe();
    i__1 = s_rsle(&io___232);
    if (i__1 != 0) {
	goto L3601;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ncomp), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L3601;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3601;
    }
    s_wsfe(&io___233);
    do_fio(&c__1, "ncomp =", (ftnlen)7);
    do_fio(&c__1, (char *)&(*ncomp), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ncomp > 50) {
	s_wsle(&io___234);
	do_lio(&c__9, &c__1, "ncomp exceeds dimension bound, ngas=", (ftnlen)
		36);
	do_lio(&c__3, &c__1, (char *)&c__50, (ftnlen)sizeof(integer));
	e_wsle();
	s_stop("", (ftnlen)0);
    }

    i__1 = *ncomp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfe(&io___236);
	do_fio(&c__1, " component types: (1) air  (2) co2  (3) n2  (4) o2", (
		ftnlen)50);
	e_wsfe();
	ic = 0;
L3621:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___237);
	do_fio(&c__1, " Enter the index of component # ", (ftnlen)32);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	i__2 = s_rsle(&io___238);
	if (i__2 != 0) {
	    goto L3621;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&icomp[i__], (ftnlen)sizeof(
		integer));
	if (i__2 != 0) {
	    goto L3621;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L3621;
	}
	s_wsfe(&io___239);
	do_fio(&c__1, "icomp =", (ftnlen)7);
	do_fio(&c__1, (char *)&icomp[i__], (ftnlen)sizeof(integer));
	e_wsfe();
	ic = 0;
L3641:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___240);
	do_fio(&c__1, " Enter the volume mixing ratio of component ", (ftnlen)
		44);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	i__2 = s_rsle(&io___241);
	if (i__2 != 0) {
	    goto L3641;
	}
	i__2 = do_lio(&c__4, &c__1, (char *)&volmix[i__], (ftnlen)sizeof(real)
		);
	if (i__2 != 0) {
	    goto L3641;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L3641;
	}
	s_wsfe(&io___242);
	do_fio(&c__1, "volmix =", (ftnlen)8);
	do_fio(&c__1, (char *)&volmix[i__], (ftnlen)sizeof(real));
	e_wsfe();
/* L3661: */
    }

/* ****    d i s c r e t e    o r d i n a t e    m e t h o d */

    ic = 0;
L3802:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___243);
    do_fio(&c__1, " Enter number of streams for discrete ordinate method (> "
	    "2): ", (ftnlen)61);
    do_fio(&c__1, " [half negative (downward) and half positive (upward)]", (
	    ftnlen)54);
    e_wsfe();
    i__1 = s_rsle(&io___244);
    if (i__1 != 0) {
	goto L3802;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*nstr), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L3802;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L3802;
    }
    s_wsfe(&io___245);
    do_fio(&c__1, "nstr =", (ftnlen)6);
    do_fio(&c__1, (char *)&(*nstr), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*nstr > 16) {
	s_wsfe(&io___246);
	do_fio(&c__1, " number of streams exceeds dimension bound: mxumu=", (
		ftnlen)50);
	do_fio(&c__1, (char *)&c__16, (ftnlen)sizeof(integer));
	e_wsfe();
	s_stop("", (ftnlen)0);
    }

/* ****   s o u r c e    f u n c t i o n s */

    ic = 0;
L4001:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___247);
    do_fio(&c__1, " Choose types of source functions: ", (ftnlen)35);
    do_fio(&c__1, "   1) direct solar beam", (ftnlen)23);
    do_fio(&c__1, "   2) internal thermal sources", (ftnlen)30);
    do_fio(&c__1, "   3) both solar and thermal sources", (ftnlen)36);
    do_fio(&c__1, " Select index of source type(s): ", (ftnlen)33);
    e_wsfe();
    i__1 = s_rsle(&io___248);
    if (i__1 != 0) {
	goto L4001;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*isource), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L4001;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L4001;
    }
    s_wsfe(&io___249);
    do_fio(&c__1, "isource =", (ftnlen)9);
    do_fio(&c__1, (char *)&(*isource), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*isource == 1 || *isource == 3) {
	if (*isource == 3) {
	    *lplanck = TRUE_;
	}

/* ****     s o l a r    f l u x e s */

	*lsolar = TRUE_;
	ic = 0;
L4041:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___250);
	do_fio(&c__1, " Enter the solar flux file format index (1 - 3):  ", (
		ftnlen)50);
	do_fio(&c__1, "  1) formatted 2) unformatted 3) list directed ", (
		ftnlen)47);
	e_wsfe();
	i__1 = s_rsle(&io___251);
	if (i__1 != 0) {
	    goto L4041;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*ifrms0), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L4041;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L4041;
	}
	s_wsfe(&io___252);
	do_fio(&c__1, "ifrms0 =", (ftnlen)8);
	do_fio(&c__1, (char *)&(*ifrms0), (ftnlen)sizeof(integer));
	e_wsfe();

	if (*ifrms0 == 1) {
	    s_wsfe(&io___253);
	    do_fio(&c__1, " Enter format for reading solar fluxes", (ftnlen)
		    38);
	    do_fio(&c__1, " [enclosed in parenthesis, ie; (2f3.5) ]: ", (
		    ftnlen)42);
	    e_wsfe();
	    s_rsfe(&io___254);
	    do_fio(&c__1, frms0, (ftnlen)40);
	    e_rsfe();
	    s_wsfe(&io___255);
	    do_fio(&c__1, "frms0 =", (ftnlen)7);
	    do_fio(&c__1, frms0, (ftnlen)40);
	    e_wsfe();
	} else {
	    s_copy(frms0, " ", (ftnlen)40, (ftnlen)1);
	}

	ic = 0;
L4062:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___256);
	do_fio(&c__1, " Enter name of solar flux file: ", (ftnlen)32);
	e_wsfe();
	s_rsfe(&io___257);
	do_fio(&c__1, solfile, (ftnlen)132);
	e_rsfe();
	s_wsfe(&io___258);
	do_fio(&c__1, "solfile =", (ftnlen)9);
	do_fio(&c__1, solfile, (ftnlen)132);
	e_wsfe();

	s_wsfe(&io___259);
	do_fio(&c__1, " opening unit, iusol1=", (ftnlen)22);
	do_fio(&c__1, (char *)&(*iusol1), (ftnlen)sizeof(integer));
	e_wsfe();
	if (*ifrms0 != 2) {
	    o__1.oerr = 1;
	    o__1.ounit = *iusol1;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = solfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = "formatted";
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L4062;
	    }
	} else {
	    o__1.oerr = 1;
	    o__1.ounit = *iusol1;
	    o__1.ofnmlen = 132;
	    o__1.ofnm = solfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = "unformatted";
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L4062;
	    }
	}

	cl__1.cerr = 0;
	cl__1.cunit = *iusol1;
	cl__1.csta = 0;
	f_clos(&cl__1);

	ic = 0;
L4081:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___260);
	do_fio(&c__1, " Enter number of lines to skip at top of file:", (
		ftnlen)46);
	e_wsfe();
	i__1 = s_rsle(&io___261);
	if (i__1 != 0) {
	    goto L4081;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*ioffs0), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L4081;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L4081;
	}
	s_wsfe(&io___262);
	do_fio(&c__1, "ioffs0 =", (ftnlen)8);
	do_fio(&c__1, (char *)&(*ioffs0), (ftnlen)sizeof(integer));
	e_wsfe();

	ic = 0;
L4201:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___263);
	do_fio(&c__1, " Enter input radiance units of input solar fluxes: ", (
		ftnlen)51);
	do_fio(&c__1, "  1) Watts/m**2/cm**-1        ", (ftnlen)30);
	do_fio(&c__1, "  2) Watts/m**2/micron        ", (ftnlen)30);
	do_fio(&c__1, "  3) Watts/m**2/nanometer     ", (ftnlen)30);
	do_fio(&c__1, "  4) Watts/m**2/Angstrom      ", (ftnlen)30);
	e_wsfe();
	i__1 = s_rsle(&io___264);
	if (i__1 != 0) {
	    goto L4201;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*iuins0), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L4201;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L4201;
	}
	s_wsfe(&io___265);
	do_fio(&c__1, "iuins0 =", (ftnlen)8);
	do_fio(&c__1, (char *)&(*iuins0), (ftnlen)sizeof(integer));
	e_wsfe();

	ic = 0;
L4221:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___266);
	do_fio(&c__1, " Specify type of abscissa for solar spectrum: ", (
		ftnlen)46);
	do_fio(&c__1, "  1) wavelength (microns)", (ftnlen)25);
	do_fio(&c__1, "  2) wavenumber", (ftnlen)15);
	do_fio(&c__1, " Enter index of choice: ", (ftnlen)24);
	e_wsfe();
	i__1 = s_rsle(&io___267);
	if (i__1 != 0) {
	    goto L4221;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*ixs0), (ftnlen)sizeof(integer))
		;
	if (i__1 != 0) {
	    goto L4221;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L4221;
	}
	s_wsfe(&io___268);
	do_fio(&c__1, "ixs0 =", (ftnlen)6);
	do_fio(&c__1, (char *)&(*ixs0), (ftnlen)sizeof(integer));
	e_wsfe();

	ic = 0;
	if (*ixs0 == 1) {
L4241:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___269);
	    do_fio(&c__1, " Enter factor to convert wavelengths to microns:", 
		    (ftnlen)48);
	    e_wsfe();
	    i__1 = s_rsle(&io___270);
	    if (i__1 != 0) {
		goto L4241;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*scwns0), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L4241;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L4241;
	    }
	    s_wsfe(&io___271);
	    do_fio(&c__1, "scwns0 =", (ftnlen)8);
	    do_fio(&c__1, (char *)&(*scwns0), (ftnlen)sizeof(real));
	    e_wsfe();
	    ic = 0;
L4261:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___272);
	    do_fio(&c__1, " Enter column numbers of wavelength and solar flu"
		    "x: ", (ftnlen)52);
	    e_wsfe();
	    i__1 = s_rsle(&io___273);
	    if (i__1 != 0) {
		goto L4261;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*icwns0), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L4261;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*icsol), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L4261;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L4261;
	    }
	    s_wsfe(&io___274);
	    do_fio(&c__1, "icwns0, icsol =", (ftnlen)15);
	    do_fio(&c__1, (char *)&(*icwns0), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*icsol), (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
L4281:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___275);
	    do_fio(&c__1, " Enter column numbers of wavenumber and solar flu"
		    "x: ", (ftnlen)52);
	    e_wsfe();
	    i__1 = s_rsle(&io___276);
	    if (i__1 != 0) {
		goto L4281;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*icwns0), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L4281;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*icsol), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L4281;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L4281;
	    }
	    s_wsfe(&io___277);
	    do_fio(&c__1, "icwns0, icsol =", (ftnlen)15);
	    do_fio(&c__1, (char *)&(*icwns0), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*icsol), (ftnlen)sizeof(integer));
	    e_wsfe();
	    *scwns0 = 1.f;
	}
	*scsol = 1.f;

	ic = 0;
L4402:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___278);
	do_fio(&c__1, " Enter number of solar incidence angles (> 0): ", (
		ftnlen)47);
	e_wsfe();
	i__1 = s_rsle(&io___279);
	if (i__1 != 0) {
	    goto L4402;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*nza), (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L4402;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L4402;
	}
	s_wsfe(&io___280);
	do_fio(&c__1, "nza =", (ftnlen)5);
	do_fio(&c__1, (char *)&(*nza), (ftnlen)sizeof(integer));
	e_wsfe();

	if (*nza > 10) {
	    s_wsfe(&io___281);
	    do_fio(&c__1, " number solar zenith angles exceeds dimension", (
		    ftnlen)45);
	    do_fio(&c__1, " bound: nsol=", (ftnlen)13);
	    do_fio(&c__1, (char *)&c__10, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_stop("", (ftnlen)0);
	}

	s_wsfe(&io___282);
	do_fio(&c__1, " Enter solar zenith and azimuth angles (degrees)", (
		ftnlen)48);
	do_fio(&c__1, " for each in incidence angle: ", (ftnlen)30);
	do_fio(&c__1, " NOTE:  Solar zenith angles less than ~2 degrees and ",
		 (ftnlen)53);
	do_fio(&c__1, "        angles between ~58 and 62 degrees often cause "
		, (ftnlen)54);
	do_fio(&c__1, "        near-singular matrices in the discrete ordina"
		"te model", (ftnlen)61);
	e_wsfe();
	i__1 = *nza;
	for (n = 1; n <= i__1; ++n) {
	    ic = 0;
L4421:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___283);
	    do_fio(&c__1, " Enter solar zenith and azimuth angles for path # "
		    , (ftnlen)50);
	    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	    e_wsfe();
	    i__2 = s_rsle(&io___284);
	    if (i__2 != 0) {
		goto L4421;
	    }
	    i__2 = do_lio(&c__4, &c__1, (char *)&sza0[n - 1], (ftnlen)sizeof(
		    real));
	    if (i__2 != 0) {
		goto L4421;
	    }
	    i__2 = do_lio(&c__4, &c__1, (char *)&phi0[n], (ftnlen)sizeof(real)
		    );
	    if (i__2 != 0) {
		goto L4421;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L4421;
	    }
	    s_wsfe(&io___286);
	    do_fio(&c__1, "sza0, phi0", (ftnlen)10);
	    do_fio(&c__1, (char *)&sza0[n - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&phi0[n], (ftnlen)sizeof(real));
	    e_wsfe();
	    umu0[n] = cos(pi * sza0[n - 1] / 180.f);
/* L4441: */
	}

	s_wsfe(&io___287);
	do_fio(&c__1, " Enter convergence criteria for azimuthal series (0.0"
		"-0.01):", (ftnlen)60);
	e_wsfe();
	s_rsle(&io___288);
	do_lio(&c__4, &c__1, (char *)&(*accur), (ftnlen)sizeof(real));
	e_rsle();
	s_wsfe(&io___289);
	do_fio(&c__1, "accur =", (ftnlen)7);
	do_fio(&c__1, (char *)&(*accur), (ftnlen)sizeof(real));
	e_wsfe();
    } else {
	*lplanck = TRUE_;
	*lsolar = FALSE_;
	*nza = 1;
    }

    ic = 0;
L4621:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___290);
    do_fio(&c__1, " Specify quantities to be included in output files: ", (
	    ftnlen)52);
    do_fio(&c__1, "   1) wavelength-dependent fluxes and radiances at the", (
	    ftnlen)54);
    do_fio(&c__1, "      computational azimuths and zenith angles at the ", (
	    ftnlen)54);
    do_fio(&c__1, "      specified output levels, and spectrally-integrated ",
	     (ftnlen)57);
    do_fio(&c__1, "      fluxes and heating rates at each computational leve"
	    "l.  ", (ftnlen)61);
    do_fio(&c__1, "   2) wavelength-dependent fluxes, radiances, and ", (
	    ftnlen)50);
    do_fio(&c__1, "      transmission values at computational zenith angles ",
	     (ftnlen)57);
    do_fio(&c__1, "      and specified output levels, and spectrally- ", (
	    ftnlen)51);
    do_fio(&c__1, "      integrated fluxes and heating rates at each ", (
	    ftnlen)50);
    do_fio(&c__1, "      computational level.  ", (ftnlen)28);
    do_fio(&c__1, "   3) wavelength-dependent fluxes and radiances at the ", (
	    ftnlen)55);
    do_fio(&c__1, "      computational azimuths and zenith angles at the ", (
	    ftnlen)54);
    do_fio(&c__1, "      specified output levels, wavelength-dependent,", (
	    ftnlen)52);
    do_fio(&c__1, "      level-dependent, pressure-weighted flux divergences"
	    ", ", (ftnlen)59);
    do_fio(&c__1, "      and spectrally-integrated fluxes and heating rates ",
	     (ftnlen)57);
    do_fio(&c__1, "      at each computational level. ", (ftnlen)35);
    e_wsfe();
    s_wsfe(&io___291);
    do_fio(&c__1, "   4) wavelength-dependent fluxes, radiances, ", (ftnlen)
	    46);
    do_fio(&c__1, "      transmission values, wavelength-dependent, ", (
	    ftnlen)49);
    do_fio(&c__1, "      level-dependent, pressure-weighted flux divergences"
	    ", ", (ftnlen)59);
    do_fio(&c__1, "      and spectrally-integrated fluxes and heating rates ",
	     (ftnlen)57);
    do_fio(&c__1, "      at each computational level. ", (ftnlen)35);
    do_fio(&c__1, "   5) same as (1) with radiances saved at arbitrary ", (
	    ftnlen)52);
    do_fio(&c__1, "      azimuths and zenith angles. ", (ftnlen)34);
    do_fio(&c__1, "   6) same as (2) with radiances saved at arbitrary ", (
	    ftnlen)52);
    do_fio(&c__1, "      azimuths and zenith angles. ", (ftnlen)34);
    do_fio(&c__1, "   7) same as (3) with radiances saved at arbitrary ", (
	    ftnlen)52);
    do_fio(&c__1, "      azimuths and zenith angles. ", (ftnlen)34);
    do_fio(&c__1, "   8) same as (4) with radiances saved at arbitrary ", (
	    ftnlen)52);
    do_fio(&c__1, "      azimuths and zenith angles.", (ftnlen)33);
    do_fio(&c__1, " Enter index of choice: ", (ftnlen)24);
    e_wsfe();
    i__1 = s_rsle(&io___292);
    if (i__1 != 0) {
	goto L4621;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*irad), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L4621;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L4621;
    }
    s_wsfe(&io___293);
    do_fio(&c__1, "irad =", (ftnlen)6);
    do_fio(&c__1, (char *)&(*irad), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*irad <= 4) {
	*usrang = FALSE_;
	*numu = *nstr;
    } else {

/* ****      set disort flag that creates output at user-specified angles */

	*usrang = TRUE_;

	ic = 0;
L4641:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___294);
	do_fio(&c__1, " Enter number of arbitrary emission zenith angles (> "
		"0): ", (ftnlen)57);
	do_fio(&c__1, " [half negative (downward) and half positive (upward)]"
		, (ftnlen)54);
	e_wsfe();
	i__1 = s_rsle(&io___295);
	if (i__1 != 0) {
	    goto L4641;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*numu), (ftnlen)sizeof(integer))
		;
	if (i__1 != 0) {
	    goto L4641;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L4641;
	}
	s_wsfe(&io___296);
	do_fio(&c__1, "numu =", (ftnlen)6);
	do_fio(&c__1, (char *)&(*numu), (ftnlen)sizeof(integer));
	e_wsfe();

	if (*numu > 16) {
	    s_wsfe(&io___297);
	    do_fio(&c__1, " number of streams exceeds dimension", (ftnlen)36);
	    do_fio(&c__1, " bound: mxumu=", (ftnlen)14);
	    do_fio(&c__1, (char *)&c__16, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_stop("", (ftnlen)0);
	}

	s_wsfe(&io___298);
	do_fio(&c__1, " Enter output emission zenith angles (180 - 0 degrees"
		"): ", (ftnlen)56);
	do_fio(&c__1, " (Note: Angles must be monotonically DECREASING) ", (
		ftnlen)49);
	e_wsfe();
	i__1 = *numu;
	for (n = 1; n <= i__1; ++n) {
	    ic = 0;
L4661:
	    ++ic;
	    if (ic > icmx) {
		s_stop("", (ftnlen)0);
	    }
	    s_wsfe(&io___299);
	    do_fio(&c__1, " Enter emission zenith angle #", (ftnlen)30);
	    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	    e_wsfe();
	    i__2 = s_rsle(&io___300);
	    if (i__2 != 0) {
		goto L4661;
	    }
	    i__2 = do_lio(&c__4, &c__1, (char *)&ang, (ftnlen)sizeof(real));
	    if (i__2 != 0) {
		goto L4661;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L4661;
	    }
	    umu[n] = cos(pi * ang / 180.f);
	    s_wsfe(&io___302);
	    do_fio(&c__1, "ang =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ang, (ftnlen)sizeof(real));
	    e_wsfe();
	    s_wsfe(&io___303);
	    do_fio(&c__1, "umu = ", (ftnlen)6);
	    do_fio(&c__1, (char *)&umu[n], (ftnlen)sizeof(real));
	    e_wsfe();
/* L4681: */
	}
    }

/* ****    specify emission azimuth angles */

    ic = 0;
L4802:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___304);
    do_fio(&c__1, " Enter number of output emission azimuth angles (> 0): ", (
	    ftnlen)55);
    e_wsfe();
    i__1 = s_rsle(&io___305);
    if (i__1 != 0) {
	goto L4802;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*nphi), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L4802;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L4802;
    }
    s_wsfe(&io___306);
    do_fio(&c__1, "nphi =", (ftnlen)6);
    do_fio(&c__1, (char *)&(*nphi), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*nphi > 16) {
	s_wsfe(&io___307);
	do_fio(&c__1, " number of output azimuths exceeds dimension", (ftnlen)
		44);
	do_fio(&c__1, " bounds: mxphi=", (ftnlen)15);
	do_fio(&c__1, (char *)&c__16, (ftnlen)sizeof(integer));
	e_wsfe();
	s_stop("", (ftnlen)0);
    }

    s_wsfe(&io___308);
    do_fio(&c__1, " Enter output emission azimuth angles (degrees): ", (
	    ftnlen)49);
    e_wsfe();
    i__1 = *nphi;
    for (n = 1; n <= i__1; ++n) {
	ic = 0;
L4821:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___309);
	do_fio(&c__1, " Enter emission azimuth angle #", (ftnlen)31);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	i__2 = s_rsle(&io___310);
	if (i__2 != 0) {
	    goto L4821;
	}
	i__2 = do_lio(&c__4, &c__1, (char *)&phi[n], (ftnlen)sizeof(real));
	if (i__2 != 0) {
	    goto L4821;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L4821;
	}
	s_wsfe(&io___311);
	do_fio(&c__1, "phi =", (ftnlen)5);
	do_fio(&c__1, (char *)&phi[n], (ftnlen)sizeof(real));
	e_wsfe();
/* L4841: */
    }

/* ****    o u t p u t    l e v e l s */

    ic = 0;
L5002:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___312);
    do_fio(&c__1, " Specify level at which output radiances are to be printe"
	    "d: ", (ftnlen)60);
    do_fio(&c__1, "   1) top of the atmosphere only, ", (ftnlen)34);
    do_fio(&c__1, "   2) surface only, ", (ftnlen)20);
    do_fio(&c__1, "   3) top of atmosphere and surface,", (ftnlen)36);
    do_fio(&c__1, "   4) one or more (monotonically-increasing) pressure lev"
	    "els", (ftnlen)60);
    do_fio(&c__1, " Choose index of output level: ", (ftnlen)31);
    e_wsfe();
    i__1 = s_rsle(&io___313);
    if (i__1 != 0) {
	goto L5002;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*levout), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L5002;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L5002;
    }
    s_wsfe(&io___314);
    do_fio(&c__1, "levout =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*levout), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*levout <= 2) {
	*nlout = 1;
	if (*levout == 1) {
	    s_copy(clev, "_toa", (ftnlen)4, (ftnlen)4);
	} else {
	    s_copy(clev, "_sur", (ftnlen)4, (ftnlen)4);
	}
    } else {
	if (*levout == 3) {
	    *nlout = 2;
	    s_copy(clev, "_toa", (ftnlen)4, (ftnlen)4);
	    s_copy(clev + 4, "_sur", (ftnlen)4, (ftnlen)4);
	} else {
	    s_wsfe(&io___316);
	    do_fio(&c__1, " Enter the number or arbitrary output levels: ", (
		    ftnlen)46);
	    do_fio(&c__1, " NOTE: the current maximum is mxlout = ", (ftnlen)
		    39);
	    do_fio(&c__1, (char *)&c__3, (ftnlen)sizeof(integer));
	    e_wsfe();
	    i__1 = s_rsle(&io___317);
	    if (i__1 != 0) {
		goto L5002;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*nlout), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L5002;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L5002;
	    }
	    s_wsle(&io___318);
	    do_lio(&c__9, &c__1, "nlout = ", (ftnlen)8);
	    do_lio(&c__3, &c__1, (char *)&(*nlout), (ftnlen)sizeof(integer));
	    e_wsle();
	    if (*nlout > 3) {
		s_wsle(&io___319);
		do_lio(&c__9, &c__1, "nlout exceeds dimension bound, mxlout=",
			 (ftnlen)38);
		do_lio(&c__3, &c__1, (char *)&c__3, (ftnlen)sizeof(integer));
		e_wsle();
		s_stop("", (ftnlen)0);
	    }

	    i__1 = *nlout;
	    for (k = 1; k <= i__1; ++k) {
		s_wsfe(&io___321);
		do_fio(&c__1, " Enter the output pressure (bars) for level #"
			": ", (ftnlen)47);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		e_wsfe();
		i__2 = s_rsle(&io___322);
		if (i__2 != 0) {
		    goto L5002;
		}
		i__2 = do_lio(&c__4, &c__1, (char *)&pout[k], (ftnlen)sizeof(
			real));
		if (i__2 != 0) {
		    goto L5002;
		}
		i__2 = e_rsle();
		if (i__2 != 0) {
		    goto L5002;
		}
		s_wsfe(&io___323);
		do_fio(&c__1, "pout =", (ftnlen)6);
		do_fio(&c__1, (char *)&pout[k], (ftnlen)sizeof(real));
		e_wsfe();
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 4;
		ici__1.iciunit = clev + (k - 1 << 2);
		ici__1.icifmt = "(1a2,i2.2)";
		s_wsfi(&ici__1);
		do_fio(&c__1, "_l", (ftnlen)2);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		e_wsfi();
/* L5011: */
	    }
	}
    }

    ic = 0;
L5021:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___324);
    do_fio(&c__1, " Enter index of desired output radiance units: ", (ftnlen)
	    47);
    do_fio(&c__1, "  1) Watts/m**2/sr/cm**-1        ", (ftnlen)33);
    do_fio(&c__1, "  2) Watts/m**2/sr/micron        ", (ftnlen)33);
    do_fio(&c__1, "  3) Watts/m**2/sr/nanometer     ", (ftnlen)33);
    do_fio(&c__1, "  4) ergs/s/cm**2/sr/cm-1        ", (ftnlen)33);
    do_fio(&c__1, "  5) photons/s/m**2/sr/cm-1      ", (ftnlen)33);
    e_wsfe();
    i__1 = s_rsle(&io___325);
    if (i__1 != 0) {
	goto L5021;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*iunits), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L5021;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L5021;
    }
    s_wsfe(&io___326);
    do_fio(&c__1, "iunits =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*iunits), (ftnlen)sizeof(integer));
    e_wsfe();

/* ****    enter characteristics of the output spectral grid. */

    ic = 0;
L6001:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___327);
    do_fio(&c__1, " Enter min and max wavenumbers of the ", (ftnlen)38);
    do_fio(&c__1, "output spectral grid (cm**-1): ", (ftnlen)31);
    e_wsfe();
    i__1 = s_rsle(&io___328);
    if (i__1 != 0) {
	goto L6001;
    }
    i__1 = do_lio(&c__5, &c__1, (char *)&(*wnmin), (ftnlen)sizeof(doublereal))
	    ;
    if (i__1 != 0) {
	goto L6001;
    }
    i__1 = do_lio(&c__5, &c__1, (char *)&(*wnmax), (ftnlen)sizeof(doublereal))
	    ;
    if (i__1 != 0) {
	goto L6001;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L6001;
    }
    s_wsfe(&io___329);
    do_fio(&c__1, "wnmin, wnmax =", (ftnlen)14);
    do_fio(&c__1, (char *)&(*wnmin), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*wnmax), (ftnlen)sizeof(doublereal));
    e_wsfe();

    ic = 0;
L6201:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___330);
    do_fio(&c__1, " Choose index of output spectrum type: ", (ftnlen)39);
    do_fio(&c__1, "  1) Full-resolution spectrum (this creates enourmous fil"
	    "es)", (ftnlen)60);
    do_fio(&c__1, "  2) Sampled spectrum convolved with a slit function ", (
	    ftnlen)53);
    e_wsfe();
/*     - '  3) SMT binned data with spectral map: ' */
    i__1 = s_rsle(&io___331);
    if (i__1 != 0) {
	goto L6201;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*isptype), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L6201;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L6201;
    }
    s_wsfe(&io___332);
    do_fio(&c__1, "isptype =", (ftnlen)9);
    do_fio(&c__1, (char *)&(*isptype), (ftnlen)sizeof(integer));
    e_wsfe();

    if (*isptype == 2) {

	ic = 0;
L6221:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___333);
	do_fio(&c__1, " Enter type of spectral response function: ", (ftnlen)
		43);
	do_fio(&c__1, "  1) boxcar", (ftnlen)11);
	do_fio(&c__1, "  2) triangular (approximate slit spectrometer)", (
		ftnlen)47);
	e_wsfe();

	i__1 = s_rsle(&io___334);
	if (i__1 != 0) {
	    goto L6221;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*islit), (ftnlen)sizeof(integer)
		);
	if (i__1 != 0) {
	    goto L6221;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L6221;
	}
	s_wsfe(&io___335);
	do_fio(&c__1, "islit =", (ftnlen)7);
	do_fio(&c__1, (char *)&(*islit), (ftnlen)sizeof(integer));
	e_wsfe();

	ic = 0;
L6241:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___336);
	do_fio(&c__1, " Enter the half-width-at-half-max: ", (ftnlen)35);
	e_wsfe();
	i__1 = s_rsle(&io___337);
	if (i__1 != 0) {
	    goto L6241;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&(*width), (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L6241;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L6241;
	}
	s_wsfe(&io___338);
	do_fio(&c__1, "width =", (ftnlen)7);
	do_fio(&c__1, (char *)&(*width), (ftnlen)sizeof(doublereal));
	e_wsfe();

	ic = 0;
L6281:
	++ic;
	if (ic > icmx) {
	    s_stop("", (ftnlen)0);
	}
	s_wsfe(&io___339);
	do_fio(&c__1, " Enter sampling resolution of output file (cm**-1):", (
		ftnlen)51);
	e_wsfe();
	i__1 = s_rsle(&io___340);
	if (i__1 != 0) {
	    goto L6281;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&(*dwn), (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L6281;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L6281;
	}
	s_wsfe(&io___341);
	do_fio(&c__1, "dwn =", (ftnlen)5);
	do_fio(&c__1, (char *)&(*dwn), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* ****   set error limits for spectral mapping method */

    ic = 0;
L6801:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___342);
    do_fio(&c__1, " Enter fractional error for tau binning: ", (ftnlen)41);
    e_wsfe();
    i__1 = s_rsle(&io___343);
    if (i__1 != 0) {
	goto L6801;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*tauerr), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L6801;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L6801;
    }
    s_wsfe(&io___344);
    do_fio(&c__1, "tauerr =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*tauerr), (ftnlen)sizeof(real));
    e_wsfe();
    ic = 0;
L6821:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___345);
    do_fio(&c__1, " Enter error for pi0 binning: ", (ftnlen)30);
    e_wsfe();
    i__1 = s_rsle(&io___346);
    if (i__1 != 0) {
	goto L6821;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*pi0err), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L6821;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L6821;
    }
    s_wsfe(&io___347);
    do_fio(&c__1, "pi0err =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*pi0err), (ftnlen)sizeof(real));
    e_wsfe();
    ic = 0;
L6841:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___348);
    do_fio(&c__1, " Enter error for <cos> binning: ", (ftnlen)32);
    e_wsfe();
    i__1 = s_rsle(&io___349);
    if (i__1 != 0) {
	goto L6841;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*phferr), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L6841;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L6841;
    }
    s_wsfe(&io___350);
    do_fio(&c__1, "phferr =", (ftnlen)8);
    do_fio(&c__1, (char *)&(*phferr), (ftnlen)sizeof(real));
    e_wsfe();
    ic = 0;
L6861:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___351);
    do_fio(&c__1, " Enter error for surface optics binning: ", (ftnlen)41);
    e_wsfe();
    i__1 = s_rsle(&io___352);
    if (i__1 != 0) {
	goto L6861;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*surferr), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L6861;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L6861;
    }
    s_wsfe(&io___353);
    do_fio(&c__1, "surferr =", (ftnlen)9);
    do_fio(&c__1, (char *)&(*surferr), (ftnlen)sizeof(real));
    e_wsfe();

/* ****   n a m e    o f    o u t p u t    f i l e */

    ic = 0;
L7001:
    ++ic;
    if (ic > icmx) {
	s_stop("", (ftnlen)0);
    }
    s_wsfe(&io___354);
    do_fio(&c__1, " Enter index of output file type: ", (ftnlen)34);
    do_fio(&c__1, "  1) ASCII,   2) Binary with header,   3) binary no header"
	    , (ftnlen)58);
    e_wsfe();
    i__1 = s_rsle(&io___355);
    if (i__1 != 0) {
	goto L7001;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ifrmout), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L7001;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L7001;
    }
    s_wsfe(&io___356);
    do_fio(&c__1, " ifrmout = ", (ftnlen)11);
    do_fio(&c__1, (char *)&(*ifrmout), (ftnlen)sizeof(integer));
    e_wsfe();

    init_spect_io__(nza, nlout, irad, ifrmout, sza0, clev, iuout, iuflx, 
	    iuheat, iustat, iutrn, name__, &len, radfile, heatfile, statfile, 
	    trnfile, flxfile, (ftnlen)4, (ftnlen)1, (ftnlen)132, (ftnlen)132, 
	    (ftnlen)132, (ftnlen)132, (ftnlen)132);

/* ****    define the first output zenith angle for radiances */
/*        (only upward radiances are saved at the top fo the atmosphere) */

    if (*usrang) {
	nzdn = 0;
	i__1 = *numu;
	for (nze = 1; nze <= i__1; ++nze) {
	    if (umu[nze] < 0.f) {
		nzdn = nze;
	    }
/* L8021: */
	}
	nzup = nzdn + 1;
    } else {
	nzdn = *nstr / 2;
	nzup = nzdn + 1;
    }

    i__1 = *nlout;
    for (nl = 1; nl <= i__1; ++nl) {
	if (*levout == 1 || *levout == 3 && nl == 1) {
	    nza_1__[nl - 1] = nzup;
	} else {
	    nza_1__[nl - 1] = 1;
	}
/* L8011: */
    }

/* ****       create an output files for each partial derivative */

    init_pd_iu__(nza, numu, nphi, nlout, ifrmout, nza_1__, radfile, lplanck, 
	    lsolar, clev, nstate, &istate[1], iu_pd__, j_ext__, &iutpd[11], &
	    iuspd[11], &iupdrad[10513], (ftnlen)132, (ftnlen)4, (ftnlen)9);

    return 0;
} /* smart_in__ */

