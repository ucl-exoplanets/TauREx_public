/* xyinterp.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int xyinterp_(real *xd, real *zd, real *xi, real *zi, 
	integer *nxmax, integer *nymax, integer *nxd, integer *nxi, integer *
	ny)
{
    /* System generated locals */
    integer zd_dim1, zd_offset, zi_dim1, zi_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, ll;
    static real xp1, xp2;
    static integer iord, iori, lstd;


/* ccccccccccccccccccccccccc  x y i n t e r p  ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e:                                                  cc */
/* c                                                                    cc */
/* c    this subroutine uses linear interpolation for an unequally-     cc */
/* c    spaced input grid to interpolate a 2-d file to a new grid.      cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c    xd    : vector of length nxd. x points where input array is     cc */
/* c            known.                                                  cc */
/* c    zd    : nxd array of values defined at points xd.               cc */
/* c    nxd   : number of points in input x-array                       cc */
/* c    xi    : vector of length nxi. x points where array is needed.   cc */
/* c    zi    : approximate value of array zd at points xi, yi.         cc */
/* c    nxi   : number of points in x-vector xi                         cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    zi    : value of input function, zd, at points xi.              cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  x y i n t e r p  ccccccccccccccccccccccccccc */




    /* Parameter adjustments */
    zi_dim1 = *nxmax;
    zi_offset = 1 + zi_dim1;
    zi -= zi_offset;
    zd_dim1 = *nxmax;
    zd_offset = 1 + zd_dim1;
    zd -= zd_offset;
    --xd;
    --xi;

    /* Function Body */
    if (*nxd == 1) {
	i__1 = *ny;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *nxi;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		zi[i__ + j * zi_dim1] = zd[j * zd_dim1 + 1];
/* L1001: */
	    }
/* L1021: */
	}
	return 0;
    }

/* ****   determine the ordering of the array xi */

    iori = 0;
    if (*nxi > 1) {
	if (xi[*nxi] - xi[1] < 0.f) {
	    iori = 1;
	}
    }

/* ****   determine the ordering of the array xd */

    iord = 0;
    if (xd[*nxd] - xd[1] < 0.f) {
	iord = 1;
    }

/*   use linear interpolation to estimate the value of */
/*   zi at the desired x values, xi. */

    if (iord == 0) {
	lstd = 2;
	xp1 = xd[1];
	xp2 = xd[2];

/* ****      enter loop over interpolated positions. */

	i__1 = *nxi;
	for (ll = 1; ll <= i__1; ++ll) {
	    l = ll;
	    if (iori == 1) {
		l = *nxi - ll + 1;
	    }

/* ****          determine if xi is within the bounds of xd */

	    if (xi[l] >= xd[1] && xi[l] <= xd[*nxd]) {

/* ****              determine if xi is between points xp1 and xp2 */

L2001:
		if (xi[l] >= xp1 && xi[l] <= xp2) {
		    i__2 = *ny;
		    for (j = 1; j <= i__2; ++j) {
			zi[l + j * zi_dim1] = zd[lstd - 1 + j * zd_dim1] + (
				zd[lstd + j * zd_dim1] - zd[lstd - 1 + j * 
				zd_dim1]) * (xi[l] - xd[lstd - 1]) / (xd[lstd]
				 - xd[lstd - 1]);
/* L2021: */
		    }
		} else {

/* ****             update counters and xp1, xp2 */

		    ++lstd;
		    xp1 = xp2;
		    xp2 = xd[lstd];
		    goto L2001;
		}
	    } else {

/* ****           xi(l) is outside the range xd(1)-xd(nxd).  extrapolate */

		if (xi[l] < xd[1]) {
		    xp1 = xd[1];
		    xp2 = xd[2];
		    lstd = 2;
		} else {
		    xp1 = xd[*nxd - 1];
		    xp2 = xd[*nxd];
		    lstd = *nxd;
		}

		i__2 = *ny;
		for (j = 1; j <= i__2; ++j) {
		    zi[l + j * zi_dim1] = zd[lstd - 1 + j * zd_dim1] + (zd[
			    lstd + j * zd_dim1] - zd[lstd - 1 + j * zd_dim1]) 
			    * (xi[l] - xd[lstd - 1]) / (xd[lstd] - xd[lstd - 
			    1]);
/* L2041: */
		}
	    }
/* L2062: */
	}
    } else {

/* ****    the input xd array is monotonically decreasing. */

	lstd = 2;
	xp1 = xd[1];
	xp2 = xd[2];

/* ****      enter loop over interpolated positions. */

	i__1 = *nxi;
	for (ll = 1; ll <= i__1; ++ll) {
	    l = ll;
	    if (iori == 0) {
		l = *nxi - ll + 1;
	    }

/* ****          determine if xi is within the bounds of xd */

	    if (xi[l] <= xd[1] && xi[l] >= xd[*nxd]) {

/* ****              determine if xi is between points xp1 and xp2 */

L3001:
		if (xi[l] <= xp1 && xi[l] >= xp2) {
		    i__2 = *ny;
		    for (j = 1; j <= i__2; ++j) {
			zi[l + j * zi_dim1] = zd[lstd - 1 + j * zd_dim1] + (
				zd[lstd + j * zd_dim1] - zd[lstd - 1 + j * 
				zd_dim1]) * (xi[l] - xd[lstd - 1]) / (xd[lstd]
				 - xd[lstd - 1]);
/* L3021: */
		    }
		} else {

/* ****             update counters and xp1, xp2 */

		    ++lstd;
		    xp1 = xp2;
		    xp2 = xd[lstd];
		    goto L3001;
		}
	    } else {

/* ****           xi(l) is outside the range xd(1)-xd(nxd).  extrapolate */

		if (xi[l] > xd[1]) {
		    xp1 = xd[1];
		    xp2 = xd[2];
		    lstd = 2;
		} else {
		    xp1 = xd[*nxd - 1];
		    xp2 = xd[*nxd];
		    lstd = *nxd;
		}

		i__2 = *ny;
		for (j = 1; j <= i__2; ++j) {
		    zi[l + j * zi_dim1] = zd[lstd - 1 + j * zd_dim1] + (zd[
			    lstd + j * zd_dim1] - zd[lstd - 1 + j * zd_dim1]) 
			    * (xi[l] - xd[lstd - 1]) / (xd[lstd] - xd[lstd - 
			    1]);
/* L3041: */
		}
	    }
/* L3062: */
	}
    }

    return 0;
} /* xyinterp_ */

