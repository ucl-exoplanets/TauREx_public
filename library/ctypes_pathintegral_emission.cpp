/* 

    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Forward model for emission

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

    Compile with:   g++ -fPIC -shared -o ctypes_pathintegral_emission.so ctypes_pathintegral_emission.cpp

 */

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

extern "C" {

    void path_integral(const double * wngrid,
                       const int nwngrid,
                       const int nlayers,
                       const int nactive,
                       const int cia,
                       const double * sigma_array,
                       const double * sigma_temp,
                       const int sigma_ntemp,
                       const int cia_npairs,
                       const double * cia_idx,
                       const int cia_nidx,
                       const double * sigma_cia,
                       const double * sigma_cia_temp,
                       const int sigma_cia_ntemp,
                       const double * density,
                       const double * z,
                       const double * active_mixratio,
                       const double * inactive_mixratio,
                       const double * temperature,
                       const double planet_radius,
                       const double star_radius,
                       const double * star_sed,
                       void * FpFsv,
                       void * tauv) {

        double * FpFs = (double *) FpFsv;
        double * tau_total = (double *) tauv;

        // setting up arrays and variables
        double* dz = new double[nlayers];
        double* sigma_interp = new double[nwngrid*nlayers*nactive];
        double* sigma_cia_interp = new double[nwngrid * nlayers * cia_npairs];
        double sigma, sigma_l, sigma_r;
        double tau, dtau;
        double p;
        int count, count2, t_idx;
        double I_total, BB_wl, exponent;
        double h, c, kb, pi;
        double x1_idx[cia_npairs][nlayers];
        double x2_idx[cia_npairs][nlayers];

        h = 6.62606957e-34;
        c = 299792458;
        kb = 1.3806488e-23;
        pi= 3.14159265359;

        //dz array
//        for (int j=0; j<(nlayers); j++) {
//            if ((j+1) == nlayers) {
//                dz[j] = z[j] - z[j-1];
//            } else {
//                dz[j] = z[j+1] - z[j];
//            }
//        }
        for (int j=0; j<(nlayers-1); j++) {
                dz[j] = z[j+1] - z[j];
        }

        // interpolate sigma array to the temperature profile
        for (int j=0; j<nlayers; j++) {
            if (temperature[j] > sigma_temp[sigma_ntemp-1]) {
                for (int wn=0; wn<nwngrid; wn++) {
                    for (int l=0;l<nactive;l++) {
                        sigma_interp[wn + nwngrid*(j + l*nlayers)] = sigma_array[wn + nwngrid*(sigma_ntemp-1 + sigma_ntemp*(j + l*nlayers))];
                    }
                }
            } else if (temperature[j] < sigma_temp[0]) {
                for (int wn=0; wn<nwngrid; wn++) {
                    for (int l=0;l<nactive;l++) {
                        sigma_interp[wn + nwngrid*(j + l*nlayers)] = sigma_array[wn + nwngrid*(sigma_ntemp*(j + l*nlayers))];
                    }
                }
            } else {
                if (sigma_ntemp == 1) { // This only happens for create_spectrum (when temperature is part of sigma_t)
                    for (int wn=0; wn<nwngrid; wn++) {
                        for (int l=0;l<nactive;l++) {
                            sigma_interp[wn + nwngrid*(j + l*nlayers)] = sigma_array[wn + nwngrid*(sigma_ntemp*(j + l*nlayers))];
                        }
                    }
                } else {
                    for (int t=1; t<sigma_ntemp; t++) {
                        if ((temperature[j] >= sigma_temp[t-1]) && (temperature[j] < sigma_temp[t])) {
                            for (int wn=0; wn<nwngrid; wn++) {
                                for (int l=0;l<nactive;l++) {
                                    sigma_l = sigma_array[wn + nwngrid*(t-1 + sigma_ntemp*(j + l*nlayers))];
                                    sigma_r = sigma_array[wn + nwngrid*(t + sigma_ntemp*(j + l*nlayers))];
                                    sigma = sigma_l + (sigma_r-sigma_l)*(temperature[j]-sigma_temp[t-1])/(sigma_temp[t]-sigma_temp[t-1]);
                                    sigma_interp[wn + nwngrid*(j + l*nlayers)] = sigma;
                                }
                            }
                        }
                    }
                }
            }
        }


        // interpolate sigma CIA array to the temperature profile
        for (int j=0; j<nlayers; j++) {
            for (int t=1; t<sigma_cia_ntemp; t++) {
                if ((temperature[j] > sigma_cia_temp[t-1]) && (temperature[j] < sigma_cia_temp[t])) {
                    for (int wn=0; wn<nwngrid; wn++) {
                        for (int l=0;l<cia_npairs;l++) {
                            sigma_l = sigma_cia[wn + nwngrid*(t-1 + sigma_cia_ntemp*l)];
                            sigma_r = sigma_cia[wn + nwngrid*(t + sigma_cia_ntemp*l)];
                            sigma = sigma_l + (sigma_r-sigma_l)*(log10(temperature[j])-log10(sigma_temp[t-1]))/(log10(sigma_temp[t])-log10(sigma_temp[t-1]));
                            sigma_cia_interp[wn +  nwngrid*(j + l*nlayers)] = sigma;
                        }
                    }
                }
            }
        }

        // get mixing ratio of individual molecules in the collision induced absorption (CIA) pairs
        for (int c=0; c<cia_npairs;c++) {
             if (int(cia_idx[c*2]) >= nactive) {
                for (int j=0;j<nlayers;j++) {

                    x1_idx[c][j] = inactive_mixratio[j+nlayers*(int(cia_idx[c*2])-nactive)];
                    x2_idx[c][j] = inactive_mixratio[j+nlayers*(int(cia_idx[c*2+1])-nactive)];
                }
             } else {
                for (int j=0;j<nlayers;j++) {
                    x1_idx[c][j] = active_mixratio[j+nlayers*int(cia_idx[c*2])];
                    x2_idx[c][j] = active_mixratio[j+nlayers*int(cia_idx[c*2+1])];
                }
            }
         }


        // calculate emission
        count2 = 0;
        for (int wn=0; wn < nwngrid; wn++) {

            tau = 0.0;
            I_total = 0.0;

            // surface layer
            for(int j = 0; j < nlayers; j++){
    			for (int l=0;l<nactive;l++) { // active gases
                    tau += (sigma_interp[wn + nwngrid*(j + l*nlayers)] * active_mixratio[j+nlayers*l] * density[j] * dz[j]);
                }
                if (cia == 1) { // cia
                    for (int c=0; c<cia_npairs;c++) {
                        tau += sigma_cia[wn + nwngrid*c] * x1_idx[c][j]*x2_idx[c][j] * density[j]*density[j] * dz[j];
                    }
                }
            }
            exponent = exp((h * c) / ((10000./wngrid[wn])*1e-6  * kb * temperature[0]));
            BB_wl = ((pi*2.0*h*pow(c,2))/pow((10000./wngrid[wn])*1e-6,5) * (1.0/(exponent - 1)))* 1e-6; // (W/m^2/micron)

            I_total += BB_wl * (exp(-1.0*tau));

            //other layers
    		for (int j=1; j<(nlayers); j++) {
    			tau = 0.0;
    			dtau = 0.0;
    			for (int k=j; k < (nlayers); k++) {
                    for (int l=0;l<nactive;l++) { // active gases
                        tau += (sigma_interp[wn + nwngrid*(k + l*nlayers)] * active_mixratio[k+nlayers*l] * density[k] * dz[k]);
                    }
                    if (cia == 1) { // cia
                        for (int c=0; c<cia_npairs;c++) {
                            tau += sigma_cia[wn + nwngrid*c] * x1_idx[c][k]*x2_idx[c][k] * density[k]*density[k] * dz[k];
                        }
                    }
    			}
    			// get dtau
                for (int l=0;l<nactive;l++) { // active gases
                    dtau += (sigma_interp[wn + nwngrid*(j + l*nlayers)] * active_mixratio[j+nlayers*l] * density[j] * dz[j]);
                }
                if (cia == 1) { // cia
                    for (int c=0; c<cia_npairs;c++) {
                        dtau += sigma_cia[wn + nwngrid*c] * x1_idx[c][j]*x2_idx[c][j] * density[j]*density[j] * dz[j];
                    }
                }

                // get blackbody
                if (temperature[j] != temperature[j-1]) {
                    exponent = exp((h * c) / ((10000./wngrid[wn])*1e-6  * kb * temperature[j]));
                    BB_wl = ((pi*2.0*h*pow(c,2))/pow((10000./wngrid[wn])*1e-6,5) * (1.0/(exponent - 1)))* 1e-6; // (W/m^2/micron)
                }
                tau_total[count2] = (BB_wl * (exp(-1.0*tau))* dtau);
                I_total += (BB_wl * (exp(-1.0*tau))* dtau);
                count2 += 1;
            }
//    		FpFs[wn] = I_total;
            FpFs[wn] = (I_total/star_sed[wn]) * pow((planet_radius/star_radius), 2);
        }

        delete dz;
        delete sigma_interp;

    }
}
