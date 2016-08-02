/* 

    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Forward model for transmission

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

    Compile with:   g++ -fPIC -shared -o ctypes_pathintegral_transmission_xsec.so ctypes_pathintegral_transmission_xsec.cpp

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

    void path_integral(const int nwngrid,
                       const int nlayers,
                       const int nactive,
                       const int ninactive,
                       const int rayleigh,
                       const int cia,
                       const int clouds,
                       const double * sigma_array,
                       const double * sigma_temp,
                       const int sigma_ntemp,
                       const double * sigma_rayleigh,
                       const int cia_npairs,
                       const double * cia_idx,
                       const int cia_nidx,
                       const double * sigma_cia,
                       const double * sigma_cia_temp,
                       const int sigma_cia_ntemp,
                       const double cloud_topP,
                       const double * pressure,
                       const double * density,
                       const double * z,
                       const double * z_levels,
                       const double * active_mixratio,
                       const double * inactive_mixratio,
                       const double * temperature,
                       const double planet_radius,
                       const double star_radius,
                       void * absorptionv,
                       void * tauv) {


        double * absorption = (double *) absorptionv;
        double * tau = (double *) tauv;

        // setting up arrays and variables
        //double* tau = new double[nlayers*nwngrid];
        double* dz = new double[nlayers];
        double* dlarray = new double[nlayers*nlayers];
        double* sigma_interp = new double[nwngrid*nlayers*nactive];
        double* sigma_cia_interp = new double[nwngrid * nlayers * cia_npairs];
        double sigma, sigma_l, sigma_r;
        double x1_idx[cia_npairs][nlayers];
        double x2_idx[cia_npairs][nlayers];
        double tautmp, exptau,  integral;
        double cld_factor, cld_rho;
        double p;
        int count, count2, t_idx;

        // dl array
        count = 0;
        for (int j=0; j<(nlayers); j++) {
            for (int k=0; k < (nlayers - j); k++) {
                if (k == 0) {
                    dlarray[count] = 2.0 * sqrt(pow((z_levels[j+1]+planet_radius),2) - pow((z[j]+planet_radius),2));
                } else {
                    p = pow((z[j]+planet_radius),2);
                    dlarray[count] = 2.0 * (sqrt(pow((z_levels[k+j+1]+planet_radius),2) - p) -  sqrt(pow((z[k+j]+planet_radius),2) - p));
                }
                count += 1;
            }
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
            if (sigma_ntemp == 1) { //
                for (int wn=0; wn<nwngrid; wn++) {
                     for (int l=0;l<cia_npairs;l++) {
                        sigma_cia_interp[wn +  nwngrid*(j + l*nlayers)] = sigma_cia[wn + nwngrid*(sigma_cia_ntemp*l)];
                     }
                }
            } else {
                if (temperature[j] > sigma_cia_temp[sigma_cia_ntemp-1]) {
                    for (int wn=0; wn<nwngrid; wn++) {
                         for (int l=0;l<cia_npairs;l++) {
                            sigma_cia_interp[wn +  nwngrid*(j + l*nlayers)] = sigma_cia[wn + nwngrid*(sigma_cia_ntemp-1 + sigma_cia_ntemp*l)];
                         }
                    }
                } else if  (temperature[j] <  sigma_cia_temp[0]) {
                    for (int wn=0; wn<nwngrid; wn++) {
                         for (int l=0;l<cia_npairs;l++) {
                            sigma_cia_interp[wn +  nwngrid*(j + l*nlayers)] = sigma_cia[wn + nwngrid*(sigma_cia_ntemp*l)];
                         }
                    }
                } else {
                    for (int t=1; t<sigma_cia_ntemp; t++) {
                        if ((temperature[j] > sigma_cia_temp[t-1]) && (temperature[j] < sigma_cia_temp[t])) {
                            for (int wn=0; wn<nwngrid; wn++) {
                                for (int l=0;l<cia_npairs;l++) {
                                    sigma_l = sigma_cia[wn + nwngrid*(t-1 + sigma_cia_ntemp*l)];
                                    sigma_r = sigma_cia[wn + nwngrid*(t + sigma_cia_ntemp*l)];
                                    sigma = sigma_l + (sigma_r-sigma_l)*(temperature[j]-sigma_cia_temp[t-1])/(sigma_cia_temp[t]-sigma_cia_temp[t-1]);
                                    sigma_cia_interp[wn +  nwngrid*(j + l*nlayers)] = sigma;
                                }
                            }
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

        // calculate absorption
        count2 = 0;
        for (int wn=0; wn < nwngrid; wn++) {
            count = 0;
            integral = 0.0;
            //cout << " integral 0 " << integral << endl;
            //cout << " count 0 " << count << endl;
    		for (int j=0; j<(nlayers); j++) { 	// loop through atmosphere layers, z[0] to z[nlayers]
    			tautmp = 0.0;
                if ((clouds == 1) && (pressure[j] >= cloud_topP)) {
                    //cout << j << " YES " << pressure[j] << endl;
                    for (int k=0; k < (nlayers-j); k++) { // loop through each layer to sum up path length
                        count += 1;
                    }
                    //exptau = exp(-tautmp);
                    integral += ((planet_radius+z[j])*(1.0)*dz[j]);
                    tau[count2] = 1.0;
                    //cout << count2 << " " << tau[count2] << endl;
                    //cout << count << " j " << j << " z " << z[j] << " dz " << dz[j] << " exptau  " << exptau << " integral " << integral << endl;
                    count2 += 1;
                } else {
                    //cout << j << " NO " << pressure[j] << endl;

                    for (int k=0; k < (nlayers-j); k++) { // loop through each layer to sum up path length


                        // calculate optical depth due to clouds
                        // calculate optical depths due to active absorbing gases (absorption + rayleigh scattering)
                        for (int l=0;l<nactive;l++) {
                            sigma = sigma_interp[wn + nwngrid*((k+j) + l*nlayers)];
                            tautmp += (sigma * active_mixratio[k+j+nlayers*l] * density[k+j] * dlarray[count]);
                            //cout << " j " << j  << " k " << k  << " count " << count << " sigma " << sigma << " active_mixratio " << active_mixratio[k+j+nlayers*l] << " density " << density[k+j] << " dlarray " << dlarray[count] << " tau " << (sigma * active_mixratio[k+j+nlayers*l] * density[k+j] * dlarray[count]) << endl;
                            //cout << " j " << j  << " k " << k  << " count " << count << " sigma_rayleigh " << sigma_rayleigh[wn + nwngrid*l] << " active_mixratio " << active_mixratio[k+j+nlayers*l] << " density " << density[k+j] << " dlarray " << dlarray[count] << endl;
                            if (rayleigh == 1) {
                                tautmp += sigma_rayleigh[wn + nwngrid*l] * active_mixratio[k+j+nlayers*l] * density[j+k] * dlarray[count];
                            }
                        }

                        // calculating optical depth due inactive gases (rayleigh scattering)
                        if (rayleigh == 1) {
                            for (int l=0; l<ninactive; l++) {
                                //cout << sigma_rayleigh[wn + nwngrid*(l+nactive)] << " " << inactive_mixratio[k+j+nlayers*l] << " " << density[j+k] << " " << dlarray[count] << endl;
                                //tautmp += sigma_rayleigh[wn + nwngrid*(l+nactive)] * inactive_mixratio[k+j+nlayers*l] * density[j+k] * dlarray[count];
                                tautmp += sigma_rayleigh[wn + nwngrid*(l+nactive)] * inactive_mixratio[k+j+nlayers*l] * density[j+k] * dlarray[count];
                            }
                        }
                        // calculating optical depth due to collision induced absorption
                        if (cia == 1) {
                            for (int c=0; c<cia_npairs;c++) {
                                tautmp += sigma_cia[wn + nwngrid*c] * x1_idx[c][k+j]*x2_idx[c][k+j] * density[j+k]*density[j+k] * dlarray[count];
                            }
                        }

                        count += 1;
                    }
                    exptau = exp(-tautmp);
                    integral += ((planet_radius+z[j])*(1.0-exptau)*dz[j]);
                    tau[count2] = 1.0 - exptau;
                    //cout << count2 << " " << tau[count2] << endl;
                    //cout << count << " j " << j << " z " << z[j] << " dz " << dz[j] << " exptau  " << exptau << " integral " << integral << endl;
                    count2 += 1;
                }
            }
            integral *= 2.0;
            //cout << integral << endl;
            absorption[wn] = ((planet_radius*planet_radius) + integral) / (star_radius*star_radius);
            //cout << wn << " " << absorption[wn] << endl;
        }
//        cout << "END" << endl;

        delete[] dlarray;
        delete[] sigma_interp;
        delete[] sigma_cia_interp;
        dlarray = NULL;
        sigma_interp = NULL;
        sigma_cia_interp = NULL;
    }
}
