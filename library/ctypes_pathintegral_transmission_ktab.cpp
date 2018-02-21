/* 

    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Forward model for transmission

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

    Compile with:
    g++ -fPIC -shared -o ctypes_pathintegral_transmission_ktab.so ctypes_pathintegral_transmission_ktab.cpp

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
					   const int mie,
                       const int cia,
                       const int clouds,
                       const double * ktab_array,
                       const double * ktab_temp,
                       const int ktab_ntemp,
                       const int ngauss,
                       const double * ktab_weights,
                       const double * sigma_rayleigh,
                       const int cia_npairs,
                       const double * cia_idx,
                       const int cia_nidx,
                       const double * sigma_cia,
                       const double * sigma_cia_temp,
                       const int sigma_cia_ntemp,
                       const double cloud_topP,
					   const double mie_topP,
					   const double mie_bottomP,
					   const double * sigma_mie,
                       const double * pressure,
                       const double * density,
                       const double * z,
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
        double* ktab_interp = new double[ngauss*nwngrid*nlayers*nactive];
        double* sigma_cia_interp = new double[nwngrid * nlayers * cia_npairs];
        double sigma, sigma_l, sigma_r;
        double x1_idx[cia_npairs][nlayers];
        double x2_idx[cia_npairs][nlayers];
        double tautmp, transtmp, transtot, exptau,  integral;
        double cld_factor, cld_rho;
        double p;
        int count, count_orig, count2, t_idx;

        //dz array
        for (int j=0; j<(nlayers); j++) {
            if ((j+1) == nlayers) {
                dz[j] = z[j] - z[j-1];
            } else {
                dz[j] = z[j+1] - z[j];
            }
        }

        // dl array
        count = 0;
        for (int j=0; j<(nlayers); j++) {
            for (int k=0; k < (nlayers - j); k++) {
                p = pow((planet_radius+dz[0]/2.+z[j]),2);
                if (k == 0) {
                    dlarray[count] = 2.0 * sqrt(pow((planet_radius + dz[0]/2. + z[j] + dz[j]/2.),2) - p);
                } else {
                    dlarray[count] = 2.0 * (sqrt(pow((planet_radius + dz[0]/2. + z[k+j] + dz[j+k]/2.),2) - p) -  sqrt(pow((planet_radius + dz[0]/2. + z[k+j-1] + dz[j+k-1]/2.  ),2) - p));
                }
                count += 1;
            }
        }


        // interpolate ktab array to the temperature profile
        for (int j=0; j<nlayers; j++) {
            if (temperature[j] > ktab_temp[ktab_ntemp-1]) { // temperature higher than ktab
                for (int wn=0; wn<nwngrid; wn++) {
                    for (int l=0;l<nactive;l++) {
                        for (int g=0; g<ngauss; g++) {
                            ktab_interp[g + ngauss*(wn + nwngrid*(j + l*nlayers))] =  ktab_array[g + ngauss*(wn + nwngrid*(ktab_ntemp-1 + ktab_ntemp*(j + l*nlayers)))];
                        }
                    }
                }
            } else if (temperature[j] < ktab_temp[0]) { // temperature lower than ktab
                for (int wn=0; wn<nwngrid; wn++) {
                    for (int l=0;l<nactive;l++) {
                        for (int g=0; g<ngauss; g++) {
                            ktab_interp[g + ngauss*(wn + nwngrid*(j + l*nlayers))] = ktab_array[g + ngauss*(wn + nwngrid*(ktab_ntemp*(j + l*nlayers)))];
                        }
                    }
                }
            } else {
                if (ktab_ntemp == 1) { // This only happens for create_spectrum (when temperature is part of sigma_t)
                    for (int wn=0; wn<nwngrid; wn++) {
                        for (int l=0;l<nactive;l++) {
                            for (int g=0; g<ngauss; g++) {
                                ktab_interp[g + ngauss * (wn + nwngrid*(j + l*nlayers))] = ktab_array[g + ngauss*(wn + nwngrid*(ktab_ntemp*(j + l*nlayers)))];
                            }
                        }
                    }
                } else {
                    for (int t=1; t<ktab_ntemp; t++) {
                        if ((temperature[j] >= ktab_temp[t-1]) && (temperature[j] < ktab_temp[t])) {
                            for (int wn=0; wn<nwngrid; wn++) {
                                for (int l=0;l<nactive;l++) {
                                    for (int g=0; g<ngauss; g++) {
                                        sigma_l = ktab_array[g + ngauss*(wn + nwngrid*(t-1 + ktab_ntemp*(j + l*nlayers)))];
                                        sigma_r = ktab_array[g + ngauss*(wn + nwngrid*(t + ktab_ntemp*(j + l*nlayers)))];
                                        sigma = sigma_l + (sigma_r-sigma_l)*(log10(temperature[j])-log10(ktab_temp[t-1]))/(log10(ktab_temp[t])-log10(ktab_temp[t-1]));
                                        ktab_interp[g + ngauss*(wn + nwngrid*(j + l*nlayers))] = sigma;
                                     }
                                }
                            }
                        }
                    }
                }
            }
        }

        // interpolate sigma CIA array to the temperature profile
        for (int j=0; j<nlayers; j++) {
            if (sigma_cia_ntemp == 1) { //
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
                        if ((temperature[j] >= sigma_cia_temp[t-1]) && (temperature[j] < sigma_cia_temp[t])) {
                            #pragma omp parallel for private(sigma_l, sigma_r, sigma)
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

                count_orig = count;

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

                    transtot = 1.;
                    for (int l=0;l<nactive;l++) {
                        transtmp = 0;
                        for (int g=0; g<ngauss; g++) {
                            count = count_orig;
                            tautmp = 0;
                            for (int k=0; k < (nlayers-j); k++) { // loop through each layer to sum up path length
                                // calculate optical depths due to active absorbing gases (absorption + rayleigh scattering)
                                sigma = ktab_interp[g + ngauss*(wn + nwngrid*((k+j) + l*nlayers))];
                                tautmp += (sigma * active_mixratio[k+j+nlayers*l] * density[k+j] * dlarray[count]);
                                count += 1;
                            }
                            //cout << tautmp << " " << exp(-tautmp) << " " << ktab_weights[g] << endl;
                            transtmp += exp(-tautmp) * ktab_weights[g];
                        }
                        transtot *= transtmp;
                    }

                    // calculate rayleigh and cia separtely from cross sections
                    // firstly get tau
                    count = count_orig;
                    tautmp = 0.;
                    for (int k=0; k < (nlayers-j); k++) { // loop through each layer to sum up path length
                        if (rayleigh == 1) {
                           for (int l=0;l<nactive;l++) {
                                tautmp += sigma_rayleigh[wn + nwngrid*l] * active_mixratio[k+j+nlayers*l] * density[j+k] * dlarray[count];
                           }
                           for (int l=0; l<ninactive; l++) {
                                tautmp += sigma_rayleigh[wn + nwngrid*(l+nactive)] * inactive_mixratio[k+j+nlayers*l] * density[j+k] * dlarray[count];
                           }
                        }
                        if (cia == 1) {
                            for (int c=0; c<cia_npairs;c++) {
                                tautmp += sigma_cia_interp[wn + nwngrid*c] * x1_idx[c][k+j]*x2_idx[c][k+j] * density[j+k]*density[j+k] * dlarray[count];
                            }
                        }
                        //calculating mie scattering model
						if ((mie == 1) && (pressure[j] >= mie_topP) && (pressure[j] <= mie_bottomP)){
							tautmp += sigma_mie[wn] * density[j+k] *dlarray[count];
						}
                        count += 1;
                    }
                    transtot *= exp(-tautmp);

                    integral += ((planet_radius+z[j])*(1.0-transtot)*dz[j]);
                    tau[count2] = 1.0 - transtot;
                    //cout << count2 << " " << transtot << endl;
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
        delete[] ktab_interp;
        delete[] sigma_cia_interp;
        dlarray = NULL;
        ktab_interp = NULL;
        sigma_cia_interp = NULL;
    }
}
