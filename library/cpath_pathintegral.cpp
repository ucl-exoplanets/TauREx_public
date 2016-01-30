/* 

    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Forward model for transmission

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

    Compile with:   g++ -fPIC -shared -o cpath_pathintegral.so cpath_pathintegral.cpp

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
                       const double * sigma_array,
                       const double * sigma_temp,
                       const int sigma_ntemp,
                       const double * sigma_rayleigh,
                       const int cia_npairs,
                       const int * cia_idx,
                       const int cia_nidx,
                       const double * sigma_cia,
                       const double * sigma_cia_temp,
                       const int sigma_cia_ntemp,
                       const int clouds,
                       const double * sigma_cloud,
                       const double * clouds_density,
                       const double * density,
                       const double * z,
                       const double * dz,
                       const double ** active_mixratio,
                       const double ** inactive_mixratio,
                       const double * temperature,
                       const double planet_radius,
                       const double star_radius,
                       void * absorptionv) {


        double * absorption = (double *) absorptionv;

        // setting up arrays and variables
        double* dlarray = new double[nlayers*nlayers];
        double* sigma_interp = new double[nwngrid*nlayers*nactive];
        double* sigma_cia_interp = new double[nwngrid * nlayers * cia_npairs];
        double sigma, sigma_l, sigma_r;
        double x1_idx[cia_npairs][nlayers];
        double x2_idx[cia_npairs][nlayers];
        double tau, exptau,  integral;
        double cld_factor, cld_rho;
        double p;
        int count, t_idx;

        // dl array
        count = 0;
        for (int j=0; j<(nlayers); j++) {
            for (int k=1; k < (nlayers - j); k++) {
                p = pow((z[j]+planet_radius),2);
                dlarray[count] = 2.0 * (sqrt(pow((z[k+j]+planet_radius),2) - p) - sqrt(pow((z[k-1+j]+planet_radius),2) - p));
                //cout << dlarray[count] << " " << p << " " << count << endl;
                count += 1;
            }
        }

        // interpolate sigma array to the temperature profile
        for (int j=0; j<nlayers; j++) {

            //cout << temperature[j] << " " << sigma_temp[sigma_ntemp-1] << endl;

            if (temperature[j] > sigma_temp[sigma_ntemp-1]) {
                for (int wn=0; wn<nwngrid; wn++) {
                    for (int l=0;l<nactive;l++) {
                        sigma_interp[wn + nwngrid*(j + l*nlayers)] = sigma_array[wn + nwngrid*(sigma_ntemp-1 + sigma_ntemp*(j + l*nlayers))];
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
                                //cout << sigma_l << " " << sigma_r << " " << temperature[j] << endl;
                                sigma_interp[wn + nwngrid*(j + l*nlayers)] = sigma;
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
                            sigma = sigma_l + (sigma_r-sigma_l)*(temperature[j]-sigma_temp[t-1])/(sigma_temp[t]-sigma_temp[t-1]);
                            sigma_cia_interp[wn +  nwngrid*(j + l*nlayers)] = sigma;
                        }
                    }
                }
            }
        }

        // get mixing ratio of individual molecules in the collision induced absorption (CIA) pairs
        for (int c=0; c<cia_npairs;c++) {
            if (cia_idx[c*2] >= nactive) {
                for (int j=0;j<nlayers;j++) { x1_idx[c][j] = inactive_mixratio[cia_idx[c*2]-nactive][j]; }
            } else {
                for (int j=0;j<nlayers;j++) { x1_idx[c][j] = active_mixratio[cia_idx[c*2]][j]; }
            }
            if (cia_idx[c*2] >= nactive) {
                for (int j=0;j<nlayers;j++) { x2_idx[c][j] = inactive_mixratio[cia_idx[c*2+1]-nactive][j]; }
            } else {
                for (int j=0;j<nlayers;j++) { x2_idx[c][j] = active_mixratio[cia_idx[c*2+1]][j]; }
            }
        }

        // calculate absorption
        for (int wn=0; wn < nwngrid; wn++) {
            count = 0;
            integral = 0.0;
    		for (int j=0; j<(nlayers); j++) { 	// loop through atmosphere layers, z[0] to z[nlayers]
    			tau = 0.0;
    			for (int k=1; k < (nlayers-j); k++) { // loop through each layer to sum up path length
                    // calculate optical depths due to active absorbing gases (absorption + rayleigh scattering)
    				for (int l=0;l<nactive;l++) {
                        sigma = sigma_interp[wn + nwngrid*((k+j) + l*nlayers)];
                        tau += (sigma * active_mixratio[l][k+j] * density[k+j] * dlarray[count]);
                        //cout << " j " << j  << " k " << k  << " count " << count << " sigma " << sigma << " active_mixratio " << active_mixratio[l][k+j] << " density " << density[k+j] << " dlarray " << dlarray[count] << endl;
                        //cout << " j " << j  << " k " << k  << " count " << count << " sigma_rayleigh " << sigma_rayleigh[wn + nwngrid*l] << " active_mixratio " << active_mixratio[l][k+j] << " density " << density[k+j] << " dlarray " << dlarray[count] << endl;
                        tau += sigma_rayleigh[wn + nwngrid*l] * active_mixratio[l][k+j] * density[j+k] * dlarray[count];
                    }

                    // calculating optical depth due inactive gases (rayleigh scattering)
                    for (int l=0; l<ninactive; l++) {
                        //cout << sigma_rayleigh[wn + nwngrid*(l+nactive)] << " " << inactive_mixratio[l][k+j] << " " << density[j+k] << " " << dlarray[count] << endl;
                        tau += sigma_rayleigh[wn + nwngrid*(l+nactive)] * inactive_mixratio[l][k+j] * density[j+k] * dlarray[count];
                    }
                    // calculating optical depth due to collision induced absorption
                    for (int c=0; c<cia_npairs;c++) {
                        tau += sigma_cia[wn + nwngrid*c] * x1_idx[c][k+j]*x1_idx[c][k+j] * density[j+k]*density[j+k] * dlarray[count];
                    }

                    // calculate optical depth due to clouds
                    if (clouds == 1) {
                        // = log(cloud density), assuming linear decrease with decreasing in log pressure
                        // following Ackerman & Marley (2001), Fig. 6
                        tau += sigma_cloud[wn] * dlarray[count]*1.0e2 * clouds_density[j+k] * 1.0e-6; // convert path lenth from m to cm, and density from g m^-3 to g cm^-3
                    }


                    count += 1;
                }
                exptau = exp(-tau);
		        integral += ((planet_radius+z[j])*(1.0-exptau)*dz[j]);
            }
            integral *= 2.0;
            absorption[wn] = ((planet_radius*planet_radius) + integral) / (star_radius*star_radius);
            //cout << absorption[wn] << endl;
        }
    }
}