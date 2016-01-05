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
                       const double * z, // altitude profile
                       const double * dz,
                       const double * density,
                       const double ** active_mixratio,
                       const double ** inactive_mixratio,
                       const double * temperature,
                       const double planet_radius,
                       const double star_radius,
                       void * absorptionv) {
            
        double * absorption = (double *) absorptionv;

        // setting up arrays and variables
        double dlarray[nlayers*nlayers];
        double sigma, sigma_l, sigma_r, sigma_interp[nwngrid*nlayers*nactive];
        double tau, exptau,  integral;
        double p;
        int count, t_idx;

        // dl array
        count = 0;
        for (int j=0; j<(nlayers); j++) {
            for (int k=1; k < (nlayers - j); k++) {
                p = pow((z[j]+planet_radius),2);
                dlarray[count] = 2.0 * (sqrt(pow((z[k+j]+planet_radius),2) - p) - sqrt(pow((z[k-1+j]+planet_radius),2) - p));
                count += 1;
            }
        }

        // interpolate sigma array to the temperature profile
        for (int j=0; j<(nlayers); j++) {
            for (int t=1; t<(sigma_ntemp); t++) {
                if ((temperature[j] > sigma_temp[t-1]) && (temperature[j] < sigma_temp[t])) {
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


        // calculate absorption
        for (int wn=0; wn < nwngrid; wn++) {
            count = 0;
            integral = 0.0;
    		for (int j=0; j<(nlayers); j++) { 	// loop through atmosphere layers, z[0] to z[nlayers]
    			tau = 0.0;
    			for (int k=0; k < (nlayers-j); k++) { // loop through each layer to sum up path length
                    // calculate optical depths due to active absorbing gases (absorption + rayleigh scattering)
    				for (int l=0;l<nactive;l++) {
                        sigma = sigma_interp[wn + nwngrid*((k+j) + l*nlayers)];
                        tau += (sigma * active_mixratio[l][k+j] * density[k+j] * dlarray[count]);
                        tau += sigma_rayleigh[wn + nwngrid*l] * active_mixratio[l][k+j] * density[j+k] * dlarray[count];
                    }
                    // calculating optical depth due inactive gases (rayleigh scattering)
                    for (int l=0; l<ninactive; l++) {
                        //cout << sigma_rayleigh[wn + nwngrid*(l+nactive)] << " " << inactive_mixratio[l][k+j] << " " << density[j+k] << " " << dlarray[count] << endl;
                        tau += sigma_rayleigh[wn + nwngrid*(l+nactive)] * inactive_mixratio[l][k+j] * density[j+k] * dlarray[count];
                    }
                }
                exptau = exp(-tau);
		        integral += ((planet_radius+z[j])*(1.0-exptau)*dz[j]);
            }
            integral *= 2.0;
            absorption[wn] = ((planet_radius*planet_radius) + integral) / (star_radius*star_radius);
        }
    }
}