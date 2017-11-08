/* 

    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Forward model for transmission

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

    For both openmp and single core versions compile with g++:

             g++ -fPIC -shared -o ctypes_pathintegral_transmission_xsec.so ctypes_pathintegral_transmission_xsec.cpp
             g++ -fPIC -shared -fopenmp -o ctypes_pathintegral_transmission_parallel_xsec.so ctypes_pathintegral_transmission_xsec.cpp

    or Intel compiler:

             icc -fPIC -shared -o ctypes_pathintegral_transmission_xsec.so ctypes_pathintegral_transmission_xsec.cpp
             icc -fPIC -shared -openmp -o ctypes_pathintegral_transmission_parallel_xsec.so ctypes_pathintegral_transmission_xsec.cpp
 
    For a PGI compiler debugging:
 
             pgc++ -fast -ta=tesla -Minfo=all ctypes_pathintegral_transmission_xsec.cpp

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

//    	double mie_topP;
//    	mie_topP = 1e-5;

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
        double p;
        int count;

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
                            #pragma acc data copyin(nwngrid,nactive,sigma_ntemp,nlayers,temperature,sigma_temp)
                            #pragma acc data copyout(sigma_interp)
                            #pragma acc kernels
                            {
                            #pragma acc loop tile(32,4) device_type(nvidia)
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
                        if ((temperature[j] > sigma_cia_temp[t-1]) && (temperature[j] < sigma_cia_temp[t])) {
                            #pragma acc data copyin(nwngrid,cia_npairs,sigma_cia_ntemp,nlayers,temperature,sigma_cia_temp)
                            #pragma acc data copyout(sigma_cia_interp)
                            #pragma acc kernels
                            {
                            #pragma acc loop tile(32,4) device_type(nvidia)
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
        //#pragma omp parallel for schedule(dynamic) private(tautmp, sigma, count, integral, exptau)
        #pragma acc data copyin(nwngrid,nlayers,clouds,pressure,cloud_topP,planet_radius,star_radius)
        #pragma acc data copyin(dz,tau,nactive,sigma_interp,active_mixratio,density,dlarray)
        #pragma acc data copyin(rayleigh,sigma_rayleigh,ninactive,inactive_mixratio,cia)
        #pragma acc data copyin(cia_npairs,sigma_cia)
//        #pragma acc data copyin(x1_idx,x2_idx)
        #pragma acc data copyin(mie,mie_topP,mie_bottomP,sigma_mie)
        #pragma acc data create(exptau,integral,tau,absorption)
        #pragma acc kernels
        {
        #pragma acc loop tile(32,4) device_type(nvidia)
            for (int wn=0; wn < nwngrid; wn++) {
                count = 0;
                integral = 0.0;
                //cout << " integral 0 " << integral << endl;
                //cout << " count 0 " << count << endl;
                for (int j=0; j<(nlayers); j++) {     // loop through atmosphere layers, z[0] to z[nlayers]
                    tautmp = 0.0;
                    if ((clouds == 1) && (pressure[j] >= cloud_topP)) {
                        //cout << j << " YES " << pressure[j] << endl;
                        for (int k=0; k < (nlayers-j); k++) { // loop through each layer to sum up path length
                            count += 1;
                        }
                        //exptau = exp(-tautmp);
                        integral += ((planet_radius+z[j])*(1.0)*dz[j]);
                        tau[wn + j*nwngrid] = 1.0;
                        //cout << count << " j " << j << " z " << z[j] << " dz " << dz[j] << " exptau  " << exptau << " integral " << integral << endl;
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
                            //calculating mie scattering model
                            if ((mie == 1) && (pressure[j] >= mie_topP) && (pressure[j] <= mie_bottomP)){
                                tautmp += sigma_mie[wn] * density[j+k] *dlarray[count];
                            }

                            count += 1;
                        }
                        exptau = exp(-tautmp);
                        integral += ((planet_radius+z[j])*(1.0-exptau)*dz[j]);
                        tau[wn + j*nwngrid] =  exptau;
                        //cout << count << " j " << j << " z " << z[j] << " dz " << dz[j] << " exptau  " << exptau << " integral " << integral << endl;
                    }
                }
                integral *= 2.0;
                //cout << integral << endl;
                absorption[wn] = ((planet_radius*planet_radius) + integral) / (star_radius*star_radius);
                //cout << wn << " " << absorption[wn] << endl;
            }
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
