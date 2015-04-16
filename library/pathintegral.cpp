// compile as shared library
// g++ -fPIC -shared -o pathintegral.so pathintegral.cpp

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

/* Function to interpolate single values */
#ifndef interpolateValue_H
#define interpolateValue_H
double interpolateValue(double *bounds, double sig1, double sig2)
{
    /* Extract bounds... */
    const double y_low = *(bounds);
    const double y_high = *(bounds+1);
    const double new_y = *(bounds+2);
    //cout<< "Interpolating between "<< *(bounds) << " and " << *(bounds+1) << "....." <<endl;

    /* ...and define a useful value */
    const double factor = (new_y - y_low) / (y_high - y_low);


    /* Calculate new values */
    double new_val = sig1 + ((sig2-sig1) * factor);
	//interpolation formula

    return(new_val);
}
#endif

extern "C" {

void cpath_length(int nlayers, const double * zRp, void * dlarrayv) {

    double * dlarray = (double *) dlarrayv;

    int count;
    double p;

    count = 0;
    for (int j=0; j<(nlayers); j++) {	// loop through atmosphere layers, z[0] to z[nlayers]
        for (int k=1; k < (nlayers - j); k++) {
            p = pow((zRp[j]),2);
            dlarray[count] = 2.0 * (sqrt(pow((zRp[k+j]),2) - p) - sqrt(pow((zRp[k-1+j]),2) - p));
            count += 1;
        }
    }
}

 void cpath_int(
    const double ** sigma_array, const double * sigma_array_flat, int const_temp, const double * dlarray,
    const double * z, const double * dz, const double * Rsig, const double * Csig, const double ** X, const double * rho,
    double Rp, double Rs, int linecount, const int nlayers, const int n_gas, const int n_sig_temp, int include_cld,
    const double cld_lowbound, const double cld_upbound, const double * p_bar, const double * cld_sig,

    const int pressure_broadening, const double * flattened_sigma_arr,
    const double * sigma_templist, const double * sigma_preslist,
    const int nsigma_templist, const int nsigma_preslist,
    const double * pressure_array, const double temperature,

    void * absorptionv) {

    //output array to be passed back to python
    double * absorption = (double *) absorptionv;


   // setting up arrays and variables
    double tau[nlayers];
    double exptau[nlayers];
    double sigma_mol[n_gas]; // cross section for each molecule

    double Rtau, Ctau, cld_tau;
    int count;

    int t1, t2, p1[nlayers], p2[nlayers]; // temperature and pressure bounds idx for 2D interpolation of sigma array
    double F11, F12, F21, F22;
    double x, y, x1, x2, y1, y2;
    double sigma;
    int t;

    double ***sigma_array_3d = new double**[n_sig_temp];

	//if T not constant with altitude, generating 3D sigma_array from flat 1D sigma_array_flat
	//if T is constant, then use sigma_array (2d array)
    if(const_temp==0){
        for(int i =0; i<n_sig_temp; i++){
            sigma_array_3d[i] = new double*[n_gas];
           for(int j =0; j<n_gas; j++){
               sigma_array_3d[i][j] = new double[linecount];
               for(int k = 0; k<linecount;k++){
                   sigma_array_3d[i][j][k] = sigma_array_flat[(i*n_gas*linecount)+(j*linecount)+k];
               }
           }
        }
    }

    // initialise cloud quantities
	double bounds[3]={0.0}, cld_log_rho=0.0;

//    count = 0;
//    for (int j=0; j<(nlayers); j++) {	// loop through atmosphere layers, z[0] to z[nlayers]
//        for (int k=2; k < (nlayers-j); k++) {
//           // if ((k > 0) && (k < (nlayers-j))) {
//                p = pow((Rp + z[j]),2);
// 		        dlarray[count] = 2.0 * (sqrt(pow((Rp + z[k+j]),2) - p) - sqrt(pow((Rp + z[k-1+j]),2) - p));
//                count += 1;
//            //}
//        }
//    }
//

    if (pressure_broadening == 1) {
        if (const_temp==1) {
            // temperature is constant through the atmosphere (doesn't change with j)
            // find closest temperature indexes in sigma_templist
            for(int t=0;t<nsigma_templist;t++) {
                if ((temperature > sigma_templist[t]) && (temperature < sigma_templist[t+1])) {
                    t1 = t;
                    t2 = t+1;
                }
            }
        }
        // precalculate closest pressure indexes in sigma_preslist for all levels
        for (int j=0; j<(nlayers); j++) {
            for(int p=0; p<nsigma_preslist-1;p++) {
                if (p == 0 && (pressure_array[j] < sigma_preslist[p]*1e5)) {
                    // if pressure is < min(preslist), assume  p1 = p2 = min(preslist)
                    p1[j] = 0;
                    p2[j] = 1;
                    break;
                } else if ((p == nsigma_preslist-2) && (pressure_array[j] > sigma_preslist[p+1]*1e5)) {
                    // if pressure is > max(preslist), assume  p1 = p2 = max(preslist)
                    p1[j] = p;
                    p2[j] = p+1;
                    break;
                } else if ((pressure_array[j] >= sigma_preslist[p]*1e5) && (pressure_array[j] <= sigma_preslist[p+1]*1e5)) {
                    p1[j] = p;
                    p2[j] = p+1;
                    break;
                }
            }
        }
    }

    //beginning calculations
    for (int wl=0;wl < linecount; wl++) {
        count = 0;
		for (int j=0; j<(nlayers); j++) { 	// loop through atmosphere layers, z[0] to z[nlayers]

			Rtau = 0.0; 					// sum of Rayleigh optical depth
			Ctau = 0.0; 					// sum of CIA optical depth
			cld_tau = 0.0;					// sum of cloud optical depth
			tau[j] = 0.0;					// total optical depth

			for (int k=1; k < (nlayers-j); k++) { // loop through each layer to sum up path length

                // pressure is pressure_array[k]
                // temperature is temperature
                //
                // pressure bounds are p1[k] and p2[k]
                // temperature bounds are t1 and t2


			    /* Sum up taus for all gases for this path */
				for(int l=0;l<n_gas;l++) {

				    if (pressure_broadening == 1) {
                        // x = temperature
                        // y = pressure

                        // interpolate sigma_array for given pressure and temperature and molecule
                        // follow wikipedia Bilinear interpolation: http://en.wikipedia.org/wiki/Bilinear_interpolation
                        // to be changed to Hills method....

                        F11 = flattened_sigma_arr[t1*nsigma_preslist*n_gas*nlayers+p1[j]*n_gas*nlayers+l*nlayers+wl];
                        F12 = flattened_sigma_arr[t1*nsigma_preslist*n_gas*nlayers+p2[j]*n_gas*nlayers+l*nlayers+wl];
                        F21 = flattened_sigma_arr[t2*nsigma_preslist*n_gas*nlayers+p1[j]*n_gas*nlayers+l*nlayers+wl];
                        F22 = flattened_sigma_arr[t2*nsigma_preslist*n_gas*nlayers+p2[j]*n_gas*nlayers+l*nlayers+wl];

                        if (temperature > sigma_templist[nsigma_templist-1]) {
                            x = sigma_templist[nsigma_templist-1]; // if T is > max(sigma_templist) set T to max(sigma_templist)
                        } else if (temperature < sigma_templist[0]) {
                            x = sigma_templist[nsigma_templist-1]; // if T is < min(sigma_templist) set T to min(sigma_templist)
                        } else {
                            x = temperature;
                        }

                        if (pressure_array[k] > sigma_preslist[nsigma_preslist-1]) {
                            y = sigma_preslist[nsigma_preslist-1]; // if P is > max(sigma_preslist) set it to max(sigma_preslist)
                        } else if (pressure_array[k] < sigma_preslist[0]) {
                            y = sigma_preslist[0]; // if P is < min(sigma_preslist) set it to to min(sigma_preslist)
                        } else {
                            y = pressure_array[k];
                        }

                        x1 = sigma_templist[t1];
                        x2 = sigma_templist[t2];
                        y1 = sigma_preslist[p1[j]];
                        y2 = sigma_preslist[p2[j]];
                        sigma =  (1/((x2-x1)*(y2-y1))) * ( F11*(x2-x)*(y2-y)+ F21*(x-x1)*(y2-y) + F12*(x2-x)*(y-y1) + F22*(x-x1)*(y-y1));

                        // cout  << sigma << " " << F11 << " " << F12 << " " << F21 << " " << F22 << " " << x << " " << y << " " << x1 << " " << x2 << " " << y1 << " " << y2 << endl;
                        if (sigma != sigma) { // check for nans
                            sigma = 0;
                        }
                    } else {
                        if(const_temp==0){
                            sigma = sigma_array_3d[j][l][wl];
                        } else {
                            sigma = sigma_array[l][wl];
                        }


                    }
                    // cout << "sigma is " << sigma << endl;
                    tau[j] += (sigma * X[l][k+j] * rho[k+j] * dlarray[count]);
				}

				Rtau += Rsig[wl] * rho[j+k] * dlarray[count]; // calculating Rayleigh scattering optical depth
				Ctau += Csig[wl] * rho[j+k] * rho[j+k] * dlarray[count]; // calculating CIA optical depth

                //cout << dlarray[count] << endl;

//              Calculating cloud opacities if requested
                if(include_cld==1){
                    if( (p_bar[k+j]<cld_upbound) && (p_bar[k+j]>cld_lowbound) ){ // then cloud exists in this layer [k+j]

                        bounds[0]=log(cld_lowbound);
						bounds[1]=log(cld_upbound);
						bounds[2]=log(p_bar[k+j]);

                        cld_log_rho = interpolateValue(bounds,-6,-1);
                        // = log(cloud density), assuming linear decrease with decreasing log pressure
                        // following Ackerman & Marley (2001), Fig. 6
                        cld_tau += ( cld_sig[wl] * (dlarray[count]*1.0e2) * (exp(cld_log_rho)*1.0e-6) );	// convert path lenth from m to cm, and density from g m^-3 to g cm^-3
                    }
                }
				count += 1;
			}
			tau[j] += Rtau;	//adding Rayleigh scattering tau to gas tau
			tau[j] += Ctau; //adding CIA tau to gas tau
            tau[j] += cld_tau; //adding cloud tau to gas tau

            exptau[j]= exp(-tau[j]);

		}

		double integral = 0.0;
		for(int j=0; j<nlayers; j++){
		  //HOTFIX TO EQUAL TAU.CPP. TAU.CPP does not calculate the upper layer correctly
//		  if (j == nlayers-1) exptau[j] = 0.0;
		  // END OF HOTFIX
		   integral += ((Rp+z[j])*(1.0-exptau[j])*dz[j]);
		   //cout << "int" << integral << " j " << j << " dz " << dz[j] << " z " << z[j] << " exptau "  << exptau[j] << endl;
		}
		integral*=2.0;

		absorption[wl] = ((Rp*Rp) + integral) / (Rs*Rs);

    }
    if(const_temp==0){
        for(int i =0; i<n_sig_temp; i++){
            for(int j =0; j<n_gas; j++){
    //			for(int k = 0; k<n_lambda;k++){
    //		    	   delete sigma_array_3d[i][j][k];
    //		       }
                delete [] sigma_array_3d[i][j];
               }
            delete [] sigma_array_3d[i];
         }
    }
	delete [] sigma_array_3d;



}

}