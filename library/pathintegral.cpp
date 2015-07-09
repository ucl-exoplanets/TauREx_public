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
    ////cout<< "Interpolating between "<< *(bounds) << " and " << *(bounds+1) << "....." <<endl;

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
    double Rp, double Rs, int nlambda, const int nlayers, const int n_gas, const int n_sig_temp, int include_cld,
    const double cld_lowbound, const double cld_upbound, const double * p_bar, const double * cld_sig,

    const int pressure_broadening, const double * flattened_sigma_arr,
    const double * sigma_templist, const double * sigma_preslist,
    const int nsigma_templist, const int nsigma_preslist,
    const double * pressure_array, const double * temperature_array,

    const double temperature,

    void * absorptionv) {


    //output array to be passed back to python
    double * absorption = (double *) absorptionv;



   // setting up arrays and variables
    double tau[nlayers];
    double exptau[nlayers];
    double sigma_mol[n_gas]; // cross section for each molecule

    double Rtau, Ctau, cld_tau;
    int count;

    int t0[nlayers], t1, t2, p1[nlayers], p2[nlayers]; // temperature and pressure bounds idx for 2D interpolation of sigma array
    double F11, F12, F21, F22, F1, F2;
    double x, y, x1, x2, y1, y2;
    double sigma;
    int t;

    double ***sigma_array_3d = new double**[n_sig_temp];

//   //cout << " n_sig_temp " << n_sig_temp << " n_gas " << n_gas << " nlambda " <<   nlambda << endl;

	//if T not constant with altitude, generating 3D sigma_array from flat 1D sigma_array_flat
	//if T is constant, then use sigma_array (2d array)
    if(const_temp==0){
        for(int i =0; i<n_sig_temp; i++){

            sigma_array_3d[i] = new double*[n_gas];
           for(int j =0; j<n_gas; j++){
               sigma_array_3d[i][j] = new double[nlambda];
               for(int k = 0; k<nlambda;k++){

//                   //cout << (i*n_gas*nlambda)+(j*nlambda)+k << endl;
                   sigma_array_3d[i][j][k] = sigma_array_flat[(i*n_gas*nlambda)+(j*nlambda)+k];
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
            // precalculate closest temperature indexes in sigma_templist
            if (nsigma_templist == 1) {
                t1 = 0;
                t2 = 0;
            } else {

                for(int t=0;t<nsigma_templist;t++) {
                    if ((temperature >= sigma_templist[t]) && (temperature < sigma_templist[t+1])) {
                        t1 = t;
                        t2 = t+1;
                    }
                }
            }
        }

        // precalculate closest pressure indexes in sigma_preslist for all levels
        for (int j=0; j<(nlayers); j++) {
            for(int p=0; p<nsigma_preslist;p++) {
                if (p == 0 && (pressure_array[j] < sigma_preslist[p]*1e5)) {
                    // if pressure is < min(preslist), assume  p1 = p2 = min(preslist)
                    p1[j] = 0;
                    p2[j] = 1;
                    break;
                } else if ((p == nsigma_preslist-1) && (pressure_array[j] > sigma_preslist[p]*1e5)) {
                    // if pressure is > max(preslist), assume  p1 = p2 = max(preslist)
                    p1[j] = p-1;
                    p2[j] = p;
                    break;
                } else if ((pressure_array[j] >= sigma_preslist[p]*1e5) && (pressure_array[j] <= sigma_preslist[p+1]*1e5)) {
                    p1[j] = p;
                    p2[j] = p+1;
                    break;
                }
            }
        }
    } else {
        // pressure broadening is off

//        //cout << "nsigma_templist" << nsigma_templist << endl;

        if (const_temp==0) {
            // variable T profile. Precalculate T indexes in sigma_array for each layer


            for (int j=0; j<(nlayers); j++) { // loop through layers
//        //cout << "test0" << endl;
                //cout << " temp at j " << j << " is " << temperature_array[j] << endl;
                for(int t=0; t<nsigma_templist;t++) { // loop through sigma T grid

                    //cout << " is T " << temperature_array[j]  << " betweem " << sigma_templist[t] << " and " << sigma_templist[t+1] << endl;
                    if ((temperature_array[j] >= sigma_templist[t]) && (temperature_array[j] <= sigma_templist[t+1])) {
                        // find closest temperatures (upper and lower) for T in layer in sigma T grid
                        if ((temperature_array[j]-sigma_templist[t]) < (sigma_templist[t+1]-temperature_array[j])) {
                            // set idx to closest temperature between upper/lower bounds
                            t0[j] = t;
                            //cout << " t0 j " << j << " " << t << endl;
                        } else {
                            t0[j] = t+1;
                        }
                    }
                }
            }
        }
    }

//    //cout << "  end " << endl;

    //beginning calculations
    for (int wl=0;wl < nlambda; wl++) {
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



                        if (temperature > sigma_templist[nsigma_templist-1]) {
                            x = sigma_templist[nsigma_templist-1]; // if T is > max(sigma_templist) set T to max(sigma_templist)
                        } else if (temperature < sigma_templist[0]) {
                            x = sigma_templist[nsigma_templist-1]; // if T is < min(sigma_templist) set T to min(sigma_templist)
                        } else {
                            x = temperature;
                        }

                        if (pressure_array[k]*1.0e-5 > sigma_preslist[nsigma_preslist-1]) {
                            y = sigma_preslist[nsigma_preslist-1]; // if P is > max(sigma_preslist) set it to max(sigma_preslist)
                        } else if (pressure_array[k]*1.0e-5 < sigma_preslist[0]) {
                            y = sigma_preslist[0]; // if P is < min(sigma_preslist) set it to to min(sigma_preslist)
                        } else {
                            y = pressure_array[k]*1.0e-5;

                        }

                        //apply linear or bilinear interpolation:

                        if (t1 == t2) {
                            // no interpolation in temperature needed
                            if (p1[k] == p2[l]) {
                                // no interpolation needed
                                sigma = flattened_sigma_arr[t1*nsigma_preslist*n_gas*nlambda+p1[k]*n_gas*nlambda+l*nlambda+wl];
                            } else {
                                // interpolate in P
                                F1 = flattened_sigma_arr[t1*nsigma_preslist*n_gas*nlambda+p1[k]*n_gas*nlambda+l*nlambda+wl];
                                F2 = flattened_sigma_arr[t1*nsigma_preslist*n_gas*nlambda+p2[k]*n_gas*nlambda+l*nlambda+wl];
                                x1 = sigma_preslist[p1[k]];
                                x2 = sigma_preslist[p2[k]];
                                x = pressure_array[k]*1.0e-5;
                                sigma = F1 + (F2-F1)*(x-x1)/(x2-x1);
//                                //cout << k << " interpolate in P " << sigma << " F1 " << F1 << " F2 " << F2 << " x1 " << x1 << " x2 " << x2 << " x " << x << endl;
//
//                                //cout << sigma <<  " k "  << k <<  " t1 "  << t1 <<  " T_t1 "  << sigma_templist[t1] <<  " t2 "  << t2 <<  " T_t2 "  << sigma_templist[t2]
//                                     <<  " p1 "  << p1[k] <<  " P_p1 "  << sigma_preslist[p1[k]] <<  " p2 "  << p2[k] <<  " P_p2 "  << sigma_preslist[p2[k]]
//                                     <<  " P "  << y <<  " T "  << x << endl;;
//
//                                //cout << "sigma "  << sigma   << " F11 "  << F11 << " F12 "  << F12 << " F21 "  << F21 <<  " F22 "  << F22 <<  endl;
                            }
                        } else {
                            if (p1[k] == p2[l]) {
                                // linear interpolation in T
                                F1 = flattened_sigma_arr[t1*nsigma_preslist*n_gas*nlambda+p1[k]*n_gas*nlambda+l*nlambda+wl];
                                F2 = flattened_sigma_arr[t2*nsigma_preslist*n_gas*nlambda+p1[k]*n_gas*nlambda+l*nlambda+wl];
                                x1 = sigma_templist[t1];
                                x2 = sigma_templist[t2];
                                x = sigma_templist[k];
                                sigma = F1 + (F2-F1)*(x-x1)/(x2-x1);
                            } else {
                                // bilinear interpolation in both T and P
                                x1 = sigma_templist[t1];
                                x2 = sigma_templist[t2];
                                y1 = sigma_preslist[p1[k]];
                                y2 = sigma_preslist[p2[k]];
                                F11 = flattened_sigma_arr[t1*nsigma_preslist*n_gas*nlambda+p1[k]*n_gas*nlambda+l*nlambda+wl];
                                F12 = flattened_sigma_arr[t1*nsigma_preslist*n_gas*nlambda+p2[k]*n_gas*nlambda+l*nlambda+wl];
                                F21 = flattened_sigma_arr[t2*nsigma_preslist*n_gas*nlambda+p1[k]*n_gas*nlambda+l*nlambda+wl];
                                F22 = flattened_sigma_arr[t2*nsigma_preslist*n_gas*nlambda+p2[k]*n_gas*nlambda+l*nlambda+wl];
                                sigma =  (1/((x2-x1)*(y2-y1))) * ( F11*(x2-x)*(y2-y)+ F21*(x-x1)*(y2-y) + F12*(x2-x)*(y-y1) + F22*(x-x1)*(y-y1));
                            }
                        }

//                        //cout << sigma <<  " k "  << k <<  " t1 "  << t1 <<  " T_t1 "  << sigma_templist[t1] <<  " t2 "  << t2 <<  " T_t2 "  << sigma_templist[t2]
//                             <<  " p1 "  << p1[k] <<  " P_p1 "  << sigma_preslist[p1[k]] <<  " p2 "  << p2[k] <<  " P_p2 "  << sigma_preslist[p2[k]]
//                             <<  " P "  << y <<  " T "  << x << endl;;

//                        //cout << "sigma "  << sigma   << " F11 "  << F11 << " F12 "  << F12 << " F21 "  << F21 <<  " F22 "  << F22 <<  endl;

                        if (sigma != sigma) { // check for nans
                            sigma = 0;
                        }
                    } else {
                        if(const_temp==0){

                            // temperature is temperature_array[k]

                            // l is the molecule number, OK --> second array
                            // need to find closest temperature id to tempgrid for j

                            // get T profile first
                            // then get T for layer j
                            // then find closest sigma array to T, using tempgrid

                            // temperature idx for layer j

                            //cout << "wl" << wl << " t " << t0[j] << endl;

                            sigma = sigma_array_3d[t0[j]][l][wl];

                        } else {
                            sigma = sigma_array[l][wl];
                        }
                    }
                    tau[j] += (sigma * X[l][k+j] * rho[k+j] * dlarray[count]);
				}

				Rtau += Rsig[wl] * rho[j+k] * dlarray[count]; // calculating Rayleigh scattering optical depth
				Ctau += Csig[wl] * rho[j+k] * rho[j+k] * dlarray[count]; // calculating CIA optical depth

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
            //cout << " exptau " << exptau[j] << endl;
		}



		double integral = 0.0;
		for(int j=0; j<nlayers; j++){
		  //HOTFIX TO EQUAL TAU.CPP. TAU.CPP does not calculate the upper layer correctly
//		  if (j == nlayers-1) exptau[j] = 0.0;
		  // END OF HOTFIX

		   integral += ((Rp+z[j])*(1.0-exptau[j])*dz[j]);
		   //cout << "int " << integral << " j " << j << " dz " << dz[j] << " z " << z[j] << " exptau "  << exptau[j] << endl;
		}
		integral*=2.0;

        //cout << "integral ole" << endl;
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