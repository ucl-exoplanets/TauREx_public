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

void cpath_int(const double ** sigma_array, const double * dlarray, const double * z,
    const double * dz, const double * Rsig, const double * Csig, const double ** X, const double * rho,
    double Rp, double Rs, int linecount, int nlayers, int n_gas, int include_cld, const double cld_lowbound,
    const double cld_upbound, const double * p_bar, const double * cld_sig, void * absorptionv) {

    //output array to be passed back to python
    double * absorption = (double *) absorptionv;

   // setting up arrays and variables
    double tau[nlayers], exptau[nlayers];
    //double  dlarray[nlayers*nlayers];
    //double p
    double Rtau, Ctau, cld_tau;
    int count;

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

    //beginning calculations
    for (int wl=0;wl < linecount; wl++)		// loop through wavelengths
	{
//		/* Calculate scattering cross-sections (wavelength dependence) */
//		Csig = 0.0;			// CIA cross-section
//		Csig += scatterCIA(cia_coeffs[wl],atmos.fraction[0]);

        count = 0;
		/* Calculating optical depth */
		for (int j=0; j<(nlayers); j++) 	// loop through atmosphere layers, z[0] to z[nlayers]
		{
		/* Calculating optical depths */
			Rtau = 0.0; 					// sum of Rayleigh optical depth
			Ctau = 0.0; 					// sum of CIA optical depth
			cld_tau = 0.0;					// sum of cloud optical depth
			tau[j] = 0.0;					// total optical depth

			for (int k=1; k < (nlayers-j); k++) // loop through each layer to sum up path length
			{
			 /* Sum up taus for all gases for this path */
				for(int l=0;l<n_gas;l++) {
    			    tau[j] += (sigma_array[l][wl] * X[l][k+j] * rho[k+j] * dlarray[count]);
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


}

}