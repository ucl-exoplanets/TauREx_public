// compile as shared library
// g++ -fPIC -shared -o pathintegral.so pathintegral.cpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#include <stdlib.h>

using namespace std;

extern "C" {
void cpath_int(const double ** sigma_array,const double * dlarray, const double * z,
    const double * dz, const double * Rsig, const double * Csig, const double ** X, const double * rho,
    double Rp, double Rs, int linecount, int nlayers, int n_gas, void * absorptionv) {



//    cout<< "C1" <<endl;
    //output array to be passed back to python
    double * absorption = (double *) absorptionv;
//   cout<< "C2" <<endl;
   
   // setting up arrays and variables
    float tau[nlayers], exptau[nlayers];
    float Rtau, Ctau;
    int count;

//    cout<< "C3" <<endl;

    //beginning calculations
    for (int wl=0;wl < linecount; wl++)		// loop through wavelengths
	{
//		/* Calculate scattering cross-sections (wavelength dependence) */
//		Csig = 0.0;			// CIA cross-section
//		Csig += scatterCIA(cia_coeffs[wl],atmos.fraction[0]);

//        cout<< "C4" <<endl;

        count = 0;
		/* Calculating optical depth */
		for (int j=0; j<(nlayers); j++) 	// loop through atmosphere layers, z[0] to z[nlayers]
		{
		/* Calculating optical depths */
			Rtau = 0.0; 					// sum of Rayleigh optical depth
			Ctau = 0.0; 					// sum of CIA optical depth
//			cld_tau = 0.0;					// sum of cloud optical depth
			tau[j] = 0.0;					// total optical depth

//            cout<< "C5" <<endl;

			for (int k=1; k < (nlayers-j); k++) // loop through each layer to sum up path length
			{
			 /* Sum up taus for all gases for this path */
				for(int l=0;l<n_gas;l++) tau[j] += (sigma_array[l][wl] * X[l][k+j] * rho[k+j] * dlarray[count]);
				Rtau += Rsig[wl] * rho[j+k] * dlarray[count]; // calculating Rayleigh scattering optical depth
				Ctau += Csig[wl] * rho[j+k] * rho[j+k] * dlarray[count]; // calculating CIA optical depth
				count += 1;

//				cout<< "C6" <<endl;
			}
			tau[j] += Rtau;	//adding Rayleigh scattering tau to gas tau
			tau[j] += Ctau; //adding CIA tau to gas tau
			exptau[j]= exp(-tau[j]);
		}
//		cout<< "C7" <<endl;
		
		double integral = 0.0;
		for(int j=0; j<nlayers; j++){
		  //HOTFIX TO EQUAL TAU.CPP. TAU.CPP does not calculate the upper layer correctly
//		  if (j == nlayers-1) exptau[j] = 0.0;
		  // END OF HOTFIX
		   integral += ((Rp+z[j])*(1.0-exptau[j])*dz[j]);
//		   cout << j << ' '<< dz[j] << ' ' << z[j] << ' ' << exptau[j] << endl;
		}
//		cout<< "C8" <<endl;
		integral*=2.0;

		absorption[wl] = ((Rp*Rp) + integral) / (Rs*Rs);

//		free(absorptionv);

    }
}

}