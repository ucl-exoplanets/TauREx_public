// compile as shared library
// g++ -fPIC -shared -o pathintegral_emission.so pathintegral_emission.cpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdlib.h>

using namespace std;




//function creating 2D pointer array
#ifndef create2Darray_H
#define create2Darray_H
double** create2Darray(unsigned height, unsigned width)
    {
      double** array2D = 0;
      array2D = new double*[height];

      for (int h = 0; h < height; h++)
      {
            array2D[h] = new double[width];
            for (int w = 0; w < width; w++)
            {
                  // fill in some initial values
                  // (filling in zeros would be more logic, but this is just for the example)
                  array2D[h][w] = 0;
            }
      }

      return array2D;
    }

#endif


/* Function to interpolate single values */
#ifndef interpolateValue_H
#define interpolateValue_H
double interpolateValue(double *bounds, double sig1, double sig2)
{
    /* Extract bounds... */
    const double y_low = *(bounds);
    const double y_high = *(bounds+1);
    const double new_y = *(bounds+2);

    /* ...and define a useful value */
    const double factor = (new_y - y_low) / (y_high - y_low);


    /* Calculate new values */
    double new_val = sig1 + ((sig2-sig1) * factor);
	//interpolation formula

    return(new_val);
}
#endif

/* Function calculating the blackbody for an array of wavelengths and a given temperature */
#ifndef black_body_H
#define black_body_H
void black_body(double *lambda, int n_lambda, double temperature, double *BB)
{
/* Small function calculating plank black body
 * input: microns, kelvin
 * output: W/m^2/micron
 */
	double h,c,k,pi, exponent;
//	double BB[n_lambda];

    h = 6.62606957e-34;
    c = 299792458;
    k = 1.3806488e-23;
    pi= 3.14159265359;

    for (int wl=0; wl < n_lambda; wl++){ //loop through wavelengths
    	exponent = exp((h * c) / (lambda[wl]*1e-6 * k * temperature));
    	BB[wl] = (pi * (2.0*h*pow(c,2))/pow((lambda[wl]*1e-6),5) * (1.0/(exponent -1)));
    	BB[wl] *= 1e-6;
    }
}
#endif

/* Function calculating the blackbody for one wavelength and a given temperature */
#ifndef black_body_2_H
#define black_body_2_H
void black_body_2(double lambda, double temperature, double &BB)
{
/* Small function calculating plank black body
 * input: microns, kelvin
 * output: W/m^2/micron
 */
	double h,c,k,pi, exponent;
//	double BB[n_lambda];

    h = 6.62606957e-34;
    c = 299792458;
    k = 1.3806488e-23;
    pi= 3.14159265359;

//    for (int wl=0; wl < n_lambda; wl++){ //loop through wavelengths
    	exponent = exp((h * c) / (lambda*1e-6 * k * temperature));
    	BB = (pi * (2.0*h*pow(c,2))/pow((lambda*1e-6),5) * (1.0/(exponent -1)));
    	BB *= 1e-6;
//    }
}
#endif




//function that finds nearest vallue in array for given value
#ifndef find_nearest_H
#define find_nearest_H
void find_nearest(double * array,int array_size,double value,int& index,double& nearest)
{
	double diff = abs(value - array[0]);
	index = 0;

	for(int i = 0; i < array_size; i++){
		if(diff > abs(value - array[i])){
			diff = abs(value - array[i]);
			nearest = array[i];
			index = i;
		}

	}
}
#endif

//function generating the correct 2D sigma_array from 3D input and temperature
#ifndef get_sigma_array_H
#define get_sigma_array_H
double** get_sigma_array(double *** sigma_input, double temperature,
		double * sig_tempgrid,int n_gas,int n_lambda,int n_sig_temp)
{
	//takes 3d sigma_array and returns 2d sigma_array for correct temperature

	double ** sigma_output;
	sigma_output = create2Darray(n_gas,n_lambda);
	int index;
	double nearest;

	find_nearest(sig_tempgrid,n_sig_temp,temperature,index,nearest); //get the correct index for temperature grid

	for(int wl=0;wl<n_lambda;wl++){
		for(int i=0;i<n_gas;i++){
			cout << index << " " << i << " " << wl << endl;
			cout << sigma_input[index][i][wl] << endl;
			sigma_output[i][wl] = sigma_input[index][i][wl]; //cast to 2D array for given temperature
		}
	}
	return sigma_output;
}
#endif


//main path integral function
extern "C" {
void cpath_int_emission(double ** X,  double * rho, double * temperature,
		double * F_star, double * specgrid, double * sigma_array_flat, double * dzarray,
		const int n_lambda,const double Rp, const double Rs, const int nlayers,
		const int n_gas, double* sig_tempgrid, const int n_sig_temp, void * FpFs_v) {
	/*function calculating the path integral for the emission. See python code for comments */

	//casting output vector
	double *FpFs = (double *)FpFs_v;

	//generating 3D sigma_array from flat 1D sigma_array_flat
	double ***sigma_array_3d = new double**[n_sig_temp];
	for(int i =0; i<n_sig_temp; i++){
		sigma_array_3d[i] = new double*[n_gas];
	   for(int j =0; j<n_gas; j++){
		   sigma_array_3d[i][j] = new double[n_lambda];
	       for(int k = 0; k<n_lambda;k++){
	    	   sigma_array_3d[i][j][k] = sigma_array_flat[(i*n_gas*n_lambda)+(j*n_lambda)+k];
	       }
	   }
	}


//	//declaring arrays and variables
	double exptau[n_lambda], BB_layer[n_lambda],tau_total[n_lambda],BB_surf[n_lambda];
//	double tau[nlayers],dtau[nlayers];
	double sig_nearest_temp=0, BB_lambda, I_total,tau, dtau;
	int sig_index=0;

	//surface layer
//	black_body(specgrid,n_lambda,temperature[0],BB_surf);
//	find_nearest(sig_tempgrid,n_sig_temp,temperature[0],sig_index,sig_nearest_temp);

	//calculating emission spectrum per lambda
	for(int wl=0;wl<n_lambda;wl++){
		tau = 0.0;
		BB_lambda = 0.0;
		I_total = 0.0;

		black_body_2(specgrid[wl],temperature[0],BB_lambda);
		find_nearest(sig_tempgrid,n_sig_temp,temperature[0],sig_index,sig_nearest_temp);
		for(int k = 0; k < nlayers; k++){
			//selecting correct sigma_array for given temperature
			if(temperature[k] != temperature[k-1]){
				find_nearest(sig_tempgrid,n_sig_temp,temperature[k],sig_index,sig_nearest_temp);
			}
			//calculating tau through layers
			for(int i=0;i<n_gas;i++){
				tau += (sigma_array_3d[sig_index][i][wl] * X[i][k] * rho[k] * dzarray[k]);
				}
		}

		I_total += BB_lambda * (exp(-1.0*tau));

		//other layers
		find_nearest(sig_tempgrid,n_sig_temp,temperature[0],sig_index,sig_nearest_temp);
//		black_body_2(specgrid[wl],temperature[0],BB_lambda);
		for(int j=1;j<nlayers;j++){
			tau = 0.0;
			dtau = 0.0;
			for(int k=j;k<nlayers;k++){
				if(temperature[k] != temperature[k-1]){
					find_nearest(sig_tempgrid,n_sig_temp,temperature[k],sig_index,sig_nearest_temp);
				}
				for(int i=0;i<n_gas;i++){
					tau += (sigma_array_3d[sig_index][i][wl] *X[i][k] *rho[k] *dzarray[k]);
				}

				if(k==j){
					dtau = tau;
				}

				if(temperature[j] != temperature[j-1]){
					black_body_2(specgrid[wl],temperature[j],BB_lambda);
				}
			}
			I_total += (BB_lambda * (exp(-1.0*tau))* dtau);
		}
		FpFs[wl] = (I_total /F_star[wl]) * pow((Rp/Rs),2);
	}

	/*--------------------------------------------------------------------------------*/
	//deallocating pointeres and pointer arrays that are not passed back to python

	for(int i =0; i<n_sig_temp; i++){
		for(int j =0; j<n_gas; j++){
//			for(int k = 0; k<n_lambda;k++){
//		    	   delete sigma_array_3d[i][j][k];
//		       }
			delete [] sigma_array_3d[i][j];
		   }
		delete [] sigma_array_3d[i];
		}
	delete [] sigma_array_3d;

//	for (int i=1;i<n_lambda;i++){
//		delete [] exptau[i]; delete [] BB_layer[i]; delete [] tau_total[i];
//		delete [] BB_surf[i];
//	}
//	delete [] exptau; delete [] BB_layer; delete [] tau_total; delete BB_surf;
//
//	for (int i=1;i<nlayers;i++){
//		delete [] tau[i]; delete [] dtau[i];
//	}
//	delete [] tau; delete dtau;
//
//	delete [] sig_nearest_temp; delete [] BB_lambda; delete [] I_total;delete [] sig_index;



}
}


