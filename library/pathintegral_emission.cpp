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
    //cout<< "Interpolating between "<< *(bounds) << " and " << *(bounds+1) << "....." <<endl;

    /* ...and define a useful value */
    const double factor = (new_y - y_low) / (y_high - y_low);


    /* Calculate new values */
    double new_val = sig1 + ((sig2-sig1) * factor);
	//interpolation formula

    return(new_val);
}
#endif

/* Function calculating the blackbody for a given temperature */
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

//function that finds nearest vallue in array for given value
#ifndef find_nearest_H
#define find_nearest_H
void find_nearest(double * array,double value,int* index,double* nearest)
{
	double diff = abs(value - array[0]);
//	double num;
//	double nearest;
//	int index;
	int arr_size = sizeof(array)/sizeof(array[0]);

	for(int i = 0; i < arr_size; i++){
		if(diff > abs(value - array[i])){
			diff = abs(value - array[0]);
			nearest = &array[i];
			index = &i;
		}
	}
}
#endif

//function generating the correct 2D sigma_array from 3D input and temperature
#ifndef get_sigma_array_H
#define get_sigma_array_H
double** get_sigma_array(double *** sigma_input, double temperature,
		double * temperature_grid,int n_gas,int n_lambda)
{
	//takes 3d sigma_array and returns 2d sigma_array for correct temperature

	double ** sigma_output;
	sigma_output = create2Darray(n_gas,n_lambda);
	int *index;
	double *nearest;

	find_nearest(temperature_grid,temperature,index,nearest); //get the correct index for temperature grid

	for(int wl=0;wl<n_lambda;wl++){
		for(int i=0;i<n_gas;i++){
			sigma_output[i][wl] = sigma_input[*index][i][wl]; //cast to 2D array for given temperature
		}
	}
	return sigma_output;
}
#endif


//main path integral function
extern "C" {
void cpath_int_emission(double ** X,  double * rho, double * temperature,
		double * F_star, double * specgrid, double *** sigma_array_full, double * dzarray,
		const int n_lambda,const double Rp, const double Rs, const int nlayers,
		const int n_gas, const int n_sig_temp, void * FpFs_v) {
	/*function calculating the path integral for the emission. See python code for comments */


	double *FpFs = (double *)FpFs_v;

	FpFs[0] = 100;
	//declaring arrays and variables
	double I_total[n_lambda],tau_total[n_lambda];
	double BB_surf[n_lambda], exptau[n_lambda],temperature_grid[n_sig_temp];
	double tau[nlayers][n_lambda], dtau[nlayers][n_lambda];
	double ** sigma_array, *BB_layer;



	//getting temperature_grid for get_sigma_array function
	for(int t=0;t<n_sig_temp;t++){
		temperature_grid[t] =  sigma_array_full[t][0][0];
	}

	//surface layer

	for(int i=0;i<n_lambda;i++){
		cout << specgrid[i] << endl;
	}
//	cout << specgrid << endl;
//	cout << n_lambda << endl;
//	cout << temperature[0] << endl;

//	black_body(specgrid,n_lambda,temperature[0],BB_surf);
//	sigma_array = get_sigma_array(sigma_array_full,temperature[0],\
//			temperature_grid,n_gas,n_lambda);


//	for(int wl=0;wl<n_lambda;wl++){
//		for(int k = 1; k < nlayers; k++){
//			//selecting correct sigma_array for given temperature
//			if(temperature[k] != temperature[k-1]){
//				sigma_array = get_sigma_array(sigma_array_full,temperature[k],\
//							temperature_grid,n_gas,n_lambda);
//			}
//			//calculating tau for all wavelengths through layers
//			for(int i=0;i<n_gas;i++){
//				tau[0][wl] += (sigma_array[i][wl] * X[i][k] * rho[k] * dzarray[k]);
//				}
//		}
//
//		I_total[wl] += BB_surf[wl] * (exp(-1.0*tau[0][wl]));
//
//		//other layers
//		BB_layer = BB_surf;
//		sigma_array = get_sigma_array(sigma_array_full,temperature[0],\
//					temperature_grid,n_gas,n_lambda);
//
//		for(int j=1;j<nlayers;j++){
//
//			for(int k=j;k<nlayers;k++){
//				if(temperature[k] != temperature[k-1]){
//					sigma_array = get_sigma_array(sigma_array_full,temperature[k],\
//								temperature_grid,n_gas,n_lambda);
//				}
//				for(int i=0;i<n_gas;i++){
//					tau[j][wl] += (sigma_array[i][wl] *X[i][k] *rho[k] *dzarray[k]);
//					}
//
//				if(k==j){
//					dtau[j][wl] = tau[j][wl];
//				}
//
//				if(temperature[j] != temperature[j-1]){
//					black_body(specgrid,n_lambda,temperature[j],BB_layer);
//				}
//				I_total[wl] += BB_layer[wl] * (exp(-1.0*tau[j][wl]))* dtau[j][wl];
//			}
//		}
//		FpFs[wl] = (I_total[wl] /F_star[wl]) * pow((Rp/Rs),2);
//	}


}
}


