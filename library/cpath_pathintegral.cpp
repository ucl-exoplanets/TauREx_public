/* 
 
 TauREx v2 - Transmission Pathintegral

 Compile with:
 
 g++ -fPIC -shared -o cpath_pathintegral.so cpath_pathintegral.cpp

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
    
    void path_length(int nlayers, const double * zRp, void * dlarrayv) {
        
        double * dlarray = (double *) dlarrayv;
        
        int count;
        double p;
        
        count = 0;
        for (int j=0; j<(nlayers); j++) {
            for (int k=1; k < (nlayers - j); k++) {
                p = pow((zRp[j]),2);
                dlarray[count] = 2.0 * (sqrt(pow((zRp[k+j]),2) - p) - sqrt(pow((zRp[k-1+j]),2) - p));
                //cout << count << " " << j << " " << zRp[j] << " " << dlarray[count] << " " << zRp[k+j] << " " <<  pow((zRp[k+j]),2) << " " << p << endl;
                count += 1;
            }
        }
    }

    void path_integral(const int nwngrid,
                       const int nlayers,
                       const int nmol,
                       const double * sigma_array,
                       const int * sigma_np,
                       const int * sigma_nt,
                       const float * sigma_t,
                       const float * sigma_p,
                       const double * dlarray,
                       const double * z,
                       const double * dz,
                       const double * rho,
                       const double ** X,
                       const double * T,
                       const double * P,
                       const double Rp,
                       const double Rs,
                       void * absorptionv) {
            
        double * absorption = (double *) absorptionv;
        int count;
        cout << "start" << endl;
        double ****sigma = new double***[nmol];
    	for(int l=0;l<nmol;l++) {
            cout << "mol " << l <<  endl;
    	    sigma[l] = new double**[sigma_np[l]];
            for(int k=0;k<sigma_np[l];k++) {
                sigma[l][k] = new double*[sigma_nt[l]];
                cout << "pres " << k << " nt for this mol is " << sigma_nt[l] << endl;
                for(int j=0;j<sigma_nt[l];j++) {
//                    cout << "temp " << j <<  endl;
                    sigma[l][k][j] = new double[nwngrid];
                    for(int z=0;z<nwngrid;z++) {
//                        count = 0;
//                        for (int p=0;p<l;p++) {
//                            count += p*sigma_np[p]*sigma_nt[p]*nwngrid;
//                        }
//                        count += k * sigma_np[l] * sigma_np[l];
//                        count += j * sigma_nt[l];
                        sigma[l][k][j][z] = 3;
                    }
                }
            }
        }
        cout << "finish" << endl;

        cout << sigma[1][3][2][100] << endl;


        for (int wn=0; wn < nwngrid; wn++) {
            count = 0;
    		for (int j=0; j<(nlayers); j++) { 	// loop through atmosphere layers, z[0] to z[nlayers]
    			for (int k=1; k < (nlayers-j); k++) { // loop through each layer to sum up path length
    				for(int l=0;l<nmol;l++) {



                    }
                }
            }
        }

        cout << sigma_array[10] << endl;
        cout << "ciao" << endl;
    }
}