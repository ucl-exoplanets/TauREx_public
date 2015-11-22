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

        cout << "ciao" << endl;
        
    }
}