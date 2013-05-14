/************************************************************************************
** TAU.CPP - Marcell Tessenyi 2011 - v0.1
**			- Morgan Hollis 2012 - v1.8c
**
**
** This code is a 1D radiative transfer code for transmission spectroscopy of
**	extrasolar planets. It uses a line-by-line integration scheme to model 
**	transmission of the radiation from a parent star through the atmosphere of an
**	orbiting planet, in order to compare to observations of the radius ratio as a
**	function of wavelength in primary transit, and hence to infer the abundances of 
**	trace absorbers present in the planetary atmosphere. 
**
** The code reads in an atmospheric profile and absorption cross-sections (filenames
**	input by user on prompt) for the required absorbers and calculates the optical 
**	path length exp(-tau) in the transit geometry, outputting the transit depth 
**	(radius ratio) as a function of wavelength.
**
** Run './tau' to display usage instructions and run modes.
**
**
** INPUTS:  - atmospheric temperature-pressure profile file, e.g. "profile.atm"
**				10-line header, then 39 atmosphere levels,
**				column 1: Pressure, in Pascal
**				column 2: Temperature, in Kelvin
**				column 3: Altitude, in kilometres
**
**			- absorption cross-section files, e.g. "molecule1.abs"
**				one per molecule/absorber, containing,
**				column 1: wavelength (in microns)
**				column 2: absorption cross-section (in cm^2)
**
**			- [OPTIONAL] stellar radius as a function of wavelength, e.g. "rad_star.rad"
**				with,
**				column 1: wavelength (in microns)
**				column 2: stellar radius (in m)
**
**				N.B. if this option is used (and sw_rad set), the user must ensure that
**					the designated file wavelength values span the entire specified 
**					wavelength range for the model, in order for the interpolation to 
**					function correctly. 
**
**			- [OPTIONAL] H2-H2 Collision-Induced Absorption coefficient file, 
**				e.g. "h2_h2_1000K.cia", with,
**				column1: wavelength (in microns)
**				column2: absorption coefficient (in cm^-1 amagat^-2)
**
**			- [OPTIONAL] opacities from cloud models if available, e.g. "cloud1.cld"
**				formatted as,
**				column 1: wavelength (in microns)
**				column 2: mass opacities (in cm^2 g^-1)
**
**
** OUTPUT: 	- "tau_output.dat" - file containing wavelengths (column 1, in microns) and 
**								absorption/radius ratio (column 2, dimensionless).
**
**
** For more details, see M. D. J. Hollis et al., "TAU: A 1D radiative transfer code
**	for transmission spectroscopy of extrasolar planet atmospheres". Comp. Phys. Comm.
**	(2013).
**
** Email mdjh@star.ucl.ac.uk with any questions or bugs, and I'll try to help!
**
*************************************************************************************/
#include "functions.h"

int main(int argc, char *argv[])
{
/***************************** Define run variables ********************************/
const double prog_start=omp_get_wtime();


   /* Run parameters */
   const float rad_fac=-0.0;		// increase planet radius by (rad_fac)% 
									//	i.e. rad_fac > 0 for input Rp lowers radius @ 1bar

   const float s_rad_fac=+0.0;		// increase stellar radius by (s_rad_fac)%

   double lambda_min=0.00;		// define wavelength range for model (in microns)
   double lambda_max=20.00;

   float lambda_res=0.01;			// resolution of new wavelength grid (o/p spectrum) in microns
   //	if(lambda_max>lambda_min) lambda_res = (lambda_max-lambda_min) / 2000.0;

   const float mixdef = 1.0e-5;			// default mixing ratio value


   /* System parameters - e.g. for HD189733b..... */
   const float Rp =  1.138 *(1. - (rad_fac/100.)) *RJUP; // planet radius @ 1bar level: *Jupiter radius (m)
   const float Rstar =  0.788 *(1. - (s_rad_fac/100.)) *RSOL; // stellar radius: *Sun radius (m)
   const float semimajor = 0.03142 *AU; // semi-major axis: *1 AU (m)

   const float grav=23.45;			// gravitational acceleration at planetary surface (m s^-2)
   const float temp=1500;			// atmospheric temperature (K)


   /* Atmosphere parameters */
   Atmos atmos;

   const Mol H2("H2",2.0,2.0e-9,1.0001384);		// define possible bulk atmosphere constituents
   const Mol He("He",4.0,1.0e-9,1.0000350);

   atmos.ADD_MOL(H2,0.85);		// add molecules to atmosphere, with corresponding mass mixing ratios (fractional abundances), such that sum = 1
   atmos.ADD_MOL(He,0.15);

   if( atmos.CHECK_ATMOS() ) atmos.GET_MMW();	// check composition adds up to 100% and get mean relative molecular weight of atmosphere
   else {
	cout<< "WARNING: Bulk atmosphere composition doesn't add up!" <<endl;
	exit(1);
   }
   		// ...to calculate atmospheric scale height (km)
   const float H = (RGAS * temp) / (atmos.mu * grav);


   /* Filenames and switches for external file inputs (0=off, 1=on) */	
   const int sw_rad=0;		// Vary stellar radius with wavelength
   const int sw_cia=0;		// Include H2-H2 CIA
   const int sw_cld=0;		// Read in extra optical depths due to clouds, set switch to equal number of files to be read in


   const char* rad_file={"./run/rad_star.rad"};		// file from which to read stellar radius R*(lambda)

   const char* cia_file={"./run/h2_h2_1500K.cia"};	// file from which to read CIA coefficients

   const char* cld_file[]={"./run/cloud1.cld"};		// files from which to read optical depths for extra opacities (e.g. clouds)


   string arg_outFile = "./out/tau_output.dat";		// file to contain final spectrum



/***************************** Sorting options ********************************/
   string arg_atmFile;
   vector<string> arg_btFile;

   char *atmFile, *outFile;
   vector<char*> btFile;
	// multiple abs file inputs possible --> vector CONTAINING some number of pointers to chars
	//vector<char> *btFile;		whereas this would be one pointer TO a vector containing chars


   int option=0;

   if (argc > 1) option = atoi (argv[1]); // arg-to-int: convert character from 1st argument to integer
   else {
    cout << "\nPlease use option '0', '1' or '9'.\n" << endl;
    instructions(argv[0]); // if no arguments provided at program execution, display usage instructions
   }

   optionSort(option,argc,argv,arg_btFile,arg_atmFile);


   /* Assign and check some run parameters */
   int n_gas=arg_btFile.size();
   	if( !n_gas ) cout<<"No molecules entered! "<<endl;

   const int n_cld=sw_cld;

   cout << endl << "Files used:\n Atm: " <<"\t"<< arg_atmFile <<"\n"<< " Abs: ";
   for(int i=0;i<arg_btFile.size();i++) cout<<"\t"<< arg_btFile[i] <<endl;

   if(sw_rad) cout<< " R*: "<<"\t"<< rad_file <<endl;
   if(sw_cia) cout<< " CIA: "<<"\t"<< cia_file <<endl;
   if(sw_cld){
   		cout<< " Cld: ";
   		for(int i=0;i<n_cld;i++) cout<<"\t"<< cld_file[i] <<endl;
   }
   
   cout<< " n_gas = "<<n_gas <<endl;
   cout<<endl<< " O/P: " <<"\t"<< arg_outFile << endl;

   float thres = 50.0 * (atmos.mol_list[0].radius * 1.0e6);		// NB converting particle radius to microns
   if(lambda_min<thres){
	cout<< "\nWARNING: Rayleigh scatter not calculated for wavelengths below "<<thres<<" microns!" <<endl;
   }

   cout<<endl;


   /* Convert string type to char array (to pass to functions) */
   atmFile=new char[arg_atmFile.size()+1];
   atmFile[arg_atmFile.size()]=0;
   memcpy(atmFile,arg_atmFile.c_str(),arg_atmFile.size());

   for(int i=0;i<arg_btFile.size();i++){
 		btFile.push_back(new char[arg_btFile[i].size()+1]);
		btFile[i][arg_btFile[i].size()]=0;
		memcpy(btFile[i],arg_btFile[i].c_str(),arg_btFile[i].size());
   }

   outFile=new char[arg_outFile.size()+1];
   outFile[arg_outFile.size()]=0;
   memcpy(outFile,arg_outFile.c_str(),arg_outFile.size());





/***************************** Get data from .abs file(s) ********************************/

   /* Create array of data vectors for wavelengths and abs coeffs for 'n_gas' gases */
   vector<vector<double> > sigma_array;

   sigma_array.resize(n_gas+1);		// resize to height=n_gas+2

   /* Interpolate cross-sections to same wavelength grid */
   if(lambda_min<TINY) lambda_min=TINY;		// avoid potential zero division errors in scattering functions
   if(lambda_max>VBIG) lambda_max=VBIG;
   interpolateAbs(btFile,sigma_array,lambda_min,lambda_max,lambda_res,n_gas);


   vector<double> &gridwl=sigma_array[0];	// i.e. top row for wavelengths, and each middle row is a different gas
   int linecount = gridwl.size();
   cout<<"\nNew linecount: "<<linecount<<endl;
	// even though this is now number of columns in sigma_array



   /* Interpolate CIA coefficients etc. to the model wavelength grid */
   vector<double> rad_star, cia_coeffs;		// vectors for stellar radius, CIA coefficients as a function of wavelength

   vector< vector<double> > cld_coeffs;		// vector for cloud optical depths, one row for each cloud file
   cld_coeffs.resize(n_cld);



   for(int i=0;i<linecount;i++) rad_star.push_back(Rstar);		// stellar radius constant with wavelength if not read in from file
	if(sw_rad) interpolateCS(rad_file,gridwl,rad_star);


   for(int i=0;i<linecount;i++) cia_coeffs.push_back(0.0);		// CIA has no effect if no file input
	if(sw_cia) interpolateCS(cia_file,gridwl,cia_coeffs);


   for(int n=0;n<n_cld;n++){									// clouds have no effect if no file input
		for(int i=0;i<linecount;i++) cld_coeffs[n].push_back(0.0);
   }
   double low_p_bound[n_cld], up_p_bound[n_cld];

	if(sw_cld){		
		for(int n=0;n<n_cld;n++){
			low_p_bound[n]=1.0e-3;		// pressure (in bar) of lower pressure/upper altitude cloud bound
			up_p_bound[n]=0.1e0;		// pressure (in bar) of upper pressure/lower altitude cloud bound
				// set cloud vertical extent

			interpolateCS(cld_file[n],gridwl,cld_coeffs[n]);
		}
	}




/***************************** Get data from .atm file ********************************/

   int nlayers = getNumberLines(atmFile);		// number of usable lines from atm file
   cout << endl << "Number of layers from file: " << nlayers << endl;


   /* For each level, get..... */
   float p[nlayers];			// pressure (in Pascal)
   float Tp[nlayers];			// temperature (in Kelvin)
   float z[nlayers];			// altitude (in kilometres)

   vector<vector<float> > X;		// mixing ratios

   float rho[nlayers], rho_prime[nlayers];
   float tau[nlayers], exptau[nlayers];


   readAtmFile(atmFile,nlayers,p,Tp,z,X);	// read the file and send reference of arrays (p,Tp,z,X) which will have contents replaced
   //cout << "Values obtained from file " << atmFile << ":" <<endl << endl;


   /* Set default mixing ratios for gases */
   if(option==9){			// 'testing' mode
	for(int n=0;n<n_gas;n++){
		for(int m=0;m<nlayers;m++) X[n][m] = mixdef;
	}
   } else {
	if(n_gas != X.size()){
		cout<< "\nEXITING: number of .abs file don't match mixing ratio columns in .atm file!" <<endl;
		exit(1);
	}
   }


   /* Calculate number density for each layer and display atm file readout */
   float rho_tot=0.0;

   for (int layer = 0; layer < nlayers; layer++){

		rho[layer] = (p[layer])/(KBOLTZ*Tp[layer]);		// convert p/T to number density, in m^-3
		
		rho_tot += rho[layer];							// to get total number density along a vertical path (dz)
   }

   cout<< "Number density at surface: " << rho[0] << " m^-3" <<endl;
   cout<< "Total number density (dz): " << rho_tot << " m^-3" <<endl;




/***************************** Calculate path length integral ********************************/

	cout<<endl<<"=========================================="<<endl;
	cout<<"Performing calculation....."<<endl;

	cout << endl << "nlayers: " << nlayers;
	cout << endl << "Rp: " << (Rp/1000.0) << " km\t\tz[nlayers] (Atm) : " << (z[nlayers-1]/1000.0) << " km" << endl;
	cout << "Rp+Atm: " << ((Rp+z[nlayers-1])/1000.0) << " km" << endl;
	cout << "Scale height: H = " << H << " km" << endl;
	cout << "MMW: mu = " << atmos.mu << " g/mol" << endl << endl;


	float dl, Rsig, Rtau, Csig, Ctau, cld_tau;
		// initialise optical quantities
	double p_bar=0.0, bounds[3]={0.0}, cld_log_rho=0.0, absorption[linecount];
		// initialise cloud parameters and absorption variables


	for (int wl=0;wl < linecount; wl++)		// loop through wavelengths
	{
		/* Calculate scattering cross-sections (wavelength dependence) */
		Rsig = 0.0;			// Rayleigh cross-section
		Csig = 0.0;			// CIA cross-section

		if(gridwl[wl]>thres){
			for(int i=0;i<(atmos.mol_list).size();i++) {
				Rsig += ( atmos.fraction[i] * scatterRayleigh(gridwl[wl],atmos.mol_list[i]) );			// Rayleigh cross-section
			}
		}

		Csig += scatterCIA(cia_coeffs[wl],atmos.fraction[0]);


		/* Calculate optical path length */
		for (int j=0; j<(nlayers-1); j++) 	// loop through atmosphere layers, z[0] to z[nlayers]
		{
		/* Calculate layer lengths, and get optical path */
			dl = 0.0;					// element of path length
			Rtau = 0.0; 					// sum of Rayleigh optical depth
			Ctau = 0.0; 					// sum of CIA optical depth
			cld_tau = 0.0;					// sum of cloud optical depth
			tau[j] = 0.0;					// total optical depth


			for (int k=1; k < (nlayers); k++) // loop through each layer to sum up path length
			{
				dl = 2.0 * (sqrt(pow((Rp + z[k+j]),2) - pow((Rp + z[j]),2)) - sqrt(pow((Rp + z[k-1+j]),2) - pow((Rp + z[j]),2)));
					// Calculate half-path length, and double (from system geometry) to get full path distance


				/* Sum up taus for all gases for this path, recall sigma_array[0][*] = wavelengths */
				for(int l=0;l<n_gas;l++) tau[j] += (sigma_array[l+1][wl] * X[l][k+j] * rho[k+j] * dl);


				/* Calculate bulk atmos Rayleigh contribution (wavelength, density, layer length dependence) for this element of path */
				Rtau += Rsig * rho[k+j] * dl;


				/* Calculate CIA contribution (wavelength, density, layer length dependence) for this element of path */
				Ctau += Csig * rho[k+j] * rho[k+j] * dl;

				/* Calculate cloud contribution (wavelength, layer length dependence) for this element of path */
				p_bar = p[k+j] * 1.0e-5;								// convert pressure from Pa to bar

				for(int n=0;n<n_cld;n++){
					if( (p_bar<up_p_bound[n]) && (p_bar>low_p_bound[n]) ){	// then cloud exists in this layer [k+j]

						bounds[0]=log(low_p_bound[n]);
						bounds[1]=log(up_p_bound[n]);
						bounds[2]=log(p_bar);

						cld_log_rho = interpolateValue(bounds,-6,-1);
							// = log(cloud density), assuming linear decrease with decreasing log pressure
							// following Ackerman & Marley (2001), Fig. 6

						cld_tau += ( cld_coeffs[n][wl] * (dl*1.0e2) * (exp(cld_log_rho)*1.0e-6) );	// convert path lenth from m to cm, and density from g m^-3 to g cm^-3

						//cout<<"Pressure = "<<p[k+j]*1.0e-5<<" bars = "<<p[k+j]<<" Pa at level "<<k<<endl;
						//cout<<"Path "<<j<<", section "<<k+j<<": kappa="<<cld_coeffs[0][wl]<<" cm^2 g^-1, dl="<<dl<<" m, rho="<<exp(cld_log_rho)<<"g m^-3"<<endl;
						//cout<< "\twl: "<<gridwl[wl]<<"\tCloud tau: " << cld_tau <<endl<<endl;
					}
				}

			}


			/* Include extra opacities (scattering etc.) */
			tau[j] += Rtau;
			tau[j] += Ctau;
			tau[j] += cld_tau;

			exptau[j] = exp(-tau[j]);
		}


		/* Calculate area of circles of atmos (mediated by e^-tau), and sum */
		double integral=0.0;
		double dz[nlayers];
		for(int j=0; j<(nlayers-1); j++) dz[j] = z[j+1] - z[j];
		for(int j=(nlayers-1); j<nlayers; j++) dz[j] = dz[j-1];
		for(int j=0; j<nlayers; j++) integral += ((Rp+z[j])*(1-exptau[j])*dz[j]);
		integral*=2.0;

		absorption[wl] = ((Rp*Rp) + integral) / (rad_star[wl]*rad_star[wl]);

	}


	/* Output to file */
	ofstream myfile (outFile);
	if (myfile.is_open()){
		for(int wl=0;wl < linecount; wl++){
			myfile << gridwl[wl] << "\t " << absorption[wl] << endl;
			//cout << gridwl[wl] << "\t " << absorption[wl] << endl;
		}
		myfile.close();

	} else {
		cout << "Unable to open file" << endl;
		exit(1);
	}
	cout << "Complete.\nData in " << outFile << ", in 2 columns:\n\n\tWL (microns) \tAbsorption\n\n";


const double prog_end=omp_get_wtime();
cout<<"Total runtime: "<< prog_end-prog_start <<endl;

return(0);

}
