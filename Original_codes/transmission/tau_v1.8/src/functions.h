/**************************** The headers ****************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <vector>
#include <omp.h>
using namespace std;



/**************************** The definitions **************************/
#define EVER ;;
#define TINY std::numeric_limits< double >::min()
#define VBIG std::numeric_limits< double >::max()
#define PI 3.14159265	// pi
#define d2s 86400		// day-to-second conversion constant
#define d2r PI/180.		// degree-to-radian conversion constant
#define r2d 180./PI		// radian-to-degree conversion constant

const double	RSOL=6.955e8;			// radius of the Sun (m)
const double	MSOL=1.9891e30;			// mass of the Sun (kg)
const double	RJUP=6.9911e7;			// radius of Jupiter (m)
const double	MJUP=1.8986e27;			// mass of Jupiter (kg)
const double	REARTH=6.371e3;			// radius of Earth (m)
const double	MEARTH=5.9736e24;		// mass of Earth (kg)
const double	AU=1.49e11;				// 1 AU (m)
const double	KBOLTZ=1.380648813e-23;	// Boltzmann's constant (J/K)
const double	AMU=1.660538921e-27;	// Atomic mass unit (kg)
const double	AVOGADRO=6.0221415e23;	// Avogadro's number
const double	RGAS=AVOGADRO*KBOLTZ;	// Universal gas constant (J/K/mol)
const double	LO=2.68676e+25;			// Loschmidt's number (m^-3)
const double	AMA=2.68676e+25;		// Amagat (molecules m^-3)



/**************************** The classes ****************************/

class Mol {
  public:
	Mol(string mol, double wt, double rad, double rdx) {
		name=mol;		// molecule name
		weight=wt;		// relative molecular weight (amu)
		radius=rad;		// molecular radius (m)
		rindx=rdx;		// refractive index
	}	// class constructor

	string name;
	double weight, radius, rindx;
};


class Atmos {
  public:
	Atmos() {mu=0.0;def_mu=2.3;}	// class constructor

	vector<Mol> mol_list;
	vector<double> fraction;
	double def_mu;	// default atmos 85% H2, 15% H2 --> mu~2.3
	double mu;

	void ADD_MOL(Mol, double);
	void GET_MMW();
	int CHECK_ATMOS();
};

void Atmos::ADD_MOL(Mol mol, double frac){
	// input mass mixing ratio as 'frac', such that e.g. if atmosphere 80% H2, frac_H2=0.8

	mol_list.push_back(mol);
	fraction.push_back(frac);
}

void Atmos::GET_MMW(){
	int nmols=mol_list.size();

	for(int i=0;i<nmols;i++) mu += (fraction[i] * mol_list[i].weight);
}

int Atmos::CHECK_ATMOS(){
	int nmols=mol_list.size();
	double tot_frac=0.0;

	for(int i=0;i<nmols;i++) tot_frac += fraction[i];

	return((tot_frac != 1.0) ? 0 : 1);
}



/**************************** The functions ****************************/

/* Usage instructions */
#ifndef instructions_H
#define instructions_H
void instructions(char argv[256])
{
	cout << "Usage: " << argv << " 0 [[atmfile]] [absfile]" << endl;
	cout << "\tor " << argv << " 1 [atmfile]" << endl;

	cout << "\n\t0,1: Only options for now. \n\t" 
	<< "atmfile: optional (if no absfile specified). If not provided, default atm file used\n\t"
	<< "absfile: optional. If not provided, default abs file used\n\t" << endl;
	exit(1);
}
#endif



/* Option sorting */
#ifndef optionSort_H
#define optionSort_H
void optionSort(int option, int argc, char* argv[], vector<string> &arg_btFile, string &arg_atmFile)
{

// What to do with entered options
	int trigger2=0,trigger3=0;

	if (argc <= 1)  	// if no arguments provided at program execution, display usage instructions
	{
		instructions(argv[0]);
	}
	else if (argc <= 3)		// prevent console giving rubbish values into arguments. Checks count of arguments
	{
		trigger2=1;
	}
	else if (argc <=4)
	{
		trigger2=1;
		trigger3=1;		// if we have 4 argv: (filename, 0, myatm, myabs) myatm and myabs MUST be present!
	}


	if (option == 0)  // what to do with arguments 0= single file read; 1= user input abs files
	{
		if (!argv[2])		// if no argument #2 given, switch to default
		{
			arg_atmFile = "./run/profile.atm";
			cout << "\nNo atm file provided in arguments, using code default...\n";
			
			
		}
		else {
			arg_atmFile = argv[2];
			cout << "You have provided atm file: " << arg_atmFile << endl;
		}
		
		if (!argv[3] || (trigger3 == 0))	// if no argument #3 given, console sends garbage as value 3 sometimes, use trigger to prevent this
		{
			arg_btFile.push_back("./run/h2o_1500K.abs");
			cout << "No abs file provided in arguments, using code default...\n";
		}
		else {
			arg_btFile.push_back(argv[3]);
			cout << "You have provided abs file: " << arg_btFile[0] << endl;
		} 

	}
	else if (option == 1)  // what to do with arguments 0=file read 1= user input abs files
	{
		if (!argv[2])		// if no argument #2 given, switch to default
		{
			arg_atmFile = "./run/profile.atm";
			cout << "\nNo atm file provided in arguments, using code default...\n";
			
			
		}
		else {
			arg_atmFile = argv[2];
			cout << "You have provided atm file: " << arg_atmFile << endl;
		}


		string input="";
		cout<< "\nEnter names of gas absorption coefficient files"<<endl;
		cout<< "\t(enter 'x' when done): " <<endl;

		for(EVER){
		   getline(cin, input);
		   if(input=="x") break;
		   else arg_btFile.push_back(input);
		}

	}
	else if (option == 9)  		// testing mode
	{
		arg_atmFile = "./run/profile.atm";

		arg_btFile.push_back("./run/h2o_1500K.abs");
	}
	else{	// END OF FILE READ
		cout << "\nPlease use option '0', '1' or '9'.\n" << endl;
		instructions(argv[0]);
	}
}
#endif



/* Get number of lines in a file, assuming header of 11 lines present at TOF */
#ifndef getNumberLines_H
#define getNumberLines_H
int getNumberLines(const char* filename)
{
	string line, line1;
	ifstream myfile (filename);	// open once to count number of lines (up to ***** line OR file end)
	int linecount=1, totalline;
	if (myfile.is_open())
	{
		while ( myfile.good() )
		{
			getline (myfile,line);
			if (line == "******") break;
			linecount++;
		}
		totalline = linecount;
		linecount=1;
	}
	else
	{
		cout << "\nUnable to open file " << filename << endl << "Exiting program...\n" << endl;
		exit(1);
	}

	myfile.close();					// close file
	
	//cout << "Total lines: " << totalline << endl;
	
	return totalline-11;

}
#endif



/* Function to read in file of absorption cross-sections */
#ifndef readAbsFile_H
#define readAbsFile_H
int readAbsFile(char *file, vector<double> &wl, vector<double> &sig)
{
	/* Read data in (abs c/s file, with wavelength in micron, sigma in cm^2) */
	ifstream the_file (file);
	vector<double> in_data;
	double	d=0.0;

	if(!the_file.is_open()){
	  cout<< "Error opening data file '" << file << "'" <<endl;
	  return(1);
	}
	else while(the_file >> d) in_data.push_back(d);	// read from file and put in in_data

	the_file.close();	// close file after read-in


	/* Re-organise data, and convert sigma units to m^2 */
	int n_cols=2;

	if(in_data[0]<in_data[2]){	//reverse order (into wavelength decreasing)
		for(int i=in_data.size()-1;i>0;i-=n_cols){
		   wl.push_back(in_data[i-1]);
		   sig.push_back(in_data[i]*1.0e-4);
		}
	} else{
		for(int i=0;i<in_data.size();i+=n_cols){
		   wl.push_back(in_data[i]);
		   sig.push_back(in_data[i+1]*1.0e-4);
		}
	}

   int btlines = getNumberLines(file);
   btlines+=11;	// due to totallines-11 in readNumberLines
   cout << "\nabs lines read: " << btlines << endl;

return(0);
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



/* Function to interpolate absorption cross-section files to same wavelength grid */
#ifndef interpolateAbs_H
#define interpolateAbs_H
int interpolateAbs(vector<char*> &files, vector<vector<double> > &sigma, double wl_min, double wl_max, float res, int &n_gas)
{
	/* Define some data vectors */
	vector<double> data_xx, data_yy;	// for input data
	vector<vector<double> > xx, yy;		// for valid input data


	/* Read in .abs files */
	for(int i=0;i<files.size();i++){
	   if(!readAbsFile(files[i],data_xx,data_yy)){		// then .abs file read succesful
	   		xx.push_back(data_xx);						// add input data to 'valid data' array s.t. each row is a different gas
	   		yy.push_back(data_yy);
	   }
	}


	/* Calculate over defined range, or largest range covered by absorption cross-sections */
	if(wl_max<wl_min){
		if(xx.size()>1){

			wl_max=max(xx[0][1],xx[1][1]);
			wl_min=min(xx[0][xx[0].size()-2],xx[1][xx[1].size()-2]);
				// get largest and smallest overall wavelength values (initial values from first file)

			for(int i=2;i<n_gas;i++){
				wl_max=max(wl_max,xx[i][1]);
				wl_min=min(wl_min,xx[i][xx[i].size()-2]);
					// NB need an extra value each end for upper/lower interpolation bounds for max/min values
			}
			// updating initial values if range is different for subsequent files

		} else{
			wl_max=xx[0][0];
			wl_min=xx[0][xx[0].size()-1];
		}

	} else if(wl_max == wl_min){
			cout<<"Zero range! Exiting....."<<endl;
			exit(1);
	}
	cout<<"\tfrom "<<wl_min<<" to "<<wl_max<<" microns"<<endl;


	if(n_gas != xx.size()) cout<<"New n_gas = "<<xx.size()<<endl;
	if(!(n_gas=xx.size())) cout<<"No molecules entered! "<<endl;

	const double startTime = omp_get_wtime();


	/* Layout wavelength grid at even intervals */
	int i=0;
	for(i=0;(wl_max-(i*res)) > wl_min;i++) sigma[0].push_back(wl_max - (i*res));
		// top row of sigma 2d array is for wavelengths,
		// and now i=number of lines=wl.size()=sigma[0][*].size()


	/* Interpolate from files */
	for(int n=0;n<n_gas;n++){			// loop through gases

		/* Initialise gas abs coeff slots for gas n */
		for(int i=0;i<sigma[0].size();i++) sigma[n+1].push_back(0.0);


		/* Create a parallel region */
		//#pragma omp parallel num_threads(1)			// specify num_threads...
		#pragma omp parallel							// ...or use default num_threads
		{
			const int thread_id=omp_get_thread_num();	// get thread id on first pass
			if(n==0){
				// #pragma omp single
				// 	cout<<endl<<"My name is Legion, for we are "<< omp_get_num_threads() <<endl<<endl;
			}

			double bounds[3]={0.0};


			/* Start parallel loop */
			#pragma omp for schedule(static) nowait
			for(int j=0;j<sigma[0].size();j++){
			/* for every (new) wavelength, find equivalent location in wl grid of original file by going down original file and checking if new value is between orig_wl[k] and orig_wl[k+1] */

				/* Get interpolation bounds */
				for(int k=0;k<xx[n].size()-1;k++){		// NB need an extra value each end for upper/lower interpolation bounds for max/min values

					if(sigma[0][j]==xx[n][k]) sigma[n+1][j] = yy[n][k];
						// no interpolation needed - e.g. endpoints of smallest input file

					if((sigma[0][j]<xx[n][k]) && (sigma[0][j]>xx[n][k+1])){
						/* NB TAKE CARE WITH EQUALITY SIGNS - wl vector is in DECREASING order, so need
							val < orig_wl[k] and
							val > orig_wl[k+1] */

						bounds[0]=xx[n][k+1];
							// TAKE CARE WITH VECTOR ORDER AGAIN
						bounds[1]=xx[n][k];
						bounds[2]=sigma[0][j];

				/* And actually assign value */
						sigma[n+1][j] = interpolateValue(bounds,yy[n][k+1],yy[n][k]);
							// TAKE CARE WITH VECTOR ORDER AGAIN
						break;
					}
				}	// end of interpolation for wavelength lambda_j

			}	// end of (parallel) loop over wavelengths
			
		} 	// end of parallel region

	}	// end of loop over gases



/* Debug - check interpolation *//*
const double endTime = omp_get_wtime();
const double totalTime = endTime - startTime;

cout<<"Interpolation time: "<< totalTime <<" seconds"<<endl;
/**/

/* Debug - check output *//*
for(int n=0;n<n_gas;n++){

cout<<"Gas "<<n<<"\n============================================="<<endl;
for(int j=0;j<sigma[0].size();j++){ 
cout<<setprecision(10)<<j<<"\t\t"<<sigma[0][j]<<"\t\t"<<sigma[n+1][j]<<endl;
}
cout<<endl;

}/**/

return(0);
}
#endif





/* Function to interpolate single files to same wavelength grid */
#ifndef interpolateCS_H
#define interpolateCS_H
int interpolateCS(const char *in_file, vector<double> &wl, vector<double> &cs)
{
/*
	VARIABLES:	wl[] = base wavelength grid
				cs[] = interpolated values at wl[]
				xx[] = original file wavelength grid
				yy[] = original file data values
*/


	/* Read data in */
	ifstream the_file (in_file);
	vector<double> in_data;
	int 	n_cols=2;
	double	d=0.0;
	vector<double>	xx,yy;		// storage vectors for input cs data

	if(!the_file.is_open()){
	  cout<<endl<< "Error opening file '" << in_file << "'" <<endl;
	  return(1);
	}
	else while(the_file >> d) in_data.push_back(d);	//read from file and put in in_data

	the_file.close();	//close file after read-in

	/* Put input cs data into data vectors */
	for(int i=0;i<in_data.size()-1;i+=n_cols){
		xx.push_back(in_data[i]);		// wavelength
		yy.push_back(in_data[i+1]);		// data value
	}


	/* Interpolate data to base wavelength grid, wl[] */

		/* Create a parallel region */
		#pragma omp parallel
		{
			const int thread_id = omp_get_thread_num();		// get thread id

			double bounds[3]={0.0};


			/* Start parallel loop */
			#pragma omp for schedule(static) nowait
			for(int j=0;j<cs.size();j++){
			/* for every (new) wavelength, find equivalent location in wl grid of original file by going down original file and checking if new value is between orig_wl[k] and orig_wl[k+1] */

				/* Get interpolation bounds */
				for(int k=0;k<xx.size()-1;k++)	// NB need an extra value each end for upper/lower interpolation bounds for max/min values
				{
					if(wl[j]==xx[k]) cs[j] = yy[k];
						//no interpolation needed (values match)

					else if((wl[j]<xx[k]) && (wl[j]>xx[k+1])){
						/* NB TAKE CARE WITH EQUALITY SIGNS - wl vector is in DECREASING order, so need
							val < orig_wl[k] and
							val > orig_wl[k+1] */

						bounds[0]=xx[k+1];
							//TAKE CARE WITH VECTOR ORDER AGAIN
						bounds[1]=xx[k];
						bounds[2]=wl[j];

				/* And actually assign value */
						cs[j] = interpolateValue(bounds,yy[k+1],yy[k]);
						break;
					}
				}	// end of interpolation for wavelength lambda_j

			}	// end of (parallel) loop over wavelengths

		} 	// end of parallel region


	if( (wl.front()>xx.front()) || (wl.back()<xx.back()) ){		//TAKE CARE WITH VECTOR ORDER AGAIN
		cout<< "WARNING: range doesn't match model wavelength grid for file " << in_file <<endl;
	}


/* Debug - check output *//*
for(int j=0;j<cs.size();j++){ 
	cout<<setprecision(6)<<j<<"\t\t"<<wl[j]<<"\t\t"<<cs[j]<<endl;
}
cout<<endl;
/**/


return(0);
}
#endif





/* Get number of gases in .atm file */
#ifndef getNumberGases_H
#define getNumberGases_H
int getNumberGases(char* filename)
{
	string line, line1;
	ifstream myfile (filename);	// open once to count number of layers (up to ***** line OR file end)
	int linecount=1, totalline;
	if (myfile.is_open())
	{
		while ( myfile.good() )
		{
			getline (myfile,line);
			if (line == "******") break;
			linecount++;
		}
		totalline = linecount;
		linecount=1;
	}
	else
	{
		cout << "\nUnable to open file " << filename << endl << "Exiting program...\n" << endl;
		exit(1);
	}

	myfile.close();					// close file
	
	//cout << "Total lines: " << totalline << endl;
	
	return totalline-11;

}
#endif



/* Function to read in .atm file */
#ifndef readAtmFile_H
#define readAtmFile_H
void readAtmFile(char* filename, int numlines, float* arrayP, float* arrayT, float* arrayZ, vector<vector<float> > &arrayX)
{
	
	string line,line1;
	int linecount = 1, totalline=numlines+11, colcount=0, chicount;

	float input1, input2, input3;

	ifstream myfile2 (filename);		// re-open for data read.
	if (myfile2.is_open())
	{

		
		while ( myfile2.good() )
		{
			if ((linecount > 10) && (linecount < totalline))
			{
				/* Count number of columns in first line of data */
				if(colcount==0){
					string buf;
					stringstream ss(line1);
					vector<string> tokens;

					while(ss >> buf) tokens.push_back(buf);

					colcount=tokens.size();
					chicount=colcount -3;

					if(colcount>3) arrayX.resize(chicount);		// resize to height=n_chi_cols
					else{
						cout<< "\nWARNING: .atm file must have columns of 'p, T, z, X1 [, X2, ...]'!" <<endl;
						exit(1);
					}
				}


				myfile2 >> input1;
				myfile2 >> input2;
				myfile2 >> input3;

				arrayP[totalline-(linecount+1)] = input1;
				arrayT[totalline-(linecount+1)] = input2;
				arrayZ[totalline-(linecount+1)] = input3*1000;	// convert from km to m


				/* Get mixing ratios in remaining columns */
				float input4[chicount];
				for(int i=0;i<chicount;i++){
					myfile2 >> input4[i];
					arrayX[i].push_back(input4[i]);
				}
				
			}
			else {
				//cout << "Ignored line: " << line1 << endl;
			}
			getline (myfile2,line1);

				linecount++;
		}
		myfile2.close();
	}
	
	
}
#endif





/* Function to calculate H2 Rayleigh scattering cross-section */
#ifndef scatterRayleigh_H
#define scatterRayleigh_H
double scatterRayleigh(double lambda, Mol species)
{
/* Formula from Liou 2002, 'An Introduction to Atmospheric Radiation', pp.92-93. Also uses 'minimum volume' approximation pg.97, 
		N_dens = 1 / V_particle .

	Optical depth given by tau = sigma * L * c ; sigma = abs cross-section (m^2), L = path length (m), c = concentration (m^-3)
	
	NB This is for bulk atmos scattering ONLY (assumptions: particles much smaller than wavelength, gas sufficiently dense), 
	cloud Rayleigh + Mie included in scatterMie function. 

	IN: Wavelength (in um), path length (in m)

	OUT: Rayleigh scattering opacity cross-section per particle (in m^2)
*/


   double sigma_R=0.0;			// Rayleigh absorption coefficient (from Liou, An Introduction to Atmospheric Radiation)

   double wl=lambda *1.0e-6;		// convert wavelengths to m

   double rad=species.radius;		// molecular radius (m)

   double r_ind=species.rindx;		// molecular refractive index
   double r_sq=r_ind*r_ind;
   double r_red = (r_sq-1) / (r_sq+2);

   double delta = 0.035;					// molecular anisotropy factor
   double f_delta = (6.0+(3.0*delta)) / (6.0-(7.0*delta));	// King correction factor
   	if(species.name == "He") f_delta = 1.0;		// no asymmetry for helium molecules

   /* Find cross-section */
   sigma_R = (128.0/3.0) * (pow(PI,5) * pow(rad,6) / pow(wl,4)) *r_red*r_red *f_delta;		// gives sigma_R in m^2


return(sigma_R);
}
#endif





/* Function to convert H2-H2 CIA coefficients from A. Borysow data into cross-sections */
#ifndef scatterCIA_H
#define scatterCIA_H
double scatterCIA(double coeff, double amount)
{
/* Optical depth given by tau = alpha * L * c_1 * c_2 ; alpha = abs coeff (cm^5 mol^-2), L = path length (cm), c_i = concentration of collider i (mol cm^-3)

	IN: CIA coeffs in (cm^-1 amagat^-2), grid wavelength (in um), path length (in m) and total number density dz (in m^-3)

	OUT: H2-H2 collision-induced absorption coefficient (in m^5 mol^-2)
*/


	/* Calculate unit conversion factor from (cm^-1 amagat^-2) to (cm^5 mol^-2), i.e. into HITRAN cia format... */
	double conv_factor = 1.0 / pow((AMA*1.0e-6),2);
		// conversion factor from absorption coefficient alpha (cm^-1 amagat^-2) to (cm^5 mol^-2)
		//	 = 1/(AMA^2), with AMA in mol cm^-3

	//double conv_factor = 1.0;


	/* ...and calculate cross-section */
	double alpha = coeff * conv_factor;			// converting from (cm^-1 amagat^-2) to (cm^5 mol^-2)...

	alpha *= (amount*amount) * 1.0e-10;			// e.g. composition 85% H2, and convert from cm^5 to m^5

return(alpha);
}
#endif
