// C++ wrapper to the subroutine sigrad_sim.f  by Ryan Zielinski
// and David Ruth "David  <David.Ruth@unh.edu> from the originla nradcor.f
// by Karl Slifer
// Carlos Ayerbe <gayoso@jlab.org>
// - 08/2023: This version. Very simple code. The user must codify
//   the readout section according to its needs. 
//
// NOTES
// Future possible improvements:
// - code output to a text file. Right now it dumps output to screen. 
// - read an input file with the different parameters inputs by the users
// - independent readout methods created by users which feed the subroutine as
//   an independet method
// - more C++-ish style
// - ROOT is not used except for the TApplication and TStopwatch
//   both of them are not even really needed. Maybe if we want
//   to add extra functionalities as plotting.
// - when it is relevant, I added a lower letter at the beggining of the
//   name of the variable indicating the type of variable, d: double, v:vector etc.


// Need to clean headers not used.

#include "Rtypes.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"
#include <TROOT.h>
#include <TStyle.h>
#include <vector>
#include <math.h>
#include "TStopwatch.h"
#include <sstream>
#include "TSystem.h"
#include "TString.h"


using namespace std;

// fortran subroutines wrapped are declared as void
//- parameters to send (in order):
//-- E beam, angle, rad len before, rad len after, Z target, A target, nu, Born Cross Section
//- parameters to receive:
//-- cross section radiated
// * Angle in radians
// * Energy, nu, in MeV
// * Cross section in nb/MeV-sr

extern"C"{
  void radiate_(double *, double *, double *, double *,double *, double *,double *,double *, double *);
}

void  read_and_radiate();

void read_and_radiate()
{
  cout<<"read_and_radiate"<<endl;
  string line;
  
  // the different columns of the file
  // the number of items will depend of how many columns are in the file
  // I believe it is not necesary to read all the columns if it is not
  // necesary
  string it1, it2, it3, it4, it5, it6, it7, it8, it9, it10, it11;

  double dAngle, dEs;

 // Vectors to store the read values
  //  vector<double> vEs;
  vector<double> vEprime;
  vector<double> vnu; // this is calculated later
  vector<double> vW;
  vector<double> vXS;

  
  // The module/code (next lines) to read the file should be modify according to content of the file
   
  ifstream  ReadFile("muchana_data_13.4.csv");

  if (ReadFile) // it controls that the file exists
    {
      while (!ReadFile.eof()) // loop until the end of the file (eof)
	{
	  while ( getline( ReadFile, line )) // take a full line from the file
	    {
	      	  stringstream ss( line );      // Set up up a stream from this line

		  // I read 11 columns since this is the data from Murchana (rc-externals)
		  // probably, I don't need to read beyond 7, except if I want to make use
		  // of a different cross-section
		  ss >> it1 >> it2 >> it3 >> it4 >> it5 >> it6 >> it7 >> it8 >> it9 >> it10 >> it11; //

		  
		  //vEs.push_back( stod(it1)*1000 ); // in MeV
		  dEs = stod(it1)*1000; //in MeV. We make use of only one Es value (Ebeam)
		  vEprime.push_back( stod(it2)*1000 ); // in MeV
		  dAngle = (stod(it3)); // in deg (one angle at the time)	  		  
		  vW.push_back(stod(it3));  // in MeV 
		  vXS.push_back(stod(it6)/1000.); // in nb/MeV-sr (Born XS)

		  // E - E'
		  vnu.push_back( (stod(it1) - stod(it2)) *1000 ); // in MeV
		  
	    }
	}
    }


  // The original radcor code needed the nu values in increasing order
  // Murchana's data has Eprime in increasing order, so when calculating nu
  // the order of all the data will be in decreasing order of nu
  // With these lines, we reverse the order of the data.  
   reverse(vW.begin(), vW.end());
   reverse(vnu.begin(), vnu.end());
   reverse(vXS.begin(), vXS.end());

   double dAnglerad = dAngle /(180./3.141592654); //Angle in radians

   //**************************************************************************
   // in a perfect life, these values should be read from an external file
   // or edited in other part of the program. The simple way is having there here.

   // Numbers from Murchana's thesis (HMS)
   // Radiation lengths  
   // BEFORE interaction
   double BeWin = 1.33E-3;
   double N2Win0 = 3.36E-4;
   double USGE180Win = 1.946E-3;
   double He3_half = 0.445E-3;
   
   // AFTER interaction
   double GE180Win = 0.0213/sin(dAnglerad); // the thickness here is a function of the scattering angle
   double N2Win1 = 1.348E-3;
   double LexanWin = 4.5896E-3;
   double AirWin = 8.09E-4;
   double KevWin = 5.107E-4;
   double MylWin = 4.425E-4;
   //**************************************************************************   

   double dTb  = BeWin+N2Win0+USGE180Win+He3_half; // Total thickness before scattering in radiation lengths
   double dTa  = He3_half+GE180Win+N2Win1+LexanWin+AirWin+KevWin+MylWin; // Total thickness after scattering in radiation lengths

   // He3 Z and A. Declared double to make compatible with the subroutine
   double dZ = 2;
   double dA = 3; 


   vector<double> vRadXS;

   //auxiliary variables
   double dnuaux, dXSaux, dRadXSaux;

   for (int kk = 0; kk < vnu.size(); kk++)
     {
       dnuaux = vnu[kk];
       dXSaux = vXS[kk];		    
       
       radiate_(&dEs, &dAnglerad, &dTb, &dTa, &dZ, &dA , &dnuaux, &dXSaux, &dRadXSaux);
       cout<<dnuaux<<" "<<dRadXSaux<<endl;
       vRadXS.push_back(dRadXSaux);
     }
   
}








int main(int argc, char *argv[])
{
    TApplication *theApp = new TApplication("app", &argc, argv); //<======

  cout << "Compiled on: " __DATE__ " " __TIME__ "." << endl;
  
  TStopwatch timer;
  timer.Start();

  // TString inifile = "options.ini";

  read_and_radiate();
  
  gStyle->SetOptFit(1);
  gStyle->SetStatFormat("6.6g");

  timer.Stop();
  
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();

  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

  cout << "-----> CTRL+C to end" << endl;
  
  theApp->Run();
  
  return 1;
}
