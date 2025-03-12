/************************************************************************
 *
 *  Filename: Scan.cpp
 *
 *  Description:
 *    Scan code to use with CoMPASS
 *
 *	Author(s):
 *     Michael T. Febbraro
 *
 *  Creation Date: 9/25/2016
 *  Last modified: 10/10/2020
 *
 *  To compile: g++ -O3 -pedantic -o Scan.exe `root-config --cflags --libs` -lSpectrum NewScan.cpp
 *      - if errors copiling in Mac OSX
 *        - remove -03 option
 *        - clang++ instead of g++
 *        - $(root-config --cflags --libs) instead of 'root-config --cflags --libs'
 *
 *
 *  If "error while loading shared libraries: libcore.so: ..." occurs, type
 *  "source `root-config --prefix`/bin/thisroot.sh"
 *
 * -----------------------------------------------------
 * 	Nuclear Astrophysics, Physics Division
 *  Oak Ridge National Laboratory
 *
 */
// C++ classes
#include <iostream>
#include <fstream>
#include <string>
#include <signal.h>
#include <iomanip>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
//#include <vector>
// custom
#include "PulseAnalysis.h"
// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"

using namespace std;

typedef struct
{
  Float_t l;               // Long integral
  Float_t s;               // Short integral
  Float_t amp;             // Amplitude
  Float_t cfd;             // Trigger time
  Float_t psd;             // PSD parameter s/l
  Float_t trg;             // Detector trigger

} Detector;
// *********
// functions
// ********* 
int Scan(char* prefix);
// *************
// main function
// *************
int main(int argc, char** argv) 
{
  if(argc<2)                                                                    
  {                                                                             
    cout << "You must provide the file number as an arguement in the form of run_X.\n";                
    cout << "This is for Mike a.k.a. The Godfather of alpha,n\n";
    return 1;                                                                   
  }                                                                             
                                                                                
  char rootname[500];                                                           
  std::string rname;                                                                 
  rname = argv[1];                                                              
  strcpy(rootname,rname.c_str());;
  return Scan(rootname);
}

int Scan(char* prefix)
{
  /** ----------------------------------------------------
   *	Variable declairation
   *   ----------------------------------------------------
   */
  Detector det;

  bool unfolding = 1;

  float	X, offset, sigma, mu;
  int	multi;

  ifstream fp[17];

  string line, fileheader;

  int i,j,k, pposition,
    Tracelength,
    eventlength;

  float pulse[2000],
    CMAtrace[2000],
    SG_pulse[2000],
    SGderv2_pulse[2000],
    baseline[2000];

  Float_t amplitude,
    risetime,
    falltime,
    width,
    CFD,
    tac,
    paraL,
    paraS,
    runtime,
    steerer_X, steerer_Y,
    temp;
	
  float energy;
  
  double timestamp;

  // For SG filtered pulse
  Float_t trace_min, trace_max;
  int trace_min_loc, trace_max_loc;
  bool zero_cross;

  char 	filename[500],
        prompt[10],
        openfile[500],
        /*prefix[500],*/
        buffer[50],
        interrputPrompt;

  Float_t trgtime, prevtime, difftime;
  Float_t prevtrgtime[10];
  long	TEvt = 0;

  uint64_t buffer64;
  uint32_t buffer32;
  uint16_t buffer16;
  
  TRandom3 r;

  TF1 *f1 = new TF1("f1","gaus",0,300);

  /** ----------------------------------------------------
   *	Calibrations and threshold
   *   ----------------------------------------------------
   */

  float cal[16] =
    {   0.0236, /*0*/
	      0.0268, /*1*/
	      0.0243, /*2*/
	      0.0301, /*3*/
	      0.0253, /*4*/
  	    0.0249, /*5*/
  	    0.0265, /*6*/
  	    1.0000, /*7*/
  	    1.0000, /*8*/
  	    1.0000, /*9*/
  	    0.0338, /*10*/
  	    0.0193, /*11*/
        0.1286, /*12*/
        0.0390, /*13*/
        0.0196, /*14*/
        1.0000, /*15*/
    }; // calibration (keVee / (bit * sample)) from manual check on calibration
  //for(int i=0; i<16; i++) cal[i] = 1.;

  float threshold[16] =
    { 15700,
	    15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700,
		  15700
    };

  /** ----------------------------------------------------
   *	Get functions
   *   ----------------------------------------------------
   */

  PulseAnalysis *Analysis = new PulseAnalysis();


  /** ----------------------------------------------------
   *	Program start...
   *   ----------------------------------------------------
   */

  cout << " ------------------------------------------------ " << endl;
  cout << " | Scan.cpp - CoMPASS binary version             |" << endl;
  cout << " |   Experiment: ND 2020                         |" << endl;
  cout << " |   Date: November 2020                         |" << endl;
  cout << " |   Calibration used: none                      |" << endl;
  cout << " |   ORNL Nuclear Astrophysics                   |" << endl;
  cout << " ------------------------------------------------ " << endl;

  //cout << "Run binary file prefix ('run_#'): ";
  //cin >> prefix;
  //char* localPointer  = new char[500];
  //localPointer = argv[1];
  //prefix = localPointer;
  TH1F *trace0 = new TH1F("trace0","Trace for channel 0",200,0,199);
  //tt->Branch("trace0","TH1F", &trace0);
  TH1F *trace1 = new TH1F("trace1","Trace for channel 1",200,0,199);
  //tt->Branch("trace1","TH1F", &trace1);

  TH1F *traceCMA = new TH1F("traceCMA","Trace for CMA",200,0,199);
  //tt->Branch("traceCMA","TH1F", &traceCMA);

  TH1F *trace0C = new TH1F("trace0C","Corrected trace for channel 0",200,0,199);
  //tt->Branch("trace0C","TH1F", &trace0C);

  int numSteps = 20;

  // Open files
  for (i = 0; i < 17; i++)
  {
    //sprintf(openfile, "/media/shahina/New Volume/data/%s/RAW/DataR_CH%d@V1725_190_%s.bin",prefix,i,prefix);
    
    sprintf(openfile, "/Volumes/New Volume/Experiments_data/2526Mgan/2526Mgan_2021/ND_2021/%s/RAW/DataR_CH%d@V1725_190_%s.bin",prefix,i,prefix);
    //sprintf(openfile, "/home/shahina/Documents/PhD/25Mg_codes_Shahina/Compass_data/%s/RAW/DataR_CH%d@V1725_190_%s.bin",prefix,i,prefix);
    
    
    
    
	  //if (i == 16) {sprintf(openfile, "../%s/RAW/Data_CH0@DT5730S_10803_%s.bin",prefix,prefix);}
    fp[i].open(openfile, std::ifstream::in | std::ifstream::binary);
    if(fp[i].is_open()) {cout << openfile << " - Open!" << endl;}
    else{ cout << "ERROR: "   << openfile << " - not open!" << endl; }
  }
  
  // Create root file
  sprintf(filename, "Extracted_rootfiles/%s.root",prefix);
  TFile *ff = new TFile(filename, "RECREATE");
		
   // Create an array of trees
   TTree *tt[17];
   for (i = 0; i < 17; i++)
   {
		sprintf(buffer, "T%d", i);
		tt[i] = new TTree(buffer, fileheader.c_str());
		tt[i]->Branch("d",&det,"l:s:amp:cfd:psd:trg"); //l:lgate ,s:sgate,amp:pulse height(max), cfd, 
		tt[i]->Branch("E",&energy,"E/f");
		tt[i]->Branch("time",&timestamp,"time/d");
   }

  // Process channels
  for (j = 0; j < 17; j++)
	{

	  if(fp[j].is_open())
	    {
        TEvt = 0;

        // Binary parsing (DPP-PSD firmware)
        while (fp[j].read((char*)&buffer16, 2)) // Board number
        {
          fp[j].read((char*)&buffer16, 2);  // Channel number
          fp[j].read((char*)&buffer64, 8);  // Time stamp (ps)
		      timestamp = 1.0E-3*(double)buffer64;
          fp[j].read((char*)&buffer16, 2);  // Energy long (ADC)
		      energy = (float)buffer16;
          fp[j].read((char*)&buffer16, 2);  // Energy short (ADC)
          fp[j].read((char*)&buffer32, 4);  // Flags
          fp[j].read((char*)&buffer32, 4);  // Tracelength
          Tracelength = (int)(buffer32);
		  
	        // Reset variables
	        CFD = -1;
	        amplitude = -1;
	        paraL = 0;
	        paraS = 0;

	        // Get traces
	        for (i = 0; i < Tracelength; i++)
		    {
				fp[j].read((char*)&buffer16, 2);
			    pulse[i] = 16383 - (float)buffer16;
		        // Added traces
		        //if (j==0) {trace0->SetBinContent(i, pulse[i]);}
		        //if (j==1) {trace1->SetBinContent(i, pulse[i]);}
		    }

	        if(Tracelength > 1 && j < 16)
		    {
				// Process trace
				Analysis->CMA_Filter(pulse, Tracelength, CMAtrace, 10, pulse[0], 3.5 );
				for (i = 0; i < Tracelength; i++) 
				{
					pulse[i] -= CMAtrace[i];
					trace0->SetBinContent(i, pulse[i]);
					trace0C->SetBinContent(i, pulse[i]);
					traceCMA->SetBinContent(i, CMAtrace[i]);
					if (pulse[i] > amplitude) {amplitude = pulse[i]; pposition = i;}
				}

				// CFD timing
				//f1->SetParameters(1.0, (double)pposition, 0.1);
				//trace0->GetXaxis()->SetRangeUser(pposition - 5, pposition + 1);
				//trace0->Fit("f1","RQ");
				//mu = (float)f1->GetParameter(1);
				//sigma = (float)f1->GetParameter(2);

				//CFD = mu - sqrtf(1.38629*sigma*sigma);
				CFD = (float)pposition;

				// PSD integration
				offset = 12.0;
				if (pposition - 10 > 0 && pposition + 100 < Tracelength) {
					for (i = (pposition - 10); i < (pposition + 100); i++) {
						paraL += pulse[i];
					    if (i > pposition + offset) { paraS += pulse[i];}
				    }
				}
				
				// Input values into data struct 
				det.l = cal[j]*paraL;
				det.s = cal[j]*paraS;
				det.amp = amplitude;
				det.cfd = CFD;
				if (paraL != 0) {det.psd = paraS/paraL;}
				else {det.psd = -1;}
			
			}
			
			tt[j]->Fill();
			TEvt++;
			if (TEvt%1000==0) {cout << "\rChannel " << j << " Event counter: " << TEvt << flush;}
		}
		fp[j].close();
	    }
    }
	
	
	ff->cd();
	for (i = 0; i < 17; i++) { tt[i]->Write(); }
    ff->Close();

  cout << "\nFinsihed!" << endl;

  return 0;
}


