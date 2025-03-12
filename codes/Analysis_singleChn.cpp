/************************************************************************
 *
 *  Filename: Analysis.cpp
 *
 *  Description:
 *    Analysis routines for the a set of experiments performed at Notre dame
 *    November 2020.  This script inputs a processed root file from the 'scan'
 *    script, generates a gated neutron spectrum, runs the MLEM unfolder, and
 *    performs peak fittings / integration.
 *
 *    This script uses "UncertainNumber.h" class which tracks uncertanties
 *    during calculations.  Please see "UncertainNumber.h" for details.
 *
 *	Author(s):
 *     Michael T. Febbraro
 *
 *  Creation Date: 11/7/2019
 *  Last modified: 11/19/2020
 *
 *  If "error while loading shared libraries: libcore.so: ..." occurs, type
 *  "source `root-config --prefix`/bin/thisroot.sh"
 *
 * -----------------------------------------------------
 * 	Nuclear Astrophysics, Physics Division
 *  Oak Ridge National Laboratory
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <signal.h>
#include<sstream>

#include "PulseAnalysis.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TF1.h"

#include "UncertainNumber.h"



static UncertainNumber ALPHA(7.297352533e-3, 2.7e-11);
static UncertainNumber HBAR_C(197.3269602, 7.7e-6); 	// (MeV-fm)
static UncertainNumber AMU(931.4940954, 5.7e-8); 		// (MeV/c2)

#include "Nucleus.h"
#include "Reaction.h"

using namespace std;

void Analysis_singleChn() {

// ----------------------------------------------
// Nuclei definition (Z, A) using NIST database
// ----------------------------------------------
Nucleus Li7 (3, 7);
Nucleus Be7 (4, 7);
Nucleus C13 (6, 13);
Nucleus O16 (8, 16);
Nucleus O17 (8, 17);
Nucleus O18 (8, 18);
Nucleus Ne20 (10, 20);
Nucleus Ne21 (10, 21);
Nucleus Ne22 (10, 22);
Nucleus Mg26 (12, 26);
Nucleus Si29 (14, 29);
Nucleus a (2, 4);
Nucleus n (0, 1);
Nucleus p (1, 1);
Nucleus Si28(14,28);
Nucleus N13(7,13);
Nucleus N14(7,14);
Nucleus Mg25(12,25);
Nucleus B10(5,10);
Nucleus B11(5,11);

// ----------------------------------------------
// Create reactions
// ----------------------------------------------
Reaction rxn_13Can0 (C13, a, n, O16);
Reaction rxn_13Can1 (C13, a, n, O16);
Reaction rxn_13Can2 (C13, a, n, O16);
Reaction rxn_13Can3 (C13, a, n, O16);

Reaction rxn_7Lipn0 (Li7, p, n, Be7);
Reaction rxn_7Lipn1 (Li7, p, n, Be7);

Reaction rxn_25Mgan0 (Mg25, a, n, Si28);
Reaction rxn_25Mgan1 (Mg25, a, n, Si28); 
Reaction rxn_25Mgan2 (Mg25, a, n, Si28);

Reaction rxn_17Oan0 (O17, a, n, Ne20);
Reaction rxn_17Oan1 (O17, a, n, Ne20);

Reaction rxn_18Oan0 (O18, a, n, Ne21);
Reaction rxn_18Oan1 (O18, a, n, Ne21);

Reaction rxn_10Ban0 (B10, a, n, N13);
Reaction rxn_11Ban0 (B11, a, n, N14);

Reaction rxn_26Mgan0 (Mg26, a, n, Si29);
Reaction rxn_26Mgan1 (Mg26, a, n, Si29); 
Reaction rxn_26Mgan2 (Mg26, a, n, Si29);






// Excited states (MeV) in 16O with uncertanties
rxn_13Can1.setExcitation(6.0494, 0.0010);
rxn_13Can2.setExcitation(6.12989, 0.0004);
rxn_13Can3.setExcitation(6.9171, 0.0006);

// Excited states in 7Be with uncertanties
rxn_7Lipn1.setExcitation(0.42908, 0.00010);


//Excited states in 28Si with uncertanties
rxn_25Mgan1.setExcitation(1.77903,0.00011);
rxn_25Mgan2.setExcitation(4.61786,0.00004);


//Excited states in 20Ne 
rxn_17Oan1.setExcitation(0.870756,0.000020);


//Excited states in 21Ne 
rxn_18Oan1.setExcitation(1.98207,0.00009);


//Excited states in 29Si with uncertanties
rxn_26Mgan1.setExcitation(1.80874,0.00004);
rxn_26Mgan2.setExcitation(2.93833,0.00004);

// Neutron detector angles with uncertanties
double angles[10] =
	{45.0,
	 15.0,
	 30.0,
	 45.0,
	 60.0,
	 75.0,
	 90.0,
	 120.0,
	 140.0,
	 160.0};

// PSD gate definition
//double psd_A[10] =
    //{0.6, 1, 1, 1, 1, 1, 1, 1, 1, 1};
//double psd_B[10] =
    //{0.01, 0.02, 0.09, 0.06, 0.05, 0.04, 0.06, 0.05, 0.05, 0.04};
    
double psd_A_neutron_upper = 0.7643544;
double psd_B_neutron_upper = -0.0000280;
double psd_C_neutron_upper = 0.1710680;

double psd_A_neutron_lower = -1.3253611;
double psd_B_neutron_lower = -0.0000362;
double psd_C_neutron_lower = 0.1825187;

double psd_A_gamma_upper = 0.4175856;
double psd_B_gamma_upper = -0.0000025;
double psd_C_gamma_upper = 0.0068125;

int low, high;
int position;
double bin;

double C13an0[10];
double uncer_Li7pn0[10];

char 	filename[250];

char buffer[200];

string line,line1, runID;
string labelX = "E_{n} (MeV)";
string labelY = "";
TH1F *s = new TH1F("s","Unfolded Incident Spectrum",400,0.05,20);
s->GetXaxis()->SetTitle(labelX.c_str());
s->GetYaxis()->SetTitle(labelY.c_str());

//cout << "Root file name to be opened: ";
//cin >> filename;

//TFile *f = new TFile(filename);
//TTree *T;
    
    
ifstream fin;
fin.open("25Mg_a_n_columns.txt");
int runNumber, charge;
double ebeam;

  

// ----------------------------------------------
// Load data set infomation
// ----------------------------------------------



while (getline(fin, line1))
{

	sscanf(line1.c_str(),"%d %d %lf",&runNumber , &charge, &ebeam);
    
   

	runID = "Extracted_rootfiles/run_" + to_string(runNumber) + ".root";
	ebeam /= 1000; // convert to MeV
	rxn_13Can0.setEbeam(ebeam,0);
	rxn_13Can0.printSummary();
	cout << rxn_13Can0.getEjectileElab_thlab(55) << endl;
	
	rxn_25Mgan0.setEbeam(ebeam,0);
	rxn_25Mgan0.printSummary();
	cout << rxn_25Mgan0.getEjectileElab_thlab(55) << endl;
	
	
	rxn_25Mgan1.setEbeam(ebeam,0);
	rxn_25Mgan1.printSummary();
	cout << rxn_25Mgan1.getEjectileElab_thlab(55) << endl;
	
	rxn_25Mgan2.setEbeam(ebeam,0);
	rxn_25Mgan2.printSummary();
	cout << rxn_25Mgan2.getEjectileElab_thlab(55) << endl;
	
	
	rxn_26Mgan0.setEbeam(ebeam,0);
	rxn_26Mgan0.printSummary();
	cout << rxn_26Mgan0.getEjectileElab_thlab(55) << endl;
	
	
	rxn_26Mgan1.setEbeam(ebeam,0);
	rxn_26Mgan1.printSummary();
	cout << rxn_26Mgan1.getEjectileElab_thlab(55) << endl;
	
	rxn_26Mgan2.setEbeam(ebeam,0);
	rxn_26Mgan2.printSummary();
	cout << rxn_26Mgan2.getEjectileElab_thlab(55) << endl;
	
	
	rxn_17Oan0.setEbeam(ebeam,0);
	rxn_17Oan0.printSummary();
	cout << rxn_17Oan0.getEjectileElab_thlab(55) << endl;
	
	rxn_17Oan1.setEbeam(ebeam,0);
	rxn_17Oan1.printSummary();
	cout << rxn_17Oan1.getEjectileElab_thlab(55) << endl;
	
	rxn_18Oan0.setEbeam(ebeam,0);
	rxn_18Oan0.printSummary();
	cout << rxn_18Oan0.getEjectileElab_thlab(55) << endl;
	
	
	rxn_18Oan1.setEbeam(ebeam,0);
	rxn_18Oan1.printSummary();
	cout << rxn_18Oan1.getEjectileElab_thlab(55) << endl;
	
	
	rxn_10Ban0.setEbeam(ebeam,0);
	rxn_10Ban0.printSummary();
	cout << rxn_10Ban0.getEjectileElab_thlab(55) << endl;
	
    rxn_11Ban0.setEbeam(ebeam,0);
	rxn_11Ban0.printSummary();
	cout << rxn_11Ban0.getEjectileElab_thlab(55) << endl;
	
	rxn_26Mgan0.setEbeam(ebeam,0);
	rxn_26Mgan0.printSummary();
	cout << rxn_26Mgan0.getEjectileElab_thlab(55) << endl;
	
	rxn_26Mgan1.setEbeam(ebeam,0);
	rxn_26Mgan1.printSummary();
	cout << rxn_26Mgan1.getEjectileElab_thlab(55) << endl;
	
	rxn_26Mgan2.setEbeam(ebeam,0);
	rxn_26Mgan2.printSummary();
	cout << rxn_26Mgan2.getEjectileElab_thlab(55) << endl;
	
	rxn_7Lipn0.setEbeam(ebeam,0);
	rxn_7Lipn0.printSummary();
	cout << rxn_7Lipn0.getEjectileElab_thlab(55) << endl;
	
	rxn_7Lipn1.setEbeam(ebeam,0);
	rxn_7Lipn1.printSummary();
	cout << rxn_7Lipn1.getEjectileElab_thlab(55) << endl;
	
	
	
	


	

	TFile *f = new TFile(runID.c_str());

	//TChain chain_p1("T");
	//chain_p1.Add(filename);

	TH1F *h_p1 = new TH1F("h_p1","h_p1",600,0,6000);
	char channel[50], condition[250];

	for (int i = 0; i < 1; i++)
	{
		sprintf(channel, "T%d",i);
		TTree *T = (TTree*)f->Get(channel);
		//sprintf(condition, "d.amp<15500 && d.amp>50 && d.psd > %f/TMath::Sqrt(d.l) + %f + (0.01/6000.)*d.l", psd_A[i], psd_B[i]);
		//sprintf(condition, "d.amp<15500 && d.amp>50 && d.psd > 0.1");
        sprintf(condition,"d.psd > %f/TMath::Sqrt(d.l) + %f*d.l + %f && d.psd > %f/TMath::Sqrt(d.l) + %f*d.l + %f && d.psd < %f/TMath::Sqrt(d.l) + %f*d.l + %f", psd_A_gamma_upper,psd_B_gamma_upper,psd_C_gamma_upper,psd_A_neutron_lower,psd_B_neutron_lower,psd_C_neutron_lower,psd_A_neutron_upper,psd_B_neutron_upper,psd_C_neutron_upper);
        
		
		T->Draw("d.l>>h_p1", condition,"");

		// ----------------------------------
		// Generate input file for unfolder
		// ----------------------------------
		ofstream fop;
		fop.open ("spectrum.spe");

		fop << "500 0.01 500\n";
		for (int i = 0; i < 500; i++) {fop << h_p1->GetBinContent(i) << "\n";}
		fop.close();

		// ----------------------------------
		// Run spectrum unfolder
		// ----------------------------------
		string cmd = "./MLEM_OU.exe ";
		//string cmd = "/home/shahina/Documents/PhD/25Mg_codes_Shahina/Unfolding/MLEM_Code/Unfolding_shahina.exe ";
		system(cmd.c_str());

		// ----------------------------------
		// Load unfolded spectrum
		// ----------------------------------

		string filenames = "Unfolded.out", line;
		ifstream fp;
		fp.open(filenames.c_str());
        
        // Open the output file for storing the unfolded counts
        string outputFilename = "Unfolded_output_files/" + to_string(runNumber) + ".txt";
        ofstream outputFile(outputFilename);
        
        // Open the output file for storing the neutron energy information
        string outputFile_counts_name = "Unfolded_count_files/" + to_string(runNumber) + ".txt";
        ofstream outputFile_counts(outputFile_counts_name);
        
        
        
		float x = 0.05;
		s->Reset();
		while (getline(fp, line))
		{
            double value = atof(line.c_str());
            outputFile << x << " " << value << endl;
            s->Fill(x, value);
    	    x = x + 0.05;
		}
		fp.close();
        
        
        

		// ----------------------------------
		// Integrate peaks
		// ----------------------------------
		
		//if (position == 1) { cout << "P=1" << endl; bin = rxn_13Can0.getEjectileElab_thlab(angles[i])/0.05;}
		bin = rxn_13Can0.getEjectileElab_thlab(angles[i] + 8.3)/0.05;
		cout<<"Neutron Energy: 13C(a,n0) " << rxn_13Can0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 25Mg(a,n0) " << rxn_25Mgan0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 25Mg(a,n1) " << rxn_25Mgan1.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 25Mg(a,n2)" << rxn_25Mgan2.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 26Mg(a,n0) " << rxn_26Mgan0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 26Mg(a,n1) " << rxn_26Mgan1.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 26Mg(a,n2)" << rxn_26Mgan2.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 17O(a,n0)" << rxn_17Oan0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 17O(a,n1)" << rxn_17Oan1.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 18O(a,n0)" << rxn_18Oan0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 18O(a,n1)" << rxn_17Oan0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 10B(a,n0)" << rxn_10Ban0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 11B(a,n0)" << rxn_11Ban0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 7Li(p,n0)" << rxn_7Lipn0.getEjectileElab_thlab(55)<<endl;
		cout<<"Neutron Energy: 7Li(p,n1)" << rxn_7Lipn1.getEjectileElab_thlab(55)<<endl;
		
		
		
		//cout<<"bin: "<<bin<<endl;
		//low = (int)(bin) - 20;
		//high = (int)(bin) + 20;
		//fop_counts << s->Integral(low,high) << "  ";
        outputFile_counts << runID << " " << ebeam << " " << "13Can0"<<"  "<<rxn_13Can0.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "25Mgan0"<<"  "<<rxn_25Mgan0.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "25Mgan1"<<"  "<<rxn_25Mgan1.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "25Mgan2"<<"  "<<rxn_25Mgan2.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "26Mgan0"<<"  "<<rxn_26Mgan0.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "26Mgan1"<<"  "<<rxn_26Mgan1.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "26Mgan2"<<"  "<<rxn_26Mgan2.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "17Oan0"<<"  "<<rxn_17Oan0.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "17Oan1"<<"  "<<rxn_17Oan1.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "18Oan0"<<"  "<<rxn_18Oan0.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "18Oan1"<<"  "<<rxn_18Oan1.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "10Ban0"<<"  "<<rxn_10Ban0.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "11Ban0"<<"  "<<rxn_11Ban0.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "7Lipn0"<<"  "<<rxn_7Lipn0.getEjectileElab_thlab(55) << endl;
        outputFile_counts << runID << " " << ebeam << " " << "7Lipn1"<<"  "<<rxn_7Lipn1.getEjectileElab_thlab(55) << endl;
		
        outputFile_counts << endl;
        outputFile_counts.close();
        outputFile.close();
		
	}
    
 }




return;
}
