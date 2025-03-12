/************************************************************************
 *  Description: MLEM matrix text file plotter
 *	Author(s):
 *     Rebecca Toomey
 *
 *  Creation Date: 1/20/2021
 * -----------------------------------------------------
 *  Rutgers University
 *
 *************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <signal.h>
//#include "PulseAnalysis.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TPolyMarker3D.h"

int matrixmlem_plot(){

  ifstream matrix;
  matrix.open("matrix.txt");
  stringstream ss (stringstream::in | stringstream::out);

  int i=0,j=0;
  string line;
  int val;

  TH2F *h = new TH2F("h","h",300,0,15,1000,0,10000);

  for(i=0;i<302;i++){

    getline(matrix,line);

    if(i==0){continue;} //skip first line

    ss.clear();
    ss << line;

    while(ss>>val){

      //if ((i-1)*0.05 < 5.50){
       h->Fill((i-1)*0.05,(j-1)*10.,val);
      //}
      // else{
      //   h->Fill((i-1)*0.05, (j-1)*10.*1.5, val);
      // }
       //cout << i*0.05 << " " << j*10. << " " << val << endl;
       j++;

    }

    j=0;

  }

  h->Draw("col");

  // ofstream fop;
  // fop.open ("matrix_stil100.txt");
  //
  // // <light bins> <energy bins> <norm (only used for MC)> <bin width in MeVee>
  // fop << "1000 300 1 0.01\n";
  //
  // for (int i = 0; i < 300; i++)
  // {
  //   for (int j = 0; j < 1000; j++)
  //   {
  //       if ((i-1)*0.05 < 5.50){
  //      fop << h->GetBinContent(i, j) << " ";
  //      }
  //      else{
  //        if ((j-1)*10.<100){
  //          fop << 0 << " ";
  //         }
  //        else{
  //          fop << h->GetBinContent(i, j) << " ";
  //         }
  //      }
  //   }
  //   fop << "\n";
  // }
  // fop.close();

return 1;
}
