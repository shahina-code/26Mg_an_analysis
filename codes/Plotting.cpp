
#include <iostream>
#include <fstream>
#include <string>
#include <signal.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"


void Plotting(){

//system("./MLEM-N.exe 400");

string labelX;
string labelY;

double xdata[200], ydata[200];
float x, y;
string line, filename;
int j = 0;

labelX = "E_{n} (MeV)";
labelY = "";
TH1F *s = new TH1F("s","Unfolded Incident Spectrum",400,0.05,20);
s->GetXaxis()->SetTitle(labelX.c_str());
s->GetYaxis()->SetTitle(labelY.c_str());

labelX = "E_{n} (MeV)";
labelY = "";
TH1F *s_tof = new TH1F("s_tof","ToF Spectrum",400,0.05,20);
s_tof->GetXaxis()->SetTitle(labelX.c_str());
s_tof->GetYaxis()->SetTitle(labelY.c_str());

TH1F *error = new TH1F("error","",400,0.05,20);
labelX = "Light Response (MeVee)";
labelY = "";
TH1F *e = new TH1F("e","Estimated Response Spectrum",1500,0.01,15);
e->GetXaxis()->SetTitle(labelX.c_str());
e->GetYaxis()->SetTitle(labelY.c_str());

labelX = "Light Response (MeVee)";
labelY = "";
TH1F *i = new TH1F("i","Response Spectrum",1500,0.01,15);
i->GetXaxis()->SetTitle(labelX.c_str());
i->GetYaxis()->SetTitle(labelY.c_str());

labelX = "Light Response (MeVee)";
labelY = "#frac{#bar{E} - #bar{S}}{#sqrt{#bar{S}}}";
TH1F *r = new TH1F("r","Residual",1500,0.01,15);
r->GetXaxis()->SetTitle(labelX.c_str());
r->GetYaxis()->SetTitle(labelY.c_str());

TH1F *upper = new TH1F("upper","",1500,0.01,15);
TH1F *lower = new TH1F("lower","",1500,0.01,15);

ifstream fp;

//cout << "File to plot: ";
//cin >> filename;

	filename = "error.out";
	fp.open(filename.c_str());
	x = 0.05;
	while (getline(fp, line))
	{
			error->Fill(x, atof(line.c_str()));
			//cout << line << endl;
			x = x + 0.05;
	}
	fp.close();

	filename = "Unfolded.out";
	fp.open(filename.c_str());
	x = 0.05;
	while (getline(fp, line))
	{
			s->Fill(x, atof(line.c_str()));
			//y = sqrt(error->GetBinContent(x))/atof(line.c_str());
			//s->SetBinError(x,y);
			//cout << line << endl;
			x = x + 0.05;
	}
	fp.close();

	// *************

	filename = "tof.spe";
	fp.open(filename.c_str());
	x = 0.05;
	while (getline(fp, line))
	{
			s_tof->Fill(x, atof(line.c_str()));
			//y = sqrt(error->GetBinContent(x))/atof(line.c_str());
			//s->SetBinError(x,y);
			//cout << line << endl;
			x = x + 0.05;
	}
	fp.close();

	// *************

	filename = "estimate.out";
	fp.open(filename.c_str());
	x = 0.01;
	while (getline(fp, line))
	{
			e->Fill(x, atof(line.c_str()));
			//cout << line << endl;
			x = x + 0.01;
	}
	fp.close();

	filename = "spectrum.in";
	fp.open(filename.c_str());
	x = 0.01;
	while (getline(fp, line))
	{
			i->Fill(x, atof(line.c_str()));
			r->Fill(x, atof(line.c_str()));
			//cout << line << endl;
			x = x + 0.01;
	}
	fp.close();



	TCanvas *c1 = new TCanvas("c1","MLEM Unfolding",700,700);
	c1->Divide(1,3);

	c1->cd(1);
	s->GetXaxis()->SetRange(0,300);
	s->Draw("hist");
	gStyle->SetOptStat(1);
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.15);
	//s->SetFillColor(17);
	s->SetLineColor(1);
	s->GetXaxis()->SetLabelSize(0.07);
	s->GetYaxis()->SetLabelSize(0.07);
	s->GetXaxis()->SetTitleSize(0.06);
	s->GetXaxis()->SetTitleSize(0.06);
	s->GetYaxis()->SetTitleSize(0.07);
	s->GetYaxis()->SetTitleOffset(0.65);

	s_tof->Draw("same hist");

	cout << "Counts: " << s->Integral(86,153) << endl;



	c1->cd(2);
	i->Draw("hist");

	gStyle->SetOptStat(0);
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.15);
	i->GetXaxis()->SetRange(0,500);
	i->GetXaxis()->SetLabelSize(0.07);
	i->GetYaxis()->SetLabelSize(0.07);
	i->GetXaxis()->SetTitleSize(0.06);
	i->GetXaxis()->SetTitleSize(0.06);
	i->GetYaxis()->SetTitleSize(0.07);
	e->Draw("same hist");
	e->SetLineColor(1);



	c1->cd(3);
	r->Add(e,-1);

	for (j = 0; j<1500; j++)
	{
		if(e->GetBinContent(j) != 0)
		{
			y = r->GetBinContent(j)/sqrt(e->GetBinContent(j));
			r->SetBinContent(j,y);
		}
		upper->SetBinContent(j,1);
		lower->SetBinContent(j,-1);
	}

	r->Draw("hist");
	gStyle->SetOptStat(0);
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.15);
	r->GetXaxis()->SetRange(0,500);
	r->GetXaxis()->SetLabelSize(0.07);
	r->GetYaxis()->SetLabelSize(0.07);
	r->GetXaxis()->SetTitleSize(0.06);
	r->GetXaxis()->SetTitleSize(0.06);
	r->GetYaxis()->SetTitleSize(0.07);
	r->GetYaxis()->SetTitleOffset(0.65);
	r->GetYaxis()->CenterTitle();
	upper->Draw("same");
	lower->Draw("same");
	upper->SetLineColor(2);
	upper->SetLineStyle(2);
	lower->SetLineColor(2);
	lower->SetLineStyle(2);

	// Print Stats
	cout << (s->Integral(10,28) - s_tof->Integral(10,28))/s_tof->Integral(10,28) << endl;
	cout << (s->Integral(28,40) - s_tof->Integral(28,38))/s_tof->Integral(28,38) << endl;
	cout << (s->Integral(40,54) - s_tof->Integral(38,48))/s_tof->Integral(38,48) << endl;
	cout << (s->Integral(56,80) - s_tof->Integral(56,76))/s_tof->Integral(56,76) << endl;

	cout << endl;
	// Print Stats
	cout << sqrt(s_tof->Integral(10,28)) << endl;
	cout << sqrt(s_tof->Integral(28,38)) << endl;
	cout << sqrt(s_tof->Integral(38,48)) << endl;
	cout << sqrt(s_tof->Integral(48,80)) << endl;
}
