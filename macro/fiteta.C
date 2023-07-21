#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include<TMultiGraph.h>
#include "TROOT.h"
#include<TLegend.h>
#include<TFile.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "/sphenix/u/shuhang98/AtlasStyle.C"
#include "/sphenix/u/bseidlitz/plotstyle/AtlasUtils.C"
#include <TGraphAsymmErrors.h>

void fiteta() {

// Open the root file
	TFile *file = TFile::Open("etamassrecoout.root");

	// Retrieve the TH1F histograms
	TH1F *hetamass45 = (TH1F *)file->Get("hetamasspt45");
	TH1F *hetamassmix45 = (TH1F *)file->Get("hetamassptmix45");
	TH1F *hetamass56 = (TH1F *)file->Get("hetamasspt56");
	TH1F *hetamassmix56 = (TH1F *)file->Get("hetamassptmix56");
	TH1F *hetamass67 = (TH1F *)file->Get("hetamasspt67");
	TH1F *hetamassmix67 = (TH1F *)file->Get("hetamassptmix67");
	TH1F *hetamass78 = (TH1F *)file->Get("hetamasspt78");
	TH1F *hetamassmix78 = (TH1F *)file->Get("hetamassptmix78");
	TH1F *hetamass8u = (TH1F *)file->Get("hetamasspt8u");
	TH1F *hetamassmix8u = (TH1F *)file->Get("hetamassptmix8u");

	//sum histograms
	TH1F *hetamass = (TH1F *)file->Get("hetamass");
	//hetamass->Add(hetamass56);
	hetamass->Add(hetamass67);
	hetamass->Add(hetamass78);
	hetamass->Add(hetamass8u);
	TH1F *hetamassmix = (TH1F *)file->Get("hetamassmix");
	//hetamassmix->Add(hetamassmix56);
	hetamassmix->Add(hetamassmix67);
	hetamassmix->Add(hetamassmix78);
	hetamassmix->Add(hetamassmix8u);
	//add titles
	hetamass->SetXTitle("#gamma #gamma invariant Mass(GeV)");


	// Call Sumw2() on both histograms
	hetamass->Sumw2();
	hetamassmix->Sumw2();
	//rebin both histograms
	hetamass->Rebin(2);
	hetamassmix->Rebin(2);


	// Divide the histograms
	TH1F *hdiv = (TH1F *)hetamass->Clone("hdiv");
	hdiv->Divide(hetamassmix);

	// Define the fit function as the sum of an exponential and a Gaussian
	TF1 *fitFunc = new TF1("fitFunc", "[0]*exp(-[1]*x) + [2]/(sqrt(2*TMath::Pi()) * [4])*exp(-0.5*((x-[3])/[4])**2)", 0.3, 1);

	// Set initial parameter values
	fitFunc->SetParameter(0, 1);
	fitFunc->SetParameter(1, 1);
	fitFunc->SetParLimits(1, 0, 100);
	fitFunc->SetParameter(2, 0.01);
	fitFunc->SetParLimits(2, 0, 1);
	fitFunc->SetParameter(3, 0.6);
	fitFunc->SetParLimits(3, 0.4, 0.8);
	fitFunc->SetParameter(4, 0.05);
	fitFunc->SetParLimits(4, 0.02, 0.2);

	// Fit the histogram using the fit function
	hdiv->Fit("fitFunc", "R"); // "R" option for fitting within the specified range (0.2, 1)

	// Draw the histogram and fit function
	//TCanvas *canvas = new TCanvas("canvas", "Eta Mass Histogram", 800, 600);
	//hdiv->Draw();
	//fitFunc->Draw("SAME");
	// Define the exponential function independently
	TF1 *exponentialFunc = new TF1("exponentialFunc", "[0]*exp(-[1]*x)", 0, 2);

	// Set the parameters to the corresponding values from the original fit function
	exponentialFunc->SetParameter(0, fitFunc->GetParameter(0));
	exponentialFunc->SetParameter(1, fitFunc->GetParameter(1));

	// Scale hetamassmix by the exponential function
	TH1F *hetamassmixScaled = (TH1F *)hetamassmix->Clone("hetamassmixScaled");
	hetamassmixScaled->Multiply(exponentialFunc);

	// Subtract the scaled hetamassmix from hetamass
	TH1F *hetamassSubtracted = (TH1F *)hetamass->Clone("hetamassSubtracted");
	hetamassSubtracted->Add(hetamassmixScaled, -1);


	// Define the fit function as the sum of an exponential and a Gaussian
	TF1 *fg = new TF1("fg", "[0]/(sqrt(2*TMath::Pi()) * [2])*exp(-0.5*((x-[1])/[2])**2)", 0.35, 0.85);

	fg->SetParameter(0, 1);
	fg->SetParameter(1, 0.55);
	fg->SetParameter(2, 0.2);
	fg->SetParLimits(2, 0, 1);

	hetamassSubtracted->Fit("fg", "R");

	// Draw the subtracted histogram
	TCanvas *canvas2 = new TCanvas("canvas2", "Subtracted Eta Mass Histogram", 800, 600);
	hdiv->Draw();



}
