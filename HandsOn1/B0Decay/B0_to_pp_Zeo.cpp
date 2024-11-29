#include <iostream>
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooBinning.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TAxis.h"
#include "TStyle.h"
using namespace RooFit;

void B0_to_pp()
 {
    RooRealVar x("x", "x", 5080, 5560);
    RooRealVar mean_B0("mean_B0", "mean of gaussian for B0", 5280, 5260,5300);
    RooRealVar sigma("sigma_B0", "width of gaussian for B0", 20, 5, 50);
    RooRealVar mean_Bs("mean_Bs", "mean of gaussian for Bs", 5366, 5346,5386);
    RooRealVar lambda("lambda", "lambda of exponential", -0.0001, -1, -0.000001);

    RooGaussian sig_B0("sig_B0", "gaussian PDF for B0", x, mean_B0, sigma);
    RooGaussian sig_Bs("sig_Bs", "gaussian PDF for B0", x, mean_Bs, sigma);
    RooExponential exp("exp", "exponential PDF", x, lambda);

    //RooRealVar nsig_B0("nsig_B0", "signal count B0", 30, 0, 50);
    //RooRealVar nsig_Bs("nsig_Bs", "signal count Bs", 5, 0, 10);
    //RooRealVar nbkg("nbkg", "background count", 100, 0, 500);

    RooRealVar B0_frac("B_frac", "fraction of B0", 0.1, 0, 1);
    RooRealVar BS_frac("bkg_frac", "fraction of BS", 0.1, 0, 1);


    RooAddPdf model("model", "full model", RooArgList(sig_B0, sig_Bs, exp), RooArgList(B0_frac, BS_frac));

    RooDataSet data = *RooDataSet::read("data/rarest_b0_decay.dat", x, "v");

    float bin_width = 20;
    int nbins = (x.getMax() - x.getMin()) / bin_width;
    x.setBinning(RooBinning(nbins, x.getMin(), x.getMax()));
    model.fitTo(data, PrintLevel(-1));

    mean_B0.Print();
    mean_Bs.Print();
    sigma.Print();
    lambda.Print();


    RooPlot *frame1 = x.frame(Title("Data + Fit"));
    RooPlot *frame2 = x.frame(Title("Residual Distribution"));
    RooPlot *frame3 = x.frame(Title("Pulls"));
    data.plotOn(frame1);
    model.plotOn(frame1, LineColor(kBlue));

    RooHist *resid = frame1->residHist(); // get residuals
    RooHist *pull = frame1->pullHist(); // get pulls
    
    frame2->addPlotable(resid, "P");
    frame3->addPlotable(pull, "P");

    model.plotOn(frame1, Components("exp"), LineStyle(kDashed), LineColor(kGreen)); // Background component
    model.plotOn(frame1, Components("sig_B0"), LineStyle(kDashed), LineColor(kRed)); // B0 component
    model.plotOn(frame1, Components("sig_Bs"), LineStyle(kDashed), LineColor(kGray)); // Bs component

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    frame1->Draw();
    c1->SaveAs("output/fit.pdf");

    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 600);
    c2->Divide(2);
    c2->cd(1);
    gPad->SetLeftMargin(0.15);
	frame2->Draw();
	c2->cd(2);
    gPad->SetLeftMargin(0.05);
	frame3->Draw();
    c2->SaveAs("output/res_and_pulls.pdf");

    //TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    //gStyle->SetPalette(1) ;
    //fit_results->correlationHist()->Draw(“colz”) ; 
}
 