#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "TCanvas.h"

 using namespace RooFit;

void rarestB0()
{
    // Define x variable
    RooRealVar x("x", "x", 5000, 5600);

    // Load the unbinned dataset from the file
    RooDataSet* data = RooDataSet::read("rarest_b0_decay.dat", x, "v");

    // Define exponential background pdf
    RooRealVar tau("tau", "tau", -0.001, -1.0, 0.0);
    RooExponential background("bkg", "Exponential Background", x, tau);

    // Define gaussian peak around B0 mass
    RooRealVar mean_B0("mean_B0", "Mean of B0", 5279, 5270, 5290);
    RooRealVar sigma_B0("sigma_B0", "Width of B0", 10, 0.1, 50);
    RooGaussian gauss_B0("gauss_B0", "Gaussian for B0", x, mean_B0, sigma_B0);

    // Define gaussian peak around B0s mass
    RooRealVar mean_Bs("mean_Bs", "Mean of Bs", 5367, 5350, 5380);
    RooRealVar sigma_Bs("sigma_Bs", "Width of Bs", 10, 0.1, 50);
    RooGaussian gauss_Bs("gauss_Bs", "Gaussian for Bs", x, mean_Bs, sigma_Bs);

    // Define yields
    RooRealVar yield_B0("yield_B0", "Yield of B0", 100, 0, 1000);
    RooRealVar yield_Bs("yield_Bs", "Yield of Bs", 50, 0, 500);
    RooRealVar yield_bkg("yield_bkg", "Yield of Background", 1000, 0, 10000);

    // Create composite model
    RooAddPdf model("model", "Composite Model", RooArgList(gauss_B0, gauss_Bs, background), RooArgList(yield_B0, yield_Bs, yield_bkg));

    // Fit model to data
    model.fitTo(*data);

    // Create Frame
    TCanvas* c = new TCanvas("c", "B^{0} Invariant Mass Distribution", 800, 600);
    RooPlot* frame = x.frame();

    // Make Plot
    data->plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components("bkg"), LineColor(kRed));
    model.plotOn(frame, Components("gauss_B0"), LineColor(kGreen));
    model.plotOn(frame, Components("gauss_Bs"), LineColor(kBlue));

    // Draw
    frame->Draw();
    c->SaveAs("InvariantMassDistribution.pdf");

    // Compute residual and pull histograms
    RooHist* residuals = frame->residHist();
    RooHist* pulls = frame->pullHist();

    // Plot residuals
    TCanvas* c_resid = new TCanvas("c_resid", "Residuals", 800, 600);
    RooPlot* frame_resid = x.frame(Title("Residuals"));
    frame_resid->addPlotable((RooPlotable*)residuals, "P"); // Cast RooHist to RooPlotable
    frame_resid->Draw();
    c_resid->SaveAs("residuals_fixed.pdf");

    // Plot pulls
    TCanvas* c_pull = new TCanvas("c_pull", "Pulls", 800, 600);
    RooPlot* frame_pull = x.frame(Title("Pulls"));
    frame_pull->addPlotable((RooPlotable*)pulls, "P"); // Cast RooHist to RooPlotable
    frame_pull->Draw();
    c_pull->SaveAs("pulls_fixed.pdf");
}