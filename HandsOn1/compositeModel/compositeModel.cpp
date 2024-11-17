#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TCanvas.h"

void compositeModel() {
    // Define the observable x
    RooRealVar x("x", "Observable", -10, 10);

    // Define the Gaussian signal component
    RooRealVar mean("mean", "Mean of Gaussian", 0);
    RooRealVar sigma("sigma", "Width of Gaussian", 3);
    RooGaussian sig("sig", "Signal Gaussian PDF", x, mean, sigma);

    // Define the exponential background component
    RooRealVar tau("tau", "Tau of exponential", 10, 0.1, 50);
    RooFormulaVar C("neg_tau", "-1/tau", "-1./tau", RooArgList(tau));
    RooExponential bkg("bkg", "Background Exponential PDF", x, C);

    // Define the signal fraction parameter
    RooRealVar fsig("fsig", "Signal fraction", 0.5, 0., 1.);

    // Composite model: fsig * sig + (1 - fsig) * bkg
    RooAddPdf model("model", "Signal + Background", RooArgList(sig, bkg), RooArgList(fsig));

    // Generate a binned dataset
    RooDataHist *data = model.generateBinned(x, 1000); // 1000 events

    // Plot the data and model
    RooPlot* xframe = x.frame();
    xframe->SetTitle("Composite Model"); 

    data->plotOn(xframe, RooFit::DataError(RooAbsData::SumW2)); // Plot data
    model.plotOn(xframe);                                      // Plot model
    model.plotOn(xframe, RooFit::Components("bkg"), RooFit::LineStyle(kDashed)); // Background component
    model.plotOn(xframe, RooFit::Components("sig"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed)); // Signal component

    // Fit the model to the data
    model.fitTo(*data);

    // Plot the fit results
    TCanvas* c = new TCanvas("c", "Composite Model", 800, 600);
    xframe->Draw();
    c->SaveAs("composite_model_fit.png");
}
