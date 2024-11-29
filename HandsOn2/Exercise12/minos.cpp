#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooHistFunc.h"
#include "RooGenericPdf.h"
#include "RooMinimizer.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit;

void minos()
{
    // S e t u p   m o d e l
    // ---------------------

    // Declare neutrino energy x in range 0.5 to 14 GeV
    RooRealVar energy("energy", "Reconstructed Neutrino Energy", 0.5, 14);

    // Load unbinned dataset of the neutrinos observed by MINOS from minos_2013_data.dat
    RooDataSet data = *RooDataSet::read("minos_2013_data.dat", energy, "v");

    // Load unbinned dataset of the Monte Carlo generated non oscillated neutrino energy distribution from minos_2013_mc.dat
    RooDataSet mc_noosc = *RooDataSet::read("minos_2013_mc.dat", energy, "v");

    // Make a histogram function out of it
    RooDataSet *dd = (RooDataSet *)mc_noosc.reduce(RooArgSet(energy));
    RooDataHist *dh_mc_noosc = dd->binnedClone();
    RooHistFunc func_noosc{"func_mc_noosc", "No oscillation", energy, *dh_mc_noosc, 2};

    // Create mixing and dm2
    RooRealVar mixing("mixing", "sin^2(2*theta)", 0.5, 0, 1);
    RooRealVar dm2("dm2", "Delta m^2", 0.0024, 0.0001, 0.01);

    // Create the oscillation probably
    RooFormulaVar prob_osc("prob_osc", "Oscillation Probability", "1 - @0*sin(1.267*730*@1/@2)^2", RooArgSet(mixing, dm2, energy));

    // Create final PDF of the energy distribution of oscillated neutrino
    RooGenericPdf model{"model", "model", "@0*@1", RooArgSet(prob_osc, func_noosc)};

    // Fit model to data
    model.fitTo(data);

    // Plot the observed data and the model
    RooPlot *frame = energy.frame();
    data.plotOn(frame);
    model.plotOn(frame);

    // Save the plot
    TCanvas c1("canvas", "MINOS Data Fit", 800, 600);
    frame->Draw();
    c1.SaveAs("minos_data.png");

    // M I N I M I Z E  L I K L I H O O D  B Y  H A N D
    // ------------------------------------------------

    // Construct function object representing negative logarithm of likelihood
    RooAbsReal *nll = model.createNLL(data);

    // Minimize nll w.r.t. its parameters with Migrad
    RooMinimizer m(*nll);
    m.setVerbose(kTRUE); // enable verbose mode
    m.migrad();

    // Print parameter values after MIGRAD
    mixing.Print();
    dm2.Print();

    // Disable verbose mode
    m.setVerbose(kFALSE);

    // HESSE error calculation
    m.hesse();

    // Print parameter values after HESSE error calculation
    mixing.Print();
    dm2.Print();

    // MINOS error calculation for dm2
    m.minos(dm2);

    // Print dm2 value after MINOS error calculation
    mixing.Print();
    dm2.Print();

    // Save a snapshot of the fit results
    RooFitResult *fitResult = m.save();
    m.Print("v");

    // Contour plot for dm2 vs mixing at 1, 2, and 3 sigma
    TCanvas c2("contour", "Contour Plot", 800, 600);
    RooPlot* contour = m.contour(dm2, mixing, 2.30, 6.18, 11.83);
    contour->SetTitle("Contour Plot of dm2 vs mixing");
    contour->GetXaxis()->SetTitle("dm2");
    contour->GetYaxis()->SetTitle("mixing");
    contour->Draw();

    // Save the contour plot
    c2.SaveAs("minos_likelihood.png");
}
