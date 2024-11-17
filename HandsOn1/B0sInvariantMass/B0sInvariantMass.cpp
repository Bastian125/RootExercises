#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TFile.h"
using namespace RooFit;

TH1 *makeTH1();

void B0sInvariantMass()
{
    // Create a macro to open the “B0sInvariantMass.root” and import the corresponding binned dataset.
    TFile *file = TFile::Open("B0sInvariantMass.root");
    TH1D *hist = (TH1D *)file->Get("massaB0");

    RooRealVar x("x", "x", -10, 10);
    RooDataHist dh("dh", "dh", x, Import(*hist));

    // Create a Breit-Wigner model.
    RooRealVar mean("mean", "mean", 5, -10, 10);
    RooRealVar sigma("sigma", "sigma", 1, 0.1, 10);
    RooBreitWigner bw("bw", "bw", x, mean, sigma);

    // Fit the model to the binned dataset.
    bw.fitTo(dh, PrintLevel(-1));

    // Create a Gaussian function and fit to the data.
    RooRealVar meanG("meanG", "meanG", 5, -10, 10);
    RooRealVar sigmaG("sigmaG", "sigmaG", 1, 0.1, 10);
    RooGaussian gauss("gauss", "gauss", x, meanG, sigmaG);
    gauss.fitTo(dh, PrintLevel(-1));

    // Plot the data, and the BW and Gaussian distribution to the same canvas.
    RooPlot *frame = x.frame(Title("Breit-Wigner / Gauss fit"));
    dh.plotOn(frame);
    bw.plotOn(frame);
    gauss.plotOn(frame, LineColor(kRed));

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    frame->Draw();
    c->SaveAs("B0_mass_dist_fit.pdf");

    // Compare the fitted value with the particle mass reported in the Particle Data Group.
    mean.Print();
    sigma.Print();
    meanG.Print();
    sigmaG.Print();
}