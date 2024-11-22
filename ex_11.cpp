#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "TAxis.h"
using namespace RooFit;

void ex_11()
{
    RooWorkspace w("w");
    w.factory("x[   -20,20]");
    w.factory("mean[0]");
    // w.factory("s1[3.3]"); // last digit Matrikelnummer Unibo: 3
    w.factory("s1[3.9]"); // TEST
    // w.factory("s1[3.3,3.3,3.3]"); // TODO
    w.factory("s2[4,3,6]");
    w.factory("Gaussian::gauss1(x,mean,s1)");
    w.factory("Gaussian::gauss2(x,mean,s2)");
    // ↓ geht
    w.factory("SUM:model(f[0.5,0,1] * gauss1, gauss2)");
    // ↓ geht nicht
    // w.factory("f[0.5,0,1]");
    // w.factory("SUM:model(f * gauss1, (1-f) * gauss2)");
    // ↑
    RooAbsPdf *model = w.pdf("model");
    RooRealVar *x = w.var("x");

    // 1. Generate a dataset with 1000 events
    RooDataHist *data = model->generateBinned(*x, 1'000);
    data->SetName("data"); // TODO: TEST

    // 2. Save the workspace (model + data) to a file
    w.import(*data);
    w.writeToFile("ex_11.root");

    // 3. Minimize the likelihood…
    // NOTE: RooMinuit is deprecated, use RooMinimizer instead
    // a)
    RooAbsReal *nll = model->createNLL(*data);
    // b)
    RooMinimizer minim(*nll);
    // c)
    minim.setVerbose(true);
    // d)
    minim.migrad();
    // e)
    w.var("f")->Print();
    w.var("mean")->Print();
    w.var("s1")->Print();
    w.var("s2")->Print();

    // 4. HESSE Error Calculation
    // a)
    minim.setVerbose(false);
    // b)
    minim.hesse();
    // c)
    w.var("f")->Print();
    w.var("mean")->Print();
    w.var("s1")->Print();
    w.var("s2")->Print();

    // 5. MINOS Error Calculation for "s2"
    // a)
    minim.minos(*w.var("s2"));
    // b)
    w.var("f")->Print();
    w.var("mean")->Print();
    w.var("s1")->Print();
    w.var("s2")->Print();

    // 6.
    // a)
    RooFitResult *fit_result = minim.save();
    // b)
    fit_result->Print("v"); // v: verbose
    // c)
    gStyle->SetPalette(1);
    fit_result->correlationHist()->Draw("colz");

    // 7. Contour Plot
    // a)
    RooPlot *contour_plot = minim.contour(*w.var("f"), *w.var("s2"), 1, 2, 3);
    // b)
    // Draw the plot frame (using the standard Draw() method) and save it to a file.
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    contour_plot->Draw();
}
