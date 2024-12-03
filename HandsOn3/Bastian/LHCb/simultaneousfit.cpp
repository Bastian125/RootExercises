#include <iostream>
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooBinning.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TAxis.h"
#include "TStyle.h"

using namespace RooFit;

void simultaneousfit()
{
    // C R E A T E  M O D E L
    //-----------------------
    RooRealVar x("x", "x", 5080, 5560);
    RooRealVar mean_B0("mean_B0", "mean of gaussian for B0", 5280, 5260, 5300);
    RooRealVar sigma("sigma_B0", "width of gaussian for B0", 20, 5, 50);
    RooRealVar mean_Bs("mean_Bs", "mean of gaussian for Bs", 5366, 5346, 5386);
    RooRealVar lambda("lambda", "lambda of exponential", -0.001, -1, -0.000001);

    RooGaussian sig_B0("sig_B0", "gaussian PDF for B0", x, mean_B0, sigma);
    RooGaussian sig_Bs("sig_Bs", "gaussian PDF for B0", x, mean_Bs, sigma);
    RooExponential exp("exp", "exponential PDF", x, lambda);

    RooRealVar B0_frac("B_frac", "fraction of B0", 0.1, 0, 1);
    RooRealVar BS_frac("bkg_frac", "fraction of BS", 0.1, 0, 1);

    RooAddPdf model("model", "full model", RooArgList(sig_B0, sig_Bs, exp), RooArgList(B0_frac, BS_frac));

    // Load data
    RooDataSet data = *RooDataSet::read("rarest_b0_decay.dat", x, "v");

    // C O N T R O L  R E G I O N
    // --------------------------

    // Define control region range
    RooRealVar y("y", "y", 4000, 5000);

    // Create exponential model
    RooExponential model_ctl("model_ctl", "Exponential control region model", y, lambda);

    // Generate data
    RooDataSet *data_ctl = model_ctl.generate(y, 10000);

    // Define category to distinguish physics and control samples events
    RooCategory sample("sample", "sample");
    sample.defineType("physics");
    sample.defineType("control");

    // Construct combined dataset in (x,sample)
    RooDataSet combData("combData", "combined data", RooArgSet(x, y), Index(sample), Import({{"physics", &data}, {"control", data_ctl}}));

    // C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
    // -----------------------------------------------------------------------------------

    // Construct a simultaneous pdf using category sample as index:
    // associate model with the physics state and model_ctl with the control state
    RooSimultaneous simPdf("simPdf", "simultaneous pdf", {{"physics", &model}, {"control", &model_ctl}}, sample);

    // P e r f o r m   a   s i m u l t a n e o u s   f i t
    // ---------------------------------------------------

    // Perform simultaneous fit of model to data and model_ctl to data_ctl
    std::unique_ptr<RooFitResult> fitResult{simPdf.fitTo(combData, PrintLevel(-1), Save(), PrintLevel(-1))};
    fitResult->Print();

    RooPlot *frame1 = x.frame(Title("Physics sample"));
    combData.plotOn(frame1, Cut("sample==sample::physics"));
    simPdf.plotOn(frame1, Slice(sample, "physics"), ProjWData(sample, combData));
    simPdf.plotOn(frame1, Slice(sample, "physics"), Components("exp"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame1, Slice(sample, "physics"), Components("sig_B0"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kRed));
    simPdf.plotOn(frame1, Slice(sample, "physics"), Components("sig_Bs"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kGray));

    RooPlot *frame2 = y.frame(Title("Control sample"));
    combData.plotOn(frame2, Cut("sample==sample::control"));
    simPdf.plotOn(frame2, Slice(sample, "control"), ProjWData(sample, combData));

    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    c->Divide(2);
    c->cd(1);
    frame1->Draw();
    c->cd(2);
    frame2->Draw();
    c->SaveAs("simultaneous_fit.pdf");
}