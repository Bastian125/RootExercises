#include <iostream>
#include "RooRealVar.h"
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
#include "TH2.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
using namespace RooFit;

void B0_to_pp_TheosEIGENTUM()
 {

    // Part 1 //

    RooRealVar x("x", "x", 5080, 5560);
    RooRealVar mean_B0("mean_B0", "mean of gaussian for B0", 5280, 5260,5300);
    RooRealVar sigma("sigma", "width of gaussian", 20, 5, 50);
    RooRealVar mean_Bs("mean_Bs", "mean of gaussian for Bs", 5366, 5346,5386);
    RooRealVar lambda("lambda", "lambda of exponential", -0.0001, -1, -0.000001);

    RooGaussian sig_B0("sig_B0", "gaussian PDF for B0", x, mean_B0, sigma);
    RooGaussian sig_Bs("sig_Bs", "gaussian PDF for B0", x, mean_Bs, sigma);
    RooExponential exp("exp", "exponential PDF", x, lambda);

    RooRealVar B0_frac("B_frac", "fraction of B0", 0.1, 0, 1);
    RooRealVar BS_frac("bkg_frac", "fraction of BS", 0.1, 0, 1);


    RooAddPdf model("model", "full model", RooArgList(sig_B0, sig_Bs, exp), RooArgList(B0_frac, BS_frac));

    RooDataSet data = *RooDataSet::read("data/rarest_b0_decay.dat", x, "v");

    float bin_width = 20;
    int nbins = (x.getMax() - x.getMin()) / bin_width;
    x.setBinning(RooBinning(nbins, x.getMin(), x.getMax()));
    RooFitResult* fit_results = model.fitTo(data, PrintLevel(-1), Save());

    mean_B0.Print();
    mean_Bs.Print();
    sigma.Print();
    lambda.Print();

    // Part 2 //

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

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    gStyle->SetPalette(1) ;
    fit_results->correlationHist()->Draw("colz");
    c3->SaveAs("output/correlation.pdf");

    // Part 3 //

    mean_B0.setVal(5280); // reset all params
    sigma.setVal(20);
    mean_Bs.setVal(5366);
    lambda.setVal(-0.0001);
    B0_frac.setVal(0.1);
    BS_frac.setVal(0.1);

    RooRealVar y("y", "y", 4000, 5000);
    RooExponential model_ctl("exp_c", "control exponential PDF", y, lambda);
    lambda.setVal(-1.0e-3);

    RooDataSet* data_ctl = model_ctl.generate(y, 10000);

    RooCategory sample("sample", "sample");
    sample.defineType("physics");
    sample.defineType("control");

    RooDataSet combData("combData","combined data", RooArgSet(x,y), Index(sample), Import({{"physics", &data}, {"control", data_ctl}}) );
    RooSimultaneous simPdf("simPdf", "simultaneous pdf", {{"physics", &model}, {"control", &model_ctl}}, sample);
    RooFitResult* fitResult{simPdf.fitTo(combData, PrintLevel(-1), Save(), PrintLevel(-1))};
    fitResult->Print();

    RooPlot *frame5 = x.frame(Title("Physics sample"));
    combData.plotOn(frame5, Cut("sample==sample::physics"));
    simPdf.plotOn(frame5, Slice(sample, "physics"), ProjWData(sample, combData));
    simPdf.plotOn(frame5, Slice(sample, "physics"), Components("exp"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame5, Slice(sample, "physics"), Components("sig_B0"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kRed));
    simPdf.plotOn(frame5, Slice(sample, "physics"), Components("sig_Bs"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kGray));

    RooPlot *frame6 = y.frame(Title("Control sample"));
    combData.plotOn(frame6, Cut("sample==sample::control"));
    simPdf.plotOn(frame6, Slice(sample, "control"), ProjWData(sample, combData));

    TCanvas *c4 = new TCanvas("c4", "c4", 1200, 600);
    c4->Divide(2);
    c4->cd(1);
    frame5->Draw();
    c4->cd(2);
    frame6->Draw();
    c4->SaveAs("output/simultaneous_fit.pdf");
}
