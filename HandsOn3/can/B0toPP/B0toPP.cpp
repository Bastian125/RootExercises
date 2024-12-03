#include <TApplication.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <TRandom.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
#include <RooHistFunc.h>
#include <RooGenericPdf.h>
#include <RooBinning.h>
#include <RooPoisson.h>
#include <RooProdPdf.h>
#include <RooDataSet.h>
#include <RooStats/ModelConfig.h>
#include <RooStats/ProfileLikelihoodCalculator.h>
#include "RooHist.h"
#include "TAxis.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

using namespace RooFit;

int main(int argc, char **argv)
{
    // /-1, -0.000001 );

    TCanvas *canvas = new TCanvas("canvas", "CB Fit", 2400, 600);
    canvas->Divide(3);
    canvas->cd(1);
    RooRealVar x("x", "x", 5080.0, 5560.0);

    Double_t binWidth = 20;
    Int_t nbins = (x.getMax() - x.getMin()) / binWidth;
    x.setBins(nbins);

    // RooRealVar meanB0("meanB0", "meanB0", 5279.65, 5270, 5290);
    // RooRealVar sigmaB0("sigmaB0", "sigmaB0", 20, 15, 25);
    // RooGaussian sigB0("sigB0", "sigB0", x, meanB0, sigmaB0);
    // RooRealVar meanB0s("meanB0s", "meanB0s", 5366, 5360, 5375);
    // RooRealVar sigmaB0s("sigmaB0s", "sigmaB0s", 7, 5, 9);
    // RooGaussian sigB0s("sigB0s", "sigB0s", x, meanB0s, sigmaB0s);
    // RooRealVar lambda("lambda", "lambda", 200.0,-1, -0.000001 );
    // // RooFormulaVar invTau("invTau", "-1/tau", RooArgList(tau));
    // RooExponential bkg("bkg", "bkg", x, lambda);
    // RooRealVar NBkg("NBkg", "NBkg", 11.5/20.0);
    // RooRealVar NsigB0("NsigB0", "NsigB0", 8.0/20.0);
    // RooRealVar NsigB0s("NsigB0s", "NsigB0s", 0.5/20.0);
    // RooAddPdf model("model", "model", RooArgList(bkg, sigB0, sigB0s), RooArgList(NBkg, NsigB0, NsigB0s));
    

    RooRealVar mean_B0("mean_B0", "mean of gaussian for B0", 5280, 5260,5300);
    RooRealVar sigma_B0("sigma_B0", "width of gaussian for B0", 20, 5, 50);
    RooRealVar mean_Bs("mean_Bs", "mean of gaussian for Bs", 5366, 5346,5386);
    RooRealVar sigma_Bs("sigma_Bs", "width of gaussian for Bs", 7, 5, 9);
    RooRealVar lambda("lambda", "lambda of exponential", -0.0001, -1, -0.000001);

    RooGaussian sig_B0("sig_B0", "gaussian PDF for B0", x, mean_B0, sigma_B0);
    RooGaussian sig_Bs("sig_Bs", "gaussian PDF for B0", x, mean_Bs, sigma_Bs);
    RooExponential exp("exp", "exponential PDF", x, lambda);

    RooRealVar B0_frac("B_frac", "fraction of B0", 0.1, 0, 1);
    RooRealVar BS_frac("bkg_frac", "fraction of BS", 0.1, 0, 1);

    RooAddPdf model("model", "model", RooArgList(sig_B0, sig_Bs, exp), RooArgList(B0_frac, BS_frac));




    RooDataSet data = *RooDataSet::read("rarest_b0_decay.dat", x, "v");

    // fit the model to the data using a maximum likelihood fit
    RooFitResult* result = model.fitTo(data, Save(true));

    RooPlot *frame = x.frame(Title("B0 to PP invariant mass"));
    data.plotOn(frame);
    
    model.plotOn(frame, Components("exp"), LineStyle(kDashed), LineColor(kGreen)); 
    model.plotOn(frame, Components("sig_B0"), LineStyle(kDashed), LineColor(kRed)); 
    model.plotOn(frame, Components("sig_Bs"), LineStyle(kDashed), LineColor(kGray)); 
    model.paramOn(frame, Format("NELU", AutoPrecision(2)), Layout(0.6, 0.95, 0.9));

    model.plotOn(frame);
    frame->Draw();
    
    RooHist* risiduals = frame->residHist();
    RooHist* pulls = frame->pullHist();

    canvas->cd(2);
    RooPlot* frame_risiduals= x.frame(Title("Residual Distribution"));
    frame_risiduals->addPlotable(risiduals, "P");
    frame_risiduals->Draw();

    canvas->cd(3);
    RooPlot* frame_pulls= x.frame(Title("Pull Distribution"));
    frame_pulls->addPlotable(pulls, "P");
    frame_pulls->Draw();


    
    // Draw the frame
    
    canvas->Update();

    canvas->SaveAs("invMass.pdf");

    // Clean up
    delete frame;
    delete canvas;

    TCanvas c;
    gStyle->SetPalette(1);
    result->correlationHist()->Draw("colz");

    c.SaveAs("correlation.pdf");



    TCanvas *canvas_sim= new TCanvas("canvas", "Simultaneous Fit", 1200, 600);
    canvas_sim->Divide(2);

    RooRealVar x_control("x_control", "x_control", 4000.0, 5000.0);


    RooExponential model_control("bkg_control", "bkg_control", x_control, lambda);
    lambda.setVal(-1.0e-3);
    //generate data
    RooDataSet *data_control = model_control.generate(x_control, 10000);

    RooCategory sample("sample", "sample");
    sample.defineType("physics");
    sample.defineType("control");
    
    RooDataSet combData("combData","combined data", RooArgSet(x, x_control), Index(sample), 
                                                Import("physics", data), Import("control", *data_control));

    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
    simPdf.addPdf(model, "physics");
    simPdf.addPdf(model_control, "control");
    

    std::unique_ptr<RooFitResult> fitResult{simPdf.fitTo(combData, PrintLevel(-1), Save(), PrintLevel(-1))};
    fitResult->Print();

    canvas_sim->cd(1);
    RooPlot *frame_physics = x.frame(Title("B0 to PP invariant mass"));
    combData.plotOn(frame_physics, Cut("sample==sample::physics"));
    simPdf.plotOn(frame_physics, Slice(sample, "physics"), ProjWData(sample, combData), LineColor(kRed));
    simPdf.plotOn(frame_physics, Slice(sample, "physics"), ProjWData(sample, combData), Components("sig_B0"), LineColor(kGreen));
    simPdf.plotOn(frame_physics, Slice(sample, "physics"), ProjWData(sample, combData), Components("sig_Bs"), LineColor(kGray));
    simPdf.plotOn(frame_physics, Slice(sample, "physics"), ProjWData(sample, combData), Components("exp"), LineColor(kBlue));
    simPdf.paramOn(frame_physics, Format("NELU", AutoPrecision(2)), Layout(0.6, 0.95, 0.9));

    frame_physics->Draw();

    canvas_sim->cd(2);
    RooPlot *frame_control = x_control.frame(Title("B0 to PP invariant mass"));
    combData.plotOn(frame_control, Cut("sample==sample::control"));
    simPdf.plotOn(frame_control, Slice(sample, "control"), ProjWData(sample, combData), LineColor(kRed));
    frame_control->Draw();

    canvas_sim->SaveAs("simultaneous.pdf");

    return 0;
}
