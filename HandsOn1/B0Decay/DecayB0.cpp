#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"

using namespace RooFit;

void DecayB0() {
	RooRealVar x("x","x",5080,5550);
	x.setBins(40);
	RooDataSet data = *RooDataSet::read("rarest_b0_decay.dat", x, "v");
    // Gaussian for B0 
    RooRealVar meanB0("meanB0", "mean of B0 gaussian", 5280, 5000,5600);
    RooRealVar sigmaB0("sigmaB0", "width of B0 gaussian", 6, 0.1, 1000);
    RooGaussian gaussianB0("gaussianB0", "B0 Gaussian", x, meanB0, sigmaB0);
	// Gaussian for Bs
    RooRealVar meanBs("meanBs", "mean of Bs gaussian", 5366, 5000,5600);
    RooRealVar sigmaBs("sigmaBs", "width of Bs gaussian", 5, 0.1, 1000);
    RooGaussian gaussianBs("gaussianBs", "Bs Gaussian", x, meanBs, sigmaBs);

    // Exponential background: tau = 10
    RooRealVar tau("tau", "exponential background tau", 1, 0.01, 20);
	//RooFormulaVar("gen", "exp(-1/tau)", tau);
	RooFormulaVar invnegtau("gen", "-1/tau",tau);
    RooExponential bkgEx("bkgEx", "Exponential background", x, invnegtau);
	RooPolynomial bkgFl("bkgFl", "Flat Background", x);
	
	RooRealVar NB0("NB0", "Yield of B0 signal", 35, 0, 1000);
    RooRealVar NBs("NBs", "Yield of Bs signal", 5, 0, 1000);
    RooRealVar Nbkg("Nbkg", "Yield of Background", 100, 0, 1000);

    // Create the composite model: signal + background
    RooAddPdf model("model", "B0 + Bs + Background",
			RooArgList(gaussianB0, gaussianBs, bkgFl),
			RooArgList(NB0, NBs, Nbkg));//, RooArgList(fsig));

    // Create a frame to plot the model and the data
    RooPlot* frame = x.frame(Title("Signal + Background Model"));

    // Plot the data
    data.plotOn(frame);
	
	model.fitTo(data);
    // Plot the model
    model.plotOn(frame, Name("Model"));
	RooHist *hresid = frame->residHist();
    RooHist *hpull = frame->pullHist();

	model.plotOn(frame, Components(gaussianB0), LineColor(kRed));
	model.plotOn(frame, Components(gaussianBs), LineColor(kBlack));
	model.plotOn(frame, Components(bkgFl), LineColor(kOrange));

    // Create a canvas and draw the frame
    TCanvas* c = new TCanvas("c", "Signal + Background Fit", 2000, 600);
	c->Divide(3);
    RooPlot *frame2 = x.frame(Title("Residual Distribution"));
    RooPlot *frame3 = x.frame(Title("Residual Distribution"));
    frame2->addPlotable(hresid, "P");
    frame3->addPlotable(hpull, "P");
	c->cd(1);
    gPad->SetLeftMargin(0.15);
	frame->Draw();
	c->cd(2);
    gPad->SetLeftMargin(0.15);
	frame2->Draw();
	c->cd(3);
    gPad->SetLeftMargin(0.15);
	frame3->Draw();
	//	TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
	//legend->SetBorderSize(0); // No border
    //legend->SetFillStyle(0);  // Transparent background
    //legend->AddEntry(frame->findObject("Model"), "True Model", "l");
    //legend->AddEntry(frame->findObject("Signal"), "True Signal (Gaussian)", "l");
    //legend->AddEntry(frame->findObject("Background"), "True Background (Exponential)", "l");
    //frame->Draw();

    // Save the plot
    //c->SaveAs("signal_background_fit.png");
    
    // Perform the fit (Maximum Likelihood Fit)
    //model.fitTo(*data);
	//model.plotOn(frame,LineColor(kOrange), Name("Fitted Model"));
	//legend->AddEntry(frame->findObject("Fitted Model"), "Model fitted to Data", "l");
    // Draw the fitted model
    //frame->Draw();
	//legend->Draw();

    // Save the fit result plot
    c->SaveAs("B0_decay_data.pdf");
}

