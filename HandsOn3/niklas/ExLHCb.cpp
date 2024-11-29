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

void ExLHCb() {
	RooRealVar m("m","m(pp)[MeV/C^2]",5080,5550);
	//m.setBins(40);
	RooDataSet *data = RooDataSet::read("rarest_b0_decay.dat", m, "v");


	// Define exponential background pdf 
	RooRealVar tau("tau", "tau", -0.001, -1.0, 0.);
	RooExponential bkg("bkg", "Exponetial Background", m, tau);

    // Gaussian for B0 
    RooRealVar meanB0("meanB0", "mean of B0 gaussian", 5280, 5000,5600);
    RooRealVar sigmaB0("sigmaB0", "width of B0 gaussian", 20, 15, 25);
    RooGaussian gaussianB0("gaussianB0", "B0 Gaussian", m, meanB0, sigmaB0);
	RooRealVar fB0("fB0", "signal fraction B0", 0.5, 0.0,1.0);
	// Gaussian for Bs
    RooRealVar meanBs("meanBs", "mean of Bs gaussian", 5366, 5000,5600);
    RooRealVar sigmaBs("sigmaBs", "width of Bs gaussian", 7, 5, 9);
    RooGaussian gaussianBs("gaussianBs", "Bs Gaussian", m, meanBs, sigmaBs);
	RooRealVar fBs("fBs", "signal fraction Bs", 0.5, 0.0,1.0);

	//RooRealVar NB0("NB0", "Yield of B0 signal", 35, 0, 1000);
    //RooRealVar NBs("NBs", "Yield of Bs signal", 5, 0, 1000);
    //RooRealVar Nbkg("Nbkg", "Yield of Background", 100, 0, 1000);
	


    // Create the composite model: signal + background
    RooAddPdf model("model", "B0 + Bs + Background",
			RooArgList(gaussianB0, gaussianBs, bkg),
			RooArgList(fB0,fBs));
			//RooArgList(NB0, NBs, Nbkg));//, RooArgList(fsig));
			//
	
	// Plot the data and model fit
	RooPlot* frame = m.frame(Title("Signal + Background Model"));
	
	// Plot the data
	data->plotOn(frame);
	
	// Fit the model to the data
	RooFitResult* fit_results = model.fitTo(*data, RooFit::Save());	
	// Create a canvas to draw the correlation matrix
	TCanvas* c = new TCanvas("c", "Correlation Matrix", 800, 600);
	
	// Set the color palette for better visualization
	gStyle->SetPalette(1);
	
	// Draw the correlation matrix as a histogram
	fit_results->correlationHist()->Draw("colz");
	
	// Save the correlation matrix plot as a PNG file
	c->SaveAs("CorrelationMatrix.png");
	// Plot the model
	model.plotOn(frame, Name("Model"));
	model.plotOn(frame, Components(gaussianB0), LineColor(kRed));
	model.plotOn(frame, Components(gaussianBs), LineColor(kBlack));
	model.plotOn(frame, Components(bkg), LineColor(kOrange));
	
	// Create residual and pull histograms
	RooHist* hresid = frame->residHist();
	RooHist* hpull = frame->pullHist();
	
	RooPlot* frame2 = m.frame(Title("Residual Distribution"));
	RooPlot* frame3 = m.frame(Title("Pull Distribution"));
	
	frame2->addPlotable(hresid, "P");
	frame3->addPlotable(hpull, "P");
	
	// First canvas: Data and fit
	TCanvas* c1 = new TCanvas("c1", "Data and Fit", 800, 600);
	c1->cd();
	gPad->SetLeftMargin(0.15);
	frame->Draw();
	c1->SaveAs("Data_Fit.png"); // Save as PNG
	
	// Second canvas: Residuals and Pull distributions
	TCanvas* c2 = new TCanvas("c2", "Residuals and Pulls", 1200, 600);
	c2->Divide(2);
	
	c2->cd(1);
	gPad->SetLeftMargin(0.15);
	frame2->Draw();
	
	c2->cd(2);
	gPad->SetLeftMargin(0.15);
	frame3->Draw();
	
	c2->SaveAs("Residuals_Pulls.png"); // Save as PNG
		
	// Part 3
	//
	//

	RooRealVar y("y", "Control Region Observable", 4000, 5000);
	tau.setVal(-1.0e-3);
	//RooRealVar tau("tau", "Exponential Coefficient", -1.0e-3);
	
	RooExponential model_ctl("model_ctl", "Control Region Model", y, tau);
	
	// Generate data for control region
	RooDataSet* data_ctl = model_ctl.generate(y, 10000);
	// Define the Category to Distinguish Events
	RooCategory sample("sample", "Sample Type");
	sample.defineType("physics");
	sample.defineType("control");
	//Combine Events
	RooDataSet combData("combData", "Combined Data", RooArgSet(m, y), Index(sample),
                    Import("physics", *data), Import("control", *data_ctl));
	//Define the Simultaneous Model
	RooSimultaneous simPdf("simPdf", "Simultaneous PDF", sample);
	simPdf.addPdf(model, "physics");
	simPdf.addPdf(model_ctl, "control");
	//Fit the Simultaneous Model
	simPdf.fitTo(combData);
	//Plot the Results

	// Plot for physics region
	RooPlot* frame_phys = m.frame(Title("Physics Region"));
	combData.plotOn(frame_phys, Cut("sample==sample::physics"));
	simPdf.plotOn(frame_phys, Slice(sample, "physics"), ProjWData(sample, combData));
	TCanvas* canvas_phys = new TCanvas("canvas_phys", "Physics Region", 800, 600);
	frame_phys->Draw();
	canvas_phys->SaveAs("PhysicsRegion.png");
	
	// Plot for control region
	RooPlot* frame_ctl = y.frame(Title("Control Region"));
	combData.plotOn(frame_ctl, Cut("sample==sample::control"));
	simPdf.plotOn(frame_ctl, Slice(sample, "control"), ProjWData(sample, combData));
	TCanvas* canvas_ctl = new TCanvas("canvas_ctl", "Control Region", 800, 600);
	frame_ctl->Draw();
	canvas_ctl->SaveAs("ControlRegion.png");
	



	return;
}

