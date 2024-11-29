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

void ExMinos()
{
	RooRealVar energy("energy", "Neutrino Energy [Gev]", 0.5,14);
	energy.setBins(25);
	RooDataSet* data = RooDataSet::read("minos_2013_data.dat", energy, "v");
	RooDataSet* mc_noosc = RooDataSet::read("minos_2013_mc.dat", energy, "v");
	
	
	RooDataSet* dd = (RooDataSet*) mc_noosc->reduce(RooArgSet(energy));
	RooDataHist* dh_mc_noosc = dd->binnedClone();
	RooHistFunc func_noosc{"func_mc_noosc", "No oscillation", energy, *dh_mc_noosc, 2};

	//RooRealVar energy("e", "Reconstructed Neutrino Energy [GeV]", 0.5, 14); 
	RooRealVar mixing("mixing", "sin squared of 2 theta", 0.8,0.5,1);
	RooRealVar dm2("dm2", "delta mass squared",2.5e-3,0, 5e-3);
	RooFormulaVar sin2_term("sin2_term", "sin(1.267*dm2*730/e)^2",
                        "sin(1.267*@0*730/@1)^2", RooArgList(dm2, energy));
	RooFormulaVar prob_osc("osc_prob", "Oscillation Probability",
                       "1 - @0 * @1", RooArgList(mixing, sin2_term));

	RooGenericPdf model("model", "model", "@0 * @1", RooArgSet(prob_osc, func_noosc));

	RooFitResult* result = model.fitTo(*data, RooFit::Save());

    // Plot data and model
    RooPlot* frame = energy.frame();
    data->plotOn(frame);              // Plot data points
    model.plotOn(frame);              // Plot model fit curve

    // Create a canvas to display the plot
    TCanvas canvas("canvas", "Fit to MINOS Data", 800, 600);
    frame->Draw();                    // Draw the plot on the canvas

    // Save the plot as a PNG file
	canvas.SaveAs("minos_data.png");

	// Part 2
	// creating the Negative likelihood object
	RooAbsReal * nll = model.createNLL(*data);
	// initializing the Minimizer
	RooMinimizer minim(*nll);
	// setting the output to Verbose
	minim.setVerbose(true);
	// calling the Minimizer 
	minim.migrad();
	//Printing
	//I believ this is done anyway when calling the Minimizer
	//mixing.Print();
	//dm2.Print();
	minim.setVerbose(false);
	// calculating errors from the second derivative at maximum
	minim.hesse();

	mixing.Print();
	dm2.Print();

	minim.minos(dm2);

	mixing.Print();
	dm2.Print();

	RooFitResult *fitResult = minim.save();
	minim.Print("v");


	TCanvas c("contour", "Contour Plot", 800, 600);
	RooPlot * contour = minim.contour(dm2, mixing, 2.30, 6.18, 11.83);
	contour->SetTitle("Contour Plot of dm2 vs mixing");
	contour->GetXaxis()->SetTitle("dm2");
	contour->GetYaxis()->SetTitle("mixing");
	contour->Draw();

	c.SaveAs("Minos_Contour.png");
	

	return ;
}
