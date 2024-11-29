#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "TFile.h"


using namespace RooFit;

void ATLAS() {
	RooRealVar x("x", "invariant mass", 110, 135,"GeV");
	x.setBins(10);
	RooDataSet *data = RooDataSet::read("higgs_4l.dat", x, "v");
	

	// Signal model (Crystal Ball)
	RooRealVar mass("mass", "Higgs Mass -> mean of CB", 125, 110, 150); // Mean of the CB
	RooRealVar width("width", "Higgs Width -> width of CB", 4.1 / 2.35); // Fixed width
	RooRealVar alpha("alpha", "Alpha", 0.6);                // Fixed alpha
	RooRealVar n("n", "Exponent", 20);                     // Fixed n
	
	RooCBShape smodel("smodel", "Signal Model (Crystal Ball)", x, mass, width, alpha, n);
	
	// Background model (Polynomial of degree 2)
	RooRealVar a1("a1", "Coefficient a1", -160, -100, -200); // Background slope parameter
	RooRealVar a2("a2", "Coefficient a2", 2.7, 2, 4);     // Background curvature parameter
	
	RooPolynomial bmodel("bmodel", "Background Model", x, RooArgList(a1, a2));

	RooRealVar nsignal("nsignal", "Number of Signal Events", 50,0,123);	
	RooRealVar nbackground("nbackground", "Number of background Events", 50,0,123);	

	// Combined model
	RooAddPdf model("model", "Signal + Background Model", RooArgList(smodel, bmodel), RooArgList(nsignal, nbackground));

	
	TCanvas* c = new TCanvas("c", "Data + Fit", 800, 600);
	// Plot the data and model fit
	RooPlot* frame = x.frame(Title("Signal + Background Model"));
	
	// Plot the data
	data->plotOn(frame);
	
	RooFitResult* fit_results = model.fitTo(*data, RooFit::Save());	
	model.plotOn(frame, Name("Total Model"), LineColor(kRed));
	//model.plotOn(frame, Components(smodel), LineColor(kBlue));
	model.plotOn(frame, Components(bmodel),LineColor(kRed), FillColor(kRed),FillStyle(3001),DrawOption("F"));
	
	frame->Draw();
	c->SaveAs("ATLAS_Data_Fit.png");


	// Create the workspace
	RooWorkspace w("w", "workspace");


	// Import dataset and model
	w.import(*data);   // Import the dataset
	w.import(model);   // Import the combined model

	w.var("width")->setConstant(true);
	w.var("alpha")->setConstant(true);
	w.var("n")->setConstant(true);
	w.var("mass")->setConstant(true);

	// CONFIG 1, analysis for nsignal
	RooStats::ModelConfig mc("ModelConfig", &w);

	// Set the PDF for the model
	mc.SetPdf(model);
	// Set the parameters of interest (e.g., signal yield)
	w.var("nsignal")->Print();
	mc.SetParametersOfInterest(RooArgSet(*w.var("nsignal")));
	
	// Set the observables (e.g., the mass variable)
	mc.SetObservables(RooArgSet(mass));
	

	// Define parameters and set for nuisance parameters
	w.defineSet("nuisParams", "nbackground,a1,a2");
	
	// Ensure the variables are in the workspace
	w.var("nbackground"); // Check if nbackground exists
	w.var("a1");          // Check if a1 exists
	w.var("a2");          // Check if a2 exists
	
	// Now you can set the nuisance parameters in ModelConfig
	mc.SetNuisanceParameters(RooArgSet(*w.set("nuisParams")));

	
	// Define the set of nuisance parameters (e.g., background parameters)
	//w.defineSet("nuisParams", "nbackground, a1, a2");  // Background parameters as nuisance
	//mc.SetNuisanceParameters(*w.set("nuisParams"));

	mc.SetSnapshot(*w.var("nsignal"));

	w.import(mc);


	w.var("mass")->setConstant(false);
	w.var("nsignal")->setConstant(true);

	// Define another ModelConfig for the Higgs mass
	w.defineSet("nuisParams2", "nbackground,a1,a2,nsignal,width");  // Include nsignal, width in nuisance parameters


	RooStats::ModelConfig mc_mass("ModelConfig_mass", &w);
	mc_mass.SetPdf(model);  // Combined model
	mc_mass.SetParametersOfInterest(RooArgSet(*w.var("mass")));  // Parameter of interest is mass
	mc_mass.SetObservables(RooArgSet(*w.var("mass")));  // Observable is mass
	mc_mass.SetNuisanceParameters(*w.set("nuisParams2"));  // Define nuisance parameters	
	
	mc_mass.SetSnapshot(RooArgSet(*w.var("mass")));

	w.import(mc_mass);

	TFile* file = new TFile("workspace.root", "RECREATE");
	w.Write();  // Write the workspace to the file
	file->Close();


	return;
}

