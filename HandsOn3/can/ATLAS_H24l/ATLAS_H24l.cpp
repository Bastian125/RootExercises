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
#include <RooHist.h>
#include <TAxis.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <RooPolynomial.h>
#include <RooCBShape.h>
#include <RooStats/FeldmanCousins.h>
#include </home/can/root/tutorials/roostats/StandardProfileLikelihoodDemo.C>
#include </home/can/root/tutorials/roostats/StandardFeldmanCousinsDemo.C>
#include </home/can/root/tutorials/roostats/StandardHypoTestDemo.C>

using namespace RooFit;
using namespace RooStats;
using namespace std;

int main(int argc, char **argv)
{
    RooRealVar x{"x","invariant mass", 110, 135, "GeV"};
    x.setBins(10);
    RooPlot *frame = x.frame(Title("Invariant mass of Higgs boson"));

    TCanvas *canvas = new TCanvas("canvas", "CB Fit", 800, 600);

    RooDataSet *data = RooDataSet::read("higgs_4l.dat", RooArgList(x));
    
    // variables for the BreitWiegner and polynomial
    RooRealVar mass("mass", "mass", 125, 110, 150);
    RooRealVar width("width", "width", 4.1/2.35);
    RooRealVar alpha("alpha", "alpha", 0.6);
    RooRealVar n("n", "n", 20);
    RooRealVar a0("a0", "a0", -160, -100, -200);
    RooRealVar a1("a1", "a1", 2.7, 2, 4);


    RooRealVar nsignal("nsignal", "number of signal events", 100, 0, 1000);
    RooRealVar nbackground("nbackground", "number of background events", 100, 0, 1000);
    
    

    RooPolynomial bmodel("poly", "polynomial PDF", x, RooArgList(a0, a1));
    RooCBShape smodel("CB", "CB PDF", x, mass, width, alpha, n);

    RooAddPdf model("model", "model", RooArgList(smodel, bmodel), RooArgList(nsignal, nbackground));

    RooFitResult *result = model.fitTo(*data, Save(true));

    data->plotOn(frame);
    model.plotOn(frame);
    frame->Draw();
    canvas->SaveAs("invMass.pdf");

    RooWorkspace w = RooWorkspace("w", "workspace");
    w.import(model);
    w.import(*data);

    // create the model configuration with nSig as the signal parameter of interest
    ModelConfig mc("ModelConfig", &w);

    w.var("width")->setConstant(true);
    w.var("alpha")->setConstant(true);
    w.var("n")->setConstant(true);

    mc.SetPdf(model);
    mc.SetParametersOfInterest(nsignal);
    mc.SetObservables(x);
    w.defineSet("nuisParams", RooArgSet(nbackground,a0,a1));
    mc.SetNuisanceParameters(*w.set("nuisParams"));
    mc.SetSnapshot(nsignal);

    w.var("mass")->setConstant(true);

    w.import(mc);

    ModelConfig mc_mass("ModelConfig_mass", &w);
    mc_mass.SetPdf(model);
    mc_mass.SetParametersOfInterest(mass);
    mc_mass.SetObservables(x);
    w.defineSet("nuisParams_mass", RooArgSet(nbackground,a0,a1,nsignal,width));
    mc_mass.SetNuisanceParameters(*w.set("nuisParams_mass"));
    mc_mass.SetSnapshot(mass);

    w.import(mc_mass);

    w.writeToFile("workspace.root");

    TCanvas *canvas2 = new TCanvas("canvas2", "Likeliehood plt nSig", 800, 600);

    ProfileLikelihoodCalculator plc(*data, mc);
    plc.SetConfidenceLevel(0.68);
    LikelihoodInterval* interval = plc.GetInterval();
    LikelihoodIntervalPlot plot(interval);
    plot.SetRange(40, 100);
    plot.Draw();
    canvas2->SaveAs("profileLikelihood.pdf");



    RooDataSet* data_mass = (RooDataSet*)data->reduce(RooArgSet(x));
    FeldmanCousins fc(*data_mass, mc_mass);

    double cl = 0.90;

    fc.UseAdaptiveSampling(true);
    fc.FluctuateNumDataEntries(false);
    fc.SetNBins(100); 
    fc.SetConfidenceLevel(cl);
    // get interval
    PointSetInterval *interval2 = (PointSetInterval*)fc.GetInterval();





    // print interval
    interval2->SetConfidenceLevel(cl);
    RooRealVar *firstPOI = (RooRealVar *)mc_mass.GetParametersOfInterest()->first();
    cout << cl*100 << "% interval on " << firstPOI->GetName() << " is : [" << interval2->LowerLimit(*firstPOI) << ", " << interval2->UpperLimit(*firstPOI) << "] " << endl;


    // Hypothesis Test


    

    RooRealVar mu("mu", "mu", 1, 0, 1);
    w.import(mu);
    w.var("mu")->setConstant(true);
    RooFormulaVar nsig_scaled("nsig_scaled", "nsignal * mu", RooArgList(nsignal, mu));
    RooAddPdf model_sb("model_sb", "model_sb", RooArgList(smodel, bmodel), RooArgList(nsig_scaled, nbackground));
    RooFitResult* result_sb = model_sb.fitTo(*data, Save(true));
    result_sb->Print("v");

    w.import(model_sb);

    // //plot and safe the fit
    // TCanvas *canvas_sb = new TCanvas("canvas_sb", "Higgs stff", 800, 600);
    // RooPlot* frame_sb = x.frame();
    // data->plotOn(frame_sb);
    // model_sb.plotOn(frame_sb);
    // frame_sb->Draw();
    // canvas->SaveAs("invMass_sb.pdf");


    w.var("mu")->setConstant(false);

    ModelConfig mc_sb("ModelConfig_sb", &w);
    mc_sb.SetPdf(model_sb);
    mc_sb.SetParametersOfInterest(mu);
    mc_sb.SetObservables(x);
    w.defineSet("nuisParams_sb", RooArgSet(nsignal,nbackground,a0,a1));
    mc_sb.SetNuisanceParameters(*w.set("nuisParams_sb"));
    mc_sb.SetSnapshot(mu);

    ProfileLikelihoodCalculator plc_sb(*data, mc_sb);



    w.var("mu")->setVal(0);
    plc_sb.SetNullParameters(*w.var("mu"));
    HypoTestResult* hypotest = plc_sb.GetHypoTest();
    double alpha_value = hypotest->NullPValue();
    double significance_value = hypotest->Significance();

    cout << "p-value: " << alpha_value << endl;
    cout << "significance: " << significance_value << endl;


    return 0;
}
