#include <iostream>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooWorkspace.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooHistFunc.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "TH2.h"
#include "TStyle.h"
#include "RooBinning.h"
using namespace RooFit;

void minos()
 {
   RooRealVar energy("x", "x", 0.5, 14);
   float bin_width = 0.5;
   int nbins = (energy.getMax() - energy.getMin()) / bin_width;
   energy.setBinning(RooBinning(nbins, energy.getMin(), energy.getMax()));
   RooDataSet data = *RooDataSet::read("data/minos_2013_data.dat", energy, "v");

   RooDataSet mc_noosc = *RooDataSet::read("data/minos_2013_mc.dat", energy, "v");

   RooDataSet* dd = (RooDataSet*) mc_noosc.reduce(RooArgSet(energy)) ;
   RooDataHist* dh_mc_noosc = dd->binnedClone();
   RooHistFunc func_noosc { "func_mc_noosc", "No oscillation", energy, *dh_mc_noosc, 2 };

   RooRealVar mixing("sin2theta", "sin2theta", 0.5, 0 ,1);
   RooRealVar dm2("deltam", "deltam", 2.5e-3, 1e-4, 1e-2);
   RooRealVar L("L", "Baseline [km]", 730);
   L.setConstant(kTRUE);

   RooFormulaVar arg("arg", "1.267 * dm2 * L / E", "1.267 * @0 * @1 / @2", RooArgList(dm2, L, energy));
   RooFormulaVar sin2_term("sin2_term", "sin^2(x)", "sin(@0)^2", RooArgList(arg));
   RooFormulaVar prob_osc("P_survival", "P(mu->mu)", "1 - @0 * @1", RooArgList(mixing, sin2_term));
   RooGenericPdf model = RooGenericPdf{ "model", "model", "@0*@1", RooArgSet(func_noosc, prob_osc) };
   model.fitTo(data, PrintLevel(-1));

   RooPlot *frame1 = energy.frame(Title("Minos Fit"));
   data.plotOn(frame1);
   model.plotOn(frame1, Binning(nbins, energy.getMin(), energy.getMax()), LineStyle(2), LineColor(kBlue));
   TCanvas* c1 = new TCanvas("c1", "Minos Fit", 800, 600);
   frame1->Draw();
   c1->SaveAs("output/minos_data.png");


   RooAbsReal* NLL = model.createNLL(data);
   RooMinimizer m(*NLL);
   m.setVerbose(kTRUE);
   m.migrad();

   std::cout << "\n";
   std::cout << "Task [A]:";
   std::cout << "\n";
   mixing.Print();
   dm2.Print();

   m.setVerbose(kFALSE);
   m.hesse();

   std::cout << "\n";
   std::cout << "Task [B]:";
   std::cout << "\n";
   mixing.Print();
   dm2.Print();

   m.minos(dm2);
   std::cout << "\n";
   std::cout << "Task [C]:";
   std::cout << "\n";
   dm2.Print();

   RooFitResult* fit_result = m.save();
   std::cout << "\n";
   std::cout << "Task [D]:";
   std::cout << "\n";
   fit_result->Print("v");

   TCanvas* c2 = new TCanvas("c2", "Contour", 800, 800);
   RooPlot* contour = m.contour(mixing, dm2, 2.30,6.18,11.83);
   contour->GetXaxis()->SetRangeUser(0.6, 1);  // Set X-axis range from 2 to 8
   contour->GetYaxis()->SetRangeUser(0.0015, 0.0035);  // Set Y-axis range from 1 to 9
   contour->Draw();
   c2->SaveAs("output/minos_likelihood.png");
 }
