#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooWorkspace.h"
#include "RooMinimizer.h"
using namespace RooFit;

void exercise11()
{
   RooWorkspace w{"w"};

   // S e t u p   m o d e l
   // ---------------------

   // Declare variables x,mean,sigma and f with associated name, title, initial value and allowed range
   w.factory("x[-20, 20]");
   w.factory("mean[0]");
   w.factory("s1[3.1]");
   w.factory("s2[4, 3, 6]");
   w.factory("f[0.5, 0.0, 1.0]");

   // Build gaussian pdf in terms of x,mean and sigma
   w.factory("Gaussian::g1(x, mean, s1)");
   w.factory("Gaussian::g2(x, mean, s2)");

   // Create composite model
   w.factory("SUM::model(f*g1, g2)");

   // Retrieve pointers to variables and PDFs for later use
   auto model = w.pdf("model");
   auto x = w.var("x");

   // Generate dataset with 1000 events
   auto data = model->generate(*x, 1000);
   w.import(*data);

   // Save workspace
   w.writeToFile("exercise11_workspace.root");

   // M i n i m i z e
   // ---------------

   // Construct function object representing negative logarithm of likelihood
   RooAbsReal *nll = model->createNLL(*data);

   // Minimize nll w.r.t. its parameters with Migrad
   RooMinimizer m(*nll);
   m.setVerbose(kTRUE); // enable verbose mode
   m.migrad();

   // Display the parameter values
   w.var("f")->Print();
   w.var("mean")->Print();
   w.var("s1")->Print();
   w.var("s2")->Print();

   // Disable verbose mode
   m.setVerbose(kFALSE);

   // Minimize nll w.r.t. its parameters with Hesse
   m.hesse();

   // Display the parameter values
   w.var("f")->Print();
   w.var("mean")->Print();
   w.var("s1")->Print();
   w.var("s2")->Print();

   // MINOS error calculation for s2
   m.minos(*w.var("s2"));

   // Display the parameter values
   w.var("f")->Print();
   w.var("mean")->Print();
   w.var("s1")->Print();
   w.var("s2")->Print();

   // Save and print result
   RooFitResult *fit_results = m.save();
   m.Print("v");

   // Visualize correlation matrix
   TCanvas* c1 = new TCanvas("c1", "Correlation Matrix", 800, 600);
   gStyle->SetPalette(1);
   fit_results->correlationHist()->Draw("colz");
   c1->SaveAs("CorrelationMatrix.pdf");


   // C o n t o u r  P l o t
   // ----------------------
   // Make contour plot of mx vs sx at 1,2,3 sigma
   TCanvas* c2 = new TCanvas("c2", "Contour Plot", 800, 600);
   RooPlot* contour = m.contour(*(w.var("f")), *(w.var("s2")), 2.30, 6.18, 11.83);
   contour->Draw();
   c2->SaveAs("ContourPlot.pdf");
}
