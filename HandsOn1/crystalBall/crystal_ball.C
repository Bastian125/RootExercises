 #include "RooRealVar.h"
 #include "RooDataSet.h"
 #include "RooCBShape.h"
 #include "TCanvas.h"
 #include "RooPlot.h"
 #include "TAxis.h"
 using namespace RooFit;
  
 void crystal_ball()
 {
    // S e t u p   m o d e l
    // ---------------------
  
    // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
    RooRealVar x("x", "x", -10, 10);
    RooRealVar mean("mean", "mean of crystal ball", 0, -10, 10);
    RooRealVar sigma("sigma", "width of crystal ball", 1, 0.1, 10);
    RooRealVar alpha("alpha", "alpha value of crystal ball", 1.5, 0, 10);
    RooRealVar n("n", "n of crystal ball", 1.5, 0, 10);
    
    // Build crystal ball pdf in terms of x,mean and sigma
    RooCBShape CB("CB", "crystal ball PDF", x, mean, sigma, alpha, n);
    
    // Construct plot frame in 'x'
    RooPlot *xframe = x.frame(Title("Crystal Ball pdf."));
    
    // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
    // ---------------------------------------------------------------------------
  
    // Plot CB in frame (i.e. in x)
    CB.plotOn(xframe);
    
    // Change the value of sigma to 0.3
    sigma.setVal(0.3);
    
    // Plot CB in frame (i.e. in x) and draw frame on canvas
    CB.plotOn(xframe, LineColor(kRed));
    
    // G e n e r a t e   e v e n t s
    // -----------------------------
  
    // Generate a dataset of 1000 events in x from CB
    std::unique_ptr<RooDataSet> data{CB.generate(x, 10000)};
    
    // Make a second plot frame in x and draw both the
    // data and the pdf in the frame
    RooPlot *xframe2 = x.frame(Title("Crystal ball pdf with data"));
    data->plotOn(xframe2);
    CB.plotOn(xframe2);

  
    // F i t   m o d e l   t o   d a t a
    // -----------------------------
  
    // Fit pdf to data
    CB.fitTo(*data, PrintLevel(-1));
    
    // Print values of mean and sigma (that now reflect fitted values and errors)
    mean.Print();
    sigma.Print();
    
    // Draw all frames on a canvas
    TCanvas *c = new TCanvas("crystal_ball", "crystal_ball", 800, 400); c->Divide(2);
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    xframe->GetYaxis()->SetTitleOffset(1.6);
    xframe->Draw();
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    xframe2->GetYaxis()->SetTitleOffset(1.6);
    xframe2->Draw();

    // Save pdf
    c->SaveAs("Crystal_Ball.pdf");
 }
 