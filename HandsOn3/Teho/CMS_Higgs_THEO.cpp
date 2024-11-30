#include <iostream>
#include <fstream>

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooProduct.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/HypoTestResult.h"
using namespace std;
using namespace RooFit;
using namespace RooStats;

RooDataHist read_file(const char* filename, RooRealVar& x)
{
    RooDataHist data{"data", "data", x};
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }
    double val, weight;
    vector<pair<double, double>> entries;
    while (file >> val >> weight)
    {
        entries.emplace_back(val, weight);
    }
    file.close();

    // Skip the first and last entries (overflow)
    for (size_t i = 1; i < entries.size() - 1; ++i)
    {
        x.setVal(entries[i].first);
        data.add(x, entries[i].second);
    }

    x.setBins(entries.size() - 2);

    return data;
}

//RooDataHist read_file(const char* filename, RooRealVar&x)
//{
//    RooDataHist data{"data", "data", x};
//    std::ifstream file(filename);
//    double val, weight;
//    while (!file.eof()) {
//        file >> val >> weight;
//        x.setVal(val);
//        data.set(x, weight);
//    }
//    return data;
//
//}

void CMS_Higgs()
 {
    // Part 1 //

    RooRealVar x("x", "x", 71.5, 182.5);
    x.setUnit("GeV");
    float bin_width = 3;
    int nbins = (x.getMax() - x.getMin()) / bin_width;
    x.setBins(nbins);

    RooDataHist data = read_file("data/_cms_higgs_data.txt",x);
    RooDataHist DYto4l = read_file("data/_cms_higgs_DYto4l.txt",x);
    RooDataHist TTbarto4l = read_file("data/_cms_higgs_TTbarto4l.txt",x);
    RooDataHist ZZto4l = read_file("data/_cms_higgs_ZZto4l.txt",x);
    //RooDataHist HZZto4l = read_file("data/_cms_higgs_HZZto4l.txt",x); // leave away?

    RooHistPdf ttbar("ttbar", "ttbar", x, TTbarto4l);
    RooHistPdf drly("drly", "drly", x, DYto4l);
    RooHistPdf zz("zz", "zz", x, ZZto4l);

    const int sum_ttbar = TTbarto4l.sum(kTRUE);
    const int sum_drly = DYto4l.sum(kTRUE);
    const int sum_zz = ZZto4l.sum(kTRUE);
    const int sum_bkg = sum_ttbar + sum_drly + sum_zz;

    RooRealVar f_ttbar("f_ttbar", "fraction of ttbar", sum_ttbar/(double)sum_bkg);
    RooRealVar f_drly("f_drly", "fraction of drly", sum_drly/(double)sum_bkg);

    RooAddPdf bkg("bkg", "sum of backgrounds", RooArgList(ttbar,drly,zz), RooArgList(f_ttbar,f_drly));

    RooRealVar mean("mean", "Higgs mass", 125, 110, 140);
    RooRealVar sigma("sigma", "width of gaussian", 3);
    sigma.setConstant();

    RooGaussian sign("sign", "gaussian PDF for signal", x, mean, sigma);

    RooRealVar f_s("f_s", "fraction of signal", 0.1, 0, 1);

    RooAddPdf model("model", "model", RooArgList(sign, bkg), RooArgList(f_s));

    RooFitResult* fit_results = model.fitTo(data, PrintLevel(-1), Save());

    fit_results->Print();

    TH1* hist = model.createHistogram("model_hist", x);
    RooDataHist hist_data("hist_data", "Histogram Data", x, hist);
    RooHistPdf model_hist("model_hist", "Binned Model", x, hist_data);

    RooPlot *frame1 = x.frame(Title("Data + Fit"));
    data.plotOn(frame1);
    model_hist.plotOn(frame1, LineColor(kRed));  // Signal component
    model.plotOn(frame1, Components("bkg"), LineColor(kBlack)); // Background component
    model.paramOn(frame1);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    frame1->Draw();
    c1->SaveAs("output/CMS_fit.pdf");

    // Part 2 //

    RooWorkspace w("w", "workspace");
    ModelConfig mc("model_config", &w);
    mc.SetPdf(model);
    mc.SetObservables(RooArgSet(x));
    mc.SetParametersOfInterest(RooArgSet(*w.var("mean")));
    mc.SetNuisanceParameters(RooArgSet(*w.var("f_s")));
    w.import(mc);
    w.import(data);
    w.writeToFile("output/CMS_workspace.root");

    // create the class using data and model config
    ProfileLikelihoodCalculator plc(data, mc);
    // set the confidence level
    plc.SetConfidenceLevel(0.95);
    // compute the interval of the parameter mu
    LikelihoodInterval* interval = plc.GetInterval();
    // print the interval (get a pointer of mu from M.C.)
    auto poi = static_cast<RooRealVar*>(mc.GetParametersOfInterest()->first());
    double lowerLimit = interval->LowerLimit(*poi);
    double upperLimit = interval->UpperLimit(*poi);
    LikelihoodIntervalPlot plot(interval);

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    plot.Draw();
    c2->SaveAs("output/CMS_interval.pdf");

    // Part 3 //
 
    FeldmanCousins fc(data, mc);
    fc.SetConfidenceLevel(0.90);
    fc.UseAdaptiveSampling(true);
    fc.FluctuateNumDataEntries(false);
    fc.SetNBins(100); // number of points to test per parameter
    //fc.SetTestSize(.1);
 
    PointSetInterval* interval_fc = fc.GetInterval(); // "PDF not extentable"?
    RooRealVar *firstPOI = (RooRealVar *)mc.GetParametersOfInterest()->first();
    cout << "\n90% interval on " << firstPOI->GetName() << " is : [" << interval_fc->LowerLimit(*firstPOI) << ", " << interval_fc->UpperLimit(*firstPOI) << "] " << endl;

    // Part 4 //

    f_s.setVal(0);
    plc.SetNullParameters(f_s);
    HypoTestResult* hypotest = plc.GetHypoTest();
    double alpha = hypotest->NullPValue();
    double significance = hypotest->Significance();
    cout << "\n Alpha of Hypothesis test: " << alpha << ", Significance of Hypothesis test: " << significance << " sigma" << endl;
}
 