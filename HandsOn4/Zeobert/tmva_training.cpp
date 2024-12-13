#include <cstdlib>
#include <iostream>
#include <map>
#include <list>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

using namespace TMVA;

void tmva_training(){

    Tools::Instance();
    auto outputFile = TFile::Open("output/TMVAOutputCV.root", "RECREATE");
    Factory factory("TMVAClassification", outputFile,
                        "!V:ROC:!Correlations:!Silent:Color:"
                        "!DrawProgressBar:AnalysisType=Classification");
    DataLoader dataloader("dataset");

    std::list<std::string> vars = {"DER_mass_MMC","DER_mass_transverse_met_lep","DER_mass_vis","DER_pt_h","DER_deltaeta_jet_jet","DER_mass_jet_jet","DER_prodeta_jet_jet","DER_deltar_tau_lep","DER_pt_tot","DER_sum_pt","DER_pt_ratio_lep_tau","DER_met_phi_centrality","DER_lep_eta_centrality","PRI_tau_pt","PRI_tau_eta","PRI_tau_phi","PRI_lep_pt","PRI_lep_eta","PRI_lep_phi","PRI_met","PRI_met_phi","PRI_met_sumet","PRI_jet_num","PRI_jet_leading_pt","PRI_jet_leading_eta","PRI_jet_leading_phi","PRI_jet_subleading_pt","PRI_jet_subleading_eta","PRI_jet_subleading_phi","PRI_jet_all_pt"};
    for(auto var : vars){
        dataloader.AddVariable(var, 'F');}

    TFile signal_file("data/atlas-higgs-challenge-2014-v2-sig.root");
    TTree *signal_tree = (TTree*)signal_file.Get("tree");

    TFile bkg_file("data/atlas-higgs-challenge-2014-v2-bkg.root");
    TTree *bkg_tree = (TTree*)bkg_file.Get("tree");

    dataloader.AddSignalTree(signal_tree, 1.0);
    dataloader.AddBackgroundTree(bkg_tree, 1.0);

    dataloader.SetSignalWeightExpression("Weight");
    dataloader.SetBackgroundWeightExpression("Weight");

    TCut mycuts, mycutb;
    dataloader.PrepareTrainingAndTestTree(mycuts, mycutb,
                                            "nTrain_Signal=10000:nTrain_Background=20000:"
                                            "SplitMode=Random:NormMode=NumEvents:!V");

    factory.BookMethod(&dataloader, TMVA::Types::kCuts, "Cuts", "!H:!V:FitMethod=MC");
    factory.BookMethod(&dataloader, TMVA::Types::kFisher, "Fisher", "!H:!V");
	factory.BookMethod(&dataloader, TMVA::Types::kMLP, "MLP",
	 		 "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=100:"
			 "HiddenLayers=N+5,3:TestRate=5:!UseRegulator");
    factory.BookMethod(&dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=200:BoostType=AdaBoost");

    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    outputFile->Close();
    TMVAGui("output/TMVAOutputCV.root");
}