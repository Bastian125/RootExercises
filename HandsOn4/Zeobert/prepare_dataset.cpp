//  allows to analyse data stored in TTrees (CSVs) with a high level interface.

#include <iostream>
#include "/opt/homebrew/Cellar/root/6.32.08/include/root/ROOT/RDataFrame.hxx"
#include "/opt/homebrew/Cellar/root/6.32.08/include/root/ROOT/RCsvDS.hxx"

using namespace std;
using namespace ROOT::RDF;

void prepare_dataset()
{
  auto df = FromCSV("data/atlas-higgs-challenge-2014-v2.csv");

  // Get the number of entries

  auto nevents = *df.Count(); 
  cout << "n. dati " << nevents << "\n"; 
  // auto h1 = df.Histo1D("totale_casi");

  // Returns the names of the available columns.

  auto colNames = df.GetColumnNames();
  for (auto &&colName : colNames)  {
    std::cout << colName << ", " ;  
  }
  std::cout << '\n';

  // auto colNames = df.GetColumnNames();
  // for (auto &&colName : colNames)  {
  //   std::cout << "  loader.AddVariable(\"" << colName << "\", 'F');  \t// "
  //     << colName << " (" << df.GetColumnType(colName) << ")\n";
  // }

  // Save selected columns to disk, in a new TTree in new file.
  df.Filter("Label == \"s\" ")
    .Snapshot("tree", "data/atlas-higgs-challenge-2014-v2-sig.root");  

  df.Filter("Label == \"b\" ")
    .Snapshot("tree", "data/atlas-higgs-challenge-2014-v2-bkg.root");

}
