#include "TSystem.h"

#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/QuickPlotter.hh"

void qp_demo() {
  // Uncomment the sample name that you would like to look at
  std::string sample_name = "mcc9";
  //std::string sample_name = "mcc9.10 Pandora";
  //std::string sample_name = "mcc9.10 Super Unified";
  
  QuickPlotter qp( 1, "NuMICC1eNp" );
  
  auto& fpm = FilePropertiesManager::Instance();
  std::string xsec_analyzer_dir = gSystem->Getenv( "XSEC_ANALYZER_DIR" );
  fpm.load_file_properties( xsec_analyzer_dir + "/configs/file_properties.txt" );
  qp.plot( "shr_energy_cali",
    "NuMICC1eNp_Selected",
    0., 1., 40, { "nuselection/NeutrinoSelectionFilter" },
    "Electron Energy", "Events", "Run 1" );
}