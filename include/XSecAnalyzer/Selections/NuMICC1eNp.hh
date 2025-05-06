#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class NuMICC1eNp : public SelectionBase {

public:

  NuMICC1eNp();

  virtual std::string categorize_event( AnalysisEvent& event ) override final;
  virtual bool is_selected( AnalysisEvent& event ) override final;

};
