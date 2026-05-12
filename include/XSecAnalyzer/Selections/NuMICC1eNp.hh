#pragma once

#include <atomic>
#include <memory>
#include <string>

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class NuMICC1eNp : public SelectionBase {

public:

  NuMICC1eNp();
  ~NuMICC1eNp();

  virtual std::string categorize_event( AnalysisEvent& event ) override final;
  virtual bool is_selected( AnalysisEvent& event ) override final;

private:

  std::string beam_mode_;
  std::string file_type_;

  std::atomic<int> n_total_{0};
  std::atomic<int> n_pass_quality_{0};
  std::atomic<int> n_pass_sig_{0};
  std::atomic<int> n_pass_pre_{0};
  std::atomic<int> n_pass_bdt_loose_{0};
  std::atomic<int> n_pass_bdt_score_{0};
  std::atomic<int> n_pass_final_{0};

};
