#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

// Constants removed from framework — define locally
constexpr int kSignalTrueBin = 0;
constexpr int kOrdinaryRecoBin = 0;

class TutorialBinScheme : public BinSchemeBase {

  public:

    TutorialBinScheme();
    virtual void DefineBlocks() override;
};
