#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

// Constants removed from framework — define locally

class TutorialBinScheme : public BinSchemeBase {

  public:

    TutorialBinScheme();
    virtual void DefineBlocks() override;
};
