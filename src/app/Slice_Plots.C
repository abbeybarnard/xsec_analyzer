// Standard library includes
#include <algorithm>
#include <cstdlib>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/MCC9SystematicsCalculator.hh"
#include "XSecAnalyzer/PlotUtils.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/SliceHistogram.hh"

using NFT = NtupleFileType;

//#define USE_FAKE_DATA ""

namespace {

  void set_mc_histogram_style( int event_category, TH1* mc_hist, int color ) {
    mc_hist->SetFillColor( color );
    mc_hist->SetLineColor( color );
    mc_hist->SetStats( false );
  }

  void set_ext_histogram_style( TH1* ext_hist ) {
    ext_hist->SetFillColor( 28 );
    ext_hist->SetLineColor( 28 );
    ext_hist->SetLineWidth( 2 );
    ext_hist->SetFillStyle( 3005 );
    ext_hist->SetStats( false );
  }

  void set_bnb_data_histogram_style( TH1* bnb_hist ) {

    bnb_hist->SetLineColor( kBlack );
    bnb_hist->SetLineWidth( 3 );
    bnb_hist->SetMarkerStyle( kFullCircle );
    bnb_hist->SetMarkerSize( 0.8 );
    bnb_hist->SetStats( false );

    bnb_hist->GetXaxis()->SetTitleOffset( 0.0 );
    bnb_hist->GetXaxis()->SetTitleSize( 0.0 );
    bnb_hist->GetYaxis()->SetTitleSize( 0.05 );
    bnb_hist->GetYaxis()->CenterTitle( true );
    bnb_hist->GetXaxis()->SetLabelSize( 0.0 );

    // This prevents the first y-axis label label (0) to be clipped by the
    // ratio plot
    bnb_hist->SetMinimum( 1e-3 );
  }

  void set_stat_err_histogram_style( TH1* stat_err_hist ) {
    stat_err_hist->SetFillColor( kBlack );
    stat_err_hist->SetLineColor( kBlack );
    stat_err_hist->SetLineWidth( 2 );
    stat_err_hist->SetFillStyle( 3004 );
  }

  // Derives the reco-bin-index groups needed by
  // SystematicsCalculator::enable_detvar_smoothing() from a SliceBinning
  // configuration: one group per Slice, listing the reco bins that make up
  // that 1D projection in bin order.
  std::vector< std::vector<size_t> > build_detvar_smoothing_groups(
    const SliceBinning& sb )
  {
    std::vector< std::vector<size_t> > groups;
    for ( const auto& slice : sb.slices_ ) {
      std::vector<size_t> group;
      // bin_map_ keys are 1-based histogram bins in order; each value is
      // the (already-sorted) set of contributing reco-bin indices.
      for ( const auto& bin_pair : slice.bin_map_ ) {
        for ( size_t reco_bin_idx : bin_pair.second ) {
          group.push_back( reco_bin_idx );
        }
      }
      if ( !group.empty() ) groups.push_back( group );
    }
    return groups;
  }

} // anonymous namespace

void tutorial_slice_plots(std::string FPM_Config, std::string SYST_Config, std::string SLICE_Config, std::string Univ_Output, std::string Plot_OutputDir) {

  // Counter to ensure plots aren't overwritten
  uint FileNameCounter = 0;
  std::string Plot_Prefix = "SlicePlots";
  std::string Plot_Suffix = ".pdf";
  std::string PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;

  std::cout << "\nRunning Slice_Plots with options:" << std::endl;
  std::cout << "\tFPM_Config: " << FPM_Config << std::endl;
  std::cout << "\tSYST_Config: " << SYST_Config << std::endl;
  std::cout << "\tSLICE_Config: " <<  SLICE_Config << std::endl;
  std::cout << "\tUniv_Output: " << Univ_Output << std::endl;
  std::cout << "\tPlot_OutputDir: " << Plot_OutputDir << std::endl;
  std::cout << "\t\tWith filename: " << PlotFileName << std::endl;
  std::cout << "\n" << std::endl;

/* Now always expects file properties file in input
#ifdef USE_FAKE_DATA
  // Initialize the FilePropertiesManager and tell it to treat the NuWro
  // MC ntuples as if they were data
  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( FPM_Config );
#endif
*/
  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( FPM_Config );

  // Check that we can read the universe output file
  TFile* temp_file = new TFile(Univ_Output.c_str(), "read");
  if (!temp_file || temp_file->IsZombie()) {
    std::cerr << "Could not read file: " << Univ_Output << std::endl;
    throw;
  }
  delete temp_file;

  auto* syst_ptr = new MCC9SystematicsCalculator(Univ_Output, SYST_Config);

  // include full systematics on signal for slice plots, rather than only on response
  syst_ptr->set_syst_mode(syst_ptr->SystMode::VaryBackgroundAndSignalDirectly);

  auto& syst = *syst_ptr;

  auto* sb_ptr = new SliceBinning( SLICE_Config );
  auto& sb = *sb_ptr;

  // Opt-in smoothing of the detector-variation systematic universes, to
  // tame low-MC-statistics noise (see SystematicsCalculator.hh). Off by
  // default; set XSEC_DETVAR_SMOOTHING to enable.
  if ( std::getenv( "XSEC_DETVAR_SMOOTHING" ) ) {
    syst.enable_detvar_smoothing( build_detvar_smoothing_groups( sb ) );
    std::cout << "DetVar smoothing ENABLED (XSEC_DETVAR_SMOOTHING set)\n";
  }

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  /*
  #ifdef USE_FAKE_DATA
    // Add the EXT to the "data" when working with fake data
    // No longer needed, this is already done in SystematicsCalculator
    reco_bnb_hist->Add( reco_ext_hist );
  #endif
  */

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  // --- Raw CV-vs-variation values for each individual DetVar, computed
  // once here (not per slice, since evaluate_observable() loops over every
  // true bin for every reco bin -- expensive to redo per slice), covering
  // both the unsmoothed and smoothed states regardless of whether
  // XSEC_DETVAR_SMOOTHING was set for this run. Unlike the covariance-based
  // breakdowns above, this comparison always shows both states so you
  // never need to rerun the script to see "before vs after".
  const std::vector< std::pair<std::string, NFT> > detvar_universe_types = {
    { "detVarLYdown", NFT::kDetVarMCLYdown },
    { "detVarLYrayl", NFT::kDetVarMCLYrayl },
    { "detVarLYatten", NFT::kDetVarMCLYatten },
    { "detVarRecomb2", NFT::kDetVarMCRecomb2 },
    { "detVarSCE", NFT::kDetVarMCSCE },
    { "detVarWMAngleXZ", NFT::kDetVarMCWMAngleXZ },
    { "detVarWMAngleYZ", NFT::kDetVarMCWMAngleYZ },
    { "detVarWMX", NFT::kDetVarMCWMX },
    { "detVarWMYZ", NFT::kDetVarMCWMYZ },
  };

  auto build_raw_array = [&]( const Universe& univ ) {
    size_t num_cm_bins = syst.get_covariance_matrix_size();
    std::vector<double> arr( num_cm_bins, 0. );
    for ( size_t rb = 0u; rb < num_cm_bins; ++rb ) {
      arr[rb] = syst.evaluate_observable( univ, rb );
    }
    return arr;
  };

  auto raw_plot_smoothing_groups = build_detvar_smoothing_groups( sb );

  std::map< std::string, std::vector<double> > raw_unsmoothed;
  std::map< std::string, std::vector<double> > raw_smoothed;

  bool have_detvar_cv = syst.detvar_universes_.count( NFT::kDetVarMCCV ) != 0;
  if ( !have_detvar_cv ) {
    std::cout << "\nWarning: no detVarCV universe found in this file --"
      " skipping the raw CV-vs-variation DetVar comparison plot.\n";
  }
  else {
    const Universe& cv_univ = *syst.detvar_universes_.at( NFT::kDetVarMCCV );
    raw_unsmoothed[ "detVarCV" ] = build_raw_array( cv_univ );
    raw_smoothed[ "detVarCV" ] = raw_unsmoothed.at( "detVarCV" );
    SystematicsCalculator::smooth_values_using_groups(
      raw_smoothed.at( "detVarCV" ), raw_plot_smoothing_groups );

    for ( const auto& dv : detvar_universe_types ) {
      if ( !syst.detvar_universes_.count( dv.second ) ) {
        std::cout << "  Skipping " << dv.first
          << " in the raw comparison plot: universe not found in this file.\n";
        continue;
      }
      const Universe& var_univ = *syst.detvar_universes_.at( dv.second );
      raw_unsmoothed[ dv.first ] = build_raw_array( var_univ );
      raw_smoothed[ dv.first ] = raw_unsmoothed.at( dv.first );
      SystematicsCalculator::smooth_values_using_groups(
        raw_smoothed.at( dv.first ), raw_plot_smoothing_groups );
    }
  }

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("PredTotal") );

    auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
      << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
      << " p-value = " << chi2_result.p_value_ << '\n';

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    const auto& sel_for_cat = syst.get_selection_for_categories();
    const auto& cat_map = sel_for_cat.category_map();

    // Go in reverse so that, if the signal is defined first in the map, it
    // ends up on top. Note that this index is one-based to match the ROOT
    // histograms
    int cat_bin_index = cat_map.size();
    for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
    {
      int cat = iter->first;
      int color = iter->second.color_;
      TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
        cat_bin_index, cat_bin_index );
      temp_mc_hist->SetDirectory( nullptr );

      SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
        *temp_mc_hist, slice  );

      set_mc_histogram_style( cat, temp_slice_mc->hist_.get(), color );

      slice_pred_stack->Add( temp_slice_mc->hist_.get() );

      std::string cat_col_prefix = "MC" + std::to_string( cat );

      --cat_bin_index;
    }

    TCanvas* c1 = new TCanvas;
    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 3 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 1.2 );
    slice_bnb->hist_->SetStats( false );
    double ymax = std::max( slice_bnb->hist_->GetMaximum(),
      slice_mc_plus_ext->hist_->GetMaximum() ) * 1.2;
    slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );

    slice_bnb->hist_->Draw( "e" );

    slice_pred_stack->Draw( "hist same" );

    slice_mc_plus_ext->hist_->SetLineWidth( 3 );
    slice_mc_plus_ext->hist_->SetLineColor(kRed);
    slice_mc_plus_ext->hist_->Draw( "same hist e" );

    slice_bnb->hist_->Draw( "same e" );

    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;
    c1->SaveAs(PlotFileName.c_str());
    FileNameCounter += 1;

    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast< TH1* >(
      slice.hist_->Clone("slice_hist") );

    slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Show fractional uncertainties computed using these covariance matrices
    // in the ROOT plot. All configured fractional uncertainties will be
    // included in the output pgfplots file regardless of whether they appear
    // in this vector.
    const std::vector< std::string > cov_mat_keys = { "PredTotal",
      "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets",
      "MCstats", "EXTstats", "BNBstats",
      "detVarLYdown", "detVarLYrayl", "detVarLYatten", "detVarRecomb2",
      "detVarSCE", "detVarWMAngleXZ", "detVarWMAngleYZ", "detVarWMX",
      "detVarWMYZ"
    };

    // Subset of the above used for the dedicated individual-DetVar
    // breakdown plot below (detVar_total is included for reference).
    const std::vector< std::string > detvar_breakdown_keys = {
      "detVar_total", "detVarLYdown", "detVarLYrayl", "detVarLYatten",
      "detVarRecomb2", "detVarSCE", "detVarWMAngleXZ", "detVarWMAngleYZ",
      "detVarWMX", "detVarWMYZ"
    };

    // Loop over the various systematic uncertainties
    int color = 0;
    for ( const auto& pair : matrix_map ) {

      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;

      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
        double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }

      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend ) continue;

      frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 10;

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );
    }

    TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );

    auto* total_frac_err_hist = frac_uncertainty_hists.at( "PredTotal" );
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );

    lg2->AddEntry( total_frac_err_hist, "PredTotal", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "PredTotal" ) continue;

      // The individual DetVar systematics are reserved for their own
      // dedicated breakdown plot below (detVar_total, their sum, still
      // belongs on this combined plot).
      bool is_individual_detvar = name != "detVar_total"
        && std::find( detvar_breakdown_keys.cbegin(),
             detvar_breakdown_keys.cend(), name )
           != detvar_breakdown_keys.cend();
      if ( is_individual_detvar ) continue;

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );

      std::cout << name << " frac err in bin #1 = "
        << hist->GetBinContent( 1 )*100. << "%\n";
    }

    lg2->Draw( "same" );

    std::cout << "Total frac error in bin #1 = "
      << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";

    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;
    c2->SaveAs(PlotFileName.c_str());
    FileNameCounter += 1;

    // Dedicated breakdown of the individual DetVar systematics (rather than
    // just their sum, detVar_total). Lets us check whether the DetVars
    // expected to matter more for a given slice actually stand out, or
    // whether they've all just been shrunk to a smaller, still-noise-
    // dominated level by the smoothing.
    TCanvas* c3 = new TCanvas;
    TLegend* lg3 = new TLegend( 0.7, 0.7, 0.9, 0.9 );

    auto* detvar_total_hist = frac_uncertainty_hists.at( "detVar_total" );
    detvar_total_hist->SetStats( false );
    detvar_total_hist->GetYaxis()->SetRangeUser( 0.,
      detvar_total_hist->GetMaximum() * 1.05 );
    detvar_total_hist->SetLineColor( kBlack );
    detvar_total_hist->SetLineWidth( 3 );
    detvar_total_hist->Draw( "hist" );

    lg3->AddEntry( detvar_total_hist, "detVar_total", "l" );

    std::cout << "detVar_total frac err in bin #1 = "
      << detvar_total_hist->GetBinContent( 1 )*100. << "%\n";

    for ( const auto& name : detvar_breakdown_keys ) {
      // Already plotted above
      if ( name == "detVar_total" ) continue;

      TH1* hist = frac_uncertainty_hists.at( name );
      lg3->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );

      std::cout << "  " << name << " frac err in bin #1 = "
        << hist->GetBinContent( 1 )*100. << "%\n";
    }

    lg3->Draw( "same" );

    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;
    c3->SaveAs(PlotFileName.c_str());
    FileNameCounter += 1;

    // Raw CV-vs-variation values (not a derived uncertainty or covariance
    // matrix) for each individual DetVar, overlaid, baseline vs. smoothed
    // side by side -- lets you check whether smoothing visibly reshapes
    // the actual distributions, not just their derived uncertainties.
    if ( have_detvar_cv ) {

      TCanvas* c4 = new TCanvas;
      c4->Divide( 2, 1 );
      TLegend* lg4 = new TLegend( 0.7, 0.7, 0.9, 0.9 );

      // Projects a raw (full-reco-space) array into this slice's binning
      // via the existing SliceHistogram machinery (nullptr covariance --
      // we only want values here, not error bars).
      auto project_to_slice = [&]( const std::vector<double>& arr,
        const std::string& hist_name ) -> TH1* {
        TH1D raw_hist( hist_name.c_str(), "",
          (int)arr.size(), 0., (double)arr.size() );
        for ( size_t rb = 0u; rb < arr.size(); ++rb ) {
          raw_hist.SetBinContent( rb + 1, arr[rb] );
        }
        SliceHistogram* sh = SliceHistogram::make_slice_histogram(
          raw_hist, slice, nullptr );
        return sh->hist_.get();
      };

      // Fixed color per DetVar (same scheme as the breakdown above) so the
      // same DetVar has the same color in both the baseline and smoothed
      // pads, and matches the rest of this file's style.
      std::map< std::string, int > detvar_colors;
      {
        int color = 0;
        for ( const auto& dv : detvar_universe_types ) {
          if ( color <= 9 ) ++color;
          if ( color == 5 ) ++color;
          if ( color >= 10 ) color += 10;
          detvar_colors[ dv.first ] = color;
        }
      }

      // Build every histogram for both pads first (without drawing), so we
      // can find the shared maximum and put both pads on identical y-axes
      // before anything is actually drawn.
      struct PadHist { TH1* hist; int pad; bool is_cv; bool is_ref; std::string dv_name; };
      std::vector< PadHist > pad_hists;
      double shared_max = 0.;

      for ( int pad = 1; pad <= 2; ++pad ) {
        const auto& raw_map = ( pad == 1 ) ? raw_unsmoothed : raw_smoothed;
        const std::string pad_label = ( pad == 1 ) ? "baseline" : "smoothed";

        TH1* cv_slice_hist = project_to_slice( raw_map.at( "detVarCV" ),
          "detvarcv_raw_" + pad_label + Form( "_%zu", sl_idx ) );
        pad_hists.push_back( { cv_slice_hist, pad, true, false, "" } );
        shared_max = std::max( shared_max, cv_slice_hist->GetMaximum() );

        // The smoothed pad's CV curve above is the internal working copy
        // used only to compute the covariance numerator -- it never feeds
        // the actual reported prediction (see reco_mc_plus_ext_hist, built
        // straight from syst.cv_universe().hist_reco_, never smoothed).
        // Overlay the always-unsmoothed CV as a dashed reference so it's
        // visually obvious the true prediction hasn't moved.
        if ( pad == 2 ) {
          TH1* cv_ref_hist = project_to_slice( raw_unsmoothed.at( "detVarCV" ),
            "detvarcv_ref_" + pad_label + Form( "_%zu", sl_idx ) );
          pad_hists.push_back( { cv_ref_hist, pad, false, true, "" } );
          shared_max = std::max( shared_max, cv_ref_hist->GetMaximum() );
        }

        for ( const auto& dv : detvar_universe_types ) {
          if ( !raw_map.count( dv.first ) ) continue;

          TH1* var_slice_hist = project_to_slice( raw_map.at( dv.first ),
            dv.first + "_raw_" + pad_label + Form( "_%zu", sl_idx ) );
          pad_hists.push_back( { var_slice_hist, pad, false, false, dv.first } );
          shared_max = std::max( shared_max, var_slice_hist->GetMaximum() );
        }
      }

      shared_max *= 1.05;

      for ( int pad = 1; pad <= 2; ++pad ) {
        const std::string pad_label = ( pad == 1 ) ? "baseline" : "smoothed";
        c4->cd( pad );

        bool first_in_pad = true;
        for ( auto& ph : pad_hists ) {
          if ( ph.pad != pad ) continue;

          if ( ph.is_cv ) {
            ph.hist->SetStats( false );
            ph.hist->SetTitle(
              ( "Raw DetVar CV/variation (" + pad_label + ")" ).c_str() );
            ph.hist->GetYaxis()->SetRangeUser( 0., shared_max );
            ph.hist->SetLineColor( kBlack );
            ph.hist->SetLineWidth( 3 );
            ph.hist->Draw( first_in_pad ? "hist" : "same hist" );
            if ( pad == 1 ) lg4->AddEntry( ph.hist, "detVarCV", "l" );
          }
          else if ( ph.is_ref ) {
            ph.hist->SetLineColor( kBlack );
            ph.hist->SetLineStyle( kDashed );
            ph.hist->SetLineWidth( 2 );
            ph.hist->Draw( first_in_pad ? "hist" : "same hist" );
            lg4->AddEntry( ph.hist, "detVarCV (unsmoothed reference)", "l" );
          }
          else {
            ph.hist->SetLineColor( detvar_colors.at( ph.dv_name ) );
            ph.hist->SetLineWidth( 2 );
            ph.hist->Draw( first_in_pad ? "hist" : "same hist" );
            if ( pad == 1 ) lg4->AddEntry( ph.hist, ph.dv_name.c_str(), "l" );
          }
          first_in_pad = false;
        }
      }

      c4->cd( 1 );
      lg4->Draw( "same" );

      PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;
      c4->SaveAs(PlotFileName.c_str());
      FileNameCounter += 1;
    }

  } // slices

}

int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    std::cout << "Usage: Slice_Plots FPM_CONFIG"
	      << " SYST_Config SLICE_Config Univ_Output Plot_OutputDir\n";
    return 1;
  }

  std::string list_file_name( argv[1] );
  std::string univmake_config_file_name( argv[2] );
  std::string output_file_name( argv[3] );

  std::string FPM_Config( argv[1] );
  std::string SYST_Config( argv[2] );
  std::string SLICE_Config( argv[3] );
  std::string Univ_Output( argv[4] );
  std::string Plot_OutputDir( argv[5] );

  tutorial_slice_plots(FPM_Config, SYST_Config, SLICE_Config, Univ_Output, Plot_OutputDir);
  return 0;
}
