// Compares the systematic covariance and correlation matrices for
// detVar_total (the sum) and each of the 9 individual DetVar systematics,
// with and without the opt-in DetVar smoothing
// (SystematicsCalculator::enable_detvar_smoothing), to help judge whether
// the smoothing is cleaning up statistical noise in the detector-variation
// MC samples, or also smoothing out real physical structure -- and whether
// the DetVars expected to matter more (e.g. SCE, Recomb2) actually stand
// out from the others, rather than everything just being uniformly-shrunk
// noise.

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"

#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/MCC9SystematicsCalculator.hh"
#include "XSecAnalyzer/SliceBinning.hh"

using NFT = NtupleFileType;

namespace {

  // Mirrors Slice_Plots.C's build_detvar_smoothing_groups(): derives the
  // reco-bin-index groups enable_detvar_smoothing() expects from a
  // SliceBinning configuration, one group per 1D slice, so smoothing never
  // mixes bins from unrelated projections.
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

  // Builds a correlation TH2D from a covariance TH2D:
  // Corr(i,j) = Cov(i,j) / sqrt(Cov(i,i) * Cov(j,j))
  // No existing correlation-matrix utility exists in the codebase (checked)
  // so this is computed directly here.
  TH2D* covariance_to_correlation( const TH2D& cov, const std::string& name )
  {
    int n = cov.GetNbinsX();
    TH2D* corr = dynamic_cast< TH2D* >( cov.Clone( name.c_str() ) );
    corr->Reset();
    corr->SetDirectory( nullptr );

    std::vector< double > diag( n );
    for ( int a = 0; a < n; ++a ) diag[a] = cov.GetBinContent( a + 1, a + 1 );

    for ( int a = 0; a < n; ++a ) {
      for ( int b = 0; b < n; ++b ) {
        double denom = std::sqrt( diag[a] * diag[b] );
        double value = ( denom > 0. )
          ? cov.GetBinContent( a + 1, b + 1 ) / denom : 0.;
        corr->SetBinContent( a + 1, b + 1, value );
      }
    }
    return corr;
  }

  // Builds a short, physically-meaningful label for each of the num_cm_bins
  // global reco bins (e.g. "T_p[0.05,0.10) @ E_e[0.42,0.62)"), so the
  // matrices are self-describing once pulled into a ROOT file (e.g. via
  // uproot) without needing to cross-reference the slice config separately.
  // Bins not covered by any slice (see the background/other-vars catch-all
  // discussed with the user) are labeled accordingly rather than left blank.
  std::vector< std::string > build_bin_labels( const SliceBinning& sb,
    size_t num_cm_bins )
  {
    std::vector< std::string > labels( num_cm_bins, "bkg/other" );

    for ( const auto& slice : sb.slices_ ) {
      std::string active_name = slice.hist_->GetXaxis()->GetTitle();

      std::string other_desc;
      for ( const auto& ovs : slice.other_vars_ ) {
        if ( !other_desc.empty() ) other_desc += ", ";
        const auto& svar = sb.slice_vars_.at( ovs.var_index_ );
        other_desc += svar.name_ + "[" + std::to_string( ovs.low_bin_edge_ )
          + "," + std::to_string( ovs.high_bin_edge_ ) + ")";
      }

      for ( const auto& bin_pair : slice.bin_map_ ) {
        int local_bin = bin_pair.first; // 1-based TH1 bin in slice.hist_
        double lo = slice.hist_->GetBinLowEdge( local_bin );
        double hi = slice.hist_->GetBinLowEdge( local_bin + 1 );

        std::ostringstream oss;
        oss << active_name << "[" << lo << "," << hi << ")";
        if ( !other_desc.empty() ) oss << " @ " << other_desc;

        for ( size_t reco_bin_idx : bin_pair.second ) {
          if ( reco_bin_idx < num_cm_bins ) labels.at( reco_bin_idx ) = oss.str();
        }
      }
    }

    return labels;
  }

  // Clones a matrix and stamps the given per-bin labels onto both axes, so
  // the copy written out for uproot carries the binning scheme with it.
  // Kept as a separate clone (rather than labeling the PDF-drawn original)
  // so today's already-validated COLZ comparison plots are unaffected --
  // ROOT switches an axis to categorical (text) mode as soon as any bin
  // label is set on it.
  TH2D* make_labeled_clone( const TH2D& src, const std::string& name,
    const std::vector< std::string >& labels )
  {
    TH2D* clone = dynamic_cast< TH2D* >( src.Clone( name.c_str() ) );
    clone->SetDirectory( nullptr );
    for ( size_t i = 0u; i < labels.size(); ++i ) {
      clone->GetXaxis()->SetBinLabel( i + 1, labels[i].c_str() );
      clone->GetYaxis()->SetBinLabel( i + 1, labels[i].c_str() );
    }
    return clone;
  }

  // Blue-white-red diverging palette, centered at zero -- appropriate here
  // since both covariance and correlation matrices can have negative
  // off-diagonal entries. Avoids ROOT's sequential built-in palettes, which
  // would misleadingly imply the data has no natural zero point.
  void set_diverging_palette() {
    const int n_stops = 3;
    double stops[ n_stops ] = { 0.00, 0.50, 1.00 };
    double red[ n_stops ]   = { 0.230, 1.000, 0.706 };
    double green[ n_stops ] = { 0.299, 1.000, 0.016 };
    double blue[ n_stops ]  = { 0.754, 1.000, 0.150 };
    TColor::CreateGradientColorTable( n_stops, stops, red, green, blue, 255 );
    gStyle->SetNumberContours( 255 );
  }

} // anonymous namespace

int main( int argc, char* argv[] ) {

  if ( argc != 6 ) {
    std::cout << "Usage: DetVarCovarianceComparison FPM_Config"
      << " SYST_Config SLICE_Config Univ_Output Plot_OutputDir\n";
    return 1;
  }

  std::string FPM_Config( argv[1] );
  std::string SYST_Config( argv[2] );
  std::string SLICE_Config( argv[3] );
  std::string Univ_Output( argv[4] );
  std::string Plot_OutputDir( argv[5] );

  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( FPM_Config );

  SliceBinning sb( SLICE_Config );
  auto smoothing_groups = build_detvar_smoothing_groups( sb );

  std::cout << "\nConstructing baseline (unsmoothed) calculator...\n";
  MCC9SystematicsCalculator baseline( Univ_Output, SYST_Config );
  // Match Slice_Plots.C exactly: include full systematics on signal for
  // slice plots, rather than only on response. Without this, evaluate_
  // observable() uses the default SystMode::ForXSec instead, which is a
  // genuinely different computation -- not merely a different scale.
  baseline.set_syst_mode( baseline.SystMode::VaryBackgroundAndSignalDirectly );

  std::cout << "Constructing smoothed calculator...\n";
  MCC9SystematicsCalculator smoothed( Univ_Output, SYST_Config );
  smoothed.set_syst_mode( smoothed.SystMode::VaryBackgroundAndSignalDirectly );
  smoothed.enable_detvar_smoothing( smoothing_groups );

  std::cout << "Computing baseline covariances...\n";
  auto baseline_cov_map = baseline.get_covariances();
  std::cout << "Computing smoothed covariances...\n";
  auto smoothed_cov_map = smoothed.get_covariances();

  // Sanity check against earlier findings: bins 25 and 67 (the first bins
  // of slices 5 and 12) were found to have their fractional DetVar
  // uncertainty *increase* after smoothing, traced to bin 25 having far
  // fewer raw events than its neighbor bin 26.
  {
    TH2D* total_before = baseline_cov_map->at( "detVar_total" ).cov_matrix_.get();
    TH2D* total_after  = smoothed_cov_map->at( "detVar_total" ).cov_matrix_.get();
    for ( int b : { 25, 67 } ) {
      double before = total_before->GetBinContent( b + 1, b + 1 );
      double after  = total_after->GetBinContent( b + 1, b + 1 );
      std::cout << "bin " << b << ": baseline cov = " << before
        << ", smoothed cov = " << after << ", change = "
        << ( ( after > before ) ? "UP" : "down" ) << " ("
        << 100. * ( after - before ) / before << "%)\n";
    }
  }

  gStyle->SetOptStat( 0 );
  set_diverging_palette();

  // One page per systematic: the sum (detVar_total) first, for reference,
  // then each of the 9 individual DetVars, then every other systematic type
  // configured in systcalc_numi.conf (matches Slice_Plots.C's cov_mat_keys
  // exactly). Each page's color scale is computed from *that page's own*
  // matrices only -- sharing one scale across all of them would wash out the
  // individual DetVars, which are much smaller than their quadrature sum.
  // Note: smoothing only ever affects DetVar-type covariances (it operates
  // on detvar_universes_ specifically) -- baseline and smoothed will be
  // identical for flux/reint/xsec_total/etc., which is expected, not a bug.
  const std::vector< std::string > all_keys = {
    "detVar_total", "detVarLYdown", "detVarLYrayl", "detVarLYatten",
    "detVarRecomb2", "detVarSCE", "detVarWMAngleXZ", "detVarWMAngleYZ",
    "detVarWMX", "detVarWMYZ",
    "PredTotal", "flux", "reint", "xsec_total", "POT", "numTargets",
    "MCstats", "EXTstats", "BNBstats", "dirtNorm"
  };
  std::vector< std::string > systematics_to_plot;
  for ( const auto& name : all_keys ) {
    if ( baseline_cov_map->count( name ) && smoothed_cov_map->count( name ) ) {
      systematics_to_plot.push_back( name );
    }
    else {
      std::cout << "Skipping \"" << name
        << "\" -- not present in get_covariances() output "
        << "(not configured in " << SYST_Config << ")\n";
    }
  }

  // Per-bin physical labels ("T_p[0.05,0.10) @ E_e[0.42,0.62)", etc.),
  // shared by every systematic since they all use the same reco binning.
  size_t num_cm_bins = baseline.get_covariance_matrix_size();
  auto bin_labels = build_bin_labels( sb, num_cm_bins );

  // The true reported prediction (EXT + CV MC), built exactly as Slice_
  // Plots.C does -- never touched by smoothing (that's local to make_cov_mat
  // and never modifies data_hists_/cv_universe()). Needed as the denominator
  // to turn a covariance diagonal into a *fractional* uncertainty; without
  // it, only absolute covariance/correlation values are available.
  TH1D* reco_ext_hist = baseline.data_hists_.at( NFT::kExtBNB ).get();
  std::unique_ptr< TH1D > reco_mc_plus_ext_hist( dynamic_cast< TH1D* >(
    reco_ext_hist->Clone( "reco_mc_plus_ext_prediction" ) ) );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );
  reco_mc_plus_ext_hist->Add( baseline.cv_universe().hist_reco_.get() );

  // Raw (not derived-uncertainty) per-bin CV/variation values for each
  // individual DetVar -- mirrors Slice_Plots.C's c4 canvas exactly. Smoothing
  // never touches evaluate_observable() itself (it's local to make_cov_mat),
  // so the same raw array is smoothed here manually via the same public
  // smooth_values_using_groups() method, using the same slice-derived groups.
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
    std::vector<double> arr( num_cm_bins, 0. );
    for ( size_t rb = 0u; rb < num_cm_bins; ++rb ) {
      arr[rb] = baseline.evaluate_observable( univ, rb );
    }
    return arr;
  };

  bool have_detvar_cv = baseline.detvar_universes_.count( NFT::kDetVarMCCV ) != 0;
  std::map< std::string, std::vector<double> > raw_unsmoothed;
  std::map< std::string, std::vector<double> > raw_smoothed;

  if ( have_detvar_cv ) {
    raw_unsmoothed[ "detVarCV" ] = build_raw_array(
      *baseline.detvar_universes_.at( NFT::kDetVarMCCV ) );
    for ( const auto& dv : detvar_universe_types ) {
      if ( !baseline.detvar_universes_.count( dv.second ) ) continue;
      raw_unsmoothed[ dv.first ] = build_raw_array(
        *baseline.detvar_universes_.at( dv.second ) );
    }
    raw_smoothed = raw_unsmoothed;
    for ( auto& pair : raw_smoothed ) {
      SystematicsCalculator::smooth_values_using_groups(
        pair.second, smoothing_groups );
    }
  }
  else {
    std::cout << "No detVarCV universe present -- skipping raw"
      " CV/variation value export.\n";
  }

  std::string plot_file_name = Plot_OutputDir + "/DetVarCovarianceComparison.pdf";
  std::string root_file_name = Plot_OutputDir + "/DetVarCovarianceMatrices.root";
  TFile out_root_file( root_file_name.c_str(), "recreate" );

  out_root_file.cd();
  std::unique_ptr< TH1D > pred_labeled( dynamic_cast< TH1D* >(
    reco_mc_plus_ext_hist->Clone( "reco_mc_plus_ext_prediction" ) ) );
  pred_labeled->SetDirectory( nullptr );
  for ( size_t i = 0u; i < bin_labels.size(); ++i ) {
    pred_labeled->GetXaxis()->SetBinLabel( i + 1, bin_labels[i].c_str() );
  }
  pred_labeled->Write();

  // Raw CV/variation arrays, written out as labeled TH1D's the same way,
  // one pair (unsmoothed/smoothed) per array.
  auto write_raw_array_as_hist = [&]( const std::string& hist_name,
    const std::vector<double>& arr ) {
    TH1D h( hist_name.c_str(), "", (int)arr.size(), 0., (double)arr.size() );
    h.SetDirectory( nullptr );
    for ( size_t i = 0u; i < arr.size(); ++i ) {
      h.SetBinContent( i + 1, arr[i] );
      h.GetXaxis()->SetBinLabel( i + 1, bin_labels[i].c_str() );
    }
    h.Write();
  };
  for ( const auto& pair : raw_unsmoothed ) {
    write_raw_array_as_hist( pair.first + "_raw_unsmoothed", pair.second );
  }
  for ( const auto& pair : raw_smoothed ) {
    write_raw_array_as_hist( pair.first + "_raw_smoothed", pair.second );
  }

  TCanvas c( "c", "", 1400, 1200 );

  for ( size_t s = 0u; s < systematics_to_plot.size(); ++s ) {
    const std::string& name = systematics_to_plot[s];

    TH2D* cov_before = baseline_cov_map->at( name ).cov_matrix_.get();
    TH2D* cov_after  = smoothed_cov_map->at( name ).cov_matrix_.get();

    std::unique_ptr< TH2D > corr_before(
      covariance_to_correlation( *cov_before, name + "_corr_before" ) );
    std::unique_ptr< TH2D > corr_after(
      covariance_to_correlation( *cov_after, name + "_corr_after" ) );

    // Labeled clones, written out for later use in Python (e.g. via
    // uproot) -- the PDF below keeps drawing the plain numeric-axis
    // originals so today's already-validated comparison view is unchanged.
    out_root_file.cd();
    std::unique_ptr< TH2D > cov_before_labeled( make_labeled_clone(
      *cov_before, name + "_cov_baseline", bin_labels ) );
    std::unique_ptr< TH2D > cov_after_labeled( make_labeled_clone(
      *cov_after, name + "_cov_smoothed", bin_labels ) );
    std::unique_ptr< TH2D > corr_before_labeled( make_labeled_clone(
      *corr_before, name + "_corr_baseline", bin_labels ) );
    std::unique_ptr< TH2D > corr_after_labeled( make_labeled_clone(
      *corr_after, name + "_corr_smoothed", bin_labels ) );
    cov_before_labeled->Write();
    cov_after_labeled->Write();
    corr_before_labeled->Write();
    corr_after_labeled->Write();

    double cov_scale = std::max(
      std::max( std::abs( cov_before->GetMaximum() ),
        std::abs( cov_before->GetMinimum() ) ),
      std::max( std::abs( cov_after->GetMaximum() ),
        std::abs( cov_after->GetMinimum() ) )
    );
    if ( cov_scale <= 0. ) cov_scale = 1.; // guard an all-zero matrix

    c.Clear();
    c.Divide( 2, 2 );

    c.cd( 1 );
    cov_before->SetTitle( ( name + " covariance -- baseline" ).c_str() );
    cov_before->GetZaxis()->SetRangeUser( -cov_scale, cov_scale );
    cov_before->Draw( "COLZ" );

    c.cd( 2 );
    cov_after->SetTitle( ( name + " covariance -- smoothed" ).c_str() );
    cov_after->GetZaxis()->SetRangeUser( -cov_scale, cov_scale );
    cov_after->Draw( "COLZ" );

    c.cd( 3 );
    corr_before->SetTitle( ( name + " correlation -- baseline" ).c_str() );
    corr_before->GetZaxis()->SetRangeUser( -1., 1. );
    corr_before->Draw( "COLZ" );

    c.cd( 4 );
    corr_after->SetTitle( ( name + " correlation -- smoothed" ).c_str() );
    corr_after->GetZaxis()->SetRangeUser( -1., 1. );
    corr_after->Draw( "COLZ" );

    // ROOT's multi-page PDF idiom: "(" opens the document on the first
    // page, ")" closes it on the last page, plain filename for every page
    // in between.
    std::string page_spec = plot_file_name;
    if ( s == 0u ) page_spec += "(";
    else if ( s == systematics_to_plot.size() - 1u ) page_spec += ")";
    c.Print( page_spec.c_str() );

    std::cout << "  Drew page " << ( s + 1 ) << "/"
      << systematics_to_plot.size() << ": " << name << '\n';
  }

  std::cout << "\nWrote comparison plot (" << systematics_to_plot.size()
    << " pages) to: " << plot_file_name << '\n';

  out_root_file.Write();
  out_root_file.Close();
  std::cout << "Wrote labeled covariance/correlation matrices ("
    << ( systematics_to_plot.size() * 4 ) << " histograms), 1 prediction"
    " histogram, and " << ( raw_unsmoothed.size() + raw_smoothed.size() )
    << " raw CV/variation histograms to: " << root_file_name << '\n';

  return 0;
}
