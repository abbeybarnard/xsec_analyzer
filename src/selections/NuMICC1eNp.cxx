// XSecAnalyzer includes
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>

#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Array.hh"
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"

#include "XSecAnalyzer/Selections/NuMICC1eNp.hh"

// Initialising WC signal definition variables
namespace {

struct WcRecoCounts {
  int reco_electron = 0;
  int reco_proton = 0;
  int reco_mu = 0;
  int reco_pi0 = 0;
  int reco_pi = 0;
};

template< typename TMap >
bool read_track_value_as_float( const TMap& tree, const std::string& br_name,
  int idx, float& out_val ) {

  auto pick_from_vector = [ idx ]( const auto* vec_ptr, float& out ) -> bool {
    if ( vec_ptr == nullptr || vec_ptr->empty() ) return false;

    if ( 0 <= idx && idx < static_cast< int >( vec_ptr->size() ) ) {
      out = static_cast< float >( vec_ptr->at( idx ) );
      return true;
    }
    else {
      return false;
    }
  };

  try {
    std::vector< float >* vec_ptr = nullptr;
    tree.at( br_name ) >> vec_ptr;
    if ( pick_from_vector( vec_ptr, out_val ) ) return true;
  }
  catch ( const std::runtime_error& ) {}

  try {
    std::vector< double >* vec_ptr = nullptr;
    tree.at( br_name ) >> vec_ptr;
    if ( pick_from_vector( vec_ptr, out_val ) ) return true;
  }
  catch ( const std::runtime_error& ) {}

  try {
    std::vector< std::int32_t >* vec_ptr = nullptr;
    tree.at( br_name ) >> vec_ptr;
    if ( pick_from_vector( vec_ptr, out_val ) ) return true;
  }
  catch ( const std::runtime_error& ) {}

  try {
    std::vector< std::uint32_t >* vec_ptr = nullptr;
    tree.at( br_name ) >> vec_ptr;
    if ( pick_from_vector( vec_ptr, out_val ) ) return true;
  }
  catch ( const std::runtime_error& ) {}

  return false;
}

/////////////////////////////////////////////////////////////////////

// Choose whether to use Pandora or WC reco variables
// for the overlapping events
constexpr bool use_pandora_overlap_e = false; // Prefer WC for electron energy
constexpr bool use_pandora_overlap_angle = true; // Prefer Pandora for opening angle
constexpr bool use_pandora_overlap_lead_p_ke = false; // Prefer WC for leading proton KE (but WC track energy is not very good, so maybe prefer Pandora here?)

// Populating those WC signal definition variables
template< typename TMap >
WcRecoCounts populate_wc_reco_counters( const TMap& pf_eval ) {
  WcRecoCounts counts;

  ArrayView< int > reco_pdg;
  ArrayView< int > reco_mother;
  ArrayView< float > reco_start_momentum;

  pf_eval.at( "reco_pdg" ) >> reco_pdg;
  pf_eval.at( "reco_mother" ) >> reco_mother;
  pf_eval.at( "reco_startMomentum" ) >> reco_start_momentum;

  size_t ntrack = reco_pdg.size();
  for ( size_t i = 0u; i < ntrack; ++i ) {
    int pdg = reco_pdg[ i ];
    int mother = reco_mother[ i ];

    // Must come from neutrino
    if ( mother != 0 ) continue;

    // Skip photons + neutrons
    if ( pdg == 22 || pdg == 2112 ) continue;

    // Momentum magnitude (index 3 stores TOTAL energy)
    float energy = reco_start_momentum[ i ][ 3 ];

    if ( std::abs( pdg ) == 11 && energy - 0.000511 > 0.07 ) {
      ++counts.reco_electron;
    }
    else if ( pdg == 2212 && energy - PROTON_MASS > 0.04 ) {
      ++counts.reco_proton;
    }
    else if ( std::abs( pdg ) == 13 ) {
      ++counts.reco_mu;
    }
    else if ( pdg == 111 ) {
      ++counts.reco_pi0;
    }
    else if ( std::abs( pdg ) == 211 && energy - PI_PLUS_MASS > 0.04 ) {
      ++counts.reco_pi;
    }
  }

  return counts;
}

} // end anonymous namespace

NuMICC1eNp::NuMICC1eNp() : SelectionBase( "NuMICC1eNp" ) {
  // FV definition as in PeLEE analysis
  // x_min, x_max, y_min, y_max, z_min, z_max
  this->define_fv( 10., 246., -101., 101., 10., 986. );

  // Beam mode is passed per file by ProcessNTuples (from files_to_process)
  // via XSEC_ANALYZER_BEAM_MODE. Default to FHC if not set.
  const char* beam_mode_env = std::getenv( "XSEC_ANALYZER_BEAM_MODE" );
  beam_mode_ = beam_mode_env ? beam_mode_env : "FHC";

  // File type is passed per file by ProcessNTuples (from files_to_process)
  // via XSEC_ANALYZER_FILE_TYPE. Default to "unknown" if not set.
  const char* file_type_env = std::getenv( "XSEC_ANALYZER_FILE_TYPE" );
  file_type_ = file_type_env ? file_type_env : "unknown";
}

// This is more where the selection is.
bool NuMICC1eNp::is_selected( AnalysisEvent& ev ) {

  // Get access to the input Pandora and WC TTrees
  const auto& pandora = ev.in( "nuselection/NeutrinoSelectionFilter" );
  const auto& bdt_vars = ev.in( "wcpselection/T_BDTvars" );
  const auto& pf_eval = ev.in( "wcpselection/T_PFeval" );
  const auto& eval = ev.in( "wcpselection/T_eval" );

  // Initialising variables for Pandora fiducial volume coordinates
  float pandora_nu_vx, pandora_nu_vy, pandora_nu_vz;
  pandora.at( "reco_nu_vtx_sce_x" ) >> pandora_nu_vx;
  pandora.at( "reco_nu_vtx_sce_y" ) >> pandora_nu_vy;
  pandora.at( "reco_nu_vtx_sce_z" ) >> pandora_nu_vz;

  // Initialising variables for WC fiducial volume coordinates
  float wc_nu_vx, wc_nu_vy, wc_nu_vz;
  pf_eval.at( "reco_nuvtxX" ) >> wc_nu_vx;
  pf_eval.at( "reco_nuvtxY" ) >> wc_nu_vy;
  pf_eval.at( "reco_nuvtxZ" ) >> wc_nu_vz;

  // Initialising some other variables we need
  float pandora_shr_energy_cali, pandora_shr_energy_tot_cali;
  float pandora_trk_energy, pandora_tksh_angle;
  pandora.at( "shr_energy_cali" ) >> pandora_shr_energy_cali;
  pandora.at( "shr_energy_tot_cali" ) >> pandora_shr_energy_tot_cali;
  pandora.at( "trk_energy" ) >> pandora_trk_energy;
  pandora.at( "tksh_angle" ) >> pandora_tksh_angle;

  // -----------------------------------------------------------------------
  // TRUTH VARIABLES FOR RESOLUTION STUDIES
  //
  // These branches are only present in MC files. All reads are wrapped in a
  // single try/catch so that missing branches (e.g. extBNB) leave all truth
  // variables at BOGUS without crashing.
  //
  // mc_px/py/pz and mc_pdg are vector branches containing all final-state
  // particles. The leading proton is identified as the highest-momentum
  // proton (PDG 2212) in the list, which matches the definition used by
  // the ntuple producer for proton_e.
  // -----------------------------------------------------------------------
  float true_elec_e   = BOGUS;
  float true_proton_e = BOGUS;
  float elec_px = 0.f, elec_py = 0.f, elec_pz = 0.f;
  float mc_px_0 = 0.f, mc_py_0 = 0.f, mc_pz_0 = 0.f;
  bool has_truth_momentum = false;

  try {
    pandora.at( "elec_e" )   >> true_elec_e;
    pandora.at( "proton_e" ) >> true_proton_e;
    pandora.at( "elec_px" )  >> elec_px;
    pandora.at( "elec_py" )  >> elec_py;
    pandora.at( "elec_pz" )  >> elec_pz;

    // mc_px/py/pz are vector branches — find the leading proton as the
    // highest-momentum proton (PDG 2212) in the final state particle list
    std::vector< float >* mc_px_vec  = nullptr;
    std::vector< float >* mc_py_vec  = nullptr;
    std::vector< float >* mc_pz_vec  = nullptr;
    std::vector< int >*   mc_pdg_vec = nullptr;
    pandora.at( "mc_px" )  >> mc_px_vec;
    pandora.at( "mc_py" )  >> mc_py_vec;
    pandora.at( "mc_pz" )  >> mc_pz_vec;
    pandora.at( "mc_pdg" ) >> mc_pdg_vec;

    if ( mc_pdg_vec && mc_px_vec && mc_py_vec && mc_pz_vec
      && !mc_pdg_vec->empty() ) {

      // Find the highest-momentum proton
      float best_pmag2 = -1.f;
      int best_idx = -1;

      for ( int i = 0; i < static_cast<int>( mc_pdg_vec->size() ); ++i ) {
        if ( mc_pdg_vec->at( i ) != 2212 ) continue;
        if ( i >= static_cast<int>( mc_px_vec->size() ) ) continue;

        float px = mc_px_vec->at( i );
        float py = mc_py_vec->at( i );
        float pz = mc_pz_vec->at( i );
        float pmag2 = px*px + py*py + pz*pz;

        if ( pmag2 > best_pmag2 ) {
          best_pmag2 = pmag2;
          best_idx   = i;
        }
      }

      if ( best_idx >= 0 ) {
        mc_px_0 = mc_px_vec->at( best_idx );
        mc_py_0 = mc_py_vec->at( best_idx );
        mc_pz_0 = mc_pz_vec->at( best_idx );
        has_truth_momentum = true;
      }
    }
  }
  catch ( const std::exception& ) {}

  // Convert true proton total energy to kinetic energy by subtracting
  // the proton mass, matching what pandora_reco_track_energy and
  // wc_reco_track_energy store (both are KE, not total energy)
  float true_proton_ke = ( true_proton_e != static_cast<float>( BOGUS ) )
    ? true_proton_e - static_cast<float>( PROTON_MASS ) : BOGUS;

  // True cos(opening angle) between the true electron and leading proton
  double true_cos_opening_angle = BOGUS;
  if ( has_truth_momentum ) {
    float emag = std::sqrt( elec_px*elec_px + elec_py*elec_py + elec_pz*elec_pz );
    float pmag = std::sqrt( mc_px_0*mc_px_0 + mc_py_0*mc_py_0 + mc_pz_0*mc_pz_0 );
    if ( emag > 0.f && pmag > 0.f ) {
      float dot = elec_px*mc_px_0 + elec_py*mc_py_0 + elec_pz*mc_pz_0;
      float cos_theta = dot / ( emag * pmag );
      if ( cos_theta > 1.f ) cos_theta = 1.f;
      else if ( cos_theta < -1.f ) cos_theta = -1.f;
      true_cos_opening_angle = cos_theta;
    }
  }

  // Is the reco neutrino vertex in the fiducial volume?
  bool pandora_in_fv = this->get_fv().is_inside(
    pandora_nu_vx, pandora_nu_vy, pandora_nu_vz );
  bool wc_in_fv = this->get_fv().is_inside(
    wc_nu_vx, wc_nu_vy, wc_nu_vz );

  ////////// PANDORA SELECTION //////////

  // Read each Pandora selection variable explicitly into a local float so it
  // can be written to the output tree unconditionally. The formula strings
  // below still drive the actual cuts — these reads are independent.
  float var_swtrig_pre       = 9999.f;
  float var_slice_id         = 9999.f;
  float var_contained_frac   = 9999.f;
  float var_n_showers        = 9999.f;
  float var_n_tracks         = 9999.f;
  float var_shr_energy_tot   = 9999.f;
  float var_shrmoliereavg    = 9999.f;
  float var_tksh_distance    = 9999.f;
  try { pandora.at( "swtrig_pre" )          >> var_swtrig_pre;     } catch (...) {}
  try { pandora.at( "slice_orig_pass_id" )  >> var_slice_id;       } catch (...) {}
  try { pandora.at( "contained_fraction" )  >> var_contained_frac; } catch (...) {}
  try { pandora.at( "n_showers_contained" ) >> var_n_showers;      } catch (...) {}
  try { pandora.at( "n_tracks_contained" )  >> var_n_tracks;       } catch (...) {}
  try { pandora.at( "shr_energy_tot_cali" ) >> var_shr_energy_tot; } catch (...) {}
  try { pandora.at( "shrmoliereavg" )       >> var_shrmoliereavg;  } catch (...) {}
  try { pandora.at( "tksh_distance" )       >> var_tksh_distance;  } catch (...) {}

  static const std::string quality_cuts(
    "swtrig_pre==1 && slice_orig_pass_id==1 && contained_fraction > 0.9" );

  static const std::string signal_definition_cuts(
    "n_showers_contained==1 && n_tracks_contained > 0 &&"
    " shr_energy_tot_cali > 0.07 && trk_energy > 0.04" );

  // Properly reject sentinel values for shower Molière angle 
  static const std::string bdt_loose_cuts_scalar(
    "shrmoliereavg < 15 && shrmoliereavg > -3.4e37 && tksh_distance < 12" );

  // BDT score is precomputed in scripts/compute_bdt_scores.py and stored
  // in the Pandora tree as "bdt_score".
  float bdt_score = 0.f;
  try {
    pandora.at( "bdt_score" ) >> bdt_score;
  }
  catch ( const std::exception& ) {
    // If preprocessing has not been run, keep default score and fail cut.
    bdt_score = 0.f;
  }

  const bool beam_mode_is_rhc = ( beam_mode_ == "RHC" );
  const float bdt_score_cut = beam_mode_is_rhc ? 0.500f : 0.475f;

  // Reproduce the Python trksemlbl logic by taking the value associated
  // with the leading-track index (trk_id is one-based in the ntuple).
  float trksemlbl = 9999.f;
  std::uint32_t trk_id_u = 0u;
  pandora.at( "trk_id" ) >> trk_id_u;
  const int trk_idx = static_cast<int>( trk_id_u ) - 1;

  if ( trk_idx < 0 ) {
    trksemlbl = 9999.f;
  } else {
    bool semlbl_ok = read_track_value_as_float( pandora, "pfng2semlabel", trk_idx, trksemlbl );
    if ( !semlbl_ok ) trksemlbl = 9999.f;
  }

  // PANDORA QUALITY CUTS
  bool pandora_pass_quality = pandora_in_fv
    && pandora.formula( quality_cuts );

  // PANDORA SIGNAL DEFINITION CONSTRAINTS
  bool pandora_pass_sig = pandora.formula( signal_definition_cuts );

  // PANDORA PRESELECTION = quality cuts + signal definition constraints
  bool pandora_pass_pre = pandora_pass_quality
    && pandora_pass_sig;

  // PANDORA BDT LOOSE CUTS
  bool pandora_pass_bdt_loose = pandora_pass_pre
    && pandora.formula( bdt_loose_cuts_scalar )
    && ( trksemlbl == 1.f );

  // Beam-mode-dependent BDT score requirement
  bool pandora_pass_bdt_score = ( bdt_score > bdt_score_cut );

  // Final Pandora selection flag
  bool pandora_sel = pandora_pass_bdt_loose
    && pandora_pass_bdt_score;

  ////////// CUTFLOW COUNTERS //////////

  ++n_total_;
  if ( pandora_pass_quality ) ++n_pass_quality_;
  if ( pandora_pass_sig ) ++n_pass_sig_;
  if ( pandora_pass_pre ) ++n_pass_pre_;
  if ( pandora_pass_bdt_loose ) ++n_pass_bdt_loose_;
  if ( pandora_pass_bdt_score ) ++n_pass_bdt_score_;
  if ( pandora_sel ) ++n_pass_final_;

  ////////// WC SELECTION //////////

  // Read WC selection variables explicitly for output.
  float var_wc_numu_cc_flag = 9999.f;
  float var_wc_match_isFC   = 9999.f;
  float var_wc_nue_score    = 9999.f;
  try { bdt_vars.at( "numu_cc_flag" ) >> var_wc_numu_cc_flag; } catch (...) {}
  try { eval.at( "match_isFC" )       >> var_wc_match_isFC;   } catch (...) {}
  try { bdt_vars.at( "nue_score" )    >> var_wc_nue_score;    } catch (...) {}

  // WC QUALITY CUTS
  static const std::string wc_quality_cuts_bdt_vars( "numu_cc_flag >= 0" );
  static const std::string wc_quality_cuts_eval( "match_isFC == 1" );

  bool wc_pass_quality_bdt
    = bdt_vars.formula( wc_quality_cuts_bdt_vars );
  bool wc_pass_quality_eval
    = eval.formula( wc_quality_cuts_eval );
  bool wc_pass_quality = wc_in_fv
    && wc_pass_quality_bdt
    && wc_pass_quality_eval;

  // WC SIGNAL DEFINITION CONSTRAINTS
  WcRecoCounts wc_reco_counts = populate_wc_reco_counters( pf_eval );

  bool wc_pass_sig_pf = ( wc_reco_counts.reco_electron == 1 )
    && ( wc_reco_counts.reco_proton > 0 ) && ( wc_reco_counts.reco_mu == 0 )
    && ( wc_reco_counts.reco_pi == 0 ) && ( wc_reco_counts.reco_pi0 == 0 );

  bool wc_pass_sig = wc_pass_sig_pf;

  // WC PRESELECTION = quality cuts + signal-definition constraints
  bool wc_pass_pre = wc_pass_quality && wc_pass_sig;

  // WC BDT SCORE
  bool wc_pass_bdt = bdt_vars.formula( "nue_score > 7" );

  bool wc_sel = wc_pass_pre && wc_pass_bdt;

  ////////// UNION SELECTION //////////

  bool sel_nue_cc_union = pandora_sel || wc_sel;
  bool sel_nue_cc_intersection = pandora_sel && wc_sel;
  bool sel_nue_cc = sel_nue_cc_union;

  ////////// VARIABLE INITIALIZATION //////////

  double reco_electron_energy = BOGUS;
  bool reco_electron_energy_is_pandora = false;

  double reco_track_energy = BOGUS;
  bool reco_track_energy_is_pandora = false;

  double reco_opening_angle = BOGUS;
  bool reco_opening_angle_is_pandora = false;

  ////////// CORRECTIONS //////////

  // Pandora electron energy: shr_energy_cali / 0.83 applied for all events
  // (matches Python notebook which applies this scaling to all MC + data)
  const double pandora_reco_electron_energy = pandora_sel
    ? pandora_shr_energy_cali / 0.83 : BOGUS;

  // WireCell shower KE: read for all events regardless of file type
  // so it is available for MC plots.
  // The 0.95 data correction is applied only for real data (onBNB/extBNB),
  // matching the Python notebook which scales data/ext but not MC.
  bool is_real_data = ( file_type_ == "onBNB" || file_type_ == "extBNB" );

  float wc_reco_shower_ke = BOGUS;
  try {
    pf_eval.at( "reco_showerKE" ) >> wc_reco_shower_ke;
  }
  catch ( const std::exception& ) {}

  double wc_reco_electron_energy = BOGUS;
  if ( wc_sel ) {
    wc_reco_electron_energy = is_real_data
      ? wc_reco_shower_ke * 0.95
      : wc_reco_shower_ke;
  }

  // Pandora track energy
  const double pandora_reco_track_energy = pandora_sel
    ? pandora_trk_energy : BOGUS;

  // WireCell proton and shower momenta — read for all events (MC and data)
  const bool has_reco_proton_momentum = pf_eval.formula(
    "reco_protonMomentum[3] > 0" );
  const bool has_reco_shower_momentum = pf_eval.formula(
    "reco_showerMomentum[3] > 0" );

  double wc_reco_track_energy = BOGUS;
  if ( wc_sel && has_reco_proton_momentum ) {
    ArrayView< float > wc_reco_proton_momentum;
    pf_eval.at( "reco_protonMomentum" ) >> wc_reco_proton_momentum;
    if ( wc_reco_proton_momentum.size() > 3 ) {
      float wc_proton_total_energy = wc_reco_proton_momentum[ 3 ];
      wc_reco_track_energy = wc_proton_total_energy - PROTON_MASS;
    }
  }

  // WireCell opening angle — read for all events (MC and data)
  double wc_reco_opening_angle = BOGUS;
  if ( wc_sel && has_reco_shower_momentum && has_reco_proton_momentum ) {
    ArrayView< float > wc_reco_shower_momentum;
    ArrayView< float > wc_reco_proton_momentum;
    pf_eval.at( "reco_showerMomentum" ) >> wc_reco_shower_momentum;
    pf_eval.at( "reco_protonMomentum" ) >> wc_reco_proton_momentum;

    if ( wc_reco_shower_momentum.size() > 3 && wc_reco_proton_momentum.size() > 3 ) {
      float shower_px = wc_reco_shower_momentum[ 0 ];
      float shower_py = wc_reco_shower_momentum[ 1 ];
      float shower_pz = wc_reco_shower_momentum[ 2 ];

      float proton_px = wc_reco_proton_momentum[ 0 ];
      float proton_py = wc_reco_proton_momentum[ 1 ];
      float proton_pz = wc_reco_proton_momentum[ 2 ];

      float dot_product = shower_px * proton_px + shower_py * proton_py
        + shower_pz * proton_pz;

      float shower_mag = std::sqrt( shower_px * shower_px
        + shower_py * shower_py + shower_pz * shower_pz );
      float proton_mag = std::sqrt( proton_px * proton_px
        + proton_py * proton_py + proton_pz * proton_pz );

      if ( shower_mag > 0 && proton_mag > 0 ) {
        float cos_theta = dot_product / ( shower_mag * proton_mag );
        if ( cos_theta > 1.f ) cos_theta = 1.f;
        else if ( cos_theta < -1.f ) cos_theta = -1.f;
        wc_reco_opening_angle = cos_theta;
      }
    }
  }

  ////////// RECO VARIABLE CHOICE //////////

  const bool use_pandora_e     = pandora_sel && ( !wc_sel || use_pandora_overlap_e );
  const bool use_pandora_p_ke  = pandora_sel && ( !wc_sel || use_pandora_overlap_lead_p_ke );
  const bool use_pandora_angle = pandora_sel && ( !wc_sel || use_pandora_overlap_angle );

  const double pandora_reco_opening_angle = pandora_sel ? pandora_tksh_angle : BOGUS;

  // ELECTRON ENERGY
  if ( use_pandora_e ) {
    reco_electron_energy = pandora_reco_electron_energy;
    reco_electron_energy_is_pandora = true;
  } else if ( wc_sel ) {
    reco_electron_energy = wc_reco_electron_energy;
  }

  // LEADING PROTON KINETIC ENERGY
  if ( use_pandora_p_ke ) {
    reco_track_energy = pandora_reco_track_energy;
    reco_track_energy_is_pandora = true;
  } else if ( wc_sel ) {
    reco_track_energy = wc_reco_track_energy;
  }

  // OPENING ANGLE
  if ( use_pandora_angle ) {
    reco_opening_angle = pandora_reco_opening_angle;
    reco_opening_angle_is_pandora = true;
  } else if ( wc_sel ) {
    reco_opening_angle = wc_reco_opening_angle;
  }

  // We're done. Store the results in the output tree and return whether the
  // event passed the full selection
  auto& out = ev.out();
  int run = BOGUS_INDEX;
  int sub = BOGUS_INDEX;
  int evt = BOGUS_INDEX;
  pandora.at( "run" ) >> run;
  pandora.at( "sub" ) >> sub;
  pandora.at( "evt" ) >> evt;
  out[ "run" ] = run;
  out[ "sub" ] = sub;
  out[ "evt" ] = evt;

  // Pandora selection outputs
  out[ "pandora_pass_quality" ] = pandora_pass_quality;
  out[ "pandora_pass_pre" ] = pandora_pass_pre;
  out[ "pandora_pass_sig" ] = pandora_pass_sig;
  out[ "pandora_pass_bdt_loose" ] = pandora_pass_bdt_loose;
  out[ "pandora_pass_bdt_score" ] = pandora_pass_bdt_score;
  out[ "beam_mode_is_rhc" ] = beam_mode_is_rhc;
  out[ "pandora_bdt_score_cut" ] = bdt_score_cut;
  out[ "pandora_sel_nue_cc" ] = pandora_sel;

  // WC selection outputs
  out[ "wc_pass_quality" ] = wc_pass_quality;
  out[ "wc_pass_sig" ] = wc_pass_sig;
  out[ "wc_pass_pre" ] = wc_pass_pre;
  out[ "wc_pass_bdt" ] = wc_pass_bdt;
  out[ "wc_sel_nue_cc" ] = wc_sel;

  // Combined selection outputs (Pandora OR WC)
  out[ "sel_nue_cc_union" ] = sel_nue_cc_union;
  out[ "sel_nue_cc_intersection" ] = sel_nue_cc_intersection;

  // Main selection variable (final selection)
  out[ "sel_nue_cc" ] = sel_nue_cc;

  // BDT score
  out[ "bdt_score" ] = bdt_score;

  // Electron energy
  out[ "pandora_reco_electron_energy" ] = pandora_reco_electron_energy;
  out[ "wc_reco_electron_energy" ] = wc_reco_electron_energy;
  out[ "wc_reco_shower_ke" ] = wc_reco_shower_ke;
  out[ "reco_electron_energy" ] = reco_electron_energy;

  // Leading proton kinetic energy
  out[ "pandora_reco_track_energy" ] = pandora_reco_track_energy;
  out[ "wc_reco_track_energy" ] = wc_reco_track_energy;
  out[ "reco_track_energy" ] = reco_track_energy;

  // Opening angle between electron and leading proton
  out[ "pandora_reco_opening_angle" ] = pandora_reco_opening_angle;
  out[ "wc_reco_opening_angle" ] = wc_reco_opening_angle;
  out[ "reco_opening_angle" ] = reco_opening_angle;

  // Flags for whether the chosen reco variable came from Pandora or WC
  out[ "reco_electron_energy_is_pandora" ] = reco_electron_energy_is_pandora;
  out[ "reco_track_energy_is_pandora" ] = reco_track_energy_is_pandora;
  out[ "reco_opening_angle_is_pandora" ] = reco_opening_angle_is_pandora;

  // Truth variables for resolution studies
  // BOGUS for extBNB and any file where truth branches are absent
  out[ "true_elec_e" ] = true_elec_e;
  out[ "true_proton_ke" ] = true_proton_ke;
  out[ "true_cos_opening_angle" ] = true_cos_opening_angle;

  std::string cat = categorize_event(ev);
  out["NuMICC1eNp_Cats"] = cat;
  out["is_true_signal"] = (cat == "#nu_{e} CC0#piNp");

  // Integer event category for UniverseMaker compatibility.
  // Values match the order in configs/selections.conf for NuMICC1eNp.
  static const std::map<std::string, int> cat_to_int = {
    { "#nu_{e} CC0#piNp",       1 },
    { "#bar{#nu}_{e} CC0#piNp", 2 },
    { "#bar{#nu}_{e} CC Other", 3 },
    { "#nu_{e} CC Other",       4 },
    { "#nu_{#mu} CC #pi^{0}",   5 },
    { "#nu_{#mu} CC Other",     6 },
    { "NC #pi^{0}",             7 },
    { "NC Other",               8 },
    { "Out FV",                 9 },
    { "Unknown",               10 },
  };
  auto cat_it = cat_to_int.find( cat );
  int cat_int = ( cat_it != cat_to_int.end() ) ? cat_it->second : 0;
  out["NuMICC1eNp_EventCategory"] = cat_int;

  // Branches required by UniverseMaker to correctly apply CV weights.
  // is_mc is read from the environment variable set by ProcessNTuples.
  const char* file_type_env = std::getenv( "XSEC_ANALYZER_FILE_TYPE" );
  std::string file_type_str = file_type_env ? file_type_env : "unknown";
  bool is_mc_flag = ( file_type_str != "onBNB" && file_type_str != "extBNB" );
  out["is_mc"] = is_mc_flag;

  // GENIE tune CV weight
  float tuned_cv_weight = 1.f;
  if ( is_mc_flag ) {
    try { pandora.at( "weightTune" ) >> tuned_cv_weight; }
    catch (...) { tuned_cv_weight = 1.f; }
  }
  out["tuned_cv_weight"] = tuned_cv_weight;

  // PPFX CV weight (NuMI-specific)
  float ppfx_cv_weight = 1.f;
  if ( is_mc_flag ) {
    try { pandora.at( "ppfx_cv" ) >> ppfx_cv_weight; }
    catch (...) { ppfx_cv_weight = 1.f; }
  }
  out["ppfx_cv_weight"] = ppfx_cv_weight;

  // Normalisation weight (0.65 for dirt MC, 1.0 otherwise)
  float normalisation_weight = 1.f;
  if ( file_type_str == "dirtMC" ) normalisation_weight = 0.65f;
  out["normalisation_weight"] = normalisation_weight;

  // Copy only the weight branches listed in systcalc_numi.conf into the
  // output tree with "weight_" prefix. Writing *all* weights from the map
  // causes UniverseMaker to throw on unrecognised names (e.g. the PPFX
  // sub-component branches such as weight_ppfx_mippk_PPFXMIPPKaon that
  // are present in the input ntuples but not defined in the config).
  //
  // Three additional branches are required by UniverseMaker regardless of
  // systcalc_numi.conf:
  //   weight_TunedCentralValue_UBGenie -- checked by is_reweightable_mc_ntuple()
  //                                       to determine whether to build universes
  //   weight_splines_general_Spline    -- spline CV correction (SPLINE_WEIGHT_NAME)
  //   weight_ppfx_cv_UBPPFXCV          -- PPFX CV weight universe (PPFX_WEIGHT_NAME)
  if ( is_mc_flag ) {
    static const std::set<std::string> allowed_weights = {
      // Required by UniverseMaker internally (TUNE/SPLINE/PPFX_WEIGHT_NAME)
      "TunedCentralValue_UBGenie",
      "splines_general_Spline",
      "ppfx_cv_UBPPFXCV",
      // Flux: PPFX total
      "ppfx_all",
      // Flux: beamline geometry
      "Horn_2kA",
      "Horn1_x_3mm",   "Horn1_y_3mm",
      "Beam_spot_1_1mm", "Beam_spot_1_5mm",
      "Horn2_x_3mm",   "Horn2_y_3mm",
      "Horns_0mm_water", "Horns_2mm_water",
      "Beam_shift_x_1mm", "Beam_shift_y_1mm",
      "Target_z_7mm",
      // Reinteraction
      "reint_all",
      // GENIE cross-section: multisim
      "All_UBGenie",
      // GENIE cross-section: unisims
      "AxFFCCQEshape_UBGenie", "DecayAngMEC_UBGenie",
      "NormCCCOH_UBGenie",     "NormNCCOH_UBGenie",
      "RPA_CCQE_UBGenie",      "ThetaDelta2NRad_UBGenie",
      "Theta_Delta2Npi_UBGenie", "VecFFCCQEshape_UBGenie",
      "XSecShape_CCMEC_UBGenie",
      // SCC
      "xsr_scc_Fa3_SCC", "xsr_scc_Fv3_SCC",
    };

    std::map<std::string, std::vector<double>>* wm = nullptr;
    try {
      pandora.at( "weights" ) >> wm;
      if ( wm ) {
        for ( const auto& wgt_pair : *wm ) {
          if ( allowed_weights.count( wgt_pair.first ) ) {
            out[ "weight_" + wgt_pair.first ] = wgt_pair.second;
          }
        }
      }
    }
    catch (...) {}
  }

  // -----------------------------------------------------------------------
  // RAW SELECTION VARIABLES
  //
  // All raw inputs to the selection, written unconditionally for every event
  // so they can be plotted at any cut stage in the notebook. They are never
  // cut on here — use pandora_pass_* / wc_pass_* as stage gates when plotting.
  // Sentinel 9999 means the branch was absent or reconstruction failed.
  //
  // In the notebook: tree.keys() filtered by startswith("var_") gives the
  // full list.
  // -----------------------------------------------------------------------

  // Pandora quality cuts
  out[ "var_pandora_swtrig_pre" ]     = var_swtrig_pre;
  out[ "var_pandora_slice_id" ]       = var_slice_id;
  out[ "var_pandora_contained_frac" ] = var_contained_frac;
  // Pandora signal definition cuts
  out[ "var_pandora_n_showers" ]      = var_n_showers;
  out[ "var_pandora_n_tracks" ]       = var_n_tracks;
  out[ "var_pandora_shr_energy_tot" ] = var_shr_energy_tot;
  out[ "var_pandora_trk_energy" ]     = pandora_trk_energy;
  // Pandora BDT loose cuts
  out[ "var_pandora_shrmoliereavg" ]  = var_shrmoliereavg;
  out[ "var_pandora_tksh_distance" ]  = var_tksh_distance;
  out[ "var_pandora_trksemlbl" ]      = trksemlbl;
  // Pandora BDT score
  out[ "var_pandora_bdt_score" ]      = bdt_score;
  // WC quality cuts
  out[ "var_wc_numu_cc_flag" ]        = var_wc_numu_cc_flag;
  out[ "var_wc_match_isFC" ]          = var_wc_match_isFC;
  // WC signal definition counts
  out[ "var_wc_reco_n_electron" ]     = wc_reco_counts.reco_electron;
  out[ "var_wc_reco_n_proton" ]       = wc_reco_counts.reco_proton;
  out[ "var_wc_reco_n_mu" ]           = wc_reco_counts.reco_mu;
  out[ "var_wc_reco_n_pi" ]           = wc_reco_counts.reco_pi;
  out[ "var_wc_reco_n_pi0" ]          = wc_reco_counts.reco_pi0;
  // WC BDT score
  out[ "var_wc_nue_score" ]           = var_wc_nue_score;

  return sel_nue_cc;
}

//////////////////////////////////////////////////////////////////////////////////

std::string NuMICC1eNp::categorize_event( AnalysisEvent& ev ) {

  int nu_pdg, ccnc, nproton, npion, npi0;
  float nu_x, nu_y, nu_z, Ee;

  const auto& in = ev.in();

  in.at( "nu_pdg" ) >> nu_pdg;
  in.at( "ccnc" ) >> ccnc;
  in.at( "nproton" ) >> nproton;
  in.at( "npion" ) >> npion;
  in.at( "npi0" ) >> npi0;
  in.at( "elec_e" ) >> Ee;
  in.at( "true_nu_vtx_x" ) >> nu_x;
  in.at( "true_nu_vtx_y" ) >> nu_y;
  in.at( "true_nu_vtx_z" ) >> nu_z;

  bool sig_inFV = this->get_fv().is_inside( nu_x, nu_y, nu_z );
  bool sig_isNuE = (nu_pdg == ELECTRON_NEUTRINO);
  bool sig_isCC = (ccnc == CHARGED_CURRENT);

  bool is_signal = sig_inFV && sig_isNuE && sig_isCC
                && (Ee > 0.070)
                && (nproton > 0)
                && (npion == 0)
                && (npi0 == 0);

  if ( !sig_inFV ) {
    return "Out FV";
  }
  else if ( !sig_isCC ) {
    if ( npi0 > 0 ) return "NC #pi^{0}";
    else return "NC Other";
  }
  else if ( nu_pdg == ELECTRON_NEUTRINO ) {
    if ( is_signal ) return "#nu_{e} CC0#piNp";
    else return "#nu_{e} CC Other";
  }
  else if ( nu_pdg == ELECTRON_ANTINEUTRINO ) {
    if ( nproton > 0 && npion == 0 && npi0 == 0 ) return "#bar{#nu}_{e} CC0#piNp";
    else return "#bar{#nu}_{e} CC Other";
  }
  else if ( std::abs(nu_pdg) == MUON_NEUTRINO ) {
    if ( npi0 > 0 ) return "#nu_{#mu} CC #pi^{0}";
    else return "#nu_{#mu} CC Other";
  }

  std::cout << "Warning: Unknown event! Check the categorization logic.\n";
  return "Unknown";
}

//////////////////////////////////////////////////////////////////////////////////

NuMICC1eNp::~NuMICC1eNp() {
  std::cout << "\n==============================\n";
  std::cout << "  PANDORA CUTFLOW (final)\n";
  std::cout << "==============================\n";
  std::cout << "  total events:    " << n_total_ << "\n"
            << "  pass quality:    " << n_pass_quality_ << "\n"
            << "  pass sig:        " << n_pass_sig_ << "\n"
            << "  pass pre:        " << n_pass_pre_ << "\n"
            << "  pass bdt_loose:  " << n_pass_bdt_loose_ << "\n"
            << "  pass bdt_score:  " << n_pass_bdt_score_ << "\n"
            << "  pass final:      " << n_pass_final_ << "\n";
}