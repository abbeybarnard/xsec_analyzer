// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"

#include "XSecAnalyzer/Selections/NuMICC1eNp.hh"

NuMICC1eNp::NuMICC1eNp() : SelectionBase( "NuMICC1eNp" ) {
  // FV definition as in PeLEE analysis
  // x_min, x_max, y_min, y_max, z_min, z_max
  // Have matched this up with what Katrina used.
  this->define_fv( 10., 246., -106., 106., 10., 1026. );
}

// This is the signal definition.
std::string NuMICC1eNp::categorize_event( AnalysisEvent& ev ) {

  // Remember to add things here! These are the variables we use.
  int nu_pdg, ccnc, npi0, nelec;
  float nu_x, nu_y, nu_z, Ee;
  const auto& in = ev.in();
  in.at( "nu_pdg" ) >> nu_pdg;
  in.at( "ccnc" ) >> ccnc;
  in.at( "npi0" ) >> npi0;
  in.at( "nelec" ) >> nelec;
  in.at( "elec_e" ) >> Ee;
  in.at( "true_nu_vtx_x" ) >> nu_x;
  in.at( "true_nu_vtx_y" ) >> nu_y;
  in.at( "true_nu_vtx_z" ) >> nu_z;

  // This is the proton selection stuff
  std::vector< int >* nu_daughter_pdg;
  std::vector< float > *nu_daughter_energy, *nu_daughter_px,
    *nu_daughter_py, *nu_daughter_pz;

  in.at( "mc_pdg" ) >> nu_daughter_pdg;
  in.at( "mc_E" ) >> nu_daughter_energy;
  in.at( "mc_px" ) >> nu_daughter_px;
  in.at( "mc_py" ) >> nu_daughter_py;
  in.at( "mc_pz" ) >> nu_daughter_pz;

  int num_p_in_energy_range = 0;
  double energy_lead_p = 0.;
  bool has_pions = false;

  for ( size_t p = 0u; p < nu_daughter_pdg->size(); ++p ) {
    int pdg = nu_daughter_pdg->at( p );
    float energy = nu_daughter_energy->at( p );

    // Requires events to have at least one proton with energy > 40 MeV
    if ( pdg == PROTON ) {
      if ( energy - PROTON_MASS > 0.040 ) {
        ++num_p_in_energy_range;
        if ( energy > energy_lead_p ) {
          energy_lead_p = energy;
        }
      }
    }

    // Requires events to have no neutral pions
    else if ( pdg == PI_ZERO ) {
      has_pions = true;
    }
    
    // Requires events to have no charged pions with energy > 40 MeV
    else if ( std::abs(pdg) == PI_MINUS ) {
      if ( energy - PI_MINUS_MASS > 0.040 ) {
        has_pions = true;
      }
    }

    // Requires events to have no charged pions with energy > 40 MeV
    else if ( std::abs(pdg) == PI_PLUS ) {
      if ( energy - PI_PLUS_MASS > 0.040 ) {
        has_pions = true;
      }
    }
  }
  
  // Require signal events to be inside the true fiducial volume
  bool sig_inFV = this->get_fv().is_inside( nu_x, nu_y, nu_z );

  // Require an incident electron neutrino - no electron antineutrinos here!
  bool sig_isNuE = (nu_pdg == ELECTRON_NEUTRINO);

  // Require a charged-current interaction
  bool sig_isCC = ( ccnc == CHARGED_CURRENT );

  // Require a final-state electron above threshold
  bool sig_has_fs_electron = ( nelec > 0 ); // 30 MeV threshold

  // This is the total signal definition with everything in it - originally missed the proton and pion cuts
  bool is_signal = sig_inFV && sig_isNuE && sig_isCC && sig_has_fs_electron && (num_p_in_energy_range > 0) && (!has_pions);

  // Evaluate the true kinematic variables of interest
  float mc_electron_energy = BOGUS;
  // Check if there is a true final-state electron in this event. If there is
  // one, then store its energy
  bool has_true_electron = ( sig_isNuE && sig_isCC );
  if ( has_true_electron ) mc_electron_energy = Ee;

  // Save TRUTH information to the output TTree
  auto& out = ev.out();
  out[ "mc_is_nue" ] = sig_isNuE;
  out[ "mc_is_CC" ] = sig_isCC;
  out[ "mc_vertex_in_FV" ] = sig_inFV;
  out[ "mc_has_fs_electron" ] = sig_has_fs_electron;
  out[ "mc_is_signal" ] = is_signal; // This is true total signal!
  out[ "num_p_in_energy_range" ] = num_p_in_energy_range;
  out[ "has_pions" ] = has_pions;
  out[ "mc_electron_energy" ] = mc_electron_energy;
  out[ "energy_lead_p" ] = energy_lead_p;

  // This is where we categorize the events!
  // All events outside of the true fiducial volume should be categorized
  // as "out of fiducial volume"
  
  // Events that are not within the FV
  if ( !sig_inFV ) return "Out FV";

  // The two NC categories
  else if ( !sig_isCC ) {
    if ( npi0 > 0 ) return "NC #pi^{0}";
    else return "NC Other";
  }

  // CC muon (anti)neutrinos
  else if ( std::abs( nu_pdg ) == MUON_NEUTRINO ) {
    if ( npi0 > 0 ) return "#nu_{#mu} CC #pi^{0}";
    else return "#nu_{#mu} CC Other";
  }

  // CC electron antineutrinos
  else if ( nu_pdg  == ELECTRON_ANTINEUTRINO ) {
    if ( npi0 > 0 ) return "#bar{#nu}_{e} CC0#piNp";
    else return "#bar{#nu}_{e} CC Other";
  }

  // CC electron neutrinos - where our signal events live
  else if ( sig_isNuE ) {
    // signal events
    if ( is_signal ) return "#nu_{e} CC0#piNp";
    // non-signal nues
    else return "#nu_{e} CC Other";
  }
 
  // We shouldn't ever get here, but return "Unknown" just in case
  std::cout << "Warning: Unknown event! Check the categorization logic.\n";
  return "Unknown";
}

// This is more where the selection is.
bool NuMICC1eNp::is_selected( AnalysisEvent& ev ) {

  // Get access to the reco information needed to apply the selection
  const auto& in = ev.in();
  unsigned int shr_id, n_tracks_contained, n_showers_contained; // This is the one that is "i" compared to "I"
  int nslice, n_showers, n_tracks, swtrig_pre;
  float nu_vx, nu_vy, nu_vz, contained_frac, topo_score,
    cosmic_ip, shr_energy_cali, shr_energy_tot_cali, shr_score, hits_ratio, shrmoliereavg,
    shr_tkfit_gap10_dedx_Y, shr_distance, trkpid, shr_tkfit_dedx_Y, tksh_distance, trk_energy;

  in.at( "shr_id" ) >> shr_id;
  in.at( "nslice" ) >> nslice;
  in.at( "n_showers" ) >> n_showers;
  in.at( "n_tracks" ) >> n_tracks;
  in.at( "n_tracks_contained" ) >> n_tracks_contained;
  in.at( "reco_nu_vtx_sce_x" ) >> nu_vx;
  in.at( "reco_nu_vtx_sce_y" ) >> nu_vy;
  in.at( "reco_nu_vtx_sce_z" ) >> nu_vz;
  in.at( "contained_fraction" ) >> contained_frac;
  in.at( "topological_score" ) >> topo_score;
  in.at( "CosmicIP" ) >> cosmic_ip;
  in.at( "shr_energy_cali" ) >> shr_energy_cali;
  in.at( "shr_energy_tot_cali" ) >> shr_energy_tot_cali;
  in.at( "shr_score" ) >> shr_score;
  in.at( "hits_ratio" ) >> hits_ratio;
  in.at( "shrmoliereavg" ) >> shrmoliereavg;
  in.at( "shr_tkfit_gap10_dedx_Y" ) >> shr_tkfit_gap10_dedx_Y;
  in.at( "shr_tkfit_dedx_Y" ) >> shr_tkfit_dedx_Y;
  in.at( "shr_distance" ) >> shr_distance;
  in.at( "tksh_distance" ) >> tksh_distance;
  in.at( "trk_energy" ) >> trk_energy;
  in.at( "trkpid" ) >> trkpid;
  in.at( "n_showers_contained" ) >> n_showers_contained;
  in.at( "swtrig_pre" ) >> swtrig_pre;

  std::vector< unsigned int >* gen_vec;
  in.at( "pfp_generation_v" ) >> gen_vec;

  // PRE-SELECTION (signal definition constraints and quality cuts)
  // passes software trigger
  bool passes_software_trigger = ( swtrig_pre==1 );
  // neutrino slice
  bool has_nu_slice = ( nslice == 1 );
  // vertex inside FV
  bool in_fv = this->get_fv().is_inside( nu_vx, nu_vy, nu_vz );
  // contained fraction
  bool contained_cut_ok = ( contained_frac > 0.9 );
  // has showee
  bool has_shower = ( n_showers_contained == 1 );
  // has contained tracks
  bool has_contained_tracks = ( n_tracks_contained > 0 );
  // track energy
  bool trk_energy_ok = ( trk_energy > 0.04 ); // GeV

  bool sel_pass_preselection = passes_software_trigger
    && has_nu_slice && in_fv && contained_cut_ok
    && has_shower && has_contained_tracks && trk_energy_ok;

  // COSMIC REJECTION
  // topological score
  bool topo_ok = ( topo_score >= 0.2 );
  // cosmic impact parameter
  bool cip_ok = ( cosmic_ip >= 10 );

  bool sel_pass_cosmic_rejection = topo_ok && cip_ok;

  // OTHER
  bool shower_score_ok = ( shr_score < 0.125 );
  bool shrmoliereavg_ok = ( shrmoliereavg < 8. );
  bool trkpid_ok = ( trkpid < 0. );
  bool shr_trkfit_dedx_Y_ok = ( shr_tkfit_dedx_Y < 4. );
  bool tksh_distance_ok = ( tksh_distance < 5. );

  bool sel_pass_other = shower_score_ok && shrmoliereavg_ok
    && trkpid_ok && shr_trkfit_dedx_Y_ok && tksh_distance_ok;

  // // SHOWER IDENTIFICATION
  // // valid shower ID
  // bool valid_ID = ( shr_id != 0 ); // zero is the default value (not filled)
  // // shower pfp generation
  // // NOTE: The second expression after the && is only evaluated if we have
  // // a valid shower ID, thus avoiding any problems with an invalid
  // // argument to the std::vector::at() function.
  // bool second_generation = valid_ID && gen_vec->at( shr_id - 1 ) == 2;
  // // shower energy (threshold matches the signal definition)
  // bool energy_ok = ( ( shr_energy_cali / 0.83 ) >= 0.03 );
  // // shower score
  // bool score_ok = ( shr_score <= 0.15 );
  // // shower hits ratio
  // bool hits_ratio_ok = ( hits_ratio >= 0.5 );

  // bool sel_pass_shower_identification = valid_ID && second_generation
  //   && energy_ok && score_ok && hits_ratio_ok;

  // // ELECTRON IDENTIFICATION
  // // moliere average angle
  // bool moliere_ok = ( shrmoliereavg <= 7. );
  // // shower distance and dE/dx (default to passing the cut)
  // bool dist_and_dEdx_ok = true;
  // if ( n_tracks > 0 ) {
  //   // track present, 2D distance-dE/dx cut
  //   if ( shr_tkfit_gap10_dedx_Y >= 0. && shr_tkfit_gap10_dedx_Y < 1.75 ) {
  //     if ( shr_distance > 3.0 ) dist_and_dEdx_ok = false;
  //   }
  //   else if ( shr_tkfit_gap10_dedx_Y >= 1.75 && shr_tkfit_gap10_dedx_Y < 2.5 ) {
  //     if ( shr_distance > 12.0 ) dist_and_dEdx_ok = false;
  //   }
  //   else if ( shr_tkfit_gap10_dedx_Y >= 2.5 && shr_tkfit_gap10_dedx_Y < 3.5 ) {
  //     if ( shr_distance > 3.0 ) dist_and_dEdx_ok = false;
  //   }
  //   else if ( shr_tkfit_gap10_dedx_Y >= 3.5 && shr_tkfit_gap10_dedx_Y < 4.7 ) {
  //     dist_and_dEdx_ok = false;
  //   }
  //   else if ( shr_tkfit_gap10_dedx_Y >= 4.7 ) {
  //     if ( shr_distance > 3.0 ) dist_and_dEdx_ok = false;
  //   }
  //   else dist_and_dEdx_ok = false;
  // }
  // else {
  //   // no track, 1D dE/dx cut
  //   if ( shr_tkfit_gap10_dedx_Y < 1.7 ) dist_and_dEdx_ok = false;
  //   if ( shr_tkfit_gap10_dedx_Y > 2.7 && shr_tkfit_gap10_dedx_Y < 5.5 ) {
  //     dist_and_dEdx_ok = false;
  //   }
  // }

  // bool sel_pass_electron_identification = moliere_ok && dist_and_dEdx_ok;

  // Flag indicating whether the event passed the full selection
  // bool sel_nu_e_cc = sel_pass_preselection && sel_pass_cosmic_rejection
  //   && sel_pass_shower_identification && sel_pass_electron_identification;

  bool sel_nu_e_cc = sel_pass_preselection && sel_pass_cosmic_rejection
    && sel_pass_other;

  // // Set the reco energy of the electron candidate if we found one
  // double reco_electron_energy = BOGUS;
  // if ( sel_pass_shower_identification ) {
  //   // Apply the shower energy correction factor
  //   reco_electron_energy = shr_energy_cali / 0.83;
  // }

  // Set the reco energy of the electron candidate if we found one
  double reco_electron_energy = BOGUS;
  if ( sel_nu_e_cc ) {
    // Apply the shower energy correction factor
    reco_electron_energy = shr_energy_cali / 0.83;
  }

  // We're done. Store the results in the output tree and return whether the
  // event passed the full selection
  auto& out = ev.out();
  out[ "sel_pass_preselection" ] = sel_pass_preselection;
  out[ "sel_pass_cosmic_rejection" ] = sel_pass_cosmic_rejection;
  out[ "sel_pass_other" ] = sel_pass_other;
  out[ "sel_nu_e_cc" ] = sel_nu_e_cc;
  out[ "reco_electron_energy" ] = reco_electron_energy;

  return sel_nu_e_cc;
}
