// This script adds beamline geometry weights to NuMI  MC files
// It should be run before ProcessNTuples.C
// The beamline geometry weights are stored in a separate ROOT file, required as input
// Note - path needs to be set

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TRotation.h"

float GetNuMIAngle(double px, double py, double pz, std::string direction); 

std::vector<float> getWeights(float nu_e, float nu_angle, std::vector<TH2F> &h_weights);

float checkWeight(float weight);

// Recursively copy all contents of a TDirectory to another TDirectory.
// Uses CloneTree for TTrees to correctly copy all basket data across files.
void copyDirectoryRecursive( TDirectory* src, TDirectory* dst ) {
    TIter nextkey( src->GetListOfKeys() );
    TKey* key;
    while ( (key = (TKey*)nextkey()) ) {
        src->cd();
        TObject* obj = key->ReadObj();
        TDirectory* subdir = dynamic_cast<TDirectory*>( obj );
        TTree*      tree   = dynamic_cast<TTree*>( obj );
        if ( subdir ) {
            TDirectory* newdir = dst->mkdir( subdir->GetName() );
            newdir->cd();
            copyDirectoryRecursive( subdir, newdir );
            dst->cd();
        } else if ( tree ) {
            dst->cd();
            TTree* clone = tree->CloneTree( -1, "fast" );
            clone->Write();
        } else {
            dst->cd();
            obj->Write();
        }
        delete obj;
    }
}

int main(int argc, char *argv[]) {

    if ( argc != 4 ) {
        std::cout << "Usage: AddBeamlineGeometryWeights INPUT_FILE HORN_CURRENT_MODE OUTPUT_FILE" << std::endl;
        return 1;
    }

    std::string input_filename( argv[1] );
    std::string horn_current_mode( argv[2] );
    std::string output_filename( argv[3] );

    std::cout << "Adding beamline geometry weights to file: " << input_filename << std::endl;
    if (horn_current_mode == "FHC" || horn_current_mode == "RHC") {
        std::cout << "Horn current mode: " << horn_current_mode << std::endl;
    } else {
        std::cout << "Error: invalid horn current mode. Valid modes: FHC, RHC" << std::endl;
        exit(1);
    }

    TFile *beamlineVariationsFile = new TFile("src/app/NuMI/NuMI_Geometry_Weights_Histograms.root");
    if (!beamlineVariationsFile || beamlineVariationsFile->IsZombie()) {
        std::cerr << "Error: Could not open file 'NuMI_Geometry_Weights_Histograms.root'." << std::endl;
        exit(1);
    }

    std::vector<TH2F> h_nue, h_nuebar, h_numu, h_numubar;

    for (int i = 1; i <= 20; i++) {
        std::stringstream name_nue_ss; name_nue_ss << "EnergyTheta2D/ratio_run" << i << "_" << horn_current_mode << "_nue_CV_AV_TPC_2D";
        std::stringstream name_nuebar_ss; name_nuebar_ss << "EnergyTheta2D/ratio_run" << i << "_" << horn_current_mode << "_nuebar_CV_AV_TPC_2D";
        std::stringstream name_numu_ss; name_numu_ss << "EnergyTheta2D/ratio_run" << i << "_" << horn_current_mode << "_numu_CV_AV_TPC_2D";
        std::stringstream name_numubar_ss; name_numubar_ss << "EnergyTheta2D/ratio_run" << i << "_" << horn_current_mode << "_numubar_CV_AV_TPC_2D";
        h_nue.push_back(*(TH2F*)beamlineVariationsFile->Get(name_nue_ss.str().c_str()));
        h_nuebar.push_back(*(TH2F*)beamlineVariationsFile->Get(name_nuebar_ss.str().c_str()));
        h_numu.push_back(*(TH2F*)beamlineVariationsFile->Get(name_numu_ss.str().c_str()));
        h_numubar.push_back(*(TH2F*)beamlineVariationsFile->Get(name_numubar_ss.str().c_str()));
    }

    beamlineVariationsFile->Close();
    delete beamlineVariationsFile;

    TFile *f = new TFile(input_filename.c_str());
    if (!f->IsOpen()) { std::cout << "Could not open input file." << std::endl; exit(1); }

    TFile *f_out = new TFile(output_filename.c_str(), "recreate");
    if (!f_out->IsOpen()) { std::cout << "Could not open output file." << std::endl; exit(1); }

    // Copy all objects except nuselection using recursive copy to preserve
    // subdirectory contents (e.g. wcpselection and its TTrees)
    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        f->cd();
        TObject    *obj  = key->ReadObj();
        TDirectory *dir  = dynamic_cast<TDirectory*>(obj);
        TTree      *tree = dynamic_cast<TTree*>(obj);
        if (dir && std::string(dir->GetName()) == "nuselection") {
            continue; // handled separately below
        } else if (dir) {
            TDirectory *newdir = f_out->mkdir(dir->GetName());
            newdir->cd();
            copyDirectoryRecursive(dir, newdir);
            f_out->cd();
        } else if (tree) {
            f_out->cd();
            TTree *clone = tree->CloneTree(-1, "fast");
            clone->Write();
        } else {
            f_out->cd();
            obj->Write();
        }
    }

    // Handle nuselection directory
    TDirectory *indir = (TDirectory*)f->Get("nuselection");
    if (indir) {
        f_out->mkdir("nuselection");
        f_out->cd("nuselection");
        TIter nextkey2(indir->GetListOfKeys());
        TKey *key2;
        while ((key2 = (TKey*)nextkey2())) {
            std::string objname = key2->GetName();
            if (objname == "NeutrinoSelectionFilter" || objname == "SubRun") continue;
            TObject *obj2 = key2->ReadObj();
            obj2->Write();
        }
    }

    // Process nuselection/NeutrinoSelectionFilter and SubRun
    TTree *t_nu  = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
    TTree *t_pot = (TTree*)f->Get("nuselection/SubRun");

    int nu_pdg;
    float nu_e, true_nu_px, true_nu_py, true_nu_pz;
    t_nu->SetBranchAddress("nu_pdg",      &nu_pdg);
    t_nu->SetBranchAddress("nu_e",        &nu_e);
    t_nu->SetBranchAddress("true_nu_px",  &true_nu_px);
    t_nu->SetBranchAddress("true_nu_py",  &true_nu_py);
    t_nu->SetBranchAddress("true_nu_pz",  &true_nu_pz);
    std::map<std::string, std::vector<double>> *mc_weights_map = nullptr;
    t_nu->SetBranchAddress("weights", &mc_weights_map);
    t_nu->GetEntry(0);

    auto Horn_2kA_index = mc_weights_map->find("Horn_2kA");
    if (Horn_2kA_index != mc_weights_map->end()) {
        std::cout << "Beamline geometry weights already present." << std::endl;
        f->Close(); f_out->Close();
        delete f; delete f_out;
        exit(0);
    }

    f_out->cd("nuselection");
    TTree *t_out_nu  = t_nu->CloneTree(0);
    TTree *t_out_pot = t_pot->CloneTree();
    int nEntries = t_nu->GetEntries();
    std::cout << "Number entries: " << nEntries << std::endl;

    for (int iEntry = 0; iEntry < nEntries; iEntry++) {
        t_nu->GetEntry(iEntry);
        if ((iEntry != 0) && (nEntries >= 10) && (iEntry % (nEntries/10) == 0))
            std::cout << Form("%i0%% Completed...\n", iEntry / (nEntries/10));

        float nu_angle = GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "beam");
        std::vector<float> weights;
        if      (nu_pdg ==  12) weights = getWeights(nu_e, nu_angle, h_nue);
        else if (nu_pdg == -12) weights = getWeights(nu_e, nu_angle, h_nuebar);
        else if (nu_pdg ==  14) weights = getWeights(nu_e, nu_angle, h_numu);
        else if (nu_pdg == -14) weights = getWeights(nu_e, nu_angle, h_numubar);
        else { std::cout << "Error: cannot get beamline variation weights" << std::endl; exit(1); }

        if (weights.size() != 20) { std::cout << "Error: missing expected beamline variation weights" << std::endl; exit(1); }

        mc_weights_map->insert({"Horn_2kA",           {checkWeight(weights[0]),  checkWeight(weights[1])}});
        mc_weights_map->insert({"Horn1_x_3mm",         {checkWeight(weights[2]),  checkWeight(weights[3])}});
        mc_weights_map->insert({"Horn1_y_3mm",         {checkWeight(weights[4]),  checkWeight(weights[5])}});
        mc_weights_map->insert({"Beam_spot_1_1mm",     {checkWeight(weights[6])}});
        mc_weights_map->insert({"Beam_spot_1_5mm",     {checkWeight(weights[7])}});
        mc_weights_map->insert({"Horn2_x_3mm",         {checkWeight(weights[8]),  checkWeight(weights[9])}});
        mc_weights_map->insert({"Horn2_y_3mm",         {checkWeight(weights[10]), checkWeight(weights[11])}});
        mc_weights_map->insert({"Horns_0mm_water",     {checkWeight(weights[12])}});
        mc_weights_map->insert({"Horns_2mm_water",     {checkWeight(weights[13])}});
        mc_weights_map->insert({"Beam_shift_x_1mm",   {checkWeight(weights[14]), checkWeight(weights[15])}});
        mc_weights_map->insert({"Beam_shift_y_1mm",   {checkWeight(weights[16]), checkWeight(weights[17])}});
        mc_weights_map->insert({"Target_z_7mm",       {checkWeight(weights[18]), checkWeight(weights[19])}});
        t_out_nu->Fill();
    }

    t_out_nu->Write("NeutrinoSelectionFilter");
    t_out_pot->Write("SubRun");
    f->Close(); f_out->Close();
    delete f; delete f_out;
}


float GetNuMIAngle(double px, double py, double pz, std::string direction) {
    TRotation RotDet2Beam;
    TVector3  detxyz, BeamCoords;
    std::vector<double> rotmatrix;
    detxyz = {px, py, pz};
    rotmatrix = {
        0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
        4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
        -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };
    TVector3 newX(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    TVector3 newY(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    TVector3 newZ(rotmatrix[6], rotmatrix[7], rotmatrix[8]);
    RotDet2Beam.RotateAxes(newX, newY, newZ);
    BeamCoords = RotDet2Beam * detxyz;
    TVector3 beamdir = {0, 0, 1};
    if (direction == "target") {
        beamdir = {5502, 7259, 67270};
        beamdir = beamdir.Unit();
    } else if (direction != "beam") {
        std::cout << "Warning unknown angle type specified, you should check this" << std::endl;
    }
    return BeamCoords.Angle(beamdir) * 180 / 3.1415926;
}

std::vector<float> getWeights(float nu_e, float nu_angle, std::vector<TH2F> &h_weights) {
    std::vector<float> weights; weights.reserve(20);
    for (int i = 0; i < (int)h_weights.size(); i++) {
        int binx = h_weights[i].GetXaxis()->FindBin(nu_e);
        int biny = h_weights[i].GetYaxis()->FindBin(nu_angle);
        weights.push_back(h_weights[i].GetBinContent(binx, biny));
    }
    return weights;
}

float checkWeight(float weight) {
    if      (std::isinf(weight))   weight = 1.0;
    else if (std::isnan(weight))   weight = 1.0;
    else if (weight > 30.0)        weight = 1.0;
    else if (weight < 0.0)         weight = 1.0;
    else if (weight < 1e-4)        weight = 0.0;
    return weight;
}