// Univmake expects all MC files to have identical systematic weights structure.
// For the original processing of NuMI Dirt MC the weights were not included.
// This script adds a fake set weights with identical structure to the NuMI 
// Dirt MC files to make them compatible.
// This script should be run after addBeamlineGeometryWeights.cpp and before ProcessNTuples.C.

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"

// Recursively copy all contents of a TDirectory to another TDirectory.
// Uses CloneTree for TTrees to correctly copy all basket data across files.
void copyDirectoryRecursive( TDirectory* src, TDirectory* dst ) {
    TIter nextkey( src->GetListOfKeys() );
    TKey* key;
    while ( (key = (TKey*)nextkey()) ) {
        src->cd();
        TObject*    obj    = key->ReadObj();
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
        std::cout << "Usage: AddFakeWeights INPUT_FILE INPUT_FULL_WEIGHTS_FILE OUTPUT_FILE" << std::endl;
        return 1;
    }

    std::string input_filename( argv[1] );
    std::string input_full_weights_filename( argv[2] );
    std::string output_filename( argv[3] );

    std::cout << "Adding fake weights to file: " << input_filename
              << ", using structure from: " << input_full_weights_filename << std::endl;

    // Load file containing real weights structure
    TFile *f1 = new TFile(input_full_weights_filename.c_str());
    if (!f1->IsOpen()) { std::cout << "Could not open input file containing weights structure." << std::endl; exit(1); }

    TTree *t1 = (TTree*)f1->Get("nuselection/NeutrinoSelectionFilter");
    std::map<std::string, std::vector<double>> *weights = nullptr;
    t1->SetBranchAddress("weights", &weights);
    t1->GetEntry(0);

    std::cout << "Number of weight pairs in map: " << weights->size() << std::endl;

    // Create fake weights (all set to 1)
    std::map<std::string, std::vector<double>> fake_weights;
    std::cout << std::endl;
    for (auto pair : *weights) {
        std::cout << "Name: " << pair.first << ", Size: " << (pair.second).size() << std::endl;
        std::fill((pair.second).begin(), (pair.second).end(), 1);
        fake_weights.insert({pair.first, pair.second});
    }
    std::cout << std::endl;

    f1->Close();
    delete f1;

    // Open input file to add weights to
    TFile *f2 = new TFile(input_filename.c_str());
    if (!f2->IsOpen()) { std::cout << "Could not open input file." << std::endl; exit(1); }

    // Open output file
    TFile *f3 = new TFile(output_filename.c_str(), "recreate");
    if (!f3->IsOpen()) { std::cout << "Could not open output file." << std::endl; exit(1); }

    // Copy all objects except nuselection using recursive copy to preserve
    // subdirectory contents (e.g. wcpselection and its TTrees)
    TIter nextkey(f2->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        f2->cd();
        TObject    *obj  = key->ReadObj();
        TDirectory *dir  = dynamic_cast<TDirectory*>(obj);
        TTree      *tree = dynamic_cast<TTree*>(obj);
        if (dir && std::string(dir->GetName()) == "nuselection") {
            continue; // handled separately below
        } else if (dir) {
            TDirectory *newdir = f3->mkdir(dir->GetName());
            newdir->cd();
            copyDirectoryRecursive(dir, newdir);
            f3->cd();
        } else if (tree) {
            f3->cd();
            TTree *clone = tree->CloneTree(-1, "fast");
            clone->Write();
        } else {
            f3->cd();
            obj->Write();
        }
    }

    // Handle nuselection directory
    TDirectory *indir = (TDirectory*)f2->Get("nuselection");
    if (indir) {
        f3->mkdir("nuselection");
        f3->cd("nuselection");
        TIter nextkey2(indir->GetListOfKeys());
        TKey *key2;
        while ((key2 = (TKey*)nextkey2())) {
            std::string objname = key2->GetName();
            if (objname == "NeutrinoSelectionFilter" || objname == "SubRun") continue;
            TObject *obj2 = key2->ReadObj();
            obj2->Write();
        }
    }

    // Process NeutrinoSelectionFilter — replace weights branch with fake weights
    TTree *t2_nu  = (TTree*)f2->Get("nuselection/NeutrinoSelectionFilter");
    TTree *t2_pot = (TTree*)f2->Get("nuselection/SubRun");
    f3->cd("nuselection");
    t2_nu->SetBranchStatus("weights", 0);
    TTree *t3_nu  = t2_nu->CloneTree(0);
    TTree *t3_pot = t2_pot->CloneTree();
    std::map<std::string, std::vector<double>> new_weights;
    t3_nu->Branch("weights", "std::map<std::string, std::vector<double>>", &new_weights);

    int n_entries = t2_nu->GetEntries();
    for (int e = 0; e < n_entries; e++) {
        t2_nu->GetEntry(e);
        if ((e != 0) && (n_entries >= 10) && (e % (n_entries/10) == 0))
            std::cout << Form("%i0%% Completed...\n", e / (n_entries/10));
        new_weights = fake_weights;
        t3_nu->Fill();
    }

    t3_nu->Write("NeutrinoSelectionFilter");
    t3_pot->Write("SubRun");
    f2->Close(); f3->Close();
    delete f2; delete f3;
}