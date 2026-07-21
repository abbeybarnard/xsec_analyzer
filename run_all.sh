# Set up the code as usual 
source setup_xsec_analyzer.sh
make clean
make 

# Activate the BDT Python environment
source /exp/uboone/data/users/abarnard/bdt-env/bin/activate

# Compute the BDT scores for the files_to_process.txt files 
python scripts/compute_bdt_scores.py

# Add the fake weights for the files_to_process_bdt.txt files
# Ignores data and detector variation files 
./scripts/AddWeights.sh

# Run the selection on the files_to_process_bdt_weighted.txt files
# Produces files with the `xsec-ana-` prefix, which needs updating in file_properties.txt
./scripts/ReprocessNTuples.sh /exp/uboone/data/users/abarnard/analysis/selection_output NuMICC1eNp configs/files_to_process_bdt_weighted.txt

# Run the systematics using the binning defined in RunUnivmake.sh
./scripts/RunUnivmake.sh v3_detvars 