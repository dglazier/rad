//THIS SCRIPT INTRODUCES QUITE A BIT OF TECH DEBT
//SHOULD FIND A WAY IN THE FUTURE TO AUTOMATE SOME/MOST
//OF THE THINGS ITS DOING AT THE RAD LEVEL?
//FOR NOW GOOD AND FAST WAY TO GET RAD OUTPUTS
//MIXED AND READY FOR BRUFIT

#include "HepMCElectro.h"
#include "ParticleCreator.h"
#include "ParticleGenerator.h"
#include "ParticleModifier.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

#include "ProcessHepMCDDVCS.C"

std::vector<Long64_t> unique_random_indices(Long64_t total, Long64_t pick, UInt_t seed) {
    std::vector<Long64_t> indices(total);
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 rng(seed);
    std::shuffle(indices.begin(), indices.end(), rng);

    indices.resize(pick);
    return indices;
}

bool checkFileExists(const std::string &filename) {
  if (gSystem->AccessPathName(filename.c_str())) {
        std::cout << "File does NOT exist: " << filename << std::endl;
	return false;
    } else {
        std::cout << "File exists: " << filename << std::endl;
	return true;
    }
}

void ShuffleTree(const std::string &infile, const std::string &outfile, const std::string &treename) {
    // Open input file & tree
    TFile inFile(infile.c_str(), "READ");
    if (inFile.IsZombie()) {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }
    TTree* tIn = nullptr;
    inFile.GetObject(treename.c_str(), tIn);
    if (!tIn) {
        std::cerr << "Tree " << treename << " not found in " << infile << std::endl;
        return;
    }

    // Create output file & clone tree structure
    TFile outFile(outfile.c_str(), "RECREATE");
    TTree* tOut = tIn->CloneTree(0); // clone structure only

    // Prepare index vector
    Long64_t nEntries = tIn->GetEntries();
    std::vector<Long64_t> indices(nEntries);
    std::iota(indices.begin(), indices.end(), 0);
    
    // ROOT RNG ? std::mt19937
    TRandom3 rootRand(0); // seed
    std::mt19937 gen(static_cast<unsigned>(rootRand.Integer(1e9)));
    
    // Shuffle indices
    std::shuffle(indices.begin(), indices.end(), gen);

    // Loop over shuffled indices
    for (auto idx : indices) {
        tIn->GetEntry(idx);
        tOut->Fill();
    }

    // Write output
    tOut->Write();
    outFile.Close();
    inFile.Close();

    std::cout << "Shuffled tree saved to " << outfile << std::endl;
}

void DFMerge(const std::string plus_file, const std::string minus_file, const std::string outfile, const std::string treename, const double pol, const int nev){
  
  gBenchmark->Start("Setup");
  
  ROOT::RDataFrame df_plus(treename, plus_file);
  ROOT::RDataFrame df_minus(treename, minus_file);
  
  // Step 1: Add helicity and thread-safe RNG
  auto df_plus_heli = df_plus.Define("egen_helicity", "1")
    .Define("rand", [] { return rad::random::Generator().Rndm(); }).Define("UID","rdfentry_");
  
  auto df_minus_heli = df_minus.Define("egen_helicity", "-1")
    .Define("rand", [] { return rad::random::Generator().Rndm(); }).Define("UID","rdfentry_");
  
  
  // Step 2: calculate asym n+ and n-
  auto nplus = (double)nev * (1.0 + pol) / 2.0; //9000 for 0.8 e.g.
  auto nminus = (double)nev * (1.0 - pol) / 2.0 ; //1000 for 0.8 e.g.
  
  auto total_plus = *df_plus_heli.Count();
  auto total_minus = *df_minus_heli.Count();
  
  // Random non-overlapping indices
  auto all_plus_indices = unique_random_indices(total_plus, nplus + nminus, 42);
  auto all_minus_indices = unique_random_indices(total_minus, nplus + nminus, 43);
  
  std::vector<Long64_t> indices_plus_for_plus(all_plus_indices.begin(), all_plus_indices.begin() + nplus);
  std::vector<Long64_t> indices_plus_for_minus(all_plus_indices.begin() + nplus, all_plus_indices.end());
  
  std::vector<Long64_t> indices_minus_for_minus(all_minus_indices.begin(), all_minus_indices.begin() + nplus);
  std::vector<Long64_t> indices_minus_for_plus(all_minus_indices.begin() + nplus, all_minus_indices.end());
  
  // Convert to unordered_set for O(1) filtering
  auto set_plus_for_plus = std::unordered_set<Long64_t>(indices_plus_for_plus.begin(), indices_plus_for_plus.end());
  auto set_plus_for_minus = std::unordered_set<Long64_t>(indices_plus_for_minus.begin(), indices_plus_for_minus.end());
  auto set_minus_for_minus = std::unordered_set<Long64_t>(indices_minus_for_minus.begin(), indices_minus_for_minus.end());
  auto set_minus_for_plus = std::unordered_set<Long64_t>(indices_minus_for_plus.begin(), indices_minus_for_plus.end());
   
  gBenchmark->Stop("Setup");
  gBenchmark->Print("Setup");
  
  gBenchmark->Start("Filter");
  
  // Filter each group
    auto df_plus_plus = df_plus_heli.Filter(
        [set_plus_for_plus](ULong64_t entry) {
            return set_plus_for_plus.count(entry);
        }, {"rdfentry_"});

    auto df_plus_minus = df_minus_heli.Filter(
        [set_minus_for_plus](ULong64_t entry) {
            return set_minus_for_plus.count(entry);
        }, {"rdfentry_"});

    auto df_minus_plus = df_plus_heli.Filter(
        [set_plus_for_minus](ULong64_t entry) {
            return set_plus_for_minus.count(entry);
        }, {"rdfentry_"});

    auto df_minus_minus = df_minus_heli.Filter(
        [set_minus_for_minus](ULong64_t entry) {
            return set_minus_for_minus.count(entry);
        }, {"rdfentry_"});

  
  gBenchmark->Stop("Filter");
  gBenchmark->Print("Filter");
  
  vector<std::string> brufit_cols = {"UID","mc_t_bot","mc_GammaPolCirc","mc_GammaE","mc_Heli_Phi","mc_Heli_Theta","mc_Heli_CosTheta","mc_CM_Phi","mc_CM_Theta"};
  
  gBenchmark->Start("Snapshots");
  auto df_plus_plus_snap = df_plus_plus.Define("pol",[pol] { return pol; }).Define("pol_helicity","1").Redefine("mc_GammaPolCirc","mc_GammaPolCirc*pol").Snapshot(treename,"plus_plus_snap_test.root",brufit_cols);
  auto df_plus_minus_snap = df_plus_minus.Define("pol",[pol] { return pol; }).Define("pol_helicity","1").Redefine("mc_GammaPolCirc","mc_GammaPolCirc*pol").Snapshot(treename,"plus_minus_snap_test.root",brufit_cols);
  auto df_minus_plus_snap = df_minus_plus.Define("pol",[pol] { return -pol; }).Define("pol_helicity","-1").Redefine("mc_GammaPolCirc","mc_GammaPolCirc*pol").Snapshot(treename,"minus_plus_snap_test.root",brufit_cols);
  auto df_minus_minus_snap = df_minus_minus.Define("pol",[pol] { return -pol; }).Define("pol_helicity","-1").Redefine("mc_GammaPolCirc","mc_GammaPolCirc*pol").Snapshot(treename,"minus_minus_snap_test.root",brufit_cols);
  gBenchmark->Stop("Snapshots");
  gBenchmark->Print("Snapshots");
  
  gBenchmark->Start("hadd");
  gSystem->Exec("hadd -f pol_mixed.root plus_plus_snap_test.root plus_minus_snap_test.root minus_plus_snap_test.root minus_minus_snap_test.root");
  
  gSystem->Exec("rm *snap_test.root");
  gBenchmark->Stop("hadd");
  gBenchmark->Print("hadd");

  gBenchmark->Start("Shuffle");
  ShuffleTree("pol_mixed.root", outfile, treename);
  gSystem->Exec("rm pol_mixed.root");
  gBenchmark->Stop("Shuffle");
  gBenchmark->Print("Shuffle");
}

void MixHelicityTrees(){
  
  // Enable implicit multi-threading
  ROOT::EnableImplicitMT(32);
  gROOT->SetBatch(kTRUE);
  gBenchmark->Start("Total");
  
  std::string plus_file = "/w/work5/home/garyp/eic/Farm/data/EpIC_ep_DDVCS_18x275/18x275_ddvcs_1M_events_plus.root";
  std::string minus_file = "/w/work5/home/garyp/eic/Farm/data/EpIC_ep_DDVCS_18x275/18x275_ddvcs_1M_events_minus.root";
  
  std::string plus_flat_file = "/w/work5/home/garyp/eic/Farm/data/EpIC_ep_DDVCS_18x275/18x275_ddvcs_1M_events_plus_flat.root";
  std::string minus_flat_file = "/w/work5/home/garyp/eic/Farm/data/EpIC_ep_DDVCS_18x275/18x275_ddvcs_1M_events_minus_flat.root";
  
  //this uses DDVCS_GenHeli.C to generate flat phase space files from original files
  //commented because i dont want to redo this every time right now
  //because i have generated the flat phase space then simulated those events 
  //DDVCS_GenHeli(plus_file,plus_flat_file);
  //DDVS_GenHeli(minus_file,minus_flat_file);
  
  std::string plus_outfile = "HepMC_ddvcs_plus.root";
  std::string minus_outfile = "HepMC_ddvcs_minus.root";
  std::string plus_flat_outfile = "HepMC_ddvcs_plus_flat.root";
  std::string minus_flat_outfile = "HepMC_ddvcs_minus_flat.root";
  
  if(!checkFileExists(plus_outfile))
    ProcessHepMCDDVCS(plus_file, plus_outfile);
  
  if(!checkFileExists(minus_outfile))
    ProcessHepMCDDVCS(minus_file, minus_outfile);
  
  if(!checkFileExists(plus_flat_outfile))
    ProcessHepMCDDVCS(plus_flat_file, plus_flat_outfile);
  
  if(!checkFileExists(minus_flat_outfile))
    ProcessHepMCDDVCS(minus_flat_file, minus_flat_outfile);
  
  
  std::string mixed_outfile = "HepMC_ddvcs_mixed.root";
  std::string mixed_flat_outfile = "HepMC_ddvcs_mixed_flat.root";
  
  std::string treename="rad_tree";
  
  double pol=0.8;
  int nev=100000;
  
  DFMerge(plus_outfile, minus_outfile, mixed_outfile, treename, pol, nev);
  DFMerge(plus_flat_outfile, minus_flat_outfile, mixed_flat_outfile, treename, pol, nev);

  gBenchmark->Stop("Total");
  gBenchmark->Print("Total");
  
  
  //here we now have flat and physics helicity angle files
  //DONE 1) convert hepmc files to rad_trees adding a helcity and including all kinemtics needed by the brufit code
  //DONE 2) use brufit toy generator to fold in physics to a subset of these events
  
  // 3) Fit the generated toy with the full data-set used as SimulatedDAta in brufit to calculate the normalsiation integral
  // 4) Now do the same thing starting from the simulated and reconstructed data
  // 5) Now fit the full EPiC generated hepmc files without generating toys and again using the flat data you created as the SimulatedData in brufit
  // 6) Finally fir the simulated and reconstructed EPIc Physics events
  
} 
