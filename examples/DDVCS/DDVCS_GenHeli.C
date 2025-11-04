//This macro generates the "flat" helicity distributions
//for a given input file, and outputs a new hepmc rootfile
//the "flat" hepmc rootfile can then be fully analysed with
//the regular rad macros and workflow

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

void DDVCS_GenHeli(const TString inputfile){
  
  TString outputfile=inputfile;
  outputfile.ReplaceAll(".root","_flat.root");
  
  // Enable implicit multi-threading
  ROOT::EnableImplicitMT(32);
  
  using namespace rad::names::data_type; //for MC()
  
  gBenchmark->Start("df");
  
  //create reaction dataframe
  rad::config::HepMCElectro hepmc{"hepmc3_tree", inputfile };
  
   hepmc.AliasMomentumComponents();
    
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  hepmc.setBeamIonIndex(3);
  hepmc.setBeamElectronIndex(0);
  hepmc.setScatElectronIndex(1);
  
  hepmc.setParticleIndex("pprime",5);
  hepmc.setBaryonParticles({"pprime"});
  
  //if using existing lepton pair
  //create gprime first from their sums
  rad::config::ParticleCreator particles{hepmc};
  hepmc.setParticleIndex("ele0",6);
  hepmc.setParticleIndex("pos0",7);
  particles.Sum("gprime",{"ele0","pos0"});
  
  //use ParticleGenerator to decay gprime to 2 M_ele particles ele,pos
  // hepmc.setParticleIndex("gprime",4);//
  rad::generator::ParticleGenerator gen{hepmc};
  gen.GenerateTwoBody({"ele","pos"},
		      {rad::constant::M_ele(),rad::constant::M_ele()},
		      "gprime");
  
  rad::config::ParticleModifier mod{hepmc};
  mod.ModifyMomentum("ele0","ele");
  mod.ModifyMomentum("pos0","pos");
  mod.Apply("NewMesons");
  hepmc.setMesonParticles({"ele0","pos0"});
  
  //must call this after all particles are configured
  hepmc.makeParticleMap();
  //rad::rdf::PrintParticles(hepmc,MC());
  
  gBenchmark->Start("snapshot");
  //hepmc.BookHepMC3Snapshot(outputfile.Data());
  hepmc.HepMC3Snapshot(outputfile.Data());
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");
  
  gBenchmark->Stop("df");
  gBenchmark->Print("df");
  
}
