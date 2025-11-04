//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "HepMCElectro.h"
//#include "ParticleCreator.h"
#include "Indicing.h"
//#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
//#include "ReactionKinematicsRDF.h"
//#include "ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>


void ProcessHepMCZCombi(){
  using namespace rad::names::data_type; //for MC()
   gBenchmark->Start("df");

  //create reaction dataframe
  rad::config::HepMCElectro hepmc{"hepmc3_tree", "/home/dglazier/Dropbox/EIC/EventGenerators/elSpectro/examples/out/jpac_z3900_10x100.hepmc.root"};
  hepmc.AliasMomentumComponents();
    
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  hepmc.setBeamIonIndex(1);
  hepmc.setBeamElectronIndex(0);
  hepmc.setScatElectronCandidates({6,2});
  //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  hepmc.setParticleCandidates("ele",{2,6});
  hepmc.setParticleIndex("pos",3);
  hepmc.setParticleIndex("pip",4);
  hepmc.setParticleIndex("n",5);

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calculating reaction kinematics
  //e.g. t distributions
  //hepmc.setBaryonParticles({"n"});

  //hepmc.setMesonParticles({"ele","pos","pip"});

  //must call this after all particles are configured
  //hepmc.makeParticleMap();
  hepmc.makeCombinations();
  
  rad::rdf::Mass(hepmc,"Whad","{ele,pos,pip,n}");

 
  gBenchmark->Start("snapshot");
  hepmc.Snapshot("HepMCZCombi.root");
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");

}
