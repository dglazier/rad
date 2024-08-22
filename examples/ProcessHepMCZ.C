//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "HepMCElectro.h"
#include "ParticleCreator.h"
#include "Indicing.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

//beam components must be defined even if 0
inline constexpr std::array<double,4>  rad::beams::InitBotComponents() {return {0.,0.,0,0};}
inline constexpr std::array<double,4>  rad::beams::InitTopComponents() {return {0.,0.,0.,0};}

void ProcessHepMCZ(){
  gBenchmark->Start("df");

  //create reaction dataframe
  rad::config::HepMCElectro hepmc{"hepmc3_tree", "/home/dglazier/Dropbox/EIC/EventGenerators/elSpectro/examples/out/jpac_z3900_10x100.hepmc.root"};
  hepmc.AliasMomentumComponents();
   
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  hepmc.setBeamIonIndex(1);
  hepmc.setBeamElectronIndex(0);
  hepmc.setScatElectronIndex(6);
  //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  hepmc.setParticleIndex("idxEl",2);
  hepmc.setParticleIndex("idxPo",3);
  hepmc.setParticleIndex("idxPi",4);
  hepmc.setParticleIndex("idxN",5);

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  hepmc.setBaryonParticles({"idxN"});
  
  rad::config::ParticleCreator particles{hepmc};
  particles.Sum("idxJ",{"idxEl","idxPo"});
  particles.Sum("idxZ",{"idxPi","idxJ"});

  hepmc.setMesonParticles({"idxJ","idxPi"});

  //must call this after all particles are configured
  hepmc.makeParticleMap();
  
  //can also add missing particles
  //but must be done after ParticleMap
  //so currently cannot use as baryon...
  particles.Miss("idxCalcN",{"scat_ele","idxZ"});

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(hepmc,"W","{scat_ele}");
  rad::rdf::Mass(hepmc,"Whad","{idxEl,idxPo,idxPi,idxN}");
  rad::rdf::Mass(hepmc,"JMass","{idxEl,idxPo}");
  rad::rdf::Mass(hepmc,"ZMass","{idxEl,idxPo,idxPi}");
  rad::rdf::Mass(hepmc,"MissMassZ","{idxCalcN}");

  //t distribution, column name
  rad::rdf::TBot(hepmc,"tb_pn");
  rad::rdf::TPrimeBot(hepmc,"tbp_pn");
  rad::rdf::TTop(hepmc,"tt_gZ");
  rad::rdf::TPrimeTop(hepmc,"ttp_gZ");

  //CM production angles
  rad::rdf::CMAngles(hepmc,"CM");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  auto df0 = hepmc.CurrFrame();

  auto hW = df0.Histo1D({"W","W",100,0,20.},"mc_W");
  auto hWhad = df0.Histo1D({"Whad","Whad",100,0,20.},"mc_Whad");
  auto hMesonMass = df0.Histo1D({"MesonMass","M(e-,e+, #pi) [GeV]",100,.3,5.},"mc_ZMass");
  auto hJMass = df0.Histo1D({"JMass","M(e-,e+) [GeV]",100,.3,5.},"mc_JMass");
  auto htpn = df0.Histo1D({"tpn","t(p,n) [GeV^{2}]",100,0,5},"mc_tb_pn");
  auto htgZ = df0.Histo1D({"tgZ","t(g,Z) [GeV^{2}]",100,0,5},"mc_tt_gZ");
  auto htprimepn = df0.Histo1D({"tprimepn","t'(p,n) [GeV^{2}]",100,0,5},"mc_tbp_pn");
  auto hthCM=df0.Histo1D({"cthCM","cos(#theta_{CM})",100,-1,1},"mc_CM_CosTheta");
  auto hphCM=df0.Histo1D({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},"mc_CM_Phi");
 
  gBenchmark->Start("processing");
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
  hWhad->DrawCopy();
  new TCanvas();
  hMesonMass->DrawCopy();
  new TCanvas();
  hJMass->DrawCopy();
  new TCanvas();
  htpn->DrawCopy();
  htgZ->DrawCopy("same");
  //  htprimepn->DrawCopy("same");
  auto canCM = new TCanvas();
  canCM->Divide(2,1);
  canCM->cd(1);
  hthCM->DrawCopy();
  canCM->cd(2);
  hphCM->DrawCopy();

   
  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");
  gBenchmark->Stop("df");
  gBenchmark->Print("df");
  
  hepmc.Snapshot("HepMCZ.root");

}
