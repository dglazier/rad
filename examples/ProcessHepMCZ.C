//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "include/HepMCElectro.h"
#include "include/Indicing.h"
#include "include/BasicKinematicsRDF.h"
#include "include/ReactionKinematicsRDF.h"
#include "include/ElectronScatterKinematicsRDF.h"
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
  hepmc.setParticleIndex("el",2);
  hepmc.setParticleIndex("po",3);
  hepmc.setParticleIndex("pi",4);
  hepmc.setParticleIndex("n",5);

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  hepmc.setBaryonParticles({"n"});
  
  hepmc.setMesonParticles({"el","po","pi"});

  //must call this after all particles are configured
  hepmc.makeParticleMap();
  
  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(hepmc,"W","{scat_ele}");
  rad::rdf::Mass(hepmc,"Whad","{el,po,pi,n}");
  rad::rdf::Mass(hepmc,"JMass","{el,po}");
  rad::rdf::Mass(hepmc,"ZMass","{el,po,pi}");

  //t distribution, column name
  rad::rdf::TBot(hepmc,"tb_pn");
  rad::rdf::TPrimeBot(hepmc,"tbp_pn");
  rad::rdf::TTop(hepmc,"tt_pn");
  rad::rdf::TPrimeTop(hepmc,"ttp_pn");

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
  htprimepn->DrawCopy("same");
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
