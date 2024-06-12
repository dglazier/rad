//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "include/ePICReaction.h"
#include "include/Indicing.h"
#include "include/BasicKinematicsRDF.h"
#include "include/ReactionKinematicsRDF.h"
#include "include/ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

inline constexpr std::array<double,4>  rad::beams::BeamIonComponents() {return {0.,0.,100.,0.93827210};}
inline constexpr std::array<double,4>  rad::beams::BeamEleComponents() {return {0.,0.,-10.,0.00051099900};}

void ProcessMCMatchedTCS(){
  gBenchmark->Start("df");

  rad::config::ePICReaction epic{"events", "/home/dglazier/EIC/data/sim/TCS_gen_ab_hiDiv_10x100m_s354_novtx.0070.eicrecon.tree.edm4eic.root"};
  epic.AliasColumnsAndMatchWithMC(false);
   
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  epic.setBeamIonIndex(rad::beams::BeamEleFix());
  epic.setBeamElectronIndex(rad::beams::BeamIonFix());
  epic.setScatElectronIndex(1);
  //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  epic.setParticleIndex("el",6,11);
  epic.setParticleIndex("po",7,-11);
  epic.setParticleIndex("p",5,2212);

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  epic.setBaryonParticles({"p"});
  
  epic.setMesonParticles({"el","po"});

  //must call this after all particles are configured
  epic.makeParticleMap();

  //option filtering of reconstructed tracks
  //  epic.Filter("el_OK==1&&po_OK==1","partFilter");
  
  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  rad::rdf::Mass(epic,"Whad","{el,po,p}");
  rad::rdf::Mass(epic,"MMass","{el,po}");

  //t distribution, column name
  rad::rdf::TBot(epic,"t_pn");
  rad::rdf::TPrime(epic,"tp_pn");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  auto df0 = epic.CurrFrame();

  auto hW = df0.Histo1D({"W","W",100,0,20.},"rec_W");
  auto hWhad = df0.Histo1D({"Whad","Whad",100,0,20.},"rec_Whad");
  auto hMesonMass = df0.Histo1D({"MesonMass","M(e-,e+) [GeV]",100,.3,5.},"rec_MMass");
  auto htpn = df0.Histo1D({"tpn","t(p,n) [GeV^{2}]",100,0,5},"rec_t_pn");
  auto htprimepn = df0.Histo1D({"tprimepn","t'(p,n) [GeV^{2}]",100,0,1},"rec_tp_pn");
  auto hthCM=df0.Histo1D({"cthCM","cos(#theta_{CM})",100,-1,1},"rec_CM_CosTheta");
  auto hphCM=df0.Histo1D({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},"rec_CM_Phi");
 
  gBenchmark->Start("processing");
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
  hWhad->DrawCopy();
  new TCanvas();
  hMesonMass->DrawCopy();
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
  
  epic.Snapshot("MCMatchedTCS.root");

}
