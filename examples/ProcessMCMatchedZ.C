//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "ePICReaction.h"
#include "ParticleCreator.h"
#include "Indicing.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

//inline constexpr std::array<double,4>  rad::beams::InitBotComponents() {return {0.,0.,100.,0.93827210};}
//inline constexpr std::array<double,4>  rad::beams::InitTopComponents() {return {0.,0.,-10.,0.00051099900};}

//with afterburner need slightly altered energies
inline constexpr std::array<double,4>  rad::beams::InitBotComponents() {return {0,0,99.9339,0.938272};}
inline constexpr std::array<double,4>  rad::beams::InitTopComponents() {return {0,0,-10.007,0.000510999};}

void ProcessMCMatchedZ(){
  gBenchmark->Start("df");

  rad::config::ePICReaction epic{"events","/home/dglazier/EIC/data/sim/jpac_z3900_10x100.root"};
  //epic.AliasColumnsAndMC();
  epic.AliasColumnsAndMatchWithMC();

  
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  epic.setBeamIonIndex(rad::beams::InitBotFix());
  epic.setBeamElectronIndex(rad::beams::InitTopFix());
  epic.setScatElectronIndex(4);
  //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  epic.setParticleIndex("idxEl",0,11);
  epic.setParticleIndex("idxPo",1,-11);
  epic.setParticleIndex("idxPi",2,211);
  epic.setParticleIndex("idxN",3,2112);

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  //but must be done after ParticleMap
  //so currently cannot use as baryon...

  rad::config::ParticleCreator particles{epic};
  particles.Sum("idxJ",{"idxEl","idxPo"});
  particles.Sum("idxZ",{"idxPi","idxJ"});
  
  epic.setMesonParticles({"idxEl","idxPo","idxPi"});
  
  //can also add missing particles
  //And use those in calculations
  //Miss(name,{other final state particles})
  particles.Miss("idxCalcN",{rad::names::ScatEle().data(),"idxZ"});
  epic.setBaryonParticles({"idxCalcN"});

  //must call this after all particles are configured
  epic.makeParticleMap();
  
 
  //option filtering of reconstructed tracks
  //  epic.Filter("el_OK==1&&po_OK==1","partFilter");
  
  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  rad::rdf::MissMass(epic,"MissMass","{scat_ele,idxN,idxZ}");
  rad::rdf::Mass(epic,"Whad","{idxPi,idxEl,idxPo,idxN}");
  //  rad::rdf::Mass(epic,"Whad","{idxZ,idxN}");
  rad::rdf::Mass(epic,"JMass","{idxJ}");
  rad::rdf::Mass(epic,"ZMass","{idxZ}");
  rad::rdf::Mass(epic,"MissNMass","{idxCalcN}");

  //t distribution, column name
  rad::rdf::TTop(epic,"t_gZ");
  rad::rdf::TBot(epic,"t_pn");
  rad::rdf::TPrimeBot(epic,"tp_pn");
  rad::rdf::TPrimeTop(epic,"tp_gZ");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");
  
  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  auto df0 = epic.CurrFrame();

  auto hW = df0.Histo1D({"W","W",100,0,20.},"tru_W");
  auto hWhad = df0.Histo1D({"Whad","Whad",100,0,20.},"tru_Whad");
  auto hMesonMass = df0.Histo1D({"MesonMass","M(e-,e+, #pi) [GeV]",100,.3,5.},"tru_ZMass");
  auto hMissMass = df0.Histo1D({"MissMass","Mmiss [GeV]",1000,-10,10},"tru_MissMass");
  auto hJMass = df0.Histo1D({"JMass","M(e-,e+) [GeV]",100,.3,5.},"tru_JMass");
  auto htpn = df0.Histo1D({"tpn","t(p,n) [GeV^{2}]",100,-2,5},"tru_t_pn");
  auto htgZ = df0.Histo1D({"tgZ","t(g,Z) [GeV^{2}]",100,-2,5},"tru_t_gZ");
  auto htprimepn = df0.Histo1D({"tprimepn","t'(p,n) [GeV^{2}]",100,-2,5},"tru_tp_pn");
  auto htprimegZ = df0.Histo1D({"tprimegZ","t'(p,n) [GeV^{2}]",100,-2,5},"tru_tp_gZ");
  auto hthCM=df0.Histo1D({"cthCM","cos(#theta_{CM})",100,-1,1},"tru_CM_CosTheta");
  auto hphCM=df0.Histo1D({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},"tru_CM_Phi");
  
  gBenchmark->Start("processing");
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
  //hp->DrawCopy();
  new TCanvas();
  hWhad->DrawCopy();
  new TCanvas();
  hMesonMass->DrawCopy();
  new TCanvas();
  hMissMass->DrawCopy();
  new TCanvas();
  hJMass->DrawCopy();
  new TCanvas();
  htpn->DrawCopy();
  htgZ->DrawCopy("same");
  htprimepn->DrawCopy("same");
  htprimegZ->DrawCopy("same");
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
  
  gBenchmark->Start("snapshot");
  // epic.Snapshot("MCMatchedZ.root");
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");

}
