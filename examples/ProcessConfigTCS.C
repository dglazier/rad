//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "include/ePICReaction.h"
#include "include/Beams.h"
#include "include/Indicing.h"
#include "include/BasicKinematicsRDF.h"
#include "include/ReactionKinematicsRDF.h"
#include "include/ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

inline constexpr std::array<double,4>  rad::beams::BeamIonComponents() {return {0.,0.,100.,0.93827210};}
inline constexpr std::array<double,4>  rad::beams::BeamEleComponents() {return {0.,0.,-10.,0.00051099900};}

void ProcessConfigTCS(){
  /* Make sure beam energies are correctly set above prior to running */
  
  gBenchmark->Start("df");
  cout<<rad::beams::BeamIonComponents()[2]<<endl;
  
  rad::config::ePICReaction epic{"events", "/home/dglazier/EIC/data/sim/TCS_gen_ab_hiDiv_10x100m_s354_novtx.0070.eicrecon.tree.edm4eic.root"};
  epic.AliasColumns();
  //Assign particles names and indices
  epic.setBeamIonIndex(rad::beams::BeamIonFix());//change to taking mcparticles fourvector by index
  epic.setBeamElectronIndex(rad::beams::BeamEleFix());
  
  //2nd occurance of 11 in rec_pid
  epic.setScatElectron( rad::indice::useNthOccurance(2,11) ,{"rec_pid"} );

  //1st occurance of 11 in rec_pid ordered by pmag
  //epic.SetScatElectronIndex( indice::useNthOccuranceSortedBy(1,11) ,{"rec_pid","pmag"} );
  
  //1st occurance of 11 in rec_pid
  epic.setParticleIndex("el",rad::indice::useNthOccurance(1,11) ,{"rec_pid"},11 );
  epic.setParticleIndex("po",rad::indice::useNthOccurance(1,-11) ,{"rec_pid"},-11 );
  epic.setParticleIndex("p",rad::indice::useNthOccurance(1,2212) ,{"rec_pid"},2212 );
  
  epic.setBaryonParticles({"p"});
  epic.setMesonParticles({"el","po"});

  
  epic.makeParticleMap();


  //masses column name, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  //masses column name, {+ve particles}, {-ve particles}  
  rad::rdf::Mass(epic,"MMass","{el,po}");

  //t distribution, df, column name
  rad::rdf::TBot(epic,"t_pn");
  rad::rdf::TPrime(epic,"tp_pn");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  auto df0 = epic.CurrFrame();
  auto hW = df0.Histo1D({"W","W",100,0,70.},"rec_W");
  auto hMesonMass = df0.Histo1D({"MesonMass","M(e-,e+) [GeV]",100,.3,5.},"rec_MMass");
  auto htpn = df0.Histo1D({"tpn","t(p,n) [GeV^{2}]",100,0,5},"rec_t_pn");
  auto htprimepn = df0.Histo1D({"tprimepn","t'(p,n) [GeV^{2}]",100,0,1},"rec_tp_pn");
  auto hthCM=df0.Histo1D({"cthCM","cos(#theta_{CM})",100,-1,1},"rec_CM_CosTheta");
  auto hphCM=df0.Histo1D({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},"rec_CM_Phi");
 
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
  hW->DrawCopy();
  new TCanvas();
  hMesonMass->DrawCopy();
  new TCanvas();
  htprimepn->DrawCopy();
  htprimepn->DrawCopy("same");
  auto canCM = new TCanvas();
  canCM->Divide(2,1);
  canCM->cd(1);
  hthCM->DrawCopy();
  canCM->cd(2);
  hphCM->DrawCopy();

  gBenchmark->Stop("df");
  gBenchmark->Print("df");
}
