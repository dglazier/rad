//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "include/ePICReaction.h"
#include "include/Indicing.h"
#include "include/BasicKinematicsRDF.h"
#include "include/ReactionKinematicsRDF.h"
#include "include/ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

void ProcessConfigDISRes(){

  // rad::config::ePICReaction epic{"events", {"/home/dglazier/EIC/data/sim/pythia8NCDIS_18x275_minQ2_10_beamEffects_xAngle_-0.025_hiDiv_5.0001.eicrecon.tree.edm4eic.root","/home/dglazier/EIC/data/sim/pythia8NCDIS_18x275_minQ2_10_beamEffects_xAngle_-0.025_hiDiv_5.0001.eicrecon.tree.edm4eic.root"}};
  rad::config::ePICReaction epic{"events", "/home/dglazier/EIC/data/sim/pythia8NCDIS_18x275_minQ2_10_beamEffects_xAngle_-0.025_hiDiv_5.0001.eicrecon.tree.edm4eic.root"};
  
  // epic.AliasColumnsAndMC();
  epic.AliasColumnsAndMatchWithMC();

  // resolution functions append variable by res_ for column name
  // we can only use predefined or existing columns
  // px,py,pz,m,pmag,eta,theta,phi
  epic.ResolutionFraction("pmag");
  epic.Resolution("eta");
  epic.Resolution("phi");
  //can now use res columns like any other to calc variables
  epic.Define("deltaR","sqrt(res_eta*res_eta+res_phi*res_phi)");
  
  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  auto df0 = epic.CurrFrame();
  auto matchedPartTrackDeltaP = df0.Histo1D({"matchedPartTrackDeltaP","Delta P Between Matching Thrown and Reconstructed Charged Particle",400,-2.,2.},"res_pmag");
  auto matchedPartTrackDeltaEta = df0.Histo1D({"matchedPartTrackDeltaEta","Delta #eta Between Matching Thrown and Reconstructed Charged Particle",100, -0.25, 0.25},"res_eta");
  auto matchedPartTrackDeltaPhi = df0.Histo1D({"matchedPartTrackDeltaP","Delta #phi Between Matching Thrown and Reconstructed Charged Particle", 200, -0.2, 0.2},"res_phi");
  auto matchedPartTrackDeltaR = df0.Histo1D({"matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3},"deltaR");

  // ///////////////////////////////////////////////////////////
  // //Draw histograms
  // ///////////////////////////////////////////////////////////
  gBenchmark->Start("df");
  auto can = new TCanvas();
  can->Divide(2,2);
  can->cd(1);
  matchedPartTrackDeltaEta->DrawCopy();
  can->cd(2);
  matchedPartTrackDeltaPhi->DrawCopy();
  can->cd(3);
  matchedPartTrackDeltaP->DrawCopy();
  can->cd(4);
  matchedPartTrackDeltaR->DrawCopy();
    
  gBenchmark->Stop("df");
  gBenchmark->Print("df");
}
