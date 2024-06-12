#include <TFile.h>
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TMath.h>
#include <TVector3.h>
#include <TVector2.h>

void ResolutionAnalysis_TTreeReader(TString infile="/home/dglazier/EIC/data/sim/pythia8NCDIS_18x275_minQ2_10_beamEffects_xAngle_-0.025_hiDiv_5.0001.eicrecon.tree.edm4eic.root"){
  gBenchmark->Start("tr");
  // Set output file for the histograms
  TFile *ofile = TFile::Open("ResolutionAnalysis_Out.root","RECREATE");

  // Analysis code will go here
  // Set up input file chain
  TChain *mychain = new TChain("events");
  mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);
  // mychain->Add(infile);

  // Initialize reader
  TTreeReader tree_reader(mychain);

  // Get Particle Information
  TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
  TTreeReaderArray<float> partMomX(tree_reader, "MCParticles.momentum.x");
  TTreeReaderArray<float> partMomY(tree_reader, "MCParticles.momentum.y");
  TTreeReaderArray<float> partMomZ(tree_reader, "MCParticles.momentum.z");
  TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

  // Get Reconstructed Track Information
  TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
  TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");

  // Get Associations Between MCParticles and ReconstructedChargedParticles
  TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
  TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");
    
  // Define Histograms
  TH1D *trackMomentumRes = new TH1D("trackMomentumRes","Track Momentum Resolution", 400, -2, 2);
 
  TH1D *matchedPartTrackDeltaEta = new TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructed Charged Particle; #Delta#eta", 100, -0.25, 0.25);
  TH1D *matchedPartTrackDeltaPhi = new TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2);
  TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3);
  TH1D *matchedPartTrackDeltaMom = new TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10);
  while(tree_reader.Next()) { // Loop over events
    for(unsigned int i=0; i<partGenStat.GetSize(); i++){ // Loop over thrown particles
	if(partGenStat[i] == 1){ // Select stable thrown particles
	    int pdg = TMath::Abs(partPdg[i]);
	    if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212){ // Look at charged particles (electrons, muons, pions, kaons, protons)
		TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);

		float trueEta = trueMom.PseudoRapidity();
		float truePhi = trueMom.Phi();
	    
		for(unsigned int j=0; j<simuAssoc.GetSize(); j++){ // Loop over associations to find matching ReconstructedChargedParticle
		    if(simuAssoc[j] == i){ // Find association index matching the index of the thrown particle we are looking at
			TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle

			// Check the distance between the thrown and reconstructed particle
			float deltaEta = trueEta - recMom.PseudoRapidity();
			float deltaPhi = TVector2::Phi_mpi_pi(truePhi - recMom.Phi());
			float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
			float deltaMom = ((trueMom.Mag()) - (recMom.Mag()));
			double momRes = (recMom.Mag()- trueMom.Mag())/trueMom.Mag();
      
			trackMomentumRes->Fill(momRes);

			matchedPartTrackDeltaEta->Fill(deltaEta);
			matchedPartTrackDeltaPhi->Fill(deltaPhi);
			matchedPartTrackDeltaR->Fill(deltaR);
			matchedPartTrackDeltaMom->Fill(deltaMom);
                    }
                } // End loop over associations 
            } // End PDG check          
        } // End stable particles condition  
    } // End loop over thrown particles
  } // End loop over events 
  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
  gBenchmark->Stop("tr");
  gBenchmark->Print("tr");
}
