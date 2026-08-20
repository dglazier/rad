#include "AnalysisManager.h"
#include "HepMCElectro.h"
#include "BasicKinematicsRDF.h"
#include "KinematicsProcElectro.h"
#include "ElectronScatterKinematics.h"
#include "gammaN_2_Spin0Spin0SpinHalf.h"
#include "DefineNames.h"
#include "TCSKinematics.h"


// #include "CommonDefines.h"
// #include "HepMCElectro.h"
// #include "KinematicsProcElectro.h"
// #include "KineCalculation.h"
// #include "Indicing.h"
// #include "ElectronScatterKinematics.h"
// #include "gammaN_2_Spin0Spin0SpinHalf.h"
// #include "TCSKinematics.h"
// #include <TBenchmark.h>
 
/**
 * @brief Example Script: TCS Analysis with Combinatorics and Missing Mass.
 * Updated to use SetMesonParticles shortcut and CloneLinked.
 */
void ProcessHepMCTCSCombi(
			  std::string outdir = "tcs",
			  const int beam_ele_idx = 0,
			  const int beam_ion_idx = 1,
			  const int Role_ScatEle  = 2, 
			  const int Role_Recoil   = 3,
			  const int Role_DecayEle = 4, 
			  const int Role_DecayPos = 5, 
			  const int lep_PDG = 11) 
{
  using namespace rad;
  using namespace rad::consts::data_type; 
  ROOT::EnableImplicitMT();

  gBenchmark->Start("df");

  using Reaction = HepMCElectro;
  using Processor = KinematicsProcElectro;
  
  AnalysisManager<Reaction,Processor>  mgr{
    "TCS",
    "hepmc3_tree",
  "~/Dropbox/elSpectroTesting/out/root/dilep_phsp_10_100.root"};
  
  gSystem->Exec(Form("mkdir -p %s",outdir.c_str()));
  mgr.SetOutputDir(outdir);
  auto& rad_df = mgr.Reaction();

  rad_df.Define("impWeight","weights[4]");
  rad_df.SetupMC();
   
  rad_df.SetBeamIonIndex(beam_ion_idx);
  rad_df.SetBeamElectronIndex(beam_ele_idx);
  rad_df.SetParticleIndex(consts::ScatEle(), Role_ScatEle);
  rad_df.SetParticleIndex("ele", Role_DecayEle); 
  rad_df.SetParticleIndex("pos", Role_DecayPos); 
  rad_df.SetParticleIndex("pprime", Role_Recoil); 

  rad_df.MakeCombinations();

  mgr.AddStream(MC(),""); 
  
  // [A] Topology
  auto topology_recipe =  [](Processor& p) {
    using namespace rad::consts;

    p.Creator().Sum("gprime", {{"ele", "pos"}});       
    p.Creator().Diff("miss",    {{BeamEle(),BeamIon()},{ScatEle(), "gprime", "pprime"}});
    p.Creator().Diff("miss",    {{BeamEle(),BeamIon()},{ScatEle(), "gprime", "pprime"}});
    p.Creator().Diff("miss_hadro",    {{BeamEle(),BeamIon()},{ScatEle()}});
    p.Creator().Diff("miss_pprime",    {{BeamEle(),BeamIon()},{ScatEle(), "gprime"}});
    p.Creator().Diff("delta",{{BeamIon()},{"pprime"}});
  
    p.SetMesonParticles({"ele","pos"});
    p.SetBaryonParticles({"pprime"});
    //p.SetBaryonParticles({"Miss"});
  
    // 2. Calculate kinematic variables
    p.Q2();         
    p.xbj();
    p.y();         
    p.nu();
    p.tau();
    p.tauprime();

    p.RegisterCalc("GammaPol",rad::physics::ElS_PolVirtPhot);
    p.RegisterCalc("GammaPolCirc",rad::physics::ElS_CircPolVirtPhot);

    p.RegisterCalc("CosThetaHel",rad::gn2s0s0s12::CosThetaHel);
    p.RegisterCalc("ThetaHel",rad::gn2s0s0s12::ThetaHel);
    p.RegisterCalc("PhiHel",rad::gn2s0s0s12::PhiHel);
  
    //need to add proton rest PR production angles
    p.CosThetaCM(); 
    //p.ThetaCM();
    p.PhiCM();       
    //p.CosThetaPR();
    //p.ThetaPR();
    //p.PhiPR();
    
    // 3. Calculate Invariant Masses
    p.Mass("Qp",{"gprime"});             
    p.Mass("GMass",{"gprime"});             
    p.Mass2("Qp2",{"gprime"});
    p.Mass2("s_photo",{VirtGamma(),BeamIon()});
    p.Mass("W",{"miss_hadro"});
    p.Mass2("W2",{"miss_hadro"});
    p.Mass("Whad",{"gprime","pprime"});
    p.Mass2("Whad2",{"gprime","pprime"});
    
    // 3b. Exclusivity
    p.Mass("MissMass",{"miss"});
    p.Mass2("MissMass2", {"miss"});             
    //p.P("MissP",{"miss"});
    p.Pt("MissPt",{"miss"});
    //p.Pz("MissPz",{"miss"});
    //p.Theta("MissTheta",{"miss"});
    
    // 3c. Missing proton (later split to own toplogy recipe)
    p.Mass("MissMass_pprime",{"miss_pprime"});
    p.Mass2("MissMass2_pprime", {"miss_pprime"}); // better to use miss or particles?      
    //p.P("MissP_pprime",{"miss_pprime"});
    p.Pt("MissPt_pprime",{"miss_pprime"});
    //p.Pz("MissPz_pprime",{"miss_pprime"});
    //p.Theta("MissTheta_pprime",{"miss_pprime"});
  
    
    //4. Calculate Mandelstam t (requires beam definition)
    p.RegisterCalc("t_top", rad::physics::TTop);
    p.RegisterCalc("t_bot", rad::physics::TBot);
    p.RegisterCalc("DeltaT", rad::physics::DeltaTBot);
  
    //5. Misc Kinematics Variables for cuts etc
    p.RegisterCalc("DeltaPhiProton", rad::DeltaPhi,{{"pprime","miss_pprime"}});
    p.Energy("GammaE", {VirtGamma()});
    
    //6. Particle Basic Observables
    p.ParticleTheta({"scat_ele","ele","pos","gprime","pprime"});
    p.ParticlePhi({"scat_ele","ele","pos","gprime","pprime"});
    p.ParticleP({"scat_ele","ele","pos","gprime","pprime"});
    p.ParticleEta({"scat_ele","ele","pos","gprime","pprime"});

  };
  
  // [B] REC-SPECIFIC CORRECTIONS
  auto correction_recipe  = [](rad::KinematicsProcessor& p) {
   using namespace rad::consts;
   p.PreModifier().FixMass(ScatEle(), M_ele());
    p.PreModifier().FixMass("ele", M_ele());
    p.PreModifier().FixMass("pos", M_ele());
    p.PreModifier().FixMass("pprime", M_pro());
  };
  
  // [C] SELECTION CUTS
  auto selection_recipe  = [](rad::PhysicsSelection& s) {
    
    // At least reconstruct each particle
    s.AddCutRange("scatele_pmag_cut", "scat_ele_pmag", 0, 18); 
    s.AddCutMin("gprime_pmag_cut", "gprime_pmag", 0);
    s.AddCutRange("pprime_pmag_cut", "pprime_pmag", 0, 275);
    
    s.AddCutRange("Qp2_cut",    "Qp2", 0.0, 20.0); 
    
  };
  
  // [D] HISTOGRAMS
  auto histogram_recipe =   [](rad::histo::Histogrammer& h) {
    //kinematics
    h.Create("Q2",";Q^2[GeV^{2}]", 100, 0, 1.0, "Q2");
    h.Create("nu",";#nu [GeV]",100,0,10000.,"nu");
    h.Create("xbj",";xbj",100,0,1.,"xbj");
    h.Create("y",";y",100,0,1.,"y");
    
    h.Create("W",";W (electro) [GeV]",100,0,200.,"W");
    h.Create("Whad","W (hadronic) [GeV]",100,0,200.,"Whad");
    
    h.Create("Qp", ";Q' = M_{e+e-} [GeV]", 100, 0.0, 5.0, "Qp");
    h.Create("Qp2", ";Q'^{2} = M_{e+e-} [GeV^2]", 100, 0.0, 25.0, "Qp2");
    
    h.Create("ttop","t(p,p^') [GeV^2]",100,0,2.0,"t_top");
    h.Create("tbot","t(p,p^') [GeV^2]",100,0,2.0,"t_bot");
    // h.Create("tptop","t' top [GeV^2]",100,0,2.0,"tp_top");
    // h.Create("tpbot","t' bot [GeV^2]",100,0,2.0,"tp_bot");
    
    h.Create2D("ttop_Qp2",";t_{eXBE} [GeV^{2}]; Q'^{2} [GeV^{2}]",100, 0.0, 2.0, 100, 0.0, 25.0,"t_top","Qp2");
    h.Create2D("tbot_Qp2",";t_{BABE} [GeV^{2}]; Q'^{2} [GeV^{2}]",100, 0.0, 2.0, 100, 0.0, 25.0,"t_bot","Qp2");
    
    h.Create2D("httop_W2",";t_{eXBE} [GeV^{2}]; W^{2} [GeV^{2}]",100, 0.0, 2.0, 100, 0.0, 200.0, "t_top","W2");
    h.Create2D("htbot_W2had",";t_{BABE} [GeV^{2}]; W^{2} [GeV^{2}]",100, 0.0, 2.0, 100, 0.0, 200.0, "t_bot", "W2");
    
    h.Create2D("W2_Qp2","",100, 0.0, 25., 100, 0.0, 200.0, "Qp2","W");
 
    //CM and PR Decay Angles
    h.Create("CosThetaCM","cos(#theta_CM)",100,-1,1,"CosThetaCM");
    h.Create("PhiCM","#phi_CM",100,-TMath::Pi(),TMath::Pi(),"PhiCM");
    //h.Create("CosThetaPR","cos(#theta_{PR})",100,-1,1,"CosThetaPR");
    //h.Create("PhiPR","#phi_{PR}",100,-TMath::Pi(),TMath::Pi(),"PhiPR");
  
    //exclusivity
    h.Create("MissMass","Mmiss [GeV]",100,-10,10,"MissMass");
    h.Create("MissMass2",     "Missing Mass squared; [GeV]", 100, -50, 50, "MissMass2");
    //h.Create("missP","p_miss(e',#gamma',p')",100,0,100,"miss_pmag");
    h.Create("missPt","p_t,miss(e',#gamma',p')",100,0,10,"MissPt");
    //h.Create("missPz","p_z,miss(e',#gamma',p')",100,0,100,"MissPz");
    //h.Create("missTheta","#theta_miss(e',#gamma',p')",100,0,1,"miss_theta");
  
    //semi-exclusivity
    h.Create("MissMass_pprime","Missing Mass of Proton; M_{miss,p'} [GeV/c^{2}]",100,-5., 5., "MissMass_pprime");
    h.Create("MissMass2_pprime", "Missing Mass of Proton squared; M^{2}_{miss,p'} [GeV^{2}/c^{4}]", 100, -25., 25., "MissMass2_pprime");
    
    //scattered electron
    h.Create("scat_ele_pmag","Momentum of Scattered Electron; p_{e'} [GeV/c]",100,0,18,"scat_ele_pmag");
    h.Create("scat_ele_eta","Pseudorapidity of Scattered Electron; #eta_{e'}",100,-10,10,"scat_ele_eta");
    h.Create("scat_ele_theta","Polar Angle of Scattered Electron; #theta_{e'} [rad]",100, 0.0,TMath::Pi(),"scat_ele_theta");
    h.Create("scat_ele_phi","Azimuthal Angle of Scattered Electron; #phi_{e'} [rad]",100,-TMath::Pi(),TMath::Pi(),"scat_ele_phi");
  
    //decay ele
    h.Create("ele_pmag","Momentum of Decay Ele; p_{e-} [GeV/c]",100,0,275,"ele_pmag");
    h.Create("ele_eta","Pseudorapidity of Decay Ele; #eta_{e-}",100,-10,10,"ele_eta");
    h.Create("ele_theta","Polar Angle of Decay Ele; #theta_{e-} [rad]",100, 0.0,TMath::Pi(),"ele_theta");
    h.Create("ele_phi","Azimuthal Angle of Decay Ele; #theta_{e-} [rad]",100,-TMath::Pi(),TMath::Pi(),"ele_phi");
    
    //decay pos
    h.Create("pos_pmag","Momentum of Decay Pos; p_{e+} [GeV/c]",100,0,275,"pos_pmag");
    h.Create("pos_eta","Pseudorapidity of Decay Pos; #eta_{e+}",100,-10,10,"pos_eta");
    h.Create("pos_theta","Polar Angle of Decay Pos; #theta_{e+} [rad]",100, 0.0,TMath::Pi(),"pos_theta");
    h.Create("pos_phi","Azimuthal Angle of Decay Pos; #phi_{e+} [rad]",100,-TMath::Pi(),TMath::Pi(),"pos_phi");
    
    //recoil proton
    h.Create("pprime_pmag","Momentum of Recoil Proton; p_{p'} [GeV/c]",100,0,275,"pprime_pmag");
    h.Create("pprime_eta","Pseudorapidity of Recoil Proton; #eta_{p'}",100,-10,10,"pprime_eta");
    h.Create("pprime_theta","Polar Angle of Recoil Proton; #theta_{p'} [rad]",100, 0.0,TMath::Pi(),"pprime_theta");
    h.Create("pprime_phi","Azimuthal Angle of Recoil Proton; #phi_{p'} [rad]",100,-TMath::Pi(),TMath::Pi(),"pprime_phi");
  
    //for brufit need
    //CM_Phi Heli_theta Heli_phi GammaPolCirc=sqrt(1-epsilon)*Pol t GammaE
    h.Create("GammaPol","Polarisation of Virtual Photon;#epsilon",100,0,1,"GammaPol");
    h.Create("GammaE","Energy of Virtual Photon;E_{#gamma} [GeV]",100,0,18,"GammaE");
    h.Create("CosThetaHel","CosTheta Decay Angle;cos(#theta_{l})",100,-1,1,"CosThetaHel");
    h.Create("ThetaHel","Theta Decay Angle;#theta_{l} [rad]",100, 0.0,TMath::Pi(),"ThetaHel");
    h.Create("PhiHel","Phi Decay Angle;#phi_{l} [rad]",100,-TMath::Pi()-1,TMath::Pi()+1,"PhiHel");
 
    //Polarisation and kinematic limits for cicular
    h.Create2D("y_W","Wele missmass vs y",100,0,1,100,0,200,"y","W");
    h.Create2D("y_Escatele","E_e' vs y",100,0,1,100,0,18,"y","scat_ele_pmag");
    h.Create2D("y_CircPol","Circular Polarisation vs y",100,0,1,100,0,1,"y","GammaPolCirc");
    h.Create2D("W_CircPol","Circular Polarisation vs Wele missmiass",100,0,200,100,0,1,"W","GammaPolCirc");

  };
  mgr.ConfigureKinematics(topology_recipe);
  //mgr.ConfigureKinematics(correction_recipe);
  //mgr.ConfigureSelection(selection_recipe);
  mgr.ConfigureHistograms(histogram_recipe);
  
  //rad::rdf::PrintParticles(rad_df, MC());
  //rad::PrintDefinedColumnNames(mgr.Reaction().CurrFrame());
  mgr.Snapshot({"impWeight"});

  //Diagnostics helpers
  // rad::rdf::PrintParticles(rad_df, MC());
  rad::PrintDefinedColumnNames(mgr.Reaction().CurrFrame());
  mgr.PrintDiagnostics();
  //PrintDefinedColumnNames will give a list of all columns we can print at this point
  //optional 3rd argument = number of events to print for
  
  
  gBenchmark->Start("analysis");
  mgr.Run();
  gBenchmark->Stop("analysis");
  gBenchmark->Print("analysis");
}


 
