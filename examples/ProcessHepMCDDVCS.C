//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "HepMCElectro.h"
#include "ParticleCreator.h"
#include "ParticleGenerator.h"
//#include "ParticleGeneratorRDF.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

void ProcessHepMCDDVCS(){
  using namespace rad::names::data_type; //for MC()
   gBenchmark->Start("df");

   //create reaction dataframe
   rad::config::HepMCElectro hepmc{"hepmc3_tree", "/w/work5/home/garyp/18x275_ddvcs_events_plus.root"};
  hepmc.AliasMomentumComponents();
    
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  hepmc.setBeamIonIndex(3);
  hepmc.setBeamElectronIndex(0);
  hepmc.setScatElectronIndex(1);
  
  hepmc.setParticleIndex("pprime",5);
  hepmc.setBaryonParticles({"pprime"});
  
  //if using existing lepton pair
  rad::config::ParticleCreator particles{hepmc};
  // hepmc.setParticleIndex("ele",6);
  // hepmc.setParticleIndex("pos",7);
  // particles.Sum("gprime",{"ele","pos"});
  
  //if regenerating lepton pair
  hepmc.setParticleIndex("gprime",4);//
  rad::generator::ParticleGenerator gen{hepmc};
  double m_e = 0.000511;
  ROOT::VecOps::RVec<double> masses;
  masses.push_back(m_e);
  masses.push_back(m_e);
  gen.GenerateTwoBody({"ele","pos"},masses,"gprime");
  hepmc.setMesonParticles({"ele","pos"});

  //must call this after all particles are configured
  hepmc.makeParticleMap();
  //rad::rdf::PrintParticles(hepmc,MC());
  
  //////////////////////////////////////////////////////////
  // Now define calculated variables
  // Note reconstructed variables will have rec_ prepended
  // truth variables will have tru_ prepended
  //////////////////////////////////////////////////////////

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(hepmc,"W","{scat_ele}");
  rad::rdf::Mass(hepmc,"Whad","{gprime,pprime}");
  rad::rdf::Mass(hepmc,"GMass","{gprime}");
  rad::rdf::Q2(hepmc,"Q2");
  
  //t distribution, column name
  rad::rdf::TTop(hepmc,"t_top");
  rad::rdf::TBot(hepmc,"t_bot");
  rad::rdf::TPrimeBot(hepmc,"tp_bot");
  rad::rdf::TPrimeTop(hepmc,"tp_top");

  //CM production angles
  rad::rdf::CMAngles(hepmc,"CM");
  //Proton Rest production angles
  rad::rdf::PRAngles(hepmc,"PR");
  
  //exclusivity
  rad::rdf::MissMass(hepmc,"MissMass","{scat_ele,pprime,gprime}");
  rad::rdf::MissP(hepmc,"MissP","{scat_ele,pprime,gprime}");
  rad::rdf::MissPt(hepmc,"MissPt","{scat_ele,pprime,gprime}");
  rad::rdf::MissPz(hepmc,"MissPz","{scat_ele,pprime,gprime}");
  rad::rdf::MissTheta(hepmc,"MissTheta","{scat_ele,pprime,gprime}");
  
  //
  rad::rdf::MissMass(hepmc,"MissMassPprime","{scat_ele,pprime}");
  
  //decay angles
  rad::rdf::gn2s0s0s12::HelicityAngles(hepmc,"Heli");
  //photon polarisation
  rad::rdf::PolGammaStar(hepmc,"GammaPol");
  rad::rdf::EGammaStar(hepmc,"GammaE");
  
  //hepmc.Filter("Q2_tru>0.15");
  ///////////////////////////////////////////////////////////////
  //Define subsets of particles and corresponing variables to plot
  ///////////////////////////////////////////////////////////////
  hepmc.Define("electrons","rad::helpers::PositionsWhere(mc_pid==11)");
  hepmc.Define(MC()+"elsP",Form("Take(%spmag,electrons)",MC().data()));
  hepmc.Define("e_helicity","1");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  rad::histo::Histogrammer histo{"set1",hepmc};
  // we can create many histograms by splitting events into
  // bins, where the same histogram is produced for the given bin
  // e.g. create 10 bins in mc_W between 4 and 54 GeV 
  //histo.Splitter().AddRegularDimension(MC()+"W", rad::histo::RegularSplits(10,4,54) );
  //can add as many split dimensions as we like
  //histo.Splitter().AddRegularDimension("xxx", rad::histo::RegularSplits(nbins,low,high) );
  histo.Init({MC()});//will create histograms for mc

  //Kinematics
  histo.Create<TH1D,double>({"Q2","Q2",100,0,2.},{"Q2"});
  histo.Create<TH1D,double>({"W","W",100,0,200.},{"W"});
  histo.Create<TH1D,double>({"Whad","Whad",100,0,200.},{"Whad"});
  
  histo.Create<TH1D,double>({"MesonMass","M(e-,e+) [GeV]",100,-1,5.},{"GMass"});
  
  histo.Create<TH1D,double>({"ttop","t(p,pprime) [GeV^{2}]",100,-0.1,1.5},{"t_top"});
  histo.Create<TH1D,double>({"tbot","t(g,gprime) [GeV^{2}]",100,-0.1,1.5},{"t_bot"});
  histo.Create<TH1D,double>({"tptop","t' top [GeV^{2}]",100,-0.1,1.5},{"tp_top"});
  histo.Create<TH1D,double>({"tpbot","t' bot [GeV^{2}]",100,-0.1,1.5},{"tp_bot"});

  histo.Create<TH2D,double,double>({"tbot_mesonmass","tbot vs meson mass",100,-0.1,1.5,100,0.3,5.},{"t_bot","GMass"});
  histo.Create<TH2D,double,double>({"tpbot_mesonmass","t'bot vs meson mass",100,-0.1,1.5,100,0.3,5.},{"tp_bot","GMass"});
  histo.Create<TH2D,double,double>({"tbot_W","tbot vs W",100,-0.1,1.5,100,0.,200.},{"t_bot","W"});
  histo.Create<TH2D,double,double>({"tpbot_W","t'bot vs W",100,-0.1,1.5,100,0.,200.},{"tp_bot","W"});
  histo.Create<TH2D,double,double>({"W_mesonmass","W vs meson mass",100,0.3,5.,100,0.,200.},{"GMass","W"});
  
  //CM and PR Decay Angles
  histo.Create<TH1D,double>({"cthCM","cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"phCM","#phi_{CM}",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  histo.Create<TH1D,double>({"cthPR","cos(#theta_{PR})",100,-1,1},{"PR_CosTheta"});
  histo.Create<TH1D,double>({"phPR","#phi_{PR}",100,-TMath::Pi(),TMath::Pi()},{"PR_Phi"});
  
  //exclusivity
  histo.Create<TH1D,double>({"MissMass","Mmiss [GeV]",1000,-10,10},{"MissMass"});
  histo.Create<TH1D,double>({"missP","p_{miss}(e',#gamma',p')",105,0,105},{"MissP"});
  histo.Create<TH1D,double>({"missPt","p_{t,miss}(e',#gamma',p')",100,0,10},{"MissPt"});
  histo.Create<TH1D,double>({"missPz","p_{z,miss}(e',#gamma',p')",105,0,105},{"MissPz"});
  histo.Create<TH1D,double>({"missTheta","#theta_{miss}(e',#gamma',p')",100,0,1},{"MissTheta"});
  
  //semi-exclusivity
  histo.Create<TH1D,double>({"MissMassPprime","Mmiss {e,p'}[GeV]",100,.3,5.},{"MissMassPprime"});
  
  histo.Create<TH1D,ROOT::RVecD>({"allP","momentum of all particles",100,0,100},{"pmag"});
  histo.Create<TH1D,ROOT::RVecD>({"eleP","momentum of electrons",100,0,100},{"elsP"});
    
  //check recoil proton azimuthal distribution
  histo.Create<TH1D,double>({"scatele_phi","Azimuthal Angle of Recoil Proton",250,-TMath::Pi(),TMath::Pi()},{"phi[pprime]"});
  histo.Create<TH1D,double>({"pprime_phi","Azimuthal Angle of Recoil Proton",250,-TMath::Pi(),TMath::Pi()},{"phi[pprime]"});
  
  //for brufit need
  //theta phi pol t eepEgamma
  histo.Create<TH1D,double>({"GammaPol","Polarisation of Virtual Photon",100,0,1},{"GammaPol"});
  histo.Create<TH1D,double>({"GammaE","Energy of Virtual Photon",100,0,18},{"GammaE"});
  histo.Create<TH1D,double>({"Heli_CosTheta","#theta decay angle",100,-TMath::Pi(),TMath::Pi()},{"Heli_CosTheta"});
  histo.Create<TH1D,double>({"Heli_Phi","#phi decay angle",100,-TMath::Pi()-1,TMath::Pi()+1},{"Heli_Phi"});
 
  
  gBenchmark->Start("snapshot");
  hepmc.BookLazySnapshot("HepMCDDVCS.root");
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");

  gBenchmark->Start("processing");
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
  
  TCanvas *c00 = new TCanvas("c00","Kinematics",800,400);
  c00->Divide(4,2);
  c00->cd(1);
  histo.DrawSame("Q2",gPad);
  c00->cd(2);
  histo.DrawSame("W",gPad);
  c00->cd(3);
  histo.DrawSame("Whad",gPad);
  c00->cd(4);
  histo.DrawSame("MesonMass",gPad);
  c00->cd(5);
  histo.DrawSame("ttop",gPad);
  c00->cd(6);
  histo.DrawSame("tbot",gPad);
  c00->cd(7);
  histo.DrawSame("tptop",gPad);
  c00->cd(8);
  histo.DrawSame("tpbot",gPad);
  
  TCanvas *c01 = new TCanvas("c01","Exclusivity Plots",800,400);
  c01->Divide(3,2);
  c01->cd(1);
  histo.DrawSame("MissMass",gPad);
  c01->cd(2);
  histo.DrawSame("missP",gPad);
  c01->cd(3);
  histo.DrawSame("missPt",gPad);
  c01->cd(4);
  histo.DrawSame("missPz",gPad);
  c01->cd(5);
  histo.DrawSame("missTheta",gPad);
  //c01->cd(6);
  //histo.DrawSame("MissMass",gPad);
  
    TCanvas *c02 = new TCanvas("c02","");
  c02->Divide(2,2);
  c02->cd(1);
  histo.DrawSame("tbot_mesonmass",gPad);
  c02->cd(2);
  histo.DrawSame("tpbot_mesonmass",gPad);
  c02->cd(3);
  histo.DrawSame("tbot_W",gPad);
  c02->cd(4);
  histo.DrawSame("tpbot_W",gPad);

  TCanvas *c03 = new TCanvas("c03","CM and PR Frame Angles");
  c03->Divide(2,2);
  c03->cd(1);
  histo.DrawSame("cthCM",gPad);
  c03->cd(2);
  histo.DrawSame("phCM",gPad);
  c03->cd(3);
  histo.DrawSame("cthPR",gPad);
  c03->cd(4);
  histo.DrawSame("phPR",gPad);
  
  TCanvas *c04 = new TCanvas("c04","Helicity and Polarisation");
  c04->Divide(2,2);
  c04->cd(1);
  histo.DrawSame("Heli_CosTheta",gPad);
  c04->cd(2);
  histo.DrawSame("Heli_Phi",gPad);
  c04->cd(3);
  histo.DrawSame("GammaPol",gPad);
  c04->cd(4);
  histo.DrawSame("GammaE",gPad);
  
  // histo.DrawSame("pprime_phi");
  // histo.DrawSame("MissMassPprime");
  // histo.DrawSame("W_mesonmass");
  
  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");

  //save all histograms to file
  histo.File("HepMCDDVCS_hists.root");

  gBenchmark->Stop("df");
  gBenchmark->Print("df");
  
}
