//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "HepMCElectro.h"
#include "ParticleCreator.h"
#include "ParticleGenerator.h"
#include "ParticleModifier.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

void ProcessHepMCDDVCS(const std::string infile="/w/work5/home/garyp/eic/Farm/data/EpIC_DDVCS_ee_18x275/rootfiles/18x275_ddvcs_edecay_hplus.root", const std::string outfile="tempout.root"){
  // Enable implicit multi-threading
  //ROOT::EnableImplicitMT(4);
 
  using namespace rad::names::data_type; //for MC()
  
  gBenchmark->Start("df");
  
  //create reaction dataframe
  rad::config::HepMCElectro hepmc{"hepmc3_tree", infile };
  
  hepmc.AliasMomentumComponents();
  
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  hepmc.setBeamIonIndex(3);
  hepmc.setBeamElectronIndex(0);
  hepmc.setScatElectronIndex(1);
  
  hepmc.setParticleIndex("pprime",5);
  hepmc.setBaryonParticles({"pprime"});
  
  //if using existing lepton pair
  //create gprime first from their sums
  rad::config::ParticleCreator particles{hepmc};
  hepmc.setParticleIndex("ele",6);
  hepmc.setParticleIndex("pos",7);
  particles.Sum("gprime",{"ele","pos"});
  
  //use ParticleGenerator to decay gprime to 2 M_ele particles ele,pos
  //hepmc.setParticleIndex("gprime",4);//
  // rad::generator::ParticleGenerator gen{hepmc};
  // gen.GenerateTwoBody({"ele0","pos0"},
  // 		      {rad::constant::M_ele(),rad::constant::M_ele()},
  // 		      "gprime");
  
  hepmc.setMesonParticles({"ele","pos"});
  //hepmc.setMesonParticles({"ele0","pos0"});

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
  
  //dis kinematics
  rad::rdf::Q2(hepmc,"Q2");
  rad::rdf::nu(hepmc,"nu");
  rad::rdf::y(hepmc,"y");
  rad::rdf::xbj(hepmc,"xbj");
  
  //q2prime tau
  rad::rdf::TauPrime(hepmc,"tau");
  
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
  rad::rdf::CircPolGammaStar(hepmc,"GammaPolCirc");
  rad::rdf::EGammaStar(hepmc,"GammaE");
  
  ///////////////////////////////////////////////////////////////
  //Define subsets of particles and corresponing variables to plot
  ///////////////////////////////////////////////////////////////
  // hepmc.Define("electrons","rad::helpers::PositionsWhere(mc_pid==11)");
  // hepmc.Define(MC()+"elsP",Form("Take(%spmag,electrons)",MC().data()));
  
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
  histo.Create<TH1D,double>({"nu","nu",100,0,10000.},{"nu"});
  histo.Create<TH1D,double>({"xbj","xbj",100,0,1.},{"xbj"});
  histo.Create<TH1D,double>({"y","y",100,0,1.},{"y"});
  
  histo.Create<TH1D,double>({"W","W",100,0,200.},{"W"});
  histo.Create<TH1D,double>({"Whad","Whad",100,0,200.},{"Whad"});
  
  histo.Create<TH1D,double>({"MesonMass","M(e^{-},e^{+}) [GeV]",100,0,5.},{"GMass"});
  
  histo.Create<TH1D,double>({"ttop","t(p,p^{'}) [GeV^{2}]",100,0,2.0},{"t_top"});
  histo.Create<TH1D,double>({"tbot","t(p,p^{'}) [GeV^{2}]",100,0,2.0},{"t_bot"});
  histo.Create<TH1D,double>({"tptop","t' top [GeV^{2}]",100,0,2.0},{"tp_top"});
  histo.Create<TH1D,double>({"tpbot","t' bot [GeV^{2}]",100,0,2.0},{"tp_bot"});

  histo.Create<TH2D,double,double>({"tbot_mesonmass","tbot vs meson mass",100,0.0,0.2,100,0.3,5.},{"t_bot","GMass"});
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
  histo.Create<TH1D,double>({"MissMass","Mmiss [GeV]",100,-10,10},{"MissMass"});
  histo.Create<TH1D,double>({"missP","p_{miss}(e',#gamma',p')",100,0,100},{"MissP"});
  histo.Create<TH1D,double>({"missPt","p_{t,miss}(e',#gamma',p')",100,0,10},{"MissPt"});
  histo.Create<TH1D,double>({"missPz","p_{z,miss}(e',#gamma',p')",100,0,100},{"MissPz"});
  histo.Create<TH1D,double>({"missTheta","#theta_{miss}(e',#gamma',p')",100,0,1},{"MissTheta"});
  
  //semi-exclusivity
  histo.Create<TH1D,double>({"MissMassPprime","Mmiss {e,p'}[GeV]",100,.3,5.},{"MissMassPprime"});
  
  histo.Create<TH1D,ROOT::RVecD>({"allP","momentum of all particles",100,0,100},{"pmag"});
  // histo.Create<TH1D,ROOT::RVecD>({"eleP","momentum of electrons",100,0,100},{"elsP"});
    
  //check recoil proton azimuthal distribution
  histo.Create<TH1D,double>({"scatele_phi","Azimuthal Angle of Scattered Electron",250,-TMath::Pi(),TMath::Pi()},{"phi[scat_ele]"});
  histo.Create<TH1D,double>({"pprime_phi","Azimuthal Angle of Recoil Proton",250,-TMath::Pi(),TMath::Pi()},{"phi[pprime]"});
  
  //for brufit need
  //CM_Phi Heli_theta Heli_phi GammaPolCirc=sqrt(1-epsilon)*Pol t GammaE
  histo.Create<TH1D,double>({"GammaPol","Polarisation of Virtual Photon",100,0,1},{"GammaPol"});
  histo.Create<TH1D,double>({"GammaE","Energy of Virtual Photon",100,0,18},{"GammaE"});
  histo.Create<TH1D,double>({"Heli_CosTheta","#theta decay angle",100,-TMath::Pi(),TMath::Pi()},{"Heli_CosTheta"});
  histo.Create<TH1D,double>({"Heli_Phi","#phi decay angle",100,-TMath::Pi()-1,TMath::Pi()+1},{"Heli_Phi"});
 
  //Polarisation and kinematic limits for cicular
  histo.Create<TH2D,double,double>({"y_W","W{ele missmass} vs y",100,0,1,100,0,200},{"y","W"});
  histo.Create<TH2D,double,double>({"y_Escatele","E_{e'} vs y",100,0,1,100,0,18},{"y","pmag[scat_ele]"});
  histo.Create<TH2D,double,double>({"y_CircPol","Circular Polarisation vs y",100,0,1,100,0,1},{"y","GammaPolCirc"});
  histo.Create<TH2D,double,double>({"W_CircPol","Circular Polarisation vs W{ele missmiass}",100,0,200,100,0,1},{"W","GammaPolCirc"});
  
  
  gBenchmark->Start("snapshot");
  hepmc.BookLazySnapshot(outfile);
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");

  gBenchmark->Start("processing");
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
  
  TCanvas *c00 = new TCanvas("c00","Kinematics");
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
  c00->Print("temp00.pdf");

  TCanvas *c001 = new TCanvas("c001","Kinematics Reduced");
  c001->Divide(2,2);
  c001->cd(1)->SetLogy();
  histo.DrawSame("Q2",gPad);
  c001->cd(2);
  histo.DrawSame("W",gPad);
  c001->cd(3);
  histo.DrawSame("MesonMass",gPad);
  c001->cd(4);
  histo.DrawSame("tbot",gPad);
  c001->Print("temp001.pdf");

  //histo.DrawSame("ttop",gPad);
  
  TCanvas *c01 = new TCanvas("c01","Exclusivity Plots");
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
  c01->Print("temp01.pdf");

  TCanvas *c02 = new TCanvas("c02","2D t Distributions");
  c02->Divide(2,2);
  c02->cd(1);
  histo.Draw2D("tbot_mesonmass",MC(),"colz",gPad);
  c02->cd(2);
  histo.Draw2D("tpbot_mesonmass",MC(),"colz",gPad);
  c02->cd(3);
  histo.Draw2D("tbot_W",MC(),"colz",gPad);
  c02->cd(4);
  histo.Draw2D("tpbot_W",MC(),"colz",gPad);
  c02->Print("temp02.pdf");

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
  c03->Print("temp03.pdf");

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
  c04->Print("temp04.pdf");

  TCanvas *c05 = new TCanvas("c05","DIS Kinematics");
  c05->Divide(2,2);
  c05->cd(1);
  histo.DrawSame("Q2",gPad);
  c05->cd(2);
  histo.DrawSame("nu",gPad);
  c05->cd(3);
  histo.DrawSame("xbj",gPad);
  c05->cd(4);
  histo.DrawSame("y",gPad);
  c05->Print("temp05.pdf");

  TCanvas *c06 = new TCanvas("c06","Kin Lim CircPol Correlations");
  c06->Divide(2,2);
  c06->cd(1);
  histo.Draw2D("y_W",MC(),"colz",gPad);
  c06->cd(2);
  histo.Draw2D("y_Escatele",MC(),"colz",gPad);
  c06->cd(3);
  histo.Draw2D("y_CircPol",MC(),"colz",gPad);
  c06->cd(4);
  histo.Draw2D("W_CircPol",MC(),"colz",gPad);
  c06->Print("temp06.pdf");

  TString pdfoutfile = outfile;
  pdfoutfile = gSystem->BaseName(pdfoutfile);
  pdfoutfile.ReplaceAll(".root",".pdf");
  
  gSystem->Exec(Form("pdfunite temp*.pdf %s",pdfoutfile.Data()));
  gSystem->Exec("rm temp*.pdf");
  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");

  //save all histograms to file
  //histo.File("HepMCDDVCS_hists.root");

  gBenchmark->Stop("df");
  gBenchmark->Print("df");
  
}
