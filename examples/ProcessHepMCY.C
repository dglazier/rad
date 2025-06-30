//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "HepMCElectro.h"
#include "ParticleCreator.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

//beam components must be defined even if 0
inline constexpr std::array<double,4>  rad::beams::InitBotComponents() {return {0.,0.,0,0};}
inline constexpr std::array<double,4>  rad::beams::InitTopComponents() {return {0.,0.,0.,0};}

void ProcessHepMCZ(){
  using namespace rad::names::data_type; //for MC()
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
  hepmc.setParticleIndex("ele",2);
  hepmc.setParticleIndex("pos",3);
  hepmc.setParticleIndex("pip",4);
  hepmc.setParticleIndex("n",5);

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calculating reaction kinematics
  //e.g. t distributions
  hepmc.setBaryonParticles({"n"});

  rad::config::ParticleCreator particles{hepmc};
  particles.Sum("Jpsi",{"ele","pos"});
  particles.Sum("Zc",{"pip","Jpsi"});
  hepmc.setMesonParticles({"Jpsi","pip"});

  //must call this after all particles are configured
  hepmc.makeParticleMap();
  

  //////////////////////////////////////////////////////////
  // Now define calculated variables
  // Note reconstructed variables will have rec_ prepended
  // truth variables will have tru_ prepended
  //////////////////////////////////////////////////////////

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(hepmc,"W","{scat_ele}");
  rad::rdf::Mass(hepmc,"Whad","{Zc,n}");
  rad::rdf::Mass(hepmc,"JMass","{Jpsi}");
  rad::rdf::Mass(hepmc,"ZMass","{Zc}");
  rad::rdf::Q2(hepmc,"Q2");

  //t distribution, column name
  rad::rdf::TTop(hepmc,"t_gZ");
  rad::rdf::TBot(hepmc,"t_pn");
  rad::rdf::TPrimeBot(hepmc,"tp_pn");
  rad::rdf::TPrimeTop(hepmc,"tp_gZ");

  //CM production angles
  rad::rdf::CMAngles(hepmc,"CM");

  //exlusivity
  rad::rdf::MissMass(hepmc,"MissMass","{scat_ele,n,Zc}");
  rad::rdf::MissP(hepmc,"MissP_Meson","{scat_ele,Zc}");
  rad::rdf::MissPt(hepmc,"MissPt_Meson","{scat_ele,Zc}");
  rad::rdf::MissPz(hepmc,"MissPz_Meson","{scat_ele,Zc}");
  rad::rdf::MissTheta(hepmc,"MissTheta_Meson","{scat_ele,Zc}");


  ///////////////////////////////////////////////////////////////
  //Define subsets of particles and corresponing variables to plot
  ///////////////////////////////////////////////////////////////
  hepmc.Define("electrons","rad::helpers::PositionsWhere(mc_pid==11)");
  hepmc.Define(MC()+"elsP",Form("Take(%spmag,electrons)",MC().data()));

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

  histo.Create<TH1D,double>({"Q2","Q2",500,0,2.},{"Q2"});
  histo.Create<TH1D,double>({"W","W",100,0,20.},{"W"});
  histo.Create<TH1D,double>({"MesonMass","M(e-,e+, #pi) [GeV]",100,.3,5.},{"ZMass"});
  histo.Create<TH1D,double>({"MissMass","Mmiss [GeV]",1000,-10,10},{"MissMass"});
  histo.Create<TH1D,double>({"JMass","M(e-,e+) [GeV]",100,.3,5.},{"JMass"});
  histo.Create<TH1D,double>({"tpn","t(p,n) [GeV^{2}]",100,-2,5},{"t_pn"});
  histo.Create<TH1D,double>({"tgZ","t(g,Z) [GeV^{2}]",100,-2,5},{"t_gZ"});
  histo.Create<TH1D,double>({"cthCM","cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  histo.Create<TH1D,double>({"missP","p_{miss}(e',Z)",105,0,105},{"MissP_Meson"});
  histo.Create<TH1D,double>({"missPt","p_{t,miss}(e',Z)",100,0,10},{"MissPt_Meson"});
  histo.Create<TH1D,double>({"missPz","p_{z,miss}(e',Z)",105,0,105},{"MissPz_Meson"});
  histo.Create<TH1D,double>({"missTheta","#theta_{miss}(e',Z)",100,0,1},{"MissTheta_Meson"});
  histo.Create<TH1D,ROOT::RVecD>({"allP","momentum of all particles",100,0,100},{"pmag"});
  histo.Create<TH1D,ROOT::RVecD>({"eleP","momentum of electrons",100,0,100},{"elsP"});

  gBenchmark->Start("processing");
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
 
  //Draw all meson masss histograms on 1 canvas
  histo.DrawAll("MesonMass");
   
  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");

  //save all histograms to file
  histo.File("HepMCZ_hists.root");

  gBenchmark->Stop("df");
  gBenchmark->Print("df");
  
  gBenchmark->Start("snapshot");
  //  hepmc.Snapshot("HepMCZ.root");
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");

}
