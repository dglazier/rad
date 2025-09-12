#include "BHxsec.C"

using ROOT::Math::PxPyPzMVector;

PxPyPzMVector GetVec(const std::vector<float>& phi,
		    const std::vector<float>& theta,
		    const std::vector<float>& eta,
		    const std::vector<float>& pmag,
		    int index,
		    double mass = 0.0
		    ){
  
  if (index < 0 || index >= phi.size()) {
    std::cerr << "Invalid particle index: " << index << std::endl;
    return ROOT::Math::PxPyPzMVector(0, 0, 0, mass);
  }
  
  double pt = pmag[index] * sin(theta[index]);
  double px = pt * cos(phi[index]);
  double py = pt * sin(phi[index]);
  double pz = pmag[index] * cos(theta[index]);

  return ROOT::Math::PxPyPzMVector(px, py, pz, mass);
  
}

void IntegrateBHEvents(){
  
  TFile *f = TFile::Open("tempout.root");
  TTree *tree = (TTree*)f->Get("rad_tree");

  std::vector<float> *mc_phi = nullptr;
  std::vector<float> *mc_theta = nullptr;
  std::vector<float> *mc_eta = nullptr;
  std::vector<float> *mc_pmag = nullptr;
  int beam_ele, beam_ion;
  int scat_ele, pprime, ele, pos;
  double heli_cos_theta;
  
  tree->SetBranchAddress("mc_phi", &mc_phi);
  tree->SetBranchAddress("mc_theta", &mc_theta);
  tree->SetBranchAddress("mc_eta", &mc_eta);
  tree->SetBranchAddress("mc_pmag", &mc_pmag);
  tree->SetBranchAddress("beam_ele", &beam_ele);
  tree->SetBranchAddress("beam_ion", &beam_ion);
  tree->SetBranchAddress("scat_ele", &scat_ele);
  tree->SetBranchAddress("pprime", &pprime);
  tree->SetBranchAddress("ele", &ele);
  tree->SetBranchAddress("pos", &pos);
  tree->SetBranchAddress("mc_Heli_CosTheta", &heli_cos_theta);
  
  const double Mp = 0.938; // Proton mass
  const double Me = 0.000511; // Electron mass
  
  const double GeV2_to_fb = 0.389379e9;
  double total_sigma = 0.0;
  Long64_t count = 0;
  
  double cos_theta_lim=0.999 ;
  
  TH1D *hDeltaT = new TH1D("hDeltaT","#Delta_T from p'-p vector",100,-1,-1);
  TH2D *h2d_DeltaT_t = new TH2D("h2d_DeltaT_t","#Delta_T vs |t|",100,-1,-1,100,-1,-1);

  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    if (i%100000==0) cout << Form("%lld / %lld",i,nentries) << endl;
    tree->GetEntry(i);
    
    //beams
    auto k = GetVec(*mc_phi, *mc_theta, *mc_eta, *mc_pmag, beam_ele, Me);
    auto p = GetVec(*mc_phi, *mc_theta, *mc_eta, *mc_pmag, beam_ion, Mp);
    
    //scattered e, recoil p
    auto kp = GetVec(*mc_phi, *mc_theta, *mc_eta, *mc_pmag, scat_ele, Me);
    auto pp = GetVec(*mc_phi, *mc_theta, *mc_eta, *mc_pmag, pprime, Mp);
    
    //virtual photon 1
    auto g = kp-k;
    
    //lepton pair
    auto eplus = GetVec(*mc_phi, *mc_theta, *mc_eta, *mc_pmag, ele, Me);
    auto eminus = GetVec(*mc_phi, *mc_theta, *mc_eta, *mc_pmag, pos, Me);
    
    //virtual photon 2, the "meson"
    auto gp = eplus+eminus;
    
    //cout << heli_cos_theta << endl;
    if(fabs(heli_cos_theta)>cos_theta_lim)continue;
    
    auto delta = pp-p;
    auto delta_T = sqrt(delta.Perp2());
    hDeltaT->Fill(delta_T);
    auto t = delta.M2();
    h2d_DeltaT_t->Fill(-t,delta_T);
    
    count++;
    total_sigma += 0.0;
    //total_sigma += dsigma_BH_vec(k,kp,p,pp,eplus,eminus); 
    //cout << "TOTAL COUNT: " << total_sigma << endl;
    
  }
  double Qp2_min = 2.0, Qp2_max = 10.0;
  double t_min = 0.0001, t_max = 2.0;
  double cos_theta_min = -cos_theta_lim, cos_theta_max = cos_theta_lim;
  double phi_min = 0.0, phi_max = 2 * TMath::Pi();
  
  double fb_to_nb = 1e-6;

  double phase_space_volume = (Qp2_max - Qp2_min) * (t_max - t_min) * (cos_theta_max - cos_theta_min) * (phi_max - phi_min);
  phase_space_volume=1;
  double sigma_fb = total_sigma * phase_space_volume * GeV2_to_fb  / count;
  double sigma_nb = sigma_fb * fb_to_nb;
  
  cout << sigma_fb << "fb" << endl;
  cout << sigma_nb << "nb" << endl;

  new TCanvas();
  hDeltaT->Draw();
  new TCanvas();
  h2d_DeltaT_t->Draw("colz");

}
