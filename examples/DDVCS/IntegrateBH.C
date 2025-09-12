#include "BHxsec.C"
#include "Constants.h"

using ROOT::Math::PxPyPzMVector;
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::RotateY;

TH1D *hQ2;
TH1D *hy;
TH1D *hxbj;
TH1D *ht;

TH2D *h2d_DeltaT_t = new TH2D("h2d_DeltaT_t","#Delta_T vs |t|",100,0,1,100,0,1);


Long64_t total_count = 0;
Long64_t good_count = 0;

double MakeEvent(double Ebeam, double Pbeam, double Q2, double y, double Qp2, double t, double phi, double costhetaL, double phiL){
  total_count++;
  //minimum value of spacelike virtuality that can be evaluated for electron mass me
  auto Q2min = y*y*me*me/(1-y);
  //cout << Q2min << endl;
  if(Q2<std::pow(y*me,2)/(1.-y)){
    //cout << "Q2 lower than allowed spacelike virtuality" << endl;
    return 0.0;
  }
  // const double Ebeam=10.6;
  // const double Pbeam=M;
  double pmag_ebeam = sqrt(Ebeam*Ebeam - me*me);
  double pmag_pbeam = sqrt(Pbeam*Pbeam - M*M);
  
  PxPyPzMVector ebeam(0,0,-pmag_ebeam,me);
  PxPyPzMVector pbeam(0,0,pmag_pbeam,M);
  
  const double s_electro = (ebeam+pbeam).M2(); 
 
  auto boost_lab_to_tar = pbeam.BoostToCM();
  auto pbeam_tar = boost(pbeam,boost_lab_to_tar);
  auto ebeam_tar = boost(ebeam,boost_lab_to_tar);
  auto Ebeam_tar = ebeam_tar.E();
  double Escat_tar = Ebeam_tar * (1. - y);
  double pr_costheta_e = 1. - Q2 / (2*Ebeam_tar*Escat_tar);
    
  if (std::abs(pr_costheta_e) > 1.){
    //cout << "Theta too large " << endl;
    return 0.0;
  }
  double pr_theta_e = std::acos(pr_costheta_e);
  auto ebeam_tar_unit = ebeam_tar.Vect().Unit();
  auto ebeam_tar_rot = RotateY(ebeam_tar_unit, pr_theta_e);
  PxPyPzMVector eprime_tar(ebeam_tar_rot.X()*Escat_tar, ebeam_tar_rot.Y()*Escat_tar, ebeam_tar_rot.Z()*Escat_tar, me);
  auto q_tar = ebeam_tar - eprime_tar;
    
  //"average photon momentum"
  //qbar = 0.5(q + qp);
  //"average nucleon momentum" 
  //P = 0.5*(pbeam + pprime)
    
  auto qbar2 = -0.5 * (Q2 - Qp2 + t/2);
  double s_photo = (q_tar + pbeam_tar).M2();
  if(s_photo<M*M){
    //cout << "s below threshold for space-like/quasireal photon" << endl;
    return 0.0;
  }
  if(s_photo>s_electro){
    //cout << "s larger than maximum centre mass energy of reaction" << endl;
    return 0.0;
  }
    
  //back to lab frame
  auto eprime = boost(eprime_tar,-boost_lab_to_tar);
  PxPyPzMVector q = ebeam - eprime;
    
  if (std::pow(s_photo - Qp2 - M*M,2) - 4*Qp2*M*M < 0.0){
    //cout << "Invalid Kinematics for s_photo" << endl;
    return 0.0;
  }
    
  auto m1_2 = -Q2;
  auto m2_2 = M*M;
  auto m3_2 = Qp2;
  auto m4_2 = M*M;
  auto s4 = s_photo;
  auto E1cm = (s4 + m1_2 - m2_2) / (2 * sqrt(s4));
  auto E3cm = (s4 + m3_2 - m4_2) / (2 * sqrt(s4));

  auto p1cm = sqrt(std::pow(E1cm, 2) - m1_2);
  auto p3cm = sqrt(std::pow(E3cm, 2) - m3_2);

  auto tmin = std::pow((m1_2 - m3_2 - m2_2 + m4_2) / (2 * sqrt(s4)), 2)
    - std::pow(p1cm - p3cm, 2);
  auto tmax = std::pow((m1_2 - m3_2 - m2_2 + m4_2) / (2 * sqrt(s4)), 2)
    - std::pow(p1cm + p3cm, 2);
  if (t > tmin || t < tmax) {
    // cout << "Invalid kinematics for t = " << t << endl;
    // cout << "tmin: " << tmin << " tmax: " << tmax << endl;
    return 0.0;
  }
    
    
  auto boost_tar_to_gp = (q_tar+pbeam_tar).BoostToCM();
  auto gamma_gp = boost(q_tar,boost_tar_to_gp);
  auto pbeam_gp = boost(pbeam_tar,boost_tar_to_gp);
    
  double pmag_qp_gp = sqrt( (std::pow(s4 - Qp2 - M*M, 2) - 4 * Qp2 * M*M) / (4 * s4) );
  auto dir_gp = gamma_gp.Vect().Unit();
    
  PxPyPzMVector qp_gp(pmag_qp_gp*dir_gp.X(),pmag_qp_gp*dir_gp.Y(),pmag_qp_gp*dir_gp.Z(),sqrt(Qp2));
    
  PxPyPzMVector pprime_gp(-qp_gp.X(),-qp_gp.Y(),-qp_gp.Z(),M);
    
  auto p_E_gp = pbeam_gp.E();
  auto p_pmag_gp = pbeam_gp.P();
  auto pprime_E_gp = pprime_gp.E();
  auto pprime_pmag_gp = pprime_gp.P();
    
    
  double costheta_gp = (t - std::pow(p_E_gp - pprime_E_gp, 2) + std::pow(p_pmag_gp, 2) + std::pow(pprime_pmag_gp, 2)) / (2 * p_pmag_gp * pprime_pmag_gp);
    
  ROOT::Math::RotationY gp_rot(std::acos(costheta_gp));
  auto pprime_gp_rot = gp_rot(pprime_gp);
    
  auto delt_gp = pprime_gp_rot - pbeam_gp;
  auto delt_T = sqrt(delt_gp.Perp2());
  auto t_rot_gp = delt_gp.M2();
  auto r = sqrt(pow(s4 - Qp2 - M*M, 2) - 4*Qp2*M*M);
  auto theta_gp = std::acos(costheta_gp);
  auto delt_T_calc = sin(theta_gp) * r / (2*sqrt(s4));
  if(TMath::IsNaN(delt_T) || TMath::IsNaN(delt_T_calc)){
    // cout << "delt_T NAN" << endl;
    // cout << delt_T << " " << delt_T_calc <<endl;
    // cout << theta_gp <<" " << sin(theta_gp) << " " << r << " " << s4 << endl;
    // cout << costheta_gp << endl;
    return 0.0;
  }
  auto qp_tar = boost(qp_gp,-boost_tar_to_gp);
  auto pprime_tar = boost(pprime_gp_rot,-boost_tar_to_gp);
    
  auto t_tar = (pprime_tar - pbeam_tar).M2();
    
  ROOT::Math::AxisAngle aa(q_tar.Vect().Unit(), phi);
  ROOT::Math::Rotation3D phi_rot(aa);
  auto pprime_tar_rot = phi_rot(pprime_tar);
  auto qp_tar_rot = phi_rot(qp_tar);
    
  auto t_rot_tar = (pprime_tar_rot - pbeam_tar).M2();
    
  //cout << "t_gen: " << t << " t_rot_gp: " << t_rot_gp  << " t_tar" << t_tar << " t_rot_tar: " << t_rot_tar << endl;
    
  auto boost_tar_to_exc = qp_tar.BoostToCM();
  auto pbeam_exc = boost(pbeam_tar,boost_tar_to_exc);
  auto pprime_exc = boost(pprime_tar,boost_tar_to_exc);
    
    
  if (sqrt(Qp2)<2*me){
    //cout << "Invalid Kinematics: Q'2 too small for decay choice!" << endl;
    return 0.0;
  }
    
  auto pprime_exc_dir = pprime_exc.Vect().Unit();
  auto lep_pair_3mom = sqrt(Qp2)/2.0  *  pprime_exc_dir;
    
  PxPyPzMVector ele(lep_pair_3mom.X(),lep_pair_3mom.Y(),lep_pair_3mom.Z(),me);
  PxPyPzMVector pos(-lep_pair_3mom.X(),-lep_pair_3mom.Y(),-lep_pair_3mom.Z(),me);
    
  auto rotate_theta = pbeam_exc.Vect().Cross(pprime_exc.Vect().Unit());
  auto thetaL = std::acos(costhetaL);
    
  ROOT::Math::AxisAngle aa_theta_cm(rotate_theta, thetaL);
  ROOT::Math::Rotation3D theta_cm_rot(aa_theta_cm);
  ele = theta_cm_rot(ele);
  pos = theta_cm_rot(pos);
  ROOT::Math::AxisAngle aa_phi_cm(pprime_exc.Vect().Unit(), phiL);
  ROOT::Math::Rotation3D phi_cm_rot(aa_phi_cm);
  ele = phi_cm_rot(ele);
  pos = phi_cm_rot(pos);
    
  auto ele_tar = boost(ele,-boost_tar_to_exc);
  auto pos_tar = boost(pos,-boost_tar_to_exc);
    
  //phiS dependence in case of trans pol target, later, not needed now

  // back to lab
  auto q_lab = boost(q_tar,-boost_lab_to_tar);
  auto qp_lab = boost(qp_tar,-boost_lab_to_tar);
  auto eprime_lab = boost(eprime_tar,-boost_lab_to_tar);
  auto pprime_lab = boost(pprime_tar,-boost_lab_to_tar);
  auto ele_lab = boost(ele_tar,-boost_lab_to_tar);
  auto pos_lab = boost(pos_tar,-boost_lab_to_tar);
    
  auto xbj = Q2/(2*pbeam.Dot(q_lab));
  auto eps = 2*xbj*M/sqrt(Q2);
  //cout << "y gen vs y eps test " << y << " " << sqrt(Q2)/(eps*ebeam_tar.E()) << endl;
  auto nu = Q2/(2*M*xbj);
  //cout << "nu calc scalar vs vector TRF " << nu << " " << pbeam_tar.Dot(q_tar)/M <<  " " << ebeam_tar.E() - eprime_tar.E() << endl;
  auto photo_flux = (alpha/(2.0*TMath::Pi()*Q2)) * (1.0 + (1.0-y)*(1.0-y)/y - 2.0*(1.0-y)*Q2min/(y*Q2) ) * nu/(ebeam_tar.E()*xbj);
  //photo_flux = (alpha/(4.0*TMath::Pi()*M*ebeam_tar.E())) * (1/(y*y) -1 - (2*me*me/Q2) );
  //cout << "nu: " << nu << " dGamma/dQ2dxb: " << photo_flux << endl;
  if(photo_flux<=0){
    //cout << "Unphysical photoflux value" << endl;
    return 0.0;
  }
  auto jacob_dx_dy = xbj / y;
  photo_flux = photo_flux * jacob_dx_dy;
    
  double d4sigma = dsigma_BH(s4,Qp2, t, thetaL, phiL, delt_T);
  if (d4sigma<=0){
    //cout << "Unphysical d4sigma value" << endl;
    return 0.0;
  }
  hQ2->Fill(Q2);
  hy->Fill(y);
  hxbj->Fill(xbj);
  ht->Fill(-t);
  h2d_DeltaT_t->Fill(-t,delt_T);
  
  //cout << "Q2, y , Qp2, t, phi, costhetaL, phiL, d7sigma" << endl;
  //cout << Q2 << " " << y << " " << Qp2 << " " << t << " " << phi << " " << costhetaL << " " << phiL << " " << d4sigma*photo_flux << endl;
  good_count++;
  return d4sigma*photo_flux;  
  //return 1.0;
  //return d4sigma;
}

struct kinsetup{
  double Ebeam, Pbeam;
  double y_min, y_max;
  double Q2_min, Q2_max;
  double Qp2_min, Qp2_max;
  double t_min, t_max;
  double phi_min, phi_max;
  double phiS_min, phiS_max;
  double phiL_min, phiL_max;
  double thetaL_min, thetaL_max;
  
  // Constructor
  kinsetup(double Ebeam, double Pbeam,
	   double y_min, double y_max,
	   double Q2_min, double Q2_max,
	   double Qp2_min, double Qp2_max,
	   double t_min, double t_max,
	   double phi_min, double phi_max,
	   double phiS_min, double phiS_max,
	   double phiL_min, double phiL_max,
	   double thetaL_min, double thetaL_max)
    : Ebeam(Ebeam), Pbeam(Pbeam), 
      y_min(y_min), y_max(y_max),
      Q2_min(Q2_min), Q2_max(Q2_max),
      Qp2_min(Qp2_min), Qp2_max(Qp2_max),
      t_min(t_min), t_max(t_max),
    phi_min(phi_min), phi_max(phi_max),
    phiS_min(phiS_min), phiS_max(phiS_max),
    phiL_min(phiL_min), phiL_max(phiL_max),
    thetaL_min(thetaL_min), thetaL_max(thetaL_max) {}
};

kinsetup jlab12{10.6, M, 0.0, 1.0, 0.15, 5.0, 2.25, 9.0, -0.8, -0.1, 0.0, 2*TMath::Pi(), 0.0, 2*TMath::Pi(), 0.1, 2*TMath::Pi()-0.1, TMath::Pi()/4.0, TMath::Pi()*3.0/4.0};
kinsetup jlab22{22, M, 0.0, 1.0, 0.15, 5.0, 2.25, 9.0, -0.8, -0.1, 0.0, 2*TMath::Pi(), 0.0, 2*TMath::Pi(), 0.1, 2*TMath::Pi()-0.1, TMath::Pi()/4.0, TMath::Pi()*3.0/4.0};
kinsetup eic5x41{5, 41, 0.0, 1.0, 0.15, 5.0, 2.25, 9.0, -1.0, -0.05, 0.0, 2*TMath::Pi(), 0.0, 2*TMath::Pi(), 0.1, 2*TMath::Pi()-0.1, TMath::Pi()/4.0, TMath::Pi()*3.0/4.0};
kinsetup eic10x100{10, 100, 0.0, 1.0, 0.15, 5.0, 2.25, 9.0, -1.0, -0.05, 0.0, 2*TMath::Pi(), 0.0, 2*TMath::Pi(), 0.1, 2*TMath::Pi()-0.1, TMath::Pi()/4.0, TMath::Pi()*3.0/4.0};

kinsetup MYeic18x275{18, 275, 0.05, 1.0, 0.0, 5.0, 2.0, 10.0, -2.0, -0.0001, 0.0, 2*TMath::Pi(), 0.0, 2*TMath::Pi(), 0.0, 2*TMath::Pi(), 0, TMath::Pi()};

void IntegrateBH(kinsetup kin=jlab12){
  
  double y_min = kin.y_min, y_max = kin.y_max;
  double Q2_min = kin.Q2_min, Q2_max = kin.Q2_max;
  double Qp2_min = kin.Qp2_min, Qp2_max = kin.Qp2_max;
  double t_min = kin.t_min, t_max = kin.t_max;
  //double xb_min = 0.000001, xb_max = 1.0;
  
  double phi_min = kin.phi_min, phi_max = kin.phi_max;
  double phiS_min = kin.phiS_min, phiS_max = kin.phiS_max;
  double phiL_min = kin.phiL_min, phiL_max = kin.phiL_max;
  double thetaL_min = kin.thetaL_min, thetaL_max = kin.thetaL_max;
  
  double costhetaL_max = std::cos(thetaL_min);
  double costhetaL_min = std::cos(thetaL_max);
  
  hQ2 = new TH1D("hQ2","Q^{2} Generated By Integrator; Q^{2} [GeV^{2}]",100,0,Q2_max);
  hy = new TH1D("hy","y Generated By Integrator; y",100,0,1);
  hxbj = new TH1D("hxbj","xbj calculated; x_{bj}",100,0,1);
  ht = new TH1D("ht","-t Generated By Integrator;|t| [GeV^{2}]",100,0,2);

  //double MakeEvent(double Q2, double y, double Qp2, double t, double phi, double costhetaL, double phiL)
  auto f = [E = kin.Ebeam, P = kin.Pbeam](const double *args) {return MakeEvent(E, P, args[0], args[1], args[2], args[3], args[4], args[5], args[6]);};
  ROOT::Math::Functor functor(f,7);

  auto flatOne = [](const double *x) { return 1.0; };
  ROOT::Math::Functor flatfunc(flatOne,7);
  
  ROOT::Math::GSLMCIntegrator mc(ROOT::Math::IntegrationMultiDim::kPLAIN,1e-12,1e-12,1e6);
  //ROOT::Math::AdaptiveIntegratorMultiDim mc(0.0,1e-9,1e6,0);
  //ROOT::Math::IntegratorMultiDim mc(ROOT::Math::IntegrationMultiDim::kVEGAS,-1,-1,1e6);
  mc.SetFunction(functor); 
  //mc.SetFunction(flatfunc);
  
  double arg_min[7] = {Q2_min, y_min, Qp2_min, t_min, phi_min, costhetaL_min, phiL_min};
  double arg_max[7] = {Q2_max, y_max, Qp2_max, t_max, phi_max, costhetaL_max, phiL_max};
  double result = mc.Integral(arg_min, arg_max);
  double error = mc.Error();
  const double GEV2_TO_PB = 38937.93719; // pb per GeV^-2
  double frac = (double) good_count / total_count;
  std::cout << "good calls, total calls, r: " << good_count << " " << total_count << " " << frac << std::endl;
  double V = (Q2_max-Q2_min)*(y_max-y_min)*(Qp2_max-Qp2_min)*(t_max-t_min)*(phi_max-phi_min)*(costhetaL_max-costhetaL_min)*(phiL_max-phiL_min);
  std::cout << "VOLUME: " << V << endl;
  std::cout << "GSLMC result [GeV^-2]: " << result << " +- " << error << std::endl;
  std::cout << "GSLMC result [pb]: " << result*GEV2_TO_PB << " +- " << error*GEV2_TO_PB << std::endl;
  std::cout << "GSLMC result*volume [pb]: " << V*result*GEV2_TO_PB << " +- " << V*error*GEV2_TO_PB << std::endl;
  
  TCanvas *c00 = new TCanvas();
  c00->Divide(2,2);
  c00->cd(1);
  hQ2->SetMinimum(0);
  hQ2->Draw("hist");
  c00->cd(2);
  hy->SetMinimum(0);
  hy->Draw("hist");
  c00->cd(3);
  hxbj->SetMinimum(0);
  hxbj->Draw("hist");
  c00->cd(4);
  ht->SetMinimum(0);
  ht->Draw("hist");
  
  new TCanvas();
  h2d_DeltaT_t->Draw("colz");

}
