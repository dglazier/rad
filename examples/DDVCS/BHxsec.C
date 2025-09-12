#include "Constants.h"

using ROOT::Math::PxPyPzMVector;

const double alpha = rad::constant::alpha();
const double M = rad::constant::M_pro();
const double me = rad::constant::M_ele(); 

double dsigma_BH_vec(const PxPyPzMVector &ebeam, const PxPyPzMVector &escat, const PxPyPzMVector &p, const PxPyPzMVector &pp,  const PxPyPzMVector &k,  const PxPyPzMVector &kp) {
  const double s = (ebeam+p).M2(); 
  
    auto q = ebeam-escat;
    //auto Q2 = -q.M2();
    auto delt = pp-p;
    auto delt_T = sqrt(delt.Perp2());
    auto t = delt.M2();
    auto qp = k+kp;
    auto Qp2 = qp.M2();
    auto tau = Qp2 / (2*p.Dot(q));
    //auto xb = Q2 / (2*p.Dot(q));

    // Scalar products (simplified placeholders)
    double a = 2.0 * (k-kp).Dot(pp); // placeholder for 2(k - k')·p'
    double b = 2.0 * (k-kp).Dot(p-pp);       // placeholder for 2(k - k')·(p - p')
      
  //return dsigma_BH(s,Qp2,t,theta,phi);
  return 0; //fix later
}

//theta,phi are l+l- cm frame, i.e. "helicity/decay" angles
//Theta is scattering angle in gamma-p c.m. frame.
//sinTheta = 2*delt_T*sqrt(s) / r, according to paper.
//but delt_T is transverse part of delt vector.
double dsigma_BH(const double s, const double Qp2, const double t, const double theta, const double phi, double delta_T){
  
  // Form factors (dipole approximation)
  double Lambda2 = 0.71;
  double GD = 1.0 / pow(1.0 - t / Lambda2, 2);
  double F1 = GD;
  double F2 = (GD - F1) / (t / (4 * M * M));
  
  //kinematics
  double beta = sqrt(1.0 - 4.0*me*me / Qp2);
  double r = sqrt((s-Qp2-M*M)*(s-Qp2-M*M)  - 4*Qp2*M*M);
  double a = beta*r*TMath::Cos(theta);
  
  double tau = Qp2/(s-M*M);
  double delt_T2 = -t*(1-tau) - tau*tau*M*M;
  double delt_T = sqrt(delt_T2);
  double sigmaterm = (Qp2*(s-M*M-Qp2) + t*(s-M*M+Qp2))/r;
  double bterm1 = 2*(s-M*M)*sqrt(Qp2)*delt_T/r;
  double bterm2 = sigmaterm * beta * TMath::Cos(theta);
  double bterm3 = -beta * bterm1 * TMath::Sin(theta) * TMath::Cos(phi);
  double b = bterm2 + bterm3;
  
  //lepton propagator
  double L = ((Qp2-t)*(Qp2-t) - b*b)/4;
  double L0 = Qp2*Qp2*sin(theta)*sin(theta)/4;
  
  // A and B terms (simplified)
  double A = pow(s-M*M, 2) * delt_T*delt_T - t*a*(a + b) - M*M*b*b - t*(4*M*M-t)*Qp2 + (me*me/L)*( pow((Qp2-t)*(a+b)-(s-M*M)*b,2) + t*(4*M*M-t)*pow(Qp2-t,2) );
  double B = pow(Qp2+t, 2) + b*b + 8*me*me*Qp2 - 4*me*me*(t+2*me*me)*pow(Qp2-t, 2)/L;

  // Differential cross section
  double prefactor = alpha * alpha * alpha / (4.0 * TMath::Pi() * pow(s - M * M, 2));
  double prefactor2 = beta/(-t*L);
  double sigma = prefactor * prefactor2 * ( (F1*F1 - (t/(4*M*M))*F2*F2)*A/-t + pow(F1 + F2, 2)*B/2.0 );
  
  if (sigma<0 || TMath::IsNaN(sigma)){
    // cout << "Bad xsec, negatives for pos definite terms!!" << endl;
    // cout << "tau " << tau << endl;
    // cout << "delt_T2 " << delt_T2 << endl;
    // cout << "delta_T_full: " << delta_T << " delta_T_approx: " << delt_T << endl;
    // cout << "Qp2, s, t, r, sigmaterm " << Qp2 << " " << " " << s << " " << t <<" "<< r << " " << sigmaterm << endl;
    // cout << "bterm1, bterm2, bterm3, b" << endl;
    // cout << bterm1 << " " << bterm2 << " " << bterm3 << " " << b << endl;
    // cout << "sintheta, costheta, cosphi, phi" << sin(theta) <<" " << cos(theta) << " " << cos(phi) << " " << phi << endl;
    // cout  <<"beta: " << beta << " t: " << t << " L: " << L <<  " L0: " << L0 << endl;
    // cout << prefactor << " " << prefactor2 <<" " << A << " " << B << " " << sigma << endl;
    // cout << endl;
  }
  return sigma; // Units: GeV^-4
}
