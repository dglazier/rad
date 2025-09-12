void DrawCrossSection(const char *fname = "tempout.root") {
  
  // Open file and get tree
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) {
        Error("make_dsigma_plots", "Cannot open file %s", fname);
        return;
    }
    TTree *rad_tree = (TTree*)f->Get("rad_tree");
    if (!rad_tree) {
        Error("make_dsigma_plots", "Tree rad_tree not found in file");
        return;
    }
    
    // ==== User settings ====
    double sigma_nb   = 0.45;   // total integrated cross section from HepMC header
    double Ngen       = 1e6;    // total number of generated events
    double Qprime2_target = 5.0; // GeV^2
    // =======================

    // Convert sigma to fb
    double sigma_fb = sigma_nb * 1e6;
    double sigma_pb = sigma_nb * 1e3;
    
    // Per-event weight
    double w_event_1  = sigma_fb / Ngen;
    double w_event_2 = sigma_pb / Ngen;
    
    
    // ------------------------------------------------
    // 1. ds/d|t|  (integrated over everything else)
    // ------------------------------------------------
    rad_tree->Draw("abs(mc_t_bot)>>h_dsigma_dt(100,0,1)", Form("%g", w_event_1), "goff");
    TH1D *h_dsigma_dt = (TH1D*)gDirectory->Get("h_dsigma_dt");

    // Divide by bin width to get per GeV^2
    h_dsigma_dt->Scale(1.0 / h_dsigma_dt->GetBinWidth(1));
    h_dsigma_dt->SetTitle("d#sigma/d|t| (integrated over Q'^2 etc.)");
    h_dsigma_dt->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    h_dsigma_dt->GetYaxis()->SetTitle("fb / GeV^{2}");

    // ------------------------------------------------
    // 2. d²s/dQ'^2 d|t| at Q'^2 = target
    // ------------------------------------------------
    rad_tree->Draw("abs(mc_t_bot):mc_GMass*mc_GMass>>h2(100,3.0,7.0,100,0,0.4)",
                   Form("%g  && mc_Heli_Theta>3.14/4.0 && mc_Heli_Theta<3*3.14/4.0", w_event_2), "goff");
    TH2D *h2 = (TH2D*)gDirectory->Get("h2");

    // Scale to per GeV^4
    double dQ = h2->GetXaxis()->GetBinWidth(1);
    double dt = h2->GetYaxis()->GetBinWidth(1);
    h2->Scale(1.0 / (dQ * dt));

    // Project at Q'^2 = target
    int binQ = h2->GetXaxis()->FindBin(Qprime2_target);
    TH1D *h_t_atQ = h2->ProjectionY("h_t_atQ", binQ, binQ);
    h_t_atQ->SetTitle(Form("d^{2}#sigma/dQ'^{2}d|t| at Q'^{2} = %.1f GeV^{2}", Qprime2_target));
    h_t_atQ->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    h_t_atQ->GetYaxis()->SetTitle("pb / GeV^{4}");

    // ------------------------------------------------
    // Draw both
    // ------------------------------------------------
    TCanvas *c1 = new TCanvas("c1","Cross section plots",1200,500);
    c1->Divide(2,1);

    c1->cd(1)->SetLogy();
    h_dsigma_dt->Draw("hist");

    c1->cd(2)->SetLogy();
    //h_t_atQ->SetMaximum(1000);
    //h_t_atQ->SetMinimum(0.1);
    h_t_atQ->Draw("l");

    c1->SaveAs("dsigma_plots.pdf");
    c1->Close();
    // Cleanup
    f->Close();
}
