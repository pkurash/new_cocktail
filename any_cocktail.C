void any_cocktail(){
    
  //Analyze cocktail
  const Int_t N = 33;

  Double_t xbins[N+1] = { 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 
                          2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 
                          4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 
                          10, 11, 12, 13, 15, 20 };


  // Double_t xbins[N+1]={0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 7, 8, 10, 12, 14, 16, 18, 20};
    
    //TFile ff("../MergedResults/MergedCocktail_Full.root");

/*---------------------------------------------------------------------------*/

    TFile ff("AnalysisResults.root");
    
    TDirectoryFile *dd= ff.Get("GammaCocktailMC");
    TList *list=dd->Get("GammaCocktailMC_0.60");

    TDirectoryFile *dd2= ff.Get("HadronicCocktailMC");
    TList *list2=dd2->Get("HadronicCocktailMC_pi0_0.60");

    TH2F *Pt_Y_Pi0=list->FindObject("Pt_Y_Pi0");
    TH1D *Pt_Pi0=Pt_Y_Pi0->ProjectionX("Pt_Y_Pi0",151,250);
   
    TH2F *Pt_Y_Gamma=list->FindObject("Pt_Y_Gamma");
    TH1D *Pt_Gamma=Pt_Y_Gamma->ProjectionX("Pt_Gamma",188,213);

  // K-meson cascade sources
    TH2F *Pt_Y_Gamma_From_Pi0_From_K0s=list2->FindObject("Pt_Y_Gamma_From_Pi0_From_K0s");    
    TH2F *Pt_Y_Gamma_From_Pi0_From_K0l=list2->FindObject("Pt_Y_Gamma_From_Pi0_From_K0l"); 

    TH1D *Pt_Gamma_From_Pi0_From_K0s=Pt_Y_Gamma_From_Pi0_From_K0s->ProjectionX("Pt_Gamma_From_Pi0_From_K0s",188,213);
    TH1D *Pt_Gamma_From_Pi0_From_K0l=Pt_Y_Gamma_From_Pi0_From_K0l->ProjectionX("Pt_Gamma_From_Pi0_From_K0l",188,213);

    Pt_Gamma->Add(Pt_Gamma_From_Pi0_From_K0s,-1);
    Pt_Gamma->Add(Pt_Gamma_From_Pi0_From_K0l,-1);

 // gamma sources
  const Int_t nInputParticles = 16;
  TString fParticleListNames[nInputParticles] = {"Pi0", "Eta", "EtaPrim", "omega", "rho0", "rho+", "rho-", "phi", "J/psi", "Delta-", "Delta0", "Delta+", "Delta++", "Sigma0", "K0l", "K0s"};  
  
    TH2F  *fHistPtYGammaSource[nInputParticles];
    TH1D  *fHistPtGammaSource[nInputParticles];
    TH1D  *fHistPtGammaSource_rebinned[nInputParticles];
    TH1D  *fHist[nInputParticles];

   for(Int_t i=0;i<16;i++)
   {
   // name=Form("Pt_Y_Gamma_From_%s",fParticleListNames[i].Data());
    fHistPtYGammaSource[i]=(TH2F*)list->FindObject(Form("Pt_Y_Gamma_From_%s",fParticleListNames[i].Data()));
    fHistPtGammaSource[i]=(TH1D*)fHistPtYGammaSource[i]->ProjectionX(Form("%s",fParticleListNames[i].Data()),188,213);

       if(i==0)
       {
         fHistPtGammaSource[i]->Add(Pt_Gamma_From_Pi0_From_K0s, -1);
         fHistPtGammaSource[i]->Add(Pt_Gamma_From_Pi0_From_K0l, -1);
       }       

       fHistPtYGammaSource[i]->SetName(Form("%s",fParticleListNames[i].Data()));
       fHistPtYGammaSource[i]->SetTitle(Form("%s",fParticleListNames[i].Data()));
       fHistPtGammaSource[i]->Divide(fHistPtGammaSource[i],Pt_Gamma,1,1,"b");
   }
 

     Pt_Gamma_From_Pi0_From_K0s->Divide(Pt_Gamma_From_Pi0_From_K0s,Pt_Gamma,1,1,"b");
     Pt_Gamma_From_Pi0_From_K0l->Divide(Pt_Gamma_From_Pi0_From_K0l,Pt_Gamma,1,1,"b");

/*---------------------------------------------------------------------------*/
     
   fHistPtGammaSource[0]->GetXaxis()->SetTitle("p_{T},GeV/c");
   fHistPtGammaSource[0]->GetYaxis()->SetTitle("#gamma_{x}/#gamma_{bckgr}");
   fHistPtGammaSource[0]->GetYaxis()->SetRangeUser(1e-6,1.8);


    TH1D *hpg=Pt_Gamma->Rebin(N,"hpg",xbins);
    TH1D *hpi0=Pt_Pi0->Rebin(N,"hpi0",xbins);

    hpg->Divide(hpg,hpi0,1,1,"b");
    hpg->Scale(100./26);
    
    TCanvas *cc_g2g=new TCanvas("cc_g2g","cc_g2g",900,600);  
    cc_g2g->SetLogy();
 
    //fHistPtGammaSource[0]->Draw();
    for(Int_t i = 0; i < nInputParticles; i++)
    {
       fHistPtGammaSource[i]->Draw(Form("%s", i == 0 ? "" : "same"));
    }

     Pt_Gamma_From_Pi0_From_K0s->Draw("same");
     Pt_Gamma_From_Pi0_From_K0l->Draw("same");

/*---------------------------------------------------------------------------*/


    TH1D *hgpi0     = fHistPtGammaSource[0]->Rebin(N, "hgpi0",     xbins);
    TH1D *hgeta     = fHistPtGammaSource[1]->Rebin(N, "hgeta",     xbins);
    TH1D *hgetaprim = fHistPtGammaSource[2]->Rebin(N, "hgetaprim", xbins);
    TH1D *hgomega   = fHistPtGammaSource[3]->Rebin(N, "hgomega",   xbins);

    TH1D *hgk0s = Pt_Gamma_From_Pi0_From_K0s->Rebin(N, "hgk0s", xbins);
    TH1D *hgk0l = Pt_Gamma_From_Pi0_From_K0l->Rebin(N, "hgk0l", xbins);

    Double_t xw;

    for(Int_t i = 0; i < N; i++)
    {
      //xw = 10*(xbins[i+1] - xbins[i]);
      xw = (xbins[i+1] - xbins[i])/0.1;

      hgpi0    ->SetBinContent(i+1, (hgpi0    ->GetBinContent(i+1))/xw);
      hgeta    ->SetBinContent(i+1, hgeta    ->GetBinContent(i+1)/xw);
      hgetaprim->SetBinContent(i+1, hgetaprim->GetBinContent(i+1)/xw);
      hgomega  ->SetBinContent(i+1, hgomega  ->GetBinContent(i+1)/xw);

      hgk0s->SetBinContent(i+1, hgk0s->GetBinContent(i+1)/xw);
      hgk0l->SetBinContent(i+1, hgk0l->GetBinContent(i+1)/xw);

      /*----------------*/

      hgpi0    ->SetBinError(i+1, hgpi0    ->GetBinError(i+1)/xw);
      hgeta    ->SetBinError(i+1, hgeta    ->GetBinError(i+1)/xw);
      hgetaprim->SetBinError(i+1, hgetaprim->GetBinError(i+1)/xw);
      hgomega  ->SetBinError(i+1, hgomega  ->GetBinError(i+1)/xw);

      hgk0s->SetBinError(i+1, hgk0s->GetBinError(i+1)/xw);
      hgk0l->SetBinError(i+1, hgk0l->GetBinError(i+1)/xw);

    }



    //hgpi0->Add(hgk0s, -0.5);
    //hgeta->Add(hgk0l, -1);


    hgpi0->SetLineWidth(2);
    hgpi0->SetLineColor(kRed);
    hgpi0->SetMarkerColor(hgpi0->GetLineColor());
    hgpi0->SetAxisRange(xbins[0], xbins[N-1], "X");
    hgpi0->SetAxisRange(3e-3, 1.3, "Y");
    hgpi0->SetTitle("#pi^{0}");

    hgeta->SetLineWidth(2);
    hgeta->SetLineColor(kBlue);
    hgeta->SetMarkerColor(hgeta->GetLineColor());
    hgeta->SetTitle("#eta");

    hgetaprim->SetLineWidth(2);
    hgetaprim->SetLineColor(kViolet);
    hgetaprim->SetMarkerColor(hgetaprim->GetLineColor());
    hgetaprim->SetTitle("#eta'");

    hgomega->SetLineWidth(2);
    hgomega->SetLineColor(kMagenta);
    hgomega->SetMarkerColor(hgomega->GetLineColor());
    hgomega->SetTitle("#omega");

    hgk0s->SetLineWidth(2);  
    hgk0s->SetLineColor(kOrange);
    hgk0s->SetMarkerColor(hgk0s->GetLineColor());
    hgk0s->SetTitle("K^{0}_{S}");

    hgk0l->SetLineWidth(2);  
    hgk0l->SetLineColor(kGreen);
    hgk0l->SetMarkerColor(hgk0l->GetLineColor());
    hgk0l->SetTitle("K^{0}_{L}");


    TCanvas *cc=new TCanvas("cc","cc",900,600); 

    cc->SetLogy();
   // cc->SetLogx();
    cc->SetGridx();
    cc->SetGridy();

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    hgpi0->Draw();
    hgeta->Draw("same");
    hgetaprim->Draw("same");
    hgomega->Draw("same");
    hgk0s->Draw("same");
    hgk0l->Draw("same");
 
    cc->BuildLegend(0.8, 0.5, 0.9, 0.9);

    cc->Print("figures/cc_cocktails.pdf");

/*---------------------------------------------------------------------------*/
    
    TFile fout("ROOT/Cocktails_13TeV.root","recreate");

     for(Int_t i=  0; i < nInputParticles; i++)
     {
         fHistPtGammaSource[i]->Write();
     }
     Pt_Gamma_From_Pi0_From_K0s->Write();
     Pt_Gamma_From_Pi0_From_K0l->Write();

     Pt_Gamma->Write();
   
     hpg->Write();

     cc_g2g->Write();
     cc->Write();
    hgpi0->Write();
    hgeta->Write();
    hgetaprim->Write();
    hgomega->Write();
    hgk0s->Write();
    hgk0l->Write();

    fout.Close();


   
} 
