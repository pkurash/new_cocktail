/****************************************************************************************************************************
******    Friederike Bock, friederike.bock@cern.ch                                                                      *****
******    Mike Sas, mike.sas@cern.ch                                                                                    *****
******    Lucas Altenkaemper, lucas.altenkaemper@cern.ch                                                                *****
*****************************************************************************************************************************/
#include "CocktailFunctions.h"
#include "CocktailPlotting.h"
#include "CocktailHEPDataPPb.h"
#include <iostream>

void  CocktailProduceCompleteInputFilePPb(  TString enableCentralities5TeV  = "11111111111111111",              //"000001110",
                                            TString enableCentralities8TeV  = "0000000000000000",
                                            TString enableParticles         = "1100011101100001100110111111",
                                            TString suffix                  = "eps",
                                            Bool_t produceAllSpectraPlots   = kTRUE,
                                            Bool_t enablePtY                = kTRUE,
                                            Bool_t thesisPlots              = kFALSE
                                            //TString convertTo               = ""                      // "dN/dydpT" or "1/2pipT*dN/dydpT"
)
{

    if (thesisPlots) isThesis                                   = kTRUE;
    //================================================================================================================
    // Enable energies and centralities to be included
    //================================================================================================================

    // This is done by changing the "enableCentralitiesXTeV" string
    // "enableCentralitiesXTeV" has 9 digits, set each digit either 0(exclude) or 1(include)
    // the centralities are sorted as follows:
    //      MB, 0-5, 5-10, 0-10, 10-20, 0-20, 20-40, 20-50, 40-60, 60-80, 80-100
    // to exclude an energy, just disable all centralities
    // Example: enableCentralities5TeV="000011110" enables pPb 5TeV and the centralities 10-20, 0-20, 20-40, 20-50
    Bool_t enable5TeV                                           = kTRUE;

    Bool_t enable8TeV                                           = kTRUE;

    Bool_t includeCentrality5TeV[nCentralities]                 = {kFALSE};
    Bool_t includeCentrality8TeV[nCentralities]                 = {kFALSE};

    Int_t tempCounter                                           = 0;
    TString tempString                                          = "";
    for (Int_t i=0; i<nCentralities; i++) {
        TString tempString                                      = enableCentralities5TeV(i,1);

        if (tempString.CompareTo("1") == 0)
            includeCentrality5TeV[i]                            = kTRUE;
        else
            tempCounter++;
    }
    if (tempCounter == nCentralities)
        enable5TeV                                              = kFALSE;

    tempCounter                                                 = 0;
    tempString                                                  = "";
    for (Int_t i=0; i<nCentralities; i++) {
        tempString                                              = enableCentralities8TeV(i,1);

        if (tempString.CompareTo("1") == 0)
            includeCentrality8TeV[i]                            = kTRUE;
        else
            tempCounter++;
    }
    if (tempCounter == nCentralities)
        enable8TeV                                              = kFALSE;

    if (!enable5TeV && !enable8TeV) {
        cout << "Warning: No energies and centralities selected, stopping!" << endl;
        return;
    }


    //================================================================================================================
    // Enable particles to be included
    //================================================================================================================

    // This is done by changing the "enableParticles" string
    // "enableParticles" has 22 digits, set each digit either 0(exclude) or 1(include)
    // the particles are sorted as follows:
    //     pi0, eta, omega, eta', gamma_dir, pi^+/-, K^+/-, p/anti-p, h^+/-, phi, K^*0, rho^0, rho^+/-, Delta^0, Delta^+/-, K^0_s, Lambda, Sigma^0, Sigma^+/-, Omega^+/-, Xi^+/-, J/psi, D0, D+, D*+
    // Example: enableParticles="0000011100000000000000000" enables pi^+/-, K^+/- and p/anti-p

    Bool_t includeParticle[nParticles]                          = {kFALSE};

    tempCounter                                                 = 0;
    tempString                                                  = "";
    for (Int_t i=0; i<nParticles; i++) {
        TString tempString                                      = enableParticles(i,1);

        if (tempString.CompareTo("1") == 0)
            includeParticle[i]                                  = kTRUE;
        else
            tempCounter++;
    }

    if (tempCounter == nParticles) {
        cout << "Warning: No particles selected, stopping!" << endl;
        return;
    }


    //================================================================================================================
    //Creating output file structure
    //================================================================================================================

    TFile *outputFile                                           = new TFile("CocktailInputPPb.root","RECREATE");
    TList *lists5TeV[nCentralities]                             = {NULL};
    TList *lists8TeV[nCentralities]                             = {NULL};

    if (enable5TeV) {
        for (Int_t i=0; i<nCentralities; i++) {
            if (includeCentrality5TeV[i]) {
                lists5TeV[i]                                    = new TList();
                lists5TeV[i]->SetName(Form("%s_%s_%s", fCollSys[1].Data(), fEnergy[2].Data(), fCentrality[i].Data()));
            }
        }
    }

    if (enable8TeV) {
        for (Int_t i=0; i<nCentralities; i++) {
            if (includeCentrality8TeV[i]) {
                lists8TeV[i]                                    = new TList();
                lists8TeV[i]->SetName(Form("%s_%s_%s", fCollSys[1].Data(), fEnergy[4].Data(), fCentrality[i].Data()));
            }
        }
    }

    //================================================================================================================
    // Set correct xSection for pPb
    //================================================================================================================
    // must be changed to the right pPb values!
//     xSection = 2.09 barn (DMeson paper)
    Double_t xSection5TeVNSD                                    = 2.09*1e12;
    Double_t scaleFacNSDToMB                                    = 1./0.964;
    Double_t scaleFacMBToNSD                                    = 0.964;
    //================================================================================================================
    // Set input files
    //================================================================================================================
    TDirectoryFile* tempList                                    = NULL;

    TFile *fileNeutralMeson5TeV                                 = NULL;
    TFile *fileNeutralMeson5TeVMB                               = NULL;
    TFile *fileNeutralMeson5TeVMBNew                            = NULL;
    if (includeParticle[0] || includeParticle[1]) {
        if (enable5TeV) {
            fileNeutralMeson5TeV                                = new TFile("pPb/5_TeV/neutralMeson_data_PCMResults_pPb_DPG2014_EtaBinning.root");
            fileNeutralMeson5TeVMB                              = new TFile("pPb/5_TeV/neutralMeson_data_PCMResults_pPb_20151111_standard_CatErrors_MB.root");
            fileNeutralMeson5TeVMBNew                           = new TFile("pPb/5_TeV/CombinedResultsPaperPPb5023GeV_2017_09_01.root");
        }
    }

    // input files particles 2 - 4

    TFile *fileChargedPionKaonProton5TeV                        = NULL;
    TFile *fileChargedPionKaonProton5TeVRatios                  = NULL;
    if (includeParticle[5] || includeParticle[6] || includeParticle[7]) {
        if (enable5TeV) {
            fileChargedPionKaonProton5TeV                       = new TFile("pPb/5_TeV/pPb502.fullpT.INEL.20151204.root");
            fileChargedPionKaonProton5TeVRatios                 = new TFile("pPb/5_TeV/pPb502.fullpT.RATIOS.20151204.root");
        }
    }

    // input files particles 8 - 14

    TFile *fileK0s5TeV                                          = NULL;
    if (includeParticle[15]) {
        if (enable5TeV)
            fileK0s5TeV                                         = new TFile("pPb/5_TeV/20130721-K0s-pPb.root");
    }

    TFile *fileLambda5TeV                                       = NULL;
    TFile *fileAntiLambda5TeV                                   = NULL;
    TFile *fileLambda5TeVRatios                                 = NULL;
    if (includeParticle[16]) {
        if (enable5TeV) {
            fileLambda5TeV                                      = new TFile("pPb/5_TeV/20130721-Lambda-pPb.root");
            fileAntiLambda5TeV                                  = new TFile("pPb/5_TeV/20130721-AntiLambda-pPb.root");
            fileLambda5TeVRatios                                = new TFile("pPb/5_TeV/20130721-LambdaOverKaon-pPb.root");
        }
    }

    // histos with published mean pt and integrated yield per particle type
    TH1D* histoMeanPtPerParticleStat[2][nCentralities];
    TH1D* histoMeanPtPerParticleSys[2][nCentralities];
    TH1D* histoIntegYieldPerParticleStat[2][nCentralities];
    TH1D* histoIntegYieldPerParticleSys[2][nCentralities];
    Bool_t bMeanPtPub[2][nCentralities];
    Bool_t bIntegYieldPub[2][nCentralities];
    for (Int_t e = 0; e < 2; e++){
        for (Int_t i = 0; i < nCentralities; i++){
            bMeanPtPub[e][i]                                    = kFALSE;
            bIntegYieldPub[e][i]                                = kFALSE;
            histoMeanPtPerParticleStat[e][i]                    = new TH1D(Form("meanPtPerParticle_pub_stat-%d-%d", e, i), "", nParticles, 0.5, nParticles+0.5);
            histoMeanPtPerParticleStat[e][i]->GetYaxis()->SetTitle("#LT p_{T} #RT");
            histoMeanPtPerParticleSys[e][i]                     = new TH1D(Form("meanPtPerParticle_pub_sys-%d-%d", e, i), "", nParticles, 0.5, nParticles+0.5);
            histoMeanPtPerParticleSys[e][i]->GetYaxis()->SetTitle("#LT p_{T} #RT");
            histoIntegYieldPerParticleStat[e][i]                = new TH1D(Form("IntegYieldPerParticle_pub_stat-%d-%d", e, i), "", nParticles, 0.5, nParticles+0.5);
            histoIntegYieldPerParticleStat[e][i]->GetYaxis()->SetTitle("d#it{N}/d#it{y}");
            histoIntegYieldPerParticleSys[e][i]                 = new TH1D(Form("IntegYieldPerParticle_pub_sys-%d-%d", e, i), "", nParticles, 0.5, nParticles+0.5);
            histoIntegYieldPerParticleSys[e][i]->GetYaxis()->SetTitle("d#it{N}/d#it{y}");
            for (Int_t particleIter = 0; particleIter < nParticles; particleIter++){
                histoMeanPtPerParticleStat[e][i]->GetXaxis()->SetBinLabel(particleIter+1,fParticleLatex[particleIter]);
                histoMeanPtPerParticleStat[e][i]->SetBinContent(particleIter+1, -1);
                histoMeanPtPerParticleSys[e][i]->GetXaxis()->SetBinLabel(particleIter+1,fParticleLatex[particleIter]);
                histoMeanPtPerParticleSys[e][i]->SetBinContent(particleIter+1, -1);
                histoIntegYieldPerParticleStat[e][i]->GetXaxis()->SetBinLabel(particleIter+1,fParticleLatex[particleIter]);
                histoIntegYieldPerParticleStat[e][i]->SetBinContent(particleIter+1, -1);
                histoIntegYieldPerParticleSys[e][i]->GetXaxis()->SetBinLabel(particleIter+1,fParticleLatex[particleIter]);
                histoIntegYieldPerParticleSys[e][i]->SetBinContent(particleIter+1, -1);
            }
        }
    }

    //================================================================================================================
    // creating histos and graphs for 5TeV
    //================================================================================================================
    if (enable5TeV) {

        for (Int_t i=0; i<nCentralities; i++) {
            if (includeCentrality5TeV[i]) {

                cout << fEnergy[2].Data() << "_" << fCentrality[i].Data() << endl;

                if (enablePtY) {
                    //================================================================================================================
                    // read in pythia pt y distributions for particles (i.e. PWGGA/GammaConv/AliAnalysisTaskGammaPureMC output)
                    //================================================================================================================
                    TFile*          fFilePtYpPb5TeV                             = new TFile("pPb/5_TeV/pPb5TeV_PtYDistributions_DPMJet.root");
                    TDirectoryFile* fDirPtYpPb5TeV                              = (TDirectoryFile*)fFilePtYpPb5TeV->Get("GammaPureMC");
                    TList*          fListPtYpPb5TeV                             = (TList*)fDirPtYpPb5TeV->Get("GammaPureMC");

                    Int_t       nBinsRebinX                                     = 5;
                    Int_t       nBinsRebinY                                     = 4;
                    Double_t    relStatErrThresh                                = 0.01;

                    TH2F* histPtYPi0                                            = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_Pi0");
                    histPtYPi0->Sumw2();
                    histPtYPi0->SetName("111_pt_y");
                    NormalizePtYHistogram(histPtYPi0, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYEta                                            = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_Eta");
                    histPtYEta->Sumw2();
                    histPtYEta->SetName("221_pt_y");
                    NormalizePtYHistogram(histPtYEta, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYEtaPrim                                        = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_EtaPrim");
                    histPtYEtaPrim->Sumw2();
                    histPtYEtaPrim->SetName("331_pt_y");
                    NormalizePtYHistogram(histPtYEtaPrim, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYOmega                                          = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_Omega");
                    histPtYOmega->Sumw2();
                    histPtYOmega->SetName("223_pt_y");
                    NormalizePtYHistogram(histPtYOmega, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYRho0                                           = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_Rho0");
                    histPtYRho0->Sumw2();
                    histPtYRho0->SetName("113_pt_y");
                    NormalizePtYHistogram(histPtYRho0, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYRhoPl                                          = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_RhoPl");
                    histPtYRhoPl->Sumw2();
                    histPtYRhoPl->SetName("213_pt_y");
                    NormalizePtYHistogram(histPtYRhoPl, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYRhoMi                                          = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_RhoMi");
                    histPtYRhoMi->Sumw2();
                    histPtYRhoMi->SetName("-213_pt_y");
                    NormalizePtYHistogram(histPtYRhoMi, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYPhi                                            = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_Phi");
                    histPtYPhi->Sumw2();
                    histPtYPhi->SetName("333_pt_y");
                    NormalizePtYHistogram(histPtYPhi, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYJPsi                                           = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_JPsi");
                    histPtYJPsi->Sumw2();
                    histPtYJPsi->SetName("443_pt_y");
                    NormalizePtYHistogram(histPtYJPsi, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYSigma0                                         = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_Sigma0");
                    histPtYSigma0->Sumw2();
                    histPtYSigma0->SetName("3212_pt_y");
                    NormalizePtYHistogram(histPtYSigma0, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYK0s                                            = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_K0s");
                    histPtYK0s->Sumw2();
                    histPtYK0s->SetName("310_pt_y");
                    NormalizePtYHistogram(histPtYK0s, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYK0l                                            = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_K0l");
                    histPtYK0l->Sumw2();
                    histPtYK0l->SetName("130_pt_y");
                    NormalizePtYHistogram(histPtYK0l, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYK0star                                         = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_K0star");
                    histPtYK0star->Sumw2();
                    histPtYK0star->SetName("313_pt_y");
                    NormalizePtYHistogram(histPtYK0star, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYDeltaPlPl                                      = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_DeltaPlPl");
                    histPtYDeltaPlPl->Sumw2();
                    histPtYDeltaPlPl->SetName("2224_pt_y");
                    NormalizePtYHistogram(histPtYDeltaPlPl, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYDeltaPl                                        = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_DeltaPl");
                    histPtYDeltaPl->Sumw2();
                    histPtYDeltaPl->SetName("2214_pt_y");
                    NormalizePtYHistogram(histPtYDeltaPl, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYDeltaMi                                        = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_DeltaMi");
                    histPtYDeltaMi->Sumw2();
                    histPtYDeltaMi->SetName("1114_pt_y");
                    NormalizePtYHistogram(histPtYDeltaMi, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYDelta0                                         = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_Delta0");
                    histPtYDelta0->Sumw2();
                    histPtYDelta0->SetName("2114_pt_y");
                    NormalizePtYHistogram(histPtYDelta0, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    TH2F* histPtYLambda                                         = (TH2F*)fListPtYpPb5TeV->FindObject("Pt_Y_Lambda");
                    histPtYLambda->Sumw2();
                    histPtYLambda->SetName("3122_pt_y");
                    NormalizePtYHistogram(histPtYLambda, nBinsRebinX, nBinsRebinY, relStatErrThresh);

                    lists5TeV[i]->Add(histPtYPi0);
                    lists5TeV[i]->Add(histPtYEta);
                    lists5TeV[i]->Add(histPtYEtaPrim);
                    lists5TeV[i]->Add(histPtYOmega);
                    lists5TeV[i]->Add(histPtYRho0);
                    lists5TeV[i]->Add(histPtYRhoPl);
                    lists5TeV[i]->Add(histPtYRhoMi);
                    lists5TeV[i]->Add(histPtYPhi);
                    lists5TeV[i]->Add(histPtYJPsi);
                    lists5TeV[i]->Add(histPtYSigma0);
                    lists5TeV[i]->Add(histPtYK0s);
                    lists5TeV[i]->Add(histPtYK0l);
                    lists5TeV[i]->Add(histPtYK0star);
                    lists5TeV[i]->Add(histPtYDeltaPlPl);
                    lists5TeV[i]->Add(histPtYDeltaPl);
                    lists5TeV[i]->Add(histPtYDeltaMi);
                    lists5TeV[i]->Add(histPtYDelta0);
                    lists5TeV[i]->Add(histPtYLambda);
                }

                //================================================================================================================
                // reading and writing pi0 and eta to 5TeV list
                // input is given as fully invariant yield
                //================================================================================================================
                if (includeParticle[0]) {
                    TH1D* histoNPion5TeVStat                            = NULL;
                    TGraphAsymmErrors* graphNPion5TeVStat               = NULL;
                    TGraphAsymmErrors* graphNPion5TeVSys                = NULL;

                    if (i==0 || i == 16) {
                        // MB
                        Bool_t useOldMBFile                             = kFALSE;
                        if (useOldMBFile) {
                            if (fileNeutralMeson5TeVMB->GetListOfKeys()->Contains("Pi0_pPb_5.023TeV_0-100%")) {

                                tempList                                = (TDirectoryFile*)fileNeutralMeson5TeVMB->Get("Pi0_pPb_5.023TeV_0-100%");

                                // spectra
                                if (tempList->Get("CorrectedYieldPi0") && tempList->Get("Pi0SystError")) {
                                    cout << " - pi^0 spectrum" << endl;

                                    histoNPion5TeVStat                  = (TH1D*)tempList->Get("CorrectedYieldPi0");
                                    graphNPion5TeVSys                   = (TGraphAsymmErrors*)tempList->Get("Pi0SystError");

                                    histoNPion5TeVStat                  = ConvertYieldHisto(histoNPion5TeVStat,kFALSE, kFALSE, kTRUE, kTRUE);
                                    graphNPion5TeVSys                   = ConvertYieldGraph(graphNPion5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                                    // Scale to MBAND
                                    if (i == 0){
                                        histoNPion5TeVStat->Scale(scaleFacNSDToMB);
                                        graphNPion5TeVSys = ScaleGraph(graphNPion5TeVSys,scaleFacNSDToMB);
                                    }

                                    SetHistoProperties(histoNPion5TeVStat, Form("%sStat", fParticle[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                                    SetGraphProperties(graphNPion5TeVSys, Form("%sSys", fParticle[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");


                                    lists5TeV[i]->Add(histoNPion5TeVStat);
                                    lists5TeV[i]->Add(graphNPion5TeVSys);
                                }
                            }
                        }

                        if (!useOldMBFile) {
                            if (fileNeutralMeson5TeVMBNew->GetListOfKeys()->Contains("Pi0pPb_5.023TeV")) {

                                tempList                                    = (TDirectoryFile*)fileNeutralMeson5TeVMBNew->Get("Pi0pPb_5.023TeV");

                                TString tempMethod[7]                       = {"PCM", "PHOS", "EMCAL", "PCMEMCAL", "Dalitz", "PCMPHOS" "Comb"};
                                TString tempReadMethod[7]                   = {"PCM", "PHOS", "EMC", "PCM-EMC", "PCM-Dal", "PCM-PHOS","Comb"};
                                // spectra
                                for (Int_t j=0; j<8; j++) {
                                    TString     tempName                    = Form("graphInvYieldPi0%spPb5023GeV", tempReadMethod[j].Data());

                                    if (tempList->Get(Form("%sStatErr", tempName.Data())) && tempList->Get(Form("%sSysErr", tempName.Data()))) {
                                        cout << " - pi^0 spectrum " << tempMethod[j];

                                        if (j==1) {
                                            cout << " -> will remove last (empty) bin" << endl;

                                            // empty bin in PHOS graph from 20.0G Gev/c to 25.0 GeV/c
                                            TGraphAsymmErrors* graphNPion5TeVStatTemp   = (TGraphAsymmErrors*)tempList->Get(Form("%sStatErr",   tempName.Data()));
                                            TGraphAsymmErrors* graphNPion5TeVSysTemp    = (TGraphAsymmErrors*)tempList->Get(Form("%sSysErr",    tempName.Data()));

                                            // stat
                                            Int_t       nPointsStat         = graphNPion5TeVStatTemp->GetN();
                                            Double_t*   xStat               = graphNPion5TeVStatTemp->GetX();
                                            Double_t*   xErrUpStat          = graphNPion5TeVStatTemp->GetEXhigh();
                                            Double_t*   xErrDownStat        = graphNPion5TeVStatTemp->GetEXlow();
                                            Double_t*   yStat               = graphNPion5TeVStatTemp->GetY();
                                            Double_t*   yErrDownStat        = graphNPion5TeVStatTemp->GetEYlow();
                                            Double_t*   yErrUpStat          = graphNPion5TeVStatTemp->GetEYhigh();

                                            Double_t*   xStatNew            = new Double_t[nPointsStat-1];
                                            Double_t*   xErrUpStatNew       = new Double_t[nPointsStat-1];
                                            Double_t*   xErrDownStatNew     = new Double_t[nPointsStat-1];
                                            Double_t*   yStatNew            = new Double_t[nPointsStat-1];
                                            Double_t*   yErrUpStatNew       = new Double_t[nPointsStat-1];
                                            Double_t*   yErrDownStatNew     = new Double_t[nPointsStat-1];

                                            for (Int_t ii=0; ii<nPointsStat-1; ii++) {
                                                xStatNew[ii]                = xStat[ii];
                                                xErrUpStatNew[ii]           = xErrUpStat[ii];
                                                xErrDownStatNew[ii]         = xErrDownStat[ii];

                                                yStatNew[ii]                = yStat[ii];
                                                yErrUpStatNew[ii]           = yErrUpStat[ii];
                                                yErrDownStatNew[ii]         = yErrDownStat[ii];
                                            }
                                            graphNPion5TeVStat              = new TGraphAsymmErrors(nPointsStat-1, xStatNew, yStatNew, xErrDownStatNew, xErrUpStatNew, yErrDownStatNew, yErrUpStatNew);

                                            // sys
                                            Int_t       nPointsSys          = graphNPion5TeVSysTemp->GetN();
                                            Double_t*   xSys                = graphNPion5TeVSysTemp->GetX();
                                            Double_t*   xErrUpSys           = graphNPion5TeVSysTemp->GetEXhigh();
                                            Double_t*   xErrDownSys         = graphNPion5TeVSysTemp->GetEXlow();
                                            Double_t*   ySys                = graphNPion5TeVSysTemp->GetY();
                                            Double_t*   yErrUpSys           = graphNPion5TeVSysTemp->GetEYhigh();
                                            Double_t*   yErrDownSys         = graphNPion5TeVSysTemp->GetEYlow();

                                            Double_t*   xSysNew             = new Double_t[nPointsSys-1];
                                            Double_t*   xErrUpSysNew        = new Double_t[nPointsSys-1];
                                            Double_t*   xErrDownSysNew      = new Double_t[nPointsSys-1];
                                            Double_t*   ySysNew             = new Double_t[nPointsSys-1];
                                            Double_t*   yErrUpSysNew        = new Double_t[nPointsSys-1];
                                            Double_t*   yErrDownSysNew      = new Double_t[nPointsSys-1];

                                            for (Int_t ii=0; ii<nPointsSys-1; ii++) {
                                                xSysNew[ii]                 = xSys[ii];
                                                xErrUpSysNew[ii]            = xErrUpSys[ii];
                                                xErrDownSysNew[ii]          = xErrDownSys[ii];

                                                ySysNew[ii]                 = ySys[ii];
                                                yErrUpSysNew[ii]            = yErrUpSys[ii];
                                                yErrDownSysNew[ii]          = yErrDownSys[ii];
                                            }
                                            graphNPion5TeVSys               = new TGraphAsymmErrors(nPointsSys-1, xSysNew, ySysNew, xErrDownSysNew, xErrUpSysNew, yErrDownSysNew, yErrUpSysNew);

                                        } else {
                                            cout << endl;

                                            graphNPion5TeVStat              = (TGraphAsymmErrors*)tempList->Get(Form("%sStatErr",   tempName.Data()));
                                            graphNPion5TeVSys               = (TGraphAsymmErrors*)tempList->Get(Form("%sSysErr",    tempName.Data()));
                                        }

                                        graphNPion5TeVStat                  = ConvertYieldGraph(graphNPion5TeVStat,kFALSE, kFALSE, kTRUE, kTRUE);
                                        graphNPion5TeVSys                   = ConvertYieldGraph(graphNPion5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                                        TString         tempMethodName      = "";
                                        if (j==0)       tempMethodName      = fMethod[2];   // PCM
                                        else if (j==1)  tempMethodName      = fMethod[3];   // PHOS
                                        else if (j==2)  tempMethodName      = fMethod[4];   // EMCal
                                        else if (j==3)  tempMethodName      = fMethod[6];   // PCMEMCal
                                        else if (j==4)  tempMethodName      = fMethod[10];  // Dalitz
                                        else if (j==5)  tempMethodName      = fMethod[5];   // PCMPHOS
                                        else if (j==6)  tempMethodName      = fMethod[1];   // Comb

                                        // Scale to MBAND
                                        if (i == 0){
                                            graphNPion5TeVStat                  = ScaleGraph(graphNPion5TeVStat,scaleFacNSDToMB);
                                            graphNPion5TeVSys                   = ScaleGraph(graphNPion5TeVSys,scaleFacNSDToMB);
                                        }

                                        SetGraphProperties(graphNPion5TeVStat, Form("%s%sStat", fParticle[0].Data(), tempMethodName.Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                                        SetGraphProperties(graphNPion5TeVSys, Form("%s%sSys", fParticle[0].Data(), tempMethodName.Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                                        lists5TeV[i]->Add(graphNPion5TeVStat);
                                        lists5TeV[i]->Add(graphNPion5TeVSys);
                                    }
                                }
                            }
                        }

                    } else {
                        // centrality classes
                        if (fileNeutralMeson5TeV->GetListOfKeys()->Contains(Form("Pi0_pPb_5.023TeV_%s", fCentralityLatex[i].Data()))) {

                            tempList                                    = (TDirectoryFile*)fileNeutralMeson5TeV->Get(Form("Pi0_pPb_5.023TeV_%s", fCentralityLatex[i].Data()));

                            // spectra
                            if (tempList->Get("CorrectedYieldPi0") && tempList->Get("Pi0SystError")) {
                                cout << " - pi^0 spectrum" << endl;

                                histoNPion5TeVStat                      = (TH1D*)tempList->Get("CorrectedYieldPi0");
                                graphNPion5TeVSys                       = (TGraphAsymmErrors*)tempList->Get("Pi0SystError");

                                histoNPion5TeVStat                      = ConvertYieldHisto(histoNPion5TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                graphNPion5TeVSys                       = ConvertYieldGraph(graphNPion5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                                SetHistoProperties(histoNPion5TeVStat, Form("%sStat", fParticle[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                                SetGraphProperties(graphNPion5TeVSys, Form("%sSys", fParticle[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                                lists5TeV[i]->Add(histoNPion5TeVStat);
                                lists5TeV[i]->Add(graphNPion5TeVSys);
                            }
                        }
                    }
                }

                if (includeParticle[1]) {
                    TH1D* histoEta5TeVStat                              = NULL;
                    TGraphAsymmErrors* graphEta5TeVStat                 = NULL;
                    TGraphAsymmErrors* graphEta5TeVSys                  = NULL;

                    TH1D* histoEtaToPi05TeVStat                         = NULL;
                    TGraphAsymmErrors* graphEtaToPi05TeVStat            = NULL;
                    TGraphAsymmErrors* graphEtaToPi05TeVSys             = NULL;

                    if (i==0 || i == 16) {
                        // MB

                        Bool_t useOldMBFile                             = kFALSE;
                        if (useOldMBFile) {
                            if (fileNeutralMeson5TeVMB->GetListOfKeys()->Contains("Eta_pPb_5.023TeV_0-100%")) {

                                tempList                                    = (TDirectoryFile*)fileNeutralMeson5TeVMB->Get("Eta_pPb_5.023TeV_0-100%");

                                // spectra
                                if (tempList->Get("CorrectedYieldEta") && tempList->Get("EtaSystError")) {
                                    cout << " - eta spectrum" << endl;

                                    histoEta5TeVStat                        = (TH1D*)tempList->Get("CorrectedYieldEta");
                                    graphEta5TeVSys                         = (TGraphAsymmErrors*)tempList->Get("EtaSystError");

                                    histoEta5TeVStat                        = ConvertYieldHisto(histoEta5TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                    graphEta5TeVSys                         = ConvertYieldGraph(graphEta5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                                    // Scale to MBAND
                                    if (i == 0){
                                        histoEta5TeVStat->Scale(scaleFacNSDToMB);
                                        graphEta5TeVSys                     = ScaleGraph(graphEta5TeVSys,scaleFacNSDToMB);
                                    }
                                    SetHistoProperties(histoEta5TeVStat, Form("%sStat", fParticle[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                                    SetGraphProperties(graphEta5TeVSys, Form("%sSys", fParticle[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");


                                    lists5TeV[i]->Add(histoEta5TeVStat);
                                    lists5TeV[i]->Add(graphEta5TeVSys);
                                }

                                // ratio
                                if (tempList->Get("EtatoPi0Ratio") && tempList->Get("EtatoPi0RatioSys")) {
                                    cout << " - eta ratio" << endl;

                                    histoEtaToPi05TeVStat                   = (TH1D*)tempList->Get("EtatoPi0Ratio");
                                    graphEtaToPi05TeVSys                    = (TGraphAsymmErrors*)tempList->Get("EtatoPi0RatioSys");

                                    SetHistoProperties(histoEtaToPi05TeVStat, Form("%sTo%sStat", fParticle[1].Data(), fParticle[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                                    SetGraphProperties(graphEtaToPi05TeVSys, Form("%sTo%sSys", fParticle[1].Data(), fParticle[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");

                                    lists5TeV[i]->Add(histoEtaToPi05TeVStat);
                                    lists5TeV[i]->Add(graphEtaToPi05TeVSys);
                                }
                            }
                        }

                        if (!useOldMBFile) {
                            if (fileNeutralMeson5TeVMBNew->GetListOfKeys()->Contains("EtapPb_5.023TeV")) {

                                tempList                                    = (TDirectoryFile*)fileNeutralMeson5TeVMBNew->Get("EtapPb_5.023TeV");

                                TString tempMethod[5]                       = {"PCM", "EMCAL", "PCMEMCAL", "PCMPHOS", "Comb"};
                                TString tempReadMethod[5]                   = {"PCM", "EMC", "PCM-EMC", "PCM-PHOS", "Comb"};

                                // spectra
                                for (Int_t j=0; j<5; j++) {

                                    TString     tempName                    = Form("graphInvYieldEta%spPb5023GeV",  tempReadMethod[j].Data());

                                    if (tempList->Get(Form("%sStatErr", tempName.Data())) && tempList->Get(Form("%sSysErr", tempName.Data()))) {
                                        cout << " - eta spectrum " << tempMethod[j] << endl;

                                        graphEta5TeVStat                    = (TGraphAsymmErrors*)tempList->Get(Form("%sStatErr",   tempName.Data()));
                                        graphEta5TeVSys                     = (TGraphAsymmErrors*)tempList->Get(Form("%sSysErr",    tempName.Data()));

                                        graphEta5TeVStat                    = ConvertYieldGraph(graphEta5TeVStat,kFALSE, kFALSE, kTRUE, kTRUE);
                                        graphEta5TeVSys                     = ConvertYieldGraph(graphEta5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                                        TString         tempMethodName      = "";
                                        if (j==0)       tempMethodName      = fMethod[2];   // PCM
                                        else if (j==1)  tempMethodName      = fMethod[4];   // EMCal
                                        else if (j==2)  tempMethodName      = fMethod[6];   // PCMEMCal
                                        else if (j==3)  tempMethodName      = fMethod[5];   // PCMPHOS
                                        else if (j==4)  tempMethodName      = fMethod[1];   // Comb


                                        // Scale to MBAND
                                        if (i == 0){
                                            graphEta5TeVStat                    = ScaleGraph(graphEta5TeVStat,scaleFacNSDToMB);
                                            graphEta5TeVSys                     = ScaleGraph(graphEta5TeVSys,scaleFacNSDToMB);
                                        }

                                        SetGraphProperties(graphEta5TeVStat, Form("%s%sStat", fParticle[1].Data(), tempMethodName.Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                                        SetGraphProperties(graphEta5TeVSys, Form("%s%sSys", fParticle[1].Data(), tempMethodName.Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");


                                        lists5TeV[i]->Add(graphEta5TeVStat);
                                        lists5TeV[i]->Add(graphEta5TeVSys);
                                    }
                                }


                                // ratio
                                for (Int_t j=0; j<5; j++) {

                                    TString     tempName                    = "histo";
                                    if (j==4)   tempName                    = "graph";

                                    if (tempList->Get(Form("%sRatioEtaToPi0%spPb5023GeVStatErr", tempName.Data(), tempMethod[j].Data())) && tempList->Get(Form("graphRatioEtaToPi0%spPb5023GeVSysErr", tempMethod[j].Data()))) {
                                        cout << " - eta ratio " << tempMethod[j] << endl;

                                        TString         tempMethodName      = "";
                                        if (j==0)       tempMethodName      = fMethod[2];   // PCM
                                        else if (j==1)  tempMethodName      = fMethod[4];   // EMCal
                                        else if (j==2)  tempMethodName      = fMethod[6];   // PCMEMCal
                                        else if (j==3)  tempMethodName      = fMethod[5];   // PCMPHOS
                                        else if (j==4)  tempMethodName      = fMethod[1];   // Comb

                                        if (j!=4) {
                                            histoEtaToPi05TeVStat                       = (TH1D*)tempList->Get(Form("histoRatioEtaToPi0%spPb5023GeVStatErr", tempMethod[j].Data()));
                                            TGraphAsymmErrors* graphEtaToPi05TeVStat    = new TGraphAsymmErrors(histoEtaToPi05TeVStat);
                                            while(graphEtaToPi05TeVStat->GetY()[0] == 0 && graphEtaToPi05TeVStat->GetN() > 0) graphEtaToPi05TeVStat->RemovePoint(0);
                                            graphEtaToPi05TeVStat->Set(graphEtaToPi05TeVStat->GetN()+1);
                                            graphEtaToPi05TeVStat->SetPoint(graphEtaToPi05TeVStat->GetN()-1, 0.01, 0.0001);
                                            graphEtaToPi05TeVStat->SetPointError(graphEtaToPi05TeVStat->GetN()-1, 0.01, 0.01, 0.01, 0.01);
                                            graphEtaToPi05TeVStat->Sort();
                                            SetGraphProperties(graphEtaToPi05TeVStat, Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), tempMethodName.Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                                            lists5TeV[i]->Add(graphEtaToPi05TeVStat);
                                        } else {
                                            TGraphAsymmErrors* graphEtaToPi05TeVStat    = (TGraphAsymmErrors*)tempList->Get(Form("graphRatioEtaToPi0%spPb5023GeVStatErr", tempMethod[j].Data()));
                                            graphEtaToPi05TeVStat->Set(graphEtaToPi05TeVStat->GetN()+1);
                                            graphEtaToPi05TeVStat->SetPoint(graphEtaToPi05TeVStat->GetN()-1, 0.01, 0.0001);
                                            graphEtaToPi05TeVStat->SetPointError(graphEtaToPi05TeVStat->GetN()-1, 0.01, 0.01, 0.01, 0.01);
                                            graphEtaToPi05TeVStat->Sort();
                                            SetGraphProperties(graphEtaToPi05TeVStat, Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), tempMethodName.Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                                            lists5TeV[i]->Add(graphEtaToPi05TeVStat);
                                        }
                                        TGraphAsymmErrors* graphEtaToPi05TeVSys     = (TGraphAsymmErrors*)tempList->Get(Form("graphRatioEtaToPi0%spPb5023GeVSysErr", tempMethod[j].Data()));
                                        graphEtaToPi05TeVSys->Set(graphEtaToPi05TeVSys->GetN()+1);
                                        graphEtaToPi05TeVSys->SetPoint(graphEtaToPi05TeVSys->GetN()-1, 0.01, 0.0001);
                                        graphEtaToPi05TeVSys->SetPointError(graphEtaToPi05TeVSys->GetN()-1, 0.01, 0.01, 0.2, 0.2);
                                        graphEtaToPi05TeVSys->Sort();
                                        SetGraphProperties(graphEtaToPi05TeVSys, Form("%sTo%s%sSys", fParticle[1].Data(), fParticle[0].Data(), tempMethodName.Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                                        lists5TeV[i]->Add(graphEtaToPi05TeVSys);
                                    }
                                }
                            }
                        }

                    } else {
                        // centrality classes
                        if (fileNeutralMeson5TeV->GetListOfKeys()->Contains(Form("Eta_pPb_5.023TeV_%s", fCentralityLatex[i].Data()))) {

                            tempList                                    = (TDirectoryFile*)fileNeutralMeson5TeV->Get(Form("Eta_pPb_5.023TeV_%s", fCentralityLatex[i].Data()));

                            // spectra
                            if (tempList->Get("CorrectedYieldEta") && tempList->Get("EtaSystError")) {
                                cout << " - eta spectrum" << endl;

                                histoEta5TeVStat                        = (TH1D*)tempList->Get("CorrectedYieldEta");
                                graphEta5TeVSys                         = (TGraphAsymmErrors*)tempList->Get("EtaSystError");

                                histoEta5TeVStat                        = ConvertYieldHisto(histoEta5TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                graphEta5TeVSys                         = ConvertYieldGraph(graphEta5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                                SetHistoProperties(histoEta5TeVStat, Form("%sStat", fParticle[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                                SetGraphProperties(graphEta5TeVSys, Form("%sSys", fParticle[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                                lists5TeV[i]->Add(histoEta5TeVStat);
                                lists5TeV[i]->Add(graphEta5TeVSys);
                            }

                            // ratio
                            if (tempList->Get("EtatoPi0Ratio") && tempList->Get("EtatoPi0RatioSys")) {
                                cout << " - eta ratio" << endl;

                                histoEtaToPi05TeVStat                   = (TH1D*)tempList->Get("EtatoPi0Ratio");
                                graphEtaToPi05TeVSys                    = (TGraphAsymmErrors*)tempList->Get("EtatoPi0RatioSys");

                                SetHistoProperties(histoEtaToPi05TeVStat, Form("%sTo%sStat", fParticle[1].Data(), fParticle[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                                SetGraphProperties(graphEtaToPi05TeVSys, Form("%sTo%sSys", fParticle[1].Data(), fParticle[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");

                                lists5TeV[i]->Add(histoEtaToPi05TeVStat);
                                lists5TeV[i]->Add(graphEtaToPi05TeVSys);
                            }
                        }
                    }
                }

                //================================================================================================================
                // reading and writing omega to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[2]) {
                    cout << " - omega" << endl;
                }

                //================================================================================================================
                // reading and writing eta! to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[3]) {
                    cout << " - eta'" << endl;
                }

                //================================================================================================================
                // reading and writing gamma_dir to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[4]) {
                    cout << " - gamma_dir" << endl;
                }

                //================================================================================================================
                // reading and writing pi+-, K+- and p/bar{p} to 5TeV list
                // input is given as fully invariant yield
                //================================================================================================================
                if (includeParticle[5]) {
                    // spectra
                    TString fCentralityCurrent                      = fCentrality[i].Data();
                    if (i == 16)
                        fCentralityCurrent                          = fCentrality[0].Data();
                    if (fileChargedPionKaonProton5TeV->GetListOfKeys()->Contains(Form("hstat_pPb502_%s_pion_sum", fCentralityCurrent.Data())) && fileChargedPionKaonProton5TeV->GetListOfKeys()->Contains(Form("hsys_pPb502_%s_pion_sum", fCentralityCurrent.Data()))) {
                        cout << " - pi^+/- spectrum" << endl;

                        TH1D* histoCPionComb5TeVStat                = (TH1D*)fileChargedPionKaonProton5TeV->Get(Form("hstat_pPb502_%s_pion_sum", fCentralityCurrent.Data()));
                        TH1D* histoCPionComb5TeVSys                 = (TH1D*)fileChargedPionKaonProton5TeV->Get(Form("hsys_pPb502_%s_pion_sum", fCentralityCurrent.Data()));

                        histoCPionComb5TeVStat                      = ConvertYieldHisto(histoCPionComb5TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        histoCPionComb5TeVSys                       = ConvertYieldHisto(histoCPionComb5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        // scale by 0.5 to get averaged single particle
                        histoCPionComb5TeVStat->Scale(0.5);
                        histoCPionComb5TeVSys->Scale(0.5);

                        SetHistoProperties(histoCPionComb5TeVStat, Form("%sStat", fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoCPionComb5TeVSys, Form("%sSys", fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        if (i == 0){
                            // Scale to MBAND
                            histoCPionComb5TeVStat->Scale(scaleFacNSDToMB);
                            histoCPionComb5TeVSys->Scale(scaleFacNSDToMB);
                        }

                        lists5TeV[i]->Add(histoCPionComb5TeVStat);
                        lists5TeV[i]->Add(histoCPionComb5TeVSys);
                    }
                }

                if (includeParticle[6]) {
                    // spectra
                    TString fCentralityCurrent                      = fCentrality[i].Data();
                    if (i == 16)
                        fCentralityCurrent                          = fCentrality[0].Data();
                    if (fileChargedPionKaonProton5TeV->GetListOfKeys()->Contains(Form("hstat_pPb502_%s_kaon_sum", fCentralityCurrent.Data())) && fileChargedPionKaonProton5TeV->GetListOfKeys()->Contains(Form("hsys_pPb502_%s_kaon_sum", fCentralityCurrent.Data()))) {
                        cout << " - K^+/- spectrum" << endl;

                        TH1D* histoCKaonComb5TeVStat                = (TH1D*)fileChargedPionKaonProton5TeV->Get(Form("hstat_pPb502_%s_kaon_sum", fCentralityCurrent.Data()));
                        TH1D* histoCKaonComb5TeVSys                 = (TH1D*)fileChargedPionKaonProton5TeV->Get(Form("hsys_pPb502_%s_kaon_sum", fCentralityCurrent.Data()));

                        histoCKaonComb5TeVStat                      = ConvertYieldHisto(histoCKaonComb5TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        histoCKaonComb5TeVSys                       = ConvertYieldHisto(histoCKaonComb5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        // scale by 0.5 to get averaged single particle
                        histoCKaonComb5TeVStat->Scale(0.5);
                        histoCKaonComb5TeVSys->Scale(0.5);

                        SetHistoProperties(histoCKaonComb5TeVStat, Form("%sStat", fParticle[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoCKaonComb5TeVSys, Form("%sSys", fParticle[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        if (i == 0){
                            // Scale to MBAND
                            histoCKaonComb5TeVStat->Scale(scaleFacNSDToMB);
                            histoCKaonComb5TeVSys->Scale(scaleFacNSDToMB);
                        }

                        lists5TeV[i]->Add(histoCKaonComb5TeVStat);
                        lists5TeV[i]->Add(histoCKaonComb5TeVSys);
                    }

                    // ratios
                    if (fileChargedPionKaonProton5TeVRatios->GetListOfKeys()->Contains(Form("hstat_pPb502_%s_kaon_to_pion_sum", fCentralityCurrent.Data())) && fileChargedPionKaonProton5TeVRatios->GetListOfKeys()->Contains(Form("hsys_pPb502_%s_kaon_to_pion_sum", fCentralityCurrent.Data()))) {
                        cout << " - K^+/- ratio" << endl;

                        TH1D* histoCKaonToCPionComb5TeVStat         = (TH1D*)fileChargedPionKaonProton5TeVRatios->Get(Form("hstat_pPb502_%s_kaon_to_pion_sum", fCentralityCurrent.Data()));
                        TH1D* histoCKaonToCPionComb5TeVSys          = (TH1D*)fileChargedPionKaonProton5TeVRatios->Get(Form("hsys_pPb502_%s_kaon_to_pion_sum", fCentralityCurrent.Data()));

                        SetHistoProperties(histoCKaonToCPionComb5TeVStat, Form("%sTo%sStat", fParticle[6].Data(), fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})", "");
                        SetHistoProperties(histoCKaonToCPionComb5TeVSys, Form("%sTo%sSys", fParticle[6].Data(), fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})", "");

                        lists5TeV[i]->Add(histoCKaonToCPionComb5TeVStat);
                        lists5TeV[i]->Add(histoCKaonToCPionComb5TeVSys);
                    }
                }

                if (includeParticle[7]) {
                    // spectra
                    TString fCentralityCurrent                      = fCentrality[i].Data();
                    if (i == 16)
                        fCentralityCurrent                          = fCentrality[0].Data();
                    if (fileChargedPionKaonProton5TeV->GetListOfKeys()->Contains(Form("hstat_pPb502_%s_proton_sum", fCentralityCurrent.Data())) && fileChargedPionKaonProton5TeV->GetListOfKeys()->Contains(Form("hsys_pPb502_%s_proton_sum", fCentralityCurrent.Data()))) {
                        cout << " - p/anti-p spectrum" << endl;

                        TH1D* histoProtonComb5TeVStat               = (TH1D*)fileChargedPionKaonProton5TeV->Get(Form("hstat_pPb502_%s_proton_sum", fCentralityCurrent.Data()));
                        TH1D* histoProtonComb5TeVSys                = (TH1D*)fileChargedPionKaonProton5TeV->Get(Form("hsys_pPb502_%s_proton_sum", fCentralityCurrent.Data()));

                        histoProtonComb5TeVStat                     = ConvertYieldHisto(histoProtonComb5TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        histoProtonComb5TeVSys                      = ConvertYieldHisto(histoProtonComb5TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        // scale by 0.5 to get averaged single particle
                        histoProtonComb5TeVStat->Scale(0.5);
                        histoProtonComb5TeVSys->Scale(0.5);

                        SetHistoProperties(histoProtonComb5TeVStat, Form("%sStat", fParticle[7].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoProtonComb5TeVSys, Form("%sSys", fParticle[7].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        if (i == 0){
                            // Scale to MBAND
                            histoProtonComb5TeVStat->Scale(scaleFacNSDToMB);
                            histoProtonComb5TeVSys->Scale(scaleFacNSDToMB);
                        }

                        lists5TeV[i]->Add(histoProtonComb5TeVStat);
                        lists5TeV[i]->Add(histoProtonComb5TeVSys);
                    }

                    // ratios
                    if (fileChargedPionKaonProton5TeVRatios->GetListOfKeys()->Contains(Form("hstat_pPb502_%s_proton_to_pion_sum", fCentralityCurrent.Data())) && fileChargedPionKaonProton5TeVRatios->GetListOfKeys()->Contains(Form("hsys_pPb502_%s_proton_to_pion_sum", fCentralityCurrent.Data()))) {
                        cout << " - p/anti-p ratio" << endl;

                        TH1D* histoProtonToCKaonComb5TeVStat        = (TH1D*)fileChargedPionKaonProton5TeVRatios->Get(Form("hstat_pPb502_%s_proton_to_pion_sum", fCentralityCurrent.Data()));
                        TH1D* histoProtonToCKaonComb5TeVSys         = (TH1D*)fileChargedPionKaonProton5TeVRatios->Get(Form("hsys_pPb502_%s_proton_to_pion_sum", fCentralityCurrent.Data()));

                        SetHistoProperties(histoProtonToCKaonComb5TeVStat, Form("%sTo%sStat", fParticle[7].Data(), fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "(p + #bar{p}) / (#pi^{+} + #pi^{-})", "");
                        SetHistoProperties(histoProtonToCKaonComb5TeVSys, Form("%sTo%sSys", fParticle[7].Data(), fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "(p + #bar{p}) / (#pi^{+} + #pi^{-})", "");

                        lists5TeV[i]->Add(histoProtonToCKaonComb5TeVStat);
                        lists5TeV[i]->Add(histoProtonToCKaonComb5TeVSys);
                    }
                }

                //================================================================================================================
                // reading and writing h+- to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[8]) {
                    cout << " - h^+/-" << endl;
                    TString chargedHadHEPDataFile                               = "";
                    if (i == 0 || i == 16)
                        chargedHadHEPDataFile                                   = "pPb/5_TeV/HEPData/ChargedHadron_HEPData-ins1295687-v1-csv/Table1.csv";

                    if (chargedHadHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphChHadYieldpPbStat               = ParseHEPData(chargedHadHEPDataFile, 10, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphChHadYieldpPbSys                = ParseHEPData(chargedHadHEPDataFile, 10, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        graphChHadYieldpPbStat                                  = ConvertYieldGraph(graphChHadYieldpPbStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        graphChHadYieldpPbSys                                   = ConvertYieldGraph(graphChHadYieldpPbSys, kFALSE, kFALSE, kTRUE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphChHadYieldpPbStat                              = ScaleGraph(graphChHadYieldpPbStat, scaleFacNSDToMB*0.5);
                            graphChHadYieldpPbSys                               = ScaleGraph(graphChHadYieldpPbSys, scaleFacNSDToMB*0.5);
                        } else {
                            graphChHadYieldpPbStat                              = ScaleGraph(graphChHadYieldpPbStat, 0.5);
                            graphChHadYieldpPbSys                               = ScaleGraph(graphChHadYieldpPbSys, 0.5);
                        }
                        SetGraphProperties(graphChHadYieldpPbStat, Form("%sStat", fParticle[8].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphChHadYieldpPbSys, Form("%sSys", fParticle[8].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphChHadYieldpPbStat);
                        lists5TeV[i]->Add(graphChHadYieldpPbSys);
                    }

                }

                //================================================================================================================
                // reading and writing phi to 5TeV list
                // input is given as invariant yield
                //================================================================================================================
                if (includeParticle[9]) {

                    cout << " - phi " << endl;
                    //      MB, 0-5,    5-10,   0-10,   10-20,  0-20,   20-40,  20-50,  40-60,  60-80,  80-100
                    //      0,  1,      2,      3,      4,      5,      6,      7,      8,      9,      10
//                     yield
//                     0 - 5
//                     0.3767±0.0037stat±0.0196sys,uncorr±0.0226sys,corr
//                     5 - 10
//                     0.2878±0.003stat±0.0138sys,uncorr±0.0173sys,corr
//                     10 - 20
//                     0.2436±0.002stat±0.0118sys,uncorr±0.0146sys,corr
//                     20 - 40
//                     0.185±0.0012stat±0.0095sys,uncorr±0.0111sys,corr
//                     40 - 60
//                     0.12285±0.00085stat±0.0064sys,uncorr±0.00737sys,corr
//                     60 - 80
//                     0.06948±0.00069stat±0.00369sys,uncorr±0.00417sys,corr
//                     80 - 100
//                     0.02973±0.00045stat±0.00228sys,uncorr±0.00178sys,corr
//                     NSD
//                     0.13443±0.00049stat±0.00692sys,uncorr±0.00807sys,corr

//                     mean pt
//                     0 - 5
//                     1.437±0.009stat±0.0282sys
//                     5 - 10
//                     1.4421±0.0093stat±0.0254sys
//                     10 - 20
//                     1.4215±0.0075stat±0.024sys
//                     20 - 40
//                     1.3574±0.0056stat±0.0252sys
//                     40 - 60
//                     1.3097±0.0064stat±0.0307sys
//                     60 - 80
//                     1.2424±0.0078stat±0.0243sys
//                     80 - 100
//                     1.0549±0.0096stat±0.0303sys
//                     NSD
//                     1.3547±0.0032stat±0.0305sys

                    TString phiHEPDataFile                                      = "";
                    if (i == 0 || i == 16)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_NSD.csv";
                    else if (i == 1)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_0005.csv";
                    else if (i == 2)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_0510.csv";
                    else if (i == 4)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_1020.csv";
                    else if (i == 5)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_0020.csv";
                    else if (i == 6)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_2040.csv";
                    else if (i == 8)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_4060.csv";
                    else if (i == 9)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_6080.csv";
                    else if (i == 11)
                        phiHEPDataFile                                          = "pPb/5_TeV/HEPData/Phi_pPb_80100.csv";

                    if (phiHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphPhiYieldpPbStat                 = ParseHEPData(phiHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphPhiYieldpPbSys                  = ParseHEPData(phiHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphPhiYieldpPbStat                                = ScaleGraph(graphPhiYieldpPbStat, scaleFacNSDToMB);
                            graphPhiYieldpPbSys                                 = ScaleGraph(graphPhiYieldpPbSys, scaleFacNSDToMB);
                        }
                        SetGraphProperties(graphPhiYieldpPbStat, Form("%sStat", fParticle[9].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphPhiYieldpPbSys, Form("%sSys", fParticle[9].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphPhiYieldpPbStat);
                        lists5TeV[i]->Add(graphPhiYieldpPbSys);
                    }

                    TString pToPhiHEPDataFile                                   = "";
                    if (i == 1)
                        pToPhiHEPDataFile                                       = "pPb/5_TeV/HEPData/PToPhi_pPb_0005.csv";
                    else if (i == 2)
                        pToPhiHEPDataFile                                       = "pPb/5_TeV/HEPData/PToPhi_pPb_0510.csv";
                    else if (i == 3)
                        pToPhiHEPDataFile                                       = "pPb/5_TeV/HEPData/PToPhi_pPb_1020.csv";
                    else if (i == 6)
                        pToPhiHEPDataFile                                       = "pPb/5_TeV/HEPData/PToPhi_pPb_2040.csv";
                    else if (i == 8)
                        pToPhiHEPDataFile                                       = "pPb/5_TeV/HEPData/PToPhi_pPb_4060.csv";
                    else if (i == 9)
                        pToPhiHEPDataFile                                       = "pPb/5_TeV/HEPData/PToPhi_pPb_6080.csv";
                    else if (i == 11)
                        pToPhiHEPDataFile                                       = "pPb/5_TeV/HEPData/PToPhi_pPb_80100.csv";

                    if (pToPhiHEPDataFile.CompareTo("") != 0){
                        cout << "   added ratio to p " << endl;
                        TGraphAsymmErrors* graphPToPhiYieldpPbStat              = ParseHEPData(pToPhiHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphPToPhiYieldpPbSys               = ParseHEPData(pToPhiHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        SetGraphProperties(graphPToPhiYieldpPbStat, Form("%sTo%sStat", fParticle[7].Data(), fParticle[9].Data()), "#it{p}_{T} (GeV/#it{c})", "p/#phi", "");
                        SetGraphProperties(graphPToPhiYieldpPbSys, Form("%sTo%sSys", fParticle[7].Data(), fParticle[9].Data()), "#it{p}_{T} (GeV/#it{c})", "p/#phi", "");

                        lists5TeV[i]->Add(graphPToPhiYieldpPbStat);
                        lists5TeV[i]->Add(graphPToPhiYieldpPbSys);
                    }
                }

                //================================================================================================================
                // reading and writing K^*0 to 5TeV list
                // input is given as invariant yield
                //================================================================================================================
                if (includeParticle[10]) {
                    cout << " - K^*0" << endl;
                    //      MB, 0-5,    5-10,   0-10,   10-20,  0-20,   20-40,  20-50,  40-60,  60-80,  80-100
                    //      0,  1,      2,      3,      4,      5,      6,      7,      8,      9,      10
//                     yield
//                     0 - 20
//                     0.6162±0.0078stat±0.0371sys,uncorr±0.037sys,corr
//                     20 - 40
//                     0.4263±0.0056stat±0.0261sys,uncorr±0.0256sys,corr
//                     40 - 60
//                     0.3016±0.0039stat±0.0194sys,uncorr±0.0181sys,corr
//                     60 - 80
//                     0.1852±0.0026stat±0.0126sys,uncorr±0.0111sys,corr
//                     80 - 100
//                     0.0832±0.0013stat±0.0054sys,uncorr±0.005sys,corr
//                     NSD
//                     0.3154±0.002stat±0.0176sys,uncorr±0.0189sys,corr

//                     mean pt
//                     0 - 20
//                     1.379±0.011stat±0.02sys
//                     20 - 40
//                     1.2998±0.0097stat±0.0186sys
//                     40 - 60
//                     1.2109±0.0094stat±0.0173sys
//                     60 - 80
//                     1.108±0.0093stat±0.0214sys
//                     80 - 100
//                     0.9433±0.0092stat±0.0164sys
//                     NSD
//                     1.2698±0.0051stat±0.017sys
                    TString k0StarHEPDataFile                                   = "";
                    if (i == 0 || i == 16)
                        k0StarHEPDataFile                                       = "pPb/5_TeV/HEPData/K0Star_pPb_NSD.csv";
                    else if (i == 5)
                        k0StarHEPDataFile                                       = "pPb/5_TeV/HEPData/K0Star_pPb_0020.csv";
                    else if (i == 6)
                        k0StarHEPDataFile                                       = "pPb/5_TeV/HEPData/K0Star_pPb_2040.csv";
                    else if (i == 8)
                        k0StarHEPDataFile                                       = "pPb/5_TeV/HEPData/K0Star_pPb_4060.csv";
                    else if (i == 9)
                        k0StarHEPDataFile                                       = "pPb/5_TeV/HEPData/K0Star_pPb_6080.csv";
                    else if (i == 11)
                        k0StarHEPDataFile                                       = "pPb/5_TeV/HEPData/K0Star_pPb_80100.csv";

                    if (k0StarHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphK0StarYieldpPbStat              = ParseHEPData(k0StarHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphK0StarYieldpPbSys               = ParseHEPData(k0StarHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphK0StarYieldpPbStat                             = ScaleGraph(graphK0StarYieldpPbStat, scaleFacNSDToMB);
                            graphK0StarYieldpPbSys                              = ScaleGraph(graphK0StarYieldpPbSys, scaleFacNSDToMB);
                        }
                        SetGraphProperties(graphK0StarYieldpPbStat, Form("%sStat", fParticle[10].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphK0StarYieldpPbSys, Form("%sSys", fParticle[10].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphK0StarYieldpPbStat);
                        lists5TeV[i]->Add(graphK0StarYieldpPbSys);
                    }
                }

                //================================================================================================================
                // reading and writing rho^0 to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[11]) {
                    cout << " - rho^0" << endl;

                }

                //================================================================================================================
                // reading and writing rho^+/- to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[12]) {
                    cout << " - rho^+/-" << endl;

                }

                //================================================================================================================
                // reading and writing Delta^0 to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[13]) {
                    cout << " - Delta^0" << endl;

                }

                //================================================================================================================
                // reading and writing Delta^+/- to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[14]) {
                    cout << " - Delta^+/-" << endl;

                }

                //================================================================================================================
                // reading and writing K^0_s to 5TeV list
                // input is given as yield (dN/dydpT)
                //================================================================================================================
                if (includeParticle[15]) {
                    // spectra
                    if (fileK0s5TeV->GetListOfKeys()->Contains(Form("fHistPtK0Short_%s_OnlyStat", fCentralityOpt2[i].Data())) && fileK0s5TeV->GetListOfKeys()->Contains(Form("fHistPtK0Short_%s_OnlySyst", fCentralityOpt2[i].Data()))) {
                        cout << " - K^0_s spectrum" << endl;

                        TH1D* histoK0s5TeVStat               = (TH1D*)fileK0s5TeV->Get(Form("fHistPtK0Short_%s_OnlyStat", fCentralityOpt2[i].Data()));
                        TH1D* histoK0s5TeVSys                = (TH1D*)fileK0s5TeV->Get(Form("fHistPtK0Short_%s_OnlySyst", fCentralityOpt2[i].Data()));

                        SetHistoProperties(histoK0s5TeVStat, Form("%sStat", fParticle[15].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoK0s5TeVSys, Form("%sSys", fParticle[15].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        if (i == 0 ){
                            // Scale to MBAND
                            histoK0s5TeVStat->Scale(scaleFacNSDToMB);
                            histoK0s5TeVSys->Scale(scaleFacNSDToMB);
                        }

                        lists5TeV[i]->Add(histoK0s5TeVStat);
                        lists5TeV[i]->Add(histoK0s5TeVSys);

                        // calculate K0s/K ratio
                        if ( lists5TeV[i]->FindObject(Form("%sStat", fParticle[7].Data())) ) {
                            TH1D* histoCKaonStat                                = (TH1D*)lists5TeV[i]->FindObject("CKaonStat");
                            TH1D* histoCKaonSys                                 = (TH1D*)lists5TeV[i]->FindObject("CKaonSys");

                            TGraphErrors* graphCKaonStatRebinnedNKaonCKaon      = NULL;
                            TGraphErrors* graphCKaonSysRebinnedNKaonCKaon       = NULL;
                            TGraphErrors* graphNKaonStatRebinnedNKaonCKaon      = NULL;
                            TGraphErrors* graphNKaonSysRebinnedNKaonCKaon       = NULL;
                            TGraphErrors* graphNKaonCKaonStatErr                = NULL;
                            TGraphErrors* graphNKaonCKaonSysErr                 = NULL;
                            TGraphErrors* graphNKaonCKaonTotErr                 = CalculateRatioBetweenSpectraWithDifferentBinning(       histoK0s5TeVStat, histoK0s5TeVSys,
                                                                                                                                          histoCKaonStat, histoCKaonSys,
                                                                                                                                          kFALSE,  kFALSE,
                                                                                                                                          &graphNKaonStatRebinnedNKaonCKaon, &graphNKaonSysRebinnedNKaonCKaon,
                                                                                                                                          &graphCKaonStatRebinnedNKaonCKaon, &graphCKaonSysRebinnedNKaonCKaon,
                                                                                                                                          kTRUE, &graphNKaonCKaonStatErr, &graphNKaonCKaonSysErr
                            )    ;

                            SetGraphProperties(graphNKaonCKaonTotErr,  Form("%sTo%s%sTot", fParticle[15].Data(), fParticle[6].Data(),""), "#it{p}_{T} (GeV/#it{c})",   "#K^{0}_{s} / (K^{+}+K^{-})/2", "");
                            SetGraphProperties(graphNKaonCKaonStatErr,  Form("%sTo%s%sStat", fParticle[15].Data(), fParticle[6].Data(),""), "#it{p}_{T} (GeV/#it{c})",  "#K^{0}_{s} / (K^{+}+K^{-})/2", "");
                            SetGraphProperties(graphNKaonCKaonSysErr,  Form("%sTo%s%sSys", fParticle[15].Data(), fParticle[6].Data(),""), "#it{p}_{T} (GeV/#it{c})",  "#K^{0}_{s} / (K^{+}+K^{-})/2", "");
                            lists5TeV[i]->Add(graphNKaonCKaonTotErr);
                            lists5TeV[i]->Add(graphNKaonCKaonStatErr);
                            lists5TeV[i]->Add(graphNKaonCKaonSysErr);

                            TGraphErrors* graphCKaonStatRebinnedCKaonNKaon      = NULL;
                            TGraphErrors* graphCKaonSysRebinnedCKaonNKaon       = NULL;
                            TGraphErrors* graphNKaonStatRebinnedCKaonNKaon      = NULL;
                            TGraphErrors* graphNKaonSysRebinnedCKaonNKaon       = NULL;
                            TGraphErrors* graphCKaonNKaonStatErr                = NULL;
                            TGraphErrors* graphCKaonNKaonSysErr                 = NULL;
                            TGraphErrors* graphCKaonNKaonTotErr                 = CalculateRatioBetweenSpectraWithDifferentBinning(       histoCKaonStat, histoCKaonSys,
                                                                                                                                          histoK0s5TeVStat, histoK0s5TeVSys,
                                                                                                                                          kFALSE,  kFALSE,
                                                                                                                                          &graphNKaonStatRebinnedCKaonNKaon, &graphNKaonSysRebinnedCKaonNKaon,
                                                                                                                                          &graphCKaonStatRebinnedCKaonNKaon, &graphCKaonSysRebinnedCKaonNKaon,
                                                                                                                                          kTRUE, &graphCKaonNKaonStatErr, &graphCKaonNKaonSysErr
                            )    ;

                            SetGraphProperties(graphCKaonNKaonTotErr,  Form("%sTo%s%sTot", fParticle[6].Data(), fParticle[15].Data(),""), "#it{p}_{T} (GeV/#it{c})", "(K^{+}+K^{-})/2 / #K^{0}_{s}", "");
                            SetGraphProperties(graphCKaonNKaonStatErr,  Form("%sTo%s%sStat", fParticle[6].Data(), fParticle[15].Data(),""), "#it{p}_{T} (GeV/#it{c})", "(K^{+}+K^{-})/2 / #K^{0}_{s}", "");
                            SetGraphProperties(graphCKaonNKaonSysErr,  Form("%sTo%s%sSys", fParticle[6].Data(), fParticle[15].Data(),""), "#it{p}_{T} (GeV/#it{c})", "(K^{+}+K^{-})/2 / #K^{0}_{s}", "");
                            lists5TeV[i]->Add(graphCKaonNKaonTotErr);
                            lists5TeV[i]->Add(graphCKaonNKaonStatErr);
                            lists5TeV[i]->Add(graphCKaonNKaonSysErr);


                        }
                    }
                }

                //================================================================================================================
                // reading and writing Lambda to 5TeV list
                // input is given as yield (dN/dydpT)
                //================================================================================================================
                if (includeParticle[16]) {
                    // spectra
                    if (fileLambda5TeV->GetListOfKeys()->Contains(Form("fHistPtLambda_%s_OnlyStat", fCentralityOpt2[i].Data())) && fileLambda5TeV->GetListOfKeys()->Contains(Form("fHistPtLambda_%s_OnlySyst", fCentralityOpt2[i].Data())) && fileAntiLambda5TeV->GetListOfKeys()->Contains(Form("fHistPtAntiLambda_%s_OnlyStat", fCentralityOpt2[i].Data())) && fileAntiLambda5TeV->GetListOfKeys()->Contains(Form("fHistPtAntiLambda_%s_OnlySyst", fCentralityOpt2[i].Data()))) {
                        cout << " - Lambda spectrum" << endl;

                        TH1D* histoLambda5TeVStat                   = (TH1D*)fileLambda5TeV->Get(Form("fHistPtLambda_%s_OnlyStat", fCentralityOpt2[i].Data()));
                        TH1D* histoLambda5TeVSys                    = (TH1D*)fileLambda5TeV->Get(Form("fHistPtLambda_%s_OnlySyst", fCentralityOpt2[i].Data()));
                        TH1D* histoAntiLambda5TeVStat               = (TH1D*)fileAntiLambda5TeV->Get(Form("fHistPtAntiLambda_%s_OnlyStat", fCentralityOpt2[i].Data()));
                        TH1D* histoAntiLambda5TeVSys                = (TH1D*)fileAntiLambda5TeV->Get(Form("fHistPtAntiLambda_%s_OnlySyst", fCentralityOpt2[i].Data()));

                        histoLambda5TeVStat->Add(histoAntiLambda5TeVStat);
                        histoLambda5TeVSys->Add(histoAntiLambda5TeVSys);

                        histoLambda5TeVStat->Scale(0.5);
                        histoLambda5TeVSys->Scale(0.5);

                        SetHistoProperties(histoLambda5TeVStat, Form("%sStat", fParticle[16].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoLambda5TeVSys, Form("%sSys", fParticle[16].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        if (i == 0){
                            // Scale to MBAND
                            histoLambda5TeVStat->Scale(scaleFacNSDToMB);
                            histoLambda5TeVSys->Scale(scaleFacNSDToMB);
                        }

                        lists5TeV[i]->Add(histoLambda5TeVStat);
                        lists5TeV[i]->Add(histoLambda5TeVSys);

                        // calculate Lambda/p ratio
                        if ( lists5TeV[i]->FindObject(Form("%sStat", fParticle[7].Data())) ) {
                            TH1D* histoProtonStat                                    = (TH1D*)lists5TeV[i]->FindObject("ProtonStat");
                            TH1D* histoProtonSys                                     = (TH1D*)lists5TeV[i]->FindObject("ProtonSys");

                            TGraphErrors* graphProtonStatRebinnedLambdaProton        = NULL;
                            TGraphErrors* graphProtonSysRebinnedLambdaProton         = NULL;
                            TGraphErrors* graphLambdaStatRebinnedLambdaProton        = NULL;
                            TGraphErrors* graphLambdaSysRebinnedLambdaProton         = NULL;
                            TGraphErrors* graphLambdaProtonStatErr                   = NULL;
                            TGraphErrors* graphLambdaProtonSysErr                    = NULL;
                            TGraphErrors* graphLambdaProtonTotErr                    = CalculateRatioBetweenSpectraWithDifferentBinning(   histoLambda5TeVStat, histoLambda5TeVSys,
                                                                                                                                           histoProtonStat, histoProtonSys,
                                                                                                                                           kFALSE,  kFALSE,
                                                                                                                                           &graphLambdaStatRebinnedLambdaProton, &graphLambdaSysRebinnedLambdaProton,
                                                                                                                                           &graphProtonStatRebinnedLambdaProton, &graphProtonSysRebinnedLambdaProton,
                                                                                                                                           kTRUE, &graphLambdaProtonStatErr, &graphLambdaProtonSysErr
                            )    ;

                            SetGraphProperties(graphLambdaProtonTotErr,  Form("%sTo%s%sTot", fParticle[16].Data(), fParticle[7].Data(),""), "#it{p}_{T} (GeV/#it{c})",   "#Lambda / p", "");
                            SetGraphProperties(graphLambdaProtonStatErr,  Form("%sTo%s%sStat", fParticle[16].Data(), fParticle[7].Data(),""), "#it{p}_{T} (GeV/#it{c})",  "#Lambda / p", "");
                            SetGraphProperties(graphLambdaProtonSysErr,  Form("%sTo%s%sSys", fParticle[16].Data(), fParticle[7].Data(),""), "#it{p}_{T} (GeV/#it{c})",  "#Lambda / p", "");
                            lists5TeV[i]->Add(graphLambdaProtonTotErr);
                            lists5TeV[i]->Add(graphLambdaProtonStatErr);
                            lists5TeV[i]->Add(graphLambdaProtonSysErr);
                        }
                    }

                    // ratios
                    if (fileLambda5TeVRatios->GetListOfKeys()->Contains(Form("fHistPtLamOverKaon_%s_OnlyStat", fCentralityOpt2[i].Data())) && fileLambda5TeVRatios->GetListOfKeys()->Contains(Form("fHistPtLamOverKaon_%s_OnlySyst", fCentralityOpt2[i].Data()))) {
                        cout << " - Lambda ratio" << endl;

                        TH1D* histoLambdaToKaon5TeVStat             = (TH1D*)fileLambda5TeVRatios->Get(Form("fHistPtLamOverKaon_%s_OnlyStat", fCentralityOpt2[i].Data()));
                        TH1D* histoLambdaToKaon5TeVSys              = (TH1D*)fileLambda5TeVRatios->Get(Form("fHistPtLamOverKaon_%s_OnlySyst", fCentralityOpt2[i].Data()));

                        SetHistoProperties(histoLambdaToKaon5TeVStat, Form("%sTo%sStat", fParticle[16].Data(), fParticle[15].Data()), "#it{p}_{T} (GeV/#it{c})", "#Lambda / K^{0}_{s}", "");
                        SetHistoProperties(histoLambdaToKaon5TeVSys, Form("%sTo%sSys", fParticle[16].Data(), fParticle[15].Data()), "#it{p}_{T} (GeV/#it{c})", "#Lambda / K^{0}_{s}", "");

                        lists5TeV[i]->Add(histoLambdaToKaon5TeVStat);
                        lists5TeV[i]->Add(histoLambdaToKaon5TeVSys);
                    }

                }


                //================================================================================================================
                // reading and writing Sigma^0 to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[17]) {
                    cout << " - Sigma^0" << endl;

                }

                //================================================================================================================
                // reading and writing Sigma^+/- to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[18]) {
                    cout << " - Sigma^+/-" << endl;

                }

                //================================================================================================================
                // reading and writing Omega^+/- to 5TeV list
                // input is given as fully invariant yield
                //================================================================================================================
                if (includeParticle[19]) {
                    // spectra
                    //      MB, 0-5,    5-10,   0-10,   10-20,  0-20,   20-40,  20-50,  40-60,  60-80,  80-100
                    //      0,  1,      2,      3,      4,      5,      6,      7,      8,      9,      10
                    TString OmegaHEPDataFile                            = "";
                    if (i == 1){
                        //   Omega- mean Pt = 1.818 0.056 -0.056 0.076 -0.076 0.025 -0.025
                        //   Omega+ mean Pt = 1.761 0.06 -0.06 0.075 -0.075 0.025 -0.025
                        //                     0–5%	45 ± 1	0.2354 ± 0.0020 ± 0.0161
                        OmegaHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Omega_pPb_0005.csv";
                        bMeanPtPub[0][i]                                = kTRUE;
                        bIntegYieldPub[0][i]                            = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(19+1, (1.818+1.761)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(19+1, (1.818+1.761)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(19+1, (0.056+0.06)/2 );
                        histoMeanPtPerParticleSys[0][i]->SetBinError(19+1, (0.076+0.075)/2);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(19+1, 0.2354);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(19+1, 0.2354);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(19+1, 0.0020);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(19+1, 0.0161);
                    } else if (i == 2){
                        //   Omega- mean Pt = 1.829 0.062 -0.062 0.079 -0.079 0.03 -0.03
                        //   Omega+ mean Pt = 1.847 0.062 -0.062 0.086 -0.086 0.032 -0.032
                        //                     5–10%	36.2 ± 0.8	0.1861 ± 0.0016 ± 0.0138
                        OmegaHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Omega_pPb_0510.csv";
                        bMeanPtPub[0][i]                                = kTRUE;
                        bIntegYieldPub[0][i]                            = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(19+1, (1.829+1.847)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(19+1, (1.829+1.847)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(19+1, (0.062+0.062)/2 );
                        histoMeanPtPerParticleSys[0][i]->SetBinError(19+1, (0.079+0.086)/2);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(19+1, 0.1861);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(19+1, 0.1861);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(19+1, 0.0016);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(19+1, 0.0138);
                    } else if (i == 4){
                        //   Omega- mean Pt = 1.771 0.058 -0.058 0.082 -0.082 0.038 -0.038
                        //   Omega+ mean Pt = 1.739 0.048 -0.048 0.073 -0.073 0.034 -0.034
                        //                     10–20%	30.5 ± 0.7	0.1500 ± 0.0010 ± 0.0112
                        OmegaHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Omega_pPb_1020.csv";
                        bMeanPtPub[0][i]                                = kTRUE;
                        bIntegYieldPub[0][i]                            = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(19+1, (1.771+1.739)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(19+1, (1.771+1.739)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(19+1, (0.058+0.048)/2 );
                        histoMeanPtPerParticleSys[0][i]->SetBinError(19+1, (0.082+0.073)/2);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(19+1, 0.1500);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(19+1, 0.1500);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(19+1, 0.0010);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(19+1, 0.0112);
                    } else if (i == 6){
                        //   Omega- mean Pt = 1.668 0.04 -0.04 0.082 -0.082 0.029 -0.029
                        //   Omega+ mean Pt = 1.707 0.036 -0.036 0.061 -0.061 0.022 -0.022
                        //                     20–40%	23.2 ± 0.5	0.1100 ± 0.0006 ± 0.0085
                        OmegaHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Omega_pPb_2040.csv";
                        bMeanPtPub[0][i]                                = kTRUE;
                        bIntegYieldPub[0][i]                            = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(19+1, (1.668+1.707)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(19+1, (1.668+1.707)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(19+1, (0.04+0.036)/2 );
                        histoMeanPtPerParticleSys[0][i]->SetBinError(19+1, (0.082+0.061)/2);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(19+1, 0.1100);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(19+1, 0.1100);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(19+1, 0.0006);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(19+1, 0.0085);
                    } else if (i == 8){
                        //   Omega- mean Pt = 1.563 0.033 -0.033 0.067 -0.067 0.035 -0.035
                        //   Omega+ mean Pt = 1.593 0.04 -0.04 0.058 -0.058 0.031 -0.031
                        //                     40–60%	16.1 ± 0.4	0.0726 ± 0.0006 ± 0.0065
                        OmegaHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Omega_pPb_4060.csv";
                        bMeanPtPub[0][i]                                = kTRUE;
                        bIntegYieldPub[0][i]                            = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(19+1, (1.563+1.593)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(19+1, (1.563+1.593)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(19+1, (0.033+0.04)/2 );
                        histoMeanPtPerParticleSys[0][i]->SetBinError(19+1, (0.067+0.058)/2);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(19+1, 0.0726);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(19+1, 0.0726);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(19+1, 0.0006);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(19+1, 0.0065);
                    } else if (i == 9){
                        //   Omega- mean Pt = 1.544 0.055 -0.055 0.084 -0.084 0.04 -0.04
                        //   Omega+ mean Pt = 1.382 0.078 -0.078 0.076 -0.076 0.036 -0.036
                        //                     60–80%	9.8 ± 0.24	0.0398 ± 0.0004 ± 0.0031
                        OmegaHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Omega_pPb_6080.csv";
                        bMeanPtPub[0][i]                                = kTRUE;
                        bIntegYieldPub[0][i]                            = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(19+1, (1.544+1.382)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(19+1, (1.544+1.382)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(19+1, (0.055+0.078)/2 );
                        histoMeanPtPerParticleSys[0][i]->SetBinError(19+1, (0.084+0.076)/2);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(19+1, 0.0398);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(19+1, 0.0398);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(19+1, 0.0004);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(19+1, 0.0031);
                    } else if (i == 11){
                        //   Omega- mean Pt = 1.35 0.11 -0.11 0.106 -0.106 0.048 -0.048
                        //   Omega+ mean Pt = 1.187 0.097 -0.097 0.096 -0.096 0.043 -0.043
                        //                     80–100%	4.3 ± 0.1	0.0143 ± 0.0003 ± 0.0015
                        OmegaHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Omega_pPb_80100.csv";
                        bMeanPtPub[0][i]                                = kTRUE;
                        bIntegYieldPub[0][i]                            = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(19+1, (1.35+1.187)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(19+1, (1.35+1.187)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(19+1, (0.11+0.097)/2 );
                        histoMeanPtPerParticleSys[0][i]->SetBinError(19+1, (0.106+0.096)/2);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(19+1, 0.0143);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(19+1, 0.0143);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(19+1, 0.0003);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(19+1, 0.0015);
                    }

                    if (OmegaHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphOmegaYieldpPbStat       = ParseHEPData(OmegaHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphOmegaYieldpPbSys        = ParseHEPData(OmegaHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphOmegaYieldpPbStat                      = ScaleGraph(graphOmegaYieldpPbStat, scaleFacNSDToMB);
                            graphOmegaYieldpPbSys                       = ScaleGraph(graphOmegaYieldpPbSys, scaleFacNSDToMB);
                        }
                        graphOmegaYieldpPbStat                          = ConvertYieldGraph(graphOmegaYieldpPbStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        graphOmegaYieldpPbSys                           = ConvertYieldGraph(graphOmegaYieldpPbSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        SetGraphProperties(graphOmegaYieldpPbStat, Form("%sStat", fParticle[19].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphOmegaYieldpPbSys, Form("%sSys", fParticle[19].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphOmegaYieldpPbStat);
                        lists5TeV[i]->Add(graphOmegaYieldpPbSys);
                    }
                }

                //================================================================================================================
                // reading and writing Xi^+/- to 5TeV list
                // input is given as fully invariant yield
                //================================================================================================================
                if (includeParticle[20]) {
                    // spectra
                    // spectra
                    //      MB, 0-5,    5-10,   0-10,   10-20,  0-20,   20-40,  20-50,  40-60,  60-80,  80-100
                    //      0,  1,      2,      3,      4,      5,      6,      7,      8,      9,      10
//                     Yields
//                     0–5%	45 ± 1	0.0260 ± 0.0011 ± 0.0034
//                     5–10%	36.2 ± 0.8	0.0215 ± 0.0008 ± 0.0029
//                     10–20%	30.5 ± 0.7	0.0167 ± 0.0006 ± 0.0022
//                     20–40%	23.2 ± 0.5	0.0120 ± 0.0005 ± 0.0016
//                     40–60%	16.1 ± 0.4	0.0072 ± 0.0003 ± 0.0010
//                     60–80%	9.8 ± 0.24	0.0042 ± 0.0002 ± 0.0006
//                     80–100%	4.3 ± 0.1	0.0013 ± 0.0003 ± 0.0003

                    TString XiHEPDataFile                            = "";
                    if (i == 1){
                        XiHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Xi_pPb_0005.csv";
                        // Xi+ meanPt = 1.574 0.01 -0.01 0.037 -0.037 0.012 -0.012
                        // Xi- meanPt = 1.549 0.01 -0.01 0.037 -0.037 0.012 -0.012
                        bMeanPtPub[0][i]                             = kTRUE;
                        bIntegYieldPub[0][i]                         = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(20+1, (1.574+1.549)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(20+1, 0.01 );
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(20+1, (1.574+1.549)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinError(20+1, 0.037);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(20+1, 0.0260);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(20+1, 0.0011);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(20+1, 0.0260);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(20+1, 0.0034);

                    } else if (i == 2){
                        XiHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Xi_pPb_0510.csv";
                        // Xi+ meanPt = 1.529 0.009 -0.009 0.036 -0.036 0.014 -0.014
                        // Xi- meanPt = 1.516 0.009 -0.009 0.036 -0.036 0.014 -0.014
                        bMeanPtPub[0][i]                             = kTRUE;
                        bIntegYieldPub[0][i]                         = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(20+1, (1.529+1.516)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(20+1, 0.009 );
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(20+1, (1.529+1.516)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinError(20+1, 0.036);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(20+1, 0.0215);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(20+1, 0.0008);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(20+1, 0.0215);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(20+1, 0.0029);
                    } else if (i == 4){
                        XiHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Xi_pPb_1020.csv";
                        // Xi+ meanPt = 1.515 0.007 -0.007 0.034 -0.034 0.016 -0.016
                        // Xi- meanPt = 1.502 0.007 -0.007 0.034 -0.034 0.016 -0.016
                        bMeanPtPub[0][i]                             = kTRUE;
                        bIntegYieldPub[0][i]                         = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(20+1, (1.515+1.502)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(20+1, 0.007 );
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(20+1, (1.515+1.502)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinError(20+1, 0.034);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(20+1, 0.0167);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(20+1, 0.0006);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(20+1, 0.0167);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(20+1, 0.0022);
                    } else if (i == 6){
                        XiHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Xi_pPb_2040.csv";
                        // Xi+ meanPt = 1.468 0.006 -0.006 0.032 -0.032 0.011 -0.011
                        // Xi- meanPt = 1.453 0.006 -0.006 0.032 -0.032 0.011 -0.011
                        bMeanPtPub[0][i]                             = kTRUE;
                        bIntegYieldPub[0][i]                         = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(20+1, (1.468+1.453)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(20+1, 0.006 );
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(20+1, (1.468+1.453)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinError(20+1, 0.032);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(20+1, 0.0120);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(20+1, 0.0005);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(20+1, 0.0120);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(20+1, 0.0016);
                    } else if (i == 8){
                        XiHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Xi_pPb_4060.csv";
                        // Xi+ meanPt = 1.394 0.012 -0.012 0.03 -0.03 0.016 -0.016
                        // Xi- meanPt = 1.358 0.007 -0.007 0.029 -0.029 0.015 -0.015
                        bMeanPtPub[0][i]                             = kTRUE;
                        bIntegYieldPub[0][i]                         = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(20+1, (1.394+1.358)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(20+1, (0.012+0.007)/2. );
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(20+1, (1.394+1.358)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinError(20+1, (0.03+0.029)/2.);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(20+1, 0.0072);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(20+1, 0.0003);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(20+1, 0.0072);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(20+1, 0.0010);
                    } else if (i == 9){
                        XiHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Xi_pPb_6080.csv";
                        // Xi+ meanPt = 1.28 0.009 -0.009 0.032 -0.032 0.015 -0.015
                        // Xi- meanPt = 1.27 0.009 -0.009 0.041 -0.041 0.019 -0.019
                        bMeanPtPub[0][i]                             = kTRUE;
                        bIntegYieldPub[0][i]                         = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(20+1, (1.28+1.27)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(20+1, 0.009 );
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(20+1, (1.28+1.27)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinError(20+1, (0.032+0.041)/2.);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(20+1, 0.0042);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(20+1, 0.0002);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(20+1, 0.0042);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(20+1, 0.0006);

                    } else if (i == 11){
                        XiHEPDataFile                                = "pPb/5_TeV/HEPData/XiAndOmega_HEPData-ins1411084-v1-csv/Xi_pPb_80100.csv";
                        // Xi+ meanPt = 1.155 0.018 -0.018 0.037 -0.037 0.017 -0.017
                        // Xi- meanPt = 1.118 0.017 -0.017 0.04 -0.04 0.018 -0.018
                        bMeanPtPub[0][i]                             = kTRUE;
                        bIntegYieldPub[0][i]                         = kTRUE;
                        histoMeanPtPerParticleStat[0][i]->SetBinContent(20+1, (1.155+1.118)/2.);
                        histoMeanPtPerParticleStat[0][i]->SetBinError(20+1, 0.0175 );
                        histoMeanPtPerParticleSys[0][i]->SetBinContent(20+1, (1.155+1.118)/2.);
                        histoMeanPtPerParticleSys[0][i]->SetBinError(20+1, (0.037+0.04)/2.);
                        histoIntegYieldPerParticleStat[0][i]->SetBinContent(20+1, 0.0013);
                        histoIntegYieldPerParticleStat[0][i]->SetBinError(20+1, 0.0003);
                        histoIntegYieldPerParticleSys[0][i]->SetBinContent(20+1, 0.0013);
                        histoIntegYieldPerParticleSys[0][i]->SetBinError(20+1, 0.0003);
                    }

                    if (XiHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphXiYieldpPbStat       = ParseHEPData(XiHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphXiYieldpPbSys        = ParseHEPData(XiHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphXiYieldpPbStat                      = ScaleGraph(graphXiYieldpPbStat, scaleFacNSDToMB);
                            graphXiYieldpPbSys                       = ScaleGraph(graphXiYieldpPbSys, scaleFacNSDToMB);
                        }
                        graphXiYieldpPbStat                          = ConvertYieldGraph(graphXiYieldpPbStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        graphXiYieldpPbSys                           = ConvertYieldGraph(graphXiYieldpPbSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        SetGraphProperties(graphXiYieldpPbStat, Form("%sStat", fParticle[20].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphXiYieldpPbSys, Form("%sSys", fParticle[20].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphXiYieldpPbStat);
                        lists5TeV[i]->Add(graphXiYieldpPbSys);
                    }


                }

                //================================================================================================================
                // reading and writing J/psi to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[21]) {

                    cout << " - J/psi" << endl;

                }

                //================================================================================================================
                // reading and writing D0 to 5TeV list
                // input is given as xsection
                //================================================================================================================

                if (includeParticle[22]) {

                    cout << " - D0" << endl;
                    TString DZeroHEPDataFile                            = "";
                    if (i == 0 || i == 16)
                        DZeroHEPDataFile                                = "pPb/5_TeV/HEPData/DZero_pPb_NSD_HEPData-ins1465513-v1-Table5.csv";

                    if (DZeroHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphDZeroYieldpPbStat       = ParseHEPData(DZeroHEPDataFile, 12, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphDZeroYieldpPbSys        = ParseHEPData(DZeroHEPDataFile, 12, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphDZeroYieldpPbStat                      = ScaleGraph(graphDZeroYieldpPbStat, 1.e3/xSection5TeVNSD*scaleFacNSDToMB);
                            graphDZeroYieldpPbSys                       = ScaleGraph(graphDZeroYieldpPbSys, 1.e3/xSection5TeVNSD*scaleFacNSDToMB);
                        } else if ( i == 16){
                            // Scale to NSD yield
                            graphDZeroYieldpPbStat                      = ScaleGraph(graphDZeroYieldpPbStat, 1.e3/xSection5TeVNSD);
                            graphDZeroYieldpPbSys                       = ScaleGraph(graphDZeroYieldpPbSys, 1.e3/xSection5TeVNSD);
                        }
                        SetGraphProperties(graphDZeroYieldpPbStat, Form("%sStat", fParticle[22].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphDZeroYieldpPbSys, Form("%sSys", fParticle[22].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphDZeroYieldpPbStat);
                        lists5TeV[i]->Add(graphDZeroYieldpPbSys);
                    }
                    //

                }
                //================================================================================================================
                // reading and writing D+ to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[23]) {

                    cout << " - D+" << endl;
                    TString DHEPDataFile                            = "";
                    if (i == 0 || i == 16)
                        DHEPDataFile                                = "pPb/5_TeV/HEPData/DMesons_HEPData-ins1296081-v1-csv/DPlus_pPb_NSD.csv";

                    if (DHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphDYieldpPbStat       = ParseHEPData(DHEPDataFile, 12, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphDYieldpPbSys        = ParseHEPData(DHEPDataFile, 12, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphDYieldpPbStat                      = ScaleGraph(graphDYieldpPbStat, 1.e3/xSection5TeVNSD*scaleFacNSDToMB);
                            graphDYieldpPbSys                       = ScaleGraph(graphDYieldpPbSys, 1.e3/xSection5TeVNSD*scaleFacNSDToMB);
                        } else if ( i == 16){
                            // Scale to NSD yield
                            graphDYieldpPbStat                      = ScaleGraph(graphDYieldpPbStat, 1.e3/xSection5TeVNSD);
                            graphDYieldpPbSys                       = ScaleGraph(graphDYieldpPbSys, 1.e3/xSection5TeVNSD);
                        }
                        SetGraphProperties(graphDYieldpPbStat, Form("%sStat", fParticle[23].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphDYieldpPbSys, Form("%sSys", fParticle[23].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphDYieldpPbStat);
                        lists5TeV[i]->Add(graphDYieldpPbSys);
                    }

                }
                //================================================================================================================
                // reading and writing D* to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[24]) {

                    cout << " - D*" << endl;
//
                    TString DHEPDataFile                            = "";
                    if (i == 0 || i == 16)
                        DHEPDataFile                                = "pPb/5_TeV/HEPData/DMesons_HEPData-ins1296081-v1-csv/DStarPlus_pPb_NSD.csv";

                    if (DHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphDYieldpPbStat       = ParseHEPData(DHEPDataFile, 12, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphDYieldpPbSys        = ParseHEPData(DHEPDataFile, 12, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphDYieldpPbStat                      = ScaleGraph(graphDYieldpPbStat, 1.e3/xSection5TeVNSD*scaleFacNSDToMB);
                            graphDYieldpPbSys                       = ScaleGraph(graphDYieldpPbSys, 1.e3/xSection5TeVNSD*scaleFacNSDToMB);
                        } else if (i == 16){
                            // Scale to NSD yield
                            graphDYieldpPbStat                      = ScaleGraph(graphDYieldpPbStat, 1.e3/xSection5TeVNSD);
                            graphDYieldpPbSys                       = ScaleGraph(graphDYieldpPbSys, 1.e3/xSection5TeVNSD);
                        }

                        SetGraphProperties(graphDYieldpPbStat, Form("%sStat", fParticle[24].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphDYieldpPbSys, Form("%sSys", fParticle[24].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphDYieldpPbStat);
                        lists5TeV[i]->Add(graphDYieldpPbSys);
                    }

                }
                //================================================================================================================
                // reading and writing Ds+ to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[25]) {

                    cout << " - Ds+" << endl;
                    TString DHEPDataFile                            = "";
                    if (i == 0 || i == 16)
                        DHEPDataFile                                = "pPb/5_TeV/HEPData/DMesons_HEPData-ins1296081-v1-csv/DSPlus_pPb_NSD.csv";

                    if (DHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphDYieldpPbStat       = ParseHEPData(DHEPDataFile, 12, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphDYieldpPbSys        = ParseHEPData(DHEPDataFile, 12, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphDYieldpPbStat                      = ScaleGraph(graphDYieldpPbStat, 1.e3/xSection5TeVNSD*scaleFacNSDToMB);
                            graphDYieldpPbSys                       = ScaleGraph(graphDYieldpPbSys, 1.e3/xSection5TeVNSD*scaleFacNSDToMB);
                        } else if (i == 16){
                            // Scale to NSD yield
                            graphDYieldpPbStat                      = ScaleGraph(graphDYieldpPbStat, 1.e3/xSection5TeVNSD);
                            graphDYieldpPbSys                       = ScaleGraph(graphDYieldpPbSys, 1.e3/xSection5TeVNSD);
                        }
                        SetGraphProperties(graphDYieldpPbStat, Form("%sStat", fParticle[25].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphDYieldpPbSys, Form("%sSys", fParticle[25].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphDYieldpPbStat);
                        lists5TeV[i]->Add(graphDYieldpPbSys);
                    }
                }
                //================================================================================================================
                // reading and writing Sigma*+ to 5TeV list
                // input is given as yield (dN/dydpT)
                //================================================================================================================
                if (includeParticle[26]) {
                    cout << " - Sigma*+" << endl;
                    // spectra
                    //      MB, 0-5,    5-10,   0-10,   10-20,  0-20,   20-40,  20-50,  40-60,  60-80,  80-100, 20-30,  30-40,  20-60   60-100
                    //      0,  1,      2,      3,      4,      5,      6,      7,      8,      9,      10,     11,     12,     13,     14
//                     Σ∗±Σ∗±
//                     NSD
//                     49.0 ± 0.6 ± 6.5
//                     1.367 ± 0.009 ± 0.061
//
//                     0–20%
//                     90.3 ± 1.4 ± 7.9
//                     1.495 ± 0.012 ± 0.046
//
//                     20–60%
//                     52.2 ± 0.8 ± 6.0
//                     1.342 ± 0.010 ± 0.055
//
//                     60–100%
//                     15.2 ± 0.4 ± 2.4
//                     1.173 ± 0.015 ± 0.067
                    TString SigmaStarPlusHEPDataFile                            = "";
                    if (i == 0 || i == 16)
                        SigmaStarPlusHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/SigmaStar_pPb_NSD.csv";
                    else if (i == 5)
                        SigmaStarPlusHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/SigmaStar_pPb_0020.csv";
                    else if (i == 14)
                        SigmaStarPlusHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/SigmaStar_pPb_2060.csv";
                    else if (i == 15)
                        SigmaStarPlusHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/SigmaStar_pPb_60100.csv";

                    if (SigmaStarPlusHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphSigmaStarPlusYieldpPbStat       = ParseHEPData(SigmaStarPlusHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphSigmaStarPlusYieldpPbSys        = ParseHEPData(SigmaStarPlusHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphSigmaStarPlusYieldpPbStat                      = ScaleGraph(graphSigmaStarPlusYieldpPbStat, scaleFacNSDToMB);
                            graphSigmaStarPlusYieldpPbSys                       = ScaleGraph(graphSigmaStarPlusYieldpPbSys, scaleFacNSDToMB);
                        }
                        SetGraphProperties(graphSigmaStarPlusYieldpPbStat, Form("%sStat", fParticle[26].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphSigmaStarPlusYieldpPbSys, Form("%sSys", fParticle[26].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphSigmaStarPlusYieldpPbStat);
                        lists5TeV[i]->Add(graphSigmaStarPlusYieldpPbSys);
                    }
                }
                //================================================================================================================
                // reading and writing Xi*0 to 5TeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[27]) {

                    cout << " - XiStar0" << endl;
                    // spectra
                    //      MB, 0-5,    5-10,   0-10,   10-20,  0-20,   20-40,  20-50,  40-60,  60-80,  80-100, 20-30,  30-40,  20-60   60-100
                    //      0,  1,      2,      3,      4,      5,      6,      7,      8,      9,      10,     11,     12,     13,     14
//                     1/2( Ξ∗0+Ξ⎯⎯⎯⎯0Ξ∗0+Ξ¯0 )
//                     NSD
//                     12.5 ± 0.3 ± 1.1
//                     1.540 ± 0.016 ± 0.071
//
//                     0–20%
//                     27.3 ± 0.6 ± 2.8
//                     1.626 ± 0.016 ± 0.068
//
//                     20–40%
//                     17.7 ± 0.5 ± 2.4
//                     1.482 ± 0.020 ± 0.100
//
//                     40–60%
//                     10.7 ± 0.3 ± 1.6
//                     1.459 ± 0.025 ± 0.114
//
//                     60–100%
//                     3.6 ± 0.1 ± 0.5
//                     1.377 ± 0.023 ± 0.089
                    TString XiStarZeroHEPDataFile                            = "";
                    if (i == 0 || i == 16)
                        XiStarZeroHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/XiStar_pPb_NSD.csv";
                    else if (i == 5)
                        XiStarZeroHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/XiStar_pPb_0020.csv";
                    else if (i == 6)
                        XiStarZeroHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/XiStar_pPb_2040.csv";
                    else if (i == 8)
                        XiStarZeroHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/XiStar_pPb_4060.csv";
                    else if (i == 15)
                        XiStarZeroHEPDataFile                                = "pPb/5_TeV/HEPData/XiStar0_SigmaStarPlus_HEPData-ins1510878-v1-csv/XiStar_pPb_60100.csv";

                    if (XiStarZeroHEPDataFile.CompareTo("") != 0){
                        cout << "   added spec" << endl;
                        TGraphAsymmErrors* graphXiStarZeroYieldpPbStat       = ParseHEPData(XiStarZeroHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
                        TGraphAsymmErrors* graphXiStarZeroYieldpPbSys        = ParseHEPData(XiStarZeroHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
                        if (i == 0){
                            // Scale to MBAND
                            graphXiStarZeroYieldpPbStat                      = ScaleGraph(graphXiStarZeroYieldpPbStat, scaleFacNSDToMB);
                            graphXiStarZeroYieldpPbSys                       = ScaleGraph(graphXiStarZeroYieldpPbSys, scaleFacNSDToMB);
                        }
                        SetGraphProperties(graphXiStarZeroYieldpPbStat, Form("%sStat", fParticle[27].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphXiStarZeroYieldpPbSys, Form("%sSys", fParticle[27].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", "");

                        lists5TeV[i]->Add(graphXiStarZeroYieldpPbStat);
                        lists5TeV[i]->Add(graphXiStarZeroYieldpPbSys);
                    }
                }
            }
        }
    }



    //================================================================================================================
    //Produce MB spectras and ratios
    //================================================================================================================
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[5], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[6], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[7], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[15], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[16], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[19], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[20], "", "");

    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[6], fParticle[5], "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[7], fParticle[5], "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[0], fParticle[16], fParticle[15], "");


    //================================================================================================================
    //Produce NSD spectras and ratios
    //================================================================================================================
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[16], fParticle[5], "", "", scaleFacMBToNSD);
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[16], fParticle[6], "", "", scaleFacMBToNSD);
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[16], fParticle[7], "", "", scaleFacMBToNSD);
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[16], fParticle[15], "", "", scaleFacMBToNSD);
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[16], fParticle[16], "", "", scaleFacMBToNSD);
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[16], fParticle[19], "", "", scaleFacMBToNSD);
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[16], fParticle[20], "", "", scaleFacMBToNSD);

//     Int_t centsToCombXiAndOmega[7]              = {1,2,4,6,8,9,11};
//     Double_t centsWeightsXiAndOmega[7]          = {0.05,0.05,0.1,0.2,0.2,0.2,0.2};
//
//     Bool_t bOmegaYield      = ProduceAveragedQuantiesInCentralityBin(  histoIntegYieldPerParticleStat[0], histoIntegYieldPerParticleSys[0], 7, centsToCombXiAndOmega,  centsWeightsXiAndOmega, 0, 19);
//     Bool_t bXiYield         = ProduceAveragedQuantiesInCentralityBin(  histoIntegYieldPerParticleStat[0], histoIntegYieldPerParticleSys[0], 7, centsToCombXiAndOmega,  centsWeightsXiAndOmega, 0, 20);
//     if (bOmegaYield || bXiYield)
//         bIntegYieldPub[0][0]                    = kTRUE;

    //================================================================================================================
    //Produce 0-10% spectra
    //================================================================================================================
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[5], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[6], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[7], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[9], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[15], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[16], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[19], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[20], "", "");

    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[6], fParticle[5], "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[7], fParticle[5], "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[3], fParticle[16], fParticle[15], "");

    //================================================================================================================
    //Produce 0-20% spectra
    //================================================================================================================
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[5], fParticle[5], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[5], fParticle[6], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[5], fParticle[7], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[5], fParticle[15], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[5], fParticle[16], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[5], fParticle[19], "", "");
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[5], fParticle[20], "", "");

    //================================================================================================================
    //Produce 60-100% spectra
    //================================================================================================================
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[5], "", ""); // ch pi
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[6], "", ""); // ch K
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[7], "", ""); // p
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[9], "", ""); // Phi
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[10], "", ""); // K0*
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[15], "", ""); // K0s
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[16], "", ""); // Lambda
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[20], "", ""); // Xi
    ProduceSpectrumInCentralityBin(lists5TeV, fCentrality[15], fParticle[19], "", ""); // COmega

    //================================================================================================================
    //Produce plots containing all particle spectra if specified
    //================================================================================================================
    if (produceAllSpectraPlots) {
        if (enable5TeV) {
            for (Int_t i=0; i<nCentralities; i++) {
                if (bMeanPtPub[0][i]){
                    histoMeanPtPerParticleStat[0][i]->SetName("meanPtPerParticle_pub_stat");
                    lists5TeV[i]->Add(histoMeanPtPerParticleStat[0][i]);
                    histoMeanPtPerParticleSys[0][i]->SetName("meanPtPerParticle_pub_sys");
                    lists5TeV[i]->Add(histoMeanPtPerParticleSys[0][i]);
                }
                if (bIntegYieldPub[0][i]){
                    histoIntegYieldPerParticleStat[0][i]->SetName("IntegYieldPerParticle_pub_stat");
                    lists5TeV[i]->Add(histoIntegYieldPerParticleStat[0][i]);
                    histoIntegYieldPerParticleSys[0][i]->SetName("IntegYieldPerParticle_pub_sys");
                    lists5TeV[i]->Add(histoIntegYieldPerParticleSys[0][i]);
                }
                if (includeCentrality5TeV[i]) {
                    ProduceParticleSpectraPlotFromList(lists5TeV[i], fCollSys[1], fEnergy[2], fCentrality[i], suffix);
                    ProduceParticleSpectraPlotFromListOnlyFinal(lists5TeV[i], fCollSys[1], fEnergy[2], fCentrality[i], suffix);
                    ProduceParticleRatioPlotFromList(lists5TeV[i], fCollSys[1], fEnergy[2], fCentrality[i], suffix);
                }
            }
        }

        if (enable8TeV) {
            for (Int_t i=0; i<nCentralities; i++) {
                if (includeCentrality8TeV[i]) {
                    ProduceParticleSpectraPlotFromList(lists8TeV[i], fCollSys[1], fEnergy[4], fCentrality[i], suffix);
                    ProduceParticleSpectraPlotFromListOnlyFinal(lists5TeV[i], fCollSys[1], fEnergy[4], fCentrality[i], suffix);
                }
            }
        }
    }

    //================================================================================================================
    //Saving the TLists to the final file
    //================================================================================================================
    cout << "writing lists" << endl;
    outputFile->cd();
    if (enable5TeV) {
        for (Int_t i=0; i<nCentralities; i++) {
            if (includeCentrality5TeV[i]) {
                lists5TeV[i]->Write(Form("%s_%s", fEnergy[2].Data(), fCentrality[i].Data()), TObject::kSingleKey);
            }
        }
    }
    outputFile->Close();
}
