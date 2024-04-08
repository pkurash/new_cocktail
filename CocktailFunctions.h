/****************************************************************************************************************************
 ******    Mike Sas, mike.sas@cern.ch                                                                                   *****
 ******    Friederike Bock, friederike.bock@cern.ch                                                                     *****
 ******    Lucas Altenkaemper, lucas.altenkamper@cern.ch                                                                *****
 ******    Lucia Leardini, lucia.leardini@cern.ch                                                                       *****
 ****************************************************************************************************************************/
#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include <string>
#include "TGaxis.h"
#include "TFractionFitter.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "THStack.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "TEllipse.h"
#include "TKey.h"
#include "Rtypes.h"
#include "TRandom.h"
#include "TSpline.h"

//================================================================================================================
//Some general definitions
//================================================================================================================
const Int_t nCentralities               = 21;
const Int_t nParticles                  = 29;
const Int_t nMethods                    = 20;
TString fCollSys[3]                     = {"pp", "pPb", "PbPb"};
TString fEnergy[6]                      = {"0.9TeV", "2.76TeV", "5.02TeV", "7TeV", "8TeV", "13TeV"};
TString fEnergyLatex[6]                 = {"#sqrt{s_{NN}} = 900 GeV", "#sqrt{s_{NN}} = 2.76 TeV", "#sqrt{s_{NN}} = 5.02 TeV", "#sqrt{s_{NN}} = 7 TeV", "#sqrt{s_{NN}} = 8 TeV", "#sqrt{s_{NN}} = 13 TeV"};
TString fEnergyLatexPP[6]               = {"#sqrt{s} = 900 GeV", "#sqrt{s} = 2.76 TeV", "#sqrt{s} = 5.02 TeV", "#sqrt{s} = 7 TeV", "#sqrt{s} = 8 TeV", "#sqrt{s} = 13 TeV"};
//                                          0     1         2         3         4         5         6         7         8         9         10        11         12        13      14      15       16     17      18      19      20   
TString fCentrality[nCentralities]      = {"MB", "0005",   "0510",   "0010",   "1020",   "0020",   "2040",   "2050",   "4060",   "6080",   "8090",   "80100",   "2030",   "3040", "2060", "60100", "NSD", "4050", "5060", "6070", "7080"};
TString fCentralityLatex[nCentralities] = {"MB", "0-5%",   "5-10%",  "0-10%",  "10-20%", "0-20%",  "20-40%", "20-50%", "40-60%", "60-80%", "80-90%", "80-100%", "20-30%", "30-40%", "20-60%", "60-100%", "NSD", "40-50%", "50-60%", "60-70%", "70-80%"};
TString fCentralityOpt2[nCentralities]  = {"MB", "00to05", "05to10", "00to10", "10to20", "00to20", "20to40", "20to50", "40to60", "60to80", "80to90", "80to100", "20to30", "30to40", "20to60", "60to100", "MB", "40to50", "50to60", "60to70", "70to80"};
TString fCentralityOpt3[nCentralities]  = {"MB", "0to5",   "5to10",  "0to10",  "10to20", "00to20", "20to40", "20to50", "40to60", "60to80", "80to90", "80to100", "20to30", "30to40", "20to60", "60to100", "MB", "40to50", "50to60", "60to70", "70to80"};
TString fCentralityOpt4[nCentralities]  = {"MB", "05",     "510",    "010",    "1020",   "020",    "2040",   "2050",   "4060",   "6080",   "8090",   "80100",   "2030",   "3040", "2060", "60100", "NSD", "4050", "5060", "6070", "7080"};
TString fCentralityOpt5[nCentralities]  = {"MB", "0.00to5.00",     "5.00to10.00",    "0.00to10.00",    "10.00to20.00",   "0.00to20.00",    "20.00to40.00",   "20.00to50.00",   "40.00to60.00",   "60.00to80.00",   "80.00to90.00",   "80.00to100.00",   "20.00to30.00",   "30.00to40.00", "20.00to60.00", "60.00to100.00", "NSD", "40.00to50.00", "50.00to60.00", "60.00to70.00", "70.00to80.00"};
TString fCentralityOpt6[nCentralities]  = {"MB", "00-05%", "05-10%", "00-10%", "10-20%", "00-20%", "20-40%", "20-50%", "40-60%", "60-80%", "80-90%", "80-100%", "20-30%", "30-40%", "20-60%", "60-100%", "MB", "40-50%", "50-60%", "60-70%", "70-80%"};
                                         //0,   1       2      3       4        5          6          7               8        9         10        11         12        13       14        15          16         17           18             19
TString fMethod[nMethods]               = {"", "Comb", "PCM", "PHOS", "EMCal", "PCMPHOS", "PCMEMCal", "EMCalMerged", "lowPt", "highPt", "Dalitz", "PCMPass2","MCParam","OldPub","EMCalMB","EMCalEMC7","EMCalEGA","PCMEMCalMB","PCMEMCalEMC7","PCMEMCalEGA"};
TString fMethodLabel[nMethods]          = {"", "Comb", "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC", "mEMC", "lowPt", "highPt", "PCM-Dal", "PCMPass2","MCParam","pub","EMC,MB","EMC,EMC7","EMC,EGA","PCM-EMC,MB","PCM-EMC,EMC7","PCM-EMC,EGA"};
TString fParticle[nParticles]           = {"NPion", "Eta", "Omega", "EtaPrime", "GammaDir", "CPion", "CKaon", "Proton", "CHadron", "Phi", "NKaonStar",
                                           "NRho", "CRho", "NDelta", "CDelta", "NKaonSubS", "Lambda", "NSigma", "CSigma", "COmega", "CXi", "JPsi",
                                           "DZero", "DPlus", "DStarPlus", "DSPlus", "CSigmaStar", "NXiStar", "NKaonSubL"};
TString fParticleLatex[nParticles]      = {"#pi^{0}", "#eta", "#omega", "#eta'", "#gamma_{dir}", "(#pi^{+}+#pi^{-})/2", "(K^{+}+K^{-})/2", "(p+#bar{p})/2", "(h^{+}+h^{-})/2",
                                           "#phi", "K^{0*}", "#rho^{0}", "(#rho^{+}+#rho^{-})/2", "#Delta^{0}", "(#Delta^{+}+#Delta^{-})/2", "K^{0}_{s}", "#Lambda", "#Sigma^{0}",
                                           "(#Sigma^{+}+#Sigma^{-})/2", "(#Omega^{-}+#bar{#Omega}^{+})/2", "(#Xi^{-}+#bar{#Xi}^{+})/ 2", "J/#psi", "D^{0}", "D^{+}", "D^{*+}",
                                           "D_{S}^{+}", "#Sigma^{*+}", "(#Xi^{*0}+#bar{#Xi}^{*0})/ 2", "K^{0}_{l}"};
Bool_t fIsMeson[nParticles]             = { kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kTRUE, kTRUE, kFALSE, kFALSE,
                                            kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kTRUE, kFALSE, kFALSE,
                                            kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kTRUE};

TString fRatioParticleLatex[nParticles] = { "#pi^{0}", "#eta", "#omega", "#eta'", "#gamma_{dir}", "(#pi^{+}+#pi^{-})", "(K^{+}+K^{-})", "(p+#bar{p})", "(h^{+}+h^{-})",
                                            "#phi", "K^{0*}", "#rho^{0}", "(#rho^{+}+#rho^{-})", "#Delta^{0}", "(#Delta^{+}+#Delta^{-})", "K^{0}_{s}", "#Lambda", "#Sigma^{0}",
                                            "(#Sigma^{+}+#Sigma^{-})/2", "(#Omega^{-}+#bar{#Omega}^{+})/2", "(#Xi^{-}+#bar{#Xi}^{+})/2", "J/#psi", "D^{0}", "D^{+}", "D^{*+}", "D_{S}^{+}",
                                            "#Sigma^{*+}", "(#Xi^{*0}+#bar{#Xi}^{*0})/ 2","K^{0}_{l}"};
TString pdgCodes[27]                    = { "111", "221", "113", "223", "331", "333", "443", "3212", "310", "2224", "2214", "1114", "2114", "213", "-213", "313", "130", "3122",
                                            "321", "-321", "-3334", "3334", "-3312", "3312", "3224", "3114", "3214"};
Int_t fParticlePDG[nParticles]          = { 111,    //#pi^{0}
                                            221,    //#eta
                                            223,    //#omega
                                            331,    //#eta'
                                            22,     //#gamma
                                            211,    //#pi^+
                                            321,    //K^+
                                            2212,   //p
                                            211,    //#pi^+
                                            333,    //#phi
                                            313,    //K^{0*}
                                            113,    //#rho^0
                                            213,    //#rho^+
                                            2114,   //#Delta^0
                                            1114,   //#Delta^0
                                            310,    //K^0_s
                                            3122,   //#Lambda
                                            3212,   //#Sigma^0
                                            3222,   //#Sigma^+
                                            3334,   //#Omega^-
                                            3312,   //#Xi^-
                                            443,    //J/#psi
                                            421,    //D^0
                                            411,    //D^+
                                            413,    //D^*+
                                            431,    //D^+_s
                                            3224,   //#Xi^{*-}
                                            3214,   //#Sigma^{*0}
                                            130     //K^0_l
                                          };
Color_t paint[29]                       = { kBlue+1, kGreen+3, kRed+1, kAzure+4, kSpring+4, kPink+8, kBlue-4, kTeal-6, kYellow+3, kOrange+2, kViolet+1,
                                            kCyan+1, kRed-3, kBlue+2, kTeal+3, kOrange+4, kPink-8, kViolet-5, kAzure+1, kSpring-7, kBlack, kGray+2,
                                            kAzure-7, kTeal+1, kViolet-7, kRed-7, kBlue+1, kGreen+3, kRed-7, };
Style_t markers[30]                     = { 20,  21,  34,  29,  24,
                                            25,  28,  33,  27,  30,
                                            20,  21,  34,  29,  24,
                                            25,  28,  33,  27,  30,
                                            20,  21,  34,  29,  24,
                                            25,  28,  33,  27,  30};
Size_t markerSize[30]                   = { 1.5, 1.4, 1.5, 2,  1.6,
                                            1.5, 1.6, 2.5, 2,  2.1,
                                            1.5, 1.4, 1.5, 2,  1.6,
                                            1.5, 1.6, 2.5, 2,  2.1,
                                            1.5, 1.4, 1.5, 2,  1.6,
                                            1.5, 1.6, 2.5, 2,  2.1};
Style_t lineStyle[10]                   = { 1,  5,  6,  7,  8,
                                            1,  5,  6,  7,  8};
Color_t lineColor[10]                   = { kBlack,  kGreen-6,  kRed-6, kGray+1,  kBlue-6,
                                            kRed-6,  kRed-4,  kGray+2,  kGreen-4,  kAzure-7};

Bool_t isThesis                         = kFALSE;
TString initFitFunctionsUsed[nParticles][nCentralities][nMethods][10];
Int_t nFitFunctionsUsed[nParticles][nCentralities][nMethods];


struct AliConvDataObject {
    Double_t valueX;
    Double_t errorXLow;
    Double_t errorXHigh;
    Double_t valueY;
    Double_t errorYStatLow;
    Double_t errorYStatHigh;
    Double_t errorYSystLow;
    Double_t errorYSystHigh;
    Double_t errorYTotLow;
    Double_t errorYTotHigh;
};

//================================================================================================================
// basic objects for TSplines
//================================================================================================================
Int_t counterGraphs = 0;
std::vector<TSpline3*> vecGraphs(0);

//================================================================================================================
// Return particle Name based on Mass
//================================================================================================================
TString GetParticleNameFromMass(Double_t mass){
    for (Int_t particleIter = 0; particleIter < nParticles; particleIter++){
        if ( mass == (Double_t)TDatabasePDG::Instance()->GetParticle(fParticlePDG[particleIter])->Mass())
            return fParticle[particleIter];
    }
    return "";
}

//================================================================================================================
// Return particle Label based on Mass
//================================================================================================================
TString GetParticleLabelFromMass(Double_t mass){
    for (Int_t particleIter = 0; particleIter < nParticles; particleIter++){
        if ( mass == (Double_t)TDatabasePDG::Instance()->GetParticle(fParticlePDG[particleIter])->Mass())
            return fParticleLatex[particleIter];
    }
    return "";
}

//================================================================================================================
// Return particle Label based on Mass
//================================================================================================================
TString GetIsMesonFromMass(Double_t mass){
    for (Int_t particleIter = 0; particleIter < nParticles; particleIter++){
        if ( mass == (Double_t)TDatabasePDG::Instance()->GetParticle(fParticlePDG[particleIter])->Mass())
            return fIsMeson[particleIter];
    }
    return "";
}
//================================================================================================================
//Return particle name from PDG code
//================================================================================================================
TString GetParticleNameFromPDG(TString particle) {
    if ( particle.Contains( "111") ){
        return "NPion";
    } else if ( particle.Contains("3312") || particle.Contains("-3312")){
        return "CXi";
    } else if ( particle.Contains("3334") || particle.Contains("-3334")){
        return "COmega";
    } else if ( particle.Contains("2212") ){
        return "Proton";
    } else if ( particle.Contains("3122") ){
        return "Lambda";
    } else if ( particle.Contains("3224") || particle.Contains("3114")){
        return "CSigma";
    } else if ( particle.Contains( "211") || particle.Contains("-211")){
        return "CPion";
    } else if ( particle.Contains( "321") || particle.Contains("-321")){
        return "CKaon";
    } else if ( particle.Contains( "221") ){
        return "Eta";
    } else if ( particle.Contains( "310") ){
        return "NKaonSubS";
    } else if ( particle.Contains( "130") ){
        return "NKaonSubL";
    } else if ( particle.Contains( "333") ){
        return "Phi";
    } else if ( particle.Contains( "313") ){
        return "NKaonStar";
    } else if ( particle.Contains( "113") ){
        return "NRho";
    } else if ( particle.Contains( "213") || particle.Contains("-213")){
        return "CRho";
    } else if ( particle.Contains( "223") ){
        return "Omega";
    } else if ( particle.Contains( "0") ){
        return "CHadron";
    } else {
        cout << "Name not found for "<< particle.Data() << " in GetParticleNameFromPDG" << endl;
        return "";
    }

}
//================================================================================================================
//Return particle color
//================================================================================================================
Color_t GetParticleColor(TString particle) {
    if (particle.Contains("NPion")){
        return kBlue+1;
    } else if (particle.Contains("Eta")){
        return kRed+1;
    } else if (particle.Contains("CPion")){
        return kBlue-6;
    } else if (particle.Contains("CKaon")){
        return kRed-6;
    } else if (particle.Contains("NKaonSubS")){
        return kPink+8;
    } else if (particle.Contains("NKaonSubL")){
        return kPink+9;
    } else if (particle.Contains("Proton")){
        return kSpring+4;
    } else if (particle.Contains("Lambda")){
        return kViolet+1;
    } else if (particle.Contains("Phi")){
        return kGreen+3;
    } else if (particle.Contains("NKaonStar")){
        return 807;
    } else if (particle.Contains("CXi")){
        return 800;
    } else if (particle.Contains("COmega")){
        return kCyan+2;
    } else if (particle.Contains("NRho")){
        return kPink-2;
    } else if (particle.Contains("CRho")){
        return kSpring+2;
    } else if (particle.Contains("CSigmaStar")){
        return 802;
    } else if (particle.Contains("NXiStar")){
        return kTeal+4;
    } else if (particle.Contains("DZero")){
        return kAzure;
    } else if (particle.Contains("DPlus")){
        return kAzure+2;
    } else if (particle.Contains("DStarPlus")){
        return kAzure-7;
    } else if (particle.Contains("DSPlus")){
        return kViolet+4;
    } else if (particle.Contains("Omega")){
        return kOrange-2;
    } else if (particle.Contains("CSigma")){
        return kPink-6;
    } else if (particle.Contains("CHadron")){
        return kViolet+2;
    } else {
        cout << "setting color for "<< particle.Data() << ", which isn't known for the color settings" << endl;
        return kGray+1;
    }

}

//================================================================================================================
//Return particle color
//================================================================================================================
Style_t GetParticleMarkerStyle(TString particle) {

    if (particle.Contains("NPion"))
        return kFullCircle;
    else if (particle.Contains("Eta"))
        return kFullSquare;
    else if (particle.Contains("CPion"))
        return kOpenCircle;
    else if (particle.Contains("CKaon"))
        return kOpenSquare;
    else if (particle.Contains("NKaonSubS"))
        return kOpenCross;
    else if (particle.Contains("NKaonSubL"))
        return kOpenStar;
    else if (particle.Contains("Proton"))
        return kFullCross;
    else if (particle.Contains("Lambda"))
        return kFullDiamond;
    else if (particle.Contains("Phi"))
        return kOpenDiamond;
    else if (particle.Contains("NKaonStar"))
        return kFullStar;
    else if (particle.Contains("CXi"))
        return kOpenStar;
    else if (particle.Contains("COmega"))
        return kFullCircle;
    else if (particle.Contains("NRho"))
        return kOpenCircle;
    else if (particle.Contains("CRho"))
        return kOpenCircle;
    else if (particle.Contains("NXiStar"))
        return kOpenSquare;
    else if (particle.Contains("CSigmaStar"))
        return kOpenCircle;
    else if (particle.Contains("DZero"))
        return kFullCross;
    else if (particle.Contains("DPlus"))
        return kFullStar;
    else if (particle.Contains("DStarPlus"))
        return kOpenCross;
    else if (particle.Contains("DSPlus"))
        return kFullDiamond;
    else if (particle.Contains("Omega"))
        return kOpenCross;
    else if (particle.Contains("CSigma"))
        return kOpenCircle;
    else
        return 1;
}

//================================================================================================================
//Return particle color
//================================================================================================================
Size_t GetParticleMarkerSize(TString particle) {

    if (particle.Contains("NPion"))
        return 1.5;
    else if (particle.Contains("Eta"))
        return 1.4;
    else if (particle.Contains("CPion"))
        return 1.6;
    else if (particle.Contains("CKaon"))
        return 1.5;
    else if (particle.Contains("NKaonSubS"))
        return 1.8;
    else if (particle.Contains("NKaonSubL"))
        return 1.8;
    else if (particle.Contains("Proton"))
        return 1.7;
    else if (particle.Contains("Lambda"))
        return 2.5;
    else if (particle.Contains("Phi"))
        return 2.6;
    else if (particle.Contains("NKaonStar"))
        return 2.3;
    else if (particle.Contains("CXi"))
        return 2.6;
    else if (particle.Contains("COmega"))
        return 1.5;
    else if (particle.Contains("NRho"))
        return 1.6;
    else if (particle.Contains("NXiStar"))
        return 1.6;
    else if (particle.Contains("CSigmaStar"))
        return 1.5;
    else if (particle.Contains("DZero"))
        return 1.7;
    else if (particle.Contains("DPlus"))
        return 2.3;
    else if (particle.Contains("DStarPlus"))
        return 1.8;
    else if (particle.Contains("DSPlus"))
        return 2.5;
    else if (particle.Contains("Omega"))
        return 1.8;
    else if (particle.Contains("CSigma"))
        return 1.6;
    else
        return 1.5;
}

//================================================================================================================
//Return fit Label
//================================================================================================================
TString GetFitLabel(TString fitName) {

    if (fitName.Contains("oHag"))
        return "mod. Hagedorn";
    else if (fitName.Contains("tcm"))
        return "two-component model";
    else if (fitName.Contains("blastwave"))
        return "blast-wave";
    else if (fitName.Contains("tmpt"))
        return "Levy-Tsallis";
    else if (fitName.Contains("maxSys"))
        return "shifted upwards";
    else if (fitName.Contains("minSys"))
        return "shifted downwards";
    else if (fitName.Contains("hardest"))
        return "hardest spectrum";
    else if (fitName.Contains("softest"))
        return "softest spectrum";
    else
        return "";
}

//================================================================================================================
//Return particle iterator
//================================================================================================================
Int_t GetParticleIterator(TString particle) {

    Int_t iterator                      = 0;
    while (particle.CompareTo(fParticle[iterator]) != 0 && iterator < nParticles)
        iterator++;
    if (iterator>=nParticles) iterator  = -1;

    return iterator;
}

//================================================================================================================
//Return centrality iterator
//================================================================================================================
Int_t GetCentralityIterator(TString centrality) {

    Int_t iterator                              = 0;
    if (centrality.CompareTo("")!=0) {
        while (centrality.CompareTo(fCentrality[iterator]) != 0 && iterator < nCentralities)
            iterator++;
        if (iterator>=nCentralities) iterator   = -1;
    }

    return iterator;
}

//================================================================================================================
//Return method iterator
//================================================================================================================
Int_t GetMethodIterator(TString method) {

    Int_t iterator                      = 0;
    while (method.CompareTo(fMethod[iterator]) != 0 && iterator < nMethods)
        iterator++;
    if (iterator>=nMethods) iterator    = -1;

    return iterator;
}

//================================================================================================================
//Return OutputString for energy
//================================================================================================================
TString GetOutputStringEnergy(TString energyInput){
    if (energyInput.CompareTo("900GeV") == 0 || energyInput.CompareTo("0.9TeV") == 0)
        return "900GeV";
    else if (energyInput.CompareTo("2760GeV") == 0 || energyInput.CompareTo("2.76TeV") == 0)
        return "2760GeV";
    else if (energyInput.CompareTo("5023GeV") == 0 || energyInput.CompareTo("5TeV") == 0 || energyInput.CompareTo("5.023TeV") == 0)
        return "5TeV";
    else if (energyInput.CompareTo("5020GeV") == 0 || energyInput.CompareTo("5.02TeV") == 0 )
        return "5020GeV";
    else if (energyInput.CompareTo("7000GeV") == 0 || energyInput.CompareTo("7TeV") == 0)
        return "7TeV";
    else if (energyInput.CompareTo("8000GeV") == 0 || energyInput.CompareTo("8TeV") == 0)
        return "8TeV";
    else if (energyInput.CompareTo("13000GeV") == 0 || energyInput.CompareTo("13TeV") == 0)
        return "13TeV";
    else
        return "";
    return "";
}

//================================================================================================================
//Initialize options for spectra parametrization for integrated yield
//================================================================================================================
Bool_t InitializeFinalCalcIntegYield(TString paramSettingsFileName = "") {

    // initialize arrays that will be read from parametrization setting file
    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t j=0; j<nCentralities; j++) {
            for (Int_t k=0; k<nMethods; k++) {
                for (Int_t l=0; l<10; l++) {
                    initFitFunctionsUsed[i][j][k][l]        = "";
                }
                nFitFunctionsUsed[i][j][k]                  = 0;
            }
        }
    }

    // read in parametrization settings file
    ifstream fileParamSettings;
    cout << "Settings-file: " << paramSettingsFileName.Data() << endl;
    fileParamSettings.open(paramSettingsFileName,ios_base::in);
    if (!fileParamSettings) {
        cout << "ERROR: Parametrization settings " << paramSettingsFileName.Data() << " not found!" << endl;
        return kFALSE;
    }

    // read settings from file
    TString particleFromFile, centralityFromFile, methodFromFile, ptConstRelSysFromFile;
    TString tempString;
    Int_t i,j,k;
    std::string line;
    for( std::string line; getline(fileParamSettings, line); ) {

        // get basic settings
        fileParamSettings >> particleFromFile >> centralityFromFile >> methodFromFile;
        if (methodFromFile.CompareTo("-") == 0)     methodFromFile                      = "";
        if (centralityFromFile.CompareTo("-") == 0) centralityFromFile                  = "MB";
        i                                                                               = GetParticleIterator(particleFromFile);
        j                                                                               = GetCentralityIterator(centralityFromFile);
        k                                                                               = GetMethodIterator(methodFromFile);
        if (i<0 || j<0 || k<0) continue;
        cout << particleFromFile.Data() << endl;

        // get parameter and parameter limits
        Int_t counter                                                                   = 0;
        fileParamSettings >> tempString;
        while (tempString.CompareTo("stop")!=0 && counter < 10) {
            if (tempString.CompareTo("-")!=0) {
                initFitFunctionsUsed[i][j][k][counter]                                  = tempString;
            }
            fileParamSettings >> tempString;
            counter++;
        }
        nFitFunctionsUsed[i][j][k]                                                      = counter;
    }

    return kTRUE;
}


//================================================================================================================
//Function that scales graphs
//================================================================================================================
TGraph* ScaleGraph (TGraph* graph, Double_t scaleFac){

    Double_t * xValue       = graph->GetX();
    Double_t * yValue       = graph->GetY();
    Int_t nPoints           = graph->GetN();

    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]           = yValue[i]*scaleFac;
    }

    graph                   = new TGraph(nPoints,xValue,yValue);
    return graph;
}


//================================================================================================================
//Function that scales graphs with asymmetric errors
//================================================================================================================
TGraphAsymmErrors* ScaleGraph(TGraphAsymmErrors* graph, Double_t scaleFac){

    Double_t* xValue                = graph->GetX();
    Double_t* yValue                = graph->GetY();
    Double_t* xErrorLow             = graph->GetEXlow();
    Double_t* xErrorHigh            = graph->GetEXhigh();
    Double_t* yErrorLow             = graph->GetEYlow();
    Double_t* yErrorHigh            = graph->GetEYhigh();
    Int_t nPoints                   = graph->GetN();

    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]                   = yValue[i]*scaleFac;
        yErrorLow[i]                = yErrorLow[i]*scaleFac;
        yErrorHigh[i]               = yErrorHigh[i]*scaleFac;
    }

    graph                           = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    return graph;
}


//================================================================================================================
//Function that scales graphs with errors
//================================================================================================================
TGraphErrors* ScaleGraph (TGraphErrors* graph, Double_t scaleFac){

    Double_t* xValue            = graph->GetX();
    Double_t* yValue            = graph->GetY();
    Double_t* xError            = graph->GetEX();
    Double_t* yError            = graph->GetEY();
    Int_t nPoints               = graph->GetN();

    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]               = yValue[i]*scaleFac;
        yError[i]               = yError[i]*scaleFac;
    }

    graph                       = new TGraphErrors(nPoints,xValue,yValue,xError,yError);
    return graph;
}

//================================================================================================================
//Function to get the xrange of a graph (there would also be a TGraph method doing this)
//================================================================================================================
Double_t GetXRangeFromGraph(TGraph* graph, Bool_t returnMax) {

    Double_t returnValue            = 0;

    if (returnMax) {
        returnValue                 = TMath::MaxElement(graph->GetN(),graph->GetX());
    } else {
        returnValue                 = TMath::MinElement(graph->GetN(),graph->GetX());
    }

    return returnValue;
}

Double_t GetXRangeFromGraph(TGraphErrors* graph, Bool_t returnMax) {

    Double_t returnValue            = 0;

    if (returnMax) {
        returnValue                 = TMath::MaxElement(graph->GetN(),graph->GetX()) + graph->GetErrorX(graph->GetN()-1);
    } else {
        returnValue                 = TMath::MinElement(graph->GetN(),graph->GetX()) - graph->GetErrorX(0);
    }

    return returnValue;
}

Double_t GetXRangeFromGraph(TGraphAsymmErrors* graph, Bool_t returnMax) {

    Double_t returnValue            = 0;

    if (returnMax) {
        returnValue                 = TMath::MaxElement(graph->GetN(),graph->GetX()) + graph->GetErrorXhigh(graph->GetN()-1);
    } else {
        returnValue                 = TMath::MinElement(graph->GetN(),graph->GetX()) - graph->GetErrorXlow(0);
    }

    return returnValue;
}

//================================================================================================================
//Function to transform TGraph to TH1D
//================================================================================================================
TH1D *GraphToHist_withErrors(TGraphAsymmErrors *graph, TString name = ""){
    Double_t* xValue        = graph->GetX();
    Double_t* yValue        = graph->GetY();
    Double_t* Exhigh        = graph->GetEXhigh();
    Double_t* Exlow         = graph->GetEXlow();
    Double_t* Eyhigh        = graph->GetEYhigh();
    Double_t* Eylow         = graph->GetEYlow();
    Int_t nPoints           = graph->GetN();

    Double_t *newBinningX   = new Double_t[nPoints+1];
    for(Int_t i = 0;i<nPoints;i++)
        newBinningX[i]      = xValue[i]-Exlow[i];
        newBinningX[nPoints] = xValue[nPoints-1]+Exhigh[nPoints-1];
    TH1D *hist              = new TH1D(name,"",nPoints,newBinningX);

    for(Int_t i = 1;i<=nPoints;i++){
      hist->SetBinContent(i,yValue[i-1]);
      if (Eyhigh[i-1]<Eylow[i-1])Eyhigh[i-1]=Eylow[i-1];
      hist->SetBinError(i,Eyhigh[i-1]);
    }

    return hist;
}

TH1D *GraphToHist_withErrors(TGraphErrors *graph, TString name = ""){
    Double_t* xValue        = graph->GetX();
    Double_t* yValue        = graph->GetY();
    Double_t* Ex            = graph->GetEX();
    Double_t* Ey            = graph->GetEY();
    Int_t nPoints           = graph->GetN();

    Double_t *newBinningX   = new Double_t[nPoints+1];
    for(Int_t i = 0;i<nPoints;i++)
        newBinningX[i]      = xValue[i]-Ex[i];
        newBinningX[nPoints] = xValue[nPoints-1]+Ex[nPoints-1];
    TH1D *hist              = new TH1D(name,"",nPoints,newBinningX);

    for(Int_t i = 1;i<=nPoints;i++){
      hist->SetBinContent(i,yValue[i-1]);
      hist->SetBinError(i,Ey[i-1]);
    }

    return hist;
}

//================================================================================================================
//Function to transform TH1D to TGraphAsymmErrors
//================================================================================================================
TGraphAsymmErrors* HistToGraph(TH1D* hist) {

    if (!hist) return NULL;

    Int_t       nBins           = hist->GetNbinsX();
    Double_t*   x               = new Double_t[nBins];
    Double_t*   xErrUp          = new Double_t[nBins];
    Double_t*   xErrDown        = new Double_t[nBins];
    Double_t*   y               = new Double_t[nBins];
    Double_t*   yErrUp          = new Double_t[nBins];
    Double_t*   yErrDown        = new Double_t[nBins];

    for (Int_t i=1; i<nBins+1; i++) {
        x[i-1]                  = hist->GetBinCenter(i);
        xErrUp[i-1]             = hist->GetBinWidth(i)/2;
        xErrDown[i-1]           = hist->GetBinWidth(i)/2;
        y[i-1]                  = hist->GetBinContent(i);
        yErrUp[i-1]             = hist->GetBinError(i);
        yErrDown[i-1]           = hist->GetBinError(i);
    }

    TGraphAsymmErrors* graph    = new TGraphAsymmErrors(nBins,x,y,xErrDown,xErrUp,yErrDown,yErrUp);
    return graph;
}


//================================================================================================================
//Function that converts any yield histogram
//================================================================================================================
TH1D* ConvertYieldHisto(TH1D* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt){

    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }

    Int_t nBins                 = input->GetNbinsX();
    Double_t newValue           = 0;
    Double_t newErrorValue      = 0;
    Double_t correctionValue    = 1;

    //correct by 2pi if specified
    if (DivideBy2pi) input->Scale(1/(2*TMath::Pi()));
    if (MultiplyBy2pi) input->Scale(2*TMath::Pi());

    for(Int_t i=0;i<nBins;i++){

        //correct by 1/Pt if specified
        if(DivideByPt)    correctionValue  = 1/(input->GetBinCenter(i+1));
        if(MultiplyByPt)  correctionValue  = input->GetBinCenter(i+1);

        //set the value and error of the bin
        input->SetBinContent(i+1,input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,input->GetBinError(i+1)*correctionValue);
    }

    return input;
}

//================================================================================================================
//Function that converts any yield histogram
//================================================================================================================
TH1F* ConvertYieldHisto(TH1F* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt){

    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }

    Int_t nBins                 = input->GetNbinsX();
    Float_t newValue           = 0;
    Float_t newErrorValue      = 0;
    Float_t correctionValue    = 1;

    //correct by 2pi if specified
    if (DivideBy2pi) input->Scale(1/(2*TMath::Pi()));
    if (MultiplyBy2pi) input->Scale(2*TMath::Pi());

    for(Int_t i=0;i<nBins;i++){

        //correct by 1/Pt if specified
        if(DivideByPt)    correctionValue  = 1/(input->GetBinCenter(i+1));
        if(MultiplyByPt)  correctionValue  = input->GetBinCenter(i+1);

        //set the value and error of the bin
        input->SetBinContent(i+1,input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,input->GetBinError(i+1)*correctionValue);
    }

    return input;
}

//================================================================================================================
//Function that converts any yield graph
//================================================================================================================
TGraph* ConvertYieldGraph(TGraph* inputGraph, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt) {

    if (!inputGraph) {
        cout << "Error: Graph is NULL" << endl;
        return NULL;
    }

    if (DivideBy2pi) inputGraph                 = ScaleGraph(inputGraph, 1/(2*TMath::Pi()));
    if (MultiplyBy2pi) inputGraph               = ScaleGraph(inputGraph, 2*TMath::Pi());

    Double_t* xValue                            = inputGraph->GetX();
    Double_t* yValue                            = inputGraph->GetY();
    Int_t nPoints                               = inputGraph->GetN();

    if (DivideByPt || MultiplyByPt) {
        Double_t correctionValue                = 1;
        for (Int_t i=0; i<nPoints; i++) {

            if (DivideByPt) correctionValue     = 1/xValue[i];
            if (MultiplyByPt) correctionValue   = xValue[i];

            yValue[i]                           = yValue[i]*correctionValue;
        }
    }

    inputGraph                                  = new TGraph(nPoints,xValue,yValue);

    return inputGraph;
}

//================================================================================================================
//Function that converts any yield graph with asymmetric errors
//================================================================================================================
TGraphAsymmErrors* ConvertYieldGraph(TGraphAsymmErrors* inputGraph, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt) {

    if (!inputGraph) {
        cout << "Error: Graph is NULL" << endl;
        return NULL;
    }

    if (DivideBy2pi) inputGraph                 = ScaleGraph(inputGraph, 1/(2*TMath::Pi()));
    if (MultiplyBy2pi) inputGraph               = ScaleGraph(inputGraph, 2*TMath::Pi());

    Double_t* xValue                            = inputGraph->GetX();
    Double_t* yValue                            = inputGraph->GetY();
    Double_t* xErrorLow                         = inputGraph->GetEXlow();
    Double_t* xErrorHigh                        = inputGraph->GetEXhigh();
    Double_t* yErrorLow                         = inputGraph->GetEYlow();
    Double_t* yErrorHigh                        = inputGraph->GetEYhigh();
    Int_t nPoints                               = inputGraph->GetN();

    if (DivideByPt || MultiplyByPt) {
        Double_t correctionValue                = 1;
        for (Int_t i=0; i<nPoints; i++) {

            if (DivideByPt) correctionValue     = 1/xValue[i];
            if (MultiplyByPt) correctionValue   = xValue[i];

            yValue[i]                           = yValue[i]*correctionValue;
            yErrorLow[i]                        = yErrorLow[i]*correctionValue;
            yErrorHigh[i]                       = yErrorHigh[i]*correctionValue;
        }
    }

    inputGraph                                  = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);

    return inputGraph;
}

//================================================================================================================
//Function that converts any yield graph with errors
//================================================================================================================
TGraphErrors* ConvertYieldGraph(TGraphErrors* inputGraph, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt) {

    if (!inputGraph) {
        cout << "Error: Graph is NULL" << endl;
        return NULL;
    }

    if (DivideBy2pi) inputGraph                 = ScaleGraph(inputGraph, 1/(2*TMath::Pi()));
    if (MultiplyBy2pi) inputGraph               = ScaleGraph(inputGraph, 2*TMath::Pi());

    Double_t* xValue                            = inputGraph->GetX();
    Double_t* yValue                            = inputGraph->GetY();
    Double_t* xError                            = inputGraph->GetEX();
    Double_t* yError                            = inputGraph->GetEY();
    Int_t nPoints                               = inputGraph->GetN();

    if (DivideByPt || MultiplyByPt) {
        Double_t correctionValue                = 1;
        for (Int_t i=0; i<nPoints; i++) {

            if (DivideByPt) correctionValue     = 1/xValue[i];
            if (MultiplyByPt) correctionValue   = xValue[i];

            yValue[i]                           = yValue[i]*correctionValue;
            yError[i]                           = yError[i]*correctionValue;
        }
    }

    inputGraph                                  = new TGraphErrors(nPoints,xValue,yValue,xError,yError);

    return inputGraph;
}

//================================================================================================================
//Function that adds up to graphs
//================================================================================================================
TGraphErrors* Add2TGraphErrorsSameBinning(TGraphErrors* inputgraph1,TGraphErrors* inputgraph2){

    Double_t* xValue        = inputgraph1->GetX();
    Double_t* xError        = inputgraph1->GetEX();

    Double_t* yValue1       = inputgraph1->GetY();
    Double_t* yError1       = inputgraph1->GetEY();
    Double_t* yValue2       = inputgraph2->GetY();
    Double_t* yError2       = inputgraph2->GetEY();

    Double_t* yValue        = inputgraph2->GetY();
    Double_t* yError        = inputgraph2->GetEY();

    Int_t nPoints           = inputgraph1->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]           = yValue1[i]+ yValue2[i];
        yError[i]           = TMath::Sqrt(TMath::Power(yError1[i],2) + TMath::Power(yError2[i],2));
    }
    TGraphErrors* returnGraph = new TGraphErrors(nPoints,xValue,yValue,xError,yError);
    return returnGraph;
}

//================================================================================================================
//Function that adds up N graphs
//================================================================================================================
TGraphErrors* AddNTGraphErrorsSameBinning( TGraphErrors** inputgraphs,
                                           Int_t ngraphs                    = 1,
                                           Int_t errorMode                  = 0         //0 errors added in quadrature, 1 find largest relative error
                                         ){

    TGraphErrors*  dummy    = (TGraphErrors*)inputgraphs[0]->Clone("dummy");
    Double_t* xValue        = dummy->GetX();
    Double_t* xError        = dummy->GetEX();

    Double_t* yValue        = dummy->GetY();
    Double_t* yError        = dummy->GetEY();

    Int_t nPoints           = dummy->GetN();
    for (Int_t k = 0; k < ngraphs; k++){
        if (inputgraphs[k]->GetN() != nPoints){
            cout << "ERROR: graphs don't have the same number of points" << endl;
            if ( inputgraphs[k]->GetN() < nPoints){
                nPoints     = inputgraphs[k]->GetN();
            }
        }
        for (Int_t i = 0; i < nPoints; i++){
            if (inputgraphs[k]->GetX()[i] != xValue[i]){
                cout << "ERROR: graphs don't have the same binning" << endl;
                return NULL;
            }
        }
    }

    if (errorMode == 1){ // finding the largest relative error
        for (Int_t i = 0; i< nPoints; i++) {
            yValue[i]                   = 0;
            yError[i]                   = 0;

            std::vector<Double_t> sysErrsVector(ngraphs);
            for (Int_t k = 0; k< ngraphs; k++){
                yValue[i]               = yValue[i] + inputgraphs[k]->GetY()[i];
                if (((TGraphErrors*)inputgraphs[k])->GetY()[i])
                    sysErrsVector[k]    = ((TGraphErrors*)inputgraphs[k])->GetEY()[i]/((TGraphErrors*)inputgraphs[k])->GetY()[i];
            }
            std::sort (sysErrsVector.begin(), sysErrsVector.end());
            yError[i]                   = sysErrsVector[sysErrsVector.size()-1] * yValue[i];
        }
    } else { // adding errors in quadrature
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]           = 0;
            yError[i]           = 0;
            for (Int_t k = 0; k< ngraphs; k++){
                yValue[i]       = yValue[i] + inputgraphs[k]->GetY()[i];
                yError[i]       = yError[i] + TMath::Power(inputgraphs[k]->GetEY()[i],2);
            }
            yError[i]           = TMath::Sqrt(yError[i]);
        }
    }
    TGraphErrors* returnGraph = new TGraphErrors(nPoints,xValue,yValue,xError,yError);
    return returnGraph;
}

//================================================================================================================
//Function that adds up to graphs
//================================================================================================================
TGraphAsymmErrors* Add2TGraphErrorsSameBinning(TGraphAsymmErrors* inputgraph1,TGraphAsymmErrors* inputgraph2){

    Double_t* xValue        = inputgraph1->GetX();
    Double_t* xErrorUp      = inputgraph1->GetEXhigh();
    Double_t* xErrorDown    = inputgraph1->GetEXlow();

    Double_t* yValue1       = inputgraph1->GetY();
    Double_t* yErrorUp1     = inputgraph1->GetEYhigh();
    Double_t* yErrorDown1   = inputgraph1->GetEYlow();
    Double_t* yValue2       = inputgraph2->GetY();
    Double_t* yErrorUp2     = inputgraph2->GetEYhigh();
    Double_t* yErrorDown2   = inputgraph2->GetEYlow();

    Double_t* yValue        = inputgraph1->GetY();
    Double_t* yErrorUp      = inputgraph1->GetEYhigh();
    Double_t* yErrorDown    = inputgraph1->GetEYlow();

    Int_t nPoints           = inputgraph1->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]           = yValue1[i] + yValue2[i];
        yErrorUp[i]         = TMath::Sqrt(yErrorUp1[i]*yErrorUp1[i] + yErrorUp2[i]*yErrorUp2[i]);
        yErrorDown[i]       = TMath::Sqrt(yErrorDown1[i]*yErrorDown1[i] + yErrorDown2[i]*yErrorDown2[i]);
    }

    TGraphAsymmErrors* returnGraph = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorDown,xErrorUp,yErrorDown,yErrorUp);
    return returnGraph;
}


//================================================================================================================
//Function that adds up N graphs
//================================================================================================================
TGraphAsymmErrors* AddNTGraphErrorsSameBinning(TGraphAsymmErrors** inputgraphs,
                                               Int_t ngraphs                    = 1,
                                               Int_t errorMode                  = 0         //0 errors added in quadrature, 1 find largest relative error
                                              ){
    TGraphAsymmErrors*  dummy    = (TGraphAsymmErrors*)inputgraphs[0]->Clone("dummy");
    Double_t* xValue        = dummy->GetX();
    Double_t* xErrorLow     = dummy->GetEXlow();
    Double_t* xErrorHigh    = dummy->GetEXhigh();

    Double_t* yValue        = dummy->GetY();
    Double_t* yErrorLow     = dummy->GetEYlow();
    Double_t* yErrorHigh    = dummy->GetEYhigh();

    Int_t nPoints           = dummy->GetN();
    for (Int_t k = 0; k < ngraphs; k++){
        if (inputgraphs[k]->GetN() != nPoints){
            cout << "ERROR: graphs don't have the same number of points" << endl;
            cout << k << "\t" << inputgraphs[k]->GetN() << "\t" << nPoints << endl;;
            if ( inputgraphs[k]->GetN() < nPoints){
                nPoints     = inputgraphs[k]->GetN();
            }
        }
        for (Int_t i = 0; i < nPoints; i++){
            if (inputgraphs[k]->GetX()[i] != xValue[i]){
                cout << "ERROR: graphs don't have the same binning" << endl;
                for (Int_t l = 0; l < ngraphs; l++){
                    cout << "input \t" << l << endl;
                    inputgraphs[l]->Print();
                }
                return NULL;
            }
        }
    }
    if (errorMode == 1){ // finding the largest relative error
        for (Int_t i = 0; i< nPoints; i++) {
            yValue[i]                   = 0;
            yErrorLow[i]                = 0;
            yErrorHigh[i]               = 0;

            std::vector<Double_t> sysErrsHVector(ngraphs);
            std::vector<Double_t> sysErrsLVector(ngraphs);
            for (Int_t k = 0; k< ngraphs; k++){
                yValue[i]               = yValue[i] + inputgraphs[k]->GetY()[i];
                if (((TGraphAsymmErrors*)inputgraphs[k])->GetY()[i]){
                    sysErrsHVector[k]   = ((TGraphAsymmErrors*)inputgraphs[k])->GetEYhigh()[i]/((TGraphAsymmErrors*)inputgraphs[k])->GetY()[i];
                    sysErrsLVector[k]   = ((TGraphAsymmErrors*)inputgraphs[k])->GetEYlow()[i]/((TGraphAsymmErrors*)inputgraphs[k])->GetY()[i];
                }
            }
            std::sort (sysErrsHVector.begin(), sysErrsHVector.end());
            std::sort (sysErrsLVector.begin(), sysErrsLVector.end());
            yErrorLow[i]                = sysErrsLVector[sysErrsLVector.size()-1] * yValue[i];
            yErrorHigh[i]               = sysErrsHVector[sysErrsHVector.size()-1] * yValue[i];
        }
    } else {
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]           = 0;
            yErrorLow[i]        = 0;
            yErrorHigh[i]       = 0;
            for (Int_t k = 0; k< ngraphs; k++){
                yValue[i]       = yValue[i] + inputgraphs[k]->GetY()[i];
                yErrorLow[i]    = yErrorLow[i] + TMath::Power(inputgraphs[k]->GetEYlow()[i],2);
                yErrorHigh[i]   = yErrorHigh[i] + TMath::Power(inputgraphs[k]->GetEYhigh()[i],2);
            }
            yErrorLow[i]        = TMath::Sqrt(yErrorLow[i]);
            yErrorHigh[i]       = TMath::Sqrt(yErrorHigh[i]);
        }
    }
    TGraphAsymmErrors* returnGraph = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh, yErrorLow, yErrorHigh);
    return returnGraph;
}


//================================================================================================================
//Function that returns a histogram from the generalized input file
//================================================================================================================
TH1D* GetHistoFromOfficialInput(TFile* filename, TH1D* histo,TString histoName,TString histoOutputName,TList* list){

    histo = (TH1D*)filename->Get(histoName.Data());
    histo->SetName(histoOutputName.Data());
    list->Add(histo);

    return histo;
}

//================================================================================================================
//Function that returns a histogram from the generalized input file
//================================================================================================================
TH1D* GetHistoFromGeneralizedInput(TFile* inputfile, TString particle, TString centrality){

  TList *list = (TList*)inputfile->Get(Form("PbPb_2760_cent_%s",centrality.Data()));
  if(list==NULL) cout << "CENTRALITY LIST NOT FOUND IN THE INPUT FILE" << endl;

  TString histoName = "";
  //here we add a lot of if statements to collect the right histograms
  if(particle.CompareTo("ChargedPion") == 0) histoName = Form("ChargedPion_HighPtStat%s",centrality.Data());

  //finally get the histogram and return it
  TH1D* histo = (TH1D*)list->FindObject(histoName.Data());
  gStyle->SetOptStat(0);
  histo->SetTitle(Form("%s cent %s", particle.Data(), centrality.Data()));
  histo->GetXaxis()->SetTitle("#it{p}_{T}(GeV/c)");
  histo->GetYaxis()->SetTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dydp_{T}} (GeV/c)^{-2}");
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleSize(0.03);
  histo->GetXaxis()->SetTitleOffset(0.9);
  histo->GetYaxis()->SetTitleOffset(1.5);
  histo->GetXaxis()->SetRangeUser(0,20);

  if(histo!=0){
    return histo;
  }else{
    cout << "HISTOGRAM NOT FOUND IN THE INPUT FILE" << endl;
    return 0;
  }
}

//================================================================================================================
// Function to add data point at zero
//================================================================================================================
void SetPointAtZero(TGraphAsymmErrors* graph, Float_t px, Float_t py, Float_t pxErr, Float_t pyErr){

    if (graph) {
        graph->Set(graph->GetN()+1);
        graph->SetPoint(graph->GetN()-1, px, py);
        graph->SetPointError(graph->GetN()-1, pxErr, pxErr, pyErr, pyErr);
        graph->Sort();
        return;
    } else {
        cout << "GRAPH IS NULL" << endl;
        return;
    }
}

//================================================================================================================
// Function to add data point at zero
//================================================================================================================
void SetPointAtZero(TH1* histo, Float_t px, Float_t py, Float_t pxErr, Float_t pyErr){

    if (histo) {
        histo->SetBinContent(histo->FindBin(px),py);
        histo->SetBinError(histo->FindBin(px),pyErr);
        return;
    } else {
        cout << "HISTO IS NULL" << endl;
        return;
    }
}

//================================================================================================================
// Function to properly set the properties of the graphs
//================================================================================================================
void SetGraphProperties(TGraph* graph, TString nameGraph, TString xAxisTitle, TString yAxisTitle, TString graphTitle = ""){

    if (graph) {
        graph->SetName(nameGraph.Data());
        graph->SetTitle(graphTitle.Data());
        graph->GetXaxis()->SetTitle(xAxisTitle.Data());
        graph->GetYaxis()->SetTitle(yAxisTitle.Data());
        return;
    } else {
        cout << "GRAPH IS NULL" << endl;
        return;
    }
}

//================================================================================================================
// Function to properly set the properties of the histo
//================================================================================================================
void SetHistoProperties(TH1* histo, TString nameGraph, TString xAxisTitle, TString yAxisTitle, TString graphTitle = ""){

    if (histo) {
        histo->SetName(nameGraph.Data());
        histo->SetTitle(graphTitle.Data());
        histo->GetXaxis()->SetTitle(xAxisTitle.Data());
        histo->GetYaxis()->SetTitle(yAxisTitle.Data());
        return;
    } else {
        cout << "HISTOGRAM IS NULL" << endl;
        return;
    }
}

//================================================================================================================
//Function that determines if the relative syst. err. given in the settings file should be used
//return values:    0) pt const. syst. err. is zero                     -> don't use
//                  1) pt const. syst. err. given and in range          -> use
//                  2) pt const. syst. err. given but exceeds limits    -> use full syst. err.
//================================================================================================================
Int_t UseRelPtConstSystErr(TH1D* systHist, Double_t relPtConstSyst) {

    if (!systHist) return 0;

    // no relative pt const syst err given, return 0
    if (!relPtConstSyst) return 0;

    // relative pt const syst err != 0, check for size
    Double_t relSystErr = 0.;
    for (Int_t i=1; i<=systHist->GetNbinsX(); i++) {
        relSystErr      = systHist->GetBinError(i) / systHist->GetBinContent(i);
        if (relPtConstSyst > relSystErr) return 2;
    }
    return 1;
}

Int_t UseRelPtConstSystErr(TGraphErrors* systGraph, Double_t relPtConstSyst) {

    if (!systGraph) return 0;

    Int_t       nPoints     = systGraph->GetN();
    Double_t*   yValue      = systGraph->GetY();
    Double_t*   yError      = systGraph->GetEY();

    // no relative pt const syst err given, return 0
    if (!relPtConstSyst) return 0;

    // relative pt const syst err != 0, check for size
    Double_t relSystErr     = 0.;
    for (Int_t i=0; i<nPoints; i++) {
        relSystErr          = yError[i] / yValue[i];
        if (relPtConstSyst > relSystErr) return 2;
    }
    return 1;
}

Int_t UseRelPtConstSystErr(TGraphAsymmErrors* systGraph, Double_t relPtConstSyst) {

    if (!systGraph) return 0;

    Int_t       nPoints     = systGraph->GetN();
    Double_t*   yValue      = systGraph->GetY();
    Double_t*   yErrorUp    = systGraph->GetEYhigh();
    Double_t*   yErrorDown  = systGraph->GetEYlow();

    // no relative pt const syst err given, return 0
    if (!relPtConstSyst) return 0;

    // relative pt const syst err != 0, check for size
    Double_t relSystErrUp   = 0.;
    Double_t relSystErrDown = 0.;
    for (Int_t i=0; i<nPoints; i++) {
        relSystErrUp        = yErrorUp[i] / yValue[i];
        relSystErrDown      = yErrorDown[i] / yValue[i];
        if ( (relPtConstSyst > relSystErrUp) || (relPtConstSyst > relSystErrDown) ) return 2;
    }
    return 1;
}

//================================================================================================================
//Function that shifts the input spectrum according to systematic errors
//================================================================================================================

// ___________ TH1 type input spectra ____________________________________________________________________________
TH1D* ShiftSpectraWithSyst(TH1D* centralValue, TH1D* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat histogram is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetNbinsX() || !sysErr->GetNbinsX()) return NULL;

    TH1D* shiftedHistogram                  = (TH1D*)centralValue->Clone(centralValue->GetName());

    Double_t sign                           = 0;
    if (shiftUp) {
        sign                                = 1;
        shiftedHistogram->SetName(Form("%s_constShiftUp", centralValue->GetName()));
    } else {
        sign                                = -1;
        shiftedHistogram->SetName(Form("%s_constShiftDown", centralValue->GetName()));
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (!UseRelPtConstSystErr(sysErr, relPtConstSyst)) return shiftedHistogram;
    for (Int_t i = 1; i <= centralValue->GetNbinsX(); i++) {
        if (!centralValue->GetBinContent(i)) continue;
        if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 1) {
            shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+sign*relPtConstSyst*centralValue->GetBinContent(i));
        } else {
            for (Int_t j = 1; j <= sysErr->GetNbinsX(); j++) {
                if (centralValue->GetBinCenter(i) == sysErr->GetBinCenter(j)) continue;
                shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+sign*sysErr->GetBinError(j));
            }
        }
    }

    return shiftedHistogram;
}

TH1D* ShiftSpectraWithSyst(TH1D* centralValue, TGraphErrors* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat histogram is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetNbinsX() || !sysErr->GetN()) return NULL;

    TH1D*       shiftedHistogram            = (TH1D*)centralValue->Clone(centralValue->GetName());
    Int_t       nPoints                     = sysErr->GetN();
    Double_t*   xValue                      = sysErr->GetX();
    Double_t*   yValue                      = sysErr->GetY();
    Double_t*   yError                      = sysErr->GetEY();

    Double_t sign                           = 0;
    if (shiftUp) {
        sign                                = 1;
        shiftedHistogram->SetName(Form("%s_constShiftUp", centralValue->GetName()));
    } else {
        sign                                = -1;
        shiftedHistogram->SetName(Form("%s_constShiftDown", centralValue->GetName()));
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (!UseRelPtConstSystErr(sysErr, relPtConstSyst)) return shiftedHistogram;
    for (Int_t i = 1; i <= centralValue->GetNbinsX(); i++) {
        if (!centralValue->GetBinContent(i)) continue;
        if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 1) {
            shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+sign*relPtConstSyst*centralValue->GetBinContent(i));
        } else {
            for (Int_t j=0; j<nPoints; j++) {
                if(centralValue->GetBinCenter(i)!=xValue[j]) continue;
                shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+sign*yError[j]);
            }
        }
    }

    return shiftedHistogram;
}

TH1D* ShiftSpectraWithSyst(TH1D* centralValue, TGraphAsymmErrors* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat histogram is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetNbinsX() || !sysErr->GetN()) return NULL;

    TH1D*       shiftedHistogram            = (TH1D*)centralValue->Clone(centralValue->GetName());
    Int_t       nPoints                     = sysErr->GetN();
    Double_t*   xValue                      = sysErr->GetX();
    Double_t*   yValue                      = sysErr->GetY();
    Double_t*   yError                      = NULL;

    Double_t sign                           = 0;
    if (shiftUp) {
        sign                                = 1;
        yError                              = sysErr->GetEYhigh();
        shiftedHistogram->SetName(Form("%s_constShiftUp", centralValue->GetName()));
    } else {
        sign                                = -1;
        yError                              = sysErr->GetEYlow();
        shiftedHistogram->SetName(Form("%s_constShiftDown", centralValue->GetName()));
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (!UseRelPtConstSystErr(sysErr, relPtConstSyst)) return shiftedHistogram;
    for (Int_t i = 1; i <= centralValue->GetNbinsX(); i++) {
        if (!centralValue->GetBinContent(i)) continue;
        if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 1) {
            shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+sign*relPtConstSyst*centralValue->GetBinContent(i));
        } else {
            for (Int_t j=0; j<nPoints; j++) {
                if(centralValue->GetBinCenter(i)!=xValue[j]) continue;
                shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+sign*yError[j]);
            }
        }
    }

    return shiftedHistogram;
}

// ___________ TGraphErrors type input spectra ___________________________________________________________________
TGraphErrors* ShiftSpectraWithSyst(TGraphErrors* centralValue, TH1D* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetNbinsX()) return NULL;

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xError                                    = centralValue->GetEX();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yError                                    = centralValue->GetEY();

    TString tempName                                    = "";
    Double_t sign                                       = 0;
    if (shiftUp) {
        sign                                            = 1;
        tempName                                        = Form("%s_constShiftUp", centralValue->GetName());
    } else {
        sign                                            = -1;
        tempName                                        = Form("%s_constShiftDown", centralValue->GetName());
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) != 0) {
        for (Int_t i = 0; i < nPoints; i++) {
            if (!yValue[i]) continue;
            if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 1) {
                yValueShifted[i]                        = yValue[i] * (1 + sign*relPtConstSyst);
            } else {
                for (Int_t j=1; j<=sysErr->GetNbinsX(); j++) {
                    if (xValue[i]!=sysErr->GetBinCenter(j)) continue;
                    yValueShifted[i]                    = yValue[i] + sign*sysErr->GetBinError(j);
                }
            }
        }
    }

    TGraphErrors* returnGraph                           = new TGraphErrors(nPoints,xValue,yValueShifted,xError,yError);
    returnGraph->SetName(tempName);
    return returnGraph;
}

TGraphErrors* ShiftSpectraWithSyst(TGraphErrors* centralValue, TGraphErrors* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetN()) return NULL;

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xError                                    = centralValue->GetEX();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yError                                    = centralValue->GetEY();

    Double_t* yErrorSys                                 = sysErr->GetEY();

    TString tempName                                    = "";
    Double_t sign                                       = 0;
    if (shiftUp) {
        sign                                            = 1;
        tempName                                        = Form("%s_constShiftUp", centralValue->GetName());
    } else {
        sign                                            = -1;
        tempName                                        = Form("%s_constShiftDown", centralValue->GetName());
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) != 0) {
        for (Int_t i = 0; i < nPoints; i++) {
            if (!yValue[i]) continue;
            if (UseRelPtConstSystErr(sysErr, relPtConstSyst)==1) {
                yValueShifted[i]                            = yValue[i] * (1 + sign*relPtConstSyst);
            } else {
                yValueShifted[i]                            = yValue[i] + sign*yErrorSys[i];
            }
        }
    }

    TGraphErrors* returnGraph                           = new TGraphErrors(nPoints,xValue,yValueShifted,xError,yError);
    returnGraph->SetName(tempName);
    return returnGraph;
}

TGraphErrors* ShiftSpectraWithSyst(TGraphErrors* centralValue, TGraphAsymmErrors* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetN()) return NULL;

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xError                                    = centralValue->GetEX();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yError                                    = centralValue->GetEY();

    Double_t* yErrorSys                                 = NULL;

    TString tempName                                    = "";
    Double_t sign                                       = 0;
    if (shiftUp) {
        yErrorSys                                       = sysErr->GetEYhigh();
        sign                                            = 1;
        tempName                                        = Form("%s_constShiftUp", centralValue->GetName());
    } else {
        yErrorSys                                       = sysErr->GetEYlow();
        sign                                            = -1;
        tempName                                        = Form("%s_constShiftDown", centralValue->GetName());
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) != 0) {
        for (Int_t i = 0; i < nPoints; i++) {
            if (!yValue[i]) continue;
            if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 1) {
                yValueShifted[i]                            = yValue[i] * (1 + sign*relPtConstSyst);
            } else {
                yValueShifted[i]                            = yValue[i] + sign*yErrorSys[i];
            }
        }
    }

    TGraphErrors* returnGraph                           = new TGraphErrors(nPoints,xValue,yValueShifted,xError,yError);
    returnGraph->SetName(tempName);
    return returnGraph;
}

// ___________ TGraphAsymmErrors type input spectra ______________________________________________________________
TGraphAsymmErrors* ShiftSpectraWithSyst(TGraphAsymmErrors* centralValue, TH1D* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetNbinsX()) return NULL;

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xErrorLow                                 = centralValue->GetEXlow();
    Double_t* xErrorHigh                                = centralValue->GetEXhigh();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yErrorLow                                 = centralValue->GetEYlow();
    Double_t* yErrorHigh                                = centralValue->GetEYhigh();

    TString tempName                                    = "";
    Double_t sign                                       = 0;
    if (shiftUp) {
        sign                                            = 1;
        tempName                                        = Form("%s_constShiftUp", centralValue->GetName());
    } else {
        sign                                            = -1;
        tempName                                        = Form("%s_constShiftDown", centralValue->GetName());
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) != 0) {
        for (Int_t i = 0; i < nPoints; i++) {
            if (!yValue[i]) continue;
            if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 1) {
                yValueShifted[i]                            = yValue[i] * (1 + sign*relPtConstSyst);
            } else {
                for (Int_t j=1; j<=sysErr->GetNbinsX(); j++) {
                    if (xValue[i]!=sysErr->GetBinCenter(j)) continue;
                    yValueShifted[i]                        = yValue[i] + sign*sysErr->GetBinError(j);
                }
            }
        }
    }

    TGraphAsymmErrors* returnGraph                      = new TGraphAsymmErrors(nPoints,xValue,yValueShifted,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    returnGraph->SetName(tempName);
    return returnGraph;
}

TGraphAsymmErrors* ShiftSpectraWithSyst(TGraphAsymmErrors* centralValue, TGraphErrors* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetN()) return NULL;

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xErrorLow                                 = centralValue->GetEXlow();
    Double_t* xErrorHigh                                = centralValue->GetEXhigh();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yErrorLow                                 = centralValue->GetEYlow();
    Double_t* yErrorHigh                                = centralValue->GetEYhigh();

    Double_t* yErrorSys                                 = sysErr->GetEY();

    TString tempName                                    = "";
    Double_t sign                                       = 0;
    if (shiftUp) {
        sign                                            = 1;
        tempName                                        = Form("%s_constShiftUp", centralValue->GetName());
    } else {
        sign                                            = -1;
        tempName                                        = Form("%s_constShiftDown", centralValue->GetName());
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) != 0) {
        for (Int_t i = 0; i < nPoints; i++) {
            if (!yValue[i]) continue;
            if (UseRelPtConstSystErr(sysErr, relPtConstSyst)==1) {
                yValueShifted[i]                            = yValue[i] * (1 + sign*relPtConstSyst);
            } else {
                yValueShifted[i]                            = yValue[i] + sign*yErrorSys[i];
            }
        }
    }

    TGraphAsymmErrors* returnGraph                      = new TGraphAsymmErrors(nPoints,xValue,yValueShifted,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    returnGraph->SetName(tempName);
    return returnGraph;
}

TGraphAsymmErrors* ShiftSpectraWithSyst(TGraphAsymmErrors* centralValue, TGraphAsymmErrors* sysErr, Bool_t shiftUp, Double_t relPtConstSyst){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetN()) return NULL;

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xErrorLow                                 = centralValue->GetEXlow();
    Double_t* xErrorHigh                                = centralValue->GetEXhigh();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yErrorLow                                 = centralValue->GetEYlow();
    Double_t* yErrorHigh                                = centralValue->GetEYhigh();

    Double_t* yErrorSys                                 = NULL;

    TString tempName                                    = "";
    Double_t sign                                       = 0;
    if (shiftUp) {
        yErrorSys                                       = sysErr->GetEYhigh();
        sign                                            = 1;
        tempName                                        = Form("%s_constShiftUp", centralValue->GetName());
    } else {
        yErrorSys                                       = sysErr->GetEYlow();
        sign                                            = -1;
        tempName                                        = Form("%s_constShiftDown", centralValue->GetName());
    }

    // constant shift if pt-const rel. syst. err. is given or larger than the syst. errs. assigned (then use full syst. err.)
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) != 0) {
        for (Int_t i = 0; i < nPoints; i++) {
            if (!yValue[i]) continue;
            if (UseRelPtConstSystErr(sysErr, relPtConstSyst)==1) {
                yValueShifted[i]                            = yValue[i] * (1 + sign*relPtConstSyst);
            } else {
                yValueShifted[i]                            = yValue[i] + sign*yErrorSys[i];
            }
        }
    }

    TGraphAsymmErrors* returnGraph                      = new TGraphAsymmErrors(nPoints,xValue,yValueShifted,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    returnGraph->SetName(tempName);
    return returnGraph;
}

//================================================================================================================
// Function that shifts spectra with pt dependent slope for systematic errors evaluation
//================================================================================================================

// ___________ TH1 type input spectra ____________________________________________________________________________
TH1D* ShiftSpectraWithSlopeSyst(TH1D* centralValue, TH1D* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat histogram is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetNbinsX() || !sysErr->GetNbinsX()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString          = "";
    if (mode == 1)      modeString          = "lin";
    else if (mode == 2) modeString          = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TH1D* shiftedHistogram                  = (TH1D*)centralValue->Clone(centralValue->GetName());
    if (shiftA) shiftedHistogram->SetName(Form("%s_%sShiftA", centralValue->GetName(), modeString.Data()));
    else        shiftedHistogram->SetName(Form("%s_%sShiftB", centralValue->GetName(), modeString.Data()));

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                      = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                      = 0.;

    Double_t factor                         = 0.;
    Double_t effectiveError                 = 0.;
    for (Int_t i = 1; i <= centralValue->GetNbinsX(); i++) {
        if (!centralValue->GetBinContent(i)) continue;
        for (Int_t j = 1; j <= sysErr->GetNbinsX(); j++) {
            if (!(centralValue->GetBinCenter(i) == sysErr->GetBinCenter(j))) continue;
            factor                          = shiftWeight->Eval(centralValue->GetBinCenter(i));
            if (relPtConstSyst) {
                effectiveError              = sysErr->GetBinError(j) - relPtConstSyst*centralValue->GetBinContent(i);
                effectiveError              = factor * effectiveError;
            } else {
                effectiveError              = factor * sysErr->GetBinError(j);
            }
            shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+effectiveError);
        }
    }

    return shiftedHistogram;
}

TH1D* ShiftSpectraWithSlopeSyst(TH1D* centralValue, TGraphErrors* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat histogram is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetNbinsX() || !sysErr->GetN()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString          = "";
    if (mode == 1)      modeString          = "lin";
    else if (mode == 2) modeString          = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TH1D* shiftedHistogram                  = (TH1D*)centralValue->Clone(centralValue->GetName());
    if (shiftA) shiftedHistogram->SetName(Form("%s_%sShiftA", centralValue->GetName(), modeString.Data()));
    else        shiftedHistogram->SetName(Form("%s_%sShiftB", centralValue->GetName(), modeString.Data()));

    Double_t* xValue                        = sysErr->GetX();
    Double_t* yValue                        = sysErr->GetY();
    Double_t* yError                        = sysErr->GetEY();

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                      = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                      = 0.;

    Double_t factor                         = 0.;
    Double_t effectiveError                 = 0.;
    for (Int_t i = 1; i <= centralValue->GetNbinsX(); i++) {
        if (!centralValue->GetBinContent(i)) continue;
        for (Int_t j = 0; j < sysErr->GetN(); j++) {
            if (!(centralValue->GetBinCenter(i) == xValue[j])) continue;
            factor                          = shiftWeight->Eval(centralValue->GetBinCenter(i));
            if (relPtConstSyst) {
                effectiveError              = yError[j] - relPtConstSyst*centralValue->GetBinContent(i);
                effectiveError              = factor * effectiveError;
            } else {
                effectiveError              = factor * yError[j];
            }
            shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+effectiveError);
        }
    }

    return shiftedHistogram;
}

TH1D* ShiftSpectraWithSlopeSyst(TH1D* centralValue, TGraphAsymmErrors* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat histogram is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetNbinsX() || !sysErr->GetN()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString          = "";
    if (mode == 1)      modeString          = "lin";
    else if (mode == 2) modeString          = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TH1D* shiftedHistogram                  = (TH1D*)centralValue->Clone(centralValue->GetName());
    if (shiftA) shiftedHistogram->SetName(Form("%s_%sShiftA", centralValue->GetName(), modeString.Data()));
    else        shiftedHistogram->SetName(Form("%s_%sShiftB", centralValue->GetName(), modeString.Data()));

    Double_t* xValue                        = sysErr->GetX();
    Double_t* yValue                        = sysErr->GetY();
    Double_t* yErrorUp                      = sysErr->GetEYhigh();
    Double_t* yErrorDown                    = sysErr->GetEYlow();

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                      = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                      = 0.;

    Double_t factor                         = 0.;
    Double_t effectiveError                 = 0.;
    for (Int_t i = 1; i <= centralValue->GetNbinsX(); i++) {
        if (!centralValue->GetBinContent(i)) continue;
        for (Int_t j = 0; j < sysErr->GetN(); j++) {
            if (!(centralValue->GetBinCenter(i) == xValue[j])) continue;
            factor                          = shiftWeight->Eval(centralValue->GetBinCenter(i));
            if (factor > 0) {
                if (relPtConstSyst) {
                    effectiveError          = yErrorUp[j] - relPtConstSyst*centralValue->GetBinContent(i);
                    effectiveError          = factor * effectiveError;
                } else {
                    effectiveError          = factor * yErrorUp[j];
                }
            } else {
                if (relPtConstSyst) {
                    effectiveError          = yErrorDown[j] - relPtConstSyst*centralValue->GetBinContent(i);
                    effectiveError          = factor * effectiveError;
                } else {
                    effectiveError          = factor * yErrorDown[j];
                }
            }
            shiftedHistogram->SetBinContent(i, centralValue->GetBinContent(i)+effectiveError);
        }
    }

    return shiftedHistogram;
}

// ___________ TGraphErrors type input spectra ___________________________________________________________________
TGraphErrors* ShiftSpectraWithSlopeSyst(TGraphErrors* centralValue, TH1D* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetNbinsX()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString                      = "";
    if (mode == 1)      modeString                      = "lin";
    else if (mode == 2) modeString                      = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TString     tempName                                = "";
    if (shiftA) tempName                                = Form("%s_%sShiftA", centralValue->GetName(), modeString.Data());
    else        tempName                                = Form("%s_%sShiftB", centralValue->GetName(), modeString.Data());

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xError                                    = centralValue->GetEX();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yError                                    = centralValue->GetEY();

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                                  = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                                  = 0.;

    Double_t factor                                     = 0.;
    Double_t effectiveError                             = 0.;
    for (Int_t i = 0; i < nPoints; i++) {
        if (!yValue[i]) continue;
        for (Int_t j=1; j<= sysErr->GetNbinsX(); j++) {
            if (!(xValue[i] == sysErr->GetBinCenter(j))) continue;
            factor                                      = shiftWeight->Eval(xValue[i]);
            if (relPtConstSyst) {
                effectiveError                          = sysErr->GetBinError(j) - relPtConstSyst*yValue[i];
                effectiveError                          = factor * effectiveError;
            } else {
                effectiveError                          = factor * sysErr->GetBinError(j);
            }
            yValueShifted[i]                            = yValue[i] + effectiveError;
        }
    }

    TGraphErrors* returnGraph                           = new TGraphErrors(nPoints,xValue,yValueShifted,xError,yError);
    returnGraph->SetName(tempName);
    return returnGraph;
}

TGraphErrors* ShiftSpectraWithSlopeSyst(TGraphErrors* centralValue, TGraphErrors* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetN()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString                      = "";
    if (mode == 1)      modeString                      = "lin";
    else if (mode == 2) modeString                      = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TString     tempName                                = "";
    if (shiftA) tempName                                = Form("%s_%sShiftA", centralValue->GetName(), modeString.Data());
    else        tempName                                = Form("%s_%sShiftB", centralValue->GetName(), modeString.Data());


    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xError                                    = centralValue->GetEX();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yError                                    = centralValue->GetEY();

    Double_t* yErrorSys                                 = sysErr->GetEY();

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                                  = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                                  = 0.;

    Double_t factor                                     = 0.;
    Double_t effectiveError                             = 0.;
    for (Int_t i = 0; i < nPoints; i++) {
        if (!yValue[i]) continue;
        factor                                          = shiftWeight->Eval(xValue[i]);
        if (relPtConstSyst) {
            effectiveError                              = yErrorSys[i] - relPtConstSyst*yValue[i];
            effectiveError                              = factor * effectiveError;
        } else {
            effectiveError                              = factor * yErrorSys[i];
        }
        yValueShifted[i]                                = yValue[i] + effectiveError;
    }

    TGraphErrors* returnGraph                           = new TGraphErrors(nPoints,xValue,yValueShifted,xError,yError);
    returnGraph->SetName(tempName);
    return returnGraph;
}

TGraphErrors* ShiftSpectraWithSlopeSyst(TGraphErrors* centralValue, TGraphAsymmErrors* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetN()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString                      = "";
    if (mode == 1)      modeString                      = "lin";
    else if (mode == 2) modeString                      = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TString     tempName                                = "";
    if (shiftA) tempName                                = Form("%s_%sShiftA", centralValue->GetName(), modeString.Data());
    else        tempName                                = Form("%s_%sShiftB", centralValue->GetName(), modeString.Data());

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xError                                    = centralValue->GetEX();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yError                                    = centralValue->GetEY();

    Double_t* yErrorSysUp                               = sysErr->GetEYhigh();
    Double_t* yErrorSysLow                              = sysErr->GetEYlow();

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                                  = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                                  = 0.;

    Double_t factor                                     = 0.;
    Double_t effectiveError                             = 0.;
    for (Int_t i = 0; i < nPoints; i++) {
        if (!yValue[i]) continue;
        factor                                          = shiftWeight->Eval(xValue[i]);
        if (factor > 0) {
            if (relPtConstSyst) {
                effectiveError                          = yErrorSysUp[i] - relPtConstSyst*yValue[i];
                effectiveError                          = factor * effectiveError;
            } else {
                effectiveError                          = factor * yErrorSysUp[i];
            }
        } else {
            if (relPtConstSyst) {
                effectiveError                          = yErrorSysLow[i] - relPtConstSyst*yValue[i];
                effectiveError                          = factor * effectiveError;
            } else {
                effectiveError                          = factor * yErrorSysLow[i];
            }
        }
        yValueShifted[i]                                = yValue[i] + effectiveError;
    }

    TGraphErrors* returnGraph                           = new TGraphErrors(nPoints,xValue,yValueShifted,xError,yError);
    returnGraph->SetName(tempName);
    return returnGraph;
}

// ___________ TGraphAsymmErrors type input spectra ______________________________________________________________
TGraphAsymmErrors* ShiftSpectraWithSlopeSyst(TGraphAsymmErrors* centralValue, TH1D* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetNbinsX()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString                      = "";
    if (mode == 1)      modeString                      = "lin";
    else if (mode == 2) modeString                      = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TString     tempName                                = "";
    if (shiftA) tempName                                = Form("%s_%sShiftA", centralValue->GetName(), modeString.Data());
    else        tempName                                = Form("%s_%sShiftB", centralValue->GetName(), modeString.Data());

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xErrorLow                                 = centralValue->GetEXlow();
    Double_t* xErrorHigh                                = centralValue->GetEXhigh();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yErrorLow                                 = centralValue->GetEYlow();
    Double_t* yErrorHigh                                = centralValue->GetEYhigh();

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                                  = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                                  = 0.;

    Double_t factor                                     = 0.;
    Double_t effectiveError                             = 0.;
    for (Int_t i = 0; i < nPoints; i++) {
        if (!yValue[i]) continue;
        for (Int_t j=1; j<= sysErr->GetNbinsX(); j++) {
            if (!(xValue[i] == sysErr->GetBinCenter(j))) continue;
            factor                                      = shiftWeight->Eval(xValue[i]);
            if (relPtConstSyst) {
                effectiveError                          = sysErr->GetBinError(j) - relPtConstSyst*yValue[i];
                effectiveError                          = factor * effectiveError;
            } else {
                effectiveError                          = factor * sysErr->GetBinError(j);
            }
            yValueShifted[i]                            = yValue[i] + effectiveError;
        }
    }

    TGraphAsymmErrors* returnGraph                      = new TGraphAsymmErrors(nPoints,xValue,yValueShifted,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    returnGraph->SetName(tempName);
    return returnGraph;
}

TGraphAsymmErrors* ShiftSpectraWithSlopeSyst(TGraphAsymmErrors* centralValue, TGraphErrors* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetN()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString                      = "";
    if (mode == 1)      modeString                      = "lin";
    else if (mode == 2) modeString                      = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TString     tempName                                = "";
    if (shiftA) tempName                                = Form("%s_%sShiftA", centralValue->GetName(), modeString.Data());
    else        tempName                                = Form("%s_%sShiftB", centralValue->GetName(), modeString.Data());

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xErrorLow                                 = centralValue->GetEXlow();
    Double_t* xErrorHigh                                = centralValue->GetEXhigh();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yErrorLow                                 = centralValue->GetEYlow();
    Double_t* yErrorHigh                                = centralValue->GetEYhigh();

    Double_t* yErrorSys                                 = sysErr->GetEY();

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                                  = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                                  = 0.;

    Double_t factor                                     = 0.;
    Double_t effectiveError                             = 0.;
    for (Int_t i = 0; i < nPoints; i++) {
        if (!yValue[i]) continue;
        factor                                          = shiftWeight->Eval(xValue[i]);
        if (relPtConstSyst) {
            effectiveError                              = yErrorSys[i] - relPtConstSyst*yValue[i];
            effectiveError                              = factor * effectiveError;
        } else {
            effectiveError                              = factor * yErrorSys[i];
        }
        yValueShifted[i]                                = yValue[i] + effectiveError;
    }

    TGraphAsymmErrors* returnGraph                      = new TGraphAsymmErrors(nPoints,xValue,yValueShifted,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    returnGraph->SetName(tempName);
    return returnGraph;
}

TGraphAsymmErrors* ShiftSpectraWithSlopeSyst(TGraphAsymmErrors* centralValue, TGraphAsymmErrors* sysErr, Bool_t shiftA, Double_t relPtConstSyst, TF1* shiftWeight, Int_t mode){

    if (!centralValue) {
        cout << "Warning: Stat graph is null!" << endl;
        return NULL;
    }
    if (!centralValue->GetN() || !sysErr->GetN()) return NULL;

    if (!shiftWeight) {
        cout << "Warning: Slope shift is null!" << endl;
        return NULL;
    }

    TString             modeString                      = "";
    if (mode == 1)      modeString                      = "lin";
    else if (mode == 2) modeString                      = "pol2";
    else {
        cout << "Warning: Mode " << mode << " not implemented for slope shift!" << endl;
        return NULL;
    }

    TString     tempName                                = "";
    if (shiftA) tempName                                = Form("%s_%sShiftA", centralValue->GetName(), modeString.Data());
    else        tempName                                = Form("%s_%sShiftB", centralValue->GetName(), modeString.Data());

    Int_t nPoints                                       = centralValue->GetN();
    Double_t* xValue                                    = centralValue->GetX();
    Double_t* xErrorLow                                 = centralValue->GetEXlow();
    Double_t* xErrorHigh                                = centralValue->GetEXhigh();
    Double_t* yValue                                    = centralValue->GetY();
    Double_t* yValueShifted                             = new Double_t[nPoints];
    for (Int_t i=0; i<nPoints; i++) yValueShifted[i]    = yValue[i];
    Double_t* yErrorLow                                 = centralValue->GetEYlow();
    Double_t* yErrorHigh                                = centralValue->GetEYhigh();

    Double_t* yErrorSysUp                               = sysErr->GetEYhigh();
    Double_t* yErrorSysLow                              = sysErr->GetEYlow();

    // check for size of relPtConstSyst
    if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 0)
        relPtConstSyst                                  = 0.;
    else if (UseRelPtConstSystErr(sysErr, relPtConstSyst) == 2)
        relPtConstSyst                                  = 0.;

    Double_t factor                                     = 0.;
    Double_t effectiveError                             = 0.;
    for (Int_t i = 0; i < nPoints; i++) {
        if (!yValue[i]) continue;
        factor                                          = shiftWeight->Eval(xValue[i]);
        if (factor > 0) {
            if (relPtConstSyst) {
                effectiveError                          = yErrorSysUp[i] - relPtConstSyst*yValue[i];
                effectiveError                          = factor * effectiveError;
            } else {
                effectiveError                          = factor * yErrorSysUp[i];
            }
        } else {
            if (relPtConstSyst) {
                effectiveError                          = yErrorSysLow[i] - relPtConstSyst*yValue[i];
                effectiveError                          = factor * effectiveError;
            } else {
                effectiveError                          = factor * yErrorSysLow[i];
            }
        }
        yValueShifted[i]                                = yValue[i] + effectiveError;
    }

    TGraphAsymmErrors* returnGraph                      = new TGraphAsymmErrors(nPoints,xValue,yValueShifted,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    returnGraph->SetName(tempName);
    return returnGraph;
}

//================================================================================================================
//Function that calculates the factors for non-constant systematic shift
//================================================================================================================
TF1* CalculateNonConstantShiftFactor(TH1D* hist, Int_t mode, Bool_t shiftA) {

    if (!hist) {
        cout << "CalculateNonConstantShiftFactor: No input!" << endl;
        return NULL;
    }

    TString     name            = hist->GetName();
    Double_t    ptMin           = 0.;
    for (Int_t i=1; i<hist->GetNbinsX()+1; i++) {
        if (hist->GetBinContent(i)) {
            ptMin               = hist->GetBinCenter(i);
            break;
        }
    }
    Double_t    ptMax           = 0.;
    for (Int_t i=hist->GetNbinsX(); i>=1; i--) {
        if (hist->GetBinContent(i)) {
            ptMax               = hist->GetBinCenter(i);
            break;
        }
    }

    TString     var             = "";
    Double_t    sign            = 1.;
    if (shiftA) {
        var                     = "A";
        sign                    = 1.;
    } else {
        var                     = "B";
        sign                    = -1.;
    }

    TF1* func                   = NULL;
    if (mode == 1) {
        // lin
        Double_t deltaPt        = ptMax - ptMin;
        Double_t slope          = 2 / deltaPt;
        Double_t intercept      = 1 + slope*ptMin;
        func                    = new TF1(Form("%s_NonConstantSystFactorLin%s", name.Data(), var.Data()), "[0]-[1]*x", ptMin, ptMax);
        func->SetParameter(0, sign*intercept);
        func->SetParameter(1, sign*slope);
    }
    else if (mode == 2) {
        // pol2
        Double_t    par1        = ptMin + (ptMax - ptMin)/2;
        Double_t    par0        = 0.;
        if (shiftA) par0        = 2/(ptMax - par1)/(ptMax - par1);
        else        par0        = -2/(ptMax - par1)/(ptMax - par1);
        Double_t    par2        = 0.;
        if (shiftA) par2        = -1.;
        else        par2        = 1.;
        func                    = new TF1(Form("%s_NonConstantSystFactorPol2%s", name.Data(), var.Data()), "[0]*(x-[1])*(x-[1])+[2]", ptMin, ptMax);
        func->SetParameter(0, par0);
        func->SetParameter(1, par1);
        func->SetParameter(2, par2);
    }
    else
    {
        cout << "Warning: Mode " << mode << " not implemented!" << endl;
        return NULL;
    }

    return func;
}

TF1* CalculateNonConstantShiftFactor(TGraphErrors* graph, Int_t mode, Bool_t shiftA) {

    if (!graph) {
        cout << "CalculateNonConstantShiftFactor: No input!" << endl;
        return NULL;
    }

    TString     name            = graph->GetName();
    Int_t       nPoints         = graph->GetN();
    Double_t*   xValue          = graph->GetX();
    Double_t*   yValue          = graph->GetY();
    Double_t    ptMin           = 0.;
    for (Int_t i=0; i<nPoints; i++) {
        if (yValue[i]) {
            ptMin               = xValue[i];
            break;
        }
    }
    Double_t    ptMax           = 0.;
    for (Int_t i=nPoints-1; i>=0; i--) {
        if (yValue[i]) {
            ptMax               = xValue[i];
            break;
        }
    }

    TString     var             = "";
    Double_t    sign            = 1.;
    if (shiftA) {
        var                     = "A";
        sign                    = 1.;
    } else {
        var                     = "B";
        sign                    = -1.;
    }

    TF1* func                   = NULL;
    if (mode == 1) {
        // lin
        Double_t deltaPt        = ptMax - ptMin;
        Double_t slope          = 2 / deltaPt;
        Double_t intercept      = 1 + slope*ptMin;
        func                    = new TF1(Form("%s_NonConstantSystFactorLin%s", name.Data(), var.Data()), "[0]-[1]*x", ptMin, ptMax);
        func->SetParameter(0, sign*intercept);
        func->SetParameter(1, sign*slope);
    }
    else if (mode == 2) {
        // pol2
        Double_t    par1        = ptMin + (ptMax - ptMin)/2;
        Double_t    par0        = 0.;
        if (shiftA) par0        = 2/(ptMax - par1)/(ptMax - par1);
        else        par0        = -2/(ptMax - par1)/(ptMax - par1);
        Double_t    par2        = 0.;
        if (shiftA) par2        = -1.;
        else        par2        = 1.;
        func                    = new TF1(Form("%s_NonConstantSystFactorPol2%s", name.Data(), var.Data()), "[0]*(x-[1])*(x-[1])+[2]", ptMin, ptMax);
        func->SetParameter(0, par0);
        func->SetParameter(1, par1);
        func->SetParameter(2, par2);
    }
    else
    {
        cout << "Warning: Mode " << mode << " not implemented!" << endl;
        return NULL;
    }

    return func;
}

TF1* CalculateNonConstantShiftFactor(TGraphAsymmErrors* graph, Int_t mode, Bool_t shiftA) {

    if (!graph) {
        cout << "CalculateNonConstantShiftFactor: No input!" << endl;
        return NULL;
    }

    TString     name            = graph->GetName();
    Int_t       nPoints         = graph->GetN();
    Double_t*   xValue          = graph->GetX();
    Double_t*   yValue          = graph->GetY();
    Double_t    ptMin           = 0.;
    for (Int_t i=0; i<nPoints; i++) {
        if (yValue[i]) {
            ptMin               = xValue[i];
            break;
        }
    }
    Double_t    ptMax           = 0.;
    for (Int_t i=nPoints-1; i>=0; i--) {
        if (yValue[i]) {
            ptMax               = xValue[i];
            break;
        }
    }

    TString     var             = "";
    Double_t    sign            = 1.;
    if (shiftA) {
        var                     = "A";
        sign                    = 1.;
    } else {
        var                     = "B";
        sign                    = -1.;
    }

    TF1* func                   = NULL;
    if (mode == 1) {
        // lin
        Double_t deltaPt        = ptMax - ptMin;
        Double_t slope          = 2 / deltaPt;
        Double_t intercept      = 1 + slope*ptMin;
        func                    = new TF1(Form("%s_NonConstantSystFactorLin%s", name.Data(), var.Data()), "[0]-[1]*x", ptMin, ptMax);
        func->SetParameter(0, sign*intercept);
        func->SetParameter(1, sign*slope);
    }
    else if (mode == 2) {
        // pol2
        Double_t    par1        = ptMin + (ptMax - ptMin)/2;
        Double_t    par0        = 0.;
        if (shiftA) par0        = 2/(ptMax - par1)/(ptMax - par1);
        else        par0        = -2/(ptMax - par1)/(ptMax - par1);
        Double_t    par2        = 0.;
        if (shiftA) par2        = -1.;
        else        par2        = 1.;
        func                    = new TF1(Form("%s_NonConstantSystFactorPol2%s", name.Data(), var.Data()), "[0]*(x-[1])*(x-[1])+[2]", ptMin, ptMax);
        func->SetParameter(0, par0);
        func->SetParameter(1, par1);
        func->SetParameter(2, par2);
    }
    else
    {
        cout << "Warning: Mode " << mode << " not implemented!" << endl;
        return NULL;
    }

    return func;
}


//================================================================================================================
//Function that produces ratio to fit histogram
//================================================================================================================
TH1D* CalculateRatioToFit(TH1D* histo, TF1* fit) {

    // set name
    TString name                        = Form("%s_RatioToFit", histo->GetName());

    // fit range
    Double_t ptMin                      = 0;
    Double_t ptMax                      = 0;
    fit->GetRange(ptMin, ptMax);

    // ratio to fit
    TH1D* histoRatioToFit               = (TH1D*)histo->Clone(name);
    histoRatioToFit->Sumw2();
    histoRatioToFit->Divide(fit);

    // set bins outside of function range to zero
    for (Int_t i=0; i<histoRatioToFit->GetNbinsX()+1; i++) {
        if (histoRatioToFit->GetXaxis()->GetBinUpEdge(i) <= ptMin || histoRatioToFit->GetXaxis()->GetBinLowEdge(i) >= ptMax) {
            histoRatioToFit->SetBinContent(i,   0);
            histoRatioToFit->SetBinError(i,     0);
        }
    }

    return histoRatioToFit;
}

TH1D* CalculateRatioToFit(TGraphErrors* graph, TF1* fit) {

    // this method integrates the fit function and and the graph in pt bins
    // to calculate the ratio to account for bin shifted spectra

    // set name
    TString name                        = Form("%s_RatioToFit", graph->GetName());

    // graph
    Int_t nPoints                       = graph->GetN();
    Double_t* xValue                    = graph->GetX();
    Double_t* xError                    = graph->GetEX();
    Double_t* yValue                    = graph->GetY();
    Double_t* yError                    = graph->GetEY();

    // fit
    Double_t ptMin                      = 0;
    Double_t ptMax                      = 0;
    fit->GetRange(ptMin, ptMax);
    ptMin                               = ptMin - 1e-12;     // dirty bugfix, selection of bins outside pt-range of fit broke sometimes for no obvious reason

    // ratio histogram
    Double_t* binning                   = new Double_t[nPoints+2];
    for (Int_t i=0; i<nPoints+2; i++) {
        if (i==0)
            binning[i]                  = 0;
        else if (i<nPoints)
            binning[i]                  = xValue[i-1] - xError[i-1];
        else
            binning[i]                  = xValue[i-2] + xError[i-2];
    }
    TH1D* dummyhistoRatioToFit               = new TH1D("dummy", "", nPoints+1, binning);
    TH1D* histoRatioToFit                    = (TH1D*)dummyhistoRatioToFit->Clone(name);

    Double_t tempGraphVal               = 0;
    Double_t tempGraphErr               = 0;
    Double_t tempFitVal                 = 0;
    Double_t tempFitErr                 = 0;
    Double_t tempBinContent             = 0;
    Double_t tempBinError               = 0;

    // for integrating the fit
    Int_t np                            = 1000;
    Double_t *x                         = new Double_t[np];
    Double_t *w                         = new Double_t[np];

    for (Int_t i=2; i<histoRatioToFit->GetNbinsX()+1; i++) {
        // integrate fit in pt bin
        if (histoRatioToFit->GetXaxis()->GetBinLowEdge(i) < ptMin ||  histoRatioToFit->GetXaxis()->GetBinUpEdge(i) > ptMax) {
            tempFitVal                  = 0;
            tempFitErr                  = 0;
        } else {
            fit->CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            tempFitVal                  = fit->IntegralFast(np, x, w, xValue[i-2] - xError[i-2], xValue[i-2] + xError[i-2]);
            //tempFitErr                  = fit->IntegralError(xValue[i-1] - xError[i-1], xValue[i-1] + xError[i-1]);  // needs covariance matrix, but maybe not even reasonable for ratio to fit
            tempFitErr                  = 0;
        }

        // integrate graph in bin
        tempGraphVal                    = yValue[i-2] / (2*xError[i-2]);
        tempGraphErr                    = yError[i-2] * (2*xError[i-2]);

        // calculate bin content and error
        if (tempFitVal) {
            tempBinContent              = tempGraphVal/tempFitVal;
            tempBinError                = TMath::Power(tempGraphErr, 2) / TMath::Power(tempFitVal, 2) + TMath::Power(tempGraphVal, 2) * TMath::Power(tempFitErr, 2) / TMath::Power(tempFitVal, 4);
        } else {
            tempBinContent              = 0;
            tempBinError                = 0;
        }

        // fill historam
        histoRatioToFit->SetBinContent(i,   tempBinContent);
        histoRatioToFit->SetBinError(i,     TMath::Sqrt(tempBinError));
    }

    delete dummyhistoRatioToFit;
    return histoRatioToFit;
}

TH1D* CalculateRatioToFit(TGraphAsymmErrors* graph, TF1* fit) {

    // this method integrates the fit function and and the graph in pt bins
    // to calculate the ratio to account for bin shifted spectra

    // set name
    TString name                        = Form("%s_RatioToFit", graph->GetName());

    // graph
    Int_t nPoints                       = graph->GetN();
    Double_t* xValue                    = graph->GetX();
    Double_t* xErrorLow                 = graph->GetEXlow();
    Double_t* xErrorHigh                = graph->GetEXhigh();
    Double_t* yValue                    = graph->GetY();
    Double_t* yErrorLow                 = graph->GetEYlow();
    Double_t* yErrorHigh                = graph->GetEYhigh();

    // fit
    Double_t ptMin                      = 0;
    Double_t ptMax                      = 0;
    fit->GetRange(ptMin, ptMax);
    ptMin                               = ptMin - 1e-12;     // dirty bugfix, selection of bins outside pt-range of fit broke sometimes for no obvious reason

    // ratio histogram
    Double_t* binning                   = new Double_t[nPoints+2];
    for (Int_t i=0; i<nPoints+2; i++) {
        if (i==0)
            binning[i]                  = 0;
        else if (i<nPoints)
            binning[i]                  = xValue[i-1] - xErrorLow[i-1];
        else
            binning[i]                  = xValue[i-2] + xErrorHigh[i-2];
    }
    TH1D* dummyhistoRatioToFit               = new TH1D("dummy", "", nPoints+1, binning);
    TH1D* histoRatioToFit                    = (TH1D*)dummyhistoRatioToFit->Clone(name);

    Double_t tempGraphVal               = 0;
    Double_t tempGraphErr               = 0;
    Double_t tempFitVal                 = 0;
    Double_t tempFitErr                 = 0;
    Double_t tempBinContent             = 0;
    Double_t tempBinError               = 0;

    // for integrating the fit
    Int_t np                            = 1000;
    Double_t *x                         = new Double_t[np];
    Double_t *w                         = new Double_t[np];

    for (Int_t i=2; i<histoRatioToFit->GetNbinsX()+1; i++) {

        // integrate fit in pt bin
        if (histoRatioToFit->GetXaxis()->GetBinLowEdge(i) < ptMin ||  histoRatioToFit->GetXaxis()->GetBinUpEdge(i) > ptMax) {
            tempFitVal                  = 0;
            tempFitErr                  = 0;
        } else {
            fit->CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            tempFitVal                  = fit->IntegralFast(np, x, w, xValue[i-2] - xErrorLow[i-2], xValue[i-2] + xErrorHigh[i-2]);
            //tempFitErr                  = fit->IntegralError(xValue[i-1] - xErrorLow[i-1], xValue[i-1] + xErrorHigh[i-1]);  // needs covariance matrix, but maybe not even reasonable for ratio to fit
            tempFitErr                  = 0;
        }

        // integrate graph in bin
        tempGraphVal                    = yValue[i-2] * (xErrorLow[i-2] + xErrorHigh[i-2]);
        tempGraphErr                    = (yErrorLow[i-2] + yErrorHigh[i-2]) / 2 * (xErrorLow[i-2] + xErrorHigh[i-2]);

        // calculate bin content and error
        if (tempFitVal) {
            tempBinContent              = tempGraphVal/tempFitVal;
            tempBinError                = TMath::Power(tempGraphErr, 2) / TMath::Power(tempFitVal, 2) + TMath::Power(tempGraphVal, 2) * TMath::Power(tempFitErr, 2) / TMath::Power(tempFitVal, 4);
        } else {
            tempBinContent              = 0;
            tempBinError                = 0;
        }

        // fill historam
        histoRatioToFit->SetBinContent(i,   tempBinContent);
        histoRatioToFit->SetBinError(i,     TMath::Sqrt(tempBinError));
    }

    delete dummyhistoRatioToFit;
    return histoRatioToFit;
}

//================================================================================================================
//Function that produces spectra in certain cent. class
//================================================================================================================
void ProduceSpectrumInCentralityBin(TList* list[nCentralities], TString desiredCentrality, TString particle, TString particle2, TString method, Double_t scaleFacNSD = 1 ) {

    // input objects
    TObject* tempInputObjStat[nCentralities];
    TObject* tempInputObjSys[nCentralities];

    cout << "trying to combine \t "<< desiredCentrality.Data() << " for \t" << particle.Data() << "\t" << method.Data() << endl;
    // get input objects
    for (Int_t i=0; i<nCentralities; i++) {

        TString tempListName = "";
        if (list[i]) tempListName = list[i]->GetName();
        if (list[i] && !(tempListName.Contains(desiredCentrality))) {

            if (particle2.CompareTo("") == 0) {
                // spectra
                if (list[i]->FindObject(Form("%s%sStat", particle.Data(), method.Data())) && list[i]->FindObject(Form("%s%sSys", particle.Data(), method.Data()))) {
                    tempInputObjStat[i]             = (TObject*)list[i]->FindObject(Form("%s%sStat", particle.Data(), method.Data()));
                    tempInputObjSys[i]              = (TObject*)list[i]->FindObject(Form("%s%sSys", particle.Data(), method.Data()));
                } else {
                    tempInputObjStat[i]             = NULL;
                    tempInputObjSys[i]              = NULL;
                }
            } else {
                // ratios
                if (list[i]->FindObject(Form("%sTo%s%sStat", particle.Data(), particle2.Data(), method.Data())) && list[i]->FindObject(Form("%sTo%s%sSys", particle.Data(), particle2.Data(), method.Data()))) {
                    tempInputObjStat[i]             = (TObject*)list[i]->FindObject(Form("%sTo%s%sStat", particle.Data(), particle2.Data(), method.Data()));
                    tempInputObjSys[i]              = (TObject*)list[i]->FindObject(Form("%sTo%s%sSys", particle.Data(), particle2.Data(), method.Data()));
                } else {
                    tempInputObjStat[i]             = NULL;
                    tempInputObjSys[i]              = NULL;
                }
            }
        } else {
            tempInputObjStat[i]                     = NULL;
            tempInputObjSys[i]                      = NULL;
        }
        cout << "next cent" << endl;
    }

    // get class name
    TString ClassName;
    for (Int_t i=0; i<nCentralities; i++) {
        if (tempInputObjStat[i] && tempInputObjSys[i]) {
            ClassName                               = (TString)tempInputObjStat[i]->ClassName();
            break;
        }
    }

    // create output object according to input type
    TH1D* tempOutputHistStat                        = NULL;
    TH1D* tempOutputHistSys                         = NULL;
    TGraphErrors* tempOutputGraphStat               = NULL;
    TGraphErrors* tempOutputGraphSys                = NULL;
    TGraphAsymmErrors* tempOutputGraphAsymStat      = NULL;
    TGraphAsymmErrors* tempOutputGraphAsymSys       = NULL;

    Bool_t writeToList                              = kTRUE;
    Double_t* sysErrsPerCent[nCentralities];

    // calculating MB
    if (desiredCentrality == fCentrality[0] || desiredCentrality == fCentrality[16]) {
        if (tempInputObjStat[1] && tempInputObjStat[2] && tempInputObjStat[4] && tempInputObjStat[6] && tempInputObjStat[8] && tempInputObjStat[9] && tempInputObjStat[11] && !tempInputObjStat[3]) {
            // 0-5 + 5-10 + 10-20 + 20-40 + 40-60 + 60-80 + 80-100

            if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[1]->Clone(tempInputObjStat[1]->GetName());
                tempOutputHistStat->Scale(0.05);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[2], 0.05);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[4], 0.1);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[6], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[8], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[9], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[11], 0.2);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[1]->Clone(tempInputObjSys[1]->GetName());
                tempOutputHistSys->Scale(0.05);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[2], 0.05);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[4], 0.1);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[6], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[8], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[9], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[11], 0.2);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(7);
                    if (((TH1D*)tempInputObjSys[1])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[1])->GetBinError(i)/((TH1D*)tempInputObjSys[1])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[2])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[2])->GetBinError(i)/((TH1D*)tempInputObjSys[2])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[4])->GetBinContent(i)) sysErrsVector[2]     = ((TH1D*)tempInputObjSys[4])->GetBinError(i)/((TH1D*)tempInputObjSys[4])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[6])->GetBinContent(i)) sysErrsVector[3]     = ((TH1D*)tempInputObjSys[6])->GetBinError(i)/((TH1D*)tempInputObjSys[6])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[8])->GetBinContent(i)) sysErrsVector[4]     = ((TH1D*)tempInputObjSys[8])->GetBinError(i)/((TH1D*)tempInputObjSys[8])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[9])->GetBinContent(i)) sysErrsVector[5]     = ((TH1D*)tempInputObjSys[9])->GetBinError(i)/((TH1D*)tempInputObjSys[9])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[11])->GetBinContent(i)) sysErrsVector[6]    = ((TH1D*)tempInputObjSys[11])->GetBinError(i)/((TH1D*)tempInputObjSys[11])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[7]                                     = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[7]                                      = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
                Int_t intOrCents[7]                                                         = {1, 2, 4, 6, 8, 9, 11};
                Double_t scaleFacs[7]                                                       = {0.05, 0.05, 0.1, 0.2, 0.2, 0.2, 0.2};
                for (Int_t k = 0; k < 7; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 7, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 7, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[7]                                = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[7]                                 = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
                Int_t intOrCents[7]                                                         = {1, 2, 4, 6, 8, 9, 11};
                Double_t scaleFacs[7]                                                       = {0.05, 0.05, 0.1, 0.2, 0.2, 0.2, 0.2};
                for (Int_t k = 0; k < 7; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 7, 0);
                if (tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 7, 1);
                if (tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else if (tempInputObjStat[3] && tempInputObjStat[4] && tempInputObjStat[6] && tempInputObjStat[8] && tempInputObjStat[9] && tempInputObjStat[11] && !tempInputObjStat[5]) {
            // 0-10 + 10-20 + 20-40 + 40-60 + 60-80 + 80-100

            if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[3]->Clone(tempInputObjStat[3]->GetName());
                tempOutputHistStat->Scale(0.1);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[4], 0.1);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[6], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[8], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[9], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[11], 0.2);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[3]->Clone(tempInputObjSys[3]->GetName());
                tempOutputHistSys->Scale(0.1);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[4], 0.1);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[6], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[8], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[9], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[11], 0.2);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(6);
                    if (((TH1D*)tempInputObjSys[3])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[3])->GetBinError(i)/((TH1D*)tempInputObjSys[3])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[4])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[4])->GetBinError(i)/((TH1D*)tempInputObjSys[4])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[6])->GetBinContent(i)) sysErrsVector[2]     = ((TH1D*)tempInputObjSys[6])->GetBinError(i)/((TH1D*)tempInputObjSys[6])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[8])->GetBinContent(i)) sysErrsVector[3]     = ((TH1D*)tempInputObjSys[8])->GetBinError(i)/((TH1D*)tempInputObjSys[8])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[9])->GetBinContent(i)) sysErrsVector[4]     = ((TH1D*)tempInputObjSys[9])->GetBinError(i)/((TH1D*)tempInputObjSys[9])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[11])->GetBinContent(i)) sysErrsVector[5]    = ((TH1D*)tempInputObjSys[11])->GetBinError(i)/((TH1D*)tempInputObjSys[11])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[6]                                     = {NULL, NULL, NULL, NULL, NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[6]                                      = {NULL, NULL, NULL, NULL, NULL, NULL};
                Int_t intOrCents[6]                                                         = {3, 4, 6, 8, 9, 11};
                Double_t scaleFacs[6]                                                       = {0.1, 0.1, 0.2, 0.2, 0.2, 0.2};
                for (Int_t k = 0; k < 6; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 6, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 6, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[6]                                = {NULL, NULL, NULL, NULL, NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[6]                                 = {NULL, NULL, NULL, NULL, NULL, NULL};
                Int_t intOrCents[6]                                                         = {3, 4, 6, 8, 9, 11};
                Double_t scaleFacs[6]                                                       = {0.1, 0.1, 0.2, 0.2, 0.2, 0.2};
                for (Int_t k = 0; k < 6; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 6, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 6, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else if (tempInputObjStat[5] && tempInputObjStat[6] && tempInputObjStat[8] && tempInputObjStat[9] && tempInputObjStat[11]) {
            // 0-20 + 20-40 + 40-60 + 60-80 + 80-100

            if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[5]->Clone(tempInputObjStat[5]->GetName());
                tempOutputHistStat->Scale(0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[6], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[8], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[9], 0.2);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[11], 0.2);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[5]->Clone(tempInputObjSys[5]->GetName());
                tempOutputHistSys->Scale(0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[6], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[8], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[9], 0.2);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[11], 0.2);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(5);
                    if (((TH1D*)tempInputObjSys[5])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[5])->GetBinError(i)/((TH1D*)tempInputObjSys[5])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[6])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[6])->GetBinError(i)/((TH1D*)tempInputObjSys[6])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[8])->GetBinContent(i)) sysErrsVector[2]     = ((TH1D*)tempInputObjSys[8])->GetBinError(i)/((TH1D*)tempInputObjSys[8])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[9])->GetBinContent(i)) sysErrsVector[3]     = ((TH1D*)tempInputObjSys[9])->GetBinError(i)/((TH1D*)tempInputObjSys[9])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[11])->GetBinContent(i)) sysErrsVector[4]    = ((TH1D*)tempInputObjSys[11])->GetBinError(i)/((TH1D*)tempInputObjSys[11])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[5]                                     = {NULL, NULL, NULL, NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[5]                                      = {NULL, NULL, NULL, NULL, NULL};
                Int_t intOrCents[5]                                                         = {5, 6, 8, 9, 11};
                Double_t scaleFacs[5]                                                       = {0.2, 0.2, 0.2, 0.2, 0.2};
                for (Int_t k = 0; k < 5; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 5, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 5, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[5]                                = {NULL, NULL, NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[5]                                 = {NULL, NULL, NULL, NULL};
                Int_t intOrCents[5]                                                         = {5, 6, 8, 9, 11};
                Double_t scaleFacs[5]                                                       = {0.2, 0.2, 0.2, 0.2, 0.2};
                for (Int_t k = 0; k < 5; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 5, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 5, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        }  else {
            cout << "\trequired spectra not found" << endl;
            writeToList                                                                     = kFALSE;
        }

        if (writeToList) {
            if (desiredCentrality == fCentrality[0]){
                if (ClassName.BeginsWith("TH1")) {
                    if (tempOutputHistStat) list[0]->Add(tempOutputHistStat);
                    if (tempOutputHistSys) list[0]->Add(tempOutputHistSys);
                }
                if (ClassName.CompareTo("TGraphErrors") == 0) {
                    if (tempOutputGraphStat) list[0]->Add(tempOutputGraphStat);
                    if (tempOutputGraphSys) list[0]->Add(tempOutputGraphSys);
                }
                if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                    if (tempOutputGraphAsymStat) list[0]->Add(tempOutputGraphAsymStat);
                    if (tempOutputGraphAsymStat) list[0]->Add(tempOutputGraphAsymSys);
                }
            } else if (desiredCentrality == fCentrality[16]){
                if (ClassName.BeginsWith("TH1")) {
                    if (tempOutputHistStat){
                        tempOutputHistStat->Scale(scaleFacNSD);
                        list[16]->Add(tempOutputHistStat);
                    }
                    if (tempOutputHistSys){
                        tempOutputHistSys->Scale(scaleFacNSD);
                        list[16]->Add(tempOutputHistSys);
                    }
                }
                if (ClassName.CompareTo("TGraphErrors") == 0) {
                    if (tempOutputGraphStat){
                        tempOutputGraphStat     = ScaleGraph(tempOutputGraphStat, scaleFacNSD);
                        list[16]->Add(tempOutputGraphStat);
                    }
                    if (tempOutputGraphSys){
                        tempOutputGraphSys      = ScaleGraph(tempOutputGraphSys, scaleFacNSD);
                        list[16]->Add(tempOutputGraphSys);
                    }
                }
                if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                    if (tempOutputGraphStat){
                        tempOutputGraphStat     = ScaleGraph(tempOutputGraphStat, scaleFacNSD);
                        list[16]->Add(tempOutputGraphAsymStat);
                    }
                    if (tempOutputGraphSys){
                        tempOutputGraphSys      = ScaleGraph(tempOutputGraphSys, scaleFacNSD);
                        list[16]->Add(tempOutputGraphAsymSys);
                    }
                }
            }
        }
    }

    // claculating 0-20%
    if (desiredCentrality == fCentrality[5]) {
        if (tempInputObjStat[1] && tempInputObjStat[2] && tempInputObjStat[4] && !tempInputObjStat[3]) {
            // 0-5 + 5-10 + 10-20

            if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[1]->Clone(tempInputObjStat[1]->GetName());
                tempOutputHistStat->Scale(0.25);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[2], 0.25);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[4], 0.5);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[1]->Clone(tempInputObjSys[1]->GetName());
                tempOutputHistSys->Scale(0.25);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[2], 0.25);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[4], 0.5);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(3);
                    if (((TH1D*)tempInputObjSys[1])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[1])->GetBinError(i)/((TH1D*)tempInputObjSys[1])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[2])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[2])->GetBinError(i)/((TH1D*)tempInputObjSys[2])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[4])->GetBinContent(i)) sysErrsVector[2]     = ((TH1D*)tempInputObjSys[4])->GetBinError(i)/((TH1D*)tempInputObjSys[4])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[3]                                     = {NULL, NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[3]                                      = {NULL, NULL, NULL};
                Int_t intOrCents[3]                                                         = {1, 2, 4};
                Double_t scaleFacs[3]                                                       = {0.25, 0.25, 0.5};
                for (Int_t k = 0; k < 3; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 3, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 3, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[3]                                = {NULL, NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[3]                                 = {NULL, NULL, NULL};
                Int_t intOrCents[3]                                                         = {1, 2, 4};
                Double_t scaleFacs[3]                                                       = {0.25, 0.25, 0.5};
                for (Int_t k = 0; k < 3; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 3, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 3, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else if (tempInputObjStat[3] && tempInputObjStat[4]) {
            // 0-10 + 10-20
            if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[3]->Clone(tempInputObjStat[3]->GetName());
                tempOutputHistStat->Scale(0.5);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[4], 0.5);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[3]->Clone(tempInputObjSys[3]->GetName());
                tempOutputHistSys->Scale(0.5);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[4], 0.5);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(2);
                    if (((TH1D*)tempInputObjSys[3])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[3])->GetBinError(i)/((TH1D*)tempInputObjSys[3])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[4])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[4])->GetBinError(i)/((TH1D*)tempInputObjSys[4])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[2]                                     = {NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[2]                                      = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {3, 4};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[2]                                = {NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[2]                                 = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {3, 4};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else {
            cout << "\trequired spectra not found" << endl;
            if (tempInputObjStat[1]) cout << "missing 0-5%" << endl;
            if (tempInputObjStat[2]) cout << "missing 5-10%" << endl;
            if (tempInputObjStat[3]) cout << "missing 0-10%" << endl;
            if (tempInputObjStat[4]) cout << "missing 10-20%" << endl;
            writeToList                                                                     = kFALSE;
        }

        if (writeToList) {
            if (ClassName.BeginsWith("TH1")) {
                list[5]->Add(tempOutputHistStat);
                list[5]->Add(tempOutputHistSys);
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                list[5]->Add(tempOutputGraphStat);
                list[5]->Add(tempOutputGraphSys);
            }
            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                list[5]->Add(tempOutputGraphAsymStat);
                list[5]->Add(tempOutputGraphAsymSys);
            }
        }
    }

    // calcualating 0-10%
    if (desiredCentrality == fCentrality[3]) {
        if (tempInputObjStat[1] && tempInputObjStat[2]) {
            // 0-5 + 5-10

            if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[1]->Clone(tempInputObjStat[1]->GetName());
                tempOutputHistStat->Scale(0.5);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[2], 0.5);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[1]->Clone(tempInputObjSys[1]->GetName());
                tempOutputHistSys->Scale(0.5);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[2], 0.5);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(2);
                    if (((TH1D*)tempInputObjSys[1])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[1])->GetBinError(i)/((TH1D*)tempInputObjSys[1])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[2])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[2])->GetBinError(i)/((TH1D*)tempInputObjSys[2])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[2]                                     = {NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[2]                                      = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {1, 2};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[2]                                = {NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[2]                                 = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {1, 2};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else {
            cout << "\trequired spectra not found" << endl;
            writeToList                                                                     = kFALSE;
        }

        if (writeToList) {
            if (ClassName.BeginsWith("TH1")) {
                list[3]->Add(tempOutputHistStat);
                list[3]->Add(tempOutputHistSys);
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                list[3]->Add(tempOutputGraphStat);
                list[3]->Add(tempOutputGraphSys);
            }
            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                list[3]->Add(tempOutputGraphAsymStat);
                list[3]->Add(tempOutputGraphAsymSys);
            }
        }
    }

    // calcualating 20-40%
    if (desiredCentrality == fCentrality[6]) {
        if (tempInputObjStat[12] && tempInputObjStat[13]) {
            // 20-30 + 30-40

          if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[12]->Clone(tempInputObjStat[12]->GetName());
                tempOutputHistStat->Scale(0.5);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[13], 0.5);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[12]->Clone(tempInputObjSys[12]->GetName());
                tempOutputHistSys->Scale(0.5);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[13], 0.5);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(2);
                    if (((TH1D*)tempInputObjSys[12])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[12])->GetBinError(i)/((TH1D*)tempInputObjSys[12])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[13])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[13])->GetBinError(i)/((TH1D*)tempInputObjSys[13])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }

            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[2]                                     = {NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[2]                                      = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {12, 13};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[2]                                = {NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[2]                                 = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {12, 13};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else {
            cout << "\trequired spectra not found" << endl;
            writeToList                                                                     = kFALSE;
        }



        if (writeToList) {
            if (ClassName.BeginsWith("TH1")) {
                list[6]->Add(tempOutputHistStat);
                list[6]->Add(tempOutputHistSys);
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                list[6]->Add(tempOutputGraphStat);
                list[6]->Add(tempOutputGraphSys);
            }
            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                list[6]->Add(tempOutputGraphAsymStat);
                list[6]->Add(tempOutputGraphAsymSys);
            }
        }
    }

        // calcualating 20-40%
    if (desiredCentrality == fCentrality[6]) {
        if (tempInputObjStat[12] && tempInputObjStat[13]) {
            // 20-30 + 30-40

          if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[12]->Clone(tempInputObjStat[12]->GetName());
                tempOutputHistStat->Scale(0.5);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[13], 0.5);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[12]->Clone(tempInputObjSys[12]->GetName());
                tempOutputHistSys->Scale(0.5);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[13], 0.5);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(2);
                    if (((TH1D*)tempInputObjSys[12])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[12])->GetBinError(i)/((TH1D*)tempInputObjSys[12])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[13])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[13])->GetBinError(i)/((TH1D*)tempInputObjSys[13])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }

            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[2]                                     = {NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[2]                                      = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {12, 13};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[2]                                = {NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[2]                                 = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {12, 13};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else {
            cout << "\trequired spectra not found" << endl;
            writeToList                                                                     = kFALSE;
        }



        if (writeToList) {
            if (ClassName.BeginsWith("TH1")) {
                list[6]->Add(tempOutputHistStat);
                list[6]->Add(tempOutputHistSys);
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                list[6]->Add(tempOutputGraphStat);
                list[6]->Add(tempOutputGraphSys);
            }
            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                list[6]->Add(tempOutputGraphAsymStat);
                list[6]->Add(tempOutputGraphAsymSys);
            }
        }
    }


    // calcualating 60-100%
    if (desiredCentrality == fCentrality[15]) {
        if (tempInputObjStat[9] && tempInputObjStat[11]) {
            // 60-80 & 80-100

            if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[9]->Clone(tempInputObjStat[9]->GetName());
                tempOutputHistStat->Scale(0.5);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[11], 0.5);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[9]->Clone(tempInputObjSys[9]->GetName());
                tempOutputHistSys->Scale(0.5);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[11], 0.5);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(2);
                    if (((TH1D*)tempInputObjSys[9])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[9])->GetBinError(i)/((TH1D*)tempInputObjSys[9])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[11])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[11])->GetBinError(i)/((TH1D*)tempInputObjSys[11])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }

            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[2]                                     = {NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[2]                                      = {NULL, NULL};

                Int_t intOrCents[2]                                                         = {9, 11};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[2]                                = {NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[2]                                 = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {9, 11};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else {
            cout << "\trequired spectra not found" << endl;
            writeToList                                                                     = kFALSE;
        }

        if (writeToList) {
            if (ClassName.BeginsWith("TH1")) {
                list[15]->Add(tempOutputHistStat);
                list[15]->Add(tempOutputHistSys);
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                list[15]->Add(tempOutputGraphStat);
                list[15]->Add(tempOutputGraphSys);
            }
            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                list[15]->Add(tempOutputGraphAsymStat);
                list[15]->Add(tempOutputGraphAsymSys);
            }
        }
	cout << "done" << endl;
    }


    // calcualating 40-60%
    if (desiredCentrality == fCentrality[8]) {
        if (tempInputObjStat[17] && tempInputObjStat[18]) {
            // 40-50 + 50-60

          if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[17]->Clone(tempInputObjStat[17]->GetName());
                tempOutputHistStat->Scale(0.5);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[18], 0.5);

                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[17]->Clone(tempInputObjSys[17]->GetName());
                tempOutputHistSys->Scale(0.5);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[18], 0.5);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(2);
                    if (((TH1D*)tempInputObjSys[17])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[17])->GetBinError(i)/((TH1D*)tempInputObjSys[17])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[18])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[18])->GetBinError(i)/((TH1D*)tempInputObjSys[18])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());
                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }

            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[2]                                     = {NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[2]                                      = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {17, 18};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[2]                                = {NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[2]                                 = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {17, 18};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else {
            cout << "\trequired spectra not found" << endl;
            writeToList                                                                     = kFALSE;
        }

        if (writeToList) {
            if (ClassName.BeginsWith("TH1")) {
                list[8]->Add(tempOutputHistStat);
                list[8]->Add(tempOutputHistSys);
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                list[8]->Add(tempOutputGraphStat);
                list[8]->Add(tempOutputGraphSys);
            }
            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                list[8]->Add(tempOutputGraphAsymStat);
                list[8]->Add(tempOutputGraphAsymSys);
            }
        }
	cout << "done" << endl;
    }

    // calcualating 60-80%
    if (desiredCentrality == fCentrality[9]) {
        if (tempInputObjStat[19] && tempInputObjStat[20]) {
            // 60-70 + 70-80

          if (ClassName.BeginsWith("TH1")) {
                tempOutputHistStat                                                          = (TH1D*)tempInputObjStat[19]->Clone(tempInputObjStat[19]->GetName());
                tempOutputHistStat->Scale(0.5);
                tempOutputHistStat->Add((TH1D*)tempInputObjStat[20], 0.5);
                tempOutputHistSys                                                           = (TH1D*)tempInputObjSys[19]->Clone(tempInputObjSys[19]->GetName());
                tempOutputHistSys->Scale(0.5);
                tempOutputHistSys->Add((TH1D*)tempInputObjSys[20], 0.5);

                // set sys error
                for (Int_t i=1; i<tempOutputHistSys->GetNbinsX()+1; i++) {
                    std::vector<Double_t> sysErrsVector(2);
                    if (((TH1D*)tempInputObjSys[19])->GetBinContent(i)) sysErrsVector[0]     = ((TH1D*)tempInputObjSys[19])->GetBinError(i)/((TH1D*)tempInputObjSys[19])->GetBinContent(i);
                    if (((TH1D*)tempInputObjSys[20])->GetBinContent(i)) sysErrsVector[1]     = ((TH1D*)tempInputObjSys[20])->GetBinError(i)/((TH1D*)tempInputObjSys[20])->GetBinContent(i);
                    std::sort (sysErrsVector.begin(), sysErrsVector.end());

                    tempOutputHistSys->SetBinError(i, sysErrsVector[sysErrsVector.size()-1] * tempOutputHistSys->GetBinContent(i));
                }
            }

            if (ClassName.CompareTo("TGraphErrors") == 0) {
                TGraphErrors* tempOutputGraphStatArr[2]                                     = {NULL, NULL};
                TGraphErrors* tempOutputGraphSysArr[2]                                      = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {19, 20};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphStat                                                         = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphStat) tempOutputGraphStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphSys                                                          = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphSys) tempOutputGraphSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }

            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                TGraphAsymmErrors* tempOutputGraphStatArr[2]                                = {NULL, NULL};
                TGraphAsymmErrors* tempOutputGraphSysArr[2]                                 = {NULL, NULL};
                Int_t intOrCents[2]                                                         = {19, 20};
                Double_t scaleFacs[2]                                                       = {0.5, 0.5};
                for (Int_t k = 0; k < 2; k++){
                    tempOutputGraphStatArr[k]                                               = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphStatArr[k]                                               = ScaleGraph(tempOutputGraphStatArr[k],scaleFacs[k]);
                    tempOutputGraphSysArr[k]                                                = (TGraphAsymmErrors*)tempInputObjStat[intOrCents[k]]->Clone(tempInputObjStat[intOrCents[k]]->GetName());
                    tempOutputGraphSysArr[k]                                                = ScaleGraph(tempOutputGraphSysArr[k],scaleFacs[k]);
                }

                tempOutputGraphAsymStat                                                     = AddNTGraphErrorsSameBinning(tempOutputGraphStatArr, 2, 0);
                if(tempOutputGraphAsymStat) tempOutputGraphAsymStat->SetName(tempInputObjStat[intOrCents[0]]->GetName());

                tempOutputGraphAsymSys                                                      = AddNTGraphErrorsSameBinning(tempOutputGraphSysArr, 2, 1);
                if(tempOutputGraphAsymSys) tempOutputGraphAsymSys->SetName(tempInputObjSys[intOrCents[0]]->GetName());
            }
        } else {
            cout << "\trequired spectra not found" << endl;
            writeToList                                                                     = kFALSE;
        }

        if (writeToList) {
            if (ClassName.BeginsWith("TH1")) {
                list[9]->Add(tempOutputHistStat);
                list[9]->Add(tempOutputHistSys);
            }
            if (ClassName.CompareTo("TGraphErrors") == 0) {
                list[9]->Add(tempOutputGraphStat);
                list[9]->Add(tempOutputGraphSys);
            }
            if (ClassName.CompareTo("TGraphAsymmErrors") == 0) {
                list[9]->Add(tempOutputGraphAsymStat);
                list[9]->Add(tempOutputGraphAsymSys);
            }
        }
	cout << "done" << endl;
    }

}

//================================================================================================================
//Function that produces spectra in certain cent. class
//================================================================================================================
Bool_t ProduceAveragedQuantiesInCentralityBin(  TH1D** histoArrayStat,
                                                TH1D** histoArraySys,
                                                Int_t nCentsToComb,
                                                Int_t* centArrToComb,
                                                Double_t* centWeights,
                                                Int_t desiredCentralityIter,
                                                Int_t particleIter
                                           ) {

    // input objects
    Bool_t bStat[nCentsToComb];
    Bool_t bSys[nCentsToComb];
    Double_t values[nCentsToComb];
    Double_t errorsStat[nCentsToComb];
    Double_t errorsSys[nCentsToComb];
    cout << "trying to combine \t "<< desiredCentralityIter << " for \t" << particleIter << endl;

    // get input hists
    for (Int_t i=0; i<nCentsToComb; i++) {
        bStat[i]                        = kFALSE;
        bSys[i]                         = kFALSE;
        if (histoArrayStat[centArrToComb[i]]->GetBinContent(particleIter+1) != -1){
            bStat[i]                    = kTRUE;
            values[i]                   = histoArrayStat[centArrToComb[i]]->GetBinContent(particleIter+1);
            errorsStat[i]               = histoArrayStat[centArrToComb[i]]->GetBinError(particleIter+1);
        }
        if (histoArraySys[centArrToComb[i]]->GetBinContent(particleIter+1) != -1){
            bSys[i]                    = kTRUE;
            errorsSys[i]               = histoArraySys[centArrToComb[i]]->GetBinError(particleIter+1);
        }
        if ( !(bStat[i] && bSys[i] )){
            cout << "desired particle didn't have values filled for cent " << centArrToComb[i] << endl;
            return kFALSE;
        }
        cout << centArrToComb[i] << "\t\t: " << values[i] << "\t +- " <<  errorsStat[i] << "\t +- " << errorsSys[i] << endl;
    }
    Double_t value                      = 0;
    Double_t errStat                    = 0;
    Double_t errSys                     = 0;

    for (Int_t i = 0; i < nCentsToComb; i++){
        value                           = value + centWeights[i]*values[i];
        errStat                         = errStat + TMath::Power(centWeights[i]*errorsStat[i],2);
        errSys                          = errSys + TMath::Power(centWeights[i]*errorsSys[i],2);
    }
    errStat                             = TMath::Sqrt(errStat);
    errSys                              = TMath::Sqrt(errSys);

    cout << value << "\t +- " <<  errStat << "\t +- " << errSys << endl;

    histoArrayStat[desiredCentralityIter]->SetBinContent(particleIter+1, value);
    histoArrayStat[desiredCentralityIter]->SetBinError(errStat+1, value);
    histoArraySys[desiredCentralityIter]->SetBinContent(particleIter+1, value);
    histoArraySys[desiredCentralityIter]->SetBinError(errSys+1, value);

    return kTRUE;
}


//================================================================================================================
//Function to process pt-y histograms
//================================================================================================================
void NormalizePtYHistogram(TH2F* hist, Int_t nRebinX, Int_t nRebinY, Double_t thresh) {

    if (!hist) return;

    // rebin histogram
    hist->Rebin2D(nRebinX, nRebinY);

    // normalize pt slices, only if certain number of bins filled, otherwise set to zero
    Int_t       emptyBins       = 0;
    Double_t    newBinContent   = 0.;
    Double_t    newBinError     = 0.;
    Double_t    norm            = 0.;
    Double_t    integral        = 0.;
    Double_t    integralErr     = 0.;
    Double_t    deltaY          = 0.;
    for (Int_t xBin=1; xBin<=hist->GetNbinsX(); xBin++) {

        // check for empty bins in pt slices
        emptyBins               = 0;
        for (Int_t yBin=1; yBin<=hist->GetNbinsY(); yBin++) {
            if (!(hist->GetBinContent(xBin, yBin))) emptyBins++;
        }

        // rapidity range
        deltaY                  = hist->GetYaxis()->GetXmax() - hist->GetYaxis()->GetXmin();

        // calculate norm
        integral                = hist->IntegralAndError(xBin, xBin, 1, hist->GetNbinsY(), integralErr);
        norm                    = deltaY/integral/hist->GetYaxis()->GetBinWidth(1);

        // only use pt slices w/o empty bins
        newBinContent           = 0.;
        newBinError             = 0.;
        if (emptyBins || integralErr/integral > thresh) {
            for (Int_t yBin=1; yBin<=hist->GetNbinsY(); yBin++) {
                hist->SetBinContent(xBin, yBin, 0);
                hist->SetBinError(  xBin, yBin, 0);
            }
        } else {
            for (Int_t yBin=1; yBin<=hist->GetNbinsY(); yBin++) {
                newBinContent   = norm * hist->GetBinContent(xBin, yBin);
                newBinError     = norm * hist->GetBinError(xBin, yBin);

                hist->SetBinContent(xBin, yBin, newBinContent);
                hist->SetBinError(  xBin, yBin, newBinError);
            }
        }
    }
}
//================================================================================================================
//Function that reads txt file for v2 input
//================================================================================================================
TGraphAsymmErrors* ReadFileforv2Input(TString particle, TString centrality, Bool_t ReturnStatGraph){

  //opening the file
  ifstream inputfile;
  inputfile.open(Form("PbPb/2760_GeV/v2/%s%s.txt",particle.Data(),centrality.Data()));
  if (!inputfile.is_open()){
   cout << "NO CORRECT v2 INPUT FILE... RETURNING!.." << endl;
   return 0;
  }

  string line;
  Int_t nlines = 0; Int_t nstart = 0;
  if (inputfile.is_open()){
    //cout << "txt file (ReadFileforv2Input) is open..." << endl;
    while ( getline (inputfile,line) ){
      //cout << line << '\n';
      nlines++;
      if(line.compare(0,5,"xdesc")==0) nstart = nlines;
    }
  }

  inputfile.clear();
  inputfile.seekg(0, inputfile.beg);
  for(Int_t i=0; i<nstart; i++){
    getline (inputfile,line);
  }

  const Int_t nPoints = nlines-nstart-1;
  Float_t x[nPoints],xlow[nPoints],xhigh[nPoints],y[nPoints],yStatHigh[nPoints],yStatLow[nPoints],ySystHigh[nPoints],ySystLow[nPoints];
  for(Int_t i=0; i<(nlines-nstart-1); i++){
      inputfile >> x[i] >> xlow[i] >> xhigh[i] >> y[i] >> yStatLow[i] >> yStatHigh[i] >> ySystLow[i] >> ySystHigh[i] ;
      yStatHigh[i] = TMath::Abs(yStatHigh[i]); yStatLow[i] = TMath::Abs(yStatLow[i]); ySystHigh[i] = TMath::Abs(ySystHigh[i]); ySystLow[i] = TMath::Abs(ySystLow[i]);
      //cout << x[i] << " " << xlow[i] << " " << xhigh[i] << " " << y[i] << " " << yStatLow[i] << " " << yStatHigh[i] << " " << ySystLow[i] << " " << ySystHigh[i] << endl;
  }
  inputfile.close();

  TString     tempName;
  if(ReturnStatGraph){
    tempName = Form("%s_%s_STAT",particle.Data(),centrality.Data());
  }else{
    tempName = Form("%s_%s_SYST",particle.Data(),centrality.Data());
  }
  TGraphAsymmErrors* returnGraphStat                      = new TGraphAsymmErrors(nPoints,x,y,0,0,yStatLow,yStatHigh);
  TGraphAsymmErrors* returnGraphSyst                      = new TGraphAsymmErrors(nPoints,x,y,0,0,ySystLow,ySystHigh);
  returnGraphStat->SetName(tempName);
  returnGraphSyst->SetName(tempName);

  returnGraphStat->SetTitle(tempName);
  returnGraphSyst->SetTitle(tempName);


  if(ReturnStatGraph){
    return returnGraphStat;
  }else{
    return returnGraphSyst;
  }

}

//================================================================================================================
//Function to parse TGraphAsymmErrors from HEP data file
// - lines that should not be used from the file must be commented out with "%"
// - totalNumberOfColumns should be set to the total number of columns in the file, which are not commented out
// - column numbering starts at 0
// - rows should be given in order of increasing pt
// - if x error column numbers are set to -1, half the distance to the neighbouring bin is taken as value (syemmtric around x)
// - this function cannot deal with only one of columnXErrLow and columnXErrHigh being set to -1
// - this function can not deal with missing y errors (i.e. columnYErrLow or columnYErrHigh set to -1)
// - for symmetric errors and errors stored in column, e.g. columnXErrLow and columnXErrHigh can be set to the same number
// - if err. columns store error values, set isXErrVal or isYErrVal to true
// - if err. columns store bin boundaires, set isXErrVal or isYErrVal to false
//================================================================================================================
TGraphAsymmErrors* ParseHEPData(TString hepDataFile,
                                Int_t   totalNumberOfColumns,
                                Int_t   columnX,
                                Int_t   columnXErrLow,
                                Int_t   columnXErrHigh,
                                Int_t   columnY,
                                Int_t   columnYErrLow,
                                Int_t   columnYErrHigh,
                                Bool_t  isXErrVal,
                                Bool_t  isYErrVal,
                                Bool_t  debugMode = kFALSE) {

    // create streamer
    ifstream file;
    if (debugMode) cout << "HEP data file: " << hepDataFile.Data() << endl;
    file.open(hepDataFile,ios_base::in);
    if (!file) {
        cout << "ERROR: HEP data file " << hepDataFile.Data() << " not found!" << endl;
        return NULL;
    }

    // check for correct column numbers
    if (columnX<0) {
        cout << "ERROR: columnX set to " << columnX << endl;
        return NULL;
    }
    if (columnY<0) {
        cout << "ERROR: columnY set to " << columnY << endl;
        return NULL;
    }
    if (columnYErrLow<0 || columnYErrHigh<0) {
        cout << "ERROR: columnYErrLow set to " << columnYErrLow << " and columnYErrHigh set to " << columnYErrHigh << endl;
        return NULL;
    }

    // initialize vectors for temporary storage of values
    std::vector<Double_t> xVal;
    std::vector<Double_t> xErrLow;
    std::vector<Double_t> xErrHigh;
    std::vector<Double_t> yVal;
    std::vector<Double_t> yErrLow;
    std::vector<Double_t> yErrHigh;

    // read from file
    TString                 tempString;
    std::vector<TString>    tempStringColumn(totalNumberOfColumns);
    std::string line;
    for( std::string line; getline(file, line); ) {
        file >> tempString;
        if (!tempString.BeginsWith("%") && !tempString.BeginsWith("%") && tempString.CompareTo("")) {
            tempStringColumn[0]     = tempString;
            if (debugMode) cout << tempStringColumn[0].Data() << "\t";
            for (Int_t i=1; i<totalNumberOfColumns; i++) {
                file >> tempStringColumn[i];
                if (debugMode) cout << tempStringColumn[i].Data() << "\t";
            }
            if (debugMode) cout << endl;

            // x value and error
            xVal.push_back(tempStringColumn[columnX].Atof());
            if (columnXErrLow>=0)   xErrLow.push_back(tempStringColumn[columnXErrLow].Atof());
            else                    xErrLow.push_back(-1);
            if (columnXErrHigh>=0)  xErrHigh.push_back(tempStringColumn[columnXErrHigh].Atof());
            else                    xErrHigh.push_back(-1);

            // y value and error
            yVal.push_back(tempStringColumn[columnY].Atof());
            yErrLow.push_back(tempStringColumn[columnYErrLow].Atof());
            yErrHigh.push_back(tempStringColumn[columnYErrHigh].Atof());
        } else
            continue;
    }

    // check for equal number of rows for each column
    Bool_t  isEqualNumberOfRows     = kTRUE;
    Int_t   nRowsTemp[6];
    nRowsTemp[0]                    = xVal.size();
    nRowsTemp[1]                    = xErrLow.size();
    nRowsTemp[2]                    = xErrHigh.size();
    nRowsTemp[3]                    = yVal.size();
    nRowsTemp[4]                    = yErrLow.size();
    nRowsTemp[5]                    = yErrHigh.size();
    for (Int_t i=0; i<5; i++) {
        if (nRowsTemp[i]!=nRowsTemp[i+1]) {
            isEqualNumberOfRows     = kFALSE;
            break;
        }
    }
    if (!isEqualNumberOfRows) {
        cout << "number of rows in " << hepDataFile.Data() << " are not equal for different columns!" << endl;
        return NULL;
    }
    Int_t nRows                     = xVal.size();

    // calculate x errors if necessary (i.e. column numbers set to -1)
    std::vector<Double_t> tempXErr(xVal.size());
    if (columnXErrLow<0 || columnXErrHigh<0) {
        for (Int_t i=0; i<nRows; i++) {

            // calculate x error
            if (i==0)               tempXErr[i] = (xVal[1]-xVal[0])/2;
            else if (i==nRows-1)    tempXErr[i] = xVal[i]-(xVal[i-1] + tempXErr[i-1]);
            else                    tempXErr[i] = (xVal[i]-xVal[i-1])/2;

            // set error
            xErrLow[i]              = tempXErr[i];
            xErrHigh[i]             = tempXErr[i];
        }
    }

    // calculate errors if bin boundaries were given
    if (!isXErrVal && columnXErrLow>=0 && columnXErrHigh>=0) {
        for (Int_t i=0; i<nRows; i++) {
            xErrLow[i]              = TMath::Abs(xVal[i]-xErrLow[i]);
            xErrHigh[i]             = TMath::Abs(xErrHigh[i]-xVal[i]);
        }
    }
    if (!isYErrVal) {
        for (Int_t i=0; i<nRows; i++) {
            yErrLow[i]              = TMath::Abs(yVal[i]-yErrLow[i]);
            yErrHigh[i]             = TMath::Abs(yErrHigh[i]-yVal[i]);
        }
    }

    // set errors to absolute values, direction is taken care of by TGraphAsymmErrors
    for (Int_t i=0; i<nRows; i++) {
        xErrLow[i]                  = TMath::Abs(xErrLow[i]);
        xErrHigh[i]                 = TMath::Abs(xErrHigh[i]);

        yErrLow[i]                  = TMath::Abs(yErrLow[i]);
        yErrHigh[i]                 = TMath::Abs(yErrHigh[i]);
    }

    // cout values (debug mode)
    if (debugMode) {
        cout << "nRows = " << nRows << endl;
        for (Int_t i=0; i<nRows; i++) {
            cout << "x = " << xVal[i] << "\t+ " << xErrHigh[i] << "\t- " << xErrLow[i] << "\t y = " << yVal[i] << "\t+ " << yErrHigh[i] << "\t- " << yErrLow[i] << endl;
        }
    }

    // create TGraphAsymmErrors
    TGraphAsymmErrors* graph        = new TGraphAsymmErrors(nRows);
    for (Int_t i=0; i<nRows; i++) {
        graph->SetPoint(        i, xVal[i], yVal[i]);
        graph->SetPointError(   i, xErrLow[i], xErrHigh[i], yErrLow[i], yErrHigh[i]);
    }
    return graph;
}


//=================================================================================================================
// Create TGraphAsymmErrors without xerrors
//=================================================================================================================
void ProduceGraphAsymmWithoutXErrors(TGraphAsymmErrors* inputgraph){
    Int_t n                 = inputgraph->GetN();
    Double_t* xValue        = inputgraph->GetX();
    Double_t* xErrorHigh    = inputgraph->GetEXhigh();
    Double_t* xErrorLow     = inputgraph->GetEXlow();
    Double_t* yValue        = inputgraph->GetY();
    Double_t* yErrorLow     = inputgraph->GetEYlow();
    Double_t* yErrorHigh    = inputgraph->GetEYhigh();
    for (Int_t i= 0; i < n; i++){
        xErrorHigh[i]       = 0.;
        xErrorLow[i]        = 0.;
    }
    inputgraph              = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    //     inputgraph->Print();
}

//=================================================================================================================
// Create TGraphErrors without xerrors
//=================================================================================================================
void ProduceGraphWithoutXErrors(TGraphErrors* inputgraph){
    Int_t n                 = inputgraph->GetN();
    Double_t* xValue        = inputgraph->GetX();
    Double_t* xError        = inputgraph->GetEX();
    Double_t* yValue        = inputgraph->GetY();
    Double_t* yError        = inputgraph->GetEY();
    for (Int_t i= 0; i < n; i++){
        xError[i]           = 0.;
    }
    inputgraph              = new TGraphErrors(n,xValue,yValue,xError,yError);
    //     inputgraph->Print();
}

//=================================================================================================================
// mT scaling with TH1
//=================================================================================================================
// MASS   0=>PIZERO, 1=>ETA, 2=>RHO0, 3=>OMEGA, 4=>ETAPRIME, 5=>PHI, 6=>JPSI, 7=>SIGMA, 8=>K0s, 9=>DELTA++, 10=>DELTA+, 11=>DELTA-, 12=>DELTA0, 13=>Rho+, 14=>Rho-, 15=>K0*, 16=>K0l, 17=>Lambda, 18=>K+, 19=>K-, 20=>Omega+, 21=>Omega-, 22=>Xi+, 23=>Xi-, 24=>Sigma+, 25=>Sigma-
//const Double_t AliGenEMlibV2::fgkHM[26] = {0.1349766, 0.547853, 0.77549, 0.78265, 0.95778, 1.019455, 3.096916, 1.192642, 0.497614, 1.2311, 1.2349, 1.2349, 1.23340, 0.77549, 0.77549, 0.896, 0.497614, 1.115683, 0.493677, 0.493677, 1.67245, 1.67245, 1.32171, 1.32171, 1.3828, 1.3872};
//const Double_t AliGenEMlibV2::fgkMtFactor[3][26] = {
//  {1., 0.476, 1.0, 0.85, 0.4, 0.13, 1., 0.49, 0.575, 1, 1, 1, 1, 1.0, 1.0, 1.0, 0.575, 0.18, 0.41, 0.41, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
//};
TH1F* mTscaleHistogramFromTF1(TF1* input, TString inputName, Double_t massInput, Double_t massOutput, Double_t ratio, Int_t changeQ){
  Double_t ptBinning[65]              = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                         1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                         2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                                         3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
                                         5.0, 5.4, 5.8, 6.2, 6.6, 7.0, 7.5, 8.0, 8.5, 9.0,
                                         10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0,
                                         30.0, 35.0, 40.0, 45.0, 50.0};
  Int_t nBinsX                        = 64;

  TH1F *mtScaled = new TH1F(Form("mtScaled_%s",inputName.Data()),Form("mtScaled_%s",inputName.Data()),nBinsX,ptBinning);
  for (Int_t i=1; i<=mtScaled->GetNbinsX(); i++) {
    Double_t ptInput = mtScaled->GetBinCenter(i);
    Double_t mtOut = TMath::Sqrt(massOutput*massOutput + ptInput*ptInput);
    Double_t R = input->Eval(TMath::Sqrt(mtOut*mtOut - massInput*massInput)) * ratio;
    if(changeQ == 1) R /= (2*TMath::Pi()*ptInput); //bin shift error not considered (yet)
    mtScaled->SetBinContent(i,R);
    mtScaled->SetBinError(i,R*0.05);
  }

//  mtScaled->Rebin(nBinsX,Form("mtScaled_%s_rebin",inputName.Data()),ptBinning);
//  mtScaled = (TH1F*)gDirectory->Get(Form("mtScaled_%s_rebin",inputName.Data()));
  return mtScaled;
}


//=================================================================================================================
// Extract Binning of spectrum as array
// -> return number of bins
// -> Reference to binning array handed
//=================================================================================================================
Int_t GetBinning(TObject *Obj_Dummy, Double_t* doubleBinningX){
    TString ClassName                           = Obj_Dummy->ClassName();
    if(ClassName.BeginsWith("TH1")){
        TH1D *histo                             = (TH1D*)Obj_Dummy;
        Int_t bin                               = 0;
        for (Int_t i = 1; i < histo->GetNbinsX()+1; i++){
            if (histo->GetBinContent(i) != 0){
                doubleBinningX[bin]             = histo->GetBinLowEdge(i);
                doubleBinningX[bin+1]           = histo->GetXaxis()->GetBinUpEdge(i);
                bin++;
            }
        }
        return bin+1;
    } else if(ClassName.CompareTo("TGraphErrors")==0){
        TGraphErrors *graph                     = (TGraphErrors*)Obj_Dummy;
        Double_t* binCenter                     = graph->GetX();
        Double_t* binContent                    = graph->GetY();
        Double_t* binError                      = graph->GetEX();
        Int_t nBins                             = graph->GetN();
        Int_t bin                               = 0;
        for (Int_t i = 0; i < nBins; i++){
            if (binContent[i] != 0){
                doubleBinningX[bin]             = binCenter[i]-binError[i];
                doubleBinningX[bin+1]           = binCenter[i]+binError[i];
                bin++;
            }
        }
        return bin+1;
    } else if(ClassName.CompareTo("TGraphAsymmErrors")==0){
        TGraphAsymmErrors *graph                = (TGraphAsymmErrors*)Obj_Dummy;
        Double_t* binCenter                     = graph->GetX();
        Double_t* binContent                    = graph->GetY();
        Double_t* binErrorDown                  = graph->GetEXlow();
        Double_t* binErrorUp                    = graph->GetEXhigh();
        Int_t nBins                             = graph->GetN();
        Int_t bin                               = 0;
        for (Int_t i = 0; i < nBins; i++){
            if (binContent[i] != 0){
                doubleBinningX[bin]             = binCenter[i]-binErrorDown[i];
                doubleBinningX[bin+1]           = binCenter[i]+binErrorUp[i];
                bin++;
            }
        }
        return bin+1;
    } else {
        cout << " class not defined" << endl;
        return 0;
    }

}


//=================================================================================================================
// Compare two different binning and return a common binning if possible
// -> return number of bins of common binning
// -> References to comb binning array, rebinnged version of older handed
//=================================================================================================================
Int_t CompareBinning(   Int_t nBinsA,
                        Double_t* binningA,
                        Int_t nBinsB,
                        Double_t* binningB,
                        Double_t* newBinning,
                        Int_t* nBinsToBeCombinedA,
                        Int_t* nBinsToBeCombinedB,
                        TString returnStr           = "",
                        Double_t decisionBoundary   = 0.0000001
                    ){
    Int_t startingBin = 0;
    //Finding startBin
    cout << "Binning A" << endl;
        for (Int_t i = 0; i < nBinsA ; i++){
            cout << binningA[i] << "\t" ;
        }
    cout << endl;
    cout << "Binning B" << endl;
        for (Int_t i = 0; i < nBinsB ; i++){
            cout << binningB[i] << "\t" ;
        }
    cout << endl;
    if (binningA[0] < binningB[0]){
        while ((binningB[0] - binningA[startingBin])>decisionBoundary && startingBin <nBinsA){
            startingBin++;
        }
        cout << "Binning A" << endl;
        for (Int_t i = startingBin; i < nBinsA ; i++){
            cout << binningA[i] << "\t" ;
        }
        cout << endl;
        cout << "Binning B" << endl;
        for (Int_t i = 0; i < nBinsB ; i++){
            cout << binningB[i] << "\t" ;
        }
        cout << endl;

        Int_t c                             = 1;
        Int_t startBinNewBins               = 1;
        Int_t binsToBeMergedB               = 1;
        Int_t binsToBeMergedA               = 1;
        cout << "Binning A starts earlier, combined binning will start at " << startingBin << " with " << binningA[startingBin] << endl;
        newBinning[0] = binningA[startingBin];
        for (Int_t i = startingBin+1; i < nBinsA; i++){
            if (c > nBinsB-1) return startBinNewBins;
            newBinning[startBinNewBins]     = binningA[i];
            binsToBeMergedB                 = 1;
            binsToBeMergedA                 = 1;
            cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;
            while(TMath::Abs(binningA[i]- binningB[c]) > decisionBoundary ){
                if (c > nBinsB-1) return startBinNewBins;
                if (i > nBinsA-1) return startBinNewBins;

                if((binningA[i] - binningB[c])>decisionBoundary){
                cout << "nächstes bin ist größer in B" << endl;
                cout << " müssen einen hoch" << endl;
                c++;
                binsToBeMergedB++;
                }else{
                cout << "nächstes bin ist größer in A" << endl;
                cout << " müssen einen hoch" << endl;
                i++;
                newBinning[startBinNewBins]     = binningA[i];
                binsToBeMergedA++;
                }
                cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;

            }
            nBinsToBeCombinedB[startBinNewBins-1]       = binsToBeMergedB;
            nBinsToBeCombinedA[startBinNewBins-1]       = binsToBeMergedA;
            if (c > nBinsB-1) return startBinNewBins;
            cout << "exiting loop " << binningA[i] << "\t" << binningB[c] << "\t bins needed to merged A :"<<nBinsToBeCombinedA[startBinNewBins-1] <<" B :"<<nBinsToBeCombinedB[startBinNewBins-1] <<endl;
            startBinNewBins++;
            c++;
        }
        return startBinNewBins;
    } else  if (binningB[0] < binningA[0]){
        while (!(TMath::Abs(binningA[0]- binningB[startingBin]) < decisionBoundary) && startingBin< nBinsB){
            cout << "deviation to first bin in A for startBin: "<< startingBin << "\t" <<TMath::Abs(binningA[0]- binningB[startingBin]) << endl;
            startingBin++;
        }
        Int_t check2 = 0;
        while (startingBin == nBinsB){
            cout << "Failed to evalute starting point in attempt " << check2    << endl;
            check2++;
            startingBin=0;
            while (!(TMath::Abs(binningA[check2]- binningB[startingBin]) < decisionBoundary) && startingBin< nBinsB){
                cout << TMath::Abs(binningA[check2]- binningB[startingBin]) << endl;
                startingBin++;
            }
        }
        cout << "Binning A" << endl;
        for (Int_t i = 0; i < nBinsA ; i++){
            cout << binningA[i] << "\t" ;
        }
        cout << endl;
        cout << "Binning B" << endl;
        for (Int_t i = startingBin; i < nBinsB ; i++){
            cout << binningB[i] << "\t" ;
        }
        cout << endl;
        cout << "Binning B starts earlier, combined binning will start at " << startingBin << " with " << binningB[startingBin] << endl;

        Int_t c                             = startingBin+1;
        Int_t startBinNewBins               = 1;
        Int_t binsToBeMergedB               = 1;
        Int_t binsToBeMergedA               = 1;
        newBinning[0]                       = binningA[0];
        for (Int_t i = 1; i < nBinsA; i++){
            if (c > nBinsB-1) return startBinNewBins;
            newBinning[startBinNewBins]     = binningA[i];
            binsToBeMergedB                 = 1;
            binsToBeMergedA                 = 1;
                cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;
                while(TMath::Abs(binningA[i]- binningB[c]) > decisionBoundary ){
                if (c > nBinsB-1) return startBinNewBins;
                if (i > nBinsA-1) return startBinNewBins;

                if((binningA[i] - binningB[c])>decisionBoundary){
                    cout << "nächstes bin ist größer in B" << endl;
                    cout << " müssen einen hoch" << endl;
                    c++;
                    binsToBeMergedB++;
                }else{
                    cout << "nächstes bin ist größer in A" << endl;
                    cout << " müssen einen hoch" << endl;
                    i++;
                    newBinning[startBinNewBins]     = binningA[i];
                    binsToBeMergedA++;
                }
                cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;

                }
            nBinsToBeCombinedB[startBinNewBins-1]   = binsToBeMergedB;
            nBinsToBeCombinedA[startBinNewBins-1]   = binsToBeMergedA;
            cout << "exiting loop " << binningA[i] << "\t" << binningB[c] << "\t bins needed to merged A :"<<nBinsToBeCombinedA[startBinNewBins-1] <<" B :"<<nBinsToBeCombinedB[startBinNewBins-1] <<endl;
            startBinNewBins++;
            c++;
        }
        return startBinNewBins;

    } else {
        cout << "Binning A" << endl;
        for (Int_t i = 0; i < nBinsA ; i++){
            cout << binningA[i] << "\t" ;
        }
        cout << endl;
        cout << "Binning B" << endl;
        for (Int_t i = startingBin; i < nBinsB ; i++){
            cout << binningB[i] << "\t" ;
        }
        cout << endl;

        cout << "Both start at the same value " << binningA[0] << endl;

        Int_t c                             = 1;
        Int_t startBinNewBins               = 1;
        Int_t binsToBeMergedB               = 1;
        Int_t binsToBeMergedA               = 1;
        newBinning[0]                       = binningA[0];
        for (Int_t i = 1; i < nBinsA; i++){
            if (c > nBinsB-1) return startBinNewBins-1;
            newBinning[startBinNewBins]     = binningA[i];
            binsToBeMergedB                 = 1;
            binsToBeMergedA                 = 1;
                cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;
                while(TMath::Abs(binningA[i]- binningB[c]) > decisionBoundary ){
                if (c > nBinsB-1) return startBinNewBins;
                if (i > nBinsA-1) return startBinNewBins;

                if((binningA[i] - binningB[c])>decisionBoundary){
//                     cout << "nächstes bin ist größer in B" << endl;
//                     cout << " müssen einen hoch" << endl;
                    c++;
                    binsToBeMergedB++;
                }else{
//                     cout << "nächstes bin ist größer in A" << endl;
//                     cout << " müssen einen hoch" << endl;
                    i++;
                    newBinning[startBinNewBins]     = binningA[i];
                    binsToBeMergedA++;
                }
                cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;

                }
            nBinsToBeCombinedB[startBinNewBins-1]   = binsToBeMergedB;
            nBinsToBeCombinedA[startBinNewBins-1]   = binsToBeMergedA;
                cout << "exiting loop " << binningA[i] << "\t" << binningB[c] << "\t bins needed to merged A :"<<nBinsToBeCombinedA[startBinNewBins-1] <<" B :"<<nBinsToBeCombinedB[startBinNewBins-1] <<endl;
            startBinNewBins++;
            c++;
        }
        return startBinNewBins;
    }
    return 0;
    if(returnStr){}
}

//=================================================================================================================
// Extract first common bin between 2 binnings
// -> return first common bin
//=================================================================================================================
Int_t FindFirstCommonBin(   Double_t* vectorNewBinning,
                            Double_t* oldBinning,
                            Int_t commonBin             = 0,
                            Double_t decisionBoundary   = 0.0000001
                        ){
    Int_t startingBin   = 0;
    //Finding startBin

    while (startingBin<commonBin && !(TMath::Abs(vectorNewBinning[0] - oldBinning[startingBin])<decisionBoundary)){
        startingBin++;
        cout << startingBin << "\t" <<!(TMath::Abs(vectorNewBinning[0] - oldBinning[startingBin])<decisionBoundary) << endl;
    }
    if(startingBin>=commonBin){
        startingBin = 0;
        while (startingBin<commonBin && !(vectorNewBinning[0] - oldBinning[startingBin]<decisionBoundary)){
            startingBin++;
            cout << startingBin << "\t" <<!(vectorNewBinning[0] - oldBinning[startingBin]<decisionBoundary) << endl;
        }
        startingBin--;
    }
    return startingBin;
}

//=================================================================================================================
// Rebin orginal object to match a given binning taking into account stat and sys errors
// -> references to new objects given
//=================================================================================================================
void RebinObjects(  TObject* Obj_DummyStat,
                    TObject* Obj_DummySyst,
                    Double_t* vectorNewBinning,
                    Int_t* vectorRebinFactors,
                    Int_t nCommonBins,
                    AliConvDataObject* outputObject,
                    TString ClassNameStat,
                    TString ClassNameSyst,
                    Bool_t scaleByBinCenter             = kFALSE
                ) { //

    TGraphAsymmErrors *graph;
    TGraphAsymmErrors *graphCopy;
    TGraphErrors      *graph2;
    TGraphErrors      *graphCopy2;

    Double_t binningOldStat[200] ;
    Double_t binningOldSyst[200] ;
    Int_t validBinsOldStat          = GetBinning(Obj_DummyStat, binningOldStat);
    Int_t validBinsOldSys           = GetBinning(Obj_DummySyst, binningOldSyst);

    Int_t firstBinStat              = FindFirstCommonBin(vectorNewBinning, binningOldStat,nCommonBins);
    Int_t firstBinSyst              = FindFirstCommonBin(vectorNewBinning, binningOldSyst,nCommonBins);

    cout << "FirstBin stat "    << firstBinStat     << "\t" << binningOldStat[firstBinStat]     << endl;
    cout << "FirstBin sys "     << firstBinSyst     << "\t" << binningOldSyst[firstBinSyst]     << endl;
    cout << "statistical Errors" << endl;
    if(ClassNameStat.BeginsWith("TH1")){
        TH1D *histo                 = (TH1D*)Obj_DummyStat;
        Int_t indBin                = 0;
        cout << "using histo as input " << endl;
        while (vectorNewBinning[0] > histo->GetBinLowEdge(indBin))    { indBin++;}
        for (Int_t commonBin                        = 0; commonBin < nCommonBins-1; commonBin++){
            outputObject[commonBin].valueY          = 0;
            outputObject[commonBin].errorYStatLow   = 0;
            outputObject[commonBin].valueX          = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2;
            outputObject[commonBin].errorXLow       = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
            outputObject[commonBin].errorXHigh      = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
            for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
                if (scaleByBinCenter){
                    cout << indBin << "\t" << histo->GetBinCenter(indBin) << "\t"<< (histo->GetBinContent(indBin)*histo->GetBinCenter(indBin)*histo->GetBinWidth(indBin)*2*TMath::Pi()) << endl;
                    outputObject[commonBin].valueY          = outputObject[commonBin].valueY+(histo->GetBinContent(indBin)*histo->GetBinCenter(indBin)*histo->GetBinWidth(indBin)*2*TMath::Pi());
                    outputObject[commonBin].errorYStatLow   = outputObject[commonBin].errorYStatLow+TMath::Power((histo->GetBinError(indBin)*histo->GetBinCenter(indBin)*histo->GetBinWidth(indBin)*2*TMath::Pi()),2);
                }    else {
                    cout << indBin << "\t" << histo->GetBinCenter(indBin) << "\t"<<(histo->GetBinContent(indBin)*histo->GetBinWidth(indBin)) << endl;
                    outputObject[commonBin].valueY          = outputObject[commonBin].valueY+(histo->GetBinContent(indBin)*histo->GetBinWidth(indBin));
                    outputObject[commonBin].errorYStatLow   = outputObject[commonBin].errorYStatLow+TMath::Power((histo->GetBinError(indBin)*histo->GetBinWidth(indBin)),2);
                }
                indBin++;
            }
            outputObject[commonBin].errorYStatLow   = TMath::Sqrt(outputObject[commonBin].errorYStatLow);
            outputObject[commonBin].errorYStatHigh  = outputObject[commonBin].errorYStatLow;
            cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYStatHigh<<"\t-"<< outputObject[commonBin].errorYStatLow<<endl;
        }
    } else if(ClassNameStat.CompareTo("TGraphErrors")==0){
        cout << "using TGraphErrors as input " << endl;
        graph2         = (TGraphErrors*)Obj_DummyStat;
        graphCopy2     = (TGraphErrors*)graph2->Clone("GraphCopy");
        Double_t* valueX            = graphCopy2->GetX();
        Double_t* valueY            = graphCopy2->GetY();
        Double_t* errorX            = graphCopy2->GetEX();
        Double_t* errorY            = graphCopy2->GetEY();
        for (Int_t i = 0; i < graphCopy2->GetN();i++){
            cout << valueX[i] << "\t" << valueY[i] << "\t+-" << errorY[i] << endl;
            if (scaleByBinCenter){
                valueY[i]           = valueY[i]*valueX[i]*(errorX[i]+errorX[i])*2*TMath::Pi();
                errorY[i]           = errorY[i]*valueX[i]*(errorX[i]+errorX[i])*2*TMath::Pi();
            }    else {
                valueY[i]           = valueY[i]*(errorX[i]+errorX[i]);
                errorY[i]           = errorY[i]*(errorX[i]+errorX[i]);
            }
        }
        Int_t indBin                                = firstBinStat;
        for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
            outputObject[commonBin].valueY          = 0;
            outputObject[commonBin].errorYStatLow   = 0;
            outputObject[commonBin].valueX          = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2;
            outputObject[commonBin].errorXLow       = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
            outputObject[commonBin].errorXHigh      = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
            for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
                cout << indBin << "\t" << valueY[indBin] << endl;
                outputObject[commonBin].valueY          = outputObject[commonBin].valueY+valueY[indBin];
                outputObject[commonBin].errorYStatLow   = outputObject[commonBin].errorYStatLow+TMath::Power(errorY[indBin],2);
                indBin++;
            }
            outputObject[commonBin].errorYStatLow   = TMath::Sqrt(outputObject[commonBin].errorYStatLow);
            outputObject[commonBin].errorYStatHigh  = outputObject[commonBin].errorYStatLow;
            cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYStatHigh<<"\t-"<< outputObject[commonBin].errorYStatLow<<endl;
        }
    } else if(ClassNameStat.CompareTo("TGraphAsymmErrors")==0){
        cout << "using TGraphAsymmErrors as input " << endl;
        graph        = (TGraphAsymmErrors*)Obj_DummyStat;
        graphCopy    = (TGraphAsymmErrors*)graph->Clone("GraphCopy");
        Double_t* valueX                = graphCopy->GetX();
        Double_t* valueY                = graphCopy->GetY();
        Double_t* errorXlow             = graphCopy->GetEXlow();
        Double_t* errorXhigh            = graphCopy->GetEXhigh();
        Double_t* errorYlow             = graphCopy->GetEYlow();
        Double_t* errorYhigh            = graphCopy->GetEYhigh();
        for (Int_t i = 0; i < graphCopy->GetN();i++){
            if (scaleByBinCenter){
                valueY[i]               = valueY[i]*valueX[i]*(errorXlow[i]+errorXhigh[i])*2*TMath::Pi();
                errorYlow[i]            = errorYlow[i]*valueX[i]*(errorXlow[i]+errorXhigh[i])*2*TMath::Pi();
                errorYhigh[i]           = errorYhigh[i]*valueX[i]*(errorXlow[i]+errorXhigh[i])*2*TMath::Pi();
            }    else {
                valueY[i]               = valueY[i]*(errorXlow[i]+errorXhigh[i]);
                errorYlow[i]            = errorYlow[i]*(errorXlow[i]+errorXhigh[i]);
                errorYhigh[i]           = errorYhigh[i]*(errorXlow[i]+errorXhigh[i]);
            }
        }
        cout<< "after modification" << endl;
//         graphCopy->Print();
        Int_t indBin                                = firstBinStat;
        for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
            outputObject[commonBin].valueY          = 0;
            outputObject[commonBin].errorYStatLow   = 0;
            outputObject[commonBin].errorYStatHigh  = 0;
            outputObject[commonBin].valueX          = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2;
            outputObject[commonBin].errorXLow       = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
            outputObject[commonBin].errorXHigh      = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];

            for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
                cout << indBin << "\t" << valueY[indBin] << endl;
                outputObject[commonBin].valueY          = outputObject[commonBin].valueY+valueY[indBin];
                outputObject[commonBin].errorYStatHigh  = outputObject[commonBin].errorYStatHigh + TMath::Power(errorYhigh[indBin],2);
                outputObject[commonBin].errorYStatLow   = outputObject[commonBin].errorYStatLow + TMath::Power(errorYlow[indBin],2);
                indBin++;
            }
            outputObject[commonBin].errorYStatLow   = TMath::Sqrt(outputObject[commonBin].errorYStatLow);
            outputObject[commonBin].errorYStatHigh  = TMath::Sqrt(outputObject[commonBin].errorYStatHigh);
            cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYStatHigh<<"\t-"<< outputObject[commonBin].errorYStatLow<<endl;
        }
    } else {
        cout << " class for Stat not defined" << endl;
        return;
    }


    cout << "systematic errors" << endl;
    if(ClassNameSyst.BeginsWith("TH1")){
        cout << "\n \n TH1" << endl << endl;
        TH1D *histo                                 = (TH1D*)Obj_DummySyst;
        Int_t indBin                                = 0;
        while (vectorNewBinning[0] > histo->GetBinLowEdge(indBin))    { indBin++;}
        for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
            outputObject[commonBin].errorYSystLow   = 0;
            for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
//                 if (scaleByBinCenter){
//                     cout << indBin << "\t" << histo->GetBinCenter(indBin) << "\t"<< histo->GetBinWidth(indBin) <<"\t" <<histo->GetBinError(indBin)/histo->GetBinContent(indBin)*100 <<"\t"<< histo->GetBinError(indBin)/histo->GetBinContent(indBin)*histo->GetBinWidth(indBin)<<endl; //
                outputObject[commonBin].errorYSystLow   = outputObject[commonBin].errorYSystLow + histo->GetBinError(indBin)/histo->GetBinContent(indBin)*histo->GetBinWidth(indBin);
//                     cout << indBin << "\t" << outputObject[commonBin].errorYSystLow<<endl; //
//                 }    else {
//                     cout << indBin << "\t" << histo->GetBinCenter(indBin) << "\t"<<(histo->GetBinContent(indBin)*histo->GetBinWidth(indBin)) << endl;
//                     outputObject[commonBin].errorYSystLow = outputObject[commonBin].errorYSystLow+(histo->GetBinContent(indBin)*histo->GetBinWidth(indBin))*histo->GetBinContent(indBin);

//                 }
                indBin++;
            }
            outputObject[commonBin].errorYSystLow   = outputObject[commonBin].errorYSystLow/(outputObject[commonBin].errorXLow*2)*outputObject[commonBin].valueY;
            outputObject[commonBin].errorYSystHigh  = outputObject[commonBin].errorYSystLow;
            outputObject[commonBin].errorYTotHigh   = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystHigh,2) + TMath::Power(outputObject[commonBin].errorYStatHigh,2));
            outputObject[commonBin].errorYTotLow    = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystLow,2) + TMath::Power(outputObject[commonBin].errorYStatLow,2));

            cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYSystHigh<<"\t-"<< outputObject[commonBin].errorYSystLow<< "\t+" << outputObject[commonBin].errorYTotHigh<<"\t-"<< outputObject[commonBin].errorYTotLow<<"\t syst \t"<< outputObject[commonBin].errorYSystLow/outputObject[commonBin].valueY*100 <<"%"<<endl;
        }
    } else if(ClassNameSyst.CompareTo("TGraphErrors")==0){
        cout << "\n \n TGraphErrors" << endl << endl;
        graph2                         = (TGraphErrors*)Obj_DummySyst;
        graphCopy2                     = (TGraphErrors*)graph2->Clone("GraphCopy");
        Double_t* valueX                            = graphCopy2->GetX();
        Double_t* valueY                            = graphCopy2->GetY();
        Double_t* errorX                            = graphCopy2->GetEX();
        Double_t* errorY                            = graphCopy2->GetEY();
        for (Int_t i = 0; i < graphCopy2->GetN();i++){
            if (scaleByBinCenter){
                cout << i <<"\t"<<"before: " <<valueY[i] << "\t" << errorY[i] << "\t" << errorY[i]/valueY[i]*100<<"\t" << (errorX[i]+errorX[i])<< endl;
                errorY[i]                           = errorY[i]*(errorX[i]+errorX[i])/valueY[i];
                valueY[i]                           = valueY[i]*valueX[i]*(errorX[i]+errorX[i])*2*TMath::Pi();
            }    else {
                valueY[i]                           = valueY[i]*(errorX[i]+errorX[i]);
                errorY[i]                           = errorY[i]*(errorX[i]+errorX[i])/valueY[i] ;
            }
        }
        Int_t indBin                                = firstBinSyst;
        for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
            outputObject[commonBin].errorYSystLow   = 0;
            for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
//                 cout << indBin << "\t" << valueY[indBin] << endl;
                outputObject[commonBin].errorYSystLow   = outputObject[commonBin].errorYSystLow+errorY[indBin];
                cout << indBin << "\t" << valueY[indBin] << "\t"<< outputObject[commonBin].errorYSystLow<<endl;
                indBin++;
            }
            outputObject[commonBin].errorYSystLow   = outputObject[commonBin].errorYSystLow*outputObject[commonBin].valueY/ (outputObject[commonBin].errorXLow+outputObject[commonBin].errorXHigh);
            outputObject[commonBin].errorYSystHigh  = outputObject[commonBin].errorYSystLow;
            outputObject[commonBin].errorYTotHigh   = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystHigh,2) + TMath::Power(outputObject[commonBin].errorYStatHigh,2));
            outputObject[commonBin].errorYTotLow    = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystLow,2) + TMath::Power(outputObject[commonBin].errorYStatLow,2));

            cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYSystHigh<<"\t-"<< outputObject[commonBin].errorYSystLow<< "\t+" << outputObject[commonBin].errorYTotHigh<<"\t-"<< outputObject[commonBin].errorYTotLow<<endl;
        }
    } else if(ClassNameSyst.CompareTo("TGraphAsymmErrors")==0){
        cout << "\n \n TGraphAsymmErrors" << endl << endl;
        graph                    = (TGraphAsymmErrors*)Obj_DummySyst;
        graphCopy                = (TGraphAsymmErrors*)graph->Clone("GraphCopy");
        graphCopy->Print();
        Double_t* valueX                            = graphCopy->GetX();
        Double_t* valueY                            = graphCopy->GetY();
        Double_t* errorXlow                         = graphCopy->GetEXlow();
        Double_t* errorXhigh                        = graphCopy->GetEXhigh();
        Double_t* errorYlow                         = graphCopy->GetEYlow();
        Double_t* errorYhigh                        = graphCopy->GetEYhigh();
        for (Int_t i = 0; i < graphCopy->GetN();i++){
            if (scaleByBinCenter){
                cout << i <<"\t"<<"before: " <<valueY[i] << "\t" << errorYlow[i] << "\t" << errorYlow[i]/valueY[i]*100<<"\t" << (errorXlow[i]+errorXhigh[i])<< endl;
                errorYlow[i]                        = errorYlow[i]*(errorXlow[i]+errorXhigh[i])/valueY[i];
                errorYhigh[i]                       = errorYhigh[i]*(errorXlow[i]+errorXhigh[i])/valueY[i];
                valueY[i]                           = valueY[i]*valueX[i]*(errorXlow[i]+errorXhigh[i])*2*TMath::Pi();
            }    else {
                cout << "entered here" << endl;
                errorYlow[i]                        = errorYlow[i]*(errorXlow[i]+errorXhigh[i])/valueY[i];
                errorYhigh[i]                       = errorYhigh[i]*(errorXlow[i]+errorXhigh[i])/valueY[i];
                valueY[i]                           = valueY[i]*(errorXlow[i]+errorXhigh[i]);

            }
        }
        graphCopy->Print();
        Int_t indBin                                = firstBinSyst;
        cout << firstBinSyst << endl;
        for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
            outputObject[commonBin].errorYSystLow   = 0;
            outputObject[commonBin].errorYSystHigh  = 0;
            for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){

                outputObject[commonBin].errorYSystHigh  = outputObject[commonBin].errorYSystHigh + errorYhigh[indBin];
                outputObject[commonBin].errorYSystLow   = outputObject[commonBin].errorYSystLow + errorYlow[indBin];
                cout << indBin << "\t" << valueY[indBin] << "\t"<< outputObject[commonBin].errorYSystHigh<<endl;
                indBin++;
            }
            outputObject[commonBin].errorYSystLow   = outputObject[commonBin].errorYSystLow*outputObject[commonBin].valueY/ (outputObject[commonBin].errorXLow+outputObject[commonBin].errorXHigh);
            outputObject[commonBin].errorYSystHigh  = outputObject[commonBin].errorYSystHigh*outputObject[commonBin].valueY/ (outputObject[commonBin].errorXLow+outputObject[commonBin].errorXHigh);
            outputObject[commonBin].errorYTotHigh   = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystHigh,2) + TMath::Power(outputObject[commonBin].errorYStatHigh,2));
            outputObject[commonBin].errorYTotLow    = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystLow,2) + TMath::Power(outputObject[commonBin].errorYStatLow,2));

            cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYSystHigh<<"\t-"<< outputObject[commonBin].errorYSystLow<< "\t+" << outputObject[commonBin].errorYTotHigh<<"\t-"<< outputObject[commonBin].errorYTotLow<<"\t syst \t"<< outputObject[commonBin].errorYSystLow/outputObject[commonBin].valueY*100 <<"%"<<endl;
        }
    } else {
        cout << " class for Syst not defined" << endl;
        return;
    }
    return;
if (validBinsOldStat || validBinsOldSys){}
}

//=================================================================================================================
// Function which compares 2 Spectra with each other and returns the ratio
// statistical and systematic errors have to be handed over
//   - Graphs should be handed over as Clones of the original object, otherwise
//     they will be modified
//=================================================================================================================
TGraphErrors* CalculateRatioBetweenSpectraWithDifferentBinning( TObject* Obj_DummyAStat,
                                                                TObject* Obj_DummyASyst,
                                                                TObject* Obj_DummyBStat,
                                                                TObject* Obj_DummyBSyst,
                                                                Bool_t scaleByBinCenterA            = kTRUE,
                                                                Bool_t scaleByBinCenterB            = kTRUE,
                                                                TGraphErrors** graphRebinnedAStat   = NULL,
                                                                TGraphErrors** graphRebinnedASyst   = NULL,
                                                                TGraphErrors** graphRebinnedBStat   = NULL,
                                                                TGraphErrors** graphRebinnedBSyst   = NULL,
                                                                Bool_t calcSepErrGraphs             = kFALSE,
                                                                TGraphErrors** graphRatioStat       = NULL,
                                                                TGraphErrors** graphRatioSys        = NULL
                                                              ){

    cout << "Reading from Object A " << endl;
    TString ClassNameA      = Obj_DummyAStat->ClassName();
    Int_t nBinsXAA           = 0;
    if(ClassNameA.BeginsWith("TH1")){
        TH1D *histo         = (TH1D*)Obj_DummyAStat;
        nBinsXAA             = histo->GetNbinsX()+1;
    } else if(ClassNameA.BeginsWith("TGraph")){
        TGraphErrors *graph = (TGraphErrors*)Obj_DummyAStat;
        nBinsXAA             = graph->GetN()+1;
    }
    const Int_t nBinsXA = nBinsXAA;
    Double_t binningXA[nBinsXA];
    Int_t validBinsA        = GetBinning(Obj_DummyAStat, binningXA);

    cout << "Reading from Object B " << endl;
    TString ClassNameB      = Obj_DummyBStat->ClassName();
    Int_t nBinsXBB           = 0 ;
    if(ClassNameB.BeginsWith("TH1")){
        TH1D *histo         = (TH1D*)Obj_DummyBStat;
        nBinsXBB             = histo->GetNbinsX()+1;
    } else if(ClassNameB.BeginsWith("TGraph")){
        TGraphErrors *graph = (TGraphErrors*)Obj_DummyBStat;
        nBinsXBB             = graph->GetN()+1;
    }
    const Int_t nBinsXB = nBinsXBB;
    Double_t binningXB[nBinsXB];
    Int_t validBinsB        = GetBinning(Obj_DummyBStat, binningXB);

    Int_t nBinsCombb;
    if (nBinsXB < nBinsXA){
        nBinsCombb           = nBinsXA;
    } else {
        nBinsCombb           = nBinsXB;
    }
//     for (Int_t i = 0; i < nBinsXB; i++){
//         cout << binningXB[i] << "\t," ;
//     }
    const Int_t nBinsComb = nBinsCombb;
    Double_t binningCombined[nBinsComb];
    Int_t binningCombinedBinsToBeMergedA[nBinsComb];
    Int_t binningCombinedBinsToBeMergedB[nBinsComb];
    const Int_t nBinsNew          = CompareBinning( validBinsA, binningXA,validBinsB, binningXB, binningCombined, binningCombinedBinsToBeMergedA,binningCombinedBinsToBeMergedB, "A");

    cout << "Object A"  << endl;
    AliConvDataObject rebinnedA[nBinsComb];
    RebinObjects(Obj_DummyAStat,Obj_DummyASyst, binningCombined, binningCombinedBinsToBeMergedA, nBinsNew,rebinnedA, Obj_DummyAStat->ClassName(), Obj_DummyASyst->ClassName(),scaleByBinCenterA);

    cout << "Object B"  << endl;
    AliConvDataObject rebinnedB[nBinsComb];
    RebinObjects(Obj_DummyBStat,Obj_DummyBSyst, binningCombined, binningCombinedBinsToBeMergedB, nBinsNew,rebinnedB, Obj_DummyBStat->ClassName(), Obj_DummyBSyst->ClassName(), scaleByBinCenterB);

    Double_t ratioX[nBinsComb];
    Double_t errorX[nBinsComb];
    Double_t ratioY[nBinsComb];
    Double_t errorY[nBinsComb];
    Double_t errorYStat[nBinsComb];
    Double_t errorYSys[nBinsComb];
//        cout << nBinsNew-1 << "\t" <<nBinsComb << endl;
    for (Int_t commonBin = 0; commonBin < nBinsNew-1; commonBin++){
            cout << "A "<<commonBin<< "\t"  <<rebinnedA[commonBin].valueX << "\t +- " << rebinnedA[commonBin].errorXLow << " \t " << rebinnedA[commonBin].valueY<< "\t+" << rebinnedA[commonBin].errorYStatHigh<<"\t-"<< rebinnedA[commonBin].errorYStatLow<< "\t+" << rebinnedA[commonBin].errorYSystHigh<<"\t-"<< rebinnedA[commonBin].errorYSystLow<< "\t+" << rebinnedA[commonBin].errorYTotHigh<<"\t-"<< rebinnedA[commonBin].errorYTotLow<< "\t" << rebinnedA[commonBin].errorYTotHigh/rebinnedA[commonBin].valueY*100 << "%"<<endl;
            cout << "B " <<commonBin<< "\t"  <<rebinnedB[commonBin].valueX << "\t +- " << rebinnedB[commonBin].errorXLow << " \t " << rebinnedB[commonBin].valueY<< "\t+" << rebinnedB[commonBin].errorYStatHigh<<"\t-"<< rebinnedB[commonBin].errorYStatLow<< "\t+" << rebinnedB[commonBin].errorYSystHigh<<"\t-"<< rebinnedB[commonBin].errorYSystLow<< "\t+" << rebinnedB[commonBin].errorYTotHigh<<"\t-"<< rebinnedB[commonBin].errorYTotLow<< "\t" << rebinnedB[commonBin].errorYTotHigh/rebinnedB[commonBin].valueY*100 << "%"<<endl;
        ratioX[commonBin]   = rebinnedA[commonBin].valueX;
        errorX[commonBin]   = rebinnedA[commonBin].errorXHigh;
        ratioY[commonBin]   = rebinnedA[commonBin].valueY/rebinnedB[commonBin].valueY;
        errorY[commonBin]   = TMath::Sqrt(TMath::Power(rebinnedA[commonBin].errorYTotHigh/rebinnedA[commonBin].valueY,2) +TMath::Power(rebinnedB[commonBin].errorYTotHigh/rebinnedB[commonBin].valueY,2))*ratioY[commonBin];
        errorYStat[commonBin]   = TMath::Sqrt(TMath::Power(rebinnedA[commonBin].errorYStatHigh/rebinnedA[commonBin].valueY,2) +TMath::Power(rebinnedB[commonBin].errorYStatHigh/rebinnedB[commonBin].valueY,2))*ratioY[commonBin];
        errorYSys[commonBin]    = (((rebinnedA[commonBin].errorYSystHigh/rebinnedA[commonBin].valueY) +(rebinnedB[commonBin].errorYSystHigh/rebinnedB[commonBin].valueY))/2)*ratioY[commonBin];
        cout << "Ratio: " << ratioX[commonBin] << "\t" <<  errorX[commonBin] << "\t" << ratioY[commonBin] << "\t" << errorY[commonBin] << "\t" << errorY[commonBin]/ratioY[commonBin]*100 << "%"<<endl;
        cout << "Ratio stat: " << ratioX[commonBin] << "\t" <<  errorX[commonBin] << "\t" << ratioY[commonBin] << "\t" << errorYStat[commonBin] << "\t" << errorYStat[commonBin]/ratioY[commonBin]*100 << "%"<<endl;
        cout << "Ratio sys: " << ratioX[commonBin] << "\t" <<  errorX[commonBin] << "\t" << ratioY[commonBin] << "\t" << errorYSys[commonBin] << "\t" << errorYStat[commonBin]/ratioY[commonBin]*100 << "%"<<endl;

    }


    Double_t rebinnedSpectrumAY[nBinsNew-1];
    Double_t rebinnedSpectrumAYStatErr[nBinsNew-1];
    Double_t rebinnedSpectrumAYSysErr[nBinsNew-1];
    Double_t rebinnedSpectrumBY[nBinsNew-1];
    Double_t rebinnedSpectrumBYStatErr[nBinsNew-1];
    Double_t rebinnedSpectrumBYSysErr[nBinsNew-1];
    for (Int_t commonBin = 0; commonBin < nBinsNew-1; commonBin++){
        rebinnedSpectrumAY[commonBin]           = rebinnedA[commonBin].valueY;
        rebinnedSpectrumBY[commonBin]           = rebinnedB[commonBin].valueY;
        rebinnedSpectrumAYStatErr[commonBin]    = rebinnedA[commonBin].errorYStatHigh;
        rebinnedSpectrumBYStatErr[commonBin]    = rebinnedB[commonBin].errorYStatHigh;
        rebinnedSpectrumAYSysErr[commonBin]     = TMath::Sqrt(rebinnedA[commonBin].errorYSystHigh*rebinnedA[commonBin].errorYSystHigh -rebinnedA[commonBin].valueY*rebinnedA[commonBin].valueY) ;
        rebinnedSpectrumBYSysErr[commonBin]     = rebinnedB[commonBin].errorYSystHigh;
    }

    (*graphRebinnedAStat)   = new TGraphErrors(nBinsNew-1,ratioX,rebinnedSpectrumAY,errorX,rebinnedSpectrumAYStatErr);
    (*graphRebinnedASyst)   = new TGraphErrors(nBinsNew-1,ratioX,rebinnedSpectrumAY,errorX,rebinnedSpectrumAYSysErr);
    (*graphRebinnedBStat)   = new TGraphErrors(nBinsNew-1,ratioX,rebinnedSpectrumBY,errorX,rebinnedSpectrumBYStatErr);
    (*graphRebinnedBSyst)   = new TGraphErrors(nBinsNew-1,ratioX,rebinnedSpectrumBY,errorX,rebinnedSpectrumBYSysErr);

    if (calcSepErrGraphs){
        (*graphRatioStat)   = new TGraphErrors(nBinsNew-1,ratioX,ratioY,errorX,errorYStat);
        (*graphRatioSys)    = new TGraphErrors(nBinsNew-1,ratioX,ratioY,errorX,errorYSys);
    }
    TGraphErrors* returnGraph                   =  new TGraphErrors(nBinsNew-1,ratioX,ratioY,errorX,errorY);
    return returnGraph;

}

