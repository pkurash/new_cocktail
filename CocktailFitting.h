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
#include <fstream>
#include <math.h>
#include <string>
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
#include "TObjString.h"
#include "TFormula.h"
#include "TSpline.h"

//================================================================================================================
// definition of the fields in the histogram returned
//================================================================================================================
enum EValue_t {
    kYield = 1,
    kYieldStat,
    kYieldSysHi,
    kYieldSysLo,
    kMean,
    kMeanStat,
    kMeanSysHi,
    kMeanSysLo,
    kExtra
};

//================================================================================================================
//Read from settings files by initializer
//================================================================================================================
Double_t initParameter[nParticles][nCentralities][nMethods][15];
Double_t initParameterLimit[nParticles][nCentralities][nMethods][15][2];
Double_t initXRange[nParticles][nCentralities][nMethods][2];
Double_t initXRangeData[nParticles][nCentralities][nMethods][2];
Double_t initXRangeInteg[nParticles][nCentralities][nMethods][2];
TString initFitOption[nParticles][nCentralities][nMethods];
TString initFitFunction[nParticles][nCentralities][nMethods];
TString initDataToBeFitted[nParticles][nCentralities][nMethods];

Double_t initRatioParameter[nParticles][nParticles][nCentralities][nMethods][30];
Double_t initRatioParameterLimit[nParticles][nParticles][nCentralities][nMethods][30][2];
Double_t initRatioXRange[nParticles][nParticles][nCentralities][nMethods][2];
TString initRatioFitOption[nParticles][nParticles][nCentralities][nMethods];
TString initRatioFitFunction[nParticles][nParticles][nCentralities][nMethods];

Double_t initPtConstRelSyst[nParticles][nCentralities][nMethods];
Double_t initRatioPtConstRelSyst[nParticles][nParticles][nCentralities][nMethods];
//Double_t initPtConstSyst[nParticles][nCentralities][nMethods][2];

//================================================================================================================
//Set parameter limits
//================================================================================================================
TF1* SetParameterLimits(TF1* fct, Double_t* parLimitsDown, Double_t* parLimitsUp) {

    if (parLimitsDown && parLimitsUp) {
        for (Int_t i=0; i<fct->GetNpar(); i++) {
            if (parLimitsDown[i]!=-9999 && parLimitsUp[i]!=-9999) fct->SetParLimits(i, parLimitsDown[i], parLimitsUp[i]);
        }
    }

    return fct;
}

//================================================================================================================
//Multiply two TF1s
//================================================================================================================
TF1* MultiplyTF1(TF1* f1, TF1* f2, TString name) {

    if (!f1 || !f2) return NULL;

    Double_t xmin, xmax;
    f1->GetRange(xmin, xmax);
    Int_t nPar1                         = f1->GetNpar();
    Int_t nPar2                         = f2->GetNpar();
    TString formula1                    = f1->GetExpFormula();
    TString formula2                    = f2->GetExpFormula();

    for (Int_t i = 0;   i<nPar2;        i++) formula2.ReplaceAll(Form("[%d]",i),    Form("[-%d-]",i+nPar1));
    for (Int_t i=nPar1; i<nPar1+nPar2;  i++) formula2.ReplaceAll(Form("[-%d-]",i),  Form("[%d]",i));

    TF1* result = new TF1(name.Data(),Form("(%s)*(%s)",formula1.Data(), formula2.Data()), xmin, xmax);
    for (Int_t i = 0; i < nPar1; i++ ){
        result->SetParameter(i, f1->GetParameter(i));
    }
    for (Int_t j = 0; j < nPar2; j++ ){
        result->SetParameter(nPar1+j, f2->GetParameter(j));
    }

    return result;
}

//================================================================================================================
//Divide two TF1s
//================================================================================================================
TF1* DivideTF1(TF1* f1, TF1* f2, TString name) {

    if (!f1 || !f2) return NULL;

    Double_t xmin, xmax;
    f1->GetRange(xmin, xmax);
    Int_t nPar1                         = f1->GetNpar();
    Int_t nPar2                         = f2->GetNpar();
    TString formula1                    = f1->GetExpFormula();
    TString formula2                    = f2->GetExpFormula();

    for (Int_t i = 0;   i<nPar2;        i++) formula2.ReplaceAll(Form("[%d]",i),    Form("[-%d-]",i+nPar1));
    for (Int_t i=nPar1; i<nPar1+nPar2;  i++) formula2.ReplaceAll(Form("[-%d-]",i),  Form("[%d]",i));

    TF1* result = new TF1(name.Data(),Form("(%s)/(%s)",formula1.Data(), formula2.Data()), xmin, xmax);
    for (Int_t i = 0; i < nPar1; i++ ){
        result->SetParameter(i, f1->GetParameter(i));
    }
    for (Int_t j = 0; j < nPar2; j++ ){
        result->SetParameter(nPar1+j, f2->GetParameter(j));
    }

    return result;
}

//================================================================================================================
//Initialize options for spectra parametrization
//================================================================================================================
Bool_t InitializeFitting(TString paramSettingsFileName = "") {

    // initialize arrays that will be read from parametrization setting file
    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t j=0; j<nCentralities; j++) {
            for (Int_t k=0; k<nMethods; k++) {
                for (Int_t l=0; l<15; l++) {
                    initParameter[i][j][k][l]               = -9999;
                    initParameterLimit[i][j][k][l][0]       = -9999;
                    initParameterLimit[i][j][k][l][1]       = -9999;
                }
                initXRange[i][j][k][0]                      = -9999;
                initXRange[i][j][k][1]                      = -9999;
                initFitOption[i][j][k]                      = "";
                initFitFunction[i][j][k]                    = "";
                initPtConstRelSyst[i][j][k]                 = 0.;
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

        // get pt const rel sys err
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initPtConstRelSyst[i][j][k]                 = std::stod(tempString.Data());

        // get x range
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initXRange[i][j][k][0]                      = std::stod(tempString.Data());
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initXRange[i][j][k][1]                      = std::stod(tempString.Data());

        cout << "fit range: " << initXRange[i][j][k][0] << "\t" << initXRange[i][j][k][1] << endl;

        // get func
        fileParamSettings >> initFitFunction[i][j][k];

        // get fit option
        fileParamSettings >> initFitOption[i][j][k];
        cout << initFitFunction[i][j][k].Data() << "\t" << initFitOption[i][j][k].Data() << endl;

        // get parameter and parameter limits
        Int_t counter                                                                   = 0;
        Int_t counterTotal                                                              = 0;
        fileParamSettings >> tempString;
        while (tempString.CompareTo("stop")!=0 && counterTotal < 15) {
            if (tempString.CompareTo("-")!=0) {
                if (counter == 0)       initParameter[i][j][k][counterTotal]            = std::stod(tempString.Data());
                else if (counter==1)    initParameterLimit[i][j][k][counterTotal][0]    = std::stod(tempString.Data());
                else if (counter==2)    initParameterLimit[i][j][k][counterTotal][1]    = std::stod(tempString.Data());
            }
            fileParamSettings >> tempString;

            counter++;
            if (counter>2) counterTotal++;
            if (counter>2)  counter = 0;
        }
    }

    return kTRUE;
}

//================================================================================================================
//Initialize options for spectra parametrization for integrated yield
//================================================================================================================
Bool_t InitializeFittingIntegYield(TString paramSettingsFileName = "") {

    // initialize arrays that will be read from parametrization setting file
    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t j=0; j<nCentralities; j++) {
            for (Int_t k=0; k<nMethods; k++) {
                for (Int_t l=0; l<15; l++) {
                    initParameter[i][j][k][l]               = -9999;
                    initParameterLimit[i][j][k][l][0]       = -9999;
                    initParameterLimit[i][j][k][l][1]       = -9999;
                }
                initXRange[i][j][k][0]                      = -9999;
                initXRange[i][j][k][1]                      = -9999;
                initXRangeData[i][j][k][0]                  = -9999;
                initXRangeData[i][j][k][1]                  = -9999;
                initXRangeInteg[i][j][k][0]                 = -9999;
                initXRangeInteg[i][j][k][1]                 = -9999;
                initFitOption[i][j][k]                      = "";
                initFitFunction[i][j][k]                    = "";
                initDataToBeFitted[i][j][k]                 = "";
                initPtConstRelSyst[i][j][k]                 = 0.;
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

        // get pt const rel sys err
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initDataToBeFitted[i][j][k]                 = tempString.Data();
        else initDataToBeFitted[i][j][k]                                                = "Stat";

        // get pt const rel sys err
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initPtConstRelSyst[i][j][k]                 = std::stod(tempString.Data());

        // get x range
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initXRange[i][j][k][0]                      = std::stod(tempString.Data());
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initXRange[i][j][k][1]                      = std::stod(tempString.Data());

        cout << "fit range: " << initXRange[i][j][k][0] << "\t" << initXRange[i][j][k][1] << endl;

        // get integ x range
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initXRangeData[i][j][k][0]                 = std::stod(tempString.Data());
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initXRangeData[i][j][k][1]                 = std::stod(tempString.Data());
        cout << "integ range: " << initXRangeData[i][j][k][0] << "\t" << initXRangeData[i][j][k][1] << endl;

        // get integ x range
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initXRangeInteg[i][j][k][0]                 = std::stod(tempString.Data());
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initXRangeInteg[i][j][k][1]                 = std::stod(tempString.Data());

        cout << "integ range: " << initXRangeInteg[i][j][k][0] << "\t" << initXRangeInteg[i][j][k][1] << endl;

        // get func
        fileParamSettings >> initFitFunction[i][j][k];

        // get fit option
        fileParamSettings >> initFitOption[i][j][k];
        cout << initFitFunction[i][j][k].Data() << "\t" << initFitOption[i][j][k].Data() << endl;

        // get parameter and parameter limits
        Int_t counter                                                                   = 0;
        Int_t counterTotal                                                              = 0;
        fileParamSettings >> tempString;
        while (tempString.CompareTo("stop")!=0 && counterTotal < 15) {
            if (tempString.CompareTo("-")!=0) {
                if (counter == 0)       initParameter[i][j][k][counterTotal]            = std::stod(tempString.Data());
                else if (counter==1)    initParameterLimit[i][j][k][counterTotal][0]    = std::stod(tempString.Data());
                else if (counter==2)    initParameterLimit[i][j][k][counterTotal][1]    = std::stod(tempString.Data());
            }
            fileParamSettings >> tempString;

            counter++;
            if (counter>2) counterTotal++;
            if (counter>2)  counter = 0;
        }
    }

    return kTRUE;
}

//================================================================================================================
//Initialize options for ratio parametrization
//================================================================================================================
Bool_t InitializeRatioFitting(TString paramSettingsFileName = "") {

    // initialize arrays that will be read from parametrization setting file
    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t j=0; j<nParticles; j++) {
            for (Int_t k=0; k<nCentralities; k++) {
                for (Int_t l=0; l<nMethods; l++) {
                    for (Int_t m=0; m<30; m++) {
                        initRatioParameter[i][j][k][l][m]           = -9999;
                        initRatioParameterLimit[i][j][k][l][m][0]   = -9999;
                        initRatioParameterLimit[i][j][k][l][m][1]   = -9999;
                    }
                    initRatioXRange[i][j][k][l][0]                  = -9999;
                    initRatioXRange[i][j][k][l][1]                  = -9999;
                    initRatioFitOption[i][j][k][l]                  = "";
                    initRatioFitFunction[i][j][k][l]                = "";
                    initRatioPtConstRelSyst[i][j][k][l]             = 0.;
                }
            }
        }
    }

    // read in parametrization settings file
    ifstream fileParamSettings;
    fileParamSettings.open(paramSettingsFileName,ios_base::in);
    if (!fileParamSettings) {
        cout << "ERROR: Ratio parametrization settings " << paramSettingsFileName.Data() << " not found!" << endl;
        return kFALSE;
    }

    // read settings from file
    TString particleRatioFromFile, particle1, particle2, centralityFromFile, methodFromFile;
    TString tempString;
    Int_t i,j,k,l;
    std::string line;
    for( std::string line; getline(fileParamSettings, line); ) {

        // get basic settings
        fileParamSettings >> particleRatioFromFile >> centralityFromFile >> methodFromFile;
        particle1                                                                               = particleRatioFromFile(0,particleRatioFromFile.Index("To"));
        particle2                                                                               = particleRatioFromFile(particleRatioFromFile.Index("To")+2,particleRatioFromFile.Length()-particleRatioFromFile.Index("To")+2);
        if (methodFromFile.CompareTo("-") == 0)     methodFromFile                              = "";
        if (centralityFromFile.CompareTo("-") == 0) centralityFromFile                          = "MB";
        i                                                                                       = GetParticleIterator(particle1);
        j                                                                                       = GetParticleIterator(particle2);
        k                                                                                       = GetCentralityIterator(centralityFromFile);
        l                                                                                       = GetMethodIterator(methodFromFile);
        if (i<0 || j<0 || k<0 || l<0) continue;

        // get pt const rel sys err
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initRatioPtConstRelSyst[i][j][k][l]                 = std::stod(tempString.Data());

        // get x range
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initRatioXRange[i][j][k][l][0]                      = std::stod(tempString.Data());
        fileParamSettings >> tempString;
        if (tempString.CompareTo("-") != 0) initRatioXRange[i][j][k][l][1]                      = std::stod(tempString.Data());

        // get func
        fileParamSettings >> initRatioFitFunction[i][j][k][l];

        // get fit option
        fileParamSettings >> initRatioFitOption[i][j][k][l];

        // get parameter and parameter limits
        Int_t counter                                                                           = 0;
        Int_t counterTotal                                                                      = 0;
        fileParamSettings >> tempString;
        while (tempString.CompareTo("stop")!=0 && counterTotal < 30) {
            if (tempString.CompareTo("-")!=0) {
                if (counter == 0)       initRatioParameter[i][j][k][l][counterTotal]            = std::stod(tempString.Data());
                else if (counter==1)    initRatioParameterLimit[i][j][k][l][counterTotal][0]    = std::stod(tempString.Data());
                else if (counter==2)    initRatioParameterLimit[i][j][k][l][counterTotal][1]    = std::stod(tempString.Data());
            }
            fileParamSettings >> tempString;

            counter++;
            if (counter>2) counterTotal++;
            if (counter>2) counter = 0;
        }
    }

    return kTRUE;
}

//================================================================================================================
// Define different Blastwave functions - BOLTZMANN-GIBBS Blast-wave
//================================================================================================================
static TF1 *fBGBlastWave_Integrand = NULL;

//**********************************************************
// define integrand
//**********************************************************
Double_t BGBlastWave_Integrand(  const Double_t *x,
                                 const Double_t *p){

    /*
     *     x[0] -> r (radius)
     *     p[0] -> mT (transverse mass)
     *     p[1] -> pT (transverse momentum)
     *     p[2] -> beta_max (surface velocity)
     *     p[3] -> T (freezout temperature)
     *     p[4] -> n (velocity profile)
     */

    Double_t r = x[0];
    Double_t mt = p[0];
    Double_t pt = p[1];
    Double_t beta_max = p[2];
    Double_t temp_1 = 1. / p[3];
    Double_t n = p[4];

    Double_t beta = beta_max * TMath::Power(r, n);
    if (beta > 0.9999999999999999) beta = 0.9999999999999999;
    Double_t rho = TMath::ATanH(beta);
    Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
    if (argI0 > 700.) argI0 = 700.;
    Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
    //  if (argI0 > 100 || argI0 < -100)
    //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
    return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);

}

//**********************************************************
// define functional form
//**********************************************************
Double_t BGBlastWave_Func( const Double_t *x,
                           const Double_t *p ){
    /* dN/dpt */

    Double_t pt         = x[0];
    Double_t mass       = p[0];
    Double_t mt         = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max   = p[1];
    Double_t temp       = p[2];
    Double_t n          = p[3];
    Double_t norm       = p[4];

    if (!fBGBlastWave_Integrand)
        fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
    cout << "here" << endl;
    //Double_t integral   = fBGBlastWave_Integrand->Integral(0., 1., (Double_t *)0, 1.e-6);
    Double_t integral   = fBGBlastWave_Integrand->Integral(0., 1., 1.e-6);
    return norm * pt * integral;
}

//================================================================================================================
// Define fit function using TSpline
//================================================================================================================
Double_t methodTSplineTF1(Double_t *x, Double_t *par)
{
  if((Int_t)vecGraphs.size()<counterGraphs){
    cout << "ERROR in methodTSplineTF1: proper TGraph has not been defined, no TSpline fit possible!" << endl;
    return 0;
  }else{
    return vecGraphs.at((Int_t)par[0])->Eval(x[0]);
  }
}
//================================================================================================================
// Define fit function using TSpline+PowerLaw
//================================================================================================================
Double_t methodTSplinePowerLawTF1(Double_t *x, Double_t *par)
{
  if((Int_t)vecGraphs.size()<counterGraphs){
    cout << "ERROR in methodTSplinePowerLawTF1: proper TGraph has not been defined, no TSpline fit possible!" << endl;
    return 0;
  }else{
    if(x[0]<par[1]) return vecGraphs.at((Int_t)par[0])->Eval(x[0]);
    else return par[2] * 1/pow(x[0],par[3]);
  }
}
//================================================================================================================
//FitObject
//================================================================================================================
TF1* FitObject( TObject *Obj_Dummy      = NULL,
                TString FunctionName    = "",
                TString mesonType       = "Pi0",
                TString centrality      = "",
                TString method          = "",
                TString type            = "",
                Bool_t convertToHist    = kFALSE
               ) {

    TString ClassName                                                   = Obj_Dummy->ClassName();
    if(convertToHist){
        if (ClassName.Contains("TH1"))
            cout << "is histo already" << endl;
        else if (ClassName.CompareTo("TGraphErrors") == 0)
            Obj_Dummy = GraphToHist_withErrors((TGraphErrors*)Obj_Dummy, "Obj_Dummy_Hist");
        else if (ClassName.CompareTo("TGraphAsymmErrors") == 0)
            Obj_Dummy = GraphToHist_withErrors((TGraphAsymmErrors*)Obj_Dummy, "Obj_Dummy_Hist");
        ClassName = Obj_Dummy->ClassName();
    }
    // xRange
    Double_t xmin                                                       = initXRange[GetParticleIterator(mesonType)][GetCentralityIterator(centrality)][GetMethodIterator(method)][0];
    Double_t xmax                                                       = initXRange[GetParticleIterator(mesonType)][GetCentralityIterator(centrality)][GetMethodIterator(method)][1];
    if (xmin == -9999) {
        if (ClassName.Contains("TH1"))                          xmin    = ((TH1D*)Obj_Dummy)->GetXaxis()->GetXmin();
        else if (ClassName.CompareTo("TGraphErrors") == 0)      xmin    = GetXRangeFromGraph((TGraphErrors*)Obj_Dummy, kFALSE);
        else if (ClassName.CompareTo("TGraphAsymmErrors") == 0) xmin    = GetXRangeFromGraph((TGraphAsymmErrors*)Obj_Dummy, kFALSE);
    }
    if (xmax == -9999) {
        if (ClassName.Contains("TH1"))                          xmax    = ((TH1D*)Obj_Dummy)->GetXaxis()->GetXmax();
        else if (ClassName.CompareTo("TGraphErrors") == 0)      xmax    = GetXRangeFromGraph((TGraphErrors*)Obj_Dummy, kTRUE);
        else if (ClassName.CompareTo("TGraphAsymmErrors") == 0) xmax    = GetXRangeFromGraph((TGraphAsymmErrors*)Obj_Dummy, kTRUE);
    }

    // fit option
    TString FitOptions                                                  = initFitOption[GetParticleIterator(mesonType)][GetCentralityIterator(centrality)][GetMethodIterator(method)];

    // start parameter
    Double_t *Parameter                                                 = NULL;
    Int_t nPars                                                         = 0;
    for (Int_t i=0; i<15; i++) {
        if (initParameter[GetParticleIterator(mesonType)][GetCentralityIterator(centrality)][GetMethodIterator(method)][i] != -9999) nPars++;
    }
    if (nPars>0) {
        Parameter                                                       = new Double_t[nPars];
        for (Int_t i=0; i<nPars; i++)
            Parameter[i]                                                = initParameter[GetParticleIterator(mesonType)][GetCentralityIterator(centrality)][GetMethodIterator(method)][i];
    }

    // start parameter limits (if -9999, don't use)
    Double_t *ParameterLimitsDown                                       = NULL;
    Double_t *ParameterLimitsUp                                         = NULL;
    if (nPars) {
        ParameterLimitsDown                                             = new Double_t[nPars];
        ParameterLimitsUp                                               = new Double_t[nPars];
        for (Int_t i=0; i<nPars; i++) {
            ParameterLimitsDown[i]                                      = initParameterLimit[GetParticleIterator(mesonType)][GetCentralityIterator(centrality)][GetMethodIterator(method)][i][0];
            ParameterLimitsUp[i]                                        = initParameterLimit[GetParticleIterator(mesonType)][GetCentralityIterator(centrality)][GetMethodIterator(method)][i][1];
        }
    }

    // mass
    Double_t mass;
    if (mesonType.CompareTo(fParticle[0]) == 0){
        // Pi0
        mass = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    } else if (mesonType.CompareTo(fParticle[1]) == 0){
        // Eta
        mass = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    } else if (mesonType.CompareTo(fParticle[2]) == 0){
        // Omega
        mass = TDatabasePDG::Instance()->GetParticle(223)->Mass();
    } else if (mesonType.CompareTo(fParticle[3]) == 0){
        // EtaPrima
        mass = TDatabasePDG::Instance()->GetParticle(331)->Mass();
    } else if (mesonType.CompareTo(fParticle[4]) == 0){
        // GammaDir
        mass = 0;
    } else if (mesonType.CompareTo(fParticle[5]) == 0){
        // CPion
        mass = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    } else if (mesonType.CompareTo(fParticle[6]) == 0){
        // CKaon
        mass = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    } else if (mesonType.CompareTo(fParticle[7]) == 0){
        // Proton
        mass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    } else if (mesonType.CompareTo(fParticle[8]) == 0){
        // CHadron
        mass = 0;
    } else if (mesonType.CompareTo(fParticle[9]) == 0){
        // Phi
        mass = TDatabasePDG::Instance()->GetParticle(333)->Mass();
    } else if (mesonType.CompareTo(fParticle[10]) == 0){
        // NKaonStar
        mass = TDatabasePDG::Instance()->GetParticle(313)->Mass();
    } else if (mesonType.CompareTo(fParticle[11]) == 0){
        // NRho
        mass = TDatabasePDG::Instance()->GetParticle(113)->Mass();
    } else if (mesonType.CompareTo(fParticle[12]) == 0){
        // CRho
        mass = TDatabasePDG::Instance()->GetParticle(213)->Mass();
    } else if (mesonType.CompareTo(fParticle[13]) == 0){
        // NDelta
        mass = TDatabasePDG::Instance()->GetParticle(2114)->Mass();
    } else if (mesonType.CompareTo(fParticle[14]) == 0){
        // CDelta
        mass = TDatabasePDG::Instance()->GetParticle(2214)->Mass();     // Delta^+
    } else if (mesonType.CompareTo(fParticle[15]) == 0){
        // NKaonSubS
        mass = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    } else if (mesonType.CompareTo(fParticle[16]) == 0){
        // Lambda
        mass = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    } else if (mesonType.CompareTo(fParticle[17]) == 0){
        // NSigma
        mass = TDatabasePDG::Instance()->GetParticle(3212)->Mass();
    } else if (mesonType.CompareTo(fParticle[18]) == 0){
        // CSigma
        mass = TDatabasePDG::Instance()->GetParticle(3222)->Mass();
    } else if (mesonType.CompareTo(fParticle[19]) == 0){
        // COmega
        mass = TDatabasePDG::Instance()->GetParticle(3334)->Mass();
    } else if (mesonType.CompareTo(fParticle[20]) == 0){
        // CXi
        mass = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
    } else if (mesonType.CompareTo(fParticle[21]) == 0){
        // JPsi
        mass = TDatabasePDG::Instance()->GetParticle(443)->Mass();
    } else {
        mass = 0;
    }
    cout <<"Meson Mass: "<< mass << endl;

    if(Obj_Dummy == NULL) {
        if(type.CompareTo("h") == 0 || type.CompareTo("H") == 0){
            TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*pow([2]/([2]+x),[1])");
            Hagedorn_Dummy->SetParNames("C_{H}","n","p_{0} (GeV/c)");
            Hagedorn_Dummy->SetParameters(19.,6.8,0.84);
            Hagedorn_Dummy->SetName(FunctionName);
            return Hagedorn_Dummy;
        }
        if(type.CompareTo("powPure") == 0 || type.CompareTo("PowPure") == 0){
            cout <<Form("fitting %s with Pure Powerlaw",FunctionName.Data()) << endl;
            TF1 *PowerLaw_Dummy = new TF1("PowerLawPure_Dummy","[0] * 1/pow(x,[1])");
            PowerLaw_Dummy->SetParNames("A_{pow}","n");
            PowerLaw_Dummy->SetParameters(2.,5.);
            PowerLaw_Dummy->SetName(FunctionName);
            return PowerLaw_Dummy;
        } else
        if(type.CompareTo("p") == 0 || type.CompareTo("P") == 0){
            TF1 *PowerLaw_Dummy = new TF1("PowerLaw_Dummy","[0] * 2 / TMath::Pi() * ([1]-1.)*([1]-2.)/pow([1]-3.,2) /x * pow(1+2*x/[2]/([1]-3),-[1])");
            PowerLaw_Dummy->SetParNames("A_{pow}","n","p_{0}");
            PowerLaw_Dummy->SetParameters(2.,5.,0.37);
            PowerLaw_Dummy->SetName(FunctionName);
            return PowerLaw_Dummy;
        }
        if(type.CompareTo("l") == 0 || type.CompareTo("L") == 0){
            TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
            Levy_Dummy->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
            Levy_Dummy->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
            Levy_Dummy->SetName(FunctionName);
            return Levy_Dummy;
        }
        if(type.CompareTo("b") == 0 || type.CompareTo("B") == 0){
            TF1 *Boltzmann_Dummy = new TF1("Boltzmann_Dummy",Form("[0] *sqrt(x*x+%.10f*%.10f)* exp(- sqrt(x*x+%.10f*%.10f)/[1])",mass,mass,mass,mass));
            Boltzmann_Dummy->SetParNames("C_{B}","T_{Boltzm.} (GeV/c)");
            Boltzmann_Dummy->SetParameters(2.5,0.3); // standard parameter optimize if necessary
            Boltzmann_Dummy->SetName(FunctionName);
            return Boltzmann_Dummy;
        }
        if(type.CompareTo("e") == 0 || type.CompareTo("E") == 0){
            TF1 *Exponential_Dummy = new TF1("Exponential_Dummy",Form("[0]/( 2 * TMath::Pi())/([1]*(%.10f+[1]))*exp(-(sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass,mass));
            Exponential_Dummy->SetParameters(1.1,0.37); // standard parameter optimize if necessary
            Exponential_Dummy->SetParNames("C_{E}","T_{exp} (GeV/c)");
            Exponential_Dummy->SetName(FunctionName);
            return Exponential_Dummy;
        }
        if(type.CompareTo("m") == 0 || type.CompareTo("M") == 0){
            TF1 *ModPowerLaw_Dummy = new TF1("ModPowerLaw_Dummy","[0]*pow((1 + (x)/[1]),-[2])");
            ModPowerLaw_Dummy->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
            ModPowerLaw_Dummy->SetParNames("A","p_{0}","n");
            ModPowerLaw_Dummy->SetName(FunctionName);
            return ModPowerLaw_Dummy;
        }
        if(type.CompareTo("6pol") == 0 || type.CompareTo("6POL") == 0){
            TF1 *Pol6_Dummy = new TF1("Pol6_Dummy","[0]+[1]*pow(x,2)+[2]*pow(x,4)+[3]*pow(x,6)");
            Pol6_Dummy->SetParameters(1.,1.,1.,1.); // standard parameter optimize if necessary
            Pol6_Dummy->SetParNames("a","b","c","d");
            Pol6_Dummy->SetName(FunctionName);
            return Pol6_Dummy;
        }
        if(type.CompareTo("doubqcd") == 0 || type.CompareTo("doubqcd") == 0){
            cout <<Form("fitting %s with doubqcd",FunctionName.Data()) << endl;
            TF1 *QCD_Dummy = new TF1("QCD_Dummy","(x<=[5])*[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))+(x>[5])*[6]*TMath::Power(x,-1*([7]+[8]/(TMath::Power(x,[9])+[10])))");
            QCD_Dummy->SetParameters(24,6.7,-6.5,1.,10,24,2,6.7,-6.5,1.,10); // standard parameter optimize if necessary
            QCD_Dummy->SetParNames("a1","b1","c1","d1","e1","pT","a2","b2","c2","d2","e2");
            QCD_Dummy->SetName(FunctionName);
            return QCD_Dummy;
        }
        if(type.CompareTo("qcdtsal") == 0 || type.CompareTo("qcdtsal") == 0){
            cout <<Form("fitting %s with qcdtsal",FunctionName.Data()) << endl;
            TF1 *QCD_Dummy = new TF1("QCD_Dummy","(x<=[5])*[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))+(x>[5])*[6] / ( 2 * TMath::Pi())*([7]-1.)*([7]-2.) / ([7]*[8]*([7]*[8]+%.10f*([7]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([7]*[8]), -[7])");
            QCD_Dummy->SetParameters(24,6.7,-6.5,1.,10,2,2.,5.,0.18); // standard parameter optimize if necessary
            QCD_Dummy->SetParNames("a1","b1","c1","d1","e1","pT","a2","b2","c2");
            QCD_Dummy->SetName(FunctionName);
            return QCD_Dummy;
        }
        if(type.CompareTo("rad") == 0 || type.CompareTo("RAD") == 0){
            TF1 *Rad_Dummy = new TF1("Rad_Dummy","(x<=[3])*x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+0.135*([1]-2.)))*pow(1.+(sqrt(x*x+0.135*0.135)-0.135)/([1]*[2]), -[1]) + (x>[3]) * x *[0] * ([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+0.135*([1]-2.))) * pow([1]*[2]/([3]+[1]*[2]-0.135),[1]) * pow([3],[4]) * pow( 1./x, [4] )");
            Rad_Dummy->SetLineWidth(1);
            Double_t par0min3 = 150.;
            Double_t par0max3 = 1500.;
            Double_t par1min3 = 4.;
            Double_t par1max3 = 40.;
            Double_t par2min3 = 0.007;
            Double_t par2max3 = .7;
            Double_t par3min3 = 4.36;
            Double_t par3max3 = 4.40;
            Double_t par4min3 = 2.;
            Double_t par4max3 = 18.;
            Rad_Dummy->SetParLimits(0, par0min3, par0max3);
            Rad_Dummy->SetParLimits(1, par1min3, par1max3);
            Rad_Dummy->SetParLimits(2, par2min3, par2max3);
            Rad_Dummy->SetParLimits(3, par3min3, par3max3);
            Rad_Dummy->SetParLimits(4, par4min3, par4max3);
            Rad_Dummy->SetParNames("a","b","c","d");
            Rad_Dummy->SetName(FunctionName);
            return Rad_Dummy;
        }
        if(type.CompareTo("tmpt") == 0 || type.CompareTo("TMPT") == 0){
            // Tsallis Dummy multiplied with pt
            TF1 *Levy_Dummy = new TF1("Tsallis_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * x* pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
            Levy_Dummy->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
            Levy_Dummy->SetParameters(2.e11,7., 0.137) ; // standard parameter optimize if necessary
            Levy_Dummy->SetName(FunctionName);
            return Levy_Dummy;
        }
        if(type.CompareTo("hmpt") == 0 || type.CompareTo("HMPT") == 0){
            // Hagedorn Dummy multiplied with pt
            TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*x*pow([2]/([2]+x),[1])");
            Hagedorn_Dummy->SetParNames("C_{H}","n","p_{0} (GeV/c)");
            Hagedorn_Dummy->SetParameters(1.,7.,0.37);
            Hagedorn_Dummy->SetName(FunctionName);
            return Hagedorn_Dummy;
        }
        if(type.CompareTo("qmpt") == 0 || type.CompareTo("QMPT") == 0){
            TF1 *QCD_Dummy = new TF1("QCD_Dummy","([0]*x*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4]))))");
            QCD_Dummy->SetParNames("a","b","c","d","e");
            QCD_Dummy->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
            QCD_Dummy->SetName(FunctionName);
            return QCD_Dummy;
        }
        if(type.CompareTo("doHag") == 0 || type.CompareTo("DOHag") == 0){
            cout <<Form("fitting %s with mHag and pow law",FunctionName.Data()) << endl;
            TF1 *HagPow_Dummy = new TF1("HagPow_Dummy","(x<=[7])*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])+(x>[7])*[5]*1/pow(x,[6])");
            HagPow_Dummy->SetParNames("a","b","c","d","e","A","n","pt");
            HagPow_Dummy->SetParameters(30.,0.37,0.07,0.68,6.1,2.,5.,5.); // standard parameter optimize if necessary
            HagPow_Dummy->SetName(FunctionName);
            return HagPow_Dummy;
        }
        if(type.CompareTo("oHag") == 0 || type.CompareTo("OHag") == 0){
            cout << "entered"<< endl;
            TF1 *ModPowerLaw_Dummy2 = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])");
            ModPowerLaw_Dummy2->SetParNames("a","b","c","d","e");
            ModPowerLaw_Dummy2->SetParameters(30.,0.37,0.07,0.68,6.1);
            ModPowerLaw_Dummy2->SetName(FunctionName);
            return ModPowerLaw_Dummy2;
        }
        if(type.CompareTo("mohag") == 0 || type.CompareTo("MOHAG") == 0){
            cout <<Form("fitting %s with ModHagedorn",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*x*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])");
            //TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x)+x/[2],-[3])");
            ModPowerLaw_Dummy->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
            ModPowerLaw_Dummy->SetParNames("a","b","c","d","e");
            ModPowerLaw_Dummy->SetName(FunctionName);
            return ModPowerLaw_Dummy;
        }
        if(type.CompareTo("tcm") == 0 || type.CompareTo("TCM") == 0){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
            cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
            TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
            //             TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(1 + x*x/TMath::Power(1+x*x/([3]*[3]*[4]),-[4]) )",mass,mass,mass));
            if (mesonType.CompareTo("Pi0") == 0){
                TwoCompModel_Dummy->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
            } else if (mesonType.CompareTo("Eta") == 0){
                TwoCompModel_Dummy->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
            }
            TwoCompModel_Dummy->SetParNames("Ae","Te","A","T","n");
            TwoCompModel_Dummy->SetName(FunctionName);
            return TwoCompModel_Dummy;
        }
        if(type.CompareTo("tcmpt") == 0 || type.CompareTo("TCMPT") == 0){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
            cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
            TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form(" x*[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + x*[2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
            //             TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("x*[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + x*[2]/(1 + x*x/TMath::Power(1+x*x/([3]*[3]*[4]),-[4]) )",mass,mass,mass));
            if (mesonType.CompareTo("Pi0") == 0){
                TwoCompModel_Dummy->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
            } else if (mesonType.CompareTo("Eta") == 0){
                TwoCompModel_Dummy->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
            }
            TwoCompModel_Dummy->SetParNames("Ae","Te","A","T","n");
            TwoCompModel_Dummy->SetName(FunctionName);
            return TwoCompModel_Dummy;
        }
        if(type.CompareTo("kfunc") == 0 || type.CompareTo("KFUNC") == 0){
            // combine exp for flow and two Hagedorn functions
            TF1 *KFunc_Dummy = new TF1("KFunc_Dummy",Form("[0]*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[3]) + [4]/(1+x/[5])^[6] + [7]/(1+x/[8])^[9]",mass,mass),0,200);
            KFunc_Dummy->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            KFunc_Dummy->SetParameters(1e3,0.6,0.14,0.1,1e3,1,8,0.,1,4) ; // standard parameter optimize if necessary
            KFunc_Dummy->SetName(FunctionName);
            return KFunc_Dummy;
        }
        if(type.CompareTo("modkfunc") == 0 || type.CompareTo("MODKFUNC") == 0){
            // combine exp for flow and two Hagedorn functions
            TF1 *modKFunc_Dummy = new TF1("modKFunc_Dummy",Form("[0]*(sqrt(x*x+%.10f*%.10f)-x*[1])/sqrt(1-[1]*[1])*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2]) + [3]/(1+x/[4])^[5] + [6]/(1+x/[7])^[8]",mass,mass,mass,mass),0,200);
            modKFunc_Dummy->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            modKFunc_Dummy->SetParameters(1e3,0.6,0.1,1e3,1,8,0.,1,4); // standard parameter optimize if necessary
            modKFunc_Dummy->SetName(FunctionName);
            return modKFunc_Dummy;
        }
        if(type.CompareTo("modkfuncpt") == 0 || type.CompareTo("MODKFUNCPT") == 0){
            // combine exp for flow and two Hagedorn functions
            TF1 *modKFunc_Dummy = new TF1("modKFunc_Dummy",Form("x*([0]*(sqrt(x*x+%.10f*%.10f)-x*[1])/sqrt(1-[1]*[1])*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2]) + [3]/(1+x/[4])^[5] + [6]/(1+x/[7])^[8])",mass,mass,mass,mass),0,200);
            modKFunc_Dummy->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            modKFunc_Dummy->SetParameters(1e3,0.6,0.1,1e3,1,8,0.,1,4); // standard parameter optimize if necessary
            modKFunc_Dummy->SetName(FunctionName);
            return modKFunc_Dummy;
        }
        if (type.CompareTo("blastWave") == 0){
            TF1 *blastWave_Dummy = new TF1("blastWave_Dummy", BGBlastWave_Func, 0., 200., 5);
            blastWave_Dummy->SetParNames("mass", "beta_max", "T", "n", "norm");
            blastWave_Dummy->SetParameters(mass, 0.9, 0.1, 1, 1.e6);
            blastWave_Dummy->FixParameter(0, mass);
            blastWave_Dummy->SetParLimits(1, 0.01, 0.99);
            blastWave_Dummy->SetParLimits(2, 0.01, 1.);
            blastWave_Dummy->SetParLimits(3, 0.01, 50.);
            return blastWave_Dummy;
        }
    }

    TF1 *FitFunction = new TF1();

    if(ClassName.BeginsWith("TH1")){
        TH1D *Obj = (TH1D*)Obj_Dummy;
        if(type.CompareTo("xqcd") == 0 || type.CompareTo("XQCD") == 0){
            cout <<Form("fitting %s with XQCD",FunctionName.Data()) << endl;
            TF1 *XQCD_Dummy = new TF1("QCD_Dummy","[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4]*TMath::Power(x,0.5))))",0,200);
            FitFunction = (TF1*)XQCD_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10,0.5); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("a","b","c","d","e");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("xqcdpt") == 0 || type.CompareTo("XQCDPT") == 0){
            cout <<Form("fitting %s with XQCD times pt",FunctionName.Data()) << endl;
            TF1 *XQCD_Dummy = new TF1("QCD_Dummy","x*[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4]*TMath::Power(x,0.5))))",0,200);
            FitFunction = (TF1*)XQCD_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10,0.5); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("a","b","c","d","e");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("doubqcd") == 0 || type.CompareTo("DOUBQCD") == 0){
            cout <<Form("fitting %s with QCD",FunctionName.Data()) << endl;
            TF1 *QCD_Dummy = new TF1("QCD_Dummy","(x<=[5])*[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))+(x>[5])*[6]*TMath::Power(x,-1*([7]+[8]/(TMath::Power(x,[9])+[10])))",0,200);
            FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10,24,2,6.7,-6.5,1.,10); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7],Parameter[8],Parameter[9],Parameter[10]);
            FitFunction->SetParNames("a1","b1","c1","d1","e1","pT","a2","b2","c2","d2","e2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("qcd") == 0 || type.CompareTo("QCD") == 0){
            cout <<Form("fitting %s with QCD",FunctionName.Data()) << endl;
            TF1 *QCD_Dummy = new TF1("QCD_Dummy","[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))",0,200);
            FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("a","b","c","d","e");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("qcdpt") == 0 || type.CompareTo("QCDPT") == 0){
            cout <<Form("fitting %s with QCD times pt",FunctionName.Data()) << endl;
            TF1 *QCD_Dummy = new TF1("QCD_Dummy","x*[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))",0,200);
            FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("a","b","c","d","e");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("h") == 0 || type.CompareTo("H") == 0){
            cout <<Form("fitting %s with Hagedorn",FunctionName.Data()) << endl;
            TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*pow([2]/([2]+x),[1])",0,200);
            FitFunction = (TF1*)Hagedorn_Dummy->Clone(FunctionName);
            if(Parameter == NULL) FitFunction->SetParameters(19.,6.8,0.84);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("C_{H}","n","p_{0} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("p") == 0 || type.CompareTo("P") == 0){
            cout <<Form("fitting %s with Powerlaw",FunctionName.Data()) << endl;
            TF1 *PowerLaw_Dummy = new TF1("PowerLaw_Dummy","[0] * 2 / TMath::Pi() * ([1]-1.)*([1]-2.)/pow([1]-3.,2) /x * pow(1+2*x/[2]/([1]-3),-[1])",0,200);
            FitFunction = (TF1*)PowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.37); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("A_{pow}","n","p_{0}");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("l") == 0 || type.CompareTo("L") == 0){
            cout <<Form("fitting %s with Levy",FunctionName.Data()) << endl;
            TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("tmpt") == 0 || type.CompareTo("TMPT") == 0){
            // Tsallis Dummy multiplied with pt
            TF1 *Levy_Dummy = new TF1("Tsallis_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * x* pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("b") == 0 || type.CompareTo("B") == 0){
            cout <<Form("fitting %s with Boltzmann",FunctionName.Data()) << endl;
            TF1 *Boltzmann_Dummy =  new TF1("Boltzmann_Dummy",Form("[0] *sqrt(x*x+%.10f*%.10f)* exp(- sqrt(x*x+%.10f*%.10f)/[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Boltzmann_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.5,0.3); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1]);
            FitFunction->SetParNames("C_{B}","T_{Boltzm.} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("e") == 0 || type.CompareTo("E") == 0){
            cout <<Form("fitting %s with Exponential",FunctionName.Data()) << endl;
            TF1 *Exponential_Dummy = new TF1("Exponential_Dummy",Form("[0]/( 2 * TMath::Pi())/([1]*(%.10f+[1]))*exp(-(sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Exponential_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1.1,0.37); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1]);
            FitFunction->SetParNames("C_{E}","T_{exp} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("m") == 0 || type.CompareTo("M") == 0){
            cout <<Form("fitting %s with ModPowerlaw",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy = new TF1("ModPowerLaw_Dummy","[0]*pow((1 + (x)/[1]),-[2])",0,200);
            FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("A","p_{0}","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("doHag") == 0 || type.CompareTo("DOHag") == 0){
            cout <<Form("fitting %s with oHag and power law",FunctionName.Data()) << endl;
            TF1 *HagPow_Dummy = new TF1("HagPow_Dummy","(x<=[7])*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])+(x>[7])*[5]*1/pow(x,[6])",0,200);
            FitFunction = (TF1*)HagPow_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.7,0.07,1.,6.1,130.,5.,5.5); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7]);
            FitFunction->SetParNames("a","b","c","d","e","A","n","pt");
            //FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            FitFunction->SetParLimits(7,4.5,6.5);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("doHagPt") == 0 || type.CompareTo("DOHagPT") == 0){
            cout <<Form("fitting %s with oHag and power law * pT",FunctionName.Data()) << endl;
            TF1 *HagPow_Dummy = new TF1("HagPow_Dummy","(x<=[7])*x*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])+(x>[7])*x*[5]*1/pow(x,[6])",0,200);
            FitFunction = (TF1*)HagPow_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.7,0.07,1.,6.1,130.,5.,5.5); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7]);
            FitFunction->SetParNames("a","b","c","d","e","A","n","pt");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            //FitFunction->SetParLimits(7,4.5,6.5);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("oHag") == 0 || type.CompareTo("OHag") == 0){
            cout <<Form("fitting %s with ModHagedorn",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])",0,200);
            FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("A (mbGeV^{-2}c^{3})","a [(GeV/c)^{-1}]","b [(GeV/c)^{-1}]","p_{0} (GeV/c)","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("oHagPt") == 0 || type.CompareTo("OHagPt") == 0){
            cout <<Form("fitting %s with ModHagedorn times pt",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","x*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])",0,200);
            FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("A (mbGeV^{-2}c^{3})","a [(GeV/c)^{-1}]","b [(GeV/c)^{-1}]","p_{0} (GeV/c)","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("6pol") == 0 || type.CompareTo("6POL") == 0){
            cout <<Form("fitting %s with Polynom order 6",FunctionName.Data()) << endl;
            TF1 *Pol6_Dummy = new TF1("Pol6_Dummy","[0]+[1]*pow(x,2)+[2]*pow(x,4)+[3]*pow(x,6)",0,200);
            FitFunction = (TF1*)Pol6_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1.,1.,1.,1.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3]);
            FitFunction->SetParNames("a","b","c","d");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("rad") == 0 || type.CompareTo("RAD") == 0){
            cout <<Form("fitting %s with Radoslav Func",FunctionName.Data()) << endl;
            TF1 *Rad_Dummy = new TF1("Rad_Dummy","(x<=[3])*x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+0.135*([1]-2.)))*pow(1.+(sqrt(x*x+0.135*0.135)-0.135)/([1]*[2]), -[1]) + (x>[3]) * x *[0] * ([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+0.135*([1]-2.))) * pow([1]*[2]/([3]+[1]*[2]-0.135),[1]) * pow([3],[4]) * pow( 1./x, [4] )",0,200);
            FitFunction = (TF1*)Rad_Dummy->Clone(FunctionName);
            Double_t par0min3 = 10.;
            Double_t par0max3 = 1500.;
            Double_t par1min3 = 4.;
            Double_t par1max3 = 40.;
            Double_t par2min3 = 0.07;
            Double_t par2max3 = .7;
            Double_t par3min3 = 3.;
            Double_t par3max3 = 5.;
            Double_t par4min3 = 2.;
            Double_t par4max3 = 18.;
            if(Parameter != NULL){
                par0min3 = Parameter[0];
                par0max3 = Parameter[1];
                par1min3 = Parameter[2];
                par1max3 = Parameter[3];
                par2min3 = Parameter[4];
                par2max3 = Parameter[5];
                par3min3 = Parameter[6];
                par3max3 = Parameter[7];
                par4min3 = Parameter[8];
                par4max3 = Parameter[9];
            }
            FitFunction->SetParLimits(0, par0min3, par0max3);
            FitFunction->SetParLimits(1, par1min3, par1max3);
            FitFunction->SetParLimits(2, par2min3, par2max3);
            FitFunction->SetParLimits(3, par3min3, par3max3);
            FitFunction->SetParLimits(4, par4min3, par4max3);
            FitFunction->SetParNames("a","b","c","d");
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("tcm") == 0 || type.CompareTo("TCM") == 0){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
            cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
            TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass),0,200);
            FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
            if (mesonType.CompareTo("Pi0") == 0){
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            } else if (mesonType.CompareTo("Eta") == 0){
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            }
            FitFunction->SetParNames("Ae","Te","A","T","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("tcmpt") == 0 || type.CompareTo("TCMPT") == 0){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
            cout <<Form("fitting %s with two component model by Bylinkin by pT",FunctionName.Data()) << endl;
            TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form(" x*([0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4])) )",mass,mass,mass),0,200);
            FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
            if (mesonType.CompareTo("Pi0") == 0){
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            } else if (mesonType.CompareTo("Eta") == 0){
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); //TGraphAsymmErrors
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            }
            FitFunction->SetParNames("Ae","Te","A","T","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("kfunc") == 0 || type.CompareTo("KFUNC") == 0){
            cout <<Form("fitting %s with combined exp for flow and two Hagedorn functions ",FunctionName.Data()) << endl;
            TF1 *KFunc_Dummy = new TF1("KFunc_Dummy",Form("[0]*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2]) + [3]/(1+x/[4])^[5] + [6]/(1+x/[7])^[8]",mass,mass),0,200);
            FitFunction = (TF1*)KFunc_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1e3,0.6,0.1,1e3,1,8,0.,1,4); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7],Parameter[8]);
            FitFunction->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("modkfunc") == 0 || type.CompareTo("MODKFUNC") == 0){
            cout <<Form("fitting %s with combined exp for flow and two Hagedorn functions ",FunctionName.Data()) << endl;
            TF1 *modKFunc_Dummy = new TF1("modKFunc_Dummy",Form("[0]*(sqrt(x*x+%.10f*%.10f)-x*[1])/sqrt(1-[1]*[1])*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2]) + [3]/(1+x/[4])^[5] + [6]/(1+x/[7])^[8]",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)modKFunc_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1e3,0.6,0.1,1e3,1,8,0.,1,4); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7],Parameter[8]);
            FitFunction->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("modkfuncpt") == 0 || type.CompareTo("MODKFUNCPT") == 0){
            cout <<Form("fitting %s with combined exp for flow and two Hagedorn functions ",FunctionName.Data()) << endl;
            TF1 *modKFunc_Dummy = new TF1("modKFunc_Dummy",Form("x*([0]*(sqrt(x*x+%.10f*%.10f)-x*[1])/sqrt(1-[1]*[1])*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2]) + [3]/(1+x/[4])^[5] + [6]/(1+x/[7])^[8])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)modKFunc_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1e3,0.6,0.1,1e3,1,8,0.,1,4); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7],Parameter[8]);
            FitFunction->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if (type.CompareTo("blastWave") == 0){
            cout <<Form("fitting %s with blast wave function",FunctionName.Data()) << endl;
            TF1 *blastWave_Dummy = new TF1("blastWave_Dummy", BGBlastWave_Func, 0., 10., 5);
            FitFunction = (TF1*)blastWave_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(mass, 0.9, 0.1, 1, 1.e6); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("mass", "beta_max", "T", "n", "norm");
            FitFunction->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        //ATTENTION: mT scaling will not work in AliPhysics task when using TSpline!! in that case, need to mT scale in cocktail framework and provide all needed spectra
        if (type.CompareTo("spline") == 0){
            cout <<Form("fitting %s with spline function",FunctionName.Data()) << endl;
            TGraphAsymmErrors* tempGraph = new TGraphAsymmErrors(Obj);
            TSpline3* spline = new TSpline3(Form("spline_%s",FunctionName.Data()), tempGraph,"b2e2",0,0);
            vecGraphs.push_back(spline);
            FitFunction = new TF1(Form("spline_TF1_%s",FunctionName.Data()), methodTSplineTF1, xmin, xmax, tempGraph->GetN()); // npars = 2*nodes+2
            FitFunction->FixParameter(0,counterGraphs+0.1);
            FitFunction->SetNpx(1000);
            counterGraphs++;
        }
        //ATTENTION: mT scaling will not work in AliPhysics task when using TSpline!! in that case, need to mT scale in cocktail framework and provide all needed spectra
        if (type.CompareTo("spline_p") == 0){
            cout <<Form("fitting %s with spline_p function",FunctionName.Data()) << endl;
            TGraphAsymmErrors* tempGraph = new TGraphAsymmErrors(Obj);
            TSpline3* spline = new TSpline3(Form("spline_p_%s",FunctionName.Data()), tempGraph,"b2e2",0,0);
            vecGraphs.push_back(spline);
            FitFunction = new TF1(FunctionName, methodTSplinePowerLawTF1, xmin, xmax, tempGraph->GetN()); // npars = 2*nodes+2
            FitFunction->SetNpx(1000);
            FitFunction->FixParameter(0,counterGraphs+0.1);
            counterGraphs++;
            if(Parameter == NULL){
              FitFunction->SetParameter(1,10.);
              FitFunction->SetParameter(2,2.);
              FitFunction->SetParameter(3,5.);
            }else{
              FitFunction->SetParameter(1,Parameter[0]);
              FitFunction->SetParLimits(1,ParameterLimitsDown[0],ParameterLimitsUp[0]);
              FitFunction->SetParameter(2,Parameter[1]);
              FitFunction->SetParLimits(2,ParameterLimitsDown[1],ParameterLimitsUp[1]);
              FitFunction->SetParameter(3,Parameter[2]);
              FitFunction->SetParLimits(3,ParameterLimitsDown[2],ParameterLimitsUp[2]);
            }
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
    }

    if(ClassName.BeginsWith("TGraph")){
        TGraphAsymmErrors *Obj = (TGraphAsymmErrors*)Obj_Dummy;
        if(type.CompareTo("h") == 0 || type.CompareTo("H") == 0){
            cout <<Form("fitting %s with Hagedorn",FunctionName.Data()) << endl;
            TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*pow([2]/([2]+x),[1])",0,200);
            FitFunction = (TF1*)Hagedorn_Dummy->Clone(FunctionName);
            if(Parameter == NULL) FitFunction->SetParameters(19.,6.8,0.84); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("C_{H}","n","p_{0} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("powPure") == 0 || type.CompareTo("PowPure") == 0){
            cout <<Form("fitting %s with Pure Powerlaw",FunctionName.Data()) << endl;
            TF1 *PowerLaw_Dummy = new TF1("PowerLawPure_Dummy","[0] * 1/pow(x,[1])",0,200);
            FitFunction = (TF1*)PowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1]);
            FitFunction->SetParNames("A_{pow}","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("p") == 0 || type.CompareTo("P") == 0){
            cout <<Form("fitting %s with Powerlaw",FunctionName.Data()) << endl;
            TF1 *PowerLaw_Dummy = new TF1("PowerLaw_Dummy","[0] * 2 / TMath::Pi() * ([1]-1.)*([1]-2.)/pow([1]-3.,2) /x * pow(1+2*x/[2]/([1]-3),-[1])",0,200);
            FitFunction = (TF1*)PowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.37); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("A_{pow}","n","p_{0}");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("l")  == 0|| type.CompareTo("L") == 0){
            cout <<Form("fitting %s with Levy",FunctionName.Data()) << endl;
            TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("tmpt") == 0 || type.CompareTo("TMPT") == 0){
            // Tsallis Dummy multiplied with pt
            TF1 *Levy_Dummy = new TF1("Tsallis_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * x* pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("qcd") == 0 || type.CompareTo("QCD") == 0){
            cout <<Form("fitting %s with QCD",FunctionName.Data()) << endl;
            TF1 *QCD_Dummy = new TF1("QCD_Dummy","[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))",0,200);
            FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("a","b","c","d","e");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("doHag") == 0 || type.CompareTo("DOHag") == 0){
            cout <<Form("fitting %s with oHag and power law",FunctionName.Data()) << endl;
            TF1 *HagPow_Dummy = new TF1("HagPow_Dummy","(x<=[7])*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])+(x>[7])*[5]*1/pow(x,[6])",0,200);
            FitFunction = (TF1*)HagPow_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.7,0.07,1.,6.1,130.,5.,5.5); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7]);
            FitFunction->SetParNames("a","b","c","d","e","A","n","pt");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("doHagPt") == 0 || type.CompareTo("DOHagPT") == 0){
            cout <<Form("fitting %s with oHag and power law * pT",FunctionName.Data()) << endl;
            TF1 *HagPow_Dummy = new TF1("HagPow_Dummy","(x<=[7])*x*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])+(x>[7])*x*[5]*1/pow(x,[6])",0,200);
            FitFunction = (TF1*)HagPow_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.7,0.07,1.,6.1,130.,5.,5.5); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7]);
            FitFunction->SetParNames("a","b","c","d","e","A","n","pt");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            //FitFunction->SetParLimits(7,4.5,6.5);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("oHag") == 0 || type.CompareTo("OHag") == 0){
            cout <<Form("fitting %s with ModHagedorn",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])",0,200);
            FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("A (mbGeV^{-2}c^{3})","a [(GeV/c)^{-1}]","b [(GeV/c)^{-1}]","p_{0} (GeV/c)","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("oHagPt") == 0 || type.CompareTo("OHagPt") == 0){
            cout <<Form("fitting %s with ModHagedorn times pt",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","x*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])",0,200);
            FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("A (mbGeV^{-2}c^{3})","a [(GeV/c)^{-1}]","b [(GeV/c)^{-1}]","p_{0} (GeV/c)","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("b") == 0|| type.CompareTo("B") == 0){
            cout <<Form("fitting %s with Boltzmann",FunctionName.Data()) << endl;
            TF1 *Boltzmann_Dummy =  new TF1("Boltzmann_Dummy",Form("[0] *sqrt(x*x+%.10f*%.10f)* exp(- sqrt(x*x+%.10f*%.10f)/[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Boltzmann_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.5,0.3); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1]);
            FitFunction->SetParNames("C_{B}","T_{Boltzm.} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("e") == 0 || type.CompareTo("E") == 0){
            cout <<Form("fitting %s with Exponential",FunctionName.Data()) << endl;
            TF1 *Exponential_Dummy = new TF1("Exponential_Dummy",Form("[0]/( 2 * TMath::Pi())/([1]*(%.10f+[1]))*exp(-(sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Exponential_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1.1,0.37); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1]);
            FitFunction->SetParNames("C_{E}","T_{exp} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("m") == 0 || type.CompareTo("M") == 0){
            cout <<Form("fitting %s with ModPowerlaw",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy = new TF1("ModPowerLaw_Dummy","[0]*pow((1 + x/[1]),-[2])",0,200);
            FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("A","p_{0}","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("tcm") == 0 || type.CompareTo("TCM") == 0){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
            cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
            TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass),0,200);
            FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
            if (mesonType.CompareTo("Pi0") == 0){
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            } else if (mesonType.CompareTo("Eta") == 0){
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); //TGraphAsymmErrors
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            }
            FitFunction->SetParNames("Ae","Te","A","T","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("tcmpt") == 0 || type.CompareTo("TCMPT") == 0){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
            cout <<Form("fitting %s with two component model by Bylinkin by pT",FunctionName.Data()) << endl;
            TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form(" x*([0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4])) )",mass,mass,mass),0,200);
            FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
            if (mesonType.CompareTo("Pi0") == 0){
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            } else if (mesonType.CompareTo("Eta") == 0){
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); //TGraphAsymmErrors
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            }
            FitFunction->SetParNames("Ae","Te","A","T","n");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("kfunc") == 0 || type.CompareTo("KFUNC") == 0){
            cout <<Form("fitting %s with combined exp for flow and two Hagedorn functions ",FunctionName.Data()) << endl;
            TF1 *KFunc_Dummy = new TF1("KFunc_Dummy",Form("[0]*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2]) + [3]/(1+x/[4])^[5] + [6]/(1+x/[7])^[8]",mass,mass),0,200);
            FitFunction = (TF1*)KFunc_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1e3,0.6,0.1,1e3,1,8,0.,1,4); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7],Parameter[8]);
            FitFunction->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("modkfunc") == 0 || type.CompareTo("MODKFUNC") == 0){
            cout <<Form("fitting %s with combined exp for flow and two Hagedorn functions ",FunctionName.Data()) << endl;
            TF1 *modKFunc_Dummy = new TF1("modKFunc_Dummy",Form("[0]*(sqrt(x*x+%.10f*%.10f)-x*[1])/sqrt(1-[1]*[1])*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2]) + [3]/(1+x/[4])^[5] + [6]/(1+x/[7])^[8]",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)modKFunc_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1e3,0.6,0.1,1e3,1,8,0.,1,4); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7],Parameter[8]);
            FitFunction->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("modkfuncpt") == 0 || type.CompareTo("MODKFUNCPT") == 0){
            cout <<Form("fitting %s with combined exp for flow and two Hagedorn functions ",FunctionName.Data()) << endl;
            TF1 *modKFunc_Dummy = new TF1("modKFunc_Dummy",Form("x*([0]*(sqrt(x*x+%.10f*%.10f)-x*[1])/sqrt(1-[1]*[1])*exp((x*[1] - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2]) + [3]/(1+x/[4])^[5] + [6]/(1+x/[7])^[8])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)modKFunc_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1e3,0.6,0.1,1e3,1,8,0.,1,4); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7],Parameter[8]);
            FitFunction->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("rad")  == 0|| type.CompareTo("RAD") == 0){
            cout <<Form("fitting %s with Radoslav Func",FunctionName.Data()) << endl;
            TF1 *Rad_Dummy = new TF1("Rad_Dummy","(x<=[3])*x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+0.135*([1]-2.)))*pow(1.+(sqrt(x*x+0.135*0.135)-0.135)/([1]*[2]), -[1]) + (x>[3]) * x *[0] * ([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+0.135*([1]-2.))) * pow([1]*[2]/([3]+[1]*[2]-0.135),[1]) * pow([3],[4]) * pow( 1./x, [4] )",0,200);
            FitFunction = (TF1*)Rad_Dummy->Clone(FunctionName);
            Double_t par0min3 = 150.;
            Double_t par0max3 = 1500.;
            Double_t par1min3 = 4.;
            Double_t par1max3 = 40.;
            Double_t par2min3 = 0.07;
            Double_t par2max3 = .7;
            Double_t par3min3 = 4.36;
            Double_t par3max3 = 4.40;
            Double_t par4min3 = 2.;
            Double_t par4max3 = 18.;
            if(Parameter != NULL){
                par0min3 = Parameter[0];
                par0max3 = Parameter[1];
                par1min3 = Parameter[2];
                par1max3 = Parameter[3];
                par2min3 = Parameter[4];
                par2max3 = Parameter[5];
                par3min3 = Parameter[6];
                par3max3 = Parameter[7];
                par4min3 = Parameter[8];
                par4max3 = Parameter[9];
            }
            FitFunction->SetParLimits(0, par0min3, par0max3);
            FitFunction->SetParLimits(1, par1min3, par1max3);
            FitFunction->SetParLimits(2, par2min3, par2max3);
            FitFunction->SetParLimits(3, par3min3, par3max3);
            FitFunction->SetParLimits(4, par4min3, par4max3);
            FitFunction->SetParNames("a","b","c","d");
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            return FitFunction;
        }
        if(type.CompareTo("ptredt") == 0 || type.CompareTo("ptredt") == 0){
            cout <<Form("fitting %s with Tsallis (*Pt)",FunctionName.Data()) << endl;
            TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass),0,200);
            FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
            FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if (type.CompareTo("blastWave") == 0){
            cout <<Form("fitting %s with blast wave function",FunctionName.Data()) << endl;
            TF1 *blastWave_Dummy = new TF1("blastWave_Dummy", BGBlastWave_Func, 0., 10., 5);
            FitFunction = (TF1*)blastWave_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(mass, 0.9, 0.1, 1, 1.e6); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("mass", "beta_max", "T", "n", "norm");
            FitFunction->SetParNames("N_{0}","v_{flow}","M","T_{kin}","N_{1}","p01","pow1","N_{2}","p02","pow2");
            FitFunction = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        //ATTENTION: mT scaling will not work in AliPhysics task when using TSpline!! in that case, need to mT scale in cocktail framework and provide all needed spectra
        if (type.CompareTo("spline") == 0){
            cout <<Form("fitting %s with spline function",FunctionName.Data()) << endl;
            TSpline3* spline = new TSpline3(Form("spline_%s",FunctionName.Data()), Obj,"b2e2",0,0);
            vecGraphs.push_back(spline);
            FitFunction = new TF1(Form("spline_TF1_%s",FunctionName.Data()), methodTSplineTF1, xmin, xmax, Obj->GetN()); // npars = 2*nodes+2
            FitFunction->FixParameter(0,counterGraphs+0.1);
            FitFunction->SetNpx(1000);
            counterGraphs++;
        }
        //ATTENTION: mT scaling will not work in AliPhysics task when using TSpline!! in that case, need to mT scale in cocktail framework and provide all needed spectra
        if (type.CompareTo("spline_p") == 0){
            cout <<Form("fitting %s with spline_p function",FunctionName.Data()) << endl;
            TSpline3* spline = new TSpline3(Form("spline_p_%s",FunctionName.Data()), Obj,"b2e2",0,0);
            vecGraphs.push_back(spline);
            FitFunction = new TF1(FunctionName, methodTSplinePowerLawTF1, xmin, xmax, Obj->GetN()); // npars = 2*nodes+2
            FitFunction->SetNpx(1000);
            FitFunction->FixParameter(0,counterGraphs+0.1);
            counterGraphs++;
            if(Parameter == NULL){
              FitFunction->SetParameter(1,10.);
              FitFunction->SetParameter(2,2.);
              FitFunction->SetParameter(3,5.);
            }else{
              FitFunction->SetParameter(1,Parameter[0]);
              FitFunction->SetParLimits(1,ParameterLimitsDown[0],ParameterLimitsUp[0]);
              FitFunction->SetParameter(2,Parameter[1]);
              FitFunction->SetParLimits(2,ParameterLimitsDown[1],ParameterLimitsUp[1]);
              FitFunction->SetParameter(3,Parameter[2]);
              FitFunction->SetParLimits(3,ParameterLimitsDown[2],ParameterLimitsUp[2]);
            }
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }

    }

    if (Parameter)              delete Parameter;
    if (ParameterLimitsDown)    delete ParameterLimitsDown;
    if (ParameterLimitsUp)      delete ParameterLimitsUp;
    if (convertToHist && Obj_Dummy)              delete Obj_Dummy;

    return FitFunction;
}


//================================================================================================================
//Function to fit the ratio of two particle spectra
//================================================================================================================
TF1* FitRatio(  TObject *Obj_Dummy      = NULL,
                TString FunctionName    = "",
                TString mesonTypeNum    = "Pi0",
                TString mesonTypeDenom  = "",
                TString centrality      = "",
                TString method          = "",
                TString type            = "",
                Bool_t convertToHist    = kFALSE
              ) {
    TString ClassName                                                   = Obj_Dummy->ClassName();
    if(convertToHist){
        if (ClassName.Contains("TH1"))
            cout << "is histo already" << endl;
        else if (ClassName.CompareTo("TGraphErrors") == 0)
            Obj_Dummy = GraphToHist_withErrors((TGraphErrors*)Obj_Dummy, "Obj_Dummy_Hist");
        else if (ClassName.CompareTo("TGraphAsymmErrors") == 0)
            Obj_Dummy = GraphToHist_withErrors((TGraphAsymmErrors*)Obj_Dummy, "Obj_Dummy_Hist");
        ClassName = Obj_Dummy->ClassName();
    }

    // xRange
    Double_t xmin                                                       = initRatioXRange[GetParticleIterator(mesonTypeNum)][GetParticleIterator(mesonTypeDenom)][GetCentralityIterator(centrality)][GetMethodIterator(method)][0];
    Double_t xmax                                                       = initRatioXRange[GetParticleIterator(mesonTypeNum)][GetParticleIterator(mesonTypeDenom)][GetCentralityIterator(centrality)][GetMethodIterator(method)][1];
    if (xmin == -9999) {
        if (ClassName.Contains("TH1"))                          xmin    = ((TH1D*)Obj_Dummy)->GetXaxis()->GetXmin();
        else if (ClassName.CompareTo("TGraphErrors") == 0)      xmin    = GetXRangeFromGraph((TGraphErrors*)Obj_Dummy, kFALSE);
        else if (ClassName.CompareTo("TGraphAsymmErrors") == 0) xmin    = GetXRangeFromGraph((TGraphAsymmErrors*)Obj_Dummy, kFALSE);
    }
    if (xmax == -9999) {
        if (ClassName.Contains("TH1"))                          xmax    = ((TH1D*)Obj_Dummy)->GetXaxis()->GetXmax();
        else if (ClassName.CompareTo("TGraphErrors") == 0)      xmax    = GetXRangeFromGraph((TGraphErrors*)Obj_Dummy, kTRUE);
        else if (ClassName.CompareTo("TGraphAsymmErrors") == 0) xmax    = GetXRangeFromGraph((TGraphAsymmErrors*)Obj_Dummy, kTRUE);
    }

    // fit option
    TString FitOptions                                                  = initRatioFitOption[GetParticleIterator(mesonTypeNum)][GetParticleIterator(mesonTypeDenom)][GetCentralityIterator(centrality)][GetMethodIterator(method)];

    // start parameter
    Double_t *Parameter                                                 = NULL;
    Int_t nPars                                                         = 0;
    for (Int_t i=0; i<30; i++) {
        if (initRatioParameter[GetParticleIterator(mesonTypeNum)][GetParticleIterator(mesonTypeDenom)][GetCentralityIterator(centrality)][GetMethodIterator(method)][i] != -9999) nPars++;
    }
    if (nPars>0) {
        Parameter                                                       = new Double_t[nPars];
        for (Int_t i=0; i<nPars; i++)
            Parameter[i]                                                = initRatioParameter[GetParticleIterator(mesonTypeNum)][GetParticleIterator(mesonTypeDenom)][GetCentralityIterator(centrality)][GetMethodIterator(method)][i];
    }

    // start parameter limits (if -9999, don't use)
    Double_t *ParameterLimitsDown                                       = NULL;
    Double_t *ParameterLimitsUp                                         = NULL;
    if (nPars) {
        ParameterLimitsDown                                             = new Double_t[nPars];
        ParameterLimitsUp                                               = new Double_t[nPars];
        for (Int_t i=0; i<nPars; i++) {
            ParameterLimitsDown[i]                                      = initRatioParameterLimit[GetParticleIterator(mesonTypeNum)][GetParticleIterator(mesonTypeDenom)][GetCentralityIterator(centrality)][GetMethodIterator(method)][i][0];
            ParameterLimitsUp[i]                                        = initRatioParameterLimit[GetParticleIterator(mesonTypeNum)][GetParticleIterator(mesonTypeDenom)][GetCentralityIterator(centrality)][GetMethodIterator(method)][i][1];
        }
    }

    TString mesonType[2]    = {mesonTypeNum, mesonTypeDenom};
    Double_t mass[2]        = {0, 0};

    for (Int_t i = 0; i < 2; i++) {
        if (mesonType[i].CompareTo(fParticle[0]) == 0){
            // Pi0
            mass[i] = TDatabasePDG::Instance()->GetParticle(111)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[1]) == 0){
            // Eta
            mass[i] = TDatabasePDG::Instance()->GetParticle(221)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[2]) == 0){
            // Omega
            mass[i] = TDatabasePDG::Instance()->GetParticle(223)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[3]) == 0){
            // EtaPrima
            mass[i] = TDatabasePDG::Instance()->GetParticle(331)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[4]) == 0){
            // GammaDir
            mass[i] = 0;
        } else if (mesonType[i].CompareTo(fParticle[5]) == 0){
            // CPion
            mass[i] = TDatabasePDG::Instance()->GetParticle(211)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[6]) == 0){
            // CKaon
            mass[i] = TDatabasePDG::Instance()->GetParticle(321)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[7]) == 0){
            // Proton
            mass[i] = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[8]) == 0){
            // CHadron
            mass[i] = 0;
        } else if (mesonType[i].CompareTo(fParticle[9]) == 0){
            // Phi
            mass[i] = TDatabasePDG::Instance()->GetParticle(333)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[10]) == 0){
            // NKaonStar
            mass[i] = TDatabasePDG::Instance()->GetParticle(313)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[11]) == 0){
            // NRho
            mass[i] = TDatabasePDG::Instance()->GetParticle(113)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[12]) == 0){
            // CRho
            mass[i] = TDatabasePDG::Instance()->GetParticle(213)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[13]) == 0){
            // NDelta
            mass[i] = TDatabasePDG::Instance()->GetParticle(2114)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[14]) == 0){
            // CDelta
            mass[i] = TDatabasePDG::Instance()->GetParticle(2214)->Mass();     // Delta^+
        } else if (mesonType[i].CompareTo(fParticle[15]) == 0){
            // NKaonSubS
            mass[i] = TDatabasePDG::Instance()->GetParticle(310)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[16]) == 0){
            // Lambda
            mass[i] = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[17]) == 0){
            // NSigma
            mass[i] = TDatabasePDG::Instance()->GetParticle(3212)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[18]) == 0){
            // CSigma
            mass[i] = TDatabasePDG::Instance()->GetParticle(3222)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[19]) == 0){
            // COmega
            mass[i] = TDatabasePDG::Instance()->GetParticle(3334)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[20]) == 0){
            // CXi
            mass[i] = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
        } else if (mesonType[i].CompareTo(fParticle[21]) == 0){
            // JPsi
            mass[i] = TDatabasePDG::Instance()->GetParticle(443)->Mass();
        } else {
            mass[i] = 0;
        }
        cout << "Meson " << i << " mass: "<< mass[i] << endl;
    }

    if(Obj_Dummy == NULL) {
        cout << "Obj_Dummy is NULL, not implemented yet, skipping" << endl;
        return NULL;
    }

    TF1 *FitFunction                = new TF1();

    if(ClassName.BeginsWith("TH1")){
        TH1D *Obj                   = (TH1D*)Obj_Dummy;
        if(type.CompareTo("const") == 0 || type.CompareTo("c") == 0){
            cout <<Form("fitting %s with const",FunctionName.Data()) << endl;
            TF1 *const_Dummy        = new TF1("const_Dummy","[0]",0,200);
            FitFunction             = (TF1*)const_Dummy->Clone(FunctionName);
            FitFunction->SetRange(xmin, xmax);
            if(Parameter == NULL)FitFunction->SetParameter(0, 1); // standard parameter optimize if necessary
            else FitFunction->SetParameter(0, Parameter[0]);
            FitFunction->SetParNames("ratio");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
        }
        if (type.CompareTo("pol1") == 0 || type.CompareTo("POL1") == 0){
            cout << Form("fitting %s with linear ratio", FunctionName.Data()) << endl;
            TF1 *lin_dummy          = new TF1("lin_dummy", "pol1",0,200);
            FitFunction             = (TF1*)lin_dummy->Clone(FunctionName);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if (type.CompareTo("qcdRatio") == 0 || type.CompareTo("QCDRATIO") == 0){
            cout << Form("fitting %s with qcd ratio", FunctionName.Data()) << endl;
            TF1 *qcd_dummy          = new TF1("qcd_dummy", "[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))/[5]/TMath::Power(x,-1*([6]+[7]/(TMath::Power(x,[8])+[9])))",0,200);
            FitFunction             = (TF1*)qcd_dummy->Clone(FunctionName);
            if (Parameter == NULL) FitFunction->SetParameters(24, 6.7, -6.5, 1., 10, 24, 6.7, -6.5, 1., 10);
            else FitFunction->SetParameters(Parameter[0], Parameter[1], Parameter[2], Parameter[3], Parameter[4], Parameter[5], Parameter[6], Parameter[7], Parameter[8], Parameter[9]);
            FitFunction->SetParNames("a_num", "b_num", "c_num", "d_num", "e_num", "a_denom", "b_denom", "c_denom", "d_denom", "e_denom");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("oHagRatio") == 0 || type.CompareTo("OHagRATIO") == 0){
            cout << Form("fitting %s with ModHagedorn ratio",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy  = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])/([5]*pow(exp(-[6]*x-abs([7])*x*x)+x/[8],-[9]))",0,200);
            FitFunction             = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.,450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4], Parameter[5], Parameter[6], Parameter[7], Parameter[8], Parameter[9]);
            FitFunction->SetParNames("A_num","a_num","b_num","p_{0}_num","n_num","A_denom","a_denom","b_denom","p_{0}_denom","n_denom");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("lTooHag") == 0){
            cout << Form("fitting %s with levy to mod. hagendorn ratio",FunctionName.Data()) << endl;
            TF1 *lTooHag_Dummy      = new TF1("lTooHag_Dummy",Form("([0]/(2*TMath::Pi())*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.)))*pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]),-[1]))/([3]*pow(exp(-[4]*x-abs([5])*x*x)+x/[6],-[7]))",mass[0],mass[0],mass[0],mass[0]),0,200);
            FitFunction             = (TF1*)lTooHag_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18,450.,0.37,0.2,0.7,8.);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4], Parameter[5], Parameter[6], Parameter[7]);
            FitFunction->SetParNames("(dN/dy)_num","n_num","(T_{Levy})_num","A_denom","a_denom","b_denom","p_{0}_denom","n_denom");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("mTooHag") == 0){
            cout << Form("fitting %s with mod. power law to mod. hagendorn ratio",FunctionName.Data()) << endl;
            TF1 *mTooHag_Dummy      = new TF1("mTooHag_Dummy","([0]*pow((1 + (x)/[1]),-[2]))/([3]*pow(exp(-[4]*x-abs([5])*x*x)+x/[6],-[7]))",0,200);
            FitFunction             = (TF1*)mTooHag_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.,450.,0.37,0.2,0.7,8.);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4], Parameter[5], Parameter[6], Parameter[7]);
            FitFunction->SetParNames("A_num","(p_{0})_num","n_num","A_denom","a_denom","b_denom","p_{0}_denom","n_denom");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("blastWave") == 0){
            // 0: normalization constant
            // 1: flow velocity beta
            // 2: kinetic freeze-out temperature
            // 3: constant particle ratio at high pt
            // 4: pt at which hard scattering becomes important
            // 5: width of the transition
            cout << Form("fitting %s with blastwave inspired fit",FunctionName.Data()) << endl;
            TString sBWRatio        = Form("[0]*exp((sqrt(x*x+%.10f*%.10f) - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2])", mass[1], mass[1], mass[0], mass[0]);
            TString sTrans          = "1./(1. + exp((x - [4])/[5]))";
            TString sComb           = sTrans + "*" + sBWRatio + " + (1 - " + sTrans + ") * [3]";
            TF1* blastWave_Dummy    = new TF1("blastWave_Dummy", sComb.Data(),0,200);
            FitFunction             = (TF1*)blastWave_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(0.5,0.05,0.1,100,3,0.5);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4], Parameter[5]);
            FitFunction->SetParNames("c_{norm}","beta_{flow}","T_{kin. freeze out}","N_{num}/N_{denom}","pt_{hard}","w_{trans}");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("softHard") == 0){
            // 0: normalization constant
            // 1: flow velocity beta
            // 2: kinetic freeze-out temperature
            // 3: constant particle ratio at high pt
            // 4: relative normalization soft/hard
            // 5: Parameter p0 of Hagedorn-type function for hard scattering part
            // 6: Power n of Hagedorn-type function for hard scattering part
            cout << Form("fitting %s with soft+hard ratio",FunctionName.Data()) << endl;
            TString sCombNum1       = Form("([0]*exp(([1]*x - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2])", mass[0], mass[0]);
            TString sCombNum2       = "+ [3]*[4]*TMath::Power(1 + x/[5]*x/[5],-[6])) / ";
            TString sCombDen1       = Form("(exp(([1]*x - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2])", mass[1], mass[1]);
            TString sCombDen2       = "+ [4]*TMath::Power(1 + x/[5]*x/[5],-[6]))";
            TString sComb           = sCombNum1 + sCombNum2 + sCombDen1 + sCombDen2;
            TF1* softHard_Dummy     = new TF1("softHard_Dummy", sComb.Data(),0,200);
            FitFunction             = (TF1*)softHard_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(0.5,0.05,0.1,100,0.5,2,5);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6]);
            FitFunction->SetParNames("c_{norm}","beta_{flow}","T_{kin. freeze out}","N_{num}/N_{denom}","c_{soft to hard}","pt_{0}","n");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("empirical") == 0){
            cout << Form("fitting %s with empirical ratio",FunctionName.Data()) << endl;
            TF1* empirical_Dummy    = new TF1("empirical_Dummy", "[0]/(1+exp(-(x-[1])/[2])) + [3]*Gaus(x,[4],[5])",0,200);
            FitFunction             = (TF1*)empirical_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1,1,0.5,0,2,1);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5]);
            FitFunction->SetParNames("c_{norm}","a","b","gaus_{amp}","gaus_{mu}","gaus_{sigma}");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
    }

    if (ClassName.BeginsWith("TGraph")) {
        TGraphErrors *Obj           = (TGraphErrors*)Obj_Dummy;
        if (type.CompareTo("pol1") == 0 || type.CompareTo("POL1") == 0){
            cout << Form("fitting %s with linear ratio", FunctionName.Data()) << endl;
            TF1 *lin_dummy          = new TF1("lin_dummy", "pol1",0,200);
            FitFunction             = (TF1*)lin_dummy->Clone(FunctionName);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if (type.CompareTo("qcdRatio") == 0 || type.CompareTo("QCDRATIO") == 0){
            cout << Form("fitting %s with qcd ratio", FunctionName.Data()) << endl;
            TF1* qcd_dummy          = new TF1("qcd_dummy", "[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))/[5]/TMath::Power(x,-1*([6]+[7]/(TMath::Power(x,[8])+[9])))",0,200);
            FitFunction             = (TF1*)qcd_dummy->Clone(FunctionName);
            if (Parameter == NULL) FitFunction->SetParameters(24, 6.7, -6.5, 1., 10, 24, 6.7, -6.5, 1., 10);
            else FitFunction->SetParameters(Parameter[0], Parameter[1], Parameter[2], Parameter[3], Parameter[4], Parameter[5], Parameter[6], Parameter[7], Parameter[8], Parameter[9]);
            FitFunction->SetParNames("a_num", "b_num", "c_num", "d_num", "e_num", "a_denom", "b_denom", "c_denom", "d_denom", "e_denom");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("oHagRatio") == 0 || type.CompareTo("OHagRATIO") == 0){
            cout <<Form("fitting %s with ModHagedorn ratio",FunctionName.Data()) << endl;
            TF1 *ModPowerLaw_Dummy  = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])/([5]*pow(exp(-[6]*x-abs([7])*x*x)+x/[8],-[9]))",0,200);
            FitFunction             = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.,450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4], Parameter[5], Parameter[6], Parameter[7], Parameter[8], Parameter[9]);
            FitFunction->SetParNames("A_num","a_num","b_num","p_{0}_num","n_num","A_denom","a_denom","b_denom","p_{0}_denom","n_denom");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("lTooHag") == 0){
            cout << Form("fitting %s with levy to mod. hagendorn ratio",FunctionName.Data()) << endl;
            TF1 *lTooHag_Dummy      = new TF1("lTooHag_Dummy",Form("([0]/(2*TMath::Pi())*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.)))*pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]),-[1]))/([3]*pow(exp(-[4]*x-abs([5])*x*x)+x/[6],-[7]))",mass[0],mass[0],mass[0],mass[0]),0,200);
            FitFunction             = (TF1*)lTooHag_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18,450.,0.37,0.2,0.7,8.);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4], Parameter[5], Parameter[6], Parameter[7]);
            FitFunction->SetParNames("(dN/dy)_num","n_num","(T_{Levy})_num","A_denom","a_denom","b_denom","p_{0}_denom","n_denom");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("mTooHag") == 0){
            cout << Form("fitting %s with mod. power law to mod. hagendorn ratio",FunctionName.Data()) << endl;
            TF1 *mTooHag_Dummy      = new TF1("mTooHag_Dummy","([0]*pow((1 + (x)/[1]),-[2]))/([3]*pow(exp(-[4]*x-abs([5])*x*x)+x/[6],-[7]))",0,200);
            FitFunction             = (TF1*)mTooHag_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.,450.,0.37,0.2,0.7,8.);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4], Parameter[5], Parameter[6], Parameter[7]);
            FitFunction->SetParNames("A_num","(p_{0})_num","n_num","A_denom","a_denom","b_denom","p_{0}_denom","n_denom");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("blastWave") == 0){
            // 0: normalization constant
            // 1: flow velocity beta
            // 2: kinetic freeze-out temperature
            // 3: constant particle ratio at high pt
            // 4: pt at which hard scattering becomes important
            // 5: width of the transition
            cout << Form("fitting %s with blastwave inspired fit",FunctionName.Data()) << endl;
            TString sBWRatio        = Form("[0]*exp((sqrt(x*x+%.10f*%.10f) - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2])", mass[1], mass[1], mass[0], mass[0]);
            TString sTrans          = "1./(1. + exp((x - [4])/[5]))";
            TString sComb           = sTrans + "*" + sBWRatio + " + (1 - " + sTrans + ") * [3]";
            TF1* blastWave_Dummy    = new TF1("blastWave_Dummy", sComb.Data(),0,200);
            FitFunction             = (TF1*)blastWave_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(0.5,0.05,0.1,100,3,0.5);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4], Parameter[5]);
            FitFunction->SetParNames("c_{norm}","beta_{flow}","T_{kin. freeze out}","N_{num}/N_{denom}","pt_{hard}","w_{trans}");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("softHard") == 0){
            // 0: normalization constant
            // 1: flow velocity beta
            // 2: kinetic freeze-out temperature
            // 3: constant particle ratio at high pt
            // 4: relative normalization soft/hard
            // 5: Parameter p0 of Hagedorn-type function for hard scattering part
            // 6: Power n of Hagedorn-type function for hard scattering part
            cout << Form("fitting %s with soft+hard ratio",FunctionName.Data()) << endl;
            TString sCombNum1       = Form("([0]*exp(([1]*x - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2])", mass[0], mass[0]);
            TString sCombNum2       = "+ [3]*[4]*TMath::Power(1 + x/[5]*x/[5],-[6])) / ";
            TString sCombDen1       = Form("(exp(([1]*x - sqrt(x*x+%.10f*%.10f))/sqrt(1-[1]*[1])/[2])", mass[1], mass[1]);
            TString sCombDen2       = "+ [4]*TMath::Power(1 + x/[5]*x/[5],-[6]))";
            TString sComb           = sCombNum1 + sCombNum2 + sCombDen1 + sCombDen2;
            TF1* softHard_Dummy     = new TF1("softHard_Dummy", sComb.Data(),0,200);
            FitFunction             = (TF1*)softHard_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(0.5,0.05,0.1,100,0.5,2,5);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6]);
            FitFunction->SetParNames("c_{norm}","beta_{flow}","T_{kin. freeze out}","N_{num}/N_{denom}","c_{soft to hard}","pt_{0}","n");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
        if(type.CompareTo("empirical") == 0){
            cout << Form("fitting %s with empirical ratio",FunctionName.Data()) << endl;
            TF1* empirical_Dummy    = new TF1("empirical_Dummy", "[0]/(1+exp(-(x-[1])/[2])) + gaus(x,[3],[4])",0,200);
            FitFunction             = (TF1*)empirical_Dummy->Clone(FunctionName);
            if(Parameter == NULL)FitFunction->SetParameters(1,1,0.5,2,1);
            else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
            FitFunction->SetParNames("c_{norm}","a","b","gaus_mu","gaus_sigma");
            FitFunction             = SetParameterLimits(FitFunction, ParameterLimitsDown, ParameterLimitsUp);
            Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
        }
    }

    if (Parameter)              delete Parameter;
    if (ParameterLimitsDown)    delete ParameterLimitsDown;
    if (ParameterLimitsUp)      delete ParameterLimitsUp;
    if (convertToHist && Obj_Dummy)              delete Obj_Dummy;

    return FitFunction;
}


void PlotCocktailParametrizationListFullRange(  TList* list,
                                                TList* inputlist,
                                                TString collSys,
                                                TString energy,
                                                TString cent,
                                                TString suffix,
                                                TString paramDirName  );
//================================================================================================================
//Write parametrizations to file
//================================================================================================================
void WriteParametrizationsFile(TList* list, TList* inputList, TString uniqueName, Int_t systSwitch = 0, TString suffix = "pdf") {

    // systematic error eval
    TString systLabel                                           = "";
    if (systSwitch == 1) systLabel                              = "constShiftUp";
    if (systSwitch == 2) systLabel                              = "constShiftDown";
    if (systSwitch == 3) systLabel                              = "linShiftA";
    if (systSwitch == 4) systLabel                              = "linShiftB";
    if (systSwitch == 5) systLabel                              = "pol2ShiftA";
    if (systSwitch == 6) systLabel                              = "pol2ShiftB";
    if (systSwitch == 7) systLabel                              = "mtScalingUp";
    if (systSwitch == 8) systLabel                              = "mtScalingDown";

    cout << "systLabel: " << systLabel << endl;

    // get coll. sys, energy and centrality from list
    TString collSys, energy, cent, paramDirName, settingsFileName;
    TString listName                                            = list->GetName();
    TObjArray *arr                                              = listName.Tokenize("_");
    if (listName.BeginsWith("pp")) {
        collSys                                                 = "pp";
        energy                                                  = ((TObjString*)arr->At(1))->GetString();
        cent                                                    = "";
        if (systLabel.CompareTo("") == 0)   paramDirName        = Form("%s_%s", energy.Data(), uniqueName.Data());
        else                                paramDirName        = Form("%s_%s_%s", energy.Data(), uniqueName.Data(), systLabel.Data());
        settingsFileName                                        = Form("parametrizationSettings/%s_%s_%s_cocktail_settings.dat", collSys.Data(), energy.Data(), uniqueName.Data());
    } else {
        collSys                                                 = ((TObjString*)arr->At(0))->GetString();
        energy                                                  = ((TObjString*)arr->At(1))->GetString();
        cent                                                    = ((TObjString*)arr->At(2))->GetString();
        if (systLabel.CompareTo("") == 0)   paramDirName        = Form("%s_%s_%s", energy.Data(), cent.Data(), uniqueName.Data());
        else                                paramDirName        = Form("%s_%s_%s_%s", energy.Data(), cent.Data(), uniqueName.Data(), systLabel.Data());
        settingsFileName                                        = Form("parametrizationSettings/%s_%s_%s_%s_cocktail_settings.dat", collSys.Data(), energy.Data(), cent.Data(), uniqueName.Data());
    }

    // get correct settings file
    ifstream settings;
    settings.open(settingsFileName,ios_base::in);
    if (!settings) {
        cout << "ERROR: Parametrization settings file " << settingsFileName.Data() << " not found!" << endl;
        return;
    }

    // mt scale factor histogram
    TH1D* histoMtScaleFactor                                    = new TH1D("histoMtScaleFactor", "", 26, 0.5, 26.5);
    for (Int_t i=1; i<27; i++) histoMtScaleFactor->SetBinContent(i, -9999);
    histoMtScaleFactor->GetYaxis()->SetTitle("mt scaling factor");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(1,"111");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(2,"221");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(3,"113");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(4,"223");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(5,"331");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(6,"333");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(7,"443");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(8,"3212");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(9,"310");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(10,"2224");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(11,"2214");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(12,"1114");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(13,"2114");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(14,"213");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(15,"-213");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(16,"313");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(17,"130");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(18,"3122");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(19,"321");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(20,"-321");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(21,"-3334");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(22,"3334");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(23,"-3312");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(24,"3312");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(25,"3224");
    histoMtScaleFactor->GetXaxis()->SetBinLabel(26,"3114");

    // parametrization
    TF1* parametrization                                        = NULL;
    TF1* parametrization2                                       = NULL;
    TF1* ratio                                                  = NULL;
    Double_t mtScaleFactor                                      = 0.;
    Double_t mtScaleFactorUp                                    = 0.;
    Double_t mtScaleFactorDown                                  = 0.;
    Double_t mtScaleFactorSet                                   = 0.;

    // output list
    TList* parametrizationList                                  = new TList();

    // read from settings file
    TString inputString[8];
    for (Int_t i=0; i<8; i++) inputString[i]                    = "";
    TString tempString                                          = "";
    std::string line;
    for( std::string line; getline(settings, line); ) {
      // get settings
      settings >> tempString;
      if(strcmp(tempString.Data(),"")!=0){
        inputString[0]                                          = tempString;
        Int_t counter                                           = 1;
        while (tempString.CompareTo("stop")!=0 && counter < 8) {
            settings >> tempString;
            inputString[counter]                                = tempString;
            counter++;
        }

        // correctly set method
        if (inputString[5].CompareTo("-") == 0) inputString[5]  = "";
        if (inputString[7].CompareTo("-") == 0) inputString[7]  = "";

        // set parametrizations
        if (inputString[6].CompareTo("stop") == 0) {

            // take specified parametrization
            if(systLabel.CompareTo("")==0 || systLabel.Contains("mtScaling")){
                cout << "trying to find : " << Form("%s%sStat_Fit", inputString[4].Data(), inputString[5].Data()) << endl;
                parametrization2                                = (TF1*)list->FindObject(Form("%s%sStat_Fit", inputString[4].Data(), inputString[5].Data()));
            } else {
                cout << "trying to find : " << Form("%s%sStat_%s_Fit", inputString[4].Data(), inputString[5].Data(), systLabel.Data()) << endl;
                parametrization2                                = (TF1*)list->FindObject(Form("%s%sStat_%s_Fit", inputString[4].Data(), inputString[5].Data(), systLabel.Data()));
            }
            if (parametrization2) {
                parametrization                                 = (TF1*)parametrization2->Clone(Form("%s_pt",inputString[0].Data()));
                parametrizationList->Add(parametrization);
            } else {
                cout << "failed" << endl;
            }
        } else {

            // get ratio and parametrization
            if(systLabel.CompareTo("")==0 || systLabel.Contains("mtScaling")) {
                cout << "trying to find : " << Form("%s%sStat_Fit", inputString[4].Data(), inputString[5].Data()) << " and " << Form("%s%sStat_Fit", inputString[6].Data(), inputString[7].Data()) << endl;
                parametrization2                                = (TF1*)list->FindObject(Form("%s%sStat_Fit", inputString[4].Data(), inputString[5].Data()));
                ratio                                           = (TF1*)list->FindObject(Form("%s%sStat_Fit", inputString[6].Data(), inputString[7].Data()));
            } else {
                cout << "trying to find : " << Form("%s%sStat_%s_Fit", inputString[4].Data(), inputString[5].Data(), systLabel.Data()) << " and " << Form("%s%sStat_%s_Fit", inputString[6].Data(), inputString[7].Data(), systLabel.Data()) << endl;
                parametrization2                                = (TF1*)list->FindObject(Form("%s%sStat_%s_Fit", inputString[4].Data(), inputString[5].Data(), systLabel.Data()));
                ratio                                           = (TF1*)list->FindObject(Form("%s%sStat_%s_Fit", inputString[6].Data(), inputString[7].Data(), systLabel.Data()));
            }

            // calculate parametrization from ratio
            if (parametrization2 && ratio) {
                parametrization                                 = MultiplyTF1(ratio, parametrization2, Form("%s_pt",inputString[0].Data()));
                parametrizationList->Add(parametrization);
            } else {
                cout << "failed" << endl;
            }
        }

        // set mt scale factor (if specified)
        if (inputString[1].CompareTo("-") == 0) mtScaleFactor       = -9999;
        else                                    mtScaleFactor       = (Double_t)inputString[1].Atof();
        if (inputString[2].CompareTo("-") == 0) mtScaleFactorUp     = -9999;
        else                                    mtScaleFactorUp     = (Double_t)inputString[2].Atof();
        if (inputString[3].CompareTo("-") == 0) mtScaleFactorDown   = -9999;
        else                                    mtScaleFactorDown   = (Double_t)inputString[3].Atof();

        if (systSwitch == 7)                    mtScaleFactorSet    = mtScaleFactorUp;
        else if (systSwitch == 8)               mtScaleFactorSet    = mtScaleFactorDown;
        else                                    mtScaleFactorSet    = mtScaleFactor;

        // fill to histogram
        Int_t pdgCode                                           = inputString[0].Atoi();
        switch (pdgCode) {
            case 111:                                                   // pi0
                histoMtScaleFactor->SetBinContent(1, mtScaleFactorSet);
                break;
            case 221:                                                   // eta
                histoMtScaleFactor->SetBinContent(2, mtScaleFactorSet);
                break;
            case 113:                                                   // rho0
                histoMtScaleFactor->SetBinContent(3, mtScaleFactorSet);
                break;
            case 223:                                                   // omega
                histoMtScaleFactor->SetBinContent(4, mtScaleFactorSet);
                break;
            case 331:                                                   // eta'
                histoMtScaleFactor->SetBinContent(5, mtScaleFactorSet);
                break;
            case 333:                                                   // phi
                histoMtScaleFactor->SetBinContent(6, mtScaleFactorSet);
                break;
            case 443:                                                   // J/Psi
                histoMtScaleFactor->SetBinContent(7, mtScaleFactorSet);
                break;
            case 3212:                                                  // Sigma0
                histoMtScaleFactor->SetBinContent(8, mtScaleFactorSet);
                break;
            case 310:                                                   // K0s
                histoMtScaleFactor->SetBinContent(9, mtScaleFactorSet);
                break;
            case 2224:                                                  // Delta++
                histoMtScaleFactor->SetBinContent(10, mtScaleFactorSet);
                break;
            case 2214:                                                  // Delta+
                histoMtScaleFactor->SetBinContent(11, mtScaleFactorSet);
                break;
            case 1114:                                                  // Delta-
                histoMtScaleFactor->SetBinContent(12, mtScaleFactorSet);
                break;
            case 2114:                                                  // Delta0
                histoMtScaleFactor->SetBinContent(13, mtScaleFactorSet);
                break;
            case 213:                                                   // rho+
                histoMtScaleFactor->SetBinContent(14, mtScaleFactorSet);
                break;
            case -213:                                                  // rho-
                histoMtScaleFactor->SetBinContent(15, mtScaleFactorSet);
                break;
            case 313:                                                   // K*0
                histoMtScaleFactor->SetBinContent(16, mtScaleFactorSet);
                break;
            case 130:                                                   // K0l
                histoMtScaleFactor->SetBinContent(17, mtScaleFactorSet);
                break;
            case 3122:                                                  // Lambda
                histoMtScaleFactor->SetBinContent(18, mtScaleFactorSet);
                break;
            case 321:                                                   // K+
                histoMtScaleFactor->SetBinContent(19, mtScaleFactorSet);
                break;
            case -321:                                                  // K-
                histoMtScaleFactor->SetBinContent(20, mtScaleFactorSet);
                break;
            case -3334:                                                 // Omega+
                histoMtScaleFactor->SetBinContent(21, mtScaleFactorSet);
                break;
            case 3334:                                                  // Omega-
                histoMtScaleFactor->SetBinContent(22, mtScaleFactorSet);
                break;
            case -3312:                                                 // Xi+
                histoMtScaleFactor->SetBinContent(23, mtScaleFactorSet);
                break;
            case 3312:                                                  // Xi-
                histoMtScaleFactor->SetBinContent(24, mtScaleFactorSet);
                break;
            case 3224:                                                  // Sigma*+
                histoMtScaleFactor->SetBinContent(25, mtScaleFactorSet);
                break;
            case 3114:                                                  // Sigma*-
                histoMtScaleFactor->SetBinContent(26, mtScaleFactorSet);
                break;
        }

        // reset
        for (Int_t i=0; i<8; i++) inputString[i]                = "";
        mtScaleFactor                                           = -9999;
        mtScaleFactorUp                                         = -9999;
        mtScaleFactorDown                                       = -9999;
        mtScaleFactorSet                                        = -9999;
      }
    }

    // pt-y distributions
    TList*  ptYList                                             = new TList();
    TString pdgCodes[26]                                        = {"111", "221", "113", "223", "331", "333", "443", "3212", "310", "2224", "2214", "1114", "2114", "213", "-213", "313", "130", "3122",
                                                                    "321", "-321", "-3334", "3334", "-3312", "3312", "3224", "3114"};
    TH2F*   tempPtYHist                                         = NULL;
    if (inputList) {
        for (Int_t i=0; i<26; i++) {
            tempPtYHist                                         = (TH2F*)inputList->FindObject(Form("%s_pt_y", pdgCodes[i].Data()));
            if (tempPtYHist) ptYList->Add(tempPtYHist);
        }
    }

    // create directory for parametrizations
    gSystem->Exec("mkdir -p parametrizations");

    // plot parametrizations together
    Int_t centIte       = GetCentralityIterator(cent);
    TString centForPlot = fCentralityLatex[centIte];
    if(systSwitch == 0)
        PlotCocktailParametrizationListFullRange(parametrizationList,inputList,collSys,energy,centForPlot,suffix,paramDirName);

    // create output file containing only parametrizations
    TString energyOutputName                                    = GetOutputStringEnergy(energy);
    TFile* paramFile                                            = new TFile(Form("parametrizations/%s_%s.root", collSys.Data(), energyOutputName.Data()), "UPDATE");
    paramFile->cd();
    TDirectory*     paramDir                                    = (TDirectory*)paramFile->Get(Form("%s", paramDirName.Data()));
    if (!paramDir)  paramDir                                    = paramFile->mkdir(Form("%s", paramDirName.Data()));
    paramDir->cd();
    cout << "---- writing to " << paramFile->GetName() << "/" <<paramDirName.Data() << endl;
    histoMtScaleFactor->Write(histoMtScaleFactor->GetName(), TObject::kOverwrite);
    TF1* tempParam                                              = NULL;
    for (Int_t i=0; i<parametrizationList->GetEntries(); i++) {
        tempParam                                               = (TF1*)parametrizationList->At(i);
        tempParam->Write(tempParam->GetName(), TObject::kOverwrite);
    }

    paramFile->cd();
    TDirectory*     paramDirY                                   = (TDirectory*)paramFile->Get(Form("%s", energy.Data()));
    if (!paramDirY)  paramDirY                                  = paramFile->mkdir(Form("%s", energy.Data()));
    paramDirY->cd();

    TH2F* tempHist                                              = NULL;
    for (Int_t i=0; i<ptYList->GetEntries(); i++) {
        tempHist                                                = (TH2F*)ptYList->At(i);
        if (tempHist) tempHist->Write(tempHist->GetName(), TObject::kOverwrite);
    }

    paramFile->Close();

    delete paramFile;
    delete parametrizationList;
    delete histoMtScaleFactor;

}

//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_ReturnExtremeHisto(TH1 *hin, Float_t sign){
    Double_t ptlow          = 0;
    Double_t pthigh         = 0;
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.) continue;
        ptlow           = hin->GetBinLowEdge(ibin + 1);
        break;
    }
    for (Int_t ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
        if (hin->GetBinError(ibin + 1) <= 0.) continue;
        pthigh          = hin->GetBinLowEdge(ibin + 2);
        break;
    }

    Double_t mean       = hin->GetMean();
    Double_t maxdiff    = 0.;
    TH1 *hmax           = NULL;
    for (Int_t inode = 0; inode < hin->GetNbinsX(); inode++) {

        Double_t ptnode = hin->GetBinCenter(inode + 1);
        TH1 *hout       = (TH1 *)hin->Clone(Form("%s_extremehard", hin->GetName()));

        for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
            if (hin->GetBinError(ibin + 1) <= 0.) continue;
            Double_t val    = hin->GetBinContent(ibin + 1);
            Double_t err    = hin->GetBinError(ibin + 1);
            Double_t cen    = hin->GetBinCenter(ibin + 1);
            if (cen < ptnode)
                err         *= -1. + (cen - ptlow) / (ptnode - ptlow);
            else
                err         *= (cen - ptnode) / (pthigh - ptnode);

            hout->SetBinContent(ibin + 1, val + sign * err);
        }

        Double_t diff       = TMath::Abs(mean - hout->GetMean());
        if (diff > maxdiff) {
            //      printf("found max at %f\n", ptnode);
            if (hmax) delete hmax;
            hmax            = (TH1 *)hout->Clone("hmax");
            maxdiff         = diff;
        }
        delete hout;
    }
    return hmax;
}

//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_ReturnExtremeSoftHisto(TH1 *hin){
    return YieldMean_ReturnExtremeHisto(hin, -1);
}

//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_ReturnExtremeHardHisto(TH1 *hin){
    return YieldMean_ReturnExtremeHisto(hin, 1);
}


//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_LowExtrapolationHisto(   TH1 *h,
                                        TF1 *f,
                                        Double_t min,
                                        Double_t binwidth,
                                        TString nameHisto   = "hlo"
){
    /* find lowest edge in histo */
    Int_t binlo         = 0;
    Double_t lo         = 0;
    for (Int_t ibin = 1; ibin < h->GetNbinsX() + 1; ibin++) {
        if (h->GetBinContent(ibin) != 0.) {
            binlo = ibin;
            lo = h->GetBinLowEdge(ibin);
            break;
        }
    }

    Int_t nbins = (lo - min) / binwidth;
    if(nbins<1)
        return 0x0;
    TH1 *hlo = new TH1F(nameHisto.Data(), nameHisto.Data(), nbins, min, lo);

    /* integrate function in histogram bins */
    Double_t cont, err, width;
    for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
        width = hlo->GetBinWidth(ibin + 1);
        //cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (Double_t *)0, 1.e-6);
        cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), 1.e-6);
        err = f->IntegralError(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
        hlo->SetBinContent(ibin + 1, cont / width);
        hlo->SetBinError(ibin + 1, err / width);
    }

    return hlo;
}

//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_HighExtrapolationHisto(  TH1 *h,
                                        TF1 *f,
                                        Double_t max,
                                        Double_t binwidth
){
    /* find highest edge in histo */
    Int_t binhi     = 0;
    Double_t hi     = 0;
    for (Int_t ibin = h->GetNbinsX(); ibin > 0; ibin--) {
        if (h->GetBinContent(ibin) != 0.) {
            binhi = ibin + 1;
            hi = h->GetBinLowEdge(ibin + 1);
            break;
        }
    }
    if(max<hi) {
        Printf("Warning! You should probably set a higher max value (Max = %f, hi = %f)", max, hi);
        return 0x0;
    }
    Int_t nbins = (max - hi) / binwidth;
    if(nbins<1)
        return 0x0;
    TH1 *hhi = new TH1F("hhi", "", nbins, hi, max);

    /* integrate function in histogram bins */
    Double_t cont, err, width;
    for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
        width = hhi->GetBinWidth(ibin + 1);
        //cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (Double_t *)0, 1.e-6);
        cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), 1.e-6);
        err = f->IntegralError(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
        hhi->SetBinContent(ibin + 1, cont / width);
        hhi->SetBinError(ibin + 1, err / width);
    }

    return hhi;
}

//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_ReturnRandom(TH1 *hin){
    TH1 *hout = (TH1 *)hin->Clone("hout");
    hout->Reset();
    Double_t cont, err;
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.) continue;
        cont = hin->GetBinContent(ibin + 1);
        err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, gRandom->Gaus(cont, err));
        hout->SetBinError(ibin + 1, err);
    }
    return hout;
}

//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_ReturnCoherentRandom(TH1 *hin){
    if(!hin)
        return 0x0;
    TH1 *hout = (TH1 *)hin->Clone("hout");
    hout->Reset();
    Double_t cont, err, cohe;
    cohe = gRandom->Gaus(0., 1.);
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.) continue;
        cont = hin->GetBinContent(ibin + 1);
        err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, cont + cohe * err);
        hout->SetBinError(ibin + 1, err);
    }
    return hout;
}

//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_ReturnExtremeHighHisto(TH1 *hin){
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehigh", hin->GetName()));
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.) continue;
        Double_t val = hin->GetBinContent(ibin + 1);
        Double_t err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, val + err);
    }
    return hout;
}

//================================================================================================================
//
//================================================================================================================
TH1* YieldMean_ReturnExtremeLowHisto(TH1 *hin){
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremelow", hin->GetName()));
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.) continue;
        Double_t val = hin->GetBinContent(ibin + 1);
        Double_t err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, val - err);
    }
    return hout;
}


//================================================================================================================
//
//================================================================================================================
void YieldMean_IntegralMean(    TH1 *hdata,
                                TH1 *hlo,
                                TH1 *hhi,
                                Double_t &integral,
                                Double_t &mean,
                                Double_t &extra,
                                Bool_t printinfo        = kFALSE
){

    //*************************************************
    //compute integrals
    //*************************************************
    Double_t cont, err, width, cent;
    // integrated yield
    Double_t I = 0., IX = 0., Ierr = 0., IXerr = 0., Ilerr = 0., IXlerr = 0.;
    // mean pt
    Double_t M = 0., Merr = 0., Mlerr = 0., C;
    // extrapolated yield
    Double_t E = 0;
    Double_t dataonly=0.0;

    // integrate the data
    for (Int_t ibin = 0; ibin < hdata->GetNbinsX(); ibin++) {
        cent    = hdata->GetBinCenter(ibin + 1);
        width   = hdata->GetBinWidth(ibin + 1);
        cont    = width * hdata->GetBinContent(ibin + 1);
        err     = width * hdata->GetBinError(ibin + 1);
        if (err <= 0.) continue;
        I       += cont;
        IX      += cont * cent;
    }

    dataonly=I;
    // integrate low
    if(hlo) {
        for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
            cent    = hlo->GetBinCenter(ibin + 1);
            width   = hlo->GetBinWidth(ibin + 1);
            cont    = width * hlo->GetBinContent(ibin + 1);
            err     = width * hlo->GetBinError(ibin + 1);
            if (err <= 0.) continue;
            I       += cont;
            IX      += cont * cent;
            E       += cont;
        }
    }
    // integrate high
    if(printinfo)
        cout<<"low part data only = "<<dataonly<<" total = "<<I<<" ratio= "<<dataonly/I<<endl;
    if(hhi) {
        for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
            cent    = hhi->GetBinCenter(ibin + 1);
            width   = hhi->GetBinWidth(ibin + 1);
            cont    = width * hhi->GetBinContent(ibin + 1);
            err     = width * hhi->GetBinError(ibin + 1);
            if (err <= 0.) continue;
            I       += cont;
            IX      += cont * cent;
            E       += cont;
        }
    }
    // set values
    integral    = I;
    mean        = IX / I;
    extra       = E;
    if(printinfo)
        cout<<"low+high data only = "<<dataonly<<" total = "<<I<<" ratio= "<<dataonly/I<<endl;
}


//================================================================================================================
// Calculate integrated yields for different particle species
//================================================================================================================
TH1* YieldMean(     TObject* objStat,
                    TObject* objSys,
                    TF1 *fitInt             = NULL,
                    Double_t min            = 0.,
                    Double_t max            = 10.,
                    Double_t loprecision    = 0.01,
                    Double_t hiprecision    = 0.1,
                    Option_t *opt           = "0q",
                    TString logfilename     ="log.root",
                    Double_t minfit         = 0.0,
                    Double_t maxfit         = 10.0,
                    TString plotNameStat    = "",
                    TString plotNameSys     = "",
                    TString histNameOut     = "hOut",
                    TF1** fitFunctionSys    = NULL

               ){

    //*****************************************************************
    // Create histograms with systematic and statistical errors
    //*****************************************************************
    TH1D* hStat             = NULL;
    TH1D* hSys              = NULL;
    TString classNameStat   = objStat->ClassName();
    TString classNameSys    = objSys->ClassName();
    // fill histogram with statistical errors
    if (classNameStat.Contains("TH1")){
        hStat               = (TH1D*)objStat->Clone("stat");
    } else if (classNameStat.Contains("TGraphErrors")){
        Double_t ptBins[((TGraphErrors*)objStat)->GetN()+1];
        for (Int_t n = 0; n < ((TGraphErrors*)objStat)->GetN(); n++){
            ptBins[n]       = ((TGraphErrors*)objStat)->GetX()[n] - ((TGraphErrors*)objStat)->GetEX()[n];
        }
        ptBins[((TGraphErrors*)objStat)->GetN()]    = ((TGraphErrors*)objStat)->GetX()[((TGraphErrors*)objStat)->GetN()-1] + ((TGraphErrors*)objStat)->GetEX()[((TGraphErrors*)objStat)->GetN()-1];

        hStat               = new TH1D("stat", "stat", ((TGraphErrors*)objStat)->GetN(), ptBins);
        for (Int_t n = 0; n < ((TGraphErrors*)objStat)->GetN(); n++){
            hStat->SetBinContent(n+1, ((TGraphErrors*)objStat)->GetY()[n]);
            hStat->SetBinError(n+1, ((TGraphErrors*)objStat)->GetEY()[n]);
        }
    } else if (classNameStat.Contains("TGraphAsymmErrors")){
        Double_t ptBins[((TGraphAsymmErrors*)objStat)->GetN()+1];
        for (Int_t n = 0; n < ((TGraphAsymmErrors*)objStat)->GetN(); n++){
            ptBins[n]       = ((TGraphAsymmErrors*)objStat)->GetX()[n] - ((TGraphAsymmErrors*)objStat)->GetEXlow()[n];
            cout << ptBins[n]<< ",";
        }
        ptBins[((TGraphAsymmErrors*)objStat)->GetN()]    = ((TGraphAsymmErrors*)objStat)->GetX()[((TGraphAsymmErrors*)objStat)->GetN()-1] + ((TGraphAsymmErrors*)objStat)->GetEXhigh()[((TGraphAsymmErrors*)objStat)->GetN()-1];
        cout << ptBins[((TGraphAsymmErrors*)objStat)->GetN()] << endl;;
        hStat               = new TH1D("stat", "stat", ((TGraphAsymmErrors*)objStat)->GetN(), ptBins);
        for (Int_t n = 0; n < ((TGraphAsymmErrors*)objStat)->GetN(); n++){
            hStat->SetBinContent(n+1, ((TGraphAsymmErrors*)objStat)->GetY()[n]);
            hStat->SetBinError(n+1, ((TGraphAsymmErrors*)objStat)->GetEYhigh()[n]);
        }
        ((TGraphAsymmErrors*)objStat)->Print();
    } else {
        return NULL;
    }

    // fill histogram with statistical errors
    if (classNameSys.Contains("TH1")){
        hSys                = (TH1D*)objSys->Clone("sys");
    } else if (classNameSys.Contains("TGraphErrors")){
        Double_t ptBins[((TGraphErrors*)objSys)->GetN()+1];
        for (Int_t n = 0; n < ((TGraphErrors*)objSys)->GetN(); n++){
            ptBins[n]       = ((TGraphErrors*)objSys)->GetX()[n] - ((TGraphErrors*)objSys)->GetEX()[n];
        }
        ptBins[((TGraphErrors*)objSys)->GetN()]    = ((TGraphErrors*)objSys)->GetX()[((TGraphErrors*)objSys)->GetN()-1] + ((TGraphErrors*)objSys)->GetEX()[((TGraphErrors*)objSys)->GetN()-1];

        hSys               = new TH1D("sys", "sys", ((TGraphErrors*)objSys)->GetN(), ptBins);
        for (Int_t n = 0; n < ((TGraphErrors*)objSys)->GetN(); n++){
            hSys->SetBinContent(n+1, ((TGraphErrors*)objSys)->GetY()[n]);
            hSys->SetBinError(n+1, ((TGraphErrors*)objSys)->GetEY()[n]);
        }
    } else if (classNameSys.Contains("TGraphAsymmErrors")){
        Double_t ptBins[((TGraphAsymmErrors*)objSys)->GetN()+1];
        for (Int_t n = 0; n < ((TGraphAsymmErrors*)objSys)->GetN(); n++){
            ptBins[n]       = ((TGraphAsymmErrors*)objSys)->GetX()[n] - ((TGraphAsymmErrors*)objSys)->GetEXlow()[n];
        }
        ptBins[((TGraphAsymmErrors*)objSys)->GetN()]    = ((TGraphAsymmErrors*)objSys)->GetX()[((TGraphAsymmErrors*)objSys)->GetN()-1] + ((TGraphAsymmErrors*)objSys)->GetEXhigh()[((TGraphAsymmErrors*)objSys)->GetN()-1];
        hSys               = new TH1D("sys", "sys", ((TGraphAsymmErrors*)objSys)->GetN(), ptBins);
        for (Int_t n = 0; n < ((TGraphAsymmErrors*)objSys)->GetN(); n++){
            hSys->SetBinContent(n+1, ((TGraphAsymmErrors*)objSys)->GetY()[n]);
            hSys->SetBinError(n+1, ((TGraphAsymmErrors*)objSys)->GetEYhigh()[n]);
        }
        ((TGraphAsymmErrors*)objSys)->Print();
    } else {
        return NULL;
    }

    // set minimum and maximum of spectrum
    if(maxfit>max)
        max         = maxfit;
    if(minfit<min)
        min         = minfit;

    // set many iterations when fitting the data so we don't
    // stop minimization with MAX_CALLS
    TVirtualFitter::SetMaxIterations(1000000);

    // create output histo
    Double_t integral, mean, extra;
    TH1 *hout               = new TH1D(histNameOut.Data(), "", 9, 0, 9);
    TString binLabels[9]    = {"dN/dy", "stat Err dN/dy", "sys Err dN/dy h.", "sys Err dN/dy l.", "<p_{T}>", "stat Err <p_{T}>", "sys Err <p_{T}> h.", "sys Err <p_{T}> l.", "extrapolated dN/dy"};
    for (Int_t i = 1; i< 10; i++){
        hout->GetXaxis()->SetBinLabel(i,binLabels[i-1]);
    }
    // create histo with stat+sys errors
    TH1 *htot       = (TH1 *)hStat->Clone(Form("%sfittedwith%s",hStat->GetName(),fitInt->GetName()));
    for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
        htot->SetBinError(ibin + 1, TMath::Sqrt(hSys->GetBinError(ibin + 1) * hSys->GetBinError(ibin + 1) + hStat->GetBinError(ibin + 1) * hStat->GetBinError(ibin + 1)));
    }

    //*****************************************************************
    // measure the central value
    //*****************************************************************
    Int_t fitres;
    Int_t trials    = 0;
    do {
        fitres      = htot->Fit(fitInt, opt,"",minfit,maxfit);
        Printf("Trial: %d", trials++);
        if(trials > 10) {
            Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
            break;
        }
    }
    while (fitres != 0);
    TFile* filewithfits     = TFile::Open(logfilename.Data(),"UPDATE");
    htot->Write();
    filewithfits->Close();
    delete filewithfits;

    cout<<" Fit sys+stat for " <<fitInt->GetName()<<endl;
    cout<<"NDF="<<fitInt->GetNDF()<<" Chi^2="<<fitInt->GetChisquare()<<" Chi^2/NDF="<<fitInt->GetChisquare()/fitInt->GetNDF()<<endl;

    TH1 *hlowTot            = YieldMean_LowExtrapolationHisto(htot, fitInt, min, loprecision,"hlowTot");
    TH1 *hhighTot           = YieldMean_HighExtrapolationHisto(htot, fitInt, max, hiprecision);
    YieldMean_IntegralMean(htot, hlowTot, hhighTot, integral, mean, extra, kTRUE);
    hout->SetBinContent(kYield, integral);
    hout->SetBinContent(kMean, mean);
    hout->SetBinContent(kExtra, extra);

    delete hlowTot;
    delete hhighTot;
    //*************************************************
    // Statistical errors
    //*************************************************
    TCanvas *cCanvasStat = new TCanvas("cCanvasStat");
    cCanvasStat->Divide(2, 1);

    // fit with stat error
    trials          = 0;
    do {
        fitres      = hStat->Fit(fitInt, opt,"",minfit,maxfit);
        Printf("Trial: %d", trials++);
        if(trials > 10) {
            Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
            break;
        }
    }
    while (fitres != 0);
    TH1 *hlowStat           = YieldMean_LowExtrapolationHisto(hStat, fitInt, min, loprecision, "hlowStat");
    TH1 *hhighStat          = YieldMean_HighExtrapolationHisto(hStat, fitInt, max, hiprecision);

    // random generation with integration (coarse)
    TH1 *hIntegral_tmp  = new TH1F("hIntegral_tmp", "", 1000, 0.75 * integral, 1.25 * integral);
    TH1 *hMean_tmp      = new TH1F("hMean_tmp", "", 1000, 0.75 * mean, 1.25 * mean);
    for (Int_t irnd = 0; irnd < 100; irnd++) {
        /* get random histogram */
        TH1 *hrnd       = YieldMean_ReturnRandom(hStat);
        /* fit */
        TH1 *hrndlo     = YieldMean_ReturnCoherentRandom(hlowStat);
        TH1 *hrndhi     = YieldMean_ReturnCoherentRandom(hhighStat);
        /* integrate */
        YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean, extra);
        hIntegral_tmp->Fill(integral);
        hMean_tmp->Fill(mean);
        delete hrnd;
        delete hrndlo;
        delete hrndhi;
    }
    // random generation with integration (fine)
    TH1 *hIntegral      = new TH1F("hIntegral", "", 100,
                                    hIntegral_tmp->GetMean() - 10. * hIntegral_tmp->GetRMS(),
                                    hIntegral_tmp->GetMean() + 10. * hIntegral_tmp->GetRMS());
    TH1 *hMean          = new TH1F("hMean", "", 100,
                                    hMean_tmp->GetMean() - 10. * hMean_tmp->GetRMS(),
                                    hMean_tmp->GetMean() + 10. * hMean_tmp->GetRMS());
    for (Int_t irnd = 0; irnd < 1000; irnd++) {
        /* get random histogram */
        TH1 *hrnd       = YieldMean_ReturnRandom(hStat);
        /* fit */
        TH1 *hrndlo     = YieldMean_ReturnCoherentRandom(hlowStat);
        TH1 *hrndhi     = YieldMean_ReturnCoherentRandom(hhighStat);
        /* integrate */
        YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean, extra);
        hIntegral->Fill(integral);
        hMean->Fill(mean);
        delete hrnd;
        delete hrndlo;
        delete hrndhi;
    }

    TF1 *gaus       = (TF1 *)gROOT->GetFunction("gaus");

    cCanvasStat->cd(1);
    hIntegral->Fit(gaus, "q");
    integral        = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kYieldStat, integral);

    cCanvasStat->cd(2);
    hMean->Fit(gaus, "q");
    mean            = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kMeanStat, mean);

    cCanvasStat->SaveAs(plotNameStat.Data());
    delete hlowStat;
    delete hhighStat;
    delete hIntegral_tmp;
    delete hMean_tmp;
    delete hIntegral;
    delete hMean;
    delete cCanvasStat;

    // Systematics

    TCanvas *cCanvasSys = new TCanvas("cCanvasYieldSys");
    cCanvasSys->Divide(2, 1);
    cCanvasSys->cd(1)->DrawFrame(min, 1.e-8, max, 1.e3);
    cCanvasSys->cd(1)->SetLogy();
    hSys->SetMarkerStyle(20);
    hSys->SetMarkerColor(1);
    hSys->SetMarkerSize(1);
    hSys->Draw("same");
    cCanvasSys->cd(2)->DrawFrame(min, 1.e-8, max, 1.e3);
    cCanvasSys->cd(2)->SetLogy();
    hSys->Draw("same");

    // Systematics error up
    TH1 *hhigh      = YieldMean_ReturnExtremeHighHisto(hSys);
    trials = 0;
    do {
        fitres = hhigh->Fit(fitInt, opt,"",minfit,maxfit);
        Printf("Trial: %d", trials++);
        if(trials > 10) {
            Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
            break;
        }
    } while (fitres != 0);
    cout << "extreme high" << endl;

    TH1 *hlowMaxSys         = YieldMean_LowExtrapolationHisto(hhigh, fitInt, min, loprecision, "hlowMaxSys");
    TH1 *hhighMaxSys        = YieldMean_HighExtrapolationHisto(hhigh, fitInt, max, hiprecision);
    YieldMean_IntegralMean(hhigh, hlowMaxSys, hhighMaxSys, integral, mean, extra);
    integral                = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysHi, integral);
    delete hlowMaxSys;
    delete hhighMaxSys;

    cCanvasSys->cd(1);
    fitInt->SetLineColor(2);
    fitInt->DrawCopy("same");

    // set return function for maximum syst shift
    fitFunctionSys[0]       = (TF1*)fitInt->Clone();
    fitFunctionSys[0]->SetName(Form("fit_%s_max",histNameOut.Data()));

    // Systematic error hard
    TH1 *hhard = YieldMean_ReturnExtremeHardHisto(hSys);
    trials = 0;
    do {
        fitres = hhard->Fit(fitInt, opt,"",minfit,maxfit);
        Printf("Trial: %d", trials++);
        if(trials > 10) {
            Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
            break;
        }
    } while (fitres != 0);
    cout << "extreme hard" << endl;
    TH1 *hlowHard           = YieldMean_LowExtrapolationHisto(hhard, fitInt, min, loprecision, "hlowHard");
    TH1 *hhighHard          = YieldMean_HighExtrapolationHisto(hhard, fitInt, max, hiprecision);
    YieldMean_IntegralMean(hhard, hlowHard, hhighHard, integral, mean, extra);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysHi, mean);
    delete hlowHard;
    delete hhighHard;

    cCanvasSys->cd(2);
    fitInt->SetLineColor(2);
    fitInt->DrawCopy("same");
    // set return function for maximum syst shift
    fitFunctionSys[2]       = (TF1*)fitInt->Clone();
    fitFunctionSys[2]->SetName(Form("fit_%s_hard",histNameOut.Data()));


    // Systematic error low
    TH1 *hlow = YieldMean_ReturnExtremeLowHisto(hSys);
    trials = 0;
    do {
        fitres = hlow->Fit(fitInt, opt,"",minfit,maxfit);
        Printf("Trial: %d", trials++);
        if(trials > 10) {
            Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
            break;
        }
    } while (fitres != 0);
    cout << "extreme low" << endl;
    TH1 *hlowMinSys         = YieldMean_LowExtrapolationHisto(hlow, fitInt, min, loprecision, "hlowMinSys");
    TH1 *hhighMinSys        = YieldMean_HighExtrapolationHisto(hlow, fitInt, max, hiprecision);
    YieldMean_IntegralMean(hlow, hlowMinSys, hhighMinSys, integral, mean, extra);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysLo, integral);
    delete hlowMinSys;
    delete hhighMinSys;

    cCanvasSys->cd(1);
    fitInt->SetLineColor(4);
    fitInt->DrawCopy("same");
    fitFunctionSys[1]       = (TF1*)fitInt->Clone();
    fitFunctionSys[1]->SetName(Form("fit_%s_low",histNameOut.Data()));

    // Systematic error soft
    TH1 *hsoft = YieldMean_ReturnExtremeSoftHisto(hSys);
    trials = 0;
    do {
        fitres = hsoft->Fit(fitInt, opt,"",minfit,maxfit);
        Printf("Trial: %d", trials++);
        if(trials > 10) {
            Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
            break;
        }
    } while (fitres != 0);
    cout << "extreme soft" << endl;

    TH1 *hlowSoft           = YieldMean_LowExtrapolationHisto(hsoft, fitInt, min, loprecision, "hlowSoft");
    TH1 *hhighSoft          = YieldMean_HighExtrapolationHisto(hsoft, fitInt, max, hiprecision);
    YieldMean_IntegralMean(hsoft, hlowSoft, hhighSoft, integral, mean, extra);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysLo, mean);
    delete hlowSoft;
    delete hhighSoft;

    cCanvasSys->cd(2);
    fitInt->SetLineColor(4);
    fitInt->DrawCopy("same");
    fitFunctionSys[3]       = (TF1*)fitInt->Clone();
    fitFunctionSys[3]->SetName(Form("fit_%s_soft",histNameOut.Data()));

    cCanvasSys->SaveAs(plotNameSys.Data());
    delete cCanvasSys;

    return hout;

}
