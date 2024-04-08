/*****************************************************************************************************************************
 ******    Friederike Bock, friederike.bock@cern.ch                                                                      *****
 ******    Mike Sas, mike.sas@cern.ch                                                                                    *****
 ******    Lucas Altenkaemper, lucas.altenkaemper@cern.ch                                                                *****
 *****************************************************************************************************************************/
#include "CocktailFunctions.h"
#include "CocktailPlotting.h"
//#include "CocktailFitting.h" // is included in CocktailPlotting.h because there a function is needed
#include <iostream>

void CocktailInputParametrization(  TString inputFileName                       = "CocktailInputPP.root",
                                    TString energy                              = "",
                                    TString suffix                              = "eps",
                                    TString parametrizationSettingsFile         = "",
                                    TString parametrizationSettingsFileRatio    = "",
                                    TString listOfUniqueNames                   = "allParams.txt",
                                    Bool_t doPtDepSysErrFactorPlotsSep          = kFALSE,
                                    Bool_t doParamRatioPlotsSep                 = kFALSE,
                                    Int_t plotXaxisRatios                       = 0,
                                    Bool_t thesisPlots                          = kFALSE,
                                    Bool_t convertToHistograms                  = kFALSE
                                 ) {

    counterGraphs = 0;
    vecGraphs.clear();

    if (thesisPlots) isThesis = kTRUE;

    std::vector<TString> vecUniqueNames(0);
    std::vector<Int_t> vecUniqueNamesEnableSys(0);

    // input file
    TFile* inputFile                                        = NULL;
    if (inputFileName.CompareTo("") == 0) {
        std::cout << "ERROR: No input file specified, returning!" << std::endl;
        return;
    } else {
        inputFile                                           = new TFile(inputFileName.Data());
        std::cout << inputFileName.Data() << std::endl;
    }

    if (!inputFile) {
        std::cout << "ERROR: Input file " << inputFileName.Data() << "not found, returning!" << std::endl;
        return;
    }

    // check parametrization settings file and initialize
    Bool_t doSpectra, doRatios;
    if (parametrizationSettingsFile.CompareTo("") != 0) {
        if (InitializeFitting(parametrizationSettingsFile)) {
            std::cout << "Will parametrize spectra." << std::endl;
            doSpectra                                       = kTRUE;
        } else {
            doSpectra                                       = kFALSE;
        }
    } else {
        std::cout << "WARNING: No parametrization settings file specified, skipping parametrization of spectra!" << std::endl;
        doSpectra                                           = kFALSE;
    }
    if (parametrizationSettingsFileRatio.CompareTo("") != 0) {
        if (InitializeRatioFitting(parametrizationSettingsFileRatio)) {
            std::cout << "Will parametrize ratios." << std::endl;
            doRatios                                        = kTRUE;
        } else {
            doRatios                                        = kFALSE;
        }
    } else {
        std::cout << "WARNING: No ratio parametrization settings file specified, skipping parametrization of ratios!" << std::endl;
        doRatios                                           = kFALSE;
    }

    // return if neither spectra nor ratios will be parametrized
    if (!doSpectra && !doRatios) {
        std::cout << "ERROR: Neither spectra nor ratios will be parametrized, returning!" << std::endl;
        return;
    }

    if (listOfUniqueNames.CompareTo("") == 0) {
        std::cout << "ERROR: No file name for the uniqueName list is given" << std::endl;
        return;
    }


    cout << "INFO: You have chosen the given the following file with the list of uniqueNames: " << listOfUniqueNames.Data() << endl;
    ifstream fileConfigQA;
    fileConfigQA.open(listOfUniqueNames,ios_base::in);
    if (!fileConfigQA) {
        cout << "ERROR: settings " << listOfUniqueNames.Data() << " not found!" << endl;
        return;
    }

    // read settings from file
    for( TString tempLine; tempLine.ReadLine(fileConfigQA, kTRUE); ) {
        // check if line should be considered
        if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
            continue;
        }
        cout << tempLine.Data() << endl;
        TObjArray *tempArr  = tempLine.Tokenize("\t");
        if(tempArr->GetEntries()<1){
            cout << "nothing to be done" << endl;
            delete tempArr;
            continue;
        } else if (tempArr->GetEntries() == 1 ){
            // Separate the string according to space
            tempArr       = tempLine.Tokenize(" ");
            if(tempArr->GetEntries() == 1){
                cout << "ERROR you are using an unventional method of separating the 2 expected arguments per line, skipping this input line, please use spaces or tabs to separate the two arguments" << endl;
                continue;
            }
        }
        cout << "writing config: " << tempLine.Data() << "\t to vectors." << endl;
        vecUniqueNames.push_back(((TObjString*)tempArr->At(0))->GetString());
        vecUniqueNamesEnableSys.push_back(((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi());
    }
    if ((Int_t)vecUniqueNames.size() == 0){
        cout << "ERROR: list of unique names is empty! Aborting!" << endl;
        return;
    }

    // output file
    TString outputFileName                                  = inputFileName.ReplaceAll(".root", 5, "_Param_", 7);
    outputFileName                                          = Form("%s%s.root", outputFileName.Data(), energy.Data());
    TFile *outputFile                                       = new TFile(outputFileName.Data(), "RECREATE");

    TString textOutputFileName                              = outputFileName.ReplaceAll(".root", 5, ".dat", 4);
    ofstream textoutput;
    textoutput.open(textOutputFileName);

    // declare input and output lists
    Int_t numberOfLists                                     = (Int_t)((TList*)inputFile->GetListOfKeys())->GetEntries();
    TList** inputList                                       = new TList*[numberOfLists];
    TList** outputList                                      = new TList*[numberOfLists];

    TString tempStatClass                                   = "";
    TString tempSysClass                                    = "";

    // declare temp spectrum histograms, graphs and fits
    TObject*    tempSpecStat                                = NULL;
    TObject*    tempSpecSys                                 = NULL;
    TObject*    tempSpecStatPtConstShiftUp                  = NULL;
    TObject*    tempSpecStatPtConstShiftDown                = NULL;
    TObject*    tempSpecStatPtLinShiftA                     = NULL;
    TObject*    tempSpecStatPtLinShiftB                     = NULL;
    TObject*    tempSpecStatPtPol2ShiftA                    = NULL;
    TObject*    tempSpecStatPtPol2ShiftB                    = NULL;
    TF1*        tempSpecFactorLinShiftA                     = NULL;
    TF1*        tempSpecFactorLinShiftB                     = NULL;
    TF1*        tempSpecFactorPol2ShiftA                    = NULL;
    TF1*        tempSpecFactorPol2ShiftB                    = NULL;
    TH1D*       tempSpecStatRatioToFit                      = NULL;
    TH1D*       tempSpecStatPtConstShiftUpRatioToFit        = NULL;
    TH1D*       tempSpecStatPtConstShiftDownRatioToFit      = NULL;
    TH1D*       tempSpecStatPtLinShiftARatioToFit           = NULL;
    TH1D*       tempSpecStatPtLinShiftBRatioToFit           = NULL;
    TH1D*       tempSpecStatPtPol2ShiftARatioToFit          = NULL;
    TH1D*       tempSpecStatPtPol2ShiftBRatioToFit          = NULL;
    TF1*        tempSpecFit                                 = NULL;
    TF1*        tempSpecFitPtConstUp                        = NULL;
    TF1*        tempSpecFitPtConstDown                      = NULL;
    TF1*        tempSpecFitPtLinA                           = NULL;
    TF1*        tempSpecFitPtLinB                           = NULL;
    TF1*        tempSpecFitPtPol2A                          = NULL;
    TF1*        tempSpecFitPtPol2B                          = NULL;

    // declare temp ratio histograms, graphs and fits
    TObject*    tempRatioStat                               = NULL;
    TObject*    tempRatioSys                                = NULL;
    TObject*    tempRatioStatPtConstShiftUp                 = NULL;
    TObject*    tempRatioStatPtConstShiftDown               = NULL;
    TObject*    tempRatioStatPtLinShiftA                    = NULL;
    TObject*    tempRatioStatPtLinShiftB                    = NULL;
    TObject*    tempRatioStatPtPol2ShiftA                   = NULL;
    TObject*    tempRatioStatPtPol2ShiftB                   = NULL;
    TF1*        tempRatioFactorLinShiftA                    = NULL;
    TF1*        tempRatioFactorLinShiftB                    = NULL;
    TF1*        tempRatioFactorPol2ShiftA                   = NULL;
    TF1*        tempRatioFactorPol2ShiftB                   = NULL;
    TH1D*       tempRatioStatRatioToFit                     = NULL;
    TH1D*       tempRatioStatPtConstShiftUpRatioToFit       = NULL;
    TH1D*       tempRatioStatPtConstShiftDownRatioToFit     = NULL;
    TH1D*       tempRatioStatPtLinShiftARatioToFit          = NULL;
    TH1D*       tempRatioStatPtLinShiftBRatioToFit          = NULL;
    TH1D*       tempRatioStatPtPol2ShiftARatioToFit         = NULL;
    TH1D*       tempRatioStatPtPol2ShiftBRatioToFit         = NULL;
    TF1*        tempRatioFit                                = NULL;
    TF1*        tempRatioFitPtConstUp                       = NULL;
    TF1*        tempRatioFitPtConstDown                     = NULL;
    TF1*        tempRatioFitPtLinA                          = NULL;
    TF1*        tempRatioFitPtLinB                          = NULL;
    TF1*        tempRatioFitPtPol2A                         = NULL;
    TF1*        tempRatioFitPtPol2B                         = NULL;

    // loop over input lists
    TIter keyIter(inputFile->GetListOfKeys());
    TKey* key                                               = NULL;
    Int_t listIter                                          = 0;
    while ( (key = (TKey*)keyIter()) ) {

        // input list
        inputList[listIter]                                 = (TList*)key->ReadObj();

        // check for energy and set output list
        if (((TString)inputList[listIter]->GetName()).Contains(energy)) {
            outputList[listIter]                            = new TList();
            outputList[listIter]->SetName(inputList[listIter]->GetName());
        } else {
            outputList[listIter]                            = NULL;
            continue;
        }

        // create directory for plots
        TString outputDir                                   = Form("plots/%s/%s",suffix.Data(),inputList[listIter]->GetName());
        gSystem->Exec("mkdir -p "+outputDir);

        // skip empty or corrupted input lists
        if (!inputList[listIter])
            continue;
        if (inputList[listIter]->IsEmpty() || inputList[listIter]->IsZombie()) {
            cout << inputList[listIter]->GetName() << " empty or zombie" << endl;
            continue;
        }
        cout << "---- " << inputList[listIter]->GetName() << endl;

        // get coll. sys, energy and centrality from list
        TString fCollSysTemp, fEnergyTemp, fCentralityTemp;
        TString inputListName                               = inputList[listIter]->GetName();
        TObjArray *arr                                      = inputListName.Tokenize("_");

        if (inputListName.BeginsWith("pp")) {
            fCollSysTemp                                    = "pp";
            fEnergyTemp                                     = ((TObjString*)arr->At(1))->GetString();
            fCentralityTemp                                 = "";
        } else {
            fCollSysTemp                                    = ((TObjString*)arr->At(0))->GetString();
            fEnergyTemp                                     = ((TObjString*)arr->At(1))->GetString();
            fCentralityTemp                                 = ((TObjString*)arr->At(2))->GetString();
        }

        // loop over particles in input list
        for (Int_t particleIter = 0; particleIter < nParticles; particleIter++) {
            for (Int_t methodIter = 0; methodIter < nMethods; methodIter++) {
                // check if stat&sys spectra for particle are contained in input list
                if ( inputList[listIter]->FindObject(Form("%s%sStat", fParticle[particleIter].Data(), fMethod[methodIter].Data())) && inputList[listIter]->FindObject(Form("%s%sSys", fParticle[particleIter].Data(), fMethod[methodIter].Data())) ) {

                    // for testing purposes
                    cout << "-------- " << fParticle[particleIter].Data() << " " << fMethod[methodIter].Data() << " spectra found" << endl;

                    // check for pt-const rel. syst. err
                    Double_t ptConstRelSyst                 = initPtConstRelSyst[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])];

                    // get type stat&sys
                    tempStatClass                           = ((TObject*)inputList[listIter]->FindObject(Form("%s%sStat", fParticle[particleIter].Data(), fMethod[methodIter].Data())))->ClassName();
                    tempSysClass                            = ((TObject*)inputList[listIter]->FindObject(Form("%s%sSys",  fParticle[particleIter].Data(), fMethod[methodIter].Data())))->ClassName();

                    // get and shift spectra using sys errs
                    tempSpecStat                            = (TObject*)inputList[listIter]->FindObject(Form("%s%sStat", fParticle[particleIter].Data(), fMethod[methodIter].Data()));
                    tempSpecSys                             = (TObject*)inputList[listIter]->FindObject(Form("%s%sSys",  fParticle[particleIter].Data(), fMethod[methodIter].Data()));
                    if (tempStatClass.Contains("TH1") && tempSysClass.Contains("TH1")) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else if (tempStatClass.Contains("TH1") && tempSysClass.CompareTo("TGraphErrors") == 0) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else if (tempStatClass.Contains("TH1") && tempSysClass.CompareTo("TGraphAsymmErrors") == 0) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TH1D*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else if (tempStatClass.CompareTo("TGraphErrors") == 0 && tempSysClass.Contains("TH1")) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else if (tempStatClass.CompareTo("TGraphErrors") == 0 && tempSysClass.CompareTo("TGraphErrors") == 0) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else if (tempStatClass.CompareTo("TGraphErrors") == 0 && tempSysClass.CompareTo("TGraphAsymmErrors") == 0) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TGraphErrors*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0 && tempSysClass.Contains("TH1")) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TH1D*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TH1D*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0 && tempSysClass.CompareTo("TGraphErrors") == 0) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0 && tempSysClass.CompareTo("TGraphAsymmErrors") == 0) {
                        // calculate pt-dep. factors for syst. err.
                        tempSpecFactorLinShiftA             = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 1, kTRUE);
                        tempSpecFactorLinShiftB             = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 1, kFALSE);
                        tempSpecFactorPol2ShiftA            = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 2, kTRUE);
                        tempSpecFactorPol2ShiftB            = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempSpecStat, 2, kFALSE);

                        // shift spectra
                        tempSpecStatPtConstShiftUp          = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst);
                        tempSpecStatPtConstShiftDown        = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst);
                        tempSpecStatPtLinShiftA             = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorLinShiftA, 1);
                        tempSpecStatPtLinShiftB             = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorLinShiftB, 1);
                        tempSpecStatPtPol2ShiftA            = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kTRUE,     ptConstRelSyst, tempSpecFactorPol2ShiftA, 2);
                        tempSpecStatPtPol2ShiftB            = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, kFALSE,    ptConstRelSyst, tempSpecFactorPol2ShiftB, 2);
                    } else {
                        cout << "stat and sys type not known, skipping" << endl;
                        continue;
                    }

                    // write spectra to list
                    if (tempSpecStat)                   outputList[listIter]->Add(tempSpecStat);
                    if (tempSpecStatPtConstShiftUp)     outputList[listIter]->Add(tempSpecStatPtConstShiftUp);
                    if (tempSpecStatPtConstShiftDown)   outputList[listIter]->Add(tempSpecStatPtConstShiftDown);
                    if (tempSpecStatPtLinShiftA)        outputList[listIter]->Add(tempSpecStatPtLinShiftA);
                    if (tempSpecFactorLinShiftA)        outputList[listIter]->Add(tempSpecFactorLinShiftA);
                    if (tempSpecFactorLinShiftB)        outputList[listIter]->Add(tempSpecFactorLinShiftB);
                    if (tempSpecStatPtLinShiftB)        outputList[listIter]->Add(tempSpecStatPtLinShiftB);
                    if (tempSpecStatPtPol2ShiftA)       outputList[listIter]->Add(tempSpecStatPtPol2ShiftA);
                    if (tempSpecStatPtPol2ShiftB)       outputList[listIter]->Add(tempSpecStatPtPol2ShiftB);
                    if (tempSpecFactorPol2ShiftA)       outputList[listIter]->Add(tempSpecFactorPol2ShiftA);
                    if (tempSpecFactorPol2ShiftB)       outputList[listIter]->Add(tempSpecFactorPol2ShiftB);

                    // plot pt dept sys err factors
                    TString nameFactorPlot                  = "";
                    if (doPtDepSysErrFactorPlotsSep) {

                        nameFactorPlot                      = Form("spectra_%s%s_ptDepFactorLin",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                        ProducePlotPtDepSystFactors(tempSpecFactorLinShiftA,tempSpecFactorLinShiftB,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameFactorPlot);

                        nameFactorPlot                      = Form("spectra_%s%s_ptDepFactorPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                        ProducePlotPtDepSystFactors(NULL,NULL,tempSpecFactorPol2ShiftA,tempSpecFactorPol2ShiftB,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameFactorPlot);

                    } else {

                        nameFactorPlot                      = Form("spectra_%s%s_ptDepFactorLinPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                        ProducePlotPtDepSystFactors(tempSpecFactorLinShiftA,tempSpecFactorLinShiftB,tempSpecFactorPol2ShiftA,tempSpecFactorPol2ShiftB,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameFactorPlot);
                    }

                    // initializing fit function
                    TString type                            = initFitFunction[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])];
                    if (doSpectra && type.CompareTo("")!=0) {

                        // fit spectra
                        tempSpecFit                         = FitObject(tempSpecStat,Form("%s_Fit",tempSpecStat->GetName()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type,convertToHistograms);
                        if (tempSpecStatPtConstShiftUp)
                            tempSpecFitPtConstUp            = FitObject(tempSpecStatPtConstShiftUp,Form("%s_Fit",tempSpecStatPtConstShiftUp->GetName()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type,convertToHistograms);
                        if (tempSpecStatPtConstShiftDown)
                            tempSpecFitPtConstDown          = FitObject(tempSpecStatPtConstShiftDown,Form("%s_Fit",tempSpecStatPtConstShiftDown->GetName()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type,convertToHistograms);
                        if (tempSpecStatPtLinShiftA)
                            tempSpecFitPtLinA               = FitObject(tempSpecStatPtLinShiftA,Form("%s_Fit",tempSpecStatPtLinShiftA->GetName()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type,convertToHistograms);
                        if (tempSpecStatPtLinShiftB)
                            tempSpecFitPtLinB               = FitObject(tempSpecStatPtLinShiftB,Form("%s_Fit",tempSpecStatPtLinShiftB->GetName()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type,convertToHistograms);
                        if (tempSpecStatPtPol2ShiftA)
                            tempSpecFitPtPol2A              = FitObject(tempSpecStatPtPol2ShiftA,Form("%s_Fit",tempSpecStatPtPol2ShiftA->GetName()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type,convertToHistograms);
                        if (tempSpecStatPtPol2ShiftB)
                            tempSpecFitPtPol2B              = FitObject(tempSpecStatPtPol2ShiftB,Form("%s_Fit",tempSpecStatPtPol2ShiftB->GetName()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type,convertToHistograms);

                        // write fitparameter to file
                        if (tempSpecFit) {
                            textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << " " << fMethod[methodIter] << " - " << type << endl;
                            textoutput << "chi^2 / ndf = " << tempSpecFit->GetChisquare() << " / " << tempSpecFit->GetNDF() << endl;
                            for(Int_t j = 0; j < tempSpecFit->GetNpar(); j++) {
                                if (tempSpecFit->GetParameter(j) != tempSpecFit->GetParameter(j))
                                    textoutput << "par" << j << " ( " << tempSpecFit->GetParName(j) << " ) " << " = nan" << endl;
                                else
                                    textoutput << "par" << j << " ( " << tempSpecFit->GetParName(j) << " ) " << tempSpecFit->GetParameter(j) << " +/- " << tempSpecFit->GetParError(j) << endl;
                            }
                            textoutput << endl;
                        }
                        if (tempSpecFitPtConstUp) {
                            textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << " " << fMethod[methodIter] << " const shift up - " << type << endl;
                            textoutput << "chi^2 / ndf = " << tempSpecFitPtConstUp->GetChisquare() << " / " << tempSpecFitPtConstUp->GetNDF() << endl;
                            for(Int_t j = 0; j < tempSpecFitPtConstUp->GetNpar(); j++) {
                                if (tempSpecFitPtConstUp->GetParameter(j)!=tempSpecFitPtConstUp->GetParameter(j))
                                    textoutput << "par" << j << " ( " << tempSpecFitPtConstUp->GetParName(j) << " ) " << " = nan" << endl;
                                else
                                    textoutput << "par" << j << " ( " << tempSpecFitPtConstUp->GetParName(j) << " ) " << " = " << tempSpecFitPtConstUp->GetParameter(j) << " +/- " << tempSpecFitPtConstUp->GetParError(j) << endl;
                            }
                            textoutput << endl;
                        }
                        if (tempSpecFitPtConstDown) {
                            textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << " " << fMethod[methodIter] << " const shift down - " << type << endl;
                            textoutput << "chi^2 / ndf = " << tempSpecFitPtConstDown->GetChisquare() << " / " << tempSpecFitPtConstDown->GetNDF() << endl;
                            for(Int_t j = 0; j < tempSpecFitPtConstDown->GetNpar(); j++) {
                                if (tempSpecFitPtConstDown->GetParameter(j)!=tempSpecFitPtConstDown->GetParameter(j))
                                    textoutput << "par" << j << " ( " << tempSpecFitPtConstDown->GetParName(j) << " ) " << " = nan" << endl;
                                else
                                    textoutput << "par" << j << " ( " << tempSpecFitPtConstDown->GetParName(j) << " ) " << " = " << tempSpecFitPtConstDown->GetParameter(j) << " +/- " << tempSpecFitPtConstDown->GetParError(j) << endl;
                            }
                            textoutput << endl;
                        }
                        if (tempSpecFitPtLinA) {
                            textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << " " << fMethod[methodIter] << " lin A - " << type << endl;
                            textoutput << "chi^2 / ndf = " << tempSpecFitPtLinA->GetChisquare() << " / " << tempSpecFitPtLinA->GetNDF() << endl;
                            for(Int_t j = 0; j < tempSpecFitPtLinA->GetNpar(); j++) {
                                if (tempSpecFitPtLinA->GetParameter(j)!=tempSpecFitPtLinA->GetParameter(j))
                                    textoutput << "par" << j << " ( " << tempSpecFitPtLinA->GetParName(j) << " ) " << " = nan" << endl;
                                else
                                    textoutput << "par" << j << " ( " << tempSpecFitPtLinA->GetParName(j) << " ) " << " = " << tempSpecFitPtLinA->GetParameter(j) << " +/- " << tempSpecFitPtLinA->GetParError(j) << endl;
                            }
                            textoutput << endl;
                        }
                        if (tempSpecFitPtLinB) {
                            textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << " " << fMethod[methodIter] << " lin B - " << type << endl;
                            textoutput << "chi^2 / ndf = " << tempSpecFitPtLinB->GetChisquare() << " / " << tempSpecFitPtLinB->GetNDF() << endl;
                            for(Int_t j = 0; j < tempSpecFitPtLinB->GetNpar(); j++) {
                                if (tempSpecFitPtLinB->GetParameter(j)!=tempSpecFitPtLinB->GetParameter(j))
                                    textoutput << "par" << j << " ( " << tempSpecFitPtLinB->GetParName(j) << " ) " << " = nan" << endl;
                                else
                                    textoutput << "par" << j << " ( " << tempSpecFitPtLinB->GetParName(j) << " ) " << " = " << tempSpecFitPtLinB->GetParameter(j) << " +/- " << tempSpecFitPtLinB->GetParError(j) << endl;
                            }
                            textoutput << endl;
                        }
                        if (tempSpecFitPtPol2A) {
                            textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << " " << fMethod[methodIter] << " pol2 A - " << type << endl;
                            textoutput << "chi^2 / ndf = " << tempSpecFitPtPol2A->GetChisquare() << " / " << tempSpecFitPtPol2A->GetNDF() << endl;
                            for(Int_t j = 0; j < tempSpecFitPtPol2A->GetNpar(); j++) {
                                if (tempSpecFitPtPol2A->GetParameter(j)!=tempSpecFitPtPol2A->GetParameter(j))
                                    textoutput << "par" << j << " ( " << tempSpecFitPtPol2A->GetParName(j) << " ) " << " = nan" << endl;
                                else
                                    textoutput << "par" << j << " ( " << tempSpecFitPtPol2A->GetParName(j) << " ) " << " = " << tempSpecFitPtPol2A->GetParameter(j) << " +/- " << tempSpecFitPtPol2A->GetParError(j) << endl;
                            }
                            textoutput << endl;
                        }
                        if (tempSpecFitPtPol2B) {
                            textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << " " << fMethod[methodIter] << " pol2 B - " << type << endl;
                            textoutput << "chi^2 / ndf = " << tempSpecFitPtPol2B->GetChisquare() << " / " << tempSpecFitPtPol2B->GetNDF() << endl;
                            for(Int_t j = 0; j < tempSpecFitPtPol2B->GetNpar(); j++) {
                                if (tempSpecFitPtPol2B->GetParameter(j)!=tempSpecFitPtPol2B->GetParameter(j))
                                    textoutput << "par" << j << " ( " << tempSpecFitPtPol2B->GetParName(j) << " ) " << " = nan" << endl;
                                else
                                    textoutput << "par" << j << " ( " << tempSpecFitPtPol2B->GetParName(j) << " ) " << " = " << tempSpecFitPtPol2B->GetParameter(j) << " +/- " << tempSpecFitPtPol2B->GetParError(j) << endl;
                            }
                            textoutput << endl;
                        }

                        // add spectra fits to output list
                        if (tempSpecFit)            outputList[listIter]->Add(tempSpecFit);
                        if (tempSpecFitPtConstUp)   outputList[listIter]->Add(tempSpecFitPtConstUp);
                        if (tempSpecFitPtConstDown) outputList[listIter]->Add(tempSpecFitPtConstDown);
                        if (tempSpecFitPtLinA)      outputList[listIter]->Add(tempSpecFitPtLinA);
                        if (tempSpecFitPtLinB)      outputList[listIter]->Add(tempSpecFitPtLinB);
                        if (tempSpecFitPtPol2A)     outputList[listIter]->Add(tempSpecFitPtPol2A);
                        if (tempSpecFitPtPol2B)     outputList[listIter]->Add(tempSpecFitPtPol2B);

                        // calculate ratios to fit
                        if (tempStatClass.Contains("TH1")) {
                            if (tempSpecFit)
                                tempSpecStatRatioToFit                  = CalculateRatioToFit((TH1D*)tempSpecStat, tempSpecFit);
                            if (tempSpecFitPtConstUp)
                                tempSpecStatPtConstShiftUpRatioToFit    = CalculateRatioToFit((TH1D*)tempSpecStatPtConstShiftUp, tempSpecFitPtConstUp);
                            if (tempSpecFitPtConstDown)
                                tempSpecStatPtConstShiftDownRatioToFit  = CalculateRatioToFit((TH1D*)tempSpecStatPtConstShiftDown, tempSpecFitPtConstDown);
                            if (tempSpecFitPtLinA)
                                tempSpecStatPtLinShiftARatioToFit       = CalculateRatioToFit((TH1D*)tempSpecStatPtLinShiftA, tempSpecFitPtLinA);
                            if (tempSpecFitPtLinB)
                                tempSpecStatPtLinShiftBRatioToFit       = CalculateRatioToFit((TH1D*)tempSpecStatPtLinShiftB, tempSpecFitPtLinB);
                            if (tempSpecFitPtPol2A)
                                tempSpecStatPtPol2ShiftARatioToFit      = CalculateRatioToFit((TH1D*)tempSpecStatPtPol2ShiftA, tempSpecFitPtPol2A);
                            if (tempSpecFitPtPol2B)
                                tempSpecStatPtPol2ShiftBRatioToFit      = CalculateRatioToFit((TH1D*)tempSpecStatPtPol2ShiftB, tempSpecFitPtPol2B);
                        } else if (tempStatClass.CompareTo("TGraphErrors") == 0) {
                            if (tempSpecFit)
                                tempSpecStatRatioToFit                  = CalculateRatioToFit((TGraphErrors*)tempSpecStat, tempSpecFit);
                            if (tempSpecFitPtConstUp)
                                tempSpecStatPtConstShiftUpRatioToFit    = CalculateRatioToFit((TGraphErrors*)tempSpecStatPtConstShiftUp, tempSpecFitPtConstUp);
                            if (tempSpecFitPtConstDown)
                                tempSpecStatPtConstShiftDownRatioToFit  = CalculateRatioToFit((TGraphErrors*)tempSpecStatPtConstShiftDown, tempSpecFitPtConstDown);
                            if (tempSpecFitPtLinA)
                                tempSpecStatPtLinShiftARatioToFit       = CalculateRatioToFit((TGraphErrors*)tempSpecStatPtLinShiftA, tempSpecFitPtLinA);
                            if (tempSpecFitPtLinB)
                                tempSpecStatPtLinShiftBRatioToFit       = CalculateRatioToFit((TGraphErrors*)tempSpecStatPtLinShiftB, tempSpecFitPtLinB);
                            if (tempSpecFitPtPol2A)
                                tempSpecStatPtPol2ShiftARatioToFit      = CalculateRatioToFit((TGraphErrors*)tempSpecStatPtPol2ShiftA, tempSpecFitPtPol2A);
                            if (tempSpecFitPtPol2B)
                                tempSpecStatPtPol2ShiftBRatioToFit      = CalculateRatioToFit((TGraphErrors*)tempSpecStatPtPol2ShiftB, tempSpecFitPtPol2B);
                        } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0) {
                            if (tempSpecFit)
                                tempSpecStatRatioToFit                  = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecStat, tempSpecFit);
                            if (tempSpecFitPtConstUp)
                                tempSpecStatPtConstShiftUpRatioToFit    = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecStatPtConstShiftUp, tempSpecFitPtConstUp);
                            if (tempSpecFitPtConstDown)
                                tempSpecStatPtConstShiftDownRatioToFit  = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecStatPtConstShiftDown, tempSpecFitPtConstDown);
                            if (tempSpecFitPtLinA)
                                tempSpecStatPtLinShiftARatioToFit       = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecStatPtLinShiftA, tempSpecFitPtLinA);
                            if (tempSpecFitPtLinB)
                                tempSpecStatPtLinShiftBRatioToFit       = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecStatPtLinShiftB, tempSpecFitPtLinB);
                            if (tempSpecFitPtPol2A)
                                tempSpecStatPtPol2ShiftARatioToFit      = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecStatPtPol2ShiftA, tempSpecFitPtPol2A);
                            if (tempSpecFitPtPol2B)
                                tempSpecStatPtPol2ShiftBRatioToFit      = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecStatPtPol2ShiftB, tempSpecFitPtPol2B);
                        }

                        // add ratios to fits to output list
                        if (tempSpecStatRatioToFit)                 outputList[listIter]->Add(tempSpecStatRatioToFit);
                        if (tempSpecStatPtConstShiftUpRatioToFit)   outputList[listIter]->Add(tempSpecStatPtConstShiftUpRatioToFit);
                        if (tempSpecStatPtConstShiftDownRatioToFit) outputList[listIter]->Add(tempSpecStatPtConstShiftDownRatioToFit);
                        if (tempSpecStatPtLinShiftARatioToFit)      outputList[listIter]->Add(tempSpecStatPtLinShiftARatioToFit);
                        if (tempSpecStatPtLinShiftBRatioToFit)      outputList[listIter]->Add(tempSpecStatPtLinShiftBRatioToFit);
                        if (tempSpecStatPtPol2ShiftARatioToFit)     outputList[listIter]->Add(tempSpecStatPtPol2ShiftARatioToFit);
                        if (tempSpecStatPtPol2ShiftBRatioToFit)     outputList[listIter]->Add(tempSpecStatPtPol2ShiftBRatioToFit);

                        // fit x range for plots
                        Double_t xMinFit = initXRange[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][0];
                        Double_t xMaxFit = initXRange[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][1];

                        // produce plot to check fit quality
                        TString namePlot                        = "";
                        if (tempStatClass.Contains("TH1")) {

                            namePlot                            = Form("spectra_%s%s_default",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TH1D*)tempSpecStat,NULL,NULL,tempSpecFit,NULL,NULL,tempSpecStatRatioToFit,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptConst",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TH1D*)tempSpecStat,(TH1D*)tempSpecStatPtConstShiftUp,(TH1D*)tempSpecStatPtConstShiftDown,tempSpecFit,tempSpecFitPtConstUp,tempSpecFitPtConstDown,tempSpecStatRatioToFit,tempSpecStatPtConstShiftUpRatioToFit,tempSpecStatPtConstShiftDownRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptLin",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TH1D*)tempSpecStat,(TH1D*)tempSpecStatPtLinShiftA,(TH1D*)tempSpecStatPtLinShiftB,tempSpecFit,tempSpecFitPtLinA,tempSpecFitPtLinB,tempSpecStatRatioToFit,tempSpecStatPtLinShiftARatioToFit,tempSpecStatPtLinShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TH1D*)tempSpecStat,(TH1D*)tempSpecStatPtPol2ShiftA,(TH1D*)tempSpecStatPtPol2ShiftB,tempSpecFit,tempSpecFitPtPol2A,tempSpecFitPtPol2B,tempSpecStatRatioToFit,tempSpecStatPtPol2ShiftARatioToFit,tempSpecStatPtPol2ShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                        } else if (tempStatClass.CompareTo("TGraphErrors") == 0) {
                            namePlot                            = Form("spectra_%s%s_default",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TGraphErrors*)tempSpecStat,NULL,NULL,tempSpecFit,NULL,NULL,tempSpecStatRatioToFit,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptConst",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TGraphErrors*)tempSpecStat,(TGraphErrors*)tempSpecStatPtConstShiftUp,(TGraphErrors*)tempSpecStatPtConstShiftDown,tempSpecFit,tempSpecFitPtConstUp,tempSpecFitPtConstDown,tempSpecStatRatioToFit,tempSpecStatPtConstShiftUpRatioToFit,tempSpecStatPtConstShiftDownRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptLin",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TGraphErrors*)tempSpecStat,(TGraphErrors*)tempSpecStatPtLinShiftA,(TGraphErrors*)tempSpecStatPtLinShiftB,tempSpecFit,tempSpecFitPtLinA,tempSpecFitPtLinB,tempSpecStatRatioToFit,tempSpecStatPtLinShiftARatioToFit,tempSpecStatPtLinShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TGraphErrors*)tempSpecStat,(TGraphErrors*)tempSpecStatPtPol2ShiftA,(TGraphErrors*)tempSpecStatPtPol2ShiftB,tempSpecFit,tempSpecFitPtPol2A,tempSpecFitPtPol2B,tempSpecStatRatioToFit,tempSpecStatPtPol2ShiftARatioToFit,tempSpecStatPtPol2ShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                        } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0) {

                            namePlot                            = Form("spectra_%s%s_default",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempSpecStat,NULL,NULL,tempSpecFit,NULL,NULL,tempSpecStatRatioToFit,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptConst",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempSpecStat,(TGraphAsymmErrors*)tempSpecStatPtConstShiftUp,(TGraphAsymmErrors*)tempSpecStatPtConstShiftDown,tempSpecFit,tempSpecFitPtConstUp,tempSpecFitPtConstDown,tempSpecStatRatioToFit,tempSpecStatPtConstShiftUpRatioToFit,tempSpecStatPtConstShiftDownRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptLin",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempSpecStat,(TGraphAsymmErrors*)tempSpecStatPtLinShiftA,(TGraphAsymmErrors*)tempSpecStatPtLinShiftB,tempSpecFit,tempSpecFitPtLinA,tempSpecFitPtLinB,tempSpecStatRatioToFit,tempSpecStatPtLinShiftARatioToFit,tempSpecStatPtLinShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            namePlot                            = Form("spectra_%s%s_ptPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempSpecStat,(TGraphAsymmErrors*)tempSpecStatPtPol2ShiftA,(TGraphAsymmErrors*)tempSpecStatPtPol2ShiftB,tempSpecFit,tempSpecFitPtPol2A,tempSpecFitPtPol2B,tempSpecStatRatioToFit,tempSpecStatPtPol2ShiftARatioToFit,tempSpecStatPtPol2ShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],type,suffix,"","",namePlot,xMinFit,xMaxFit,kFALSE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);
                        }

                        // plot ratios of syst. fits to central value
                        TString nameParamRatioPlot         = "";
                        if (doParamRatioPlotsSep) {

                            nameParamRatioPlot                 = Form("spectra_%s%s_paramConstToCentral",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotParamRatios(tempSpecFit,tempSpecFitPtConstUp,tempSpecFitPtConstDown,NULL,NULL,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameParamRatioPlot);

                            nameParamRatioPlot                 = Form("spectra_%s%s_paramLinToCentral",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotParamRatios(tempSpecFit,NULL,NULL,tempSpecFitPtLinA,tempSpecFitPtLinB,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameParamRatioPlot);

                            nameParamRatioPlot                 = Form("spectra_%s%s_paramPol2ToCentral",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotParamRatios(tempSpecFit,NULL,NULL,NULL,NULL,tempSpecFitPtPol2A,tempSpecFitPtPol2B,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameParamRatioPlot);

                        } else {

                            nameParamRatioPlot                 = Form("spectra_%s%s_paramAllToCentral",fParticle[particleIter].Data(),fMethod[methodIter].Data());
                            ProducePlotParamRatios(tempSpecFit,tempSpecFitPtConstUp,tempSpecFitPtConstDown,tempSpecFitPtLinA,tempSpecFitPtLinB,tempSpecFitPtPol2A,tempSpecFitPtPol2B,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameParamRatioPlot);
                        }
                    }
                }

                // check if ratio is in input list
                for (Int_t i = 0; i < nParticles; i++) {
                    if ( inputList[listIter]->FindObject(Form("%sTo%s%sStat", fParticle[particleIter].Data(), fParticle[i].Data(), fMethod[methodIter].Data())) && inputList[listIter]->FindObject(Form("%sTo%s%sSys", fParticle[particleIter].Data(), fParticle[i].Data(), fMethod[methodIter].Data())) ) {

                        // for testing purposes
                        cout << "-------- " << fParticle[particleIter].Data() << " to " << fParticle[i].Data() << " " << fMethod[methodIter].Data() << " ratio found" << endl;

                        // skip if ratio was already found and written to output list
                        if (outputList[listIter]->FindObject(Form("%sTo%s%sStat", fParticle[particleIter].Data(), fParticle[i].Data(), fMethod[methodIter].Data()))) continue;

                        // check for pt-const rel. syst. err
                        Double_t ptConstRelSyst                     = initRatioPtConstRelSyst[GetParticleIterator(fParticle[particleIter])][GetParticleIterator(fParticle[i])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])];

                        // get type of sys errs
                        tempStatClass                               = ((TObject*)inputList[listIter]->FindObject(Form("%sTo%s%sStat", fParticle[particleIter].Data(), fParticle[i].Data(), fMethod[methodIter].Data())))->ClassName();
                        tempSysClass                                = ((TObject*)inputList[listIter]->FindObject(Form("%sTo%s%sSys", fParticle[particleIter].Data(), fParticle[i].Data(), fMethod[methodIter].Data())))->ClassName();

                        // get ratio and shift using sys errs
                        tempRatioStat                               = (TObject*)inputList[listIter]->FindObject(Form("%sTo%s%sStat", fParticle[particleIter].Data(), fParticle[i].Data(), fMethod[methodIter].Data()));
                        tempRatioSys                                = (TObject*)inputList[listIter]->FindObject(Form("%sTo%s%sSys", fParticle[particleIter].Data(), fParticle[i].Data(), fMethod[methodIter].Data()));
                        if (tempStatClass.Contains("TH1") && tempSysClass.Contains("TH1")) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else if (tempStatClass.Contains("TH1") && tempSysClass.CompareTo("TGraphErrors") == 0) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else if (tempStatClass.Contains("TH1") && tempSysClass.CompareTo("TGraphAsymmErrors") == 0) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TH1D*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TH1D*)ShiftSpectraWithSyst(      (TH1D*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TH1D*)ShiftSpectraWithSlopeSyst( (TH1D*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else if (tempStatClass.CompareTo("TGraphErrors") == 0 && tempSysClass.Contains("TH1")) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else if (tempStatClass.CompareTo("TGraphErrors") == 0 && tempSysClass.CompareTo("TGraphErrors") == 0) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else if (tempStatClass.CompareTo("TGraphErrors") == 0 && tempSysClass.CompareTo("TGraphAsymmErrors") == 0) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TGraphErrors*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TGraphErrors*)ShiftSpectraWithSyst(      (TGraphErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TGraphErrors*)ShiftSpectraWithSlopeSyst( (TGraphErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0 && tempSysClass.Contains("TH1")) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TH1D*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0 && tempSysClass.CompareTo("TGraphErrors") == 0) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TGraphErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0 && tempSysClass.CompareTo("TGraphAsymmErrors") == 0) {
                            // calculate pt-dep. factors for syst. err.
                            tempRatioFactorLinShiftA                = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 1, kTRUE);
                            tempRatioFactorLinShiftB                = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 1, kFALSE);
                            tempRatioFactorPol2ShiftA               = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 2, kTRUE);
                            tempRatioFactorPol2ShiftB               = CalculateNonConstantShiftFactor((TGraphAsymmErrors*)tempRatioStat, 2, kFALSE);

                            // shift ratios
                            tempRatioStatPtConstShiftUp             = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst);
                            tempRatioStatPtConstShiftDown           = (TGraphAsymmErrors*)ShiftSpectraWithSyst(      (TGraphAsymmErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE, ptConstRelSyst);
                            tempRatioStatPtLinShiftA                = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorLinShiftA, 1);
                            tempRatioStatPtLinShiftB                = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorLinShiftB, 1);
                            tempRatioStatPtPol2ShiftA               = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kTRUE,  ptConstRelSyst, tempRatioFactorPol2ShiftA, 2);
                            tempRatioStatPtPol2ShiftB               = (TGraphAsymmErrors*)ShiftSpectraWithSlopeSyst( (TGraphAsymmErrors*)tempRatioStat, (TGraphAsymmErrors*)tempRatioSys,  kFALSE,  ptConstRelSyst, tempRatioFactorPol2ShiftB, 2);
                        } else {
                            cout << "stat and sys type not knwon, skipping" << endl;
                            continue;
                        }

                        // write ratios to list
                        if (tempRatioStat)                  outputList[listIter]->Add(tempRatioStat);
                        if (tempRatioStatPtConstShiftUp)    outputList[listIter]->Add(tempRatioStatPtConstShiftUp);
                        if (tempRatioStatPtConstShiftDown)  outputList[listIter]->Add(tempRatioStatPtConstShiftDown);
                        if (tempRatioStatPtLinShiftA)       outputList[listIter]->Add(tempRatioStatPtLinShiftA);
                        if (tempRatioStatPtLinShiftB)       outputList[listIter]->Add(tempRatioStatPtLinShiftB);
                        if (tempRatioFactorLinShiftA)       outputList[listIter]->Add(tempRatioFactorLinShiftA);
                        if (tempRatioFactorLinShiftB)       outputList[listIter]->Add(tempRatioFactorLinShiftB);
                        if (tempRatioStatPtPol2ShiftA)      outputList[listIter]->Add(tempRatioStatPtPol2ShiftA);
                        if (tempRatioStatPtPol2ShiftB)      outputList[listIter]->Add(tempRatioStatPtPol2ShiftB);
                        if (tempRatioFactorPol2ShiftA)      outputList[listIter]->Add(tempRatioFactorPol2ShiftA);
                        if (tempRatioFactorPol2ShiftB)      outputList[listIter]->Add(tempRatioFactorPol2ShiftB);

                        // plot pt dept sys err factors
                        TString nameRatioFactorPlot                 = "";
                        if (doPtDepSysErrFactorPlotsSep) {

                            nameRatioFactorPlot                     = Form("ratio_%s%s_%s_ptDepFactorLin",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                            ProducePlotPtDepSystFactors(tempRatioFactorLinShiftA,tempRatioFactorLinShiftB,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameRatioFactorPlot);

                            nameRatioFactorPlot                     = Form("ratio_%s%s_%s_ptDepFactorPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                            ProducePlotPtDepSystFactors(NULL,NULL,tempRatioFactorPol2ShiftA,tempRatioFactorPol2ShiftB,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameRatioFactorPlot);

                        } else {

                            nameRatioFactorPlot                     = Form("ratio_%s%s_%s_ptDepFactorLinPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                            TString labelRatio                      = Form ("%s/%s", fParticleLatex[particleIter].Data(), fParticleLatex[i].Data() );
                            ProducePlotPtDepSystFactors(tempRatioFactorLinShiftA,tempRatioFactorLinShiftB,tempRatioFactorPol2ShiftA,tempRatioFactorPol2ShiftB,fCollSysTemp,fEnergyTemp,fCentralityTemp,labelRatio,fMethodLabel[methodIter],suffix,nameRatioFactorPlot);
                        }

                        // initializing fit function
                        TString typeRatio                           = initRatioFitFunction[GetParticleIterator(fParticle[particleIter])][GetParticleIterator(fParticle[i])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])];
                        if (doRatios && typeRatio.CompareTo("")!=0) {

                            // fit ratios
                            tempRatioFit                            = FitRatio(tempRatioStat,Form("%s_Fit",tempRatioStat->GetName()),fParticle[particleIter],fParticle[i],fCentralityTemp,fMethod[methodIter],typeRatio,convertToHistograms);
                            if (tempRatioStatPtConstShiftUp)
                                tempRatioFitPtConstUp               = FitRatio(tempRatioStatPtConstShiftUp,Form("%s_Fit",tempRatioStatPtConstShiftUp->GetName()),fParticle[particleIter],fParticle[i],fCentralityTemp,fMethod[methodIter],typeRatio,convertToHistograms);
                            if (tempRatioStatPtConstShiftDown)
                                tempRatioFitPtConstDown             = FitRatio(tempRatioStatPtConstShiftDown,Form("%s_Fit",tempRatioStatPtConstShiftDown->GetName()),fParticle[particleIter],fParticle[i],fCentralityTemp,fMethod[methodIter],typeRatio,convertToHistograms);
                            if (tempRatioStatPtLinShiftA)
                                tempRatioFitPtLinA                  = FitRatio(tempRatioStatPtLinShiftA,Form("%s_Fit",tempRatioStatPtLinShiftA->GetName()),fParticle[particleIter],fParticle[i],fCentralityTemp,fMethod[methodIter],typeRatio,convertToHistograms);
                            if (tempRatioStatPtLinShiftB)
                                tempRatioFitPtLinB                  = FitRatio(tempRatioStatPtLinShiftB,Form("%s_Fit",tempRatioStatPtLinShiftB->GetName()),fParticle[particleIter],fParticle[i],fCentralityTemp,fMethod[methodIter],typeRatio,convertToHistograms);
                            if (tempRatioStatPtPol2ShiftA)
                                tempRatioFitPtPol2A                 = FitRatio(tempRatioStatPtPol2ShiftA,Form("%s_Fit",tempRatioStatPtPol2ShiftA->GetName()),fParticle[particleIter],fParticle[i],fCentralityTemp,fMethod[methodIter],typeRatio,convertToHistograms);
                            if (tempRatioStatPtPol2ShiftB)
                                tempRatioFitPtPol2B                 = FitRatio(tempRatioStatPtPol2ShiftB,Form("%s_Fit",tempRatioStatPtPol2ShiftB->GetName()),fParticle[particleIter],fParticle[i],fCentralityTemp,fMethod[methodIter],typeRatio,convertToHistograms);

                            // write fitparameter to file
                            if (tempRatioFit) {
                                textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << "/" << fParticle[i] << " " << fMethod[methodIter].Data() << " - " << typeRatio << endl;
                                textoutput << "chi^2 / ndf = " << tempRatioFit->GetChisquare() << " / " << tempRatioFit->GetNDF() << endl;
                                for(Int_t j = 0; j < tempRatioFit->GetNpar(); j++) {
                                    if (tempRatioFit->GetParameter(j)!=tempRatioFit->GetParameter(j))
                                        textoutput << "par" << j << " ( " << tempRatioFit->GetParName(j) << " ) " << " = nan" << endl;
                                    else
                                        textoutput << "par" << j << " ( " << tempRatioFit->GetParName(j) << " ) " << " = " << tempRatioFit->GetParameter(j) << " +/- " << tempRatioFit->GetParError(j) << endl;
                                }
                                textoutput << endl;
                            }
                            if (tempRatioFitPtConstUp) {
                                textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << "/" << fParticle[i] << " " << fMethod[methodIter].Data() << " const shift up - " << typeRatio << endl;
                                textoutput << "chi^2 / ndf = " << tempRatioFitPtConstUp->GetChisquare() << " / " << tempRatioFitPtConstUp->GetNDF() << endl;
                                for(Int_t j = 0; j < tempRatioFitPtConstUp->GetNpar(); j++) {
                                    if (tempRatioFitPtConstUp->GetParameter(j)!=tempRatioFitPtConstUp->GetParameter(j))
                                        textoutput << "par" << j << " ( " << tempRatioFitPtConstUp->GetParName(j) << " ) " << " = nan" << endl;
                                    else
                                        textoutput << "par" << j << " ( " << tempRatioFitPtConstUp->GetParName(j) << " ) " << " = " << tempRatioFitPtConstUp->GetParameter(j) << " +/- " << tempRatioFitPtConstUp->GetParError(j) << endl;
                                }
                                textoutput << endl;
                            }
                            if (tempRatioFitPtConstDown) {
                                textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << "/" << fParticle[i] << " " << fMethod[methodIter].Data() << " const shift down - " << typeRatio << endl;
                                textoutput << "chi^2 / ndf = " << tempRatioFitPtConstDown->GetChisquare() << " / " << tempRatioFitPtConstDown->GetNDF() << endl;
                                for(Int_t j = 0; j < tempRatioFitPtConstDown->GetNpar(); j++) {
                                    if (tempRatioFitPtConstDown->GetParameter(j)!=tempRatioFitPtConstDown->GetParameter(j))
                                        textoutput << "par" << j << " ( " << tempRatioFitPtConstDown->GetParName(j) << " ) " << " = nan" << endl;
                                    else
                                        textoutput << "par" << j << " ( " << tempRatioFitPtConstDown->GetParName(j) << " ) " << " = " << tempRatioFitPtConstDown->GetParameter(j) << " +/- " << tempRatioFitPtConstDown->GetParError(j) << endl;
                                }
                                textoutput << endl;
                            }
                            if (tempRatioFitPtLinA) {
                                textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << "/" << fParticle[i] << " " << fMethod[methodIter].Data() << " lin A - " << typeRatio << endl;
                                textoutput << "chi^2 / ndf = " << tempRatioFitPtLinA->GetChisquare() << " / " << tempRatioFitPtLinA->GetNDF() << endl;
                                for(Int_t j = 0; j < tempRatioFitPtLinA->GetNpar(); j++) {
                                    if (tempRatioFitPtLinA->GetParameter(j)!=tempRatioFitPtLinA->GetParameter(j))
                                        textoutput << "par" << j << " ( " << tempRatioFitPtLinA->GetParName(j) << " ) " << " = nan" << endl;
                                    else
                                        textoutput << "par" << j << " ( " << tempRatioFitPtLinA->GetParName(j) << " ) " << " = " << tempRatioFitPtLinA->GetParameter(j) << " +/- " << tempRatioFitPtLinA->GetParError(j) << endl;
                                }
                                textoutput << endl;
                            }
                            if (tempRatioFitPtLinB) {
                                textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << "/" << fParticle[i] << " " << fMethod[methodIter].Data() << " lin B - " << typeRatio << endl;
                                textoutput << "chi^2 / ndf = " << tempRatioFitPtLinB->GetChisquare() << " / " << tempRatioFitPtLinB->GetNDF() << endl;
                                for(Int_t j = 0; j < tempRatioFitPtLinB->GetNpar(); j++) {
                                    if (tempRatioFitPtLinB->GetParameter(j)!=tempRatioFitPtLinB->GetParameter(j))
                                        textoutput << "par" << j << " ( " << tempRatioFitPtLinB->GetParName(j) << " ) " << " = nan" << endl;
                                    else
                                        textoutput << "par" << j << " ( " << tempRatioFitPtLinB->GetParName(j) << " ) " << " = " << tempRatioFitPtLinB->GetParameter(j) << " +/- " << tempRatioFitPtLinB->GetParError(j) << endl;
                                }
                                textoutput << endl;
                            }
                            if (tempRatioFitPtPol2A) {
                                textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << "/" << fParticle[i] << " " << fMethod[methodIter].Data() << " pol2 A - " << typeRatio << endl;
                                textoutput << "chi^2 / ndf = " << tempRatioFitPtPol2A->GetChisquare() << " / " << tempRatioFitPtPol2A->GetNDF() << endl;
                                for(Int_t j = 0; j < tempRatioFitPtPol2A->GetNpar(); j++) {
                                    if (tempRatioFitPtPol2A->GetParameter(j)!=tempRatioFitPtPol2A->GetParameter(j))
                                        textoutput << "par" << j << " ( " << tempRatioFitPtPol2A->GetParName(j) << " ) " << " = nan" << endl;
                                    else
                                        textoutput << "par" << j << " ( " << tempRatioFitPtPol2A->GetParName(j) << " ) " << " = " << tempRatioFitPtPol2A->GetParameter(j) << " +/- " << tempRatioFitPtPol2A->GetParError(j) << endl;
                                }
                                textoutput << endl;
                            }
                            if (tempRatioFitPtPol2B) {
                                textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << "/" << fParticle[i] << " " << fMethod[methodIter].Data() << " pol2 B - " << typeRatio << endl;
                                textoutput << "chi^2 / ndf = " << tempRatioFitPtPol2B->GetChisquare() << " / " << tempRatioFitPtPol2B->GetNDF() << endl;
                                for(Int_t j = 0; j < tempRatioFitPtPol2B->GetNpar(); j++) {
                                    if (tempRatioFitPtPol2B->GetParameter(j)!=tempRatioFitPtPol2B->GetParameter(j))
                                        textoutput << "par" << j << " ( " << tempRatioFitPtPol2B->GetParName(j) << " ) " << " = nan" << endl;
                                    else
                                        textoutput << "par" << j << " ( " << tempRatioFitPtPol2B->GetParName(j) << " ) " << " = " << tempRatioFitPtPol2B->GetParameter(j) << " +/- " << tempRatioFitPtPol2B->GetParError(j) << endl;
                                }
                                textoutput << endl;
                            }

                            // add ratio fits to output list
                            if (tempRatioFit)               outputList[listIter]->Add(tempRatioFit);
                            if (tempRatioFitPtConstUp)      outputList[listIter]->Add(tempRatioFitPtConstUp);
                            if (tempRatioFitPtConstDown)    outputList[listIter]->Add(tempRatioFitPtConstDown);
                            if (tempRatioFitPtLinA)         outputList[listIter]->Add(tempRatioFitPtLinA);
                            if (tempRatioFitPtLinB)         outputList[listIter]->Add(tempRatioFitPtLinB);
                            if (tempRatioFitPtPol2A)        outputList[listIter]->Add(tempRatioFitPtPol2A);
                            if (tempRatioFitPtPol2B)        outputList[listIter]->Add(tempRatioFitPtPol2B);

                            // calculate ratios to fit
                            if (tempStatClass.Contains("TH1")) {
                                if (tempRatioFit)
                                    tempRatioStatRatioToFit                 = CalculateRatioToFit((TH1D*)tempRatioStat, tempRatioFit);
                                if (tempRatioFitPtConstUp)
                                    tempRatioStatPtConstShiftUpRatioToFit   = CalculateRatioToFit((TH1D*)tempRatioStatPtConstShiftUp, tempRatioFitPtConstUp);
                                if (tempRatioFitPtConstDown)
                                    tempRatioStatPtConstShiftDownRatioToFit = CalculateRatioToFit((TH1D*)tempRatioStatPtConstShiftDown, tempRatioFitPtConstDown);
                                if (tempRatioFitPtLinA)
                                    tempRatioStatPtLinShiftARatioToFit      = CalculateRatioToFit((TH1D*)tempRatioStatPtLinShiftA, tempRatioFitPtLinA);
                                if (tempRatioFitPtLinB)
                                    tempRatioStatPtLinShiftBRatioToFit      = CalculateRatioToFit((TH1D*)tempRatioStatPtLinShiftB, tempRatioFitPtLinB);
                                if (tempRatioFitPtPol2A)
                                    tempRatioStatPtPol2ShiftARatioToFit     = CalculateRatioToFit((TH1D*)tempRatioStatPtPol2ShiftA, tempRatioFitPtPol2A);
                                if (tempRatioFitPtPol2B)
                                    tempRatioStatPtPol2ShiftBRatioToFit     = CalculateRatioToFit((TH1D*)tempRatioStatPtPol2ShiftB, tempRatioFitPtPol2B);
                            } else if (tempStatClass.CompareTo("TGraphErrors") == 0) {
                                if (tempRatioFit)
                                    tempRatioStatRatioToFit                 = CalculateRatioToFit((TGraphErrors*)tempRatioStat, tempRatioFit);
                                if (tempRatioFitPtConstUp)
                                    tempRatioStatPtConstShiftUpRatioToFit   = CalculateRatioToFit((TGraphErrors*)tempRatioStatPtConstShiftUp, tempRatioFitPtConstUp);
                                if (tempRatioFitPtConstDown)
                                    tempRatioStatPtConstShiftDownRatioToFit = CalculateRatioToFit((TGraphErrors*)tempRatioStatPtConstShiftDown, tempRatioFitPtConstDown);
                                if (tempRatioFitPtLinA)
                                    tempRatioStatPtLinShiftARatioToFit      = CalculateRatioToFit((TGraphErrors*)tempRatioStatPtLinShiftA, tempRatioFitPtLinA);
                                if (tempRatioFitPtLinB)
                                    tempRatioStatPtLinShiftBRatioToFit      = CalculateRatioToFit((TGraphErrors*)tempRatioStatPtLinShiftB, tempRatioFitPtLinB);
                                if (tempRatioFitPtPol2A)
                                    tempRatioStatPtPol2ShiftARatioToFit     = CalculateRatioToFit((TGraphErrors*)tempRatioStatPtPol2ShiftA, tempRatioFitPtPol2A);
                                if (tempRatioFitPtPol2B)
                                    tempRatioStatPtPol2ShiftBRatioToFit     = CalculateRatioToFit((TGraphErrors*)tempRatioStatPtPol2ShiftB, tempRatioFitPtPol2B);
                            } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0) {
                                if (tempRatioFit)
                                    tempRatioStatRatioToFit                 = CalculateRatioToFit((TGraphAsymmErrors*)tempRatioStat, tempRatioFit);
                                if (tempRatioFitPtConstUp)
                                    tempRatioStatPtConstShiftUpRatioToFit   = CalculateRatioToFit((TGraphAsymmErrors*)tempRatioStatPtConstShiftUp, tempRatioFitPtConstUp);
                                if (tempRatioFitPtConstDown)
                                    tempRatioStatPtConstShiftDownRatioToFit = CalculateRatioToFit((TGraphAsymmErrors*)tempRatioStatPtConstShiftDown, tempRatioFitPtConstDown);
                                if (tempRatioFitPtLinA)
                                    tempRatioStatPtLinShiftARatioToFit      = CalculateRatioToFit((TGraphAsymmErrors*)tempRatioStatPtLinShiftA, tempRatioFitPtLinA);
                                if (tempRatioFitPtLinB)
                                    tempRatioStatPtLinShiftBRatioToFit      = CalculateRatioToFit((TGraphAsymmErrors*)tempRatioStatPtLinShiftB, tempRatioFitPtLinB);
                                if (tempRatioFitPtPol2A)
                                    tempRatioStatPtPol2ShiftARatioToFit     = CalculateRatioToFit((TGraphAsymmErrors*)tempRatioStatPtPol2ShiftA, tempRatioFitPtPol2A);
                                if (tempRatioFitPtPol2B)
                                    tempRatioStatPtPol2ShiftBRatioToFit     = CalculateRatioToFit((TGraphAsymmErrors*)tempRatioStatPtPol2ShiftB, tempRatioFitPtPol2B);
                            }

                            // add ratios to fits to output list
                            if (tempRatioStatRatioToFit)                    outputList[listIter]->Add(tempRatioStatRatioToFit);
                            if (tempRatioStatPtConstShiftUpRatioToFit)      outputList[listIter]->Add(tempRatioStatPtConstShiftUpRatioToFit);
                            if (tempRatioStatPtConstShiftDownRatioToFit)    outputList[listIter]->Add(tempRatioStatPtConstShiftDownRatioToFit);
                            if (tempRatioStatPtLinShiftARatioToFit)         outputList[listIter]->Add(tempRatioStatPtLinShiftARatioToFit);
                            if (tempRatioStatPtLinShiftBRatioToFit)         outputList[listIter]->Add(tempRatioStatPtLinShiftBRatioToFit);
                            if (tempRatioStatPtPol2ShiftARatioToFit)        outputList[listIter]->Add(tempRatioStatPtPol2ShiftARatioToFit);
                            if (tempRatioStatPtPol2ShiftBRatioToFit)        outputList[listIter]->Add(tempRatioStatPtPol2ShiftBRatioToFit);

                            // fit x range for plots
                            Double_t xMinFit = initRatioXRange[GetParticleIterator(fParticle[particleIter])][GetParticleIterator(fParticle[i])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][0];
                            Double_t xMaxFit = initRatioXRange[GetParticleIterator(fParticle[particleIter])][GetParticleIterator(fParticle[i])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][1];

                            // produce plot to check fit quality
                            TString namePlot                        = "";
                            if (tempStatClass.Contains("TH1")) {
                                namePlot                            = Form("ratio_%s%s_%s_default",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TH1D*)tempRatioStat,NULL,NULL,tempRatioFit,NULL,NULL,tempRatioStatRatioToFit,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptConst",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TH1D*)tempRatioStat,(TH1D*)tempRatioStatPtConstShiftUp,(TH1D*)tempRatioStatPtConstShiftDown,tempRatioFit,tempRatioFitPtConstUp,tempRatioFitPtConstDown,tempRatioStatRatioToFit,tempRatioStatPtConstShiftUpRatioToFit,tempRatioStatPtConstShiftDownRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptLin",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TH1D*)tempRatioStat,(TH1D*)tempRatioStatPtLinShiftA,(TH1D*)tempRatioStatPtLinShiftB,tempRatioFit,tempRatioFitPtLinA,tempRatioFitPtLinB,tempRatioStatRatioToFit,tempRatioStatPtLinShiftARatioToFit,tempRatioStatPtLinShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TH1D*)tempRatioStat,(TH1D*)tempRatioStatPtPol2ShiftA,(TH1D*)tempRatioStatPtPol2ShiftB,tempRatioFit,tempRatioFitPtPol2A,tempRatioFitPtPol2B,tempRatioStatRatioToFit,tempRatioStatPtPol2ShiftARatioToFit,tempRatioStatPtPol2ShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            } else if (tempStatClass.CompareTo("TGraphErrors") == 0) {
                                namePlot                            = Form("ratio_%s%s_%s_default",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TGraphErrors*)tempRatioStat,NULL,NULL,tempRatioFit,NULL,NULL,tempRatioStatRatioToFit,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptConst",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TGraphErrors*)tempRatioStat,(TGraphErrors*)tempRatioStatPtConstShiftUp,(TGraphErrors*)tempRatioStatPtConstShiftDown,tempRatioFit,tempRatioFitPtConstUp,tempRatioFitPtConstDown,tempRatioStatRatioToFit,tempRatioStatPtConstShiftUpRatioToFit,tempRatioStatPtConstShiftDownRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptLin",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TGraphErrors*)tempRatioStat,(TGraphErrors*)tempRatioStatPtLinShiftA,(TGraphErrors*)tempRatioStatPtLinShiftB,tempRatioFit,tempRatioFitPtLinA,tempRatioFitPtLinB,tempRatioStatRatioToFit,tempRatioStatPtLinShiftARatioToFit,tempRatioStatPtLinShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TGraphErrors*)tempRatioStat,(TGraphErrors*)tempRatioStatPtPol2ShiftA,(TGraphErrors*)tempRatioStatPtPol2ShiftB,tempRatioFit,tempRatioFitPtPol2A,tempRatioFitPtPol2B,tempRatioStatRatioToFit,tempRatioStatPtPol2ShiftARatioToFit,tempRatioStatPtPol2ShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                            } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0) {
                                namePlot                            = Form("ratio_%s%s_%s_default",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempRatioStat,NULL,NULL,tempRatioFit,NULL,NULL,tempRatioStatRatioToFit,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptConst",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempRatioStat,(TGraphAsymmErrors*)tempRatioStatPtConstShiftUp,(TGraphAsymmErrors*)tempRatioStatPtConstShiftDown,tempRatioFit,tempRatioFitPtConstUp,tempRatioFitPtConstDown,tempRatioStatRatioToFit,tempRatioStatPtConstShiftUpRatioToFit,tempRatioStatPtConstShiftDownRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptLin",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempRatioStat,(TGraphAsymmErrors*)tempRatioStatPtLinShiftA,(TGraphAsymmErrors*)tempRatioStatPtLinShiftB,tempRatioFit,tempRatioFitPtLinA,tempRatioFitPtLinB,tempRatioStatRatioToFit,tempRatioStatPtLinShiftARatioToFit,tempRatioStatPtLinShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);

                                namePlot                            = Form("ratio_%s%s_%s_ptPol2",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempRatioStat,(TGraphAsymmErrors*)tempRatioStatPtPol2ShiftA,(TGraphAsymmErrors*)tempRatioStatPtPol2ShiftB,tempRatioFit,tempRatioFitPtPol2A,tempRatioFitPtPol2B,tempRatioStatRatioToFit,tempRatioStatPtPol2ShiftARatioToFit,tempRatioStatPtPol2ShiftBRatioToFit,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],typeRatio,suffix,"","",namePlot,xMinFit,xMaxFit,kTRUE,kFALSE,-1000,-1000,0.5,1.5,"",plotXaxisRatios);
                            }

                            // plot ratios of syst. fits to central value
                            TString nameRatioParamRatioPlot             = "";
                            if (doParamRatioPlotsSep) {

                                nameRatioParamRatioPlot                 = Form("ratio_%s%s_%s_paramConstToCentral",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotParamRatios(tempRatioFit,tempRatioFitPtConstUp,tempRatioFitPtConstDown,NULL,NULL,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameRatioParamRatioPlot);

                                nameRatioParamRatioPlot                 = Form("ratio_%s%s_%s_paramLinToCentral",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotParamRatios(tempRatioFit,NULL,NULL,tempRatioFitPtLinA,tempRatioFitPtLinB,NULL,NULL,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameRatioParamRatioPlot);

                                nameRatioParamRatioPlot                 = Form("ratio_%s%s_%s_paramPol2ToCentral",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotParamRatios(tempRatioFit,NULL,NULL,NULL,NULL,tempRatioFitPtPol2A,tempRatioFitPtPol2B,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameRatioParamRatioPlot);

                            } else {

                                nameRatioParamRatioPlot                 = Form("ratio_%s%s_%s_paramAllToCentral",fParticle[particleIter].Data(),fMethod[methodIter].Data(),fParticle[i].Data());
                                ProducePlotParamRatios(tempRatioFit,tempRatioFitPtConstUp,tempRatioFitPtConstDown,tempRatioFitPtLinA,tempRatioFitPtLinB,tempRatioFitPtPol2A,tempRatioFitPtPol2B,fCollSysTemp,fEnergyTemp,fCentralityTemp,fParticle[particleIter],fMethodLabel[methodIter],suffix,nameRatioParamRatioPlot);
                            }
                        }
                    }
                }
            }
        }

        listIter++;
    }

    // write parametrizations file
    for (Int_t i=0; i<listIter; i++) {
        for (Int_t k = 0; k<(Int_t)vecUniqueNames.size();k++){
            for(Int_t j=0; j<9; j++) {
                if (j == 0){
                    if (outputList[i]) WriteParametrizationsFile(outputList[i], inputList[i], ((TString)vecUniqueNames.at(k)).Data(), j, suffix);
                } else if ((Int_t)vecUniqueNamesEnableSys.at(k) > 0 && j > 1 ){
                    if (outputList[i]) WriteParametrizationsFile(outputList[i], inputList[i], ((TString)vecUniqueNames.at(k)).Data(), j, suffix);
                }
            }
        }
    }

    // write output lists to file
    cout << "---- writing to " << outputFile->GetName() << endl;
    outputFile->cd();
    for (Int_t i=0; i<listIter; i++) {
        if (outputList[i]) outputList[i]->Write(outputList[i]->GetName(), TObject::kSingleKey);
    }
    outputFile->Close();
    textoutput.close();

    // free pointer
    delete inputFile;
    delete outputFile;
    delete[] inputList;
    delete[] outputList;
}
