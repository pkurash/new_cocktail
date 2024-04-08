/*****************************************************************************************************************************
 ******    Friederike Bock, friederike.bock@cern.ch                                                                      *****
 *****************************************************************************************************************************/
#include "CocktailFunctions.h"
#include "CocktailPlotting.h"
#include <iostream>

void CalculateIntegratedYield(  TString inputFileName                       = "CocktailInputPP.root",
                                TString energy                              = "",
                                TString suffix                              = "eps",
                                TString parametrizationSettingsFile         = "",
                                Bool_t thesisPlots                          = kFALSE
                             ) {

    if (thesisPlots){
        isThesis                                            = kTRUE;
        cout << "setting thesis mode" << endl;
    }

    // input file
    TFile* inputFile                                        = NULL;
    if (inputFileName.CompareTo("") == 0) {
        std::cout << "ERROR: No inp ut file specified, returning!" << std::endl;
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
    Bool_t doSpectra;
    if (parametrizationSettingsFile.CompareTo("") != 0) {
        if (InitializeFittingIntegYield(parametrizationSettingsFile)) {
            std::cout << "Will parametrize spectra." << std::endl;
            doSpectra                                       = kTRUE;
        } else {
            doSpectra                                       = kFALSE;
        }
    } else {
        std::cout << "WARNING: No parametrization settings file specified, skipping parametrization of spectra!" << std::endl;
        doSpectra                                           = kFALSE;
    }

    TString energyOut                                       = energy.Data();
    energyOut.ReplaceAll(".","_");
    // output file
    TString outputFileName                                  = inputFileName.ReplaceAll(".root", 5, "_IntegYield_", 12);
    outputFileName                                          = Form("%s%s.root", outputFileName.Data(), energyOut.Data());
    TFile *outputFile                                       = new TFile(outputFileName.Data(), "UPDATE");

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
    TF1*        tempSpecFit                                 = NULL;
    TH1D*       tempSpecStatRatioToFit                      = NULL;
    TH1D*       tempSpecSysRatioToFit                       = NULL;
    TF1* tempSysFits[4]                                     = {NULL, NULL, NULL, NULL};
    TString tempFitNames[4]                                 = {"maxSys", "minSys", "hardest", "softest"};

    // loop over input lists
    TIter keyIter(inputFile->GetListOfKeys());
    TKey* key                                               = NULL;
    Int_t listIter                                          = 0;
    while ( (key = (TKey*)keyIter()) ) {

        // input list
        inputList[listIter]                                 = (TList*)key->ReadObj();

        // check for energy and set output list
        if (((TString)inputList[listIter]->GetName()).Contains(energy)) {
            outputList[listIter]                             = (TList*)outputFile->Get(inputList[listIter]->GetName());
            if (!outputList[listIter]){
                outputList[listIter]                         = new TList();
                outputList[listIter]->SetName(inputList[listIter]->GetName());
            }
        } else {
            outputList[listIter]                            = NULL;
            continue;
        }

        // create directory for plots
        TString outputDir                                   = Form("plots/%s/%s/IntegYields/base",suffix.Data(),inputList[listIter]->GetName());
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
            for (Int_t methodIter = 0; methodIter < 2; methodIter++) {

                // check if stat&sys spectra for particle are contained in input list
                if ( inputList[listIter]->FindObject(Form("%s%sStat", fParticle[particleIter].Data(), fMethod[methodIter].Data())) &&
                     inputList[listIter]->FindObject(Form("%s%sSys", fParticle[particleIter].Data(), fMethod[methodIter].Data())) ) {

                    // for testing purposes
                    cout << "-------- " << fParticle[particleIter].Data() << " " << fMethod[methodIter].Data() << " spectra found" << endl;

                    // get type stat&sys
                    tempStatClass                           = ((TObject*)inputList[listIter]->FindObject(Form("%s%sStat", fParticle[particleIter].Data(), fMethod[methodIter].Data())))->ClassName();
                    tempSysClass                            = ((TObject*)inputList[listIter]->FindObject(Form("%s%sSys",  fParticle[particleIter].Data(), fMethod[methodIter].Data())))->ClassName();

                    TList* currentParticleList              = new TList();

                    // get spectra with stat and sys errs
                    tempSpecStat                            = (TObject*)inputList[listIter]->FindObject(Form("%s%sStat", fParticle[particleIter].Data(), fMethod[methodIter].Data()));
                    tempSpecSys                             = (TObject*)inputList[listIter]->FindObject(Form("%s%sSys",  fParticle[particleIter].Data(), fMethod[methodIter].Data()));
                    // initializing fit function
                    TString type                            = initFitFunction[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])];
                    TString namePlotStat                    = Form("%s/SpectraWithFitStat_%s%s_%s.%s", outputDir.Data(), fParticle[particleIter].Data(), fMethod[methodIter].Data(), type.Data(), suffix.Data());
                    TString namePlotSys                     = Form("%s/SpectraWithFitSys_%s%s_%s.%s", outputDir.Data(), fParticle[particleIter].Data(), fMethod[methodIter].Data(), type.Data(), suffix.Data());


                    if (doSpectra && type.CompareTo("")!=0) {
                        // fit spectra
                        if (initDataToBeFitted[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])].CompareTo("Stat") == 0){
                            cout << "fitting stat errors" <<endl;
                            tempSpecFit                         = FitObject(tempSpecStat,Form("%s_Fit_%s",tempSpecStat->GetName(),type.Data()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type);
                        } else{
                            cout << "fitting sys errors" <<endl;
                            tempSpecFit                         = FitObject(tempSpecSys,Form("%s_Fit_%s",tempSpecStat->GetName(),type.Data()),fParticle[particleIter],fCentralityTemp,fMethod[methodIter],type);
                        }
                        for (Int_t k = 0; k< 4; k++)
                            tempSysFits[k]                  = NULL;

                        // write fitparameter to file
                        if (tempSpecFit) {
                            tempSpecFit->SetNpx(10000);
                            Double_t chi2                   = tempSpecFit->GetChisquare()/tempSpecFit->GetNDF();
//                             if (chi2 < 20){
                                textoutput << inputList[listIter]->GetName() << " " << fParticle[particleIter] << " " << fMethod[methodIter] << " - " << type << endl;
                                textoutput << "chi^2 / ndf = " << tempSpecFit->GetChisquare() << " / " << tempSpecFit->GetNDF() << endl;
                                for(Int_t j = 0; j < tempSpecFit->GetNpar(); j++) {
                                    if (tempSpecFit->GetParameter(j) != tempSpecFit->GetParameter(j))
                                        textoutput << "par" << j << " ( " << tempSpecFit->GetParName(j) << " ) " << " = nan" << endl;
                                    else
                                        textoutput << "par" << j << " ( " << tempSpecFit->GetParName(j) << " ) " << tempSpecFit->GetParameter(j) << " +/- " << tempSpecFit->GetParError(j) << endl;
                                }
                                textoutput << endl;

                                cout << initXRange[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][0] << "\t" << initXRange[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][1] << endl;
                                Double_t ptMinMeas              = initXRangeData[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][0];
                                Double_t ptMaxMeas              = initXRangeData[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][1];
                                Double_t ptMinInteg             = initXRangeInteg[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][0];
                                Double_t ptMaxInteg             = initXRangeInteg[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][1];
                                TH1D* histIntYield              = (TH1D*)YieldMean( tempSpecStat, tempSpecSys, tempSpecFit,
                                                                                    ptMinMeas, ptMaxMeas, 0.01, 0.1, "0q", "log.out", ptMinInteg, ptMaxInteg,
                                                                                    namePlotStat.Data(), namePlotSys.Data(),
                                                                                    Form("Yield_%s_%s",fParticle[particleIter].Data(), type.Data()), tempSysFits);
                                cout << "trying to add Yield hist" << endl;
                                if (histIntYield)               outputList[listIter]->Add(histIntYield);
                                cout << "done" << endl;

//                             }
                        }

                        // add spectra fits to output list


                        // calculate ratios to fit
                        if (tempStatClass.Contains("TH1")) {
                            // write spectra to list
                            ((TH1*)tempSpecStat)->SetName(Form("%sStat", fParticle[particleIter].Data()));
                            outputList[listIter]->Add(tempSpecStat);
                            currentParticleList->Add(tempSpecStat);
                            if (tempSpecFit){
                                tempSpecStatRatioToFit                  = CalculateRatioToFit((TH1D*)tempSpecStat, tempSpecFit);
                                tempSpecFit->SetName(Form("%s_Fit_%s",tempSpecStat->GetName(),type.Data()));
                                outputList[listIter]->Add(tempSpecFit);
                            }
                        } else if (tempStatClass.CompareTo("TGraphErrors") == 0) {
                            ((TGraphErrors*)tempSpecStat)->SetName(Form("%sStat", fParticle[particleIter].Data()));
                            outputList[listIter]->Add(tempSpecStat);
                            currentParticleList->Add(tempSpecStat);
                            if (tempSpecFit){
                                tempSpecStatRatioToFit                  = CalculateRatioToFit((TGraphErrors*)tempSpecStat, tempSpecFit);
                                tempSpecFit->SetName(Form("%s_Fit_%s",tempSpecStat->GetName(),type.Data()));
                                outputList[listIter]->Add(tempSpecFit);
                                currentParticleList->Add(tempSpecFit);
                            }
                        } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0) {
                            ((TGraphAsymmErrors*)tempSpecStat)->SetName(Form("%sStat", fParticle[particleIter].Data()));
                            outputList[listIter]->Add(tempSpecStat);
                            currentParticleList->Add(tempSpecStat);
                            if (tempSpecFit){
                                tempSpecStatRatioToFit                  = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecStat, tempSpecFit);
                                tempSpecFit->SetName(Form("%s_Fit_%s",tempSpecStat->GetName(),type.Data()));
                                outputList[listIter]->Add(tempSpecFit);
                                currentParticleList->Add(tempSpecFit);
                            }
                        }
                        // calculate ratios to fit
                        if (tempSysClass.Contains("TH1")) {
                            ((TH1*)tempSpecSys)->SetName(Form("%sSys", fParticle[particleIter].Data()));
                            outputList[listIter]->Add(tempSpecSys);
                            currentParticleList->Add(tempSpecSys);
                            if (tempSpecFit)
                                tempSpecSysRatioToFit                   = CalculateRatioToFit((TH1D*)tempSpecSys, tempSpecFit);
                        } else if (tempSysClass.CompareTo("TGraphErrors") == 0) {
                            ((TGraphErrors*)tempSpecSys)->SetName(Form("%sSys", fParticle[particleIter].Data()));
                            outputList[listIter]->Add(tempSpecSys);
                            currentParticleList->Add(tempSpecSys);
                            if (tempSpecFit)
                                tempSpecSysRatioToFit                   = CalculateRatioToFit((TGraphErrors*)tempSpecSys, tempSpecFit);
                        } else if (tempSysClass.CompareTo("TGraphAsymmErrors") == 0) {
                            ((TGraphAsymmErrors*)tempSpecSys)->SetName(Form("%sSys", fParticle[particleIter].Data()));
                            outputList[listIter]->Add(tempSpecSys);
                            currentParticleList->Add(tempSpecSys);
                            if (tempSpecFit)
                                tempSpecSysRatioToFit                   = CalculateRatioToFit((TGraphAsymmErrors*)tempSpecSys, tempSpecFit);
                        }

                        // add ratios to fits to output list
                        if (tempSpecStatRatioToFit)                 outputList[listIter]->Add(tempSpecStatRatioToFit);
                        if (tempSpecSysRatioToFit)                  outputList[listIter]->Add(tempSpecSysRatioToFit);
                        // fit x range for plots
                        Double_t xMinFit = initXRange[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][0];
                        Double_t xMaxFit = initXRange[GetParticleIterator(fParticle[particleIter])][GetCentralityIterator(fCentralityTemp)][GetMethodIterator(fMethod[methodIter])][1];

                        // produce plot to check fit quality
                        TString namePlot                        = "";
                        if (tempStatClass.Contains("TH1") && tempSysClass.Contains("TH1")) {
                            namePlot                            = Form("spectra_%s%s_DefaultPlusSys_%s",fParticle[particleIter].Data(),fMethod[methodIter].Data(),type.Data());
                            ProducePlotWithRatioToFit((TH1D*)tempSpecStat, (TH1D*)tempSpecSys, NULL,
                                                      tempSpecFit, NULL, NULL,
                                                      tempSpecStatRatioToFit, tempSpecSysRatioToFit, NULL,
                                                      fCollSysTemp, fEnergyTemp, fCentralityTemp, fParticle[particleIter], fMethod[methodIter], type, suffix, "", "", namePlot, xMinFit, xMaxFit, kFALSE, kTRUE, -1000, -1000, 0.7, 1.3, outputDir);
                            namePlot                            = Form("spectra_%s%s_DefaultPlusSysZoomLowPt_%s",fParticle[particleIter].Data(),fMethod[methodIter].Data(),type.Data());
                            ProducePlotWithRatioToFit((TH1D*)tempSpecStat, (TH1D*)tempSpecSys, NULL,
                                                      tempSpecFit, NULL, NULL,
                                                      tempSpecStatRatioToFit, tempSpecSysRatioToFit, NULL,
                                                      fCollSysTemp, fEnergyTemp, fCentralityTemp, fParticle[particleIter], fMethod[methodIter], type, suffix, "", "", namePlot, xMinFit, xMaxFit, kFALSE, kTRUE, 0, 2., 0.8, 1.2, outputDir);
                        } else if (tempStatClass.CompareTo("TGraphErrors") == 0) {
                            namePlot                            = Form("spectra_%s%s_DefaultPlusSys_%s",fParticle[particleIter].Data(),fMethod[methodIter].Data(),type.Data());
                            ProducePlotWithRatioToFit((TGraphErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, NULL,
                                                      tempSpecFit, NULL, NULL,
                                                      tempSpecStatRatioToFit, tempSpecSysRatioToFit, NULL,
                                                      fCollSysTemp, fEnergyTemp, fCentralityTemp, fParticle[particleIter], fMethod[methodIter], type, suffix, "", "", namePlot, xMinFit, xMaxFit, kFALSE, kTRUE, -1000, -1000, 0.7, 1.3, outputDir);
                            namePlot                            = Form("spectra_%s%s_DefaultPlusSysZoomLowPt_%s",fParticle[particleIter].Data(),fMethod[methodIter].Data(),type.Data());
                            ProducePlotWithRatioToFit((TGraphErrors*)tempSpecStat, (TGraphErrors*)tempSpecSys, NULL,
                                                      tempSpecFit, NULL, NULL,
                                                      tempSpecStatRatioToFit, tempSpecSysRatioToFit, NULL,
                                                      fCollSysTemp, fEnergyTemp, fCentralityTemp, fParticle[particleIter], fMethod[methodIter], type, suffix, "", "", namePlot, xMinFit, xMaxFit, kFALSE, kTRUE, 0, 2., 0.8, 1.2, outputDir);
                        } else if (tempStatClass.CompareTo("TGraphAsymmErrors") == 0) {
                            namePlot                            = Form("spectra_%s%s_DefaultPlusSys_%s",fParticle[particleIter].Data(),fMethod[methodIter].Data(),type.Data());
                            ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, NULL,
                                                       tempSpecFit, NULL, NULL,
                                                       tempSpecStatRatioToFit, tempSpecSysRatioToFit, NULL,
                                                       fCollSysTemp, fEnergyTemp, fCentralityTemp, fParticle[particleIter], fMethod[methodIter], type, suffix, "", "", namePlot, xMinFit, xMaxFit, kFALSE, kTRUE, -1000, -1000, 0.7, 1.3, outputDir);
                            namePlot                            = Form("spectra_%s%s_DefaultPlusSysZoomLowPt_%s",fParticle[particleIter].Data(),fMethod[methodIter].Data(),type.Data());
                            ProducePlotWithRatioToFit((TGraphAsymmErrors*)tempSpecStat, (TGraphAsymmErrors*)tempSpecSys, NULL,
                                                      tempSpecFit, NULL, NULL,
                                                      tempSpecStatRatioToFit, tempSpecSysRatioToFit, NULL,
                                                      fCollSysTemp, fEnergyTemp, fCentralityTemp, fParticle[particleIter], fMethod[methodIter], type, suffix, "", "", namePlot, xMinFit, xMaxFit, kFALSE, kTRUE, 0, 2., 0.8, 1.2, outputDir);
                        }

                        cout << "checking if fits for sys converged" << endl;

                        for (Int_t k = 0; k< 4; k++){
                            cout << k << "\t" << tempSysFits[k] << endl;
                            if (tempSysFits[k]){
                                tempSysFits[k]->SetName(Form("%sStat_Fit_%s", fParticle[particleIter].Data(), tempFitNames[k].Data()));
                                currentParticleList->Add(tempSysFits[k]);
                            }
                        }

                        ProduceSingleParticleSpectraPlotFromListOnlyFinalWithFits(currentParticleList, fCollSysTemp, fEnergyTemp, fCentralityTemp, suffix, outputDir, fParticle[particleIter]);
                        ProduceSingleParticleSpectraPlotFromListOnlyFinalWithFits(currentParticleList, fCollSysTemp, fEnergyTemp, fCentralityTemp, suffix, outputDir, fParticle[particleIter], kFALSE);
                        ProduceSingleParticleSpectraPlotFromListOnlyFinalWithFits(currentParticleList, fCollSysTemp, fEnergyTemp, fCentralityTemp, suffix, outputDir, fParticle[particleIter], kFALSE, 2.);

                        delete currentParticleList;
                    }
                }
            }
        }

        listIter++;
    }

    // write output lists to file
    cout << "---- writing to " << outputFile->GetName() << endl;
    outputFile->cd();
    for (Int_t i=0; i<listIter; i++) {
        if (outputList[i]) outputList[i]->Write(outputList[i]->GetName(), TObject::kOverwrite+TObject::kSingleKey);
    }

    outputFile->Close();
    textoutput.close();

    // free pointer
    delete inputFile;
    delete outputFile;
    delete[] inputList;
    delete[] outputList;
}
