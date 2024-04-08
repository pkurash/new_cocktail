/*****************************************************************************************************************************
 ******    Friederike Bock, friederike.bock@cern.ch                                                                      *****
 *****************************************************************************************************************************/
#include "CocktailFunctions.h"
#include "CocktailPlotting.h"
#include <iostream>

void PlotIntegratedYieldsAndMeanPt( TString inputFileName                       = "CocktailInputPP.root",
                                    TString energy                              = "",
                                    TString centrality                          = "",
                                    TString suffix                              = "eps",
                                    TString parametrizationSettingsFile         = "",
                                    Bool_t thesisPlots                          = kFALSE
                                   ) {

    if (thesisPlots){
        isThesis                                            = kTRUE;
        cout << "set thesis mode" << endl;
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
        doSpectra                                           = InitializeFinalCalcIntegYield(parametrizationSettingsFile);
    } else {
        std::cout << "WARNING: No settings file specified, stopping here!" << std::endl;
    }

    cout << "here" << endl;
    TString energyOut                                       = energy.Data();
    energyOut.ReplaceAll(".","_");
    // output file
    TString outputFileName                                  = inputFileName.ReplaceAll("IntegYield", 10, "IntegYieldComp", 16);
    outputFileName                                          = Form("%s%s.root", outputFileName.Data(), energyOut.Data());
    TFile *outputFile                                       = new TFile(outputFileName.Data(), "UPDATE");

    TString textOutputFileName                              = outputFileName.ReplaceAll(".root", 5, ".dat", 4);
    ofstream textoutput;
    textoutput.open(textOutputFileName);

    // declare input and output lists
    Int_t numberOfLists                                     = (Int_t)((TList*)inputFile->GetListOfKeys())->GetEntries();
    TList** inputList                                       = new TList*[numberOfLists];
    TList* outputList                                       = new TList();

    TString tempStatClass                                   = "";
    TString tempSysClass                                    = "";

    // declare temp spectrum histograms, graphs and fits
    TObject*    tempSpecStat                                = NULL;
    TObject*    tempSpecSys                                 = NULL;
    TF1*        tempSpecFit                                 = NULL;
    TH1D*       tempSpecStatRatioToFit                      = NULL;
    TH1D*       tempSpecSysRatioToFit                       = NULL;

    TH1D*       particleYieldStat                           = new TH1D("IntegratedYieldStat", "", nParticles, 0.5, nParticles+0.5);
    particleYieldStat->GetYaxis()->SetTitle("dN/dy");
    TH1D*       particleYieldSys                            = new TH1D("IntegratedYieldSys", "", nParticles, 0.5, nParticles+0.5);
    particleYieldSys->GetYaxis()->SetTitle("dN/dy");
    TH1D*       particleYieldSysFunc                        = new TH1D("IntegratedYieldSysFunc", "", nParticles, 0.5, nParticles+0.5);
    particleYieldSysFunc->GetYaxis()->SetTitle("dN/dy");
    TH1D*       meanPtStat                                  = new TH1D("meanPtStat", "", nParticles, 0.5, nParticles+0.5);
    meanPtStat->GetYaxis()->SetTitle("<p_{T}>");
    TH1D*       meanPtSys                                   = new TH1D("meanPtSys", "", nParticles, 0.5, nParticles+0.5);
    meanPtSys->GetYaxis()->SetTitle("<p_{T}>");
    TH1D*       meanPtSysFunc                               = new TH1D("meanPtSysFunc", "", nParticles, 0.5, nParticles+0.5);
    meanPtSysFunc->GetYaxis()->SetTitle("<p_{T}>");
    TH1D*       particleYieldExtra                          = new TH1D("IntegratedYieldExtrapolated", "", nParticles, 0.5, nParticles+0.5);
    particleYieldExtra->GetYaxis()->SetTitle("dN/dy");
    for (Int_t particleIter = 0; particleIter < nParticles; particleIter++){
        particleYieldStat->GetXaxis()->SetBinLabel(particleIter+1,fParticle[particleIter]);
        particleYieldSys->GetXaxis()->SetBinLabel(particleIter+1,fParticle[particleIter]);
        particleYieldSysFunc->GetXaxis()->SetBinLabel(particleIter+1,fParticle[particleIter]);
        meanPtStat->GetXaxis()->SetBinLabel(particleIter+1,fParticle[particleIter]);
        meanPtSys->GetXaxis()->SetBinLabel(particleIter+1,fParticle[particleIter]);
        meanPtSysFunc->GetXaxis()->SetBinLabel(particleIter+1,fParticle[particleIter]);
        particleYieldExtra->GetXaxis()->SetBinLabel(particleIter+1,fParticle[particleIter]);
    }

    TGraphAsymmErrors* graphParticleYieldStat[3]            = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphParticleYieldSys[3]             = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphParticleYieldSysFunc[3]         = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphMeanPtStat[3]                   = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphMeanPtSys[3]                    = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphMeanPtSysFunc[3]                = {NULL, NULL, NULL};
    TString addNameGraphs[3]                                = {"", "Meson", "Baryon"};

    for (Int_t i = 0; i < 3; i++){
        graphParticleYieldStat[i]                           = new TGraphAsymmErrors(nParticles);
        graphParticleYieldStat[i]->SetName(Form("IntegratedYieldStatVsMass%s", addNameGraphs[i].Data()));
        graphParticleYieldSys[i]                            = new TGraphAsymmErrors(nParticles);
        graphParticleYieldSys[i]->SetName(Form("IntegratedYieldSysVsMass%s", addNameGraphs[i].Data()));
        graphParticleYieldSysFunc[i]                        = new TGraphAsymmErrors(nParticles);
        graphParticleYieldSysFunc[i]->SetName(Form("IntegratedYieldSysFuncVsMass%s", addNameGraphs[i].Data()));
        graphMeanPtStat[i]                                  = new TGraphAsymmErrors(nParticles);
        graphMeanPtStat[i]->SetName(Form("meanPtStatVsMass%s", addNameGraphs[i].Data()));
        graphMeanPtSys[i]                                   = new TGraphAsymmErrors(nParticles);
        graphMeanPtSys[i]->SetName(Form("meanPtSysVsMass%s", addNameGraphs[i].Data()));
        graphMeanPtSysFunc[i]                               = new TGraphAsymmErrors(nParticles);
        graphMeanPtSysFunc[i]->SetName(Form("meanPtSysFuncVsMass%s", addNameGraphs[i].Data()));
    }

    TGraphAsymmErrors* graphExtraPolationFrac               = new TGraphAsymmErrors(nParticles);
    graphExtraPolationFrac->SetName("extraPolationFracVsMass");

    // loop over input lists
    TIter keyIter(inputFile->GetListOfKeys());
    TKey* key                                               = NULL;
    Int_t listIter                                          = 0;
    while ( (key = (TKey*)keyIter()) ) {

        // input list
        inputList[listIter]                                 = (TList*)key->ReadObj();

        // create directory for plots
        TString outputDir                                   = Form("plots/%s/%s/IntegYields",suffix.Data(),inputList[listIter]->GetName());
        gSystem->Exec("mkdir -p "+outputDir);

        // skip empty or corrupted input lists
        if (!inputList[listIter])
            continue;
        if (inputList[listIter]->IsEmpty() || inputList[listIter]->IsZombie()) {
            cout << inputList[listIter]->GetName() << " empty or zombie" << endl;
            continue;
        }
        cout << "---- " << inputList[listIter]->GetName() << endl;

        // set output list
        outputList->SetName(inputList[listIter]->GetName());

        // get coll. sys, energy and centrality from list
        TString fCollSysTemp, fEnergyTemp, fCentralityTemp;
        TString inputListName                               = inputList[listIter]->GetName();
        TList* currentList                                  = new TList();
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

            // check if stat&sys spectra for particle are contained in input list
            if ( inputList[listIter]->FindObject(Form("%sStat", fParticle[particleIter].Data())) && inputList[listIter]->FindObject(Form("%sSys", fParticle[particleIter].Data())) ) {

                // for testing purposes
                cout << "-------- " << fParticle[particleIter].Data()  << " spectra found" << endl;

                TList* currentParticleList              = new TList();

                // get type stat&sys
                tempStatClass                           = ((TObject*)inputList[listIter]->FindObject(Form("%sStat", fParticle[particleIter].Data())))->ClassName();
                tempSysClass                            = ((TObject*)inputList[listIter]->FindObject(Form("%sSys",  fParticle[particleIter].Data())))->ClassName();

                // get and shift spectra using sys errs
                tempSpecStat                            = (TObject*)inputList[listIter]->FindObject(Form("%sStat", fParticle[particleIter].Data()));
                tempSpecSys                             = (TObject*)inputList[listIter]->FindObject(Form("%sSys",  fParticle[particleIter].Data()));

                if (tempSpecStat){
                    currentList->Add(tempSpecStat);
                    currentParticleList->Add(tempSpecStat);
                }
                if (tempSpecSys){
                    currentList->Add(tempSpecSys);
                    currentParticleList->Add(tempSpecSys);
                }

                Double_t extraYield[3]                  = {0, 0, 0};
                Double_t meanPt[3]                      = {0, 0, 0};

                for (Int_t ff = 0; ff < nFitFunctionsUsed[particleIter][GetCentralityIterator(fCentralityTemp)][0]; ff++){
                    tempSpecFit                         = NULL;
                    TString currFitName                 = initFitFunctionsUsed[particleIter][GetCentralityIterator(fCentralityTemp)][0][ff].Data();
                    cout << currFitName.Data() << endl;
                    if (inputList[listIter]->FindObject(Form("%sStat_Fit_%s", fParticle[particleIter].Data(), currFitName.Data()))) {
                        tempSpecFit       = (TF1*)inputList[listIter]->FindObject(Form("%sStat_Fit_%s", fParticle[particleIter].Data(), currFitName.Data()));
                        if (tempSpecFit){
                            currentList->Add(tempSpecFit);
                            currentParticleList->Add(tempSpecFit);
                        }
                    }
                    cout << Form("Yield_%s_%s", fParticle[particleIter].Data(), currFitName.Data()) << endl;
                    if (inputList[listIter]->FindObject(Form("Yield_%s_%s", fParticle[particleIter].Data(), currFitName.Data()))) {
//                         cout << "found" << endl;
                        TH1D* tempYield       = (TH1D*)inputList[listIter]->FindObject(Form("Yield_%s_%s", fParticle[particleIter].Data(), currFitName.Data()));
                        if (ff == 0){
                            particleYieldStat->SetBinContent(particleIter+1,tempYield->GetBinContent(kYield));
                            particleYieldSys->SetBinContent(particleIter+1,tempYield->GetBinContent(kYield));
                            particleYieldSysFunc->SetBinContent(particleIter+1,tempYield->GetBinContent(kYield));
                            particleYieldStat->SetBinError(particleIter+1,tempYield->GetBinContent(kYieldStat));
                            particleYieldSys->SetBinError(particleIter+1,TMath::Sqrt(TMath::Power(tempYield->GetBinContent(kYieldSysHi),2) + TMath::Power(tempYield->GetBinContent(kYieldSysLo),2)));
                            meanPtStat->SetBinContent(particleIter+1,tempYield->GetBinContent(kMean));
                            meanPtSys->SetBinContent(particleIter+1,tempYield->GetBinContent(kMean));
                            meanPtSysFunc->SetBinContent(particleIter+1,tempYield->GetBinContent(kMean));
                            meanPtStat->SetBinError(particleIter+1,tempYield->GetBinContent(kMeanStat));
                            meanPtSys->SetBinError(particleIter+1,TMath::Sqrt(TMath::Power(tempYield->GetBinContent(kMeanSysHi),2) + TMath::Power(tempYield->GetBinContent(kMeanSysLo),2)));
                            particleYieldExtra->SetBinContent(particleIter+1,tempYield->GetBinContent(kExtra));
                            extraYield[1]               = tempYield->GetBinContent(kExtra);
                            extraYield[0]               = tempYield->GetBinContent(kExtra);
                            extraYield[2]               = tempYield->GetBinContent(kExtra);
                            meanPt[1]                   = tempYield->GetBinContent(kMean);
                            meanPt[0]                   = tempYield->GetBinContent(kMean);
                            meanPt[2]                   = tempYield->GetBinContent(kMean);
                            cout << ff << "\t" << meanPt[1] << endl;
                        } else {
                            if (tempYield->GetBinContent(kExtra) < extraYield[1]){
                                if (tempYield->GetBinContent(kExtra) < extraYield[0])
                                    extraYield[0]       = tempYield->GetBinContent(kExtra);
                            } else {
                                if (tempYield->GetBinContent(kExtra) > extraYield[2])
                                    extraYield[2]       = tempYield->GetBinContent(kExtra);
                            }
                            cout << ff << "\t" << tempYield->GetBinContent(kMean) << endl;
                            if (tempYield->GetBinContent(kMean) < meanPt[1]){
                                if (tempYield->GetBinContent(kMean) < meanPt[0])
                                    meanPt[0]           = tempYield->GetBinContent(kMean);
                            } else {
                                if (tempYield->GetBinContent(kMean) > meanPt[2])
                                    meanPt[2]           = tempYield->GetBinContent(kMean);
                            }
                        }
                    }
                }
                Double_t maxErrYield                    = TMath::Abs(extraYield[1]-extraYield[0]);
                if (TMath::Abs(extraYield[2]-extraYield[1]) > maxErrYield )
                    maxErrYield                         = TMath::Abs(extraYield[2]-extraYield[1]);
                particleYieldSysFunc->SetBinError(particleIter+1,maxErrYield);

                cout << "MEAN PT values " << endl;
                cout << meanPt[1] << "\t" << meanPt[0] << "\t" << meanPt[2] << endl;
                Double_t maxErrMeanPt                   = TMath::Abs(meanPt[1]-meanPt[0]);
                if (TMath::Abs(meanPt[2]-meanPt[1]) > maxErrMeanPt )
                    maxErrMeanPt                        = TMath::Abs(meanPt[2]-meanPt[1]);
                cout << "err:  " << maxErrMeanPt << endl;
                meanPtSysFunc->SetBinError(particleIter+1,maxErrMeanPt);

                Double_t mass                           = TDatabasePDG::Instance()->GetParticle(fParticlePDG[particleIter])->Mass();
                // fill vs mass yields
                graphParticleYieldStat[0]->SetPoint(particleIter, mass, particleYieldStat->GetBinContent(particleIter+1));
                graphParticleYieldStat[0]->SetPointError(particleIter, 0, 0, particleYieldStat->GetBinError(particleIter+1), particleYieldStat->GetBinError(particleIter+1) );
                graphParticleYieldSys[0]->SetPoint(particleIter, mass, particleYieldSys->GetBinContent(particleIter+1));
                graphParticleYieldSys[0]->SetPointError(particleIter, 0, 0, particleYieldSys->GetBinError(particleIter+1), particleYieldSys->GetBinError(particleIter+1) );
                graphParticleYieldSysFunc[0]->SetPoint(particleIter, mass, particleYieldSysFunc->GetBinContent(particleIter+1));
                graphParticleYieldSysFunc[0]->SetPointError(particleIter, 0, 0, particleYieldSysFunc->GetBinError(particleIter+1), particleYieldSysFunc->GetBinError(particleIter+1) );
                // fill vs mass mean Pt
                graphMeanPtStat[0]->SetPoint(particleIter, mass, meanPtStat->GetBinContent(particleIter+1));
                graphMeanPtStat[0]->SetPointError(particleIter, 0, 0, meanPtStat->GetBinError(particleIter+1), meanPtStat->GetBinError(particleIter+1) );
                graphMeanPtSys[0]->SetPoint(particleIter, mass, meanPtSys->GetBinContent(particleIter+1));
                graphMeanPtSys[0]->SetPointError(particleIter, 0, 0, meanPtSys->GetBinError(particleIter+1), meanPtSys->GetBinError(particleIter+1) );
                graphMeanPtSysFunc[0]->SetPoint(particleIter, mass, meanPtSysFunc->GetBinContent(particleIter+1));
                graphMeanPtSysFunc[0]->SetPointError(particleIter, 0, 0, meanPtSysFunc->GetBinError(particleIter+1), meanPtSysFunc->GetBinError(particleIter+1) );

                if (fIsMeson[particleIter]){
                    // fill vs mass yields
                    graphParticleYieldStat[1]->SetPoint(particleIter, mass, particleYieldStat->GetBinContent(particleIter+1));
                    graphParticleYieldStat[1]->SetPointError(particleIter, 0, 0, particleYieldStat->GetBinError(particleIter+1), particleYieldStat->GetBinError(particleIter+1) );
                    graphParticleYieldSys[1]->SetPoint(particleIter, mass, particleYieldSys->GetBinContent(particleIter+1));
                    graphParticleYieldSys[1]->SetPointError(particleIter, 0, 0, particleYieldSys->GetBinError(particleIter+1), particleYieldSys->GetBinError(particleIter+1) );
                    graphParticleYieldSysFunc[1]->SetPoint(particleIter, mass, particleYieldSysFunc->GetBinContent(particleIter+1));
                    graphParticleYieldSysFunc[1]->SetPointError(particleIter, 0, 0, particleYieldSysFunc->GetBinError(particleIter+1), particleYieldSysFunc->GetBinError(particleIter+1) );
                    // fill vs mass mean Pt
                    graphMeanPtStat[1]->SetPoint(particleIter, mass, meanPtStat->GetBinContent(particleIter+1));
                    graphMeanPtStat[1]->SetPointError(particleIter, 0, 0, meanPtStat->GetBinError(particleIter+1), meanPtStat->GetBinError(particleIter+1) );
                    graphMeanPtSys[1]->SetPoint(particleIter, mass, meanPtSys->GetBinContent(particleIter+1));
                    graphMeanPtSys[1]->SetPointError(particleIter, 0, 0, meanPtSys->GetBinError(particleIter+1), meanPtSys->GetBinError(particleIter+1) );
                    graphMeanPtSysFunc[1]->SetPoint(particleIter, mass, meanPtSysFunc->GetBinContent(particleIter+1));
                    graphMeanPtSysFunc[1]->SetPointError(particleIter, 0, 0, meanPtSysFunc->GetBinError(particleIter+1), meanPtSysFunc->GetBinError(particleIter+1) );

                    graphParticleYieldStat[2]->SetPoint(particleIter, -1, 0);
                    graphParticleYieldSys[2]->SetPoint(particleIter, -1, 0);
                    graphParticleYieldSysFunc[2]->SetPoint(particleIter, -1, 0);
                    // fill vs mass mean Pt
                    graphMeanPtStat[2]->SetPoint(particleIter, -1, 0);
                    graphMeanPtSys[2]->SetPoint(particleIter, -1, 0);
                    graphMeanPtSysFunc[2]->SetPoint(particleIter, -1, 0);

                } else {
                    // fill vs mass yields
                    graphParticleYieldStat[2]->SetPoint(particleIter, mass, particleYieldStat->GetBinContent(particleIter+1));
                    graphParticleYieldStat[2]->SetPointError(particleIter, 0, 0, particleYieldStat->GetBinError(particleIter+1), particleYieldStat->GetBinError(particleIter+1) );
                    graphParticleYieldSys[2]->SetPoint(particleIter, mass, particleYieldSys->GetBinContent(particleIter+1));
                    graphParticleYieldSys[2]->SetPointError(particleIter, 0, 0, particleYieldSys->GetBinError(particleIter+1), particleYieldSys->GetBinError(particleIter+1) );
                    graphParticleYieldSysFunc[2]->SetPoint(particleIter, mass, particleYieldSysFunc->GetBinContent(particleIter+1));
                    graphParticleYieldSysFunc[2]->SetPointError(particleIter, 0, 0, particleYieldSysFunc->GetBinError(particleIter+1), particleYieldSysFunc->GetBinError(particleIter+1) );
                    // fill vs mass mean Pt
                    graphMeanPtStat[2]->SetPoint(particleIter, mass, meanPtStat->GetBinContent(particleIter+1));
                    graphMeanPtStat[2]->SetPointError(particleIter, 0, 0, meanPtStat->GetBinError(particleIter+1), meanPtStat->GetBinError(particleIter+1) );
                    graphMeanPtSys[2]->SetPoint(particleIter, mass, meanPtSys->GetBinContent(particleIter+1));
                    graphMeanPtSys[2]->SetPointError(particleIter, 0, 0, meanPtSys->GetBinError(particleIter+1), meanPtSys->GetBinError(particleIter+1) );
                    graphMeanPtSysFunc[2]->SetPoint(particleIter, mass, meanPtSysFunc->GetBinContent(particleIter+1));
                    graphMeanPtSysFunc[2]->SetPointError(particleIter, 0, 0, meanPtSysFunc->GetBinError(particleIter+1), meanPtSysFunc->GetBinError(particleIter+1) );

                    graphParticleYieldStat[1]->SetPoint(particleIter, -1, 0);
                    graphParticleYieldSys[1]->SetPoint(particleIter, -1, 0);
                    graphParticleYieldSysFunc[1]->SetPoint(particleIter, -1, 0);
                    // fill vs mass mean Pt
                    graphMeanPtStat[1]->SetPoint(particleIter, -1, 0);
                    graphMeanPtSys[1]->SetPoint(particleIter, -1, 0);
                    graphMeanPtSysFunc[1]->SetPoint(particleIter, -1, 0);

                }

                graphExtraPolationFrac->SetPoint(particleIter, mass, particleYieldExtra->GetBinContent(particleIter+1)/particleYieldStat->GetBinContent(particleIter+1) );

                ProduceSingleParticleSpectraPlotFromListOnlyFinalWithFits(currentParticleList, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, fParticle[particleIter]);
                ProduceSingleParticleSpectraPlotFromListOnlyFinalWithFits(currentParticleList, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, fParticle[particleIter], kFALSE);
                ProduceSingleParticleSpectraPlotFromListOnlyFinalWithFits(currentParticleList, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, fParticle[particleIter], kFALSE, 2.);
                delete currentParticleList;
            } else {
              cout <<  fParticle[particleIter].Data() << endl;
                // fill vs mass yields
                for (Int_t i = 0; i < 3; i++){
                    graphParticleYieldStat[i]->SetPoint(particleIter, -1, 0);
                    graphParticleYieldSys[i]->SetPoint(particleIter, -1, 0);
                    graphParticleYieldSysFunc[i]->SetPoint(particleIter, -1, 0);
                    // fill vs mass mean Pt
                    graphMeanPtStat[i]->SetPoint(particleIter, -1, 0);
                    graphMeanPtSys[i]->SetPoint(particleIter, -1, 0);
                    graphMeanPtSysFunc[i]->SetPoint(particleIter, -1, 0);
                    graphExtraPolationFrac->SetPoint(particleIter, -1, 0);
                }
            }
        }

        ProduceParticleSpectraPlotFromListOnlyFinalWithFits(currentList, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir);
        ProduceParticleSpectraPlotFromListOnlyFinalWithFits(currentList, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, kFALSE);
        ProduceParticleSpectraPlotFromListOnlyFinalWithFits(currentList, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, kFALSE, 3.);

        graphExtraPolationFrac->Sort();
        while (graphExtraPolationFrac->GetX()[0] < 0 && graphExtraPolationFrac->GetN() > 0) {
            graphExtraPolationFrac->RemovePoint(0);
        }
        if(graphExtraPolationFrac->GetN() > 0) outputList->Add(graphExtraPolationFrac);
        else delete graphExtraPolationFrac;

        // sorting of Yield graphs vs mass
        for (Int_t i = 0; i< 2; i++){
            graphParticleYieldStat[i]->Sort();
            graphParticleYieldSys[i]->Sort();
            graphParticleYieldSysFunc[i]->Sort();
            while (graphParticleYieldStat[i]->GetX()[0] < 0 && graphParticleYieldStat[i]->GetN() > 0) {
                graphParticleYieldStat[i]->RemovePoint(0);
                graphParticleYieldSys[i]->RemovePoint(0);
                graphParticleYieldSysFunc[i]->RemovePoint(0);
            }
            if(graphParticleYieldStat[i]->GetN() > 0) outputList->Add(graphParticleYieldStat[i]);
            else delete graphParticleYieldStat[i];
            if(graphParticleYieldSys[i]->GetN() > 0) outputList->Add(graphParticleYieldSys[i]);
            else delete graphParticleYieldSys[i];
            if(graphParticleYieldSysFunc[i]->GetN() > 0) outputList->Add(graphParticleYieldSysFunc[i]);
            else delete graphParticleYieldSysFunc[i];
        }

        if (graphParticleYieldStat[0] && graphParticleYieldSys[0] && graphParticleYieldSysFunc[0]){
            // write output into log file
            textoutput << "=====================================================================" << endl;
            textoutput << "====================== integrated Yield vs Mass =====================" << endl;
            textoutput << "=====================================================================" << endl;
            textoutput << " mass \t Yield \t stat Err \t sys Err \t sys Err func \t extraPolationFraction" << endl;
            for (Int_t i = 0; i < graphParticleYieldStat[0]->GetN(); i++){
                textoutput << "" << graphParticleYieldStat[0]->GetX()[i] << "\t" << graphParticleYieldStat[0]->GetY()[i] << "\t" << graphParticleYieldStat[0]->GetEYhigh()[i] << "\t" << graphParticleYieldSys[0]->GetEYhigh()[i] << "\t" << graphParticleYieldSysFunc[0]->GetEYhigh()[i] << "\t" << graphExtraPolationFrac->GetY()[i] << endl;
            }
            // plot as function of the mass
            ProduceMassDependentPlot(  graphParticleYieldStat, graphParticleYieldSys, graphParticleYieldSysFunc, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, "Yield", "d#it{N}/d#it{y}", kTRUE, kFALSE);
            ProduceMassDependentPlot(  graphParticleYieldStat, graphParticleYieldSys, graphParticleYieldSysFunc, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, "Yield", "d#it{N}/d#it{y}", kTRUE, kTRUE);

            ProduceMassDependentPlotWithParticleOnAxis(  graphParticleYieldStat, graphParticleYieldSys, graphParticleYieldSysFunc, fCollSysTemp, fEnergyTemp, "", suffix, outputDir, "Yield", "d#it{N}/d#it{y}", kTRUE, kFALSE);
            ProduceMassDependentPlotWithParticleOnAxis(  graphParticleYieldStat, graphParticleYieldSys, graphParticleYieldSysFunc, fCollSysTemp, fEnergyTemp, "", suffix, outputDir, "Yield", "d#it{N}/d#it{y}", kTRUE, kTRUE);

            textoutput << "=====================================================================" << endl;
            textoutput << "====================== integrated Ratio =====================" << endl;
            textoutput << "=====================================================================" << endl;
            textoutput << " particles \t ratio \t stat Err \t sys Err \t sys Err func " << endl;
            for (Int_t i = 0; i < graphParticleYieldStat[0]->GetN(); i++){
                for (Int_t j = 0; j < graphParticleYieldStat[0]->GetN(); j++){
                    if (i == j) continue;

                    Double_t ratio              = graphParticleYieldStat[0]->GetY()[i]/graphParticleYieldStat[0]->GetY()[j];
                    Double_t ratioErrStat       = ratio*TMath::Sqrt(  TMath::Power(graphParticleYieldStat[0]->GetEYhigh()[i]/graphParticleYieldStat[0]->GetY()[i],2)
                                                                    + TMath::Power(graphParticleYieldStat[0]->GetEYhigh()[j]/graphParticleYieldStat[0]->GetY()[j],2));
                    Double_t ratioErrSys        = ratio*TMath::Sqrt(  TMath::Power(graphParticleYieldSys[0]->GetEYhigh()[i]/graphParticleYieldSys[0]->GetY()[i],2)
                                                                    + TMath::Power(graphParticleYieldSys[0]->GetEYhigh()[j]/graphParticleYieldSys[0]->GetY()[j],2));
                    Double_t ratioErrSysFunc    = ratio*TMath::Sqrt(  TMath::Power(graphParticleYieldSysFunc[0]->GetEYhigh()[i]/graphParticleYieldSysFunc[0]->GetY()[i],2)
                                                                    + TMath::Power(graphParticleYieldSysFunc[0]->GetEYhigh()[j]/graphParticleYieldSysFunc[0]->GetY()[j],2));
                    textoutput  << GetParticleLabelFromMass(graphParticleYieldStat[0]->GetX()[i])<< "/" << GetParticleLabelFromMass(graphParticleYieldStat[0]->GetX()[j]) << "\t:\t"
                                << ratio << "\t" << ratioErrStat << "\t" << ratioErrSys << "\t" << ratioErrSysFunc << endl;
                }
            }
        }

        // sorting of mean pt graphs vs mass
        for (Int_t i = 0; i< 2; i++){
            graphMeanPtStat[i]->Sort();
            graphMeanPtSys[i]->Sort();
            graphMeanPtSysFunc[i]->Sort();
            while (graphMeanPtStat[i]->GetX()[0] < 0 && graphMeanPtStat[i]->GetN() > 0){
                graphMeanPtStat[i]->RemovePoint(0);
                graphMeanPtSys[i]->RemovePoint(0);
                graphMeanPtSysFunc[i]->RemovePoint(0);
            }
            if(graphMeanPtStat[i]->GetN() > 0) outputList->Add(graphMeanPtStat[i]);
            else delete graphMeanPtStat[i];
            if(graphMeanPtSys[i]->GetN() > 0) outputList->Add(graphMeanPtSys[i]);
            else delete graphMeanPtSys[i];
            if(graphMeanPtSysFunc[i]->GetN() > 0) outputList->Add(graphMeanPtSysFunc[i]);
            else delete graphMeanPtSysFunc[i];
        }
        if (graphMeanPtStat[0] && graphMeanPtSys[0] && graphMeanPtSysFunc[0]){
            // write output into log file
            textoutput << "=====================================================================" << endl;
            textoutput << "====================== mean Pt vs Mass ==============================" << endl;
            textoutput << "=====================================================================" << endl;
            textoutput << " mass \t mean pt \t stat Err \t sys Err \t sys Err func " << endl;
            for (Int_t i = 0; i < graphMeanPtStat[0]->GetN(); i++){
                textoutput << "" << graphMeanPtStat[0]->GetX()[i] << "\t" << graphMeanPtStat[0]->GetY()[i] << "\t" << graphMeanPtStat[0]->GetEYhigh()[i] << "\t" << graphMeanPtSys[0]->GetEYhigh()[i] << "\t" << graphMeanPtSysFunc[0]->GetEYhigh()[i] << endl;
            }
            // plot as function of the mass
            ProduceMassDependentPlot(  graphMeanPtStat, graphMeanPtSys, graphMeanPtSysFunc, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, "MeanPt", "#LT #it{p}_{T} #GT (GeV/#it{c})", kFALSE, kFALSE);
            ProduceMassDependentPlot(  graphMeanPtStat, graphMeanPtSys, graphMeanPtSysFunc, fCollSysTemp, fEnergyTemp, centrality, suffix, outputDir, "MeanPt", "#LT #it{p}_{T} #GT (GeV/#it{c})", kFALSE, kTRUE);
        }

        outputList->Add(particleYieldStat);
        outputList->Add(particleYieldSys);
        outputList->Add(particleYieldSysFunc);
        outputList->Add(meanPtStat);
        outputList->Add(meanPtSys);
        outputList->Add(meanPtSysFunc);
        outputList->Add(particleYieldExtra);

        listIter++;
    }
    textoutput.close();

    // write output lists to file
    cout << "---- writing to " << outputFile->GetName() << endl;
    outputFile->cd();
        if (outputList) outputList->Write(outputList->GetName(), TObject::kOverwrite+TObject::kSingleKey);
    outputFile->Close();


    // free pointer
    delete inputFile;
    delete outputFile;
    delete[] inputList;
    delete outputList;
}
