#if !defined(__CINT__) || defined(__CLING__)
  R__LOAD_LIBRARY(libAMPT)
  R__LOAD_LIBRARY(libTAmpt)
  R__LOAD_LIBRARY(libEGPythia6)
  R__LOAD_LIBRARY(libpythia6)
  R__LOAD_LIBRARY(libAliPythia6)
  #include "AliGenerator.h"
  #include "AliGenAmpt.h"
  #include "AliDecayerPythia.h"
  #include "AliGenEMCocktailV2.h"
  #include "AliRun.h"
#endif
// runStandaloneCocktail.C
// =====================
// This macro runs cocktail generation with AliMCGen framework
//
// orginal author: M. Verweij
// modified by: F.Bock, L.Altenkaemper, N.Schmidt

#include <ctime>
#include "TGrid.h"

//const Bool_t saveManager = kFALSE;
const Bool_t saveManager = kTRUE;
void runUseSavedManager(TString, Int_t, Bool_t); 
void LoadLibs();
void runStandaloneCocktail( Long64_t nEvents        = 200,
                            Int_t collisionsSystem  = 0,            // 0: pp 7TeV, 1: pPb 5TeV
                            Int_t motherSelect      = 262143, /*257,67108863, 262143*/// includes all particles in generator (see end of macro for further examples)
                            Int_t decayMode         = 0,             // 0: kAll, 1: kGammaEM, 2: kElectronEM, 3: kDiElectronEM
			    TString paramFile = "parametrizations/Parametrizations_13TeV_mt_full_Andreas.root",
			    TString paramFileDir = "param"
){

    AliAnalysisManager* mgr         = new AliAnalysisManager("MCGenEMCocktail");
    AliDummyHandler*    dumH        = new AliDummyHandler();

    // create blank ESD event
    AliESDEvent *esdE               = new AliESDEvent();
    esdE->CreateStdContent();
    AliESDVertex *vtx               = new AliESDVertex(0.,0.,100);
    vtx->SetName("VertexTracks");
    vtx->SetTitle("VertexTracks");
    esdE->SetPrimaryVertexTracks(vtx);
    if (esdE->GetPrimaryVertex()) Printf("vtx set");
    dumH->SetEvent(esdE);
    mgr->SetInputEventHandler(dumH);

    // greate MC input Handler
    AliMCGenHandler* mcInputHandler = new AliMCGenHandler();
    mgr->SetMCtruthEventHandler(mcInputHandler);

    // cocktail generV2ator
    AliGenerator* gener = 0x0;

//    if (collisionsSystem == 0) {
    AliGenEMCocktailV2 *generV2     = new AliGenEMCocktailV2();

    AliDecayerPythia *decayer  = new AliDecayerPythia();
    decayer->SetDecayerExodus();
    decayer->DecayLongLivedParticles();
  
    generV2->SetParametrizationFile(paramFile);
    generV2->SetParametrizationFileDirectory(paramFileDir);
    generV2->SetNPart(10);                       // source multiplicity per event
    generV2->SetPtRange(0.5,100.);
    generV2->SetFixedEventPlane(0) ;
    generV2->SetDynamicalPtRange(0);
    generV2->SetUseYWeighting(0);
    generV2->SetYRange(-1., 1.);
    generV2->SetPhiRange(0., 360.);
    generV2->SetOrigin(0.,0.,0.); 
    generV2->SetSigma(0.,0.,0.);
    generV2->SetVertexSmear(kPerEvent);
    generV2->SetTrackingFlag(0);
    generV2->SelectMotherParticles(motherSelect);
    generV2->SetCollisionSystem((AliGenEMlibV2::CollisionSystem_t)collisionsSystem);
    Int_t centrality = 0;
    generV2->SetCentrality((AliGenEMlibV2::Centrality_t)centrality);
    generV2->SetDecayer(decayer);

    generV2->SetDecayMode(kGammaEM);    	// kGammaEM      => single photon
    generV2->SetWeightingMode(kNonAnalog); 	// select weighting:
                      // kNonAnalog => weight ~ dN/dp_T
                      // kAnalog    => weight ~ 1
    Printf("I-AddMCEMCocktailV2 calls AliGenEMCocktailV2::CreateCocktail()");
    generV2->CreateCocktail();
    Printf("I-AddMCEMCocktailV2 calls AliGenEMCocktailV2::Init()");
    generV2->Init();

 /*
    gener = reinterpret_cast<AliGenerator*>
        (gInterpreter->ExecuteMacro(
	 Form("$ALICE_PHYSICS/PWG/Cocktail/macros/AddMCEMCocktailV2.C(%d, %d, %d, %d, \"%s\", \"%s\", %d, %f, %f, %d, %d, %d, %d, %d",
	 collisionsSystem, 0, decayMode, motherSelect, paramFile.Data(), paramFileDir.Data(), 100, 0.5, 100., 2000, 1, 1, 0, 0)));
*/
    gener = reinterpret_cast<AliGenerator*>(generV2);
  //  }

    mcInputHandler->SetGenerator(gener);
    mcInputHandler->SetSeedMode(2);

    AliAnalysisTask *taskA = reinterpret_cast<AliAnalysisTask*>(
          gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_GammaCocktailMC.C(kTRUE, 50, 0.05, \"60\")"));

    AliAnalysisTask *taskB = reinterpret_cast<AliAnalysisTask*>(
                     gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_HadronicCocktailMC.C(0, kFALSE, 50, 0.05, \"60\")"));
    AliAnalysisTask *taskC = reinterpret_cast<AliAnalysisTask*>(
                     gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_HadronicCocktailMC.C(1, kFALSE, 50, 0.05, \"60\")"));
    AliAnalysisTask *taskD = reinterpret_cast<AliAnalysisTask*>(
                     gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_HadronicCocktailMC.C(2, kFALSE, 50, 0.05, \"60\")"));

    // give some indication everything is running fine
    mgr->SetUseProgressBar(1, 10);
    gAlice = new AliRun();

    // Run analysis
    if(saveManager) {
        TFile *fileMgr = new TFile("AnalysisManager.root","RECREATE");
        mgr->Write();
        std::cout << "Written" << endl;
        fileMgr->Write();
        fileMgr->Close();

        runUseSavedManager("AnalysisManager.root",nEvents, 0);
    } else {
        mgr->SetDebugLevel(11);
        mgr->InitAnalysis();
        mgr->PrintStatus();
        mgr->EventLoop(nEvents);
    }

}

void runUseSavedManager(TString fileMgr = "AnalysisManager.root", Int_t nEvents = 1000, Bool_t useGrid = kFALSE) 
{

    Printf("run saved");
    //LoadLibs();

    TFile *f = new TFile(fileMgr.Data());
//    AliAnalysisManager *mgr = dynamic_cast<AliAnalysisManager*>(f->Get("MCGenThermalModel"));
    AliAnalysisManager *mgr = dynamic_cast<AliAnalysisManager*>(f->Get("MCGenEMCocktail"));
     if(!mgr) {
        Printf("Did not find manager");
        return;
    }
    gRandom->Print();
    mgr->SetDebugLevel(11);
    mgr->InitAnalysis();
    mgr->PrintStatus();
    gRandom->Print();
    mgr->EventLoop(nEvents);
}

void LoadLibs()
{

  // Load common libraries (better too many than too few)
	gSystem->Load("libCore.so");  
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libGui.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libMinuit");
	gSystem->Load("libMinuit2");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libOADB");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");  
	gSystem->Load("libPWGGAGammaConv.so");
	gSystem->Load("libEve.so");   
	gSystem->Load("libCDB.so");
	gSystem->Load("libProof.so");
	gSystem->Load("libRAWDatabase.so");
	gSystem->Load("libSTEER.so");
	gSystem->Load("libEVGEN.so");
	gSystem->Load("libESDfilter.so");	
	gSystem->Load("libTender.so");
	gSystem->Load("libTOFbase.so");
	gSystem->Load("libTOFsim.so");
	gSystem->Load("libTOFrec.so");
	gSystem->Load("libTRDbase.so");
	gSystem->Load("libVZERObase.so");
	gSystem->Load("libVZEROrec.so");
	//gSystem->Load("libTenderSupplies.so");
	gSystem->Load("libXMLParser.so");
	gSystem->Load("libSTAT.so");
	gSystem->Load("libRAWDatarec.so");
	gSystem->Load("libANALYSIScalib.so");
	gSystem->Load("libCORRFW.so");
	gSystem->Load("libPWGUDbase.so");
	gSystem->Load("libTPCbase.so");
	gSystem->Load("libTPCrec.so");
	gSystem->Load("libTPCcalib.so");
	gSystem->Load("libTRDrec.so");
	gSystem->Load("libITSbase.so");
	gSystem->Load("libITSrec.so");
	gSystem->Load("libHMPIDbase.so");
	gSystem->Load("libPWGmuon.so");
	gSystem->Load("libPWGPP.so");
	gSystem->Load("libPWGHFbase.so");
	gSystem->Load("libPWGDQdielectron.so");
	gSystem->Load("libPWGHFhfe.so");
	gSystem->Load("libPHOSUtils.so");
	gSystem->Load("libPHOSbase.so");
	gSystem->Load("libPHOSpi0Calib.so");
	gSystem->Load("libPHOSrec.so");
	gSystem->Load("libPWGGAPHOSTasks.so");
	gSystem->Load("libPHOSshuttle.so");
	gSystem->Load("libPHOSsim.so");
	gSystem->Load("libPWGCaloTrackCorrBase.so");
	gSystem->Load("libEMCALUtils.so");
	gSystem->Load("libEMCALraw.so");
	gSystem->Load("libEMCALbase.so");
	gSystem->Load("libEMCALrec.so");
	gSystem->Load("libPWGTools.so");
	//gSystem->Load("libPWGEMCAL.so");
	gSystem->Load("libPWGGAUtils.so");
	gSystem->Load("libPWGGAEMCALTasks.so");
	gSystem->Load("libPWGGAGammaConv.so");
	gSystem->Load("libPWGGAPHOSTasks.so");
	gSystem->Load("libPWGCFCorrelationsBase.so");
	gSystem->Load("libPWGCFCorrelationsDPhi.so");
  
}
//================= examples for inclusion of mother particles ================================
// (62591)_10 = (1111010001111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, J/psi, Delta+, Delta0, Rho+, Rho-, K0*

// (29887)_10 = (0111010010111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, Sigma0, Delta+, Delta0, Rho+, Rho-, K0*

// (262143)_10 = (111111111111111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, Jpsi, sigma0, k0s, delta++, delta+, delta-, delta0, rho+, rho-, k0*, k0l, lambda

// (229375)_10 = (110111111111111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, Jpsi, sigma0, k0s, delta++, delta+, delta-, delta0, rho+, rho-, k0l, lambda

// (262079)_10 = (111111111110111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, sigma0, k0s, delta++, delta+, delta-, delta0, rho+, rho-, k0*, k0l, lambda

// (131072)_10 = (100000000000000000)_2
// includes: lambda

// (65536)_10 = (010000000000000000)_2
// includes: k0l

// (256)_10 = (000000000100000000)_2
// includes: k0s

// (257)_10 = (000000000100000001)_2
// includes: k0s, pi0

// (196864)_10 = (110000000100000000)_2
// includes: k0s, k0l, lambda

// (196875)_10 = (110000000100001011)_2
// includes: pi0, eta, omega, k0s, k0l, lambda

// (67108863)_10 = (11111111111111111111111111)_2
// includes: pi0, eta, rho?, omega, eta', phi, JPsi, Sigma0, K0s, Delta++, Delta+, Delta-, Delta0, rho+, rho-, K0*, K0l, Lambda, K+, K-, Omega+, Omega-, Xi+, Xi-, Sigma+, Sigma-
// (67108863)_10 = (11011000100000001000000000)_2
// includes: pi0, eta, rho?, omega, eta', phi, JPsi, Sigma0, K0s, Delta++, Delta+, Delta-, Delta0, rho+, rho-, K0*, K0l, Lambda, K+, K-, Omega+, Omega-, Xi+, Xi-, Sigma+, Sigma-
//=============================================================================================

