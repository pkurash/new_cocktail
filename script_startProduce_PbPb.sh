#! /bin/bash

# Produce new inputsFile for PbPb 2.76 TeV
if [ $1 = "newInputPbPb" ]; then
    # run pp input file compilation
    root -b -x -q -l 'CocktailProduceCompleteInputFilePbPb.C+()'

# Produce new inputsFile for PbPb 5.02 TeV
elif [ $1 = "newInputPbPb5TeV" ]; then

    rm PbPb/5020_GeV/Pi0*
    rm PbPb/5020_GeV/Eta*
    basedir=/home/meike/analysis/results/photonconvResults/PbPb/cocktailInput/
    cp $basedir/60110613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1Correction_60110613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/61210613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1Correction_61210613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/51210613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1Correction_51210613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/52310613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1Correction_52310613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/53410613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1Correction_53410613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/54610613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1Correction_54610613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/56810613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1Correction_56810613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/58910613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1Correction_58910613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/

    cp $basedir/60110613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Eta_data_GammaConvV1Correction_60110613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/61210613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Eta_data_GammaConvV1Correction_61210613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/51210613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Eta_data_GammaConvV1Correction_51210613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/52310613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Eta_data_GammaConvV1Correction_52310613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/53410613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Eta_data_GammaConvV1Correction_53410613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/54610613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Eta_data_GammaConvV1Correction_54610613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/56810613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Eta_data_GammaConvV1Correction_56810613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/
    cp $basedir/58910613_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Eta_data_GammaConvV1Correction_58910613_00200009247602008250404000_0652501500000000.root PbPb/5020_GeV/

    # 0-5-10-20-30-40-60-80-90%  Neutral pion, eta, charged pion, charged kaon, (anti-)proton, phi, K0*, K0s, Lambda, charged Omega, charged Xi spectra
    # root -b -x -q -l 'CocktailProduceCompleteInputFilePbPb.C+("","011010001110110","1100011101100001100110","eps",kTRUE,1,kTRUE,kFALSE,kFALSE)'
    # 0-5-10-20-30-40-60-80-90% Neutral pion and eta, charged pion and kaon, (anti-)proton, K0s
    root -b -x -q -l 'CocktailProduceCompleteInputFilePbPb.C+("","011010001110110001111","1100011100000001","eps",kTRUE,1,kTRUE,kFALSE,kFALSE)'
    #root -b -x -q -l 'CocktailProduceCompleteInputFilePbPb.C+("","01","0000000000000001","eps",kTRUE,1,kTRUE,kFALSE,kFALSE)' # 0-5% K0S

# parametrize inputs for 2.76 TeV
elif [ $1 = "parameterPbPb276" ]; then
    root -l -b -x -q 'CocktailInputParametrization.C+("CocktailInputPbPb.root","PbPb_2.76TeV","eps","parametrizationSettings/PbPb_2.76TeV_standard.dat","parametrizationSettings/PbPb_2.76TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_PbPb_2.76TeV.txt")'

# parametrize inputs for 5.02 TeV
elif [ $1 = "parameterPbPb5TeV" ]; then
    root -l -b -x -q 'CocktailInputParametrization.C+("CocktailInputPbPb.root","PbPb_5.02TeV","eps","parametrizationSettings/PbPb_5.02TeV_standard.dat","","parametrizationSettings/listAllUniqueNames_PbPb_5.02TeV.txt",kFALSE,kFALSE,1)'


#Integrated yields parametrizations
elif [ $1 = "calcIntegYieldPbPb276" ]; then
    # parametrize inputs
    rm CocktailInputPbPb_IntegYield_PbPb_2_76TeV.root
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPbPb.root","PbPb_2.76TeV","eps","paramSettingsIntegYield/PbPb_2.76TeV_Comb_Fit_modKPt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPbPb.root","PbPb_2.76TeV","eps","paramSettingsIntegYield/PbPb_2.76TeV_Comb_Fit_modK.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPbPb.root","PbPb_2.76TeV","eps","paramSettingsIntegYield/PbPb_2.76TeV_Comb_Fit_TCMPt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPbPb.root","PbPb_2.76TeV","eps","paramSettingsIntegYield/PbPb_2.76TeV_Comb_Fit_TCM.dat")'

#Plot integrated yields nicely together
elif [ $1 = "compIntegYieldPbPb276" ]; then
    root -x -l -b -q 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPbPb_IntegYield_PbPb_2_76TeV.root","PbPb_2.76TeV","","eps","paramSettingsIntegYield/PbPb_2.76TeV_YieldComp.dat")'
fi
