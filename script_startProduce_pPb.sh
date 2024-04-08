#/bin/bash

#compose new PPb input
if [ $1 = "newInputPPb" ]; then
    # run pPb input file compilation
    root -b -x -q -l 'CocktailProduceCompleteInputFilePPb.C+("11111111111111111","0000000000000000","1100011101100001100110111111","pdf",kTRUE,kTRUE,kFALSE)'
# parametrize inputs for 5 TeV all cents
elif [ $1 = "parameterPPb5All" ]; then
    # parametrize inputs
    root -b -x -q -l 'CocktailInputParametrization.C+("CocktailInputPPb.root","5TeV","eps","parametrizationSettings/pPb_5TeV_standard.dat","parametrizationSettings/pPb_5TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_pPb_5TeV_MB.txt",kFALSE,kFALSE,0,kFALSE)'

# parametrize inputs for 5 TeV MB
elif [ $1 = "parameterPPb5MB" ]; then
    root -b -x -q -l 'CocktailInputParametrization.C+("CocktailInputPPb.root","5TeV","eps","parametrizationSettings/pPb_5TeV_standard_onlyMB.dat","parametrizationSettings/pPb_5TeV_ratio_standard_onlyMB.dat","parametrizationSettings/listAllUniqueNames_pPb_5TeV_MB.txt",kFALSE,kFALSE,0,kFLASE)'

# parametrize inputs for 5 TeV MB
elif [ $1 = "parameterPPb5Cent" ]; then
    # parametrize inputs
    root -b -x -q -l 'CocktailInputParametrization.C+("CocktailInputPPb.root","5TeV","eps","parametrizationSettings/pPb_5TeV_standard_cent.dat","parametrizationSettings/pPb_5TeV_ratio_standard_cent.dat","parametrizationSettings/listAllUniqueNames_pPb_5TeV_MB.txt",kFALSE,kFALSE,0,kTFALSE)'


#Integrated yields parametrizations
elif [ $1 = "calcIntegYieldpPbMB" ]; then
    # parametrize inputs
    rm CocktailInputPPb_IntegYield_5TeV.root
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPPb.root","5TeV","eps","paramSettingsIntegYield/pPb_5TeV_MB_Comb_Fit_oHagPt.dat",kTRUE)'
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPPb.root","5TeV","eps","paramSettingsIntegYield/pPb_5TeV_MB_Comb_Fit_TsallisPt.dat",kTRUE)'
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPPb.root","5TeV","eps","paramSettingsIntegYield/pPb_5TeV_MB_Comb_Fit_TsallisPt2.dat",kTRUE)'
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPPb.root","5TeV","eps","paramSettingsIntegYield/pPb_5TeV_MB_Comb_Fit_TCMPt.dat",kTRUE)'

#Plot integrated yields nicely together
elif [ $1 = "compIntegYieldpPbMB" ]; then
#     root -b -x -q -l 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPPb_IntegYield_5TeV.root","5TeV","NSD","eps","paramSettingsIntegYield/pPb_5TeV_MB_YieldComp.dat",kFALSE)'
    root -b -x -q -l 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPPb_IntegYield_5TeV.root","5TeV","NSD","pdf","paramSettingsIntegYield/pPb_5TeV_MB_YieldComp.dat",kFALSE)'
fi
