#! /bin/bash

# produce new input file for PP
if [ $1 = "newInputPP" ]; then
    root -b -x -q -l 'CocktailProduceCompleteInputFilePP.C+()'

# add MC input for 8TeV PP
elif [ $1 = "newInputPPMC" ]; then
    root -l -b -q -x 'CocktailProduceCompleteInputFilePP.C+("000010",1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,"dN/dydpT","eps",kTRUE)'
    mv CocktailInputPP.root CocktailInputPP_MC.root
    root -b -x -q -l 'CocktailProduceCompleteInputFilePP.C+()'

# parametrize inputs for cockail calc
elif [ $1 = "parameter900" ]; then                  #900 GeV
    rm parametrizations/pp_900GeV.root
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_0.9TeV","eps","parametrizationSettings/pp_0.9TeV_standard.dat","parametrizationSettings/pp_0.9TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_900GeV.txt")'
elif [ $1 = "parameter2760" ]; then                 #2.76 TeV
    rm parametrizations/pp_2760GeV.root
    root -b -x -q -l 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_2.76TeV","eps","parametrizationSettings/pp_2.76TeV_standard.dat","parametrizationSettings/pp_2.76TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_2.76TeV.txt",kFALSE,kFALSE,0,kTRUE)'
elif [ $1 = "parameter5" ]; then                    #8 TeV
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_5TeV","eps","parametrizationSettings/pp_5TeV_standard.dat","parametrizationSettings/pp_5TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_5TeV.txt",kFALSE,kFALSE,0,kFALSE)'
elif [ $1 = "parameter7" ]; then                    #7 TeV
    rm parametrizations/pp_7TeV.root
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","eps","parametrizationSettings/pp_7TeV_standard.dat","parametrizationSettings/pp_7TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_7TeV.txt",kFALSE,kFALSE,0,kFALSE,kTRUE)'
    elif [ $1 = "parameter7nc" ]; then                    #7 TeV
        rm parametrizations/pp_7TeV.root
        root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","eps","parametrizationSettings/pp_7TeV_standard.dat","parametrizationSettings/pp_7TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_7TeV.txt",kFALSE,kFALSE,0,kFALSE)'
elif [ $1 = "parameter7graphs" ]; then                    #7 TeV
    rm parametrizations/pp_7TeV.root
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","pdf","parametrizationSettings/pp_7TeV_standard_Graph.dat","parametrizationSettings/pp_7TeV_ratio_standard_Graph.dat","parametrizationSettings/listAllUniqueNames_7TeV.txt")'
elif [ $1 = "parameter8" ]; then                    #8 TeV
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","eps","parametrizationSettings/pp_8TeV_standard.dat","parametrizationSettings/pp_8TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_8TeV.txt",kFALSE,kFALSE,0,kFALSE)'
elif [ $1 = "parameter8MC" ]; then                  #8 TeV MC
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP_MC.root","pp_8TeV","eps","parametrizationSettings/pp_8TeV_standardMC.dat","","parametrizationSettings/listAllUniqueNames_8TeV_MC.txt",kFALSE,kFALSE,2)'
elif [ $1 = "parameter8Full" ]; then                    #8 TeV
    rm parametrizations/pp_8TeV.root
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","eps","parametrizationSettings/pp_8TeV_standard.dat","parametrizationSettings/pp_8TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_8TeV.txt")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP_MC.root","pp_8TeV","eps","parametrizationSettings/pp_8TeV_standardMC.dat","","parametrizationSettings/listAllUniqueNames_8TeV_MC.txt",kFALSE,kFALSE,2)'

elif [ $1 = "parameterPPRun1" ]; then                  #all run1
    rm parametrizations/pp*.root
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_0.9TeV","pdf","parametrizationSettings/pp_0.9TeV_standard.dat","","parametrizationSettings/listAllUniqueNames_900GeV.txt")'
    root -b -x -q -l 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_2.76TeV","eps","parametrizationSettings/pp_2.76TeV_standard.dat","parametrizationSettings/pp_2.76TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_2.76TeV.txt")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","pdf","parametrizationSettings/pp_7TeV_standard.dat","parametrizationSettings/pp_7TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_7TeV.txt")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","parametrizationSettings/pp_8TeV_ratio_standard.dat","parametrizationSettings/listAllUniqueNames_8TeV.txt")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP_MC.root","pp_8TeV","eps","parametrizationSettings/pp_8TeV_standardMC.dat","","parametrizationSettings/listAllUniqueNames_8TeV_MC.txt",kFALSE,kFALSE,2)'

#Integrated yields parametrizations
elif [ $1 = "parameterIntYield900" ]; then         #0.9 TeV
    rm CocktailInputPP_IntegYield_0_9TeV.*         #can't write same fits into list twice
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPP.root","0.9TeV","eps","integratedyieldParametrizations/pp_0.9TeV_Comb_Fit_tcmpt.dat",kTRUE)'
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPP.root","0.9TeV","eps","integratedyieldParametrizations/pp_0.9TeV_Comb_Fit_tmpt.dat",kTRUE)'
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPP.root","0.9TeV","eps","integratedyieldParametrizations/pp_0.9TeV_Comb_Fit_oHagPt.dat",kTRUE)'
elif [ $1 = "calcIntegYieldPP2760" ]; then          #2.76 TeV
    rm CocktailInputPP_IntegYield_2_76TeV.*         #can't write same fits into list twice
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPP.root","2.76TeV","eps","paramSettingsIntegYield/pp_2.76TeV_Comb_Fit_TCMPt.dat",kTRUE)'
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPP.root","2.76TeV","eps","paramSettingsIntegYield/pp_2.76TeV_Comb_Fit_TsallisPt.dat",kTRUE)'
    root -b -x -q -l 'CalculateIntegratedYield.C+("CocktailInputPP.root","2.76TeV","eps","paramSettingsIntegYield/pp_2.76TeV_Comb_Fit_oHagPt.dat",kTRUE)'
elif [ $1 = "parameterIntYield7" ]; then            #7 TeV
    rm CocktailInputPP_IntegYield_7TeV.*            #can't write same fits into list twice
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","7TeV","eps","integratedyieldParametrizations/pp_7TeV_Comb_Fit_oHagPt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","7TeV","eps","integratedyieldParametrizations/pp_7TeV_Comb_Fit_oHagPt_2.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","7TeV","eps","integratedyieldParametrizations/pp_7TeV_Comb_Fit_tmpt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","7TeV","eps","integratedyieldParametrizations/pp_7TeV_Comb_Fit_tcmpt.dat")'
elif [ $1 = "parameterIntYield8" ]; then            #8 TeV
    rm CocktailInputPP_IntegYield_8TeV.*            #can't write same fits into list twice
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","8TeV","eps","integratedyieldParametrizations/pp_8TeV_Comb_Fit_oHagPt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","8TeV","eps","integratedyieldParametrizations/pp_8TeV_Comb_Fit_tmpt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","8TeV","eps","integratedyieldParametrizations/pp_8TeV_Comb_Fit_tcmpt.dat")'

#Plot integrated yields nicely together
elif [ $1 = "plotIntYield900" ]; then               #0.9TeV
    root -b -x -q -l 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_0_9TeV.root","0.9TeV","","eps","paramSettingsIntegYield/pp_0.9TeV_YieldComp.dat",kTRUE)'
    root -b -x -q -l 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_0_9TeV.root","0.9TeV","","pdf","paramSettingsIntegYield/pp_0.9TeV_YieldComp.dat",kTRUE)'
elif [ $1 = "plotIntYield2760" ]; then              #2.76 TeV
    root -b -x -q -l 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_2_76TeV.root","2.76TeV","","eps","paramSettingsIntegYield/pp_2.76TeV_YieldComp.dat",kTRUE)'
    root -b -x -q -l 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_2_76TeV.root","2.76TeV","","pdf","paramSettingsIntegYield/pp_2.76TeV_YieldComp.dat",kTRUE)'
elif [ $1 = "plotIntYield7" ]; then                 #7 TeV
    root -x -l -b -q 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_7TeV.root","7TeV","","eps","paramSettingsIntegYield/pp_7TeV_YieldComp.dat")'
    root -x -l -b -q 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_7TeV.root","7TeV","","pdf","paramSettingsIntegYield/pp_7TeV_YieldComp.dat")'
elif [ $1 = "plotIntYield8" ]; then                 #8 TeV
    root -x -l -b -q 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_8TeV.root","8TeV","","eps","paramSettingsIntegYield/pp_8TeV_YieldComp.dat")'
    root -x -l -b -q 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_8TeV.root","8TeV","","pdf","paramSettingsIntegYield/pp_8TeV_YieldComp.dat")'
fi
