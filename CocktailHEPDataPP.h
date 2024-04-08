Double_t        xSectionpp7TeV          = 73.2*1e-3;    // INEL!
Double_t        xSectionpp7TeVV0OR      = 62.22*1e-3;
Double_t        xSectionpp7TeVV0AND     = 54.31*1e-3;
Double_t        xSectionpp7TeVErrUp     = 2.18;
Double_t        xSectionpp7TeVErrDown   = 2.18;
Double_t        xSectionpp900GeV        = 47.78*1e-3;
Double_t        xSectionpp900GeVV0AND   = 40.06*1e-3;
Double_t        xSectionpp900GeVErrUp   = 2.39;
Double_t        xSectionpp900GeVErrDown = 1.86;
Double_t        xSectionpp2760GeV       = 55.416*1e-3;
Double_t        xSectionpp2760GeVV0AND  = 47.73*1e-3;
Double_t        xSectionpp2760GeVErr    = 3.9;
Double_t        recalcBarn              = 1e12;         //NLO in pbarn!!!!
Double_t        xSection7TeVppINEL      = 73.2*1e9;
Double_t        xSection2760GeVppINEL   = 62.8*1e9;
Double_t        xSection900GeVppINEL    = 52.5*1e9;
Double_t        mutopico                = 1e6;

TGraphAsymmErrors* CalculateParticleRatioWithFit(TGraphAsymmErrors* spec, TF1* fit) {

    if (!fit) {
        cout << "ERROR: Fit is NULL!" << endl;
        return NULL;
    }

    Int_t numPoints                     = spec->GetN();
    Double_t* xVal                      = spec->GetX();
    Double_t* xErrPlus                  = spec->GetEXhigh();
    Double_t* xErrMinus                 = spec->GetEXlow();
    Double_t* yVal                      = spec->GetY();
    Double_t* yErrPlus                  = spec->GetEYhigh();
    Double_t* yErrMinus                 = spec->GetEYlow();

    Double_t* yValRatio                 = new Double_t[numPoints];
    Double_t* yErrPlusRatio             = new Double_t[numPoints];
    Double_t* yErrMinusRatio            = new Double_t[numPoints];

    for (Int_t i=0; i<numPoints; i++) {
        Double_t integralFit            = fit->Integral(xVal[i]-xErrMinus[i], xVal[i]+xErrPlus[i])/(xErrMinus[i]+xErrPlus[i]);

        yValRatio[i]                    = yVal[i]/integralFit;
        yErrMinusRatio[i]               = yErrMinus[i]/integralFit;
        yErrPlusRatio[i]                = yErrPlus[i]/integralFit;
    }

    TGraphAsymmErrors* returnGraph      = new TGraphAsymmErrors(numPoints, xVal, yValRatio, xErrMinus, xErrPlus, yErrMinusRatio, yErrPlusRatio);
    return returnGraph;
}

//================================================================================================================
// phi to pi0 in pp 7TeV
//================================================================================================================
TGraphAsymmErrors* GetPhiToChargedPionpp7TeV(TF1* fitPi0, Bool_t returnStat){
    
    // Plot: p8208_d2x1y1
    double p8208_d2x1y1_xval[]          = { 0.45, 0.55, 0.65, 0.75, 0.85,
        0.95, 1.05, 1.15, 1.25, 1.35,
        1.45, 1.55, 1.65, 1.75, 1.85,
        1.95, 2.1, 2.3, 2.5, 2.7,
        2.9, 3.25, 3.75, 4.25, 4.75,
        5.5 };
    double p8208_d2x1y1_xerrminus[]     = { 0.04999999999999999, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999993,
        0.04999999999999993, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044,
        0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044,
        0.050000000000000044, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009,
        0.10000000000000009, 0.25, 0.25, 0.25, 0.25,
        0.5 };
    double p8208_d2x1y1_xerrplus[]      = { 0.04999999999999999, 0.04999999999999993, 0.04999999999999993, 0.050000000000000044, 0.050000000000000044,
        0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982,
        0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982,
        0.050000000000000044, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964,
        0.10000000000000009, 0.25, 0.25, 0.25, 0.25,
        0.5 };
    double p8208_d2x1y1_yval[]          = { 0.02464, 0.02461, 0.02284, 0.02289, 0.0202,
        0.01873, 0.01653, 0.01502, 0.01333, 0.0117,
        0.01003, 0.00903, 0.007908, 0.006922, 0.00577,
        0.004736, 0.003904, 0.002982, 0.002141, 0.001679,
        0.001309, 8.073E-4, 4.712E-4, 2.835E-4, 1.751E-4,
        8.503E-5 };
    double p8208_d2x1y1_yerrminus[]     = { 0.00291547594742265, 0.0023482972554597936, 0.0019859506539690254, 0.002015589243868899, 0.0017705366418123065,
        0.001586190404711868, 0.0014236923825040294, 0.0013800362314084368, 0.0012419339757008018, 0.0011884864324004712,
        9.646242791885347E-4, 9.965550662156106E-4, 8.821655173492103E-4, 7.333437120477682E-4, 5.956584591861345E-4,
        5.382434393469185E-4, 3.820013088982811E-4, 3.1072174046886387E-4, 2.490401574043833E-4, 1.7369225659193908E-4,
        1.3261975720080322E-4, 7.486848469149086E-5, 4.7309724158992934E-5, 2.9924237667817037E-5, 2.0999523804124702E-5,
        1.2698952712723992E-5 };
    double p8208_d2x1y1_yerrplus[]      = { 0.00291547594742265, 0.0023482972554597936, 0.0019859506539690254, 0.002015589243868899, 0.0017705366418123065,
        0.001586190404711868, 0.0014236923825040294, 0.0013800362314084368, 0.0012419339757008018, 0.0011884864324004712,
        9.646242791885347E-4, 9.965550662156106E-4, 8.821655173492103E-4, 7.333437120477682E-4, 5.956584591861345E-4,
        5.382434393469185E-4, 3.820013088982811E-4, 3.1072174046886387E-4, 2.490401574043833E-4, 1.7369225659193908E-4,
        1.3261975720080322E-4, 7.486848469149086E-5, 4.7309724158992934E-5, 2.9924237667817037E-5, 2.0999523804124702E-5,
        1.2698952712723992E-5 };
    double p8208_d2x1y1_ystatminus[]    = { 0.00198, 0.00101, 6.4E-4, 5.1E-4, 4.2E-4,
        3.8E-4, 3.5E-4, 3.3E-4, 3.2E-4, 3.0E-4,
        2.9E-4, 2.71E-4, 2.5E-4, 2.28E-4, 2.03E-4,
        1.75E-4, 1.06E-4, 8.8E-5, 7.0E-5, 6.0E-5,
        5.2E-5, 2.52E-5, 1.86E-5, 1.39E-5, 1.13E-5,
        5.45E-6 };
    double p8208_d2x1y1_ystatplus[]     = { 0.00198, 0.00101, 6.4E-4, 5.1E-4, 4.2E-4,
        3.8E-4, 3.5E-4, 3.3E-4, 3.2E-4, 3.0E-4,
        2.9E-4, 2.71E-4, 2.5E-4, 2.28E-4, 2.03E-4,
        1.75E-4, 1.06E-4, 8.8E-5, 7.0E-5, 6.0E-5,
        5.2E-5, 2.52E-5, 1.86E-5, 1.39E-5, 1.13E-5,
        5.45E-6 };
    //     double p8208_d2x1y1_ysystminus[]    = { 0.00214, 0.00212, 0.00188, 0.00195, 0.00172,
    //                                             0.00154, 0.00138, 0.00134, 0.0012, 0.00115,
    //                                             9.2E-4, 9.59E-4, 8.46E-4, 6.97E-4, 5.6E-4,
    //                                             5.09E-4, 3.67E-4, 2.98E-4, 2.39E-4, 1.63E-4,
    //                                             1.22E-4, 7.05E-5, 4.35E-5, 2.65E-5, 1.77E-5,
    //                                             1.147E-5 };
    
    int p8208_d2x1y1_numpoints          = 26;
    double p8208_d2x1y1_ysystminus[p8208_d2x1y1_numpoints];
    double p8208_d2x1y1_ysystplus[p8208_d2x1y1_numpoints];
    
    for (Int_t i = 0; i < p8208_d2x1y1_numpoints; i++){
        Double_t integralPi0            = fitPi0->Integral(p8208_d2x1y1_xval[i]-p8208_d2x1y1_xerrminus[i], p8208_d2x1y1_xval[i]+p8208_d2x1y1_xerrplus[i])/(p8208_d2x1y1_xerrminus[i]+p8208_d2x1y1_xerrplus[i]);
        p8208_d2x1y1_ysystminus[i]      = TMath::Sqrt(p8208_d2x1y1_yerrminus[i]*p8208_d2x1y1_yerrminus[i]-p8208_d2x1y1_ystatminus[i]*p8208_d2x1y1_ystatminus[i]);
        p8208_d2x1y1_ysystplus[i]       = TMath::Sqrt(p8208_d2x1y1_yerrplus[i]*p8208_d2x1y1_yerrplus[i]-p8208_d2x1y1_ystatplus[i]*p8208_d2x1y1_ystatplus[i]);
        
        p8208_d2x1y1_yval[i]            = p8208_d2x1y1_yval[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
        p8208_d2x1y1_yerrminus[i]       = p8208_d2x1y1_yerrminus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
        p8208_d2x1y1_yerrplus[i]        = p8208_d2x1y1_yerrplus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
        p8208_d2x1y1_ystatminus[i]      = p8208_d2x1y1_ystatminus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
        p8208_d2x1y1_ystatplus[i]       = p8208_d2x1y1_ystatplus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
        p8208_d2x1y1_ysystminus[i]      = p8208_d2x1y1_ysystminus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
        p8208_d2x1y1_ysystplus[i]       = p8208_d2x1y1_ysystplus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
    }
    
    TGraphAsymmErrors* p8208_d2x1y1     = NULL;
    if (returnStat) {
        p8208_d2x1y1                    = new TGraphAsymmErrors(p8208_d2x1y1_numpoints, p8208_d2x1y1_xval, p8208_d2x1y1_yval, p8208_d2x1y1_xerrminus, p8208_d2x1y1_xerrplus, p8208_d2x1y1_ystatminus, p8208_d2x1y1_ystatplus);
    } else {
        p8208_d2x1y1                    = new TGraphAsymmErrors(p8208_d2x1y1_numpoints, p8208_d2x1y1_xval, p8208_d2x1y1_yval, p8208_d2x1y1_xerrminus, p8208_d2x1y1_xerrplus, p8208_d2x1y1_ysystminus, p8208_d2x1y1_ysystplus);
    }
    return p8208_d2x1y1;
}

// @article{Adam:2016ich,
//       author         = "Adam, Jaroslav and others",
//       title          = "{$D$-meson production in $p$-Pb collisions at
//                         $\sqrt{s_{\rm NN}}=$5.02 TeV and in pp collisions at
//                         $\sqrt{s}=$7 TeV}",
//       collaboration  = "ALICE",
//       journal        = "Phys. Rev.",
//       volume         = "C94",
//       year           = "2016",
//       number         = "5",
//       pages          = "054908",
//       doi            = "10.1103/PhysRevC.94.054908",
//       eprint         = "1605.07569",
//       archivePrefix  = "arXiv",
//       primaryClass   = "nucl-ex",
//       reportNumber   = "CERN-EP-2016-127",
//       SLACcitation   = "%%CITATION = ARXIV:1605.07569;%%"
// }

//================================================================================================================
// D^0 meson in pp, 7 TeV collisions dsigma/dpt
//================================================================================================================
TGraphAsymmErrors* GetD0Mesonpp7TeV(Bool_t returnStat){
    
    // Plot: p9170_d1x1y1
    double p9170_d1x1y1_xval[]          = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 10.0 };
    double p9170_d1x1y1_xerrminus[]     = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 2.0 };
    double p9170_d1x1y1_xerrplus[]      = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 2.0 };
    double p9170_d1x1y1_yval[]          = { 103.0, 194.0, 113.0, 66.1, 19.9, 11.7, 4.56, 1.06 };
    double p9170_d1x1y1_yerrminus[]     = { 32.76492026543022, 44.51797389819083, 22.615260334561707, 11.66460029319479, 4.739926159762407, 2.8974644087546615, 1.2094895617573558, 0.3245689449100145 };
    double p9170_d1x1y1_yerrplus[]      = { 32.76492026543022, 44.51797389819083, 22.615260334561707, 11.66460029319479, 4.739926159762407, 2.8974644087546615, 1.2094895617573558, 0.3245689449100145 };
    double p9170_d1x1y1_ystatminus[]    = { 27.7, 30.1, 17.1, 7.77, 3.95, 2.17, 0.799, 0.228 };
    double p9170_d1x1y1_ystatplus[]     = { 27.7, 30.1, 17.1, 7.77, 3.95, 2.17, 0.799, 0.228 };
    int p9170_d1x1y1_numpoints          = 8;
    double p9170_d1x1y1_ysystminus[p9170_d1x1y1_numpoints];
    double p9170_d1x1y1_ysystplus[p9170_d1x1y1_numpoints];
    double branching                    = 1; // 0.0388;
    
    for (Int_t i = 0; i < p9170_d1x1y1_numpoints; i++){
        p9170_d1x1y1_ysystminus[i]      = TMath::Sqrt(p9170_d1x1y1_yerrminus[i]*p9170_d1x1y1_yerrminus[i]-p9170_d1x1y1_ystatminus[i]*p9170_d1x1y1_ystatminus[i]);
        p9170_d1x1y1_ysystplus[i]       = TMath::Sqrt(p9170_d1x1y1_yerrplus[i]*p9170_d1x1y1_yerrplus[i]-p9170_d1x1y1_ystatplus[i]*p9170_d1x1y1_ystatplus[i]);
        p9170_d1x1y1_yval[i]            = p9170_d1x1y1_yval[i]/(2*TMath::Pi()*p9170_d1x1y1_xval[i])*xSectionpp7TeVV0OR*recalcBarn/(branching);
        p9170_d1x1y1_yerrminus[i]       = p9170_d1x1y1_yerrminus[i]/(2*TMath::Pi()*p9170_d1x1y1_xval[i])*mutopico/(branching);
        p9170_d1x1y1_yerrplus[i]        = p9170_d1x1y1_yerrplus[i]/(2*TMath::Pi()*p9170_d1x1y1_xval[i])*mutopico/(branching);
        p9170_d1x1y1_ystatminus[i]      = p9170_d1x1y1_ystatminus[i]/(2*TMath::Pi()*p9170_d1x1y1_xval[i])*mutopico/(branching);
        p9170_d1x1y1_ystatplus[i]       = p9170_d1x1y1_ystatplus[i]/(2*TMath::Pi()*p9170_d1x1y1_xval[i])*mutopico/(branching);
        p9170_d1x1y1_ysystminus[i]      = p9170_d1x1y1_ysystminus[i]/(2*TMath::Pi()*p9170_d1x1y1_xval[i])*mutopico/(branching);
        p9170_d1x1y1_ysystplus[i]       = p9170_d1x1y1_ysystplus[i]/(2*TMath::Pi()*p9170_d1x1y1_xval[i])*mutopico/(branching);
    }
    
    TGraphAsymmErrors* p9170_d1x1y1     = NULL;
    if (returnStat) {
        p9170_d1x1y1                    = new TGraphAsymmErrors(p9170_d1x1y1_numpoints, p9170_d1x1y1_xval, p9170_d1x1y1_yval, p9170_d1x1y1_xerrminus, p9170_d1x1y1_xerrplus,p9170_d1x1y1_ystatminus, p9170_d1x1y1_ystatplus);
    } else {
        p9170_d1x1y1                    = new TGraphAsymmErrors(p9170_d1x1y1_numpoints, p9170_d1x1y1_xval, p9170_d1x1y1_yval, p9170_d1x1y1_xerrminus, p9170_d1x1y1_xerrplus,p9170_d1x1y1_ysystminus, p9170_d1x1y1_ysystplus);
    }
    return p9170_d1x1y1;
}

//================================================================================================================
// prompt D^0 meson in pp, 7 TeV collisions dsigma/dpt
//================================================================================================================
TGraphAsymmErrors* GetD0MesonPromptpp7TeV(Bool_t returnStat){
    
    // Plot: p9170_d2x1y1
    double p9170_d2x1y1_xval[]          = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 10.0, 14.0 };
    double p9170_d2x1y1_xerrminus[]     = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2.0, 2.0 };
    double p9170_d2x1y1_xerrplus[]      = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2.0, 2.0 };
    double p9170_d2x1y1_yval[]          = { 99.0, 187.0, 114.0, 59.3, 28.9, 12.4, 6.32, 3.05, 1.22, 0.213 };
    double p9170_d2x1y1_yerrminus[]     = { 34.52028389222777, 61.12487218800543, 34.69135338956957, 13.215676297488525, 6.113084327898642,
        2.5134040662018515, 1.2737272078431865, 0.6999957142725947, 0.2448999795835026, 0.06226371334894828 };
    double p9170_d2x1y1_yerrplus[]      = { 32.39645042284725, 43.713842201298206, 22.330472453577872, 9.494472075897638, 4.6558028308767545,
        2.139088590965788, 1.1578091379843225, 0.6770561276585568, 0.22965191050805564, 0.06172001296176143 };
    double p9170_d2x1y1_ystatminus[]    = { 27.7, 30.1, 10.7, 4.29, 2.13, 1.14, 0.691, 0.463, 0.126, 0.0494 };
    double p9170_d2x1y1_ystatplus[]     = { 27.7, 30.1, 10.7, 4.29, 2.13, 1.14, 0.691, 0.463, 0.126, 0.0494 };
    int p9170_d2x1y1_numpoints          = 10;
    double p9170_d2x1y1_ysystminus[p9170_d2x1y1_numpoints];
    double p9170_d2x1y1_ysystplus[p9170_d2x1y1_numpoints];
    double branching                    = 1; // 0.0388;
    
    for (Int_t i = 0; i < p9170_d2x1y1_numpoints; i++){
        p9170_d2x1y1_ysystminus[i]      = TMath::Sqrt(p9170_d2x1y1_yerrminus[i]*p9170_d2x1y1_yerrminus[i]-p9170_d2x1y1_ystatminus[i]*p9170_d2x1y1_ystatminus[i]);
        p9170_d2x1y1_ysystplus[i]       = TMath::Sqrt(p9170_d2x1y1_yerrplus[i]*p9170_d2x1y1_yerrplus[i]-p9170_d2x1y1_ystatplus[i]*p9170_d2x1y1_ystatplus[i]);
        p9170_d2x1y1_yval[i]            = p9170_d2x1y1_yval[i]/(2*TMath::Pi()*p9170_d2x1y1_xval[i])*xSectionpp7TeVV0OR*recalcBarn/(branching);
        p9170_d2x1y1_yerrminus[i]       = p9170_d2x1y1_yerrminus[i]/(2*TMath::Pi()*p9170_d2x1y1_xval[i])*mutopico/(branching);
        p9170_d2x1y1_yerrplus[i]        = p9170_d2x1y1_yerrplus[i]/(2*TMath::Pi()*p9170_d2x1y1_xval[i])*mutopico/(branching);
        p9170_d2x1y1_ystatminus[i]      = p9170_d2x1y1_ystatminus[i]/(2*TMath::Pi()*p9170_d2x1y1_xval[i])*mutopico/(branching);
        p9170_d2x1y1_ystatplus[i]       = p9170_d2x1y1_ystatplus[i]/(2*TMath::Pi()*p9170_d2x1y1_xval[i])*mutopico/(branching);
        p9170_d2x1y1_ysystminus[i]      = p9170_d2x1y1_ysystminus[i]/(2*TMath::Pi()*p9170_d2x1y1_xval[i])*mutopico/(branching);
        p9170_d2x1y1_ysystplus[i]       = p9170_d2x1y1_ysystplus[i]/(2*TMath::Pi()*p9170_d2x1y1_xval[i])*mutopico/(branching);
    }
    
    TGraphAsymmErrors* p9170_d2x1y1     = NULL;
    if (returnStat) {
        p9170_d2x1y1                    = new TGraphAsymmErrors(p9170_d2x1y1_numpoints, p9170_d2x1y1_xval, p9170_d2x1y1_yval, p9170_d2x1y1_xerrminus, p9170_d2x1y1_xerrplus,p9170_d2x1y1_ystatminus, p9170_d2x1y1_ystatplus);
    } else {
        p9170_d2x1y1                    = new TGraphAsymmErrors(p9170_d2x1y1_numpoints, p9170_d2x1y1_xval, p9170_d2x1y1_yval, p9170_d2x1y1_xerrminus, p9170_d2x1y1_xerrplus,p9170_d2x1y1_ysystminus, p9170_d2x1y1_ysystplus);
    }
    return p9170_d2x1y1;
}

//================================================================================================================
// dN_ch/dEta in pp, 0.9 TeV, http://aliceinfo.cern.ch/ArtSubmission/node/1885
//================================================================================================================
TGraphAsymmErrors* GetChargedParticlePseudorapidtypp900GeV() {
    
    // http://www.hepdata.net/record/ins1394854?version=1&table=Table1
    Int_t       nPoints         = 20;
    Double_t    xVal[]          = {-1.90, -1.70, -1.50, -1.30, -1.10, -0.90, -0.70, -0.50, -0.30, -0.10, 0.10, 0.30, 0.50, 0.70, 0.90, 1.10, 1.30, 1.50, 1.70, 1.90 };
    Double_t    xErrPlus[]      = { 0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    Double_t    xErrMinus[]     = { 0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    Double_t    yVal[]          = { 3.08,  3.16,  3.15,  3.16,  3.13,  3.08,  3.02,  2.97,  2.93,  2.90, 2.91, 2.95, 2.99, 3.02, 3.10, 3.14, 3.15, 3.17, 3.16, 3.07 };
    Double_t    yErrPlus[]      = { 0.14,  0.14,  0.13,  0.13,  0.13,  0.12,  0.11,  0.11,  0.11,  0.11, 0.11, 0.11, 0.11, 0.12, 0.12, 0.12, 0.13, 0.13, 0.13, 0.13 };
    Double_t    yErrMinus[]     = { 0.06,  0.06,  0.06,  0.06,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.06, 0.06, 0.06 };
    
    TGraphAsymmErrors* dNchdEta = new TGraphAsymmErrors(nPoints, xVal, yVal, xErrMinus, xErrPlus, yErrMinus, yErrPlus);
    return dNchdEta;
}

//================================================================================================================
// dN_ch/dEta in pp, 2.76 TeV, http://aliceinfo.cern.ch/ArtSubmission/node/1885
//================================================================================================================
TGraphAsymmErrors* GetChargedParticlePseudorapidtypp2760GeV() {
    
    // http://www.hepdata.net/record/ins1394854?version=1&table=Table4
    Int_t       nPoints         = 20;
    Double_t    xVal[]          = {-1.90, -1.70, -1.50, -1.30, -1.10, -0.90, -0.70, -0.50, -0.30, -0.10, 0.10, 0.30, 0.50, 0.70, 0.90, 1.10, 1.30, 1.50, 1.70, 1.90 };
    Double_t    xErrPlus[]      = { 0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    Double_t    xErrMinus[]     = { 0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    Double_t    yVal[]          = { 4.00,  4.08,  4.08,  4.06,  4.01,  3.94,  3.86,  3.80,  3.72,  3.71, 3.71, 3.74, 3.81, 3.87, 3.95, 4.02, 4.06, 4.08, 4.09, 4.02 };
    Double_t    yErrPlus[]      = { 0.29,  0.29,  0.29,  0.29,  0.28,  0.27,  0.27,  0.26,  0.26,  0.26, 0.26, 0.26, 0.26, 0.27, 0.27, 0.28, 0.28, 0.28, 0.28, 0.29 };
    Double_t    yErrMinus[]     = { 0.17,  0.18,  0.17,  0.17,  0.17,  0.17,  0.16,  0.16,  0.16,  0.16, 0.16, 0.16, 0.16, 0.16, 0.17, 0.17, 0.17, 0.17, 0.18, 0.17 };
    
    TGraphAsymmErrors* dNchdEta = new TGraphAsymmErrors(nPoints, xVal, yVal, xErrMinus, xErrPlus, yErrMinus, yErrPlus);
    return dNchdEta;
}

//================================================================================================================
// dN_ch/dEta in pp, 7 TeV, http://aliceinfo.cern.ch/ArtSubmission/node/1885
//================================================================================================================
TGraphAsymmErrors* GetChargedParticlePseudorapidtypp7TeV() {
    
    // http://www.hepdata.net/record/ins1394854?version=1&table=Table7
    Int_t       nPoints         = 20;
    Double_t    xVal[]          = {-1.90, -1.70, -1.50, -1.30, -1.10, -0.90, -0.70, -0.50, -0.30, -0.10, 0.10, 0.30, 0.50, 0.70, 0.90, 1.10, 1.30, 1.50, 1.70, 1.90 };
    Double_t    xErrPlus[]      = { 0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    Double_t    xErrMinus[]     = { 0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    Double_t    yVal[]          = { 4.94,  5.00,  5.01,  4.98,  4.92,  4.84,  4.74,  4.66,  4.59,  4.55, 4.55, 4.59, 4.66, 4.74, 4.84, 4.91, 4.97, 5.00, 4.99, 4.94 };
    Double_t    yErrPlus[]      = { 0.38,  0.38,  0.38,  0.37,  0.37,  0.36,  0.35,  0.35,  0.34,  0.34, 0.34, 0.34, 0.35, 0.35, 0.36, 0.37, 0.37, 0.38, 0.38, 0.38 };
    Double_t    yErrMinus[]     = { 0.19,  0.19,  0.19,  0.19,  0.18,  0.18,  0.18,  0.17,  0.17,  0.17, 0.17, 0.17, 0.17, 0.18, 0.18, 0.18, 0.19, 0.19, 0.19, 0.20 };
    
    TGraphAsymmErrors* dNchdEta = new TGraphAsymmErrors(nPoints, xVal, yVal, xErrMinus, xErrPlus, yErrMinus, yErrPlus);
    return dNchdEta;
}

//================================================================================================================
// dN_ch/dEta in pp, 8 TeV, http://aliceinfo.cern.ch/ArtSubmission/node/1885
//================================================================================================================
TGraphAsymmErrors* GetChargedParticlePseudorapidtypp8TeV() {
    
    // http://www.hepdata.net/record/ins1394854?version=1&table=Table10
    Int_t       nPoints         = 20;
    Double_t    xVal[]          = {-1.90, -1.70, -1.50, -1.30, -1.10, -0.90, -0.70, -0.50, -0.30, -0.10, 0.10, 0.30, 0.50, 0.70, 0.90, 1.10, 1.30, 1.50, 1.70, 1.90 };
    Double_t    xErrPlus[]      = { 0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    Double_t    xErrMinus[]     = { 0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    Double_t    yVal[]          = { 5.01,  5.05,  5.07,  5.02,  4.96,  4.87,  4.79,  4.69,  4.61,  4.59, 4.59, 4.62, 4.70, 4.77, 4.87, 4.95, 5.03, 5.06, 5.05, 5.01 };
    Double_t    yErrPlus[]      = { 0.40,  0.39,  0.38,  0.38,  0.37,  0.36,  0.36,  0.35,  0.34,  0.34, 0.34, 0.34, 0.35, 0.35, 0.36, 0.37, 0.38, 0.38, 0.38, 0.40 };
    Double_t    yErrMinus[]     = { 0.21,  0.19,  0.19,  0.19,  0.18,  0.18,  0.18,  0.17,  0.17,  0.17, 0.17, 0.17, 0.17, 0.18, 0.18, 0.18, 0.19, 0.19, 0.19, 0.22 };
    
    TGraphAsymmErrors* dNchdEta = new TGraphAsymmErrors(nPoints, xVal, yVal, xErrMinus, xErrPlus, yErrMinus, yErrPlus);
    return dNchdEta;
}






