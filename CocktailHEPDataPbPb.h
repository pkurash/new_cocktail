// @article{Adam:2015sza,
//       author         = "Adam, Jaroslav and others",
//       title          = "{Transverse momentum dependence of D-meson production in
//                         Pb-Pb collisions at $ \sqrt{{\mathrm{s}}_{\mathrm{NN}}}=$
//                         2.76  TeV}",
//       collaboration  = "ALICE",
//       journal        = "JHEP",
//       volume         = "03",
//       year           = "2016",
//       pages          = "081",
//       doi            = "10.1007/JHEP03(2016)081",
//       eprint         = "1509.06888",
//       archivePrefix  = "arXiv",
//       primaryClass   = "nucl-ex",
//       reportNumber   = "CERN-PH-EP-2015-252",
//       SLACcitation   = "%%CITATION = ARXIV:1509.06888;%%"
// }

//================================================================================================================
// D0 meson in 0-10% PbPb collisions dN/dy
//================================================================================================================
TGraphAsymmErrors* GetD0MesonPbPb0010(Bool_t returnStat){
    
    // Plot: p9059_d1x1y1
    double p9059_d1x1y1_xval[]          = { 1.5, 2.5, 3.5, 4.5, 5.5,
        7.0, 10.0, 14.0, 20.0 };
    double p9059_d1x1y1_xerrminus[]     = { 0.5, 0.5, 0.5, 0.5, 0.5,
        1.0, 2.0, 2.0, 4.0 };
    double p9059_d1x1y1_xerrplus[]      = { 0.5, 0.5, 0.5, 0.5, 0.5,
        1.0, 2.0, 2.0, 4.0 };
    double p9059_d1x1y1_yval[]          = { 2.65, 1.03, 0.262, 0.0744, 0.0226,
        0.00646, 0.00152, 2.58E-4, 7.16E-5 };
    double p9059_d1x1y1_yerrminus[]     = { 1.4512649654697793, 0.37099897573982593, 0.06732258164984464, 0.019246373684411306, 0.005212341124677087,
        0.001464650128870373, 3.255733404319217E-4, 7.138879463893476E-5, 3.3182224156918716E-5 };
    double p9059_d1x1y1_yerrplus[]      = { 0.8681992858785361, 0.25664894311101305, 0.045448652345256614, 0.014029942979214135, 0.00453230625620114,
        0.001375988372043892, 3.132379287378845E-4, 7.339427770609913E-5, 3.2840523747346054E-5 };
    double p9059_d1x1y1_ystatminus[]    = { 0.387, 0.0682, 0.0162, 0.00577, 0.0022,
        5.56E-4, 1.25E-4, 3.8E-5, 1.84E-5 };
    double p9059_d1x1y1_ystatplus[]     = { 0.387, 0.0682, 0.0162, 0.00577, 0.0022,
        5.56E-4, 1.25E-4, 3.8E-5, 1.84E-5 };
    int p9059_d1x1y1_numpoints          = 9;
    
    double p9059_d1x1y1_ysystminus[p9059_d1x1y1_numpoints];
    double p9059_d1x1y1_ysystplus[p9059_d1x1y1_numpoints];
    double branching                    = 1; //0.0388;
    
    for (Int_t i = 0; i < p9059_d1x1y1_numpoints; i++){
        p9059_d1x1y1_ysystminus[i]      = TMath::Sqrt(p9059_d1x1y1_yerrminus[i]*p9059_d1x1y1_yerrminus[i]-p9059_d1x1y1_ystatminus[i]*p9059_d1x1y1_ystatminus[i]);
        p9059_d1x1y1_ysystplus[i]       = TMath::Sqrt(p9059_d1x1y1_yerrplus[i]*p9059_d1x1y1_yerrplus[i]-p9059_d1x1y1_ystatplus[i]*p9059_d1x1y1_ystatplus[i]);
        
        p9059_d1x1y1_yval[i]            = p9059_d1x1y1_yval[i]/(branching);
        p9059_d1x1y1_yerrminus[i]       = p9059_d1x1y1_yerrminus[i]/(branching);
        p9059_d1x1y1_yerrplus[i]        = p9059_d1x1y1_yerrplus[i]/(branching);
        p9059_d1x1y1_ystatminus[i]      = p9059_d1x1y1_ystatminus[i]/(branching);
        p9059_d1x1y1_ystatplus[i]       = p9059_d1x1y1_ystatplus[i]/(branching);
        p9059_d1x1y1_ysystminus[i]      = p9059_d1x1y1_ysystminus[i]/(branching);
        p9059_d1x1y1_ysystplus[i]       = p9059_d1x1y1_ysystplus[i]/(branching);
    }
    
    TGraphAsymmErrors* p9059_d1x1y1     = NULL;
    if (returnStat) {
        p9059_d1x1y1                    = new TGraphAsymmErrors(p9059_d1x1y1_numpoints, p9059_d1x1y1_xval, p9059_d1x1y1_yval, p9059_d1x1y1_xerrminus, p9059_d1x1y1_xerrplus,p9059_d1x1y1_ystatminus, p9059_d1x1y1_ystatplus);
    } else {
        p9059_d1x1y1                    = new TGraphAsymmErrors(p9059_d1x1y1_numpoints, p9059_d1x1y1_xval, p9059_d1x1y1_yval, p9059_d1x1y1_xerrminus, p9059_d1x1y1_xerrplus,p9059_d1x1y1_ysystminus, p9059_d1x1y1_ysystplus);
    }
    return p9059_d1x1y1;
}

//================================================================================================================
// DPlus meson in 0-10% PbPb collisions dN/dy
//================================================================================================================
TGraphAsymmErrors* GetDPlusMesonPbPb0010(Bool_t returnStat){
    
    // Plot: p9059_d2x1y1
    double p9059_d2x1y1_xval[]          = { 3.5, 4.5, 5.5, 7.0, 10.0,
        14.0, 20.0, 30.0 };
    double p9059_d2x1y1_xerrminus[]     = { 0.5, 0.5, 0.5, 1.0, 2.0,
        2.0, 4.0, 6.0 };
    double p9059_d2x1y1_xerrplus[]      = { 0.5, 0.5, 0.5, 1.0, 2.0,
        2.0, 4.0, 6.0 };
    double p9059_d2x1y1_yval[]          = { 0.0861, 0.0301, 0.0105, 0.00295, 6.2E-4,
        1.14E-4, 3.16E-5, 6.49E-6 };
    double p9059_d2x1y1_yerrminus[]     = { 0.02883851591188423, 0.00877494159524723, 0.002970488175367813, 9.661785549265726E-4, 2.0411271395971394E-4,
        4.3373840042126774E-5, 1.2939319920304931E-5, 2.8039436513596345E-6 };
    double p9059_d2x1y1_yerrplus[]      = { 0.0271504622428422, 0.008155593172786392, 0.002860004370626031, 9.344137199335207E-4, 1.924027286708793E-4,
        4.427335541835518E-5, 1.3082262801212947E-5, 2.7458332068791067E-6 };
    double p9059_d2x1y1_ystatminus[]    = { 0.0181, 0.00486, 0.00165, 6.43E-4, 1.23E-4,
        2.98E-5, 9.84E-6, 2.14E-6 };
    double p9059_d2x1y1_ystatplus[]     = { 0.0181, 0.00486, 0.00165, 6.43E-4, 1.23E-4,
        2.98E-5, 9.84E-6, 2.14E-6 };
    int p9059_d2x1y1_numpoints          = 8;
    
    double p9059_d2x1y1_ysystminus[p9059_d2x1y1_numpoints];
    double p9059_d2x1y1_ysystplus[p9059_d2x1y1_numpoints];
    double branching                    = 1; //0.0913;
    
    for (Int_t i = 0; i < p9059_d2x1y1_numpoints; i++){
        p9059_d2x1y1_ysystminus[i]      = TMath::Sqrt(p9059_d2x1y1_yerrminus[i]*p9059_d2x1y1_yerrminus[i]-p9059_d2x1y1_ystatminus[i]*p9059_d2x1y1_ystatminus[i]);
        p9059_d2x1y1_ysystplus[i]       = TMath::Sqrt(p9059_d2x1y1_yerrplus[i]*p9059_d2x1y1_yerrplus[i]-p9059_d2x1y1_ystatplus[i]*p9059_d2x1y1_ystatplus[i]);
        
        p9059_d2x1y1_yval[i]            = p9059_d2x1y1_yval[i]/(branching);
        p9059_d2x1y1_yerrminus[i]       = p9059_d2x1y1_yerrminus[i]/(branching);
        p9059_d2x1y1_yerrplus[i]        = p9059_d2x1y1_yerrplus[i]/(branching);
        p9059_d2x1y1_ystatminus[i]      = p9059_d2x1y1_ystatminus[i]/(branching);
        p9059_d2x1y1_ystatplus[i]       = p9059_d2x1y1_ystatplus[i]/(branching);
        p9059_d2x1y1_ysystminus[i]      = p9059_d2x1y1_ysystminus[i]/(branching);
        p9059_d2x1y1_ysystplus[i]       = p9059_d2x1y1_ysystplus[i]/(branching);
    }
    
    TGraphAsymmErrors* p9059_d2x1y1     = NULL;
    if (returnStat) {
        p9059_d2x1y1                    = new TGraphAsymmErrors(p9059_d2x1y1_numpoints, p9059_d2x1y1_xval, p9059_d2x1y1_yval, p9059_d2x1y1_xerrminus, p9059_d2x1y1_xerrplus,p9059_d2x1y1_ystatminus, p9059_d2x1y1_ystatplus);
    } else {
        p9059_d2x1y1                    = new TGraphAsymmErrors(p9059_d2x1y1_numpoints, p9059_d2x1y1_xval, p9059_d2x1y1_yval, p9059_d2x1y1_xerrminus, p9059_d2x1y1_xerrplus,p9059_d2x1y1_ysystminus, p9059_d2x1y1_ysystplus);
    }
    return p9059_d2x1y1;
}

//================================================================================================================
// D*+ meson in 0-10% PbPb collisions dN/dy
//================================================================================================================
TGraphAsymmErrors* GetDStarPlusMesonPbPb0010(Bool_t returnStat){
    
    // Plot: p9059_d3x1y1
    double p9059_d3x1y1_xval[]          = { 3.5, 4.5, 5.5, 7.0, 10.0,
        14.0, 20.0, 30.0 };
    double p9059_d3x1y1_xerrminus[]     = { 0.5, 0.5, 0.5, 1.0, 2.0,
        2.0, 4.0, 6.0 };
    double p9059_d3x1y1_xerrplus[]      = { 0.5, 0.5, 0.5, 1.0, 2.0,
        2.0, 4.0, 6.0 };
    double p9059_d3x1y1_yval[]          = { 0.179, 0.0347, 0.0153, 0.00474, 7.92E-4,
        2.02E-4, 3.35E-5, 6.13E-6 };
    double p9059_d3x1y1_yerrminus[]     = { 0.07388321866296839, 0.0164251179600026, 0.0043068550010419435, 0.0012138027022543655, 2.2145753543286804E-4,
        6.434423051059047E-5, 1.268007097771933E-5, 2.3890927148187448E-6 };
    double p9059_d3x1y1_yerrplus[]      = { 0.06921972262296347, 0.01618034919277084, 0.004226310447659992, 0.0011897230770225482, 2.168392261561547E-4,
        6.214225615472937E-5, 1.271031864274063E-5, 2.35617147083993E-6 };
    double p9059_d3x1y1_ystatminus[]    = { 0.0513, 0.0142, 0.0028, 6.76E-4, 1.41E-4,
        4.49E-5, 1.02E-5, 1.91E-6 };
    double p9059_d3x1y1_ystatplus[]     = { 0.0513, 0.0142, 0.0028, 6.76E-4, 1.41E-4,
        4.49E-5, 1.02E-5, 1.91E-6 };
    int p9059_d3x1y1_numpoints          = 8;
    
    double p9059_d3x1y1_ysystminus[p9059_d3x1y1_numpoints];
    double p9059_d3x1y1_ysystplus[p9059_d3x1y1_numpoints];
    double branching                    = 1; //0.0388*0.677;
    
    for (Int_t i = 0; i < p9059_d3x1y1_numpoints; i++){
        p9059_d3x1y1_ysystminus[i]      = TMath::Sqrt(p9059_d3x1y1_yerrminus[i]*p9059_d3x1y1_yerrminus[i]-p9059_d3x1y1_ystatminus[i]*p9059_d3x1y1_ystatminus[i]);
        p9059_d3x1y1_ysystplus[i]       = TMath::Sqrt(p9059_d3x1y1_yerrplus[i]*p9059_d3x1y1_yerrplus[i]-p9059_d3x1y1_ystatplus[i]*p9059_d3x1y1_ystatplus[i]);
        
        p9059_d3x1y1_yval[i]            = p9059_d3x1y1_yval[i]/(branching);
        p9059_d3x1y1_yerrminus[i]       = p9059_d3x1y1_yerrminus[i]/(branching);
        p9059_d3x1y1_yerrplus[i]        = p9059_d3x1y1_yerrplus[i]/(branching);
        p9059_d3x1y1_ystatminus[i]      = p9059_d3x1y1_ystatminus[i]/(branching);
        p9059_d3x1y1_ystatplus[i]       = p9059_d3x1y1_ystatplus[i]/(branching);
        p9059_d3x1y1_ysystminus[i]      = p9059_d3x1y1_ysystminus[i]/(branching);
        p9059_d3x1y1_ysystplus[i]       = p9059_d3x1y1_ysystplus[i]/(branching);
    }
    
    TGraphAsymmErrors* p9059_d3x1y1     = NULL;
    if (returnStat) {
        p9059_d3x1y1                    = new TGraphAsymmErrors(p9059_d3x1y1_numpoints, p9059_d3x1y1_xval, p9059_d3x1y1_yval, p9059_d3x1y1_xerrminus, p9059_d3x1y1_xerrplus,p9059_d3x1y1_ystatminus, p9059_d3x1y1_ystatplus);
    } else {
        p9059_d3x1y1                    = new TGraphAsymmErrors(p9059_d3x1y1_numpoints, p9059_d3x1y1_xval, p9059_d3x1y1_yval, p9059_d3x1y1_xerrminus, p9059_d3x1y1_xerrplus,p9059_d3x1y1_ysystminus, p9059_d3x1y1_ysystplus);
    }
    return p9059_d3x1y1;
}

//================================================================================================================
// D0 meson in 30-50% PbPb collisions dN/dy
//================================================================================================================
TGraphAsymmErrors* GetD0MesonPbPb3050(Bool_t returnStat){
    
    // Plot: p9059_d1x1y1
    double p9059_d4x1y1_xval[]          = { 1.5, 2.5, 3.5, 4.5, 5.5,
        7.0, 10.0, 14.0 };
    double p9059_d4x1y1_xerrminus[]     = { 0.5, 0.5, 0.5, 0.5, 0.5,
        1.0, 2.0, 2.0 };
    double p9059_d4x1y1_xerrplus[]      = { 0.5, 0.5, 0.5, 0.5, 0.5,
        1.0, 2.0, 2.0 };
    double p9059_d4x1y1_yval[]          = { 0.539, 0.204, 0.0585, 0.019, 0.00793,
        0.00238, 5.23E-4, 1.14E-4 };
    double p9059_d4x1y1_yerrminus[]     = { 0.2844819853699, 0.06930036074942178, 0.01740109766652667, 0.005316822359266858, 0.002033447565097266,
        6.654990608558363E-4, 1.5211945963616885E-4, 3.9297582622853533E-5 };
    double p9059_d4x1y1_yerrplus[]      = { 0.15366053494635504, 0.04349505719044407, 0.012196237124621676, 0.0041628235609979915, 0.001813943218515949,
        6.491910350582485E-4, 1.489990939569768E-4, 3.97940950393397E-5 };
    double p9059_d4x1y1_ystatminus[]    = { 0.109, 0.0193, 0.00519, 0.00189, 8.53E-4,
        3.02E-4, 7.42E-5, 2.49E-5 };
    double p9059_d4x1y1_ystatplus[]     = { 0.109, 0.0193, 0.00519, 0.00189, 8.53E-4,
        3.02E-4, 7.42E-5, 2.49E-5 };
    int p9059_d4x1y1_numpoints          = 8;
    double p9059_d4x1y1_ysystminus[p9059_d4x1y1_numpoints];
    double p9059_d4x1y1_ysystplus[p9059_d4x1y1_numpoints];
    double branching                    = 1; //0.0388;
    
    for (Int_t i = 0; i < p9059_d4x1y1_numpoints; i++){
        p9059_d4x1y1_ysystminus[i]      = TMath::Sqrt(p9059_d4x1y1_yerrminus[i]*p9059_d4x1y1_yerrminus[i]-p9059_d4x1y1_ystatminus[i]*p9059_d4x1y1_ystatminus[i]);
        p9059_d4x1y1_ysystplus[i]       = TMath::Sqrt(p9059_d4x1y1_yerrplus[i]*p9059_d4x1y1_yerrplus[i]-p9059_d4x1y1_ystatplus[i]*p9059_d4x1y1_ystatplus[i]);
        
        p9059_d4x1y1_yval[i]            = p9059_d4x1y1_yval[i]/(branching);
        p9059_d4x1y1_yerrminus[i]       = p9059_d4x1y1_yerrminus[i]/(branching);
        p9059_d4x1y1_yerrplus[i]        = p9059_d4x1y1_yerrplus[i]/(branching);
        p9059_d4x1y1_ystatminus[i]      = p9059_d4x1y1_ystatminus[i]/(branching);
        p9059_d4x1y1_ystatplus[i]       = p9059_d4x1y1_ystatplus[i]/(branching);
        p9059_d4x1y1_ysystminus[i]      = p9059_d4x1y1_ysystminus[i]/(branching);
        p9059_d4x1y1_ysystplus[i]       = p9059_d4x1y1_ysystplus[i]/(branching);
    }
    
    TGraphAsymmErrors* p9059_d4x1y1     = NULL;
    if (returnStat) {
        p9059_d4x1y1                    = new TGraphAsymmErrors(p9059_d4x1y1_numpoints, p9059_d4x1y1_xval, p9059_d4x1y1_yval, p9059_d4x1y1_xerrminus, p9059_d4x1y1_xerrplus,p9059_d4x1y1_ystatminus, p9059_d4x1y1_ystatplus);
    } else {
        p9059_d4x1y1                    = new TGraphAsymmErrors(p9059_d4x1y1_numpoints, p9059_d4x1y1_xval, p9059_d4x1y1_yval, p9059_d4x1y1_xerrminus, p9059_d4x1y1_xerrplus,p9059_d4x1y1_ysystminus, p9059_d4x1y1_ysystplus);
    }
    return p9059_d4x1y1;
}

//================================================================================================================
// DPlus meson in 30-50% PbPb collisions dN/dy
//================================================================================================================
TGraphAsymmErrors* GetDPlusMesonPbPb3050(Bool_t returnStat){
    
    // Plot: p9059_d5x1y1
    double p9059_d5x1y1_xval[]          = { 2.5, 3.5, 4.5, 5.5, 7.0,
        10.0, 14.0 };
    double p9059_d5x1y1_xerrminus[]     = { 0.5, 0.5, 0.5, 0.5, 1.0,
        2.0, 2.0 };
    double p9059_d5x1y1_xerrplus[]      = { 0.5, 0.5, 0.5, 0.5, 1.0,
        2.0, 2.0 };
    double p9059_d5x1y1_yval[]          = { 0.0737, 0.028, 0.00998, 0.00381, 0.00106,
        2.25E-4, 3.85E-5 };
    double p9059_d5x1y1_yerrminus[]     = { 0.03259831283977746, 0.009164829512871476, 0.003259048941025587, 0.0012885674215965574, 3.4863447907514827E-4,
        6.640143070747799E-5, 1.878333303756285E-5 };
    double p9059_d5x1y1_yerrplus[]      = { 0.029043181643890185, 0.008429549216891732, 0.0029993872707604796, 0.001244991967845576, 3.340553846295551E-4,
        6.31795853104466E-5, 1.8992264214674354E-5 };
    double p9059_d5x1y1_ystatminus[]    = { 0.0214, 0.00526, 0.00197, 8.41E-4, 2.19E-4,
        3.65E-5, 1.47E-5 };
    double p9059_d5x1y1_ystatplus[]     = { 0.0214, 0.00526, 0.00197, 8.41E-4, 2.19E-4,
        3.65E-5, 1.47E-5 };
    int p9059_d5x1y1_numpoints          = 7;
    double p9059_d5x1y1_ysystminus[p9059_d5x1y1_numpoints];
    double p9059_d5x1y1_ysystplus[p9059_d5x1y1_numpoints];
    double branching                    = 1; //0.0913;
    
    for (Int_t i = 0; i < p9059_d5x1y1_numpoints; i++){
        p9059_d5x1y1_ysystminus[i]      = TMath::Sqrt(p9059_d5x1y1_yerrminus[i]*p9059_d5x1y1_yerrminus[i]-p9059_d5x1y1_ystatminus[i]*p9059_d5x1y1_ystatminus[i]);
        p9059_d5x1y1_ysystplus[i]       = TMath::Sqrt(p9059_d5x1y1_yerrplus[i]*p9059_d5x1y1_yerrplus[i]-p9059_d5x1y1_ystatplus[i]*p9059_d5x1y1_ystatplus[i]);
        
        p9059_d5x1y1_yval[i]            = p9059_d5x1y1_yval[i]/(branching);
        p9059_d5x1y1_yerrminus[i]       = p9059_d5x1y1_yerrminus[i]/(branching);
        p9059_d5x1y1_yerrplus[i]        = p9059_d5x1y1_yerrplus[i]/(branching);
        p9059_d5x1y1_ystatminus[i]      = p9059_d5x1y1_ystatminus[i]/(branching);
        p9059_d5x1y1_ystatplus[i]       = p9059_d5x1y1_ystatplus[i]/(branching);
        p9059_d5x1y1_ysystminus[i]      = p9059_d5x1y1_ysystminus[i]/(branching);
        p9059_d5x1y1_ysystplus[i]       = p9059_d5x1y1_ysystplus[i]/(branching);
    }
    
    TGraphAsymmErrors* p9059_d5x1y1     = NULL;
    if (returnStat) {
        p9059_d5x1y1                    = new TGraphAsymmErrors(p9059_d5x1y1_numpoints, p9059_d5x1y1_xval, p9059_d5x1y1_yval, p9059_d5x1y1_xerrminus, p9059_d5x1y1_xerrplus,p9059_d5x1y1_ystatminus, p9059_d5x1y1_ystatplus);
    } else {
        p9059_d5x1y1                    = new TGraphAsymmErrors(p9059_d5x1y1_numpoints, p9059_d5x1y1_xval, p9059_d5x1y1_yval, p9059_d5x1y1_xerrminus, p9059_d5x1y1_xerrplus,p9059_d5x1y1_ysystminus, p9059_d5x1y1_ysystplus);
    }
    return p9059_d5x1y1;
}

//================================================================================================================
// D*+ meson in 30-50% PbPb collisions dN/dy
//================================================================================================================
TGraphAsymmErrors* GetDStarPlusMesonPbPb3050(Bool_t returnStat){
    
    // Plot: p9059_d6x1y1
    double p9059_d6x1y1_xval[]          = { 2.5, 3.5, 4.5, 5.5, 7.0,
        10.0, 14.0 };
    double p9059_d6x1y1_xerrminus[]     = { 0.5, 0.5, 0.5, 0.5, 1.0,
        2.0, 2.0 };
    double p9059_d6x1y1_xerrplus[]      = { 0.5, 0.5, 0.5, 0.5, 1.0,
        2.0, 2.0 };
    double p9059_d6x1y1_yval[]          = { 0.0774, 0.0281, 0.00926, 0.00429, 0.00141,
        3.62E-4, 9.21E-5 };
    double p9059_d6x1y1_yerrminus[]     = { 0.03359062369173874, 0.009955506014261655, 0.003104323436757194, 0.0013598441087124657, 4.3925049800768584E-4,
        9.472544536712404E-5, 2.7733979159146996E-5 };
    double p9059_d6x1y1_yerrplus[]      = { 0.030242281659954164, 0.009215139716792144, 0.0030235285677499392, 0.0013310176557807189, 4.3141163637528367E-4,
        9.24118498895028E-5, 2.6855317909121835E-5 };
    double p9059_d6x1y1_ystatminus[]    = { 0.02, 0.00696, 0.00225, 9.46E-4, 3.32E-4,
        6.01E-5, 2.03E-5 };
    double p9059_d6x1y1_ystatplus[]     = { 0.02, 0.00696, 0.00225, 9.46E-4, 3.32E-4,
        6.01E-5, 2.03E-5 };
    int p9059_d6x1y1_numpoints          = 7;
    
    double p9059_d6x1y1_ysystminus[p9059_d6x1y1_numpoints];
    double p9059_d6x1y1_ysystplus[p9059_d6x1y1_numpoints];
    double branching                    = 1; //0.0388*0.677;
    
    for (Int_t i = 0; i < p9059_d6x1y1_numpoints; i++){
        p9059_d6x1y1_ysystminus[i]      = TMath::Sqrt(p9059_d6x1y1_yerrminus[i]*p9059_d6x1y1_yerrminus[i]-p9059_d6x1y1_ystatminus[i]*p9059_d6x1y1_ystatminus[i]);
        p9059_d6x1y1_ysystplus[i]       = TMath::Sqrt(p9059_d6x1y1_yerrplus[i]*p9059_d6x1y1_yerrplus[i]-p9059_d6x1y1_ystatplus[i]*p9059_d6x1y1_ystatplus[i]);
        
        p9059_d6x1y1_yval[i]            = p9059_d6x1y1_yval[i]/(branching);
        p9059_d6x1y1_yerrminus[i]       = p9059_d6x1y1_yerrminus[i]/(branching);
        p9059_d6x1y1_yerrplus[i]        = p9059_d6x1y1_yerrplus[i]/(branching);
        p9059_d6x1y1_ystatminus[i]      = p9059_d6x1y1_ystatminus[i]/(branching);
        p9059_d6x1y1_ystatplus[i]       = p9059_d6x1y1_ystatplus[i]/(branching);
        p9059_d6x1y1_ysystminus[i]      = p9059_d6x1y1_ysystminus[i]/(branching);
        p9059_d6x1y1_ysystplus[i]       = p9059_d6x1y1_ysystplus[i]/(branching);
    }
    
    TGraphAsymmErrors* p9059_d6x1y1     = NULL;
    if (returnStat) {
        p9059_d6x1y1                    = new TGraphAsymmErrors(p9059_d6x1y1_numpoints, p9059_d6x1y1_xval, p9059_d6x1y1_yval, p9059_d6x1y1_xerrminus, p9059_d6x1y1_xerrplus,p9059_d6x1y1_ystatminus, p9059_d6x1y1_ystatplus);
    } else {
        p9059_d6x1y1                    = new TGraphAsymmErrors(p9059_d6x1y1_numpoints, p9059_d6x1y1_xval, p9059_d6x1y1_yval, p9059_d6x1y1_xerrminus, p9059_d6x1y1_xerrplus,p9059_d6x1y1_ysystminus, p9059_d6x1y1_ysystplus);
    }
    return p9059_d6x1y1;
}


// @article{Adam:2015jda,
//       author         = "Adam, Jaroslav and others",
//       title          = "{Measurement of D$_{s}^{+}$ production and nuclear
//                         modification factor in Pb-Pb collisions at $
//                         \sqrt{{\mathrm{s}}_{\mathrm{NN}}}=$ 2.76 TeV}",
//       collaboration  = "ALICE",
//       journal        = "JHEP",
//       volume         = "03",
//       year           = "2016",
//       pages          = "082",
//       doi            = "10.1007/JHEP03(2016)082",
//       eprint         = "1509.07287",
//       archivePrefix  = "arXiv",
//       primaryClass   = "nucl-ex",
//       reportNumber   = "CERN-PH-EP-2015-253",
//       SLACcitation   = "%%CITATION = ARXIV:1509.07287;%%"
// }

//================================================================================================================
// D_s^+ meson in 0-10% PbPb collisions dN/dy
//================================================================================================================
TGraphAsymmErrors* GetDSPlusMesonPbPb0010(Bool_t returnStat){
    
    // Plot: p9019_d1x1y1
    double p9019_d1x1y1_xval[]          = { 5.0, 7.0, 10.0 };
    double p9019_d1x1y1_xerrminus[]     = { 1.0, 1.0, 2.0 };
    double p9019_d1x1y1_xerrplus[]      = { 1.0, 1.0, 2.0 };
    double p9019_d1x1y1_yval[]          = { 0.0373, 0.00485, 7.34E-4 };
    double p9019_d1x1y1_yerrminus[]     = { 0.022061051652176512, 0.0030749634144165033, 4.0388736053508775E-4 };
    double p9019_d1x1y1_yerrplus[]      = { 0.018162334651690572, 0.0025659844114881134, 3.140270688969344E-4 };
    double p9019_d1x1y1_ystatminus[]    = { 0.0142, 0.00197, 2.1E-4 };
    double p9019_d1x1y1_ystatplus[]     = { 0.0142, 0.00197, 2.1E-4 };
    int p9019_d1x1y1_numpoints          = 3;
    double p9019_d1x1y1_ysystminus[p9019_d1x1y1_numpoints];
    double p9019_d1x1y1_ysystplus[p9019_d1x1y1_numpoints];
    double branching                    = 1; //0.0224;
    
    for (Int_t i = 0; i < p9019_d1x1y1_numpoints; i++){
        p9019_d1x1y1_ysystminus[i]      = TMath::Sqrt(p9019_d1x1y1_yerrminus[i]*p9019_d1x1y1_yerrminus[i]-p9019_d1x1y1_ystatminus[i]*p9019_d1x1y1_ystatminus[i]);
        p9019_d1x1y1_ysystplus[i]       = TMath::Sqrt(p9019_d1x1y1_yerrplus[i]*p9019_d1x1y1_yerrplus[i]-p9019_d1x1y1_ystatplus[i]*p9019_d1x1y1_ystatplus[i]);
        
        p9019_d1x1y1_yval[i]            = p9019_d1x1y1_yval[i]/(branching);
        p9019_d1x1y1_yerrminus[i]       = p9019_d1x1y1_yerrminus[i]/(branching);
        p9019_d1x1y1_yerrplus[i]        = p9019_d1x1y1_yerrplus[i]/(branching);
        p9019_d1x1y1_ystatminus[i]      = p9019_d1x1y1_ystatminus[i]/(branching);
        p9019_d1x1y1_ystatplus[i]       = p9019_d1x1y1_ystatplus[i]/(branching);
        p9019_d1x1y1_ysystminus[i]      = p9019_d1x1y1_ysystminus[i]/(branching);
        p9019_d1x1y1_ysystplus[i]       = p9019_d1x1y1_ysystplus[i]/(branching);
    }
    
    TGraphAsymmErrors* p9019_d1x1y1     = NULL;
    if (returnStat) {
        p9019_d1x1y1                    = new TGraphAsymmErrors(p9019_d1x1y1_numpoints, p9019_d1x1y1_xval, p9019_d1x1y1_yval, p9019_d1x1y1_xerrminus, p9019_d1x1y1_xerrplus,p9019_d1x1y1_ystatminus, p9019_d1x1y1_ystatplus);
    } else {
        p9019_d1x1y1                    = new TGraphAsymmErrors(p9019_d1x1y1_numpoints, p9019_d1x1y1_xval, p9019_d1x1y1_yval, p9019_d1x1y1_xerrminus, p9019_d1x1y1_xerrplus,p9019_d1x1y1_ysystminus, p9019_d1x1y1_ysystplus);
    }
    return p9019_d1x1y1;
}

//================================================================================================================
// D_s^+ meson in 20-50% PbPb collisions dN/dy
//================================================================================================================
TGraphAsymmErrors* GetDSPlusMesonPbPb2050(Bool_t returnStat){
    
    // Plot: p9019_d2x1y1
    double p9019_d2x1y1_xval[]          = { 7.0, 10.0 };
    double p9019_d2x1y1_xerrminus[]     = { 1.0, 2.0 };
    double p9019_d2x1y1_xerrplus[]      = { 1.0, 2.0 };
    double p9019_d2x1y1_yval[]          = { 0.00324, 2.54E-4 };
    double p9019_d2x1y1_yerrminus[]     = { 0.001762971922635185, 1.6384178343755904E-4 };
    double p9019_d2x1y1_yerrplus[]      = { 0.0015299349659380953, 1.4224812828294087E-4 };
    double p9019_d2x1y1_ystatminus[]    = { 0.00117, 1.18E-4 };
    double p9019_d2x1y1_ystatplus[]     = { 0.00117, 1.18E-4 };
    int p9019_d2x1y1_numpoints          = 2;
    double p9019_d2x1y1_ysystminus[p9019_d2x1y1_numpoints];
    double p9019_d2x1y1_ysystplus[p9019_d2x1y1_numpoints];
    double branching                    = 1; //0.0224;
    
    for (Int_t i = 0; i < p9019_d2x1y1_numpoints; i++){
        p9019_d2x1y1_ysystminus[i]      = TMath::Sqrt(p9019_d2x1y1_yerrminus[i]*p9019_d2x1y1_yerrminus[i]-p9019_d2x1y1_ystatminus[i]*p9019_d2x1y1_ystatminus[i]);
        p9019_d2x1y1_ysystplus[i]       = TMath::Sqrt(p9019_d2x1y1_yerrplus[i]*p9019_d2x1y1_yerrplus[i]-p9019_d2x1y1_ystatplus[i]*p9019_d2x1y1_ystatplus[i]);
        
        p9019_d2x1y1_yval[i]            = p9019_d2x1y1_yval[i]/(branching);
        p9019_d2x1y1_yerrminus[i]       = p9019_d2x1y1_yerrminus[i]/(branching);
        p9019_d2x1y1_yerrplus[i]        = p9019_d2x1y1_yerrplus[i]/(branching);
        p9019_d2x1y1_ystatminus[i]      = p9019_d2x1y1_ystatminus[i]/(branching);
        p9019_d2x1y1_ystatplus[i]       = p9019_d2x1y1_ystatplus[i]/(branching);
        p9019_d2x1y1_ysystminus[i]      = p9019_d2x1y1_ysystminus[i]/(branching);
        p9019_d2x1y1_ysystplus[i]       = p9019_d2x1y1_ysystplus[i]/(branching);
    }
    
    TGraphAsymmErrors* p9019_d2x1y1     = NULL;
    if (returnStat) {
        p9019_d2x1y1                    = new TGraphAsymmErrors(p9019_d2x1y1_numpoints, p9019_d2x1y1_xval, p9019_d2x1y1_yval, p9019_d2x1y1_xerrminus, p9019_d2x1y1_xerrplus,p9019_d2x1y1_ystatminus, p9019_d2x1y1_ystatplus);
    } else {
        p9019_d2x1y1                    = new TGraphAsymmErrors(p9019_d2x1y1_numpoints, p9019_d2x1y1_xval, p9019_d2x1y1_yval, p9019_d2x1y1_xerrminus, p9019_d2x1y1_xerrplus,p9019_d2x1y1_ysystminus, p9019_d2x1y1_ysystplus);
    }
    return p9019_d2x1y1;
}
