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
#include "CocktailFitting.h"

//================================================================================================================
//Function to create canvas
//================================================================================================================
TCanvas *GetAndSetCanvas( TString name,
                         Double_t leftmargin = 0.11,
                         Double_t bottommargin = 0.1,
                         Double_t x = 1400,
                         Double_t y = 1000){

    TCanvas *canvas =  new TCanvas(name,name,x,y);
    canvas->SetLeftMargin(leftmargin);
    canvas->SetRightMargin(0.015);
    canvas->SetTopMargin(0.03);
    canvas->SetBottomMargin(bottommargin);
    canvas->SetFillColor(0);

    return canvas;
}

//================================================================================================================
//Function to set canvas settings
//================================================================================================================
void DrawCanvasSettings( TCanvas* c1,
                        Double_t leftMargin,
                        Double_t rightMargin,
                        Double_t topMargin,
                        Double_t bottomMargin){

    c1->SetTickx();
    c1->SetTicky();
    c1->SetGridx(0);
    c1->SetGridy(0);
    c1->SetLogy(0);
    c1->SetLeftMargin(leftMargin);
    c1->SetRightMargin(rightMargin);
    c1->SetTopMargin(topMargin);
    c1->SetBottomMargin(bottomMargin);
    c1->SetFillColor(0);
}

//================================================================================================================
//Function to set pad settings
//================================================================================================================
void DrawPadSettings( TPad* pad1,
                     Double_t leftMargin,
                     Double_t rightMargin,
                     Double_t topMargin,
                     Double_t bottomMargin){

    pad1->SetFillColor(0);
    pad1->GetFrame()->SetFillColor(0);
    pad1->SetBorderMode(0);
    pad1->SetLeftMargin(leftMargin);
    pad1->SetBottomMargin(bottomMargin);
    pad1->SetRightMargin(rightMargin);
    pad1->SetTopMargin(topMargin);
    pad1->SetTickx();
    pad1->SetTicky();
}

//================================================================================================================
//Function to create legend
//================================================================================================================
TLegend *GetAndSetLegend( Double_t positionX,
                          Double_t positionY,
                          Double_t entries,
                          Int_t nColumns                = 1,
                          TString header                =""
                        ){

    if(header.CompareTo("") != 0) entries++;
    Double_t positionYPlus = 0.04*1.1*(Double_t)entries;
    TLegend *legend = new TLegend(positionX,positionY,positionX+(0.25*nColumns),positionY+positionYPlus);
    legend->SetNColumns(nColumns);
    legend->SetLineColor(0);
    legend->SetLineWidth(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetLineStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    if(header.CompareTo("") != 0)legend->SetHeader(header);

    return legend;
}

//================================================================================================================
//Function to create legend
//================================================================================================================
TLegend *GetAndSetLegend( Double_t positionX,
                          Double_t positionY,
                          Double_t positionXRight,
                          Double_t positionYUp,
                          TString header              = "",
                          Int_t columns               = 2,
                          Double_t margin             = 0,
                          Double_t textSize           = 0.04,
                          Style_t textFont            = 42
                         ){

    TLegend *legend = new TLegend(positionX,positionY,positionXRight,positionYUp);
    legend->SetNColumns(columns);
    legend->SetLineColor(0);
    legend->SetLineWidth(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetLineStyle(0);
    legend->SetTextFont(textFont);
    legend->SetTextSize(textSize);
    legend->SetBorderSize(0);
    if (margin != 0) legend->SetMargin(margin);
    if (header.CompareTo("")!= 0) legend->SetHeader(header);
    return legend;
}

//================================================================================================================
//Function to set histogram settings
//================================================================================================================
void SetHistogramm( TH1 *hist,
                    TString xLabel,
                    TString yLabel,
                    Double_t rangeYlow   = -99.,
                    Double_t rangeYhigh  = -99.,
                    Double_t xOffset     = 1.0,
                    Double_t yOffset     = 1.15,
                    Double_t titleSize   = 0.04,
                    Bool_t bScale        = kTRUE,
                    Bool_t bTitleBold    = kFALSE,
                    Style_t font        = 42,
                    Int_t xNDivisions    = 510,
                    Int_t yNDivisions    = 510
                  ){

    Double_t scale = 1./gPad->GetAbsHNDC();
    //hist->GetXaxis()->SetRangeUser(rangeX[0],rangeX[1]);
    if(rangeYlow != -99.) hist->GetYaxis()->SetRangeUser(rangeYlow,rangeYhigh);
    hist->SetStats(kFALSE);
    hist->SetTitle("");
    hist->SetXTitle(xLabel);
    hist->SetYTitle(yLabel);
    hist->GetYaxis()->SetDecimals();
    hist->GetYaxis()->SetTitleOffset(yOffset/scale);
    hist->GetXaxis()->SetTitleOffset(xOffset);
    if (bScale){
        hist->GetXaxis()->SetTitleSize(0.04*scale);
        hist->GetYaxis()->SetTitleSize(0.04*scale);
        hist->GetXaxis()->SetLabelSize(0.035*scale);
        hist->GetYaxis()->SetLabelSize(0.035*scale);
    } else {
        hist->GetXaxis()->SetTitleSize(titleSize);
        hist->GetYaxis()->SetTitleSize(titleSize);
        hist->GetXaxis()->SetLabelSize(0.85*titleSize);
        hist->GetYaxis()->SetLabelSize(0.85*titleSize);
    }
    hist->GetXaxis()->SetLabelFont(font);
    hist->GetYaxis()->SetLabelFont(font);
    hist->GetXaxis()->SetTitleFont(font);
    hist->GetYaxis()->SetTitleFont(font);

    if (bTitleBold){
        hist->GetXaxis()->SetTitleFont(font+20);
        hist->GetYaxis()->SetTitleFont(font+20);
    }
    hist->SetMarkerSize(1.);
    hist->SetMarkerStyle(20);

    hist->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);
    hist->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}


//================================================================================================================
//Function to set graph settings
//================================================================================================================
void SetGraph( TGraph *graph,
              TString xLabel,
              TString yLabel,
              Double_t rangeYlow = -99.,
              Double_t rangeYhigh = -99.,
              Double_t xOffset = 1.0,
              Double_t yOffset = 1.15){

    Double_t scale = 1./gPad->GetAbsHNDC();
    if(rangeYlow != -99.) graph->GetYaxis()->SetRangeUser(rangeYlow,rangeYhigh);
    graph->GetXaxis()->SetTitle(xLabel);
    graph->GetYaxis()->SetTitle(yLabel);
    graph->GetYaxis()->SetDecimals();
    graph->GetYaxis()->SetTitleOffset(yOffset/scale);
    graph->GetXaxis()->SetTitleOffset(xOffset);
    graph->GetXaxis()->SetTitleSize(0.04*scale);
    graph->GetYaxis()->SetTitleSize(0.04*scale);
    graph->GetXaxis()->SetLabelSize(0.035*scale);
    graph->GetYaxis()->SetLabelSize(0.035*scale);
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetLabelFont(42);
    graph->SetMarkerSize(1.);
    graph->SetMarkerStyle(20);
}

//================================================================================================================
//Function to set histogram settings
//================================================================================================================
void SetStyleHistoTH1ForGraphs( TH1* histo,
                                TString XTitle,
                                TString YTitle,
                                Size_t xLableSize,
                                Size_t xTitleSize,
                                Size_t yLableSize,
                                Size_t yTitleSize,
                                Float_t xTitleOffset = 1,
                                Float_t yTitleOffset = 1,
                                Int_t xNDivisions    = 510,
                                Int_t yNDivisions    = 510,
                                Style_t textFont     = 42
                              ){

    histo->SetStats(kFALSE);

    histo->SetXTitle(XTitle);
    histo->SetYTitle(YTitle);
    histo->SetTitle("");

    histo->GetYaxis()->SetLabelFont(textFont);
    histo->GetXaxis()->SetLabelFont(textFont);
    histo->GetYaxis()->SetTitleFont(textFont+20);
    histo->GetXaxis()->SetTitleFont(textFont+20);

    histo->GetXaxis()->SetLabelSize(xLableSize);
    histo->GetXaxis()->SetTitleSize(xTitleSize);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(yLableSize);
    histo->GetYaxis()->SetTitleSize(yTitleSize);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}

//================================================================================================================
//Function to set histogram marker
//================================================================================================================
void DrawMarker(    TH1* histo1,
                    Style_t markerStyle,
                    Size_t markerSize,
                    Color_t markerColor,
                    Color_t lineColor,
                    Style_t fillStyle           = -1,
                    Color_t fillColor           = -10000
               ) {

    histo1->SetMarkerStyle(markerStyle);
    histo1->SetMarkerSize(markerSize);
    histo1->SetMarkerColor(markerColor);
    histo1->SetLineColor(lineColor);
    if (fillStyle != -1)
        histo1->SetFillStyle(fillStyle);
    if (fillStyle == 0)
        histo1->SetFillColor(0);
    if (fillStyle == 1001 && fillColor != -10000)
        histo1->SetFillColor(fillColor);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
}

void DrawMarker(    TGraph* graph,
                    Style_t markerStyle,
                    Size_t markerSize,
                    Color_t markerColor,
                    Color_t lineColor,
                    Style_t fillStyle           = -1,
                    Color_t fillColor           = -10000
               ) {

    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    if (fillStyle != -1)
        graph->SetFillStyle(fillStyle);
    graph->SetLineColor(lineColor);
    if (fillStyle == 0)
        graph->SetFillColor(0);
    if (fillStyle == 1001 && fillColor != -10000)
        graph->SetFillColor(fillColor);
    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetTitleFont(62);
    graph->GetXaxis()->SetTitleFont(62);
}

void DrawMarker(    TGraphErrors* graph,
                    Style_t markerStyle,
                    Size_t markerSize,
                    Color_t markerColor,
                    Color_t lineColor,
                    Style_t fillStyle           = -1,
                    Color_t fillColor           = -10000
               ) {

    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    if (fillStyle != -1)
        graph->SetFillStyle(fillStyle);
    graph->SetLineColor(lineColor);
    if (fillStyle == 0)
        graph->SetFillColor(0);
    if (fillStyle == 1001 && fillColor != -10000)
        graph->SetFillColor(fillColor);
    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetTitleFont(62);
    graph->GetXaxis()->SetTitleFont(62);
}

void DrawMarker(    TGraphAsymmErrors* graph,
                    Style_t markerStyle,
                    Size_t markerSize,
                    Color_t markerColor,
                    Color_t lineColor,
                    Style_t fillStyle           = -1,
                    Color_t fillColor           = -10000
               ) {

    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    if (fillStyle != -1)
        graph->SetFillStyle(fillStyle);
    graph->SetLineColor(lineColor);
    if (fillStyle == 0)
        graph->SetFillColor(0);
    if (fillStyle == 1001 && fillColor != -10000)
        graph->SetFillColor(fillColor);

    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetTitleFont(62);
    graph->GetXaxis()->SetTitleFont(62);
}

//================================================================================================================
//Function to draw a TF1
//================================================================================================================
void DrawFit(   TF1* fit1,
                Style_t lineStyle,
                Size_t lineWidth,
                Color_t lineColor ) {
    fit1->SetLineColor(lineColor);
    fit1->SetLineStyle(lineStyle);
    fit1->SetLineWidth(lineWidth);
}


//================================================================================================================
//Function to draw a line
//================================================================================================================
void DrawLine(  Float_t startX,
              Float_t endX,
              Float_t startY,
              Float_t endY,
              Float_t linew,
              Float_t lineColor = 4,
              Style_t lineStyle = 1){

    TLine * l1 = new TLine (startX,startY,endX,endY);
    l1->SetLineColor(lineColor);
    l1->SetLineWidth(linew);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
}

//================================================================================================================
//Function to set style setting for TLatex
//================================================================================================================
void SetStyleTLatex( TLatex* text,
                     Size_t textSize,
                     Width_t lineWidth,
                     Color_t textColor = 1,
                     Font_t textFont = 42,
                     Bool_t kNDC = kTRUE,
                     Short_t align = 11
                   ){
    if (kNDC) {text->SetNDC();}
    text->SetTextFont(textFont);
    text->SetTextColor(textColor);
    text->SetTextSize(textSize);
    text->SetLineWidth(lineWidth);
    text->SetTextAlign(align);
}


//================================================================================================================
//Function to obtain best choice of xRange for all histograms in list
//================================================================================================================
Double_t GetXRangeExtremaFromList(TList* list, TString type, Bool_t returnMax) {

    // get number of histograms in list
    Int_t numberOfObjects                       = list->GetEntries();

    // search list for particle spectra and fill maxima in array
    TString tempClass                           = "";
    TString tempName                            = "";
    TObject** tempObject                        = new TObject*[numberOfObjects];
    for (Int_t i=0; i<numberOfObjects; i++)
        tempObject[i]                           = new TObject;

    Double_t* tempMax                           = new Double_t[numberOfObjects];
    Double_t* tempMin                           = new Double_t[numberOfObjects];
    for (Int_t i=0; i<numberOfObjects; i++) {
        tempMax[i]                              = 0;
        tempMin[i]                              = 10000;
    }
    Double_t* tempMaxSys                        = new Double_t[numberOfObjects];
    Double_t* tempMinSys                        = new Double_t[numberOfObjects];
    for (Int_t i=0; i<numberOfObjects; i++) {
        tempMaxSys[i]                           = 0;
        tempMinSys[i]                           = 10000;
    }

    if (type.CompareTo("spectra") == 0) {
        for (Int_t i=0; i<numberOfObjects; i++) {
            tempName                            = ((TObject*)list->At(i))->GetName();
            tempClass                           = ((TObject*)list->At(i))->ClassName();

            for (Int_t particle=0; particle<nParticles; particle++) {
                for (Int_t method=0; method<9; method++) {

                    // check if object at position i is particle spectrum (stat. err.)
                    if (tempName.CompareTo(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())) == 0) {

                        // check for type of object and fill extrema in arrays
                        if (tempClass.Contains("TH1")) {

                            tempObject[i]       = (TH1D*)((TH1D*)list->FindObject(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMax[i]          = ((TH1D*)tempObject[i])->GetXaxis()->GetXmax();
                            if (((TH1D*)tempObject[i])->GetBinContent(1) > 0.){
                                tempMin[i]      = ((TH1D*)tempObject[i])->GetXaxis()->GetXmin();
                            } else {
                                Int_t bin       = 1;
                                while ( !(((TH1D*)tempObject[i])->GetBinContent(bin) > 1e-20))
                                    bin++;
                                tempMin[i]      = ((TH1D*)tempObject[i])->GetBinCenter(bin) - ((TH1D*)tempObject[i])->GetBinWidth(bin)/2;
                            }
                        } else if (tempClass.CompareTo("TGraph") == 0) {

                            tempObject[i]       = (TGraph*)((TGraph*)list->FindObject(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMax[i]          = GetXRangeFromGraph((TGraph*)tempObject[i], kTRUE);
                            tempMin[i]          = GetXRangeFromGraph((TGraph*)tempObject[i], kFALSE);

                        } else if (tempClass.CompareTo("TGraphErrors") == 0) {

                            tempObject[i]       = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMax[i]          = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kTRUE);
                            tempMin[i]          = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kFALSE);

                        } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {

                            tempObject[i]       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMax[i]          = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kTRUE);
                            tempMin[i]          = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kFALSE);

                        }
                    }
                }
            }
            for (Int_t particle=0; particle<nParticles; particle++) {
                for (Int_t method=0; method<9; method++) {

                    // check if object at position i is particle spectrum (stat. err.)
                    if (tempName.CompareTo(Form("%s%sSys", fParticle[particle].Data(), fMethod[method].Data())) == 0) {
                        // check for type of object and fill extrema in arrays
                        if (tempClass.Contains("TH1")) {

                            tempObject[i]       = (TH1D*)((TH1D*)list->FindObject(Form("%s%sSys", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMaxSys[i]       = ((TH1D*)tempObject[i])->GetXaxis()->GetXmax();
                            if (((TH1D*)tempObject[i])->GetBinContent(1) > 0.){
                                tempMaxSys[i]   = ((TH1D*)tempObject[i])->GetXaxis()->GetXmin();
                            } else {
                                Int_t bin       = 1;
                                while ( !(((TH1D*)tempObject[i])->GetBinContent(bin) > 1e-20))
                                    bin++;
                                tempMaxSys[i]   = ((TH1D*)tempObject[i])->GetBinCenter(bin) - ((TH1D*)tempObject[i])->GetBinWidth(bin)/2;
                            }

                        } else if (tempClass.CompareTo("TGraph") == 0) {

                            tempObject[i]       = (TGraph*)((TGraph*)list->FindObject(Form("%s%sSys", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMaxSys[i]       = GetXRangeFromGraph((TGraph*)tempObject[i], kTRUE);
                            tempMinSys[i]       = GetXRangeFromGraph((TGraph*)tempObject[i], kFALSE);

                        } else if (tempClass.CompareTo("TGraphErrors") == 0) {

                            tempObject[i]       = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%s%sSys", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMaxSys[i]       = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kTRUE);
                            tempMinSys[i]       = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kFALSE);

                        } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {

                            tempObject[i]       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%s%sSys", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMaxSys[i]       = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kTRUE);
                            tempMinSys[i]       = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kFALSE);

                        }
                    }
                }
            }
        }
    } else if (type.CompareTo("ratio") == 0) {
        for (Int_t i=0; i<numberOfObjects; i++) {
            tempName                            = ((TObject*)list->At(i))->GetName();
            tempClass                           = ((TObject*)list->At(i))->ClassName();

            for (Int_t particle1=0; particle1<nParticles; particle1++) {
                for (Int_t particle2=0; particle2<nParticles; particle2++) {
                    for (Int_t method=0; method<9; method++) {

                        // check if object at position i is particle spectrum (stat. err.)
                        if (tempName.CompareTo(Form("%sTo%s%sStat", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())) == 0) {

                            // check for type of object and fill extrema in arrays
                            if (tempClass.Contains("TH1")) {

                                tempObject[i]       = (TH1D*)((TH1D*)list->FindObject(Form("%sTo%s%sStat", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())))->Clone("tempObject");
                                if (((TH1D*)tempObject[i])->GetBinContent(((TH1D*)tempObject[i])->GetNbinsX()-1) > 0.){
                                    tempMax[i]          = ((TH1D*)tempObject[i])->GetXaxis()->GetXmax();
                                } else {
                                    Int_t bin       = ((TH1D*)tempObject[i])->GetNbinsX()-1;
                                    while ( !(((TH1D*)tempObject[i])->GetBinContent(bin) > 1e-20))
                                        bin--;
                                    tempMax[i]      = ((TH1D*)tempObject[i])->GetBinCenter(bin) + ((TH1D*)tempObject[i])->GetBinWidth(bin)/2;
                                }
                                if (((TH1D*)tempObject[i])->GetBinContent(1) > 0.){
                                    tempMin[i]      = ((TH1D*)tempObject[i])->GetXaxis()->GetXmin();
                                } else {
                                    Int_t bin       = 1;
                                    while ( !(((TH1D*)tempObject[i])->GetBinContent(bin) > 1e-20))
                                        bin++;
                                    tempMin[i]      = ((TH1D*)tempObject[i])->GetBinCenter(bin) - ((TH1D*)tempObject[i])->GetBinWidth(bin)/2;
                                }

                            } else if (tempClass.CompareTo("TGraph") == 0) {

                                tempObject[i]       = (TGraph*)((TGraph*)list->FindObject(Form("%sTo%s%sStat", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())))->Clone("tempObject");
                                tempMax[i]          = GetXRangeFromGraph((TGraph*)tempObject[i], kTRUE);
                                tempMin[i]          = GetXRangeFromGraph((TGraph*)tempObject[i], kFALSE);

                            } else if (tempClass.CompareTo("TGraphErrors") == 0) {

                                tempObject[i]       = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%sTo%s%sStat", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())))->Clone("tempObject");
                                tempMax[i]          = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kTRUE);
                                tempMin[i]          = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kFALSE);

                            } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {

                                tempObject[i]       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%sTo%s%sStat", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())))->Clone("tempObject");
                                tempMax[i]          = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kTRUE);
                                tempMin[i]          = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kFALSE);

                            }
                        }
                    }
                }
            }
            for (Int_t particle1=0; particle1<nParticles; particle1++) {
                for (Int_t particle2=0; particle2<nParticles; particle2++) {
                    for (Int_t method=0; method<9; method++) {

                        // check if object at position i is particle spectrum (stat. err.)
                        if (tempName.CompareTo(Form("%sTo%s%sSys", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())) == 0) {

                            // check for type of object and fill extrema in arrays
                            if (tempClass.Contains("TH1")) {

                                tempObject[i]       = (TH1D*)((TH1D*)list->FindObject(Form("%sTo%s%sSys", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())))->Clone("tempObject");
                                if (((TH1D*)tempObject[i])->GetBinContent(((TH1D*)tempObject[i])->GetNbinsX()-1) > 0.){
                                    tempMax[i]          = ((TH1D*)tempObject[i])->GetXaxis()->GetXmax();
                                } else {
                                    Int_t bin       = ((TH1D*)tempObject[i])->GetNbinsX()-1;
                                    while ( !(((TH1D*)tempObject[i])->GetBinContent(bin) > 1e-20) && bin > 1)
                                        bin--;
                                    tempMax[i]      = ((TH1D*)tempObject[i])->GetBinCenter(bin) + ((TH1D*)tempObject[i])->GetBinWidth(bin)/2;
                                }
                                if (((TH1D*)tempObject[i])->GetBinContent(1) > 0.){
                                    tempMinSys[i]   = ((TH1D*)tempObject[i])->GetXaxis()->GetXmin();
                                } else {
                                    Int_t bin       = 1;
                                    while ( !(((TH1D*)tempObject[i])->GetBinContent(bin) > 1e-20))
                                        bin++;
                                    tempMinSys[i]   = ((TH1D*)tempObject[i])->GetBinCenter(bin) - ((TH1D*)tempObject[i])->GetBinWidth(bin)/2;
                                }

                            } else if (tempClass.CompareTo("TGraph") == 0) {

                                tempObject[i]       = (TGraph*)((TGraph*)list->FindObject(Form("%sTo%s%sSys", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())))->Clone("tempObject");
                                tempMaxSys[i]       = GetXRangeFromGraph((TGraph*)tempObject[i], kTRUE);
                                tempMinSys[i]       = GetXRangeFromGraph((TGraph*)tempObject[i], kFALSE);

                            } else if (tempClass.CompareTo("TGraphErrors") == 0) {

                                tempObject[i]       = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%sTo%s%sSys", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())))->Clone("tempObject");
                                tempMaxSys[i]       = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kTRUE);
                                tempMinSys[i]       = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kFALSE);

                            } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {

                                tempObject[i]       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%sTo%s%sSys", fParticle[particle1].Data(), fParticle[particle2].Data(), fMethod[method].Data())))->Clone("tempObject");
                                tempMaxSys[i]       = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kTRUE);
                                tempMinSys[i]       = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kFALSE);

                            }
                        }
                    }
                }
            }
        }
    } else {
        for (Int_t i=0; i<numberOfObjects; i++) {
            tempClass                           = ((TObject*)list->At(i))->ClassName();

            if (tempClass.Contains("TH1")) {
                tempObject[i]                   = (TH1D*)list->At(i);
                tempMax[i]                      = ((TH1D*)tempObject[i])->GetXaxis()->GetXmax();
                tempMin[i]                      = ((TH1D*)tempObject[i])->GetXaxis()->GetXmin();
            } else if (tempClass.CompareTo("TGraph") == 0) {
                tempObject[i]                   = (TGraph*)list->At(i);
                tempMax[i]                      = GetXRangeFromGraph((TGraph*)tempObject[i], kTRUE);
                tempMin[i]                      = GetXRangeFromGraph((TGraph*)tempObject[i], kFALSE);
            } else if (tempClass.CompareTo("TGraphErrors") == 0) {
                tempObject[i]                   = (TGraphErrors*)list->At(i);
                tempMax[i]                      = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kTRUE);
                tempMin[i]                      = GetXRangeFromGraph((TGraphErrors*)tempObject[i], kFALSE);
            } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                tempObject[i]                   = (TGraphAsymmErrors*)list->At(i);
                tempMax[i]                      = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kTRUE);
                tempMin[i]                      = GetXRangeFromGraph((TGraphAsymmErrors*)tempObject[i], kFALSE);
            }
        }
    }

    // sort arrays using sort() - lowest to highest
    sort(tempMax, tempMax + numberOfObjects);
    sort(tempMin, tempMin + numberOfObjects);
    sort(tempMaxSys, tempMaxSys + numberOfObjects);
    sort(tempMinSys, tempMinSys + numberOfObjects);

    // set total min and max
    Double_t totalMax                           = tempMax[numberOfObjects-1];
    Double_t totalMin                           = tempMin[0];
    Double_t totalMaxSys                        = tempMaxSys[numberOfObjects-1];
    Double_t totalMinSys                        = tempMinSys[0];
    if (totalMin > totalMinSys)
        totalMin                                = totalMinSys;
    if (totalMax < totalMaxSys)
        totalMax                                = totalMaxSys;
    // free pointer
    delete[] tempObject;
    delete[] tempMax;
    delete[] tempMin;
    delete[] tempMaxSys;
    delete[] tempMinSys;

    // return either total maximum or total minimum
    if (returnMax) {
        return totalMax;
    } else {
        return totalMin;
    }
}

//================================================================================================================
//Function to obtain best choice of yRange for all histograms in list
//================================================================================================================
Double_t GetYRangeExtremaFromList(TList* list, Bool_t doParticleSpectra, Bool_t returnMax, TString collSystem = "") {

    // get number of histograms in list
    Int_t numberOfObjects                       = list->GetEntries();

    // search list for particle spectra and fill maxima in array
    TString tempClass                           = "";
    TString tempName                            = "";
    TObject** tempObject                        = new TObject*[numberOfObjects];
    for (Int_t i=0; i<numberOfObjects; i++)
        tempObject[i]                           = new TObject;

    Double_t* tempMax                           = new Double_t[numberOfObjects];
    Double_t* tempMin                           = new Double_t[numberOfObjects];
    for (Int_t i=0; i<numberOfObjects; i++) {
        tempMax[i]                              = 0;
        tempMin[i]                              = 0;
    }

    if (doParticleSpectra) {
        for (Int_t i=0; i<numberOfObjects; i++) {
            tempName                            = ((TObject*)list->At(i))->GetName();
            tempClass                           = ((TObject*)list->At(i))->ClassName();

            for (Int_t particle=0; particle<nParticles; particle++) {
                for (Int_t method=0; method<9; method++) {

                    // check if object at position i is particle spectrum (stat. err.)
                    if (tempName.CompareTo(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())) == 0) {

                        // check for type of object and fill extrema in arrays
                        if (tempClass.Contains("TH1")) {
                            tempObject[i]       = (TH1D*)((TH1D*)list->FindObject(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMax[i]          = 0;
                            tempMin[i]          = 1e11;
                            for (Int_t l = 0; l < ((TH1D*)tempObject[i])->GetNbinsX()+1; l++){
                                if (tempMax[i] < ((TH1D*)tempObject[i])->GetBinContent(l))
                                    tempMax[i]  = ((TH1D*)tempObject[i])->GetBinContent(l);
                                if (tempMin[i] > ((TH1D*)tempObject[i])->GetBinContent(l) && ((TH1D*)tempObject[i])->GetBinContent(l) > 0)
                                    tempMin[i]  = ((TH1D*)tempObject[i])->GetBinContent(l);
                            }
//                             tempMax[i]          = ((TH1D*)tempObject[i])->GetMaximum();
//                             tempMin[i]          = ((TH1D*)tempObject[i])->GetMinimum(0);
                        } else if (tempClass.CompareTo("TGraph") == 0) {
                            tempObject[i]       = (TGraph*)((TGraph*)list->FindObject(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMax[i]          = TMath::MaxElement(((TGraph*)tempObject[i])->GetN(),((TGraph*)tempObject[i])->GetY());
                            tempMin[i]          = TMath::MinElement(((TGraph*)tempObject[i])->GetN(),((TGraph*)tempObject[i])->GetY());
                        } else if (tempClass.CompareTo("TGraphErrors") == 0) {
                            tempObject[i]       = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMax[i]          = TMath::MaxElement(((TGraphErrors*)tempObject[i])->GetN(),((TGraphErrors*)tempObject[i])->GetY());
                            tempMin[i]          = TMath::MinElement(((TGraphErrors*)tempObject[i])->GetN(),((TGraphErrors*)tempObject[i])->GetY());
                        } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                            tempObject[i]       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%s%sStat", fParticle[particle].Data(), fMethod[method].Data())))->Clone("tempObject");
                            tempMax[i]          = TMath::MaxElement(((TGraphAsymmErrors*)tempObject[i])->GetN(),((TGraphAsymmErrors*)tempObject[i])->GetY());
                            tempMin[i]          = TMath::MinElement(((TGraphAsymmErrors*)tempObject[i])->GetN(),((TGraphAsymmErrors*)tempObject[i])->GetY());
                        }
                    }
                }
            }
        }
    } else {
        for (Int_t i=0; i<numberOfObjects; i++) {
            tempClass                           = ((TObject*)list->At(i))->ClassName();

            if (tempClass.Contains("TH1")) {
                tempObject[i]                   = (TH1D*)list->At(i);
                tempMax[i]                      = 0;
                tempMin[i]                      = 1e11;
                for (Int_t l = 0; l < ((TH1D*)tempObject[i])->GetNbinsX()+1; l++){
                    if (tempMax[i] < ((TH1D*)tempObject[i])->GetBinContent(l))
                        tempMax[i]  = ((TH1D*)tempObject[i])->GetBinContent(l);
                    if (tempMin[i] > ((TH1D*)tempObject[i])->GetBinContent(l) && ((TH1D*)tempObject[i])->GetBinContent(l) > 0)
                        tempMin[i]  = ((TH1D*)tempObject[i])->GetBinContent(l);
                }
            } else if (tempClass.CompareTo("TGraph") == 0) {
                tempObject[i]                   = (TGraph*)list->At(i);
                tempMax[i]                      = TMath::MaxElement(((TGraph*)tempObject[i])->GetN(),((TGraph*)tempObject[i])->GetY());
                tempMin[i]                      = TMath::MinElement(((TGraph*)tempObject[i])->GetN(),((TGraph*)tempObject[i])->GetY());
            } else if (tempClass.CompareTo("TGraphErrors") == 0) {
                tempObject[i]                   = (TGraphErrors*)list->At(i);
                tempMax[i]                      = TMath::MaxElement(((TGraphErrors*)tempObject[i])->GetN(),((TGraphErrors*)tempObject[i])->GetY());
                tempMin[i]                      = TMath::MinElement(((TGraphErrors*)tempObject[i])->GetN(),((TGraphErrors*)tempObject[i])->GetY());
            } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                tempObject[i]                   = (TGraphAsymmErrors*)list->At(i);
                tempMax[i]                      = TMath::MaxElement(((TGraphAsymmErrors*)tempObject[i])->GetN(),((TGraphAsymmErrors*)tempObject[i])->GetY());
                tempMin[i]                      = TMath::MinElement(((TGraphAsymmErrors*)tempObject[i])->GetN(),((TGraphAsymmErrors*)tempObject[i])->GetY());
            }
        }
    }

    // sort arrays using sort() - lowest to highest
    sort(tempMax, tempMax + numberOfObjects);
    sort(tempMin, tempMin + numberOfObjects);

    // set total min and max
    Double_t totalMax                           = tempMax[numberOfObjects-1];
    Double_t totalMin                           = tempMin[0];
    if (totalMin == 0) {
        for (Int_t i=1; i<numberOfObjects; i++) {
            if (tempMin[i] != 0 && tempMin[i]<totalMin && tempMin[i]>1e-60) {
                totalMin                        = tempMin[i];
                break;
            } else {
                continue;
            }
        }
    }
    if (totalMin == 0 && !collSystem.Contains("Pb") )
        totalMin                 = 1e-9;
    else
        totalMin                 = 1e-7;

    // free pointer
    delete[] tempObject;
    delete[] tempMax;
    delete[] tempMin;

    // return either total maximum or total minimum
    if (returnMax) {
        return totalMax;
    } else {
        return totalMin;
    }
}

//================================================================================================================
//Function to produce plot containing all particle spectra in list
//================================================================================================================
void ProduceParticleSpectraPlotFromList(TList* list, TString collSys, TString energy, TString cent, TString suffix) {

    // create output dir
    TString outputDir                       = Form("plots/%s/%s",suffix.Data(),list->GetName());
    gSystem->Exec("mkdir -p "+outputDir);

    // get xRange
    Double_t xMin                           = 0;
//     Double_t xMin                           = GetXRangeExtremaFromList(list, "spectra", kFALSE)*0.2;
    Double_t xMax                           = GetXRangeExtremaFromList(list, "spectra", kTRUE);
//     Double_t xMax                           = GetXRangeExtremaFromList(list, "spectra", kTRUE)*1.5;
//     if (xMin < 0 || xMin == 0)
//       xMin                                  = 0.01;

    // get yRange
    Double_t yMin                           = GetYRangeExtremaFromList(list, kTRUE, kFALSE, collSys) * 0.5;
    Double_t yMax                           = GetYRangeExtremaFromList(list, kTRUE, kTRUE, collSys) * 2;

    // create canvas
//     TCanvas* canvas                         = GetAndSetCanvas("canvas");
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0.11, 0.1, 1400, 2000);

    DrawCanvasSettings(canvas, 0.152, 0.015, 0.015, 0.068);
    canvas->SetLogy();
//     canvas->SetLogx();

    // count number of stat histograms in list
    Int_t histCounter                       = 0;
    TString tempName                        = "";
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                            = ((TObject*)list->At(i))->GetName();
        if (tempName.Contains("Stat") && !tempName.Contains("To"))
            histCounter++;
        else
            continue;
    }

    // get latex for enery and centrality
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }

    // get number of rows for legend
    Int_t nRows                             = 0;
    if (histCounter%2 == 0) nRows           = histCounter/2 + 1;
    else nRows                              = (histCounter+1)/2 + 1;

    // create legend
    TLegend* legend                         = GetAndSetLegend(0.4, 0.97-(0.035*nRows), 0.95, 0.97, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()),2,0.12,0.035);

    // dummy histogram
    TH1D* dummyHisto                        = new TH1D("dummyHisto", "", 100000, xMin, xMax);
    SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", yMin, yMax, 0.7, 1.65, 0.04, kFALSE, kTRUE);
    dummyHisto->GetXaxis()->SetTickLength(0.02);
    dummyHisto->GetXaxis()->SetLabelOffset(-0.001);
    dummyHisto->Draw();

    // search list for particle spectra
    TH1D* tempHist                          = NULL;
    TGraph* tempGraph                       = NULL;
    TGraphErrors* tempGraphErrs             = NULL;
    TGraphAsymmErrors* tempGraphAsymErrs    = NULL;
    TString tempClass                       = "";

    Int_t j                                 = 0;
    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t ii=0; ii<9; ii++) {
            if (list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data()))) {
                tempClass                   = ((TObject*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->ClassName();

                if (tempClass.Contains("TH1")) {

                    tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempHist");
                    DrawMarker(tempHist, markers[j], markerSize[j], paint[j],paint[j]);
                    tempHist->Draw("e1,same");
                    legend->AddEntry(tempHist, Form("%s %s", fParticleLatex[i].Data(), fMethod[ii].Data()), "lp");

                    j++;
                } else if (tempClass.CompareTo("TGraph") == 0) {

                    tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraph");
                    DrawMarker(tempGraph, markers[j], markerSize[j], paint[j],paint[j]);
                    tempGraph->DrawClone("p");
                    legend->AddEntry(tempGraph, Form("%s %s", fParticleLatex[i].Data(), fMethod[ii].Data()), "lp");

                    j++;
                } else if (tempClass.CompareTo("TGraphErrors") == 0) {

                    tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraphErrs");
                    DrawMarker(tempGraphErrs, markers[j], markerSize[j], paint[j],paint[j]);
                    tempGraphErrs->DrawClone("p");

                    legend->AddEntry(tempGraphErrs, Form("%s %s", fParticleLatex[i].Data(), fMethod[ii].Data()), "lp");

                    j++;
                } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {

                    tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data()));
                    DrawMarker(tempGraphAsymErrs, markers[j], markerSize[j], paint[j],paint[j]);
                    tempGraphAsymErrs->DrawClone("p");

                    legend->AddEntry(tempGraphAsymErrs, Form("%s %s", fParticleLatex[i].Data(), fMethod[ii].Data()), "lp");

                    j++;
                    if (j>9) j=0;
                }
            }
        }
    }

    // draw legend
    legend->Draw("same");

    // save canvas
    if (histCounter) canvas->SaveAs(Form("%s/ParticleSpectra.%s", outputDir.Data(), suffix.Data()));

    // free pointer
    delete tempHist;
    delete legend;
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot containing all particle spectra in list
//================================================================================================================
void ProduceParticleSpectraPlotFromListOnlyFinal(   TList* list,
                                                    TString collSys,
                                                    TString energy,
                                                    TString cent,
                                                    TString suffix
                                                ) {

    // create output dir
    TString outputDir                       = Form("plots/%s/%s",suffix.Data(),list->GetName());
    gSystem->Exec("mkdir -p "+outputDir);

    // get xRange
    Double_t xMin                           = GetXRangeExtremaFromList(list, "spectra", kFALSE);
    Double_t xMax                           = GetXRangeExtremaFromList(list, "spectra", kTRUE);
    cout << "Plotting from: " << xMin << "\t" << xMax << endl;
    Bool_t xLogPlot                         = kTRUE;
    if (xMin < 0.01)
        xMin = 0.035;
    //         xLogPlot                            = kFALSE;

    // get yRange
    Double_t yMin                           = GetYRangeExtremaFromList(list, kTRUE, kFALSE, collSys) * 0.5;
    Double_t yMax                           = GetYRangeExtremaFromList(list, kTRUE, kTRUE, collSys) * 2;
    if (xLogPlot && collSys.Contains("Pb"))
        yMax                                = yMax*2;


    // create canvas
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 2000);
    DrawCanvasSettings(canvas, 0.152, 0.015, 0.015, 0.068);
    canvas->SetLogy();
    if (xLogPlot)canvas->SetLogx();

    // count number of stat histograms in list
    Int_t histCounter                       = 0;
    TString tempName                        = "";
    Bool_t hasHeavyParticle                 = kFALSE;
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                            = ((TObject*)list->At(i))->GetName();
        if (tempName.Contains("Stat") && !tempName.Contains("To")){
            if ( tempName.Contains("DPlus") || tempName.Contains("DZero") ||  tempName.Contains("DStarPlus") ||  tempName.Contains("DSPlus"))
                hasHeavyParticle            = kTRUE;

            if (tempName.Contains("NPion") || tempName.Contains("Eta")){
                if (tempName.Contains("Comb"))
                    histCounter++;
                else
                    continue;
            } else {
                histCounter++;
            }
        } else {
            continue;
        }
    }
    cout << "counted: " << histCounter << " inputs"<< endl;
    if (hasHeavyParticle ){
        cout << "resetting yMin" << endl;
        yMin                                = yMin*0.005;
    }

    // get latex for enery and centrality
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }

    // get number of rows for legend
    Int_t nRows                             = 0;
    if (histCounter%2 == 0) nRows           = histCounter/2;
    else nRows                              = (histCounter+1)/2;

    // create legend
    TLegend* legend                         = GetAndSetLegend(0.20, 0.09, 0.62, 0.09+(1.05*0.04*nRows), "" ,2,0.14,0.04);
    TLatex *labelEnergy                     = new TLatex(0.94, 0.925, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);

    TString textALICE                       = "ALICE";
    if (isThesis) textALICE                 = "ALICE this thesis";

    TLatex *labelALICE                     = new TLatex(0.94, 0.885, textALICE);
    SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);

//     Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data())
    // dummy histogram
    TH1D* dummyHisto                        = new TH1D("dummyHisto", "", 1000000, xMin*0.8, xMax*1.5);
    if (xLogPlot){
        SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", yMin, yMax, 0.7, 1.65, 0.04, kFALSE, kTRUE);
        dummyHisto->GetXaxis()->SetLabelOffset(-0.01);
    } else {
        SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", yMin, yMax, 0.7, 1.65, 0.04, kFALSE, kTRUE);
        dummyHisto->GetXaxis()->SetLabelOffset(-0.001);
        dummyHisto->GetXaxis()->SetRangeUser(0,xMax);
    }
    dummyHisto->GetXaxis()->SetTickLength(0.02);

    dummyHisto->Draw();

    // search list for particle spectra
    TH1D* tempHist                          = NULL;
    TGraph* tempGraph                       = NULL;
    TGraphErrors* tempGraphErrs             = NULL;
    TGraphAsymmErrors* tempGraphAsymErrs    = NULL;
    TString tempClass                       = "";

    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t ii=0; ii<2; ii++) {
            if (list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data()))) {
                tempClass                   = ((TObject*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->ClassName();
//                 cout << Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data()) << "\t"<< tempClass.Data() << endl;
                if (tempClass.Contains("TH1")) {
                    tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempHist");
                    TGraphErrors* tempGrC   = new TGraphErrors(tempHist);
                    DrawMarker(tempGrC, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]),0);
                    tempGrC->Draw("e2");
                } else if (tempClass.CompareTo("TGraph") == 0) {
                    tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraph");
                    DrawMarker(tempGraph, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]),0);
                    tempGraph->DrawClone("e2");
                } else if (tempClass.CompareTo("TGraphErrors") == 0) {

                    tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraphErrs");
                    DrawMarker(tempGraphErrs, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]),0);
                    tempGraphErrs->DrawClone("e2");
               } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                    tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->Clone(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data()));
                    DrawMarker(tempGraphAsymErrs, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]),0);
                    tempGraphAsymErrs->DrawClone("e2");
                }
            }
            if (list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data()))) {
                tempClass                   = ((TObject*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->ClassName();
//                 cout << Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data()) << "\t"<< tempClass.Data() << endl;
                if (tempClass.Contains("TH1")) {
                    tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempHist");
                    DrawMarker(tempHist, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]));
                    tempHist->Draw("e1x0,same");
                    legend->AddEntry(tempHist, Form("%s", fParticleLatex[i].Data()), "lp");
                } else if (tempClass.CompareTo("TGraph") == 0) {
                    tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraph");
                    DrawMarker(tempGraph, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]));
                    tempGraph->DrawClone("p");
                    legend->AddEntry(tempGraph, Form("%s", fParticleLatex[i].Data()), "lp");
                } else if (tempClass.CompareTo("TGraphErrors") == 0) {
                    tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraphErrs");
                    ProduceGraphWithoutXErrors(tempGraphErrs);
                    DrawMarker(tempGraphErrs, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]));
                    tempGraphErrs->DrawClone("p");

                    legend->AddEntry(tempGraphErrs, Form("%s", fParticleLatex[i].Data()), "lp");
                } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                    tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data()));
                    ProduceGraphAsymmWithoutXErrors(tempGraphAsymErrs);
                    DrawMarker(tempGraphAsymErrs, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]));
                    tempGraphAsymErrs->DrawClone("p");

                    legend->AddEntry(tempGraphAsymErrs, Form("%s", fParticleLatex[i].Data()), "lp");
                }
            }

        }
    }
    dummyHisto->Draw("same,axis");

    // draw legend
    legend->Draw("same");
    labelEnergy->Draw();
    labelALICE->Draw();

    // save canvas
    if (histCounter) canvas->SaveAs(Form("%s/ParticleSpectra_Final.%s", outputDir.Data(), suffix.Data()));

    // free pointer
    delete tempHist;
    delete legend;
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot containing all particle spectra in list
//================================================================================================================
void ProduceParticleSpectraPlotFromListOnlyFinalWithFits(   TList* list,
                                                            TString collSys,
                                                            TString energy,
                                                            TString cent,
                                                            TString suffix,
                                                            TString outputDir,
                                                            Bool_t  xLogPlot        = kTRUE,
                                                            Double_t xMaxPlot       = -1000,
                                                            TString defaultFit      = "oHagPt"
) {

    cout << "entered" << endl;
    // create output dir
    gSystem->Exec("mkdir -p "+outputDir);
    cout << collSys.Data() << "\t"<< energy.Data() << "\t" << outputDir.Data() << endl;

    // get xRange
    Double_t xMin                           = GetXRangeExtremaFromList(list, "spectra", kFALSE);
    Double_t xMax                           = GetXRangeExtremaFromList(list, "spectra", kTRUE);
    cout << "Plotting from: " << xMin << "\t" << xMax << endl;
    if (xMin < 0.01)
        xMin                                = 0.035;
    if (xMaxPlot != -1000)
        xMax                                = xMaxPlot;

    // get yRange
    Double_t yMin                           = GetYRangeExtremaFromList(list, kTRUE, kFALSE, collSys) * 0.5;
    Double_t yMax                           = GetYRangeExtremaFromList(list, kTRUE, kTRUE, collSys) * 2;
    if (xLogPlot && collSys.Contains("Pb"))
        yMax                                = yMax*2;

    if (xMaxPlot != -1000){
        Double_t yMinFit                    = 1e8;

        for (Int_t i=0; i<nParticles; i++) {
            if (list->FindObject(Form("%sStat_Fit_%s", fParticle[i].Data(), defaultFit.Data()))) {
                TF1* tempFit2               = (TF1*)list->FindObject(Form("%sStat_Fit_%s", fParticle[i].Data(), defaultFit.Data()));
                Double_t cYM                = tempFit2->Eval(xMaxPlot);
                if (cYM < yMinFit)
                    yMinFit                 = cYM;
            }
        }
        if (yMinFit != 1e8)
            yMin                            = yMinFit*0.5;
    }

    // create canvas
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 2000);
    DrawCanvasSettings(canvas, 0.152, 0.015, 0.015, 0.068);
    canvas->SetLogy();
    if (xLogPlot)canvas->SetLogx();

    // count number of stat histograms in list
    Int_t histCounter                       = 0;
    TString tempName                        = "";
    Bool_t hasHeavyParticle                 = kFALSE;
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                            = ((TObject*)list->At(i))->GetName();
//         cout << tempName.Data() << endl;
        if (tempName.Contains("Stat") ){
            if (!tempName.Contains("To") && !tempName.Contains("Fit")) {
                if ( tempName.Contains("COmega") || tempName.Contains("Xi") ||  tempName.Contains("Phi")  ||  tempName.Contains("Sigma")  )
                    hasHeavyParticle                = kTRUE;
//                 cout << "accepted: "<< tempName.Data() << endl;
                histCounter++;
            } else {
                continue;
            }
        } else {
            continue;
        }
    }
    cout << "counted: " << histCounter << " inputs"<< endl;

    // get latex for enery and centrality
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }

    // get number of rows for legend
    Int_t nRows                             = 0;
    if (histCounter%2 == 0) nRows           = histCounter/2;
    else nRows                              = (histCounter+1)/2;

    // create legend
    TLegend* legend                         = NULL;
    if (xLogPlot ){
        legend                              = GetAndSetLegend(0.2, 0.09, 0.62, 0.08+(1.05*0.04*nRows), "" ,2,0.14,0.04);
    } else if (xMax > 10 ){
        legend                              = GetAndSetLegend(0.53, 0.885-(1.05*0.04*nRows), 0.95, 0.885, "" ,2,0.14,0.04);
    } else {
        legend                              = GetAndSetLegend(0.2, 0.09, 0.62, 0.08+(1.05*0.04*nRows), "" ,2,0.14,0.04);
    }
    TLatex *labelEnergy                     = new TLatex(0.94, 0.925, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);
    TString textALICE                       = "ALICE";
    if (isThesis) textALICE                 = "ALICE this thesis";

    TLatex *labelALICE                     = new TLatex(0.94, 0.885, textALICE);
    SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);

    if (hasHeavyParticle && xMax < 5){
        cout << "resetting yMin" << endl;
        yMin                                = yMin*0.1;
    }
    // dummy histogram
    TH1D* dummyHisto                        = NULL;

    if (xLogPlot){
        dummyHisto                          = new TH1D("dummyHisto", "", 100000, xMin*0.8, xMax*1.5);
        SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", yMin, yMax, 0.7, 1.65, 0.04, kFALSE, kTRUE);
        dummyHisto->GetXaxis()->SetLabelOffset(-0.01);
    } else {
        dummyHisto                          = new TH1D("dummyHisto", "", 100000, 0, xMax*1.5);
        SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", yMin, yMax, 0.73, 1.65, 0.04, kFALSE, kTRUE);
        dummyHisto->GetXaxis()->SetLabelOffset(-0.001);
        dummyHisto->GetXaxis()->SetRangeUser(0,xMax);
    }
    dummyHisto->GetXaxis()->SetTickLength(0.02);

    dummyHisto->Draw();

    // search list for particle spectra
    TH1D* tempHist                          = NULL;
    TGraph* tempGraph                       = NULL;
    TGraphErrors* tempGraphErrs             = NULL;
    TGraphAsymmErrors* tempGraphAsymErrs    = NULL;
    TString tempClass                       = "";
    TF1* tempFit                            = NULL;

    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t ii=0; ii<1; ii++) {
            if (list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data()))) {
                tempClass                   = ((TObject*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->ClassName();
//                 cout << Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data()) << "\t"<< tempClass.Data() << endl;
//                 cout << GetParticleMarkerStyle(fParticle[i]) << "\t" << GetParticleMarkerSize(fParticle[i]) << "\t" << GetParticleColor(fParticle[i]) << endl;
                if (tempClass.Contains("TH1")) {
                    tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempHist");
                    TGraphErrors* tempGrC   = new TGraphErrors(tempHist);
                    DrawMarker(tempGrC, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]),0);
                    tempGrC->Draw("e2");
                } else if (tempClass.CompareTo("TGraph") == 0) {
                    tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraph");
                    DrawMarker(tempGraph, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]),0);
                    tempGraph->DrawClone("e2");
                } else if (tempClass.CompareTo("TGraphErrors") == 0) {

                    tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraphErrs");
                    DrawMarker(tempGraphErrs, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]),0);
                    tempGraphErrs->DrawClone("e2");
                } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                    tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data())))->Clone(Form("%s%sSys", fParticle[i].Data(), fMethod[ii].Data()));
                    DrawMarker(tempGraphAsymErrs, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]),0);
                    tempGraphAsymErrs->DrawClone("e2");
                }
            }
            if (list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data()))) {
                tempClass                   = ((TObject*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->ClassName();
//                 cout << Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data()) << "\t"<< tempClass.Data() << endl;
                if (tempClass.Contains("TH1")) {
                    tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempHist");
                    DrawMarker(tempHist, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]));
                    tempHist->Draw("e1x0,same");
                    legend->AddEntry(tempHist, Form("%s", fParticleLatex[i].Data()), "lp");
                } else if (tempClass.CompareTo("TGraph") == 0) {
                    tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraph");
                    DrawMarker(tempGraph, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]));
                    tempGraph->DrawClone("p");
                    legend->AddEntry(tempGraph, Form("%s", fParticleLatex[i].Data()), "lp");
                } else if (tempClass.CompareTo("TGraphErrors") == 0) {
                    tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone("tempGraphErrs");
                    ProduceGraphWithoutXErrors(tempGraphErrs);
                    DrawMarker(tempGraphErrs, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]));
                    tempGraphErrs->DrawClone("p");

                    legend->AddEntry(tempGraphErrs, Form("%s", fParticleLatex[i].Data()), "lp");
                } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                    tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data())))->Clone(Form("%s%sStat", fParticle[i].Data(), fMethod[ii].Data()));
                    ProduceGraphAsymmWithoutXErrors(tempGraphAsymErrs);
                    DrawMarker(tempGraphAsymErrs, GetParticleMarkerStyle(fParticle[i]), GetParticleMarkerSize(fParticle[i]), GetParticleColor(fParticle[i]),GetParticleColor(fParticle[i]));
                    tempGraphAsymErrs->DrawClone("p");

                    legend->AddEntry(tempGraphAsymErrs, Form("%s", fParticleLatex[i].Data()), "lp");
                }
            }
            if (list->FindObject(Form("%sStat_Fit_%s", fParticle[i].Data(), defaultFit.Data()))) {
                tempFit       = (TF1*)list->FindObject(Form("%sStat_Fit_%s", fParticle[i].Data(), defaultFit.Data()));
                DrawFit(tempFit, 1, 1, GetParticleColor(fParticle[i]));
                tempFit->Draw("l,same");
            }
        }
    }
    dummyHisto->Draw("same,axis");

    // draw legend
    legend->Draw("same");
    labelEnergy->Draw();
    labelALICE->Draw();

    TString addName     = "";
    if (!xLogPlot)
        addName         = addName+"_LinX";
    if (xMaxPlot > -1000)
        addName         = addName+"_LowPt";
    // save canvas
    if (histCounter) canvas->SaveAs(Form("%s/ParticleSpectra_FinalWithFits%s.%s", outputDir.Data(), addName.Data(), suffix.Data()));

    // free pointer
    delete tempHist;
    delete legend;
    delete canvas;
    delete dummyHisto;
}


//================================================================================================================
//Function to produce plot containing all particle spectra in list
//================================================================================================================
void ProduceSingleParticleSpectraPlotFromListOnlyFinalWithFits( TList* list,
                                                                TString collSys,
                                                                TString energy,
                                                                TString cent,
                                                                TString suffix,
                                                                TString outputDir,
                                                                TString particleName,
                                                                Bool_t  xLogPlot        = kTRUE,
                                                                Double_t xMaxPlot       = -1000
) {

    // create output dir
    gSystem->Exec("mkdir -p "+outputDir);
    cout << collSys.Data() << "\t"<< energy.Data() << "\t" << outputDir.Data() << endl;

    // get xRange
    Double_t xMin                           = GetXRangeExtremaFromList(list, "spectra", kFALSE);
    Double_t xMax                           = GetXRangeExtremaFromList(list, "spectra", kTRUE);
    cout << "Plotting from: " << xMin << "\t" << xMax << endl;
    if (xMin < 0.01)
        xMin                                = 0.035;
    if (xMaxPlot != -1000)
        xMax                                = xMaxPlot;

    // get yRange
    Double_t yMin                           = GetYRangeExtremaFromList(list, kTRUE, kFALSE, collSys) * 0.5;
    Double_t yMax                           = GetYRangeExtremaFromList(list, kTRUE, kTRUE, collSys) * 4;
    if (xLogPlot && collSys.Contains("Pb"))
        yMax                                = yMax*2;

    // create canvas
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 2000);
    DrawCanvasSettings(canvas, 0.152, 0.015, 0.015, 0.068);
    canvas->SetLogy();
    if (xLogPlot)canvas->SetLogx();

    // count number of stat histograms in list
    Int_t histCounter                       = 0;
    TString tempName                        = "";
    Double_t yMinFit                        = 1e8;
    cout << "trying to find plots for " << particleName.Data() << endl;
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                            = ((TObject*)list->At(i))->GetName();
//         cout << tempName.Data() << endl;
        if ( tempName.Contains( Form("%sStat",particleName.Data()) ) ) {
            if (!tempName.Contains("To")){
                histCounter++;
//                 cout << "accepted: "<< tempName.Data() << endl;
                if (tempName.Contains("Fit")) {
                    if (xMaxPlot != -1000){
                        TF1* tempFit2               = (TF1*)list->FindObject(tempName.Data());
                        Double_t cYM                = tempFit2->Eval(xMaxPlot);
                        if (cYM < yMinFit)
                            yMinFit                 = cYM;
                        if (yMinFit != 1e8)
                            yMin                    = yMinFit*0.5;
                    }
                }
            } else {
                continue;
            }
        } else {
            continue;
        }
    }
    cout << "counted: " << histCounter << " inputs"<< endl;

    // get latex for enery and centrality
    Int_t energyIterator                    = -1;
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            energyIterator                  = i;
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    Int_t centIterator                      = -1;
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0){
            centIterator                    = i;
            cent                            = fCentralityLatex[i];
        }
    }

    // get number of rows for legend
    Int_t nRows                             = histCounter;
    // create legend
    TLegend* legend                         = NULL;
    if (xLogPlot ){
        legend                              = GetAndSetLegend(0.20, 0.09, 0.41, 0.09+(1.05*0.04*nRows), "" ,1,0.28,0.04);
    } else if (xMax > 10 ){
        legend                              = GetAndSetLegend(0.54, 0.885-(1.05*0.04*nRows), 0.95, 0.885, "" ,1,0.14,0.04);
    } else {
        legend                              = GetAndSetLegend(0.20, 0.09, 0.41, 0.09+(1.05*0.04*nRows), "" ,1,0.28,0.04);
    }
    TLatex *labelEnergy                     = new TLatex(0.94, 0.925, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);

    TString textALICE                       = "ALICE";
    if (isThesis) textALICE                 = "ALICE this thesis";

    TLatex *labelALICE                     = new TLatex(0.94, 0.885, textALICE);
    SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);

    //     Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data())
    // dummy histogram
    TH1D* dummyHisto                        = NULL;

    if (xLogPlot){
        dummyHisto                          = new TH1D("dummyHisto", "", 100000, xMin*0.8, xMax*1.5);
        SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", yMin, yMax, 0.7, 1.65, 0.04, kFALSE, kTRUE);
        dummyHisto->GetXaxis()->SetLabelOffset(-0.015);
    } else {
        dummyHisto                          = new TH1D("dummyHisto", "", 100000, 0, xMax*1.5);
        SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", yMin, yMax, 0.73, 1.65, 0.04, kFALSE, kTRUE);
        dummyHisto->GetXaxis()->SetLabelOffset(-0.001);
        dummyHisto->GetXaxis()->SetRangeUser(0,xMax);
    }
    dummyHisto->GetXaxis()->SetTickLength(0.02);

    dummyHisto->Draw();

    // search list for particle spectra
    TH1D* tempHist                          = NULL;
    TGraph* tempGraph                       = NULL;
    TGraphErrors* tempGraphErrs             = NULL;
    TGraphAsymmErrors* tempGraphAsymErrs    = NULL;
    TString tempClass                       = "";
    TF1* tempFit                            = NULL;

    if (list->FindObject(Form("%sSys", particleName.Data()))) {
        tempClass                   = ((TObject*)list->FindObject(Form("%sSys", particleName.Data())))->ClassName();
//         cout << Form("%sSys", particleName.Data()) << "\t"<< tempClass.Data() << endl;
//         cout << GetParticleMarkerStyle(particleName) << "\t" << GetParticleMarkerSize(particleName) << "\t" << GetParticleColor(particleName) << endl;
        if (tempClass.Contains("TH1")) {
            tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%sSys", particleName.Data())))->Clone("tempHist");
            TGraphErrors* tempGrC   = new TGraphErrors(tempHist);
            DrawMarker(tempGrC, GetParticleMarkerStyle(particleName), GetParticleMarkerSize(particleName), GetParticleColor(particleName),GetParticleColor(particleName),0);
            tempGrC->Draw("e2");
        } else if (tempClass.CompareTo("TGraph") == 0) {
            tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%sSys", particleName.Data())))->Clone("tempGraph");
            DrawMarker(tempGraph, GetParticleMarkerStyle(particleName), GetParticleMarkerSize(particleName), GetParticleColor(particleName),GetParticleColor(particleName),0);
            tempGraph->DrawClone("e2");
        } else if (tempClass.CompareTo("TGraphErrors") == 0) {

            tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%sSys", particleName.Data())))->Clone("tempGraphErrs");
            DrawMarker(tempGraphErrs, GetParticleMarkerStyle(particleName), GetParticleMarkerSize(particleName), GetParticleColor(particleName),GetParticleColor(particleName),0);
            tempGraphErrs->DrawClone("e2");
        } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
            tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%sSys", particleName.Data())))->Clone(Form("%sSys", particleName.Data()));
            DrawMarker(tempGraphAsymErrs, GetParticleMarkerStyle(particleName), GetParticleMarkerSize(particleName), GetParticleColor(particleName),GetParticleColor(particleName),0);
            tempGraphAsymErrs->DrawClone("e2");
        }
    }
    if (list->FindObject(Form("%sStat", particleName.Data()))) {
        tempClass                   = ((TObject*)list->FindObject(Form("%sStat", particleName.Data())))->ClassName();
//         cout << Form("%sStat", particleName.Data()) << "\t"<< tempClass.Data() << endl;
        if (tempClass.Contains("TH1")) {
            tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%sStat", particleName.Data())))->Clone("tempHist");
            DrawMarker(tempHist, GetParticleMarkerStyle(particleName), GetParticleMarkerSize(particleName), GetParticleColor(particleName),GetParticleColor(particleName));
            tempHist->Draw("e1,same");
            legend->AddEntry(tempHist, Form("%s", fParticleLatex[GetParticleIterator(particleName)].Data()), "lp");
        } else if (tempClass.CompareTo("TGraph") == 0) {
            tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%sStat", particleName.Data())))->Clone("tempGraph");
            DrawMarker(tempGraph, GetParticleMarkerStyle(particleName), GetParticleMarkerSize(particleName), GetParticleColor(particleName),GetParticleColor(particleName));
            tempGraph->DrawClone("p");
            legend->AddEntry(tempGraph, Form("%s", fParticleLatex[GetParticleIterator(particleName)].Data()), "lp");
        } else if (tempClass.CompareTo("TGraphErrors") == 0) {
            tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%sStat", particleName.Data())))->Clone("tempGraphErrs");
            ProduceGraphWithoutXErrors(tempGraphErrs);
            DrawMarker(tempGraphErrs, GetParticleMarkerStyle(particleName), GetParticleMarkerSize(particleName), GetParticleColor(particleName),GetParticleColor(particleName));
            tempGraphErrs->DrawClone("p");
            legend->AddEntry(tempGraphErrs, Form("%s", fParticleLatex[GetParticleIterator(particleName)].Data()), "lp");
        } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
            tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%sStat", particleName.Data())))->Clone(Form("%sStat", particleName.Data()));
            ProduceGraphAsymmWithoutXErrors(tempGraphAsymErrs);
            DrawMarker(tempGraphAsymErrs, GetParticleMarkerStyle(particleName), GetParticleMarkerSize(particleName), GetParticleColor(particleName),GetParticleColor(particleName));
            tempGraphAsymErrs->DrawClone("p");
            legend->AddEntry(tempGraphAsymErrs, Form("%s", fParticleLatex[GetParticleIterator(particleName)].Data()), "lp");
        }
    }

    Int_t nFits                             = 0;
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                            = ((TObject*)list->At(i))->GetName();
        if ( tempName.Contains( Form("%sStat",particleName.Data()) ) && tempName.Contains("Fit") ) {
            TObjArray *arr                  = tempName.Tokenize("_");
            TString nameFit                 = ((TObjString*)arr->At(2))->GetString();
            tempFit                         = (TF1*)list->FindObject(tempName.Data());
            DrawFit(tempFit, lineStyle[nFits], 1, lineColor[nFits]);
            tempFit->Draw("l,same");
            legend->AddEntry(tempFit, GetFitLabel(nameFit), "l");
            nFits++;
        }
    }

    dummyHisto->Draw("same,axis");

    // draw legend
    legend->Draw("same");
    labelEnergy->Draw();
    labelALICE->Draw();

    TString addName     = "";
    if (!xLogPlot)
        addName         = addName+"_LinX";
    if (xMaxPlot > -1000)
        addName         = addName+"_LowPt";
    // save canvas
    if (histCounter) canvas->SaveAs(Form("%s/Spectra_%s%s.%s", outputDir.Data(), particleName.Data(), addName.Data(), suffix.Data()));

    // free pointer
    delete tempHist;
    delete legend;
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot containing all particle ratio in list
//================================================================================================================
void ProduceParticleRatioPlotFromList(TList* list, TString collSys, TString energy, TString cent, TString suffix) {

    // create output dir
    TString outputDir                                   = Form("plots/%s/%s",suffix.Data(),list->GetName());
    gSystem->Exec("mkdir -p "+outputDir);

    // get xRange
    Double_t xMin                                       = 0;
    Double_t xMax                                       = GetXRangeExtremaFromList(list, "ratio", kTRUE);

    // get yRange
    Double_t yMin                                       = 0.;
    Double_t yMax                                       = 2.0;

    // create canvas
    TCanvas* canvas                                     = GetAndSetCanvas("canvas");
    DrawCanvasSettings(canvas, 0.11, 0.015, 0.02, 0.09);
    //     canvas->SetLogx();

    // count number of stat histograms in list
    Int_t histCounter                                   = 0;
    Bool_t hasNPionToCPion                              = kFALSE;
    TString tempName                                    = "";
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                                        = ((TObject*)list->At(i))->GetName();
        if (tempName.Contains("Stat") && tempName.Contains("NPionToCPion"))
            hasNPionToCPion                             = kTRUE;
        if (tempName.Contains("Stat") && tempName.Contains("To"))
            histCounter++;
        else
            continue;
    }
    // adjust range if pi0/pi+- is available
    if (hasNPionToCPion)
        yMax                                            = 2.3;
    // get number of rows for legend
    Int_t nRows                                         = 0;
    Int_t nColumns                                      = 2;
    if (histCounter > 9)
        nColumns++;

    if (histCounter%nColumns == 0) nRows                = histCounter/nColumns;
    else nRows                                          = (histCounter+1)/nColumns;
    // get latex for enery and centrality
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                                  = fEnergyLatexPP[i];
            else
                energy                                  = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                                        = fCentralityLatex[i];
    }

    // create legend
    TLegend* legend                                     = GetAndSetLegend(0.7-0.2*(nColumns-1), 0.94-(0.052*nRows), 0.9, 0.94, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()),
                                                                          nColumns, 0.14 ,40,43);

    // dummy histogram
    TH1D* dummyHisto                                    = new TH1D("dummyHisto", "", 1000, xMin, xMax);
    SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","Particle ratio", yMin, yMax);
    dummyHisto->Draw();

    // search list for particle spectra
    TH1D* tempHist                                      = NULL;
    TGraph* tempGraph                                   = NULL;
    TGraphErrors* tempGraphErrs                         = NULL;
    TGraphAsymmErrors* tempGraphAsymErrs                = NULL;
    TString tempClass                                   = "";

    Int_t j                                             = 0;
    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t k=0; k<nParticles; k++) {
            for (Int_t ii=0; ii<9; ii++) {
                if(i == k ) continue;
                if (list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[k].Data(), fMethod[ii].Data()))) {
                    tempClass                           = ((TObject*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[k].Data(), fMethod[ii].Data())))->ClassName();
//                     cout << "found " << Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[k].Data(), fMethod[ii].Data()) << endl;
                    if (tempClass.Contains("TH1")) {
                        tempHist                        = (TH1D*)((TH1D*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[k].Data(), fMethod[ii].Data())))->Clone("tempHist");
                        if (tempHist) {
                            DrawMarker(tempHist, markers[j], markerSize[j], paint[j],paint[j]);
                            tempHist->Draw("e1,same");
                            legend->AddEntry(tempHist, Form("%s/%s %s", fRatioParticleLatex[i].Data(), fRatioParticleLatex[k].Data(), fMethodLabel[ii].Data()), "lp");
                            j++;
                        }
                    } else if (tempClass.CompareTo("TGraph") == 0) {
                        tempGraph                       = (TGraph*)((TGraph*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[k].Data(), fMethod[ii].Data())))->Clone("tempGraph");
                        if (tempGraph) {
                            DrawMarker(tempGraph, markers[j], markerSize[j], paint[j],paint[j]);
                            tempGraph->Draw("p");
                            legend->AddEntry(tempGraph, Form("%s/%s %s", fRatioParticleLatex[i].Data(), fRatioParticleLatex[k].Data(), fMethodLabel[ii].Data()), "lp");
                            j++;
                        }
                    } else if (tempClass.CompareTo("TGraphErrors") == 0) {
                        tempGraphErrs                   = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[k].Data(), fMethod[ii].Data())))->Clone("tempGraphErrs");
                        if (tempGraphErrs) {
                            DrawMarker(tempGraphErrs, markers[j], markerSize[j], paint[j],paint[j]);
                            tempGraphErrs->Draw("p");
                            legend->AddEntry(tempGraphErrs, Form("%s/%s %s", fRatioParticleLatex[i].Data(), fRatioParticleLatex[k].Data(), fMethodLabel[ii].Data()), "lp");
                            j++;
                        }
                    } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                        tempGraphAsymErrs               = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[k].Data(), fMethod[ii].Data())))->Clone("tempGraphAsymErrs");
                        if (tempGraphAsymErrs) {
                            DrawMarker(tempGraphAsymErrs, markers[j], markerSize[j], paint[j],paint[j]);
                            tempGraphAsymErrs->Draw("p");
                            legend->AddEntry(tempGraphAsymErrs, Form("%s/%s %s", fRatioParticleLatex[i].Data(), fRatioParticleLatex[k].Data(), fMethodLabel[ii].Data()), "lp");
                            j++;
                        }
                    }
                }
            }
        }
    }

    // draw legend
    legend->Draw("same");

    // save canvas
    if (histCounter) canvas->SaveAs(Form("%s/ParticleRatio.%s", outputDir.Data(), suffix.Data()));

    // free pointer
    delete tempHist;
    delete legend;
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot containing all particle ratio in list
//================================================================================================================
void ProduceParticleRatioPlotFromListOnlyFinal(TList* list, TString collSys, TString energy, TString cent, TString suffix, Bool_t xyLogPlot = kFALSE, Bool_t onlyMesons = kFALSE) {

    // create output dir
    TString outputDir                                   = Form("plots/%s/%s",suffix.Data(),list->GetName());
    gSystem->Exec("mkdir -p "+outputDir);

    // get xRange
    Double_t xMin                                       = 0;
    Double_t xMax                                       = GetXRangeExtremaFromList(list, "ratio", kTRUE);
    if(xyLogPlot){
        xMin                                            = GetXRangeExtremaFromList(list, "ratio", kFALSE) * 0.5;
        xMax                                            = xMax * 1.3;
    }

    // get yRange
    Double_t yMin                                       = 0;
    Double_t yMax                                       = 2;
    if(xyLogPlot){
        yMin                                            = GetYRangeExtremaFromList(list, kTRUE, kFALSE, collSys) * 0.5;
        if(yMin<0.001)
            yMin                                        = 0.0015;
        yMax                                            = GetYRangeExtremaFromList(list, kTRUE, kTRUE, collSys) * 2;
    }

    // create canvas
    TCanvas* canvas                                     = GetAndSetCanvas("canvas");
    DrawCanvasSettings(canvas, 0.07, 0.015, 0.02, 0.08);
    if(xyLogPlot){
        canvas->SetLogx();
        canvas->SetLogy();
    }

    // count number of stat histograms in list
    Int_t histCounter                                   = 0;
    TString tempName                                    = "";
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                                        = ((TObject*)list->At(i))->GetName();
        if (tempName.Contains("Stat") && tempName.Contains("To")){
//             cout << "testing: " << tempName.Data() << endl;
            if (tempName.Contains("EtaToNPion")){
                if (tempName.Contains("Comb")){
//                     cout << tempName.Data() << endl;
                    histCounter++;
                } else {
                    continue;
                }
            } else if (tempName.Contains("NPion") || tempName.Contains("Eta")){
                if (!tempName.Contains("Comb")){
//                     cout << tempName.Data() << endl;
                    histCounter++;
                } else {
                    continue;
                }
            } else {
//                 cout << tempName.Data() << endl;
                histCounter++;
            }
        } else {
            continue;
        }
    }

//     cout << "number of ratio hists: " << histCounter << "\t"<<  Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data() ) << endl;

    // get number of rows for legend
    Int_t nRows                                         = 0;
    Int_t nColumns                                      = 2;
    if (histCounter > 12){
        nColumns                                        = 5;
    } else if (histCounter > 10){
        nColumns                                        = 4;
    } else if (histCounter > 6){
        nColumns                                        = 3;
    }
    if (histCounter%nColumns == 0) nRows                = histCounter/nColumns;
    else nRows                                          = (histCounter + nColumns-(histCounter%nColumns) )/nColumns;
//     cout << histCounter%nColumns << "\t" << histCounter/nColumns << "\t" << (histCounter + nColumns-(histCounter%nColumns) )/nColumns  << "\t" << nRows << endl;

    // get latex for enery and centrality
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                                  = fEnergyLatexPP[i];
            else
                energy                                  = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                                        = fCentralityLatex[i];
    }

    // create legend
    TLegend* legend                         = GetAndSetLegend(0.11, 0.95, 0.59+0.1*(nColumns-2), 0.95-(1.05*0.04*nRows), "" ,nColumns,(0.59+0.1*(nColumns-2))*0.2,0.04);
    TLatex *labelEnergy                     = new TLatex(0.94, 0.152, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);
    TString textALICE                       = "ALICE";
    if (isThesis) textALICE                 = "ALICE this thesis";

    TLatex *labelALICE                     = new TLatex(0.94, 0.112, textALICE);
    SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);

    // dummy histogram
    TH1D* dummyHisto                                    = new TH1D("dummyHisto", "", 1000, xMin, xMax);
    SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","Particle ratio", yMin, yMax, 0.9, 0.9, 0.04, kFALSE, kTRUE);
    dummyHisto->Draw();

    // search list for particle spectra
    TH1D* tempHist                          = NULL;
    TGraph* tempGraph                       = NULL;
    TGraphErrors* tempGraphErrs             = NULL;
    TGraphAsymmErrors* tempGraphAsymErrs    = NULL;
    TString tempClass                       = "";

    Int_t k                                 = 0;
    for (Int_t i=0; i<nParticles; i++) {
        for (Int_t j=0; j<nParticles; j++) {
            for (Int_t ii=0; ii<2; ii++) {
                if (ii == 1 && i > 6) continue;
                if (onlyMesons){
                    if (!(fIsMeson[i] && fIsMeson[j]))
                        continue;
                }
                if (list->FindObject(Form("%sTo%s%sSys", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data()))) {
                    tempClass                   = ((TObject*)list->FindObject(Form("%sTo%s%sSys", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->ClassName();
//                     cout << Form("%sTo%s%sSys", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data()) << "\t"<< tempClass.Data() << endl;
//                     cout << GetParticleMarkerStyle(fParticle[i]) << "\t" << GetParticleMarkerSize(fParticle[j]) << "\t" << GetParticleColor(fParticle[i]) << endl;
                    if (tempClass.Contains("TH1")) {
                        tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%sTo%s%sSys", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->Clone("tempHist");
                        if(tempHist){
                            TGraphErrors* tempGrC   = new TGraphErrors(tempHist);
                            DrawMarker(tempGrC, markers[k], markerSize[k], paint[k],paint[k],0);
                            tempGrC->Draw("e2");
                        }
                    } else if (tempClass.CompareTo("TGraph") == 0) {
                        tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%sTo%s%sSys", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->Clone("tempGraph");
                        if(tempGraph){
                            DrawMarker(tempGraph, markers[k], markerSize[k], paint[k],paint[k],0);
                            tempGraph->DrawClone("e2");
                        }
                    } else if (tempClass.CompareTo("TGraphErrors") == 0) {
                        tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%sTo%s%sSys", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->Clone("tempGraphErrs");
                        if(tempGraphErrs){
                            DrawMarker(tempGraphErrs, markers[k], markerSize[k], paint[k],paint[k],0);
                            tempGraphErrs->DrawClone("e2");
                        }
                   } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                        tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%sTo%s%sSys", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->Clone(Form("%sTo%s%sSys", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data()));
                        if(tempGraphAsymErrs){
                            DrawMarker(tempGraphAsymErrs, markers[k], markerSize[k], paint[k],paint[k],0);
                            tempGraphAsymErrs->DrawClone("e2");
                        }
                    }
                }
                if (list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data()))) {
                    tempClass                   = ((TObject*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->ClassName();
//                     cout << Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data()) << "\t"<< tempClass.Data() << endl;
//                     cout << GetParticleMarkerStyle(fParticle[i]) << "\t" << GetParticleMarkerSize(fParticle[j]) << "\t" << GetParticleColor(fParticle[i]) << endl;
                    if (tempClass.Contains("TH1")) {
                        tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->Clone("tempHist");
                        if(tempHist){
                            TGraphErrors* tempGrC   = new TGraphErrors(tempHist);
                            DrawMarker(tempGrC, markers[k], markerSize[k], paint[k],paint[k]);
                            tempGrC->Draw("pX0Z");
                            k++;
                            legend->AddEntry(tempGrC, Form("%s/%s", fRatioParticleLatex[i].Data(),fRatioParticleLatex[j].Data()), "lp");
                        }
                    } else if (tempClass.CompareTo("TGraph") == 0) {
                        tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->Clone("tempGraph");
                        if(tempGraph){
                            DrawMarker(tempGraph, markers[k], markerSize[k], paint[k],paint[k]);
                            tempGraph->DrawClone("pZ");
                            k++;
                            legend->AddEntry(tempGraph, Form("%s/%s", fRatioParticleLatex[i].Data(),fRatioParticleLatex[j].Data()), "lp");
                        }
                    } else if (tempClass.CompareTo("TGraphErrors") == 0) {
                        tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->Clone("tempGraphErrs");
                        ProduceGraphWithoutXErrors(tempGraphErrs);
                        if(tempGraphErrs){
                            DrawMarker(tempGraphErrs, markers[k], markerSize[k], paint[k],paint[k]);
                            tempGraphErrs->DrawClone("pZ");
                            k++;
                            legend->AddEntry(tempGraphErrs, Form("%s/%s", fRatioParticleLatex[i].Data(),fRatioParticleLatex[j].Data()), "lp");
                        }
                   } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {
                        tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data())))->Clone(Form("%sTo%s%sStat", fParticle[i].Data(), fParticle[j].Data(), fMethod[ii].Data()));
                        ProduceGraphAsymmWithoutXErrors(tempGraphAsymErrs);
                        if(tempGraphAsymErrs){
                            DrawMarker(tempGraphAsymErrs, markers[k], markerSize[k], paint[k],paint[k]);
                            tempGraphAsymErrs->DrawClone("pZ");
                            k++;
                            legend->AddEntry(tempGraphAsymErrs, Form("%s/%s", fRatioParticleLatex[i].Data(),fRatioParticleLatex[j].Data()), "lp");
                        }
                    }
                }
            }
        }
    }
    dummyHisto->Draw("same,axis");

    // draw legend
    legend->Draw("same");
    labelEnergy->Draw();
    labelALICE->Draw();

    // save canvas
    TString nameAdd     = "";
    if (onlyMesons)
        nameAdd         = nameAdd+"Mesons";
    if (xyLogPlot)
        nameAdd         = nameAdd+"LogLog";
    if (histCounter) canvas->SaveAs(Form("%s/ParticleRatios_Final%s.%s", outputDir.Data(), nameAdd.Data(), suffix.Data()));
    // free pointer
    delete tempHist;
    delete legend;
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot comparing to fit
//================================================================================================================
void ProducePlotWithRatioToFit( TH1D* histo1,
                                TH1D* histo2                = NULL,
                                TH1D* histo3                = NULL,
                                TF1* histo1Fit              = NULL,
                                TF1* histo2Fit              = NULL,
                                TF1* histo3Fit              = NULL,
                                TH1D* histoRatioToFit1      = NULL,
                                TH1D* histoRatioToFit2      = NULL,
                                TH1D* histoRatioToFit3      = NULL,
                                TString collSys             = "",
                                TString energy              = "",
                                TString cent                = "",
                                TString particle            = "",
                                TString method              = "",
                                TString fitLabel            = "",
                                TString suffix              = "",
                                TString xTitle              = "",
                                TString yTitle              = "",
                                TString namePlot            = "",
                                Double_t xMinFit            = 0.,
                                Double_t xMaxFit            = 15.,
                                Bool_t isRatioPlot          = kFALSE,
                                Bool_t secondIsSys          = kFALSE,
                                Double_t xMinPlot           = -1000,
                                Double_t xMaxPlot           = -1000,
                                Double_t yMinPlotRatio      = 0.5,
                                Double_t yMaxPlotRatio      = 1.5,
                                TString fExternalOutputDir  = "",
                                Int_t plotLogXaxisRatios    = 0
) {

    // create output dir
    TString outputDir;
    if (fExternalOutputDir.CompareTo("") != 0){
        outputDir                           = fExternalOutputDir;
    } else {
        if (namePlot.Contains("default")){
            if (cent.CompareTo("") == 0) outputDir  = Form("plots/%s/%s_%s",suffix.Data(),collSys.Data(), energy.Data());
            else outputDir                          = Form("plots/%s/%s_%s_%s",suffix.Data(),collSys.Data(), energy.Data(), cent.Data());
        } else {
            if (cent.CompareTo("") == 0) outputDir  = Form("plots/%s/%s_%s/SysParam",suffix.Data(),collSys.Data(), energy.Data());
            else outputDir                          = Form("plots/%s/%s_%s_%s/SysParam",suffix.Data(),collSys.Data(), energy.Data(), cent.Data());
        }
    }
    gSystem->Exec("mkdir -p "+outputDir);

    // get yRange
    TList* tempList                         = new TList;
    tempList->Add(histo1);
    tempList->Add(histo2);
    tempList->Add(histo3);
    Double_t yMin                           = GetYRangeExtremaFromList(tempList, kFALSE, kFALSE, collSys) * 0.2;
    Double_t yMax                           = GetYRangeExtremaFromList(tempList, kFALSE, kTRUE, collSys) * 2.5;
    // get xRange
    Double_t xMin                           = 0;
    Double_t xMax                           = GetXRangeExtremaFromList(tempList, kFALSE, kTRUE);
    if (xMinPlot != -1000 && xMaxPlot != -1000){
        xMin                                = xMinPlot;
        xMax                                = xMaxPlot;
    }

    if (isRatioPlot) {
        yMin                                = yMin * 0.8;
        yMax                                = yMax * 1.05;
    } else {
        yMin                                = yMin * 0.8;
        yMax                                = yMax * 1.2;
        if (xMinPlot != -1000 && xMaxPlot != -1000 && histo1Fit)
            yMin                            = histo1Fit->Eval(xMaxPlot)*0.8;
    }

    // get latex for energy, centrality and particle
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }
    for (Int_t i=0; i<nParticles; i++) {
        if (particle.CompareTo(fParticle[i]) == 0)
            particle                        = fParticleLatex[i];
    }

    // specific setting for different ratio plots
    if (isRatioPlot){
        if (((TString)histo1->GetYaxis()->GetTitle()).Contains("#eta/#pi^{0}")) {
            yMin                                = 0.01;
            yMax                                = 0.9;
        }
    }

    // sys label
    TString sysLabelA                       = "";
    TString sysLabelB                       = "";
    if (namePlot.Contains("ptConst")) {
        sysLabelA                           = "const up sys";
        sysLabelB                           = "const down sys";
    }
    if (namePlot.Contains("ptLin")) {
        sysLabelA                           = "lin A sys";
        sysLabelB                           = "lin B sys";
    }
    if (namePlot.Contains("ptPol2")) {
        sysLabelA                           = "pol2 A sys";
        sysLabelB                           = "pol2 B sys";
    }

    // fit label
    if (!fitLabel.CompareTo("oHag"))
        fitLabel                            = "mod. Hagedorn";
    else if (!fitLabel.CompareTo("oHagPt"))
        fitLabel                            = "mod. Hagedorn * #it{p}_{T}";
    else if (!fitLabel.CompareTo("h"))
        fitLabel                            = "Hagedorn";
    else if (!fitLabel.CompareTo("l"))
        fitLabel                            = "Tsallis";
    else if (!fitLabel.CompareTo("softHard"))
        fitLabel                            = "soft + hard";

    // label for histo legend
    TString label                           = "";
    if (isRatioPlot) label                  = histo1->GetYaxis()->GetTitle();
    else label                              = particle;
    label                                   = Form("%s %s", label.Data(), method.Data());

    // chi2 red of fit
    Double_t chi2red1, chi2red2, chi2red3;
    chi2red1                                = histo1Fit->GetChisquare() / histo1Fit->GetNDF();
    chi2red2                                = 0;
    if(histo2Fit) chi2red2                  = histo2Fit->GetChisquare() / histo2Fit->GetNDF();
    chi2red3                                = 0.;
    if(histo3Fit) chi2red3                  = histo3Fit->GetChisquare() / histo3Fit->GetNDF();

    // create canvas and pads
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 2000);
    DrawCanvasSettings(canvas, 0.0, 0.0, 0.0, 0.0);
    TPad* pad                               = new TPad("pad", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawPadSettings(pad, 0.13, 0.015, 0.02, 0.);
    if(plotLogXaxisRatios) pad->SetLogx();
    pad->Draw();
    TPad* padRatioToFit                     = new TPad("padRatioToFit", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawPadSettings(padRatioToFit, 0.13, 0.015, 0.0, 0.25);
    if(plotLogXaxisRatios) padRatioToFit->SetLogx();
    padRatioToFit->Draw();


    if (xTitle.CompareTo("") == 0) xTitle   = histo1->GetXaxis()->GetTitle();
    if (yTitle.CompareTo("") == 0) yTitle   = histo1->GetYaxis()->GetTitle();
    if (xMinPlot != -1000 && xMaxPlot != -1000)
        histo1->GetXaxis()->SetRangeUser(xMinPlot, xMaxPlot);

    DrawMarker(histo1, 20, 1.5, kBlack, kBlack);
    if (histo2 && secondIsSys){
        SetHistogramm(histo2, xTitle, yTitle, yMin, yMax);
        DrawMarker(histo2, 20, 1.5, kBlack, kBlack, 0);
    } else if (histo2){
        SetHistogramm(histo2, xTitle, yTitle, yMin, yMax);
        DrawMarker(histo2, 24, 1.5, kRed-4, kRed-4);
    }
    if (histo3) {
        SetHistogramm(histo2, xTitle, yTitle, yMin, yMax);
        DrawMarker(histo3, 24, 1.5, kBlue-4, kBlue-4);
    }

    // plot histos + fits
    pad->cd();
    if (!isRatioPlot) pad->SetLogy();

    TH1D* dummyHisto                        = new TH1D("dummyHisto", "", 100000, xMin, xMax);
    SetHistogramm(dummyHisto, xTitle, yTitle, 0.8*yMin, 1.2*yMax, 0.7, 2.4, 60, kFALSE, kTRUE, 43);
    dummyHisto->GetYaxis()->SetTickLength(0.02);
    dummyHisto->GetXaxis()->SetTickLength(0.02);
    dummyHisto->DrawCopy();

    if (histo1Fit) {
        histo1Fit->SetNpx(10000);
        histo1Fit->SetLineColor(kBlack);
        histo1Fit->SetLineWidth(1);
    }
    if (histo2Fit) {
        histo2Fit->SetNpx(10000);
        histo2Fit->SetLineColor(kRed-4);
        histo2Fit->SetLineWidth(1);
    }
    if (histo2Fit) {
        histo3Fit->SetNpx(10000);
        histo3Fit->SetLineColor(kBlue-4);
        histo3Fit->SetLineWidth(1);
    }

    histo1->Draw("same,e1");
    if (secondIsSys && histo2) histo2->Draw("e2,same");
    else if (histo2) histo2->Draw("e1,same");
    if (histo3) histo3->Draw("e1,same");

    if (histo1Fit) histo1Fit->Draw("same");
    if (histo2Fit) histo2Fit->Draw("same");
    if (histo3Fit) histo3Fit->Draw("same");

    // legend histo
    TLegend* legendHisto                    = NULL;
    if (isRatioPlot){
        if (((TString)histo1->GetYaxis()->GetTitle()).Contains("#eta/#pi^{0}"))
            legendHisto                     = GetAndSetLegend(0.45, 0.02, 0.45+0.25, 0.02+4*0.05, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()), 1, 0.14 ,60,43);
        else
            legendHisto                     = GetAndSetLegend(0.16, 0.95-4*0.05, 0.16+0.25, 0.95, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()), 1, 0.14 ,50,43);
    } else {
        legendHisto                         = GetAndSetLegend(0.16, 0.02, 0.16+0.25, 0.02+4*0.05, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()), 1, 0.14 ,60,43);
    }
    legendHisto->AddEntry(histo1, Form("%s", label.Data()), "p");
    if (histo2 && !secondIsSys) legendHisto->AddEntry(histo2, Form("%s %s", label.Data(), sysLabelA.Data()), "p");
    if (histo3) legendHisto->AddEntry(histo3, Form("%s %s", label.Data(), sysLabelB.Data()), "p");
    legendHisto->Draw();

    // legend fit
    TLegend* legendFit                      = NULL;
    if (isRatioPlot){
        if (label.Contains("#eta/#pi^{0}"))
            legendFit                           = GetAndSetLegend(0.16, 0.95-5*0.05, 0.16+0.25, 0.95, "", 1, 0.14, 60, 43);
        else
            legendFit                           = GetAndSetLegend(0.55, 0.95-5*0.05, 0.55+0.25, 0.95, "", 1, 0.14, 50, 43);
    } else {
        legendFit                               = GetAndSetLegend(0.45, 0.95-5*0.05, 0.45+0.25, 0.95, "", 1, 0.14, 60, 43);
    }

    legendFit->SetBorderSize(0);
    legendFit->AddEntry((TObject*)0, Form("%s", fitLabel.Data()), "");
    legendFit->AddEntry((TObject*)0, Form("#it{p}_{T} = %.2f - %.2f GeV/#it{c}", xMinFit, xMaxFit), "");
    if (!isThesis){
        legendFit->AddEntry(histo1Fit, Form("standard, #chi^{2}/ndf = %.2f", chi2red1), "l");
        if (histo2Fit) legendFit->AddEntry(histo2Fit, Form("%s, #chi^{2}/ndf = %.2f", sysLabelA.Data(), chi2red2), "l");
        if (histo3Fit) legendFit->AddEntry(histo3Fit, Form("%s, #chi^{2}/ndf = %.2f", sysLabelB.Data(), chi2red3), "l");
    } else {
        legendFit->AddEntry(histo1Fit, Form("standard"), "l");
        if (histo2Fit) legendFit->AddEntry(histo2Fit, Form("%s", sysLabelA.Data()), "l");
        if (histo3Fit) legendFit->AddEntry(histo3Fit, Form("%s", sysLabelB.Data()), "l");
    }
    legendFit->Draw();

    // plot ratio histo/fit
    padRatioToFit->cd();

    TH1D* dummyHistoRatio                        = new TH1D("dummyHistoRatio", "", 100000, xMin, xMax);
    SetHistogramm(dummyHistoRatio, "#it{p}_{T} (GeV/#it{c})","#frac{data}{fit}", yMinPlotRatio, yMaxPlotRatio, 3.4, 8.1, 60, kFALSE, kTRUE, 43, 510, 505);
    if(plotLogXaxisRatios>1) dummyHistoRatio->GetYaxis()->SetRangeUser(0.81,1.19);
    if (namePlot.Contains("ptConst") || namePlot.Contains("default") ) dummyHistoRatio->GetYaxis()->SetRangeUser(0.79,1.22);
    dummyHistoRatio->DrawCopy();

    DrawMarker(histoRatioToFit1, 20, 1.5, kBlack, kBlack, 43, 43);

    if (xMinPlot != -1000 && xMaxPlot != -1000)
        histoRatioToFit1->GetXaxis()->SetRangeUser(xMinPlot, xMaxPlot);


    if (histoRatioToFit2 && secondIsSys)
        DrawMarker(histoRatioToFit2, 20, 1.5, kBlack, kBlack, 0);
    else if (histoRatioToFit2)
        DrawMarker(histoRatioToFit2, 24, 1.5, kRed-4, kRed-4);
    if (histoRatioToFit3) DrawMarker(histoRatioToFit3, 24, 1.5, kBlue-4, kBlue-4);

    histoRatioToFit1->Draw("same,e1");

    DrawLine(xMin, xMax, 1, 1, 1, kGray+2, 2);
    DrawLine(xMin, xMax, 1.1, 1.1, 1, kGray+2, 3);
    DrawLine(xMin, xMax, 0.9, 0.9, 1, kGray+2, 3);

    if (secondIsSys && histoRatioToFit2){
        histoRatioToFit2->Draw("e2,same");
        histoRatioToFit1->Draw("e1,same");
    } else if (histoRatioToFit2) {
        histoRatioToFit2->Draw("e1,same");
    }
    if (histoRatioToFit3) histoRatioToFit3->Draw("e1,same");

    canvas->SaveAs(Form("%s/%s.%s", outputDir.Data(), namePlot.Data(), suffix.Data()));

    histo1->GetXaxis()->UnZoom();
    histo1->GetYaxis()->UnZoom();
    histoRatioToFit1->GetXaxis()->UnZoom();

    // free pointer
    delete legendHisto;
    delete legendFit;
    delete dummyHisto;
    delete dummyHistoRatio;
    delete pad;
    delete padRatioToFit;
    delete canvas;
}

void ProducePlotWithRatioToFit( TGraph* graph1,
                                TGraph* graph2               = NULL,
                                TGraph* graph3               = NULL,
                                TF1* graph1Fit               = NULL,
                                TF1* graph2Fit               = NULL,
                                TF1* graph3Fit               = NULL,
                                TH1D* histoRatioToFit1       = NULL,
                                TH1D* histoRatioToFit2       = NULL,
                                TH1D* histoRatioToFit3       = NULL,
                                TString collSys              = "",
                                TString energy               = "",
                                TString cent                 = "",
                                TString particle             = "",
                                TString method               = "",
                                TString fitLabel             = "",
                                TString suffix               = "",
                                TString xTitle               = "",
                                TString yTitle               = "",
                                TString namePlot             = "",
                                Double_t xMinFit             = 0.,
                                Double_t xMaxFit             = 15.,
                                Bool_t isRatioPlot           = kFALSE,
                                Bool_t secondIsSys           = kFALSE,
                                Double_t xMinPlot            = -1000,
                                Double_t xMaxPlot           = -1000,
                                Double_t yMinPlotRatio      = 0.5,
                                Double_t yMaxPlotRatio      = 1.5,
                                TString fExternalOutputDir  = "",
                                Int_t plotLogXaxisRatios    = 0
) {

    // create output dir
    TString outputDir;
    if (fExternalOutputDir.CompareTo("") != 0){
        outputDir                           = fExternalOutputDir;
    } else {
        if (namePlot.Contains("default")){
            if (cent.CompareTo("") == 0) outputDir  = Form("plots/%s/%s_%s",suffix.Data(),collSys.Data(), energy.Data());
            else outputDir                          = Form("plots/%s/%s_%s_%s",suffix.Data(),collSys.Data(), energy.Data(), cent.Data());
        } else {
            if (cent.CompareTo("") == 0) outputDir  = Form("plots/%s/%s_%s/SysParam",suffix.Data(),collSys.Data(), energy.Data());
            else outputDir                          = Form("plots/%s/%s_%s_%s/SysParam",suffix.Data(),collSys.Data(), energy.Data(), cent.Data());
        }
    }
    gSystem->Exec("mkdir -p "+outputDir);

    // get yRange
    TList* tempList                         = new TList;
    tempList->Add(graph1);
    tempList->Add(graph2);
    tempList->Add(graph3);

    // get xRange
    Double_t xMin                           = 0;
    Double_t xMax                           = GetXRangeExtremaFromList(tempList, kFALSE, kTRUE);
    if (xMinPlot != -1000 && xMaxPlot != -1000){
        xMin                                = xMinPlot;
        xMax                                = xMaxPlot;
    }
    // get yRange
    Double_t yMin                           = GetYRangeExtremaFromList(tempList, kFALSE, kFALSE, collSys) * 0.2;
    Double_t yMax                           = GetYRangeExtremaFromList(tempList, kFALSE, kTRUE, collSys) * 2.5;

    if (isRatioPlot) {
        yMin                                = yMin * 0.8;
        yMax                                = yMax * 1.05;
    } else {
        yMin                                = yMin * 0.8;
        yMax                                = yMax * 1.2;
        if (xMinPlot != -1000 && xMaxPlot != -1000 && graph1Fit)
            yMin                            = graph1Fit->Eval(xMax)*0.8;
    }

    // specific setting for different ratio plots
    if (isRatioPlot){
        if (((TString)graph1->GetYaxis()->GetTitle()).Contains("#eta/#pi^{0}")) {
            yMin                                = 0.01;
            yMax                                = 0.9;
        }
    }

    // get latex for energy, centrality and particle
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }
    for (Int_t i=0; i<nParticles; i++) {
        if (particle.CompareTo(fParticle[i]) == 0)
            particle                        = fParticleLatex[i];
    }

    // sys label
    TString sysLabelA                       = "";
    TString sysLabelB                       = "";
    if (namePlot.Contains("ptConst")) {
        sysLabelA                           = "const up sys";
        sysLabelB                           = "const down sys";
    }
    if (namePlot.Contains("ptLin")) {
        sysLabelA                           = "lin A sys";
        sysLabelB                           = "lin B sys";
    }
    if (namePlot.Contains("ptPol2")) {
        sysLabelA                           = "pol2 A sys";
        sysLabelB                           = "pol2 B sys";
    }

    // fit label
    if (!fitLabel.CompareTo("oHag"))
        fitLabel                            = "mod. Hagedorn";
    else if (!fitLabel.CompareTo("oHagPt"))
        fitLabel                            = "mod. Hagedorn * #it{p}_{T}";
    else if (!fitLabel.CompareTo("h"))
        fitLabel                            = "Hagedorn";
    else if (!fitLabel.CompareTo("l"))
        fitLabel                            = "Tsallis";
    else if (!fitLabel.CompareTo("softHard"))
        fitLabel                            = "soft + hard";

    // label for graph legend
    TString label                           = "";
    if (isRatioPlot) label                  = graph1->GetYaxis()->GetTitle();
    else label                              = particle;
    label                                   = Form("%s %s", label.Data(), method.Data());

    // chi2 red of fit
    Double_t chi2red1, chi2red2, chi2red3;
    chi2red1                                = graph1Fit->GetChisquare() / graph1Fit->GetNDF();
    if (graph2Fit) chi2red2                 = graph2Fit->GetChisquare() / graph2Fit->GetNDF();
    else chi2red2                           = 0;
    if (graph3Fit) chi2red3                 = graph3Fit->GetChisquare() / graph3Fit->GetNDF();
    else chi2red3                           = 0;

    // create canvas and pads
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 2000);
    DrawCanvasSettings(canvas, 0.0, 0.0, 0.0, 0.0);
    TPad* pad                               = new TPad("pad", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawPadSettings(pad, 0.13, 0.015, 0.02, 0.);
    if(plotLogXaxisRatios) pad->SetLogx();
    pad->Draw();
    TPad* padRatioToFit                     = new TPad("padRatioToFit", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawPadSettings(padRatioToFit, 0.13, 0.015, 0.0, 0.25);
    if(plotLogXaxisRatios) padRatioToFit->SetLogx();
    padRatioToFit->Draw();

    if (xTitle.CompareTo("") == 0) xTitle   = graph1->GetXaxis()->GetTitle();
    if (yTitle.CompareTo("") == 0) yTitle   = graph1->GetYaxis()->GetTitle();

    TH1D* dummyHisto                        = new TH1D("dummyHisto", "", 100000, xMin, xMax);
    SetHistogramm(dummyHisto, xTitle, yTitle, 0.8*yMin, 1.2*yMax, 0.7, 1.82, 60, kFALSE, kTRUE, 43);
    dummyHisto->GetYaxis()->SetTickLength(0.02);
    dummyHisto->GetXaxis()->SetTickLength(0.02);

    // plot graphs + fits
    pad->cd();
    if (!isRatioPlot) pad->SetLogy();
    dummyHisto->Draw();


    SetGraph(graph1, xTitle, yTitle, yMin, yMax);
    DrawMarker(graph1, 20, 1.0, kBlack, kBlack);

    if (graph2 && secondIsSys){
        SetGraph(graph2, xTitle, yTitle, yMin, yMax);
        DrawMarker(graph2, 20, 1.0, kBlack, kBlack, 0);
    } else if (graph2){
        SetGraph(graph2, xTitle, yTitle, yMin, yMax);
        DrawMarker(graph2, 24, 1.0, kRed-4, kRed-4);
    }
    if (graph3) {
        SetGraph(graph2, xTitle, yTitle, yMin, yMax);
        DrawMarker(graph3, 24, 1.0, kBlue-4, kBlue-4);
    }

    if (graph1Fit) {
        graph1Fit->SetNpx(10000);
        graph1Fit->SetLineColor(kBlack);
        graph1Fit->SetLineWidth(1);
    }
    if (graph2Fit) {
        graph2Fit->SetNpx(10000);
        graph2Fit->SetLineColor(kRed-4);
        graph2Fit->SetLineWidth(1);
    }
    if (graph3Fit) {
        graph3Fit->SetNpx(10000);
        graph3Fit->SetLineColor(kBlue-4);
        graph3Fit->SetLineWidth(1);
    }

    graph1->Draw("p");
    if (secondIsSys && graph2) graph2->Draw("e2p");
    else if (graph2) graph2->Draw("p");
    if (graph3) graph3->Draw("p");

    if (graph1Fit) graph1Fit->Draw("same");
    if (graph2Fit) graph2Fit->Draw("same");
    if (graph3Fit) graph3Fit->Draw("same");

    // legend graph
    TLegend* legendGraph                    = NULL;
    if (isRatioPlot){
        if (((TString)graph1->GetYaxis()->GetTitle()).Contains("#eta/#pi^{0}"))
            legendGraph                     = GetAndSetLegend(0.45, 0.02, 0.45+0.25, 0.02+4*0.05, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()), 1, 0.14 ,60,43);
        else
            legendGraph                     = GetAndSetLegend(0.16, 0.95-4*0.05, 0.16+0.25, 0.95, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()), 1, 0.14 ,50,43);
    } else {
        legendGraph                         = GetAndSetLegend(0.16, 0.02, 0.16+0.25, 0.02+4*0.05, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()), 1, 0.14 ,60,43);
    }
    legendGraph->SetBorderSize(0);
    legendGraph->AddEntry(graph1, Form("%s", label.Data()), "p");
    if (graph2 && !secondIsSys) legendGraph->AddEntry(graph2, Form("%s %s", label.Data(), sysLabelA.Data()), "p");
    if (graph3 ) legendGraph->AddEntry(graph3, Form("%s %s", label.Data(), sysLabelB.Data()), "p");
    legendGraph->Draw();

    // legend fit
    TLegend* legendFit                      = NULL;
    if (isRatioPlot){
        if (label.Contains("#eta/#pi^{0}"))
            legendFit                           = GetAndSetLegend(0.16, 0.95-5*0.05, 0.16+0.25, 0.95, "", 1, 0.14, 60, 43);
        else
            legendFit                           = GetAndSetLegend(0.55, 0.95-5*0.05, 0.55+0.25, 0.95, "", 1, 0.14, 50, 43);
    } else {
        legendFit                               = GetAndSetLegend(0.45, 0.95-5*0.05, 0.45+0.25, 0.95, "", 1, 0.14, 60, 43);
    }
    legendFit->SetBorderSize(0);
    legendFit->AddEntry((TObject*)0, Form("%s", fitLabel.Data()), "");
    legendFit->AddEntry((TObject*)0, Form("#it{p}_{T} = %.2f - %.2f GeV/#it{c}", xMinFit, xMaxFit), "");
    if (!isThesis){
        legendFit->AddEntry(graph1Fit, Form("standard, #chi^{2}/ndf = %.2f", chi2red1), "l");
        if (graph2Fit) legendFit->AddEntry(graph2Fit, Form("%s, #chi^{2}/ndf = %.2f", sysLabelA.Data(), chi2red2), "l");
        if (graph3Fit) legendFit->AddEntry(graph3Fit, Form("%s, #chi^{2}/ndf = %.2f", sysLabelB.Data(), chi2red3), "l");
    } else {
        legendFit->AddEntry(graph1Fit, Form("standard"), "l");
        if (graph2Fit) legendFit->AddEntry(graph2Fit, Form("%s", sysLabelA.Data()), "l");
        if (graph3Fit) legendFit->AddEntry(graph3Fit, Form("%s", sysLabelB.Data()), "l");
    }
    legendFit->Draw();

    // plot ratio graph/fit
    padRatioToFit->cd();
    TH1D* dummyHistoRatio                        = new TH1D("dummyHistoRatio", "", 100000, xMin, xMax);
    SetHistogramm(dummyHistoRatio, "#it{p}_{T} (GeV/#it{c})","#frac{data}{fit}", yMinPlotRatio, yMaxPlotRatio, 3.4, 8.1, 60, kFALSE, kTRUE, 43, 510, 505);
    if(plotLogXaxisRatios>1) dummyHistoRatio->GetYaxis()->SetRangeUser(0.81,1.19);
    if (namePlot.Contains("ptConst") || namePlot.Contains("default") ) dummyHistoRatio->GetYaxis()->SetRangeUser(0.79,1.22);
    dummyHistoRatio->DrawCopy();

    DrawMarker(histoRatioToFit1, 20, 1.0, kBlack, kBlack);
    if (histoRatioToFit2 && secondIsSys)
        DrawMarker(histoRatioToFit2, 20, 1.5, kBlack, kBlack, 0);
    else if (histoRatioToFit2)
        DrawMarker(histoRatioToFit2, 24, 1.5, kRed-4, kRed-4);
    if (histoRatioToFit3) DrawMarker(histoRatioToFit3, 24, 1.0, kBlue-4, kBlue-4);

    histoRatioToFit1->Draw("same,e1");
    DrawLine(xMin, xMax, 1, 1, 1, kGray+2, 2);
    DrawLine(xMin, xMax, 1.1, 1.1, 1, kGray+2, 3);
    DrawLine(xMin, xMax, 0.9, 0.9, 1, kGray+2, 3);

    if (secondIsSys && histoRatioToFit2){
        histoRatioToFit2->Draw("e2,same");
        histoRatioToFit1->Draw("e1,same");
    } else if (histoRatioToFit2) {
        histoRatioToFit2->Draw("e1,same");
    }
    if (histoRatioToFit3) histoRatioToFit3->Draw("e1,same");

    canvas->SaveAs(Form("%s/%s.%s", outputDir.Data(), namePlot.Data(), suffix.Data()));

    graph1->GetXaxis()->UnZoom();
    graph1->GetYaxis()->UnZoom();
    histoRatioToFit1->GetXaxis()->UnZoom();

    // free pointer
    delete legendGraph;
    delete legendFit;
    delete dummyHisto;
    delete dummyHistoRatio;
    delete pad;
    delete padRatioToFit;
    delete canvas;
}

void ProducePlotWithRatioToFit( TH1D* histo1,
                                TH1D* histo2                = NULL,
                                TH1D* histo3                = NULL,
                                TH1D* histo4                = NULL,
                                TH1D* histo5                = NULL,
                                TF1* histo1Fit              = NULL,
                                TF1* histo2Fit              = NULL,
                                TF1* histo3Fit              = NULL,
                                TF1* histo4Fit              = NULL,
                                TF1* histo5Fit              = NULL,
                                TH1D* histoRatioToFit1      = NULL,
                                TH1D* histoRatioToFit2      = NULL,
                                TH1D* histoRatioToFit3      = NULL,
                                TH1D* histoRatioToFit4      = NULL,
                                TH1D* histoRatioToFit5      = NULL,
                                TString collSys             = "",
                                TString energy              = "",
                                TString cent                = "",
                                TString particle            = "",
                                TString method              = "",
                                TString fitLabel            = "",
                                TString suffix              = "",
                                TString xTitle              = "",
                                TString yTitle              = "",
                                TString namePlot            = "",
                                Double_t xMinFit            = 0.,
                                Double_t xMaxFit            = 15.,
                                Bool_t isRatioPlot          = kFALSE
                               ) {

    // create output dir
    TString outputDir;
    if (cent.CompareTo("") == 0) outputDir  = Form("plots/%s/%s_%s",suffix.Data(),collSys.Data(), energy.Data());
    else outputDir                          = Form("plots/%s/%s_%s_%s",suffix.Data(),collSys.Data(), energy.Data(), cent.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // get yRange
    TList* tempList                         = new TList;
    tempList->Add(histo1);
    tempList->Add(histo2);
    tempList->Add(histo3);
    tempList->Add(histo4);
    tempList->Add(histo5);
    Double_t yMin                           = GetYRangeExtremaFromList(tempList, kFALSE, kFALSE, collSys) * 0.2;
    Double_t yMax                           = GetYRangeExtremaFromList(tempList, kFALSE, kTRUE, collSys) * 2.5;

    if (isRatioPlot) {
        yMin                                = yMin * 0.9;
        yMax                                = yMax * 1.4;
    } else {
        yMin                                = yMin * 0.8;
        yMax                                = yMax * 1.2;
    }

    // get latex for energy, centrality and particle
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatex[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }
    for (Int_t i=0; i<nParticles; i++) {
        if (particle.CompareTo(fParticle[i]) == 0)
            particle                        = fParticleLatex[i];
    }

    // label for histo legend
    TString label                           = "";
    if (isRatioPlot) label                  = histo1->GetYaxis()->GetTitle();
    else label                              = particle;
    label                                   = Form("%s %s", label.Data(), method.Data());

    // pt range of fit
    //Double_t ptMin1, ptMax1, ptMin2, ptMax2, ptMin3, ptMax3;
    //histo1Fit->GetRange(ptMin1, ptMax1);
    //histo2Fit->GetRange(ptMin2, ptMax2);
    //histo3Fit->GetRange(ptMin3, ptMax3);

    // chi2 red of fit
    Double_t chi2red1, chi2red2, chi2red3, chi2red4, chi2red5;
    chi2red1                                = histo1Fit->GetChisquare() / histo1Fit->GetNDF();
    chi2red2                                = histo2Fit->GetChisquare() / histo2Fit->GetNDF();
    chi2red3                                = histo3Fit->GetChisquare() / histo3Fit->GetNDF();
    chi2red4                                = histo4Fit->GetChisquare() / histo4Fit->GetNDF();
    chi2red5                                = histo5Fit->GetChisquare() / histo5Fit->GetNDF();

    // create canvas and pads
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 2000);
    DrawCanvasSettings(canvas, 0.11, 0.015, 0.02, 0.09);
    TPad* pad                               = new TPad("pad", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawPadSettings(pad, 0.11, 0.015, 0.02, 0.);
    pad->Draw();
    TPad* padRatioToFit                     = new TPad("padRatioToFit", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawPadSettings(padRatioToFit, 0.11, 0.015, 0.0, 0.3);
    padRatioToFit->Draw();

    // plot histos + fits
    pad->cd();
    if (!isRatioPlot) pad->SetLogy();

    if (xTitle.CompareTo("") == 0) xTitle   = histo1->GetXaxis()->GetTitle();
    if (yTitle.CompareTo("") == 0) yTitle   = histo1->GetYaxis()->GetTitle();

    SetHistogramm(histo1, xTitle, yTitle, yMin, yMax,0.73, 1.15, 50, kFALSE, kTRUE, 43);
    if (histo2) SetHistogramm(histo2, xTitle, yTitle, yMin, yMax);
    if (histo3) SetHistogramm(histo3, xTitle, yTitle, yMin, yMax);
    if (histo4) SetHistogramm(histo4, xTitle, yTitle, yMin, yMax);
    if (histo5) SetHistogramm(histo5, xTitle, yTitle, yMin, yMax);

    DrawMarker(histo1, 20, 1.0, kBlack, kBlack);
    if (histo2) DrawMarker(histo2, 24, 1.0, kRed, kRed);
    if (histo3) DrawMarker(histo3, 24, 1.0, kBlue, kBlue);
    if (histo4) DrawMarker(histo4, 25, 1.0, kOrange-3, kOrange-3);
    if (histo5) DrawMarker(histo5, 25, 1.0, kAzure+1, kAzure+1);

    if (histo1Fit) {
        histo1Fit->SetLineColor(kBlack);
        histo1Fit->SetLineWidth(1);
    }
    if (histo2Fit) {
        histo2Fit->SetLineColor(kRed);
        histo2Fit->SetLineWidth(1);
    }
    if (histo3Fit) {
        histo3Fit->SetLineColor(kBlue);
        histo3Fit->SetLineWidth(1);
    }
    if (histo4Fit) {
        histo4Fit->SetLineColor(kOrange-3);
        histo4Fit->SetLineWidth(1);
    }
    if (histo5Fit) {
        histo5Fit->SetLineColor(kAzure+1);
        histo5Fit->SetLineWidth(1);
    }

    histo1->Draw("e1");
    if (histo2) histo2->Draw("e1,same");
    if (histo3) histo3->Draw("e1,same");
    if (histo4) histo4->Draw("e1,same");
    if (histo5) histo5->Draw("e1,same");

    if (histo1Fit) histo1Fit->Draw("same");
    if (histo2Fit) histo2Fit->Draw("same");
    if (histo3Fit) histo3Fit->Draw("same");
    if (histo4Fit) histo4Fit->Draw("same");
    if (histo5Fit) histo5Fit->Draw("same");

    // legend histo
    TLegend* legendHisto                    = NULL;
    if (isRatioPlot) legendHisto            = GetAndSetLegend(0.15, 0.7, 4, 1, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    else legendHisto                        = GetAndSetLegend(0.6, 0.7, 4, 1, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    legendHisto->AddEntry(histo1, Form("%s", label.Data()), "p");
    legendHisto->AddEntry(histo2, Form("%s upper sys", label.Data()), "p");
    legendHisto->AddEntry(histo3, Form("%s lower sys", label.Data()), "p");
    legendHisto->AddEntry(histo4, Form("%s slope up sys", label.Data()), "p");
    legendHisto->AddEntry(histo5, Form("%s slope low sys", label.Data()), "p");
    legendHisto->Draw();

    // legend fit
    TLegend* legendFit                      = NULL;
    if (isRatioPlot) legendFit              = GetAndSetLegend(0.6, 0.6, 7, 1);
    else legendFit                          = GetAndSetLegend(0.6, 0.35, 7, 1);
    legendFit->AddEntry((TObject*)0, Form("#it{p}_{T} = %.2f - %.2f GeV/#it{c}", xMinFit, xMaxFit), "");
    if(!isThesis){
        legendFit->AddEntry(histo1Fit, Form("%s standard, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red1), "l");
        legendFit->AddEntry(histo2Fit, Form("%s upper sys, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red2), "l");
        legendFit->AddEntry(histo3Fit, Form("%s lower sys, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red3), "l");
        legendFit->AddEntry(histo4Fit, Form("%s slope up sys, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red4), "l");
        legendFit->AddEntry(histo5Fit, Form("%s slope low sys, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red5), "l");
    } else {
        legendFit->AddEntry(histo1Fit, Form("%s standard", fitLabel.Data()), "l");
        legendFit->AddEntry(histo2Fit, Form("%s upper sys", fitLabel.Data()), "l");
        legendFit->AddEntry(histo3Fit, Form("%s lower sys", fitLabel.Data()), "l");
        legendFit->AddEntry(histo4Fit, Form("%s slope up sys", fitLabel.Data()), "l");
        legendFit->AddEntry(histo5Fit, Form("%s slope low sys", fitLabel.Data()), "l");
    }
    legendFit->Draw();

    // plot ratio histo/fit
    padRatioToFit->cd();
    SetStyleHistoTH1ForGraphs(histoRatioToFit1, "#it{p}_{T} (GeV/#it{c})","#frac{data}{fit}",  0.05, 0.045, 0.05, 0.045,  0.85, 0.4, 510, 505);
    histoRatioToFit1->GetYaxis()->SetRangeUser(0.5,1.5);

    DrawMarker(histoRatioToFit1, 20, 1.0, kBlack, kBlack);
    if (histoRatioToFit2) DrawMarker(histoRatioToFit2, 24, 1.0, kRed, kRed);
    if (histoRatioToFit3) DrawMarker(histoRatioToFit3, 24, 1.0, kBlue, kBlue);
    if (histoRatioToFit4) DrawMarker(histoRatioToFit4, 25, 1.0, kOrange-3, kOrange-3);
    if (histoRatioToFit5) DrawMarker(histoRatioToFit5, 25, 1.0, kAzure+1, kAzure+1);

    histoRatioToFit1->Draw("e1");
    if (histoRatioToFit2) histoRatioToFit2->Draw("e1,same");
    if (histoRatioToFit3) histoRatioToFit3->Draw("e1,same");
    if (histoRatioToFit4) histoRatioToFit4->Draw("e1,same");
    if (histoRatioToFit5) histoRatioToFit5->Draw("e1,same");

    DrawLine(histoRatioToFit1->GetXaxis()->GetXmin(), histoRatioToFit1->GetXaxis()->GetXmax(), 1, 1, 1, kGray+2, 2);
    DrawLine(histoRatioToFit1->GetXaxis()->GetXmin(), histoRatioToFit1->GetXaxis()->GetXmax(), 1.1, 1.1, 1, kGray+2, 3);
    DrawLine(histoRatioToFit1->GetXaxis()->GetXmin(), histoRatioToFit1->GetXaxis()->GetXmax(), 0.9, 0.9, 1, kGray+2, 3);

    canvas->SaveAs(Form("%s/%s.%s", outputDir.Data(), namePlot.Data(), suffix.Data()));

    // free pointer
    delete legendHisto;
    delete legendFit;
    delete pad;
    delete padRatioToFit;
    delete canvas;
}

void ProducePlotWithRatioToFit( TGraph* graph1,
                                TGraph* graph2               = NULL,
                                TGraph* graph3               = NULL,
                                TGraph* graph4               = NULL,
                                TGraph* graph5               = NULL,
                                TF1* graph1Fit               = NULL,
                                TF1* graph2Fit               = NULL,
                                TF1* graph3Fit               = NULL,
                                TF1* graph4Fit               = NULL,
                                TF1* graph5Fit               = NULL,
                                TH1D* histoRatioToFit1       = NULL,
                                TH1D* histoRatioToFit2       = NULL,
                                TH1D* histoRatioToFit3       = NULL,
                                TH1D* histoRatioToFit4       = NULL,
                                TH1D* histoRatioToFit5       = NULL,
                                TString collSys              = "",
                                TString energy               = "",
                                TString cent                 = "",
                                TString particle             = "",
                                TString method               = "",
                                TString fitLabel             = "",
                                TString suffix               = "",
                                TString xTitle               = "",
                                TString yTitle               = "",
                                TString namePlot             = "",
                                Double_t xMinFit             = 0.,
                                Double_t xMaxFit             = 15.,
                                Bool_t isRatioPlot           = kFALSE
                               ) {

    // create output dir
    TString outputDir;
    if (cent.CompareTo("") == 0) outputDir  = Form("plots/%s/%s_%s",suffix.Data(),collSys.Data(), energy.Data());
    else outputDir                          = Form("plots/%s/%s_%s_%s",suffix.Data(),collSys.Data(), energy.Data(), cent.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // get yRange
    TList* tempList                         = new TList;
    tempList->Add(graph1);
    tempList->Add(graph2);
    tempList->Add(graph3);
    tempList->Add(graph4);
    tempList->Add(graph5);

    // get xRange
    Double_t xMin                           = 0;
    Double_t xMax                           = GetXRangeExtremaFromList(tempList, kFALSE, kTRUE);

    // get yRange
    Double_t yMin                           = GetYRangeExtremaFromList(tempList, kFALSE, kFALSE, collSys) * 0.2;
    Double_t yMax                           = GetYRangeExtremaFromList(tempList, kFALSE, kTRUE, collSys) * 2.5;

    if (isRatioPlot) {
        yMin                                = yMin * 0.9;
        yMax                                = yMax * 1.4;
    } else {
        yMin                                = yMin * 0.8;
        yMax                                = yMax * 1.2;
    }

    // get latex for energy, centrality and particle
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }
    for (Int_t i=0; i<nParticles; i++) {
        if (particle.CompareTo(fParticle[i]) == 0)
            particle                        = fParticleLatex[i];
    }

    // label for graph legend
    TString label                           = "";
    if (isRatioPlot) label                  = graph1->GetYaxis()->GetTitle();
    else label                              = particle;
    label                                   = Form("%s %s", label.Data(), method.Data());

    // pt range of fit
    //Double_t ptMin1, ptMax1, ptMin2, ptMax2, ptMin3, ptMax3;
    //graph1Fit->GetRange(ptMin1, ptMax1);
    //graph2Fit->GetRange(ptMin2, ptMax2);
    //graph3Fit->GetRange(ptMin3, ptMax3);

    // chi2 red of fit
    Double_t chi2red1, chi2red2, chi2red3, chi2red4, chi2red5;
    chi2red1                                = graph1Fit->GetChisquare() / graph1Fit->GetNDF();
    chi2red2                                = graph2Fit->GetChisquare() / graph2Fit->GetNDF();
    chi2red3                                = graph3Fit->GetChisquare() / graph3Fit->GetNDF();
    chi2red4                                = graph4Fit->GetChisquare() / graph4Fit->GetNDF();
    chi2red5                                = graph5Fit->GetChisquare() / graph5Fit->GetNDF();

    // create canvas and pads
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 2000);
    DrawCanvasSettings(canvas, 0.11, 0.015, 0.02, 0.09);
    TPad* pad                               = new TPad("pad", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawPadSettings(pad, 0.11, 0.015, 0.02, 0.);
    pad->Draw();
    TPad* padRatioToFit                     = new TPad("padRatioToFit", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawPadSettings(padRatioToFit, 0.11, 0.015, 0.0, 0.3);
    padRatioToFit->Draw();

    // plot graphs + fits
    pad->cd();
    if (!isRatioPlot) pad->SetLogy();

    if (xTitle.CompareTo("") == 0) xTitle   = graph1->GetXaxis()->GetTitle();
    if (yTitle.CompareTo("") == 0) yTitle   = graph1->GetYaxis()->GetTitle();

    SetGraph(graph1, xTitle, yTitle, yMin, yMax);
    if (graph2) SetGraph(graph2, xTitle, yTitle, yMin, yMax);
    if (graph3) SetGraph(graph3, xTitle, yTitle, yMin, yMax);
    if (graph4) SetGraph(graph4, xTitle, yTitle, yMin, yMax);
    if (graph5) SetGraph(graph5, xTitle, yTitle, yMin, yMax);

    graph1->GetXaxis()->SetLimits(xMin, xMax);
    graph2->GetXaxis()->SetLimits(xMin, xMax);
    graph3->GetXaxis()->SetLimits(xMin, xMax);
    graph4->GetXaxis()->SetLimits(xMin, xMax);
    graph5->GetXaxis()->SetLimits(xMin, xMax);

    DrawMarker(graph1, 20, 1.0, kBlack, kBlack);
    if (graph2) DrawMarker(graph2, 24, 1.0, kRed, kRed);
    if (graph3) DrawMarker(graph3, 24, 1.0, kBlue, kBlue);
    if (graph4) DrawMarker(graph4, 25, 1.0, kOrange-3, kOrange-3);
    if (graph5) DrawMarker(graph5, 25, 1.0, kAzure+1, kAzure+1);

    if (graph1Fit) {
        graph1Fit->SetLineColor(kBlack);
        graph1Fit->SetLineWidth(1);
    }
    if (graph2Fit) {
        graph2Fit->SetLineColor(kRed);
        graph2Fit->SetLineWidth(1);
    }
    if (graph3Fit) {
        graph3Fit->SetLineColor(kBlue);
        graph3Fit->SetLineWidth(1);
    }
    if (graph4Fit) {
        graph4Fit->SetLineColor(kOrange-3);
        graph4Fit->SetLineWidth(1);
    }
    if (graph5Fit) {
        graph5Fit->SetLineColor(kAzure+1);
        graph5Fit->SetLineWidth(1);
    }

    graph1->Draw("pa");
    if (graph2) graph2->Draw("p");
    if (graph3) graph3->Draw("p");
    if (graph4) graph4->Draw("p");
    if (graph5) graph5->Draw("p");

    if (graph1Fit) graph1Fit->Draw("same");
    if (graph2Fit) graph2Fit->Draw("same");
    if (graph3Fit) graph3Fit->Draw("same");
    if (graph4Fit) graph4Fit->Draw("same");
    if (graph5Fit) graph5Fit->Draw("same");

    // legend graph
    TLegend* legendGraph                    = NULL;
    if (isRatioPlot) legendGraph            = GetAndSetLegend(0.15, 0.7, 4, 1, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    else legendGraph                        = GetAndSetLegend(0.6, 0.7, 4, 1, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    legendGraph->AddEntry(graph1, Form("%s", label.Data()), "p");
    legendGraph->AddEntry(graph2, Form("%s upper sys", label.Data()), "p");
    legendGraph->AddEntry(graph3, Form("%s lower sys", label.Data()), "p");
    legendGraph->AddEntry(graph4, Form("%s slope up sys", label.Data()), "p");
    legendGraph->AddEntry(graph5, Form("%s slope low sys", label.Data()), "p");
    legendGraph->Draw();

    // legend fit
    TLegend* legendFit                      = NULL;
    if (isRatioPlot) legendFit              = GetAndSetLegend(0.6, 0.6, 7, 1);
    else legendFit                          = GetAndSetLegend(0.6, 0.35, 7, 1);
    legendFit->SetMargin(0.13);
    legendFit->AddEntry((TObject*)0, Form("#it{p}_{T} = %.2f - %.2f GeV/#it{c}", xMinFit, xMaxFit), "");
    if (!isThesis){
        legendFit->AddEntry(graph1Fit, Form("%s standard, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red1), "l");
        legendFit->AddEntry(graph2Fit, Form("%s upper sys, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red2), "l");
        legendFit->AddEntry(graph3Fit, Form("%s lower sys, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red3), "l");
        legendFit->AddEntry(graph4Fit, Form("%s slope up sys, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red4), "l");
        legendFit->AddEntry(graph5Fit, Form("%s slope low sys, #chi^{2}/ndf = %.2f", fitLabel.Data(), chi2red5), "l");
    } else {
        legendFit->AddEntry(graph1Fit, Form("%s standard", fitLabel.Data()), "l");
        legendFit->AddEntry(graph2Fit, Form("%s upper sys", fitLabel.Data()), "l");
        legendFit->AddEntry(graph3Fit, Form("%s lower sys", fitLabel.Data()), "l");
        legendFit->AddEntry(graph4Fit, Form("%s slope up sys", fitLabel.Data()), "l");
        legendFit->AddEntry(graph5Fit, Form("%s slope low sys", fitLabel.Data()), "l");
    }
    legendFit->Draw();

    // plot ratio graph/fit
    padRatioToFit->cd();
    SetStyleHistoTH1ForGraphs(histoRatioToFit1, "#it{p}_{T} (GeV/#it{c})","#frac{data}{fit}", 0.14, 0.15, 0.12, 0.10,  0.85, 0.4, 510, 505);
    histoRatioToFit1->GetYaxis()->SetRangeUser(0.5,1.5);
    histoRatioToFit1->GetXaxis()->SetRangeUser(xMin,xMax);

    DrawMarker(histoRatioToFit1, 20, 1.0, kBlack, kBlack);
    if (histoRatioToFit2) DrawMarker(histoRatioToFit2, 24, 1.0, kRed, kRed);
    if (histoRatioToFit3) DrawMarker(histoRatioToFit3, 24, 1.0, kBlue, kBlue);
    if (histoRatioToFit4) DrawMarker(histoRatioToFit4, 25, 1.0, kOrange-3, kOrange-3);
    if (histoRatioToFit5) DrawMarker(histoRatioToFit5, 25, 1.0, kAzure+1, kAzure+1);

    histoRatioToFit1->Draw("e1");
    if (histoRatioToFit2) histoRatioToFit2->Draw("e1,same");
    if (histoRatioToFit3) histoRatioToFit3->Draw("e1,same");
    if (histoRatioToFit4) histoRatioToFit4->Draw("e1,same");
    if (histoRatioToFit5) histoRatioToFit5->Draw("e1,same");

    DrawLine(histoRatioToFit1->GetXaxis()->GetXmin(), histoRatioToFit1->GetXaxis()->GetXmax(), 1, 1, 1, kGray+2, 2);
    DrawLine(histoRatioToFit1->GetXaxis()->GetXmin(), histoRatioToFit1->GetXaxis()->GetXmax(), 1.1, 1.1, 1, kGray+2, 3);
    DrawLine(histoRatioToFit1->GetXaxis()->GetXmin(), histoRatioToFit1->GetXaxis()->GetXmax(), 0.9, 0.9, 1, kGray+2, 3);

    canvas->SaveAs(Form("%s/%s.%s", outputDir.Data(), namePlot.Data(), suffix.Data()));

    // free pointer
    delete legendGraph;
    delete legendFit;
    delete pad;
    delete padRatioToFit;
    delete canvas;
}

//================================================================================================================
//Function to produce plot comparing to fit
//================================================================================================================
void ProducePlotPtDepSystFactors(TF1*       funcLinA,
                                 TF1*       funcLinB,
                                 TF1*       funcPol2A,
                                 TF1*       funcPol2B,
                                 TString    collSys     = "",
                                 TString    energy      = "",
                                 TString    cent        = "",
                                 TString    particle    = "",
                                 TString    method      = "",
                                 TString    suffix      = "",
                                 TString    namePlot    = ""
) {

    // create output dir
    TString outputDir;
    if (cent.CompareTo("") == 0) outputDir  = Form("plots/%s/%s_%s",suffix.Data(),collSys.Data(), energy.Data());
    else outputDir                          = Form("plots/%s/%s_%s_%s",suffix.Data(),collSys.Data(), energy.Data(), cent.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // x range
    Double_t xMin                           = 0.;
    Double_t xMax                           = 20.;
    if (funcLinA) {
        xMin                                = funcLinA->GetXmin();
        xMax                                = funcLinA->GetXmax();
    }
    if (funcPol2A) {
        xMin                                = funcPol2A->GetXmin();
        xMax                                = funcPol2A->GetXmax();
    }

    // y range
    Double_t yMin                           = -1.5;
    Double_t yMax                           = 1.5;

    // get latex for energy, centrality and particle
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }
    for (Int_t i=0; i<nParticles; i++) {
        if (particle.CompareTo(fParticle[i]) == 0)
            particle                        = fParticleLatex[i];
    }

    // create canvas and pads
    TCanvas* canvas                         = GetAndSetCanvas("canvas");
    DrawCanvasSettings(canvas, 0.08, 0.015, 0.02, 0.09);

    // dummy histogram
    TH1D* dummyHisto                        = new TH1D("dummyHisto", "", 1000, xMin, xMax);
    SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","C(#it{p}_{T})", yMin, yMax);
    dummyHisto->GetYaxis()->SetTitleOffset(0.9);
    dummyHisto->Draw("x0");

    DrawLine(xMin,          xMax,           0,      0,      1, kGray+2, 2);
    DrawLine(xMin,          xMax,           1,      1,      1, kGray+2, 3);
    DrawLine(xMin,          xMax,           -1,     -1,     1, kGray+2, 3);
    DrawLine((xMax+xMin)/2, (xMax+xMin)/2,  -1.5,   1.5,    1, kGray+2, 3);

    // particle label
    TString particleLabel                   = "";
    if (method.CompareTo("") == 0)
        particleLabel                       = particle;
    else
        particleLabel                       = Form("%s (%s)", particle.Data(), method.Data());



    // legend
    TLegend* legendFit                      = NULL;
    if (funcLinA && funcLinB && funcPol2A && funcPol2B)
        legendFit                           = GetAndSetLegend(0.75, 0.90-(4*0.04), 0.95, 0.90, "", 1, 0.2 );
    else legendFit                          = GetAndSetLegend(0.75, 0.90-(2*0.04), 0.95, 0.90, "", 1, 0.2 );

    TLatex *labelEnergy = new TLatex(0.95, 0.915, Form("%s %s, %s %s", cent.Data(), collSys.Data(), energy.Data(), particleLabel.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);


    // lin factors
    if (funcLinA && funcLinB) {
        funcLinA->SetLineWidth(1);
        funcLinA->SetLineColor(kBlue+3);
        funcLinA->SetLineStyle(1);
        funcLinB->SetLineWidth(1);
        funcLinB->SetLineColor(kBlue-6);
        funcLinB->SetLineStyle(1);

        funcLinA->Draw("same");
        funcLinB->Draw("same");

        legendFit->AddEntry(funcLinA, "lin A", "l");
        legendFit->AddEntry(funcLinB, "lin B", "l");
    }

    // pol2 factors
    if (funcPol2A && funcPol2B) {
        funcPol2A->SetLineWidth(1);
        funcPol2A->SetLineColor(kGreen+3);
        funcPol2A->SetLineStyle(1);
        funcPol2B->SetLineWidth(1);
        funcPol2B->SetLineColor(kGreen-6);
        funcPol2B->SetLineStyle(1);

        funcPol2A->Draw("same");
        funcPol2B->Draw("same");

        legendFit->AddEntry(funcPol2A, "pol2 A", "l");
        legendFit->AddEntry(funcPol2B, "pol2 B", "l");
    }

    legendFit->Draw();
    labelEnergy->Draw();
    canvas->SaveAs(Form("%s/%s.%s", outputDir.Data(), namePlot.Data(), suffix.Data()));

    // free pointer
    delete legendFit;
    delete canvas;
    delete dummyHisto;
}


//================================================================================================================
//Function to produce plot comparing to fit
//================================================================================================================
void ProducePlotParamRatios(TF1*    funcCentral,
                            TF1*    funcConstUp,
                            TF1*    funcConstDown,
                            TF1*    funcLinA,
                            TF1*    funcLinB,
                            TF1*    funcPol2A,
                            TF1*    funcPol2B,
                            TString collSys     = "",
                            TString energy      = "",
                            TString cent        = "",
                            TString particle    = "",
                            TString method      = "",
                            TString suffix      = "",
                            TString namePlot    = ""
) {

    // create output dir
    TString outputDir;
    if (cent.CompareTo("") == 0) outputDir  = Form("plots/%s/%s_%s",suffix.Data(),collSys.Data(), energy.Data());
    else outputDir                          = Form("plots/%s/%s_%s_%s",suffix.Data(),collSys.Data(), energy.Data(), cent.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // x range
    Double_t xMin                           = 0.;
    Double_t xMax                           = 50.;

    // y range
    Double_t yMin                           = 0.5;
    Double_t yMax                           = 1.5;

    // get latex for energy, centrality and particle
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }
    for (Int_t i=0; i<nParticles; i++) {
        if (particle.CompareTo(fParticle[i]) == 0)
            particle                        = fParticleLatex[i];
    }

    // has param x
    Bool_t hasConst                         = kFALSE;
    Bool_t hasLin                           = kFALSE;
    Bool_t hasPol2                          = kFALSE;
    Int_t  hasTotal                         = 0;

    // calculate ratios to central value
    TF1* funcConstUpToCentral               = NULL;
    TF1* funcConstDownToCentral             = NULL;
    if (funcConstUp && funcConstDown) {
        hasConst                            = kTRUE;
        funcConstUpToCentral                = DivideTF1(funcConstUp,    funcCentral, "funcConstUpToCentral");
        funcConstDownToCentral              = DivideTF1(funcConstDown,  funcCentral, "funcConstDownToCentral");
        hasTotal++;
    }

    TF1* funcLinAToCentral                  = NULL;
    TF1* funcLinBToCentral                  = NULL;
    if (funcLinA && funcLinB) {
        hasLin                              = kTRUE;
        funcLinAToCentral                   = DivideTF1(funcLinA, funcCentral, "funcLinAToCentral");
        funcLinBToCentral                   = DivideTF1(funcLinB, funcCentral, "funcLinBToCentral");
        hasTotal++;
    }

    TF1* funcPol2AToCentral                 = NULL;
    TF1* funcPol2BToCentral                 = NULL;
    if (funcPol2A && funcPol2B) {
        hasPol2                             = kTRUE;
        funcPol2AToCentral                  = DivideTF1(funcPol2A, funcCentral, "funcPol2AToCentral");
        funcPol2BToCentral                  = DivideTF1(funcPol2B, funcCentral, "funcPol2BToCentral");
        hasTotal++;
    }

    // create canvas and pads
    TCanvas* canvas                         = GetAndSetCanvas("canvas");
    DrawCanvasSettings(canvas, 0.08, 0.015, 0.02, 0.09);

    // dummy histogram
    TH1D* dummyHisto                        = new TH1D("dummyHisto", "", 1000, xMin, xMax);
    SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","ratio", yMin, yMax);
    dummyHisto->GetYaxis()->SetTitleOffset(0.9);
    dummyHisto->Draw("x0");

    DrawLine(xMin, xMax, 1.0, 1.0, 1, kGray+2, 2);
    DrawLine(xMin, xMax, 1.1, 1.1, 1, kGray+2, 3);
    DrawLine(xMin, xMax, 0.9, 0.9, 1, kGray+2, 3);

    // particle label
    TString particleLabel                   = "";
    if (method.CompareTo("") == 0)
        particleLabel                       = particle;
    else
        particleLabel                       = Form("%s (%s)", particle.Data(), method.Data());

    // legend
    TLegend*                legendFit       = GetAndSetLegend(0.65, 0.90-((hasTotal*2-1)*0.04), 0.95, 0.90, "", 1, 0.2 );

    TLatex *labelEnergy = new TLatex(0.95, 0.915, Form("%s %s, %s %s", cent.Data(), collSys.Data(), energy.Data(), particleLabel.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);


    if (hasConst) {
        funcConstUpToCentral->SetLineWidth(1);
        funcConstUpToCentral->SetLineColor(kBlue+3);
        funcConstUpToCentral->SetLineStyle(1);
        funcConstDownToCentral->SetLineWidth(1);
        funcConstDownToCentral->SetLineColor(kBlue-6);
        funcConstDownToCentral->SetLineStyle(1);

        funcConstUpToCentral->Draw("same");
        funcConstDownToCentral->Draw("same");

        legendFit->AddEntry(funcConstUpToCentral, "const. up / central", "l");
        legendFit->AddEntry(funcConstDownToCentral, "const. down / central", "l");
    }

    if (hasLin) {
        funcLinAToCentral->SetLineWidth(1);
        funcLinAToCentral->SetLineColor(kRed+2);
        funcLinAToCentral->SetLineStyle(1);
        funcLinBToCentral->SetLineWidth(1);
        funcLinBToCentral->SetLineColor(kRed-4);
        funcLinBToCentral->SetLineStyle(1);

        funcLinAToCentral->Draw("same");
        funcLinBToCentral->Draw("same");

        legendFit->AddEntry(funcLinAToCentral, "lin A / central", "l");
        legendFit->AddEntry(funcLinBToCentral, "lin B / central", "l");
    }

    if (hasPol2) {
        funcPol2AToCentral->SetLineWidth(1);
        funcPol2AToCentral->SetLineColor(kGreen+3);
        funcPol2AToCentral->SetLineStyle(1);
        funcPol2BToCentral->SetLineWidth(1);
        funcPol2BToCentral->SetLineColor(kGreen-6);
        funcPol2BToCentral->SetLineStyle(1);

        funcPol2AToCentral->Draw("same");
        funcPol2BToCentral->Draw("same");

        legendFit->AddEntry(funcPol2AToCentral, "pol2 A / central", "l");
        legendFit->AddEntry(funcPol2BToCentral, "pol2 B / central", "l");
    }

    legendFit->Draw();
    labelEnergy->Draw();

    canvas->SaveAs(Form("%s/%s.%s", outputDir.Data(), namePlot.Data(), suffix.Data()));

    // free pointer
    delete legendFit;
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot with Mass dependent values
//================================================================================================================
void ProduceMassDependentPlot(  TGraphAsymmErrors** graphvsMassStat,
                                TGraphAsymmErrors** graphvsMassSys,
                                TGraphAsymmErrors** graphvsMassSysFunc,
                                TString collSys,
                                TString energy,
                                TString cent,
                                TString suffix,
                                TString outputDir,
                                TString quantityOut,
                                TString yAxisLabel,
                                Bool_t plotLogy                         = kFALSE,
                                Bool_t separateMesonsAndBaryons         = kFALSE
) {

    // create output dir
    gSystem->Exec("mkdir -p "+outputDir);
    cout << collSys.Data() << "\t"<< energy.Data() << "\t" << outputDir.Data() << endl;

    // get xRange
    Double_t xMin                           = 0;
    Double_t xMax                           = TMath::MaxElement(graphvsMassStat[0]->GetN(),graphvsMassStat[0]->GetX())*1.1;
    cout << "Plotting from: " << xMin << "\t" << xMax << endl;

    // get yRange
    Double_t yMin                           = TMath::MinElement(graphvsMassStat[0]->GetN(),graphvsMassStat[0]->GetY())*0.1;
    Double_t yMax                           = TMath::MaxElement(graphvsMassStat[0]->GetN(),graphvsMassStat[0]->GetY())*2;

    if (plotLogy){
        yMax                                = yMax*2.;
    }

    // create canvas
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 1800);
    DrawCanvasSettings(canvas, 0.10, 0.015, 0.015, 0.07);
    if (plotLogy)canvas->SetLogy();

    // get latex for enery and centrality
    Int_t energyIterator                    = -1;
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            energyIterator                  = i;
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    Int_t centIterator                      = -1;
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0){
            centIterator                    = i;
            cent                            = fCentralityLatex[i];
        }
    }

    TLatex *labelEnergy                     = new TLatex(0.94, 0.925, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);
    TString textALICE                       = "ALICE";
    if (isThesis){
        textALICE                           = "ALICE this thesis";
        cout << "thesis mode" << endl;
    }

    TLatex *labelALICE                     = new TLatex(0.94, 0.885, textALICE);
    SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);

    //     Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data())
    // dummy histogram
    TH1D* dummyHisto                        = NULL;

    dummyHisto                          = new TH1D("dummyHisto", "", 100000, 0, xMax*1.5);
    SetHistogramm(dummyHisto,"#it{M} (GeV/#it{c}^{2})",yAxisLabel, yMin, yMax, 0.73, 1.15, 0.04, kFALSE, kTRUE);
    dummyHisto->GetXaxis()->SetLabelOffset(-0.001);
    dummyHisto->GetXaxis()->SetRangeUser(xMin,xMax);
    dummyHisto->GetXaxis()->SetTickLength(0.02);
    dummyHisto->Draw();

    Int_t nErrors               =   1;

    if (!separateMesonsAndBaryons){
        if (graphvsMassSys[0]){
            DrawMarker(graphvsMassSys[0], 20, 2, kBlue+1, kBlue-6);
            gStyle->SetEndErrorSize(10);
            graphvsMassSys[0]->DrawClone("e||");
            nErrors++;
        }
        if (graphvsMassSysFunc[0]) {
            DrawMarker(graphvsMassSysFunc[0], 20, 2, kBlue+1, kBlack);
            gStyle->SetEndErrorSize(15);
            graphvsMassSysFunc[0]->DrawClone("e[]");
            nErrors++;
        }
        DrawMarker(graphvsMassStat[0], 20, 2, kBlue+1, kBlue+1);
        graphvsMassStat[0]->SetLineWidth(2);
        graphvsMassStat[0]->DrawClone("pZ");

        TLegend* legend                         = NULL;
        if (quantityOut.Contains("Yield"))
            legend                              = GetAndSetLegend(0.14, 0.10, 0.35, 0.10+(1.05*0.042*nErrors), "" ,1,0.27,0.04);
        else
            legend                              = GetAndSetLegend(0.14, 0.95-(1.05*0.04*nErrors), 0.35, 0.95, "" ,1,0.27,0.04);
        legend->AddEntry(graphvsMassStat[0], "stat. err.", "pe");
        if (graphvsMassSys[0])legend->AddEntry(graphvsMassSys[0], "sys. err.", "l");
        if (graphvsMassSysFunc[0])legend->AddEntry(graphvsMassSysFunc[0], "sys. err. functional form", "l");
        legend->Draw();
    } else {
        Int_t nColumns                          = 1;
        if (graphvsMassSys[1]){
            DrawMarker(graphvsMassSys[1], 20, 2, kBlue+1, kBlue-6);
            gStyle->SetEndErrorSize(10);
            graphvsMassSys[1]->DrawClone("e||");
            nErrors++;
        }
        if (graphvsMassSysFunc[1]) {
            DrawMarker(graphvsMassSysFunc[1], 20, 2, kBlue+1, kBlue+3);
            gStyle->SetEndErrorSize(15);
            graphvsMassSysFunc[1]->DrawClone("e[]");
            nErrors++;
        }
        DrawMarker(graphvsMassStat[1], 20, 2, kBlue+1, kBlue+1);
        graphvsMassStat[1]->SetLineWidth(2);
        graphvsMassStat[1]->DrawClone("pZ");

        if (graphvsMassSys[2]){
            DrawMarker(graphvsMassSys[2], 20, 2, kRed+1, kRed-6);
            gStyle->SetEndErrorSize(10);
            graphvsMassSys[2]->DrawClone("e||");
        }
        if (graphvsMassSysFunc[2]) {
            DrawMarker(graphvsMassSysFunc[2], 20, 2, kRed+1, kRed+3);
            gStyle->SetEndErrorSize(15);
            graphvsMassSysFunc[2]->DrawClone("e[]");
        }
        if (graphvsMassStat[2]){
            DrawMarker(graphvsMassStat[2], 20, 2, kRed+1, kRed+1);
            graphvsMassStat[2]->SetLineWidth(2);
            graphvsMassStat[2]->DrawClone("pZ");
            nColumns                            = 2;
        }
        TLegend* legend                         = NULL;
        TLegend* legend2                        = NULL;
        Double_t widthLegend                    = 0.21;
        Double_t yLegend                        = 0.95;
        Double_t xLegend                        = 0.14;
        Double_t xLegend2                       = 0.14;
        Double_t margin                         = 0.27;
        if (nColumns > 1){
            if (quantityOut.Contains("Yield"))
                yLegend                         = 0.10+(1.05*0.042*nErrors);
            else
                yLegend                         = 0.92;
            margin                              = 0.13;
            xLegend                             = 0.2;
            xLegend2                            = 0.38;
        }

        if (quantityOut.Contains("Yield")){
            legend                              = GetAndSetLegend(xLegend, 0.10, xLegend+widthLegend*nColumns, 0.10+(1.05*0.042*nErrors), "" ,nColumns,margin,0.04);
            legend2                              = GetAndSetLegend(xLegend2, 0.10, xLegend2+widthLegend, 0.10+(1.05*0.042*nErrors), "" ,1,0,0.04);
        } else {
            legend                              = GetAndSetLegend(xLegend, yLegend-(1.05*0.04*nErrors), xLegend+widthLegend*nColumns, yLegend, "" ,nColumns,margin,0.04);
            legend2                             = GetAndSetLegend(xLegend2, yLegend-(1.05*0.04*nErrors), xLegend2+widthLegend, yLegend, "" ,1,0,0.04);
        }
        if (nColumns == 1){
            legend->AddEntry(graphvsMassStat[1], "stat. err.", "pe");
            if (graphvsMassSys)legend->AddEntry(graphvsMassSys[1], "sys. err.", "l");
            if (graphvsMassSysFunc)legend->AddEntry(graphvsMassSysFunc[1], "sys. err. functional form", "l");
        } else {
            if(graphvsMassStat[1])legend->AddEntry(graphvsMassStat[1], "   ", "pe");
            if(graphvsMassStat[2])legend->AddEntry(graphvsMassStat[2], "        ", "pe");
            legend2->AddEntry((TObject*)0, "stat. err.", "");
            if(graphvsMassSys[1]){
                graphvsMassSys[1]->SetTitle("");
                legend->AddEntry(graphvsMassSys[1], "", "l");
            }
            if(graphvsMassSys[2]){
                graphvsMassSys[2]->SetTitle("");
                legend->AddEntry(graphvsMassSys[2], "", "l");
            }
            legend2->AddEntry((TObject*)0, "sys. err.", "");
            if (graphvsMassSysFunc[1])legend->AddEntry(graphvsMassSysFunc[1], " ", "l");
            if (graphvsMassSysFunc[2])legend->AddEntry(graphvsMassSysFunc[2], " ", "l");
            legend2->AddEntry((TObject*)0, "sys. err. func. form", "");

            legend2->Draw();
            TLatex *labelMeson                  = new TLatex(0.14, yLegend, "Mesons");
            SetStyleTLatex( labelMeson, 0.04,4, 1, 42, kTRUE, 11);
            labelMeson->Draw();
            TLatex *labelBaryon                  = new TLatex(0.28, yLegend, "Baryons");
            SetStyleTLatex( labelBaryon, 0.04,4, 1, 42, kTRUE, 11);
            labelBaryon->Draw();
        }
        legend->Draw();
    }
    dummyHisto->Draw("same,axis");

    // draw legend
    labelEnergy->Draw();
    labelALICE->Draw();

    // save canvas
    if (!separateMesonsAndBaryons)
        canvas->SaveAs(Form("%s/%sVsParticleMass.%s", outputDir.Data(),quantityOut.Data(), suffix.Data()));
    else
        canvas->SaveAs(Form("%s/%sVsParticleMass_SepMesonAndBaryon.%s", outputDir.Data(),quantityOut.Data(), suffix.Data()));
    // free pointer
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot with Mass dependent values
//================================================================================================================
void ProduceMassDependentPlotWithParticleOnAxis(    TGraphAsymmErrors** graphvsMassStat,
                                                    TGraphAsymmErrors** graphvsMassSys,
                                                    TGraphAsymmErrors** graphvsMassSysFunc,
                                                    TString collSys,
                                                    TString energy,
                                                    TString cent,
                                                    TString suffix,
                                                    TString outputDir,
                                                    TString quantityOut,
                                                    TString yAxisLabel,
                                                    Bool_t plotLogy                         = kFALSE,
                                                    Bool_t separateMesonsAndBaryons         = kFALSE
) {

    // create output dir
    gSystem->Exec("mkdir -p "+outputDir);
    cout << collSys.Data() << "\t"<< energy.Data() << "\t" << outputDir.Data() << endl;

    // get xRange
    Double_t xMin                           = 0;
    Double_t xMax                           = TMath::MaxElement(graphvsMassStat[0]->GetN(),graphvsMassStat[0]->GetX())*1.1;
    cout << "Plotting from: " << xMin << "\t" << xMax << endl;

    // get yRange
    Double_t yMin                           = TMath::MinElement(graphvsMassStat[0]->GetN(),graphvsMassStat[0]->GetY())*0.1;
    Double_t yMax                           = TMath::MaxElement(graphvsMassStat[0]->GetN(),graphvsMassStat[0]->GetY())*2;

    if (plotLogy){
        yMax                                = yMax*2.;
    }

    // create canvas
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 1800);
    DrawCanvasSettings(canvas, 0.10, 0.015, 0.015, 0.07);
    if (plotLogy)canvas->SetLogy();

    // get latex for enery and centrality
    Int_t energyIterator                    = -1;
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            energyIterator                  = i;
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    Int_t centIterator                      = -1;
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0){
            centIterator                    = i;
            cent                            = fCentralityLatex[i];
        }
    }

    TLatex *labelEnergy                     = new TLatex(0.94, 0.925, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);
    TString textALICE                       = "ALICE";
    if (isThesis) textALICE                 = "ALICE this thesis";

    TLatex *labelALICE                     = new TLatex(0.94, 0.885, textALICE);
    SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);

    //     Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data())
    // dummy histogram
    TH1D* dummyHisto                        = NULL;

    dummyHisto                          = new TH1D("dummyHisto", "", 100000, 0, xMax*1.5);
    SetHistogramm(dummyHisto,"#it{M} (GeV/#it{c}^{2})",yAxisLabel, yMin, yMax, 0.73, 1.15, 0.04, kFALSE, kTRUE);
    dummyHisto->GetXaxis()->SetLabelOffset(-0.001);
    dummyHisto->GetXaxis()->SetRangeUser(xMin,xMax);
    dummyHisto->GetXaxis()->SetTickLength(0.02);
    dummyHisto->Draw();

    Int_t nErrors               =   1;

    if (!separateMesonsAndBaryons){
        if (graphvsMassSys[0]){
            DrawMarker(graphvsMassSys[0], 20, 2, kBlue+1, kBlue-6);
            gStyle->SetEndErrorSize(10);
            graphvsMassSys[0]->DrawClone("e||");
            nErrors++;
        }
        if (graphvsMassSysFunc[0]) {
            DrawMarker(graphvsMassSysFunc[0], 20, 2, kBlue+1, kBlack);
            gStyle->SetEndErrorSize(15);
            graphvsMassSysFunc[0]->DrawClone("e[]");
            nErrors++;
        }
        DrawMarker(graphvsMassStat[0], 20, 2, kBlue+1, kBlue+1);
        graphvsMassStat[0]->SetLineWidth(2);
        graphvsMassStat[0]->DrawClone("pZ");

        TLegend* legend                         = NULL;
        if (quantityOut.Contains("Yield"))
            legend                              = GetAndSetLegend(0.14, 0.10, 0.35, 0.10+(1.05*0.042*nErrors), "" ,1,0.27,0.04);
        else
            legend                              = GetAndSetLegend(0.14, 0.95-(1.05*0.04*nErrors), 0.35, 0.95, "" ,1,0.27,0.04);
        legend->AddEntry(graphvsMassStat[0], "stat. err.", "pe");
        if (graphvsMassSys[0])legend->AddEntry(graphvsMassSys[0], "sys. err.", "l");
        if (graphvsMassSysFunc[0])legend->AddEntry(graphvsMassSysFunc[0], "sys. err. functional form", "l");
        legend->Draw();
    } else {
        Int_t nColumns                          = 1;
        if (graphvsMassSys[1]){
            DrawMarker(graphvsMassSys[1], 20, 2, kBlue+1, kBlue-6);
            gStyle->SetEndErrorSize(10);
            graphvsMassSys[1]->DrawClone("e||");
            nErrors++;
        }
        if (graphvsMassSysFunc[1]) {
            DrawMarker(graphvsMassSysFunc[1], 20, 2, kBlue+1, kBlue+3);
            gStyle->SetEndErrorSize(15);
            graphvsMassSysFunc[1]->DrawClone("e[]");
            nErrors++;
        }
        DrawMarker(graphvsMassStat[1], 20, 2, kBlue+1, kBlue+1);
        graphvsMassStat[1]->SetLineWidth(2);
        graphvsMassStat[1]->DrawClone("pZ");

        if (graphvsMassSys[2]){
            DrawMarker(graphvsMassSys[2], 20, 2, kRed+1, kRed-6);
            gStyle->SetEndErrorSize(10);
            graphvsMassSys[2]->DrawClone("e||");
        }
        if (graphvsMassSysFunc[2]) {
            DrawMarker(graphvsMassSysFunc[2], 20, 2, kRed+1, kRed+3);
            gStyle->SetEndErrorSize(15);
            graphvsMassSysFunc[2]->DrawClone("e[]");
        }
        if (graphvsMassStat[2]){
            DrawMarker(graphvsMassStat[2], 20, 2, kRed+1, kRed+1);
            graphvsMassStat[2]->SetLineWidth(2);
            graphvsMassStat[2]->DrawClone("pZ");
            nColumns                            = 2;
        }
        TLegend* legend                         = NULL;
        TLegend* legend2                        = NULL;
        Double_t widthLegend                    = 0.21;
        Double_t yLegend                        = 0.95;
        Double_t xLegend                        = 0.14;
        Double_t xLegend2                       = 0.14;
        Double_t margin                         = 0.27;
        if (nColumns > 1){
            if (quantityOut.Contains("Yield"))
                yLegend                         = 0.10+(1.05*0.042*nErrors);
            else
                yLegend                         = 0.92;
            margin                              = 0.13;
            xLegend                             = 0.2;
            xLegend2                            = 0.38;
        }

        if (quantityOut.Contains("Yield")){
            legend                              = GetAndSetLegend(xLegend, 0.10, xLegend+widthLegend*nColumns, 0.10+(1.05*0.042*nErrors), "" ,nColumns,margin,0.04);
            legend2                              = GetAndSetLegend(xLegend2, 0.10, xLegend2+widthLegend, 0.10+(1.05*0.042*nErrors), "" ,1,0,0.04);
        } else {
            legend                              = GetAndSetLegend(xLegend, yLegend-(1.05*0.04*nErrors), xLegend+widthLegend*nColumns, yLegend, "" ,nColumns,margin,0.04);
            legend2                             = GetAndSetLegend(xLegend2, yLegend-(1.05*0.04*nErrors), xLegend2+widthLegend, yLegend, "" ,1,0,0.04);
        }
        if (nColumns == 1){
            legend->AddEntry(graphvsMassStat[1], "stat. err.", "pe");
            if (graphvsMassSys)legend->AddEntry(graphvsMassSys[1], "sys. err.", "l");
            if (graphvsMassSysFunc)legend->AddEntry(graphvsMassSysFunc[1], "sys. err. functional form", "l");
        } else {
            if(graphvsMassStat[1])legend->AddEntry(graphvsMassStat[1], "   ", "pe");
            if(graphvsMassStat[2])legend->AddEntry(graphvsMassStat[2], "        ", "pe");
            legend2->AddEntry((TObject*)0, "stat. err.", "");
            if(graphvsMassSys[1]){
                graphvsMassSys[1]->SetTitle("");
                legend->AddEntry(graphvsMassSys[1], "", "l");
            }
            if(graphvsMassSys[2]){
                graphvsMassSys[2]->SetTitle("");
                legend->AddEntry(graphvsMassSys[2], "", "l");
            }
            legend2->AddEntry((TObject*)0, "sys. err.", "");
            if (graphvsMassSysFunc[1])legend->AddEntry(graphvsMassSysFunc[1], " ", "l");
            if (graphvsMassSysFunc[2])legend->AddEntry(graphvsMassSysFunc[2], " ", "l");
            legend2->AddEntry((TObject*)0, "sys. err. func. form", "");

            legend2->Draw();
            TLatex *labelMeson                  = new TLatex(0.14, yLegend, "Mesons");
            SetStyleTLatex( labelMeson, 0.04,4, 1, 42, kTRUE, 11);
            labelMeson->Draw();
            TLatex *labelBaryon                  = new TLatex(0.28, yLegend, "Baryons");
            SetStyleTLatex( labelBaryon, 0.04,4, 1, 42, kTRUE, 11);
            labelBaryon->Draw();
        }
        legend->Draw();
    }
    dummyHisto->Draw("same,axis");

    // draw legend
    labelEnergy->Draw();
    labelALICE->Draw();

    // save canvas
    if (!separateMesonsAndBaryons)
        canvas->SaveAs(Form("%s/%sVsParticle.%s", outputDir.Data(),quantityOut.Data(), suffix.Data()));
    else
        canvas->SaveAs(Form("%s/%sVsParticle_SepMesonAndBaryon.%s", outputDir.Data(),quantityOut.Data(), suffix.Data()));
    // free pointer
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot with cocktail parametrizations from 0--200 GeV/c
//================================================================================================================
void PlotCocktailParametrizationListFullRange(  TList* list,
                                                TList* inputlist,
                                                TString collSys,
                                                TString energy,
                                                TString cent,
                                                TString suffix,
                                                TString paramDirName  )
{
    // create output dir
    TString outputDir                       = "parametrizations";
    gSystem->Exec("mkdir -p "+outputDir);
    cout << collSys.Data() << "\t"<< energy.Data() << "\t" << outputDir.Data() << endl;

    // get xRange
    Double_t xMin                           = 0.01;
    Double_t xMax                           = 200;

    // get yRange
    //~ Double_t yMin                           = GetYRangeExtremaFromList(inputlist, kTRUE, kFALSE, collSys) * 0.5;
    Double_t yMin                           = 4e-14;
    Double_t yMax                           = GetYRangeExtremaFromList(inputlist, kTRUE, kTRUE, collSys) * 4;

    // create canvas
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0, 0, 1700, 2000);
    DrawCanvasSettings(canvas, 0.152, 0.015, 0.015, 0.068);
    canvas->SetLogy();
    canvas->SetLogx();

    Int_t fitCounter                        = 0;
    TString tempName                        = "";
    // count inputs
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                            = ((TObject*)list->At(i))->GetName();
        if ( tempName.Contains( "_pt" ) )
            fitCounter++;
        else
            continue;
    }
    cout << "counted: " << fitCounter << " inputs"<< endl;

    // get number of rows for legend
    Int_t nRows                             = fitCounter;
    // create legend
    TLegend* legend                         = NULL;
        legend                              = GetAndSetLegend(0.20, 0.09, 0.41, 0.09+(1.05*0.03*nRows), "" ,1,0.28,0.03);

    // set labels
    TLatex *labelEnergy                     = new TLatex(0.94, 0.925, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()));
    SetStyleTLatex( labelEnergy, 0.04,4, 1, 42, kTRUE, 31);
    TString textALICE                       = "ALICE";
    if (isThesis) textALICE                 = "ALICE this thesis";

    TLatex *labelALICE                     = new TLatex(0.94, 0.885, textALICE);
    SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);

    // dummy histogram
    TH1D* dummyHisto                        = NULL;
    dummyHisto                          = new TH1D("dummyHisto", "", 100000, xMin*0.8, xMax*1.5);
        SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", yMin, yMax, 0.7, 1.65, 0.04, kFALSE, kTRUE);
        dummyHisto->GetXaxis()->SetLabelOffset(-0.015);
    dummyHisto->GetXaxis()->SetTickLength(0.02);

    dummyHisto->Draw();

    // search list for particle spectra

    TF1* tempParam                                              = NULL;
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempParam                                               = (TF1*)list->At(i);
        DrawFit(tempParam,1,1.5,GetParticleColor(GetParticleNameFromPDG(tempParam->GetName())) );
        tempParam->Draw("same");
        cout << ((TString)tempParam->GetName()).Data() << "\t" << ((TString)GetParticleNameFromPDG(tempParam->GetName())).Data() << "\t" <<
                GetParticleIterator(GetParticleNameFromPDG(tempParam->GetName())) << "\t" << fParticleLatex[GetParticleIterator(GetParticleNameFromPDG(tempParam->GetName()))].Data()<< endl;
        legend->AddEntry(tempParam,Form("%s (%s)",fParticleLatex[GetParticleIterator(GetParticleNameFromPDG(tempParam->GetName()))].Data(),tempParam->GetName()),"l");
    }
    dummyHisto->Draw("same,axis");

    // draw legend
    legend->Draw();
    labelEnergy->Draw();
    labelALICE->Draw();

    // save canvas
    canvas->SaveAs(Form("%s/Parametrizations_%s_FullpT.%s", outputDir.Data(),paramDirName.Data(), suffix.Data()));

    // free pointer
    delete canvas;
    delete dummyHisto;
}

//================================================================================================================
//Function to produce plot containing all particle v2 in list
//================================================================================================================
void ProduceParticlev2PlotFromList(TList* list, TString collSys, TString energy, TString cent, TString suffix) {

    // create output dir
    TString outputDir                       = Form("plots/%s/%s",suffix.Data(),list->GetName());
    gSystem->Exec("mkdir -p "+outputDir);

    // get xRange
    Double_t xMin                           = 0;
    //     Double_t xMin                           = GetXRangeExtremaFromList(list, "spectra", kFALSE)*0.2;
    Double_t xMax                           = 7.99;
    //     Double_t xMax                           = GetXRangeExtremaFromList(list, "spectra", kTRUE)*1.5;
    //     if (xMin < 0 || xMin == 0)
    //       xMin                                  = 0.01;

    // get yRange
    Double_t yMin                           = 0.0;
    Double_t yMax                           = 0.39;

    // create canvas
    //     TCanvas* canvas                         = GetAndSetCanvas("canvas");
    TCanvas* canvas                         = GetAndSetCanvas("canvas", 0.11, 0.1, 800, 800);

    DrawCanvasSettings(canvas, 0.10, 0.015, 0.015, 0.068);
    //     canvas->SetLogy();
    //     canvas->SetLogx();

    // count number of stat histograms in list
    Int_t histCounter                       = 0;
    TString tempName                        = "";
    for (Int_t i=0; i<list->GetEntries(); i++) {
        tempName                            = ((TObject*)list->At(i))->GetName();
        if (tempName.Contains("v2") && tempName.Contains("Stat"))
            histCounter++;
        else
            continue;
    }

    // get latex for enery and centrality
    for (Int_t i=0; i<6; i++) {
        if (energy.CompareTo(fEnergy[i]) == 0) {
            if (cent.CompareTo("") == 0)
                energy                      = fEnergyLatexPP[i];
            else
                energy                      = fEnergyLatex[i];
        }
    }
    for (Int_t i=0; i<nCentralities; i++) {
        if (cent.CompareTo(fCentrality[i]) == 0)
            cent                            = fCentralityLatex[i];
    }

    // get number of rows for legend
    Int_t nRows                             = 0;
    if (histCounter%2 == 0) nRows           = histCounter/2 + 1;
    else nRows                              = (histCounter+1)/2 + 1;

    // create legend
    TLegend* legend                         = GetAndSetLegend(0.4, 0.97-(0.035*nRows), 0.95, 0.97, Form("%s %s, %s", cent.Data(), collSys.Data(), energy.Data()),2,0.12,0.035);

    // dummy histogram
    TH1D* dummyHisto                        = new TH1D("dummyHisto", "", 100000, xMin, xMax);
    SetHistogramm(dummyHisto,"#it{p}_{T}(GeV/#it{c})","v_{2}", yMin, yMax, 0.7, 1.65, 0.04, kFALSE, kTRUE);
    dummyHisto->GetXaxis()->SetTickLength(0.02);
    dummyHisto->GetXaxis()->SetLabelOffset(-0.001);
    dummyHisto->Draw();

    // search list for particle spectra
    TH1D* tempHist                          = NULL;
    TGraph* tempGraph                       = NULL;
    TGraphErrors* tempGraphErrs             = NULL;
    TGraphAsymmErrors* tempGraphAsymErrs    = NULL;
    TString tempClass                       = "";

    Int_t j                                 = 0;
    for (Int_t i=0; i<nParticles; i++) {
        if (list->FindObject(Form("v2_%sStat", fParticle[i].Data()))) {
            tempClass                   = ((TObject*)list->FindObject(Form("v2_%sStat", fParticle[i].Data())))->ClassName();

            if (tempClass.Contains("TH1")) {

                tempHist                = (TH1D*)((TH1D*)list->FindObject(Form("v2_%sStat", fParticle[i].Data())))->Clone("tempHist");
                DrawMarker(tempHist, markers[j], markerSize[j], paint[j],paint[j]);
                tempHist->Draw("e1,same");
                legend->AddEntry(tempHist, Form("%s", fParticleLatex[i].Data()), "lp");
                j++;
            } else if (tempClass.CompareTo("TGraph") == 0) {

                tempGraph               = (TGraph*)((TGraph*)list->FindObject(Form("v2_%sStat", fParticle[i].Data())))->Clone("tempGraph");
                DrawMarker(tempGraph, markers[j], markerSize[j], paint[j],paint[j]);
                tempGraph->DrawClone("p");
                legend->AddEntry(tempGraph, Form("%s", fParticleLatex[i].Data()), "lp");
                j++;
            } else if (tempClass.CompareTo("TGraphErrors") == 0) {

                tempGraphErrs           = (TGraphErrors*)((TGraphErrors*)list->FindObject(Form("v2_%sStat", fParticle[i].Data())))->Clone("tempGraphErrs");
                DrawMarker(tempGraphErrs, markers[j], markerSize[j], paint[j],paint[j]);
                tempGraphErrs->DrawClone("p");
                legend->AddEntry(tempGraphErrs, Form("%s", fParticleLatex[i].Data()), "lp");

                j++;
            } else if (tempClass.CompareTo("TGraphAsymmErrors") == 0) {

                tempGraphAsymErrs       = (TGraphAsymmErrors*)((TGraphAsymmErrors*)list->FindObject(Form("v2_%sStat", fParticle[i].Data())))->Clone(Form("v2_%sStat", fParticle[i].Data()));
                DrawMarker(tempGraphAsymErrs, markers[j], markerSize[j], paint[j],paint[j]);
                tempGraphAsymErrs->DrawClone("p");
                legend->AddEntry(tempGraphAsymErrs, Form("%s", fParticleLatex[i].Data()), "lp");

                j++;
                if (j>9) j=0;
            }
        }
    }

    // draw legend
    legend->Draw("same");

    // save canvas
    if (histCounter) canvas->SaveAs(Form("%s/Particlev2.%s", outputDir.Data(), suffix.Data()));

    // free pointer
    delete tempHist;
    delete legend;
    delete canvas;
    delete dummyHisto;
}
