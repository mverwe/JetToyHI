//Quick dirty plotting macro for jet response from trees of JetToyHI framework

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

using namespace std;

int markerColor[11] = {1,kRed+1,4,kGreen+3,kOrange+7,kGray+2,kRed+5,kAzure+10,kGreen+4,kYellow+2,kGreen};
int markerStyle[8] = {20,21,33,34,24,25,27,28};

TH1F *DrawFrame(double xmin = 0., double xmax = 1., double ymin = 0., double ymax = 1., TString xTitle = "x", TString yTitle = "y", bool setMargins = true);

void plotJetEnergyScaleCS(TString str = "JetToyHIResult.root") {

  TFile *f = new TFile(str.Data());
  TTree *tr = dynamic_cast<TTree*>(f->Get("jetTree"));

  tr->SetLineWidth(2);

  int nEvt = tr->GetEntriesFast();

  const int ncs = 5;
 
  TString strLeg[ncs] = {"CS #alpha=-2",
                        "CS #alpha=-1",
                        "CS #alpha=0",
                        "CS #alpha=1",
                        "CS #alpha=2",
  };


  TH1 *hResponsePt[ncs] = {0};
  TH1 *hResponseM[ncs]  = {0};

  //----------------------------------------------------------------
  // plot the jet pT response
  //----------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1","c1",450,400);
  gPad->SetLogy();
  tr->SetLineColor(1);
  //tr->Draw("csJetPt/sigJetPt>>hCSJetResponsePt(100,0.,2.)","sigJetPt>120. && abs(sigJetEta)<2.3");
  //hResponsePt[0] = dynamic_cast<TH1*>(gDirectory->Get("hCSJetResponsePt"));

  for(int i = 0; i<ncs; ++i) {
    tr->Draw(Form("csJet_%dPt/sigJetPt>>hCS%dJetResponsePt(100,0.,2.)",i,i),"sigJetPt>120. && abs(sigJetEta)<2.3","same");
    hResponsePt[i] = dynamic_cast<TH1*>(gDirectory->Get(Form("hCS%dJetResponsePt",i)));
  }

  for(int i = 0; i<ncs; ++i) {
    if(!hResponsePt[i]) continue;
    hResponsePt[i]->Sumw2();
    if(nEvt>0) hResponsePt[i]->Scale(1./(double)nEvt);
  }

  TH1F *fr1 = DrawFrame(0.,2.,5e-5,5e-1,"p_{T,sub}/p_{T,gen}","N/N_{evt}");
  TLegend *leg1 = new TLegend(0.23,0.6,0.45,0.93);
  leg1->SetFillColor(10);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.05);
  for(int i = 0; i<ncs; ++i) {
    if(!hResponsePt[i]) continue;
    hResponsePt[i]->SetLineColor(markerColor[i]);
    hResponsePt[i]->SetMarkerColor(markerColor[i]);
    hResponsePt[i]->SetMarkerStyle(markerStyle[i]);
    hResponsePt[i]->Draw("same");
    leg1->AddEntry(hResponsePt[i],strLeg[i].Data(),"p");
  }
  leg1->Draw();

  //----------------------------------------------------------------
  // plot the jet mass response
  //----------------------------------------------------------------
  TCanvas *c4 = new TCanvas("c4","c4",450,400);
  gPad->SetLogy();
  //tr->Draw("csJetM/sigJetM>>hCSJetResponseM(100,0.,5.)","sigJetPt>120. && abs(sigJetEta)<2.3");
  //hResponseM[0] = dynamic_cast<TH1*>(gDirectory->Get("hCSJetResponseM"));

  for(int i =0; i<ncs; ++i) {
    tr->Draw(Form("csJet_%dM/sigJetM>>hCS%dJetResponseM(100,0.,2.)",i,i),"sigJetPt>120. && abs(sigJetEta)<2.3","same");
    hResponseM[i] = dynamic_cast<TH1*>(gDirectory->Get(Form("hCS%dJetResponseM",i)));
  }

  for(int i = 0; i<ncs; ++i) {
    if(!hResponseM[i]) continue;
    hResponseM[i]->Sumw2();
    if(nEvt>0) hResponseM[i]->Scale(1./(double)nEvt);
  }
  
  TH1F *fr2 = DrawFrame(0.,5.,5e-5,5e-1,"M_{sub}/M_{gen}","N/N_{evt}");
  TLegend *leg2 = new TLegend(0.5,0.6,0.93,0.93);
  leg2->SetFillColor(10);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.05);
  for(int i = 0; i<ncs; ++i) {
    if(!hResponseM[i]) continue;
    hResponseM[i]->SetLineColor(markerColor[i]);
    hResponseM[i]->SetMarkerColor(markerColor[i]);
    hResponseM[i]->SetMarkerStyle(markerStyle[i]);
    hResponseM[i]->Draw("same");
    leg2->AddEntry(hResponseM[i],strLeg[i].Data(),"p");
  }
  leg2->Draw();
   
}

TH1F *DrawFrame(double xmin, double xmax, double ymin, double ymax, TString xTitle, TString yTitle, bool setMargins) {

  if(setMargins) {
    gPad->SetLeftMargin(0.22);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);//0.05);
    gPad->SetTopMargin(0.05);
  }

  TH1F *frame = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  frame->SetXTitle(xTitle.Data());
  frame->SetYTitle(yTitle.Data());
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->GetXaxis()->CenterTitle(true);
  frame->GetYaxis()->CenterTitle(true);

  gPad->SetTicks(1,1);

  return frame;
}
