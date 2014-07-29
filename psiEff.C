#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TString.h"
#include <iostream>
#define Pi 3.14159265
#define twoPi 6.28318530
#define EMASS 0.000511

#include "cuts.h"

Double_t InvPt(Double_t *x,Double_t *par)
{
  return par[0]*TMath::Power(TMath::Exp(par[1]*x[0])+x[0]/par[2],par[3]);
}

Double_t InvPt_test(Double_t *x,Double_t *par)
{
  return par[0]*TMath::Power(1+TMath::Power(x[0]*par[1],2),-par[2]);
}
Double_t sigmaFit(Double_t *x,Double_t *par)
{
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];

}
class StMyElectronEvent;
StMyElectronEvent *mElectronEvent;
class StMyElectron;
StMyElectron *mElectron;
StMyElectron *mElectron2;

TH1F *hMcVertexZ = new TH1F("mcVertexZ","mcVertexZ;Vz^{mc} (cm)",400,-200,200);
TH2F *hMcVertexXY = new TH2F("mcVertexXY","mcVertexXY;Vx^{mc} (cm);Vy^{mc} (cm)", 40, -2, 2, 40, -2, 2);
TH1F *hRcVertexZ = new TH1F("rcVertexZ","rcVertexZ;Vz^{rc} (cm)",400,-200,200);
TH1F *hVertexZdiff = new TH1F("vertexZdiff", "vertexZdiff;(Vz^{rc}-Vz^{mc} (cm)", 100, -5, 5);
TH1F *hRefMult = new TH1F ("refMult","refMult;Reference Multiplicity",1000,0,1000);
TH1F *hRefMultCut = new TH1F ("refMultCut","refMultCut;Reference Multiplicity after Cut",1000,0,1000);

TH1F *hNJpsi = new TH1F("hNJpsi","#Jpsi;#Jpsi;#Event", 10, 0, 10);
TH1D *hMcJpsiPt   = new TH1D("mcJpsiPt",  "mcJpsiPt;p_{T}^{mc} (GeV/c)",  300, 0, 30.);
TH1D *hMcJpsiPt_Or   = new TH1D("mcJpsiPt_Or",  "mcJpsiPt;p_{T}^{mc} (GeV/c)",  300, 0, 30.);

TH1F *hMcJpsiY  = new TH1F("mcJpsiY", "mcJpsiY; Y^{mc};", 100, -5,  5.);
TH2D *hMcJpsiPtY = new  TH2D("mcJpsiPtY","mcJpsiPtY;Y^{mc};p_{T}^{mc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH2D *hMcJpsiPtY_Or = new  TH2D("mcJpsiPtY_Or","mcJpsiPtY;Y^{mc};p_{T}^{mc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH1F *hMcJpsiPhi  = new TH1F("mcJpsiPhi", "mcJpsiPhi,#phi^{mc}", 360, -TMath::Pi(), TMath::Pi());
TH2D *hMcJpsiPtM = new  TH2D("mcJpsiPtM","mcJpsiPtM;M_{inv}^{mc}(e^{+}e^{-}) [GeV/c^{2}];p_{T}^{mc} (GeV/c)", 300, 0, 30, 1000,0,5);
TH2D *hMcJpsiPtM_1 = new  TH2D("mcJpsiPtM_1","mcJpsiPtM;M_{inv}^{mc}(e^{+}e^{-}) [GeV/c^{2}];p_{T}^{mc} (GeV/c)", 300, 0, 30, 1000,0,5);

TH1D *hMCElectronPt = new TH1D("mcElectronPt","input electron pt",300,0,30);
TH1D *hRecElectronPt_tpc = new TH1D("recElectronPt_tpc","Rec Tpc electron pt",300,0,30);
TH1D *hRecElectronPt_EMC = new TH1D("recElectronPt_emc","Rec EMC electron pt",300,0,30);

TH2D *hmcPtvsrcPt = new TH2D("mcPtvsrcPt","mcPt vs rcPt;p_{T}^{rc};p_{T}^{mc}",3500,0,35,4000,0,40);
TH2D *hmcPtvsrcPt_1 = new TH2D("mcPtvsrcPt_1","mcPt vs rcPt;p_{T}^{rc};p_{T}^{mc}",3500,0,35,4000,0,40);


TH1D *hRcJpsiPt   = new TH1D("rcJpsiPt",  "rcJpsiPt;p_{T}^{rc} (GeV/c)",  300, 0, 30.);

TH1D *hRcJpsiPt_Or   = new TH1D("rcJpsiPt_Or",  "rcJpsiPt;p_{T}^{rc} (GeV/c)",  300, 0, 30.);
TH1F *hRcJpsiY  = new TH1F("rcJpsiY", "rcJpsiY;y^{rc};", 100, -5,  5.);
TH2D *hRcJpsiPtY = new TH2D("rcJpsiPtY", "rcJpsiPtY;y^{rc};p_{T}^{rc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH2D *hRcJpsiPtY_Or = new TH2D("rcJpsiPtY_Or", "rcJpsiPtY;y^{rc};p_{T}^{rc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH1F *hRcJpsiPhi  = new TH1F("rcJpsiPhi", "rcJpsiPhi;#phi^{rc}", 360, -TMath::Pi(), TMath::Pi());
TH2D *hRcJpsiPtM = new TH2D("rcJpsiPtM", "rcJpsiPtM;M_{inv}^{rc}(e^{+}e^{-}) [GeV/c^{2}];p_{T}^{rc} (GeV/c)", 300, 0, 30, 1000,0,5);
TH2D *hRcJpsiPtDiff = new TH2D("rcJpsiPtDiff", "rcJpsiPtDiff;p_{T}^{mc}(GeV/c);p_{T}^{rc}/p_{T}^{mc}", 300, 0, 30,1000, 0, 2);
TH2D *hRcJpsiPtDiff_rc = new TH2D("rcJpsiPtDiff_rc", "rcJpsiPtDiff;p_{T}^{rc}(GeV/c);p_{T}^{rc}/p_{T}^{mc}", 300, 0, 30,1000, 0, 2);

TH3F *hRcJpsiYDiff = new TH3F("rcJpsiYDiff", "rcJpsiYDiff;p_{T}^{mc}(GeV/c);y^{mc};y^{rc}-y^{mc}", 250,0,25,40, -1, 1, 1000, -0.1, 0.1);
TH3F *hRcJpsiPhiDiff = new TH3F("rcJpsiPhiDiff", "rcJpsiPhiDiff;p_{T}^{mc}(GeV/c);y^{mc};#phi^{rc}-#phi^{mc}", 250,0,25,40, -1, 1, 1000, -0.1, 0.1);

TH2D *hJpsiMc = new TH2D("hJpsiMc","hJpsiMc;y^{mc};p_{T}^{mc} (GeV/c)",80,-2,2,300,0,30);

TH2D *hHt0JpsiTrg = new TH2D("hHt0JpsiTrg","hHt0JpsiTrg;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt0JpsiAdc0 = new TH2D("hHt0JpsiAdc0","hHt0JpsiAdc0;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt0JpsiPE = new TH2D("hHt0JpsiPE","hHt0JpsiPE;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt0JpsiNSMD = new TH2D("hHt0JpsiNSMD","hHt0JpsiNSMD;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt0JpsiDist = new TH2D("hHt0JpsiDist","hHt0JpsiDist;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);


TH2D *hHt1JpsiTrg = new TH2D("hHt1JpsiTrg","hHt1JpsiTrg;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt1JpsiAdc0 = new TH2D("hHt1JpsiAdc0","hHt1JpsiAdc0;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt1JpsiPE = new TH2D("hHt1JpsiPE","hHt1JpsiPE;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt1JpsiNSMD = new TH2D("hHt1JpsiNSMD","hHt1JpsiNSMD;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt1JpsiDist = new TH2D("hHt1JpsiDist","hHt1JpsiDist;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);

TH2D *hHt2JpsiTrg = new TH2D("hHt2JpsiTrg","hHt2JpsiTrg;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt2JpsiAdc0 = new TH2D("hHt2JpsiAdc0","hHt2JpsiAdc0;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt2JpsiPE = new TH2D("hHt2JpsiPE","hHt2JpsiPE;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt2JpsiPE_Or = new TH2D("hHt2JpsiPE_Or","hHt2JpsiPE;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);

TH2D *hHt2JpsiNSMD = new TH2D("hHt2JpsiNSMD","hHt2JpsiNSMD;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt2JpsiDist = new TH2D("hHt2JpsiDist","hHt2JpsiDist;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);

TH2F *hHt3JpsiTrg = new TH2F("hHt3JpsiTrg","hHt3JpsiTrg;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2F *hHt3JpsiAdc0 = new TH2F("hHt3JpsiAdc0","hHt3JpsiAdc0;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2F *hHt3JpsiPE = new TH2F("hHt3JpsiPE","hHt3JpsiPE;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2F *hHt3JpsiNSMD = new TH2F("hHt3JpsiNSMD","hHt3JpsiNSMD;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2F *hHt3JpsiDist = new TH2F("hHt3JpsiDist","hHt3JpsiDist;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);

TH3F *hHt0JpsiMassPE = new TH3F("hHt0JpsiMassPE","hHt0JpsiMassPE;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt0JpsiMassDist = new TH3F("hHt0JpsiMassDist","hHt0JpsiMassDist;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt1JpsiMassPE = new TH3F("hHt1JpsiMassPE","hHt1JpsiMassPE;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt1JpsiMassDist = new TH3F("hHt1JpsiMassDist","hHt1JpsiMassDist;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt2JpsiMassPE_1 = new TH3F("hHt2JpsiMassPE_1","hHt2JpsiMassPE;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt2JpsiMassPE = new TH3F("hHt2JpsiMassPE","hHt2JpsiMassPE;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt2JpsiMassDist = new TH3F("hHt2JpsiMassDist","hHt2JpsiMassDist;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt3JpsiMassPE = new TH3F("hHt3JpsiMassPE","hHt3JpsiMassPE;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt3JpsiMassDist = new TH3F("hHt3JpsiMassDist","hHt3JpsiMassDist;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);


TH3F *hJpsi3DMc = new TH3F("hJpsi3DMc","hJpsi3DMc;y^{mc};Vz(cm);p_{T}^{mc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);

TH3F *hHt0Jpsi3DAdc0 = new TH3F("hHt0Jpsi3DAdc0","hHt0Jpsi3DAdc0;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt0Jpsi3DPE = new TH3F("hHt0Jpsi3DPE","hHt0Jpsi3DPE;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt0Jpsi3DNSMD = new TH3F("hHt0Jpsi3DNSMD","hHt0Jpsi3DNSMD;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt0Jpsi3DDist = new TH3F("hHt0Jpsi3DDist","hHt0Jpsi3DDist;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);

TH3F *hHt1Jpsi3DAdc0 = new TH3F("hHt1Jpsi3DAdc0","hHt1Jpsi3DAdc0;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt1Jpsi3DPE = new TH3F("hHt1Jpsi3DPE","hHt1Jpsi3DPE;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt1Jpsi3DNSMD = new TH3F("hHt1Jpsi3DNSMD","hHt1Jpsi3DNSMD;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt1Jpsi3DDist = new TH3F("hHt1Jpsi3DDist","hHt1Jpsi3DDist;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);

TH3F *hHt2Jpsi3DAdc0 = new TH3F("hHt2Jpsi3DAdc0","hHt2Jpsi3DAdc0;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt2Jpsi3DPE = new TH3F("hHt2Jpsi3DPE","hHt2Jpsi3DPE;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt2Jpsi3DNSMD = new TH3F("hHt2Jpsi3DNSMD","hHt2Jpsi3DNSMD;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt2Jpsi3DDist = new TH3F("hHt2Jpsi3DDist","hHt2Jpsi3DDist;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);

TH3F *hHt3Jpsi3DAdc0 = new TH3F("hHt3Jpsi3DAdc0","hHt3Jpsi3DAdc0;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt3Jpsi3DPE = new TH3F("hHt3Jpsi3DPE","hHt3Jpsi3DPE;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt3Jpsi3DNSMD = new TH3F("hHt3Jpsi3DNSMD","hHt3Jpsi3DNSMD;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt3Jpsi3DDist = new TH3F("hHt3Jpsi3DDist","hHt3Jpsi3DDist;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);

void psiEff(const char* fileInDir="./out_Jpsi_200", const char* OutFile="jpsiEff_200", const int mCentCut1=0, const int mCentCut2=9){
 // cout<<"I"<<endl;
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StMyElectronMaker");
  const Int_t nFilesMax = 10000;
  //void SetDefaultSumw2(Bool_t sumw2=kTRUE);
 // TH1::SetDefaultSumw2(Bool_t sumw2=kTRUE);
//  TH2::SetDefaultSumw2(Bool_t sumw2=kTRUE);
//  TH3::SetDefaultSumw2(Bool_t sumw2=kTRUE);

  TF1 *myGaus = new TF1("myGaus","gaus",-10,10);
  myGaus->SetParameters(1,0,0.9);
  TF1 *myGaus_1 = new TF1("myGaus_1","gaus",-10,10);
  myGaus_1->SetParameters(1,0,0.9);

  mElectronEvent = new StMyElectronEvent();
  TChain chMc("mcT");
  chMc.SetBranchAddress("mcE",&mElectronEvent);

  const Double_t mDsmAdcCut[4] = {11, 15, 18, 25};
  const float pi = 3.1415926;
  void* dir = gSystem->OpenDirectory(gSystem->ExpandPathName(fileInDir));
  int nruns=0;
  char *file_name;
  TString Tname;
  char file_list[nFilesMax][256];
  do {
    file_name = (char*)gSystem->GetDirEntry(dir);
    Tname=file_name;
      if(file_name && Tname.Contains("myminimc.root")&&Tname.Contains("emb")) {

      sprintf(file_list[nruns],"%s/%s",fileInDir,file_name);
      chMc.Add(file_list[nruns]);
      cout << " read in " << file_list[nruns] << endl;
      nruns++;
    }
  } while (file_name && nruns<=nFilesMax);

//  ch.Add(fileIn);
  int iret = 0;
  int nb=0;
  cout << chMc.GetEntries() << " events in chain" << endl;
  int nevents = chMc.GetEntries();
  //nevents/=200;

  //+--------------------------------+
  //| initialize wegiht function     |
  //+--------------------------------+
  TF1 *function_W = new TF1("function_W",InvPt,0,30,4);
  function_W->SetParameters(0.4,-0.4042,4.053,-6);
//  function_W->SetParameters(0.4,-0.5905,4.034,-7.255);
  //+---------------------------------+
  //|   loop the events               |
  //+---------------------------------+
  TF1 *function_sigma = new TF1("function_sigma",sigmaFit,0,30,3);
  function_sigma->SetParameters(1.77250e-2,3.18836e-3,1.68829e-3);


  for (int i = 0; i < nevents; i++) {
    nb +=chMc.GetEvent(i);
    chMc.GetEvent(i);
    if(i%1000==0) cout<<"event# "<<i<<endl;
    int n;
    if (!mElectronEvent) continue;
    if (mElectronEvent->eventID()<=0) continue;

    hMcVertexZ->Fill(mElectronEvent->mcVertexZ());
    hMcVertexXY->Fill(mElectronEvent->mcVertexX(), mElectronEvent->mcVertexY());
    hRcVertexZ->Fill(mElectronEvent->vertexZ());
    hVertexZdiff->Fill(mElectronEvent->vertexZ() - mElectronEvent->mcVertexZ());

    Double_t vz = mElectronEvent->vertexZ();

    int Mult = mElectronEvent->refMultPos() + mElectronEvent->refMultNeg();
    hRefMult->Fill(Mult);
    //cout<<"vertexZ = "<<mElectronEvent->vertexZ()<<endl;
    int cent = getCentrality(Mult);
    if(cent<mCentCut1||cent>mCentCut2) continue;
    hRefMultCut->Fill(Mult);

    TLorentzVector jpsiMc(0.,0.,0.,0.), ePosMc(0.,0.,0.,0.), eNegMc(0.,0.,0.,0.);
    TLorentzVector jpsiRc(0.,0.,0.,0.), ePosRc(0.,0.,0.,0.), eNegRc(0.,0.,0.,0.);

    Int_t nJpsi = 0;
    for(int j=0;j<mElectronEvent->nReal();j++){
      mElectron = (StMyElectron) mElectronEvent->real()->UncheckedAt(j);
      //cout<<j<<"  "<<mElectron->geantId<<"  "<<mElectron->pGeantId<<endl;
      if(mElectron->pGeantId!=160) continue;//not from Jpsi electron
      //if(mElectron->pGeantId>0) continue;//not original electron
      if(mElectron->mcId<0) continue;
    //  if(mElectron->tpcCommonHits<10)continue;

      bool tag = kFALSE;
      for(int k=0;k<mElectronEvent->nReal();k++) {
	mElectron2 = (StMyElectron) mElectronEvent->real()->UncheckedAt(k);
	//cout<<k<<"  "<<mElectron2->geantId<<"  "<<mElectron2->pGeantId<<endl;
	if(mElectron2->pGeantId!=160) continue;
	if(mElectron2->mcId<0) continue;
	if(mElectron2->mcId==mElectron->mcId) continue;
  //     if(mElectron2->tpcCommonHits<10)continue;
	Double_t deta = mElectron->mcY - mElectron2->mcY;
	Double_t dphi = mElectron->mcPhi - mElectron2->mcPhi;
	while(dphi>2*TMath::Pi()) dphi -= 2.*TMath::Pi();
	while(dphi<0) dphi += 2.*TMath::Pi();
	while(dphi>TMath::Pi()) dphi = dphi - 2*TMath::Pi();
	if(TMath::Abs(deta)<0.1&&TMath::Abs(dphi)<0.5&&mElectron2->pId!=mElectron->pId) tag = kTRUE;
	//cout<<deta<<"  "<<dphi<<"  "<<tag<<endl;
      }
  //  cout<<tag<<endl;
     if(tag) continue;
	Double_t pt_tem = mElectron->mcPt;
	hMCElectronPt->Fill(pt_tem);
	Int_t count1 =0;
      for(int k=j+1; k<mElectronEvent->nReal();k++) {
	mElectron2 =  (StMyElectron) mElectronEvent->real()->UncheckedAt(k);
	//cout<<"xxx"<<k<<"  "<<mElectron2->geantId<<"  "<<mElectron2->pGeantId<<endl;
	if(mElectron2->pGeantId!=160) continue;
	if(mElectron2->mcId<0) continue;
//	if(mElectron2->tpcCommonHits<10)continue;
	if(mElectron2->mcId==mElectron->mcId) continue;
	if(mElectron2->pId!=mElectron->pId) continue;//from same parent
	if(mElectron->geantId==2&&mElectron2->geantId==3) {
	  ePosMc.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
	  eNegMc.SetPtEtaPhiM(mElectron2->mcPt, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
	  ePosRc.SetPtEtaPhiM(mElectron->pt, mElectron->eta, mElectron->phi, EMASS);
	  eNegRc.SetPtEtaPhiM(mElectron2->pt, mElectron2->eta, mElectron2->phi, EMASS);
	} else if(mElectron->geantId==3&&mElectron2->geantId==2) {
	  eNegMc.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
	  ePosMc.SetPtEtaPhiM(mElectron2->mcPt, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
	  eNegRc.SetPtEtaPhiM(mElectron->pt, mElectron->eta, mElectron->phi, EMASS);
	  ePosRc.SetPtEtaPhiM(mElectron2->pt, mElectron2->eta, mElectron2->phi, EMASS);
	} else {
	  cout<<"!!!! something wrong with the J/psi decay !!!!!!"<<endl;
	  continue;
	}
	nJpsi++;
	jpsiMc = ePosMc + eNegMc;
	//cout<<jpsiMc.M()<<endl;


      Double_t weight1 = 1; 
     // Double_t weight1 = function_W->Eval(jpsiMc.Pt());
     // Double_t weight = function_W->Eval(jpsiRc.Pt());
	if(mElectron->mcId>=0 && mElectron2->mcId>=0) {
      hMcJpsiPt_Or->Fill(jpsiMc.Pt());
     // hweightPt->Fill(jpsiMc.Pt(),weight1);
	  hMcJpsiPt->Fill(jpsiMc.Pt(),weight1);
	  hMcJpsiY->Fill(jpsiMc.Rapidity());
	  hMcJpsiPtY->Fill(jpsiMc.Rapidity(), jpsiMc.Pt(),weight1);
	  hMcJpsiPtY_Or->Fill(jpsiMc.Rapidity(), jpsiMc.Pt());
	  hMcJpsiPhi->Fill(jpsiMc.Phi());
	  hMcJpsiPtM->Fill(jpsiMc.Pt(), jpsiMc.M());
	  hMcJpsiPtM_1->Fill(jpsiMc.Pt(), jpsiMc.M(),weight1);

	} else {
	  cout<<" !!! something wrong with MC tracks !!!!!!"<<endl;
	  continue;
	}//all of the MC J/psi below are good!!!!!!!!



     if(mElectron->tpcCommonHits<10)continue;
     if(mElectron2->tpcCommonHits<10)continue;

	if(mElectron->id>=0 && mElectron2->id>=0) jpsiRc = ePosRc + eNegRc;
     // Double_t weight1 = function_W->Eval(jpsiMc.Pt());
    //  Double_t weight = function_W->Eval(jpsiRc.Pt());
    //////////////////////////////////////////
  // Double_t weight1 = 1.;
  // Double_t weight = 1.;
    Double_t mcPt_sigma = function_sigma->Eval(jpsiMc.Pt());
    //////////////////////////////////////////
    hmcPtvsrcPt->Fill(jpsiRc.Pt(),jpsiMc.Pt());


	if(mElectron->id>=0&&mElectron2->id>=0)
       	{
	      if(jpsiRc.Pt()>(jpsiMc.Pt()*(1.+3*mcPt_sigma)))continue;
	}
     hmcPtvsrcPt_1->Fill(jpsiRc.Pt(),jpsiMc.Pt());
	if(mElectron->id>=0&&mElectron2->id>=0){

	  //cout<<mElectron->id<<endl;
	  Double_t rcPt = jpsiRc.Pt();
	  Double_t rcY = jpsiRc.Rapidity();
	  Double_t rcPhi = jpsiRc.Phi();
	  hRcJpsiPt_Or->Fill(rcPt);
	 // hweightPt_1->Fill(rcPt,weight1);
	  hRcJpsiPt->Fill(rcPt,weight1);

	  hRcJpsiY->Fill(rcY);
	  hRcJpsiPtY->Fill(rcY, rcPt,weight1);
	  hRcJpsiPtY_Or->Fill(rcY, rcPt);
	  hRcJpsiPhi->Fill(rcPhi,weight1);
	  hRcJpsiPtM->Fill(rcPt, jpsiRc.M(),weight1);

       hRcJpsiPtDiff_rc->Fill(rcPt,rcPt/jpsiMc.Pt());

	  if(fabs(rcY)<1){
	    Double_t mcPt = jpsiMc.Pt();
	    Double_t mcY = jpsiMc.Rapidity();
	    Double_t mcPhi = jpsiMc.Phi();
	    hRcJpsiPtDiff->Fill(mcPt, rcPt/mcPt);

	    hRcJpsiYDiff->Fill(mcPt, mcY, rcY-mcY);
	    hRcJpsiPhiDiff->Fill(mcPt, mcY, rcPhi-mcPhi);

	  }//end of |eta|<1 cut
	}//end of RC electron loop

	hJpsiMc->Fill(jpsiMc.Rapidity(), jpsiMc.Pt());
	hJpsi3DMc->Fill(jpsiMc.Rapidity(), vz, jpsiMc.Pt());

	if(mElectron->id>=0&&mElectron2->id>=0) {//good reconstructed electrons
	  Double_t eta1 = mElectron->eta;
	  Double_t dca1 = mElectron->dca;
	  Double_t nsigma1 = myGaus_1->GetRandom();
	  Double_t nHitsFit1 = mElectron->nFitPts;
	  Double_t e1 = mElectron->energy;
	  Double_t adc01 = mElectron->adc0;
	  Double_t dsmAdc01 = mElectron->dsmAdc0;
	  Double_t p1 = mElectron->p;
	  Double_t pe1 = (e1>0)? p1/e1:9999;
	  Double_t pt1 = mElectron->pt;
	  Double_t nEta1 = mElectron->nEta;
	  Double_t nPhi1 = mElectron->nPhi;
	  Double_t zDist1 = mElectron->zDist;
	  Double_t phiDist1 = mElectron->phiDist;
	  bool isTPC1[4], isTrg1[4], isAdc01[4], isPE1[4], isNSMD1[4], isDist1[4];
	  for(int iht=0;iht<4;iht++) {
	    isTPC1[iht] = kFALSE;
	    isTrg1[iht] = kFALSE;
	    isAdc01[iht] = kFALSE;
	    isPE1[iht] = kFALSE;
	    isNSMD1[iht] = kFALSE;
	    isDist1[iht] = kFALSE;
	  }
	  //TPC track
	  Int_t count2=0; 
	  if(eta1>=mTpceEtaCut[0]&&eta1<=mTpceEtaCut[1] &&
	      nHitsFit1>=mTpceHitsFitCut &&
	      dca1<=mTpceDcaCut &&
	      nsigma1>mTpcenSigmaElectronCut[0]&&nsigma1<mTpcenSigmaElectronCut[1] &&
	      pt1<30.&&p1<30.) {
	    for(int iht=0;iht<4;iht++) {
	      if(pt1>=mTpcePtCut[iht]&&p1>=mTpcePCut[iht]) isTPC1[iht] = kTRUE;
	      if(count1 ==0&&count2==0)
	      {
	     	 hRecElectronPt_tpc->Fill(mElectron->pt);
		 count2++;
	      }
	    }
	  }
	  //EMC Track
	  count2=0;
	  if(nHitsFit1>=mTpceHitsFitCut &&
	      dca1<=mEmceDcaCut &&
	      eta1>=mEmceEtaCut[0]&&eta1<=mEmceEtaCut[1] &&
	      nsigma1>mEmcenSigmaElectronCut[0]&&nsigma1<mEmcenSigmaElectronCut[1] &&
	      pt1<30) {
	    for(int iht=0;iht<4;iht++) {
	      if(pt1>mEmcePtCut[iht]&&e1>0&&dsmAdc01>mDsmAdcCut[iht]) {
		isTrg1[iht] = kTRUE;
		if(adc01>mEmceAdcCut[iht]) {
		  isAdc01[iht] = kTRUE;
		  if(pe1>mEmcePECut[0]&&pe1<mEmcePECut[1]) {
	      		if(iht==2&&count1 ==0&&count2==0)
	      			{
	     	 			hRecElectronPt_EMC->Fill(mElectron->pt);
		 			count1++;
		 			count2++;
	      			}
		    isPE1[iht] = kTRUE;
		    if(nEta1>=mEmcenEtaCut&&nPhi1>=mEmcenPhiCut) {
		      isNSMD1[iht] = kTRUE;
		      if(zDist1>mEmceZDistCut[0]&&zDist1<mEmceZDistCut[1]&&
			  phiDist1>mEmcePhiDistCut[0]&&phiDist1<mEmcePhiDistCut[1]) {
			isDist1[iht] = kTRUE;
		      }
		    }
		  }
		}
	      }
	    }//for loop
	  }


	  Double_t eta2 = mElectron2->eta;
	  Double_t dca2 = mElectron2->dca;
	  Double_t nsigma2 = myGaus->GetRandom();
	  Double_t nHitsFit2 = mElectron2->nFitPts;
	  Double_t e2 = mElectron2->energy;
	  Double_t adc02 = mElectron2->adc0;
	  Double_t dsmAdc02 = mElectron2->dsmAdc0;
	  Double_t p2 = mElectron2->p;
	  Double_t pe2 = (e2>0)? p2/e2:9999;
	  Double_t pt2 = mElectron2->pt;
	  Double_t nEta2 = mElectron2->nEta;
	  Double_t nPhi2 = mElectron2->nPhi;
	  Double_t zDist2 = mElectron2->zDist;
	  Double_t phiDist2 = mElectron2->phiDist;
	  bool isTPC2[4], isTrg2[4], isAdc02[4], isPE2[4], isNSMD2[4], isDist2[4];
	  for(int iht=0;iht<4;iht++) {
	    isTPC2[iht] = kFALSE;
	    isTrg2[iht] = kFALSE;
	    isAdc02[iht] = kFALSE;
	    isPE2[iht] = kFALSE;
	    isNSMD2[iht] = kFALSE;
	    isDist2[iht] = kFALSE;
	  }
	  //TPC track
	  if(eta2>=mTpceEtaCut[0]&&eta2<=mTpceEtaCut[1] &&
	      nHitsFit2>=mTpceHitsFitCut &&
	      dca2<=mTpceDcaCut &&
	      nsigma2>mTpcenSigmaElectronCut[0]&&nsigma2<mTpcenSigmaElectronCut[1] &&
	      pt2<30.&&p2<30.) {
	    for(int iht=0;iht<4;iht++) {
	      if(pt2>=mTpcePtCut[iht]&&p2>=mTpcePCut[iht]) isTPC2[iht] = kTRUE;
	    }
	  }
	  //EMC Track
	  if(nHitsFit2>=mTpceHitsFitCut &&
	      dca2<=mEmceDcaCut &&
	      eta2>=mEmceEtaCut[0]&&eta2<=mEmceEtaCut[1] &&
	      nsigma2>mEmcenSigmaElectronCut[0]&&nsigma2<mEmcenSigmaElectronCut[1] &&
	      pt2<30) {
	    for(int iht=0;iht<4;iht++) {
	      if(pt2>mEmcePtCut[iht]&&e2>0&&dsmAdc02>mDsmAdcCut[iht]) {
		isTrg2[iht] = kTRUE;
		if(adc02>mEmceAdcCut[iht]) {
		  isAdc02[iht] = kTRUE;
		  if(pe2>mEmcePECut[0]&&pe2<mEmcePECut[1]) {
		    isPE2[iht] = kTRUE;
		    if(nEta2>=mEmcenEtaCut&&nPhi2>=mEmcenPhiCut) {
		      isNSMD2[iht] = kTRUE;
		      if(zDist2>mEmceZDistCut[0]&&zDist2<mEmceZDistCut[1]&&
			  phiDist2>mEmcePhiDistCut[0]&&phiDist2<mEmcePhiDistCut[1]) {
			isDist2[iht] = kTRUE;
		      }
		    }
		  }
		}
	      }
	    }//for loop
	  }



	  Double_t rcPt = jpsiRc.Pt();
	  Double_t rcY = jpsiRc.Rapidity();
	  Double_t rcM = jpsiRc.M();
	  //ht0
	  if((isTrg1[0]&&isTPC2[0])||(isTrg2[0]&&isTPC1[0])||(isTrg1[0]&&isTrg2[0])) {
	    hHt0JpsiTrg->Fill(rcY,rcPt,weight1);
	    //hHt0Jpsi3DTrg->Fill(rcY,vz,rcPt);
	  }
	  if((isAdc01[0]&&isTPC2[0])||(isAdc02[0]&&isTPC1[0])||(isAdc01[0]&&isAdc02[0])) {
	    hHt0JpsiAdc0->Fill(rcY,rcPt,weight1);
	    hHt0Jpsi3DAdc0->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isPE1[0]&&isTPC2[0])||(isPE2[0]&&isTPC1[0])||(isPE1[0]&&isPE2[0])) {
	    hHt0JpsiPE->Fill(rcY,rcPt,weight1);
	    hHt0Jpsi3DPE->Fill(rcY,vz,rcPt,weight1);
	    hHt0JpsiMassPE->Fill(rcY,rcPt,rcM,weight1);
	  }
	  if((isNSMD1[0]&&isTPC2[0])||(isNSMD2[0]&&isTPC1[0])||(isNSMD1[0]&&isNSMD2[0])) {
	    hHt0JpsiNSMD->Fill(rcY,rcPt,weight1);
	    hHt0Jpsi3DNSMD->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isDist1[0]&&isTPC2[0])||(isDist2[0]&&isTPC1[0])||(isDist1[0]&&isDist2[0])) {
	    hHt0JpsiDist->Fill(rcY,rcPt,weight1);
	    hHt0JpsiMassDist->Fill(rcY,rcPt,rcM,weight1);
	    hHt0Jpsi3DDist->Fill(rcY,vz,rcPt,weight1);
	  }

	  //ht1
	  if((isTrg1[1]&&isTPC2[1])||(isTrg2[1]&&isTPC1[1])||(isTrg1[1]&&isTrg2[1])) {
	    hHt1JpsiTrg->Fill(rcY,rcPt,weight1);
	    //hHt1Jpsi3DTrg->Fill(rcY,vz,rcPt);
	  }
	  if((isAdc01[1]&&isTPC2[1])||(isAdc02[1]&&isTPC1[1])||(isAdc01[1]&&isAdc02[1])) {
	    hHt1JpsiAdc0->Fill(rcY,rcPt,weight1);
	    hHt1Jpsi3DAdc0->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isPE1[1]&&isTPC2[1])||(isPE2[1]&&isTPC1[1])||(isPE1[1]&&isPE2[1])) {
	    hHt1JpsiPE->Fill(rcY,rcPt,weight1);
	    hHt1JpsiMassPE->Fill(rcY,rcPt,rcM,weight1);
	    hHt1Jpsi3DPE->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isNSMD1[1]&&isTPC2[1])||(isNSMD2[1]&&isTPC1[1])||(isNSMD1[1]&&isNSMD2[1])) {
	    hHt1JpsiNSMD->Fill(rcY,rcPt,weight1);
	    hHt1Jpsi3DNSMD->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isDist1[1]&&isTPC2[1])||(isDist2[1]&&isTPC1[1])||(isDist1[1]&&isDist2[1])) {
	    hHt1JpsiDist->Fill(rcY,rcPt,weight1);
	    hHt1JpsiMassDist->Fill(rcY,rcPt,rcM,weight1);
	    hHt1Jpsi3DDist->Fill(rcY,vz,rcPt,weight1);
	  }

	  //ht2
	  if((isTrg1[2]&&isTPC2[2])||(isTrg2[2]&&isTPC1[2])||(isTrg1[2]&&isTrg2[2])) {
	    hHt2JpsiTrg->Fill(rcY,rcPt,weight1);
	    //hHt2Jpsi3DTrg->Fill(rcY,vz,rcPt);
	  }
	  if((isAdc01[2]&&isTPC2[2])||(isAdc02[2]&&isTPC1[2])||(isAdc01[2]&&isAdc02[2])) {
	    hHt2JpsiAdc0->Fill(rcY,rcPt,weight1);
	    hHt2Jpsi3DAdc0->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isPE1[2]&&isTPC2[2])||(isPE2[2]&&isTPC1[2])||(isPE1[2]&&isPE2[2])) {
          //  hHt2JpsiPE_Or->Sumw2();
	  //  hHt2JpsiPE->Sumw2();
	  //  hHt2JpsiMassPE_1->Sumw2();
	  //  hHt2JpsiMassPE->Sumw2();
	  //  hHt2Jpsi3DPE->Sumw2();

            hHt2JpsiPE_Or->Fill(rcY,rcPt);
         // hweightPt_2->Fill(rcPt,weight1);
         // hweightPt_0->Fill(mcPt,weight1);
            hHt2JpsiPE->Fill(rcY,rcPt,weight1);
	    hHt2JpsiMassPE_1->Fill(rcY,rcPt,rcM,weight1);
	    hHt2JpsiMassPE->Fill(rcY,rcPt,rcM);
            hHt2Jpsi3DPE->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isNSMD1[2]&&isTPC2[2])||(isNSMD2[2]&&isTPC1[2])||(isNSMD1[2]&&isNSMD2[2])) {
	    hHt2JpsiNSMD->Fill(rcY,rcPt,weight1);
	    hHt2Jpsi3DNSMD->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isDist1[2]&&isTPC2[2])||(isDist2[2]&&isTPC1[2])||(isDist1[2]&&isDist2[2])) {
	    hHt2JpsiDist->Fill(rcY,rcPt,weight1);
	    hHt2JpsiMassDist->Fill(rcY,rcPt,rcM,weight1);
	    hHt2Jpsi3DDist->Fill(rcY,vz,rcPt,weight1);
	  }

	  //ht3
	  if((isTrg1[3]&&isTPC2[3])||(isTrg2[3]&&isTPC1[3])||(isTrg1[3]&&isTrg2[3])) {
	    hHt3JpsiTrg->Fill(rcY,rcPt,weight1);
	    //hHt3Jpsi3DTrg->Fill(rcY,vz,rcPt);
	  }
	  if((isAdc01[3]&&isTPC2[3])||(isAdc02[3]&&isTPC1[3])||(isAdc01[3]&&isAdc02[3])) {
	    hHt3JpsiAdc0->Fill(rcY,rcPt,weight1);
	    hHt3Jpsi3DAdc0->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isPE1[3]&&isTPC2[3])||(isPE2[3]&&isTPC1[3])||(isPE1[3]&&isPE2[3])) {
	    hHt3JpsiPE->Fill(rcY,rcPt,weight1);
	    hHt3JpsiMassPE->Fill(rcY,rcPt,rcM,weight1);
	    hHt3Jpsi3DPE->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isNSMD1[3]&&isTPC2[3])||(isNSMD2[3]&&isTPC1[3])||(isNSMD1[3]&&isNSMD2[3])) {
	    hHt3JpsiNSMD->Fill(rcY,rcPt,weight1);
	    hHt3Jpsi3DNSMD->Fill(rcY,vz,rcPt,weight1);
	  }
	  if((isDist1[3]&&isTPC2[3])||(isDist2[3]&&isTPC1[3])||(isDist1[3]&&isDist2[3])) {
	    hHt3JpsiDist->Fill(rcY,rcPt,weight1);
	    hHt3JpsiMassDist->Fill(rcY,rcPt,rcM,weight1);
	    hHt3Jpsi3DDist->Fill(rcY,vz,rcPt,weight1);
	  }
	}

	/*
	   if(mElectron->energy>0) {
	   cout<<"(pt, eta, phi)- = "<<mMcJpsi->eminusPt<<"  "<<mMcJpsi->eminusEta<<"  "<<mMcJpsi->eminusPhi<<endl;
	   cout<<"(EmcE, neta, nphi, zdist, phidist, mindist) = "<<mMcJpsi->eminusEmcE<<"  "<<mMcJpsi->eminusnEta<<"  "<<mMcJpsi->eminusnPhi<<"  "<<mMcJpsi->eminusZDist<<"  "<<mMcJpsi->eminusPhiDist<<"  "<<mMcJpsi->eminusMinDist<<endl;
	   }
	 */
      }//partner loop
    }//electron/positron loop
    hNJpsi->Fill(nJpsi);
  }//event loop

  char buf[1024];
  sprintf(buf,"%s_cent_%d_%d.root", OutFile, mCentCut1, mCentCut2);
  TFile *f = new TFile(buf,"recreate");
  f->cd();
  cout<<"writing to "<<buf<<endl;
  hMcVertexZ->Write();
  hMcVertexXY->Write();
  hRcVertexZ->Write();
  hVertexZdiff->Write();
  hRefMult->Write();
  hRefMultCut->Write();

  hNJpsi->Write();
 // hWeight->Write();
 // hmcPtweight->Write();
  hMcJpsiPt->Write();
  hMcJpsiPt_Or->Write();
  hMcJpsiY->Write();
  hMcJpsiPtY->Write();
  hMcJpsiPtY_Or->Write();
  hMcJpsiPhi->Write();
  hMcJpsiPtM_1->Write();
  hMcJpsiPtM->Write();

  hmcPtvsrcPt->Write();
   hmcPtvsrcPt_1->Write();
  hRcJpsiPt->Write();
  hRcJpsiPt_Or->Write();
  hRcJpsiY->Write();
  hRcJpsiPtY->Write();
  hRcJpsiPtY_Or->Write();
  hRcJpsiPhi->Write();
  hRcJpsiPtM->Write();


  hRcJpsiPtDiff_rc->Write();
  hRcJpsiPtDiff->Write();
  hRcJpsiYDiff->Write();
  hRcJpsiPhiDiff->Write();

  hJpsiMc->Write();

  hHt0JpsiTrg->Write();
  hHt0JpsiAdc0->Write();
  hHt0JpsiPE->Write();
  hHt0JpsiNSMD->Write();
  hHt0JpsiDist->Write();

  hHt1JpsiTrg->Write();
  hHt1JpsiAdc0->Write();
  hHt1JpsiPE->Write();
  hHt1JpsiNSMD->Write();
  hHt1JpsiDist->Write();

  hHt2JpsiTrg->Write();
  hHt2JpsiAdc0->Write();
  hHt2JpsiPE->Write();
  hHt2JpsiPE_Or->Write();
  hHt2JpsiNSMD->Write();
  hHt2JpsiDist->Write();

  hHt3JpsiTrg->Write();
  hHt3JpsiAdc0->Write();
  hHt3JpsiPE->Write();
  hHt3JpsiNSMD->Write();
  hHt3JpsiDist->Write();

  hHt0JpsiMassPE->Write();
  hHt0JpsiMassDist->Write();
  hHt1JpsiMassPE->Write();
  hHt1JpsiMassDist->Write();
  hHt2JpsiMassPE_1->Write();
  hHt2JpsiMassPE->Write();
  hHt2JpsiMassDist->Write();
  hHt3JpsiMassPE->Write();
  hHt3JpsiMassDist->Write();

  hJpsi3DMc->Write();

  //hHt0Jpsi3DTrg->Write();
  hHt0Jpsi3DAdc0->Write();
  hHt0Jpsi3DPE->Write();
  hHt0Jpsi3DNSMD->Write();
  hHt0Jpsi3DDist->Write();

  //hHt1Jpsi3DTrg->Write();
  hHt1Jpsi3DAdc0->Write();
  hHt1Jpsi3DPE->Write();
  hHt1Jpsi3DNSMD->Write();
  hHt1Jpsi3DDist->Write();

  //hHt2Jpsi3DTrg->Write();
  hHt2Jpsi3DAdc0->Write();
  hHt2Jpsi3DPE->Write();
  hHt2Jpsi3DNSMD->Write();
  hHt2Jpsi3DDist->Write();

  //hHt3Jpsi3DTrg->Write();
  hHt3Jpsi3DAdc0->Write();
  hHt3Jpsi3DPE->Write();
  hHt3Jpsi3DNSMD->Write();
  hHt3Jpsi3DDist->Write();
	hMCElectronPt->Write();
	hRecElectronPt_tpc->Write();
	hRecElectronPt_EMC->Write();
  f->Close();
  cout<<"done"<<endl;
}

//=======================================
int getCentrality(int refmult) {
  int mCentrality = -1;
  int   cent[] = { 7,9,11,12,14,15,17,19,23};//run10 AuAu200
  if(     refmult <= cent[0]) mCentrality = 0;
  else if(refmult <= cent[1]) mCentrality = 1;
  else if(refmult <= cent[2]) mCentrality = 2;
  else if(refmult <= cent[3]) mCentrality = 3;
  else if(refmult <= cent[4]) mCentrality = 4;
  else if(refmult <= cent[5]) mCentrality = 5;
  else if(refmult <= cent[6]) mCentrality = 6;
  else if(refmult <= cent[7]) mCentrality = 7;
  else if(refmult <= cent[8]) mCentrality = 8;
  else                        mCentrality = 9;

  return mCentrality;
}

