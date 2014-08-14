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

TH1F *hMcElectronPt   = new TH1F("mcElectronPt",  "mcElectronPt;p_{T}^{mc} (GeV/c)",  250, 0, 25.);
TH1F *hMcElectronY  = new TH1F("mcElectronY", "mcElectronY; Y^{mc};", 100, -5,  5.);
TH2F *hMcElectronPtY = new  TH2F("mcElectronPtY","mcElectronPtY;Y^{mc};p_{T}^{mc} (GeV/c)", 40, -2, 2, 250, 0, 25);
TH1F *hMcElectronPhi  = new TH1F("mcElectronPhi", "mcElectronPhi,#phi^{mc}", 360, -TMath::Pi(), TMath::Pi());

TH1F *hRcElectronPt   = new TH1F("rcElectronPt",  "rcElectronPt;p_{T}^{rc} (GeV/c)",  250, 0, 25.);
TH1F *hRcElectronEta  = new TH1F("rcElectronEta", "rcElectronEta;#eta^{rc};", 100, -5,  5.);
TH2F *hRcElectronPtEta = new TH2F("rcElectronPtEta", "rcElectronPtEta;#eta^{rc};p_{T}^{rc} (GeV/c)", 40, -2, 2, 250, 0, 25);
TH1F *hRcElectronPhi  = new TH1F("rcElectronPhi", "rcElectronPhi;#phi^{rc}", 360, -TMath::Pi(), TMath::Pi());
TH2F *hRcElectronPhiptNeg  = new TH2F("rcElectronPhiptNeg", "rcElectronPhi;p_{T} (GeV/c);#phi^{rc}",250,0,25., 360, -TMath::Pi(), TMath::Pi());
TH2F *hRcElectronPhiptPos  = new TH2F("rcElectronPhiptPos", "rcElectronPhi;p_{T} (GeV/c);#phi^{rc}",250,0,25., 360, -TMath::Pi(), TMath::Pi());

TH2F *hRcElectronnFitPtsPt  = new TH2F("rcElectronnFitPtsPt",  "rcElectronnFitPtsPt;p_{T}^{rc} (GeV/c);#HitsFit",  250, 0, 25., 56, -0.5, 55.5);
TH2F *hRcElectronnHitPtsPt = new TH2F("rcElectronnHitPtsPt","rcElectronnHitPtsPt;p_{T}^{rc} (GeV/c); #Hits", 250, 0, 25., 56, -0.5, 55.5);
TH2F *hRcElectronnSigEP    = new TH2F("rcElectronSigEP",  "rcElectronnSigEP;p^{rc} (GeV/c);n#sigma_{E}", 250, 0, 25.,  1000, -5, 5);
TH2F *hRcElectronDcaPt      = new TH2F("rcElectronDcaPt", "rcElectronDcaPt;p_{T}^{rc} (GeV/c);dca (cm)", 250, 0, 25,  100, 0, 5);
TH3F *hRcElectronPtDiff = new TH3F("rcElectronPtDiff", "rcElectronPtDiff;p_{T}^{mc}(GeV/c);#eta^{mc};p_{T}^{rc}/p_{T}^{mc}", 100, 0, 10, 40,-1,1, 1000, 0, 2);
TH3F *hRcElectronEtaDiff = new TH3F("rcElectronEtaDiff", "rcElectronEtaDiff;p_{T}^{mc}(GeV/c);#eta^{mc};#eta^{rc}-#eta^{mc}", 250,0,25,40, -1, 1, 1000, -0.1, 0.1);
TH3F *hRcElectronPhiDiff = new TH3F("rcElectronPhiDiff", "rcElectronPhiDiff;p_{T}^{mc}(GeV/c);#eta^{mc};#phi^{rc}-#phi^{mc}", 250,0,25,40, -1, 1, 1000, -0.1, 0.1);
TH2F *hDsmAdcvsPt  = new TH2F("hDsmAdcvsPt","DsmAdcvsPt;p_{T}^{mc} (GeV/c);DSM Adc",250,0,25,65,0,65);
TH2F *hAdc0vsPt    = new TH2F("hAdc0vsPt","Adc0vsPt;p_{T}^{mc} (GeV/c);Adc0",250,0,25,1000,0,1000);
TH2F *hHt0Adc0vsPt = new TH2F("hHt0Adc0vsPt","Ht0Adc0vsPt;p_{T}^{mc} (GeV/c);Adc0 (HT0)",250,0,25,1000,0,1000);
TH2F *hHt1Adc0vsPt = new TH2F("hHt1Adc0vsPt","Ht1Adc0vsPt;p_{T}^{mc} (GeV/c);Adc0 (HT1)",250,0,25,1000,0,1000);
TH2F *hHt2Adc0vsPt = new TH2F("hHt2Adc0vsPt","Ht2Adc0vsPt;p_{T}^{mc} (GeV/c);Adc0 (HT2)",250,0,25,1000,0,1000);
TH2F *hHt3Adc0vsPt = new TH2F("hHt3Adc0vsPt","Ht3Adc0vsPt;p_{T}^{mc} (GeV/c);Adc0 (HT3)",250,0,25,1000,0,1000);

TH2F *hHt0Adc0vsPt_1 = new TH2F("hHt0Adc0vsPt_1","Ht0Adc0vsPt_1;p_{T}^{mc} (GeV/c);Adc0(HT0)",250,0,25,1000,0,1000);
TH2F *hHt1Adc0vsPt_1 = new TH2F("hHt1Adc0vsPt_1","Ht1Adc0vsPt_1;p_{T}^{mc} (GeV/c);Adc0 (HT1)",250,0,25,1000,0,1000);
TH2F *hHt2Adc0vsPt_1 = new TH2F("hHt2Adc0vsPt_1","Ht2Adc0vsPt_1;p_{T}^{mc} (GeV/c);Adc0 (HT2)",250,0,25,1000,0,1000);
TH2F *hHt3Adc0vsPt_1 = new TH2F("hHt3Adc0vsPt_1","Ht3Adc0vsPt_1;p_{T}^{mc} (GeV/c);Adc0 (HT3)",250,0,25,1000,0,1000);

TH2F *hHt0Adc0vsPt_2 = new TH2F("hHt0Adc0vsPt_2","Ht0Adc0vsPt_2;p_{T}^{mc} (GeV/c);Adc0(HT0)",250,0,25,1000,0,1000);
TH2F *hHt1Adc0vsPt_2 = new TH2F("hHt1Adc0vsPt_2","Ht1Adc0vsPt_2;p_{T}^{mc} (GeV/c);Adc0 (HT1)",250,0,25,1000,0,1000);
TH2F *hHt2Adc0vsPt_2 = new TH2F("hHt2Adc0vsPt_2","Ht2Adc0vsPt_2;p_{T}^{mc} (GeV/c);Adc0 (HT2)",250,0,25,1000,0,1000);
TH2F *hHt3Adc0vsPt_2 = new TH2F("hHt3Adc0vsPt_2","Ht3Adc0vsPt_2;p_{T}^{mc} (GeV/c);Adc0 (HT3)",250,0,25,1000,0,1000);

TH3F *hHt0Adc0vsPtY = new TH3F("hHt0Adc0vsPtY","Ht0Adc0vsPtY;p_{T}^{mc} (GeV/c);#phi;Adc0 (HT0)",250,0,25,40,-1,1,1000,0,1000);
TH3F *hHt1Adc0vsPtY = new TH3F("hHt1Adc0vsPtY","Ht1Adc0vsPtY;p_{T}^{mc} (GeV/c);#phi;Adc0 (HT0)",250,0,25,40,-1,1,1000,0,1000);
TH3F *hHt2Adc0vsPtY = new TH3F("hHt2Adc0vsPtY","Ht2Adc0vsPtY;p_{T}^{mc} (GeV/c);#phi;Adc0 (HT0)",250,0,25,40,-1,1,1000,0,1000);
TH3F *hHt3Adc0vsPtY = new TH3F("hHt3Adc0vsPtY","Ht3Adc0vsPtY;p_{T}^{mc} (GeV/c);#phi;Adc0 (HT0)",250,0,25,40,-1,1,1000,0,1000);

TH2F *hPvsE = new TH2F("hPvsE","hpvsE;p^{rc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPmcvsE = new TH2F("hPmcvsE","hpmcvsE;p^{mc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPEvsPt = new TH2F("hPEvsPt","hPEvsPt;p_{T}^{rc} (Gev/c);p/E",250, 0, 25, 300, 0, 3);
TH2F *hnEtavsPt = new TH2F("hnEtavsPt","hnEtavsPt;p_{T}^{rc} (Gev/c);##eta",250,0,25,10,0,10);
TH2F *hnPhivsPt = new TH2F("hnPhivsPt","hnPhivsPt;p_{T}^{rc} (Gev/c);##phi",250,0,25,10,0,10);
TH2F *hZDistvsPt = new TH2F("hZDistvsPt","hZDistvsPt;p_{T}^{rc} (Gev/c);Z_{Dist}",250,0,25,200,-10,10);
TH2F *hPhiDistvsPt = new TH2F("hPhiDistvsPt","hPhiDistvsPt;p_{T}^{rc} (Gev/c);#phi_{Dist}",250,0,25,200,-0.1,0.1);

TH2F *hDsmAdc0vsAdc0 = new TH2F("hDsmAdc0vsAdc0","hDsmAdc0vsAdc0;Adc0;dsm Adc0",1000,0,1000,70,0,70);
TH2F *hPvsEHT0 = new TH2F("hPvsEHT0","hpvsEHT0;p^{rc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPmcvsEHT0 = new TH2F("hPmcvsEHT0","hpmcvsEHT0;p^{mc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPEvsPtHT0 = new TH2F("hPEvsPtHT0","hPEvsPtHT0;p_{T}^{rc} (Gev/c);p/E",250, 0, 25, 300, 0, 3);
TH2F *hnEtavsPtHT0 = new TH2F("hnEtavsPtHT0","hnEtavsPtHT0;p_{T}^{rc} (Gev/c);##eta",250,0,25,10,0,10);
TH2F *hnPhivsPtHT0 = new TH2F("hnPhivsPtHT0","hnPhivsPtHT0;p_{T}^{rc} (Gev/c);##phi",250,0,25,10,0,10);
TH2F *hZDistvsPtHT0 = new TH2F("hZDistvsPtHT0","hZDistvsPtHT0;p_{T}^{rc} (Gev/c);Z_{Dist}",250,0,25,200,-10,10);
TH2F *hPhiDistvsPtHT0 = new TH2F("hPhiDistvsPtHT0","hPhiDistvsPtHT0;p_{T}^{rc} (Gev/c);#phi_{Dist}",250,0,25,200,-0.1,0.1);

TH2F *hPvsEHT1 = new TH2F("hPvsEHT1","hpvsEHT1;p^{rc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPmcvsEHT1 = new TH2F("hPmcvsEHT1","hpmcvsEHT1;p^{mc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPEvsPtHT1 = new TH2F("hPEvsPtHT1","hPEvsPtHT1;p_{T}^{rc} (Gev/c);p/E",250, 0, 25, 300, 0, 3);
TH2F *hnEtavsPtHT1 = new TH2F("hnEtavsPtHT1","hnEtavsPtHT1;p_{T}^{rc} (Gev/c);##eta",250,0,25,10,0,10);
TH2F *hnPhivsPtHT1 = new TH2F("hnPhivsPtHT1","hnPhivsPtHT1;p_{T}^{rc} (Gev/c);##phi",250,0,25,10,0,10);
TH2F *hZDistvsPtHT1 = new TH2F("hZDistvsPtHT1","hZDistvsPtHT1;p_{T}^{rc} (Gev/c);Z_{Dist}",250,0,25,200,-10,10);
TH2F *hPhiDistvsPtHT1 = new TH2F("hPhiDistvsPtHT1","hPhiDistvsPtHT1;p_{T}^{rc} (Gev/c);#phi_{Dist}",250,0,25,200,-0.1,0.1);

TH2F *hPvsEHT2 = new TH2F("hPvsEHT2","hpvsEHT2;p^{rc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPmcvsEHT2 = new TH2F("hPmcvsEHT2","hpmcvsEHT2;p^{mc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPEvsPtHT2 = new TH2F("hPEvsPtHT2","hPEvsPtHT2;p_{T}^{rc} (Gev/c);p/E",250, 0, 25, 300, 0, 3);
TH2F *hnEtavsPtHT2 = new TH2F("hnEtavsPtHT2","hnEtavsPtHT2;p_{T}^{rc} (Gev/c);##eta",250,0,25,10,0,10);
TH2F *hnPhivsPtHT2 = new TH2F("hnPhivsPtHT2","hnPhivsPtHT2;p_{T}^{rc} (Gev/c);##phi",250,0,25,10,0,10);
TH2F *hZDistvsPtHT2 = new TH2F("hZDistvsPtHT2","hZDistvsPtHT2;p_{T}^{rc} (Gev/c);Z_{Dist}",250,0,25,200,-10,10);
TH2F *hPhiDistvsPtHT2 = new TH2F("hPhiDistvsPtHT2","hPhiDistvsPtHT2;p_{T}^{rc} (Gev/c);#phi_{Dist}",250,0,25,200,-0.1,0.1);

TH2F *hPvsEHT3 = new TH2F("hPvsEHT3","hpvsEHT3;p^{rc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPmcvsEHT3 = new TH2F("hPmcvsEHT3","hpmcvsEHT3;p^{mc} (GeV/c);E (GeV)",2500,0,25,2500,0,25);
TH2F *hPEvsPtHT3 = new TH2F("hPEvsPtHT3","hPEvsPtHT3;p_{T}^{rc} (Gev/c);p/E",250, 0, 25, 300, 0, 3);
TH2F *hnEtavsPtHT3 = new TH2F("hnEtavsPtHT3","hnEtavsPtHT3;p_{T}^{rc} (Gev/c);##eta",250,0,25,10,0,10);
TH2F *hnPhivsPtHT3 = new TH2F("hnPhivsPtHT3","hnPhivsPtHT3;p_{T}^{rc} (Gev/c);##phi",250,0,25,10,0,10);
TH2F *hZDistvsPtHT3 = new TH2F("hZDistvsPtHT3","hZDistvsPtHT3;p_{T}^{rc} (Gev/c);Z_{Dist}",250,0,25,200,-10,10);
TH2F *hPhiDistvsPtHT3 = new TH2F("hPhiDistvsPtHT3","hPhiDistvsPtHT3;p_{T}^{rc} (Gev/c);#phi_{Dist}",250,0,25,200,-0.1,0.1);

TH3F *hElectronMc = new TH3F("hElectronMc","hElectronMc;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hElectronTrk = new TH3F("hElectronTrk","hElectronTrk;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hElectronMatch = new TH3F("hElectronMatch","hElectronMatch;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);

TH3F *hHt0ElectronTrg = new TH3F("hHt0ElectronTrg","hHt0ElectronTrg;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt0ElectronAdc0 = new TH3F("hHt0ElectronAdc0","hHt0ElectronAdc0;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt0ElectronPE = new TH3F("hHt0ElectronPE","hHt0ElectronPE;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt0ElectronNSMD = new TH3F("hHt0ElectronNSMD","hHt0ElectronNSMD;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt0ElectronDist = new TH3F("hHt0ElectronDist","hHt0ElectronDist;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);

TH3F *hHt1ElectronTrg = new TH3F("hHt1ElectronTrg","hHt1ElectronTrg;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt1ElectronAdc0 = new TH3F("hHt1ElectronAdc0","hHt1ElectronAdc0;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt1ElectronPE = new TH3F("hHt1ElectronPE","hHt1ElectronPE;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt1ElectronNSMD = new TH3F("hHt1ElectronNSMD","hHt1ElectronNSMD;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt1ElectronDist = new TH3F("hHt1ElectronDist","hHt1ElectronDist;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);

TH3F *hHt2ElectronTrg = new TH3F("hHt2ElectronTrg","hHt2ElectronTrg;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt2ElectronAdc0 = new TH3F("hHt2ElectronAdc0","hHt2ElectronAdc0;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt2ElectronPE = new TH3F("hHt2ElectronPE","hHt2ElectronPE;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt2ElectronNSMD = new TH3F("hHt2ElectronNSMD","hHt2ElectronNSMD;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt2ElectronDist = new TH3F("hHt2ElectronDist","hHt2ElectronDist;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);

TH3F *hHt3ElectronTrg = new TH3F("hHt3ElectronTrg","hHt3ElectronTrg;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt3ElectronAdc0 = new TH3F("hHt3ElectronAdc0","hHt3ElectronAdc0;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt3ElectronPE = new TH3F("hHt3ElectronPE","hHt3ElectronPE;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt3ElectronNSMD = new TH3F("hHt3ElectronNSMD","hHt3ElectronNSMD;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);
TH3F *hHt3ElectronDist = new TH3F("hHt3ElectronDist","hHt3ElectronDist;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);

TH3F *hElectronTrk2 = new TH3F("hElectronTrk2","hElectronTrk2;#eta^{mc};#phi^{mc};p_{T}^{mc} (GeV/c)",40,-1.,1.,120,-TMath::Pi,TMath::Pi(),250,0,25);

TH2F *hdEtadPhi = new TH2F("hdEtadPhi","hdEtadPhi;#DeltaY;#Delta#phi",500,-2.5,2.5,628,-TMath::Pi(),TMath::Pi());
TH2F *hdEtadPhiHT1 = new TH2F("hdEtadPhiHT1","hdEtadPhiHT1;#DeltaY;#Delta#phi",500,-2.5,2.5,628,-TMath::Pi(),TMath::Pi());
TH2F *hdEtadPhiCut = new TH2F("hdEtadPhiCut","hdEtadPhiCut;#DeltaY;#Delta#phi",500,-2.5,2.5,628,-TMath::Pi(),TMath::Pi());
TH2F *hdEtadPhiHT1Cut = new TH2F("hdEtadPhiHT1Cut","hdEtadPhiHT0Cut;#DeltaY;#Delta#phi",500,-2.5,2.5,628,-TMath::Pi(),TMath::Pi());
TH2F *hDeltaRHT1vspt = new TH2F("hDeltaRHT1vspt","hDeltaRHT1vspt;p_{T} (GeV/c);dR",250,0,25,500,0,5);
TH2F *hDeltaRvspt = new TH2F("hDeltaRvspt","hDeltaRvspt;p_{T} (GeV/c);dR",250,0,25,500,0,5);
TH2F *hDeltaRHT1vsptCut = new TH2F("hDeltaRHT1vsptCut","hDeltaRHT1vsptCut;p_{T} (GeV/c);dR",250,0,25,500,0,5);
TH2F *hDeltaRvsptCut = new TH2F("hDeltaRvsptCut","hDeltaRvsptCut;p_{T} (GeV/c);dR",250,0,25,500,0,5);
TH3D *hCommonhitsvsRCPt = new TH3D("hCommonhitsvsRCPt","commonhits vs RC pT;tpc commonHits;RC p_{T} (GeV/c);#eta",50,0,50,300,0,30,120,-3,3);
TH3D *hCommonhitsvsMCPt = new TH3D("hCommonhitsvsMCPt","commonhits vs MC pT;tpc commonHits;MC p_{T} (GeV/c);#eta",50,0,50,300,0,30,120,-3,3);

void qa(const char* fileInDir="./", const char* OutFile="jpsi_qa_100", const int mCentCut1=0, const int mCentCut2=9){
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StMyElectronMaker");

	const Int_t nFilesMax = 10000;


	mElectronEvent = new StMyElectronEvent();

	TChain chMc("mcT");
	chMc.SetBranchAddress("mcE",&mElectronEvent);


	//	chMc.Add("mcJpsi.root");
	//char filename[128];
	//sprintf(filename,"/auto/stardata/starspec2/xzb/%s",file);
	const float pi = 3.1415926;
	void* dir = gSystem->OpenDirectory(gSystem->ExpandPathName(fileInDir));
	int nruns=0;
	char *file_name;
	TString Tname;
	char file_list[nFilesMax][256];
	do {
		file_name = (char*)gSystem->GetDirEntry(dir);
		//	cout<<" file_name "<<file_name<<endl;
		Tname=file_name;
		//    if (file_name && Tname.Contains("geant.root")) {
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
		for (int i = 0; i < nevents; i++) {
			nb +=chMc.GetEvent(i);
			chMc.GetEvent(i);
			//cout<<"haha1"<<endl;
			if(i%10000==0) cout<<"event# "<<i<<endl;
			int n;
			if (!mElectronEvent) continue;
			if (mElectronEvent->eventID()<=0) continue;
		//	cout<<"haha"<<endl;
			hMcVertexZ->Fill(mElectronEvent->mcVertexZ());
			hMcVertexXY->Fill(mElectronEvent->mcVertexX(), mElectronEvent->mcVertexY());
			hRcVertexZ->Fill(mElectronEvent->vertexZ());
			hVertexZdiff->Fill(mElectronEvent->vertexZ() - mElectronEvent->mcVertexZ());
			int Mult = mElectronEvent->refMultPos() + mElectronEvent->refMultNeg();
			//if(mElectron->pGeantId>0) continue;//not original electron
			hRefMult->Fill(Mult);
			int cent = getCentrality(Mult);
			if(cent<mCentCut1||cent>mCentCut2)continue;
			hRefMultCut->Fill(Mult);

			for(int j=0;j<mElectronEvent->nReal();j++){
				mElectron = (StMyElectron) mElectronEvent->real()->UncheckedAt(j);
		//		cout<<j<<"  "<<mElectron->geantId<<"  "<<mElectron->pGeantId<<endl;
				if(mElectron->pGeantId!=167) continue;//not from Jpsi electron
				//if(mElectron->pGeantId>0) continue;//not original electron
				if(mElectron->mcId<0) continue;
				hCommonhitsvsMCPt->Fill(mElectron->tpcCommonHits,mElectron->mcPt,mElectron->mcY);
				hCommonhitsvsRCPt->Fill(mElectron->tpcCommonHits,mElectron->pt,mElectron->mcY);
		//		cout<<" tpcCommonHits "<<mElectron->mElectron->tpcCommonHits<<endl;

				bool tag = kFALSE;
				for(int k=0;k<mElectronEvent->nReal();k++) {
					mElectron2 = (StMyElectron) mElectronEvent->real()->UncheckedAt(k);
					if(mElectron2->pGeantId!=167) continue;
					if(mElectron2->mcId<0) continue;
					if(mElectron2->mcId==mElectron->mcId) continue;
					if(mElectron2->pId==mElectron->pId) continue;

					Double_t deta = mElectron->mcY - mElectron2->mcY;
					Double_t dphi = mElectron->mcPhi - mElectron2->mcPhi;
					while(dphi>2*TMath::Pi()) dphi -= 2.*TMath::Pi();
					while(dphi<0) dphi += 2.*TMath::Pi();
					while(dphi>TMath::Pi()) dphi = dphi - 2*TMath::Pi();
					if(TMath::Abs(deta)<0.1&&TMath::Abs(dphi)<0.5) tag = kTRUE;
					if(k>j) {
						hdEtadPhi->Fill(deta,dphi);
						Double_t dr = sqrt(pow(deta,2)+pow(dphi,2));
						Double_t ptlow = (mElectron->pt<mElectron2->pt)? mElectron->pt:mElectron2->pt;
						//Double_t ptlow = mElectron2->pt;
						hDeltaRvspt->Fill(ptlow,dr);
						if(!tag) {
							hDeltaRvsptCut->Fill(ptlow,dr);
							hdEtadPhiCut->Fill(deta,dphi);
						}

						if(mElectron->dsmAdc0>18&&mElectron->adc0>290&&mElectron2->dsmAdc0>18&&mElectron2->adc0>290) {
							hDeltaRHT1vspt->Fill(ptlow,dr);
							hdEtadPhiHT1->Fill(deta,dphi);
							if(!tag)  {
								hDeltaRHT1vsptCut->Fill(ptlow,dr);
								hdEtadPhiHT1Cut->Fill(deta,dphi);
							}
						}
					}
				}
				if(tag) continue;

				if(mElectron->mcId>=0) {
					hMcElectronPt->Fill(mElectron->mcPt);
					hMcElectronY->Fill(mElectron->mcY);
					hMcElectronPtY->Fill(mElectron->mcY, mElectron->mcPt);
					hMcElectronPhi->Fill(mElectron->mcPhi);
				}


				if(mElectron->mcId>=0&&mElectron->id>=0){
					if(mElectron->tpcCommonHits>10){
						//cout<<mElectron->id<<endl;
						float pt = mElectron->pt;
						hRcElectronPt->Fill(pt);
						hRcElectronEta->Fill(mElectron->eta);
						hRcElectronPtEta->Fill(mElectron->eta, pt);
						hRcElectronPhi->Fill(mElectron->phi);
						//cout<<"haha"<<endl;
						if(fabs(mElectron->eta)<1){
							hRcElectronnFitPtsPt->Fill(pt, mElectron->nFitPts+1);
							hRcElectronnHitPtsPt->Fill(pt, mElectron->nHitPts+1);
							hRcElectronnSigEP->Fill(mElectron->p, mElectron->nSigE);
							hRcElectronDcaPt->Fill(pt, mElectron->dca);

							Float_t p = mElectron->p;
							Float_t poE = (mElectron->energy>0)? p/mElectron->energy:9999.;
							hPvsE->Fill(p,mElectron->energy);
							hPmcvsE->Fill(mElectron->mcPt*TMath::CosH(mElectron->mcEta),mElectron->energy);
							hPEvsPt->Fill(pt, poE);
							hnEtavsPt->Fill(pt, mElectron->nEta);
							hnPhivsPt->Fill(pt, mElectron->nPhi);
							hZDistvsPt->Fill(pt, mElectron->zDist);
							hPhiDistvsPt->Fill(pt, mElectron->phiDist);

							if(mElectron->geantId==2)
							{
								hRcElectronPhiptPos->Fill(pt,mElectron->phi);
							}
							if(mElectron->geantId==3)
							{
								hRcElectronPhiptNeg->Fill(pt,mElectron->phi);
							}

							Float_t ptMc = mElectron->mcPt;
							Float_t mcY = mElectron->mcY;
							hDsmAdcvsPt->Fill(ptMc, mElectron->dsmAdc0);
							hAdc0vsPt->Fill(ptMc, mElectron->adc0);
							hDsmAdc0vsAdc0->Fill(mElectron->adc0, mElectron->dsmAdc0);

							if(mElectron->dsmAdc0>11) {
								hHt0Adc0vsPt->Fill(ptMc, mElectron->adc0);
								hHt0Adc0vsPtY->Fill(ptMc, mcY, mElectron->adc0);
								hPvsEHT0->Fill(p,mElectron->energy);
								if(mElectron->nFitPts>=25&&mElectron->dca<1)
								{
									hHt0Adc0vsPt_1->Fill(ptMc,mElectron->adc0);
									if(mElectron->nEta>=2&&mElectron->nPhi>=2){
										if(fabs(mElectron->zDist)<2&&fabs(mElectron->phiDist)<0.01) {
											hHt0Adc0vsPt_2->Fill(ptMc,mElectron->adc0);
										}
									}
								}
								hPvsEHT0->Fill(p,mElectron->energy);
								hPmcvsEHT0->Fill(mElectron->mcPt*TMath::CosH(mElectron->mcEta),mElectron->energy);
								hPEvsPtHT0->Fill(pt, poE);
								hnEtavsPtHT0->Fill(pt, mElectron->nEta);
								hnPhivsPtHT0->Fill(pt, mElectron->nPhi);
								hZDistvsPtHT0->Fill(pt, mElectron->zDist);
								hPhiDistvsPtHT0->Fill(pt, mElectron->phiDist);
							}
							if(mElectron->dsmAdc0>15) {
								hHt1Adc0vsPt->Fill(ptMc, mElectron->adc0);
								hHt1Adc0vsPtY->Fill(ptMc, mcY, mElectron->adc0);
								hPvsEHT1->Fill(p,mElectron->energy);
								hPmcvsEHT1->Fill(mElectron->mcPt*TMath::CosH(mElectron->mcEta),mElectron->energy);
								hPEvsPtHT1->Fill(pt, poE);
								hnEtavsPtHT1->Fill(pt, mElectron->nEta);
								hnPhivsPtHT1->Fill(pt, mElectron->nPhi);
								hZDistvsPtHT1->Fill(pt, mElectron->zDist);
								hPhiDistvsPtHT1->Fill(pt, mElectron->phiDist);
								if(mElectron->nFitPts>=25&&mElectron->dca<1)
								{
									hHt1Adc0vsPt_1->Fill(ptMc,mElectron->adc0);
									if(mElectron->nEta>=2&&mElectron->nPhi>2){
										if(fabs(mElectron->zDist)<2&&fabs(mElectron->phiDist)<0.01) {
											hHt1Adc0vsPt_2->Fill(ptMc,mElectron->adc0);
										}
									}
								}
							}
							if(mElectron->dsmAdc0>18) {
								hHt2Adc0vsPt->Fill(ptMc, mElectron->adc0);
								hHt2Adc0vsPtY->Fill(ptMc, mcY, mElectron->adc0);
								hPvsEHT2->Fill(p,mElectron->energy);
								hPmcvsEHT2->Fill(mElectron->mcPt*TMath::CosH(mElectron->mcEta),mElectron->energy);
								hPEvsPtHT2->Fill(pt, poE);
								hnEtavsPtHT2->Fill(pt, mElectron->nEta);
								hnPhivsPtHT2->Fill(pt, mElectron->nPhi);
								hZDistvsPtHT2->Fill(pt, mElectron->zDist);
								hPhiDistvsPtHT2->Fill(pt, mElectron->phiDist);
								if(mElectron->nFitPts>=25&&mElectron->dca<1)
								{
									hHt2Adc0vsPt_1->Fill(ptMc,mElectron->adc0);
									if(mElectron->nEta>=2&&mElectron->nPhi>2){
										if(fabs(mElectron->zDist)<2&&fabs(mElectron->phiDist)<0.01){
											hHt2Adc0vsPt_2->Fill(ptMc,mElectron->adc0);
										}
									}
								}
							}
							if(mElectron->dsmAdc0>25) {
								hHt3Adc0vsPt->Fill(ptMc, mElectron->adc0);
								hHt3Adc0vsPtY->Fill(ptMc, mcY, mElectron->adc0);
								hPvsEHT3->Fill(p,mElectron->energy);
								hPmcvsEHT3->Fill(mElectron->mcPt*TMath::CosH(mElectron->mcEta),mElectron->energy);
								hPEvsPtHT3->Fill(pt, poE);
								hnEtavsPtHT3->Fill(pt, mElectron->nEta);
								hnPhivsPtHT3->Fill(pt, mElectron->nPhi);
								hZDistvsPtHT3->Fill(pt, mElectron->zDist);
								hPhiDistvsPtHT3->Fill(pt, mElectron->phiDist);
								if(mElectron->nFitPts>=25&&mElectron->dca<1)
								{
									hHt3Adc0vsPt_1->Fill(ptMc,mElectron->adc0);
									if(mElectron->nEta>=2&&mElectron->nPhi>2){
										if(fabs(mElectron->zDist)<2&&fabs(mElectron->phiDist)<0.01){
											hHt3Adc0vsPt_2->Fill(ptMc,mElectron->adc0);
										}
									}
								}
							}
						}//end of |eta|<1 cut
					}
				}//end of RC electron loop

				//calculate efficiency
				if(mElectron->mcId>=0&&fabs(mElectron->mcY)<1) {
					hElectronMc->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
					if(mElectron->id>=0 &&
							fabs(mElectron->eta)<1 &&
							mElectron->nFitPts>=25 &&
							mElectron->dca<1) {//Tracking
						hRcElectronPtDiff->Fill(mElectron->mcPt, mElectron->mcY, pt/mElectron->mcPt);
						hRcElectronEtaDiff->Fill(mElectron->mcPt, mElectron->mcY, mElectron->eta - mElectron->mcY);
						hRcElectronPhiDiff->Fill(mElectron->mcPt, mElectron->mcY, mElectron->phi - mElectron->mcPhi);

						hElectronTrk->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
						if(mElectron->energy>0) {
							hElectronMatch->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);

							if(mElectron->dsmAdc0>11) {
								hHt0ElectronTrg->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
								if(mElectron->adc0>180) {
									hHt0ElectronAdc0->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
									if(mElectron->p/mElectron->energy>0.3&&mElectron->p/mElectron->energy<1.5) {
										hHt0ElectronPE->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
										if(mElectron->nEta>=2&&mElectron->nPhi>=2){
											hHt0ElectronNSMD->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
											if(fabs(mElectron->zDist)<2&&fabs(mElectron->phiDist)<0.01) {
												hHt0ElectronDist->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
											}//distance
										}//#SMD strips
									}//PE
								}//Adc0
							}//Trigger


							if(mElectron->dsmAdc0>15) {
								hHt1ElectronTrg->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
								if(mElectron->adc0>250) {
									hHt1ElectronAdc0->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
									if(mElectron->p/mElectron->energy>0.3&&mElectron->p/mElectron->energy<1.5) {
										hHt1ElectronPE->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
										if(mElectron->nEta>=2&&mElectron->nPhi>=2) {
											hHt1ElectronNSMD->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
											if(fabs(mElectron->zDist)<2&&fabs(mElectron->phiDist)<0.01) {
												hHt1ElectronDist->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
											}//distance
										}//#SMD strips
									}//PE
								}//Adc0
							}//Trigger

							if(mElectron->dsmAdc0>18) {
								hHt2ElectronTrg->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
								if(mElectron->adc0>300) {
									hHt2ElectronAdc0->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
									if(mElectron->p/mElectron->energy>0.3&&mElectron->p/mElectron->energy<1.5) {
										hHt2ElectronPE->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
										if(mElectron->nEta>=2&&mElectron->nPhi>=2) {
											hHt2ElectronNSMD->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
											if(fabs(mElectron->zDist)<2&&fabs(mElectron->phiDist)<0.01) {
												hHt2ElectronDist->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
											}//distance
										}//#SMD strips
									}//PE
								}//Adc0
							}//Trigger

							if(mElectron->dsmAdc0>25) {
								hHt3ElectronTrg->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
								if(mElectron->adc0>400) {
									hHt3ElectronAdc0->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
									if(mElectron->p/mElectron->energy>0.3&&mElectron->p/mElectron->energy<1.5) {
										hHt3ElectronPE->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
										if(mElectron->nEta>=2&&mElectron->nPhi>=2) {
											hHt3ElectronNSMD->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
											if(fabs(mElectron->zDist)<2&&fabs(mElectron->phiDist)<0.01) {
												hHt3ElectronDist->Fill(mElectron->mcY,mElectron->mcPhi, mElectron->mcPt);
											}//distance
										}//#SMD strips
									}//PE
								}//Adc0
							}//Trigger


						}//EMC match
					}//Tracking

					if(mElectron->id>=0 &&
							fabs(mElectron->eta)<1 &&
							mElectron->nFitPts>=25 &&
							mElectron->dca<3) {
						hElectronTrk2->Fill(mElectron->mcY,mElectron->mcPhi,mElectron->mcPt);
					}//Tracking 2

				}
				/*
				   if(mElectron->energy>0) {
				   cout<<"(pt, eta, phi)- = "<<mMcJpsi->eminusPt<<"  "<<mMcJpsi->eminusEta<<"  "<<mMcJpsi->eminusPhi<<endl;
				   cout<<"(EmcE, neta, nphi, zdist, phidist, mindist) = "<<mMcJpsi->eminusEmcE<<"  "<<mMcJpsi->eminusnEta<<"  "<<mMcJpsi->eminusnPhi<<"  "<<mMcJpsi->eminusZDist<<"  "<<mMcJpsi->eminusPhiDist<<"  "<<mMcJpsi->eminusMinDist<<endl;
				   }
				   */
			}

		}
		char buf[1024];
		sprintf(buf,"%s_cent_%d_%d.root", OutFile, mCentCut1, mCentCut2);
		TFile *f = new TFile(buf,"recreate");
		f->cd();
		hMcVertexZ->Write();
		hMcVertexXY->Write();
		hRcVertexZ->Write();
		hVertexZdiff->Write();
		hRefMult->Write();
		hRefMultCut->Write();

		hMcElectronPt->Write();
		hMcElectronY->Write();
		hMcElectronPtY->Write();
		hMcElectronPhi->Write();

		hdEtadPhi->Write();
		hdEtadPhiHT1->Write();
		hdEtadPhiCut->Write();
		hdEtadPhiHT1Cut->Write();

		hRcElectronPt->Write();
		hRcElectronEta->Write();
		hRcElectronPtEta->Write();
		hRcElectronPhi->Write();
		hZDistvsPtHT2->Write();
		hPhiDistvsPtHT2->Write();

		hPvsEHT3->Write();
		hPmcvsEHT3->Write();
		hPEvsPtHT3->Write();
		hnEtavsPtHT3->Write();
		hnPhivsPtHT3->Write();
		hZDistvsPtHT3->Write();
		hPhiDistvsPtHT3->Write();

		hPvsEHT0->Write();
		hPmcvsEHT0->Write();
		hPEvsPtHT0->Write();
		hnEtavsPtHT0->Write();
		hnPhivsPtHT0->Write();
		hZDistvsPtHT0->Write();
		hPhiDistvsPtHT0->Write();

		hElectronMc->Write();
		hElectronTrk->Write();
		hElectronMatch->Write();

		hHt0ElectronTrg->Write();
		hHt0ElectronAdc0->Write();
		hHt0ElectronPE->Write();
		hHt0ElectronNSMD->Write();
		hHt0ElectronDist->Write();

		hHt1ElectronTrg->Write();
		hHt1ElectronAdc0->Write();
		hHt1ElectronPE->Write();
		hHt1ElectronNSMD->Write();
		hHt1ElectronDist->Write();

		hHt2ElectronTrg->Write();
		hHt2ElectronAdc0->Write();
		hHt2ElectronPE->Write();
		hHt2ElectronNSMD->Write();
		hHt2ElectronDist->Write();

		hHt3ElectronTrg->Write();
		hHt3ElectronAdc0->Write();
		hHt3ElectronPE->Write();
		hHt3ElectronNSMD->Write();
		hHt3ElectronDist->Write();

		hElectronTrk2->Write();
		hDeltaRHT1vspt->Write();
		hDeltaRvspt->Write();
		hDeltaRHT1vsptCut->Write();
		hDeltaRvsptCut->Write();

		hDsmAdcvsPt->Write();
		hDsmAdc0vsAdc0->Write();
		hAdc0vsPt->Write();
		hHt0Adc0vsPt->Write();
		hHt1Adc0vsPt->Write();
		hHt2Adc0vsPt->Write();
		hHt3Adc0vsPt->Write();

		hHt2Adc0vsPtY->Write();
		hPvsEHT2->Write();
		hPmcvsEHT2->Write();
		hPEvsPtHT2->Write();
		hRcElectronPtDiff->Write();
		hRcElectronEtaDiff->Write();
		hRcElectronPhiDiff->Write();

		hHt0Adc0vsPt_1->Write();
		hHt1Adc0vsPt_1->Write();
		hHt2Adc0vsPt_1->Write();
		hHt3Adc0vsPt_1->Write();

		hHt0Adc0vsPt_2->Write();
		hHt1Adc0vsPt_2->Write();
		hHt2Adc0vsPt_2->Write();
		hHt3Adc0vsPt_2->Write();

		hRcElectronnFitPtsPt->Write();
		hRcElectronnHitPtsPt->Write();
		hRcElectronnSigEP->Write();
		hRcElectronDcaPt->Write();
		hRcElectronPhiptNeg->Write();
		hRcElectronPhiptPos->Write();
  		hCommonhitsvsRCPt->Write();
  		hCommonhitsvsMCPt->Write();


		f->Close();
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

