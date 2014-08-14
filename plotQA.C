#define PI 3.1415927
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
#define twoPi 6.28318530
#define EMASS 0.000511


TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex){
	TLatex *latex = new TLatex(x,y,text);
	latex->SetNDC();
	latex->SetTextFont(textFont);
	latex->SetTextSize(textSize);
	latex->SetTextColor(colorIndex);
	latex->Draw("same");
	return latex;}  

Double_t PionDis(Double_t *x,Double_t *par)
{
	return par[0]/(sqrt(2*3.1415927)*par[2])*exp(-TMath::Power(x[0]-par[1],2)/(2.*TMath::Power(par[2],2)));
}
Double_t EleDis(Double_t *x,Double_t *par)
{
	return par[0]/(sqrt(2*3.1415927)*par[2])*exp(-TMath::Power(x[0]-par[1],2)/(2.*TMath::Power(par[2],2)));
}
Double_t KionDis(Double_t *x,Double_t *par)
{
	return par[0]/(sqrt(2*3.1415927)*par[2])*exp(-TMath::Power(x[0]-par[1],2)/(2.*TMath::Power(par[2],2)));
}
Double_t Disbut(Double_t *x,Double_t *par)
{
	return KionDis(x,par)+PionDis(x,&par[3])+EleDis(x,&par[6]);
}


void plotQA()
{ 
	
	//gSystem->Clear();
	gROOT->Reset();
	gStyle->Reset("plain");
	gStyle->SetOptStat(00000);
	gStyle->SetTitleX(1.5);
	gStyle->SetStatX(0.92);
	
	gStyle->SetMarkerColor(2);
	gStyle->SetPalette(1);
	gStyle->SetOptFit(0100);
	TCanvas *c1 = new TCanvas("c1","c1", 5, 5, 1024, 768);
	c1->SetFillColor(10);
	c1->SetBorderMode(0);
	c1->SetBorderSize(0);
	c1->SetFrameFillColor(10);
	c1->SetFrameBorderMode(0);
	c1->SetFrameBorderSize(0);
	c1->SetTopMargin(0.1);
	c1->SetRightMargin(0.12);
	c1->SetLeftMargin(0.1);
	c1->SetBottomMargin(0.1);
	//  c1->SetGrid(1,1);
	Float_t  num = 4.9;
	char buff[512];
	sprintf(buff,"Jpsi ~ %2.1f M",num);
	float max =0;

	TFile f("psiEff_2011_psi_cent_0_9.root");
	mcJpsiPt->Draw();
	mcJpsiPt->SetXTitle("MC psi(2s) p_{T}(GeV/c)");
	max = mcJpsiPt->GetMaximum();
	mcJpsiPt->GetYaxis()->SetRangeUser(0,1.2*max);
	c1->SaveAs("mcJpsiPt.gif");  

	mcJpsiY->GetXaxis()->SetRangeUser(-3,3);
	mcJpsiY->Draw();
	mcJpsiY->SetXTitle("MC psi(2s) Y^{mc}");
	TLatex l;
	l.SetTextSize(0.04);
	l.SetTextColor(kRed);
	//  sprintf(inf1,"Number of #psi(2s) ~ 260k");
	l.DrawLatex(-3.,10e4,buff);
	c1->SaveAs("mcJpsiY.gif");  

	mcVertexZ->Draw();
	mcVertexZ->SetXTitle("MC VertexZ (cm)");
	c1->SaveAs("mcVertexZ.gif");

	mcJpsiPhi->Draw();
	mcJpsiPhi->SetXTitle("MC psi(2s) #phi");
	max = mcJpsiPhi->GetMaximum();
	mcJpsiPhi->GetYaxis()->SetRangeUser(0,1.2*max);
	c1->SaveAs("JpsiPhi.gif");

	mcJpsiPtM->Draw("colz");
	mcJpsiPtM->GetYaxis()->SetRangeUser(3,4.5);
	mcJpsiPtM->SetXTitle("MC p_{T} GeV/c");
	mcJpsiPtM->SetYTitle("MC M_{inv}(e^{+}e^{-}) GeV/c^{2}");
	c1->SetLogz();
	c1->SaveAs("mcjpsiPtM.gif");

	mcVertexXY->Draw("colz");
	mcVertexXY->SetXTitle("MC VertexX (cm)");
	mcVertexXY->SetYTitle("MC VertexY (cm)");
	c1->SetLogz();
	c1->SaveAs("mcVertexXY.gif");

	mcJpsiPtY->Draw("colz");
	mcJpsiPtY->GetYaxis()->SetRangeUser(0,32);
	mcJpsiPtY->SetXTitle("MC psi(2s) rapidity");
	mcJpsiPtY->SetYTitle("MC psi(2s) p_{T}(GeV/c)");
	c1->SetLogz();
	c1->SaveAs("mcJpsiYvsPt.gif");

	hNJpsi->Draw();
	hNJpsi->SetXTitle("# of psi(2s)/event");
	c1->SaveAs("number.gif");

	refMult->SetXTitle("RefMult");
	refMult->GetXaxis()->SetRangeUser(0,60);
	refMult->Draw();
	c1->SetLogy(1);
	c1->SaveAs("refMult.gif");
	f.Close();

	c1->SetLogy(0);
	////////////////////////////////////////////////////////////////////////////////draw pt bins

	//  TFile *w = new TFile("jpsiqa_2011_Jpsi_100_cent_0_9.root");
	// hRcElectronnFitPtsPt->Draw("COLZ");

	// TH1F *h1 = (TH1F*)hRcElectronnFitPtsPt->ProjectionY("l1",0,300);
	// h1->Draw();

	//  TFile *H = new TFile("Jpsiqa_0.02_nSigE_histo.root");
	TFile *H = new TFile("psiqa_2011_psi_cent_0_9.root","READ");
	TFile *R = new TFile("psi_embedding_2.root","READ");
	TH2F *hRcElectronnFitPtsPt = (TH2F*)H->Get("rcElectronnFitPtsPt");
	TH2F *hRcElectronnHitPtsPt = (TH2F*)H->Get("rcElectronnHitPtsPt");
	TH2F *hRcElectronnSigEP    = (TH2F*)H->Get("rcElectronSigEP");
	TH2F *hRcElectronDcaPt     = (TH2F*)H->Get("rcElectronDcaPt");
//	TH2F *hpEvsPt              = (TH2F*)hPEvsPt->Clone();
	TH2F *hRcElectronPhiptNeg  = (TH2F*)H->Get("rcElectronPhiptNeg");
	TH2F *hRcElectronPhiptPos  = (TH2F*)H->Get("rcElectronPhiptPos");
	TH2F *hRcElectronPtEta     = (TH2F*)H->Get("rcElectronPtEta");

	TH2F *hHt0Adc0vsPt = (TH2F*)H->Get("hHt0Adc0vsPt");
	TH2F *hHt1Adc0vsPt = (TH2F*)H->Get("hHt1Adc0vsPt");
	TH2F *hHt2Adc0vsPt = (TH2F*)H->Get("hHt2Adc0vsPt");
	TH2F *hHt3Adc0vsPt = (TH2F*)H->Get("hHt3Adc0vsPt");
	
	TH2F *hDcavsPt = (TH2F*)R->Get("hDcavsPt");
	TH2F *hEtavsPt = (TH2F*)R->Get("hEtavsPt"); 
	TH2F *hPhivsPt = (TH2F*)R->Get("hPhivsPt"); 
	TH2F *hNegPhivsPt = (TH2F*)R->Get("hNegPhivsPt"); 
	TH2F *hPosPhivsPt = (TH2F*)R->Get("hPosPhivsPt"); 
	TH2F *hHitsFitvsPt = (TH2F*)R->Get("hHitsFitvsPt"); 
	TH2F *hDataHt1Adc0vsPt = (TH2F*)R->Get("hDataHt1Adc0vsPt");
	TH2F *hDataHt0Adc0vsPt = (TH2F*)R->Get("hDataHt0Adc0vsPt");
	TH2F *hDataHt2Adc0vsPt = (TH2F*)R->Get("hDataHt2Adc0vsPt");

	TH2F *hHt1PuritynSigEvsPt = (TH2F*)R->Get("hHt1PuritynSigEvsPt");
	TH2F *hHt0PuritynSigEvsPt = (TH2F*)R->Get("hHt0PuritynSigEvsPt");
	TH2F *hHt2PuritynSigEvsPt = (TH2F*)R->Get("hHt2PuritynSigEvsPt");

	TH1F *MCnFitPt;
	TH1F *DatanFitPt;

	TH1F *MCPhiPos;
	TH1F *DataPhiPos;
	TH1F *MCPhiNeg;
	TH1F *DataPhiNeg;  

	TH1F *MCEtaPos;
	TH1F *DataEtaPos;
	TH1F *MCEtaNeg;
	TH1F *DataEtaNeg;

	TH1F *MCnDca;
	TH1F *DataDca;
	TH1F *MCpoE;
	TH1F *DatapoE;

	TH1F *MCAdc;
	TH1F *DataAdc;


	Float_t ptMin = 3.0;
	Float_t ptMax = 3.5; 
	Int_t binMin = 0;
	Int_t binMax = 0;
	Int_t n = 0;
	char inf[512];
	char inf1[512];
	char buff[512];

	TH1F *h1 = (TH1F*)hRcElectronPhiptNeg->ProjectionY("o1",0,100);
	TH1F *h2 = (TH1F*)hRcElectronPhiptPos->ProjectionY("o2",0,100);
	TH1F *h3 = (TH1F*)hNegPhivsPt->ProjectionY("o3",0,100);
	TH1F *h4 = (TH1F*)hPosPhivsPt->ProjectionY("o4",0,100);

	float MCPhiMax = h1->GetMaximum();
	float MCPhiMax1 = h2->GetMaximum();
	float DataPhiMax = h3->GetMaximum();
	float DataPhiMax1 = h4->GetMaximum();

	Int_t mss = h1->GetMaximumBin();
	Int_t counts = h1->GetBinContent(mss);
	h1->GetYaxis()->SetRangeUser(0.,counts*1.4);
	h1->Draw();

	h2->SetLineColor(2);
	h2->Draw("same");

	h3->Scale(MCPhiMax/DataPhiMax);
	h3->SetLineColor(3);
	h3->Draw("same");

	h4->Scale(MCPhiMax1/DataPhiMax1);
	h4->SetLineColor(4);
	h4->Draw("same");

	sprintf(buff,"0<p_{T}<30 GeV/c)");

	drawLatex(0.35,0.90,inf1,42,0.05,1);

	TLegend *legend1 = new TLegend(.15,.73,.26,.90);
	legend1->SetTextFont(62);
	legend1->SetTextSize(0.04);
	legend1->SetBorderSize(0);
	legend1->SetFillStyle(0);
	legend1->AddEntry(h1,"MC e^{-}");
	legend1->AddEntry(h2,"MC e^{+}");
	legend1->AddEntry(h3,"Data e^{-}");
	legend1->AddEntry(h4,"Data e^{+}");
	legend1->Draw();
	c1->SaveAs("MCphiPos.gif");
/*
	for(Int_t i=0;i<24;i++)
	{
		binMin = hRcElectronnFitPtsPt->GetXaxis()->FindBin(ptMin+10e-5);
		binMax = hRcElectronnFitPtsPt->GetXaxis()->FindBin(ptMax-10e-5);

		cout<<"binMin="<<binMin<<endl;
		cout<<"binMax="<<binMax<<endl;

		//nHitFit QA;
		sprintf(inf,"nft1Mc%d",i);
		cout<<inf<<endl;
		MCnFitPt = (TH1F*)hRcElectronnFitPtsPt->ProjectionY(inf,binMin,binMax); 
		sprintf(buff,"nft1Data%d",i);
		DatanFitPt = (TH1F*)hHitsFitvsPt->ProjectionX(buff,binMin,binMax);

		Int_t Intg_low = MCnFitPt->FindBin(20);
		cout<<"Intg_low="<<Intg_low<<endl;
		Int_t Intg_hig = MCnFitPt->FindBin(50);

		Float_t DataMax = DatanFitPt->Integral(Intg_low,Intg_hig);
		cout<<"DataMax="<<DataMax<<endl;
		Float_t McMax = MCnFitPt->Integral(Intg_low,Intg_hig);
		cout<<"McMax="<<McMax<<endl;
		if(McMax<1)DataMax = 1;
		//DatanFitPt->Scale(McMax/DataMax);
		MCnFitPt->Scale(DataMax/McMax);
		MCnFitPt->Draw();
		cout<<"Intg_hig="<<Intg_hig<<endl;
		DatanFitPt->SetMarkerStyle(20);
		DatanFitPt->SetMarkerColor(2);
		DatanFitPt->SetLineColor(2);
		DatanFitPt->Draw("epsame");

		sprintf(inf1,"%2.1f<p_{T}(e)<%2.1f GeV/c",ptMin,ptMax);
		drawLatex(0.35,0.90,inf1,42,0.05,1);

		TLegend *legend = new TLegend(.15,.75,.26,.85);

		legend->SetTextFont(62);
		legend->SetTextSize(0.03);
		legend->SetBorderSize(0);
		legend->SetFillStyle(0);
		legend->AddEntry(MCnFitPt,"MC electron");
		legend->AddEntry(DatanFitPt,"Data :-2.0<n#sigma_{#pi}<2.0");
		legend->Draw();
		sprintf(inf,"McnFitPt%1d.gif",i);  
		c1->SaveAs(inf);

		//poE QA
		sprintf(buff,"poE1Mc%d",i);
		MCpoE = (TH1F*)hPEvsPtHT2->ProjectionY(buff,binMin,binMax);  
		sprintf(buff,"poE1Data%d",i);
		DatapoE = (TH1F*)hDataPEvsPtHT2->ProjectionY(buff,binMin,binMax);
		// MCpoE->Sumw2();
		DatapoE->Sumw2();
		MCpoE->Rebin(3);
		DatapoE->Rebin(3);
		Float_t McMax = MCpoE->GetMaximum();
		Float_t DataMax   = DatapoE->GetMaximum();
		if(DataMax == 0){
			cout<<"DataMax="<<DataMax<<endl;
			DataMax = 1;
		}

		TLegend *legendl = new TLegend(.15,.75,.26,.85);
		legendl->SetTextFont(62);
		legendl->SetTextSize(0.03);
		legendl->SetBorderSize(0);
		legendl->SetFillStyle(0);
		legendl->AddEntry(MCpoE,"MC electron");
		legendl->AddEntry(DatapoE,"Data :-2.0<n#sigma_{#pi}<2.0");
		legendl->Draw();

		DatapoE->Scale(McMax/DataMax);
		MCpoE->Draw();
		DatapoE->SetMarkerStyle(20);
		DatapoE->SetMarkerColor(2);
		DatapoE->SetLineColor(2);
		DatapoE->Draw("epsame");
		sprintf(inf1,"%2.1f<p_{T}(e)<%2.1f GeV/c",ptMin,ptMax);
		drawLatex(0.35,0.90,inf1,42,0.05,1);
		legendl->Draw();
		sprintf(inf,"HT2poEn%d.png",i);  
		c1->SaveAs(inf);
		//Dca
		sprintf(buff,"Dca1Mc%d",i);
		MCnDca = (TH1F*)hRcElectronDcaPt->ProjectionY(buff,binMin,binMax);  
		sprintf(buff,"Dca1Data%d",i);
		DataDca = (TH1F*)hDcavsPt->ProjectionX(buff,binMin,binMax);
		Int_t DcaIntg_low = MCnDca->FindBin(0);
		Int_t DcaIntg_hig = MCnDca->FindBin(5);
		Float_t McMax = MCnDca->Integral(DcaIntg_low,DcaIntg_hig);
		Float_t DataMax  = DataDca->Integral(DcaIntg_low,DcaIntg_hig);
		if(McMax<1)DataMax=1;
		MCnDca->Scale(DataMax/McMax);
		MCnDca->Draw();
		DataDca->SetMarkerStyle(20);
		DataDca->SetMarkerColor(2);
		DataDca->SetLineColor(2);
		DataDca->Draw("epsame");
		sprintf(inf1,"%2.1f<p_{T}(e)<%2.1f GeV/c",ptMin,ptMax);
		drawLatex(0.35,0.90,inf1,42,0.05,1);
		legend->Draw();
		sprintf(inf,"Dcan%d.gif",i);  
		c1->SaveAs(inf);


		//phi 
		sprintf(buff,"MCphi1Neg%d",i);
		MCPhiNeg = (TH1F*)hRcElectronPhiptNeg->ProjectionY(buff,binMin,binMax);  
		sprintf(buff,"MCphi1Pos%d",i);
		MCPhiPos = (TH1F*)hRcElectronPhiptPos->ProjectionY(buff,binMin,binMax);
		sprintf(buff,"Dataphi1Neg%d",i);
		DataPhiNeg = (TH1F*)hNegPhivsPt->ProjectionX(buff,binMin,binMax);
		sprintf(buff,"Dataphi1Pos%d",i);
		DataPhiPos = (TH1F*)hPosPhivsPt->ProjectionX(buff,binMin,binMax);
		MCPhiNeg->Rebin(3);
		MCPhiPos->Rebin(3);
		DataPhiNeg->Rebin(3);
		DataPhiPos->Rebin(3);


		Int_t PhiIntg_low = MCPhiNeg->FindBin(-3.14);
		Int_t PhiIntg_hig = MCPhiNeg->FindBin(3.14);
		Float_t McPhiMax = MCPhiNeg->Integral(PhiIntg_low,PhiIntg_hig);
		Float_t DataPhiMax  = DataPhiNeg->Integral(PhiIntg_low,PhiIntg_hig);
		Float_t McPhiMax1 = MCPhiPos->Integral(PhiIntg_low,PhiIntg_hig);
		Float_t DataPhiMax1  = DataPhiPos->Integral(PhiIntg_low,PhiIntg_hig);

		Int_t Datamss = MCPhiPos->GetMaximumBin();
		Int_t Datacounts = MCPhiPos->GetBinContent(Datamss);
		cout<<"MCcounts="<<Datacounts<<endl;
		MCPhiPos->GetYaxis()->SetRangeUser(0.,Datacounts*1.4);
		if(DataPhiMax<1)DataPhiMax=1;
		DataPhiPos->Scale(McPhiMax1/DataPhiMax1);
		MCPhiPos->SetLineColor(4);
		MCPhiPos->Draw();

		DataPhiPos->SetLineStyle(9);
		DataPhiPos->SetLineColor(3);
		DataPhiPos->Draw("same");

		if(DataPhiMax<1)DataPhiMax=1;
		DataPhiNeg->Scale(McPhiMax/(DataPhiMax));
		MCPhiNeg->Draw("same");

		DataPhiNeg->SetLineColor(2);
		DataPhiNeg->SetLineStyle(3);
		DataPhiNeg->Draw("same");

		sprintf(inf1,"%2.1f<p_{T}(e)<%2.1f GeV/c",ptMin,ptMax);
		drawLatex(0.35,0.90,inf1,42,0.05,1);

		TLegend *legend2 = new TLegend(.15,.73,.26,.90);
		legend2->SetTextFont(62);
		legend2->SetTextSize(0.04);
		legend2->SetBorderSize(0);
		legend2->SetFillStyle(0);
		legend2->AddEntry(MCPhiPos,"MC e^{+}");
		legend2->AddEntry(MCPhiNeg,"MC e^{-}");
		legend2->AddEntry(DataPhiPos,"Data #pi^{+}");
		legend2->AddEntry(DataPhiNeg,"Data #pi^{-}");
		legend2->Draw();
		sprintf(inf,"Phin%d.gif",i);  
		c1->SaveAs(inf);

		//eta
		sprintf(buff,"MCEtaNeg%d",i);
		MCEta = (TH1F*)hRcElectronPtEta->ProjectionX(buff,binMin,binMax);  
		sprintf(buff,"DataEtaNeg%d",i);
		DataEtaNeg = (TH1F*)hEtavsPt->ProjectionX(buff,binMin,binMax);

		int EtaIntg_low = MCEta->FindBin(-1.);
		int EtaIntg_hig = MCEta->FindBin(1.);

		float McEtaMax = MCEta->Integral(EtaIntg_low,EtaIntg_hig);
		cout<<"McEtaMax="<<McEtaMax<<endl;
		EtaIntg_low = DataEtaNeg->FindBin(-1.);
		EtaIntg_hig = DataEtaNeg->FindBin(1.);
		float DataEtaMax  = DataEtaNeg->Integral(EtaIntg_low,EtaIntg_hig);
		cout<<"DataEtaMax="<<DataEtaMax<<endl;

		Int_t Datamss = DataEtaNeg->GetMaximumBin();
		Int_t Datacounts = DataEtaNeg->GetBinContent(Datamss);
		cout<<"Datacounts="<<Datacounts<<endl;
		DataEtaNeg->GetYaxis()->SetRangeUser(0.,Datacounts*1.4);

		if(DataEtaMax<1)DataEtaMax=1;

		DataEtaNeg->SetMarkerStyle(20);
		DataEtaNeg->SetMarkerColor(2);
		DataEtaNeg->SetLineColor(2);
		DataEtaNeg->Draw("ep");

		MCEta->Scale(DataEtaMax/(McEtaMax*1.5));
		MCEta->Draw("same");

		sprintf(inf1,"%2.1f<p_{T}(e)<%2.1f GeV/c",ptMin,ptMax);
		drawLatex(0.35,0.90,inf1,42,0.05,1);

		TLegend *legend3 = new TLegend(.15,.73,.26,.90);

		legend3->SetTextFont(62);
		legend3->SetTextSize(0.04);
		legend3->SetBorderSize(0);
		legend3->SetFillStyle(0);
		legend3->AddEntry(MCEta,"MC e^{+}+e^{-}");
		legend3->AddEntry(DataEtaNeg,"Data #pi^{+}+#pi^{-}");
		legend3->Draw();
		sprintf(inf,"Etan%d.gif",i);  
		c1->SaveAs(inf);
		ptMin = ptMin + 0.5;
		ptMax = ptMax + 0.5;

	};
*/
	///////////////////////////////////////////////////////////////////
	Int_t i = 7;
	Float_t ptMin = 2.;
	Float_t ptMax = 2.5;

	TH2F *hMCAdc0vsPt = new TH2F("hMCAdc0vsPt","Ht0Adc0vsPt;p_{T}^{mc} (GeV/c);Adc0 (HT0)",250,0,25,1000,0,1000);
	TH1F *hRationPt   = new TH1F("hRationPt","hRationPt;p_{T} (GeV/c);Ratio",250,0,25);
	TH1F *hDataElectronPt = (TH1F*)hDataHt0Adc0vsPt->ProjectionX("DataElectronPt",0,1000);
	TH1F *hMCElectronPt   = (TH1F*)hHt0Adc0vsPt->ProjectionX("MCElectronP1t",0,1000);

	char *ch = "HT1: p_{T}data/embedding";
	for(Int_t j=0;j<250;j++)
	{
		Int_t MCPtCount = 0;
		Int_t DataPtCount = 0;
		Float_t ratio = 0.;
		MCPtCount = hMCElectronPt->GetBinContent(j);
		DataPtCount = hDataElectronPt->GetBinContent(j);
		if(MCPtCount==0)ratio = 1.;
		else
			ratio = (Float_t)DataPtCount/(Float_t)MCPtCount;
		cout<<DataPtCount<<" "<<MCPtCount<<" "<<ratio<<endl;
		hRationPt->SetBinContent(j,ratio);
	}

	TLatex *latex = new TLatex(8,30,ch);

	c1->SetLogy();
	hMCElectronPt->GetYaxis()->SetRangeUser(1e-1,10e+5);
	hMCElectronPt->SetLineColor(kRed);
	hMCElectronPt->Draw();
	hDataElectronPt->Draw("same");
	TLegend *legend3 = new TLegend(.15,.73,.26,.90);

	legend3->SetTextFont(62);
	legend3->SetTextSize(0.04);
	legend3->SetBorderSize(0);
	legend3->SetFillStyle(0);
	legend3->AddEntry(hMCElectronPt,"BHT0:MC electron p_{T}");
	legend3->AddEntry(hDataElectronPt,"BHT0:Data electron p_{T}");
	legend3->Draw();
	c1->SaveAs("HT1_Pt.gif");
	hRationPt->Draw();
	latex->Draw("same");
	c1->SaveAs("HT1PtRatio.gif");

	for(Int_t i=0;i<250;i++)
	{
		for(Int_t j=0;j<1000;j++)
		{
			Float_t MCElectron = 0;
			Float_t weight = 0;
			MCElectron = hHt0Adc0vsPt->GetBinContent(i,j);
			weight = hRationPt->GetBinContent(i);
			hMCAdc0vsPt->SetBinContent(i,j,MCElectron*weight);
		}
	}
	///////////////////////////////////////////////////////////////////check
	TH1F *h1 = (TH1F*)hMCAdc0vsPt->ProjectionX("sa",0,1000);
	TH1F *h2 = (TH1F*)hDataHt0Adc0vsPt->ProjectionX("sa1",0,1000);
	h1->Draw();
	h2->SetLineColor(kRed);
	h2->Draw("same");
	//////////////////////////////////////////////////////////////////



	c1->cd();
	c1->Clear();
	c1->Divide(2,2);
	Int_t k = 0;
	for(Int_t i=0;i<1;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hHt0Adc0vsPt->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hHt0Adc0vsPt->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"MCAdc%d",k);
			MCAdc = (TH1F*)hMCAdc0vsPt->ProjectionY(buff,PtbinMin,PtbinMax);
			sprintf(buff,"DataAdc%d",k);
			DataAdc = (TH1F*)hDataHt0Adc0vsPt->ProjectionY(buff,PtbinMin,PtbinMax);

			Float_t adcCut = 290;
			Float_t adcMax = 1000;
			Int_t AdcCut = MCAdc->FindBin(adcCut);
			Int_t AdcMax = MCAdc->FindBin(adcMax);

			Float_t MCadcMax = MCAdc->Integral(adcCut,adcMax);
			Float_t DataadcMax = DataAdc->Integral(adcCut,adcMax);
			int mcadcmax = MCAdc->GetMaximum();
			int dataadcmax = DataAdc->GetMaximum();
			cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
			cout<<"MCadcMax:DataadcMax="<<MCadcMax<<" "<<DataadcMax<<endl;
			cout<<"mcadcMax:dataadcMax="<<mcadcmax<<" "<<dataadcmax<<endl;
			cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
			if(j==3)
			{
				MCAdc->Scale(DataadcMax/MCadcMax);
			}
			MCAdc->Rebin(6);
			DataAdc->Rebin(6);
			Int_t MCmss = MCAdc->GetMaximumBin();
			Int_t MCcounts = MCAdc->GetBinContent(MCmss);
			MCAdc->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			MCAdc->GetXaxis()->SetTitle("Adc0");
			MCAdc->Draw();

			DataAdc->SetMarkerStyle(20);
			DataAdc->SetMarkerColor(2);
			DataAdc->SetLineColor(2);
			DataAdc->Draw("epsame");

			sprintf(buff,"BHT0:%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			TLegend *legend3 = new TLegend(.15,.73,.26,.90);
			legend3->SetTextFont(62);
			legend3->SetTextSize(0.04);
			legend3->SetBorderSize(0);
			legend3->SetFillStyle(0);
			legend3->AddEntry(MCAdc,"MC electron");
			legend3->AddEntry(DataAdc,"Data electron");
			legend3->Draw();

			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		sprintf(buff,"BHT0_Adc_pic_10_%d.gif",i); 
		c1->SaveAs(buff);

	}
	/////////////////////////////////////////////////////////////////

	//+-----------------------------------------------------+
	//| purity estimation                                   |
	//+-----------------------------------------------------+


	c1->cd();
	c1->Clear();
	c1->Divide(2,2);
	ptMin = 2.;
	ptMax = 2.5;
	Double_t purity = 0.;
	Double_t Totcounts = 0.;
	Double_t ecounts   = 0.;
	Int_t k = 0;


	TF1 *fit1 = new TF1("fit1",Disbut,-10,4,9);
	TF1 *fit2[16];
	TF1 *fit3[16];
	TF1 *fit4[16];

	for(Int_t i=0;i<0;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);

		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			c1->cd(j+1)->SetLogy();

			Int_t PtbinMin = hHt0PuritynSigEvsPt->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hHt0PuritynSigEvsPt->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"MCAdc%d",k);
			MCAdc = (TH1F*)hHt1PuritynSigEvsPt->ProjectionY(buff,PtbinMin,PtbinMax);

			Int_t MCmss = MCAdc->GetMaximumBin();
			Int_t MCcounts = MCAdc->GetBinContent(MCmss);
			MCAdc->GetYaxis()->SetRangeUser(1.,MCcounts*1.3);
			MCAdc->GetXaxis()->SetTitle("nSigE");

			MCAdc->SetMarkerStyle(21);
			MCAdc->SetMarkerSize(0.8); 
			MCAdc->SetStats(0);


			fit1->SetNpx(500);
			fit1->SetLineWidth(3);
			fit1->SetLineColor(kMagenta);
			fit1->SetParameters(500.,-7,-1.,500,-4,-.5,50,0,-.5);
			if(i>1)
			{
				fit1->SetParLimits(7,-0.01,0.02);
				if(i==2&&j==1)
				{fit1->SetParameters(200,-6.,-1.5,500,-4.,-1.0,50,0,-0.5);} 
				if(i==2&&j==2)
				{fit1->SetParameters(400,-5.,-1.5,1000,-4.,-1.0,100,0.5,-0.5);}                    
				if(i==2&&j==3)
				{fit1->SetParLimits(8,0.3,2);}
				if(i>2&&j>2)
				{fit1->SetParLimits(4,-2.5,-1.5);}
				if(i==3&&j==3)
				{fit1->SetParameters(500,-6.,-2.,500,-4.,-0.5,50,0,-0.5);}
				if(i==3&&j==2)
				{fit1->SetParameters(200,-6.,-1.5,500,-4.,-1.0,50,0,-0.5);}
			}
			MCAdc->Fit("fit1","RV+","ep");
			Int_t NDF = fit1->GetNDF();
			Float_t chi2 = fit1->GetChisquare();
			sprintf(buff,"#chi^{2}/NDF:%2.1f/%d",chi2,NDF);
			drawLatex(0.65,0.53,buff,62,0.04,1);
			Double_t par[9];
			fit1->GetParameters(par);
			cout<<ptMin<<"<pt<"<<ptMax<<endl;
			cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<endl;
			cout<<par[3]<<" "<<par[4]<<" "<<par[5]<<endl;
			cout<<par[6]<<" "<<par[7]<<" "<<par[8]<<endl;
			sprintf(buff,"fit2%d",k);
			fit2[k] = new TF1(buff,KionDis,-10,4,3);
			fit2[k]->SetParameters(par);
			fit2[k]->SetLineColor(kBlue);
			fit2[k]->SetNpx(500);
			fit2[k]->Draw("same");

			sprintf(buff,"fit3%d",k);
			fit3[k] = new TF1(buff,PionDis,-10,4,3);
			fit3[k]->SetParameters(&par[3]);
			fit3[k]->SetLineColor(kBlack);
			fit3[k]->SetNpx(500);
			fit3[k]->Draw("same");

			sprintf(buff,"fit4%d",k);
			fit4[k] = new TF1(buff,EleDis,-10,4,3);
			fit4[k]->SetParameters(&par[6]);
			fit4[k]->SetLineColor(kRed);
			fit4[k]->SetNpx(500);
			fit4[k]->Draw("same");

			TLine *a = new TLine(3.0,1.0,3.0,6.5e3);
			a->SetLineWidth(4);
			a->SetLineColor(7);
			a->Draw(); 

			TLine *b = new TLine(-0.5,1.0,-0.5,6.5e3);
			b->SetLineWidth(4);
			b->SetLineColor(7);

			TLine *c = new TLine(0.,1.0,0.,6.5e3);
			c->SetLineWidth(4);
			c->SetLineColor(7);

			TLine *d = new TLine(0.5,1.0,0.5,6.5e3);
			d->SetLineWidth(4);
			d->SetLineColor(7);
			if(ptMax<3.0)
			{
				b->Draw();
				ecounts = fit4[k]->Integral(-0.5,3.0);
				Totcounts = fit1->Integral(-0.5,3.0);
				purity = ecounts/Totcounts;
			}
			else if(ptMax<6.0)
			{
				c->Draw();
				ecounts = fit4[k]->Integral(0.,3.0);
				Totcounts = fit1->Integral(0.,3.0);
				purity = ecounts/Totcounts;
			}
			else 
			{
				d->Draw();
				ecounts = fit4[k]->Integral(0.5,3.0);
				Totcounts = fit1->Integral(0.5,3.0);
				purity = ecounts/Totcounts;
			}
			TLegend *legend3 = new TLegend(.70,.68,.80,.85);
			legend3->SetTextFont(62);
			legend3->SetTextSize(0.04);
			legend3->SetBorderSize(0);
			legend3->SetFillStyle(0);
			legend3->AddEntry(MCAdc,"data");
			legend3->AddEntry(fit2[k],"K+p");
			legend3->AddEntry(fit3[k],"#pi");
			legend3->AddEntry(fit4[k],"e");
			legend3->AddEntry(fit1,"3-Gaus fit");
			legend3->Draw();


			sprintf(buff,"purity:%3.3f",purity);
			drawLatex(0.15,0.80,buff,62,0.04,1);

			sprintf(buff,"BHT0:%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);
			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		sprintf(buff,"BHT0_nisgE_pic_c_%d.gif",i); 
		c1->SaveAs(buff);

	}
	/////////////////////////////////////////////////////////////////      
}

