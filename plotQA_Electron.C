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
#include "TPDF.h"
#include "TBrowser.h"
#include "TLegend.h"
#include "TSpecturm.h"
#include "TCanvas.h"
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

void plotQA()
{ 
	gROOT->Reset();
	gStyle->Reset("plain");
	gStyle->SetOptStat(00000);
	gStyle->SetTitleX(1.5);
	gStyle->SetStatX(0.92);
	gStyle->SetMarkerColor(2);
	gStyle->SetPalette(1);
	gStyle->SetOptFit(0100);
	TCanvas *c1 = new TCavnas("c1","c1",0,0,800,600);
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
	TString inFile = "psiElectonQA";
	char buff[512];
	float ptMin = 1.0;
	float ptMax = 1.5; 

	TPDf *ps = new TPDF(Form("%s.pdf",inFile.Data()),111);
	ps->Off():

	TFile f("psiEff_2011_psi_cent_0_9.root");
	c1->SetLogz(1);
	TH2F *hdEtadPhiHT1 = (TH2F*)f->Get("hdEtadPhiHT1");
	hdEtadPhiHT1->Draw("colz");
	pdfAction(c1,ps);

	c1->SetLogz(1);
	TH2F *hdEtadPhiHT1Cut = (TH2F*)f->Get("hdEtadPhiHT1Cut");
	hdEtadPhiHT1iCut->Draw("colz");
	pdfAction(c1,ps);

	c1->SetLogz(1);
	TH2F *hDeltaRHT1vspt = (TH2F*)f->Get("hDeltaRHT1vspt");
	hDeltaRHT1vspt->Draw("colz");
	pdfAction(c1,ps);

	c1->SetLogz(1);
	TH2F *hDeltaRHT1vsptCut = (TH2F*)f->Get("hDeltaRHT1vsptCut");
	hDeltaRHT1vsptCut->Draw("colz");
	pdfAction(c1,ps);

	c1->SetLogy(1);
	TH1F *hMcElectronPt = (TH1F*)f->Get("hMcElectronPt");
	TH1F *hRcElectronPt = (TH1F*)f->Get("hRcElectronPt");
	hRcElectronPt->SetLineColor(kRed);
	hMcElectronPt->Draw();
	hRcElectronPt->Draw("same");
	sprintf(buff,"Rc elctron pT");
	drawLatex(0.2,0.85,buff,62,0.06,1);
	sprintf(buff,"MC electron pT");
	drawLatex(0.2,0.65,buff,62,0.06,2);
	pdfAction(c1,ps);

	c1->SetLogy(0);
	TH1F *hMcElectronY = (TH1F*)f->Get("hMcElectronY");
	TH1F *hRcElectronEta = (TH1F*)f->Get("hRcElectronEta");
	hRcElectronEta->SetLineColor(kRed);
	hMcElectronY->Draw();
	hRcElectronEta->Draw("same");
	sprintf(buff,"Rc elctron eta");
	drawLatex(0.2,0.85,buff,62,0.06,1);
	sprintf(buff,"MC electron eta");
	drawLatex(0.2,0.65,buff,62,0.06,2);
	pdfAction(c1,ps);

	c1->SetLogy(0);
	TH1F *hMcElectronPhi = (TH1F*)f->Get("hMcElectronPhi");
	TH1F *hRcElectronPhi = (TH1F*)f->Get("hRcElectronPhi");
	hRcElectronPhi->SetLineColor(kRed);
	hMcElectronPhi->Draw();
	hRcElectronPhi->Draw("same");
	sprintf(buff,"Rc elctron phi");
	drawLatex(0.2,0.85,buff,62,0.06,1);
	sprintf(buff,"MC electron phi");
	drawLatex(0.2,0.65,buff,62,0.06,2);
	pdfAction(c1,ps);

	TH2F *hMcElectronPtY = (TH2F*)f->Get("hMcElectronPtY");
	TH2F *hRcElectronPtEta = (TH2F*)f->Get("hRcElectronPtEta");
	c1->Clear();
	c1->Divide(2,2);
	int k = 0;
	int maxValue =-999;
	TH1F *MCPhi;
	TH1F *RcPhi;
	for(Int_t i=0;i<5;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hMcElectronPtY->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hMcElectronPtY->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"McPhi%d",k);
			MCPhi = (TH1F*)hMcElectronPtY->ProjectionY(buff,PtbinMin,PtbinMax);
			sprintf(buff,"RcPhi%d",k);
			RcPhi = (TH1F*)hRcElectronPtEta->ProjectionY(buff,PtbinMin,PtbinMax);

			RcPhi->SetLineColor(kRed);
			maxValue = MCPhi->GetMaximum();
			MCPhi->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			MCPhi->Draw();
			RcPhi->Draw("same");
			sprintf(buff,"%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		c1->cd();
		sprintf(buff,"Rc elctron eta");
		drawLatex(0.2,0.85,buff,62,0.06,1);
		sprintf(buff,"MC electron eta");
		drawLatex(0.2,0.65,buff,62,0.06,2);
		pdfAction(c1,ps);

		//sprintf(buff,"electron_phi_%d.gf",i); 
		//c1->SaveAs(buff);
	}

	TH2F *hRcElectronnFitPtsPt = (TH2F*)f->Get("hRcElectronnFitPtsPt");
	TH2F *hRcElectronnHitPtsPt = (TH2F*)f->Get("hRcElectronnHitPtsPt");
	c1->Clear();
	c1->Divide(2,2);
	k = 0;
	maxValue =-999;
	TH1F *MCFit;
	TH1F *RcHit;
	for(Int_t i=0;i<5;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hRcElectronnFitPtsPt->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hRcElectronnFitPtsPt->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"McFit%d",k);
			MCFit = (TH1F*)hRcElectronnFitPtsPt->ProjectionY(buff,PtbinMin,PtbinMax);
			sprintf(buff,"RcHit%d",k);
			RcHit = (TH1F*)hRcElectronnHitPtsPt->ProjectionY(buff,PtbinMin,PtbinMax);

			RcHit->SetLineColor(kRed);
			maxValue = MCFit->GetMaximum();
			MCFit->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			MCFit->Draw();
			RcHit->Draw("same");
			sprintf(buff,"%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		c1->cd();
		sprintf(buff,"Rc Hit pT");
		drawLatex(0.2,0.85,buff,62,0.06,1);
		sprintf(buff,"MC Fit pT");
		drawLatex(0.2,0.65,buff,62,0.06,2);
		pdfAction(c1,ps);

		//sprintf(buff,"electron_phi_%d.gf",i); 
		//c1->SaveAs(buff);
	}

	TH2F *hRcElectronDcaPt = (TH2F*)f->Get("hRcElectronDcaPt");
	c1->Clear();
	c1->Divide(2,2);
	k = 0;
	maxValue =-999;
	TH1F *MCDca;
	for(Int_t i=0;i<5;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hRcElectronDcaPt->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hRcElectronDcaPt->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"Dca%d",k);
			MCDca = (TH1F*)hRcElectronDcaPt->ProjectionY(buff,PtbinMin,PtbinMax);

			maxValue = MCDca->GetMaximum();
			MCDca->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			MCDca->Draw();
			sprintf(buff,"%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		c1->cd();
		sprintf(buff,"Rc Dca");
		drawLatex(0.2,0.85,buff,62,0.06,1);
		pdfAction(c1,ps);
	}

	c1->Clear();
	TH2F *hDsmAdcvsPt = (TH2F*)f->Get("hDsmAdcvsPt");
	hDsmAdcvsPt->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->Divide(2,2);
	k = 0;
	maxValue =-999;
	TH1F *DsmAdc;
	for(Int_t i=0;i<5;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hDsmAdcvsPt->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hDsmAdcvsPt->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"DsmAdcvspT%d",k);
			DsmAdc = (TH1F*)hDsmAdcvsPt->ProjectionY(buff,PtbinMin,PtbinMax);

			maxValue = DsmAdc->GetMaximum();
			DsmAdc->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			DsmAdc->Draw();
			sprintf(buff,"MC:%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		c1->cd();
		sprintf(buff,"DsmAdc");
		drawLatex(0.2,0.85,buff,62,0.06,1);
		pdfAction(c1,ps);
	}

	c1->Clear();
	TH2F *hDsmAdc0vsAdc0 = (TH2F*)f->Get("hDsmAdc0vsAdc0");
	hDsmAdc0vsAdc0->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->Divide(2,2);
	k = 0;
	maxValue =-999;
	dsmadc_Min = 5;
	dsmadc_Max = 10;
	TH1F *Adc0vsDsm;
	for(Int_t i=0;i<3;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hDsmAdc0vsAdc0->GetXaxis()->FindBin(dsmadc_Min);
			Int_t PtbinMax = hDsmAdc0vsAdc0->GetXaxis()->FindBin(dsmadc_Max);
			sprintf(buff,"Adc0vsDsmAdc0%d",k);
			Adc0vsDsm = (TH1F*)hDsmAdc0vsAdc0->ProjectionX(buff,PtbinMin,PtbinMax);

			maxValue = Adc0vsDsm->GetMaximum();
			Adc0vsDsm->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			Adc0vsDsm->Draw();
			sprintf(buff,"%2.1f<DsmAdc0<%2.1f GeV/c",dsmadc_Min,dsmadc_Max);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			dsmadc_Min = dsmadc_Min + 5;
			dsmadc_Max = dsmadc_Max + 5;
			k++;
		}
		c1->cd();
		sprintf(buff,"Adc0");
		drawLatex(0.2,0.85,buff,62,0.06,1);
		pdfAction(c1,ps);
	}


	TH2F *hPvsEHT2 = (TH2D*)f->Get("hPvsEHT2");
	TH2F *hPmcvsEHT2 = (TH2D*)f->Get("hPmcvsEHT2");
	c1->Clear();
	c1->Divide(2,2);
	k = 0;
	maxValue =-999;
	TH1F *RCPvsE;
	TH1F *RCPmvsE;
	for(Int_t i=0;i<5;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hPvsEHT2->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hPvsEHT2->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"PvsE%d",k);
			RCPvsE = (TH1F*)hPvsEHT2->ProjectionY(buff,PtbinMin,PtbinMax);
			sprintf(buff,"PmvsE%d",k);
			RCPmvsE = (TH1F*)hPmcvsEHT2->ProjectionY(buff,PtbinMin,PtbinMax);

			maxValue = RCPvsE->GetMaximum();
			RCPvsE->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			RCPvsE->Draw();
			RCPmvsE->SetLineColor(kRed);
			RCPmvsE->Draw("same");
			sprintf(buff,"%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		c1->cd();
		sprintf(buff,"HT1 :RC p vs E");
		drawLatex(0.6,0.85,buff,62,0.06,1);
		sprintf(buff,"HT1: RC pm vs E");
		drawLatex(0.6,0.65,buff,62,0.06,2);
		pdfAction(c1,ps);
	}

	TH2F *hPEvsPtHT2 = (TH2F*)f->Get("hPEvsPtHT2");
	c1->Clear();
	c1->Divide(2,2);
	k = 0;
	maxValue =-999;
	TH1F *RCPE; 
	for(Int_t i=0;i<5;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hPEvsPtHT2->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hPEvsPtHT2->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"PEvspT%d",k);
			RCPE = (TH1F*)hPEvsPtHT2->ProjectionY(buff,PtbinMin,PtbinMax);

			maxValue = RCPE->GetMaximum();
			RCPE->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			RCPE->Draw();
			sprintf(buff,"%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		c1->cd();
		sprintf(buff,"HT1: RC pE");
		drawLatex(0.6,0.85,buff,62,0.06,1);
		pdfAction(c1,ps);
	}

	c1->Clear();
	TH2F *hHt2Adc0vsPt = (TH2F*)f->Get("hHt2Adc0vsPt");
	hHt2Adc0vsPt->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->Divide(2,2);
	k = 0;
	maxValue =-999;
	TH1F *Adc0;
	for(Int_t i=0;i<5;i++)
	{
		c1->cd();
		c1->Clear();
		c1->Divide(2,2);
		for(Int_t j=0;j<4;j++)
		{
			c1->cd(j+1);
			Int_t PtbinMin = hHt2Adc0vsPt->GetXaxis()->FindBin(ptMin);
			Int_t PtbinMax = hHt2Adc0vsPt->GetXaxis()->FindBin(ptMax);
			sprintf(buff,"Adc0vspT%d",k);
			Adc0 = (TH1F*)hHt2Adc0vsPt->ProjectionY(buff,PtbinMin,PtbinMax);

			maxValue = Adc0->GetMaximum();
			Adc0->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
			Adc0->Draw();
			sprintf(buff,"MC:%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
			drawLatex(0.25,0.90,buff,42,0.05,1);

			ptMin = ptMin + 0.5;
			ptMax = ptMax + 0.5;
			k++;
		}
		c1->cd();
		sprintf(buff,"Adc0");
		drawLatex(0.6,0.85,buff,62,0.06,1);
		pdfAction(c1,ps);
	}


	TH3F *hRcElectronPtDiff = (TH3F*)f->Get("hRcElectronPtDiff");
	TH2F *hRcPtDiff = (TH2F*)hRcElectronPtDiff->Project3D("xz");
	hRcPtDiff->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->cd();
	c1->SetLogy(1);
	k = 0;
	maxValue =-999;
	ptMin =1.;
	ptMax =1.5;
	TH1F *ptDiff;
	for(Int_t i=0;i<15;i++)
	{
		Int_t PtbinMin = hRcPtDiff->GetXaxis()->FindBin(ptMin);
		Int_t PtbinMax = hRcPtDiff->GetXaxis()->FindBin(ptMax);
		sprintf(buff,"ptDiff%d",k);
		ptDiff = (TH1F*)hRcPtDiff->ProjectionY(buff,PtbinMin,PtbinMax);

		maxValue = ptDiff->GetMaximum();
		ptDiff->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
		if(k>0)
		{
			ptDiff->SetLineColor(KGreen-9+k)
				ptDiff->Draw("same");
		}
		else
		{
			ptDiff->SetLineColor(KGreen-9)
				ptDiff->Draw();
		}
		sprintf(buff,"RC:%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
		drawLatex(0.15,0.90-k*0.05,buff,42,0.05,1);

		ptMin = ptMin + 0.5;
		ptMax = ptMax + 0.5;
		k++;
	}
	c1->cd();
	sprintf(buff,"pT diff");
	drawLatex(0.6,0.85,buff,62,0.06,1);
	pdfAction(c1,ps);

	TH3F *hRcElectronEtaDiff = (TH3F*)f->Get("hRcElectronEtaDiff");
	TH2F *hRcEtaDiff = (TH2F*)hRcElectronEtaDiff->Project3D("xz");
	hRcEtaDiff->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->cd();
	c1->SetLogy(1);
	k = 0;
	maxValue =-999;
	ptMin =1.;
	ptMax =1.5;
	TH1F *etaDiff;
	for(Int_t i=0;i<15;i++)
	{
		Int_t PtbinMin = hRcEtaDiff->GetXaxis()->FindBin(ptMin);
		Int_t PtbinMax = hRcEtaDiff->GetXaxis()->FindBin(ptMax);
		sprintf(buff,"etaDiff%d",k);
		etaDiff = (TH1F*)hRcEtaDiff->ProjectionY(buff,PtbinMin,PtbinMax);

		maxValue = etaDiff->GetMaximum();
		etaDiff->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
		if(k>0)
		{
			etaDiff->SetLineColor(KGreen-9+k)
				etaDiff->Draw("same");
		}
		else
		{
			etaDiff->SetLineColor(KGreen-9)
			etaDiff->Draw();
		}
		sprintf(buff,"RC:%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
		drawLatex(0.15,0.90-k*0.05,buff,42,0.05,1);

		ptMin = ptMin + 0.5;
		ptMax = ptMax + 0.5;
		k++;
	}
	sprintf(buff,"Eta diff");
	drawLatex(0.6,0.85,buff,62,0.06,1);
	pdfAction(c1,ps);

	TH3F *hRcElectronPhiDiff = (TH3F*)f->Get("hRcElectronPhiDiff");
	TH2F *hRcPhiDiff = (TH2F*)hRcElectronPhiDiff->Project3D("xz");
	hRcPhiDiff->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->cd();
	c1->SetLogy(1);
	k = 0;
	maxValue =-999;
	ptMin =1.;
	ptMax =1.5;
	TH1F *phiDiff;
	for(Int_t i=0;i<15;i++)
	{
		int PtbinMin = hRcPhiDiff->GetXaxis()->FindBin(ptMin);
		int PtbinMax = hRcPhiDiff->GetXaxis()->FindBin(ptMax);
		sprintf(buff,"phiDiff%d",k);
		phiDiff = (TH1F*)hRcPhiDiff->ProjectionY(buff,PtbinMin,PtbinMax);

		maxValue = phiDiff->GetMaximum();
		phiDiff->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
		if(k>0)
		{
			phiDiff->SetLineColor(KGreen-9+k)
				phiDiff->Draw("same");
		}
		else
		{
			phiDiff->SetLineColor(KGreen-9)
				phiDiff->Draw();
		}
		sprintf(buff,"RC:%2.1f<p_{T}<%2.1f GeV/c",ptMin,ptMax);
		drawLatex(0.15,0.90-k*0.05,buff,42,0.05,1);

		ptMin = ptMin + 0.5;
		ptMax = ptMax + 0.5;
		k++;
	}
	sprintf(buff,"pT diff");
	drawLatex(0.6,0.85,buff,62,0.06,1);
	pdfAction(c1,ps);
	
	TH2F *hptDiffvsEta = new TH2F("hptDiffvsEta","rcElectronPt Pt Diff;#eta^{mc};p_{T}^{rc}/p_{T}^{mc}",40,-1,1,1000,0,2);
	TH2F *hetaDiffvsEta = new TH2F("hetaDiffvsEta","rcElectronPt eta Diff;#eta^{mc};#eta^{rc}-#eta^{mc}",40,-1,1,1000,-0.1,0.1);
	TH2F *hphiDiffvsEta = new TH2F("hphiDiffvsEta","rcElectronPt eta Diff;#eta^{mc};#phi^{rc}-#phi^{mc}",40,-1,1,1000,-0.1,0.1);
	int temPass;
	for(int i=1; i<41; i++)
	{
		for(int j=1; j<1001; j++)
		{
			temPass = hRcElectronPtDiff->Integral(11,250,i,i,j,j);//pT > 1 GeV/c
			hptDiffvsEta->SetBinContent(temPass,i,j);

			temPass = hRcElectronEtaDiff->Integral(11,250,i,i,j,j);//pT > 1 GeV/c
			hetaDiffvsEta->SetBinContent(temPass,i,j);

			temPass = hRcElectronPhiDiff->Integral(11,250,i,i,j,j);//pT > 1 GeV/c
			hphiDiffvsEta->SetBinContent(temPass,i,j);
		}
	}

	hptDiffvsEta->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->cd();
	c1->SetLogy(1);
	k = 0;
	maxValue =-999;
	EtaMin =-1.;
	EtaMax =-0.9;
	TH1F *PtDiff;
	for(Int_t i=0;i<10;i++)
	{
		int PtbinMin = hptDiffvsEta->GetXaxis()->FindBin(EtaMin);
		int PtbinMax = hptDiffvsEta->GetXaxis()->FindBin(EtaMax);
		sprintf(buff,"PtDiff%d",k);
		PtDiff = (TH1F*)hptDiffvsEta->ProjectionY(buff,PtbinMin,PtbinMax);

		maxValue = PtDiff->GetMaximum();
		PtDiff->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
		if(k>0)
		{
			PtDiff->SetLineColor(KGreen-9+k)
				PtDiff->Draw("same");
		}
		else
		{
			PtDiff->SetLineColor(KGreen-9)
				PtDiff->Draw();
		}
		sprintf(buff,"RC:%2.1f<#eta<%2.1f GeV/c",EtaMin,EtaMax);
		drawLatex(0.15,0.90-k*0.05,buff,42,0.05,1);

		EtaMin = EtaMin + 0.1;
		EtaMax = EtaMax + 0.1;
		k++;
	}
	sprintf(buff,"pT diff");
	drawLatex(0.6,0.85,buff,62,0.06,1);
	pdfAction(c1,ps);

	hetaDiffvsEta->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->cd();
	c1->SetLogy(1);
	k = 0;
	maxValue =-999;
	EtaMin =-1.;
	EtaMax =-0.9;
	TH1F *EtaDiff;
	for(Int_t i=0;i<10;i++)
	{
		int PtbinMin = hetaDiffvsEta->GetXaxis()->FindBin(EtaMin);
		int PtbinMax = hetaDiffvsEta->GetXaxis()->FindBin(EtaMax);
		sprintf(buff,"EtaDiff%d",k);
		EtaDiff = (TH1F*)hetaDiffvsEta->ProjectionY(buff,PtbinMin,PtbinMax);

		maxValue = EtaDiff->GetMaximum();
		EtaDiff->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
		if(k>0)
		{
			EtaDiff->SetLineColor(KGreen-9+k)
			EtaDiff->Draw("same");
		}
		else
		{
			EtaDiff->SetLineColor(KGreen-9)
			EtaDiff->Draw();
		}
		sprintf(buff,"RC:%2.1f<#eta<%2.1f GeV/c",EtaMin,EtaMax);
		drawLatex(0.15,0.90-k*0.05,buff,42,0.05,1);

		EtaMin = EtaMin + 0.1;
		EtaMax = EtaMax + 0.1;
		k++;
	}
	sprintf(buff,"Eta diff");
	drawLatex(0.6,0.85,buff,62,0.06,1);
	pdfAction(c1,ps);
	
	
	hphiDiffvsEta->Draw("colz");
	pdfAction(c1,ps);
	c1->Clear();
	c1->cd();
	c1->SetLogy(1);
	k = 0;
	maxValue =-999;
	EtaMin =-1.;
	EtaMax =-0.9;
	TH1F *PhiDiff;
	for(Int_t i=0;i<10;i++)
	{
		int PtbinMin = hphiDiffvsEta->GetXaxis()->FindBin(EtaMin);
		int PtbinMax = hphiDiffvsEta->GetXaxis()->FindBin(EtaMax);
		sprintf(buff,"PhiDiff%d",k);
		PhiDiff = (TH1F*)hphiDiffvsEta->ProjectionY(buff,PtbinMin,PtbinMax);

		maxValue = PhiDiff->GetMaximum();
		PhiDiff->GetYaxis()->SetRangeUser(0.,MCcounts*1.3);
		if(k>0)
		{
			PhiDiff->SetLineColor(KGreen-9+k)
			PhiDiff->Draw("same");
		}
		else
		{
			PhiDiff->SetLineColor(KGreen-9)
			PhiDiff->Draw();
		}
		sprintf(buff,"RC:%2.1f<#eta<%2.1f GeV/c",EtaMin,EtaMax);
		drawLatex(0.15,0.90-k*0.05,buff,42,0.05,1);

		EtaMin = EtaMin + 0.1;
		EtaMax = EtaMax + 0.1;
		k++;
	}
	sprintf(buff,"Eta diff");
	drawLatex(0.6,0.85,buff,62,0.06,1);
	pdfAction(c1,ps);


	TH3F *hElectronMc = (TH3F*)f->Get("hElectronMc");
	TH3F *hHt2ElectronTrg = (TH3F*)f->Get("hHt2ElectronTrg");
	TH3F *hHt2ElectronAdc0 = (TH3F*)f->Get("hHt2ElectronAdc0");
	TH3F *hHt2ElectronPE = (TH3F*)f->Get("hHt2ElectronPE");

	TH1F *hMcPt = (TH1F*)hElectronMc->ProjectionZ("inputMcPt",1,40,1,120);
	TH1F *hRcPt = (TH1F*)hRcElectronPtDiff->ProjectionX("RcPt",1,40,0,1000);
	TH1F *hTrg = (TH1F*)hHt2ElectronTrg->ProjectionZ("trg",1,40,1,120);
	TH1F *hAdc0 = (TH1F*)hHt2ElectronAdc0->ProjectionZ("hAdc0",1,40,1,120);
	TH1F *hPE = (TH1F*)hHt2ElectronPE->ProjectionZ("hPE",1,40,1,120);
	
	hMcPt->Sumw2();
	hRcPt->Sumw2();
	hTrg->Sumw2();
	hAdc0->Sumw2();
	hPE->Sumw2();

	hRcPt->Divide(hRcPt,hMcPt,1,1);
	hTrg->Divide(hTrg,hMcPt,1,1);
	hAdc0->Divide(hAdc0,hMcPt,1,1);
	hPE->Divide(hPE,hMcPt,1,1);

	hRcPt->SetLineColor(1);
	hTrg->SetLineColor(2);
	hAdc0->SetLineColor(3);
	hRcPt->SetLineColor(4);

	hRcPt->GetXaxis()->SetTitle("p_{T}^{mc}");
	hRcPt->GetYaxis()->SetTitle("Efficiecy * Accptance");
	hRcPt->GetXaxis()->SetRangeUser(0,15);
	hRcPt->GetYaxis()->SetRangeUser(0,1.);
	hRcPt->Draw();
	hTrg->Draw("same");
	hAdc0->Draw("same");
	hPE->Draw("same");
	sprintf(buff,"Rc efficiency");
	drawLatex(0.2,0.85,buff,62,0.06,1);
	sprintf(buff,"+ Trg(dsmadc>18)");
	drawLatex(0.2,0.75,buff,62,0.06,2);
	sprintf(buff,"+ adc0(adc0>290)");
	drawLatex(0.2,0.65,buff,62,0.06,3);
	sprintf(buff,"+ pE(0.3<pE<1.5)");
	drawLatex(0.2,0.55,buff,62,0.06,4);
	pdfAction(c1,ps);
	
	ps->On();
	ps->Close();

}	
