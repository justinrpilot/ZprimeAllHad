#include <cstdlib>
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <ostream>
#include "TChain.h"
#include "TLatex.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include <vector>	
#include <algorithm>

#include "names.C"

using namespace std;
using namespace names;

float ptMap(float pt);
float bMap(float b);
float tauMap(float t);
void convertSyst(TH1F *upH, TH1F *dnH, TH1F *nomH, TH1F *targetH);
void symmetrize( TH1F *upH, TH1F *dnH, TH1F *nomH);


int makeTemplates(bool ZPNflag, bool ZPWflag, bool RSGflag){



	TString labels[names::NUM_PROCS];
	labels[names::DATA] = "data";
	labels[names::QCD]  = "qcd";
	labels[names::TT7] = "ttbar7";
	labels[names::TT7_BTAGUP] = "ttbar7_btag_up";
	labels[names::TT7_BTAGDN] = "ttbar7_btag_dn";
	labels[names::TTLOW] = "ttbarLOW";
	labels[names::TT7_SCALEUP] = "ttbar7_jes_up";
	labels[names::TT7_SCALEDN] = "ttbar7_jes_dn";
	labels[names::TT7_PUUP] = "ttbar7_pu_up";
	labels[names::TT7_PUDN] = "ttbar7_pu_dn";
	labels[names::TT7_JERUP] = "ttbar7_jer_up";
	labels[names::TT7_JERDN] = "ttbar7_jer_dn";
	labels[names::TT7_Q2UP] = "ttbar7_q2_up";
	labels[names::TT7_Q2DN] = "ttbar7_q2_dn";
	labels[names::TT7_TOPPTUP] = "ttbar7_toppt_up";
	labels[names::TT7_TOPPTDN] = "ttbar7_toppt_dn";
	labels[names::TT7_PDFUP] = "ttbar7_pdf_up";
	labels[names::TT7_PDFDN] = "ttbar7_pdf_dn";
	labels[names::TT7_PDFNOM] = "ttbar7_pdf_nom";
	labels[names::TT7_MIS] = "ttbar7_mis";
	labels[names::TT10] = "ttbar10";
	labels[names::TT10_BTAGUP] = "ttbar10_btag_up";
	labels[names::TT10_BTAGDN] = "ttbar10_btag_dn";
	labels[names::TT10_SCALEUP] = "ttbar10_jes_up";
	labels[names::TT10_SCALEDN] = "ttbar10_jes_dn";
	labels[names::TT10_JERUP] = "ttbar10_jer_up";
	labels[names::TT10_JERDN] = "ttbar10_jer_dn";
	labels[names::TT10_Q2UP] = "ttbar10_q2_up";
	labels[names::TT10_Q2DN] = "ttbar10_q2_dn";
	labels[names::TT10_PUUP] = "ttbar10_pu_up";
	labels[names::TT10_PUDN] = "ttbar10_pu_dn";
	labels[names::TT10_TOPPTUP] = "ttbar10_toppt_up";
	labels[names::TT10_TOPPTDN] = "ttbar10_toppt_dn";
	labels[names::TT10_PDFUP] = "ttbar10_pdf_up";
	labels[names::TT10_PDFDN] = "ttbar10_pdf_dn";
	labels[names::TT10_PDFNOM] = "ttbar10_pdf_nom";
	labels[names::TT10_MIS] = "ttbar10_mis";
	labels[names::ZPN10_MIS] = "ZpN10_mis";
	labels[names::ZPN15_MIS] = "ZpN15_mis";
	labels[names::ZPN20_MIS] = "ZpN20_mis";
	labels[names::ZPN30_MIS] = "ZpN30_mis";
	labels[names::ZPN10] = "ZpN10";
	labels[names::ZPN12] = "ZpN12";
	labels[names::ZPN15] = "ZpN15";
	labels[names::ZPN20] = "ZpN20";
	labels[names::ZPN30] = "ZpN30";
	labels[names::ZPN10_BTAGUP] = "ZpN10_btag_up";
	labels[names::ZPN12_BTAGUP] = "ZpN12_btag_up";
	labels[names::ZPN15_BTAGUP] = "ZpN15_btag_up";
	labels[names::ZPN20_BTAGUP] = "ZpN20_btag_up";
	labels[names::ZPN30_BTAGUP] = "ZpN30_btag_up";
	labels[names::ZPN10_BTAGDN] = "ZpN10_btag_dn";
	labels[names::ZPN12_BTAGDN] = "ZpN12_btag_dn";
	labels[names::ZPN15_BTAGDN] = "ZpN15_btag_dn";
	labels[names::ZPN20_BTAGDN] = "ZpN20_btag_dn";
	labels[names::ZPN30_BTAGDN] = "ZpN30_btag_dn";
	labels[names::ZPN10_SCALEUP] = "ZpN10_jes_up";
	labels[names::ZPN12_SCALEUP] = "ZpN12_jes_up";
	labels[names::ZPN15_SCALEUP] = "ZpN15_jes_up";
	labels[names::ZPN20_SCALEUP] = "ZpN20_jes_up";
	labels[names::ZPN30_SCALEUP] = "ZpN30_jes_up";
	labels[names::ZPN10_SCALEDN] = "ZpN10_jes_dn";
	labels[names::ZPN12_SCALEDN] = "ZpN12_jes_dn";
	labels[names::ZPN15_SCALEDN] = "ZpN15_jes_dn";
	labels[names::ZPN20_SCALEDN] = "ZpN20_jes_dn";
	labels[names::ZPN30_SCALEDN] = "ZpN30_jes_dn";
	labels[names::ZPN10_PUUP] = "ZpN10_pileup_up";
	labels[names::ZPN12_PUUP] = "ZpN12_pileup_up";
	labels[names::ZPN15_PUUP] = "ZpN15_pileup_up";
	labels[names::ZPN20_PUUP] = "ZpN20_pileup_up";
	labels[names::ZPN30_PUUP] = "ZpN30_pileup_up";
	labels[names::ZPN10_PUDN] = "ZpN10_pileup_dn";
	labels[names::ZPN12_PUDN] = "ZpN12_pileup_dn";
	labels[names::ZPN15_PUDN] = "ZpN15_pileup_dn";
	labels[names::ZPN20_PUDN] = "ZpN20_pileup_dn";
	labels[names::ZPN30_PUDN] = "ZpN30_pileup_dn";
	labels[names::ZPN10_PDFUP] = "ZpN10_pdf_up";
	labels[names::ZPN12_PDFUP] = "ZpN12_pdf_up";
	labels[names::ZPN15_PDFUP] = "ZpN15_pdf_up";
	labels[names::ZPN20_PDFUP] = "ZpN20_pdf_up";
	labels[names::ZPN30_PDFUP] = "ZpN30_pdf_up";
	labels[names::ZPN10_PDFDN] = "ZpN10_pdf_dn";
	labels[names::ZPN12_PDFDN] = "ZpN12_pdf_dn";
	labels[names::ZPN15_PDFDN] = "ZpN15_pdf_dn";
	labels[names::ZPN20_PDFDN] = "ZpN20_pdf_dn";
	labels[names::ZPN30_PDFDN] = "ZpN30_pdf_dn";
	labels[names::ZPN10_PDFNOM] = "ZpN10_pdf_nom";
	labels[names::ZPN12_PDFNOM] = "ZpN12_pdf_nom";
	labels[names::ZPN15_PDFNOM] = "ZpN15_pdf_nom";
	labels[names::ZPN20_PDFNOM] = "ZpN20_pdf_nom";
	labels[names::ZPN30_PDFNOM] = "ZpN30_pdf_nom";
	labels[names::ZPN10_JERUP] = "ZpN10_jer_up";
	labels[names::ZPN12_JERUP] = "ZpN12_jer_up";
	labels[names::ZPN15_JERUP] = "ZpN15_jer_up";
	labels[names::ZPN20_JERUP] = "ZpN20_jer_up";
	labels[names::ZPN30_JERUP] = "ZpN30_jer_up";
	labels[names::ZPN10_JERDN] = "ZpN10_jer_dn";
	labels[names::ZPN12_JERDN] = "ZpN12_jer_dn";
	labels[names::ZPN15_JERDN] = "ZpN15_jer_dn";
	labels[names::ZPN20_JERDN] = "ZpN20_jer_dn";
	labels[names::ZPN30_JERDN] = "ZpN30_jer_dn";
	labels[names::ZPW10] = "ZpW10";
	labels[names::ZPW12] = "ZpW12";
	labels[names::ZPW15] = "ZpW15";
	labels[names::ZPW20] = "ZpW20";
	labels[names::ZPW30] = "ZpW30";
	labels[names::ZPW10_BTAGUP] = "ZpW10_btag_up";
	labels[names::ZPW12_BTAGUP] = "ZpW12_btag_up";
	labels[names::ZPW15_BTAGUP] = "ZpW15_btag_up";
	labels[names::ZPW20_BTAGUP] = "ZpW20_btag_up";
	labels[names::ZPW30_BTAGUP] = "ZpW30_btag_up";
	labels[names::ZPW10_BTAGDN] = "ZpW10_btag_dn";
	labels[names::ZPW12_BTAGDN] = "ZpW12_btag_dn";
	labels[names::ZPW15_BTAGDN] = "ZpW15_btag_dn";
	labels[names::ZPW20_BTAGDN] = "ZpW20_btag_dn";
	labels[names::ZPW30_BTAGDN] = "ZpW30_btag_dn";
	labels[names::ZPW10_SCALEUP] = "ZpW10_jes_up";
	labels[names::ZPW12_SCALEUP] = "ZpW12_jes_up";
	labels[names::ZPW15_SCALEUP] = "ZpW15_jes_up";
	labels[names::ZPW20_SCALEUP] = "ZpW20_jes_up";
	labels[names::ZPW30_SCALEUP] = "ZpW30_jes_up";
	labels[names::ZPW10_SCALEDN] = "ZpW10_jes_dn";
	labels[names::ZPW12_SCALEDN] = "ZpW12_jes_dn";
	labels[names::ZPW15_SCALEDN] = "ZpW15_jes_dn";
	labels[names::ZPW20_SCALEDN] = "ZpW20_jes_dn";
	labels[names::ZPW30_SCALEDN] = "ZpW30_jes_dn";
	labels[names::ZPW10_PUUP] = "ZpW10_pileup_up";
	labels[names::ZPW12_PUUP] = "ZpW12_pileup_up";
	labels[names::ZPW15_PUUP] = "ZpW15_pileup_up";
	labels[names::ZPW20_PUUP] = "ZpW20_pileup_up";
	labels[names::ZPW30_PUUP] = "ZpW30_pileup_up";
	labels[names::ZPW10_PUDN] = "ZpW10_pileup_dn";
	labels[names::ZPW12_PUDN] = "ZpW12_pileup_dn";
	labels[names::ZPW15_PUDN] = "ZpW15_pileup_dn";
	labels[names::ZPW20_PUDN] = "ZpW20_pileup_dn";
	labels[names::ZPW30_PUDN] = "ZpW30_pileup_dn";
	labels[names::ZPW10_PDFUP] = "ZpW10_pdf_up";
	labels[names::ZPW12_PDFUP] = "ZpW12_pdf_up";
	labels[names::ZPW15_PDFUP] = "ZpW15_pdf_up";
	labels[names::ZPW20_PDFUP] = "ZpW20_pdf_up";
	labels[names::ZPW30_PDFUP] = "ZpW30_pdf_up";
	labels[names::ZPW10_PDFDN] = "ZpW10_pdf_dn";
	labels[names::ZPW12_PDFDN] = "ZpW12_pdf_dn";
	labels[names::ZPW15_PDFDN] = "ZpW15_pdf_dn";
	labels[names::ZPW20_PDFDN] = "ZpW20_pdf_dn";
	labels[names::ZPW30_PDFDN] = "ZpW30_pdf_dn";
	labels[names::ZPW10_PDFNOM] = "ZpW10_pdf_nom";
	labels[names::ZPW12_PDFNOM] = "ZpW12_pdf_nom";
	labels[names::ZPW15_PDFNOM] = "ZpW15_pdf_nom";
	labels[names::ZPW20_PDFNOM] = "ZpW20_pdf_nom";
	labels[names::ZPW30_PDFNOM] = "ZpW30_pdf_nom";
	labels[names::ZPW10_JERUP] = "ZpW10_jer_up";
	labels[names::ZPW12_JERUP] = "ZpW12_jer_up";
	labels[names::ZPW15_JERUP] = "ZpW15_jer_up";
	labels[names::ZPW20_JERUP] = "ZpW20_jer_up";
	labels[names::ZPW30_JERUP] = "ZpW30_jer_up";
	labels[names::ZPW10_JERDN] = "ZpW10_jer_dn";
	labels[names::ZPW12_JERDN] = "ZpW12_jer_dn";
	labels[names::ZPW15_JERDN] = "ZpW15_jer_dn";
	labels[names::ZPW20_JERDN] = "ZpW20_jer_dn";
	labels[names::ZPW30_JERDN] = "ZpW30_jer_dn";
	labels[names::RSG10] = "RSG10";
	labels[names::RSG15] = "RSG15";
	labels[names::RSG20] = "RSG20";
	labels[names::RSG30] = "RSG30";
	labels[names::RSG14] = "RSG14";
	labels[names::RSG18] = "RSG18";
	labels[names::RSG25] = "RSG25";
	labels[names::RSG10_BTAGUP] = "RSG10_btag_up";
	labels[names::RSG15_BTAGUP] = "RSG15_btag_up";
	labels[names::RSG20_BTAGUP] = "RSG20_btag_up";
	labels[names::RSG30_BTAGUP] = "RSG30_btag_up";
	labels[names::RSG14_BTAGUP] = "RSG14_btag_up";
	labels[names::RSG18_BTAGUP] = "RSG18_btag_up";
	labels[names::RSG25_BTAGUP] = "RSG25_btag_up";
	labels[names::RSG10_BTAGDN] = "RSG10_btag_dn";
	labels[names::RSG15_BTAGDN] = "RSG15_btag_dn";
	labels[names::RSG20_BTAGDN] = "RSG20_btag_dn";
	labels[names::RSG30_BTAGDN] = "RSG30_btag_dn";
	labels[names::RSG14_BTAGDN] = "RSG14_btag_dn";
	labels[names::RSG18_BTAGDN] = "RSG18_btag_dn";
	labels[names::RSG25_BTAGDN] = "RSG25_btag_dn";
	labels[names::RSG10_PDFUP] = "RSG10_pdf_up";
	labels[names::RSG15_PDFUP] = "RSG15_pdf_up";
	labels[names::RSG20_PDFUP] = "RSG20_pdf_up";
	labels[names::RSG30_PDFUP] = "RSG30_pdf_up";
	labels[names::RSG14_PDFUP] = "RSG14_pdf_up";
	labels[names::RSG18_PDFUP] = "RSG18_pdf_up";
	labels[names::RSG25_PDFUP] = "RSG25_pdf_up";
	labels[names::RSG10_PDFDN] = "RSG10_pdf_dn";
	labels[names::RSG15_PDFDN] = "RSG15_pdf_dn";
	labels[names::RSG20_PDFDN] = "RSG20_pdf_dn";
	labels[names::RSG30_PDFDN] = "RSG30_pdf_dn";
	labels[names::RSG14_PDFDN] = "RSG14_pdf_dn";
	labels[names::RSG18_PDFDN] = "RSG18_pdf_dn";
	labels[names::RSG25_PDFDN] = "RSG25_pdf_dn";
	labels[names::RSG10_PDFNOM] = "RSG10_pdf_nom";
	labels[names::RSG15_PDFNOM] = "RSG15_pdf_nom";
	labels[names::RSG20_PDFNOM] = "RSG20_pdf_nom";
	labels[names::RSG30_PDFNOM] = "RSG30_pdf_nom";
	labels[names::RSG14_PDFNOM] = "RSG14_pdf_nom";
	labels[names::RSG18_PDFNOM] = "RSG18_pdf_nom";
	labels[names::RSG25_PDFNOM] = "RSG25_pdf_nom";
	labels[names::RSG10_SCALEUP] = "RSG10_jes_up";
	labels[names::RSG15_SCALEUP] = "RSG15_jes_up";
	labels[names::RSG20_SCALEUP] = "RSG20_jes_up";
	labels[names::RSG30_SCALEUP] = "RSG30_jes_up";
	labels[names::RSG14_SCALEUP] = "RSG14_jes_up";
	labels[names::RSG18_SCALEUP] = "RSG18_jes_up";
	labels[names::RSG25_SCALEUP] = "RSG25_jes_up";
	labels[names::RSG10_SCALEDN] = "RSG10_jes_dn";
	labels[names::RSG15_SCALEDN] = "RSG15_jes_dn";
	labels[names::RSG20_SCALEDN] = "RSG20_jes_dn";
	labels[names::RSG30_SCALEDN] = "RSG30_jes_dn";
	labels[names::RSG14_SCALEDN] = "RSG14_jes_dn";
	labels[names::RSG18_SCALEDN] = "RSG18_jes_dn";
	labels[names::RSG25_SCALEDN] = "RSG25_jes_dn";
	labels[names::RSG10_PUUP] = "RSG10_pileup_up";
	labels[names::RSG15_PUUP] = "RSG15_pileup_up";
	labels[names::RSG20_PUUP] = "RSG20_pileup_up";
	labels[names::RSG30_PUUP] = "RSG30_pileup_up";
	labels[names::RSG14_PUUP] = "RSG14_pileup_up";
	labels[names::RSG18_PUUP] = "RSG18_pileup_up";
	labels[names::RSG25_PUUP] = "RSG25_pileup_up";
	labels[names::RSG10_PUDN] = "RSG10_pileup_dn";
	labels[names::RSG15_PUDN] = "RSG15_pileup_dn";
	labels[names::RSG20_PUDN] = "RSG20_pileup_dn";
	labels[names::RSG30_PUDN] = "RSG30_pileup_dn";
	labels[names::RSG14_PUDN] = "RSG14_pileup_dn";
	labels[names::RSG18_PUDN] = "RSG18_pileup_dn";
	labels[names::RSG25_PUDN] = "RSG25_pileup_dn";
	labels[names::RSG10_JERUP] = "RSG10_jer_up";
	labels[names::RSG15_JERUP] = "RSG15_jer_up";
	labels[names::RSG20_JERUP] = "RSG20_jer_up";
	labels[names::RSG30_JERUP] = "RSG30_jer_up";
	labels[names::RSG14_JERUP] = "RSG14_jer_up";
	labels[names::RSG18_JERUP] = "RSG18_jer_up";
	labels[names::RSG25_JERUP] = "RSG25_jer_up";
	labels[names::RSG10_JERDN] = "RSG10_jer_dn";
	labels[names::RSG15_JERDN] = "RSG15_jer_dn";
	labels[names::RSG20_JERDN] = "RSG20_jer_dn";
	labels[names::RSG30_JERDN] = "RSG30_jer_dn";
	labels[names::RSG14_JERDN] = "RSG14_jer_dn";
	labels[names::RSG18_JERDN] = "RSG18_jer_dn";
	labels[names::RSG25_JERDN] = "RSG25_jer_dn";

	ofstream txtfile("bins.txt");	
	TFile *outFile = new TFile("bins.root", "RECREATE");


	int nBins = 50;
	int min = 50;
	int max = 5000;
	int nTagCats = 7;

	TString tagLabels[nTagCats];
	tagLabels[0] = "btag0";
	tagLabels[1] = "btag1";
	tagLabels[2] = "btag2";
	tagLabels[3] = "btag0_hiDY";	
	tagLabels[4] = "btag1_hiDY";	
	tagLabels[5] = "btag2_hiDY";	
	tagLabels[6] = "all";	


	TH1F *histos[names::NUM_PROCS][nTagCats];
	for (int h = 0; h < names::NUM_PROCS; h++){
		histos[h][0] = new TH1F(labels[h]+tagLabels[0], labels[h]+tagLabels[0], nBins, min, max);
		histos[h][0]->Sumw2();
		histos[h][1] = new TH1F(labels[h]+tagLabels[1], labels[h]+tagLabels[1], nBins, min, max);
		histos[h][1]->Sumw2();
		histos[h][2] = new TH1F(labels[h]+tagLabels[2], labels[h]+tagLabels[2], nBins, min, max);
		histos[h][2]->Sumw2();
		histos[h][3] = new TH1F(labels[h]+tagLabels[3], labels[h]+tagLabels[3], nBins, min, max);
		histos[h][3]->Sumw2();
		histos[h][4] = new TH1F(labels[h]+tagLabels[4], labels[h]+tagLabels[4], nBins, min, max);
		histos[h][4]->Sumw2();
		histos[h][5] = new TH1F(labels[h]+tagLabels[5], labels[h]+tagLabels[5], nBins, min, max);
		histos[h][5]->Sumw2();
		histos[h][6] = new TH1F(labels[h]+tagLabels[6], labels[h]+tagLabels[6], nBins, min, max);
		histos[h][6]->Sumw2();
	}

	
	
	

	TFile *mistagFile = new TFile("Apr17_mistag.root", "READ");
	TH3D *mistagH = (TH3D *) mistagFile->Get("MISTAG_RATE_SUB_TTBAR_Inclusive");

	TH2F *binFillHist = new TH2F("binFill", "binFill", 50, 0, 50, 2*29*4*9, 0, 2*29*4*9);

	int mistagBins[3];
	mistagBins[0] = mistagH->GetNbinsX();
	mistagBins[1] = mistagH->GetNbinsY();
	mistagBins[2] = mistagH->GetNbinsZ();

	cout << nBins << " " << mistagBins[0] << " " << mistagBins[1] << " " << mistagBins[2] << endl;
	int arrayForErrors[nTagCats][nBins][mistagBins[0]][mistagBins[1]][mistagBins[2]];


	for (int x = 0; x < nBins; x++){
		for (int y = 0; y < mistagBins[0]; y++){
			for (int z = 0; z < mistagBins[1]; z++){
				for (int t = 0; t < mistagBins[2]; t++){
					arrayForErrors[0][x][y][z][t] = 0;
					arrayForErrors[1][x][y][z][t] = 0;
					arrayForErrors[2][x][y][z][t] = 0;
					arrayForErrors[3][x][y][z][t] = 0;
					arrayForErrors[4][x][y][z][t] = 0;
					arrayForErrors[5][x][y][z][t] = 0;
					arrayForErrors[6][x][y][z][t] = 0;
				}
			}
		}
	}


	TString files[names::NUM_PROCS];
	files[names::DATA] = "FINAL/data.root";
	files[names::QCD] = files[names::DATA];
	files[names::TT7] = "FINAL/ttjets7.root";
	files[names::TT7_BTAGUP] = "FINAL/ttjets7.root";
	files[names::TT7_BTAGDN] = "FINAL/ttjets7.root";
	files[names::TT7_TOPPTUP] = "FINAL/ttjets7.root";
	files[names::TT7_TOPPTDN] = "FINAL/ttjets7.root";
	files[names::TTLOW] = "Aug19_ttjetsFULL_sec0_ttpair__TriggernoWeight.root";
	files[names::TT10] = "FINAL/ttjets10.root";
	files[names::TT10_BTAGUP] = "FINAL/ttjets10.root";
	files[names::TT10_BTAGDN] = "FINAL/ttjets10.root";
	files[names::TT10_TOPPTUP] = "FINAL/ttjets10.root";
	files[names::TT10_TOPPTDN] = "FINAL/ttjets10.root";
	files[names::TT7_SCALEUP] = "FINAL/ttjets7_scaleup.root";
	files[names::TT10_SCALEUP] = "FINAL/ttjets10_scaleup.root";
	files[names::TT7_PUUP] = "FINAL/ttjets7.root";
	files[names::TT10_PUUP] = "FINAL/ttjets10.root";
	files[names::TT7_JERUP] = "FINAL/ttjets7_ptsmearup.root";
	files[names::TT10_JERUP] = "FINAL/ttjets10_ptsmearup.root";
	files[names::TT7_Q2UP] = "FINAL/ttjets7_q2up.root";
	files[names::TT7_PDFUP] = "FINAL/ttjets7_pdf.root";
	files[names::TT7_PDFNOM] = "FINAL/ttjets7_pdf.root";
	files[names::TT10_Q2UP] = "FINAL/ttjets10_q2up.root";
	files[names::TT10_PDFUP] = "FINAL/ttjets10_pdf.root";
	files[names::TT10_PDFNOM] = "FINAL/ttjets10_pdf.root";
	files[names::TT7_SCALEDN] = "FINAL/ttjets7_scaledown.root";
	files[names::TT10_SCALEDN] = "FINAL/ttjets10_scaledown.root";
	files[names::TT7_JERDN] = "FINAL/ttjets7_ptsmeardown.root";
	files[names::TT10_JERDN] = "FINAL/ttjets10_ptsmeardown.root";
	files[names::TT7_Q2DN] = "FINAL/ttjets7_q2down.root";
	files[names::TT10_Q2DN] = "FINAL/ttjets10_q2down.root";
	files[names::TT7_PUDN] = "FINAL/ttjets7.root";
	files[names::TT10_PUDN] = "FINAL/ttjets10.root";
	files[names::TT10_PDFDN] = "FINAL/ttjets10_pdf.root";
	files[names::TT7_PDFDN] = "FINAL/ttjets7_pdf.root";
	files[names::TT7_MIS] = files[names::TT7];
	files[names::TT10_MIS] = files[names::TT10];
	files[names::ZPN20_MIS] = "Feb12_Zp20_mistag.root";
	files[names::ZPN10_MIS] = "Feb12_Zp10_mistag.root";
	files[names::ZPN15_MIS] = "Feb12_Zp15_mistag.root";
	files[names::ZPN30_MIS] = "Feb12_Zp30_mistag.root";
	files[names::ZPN20] = "FINAL/Zp20.root";
	files[names::ZPN10] = "FINAL/Zp10.root";
	files[names::ZPN12] = "FINAL/Zp12.root";
	files[names::ZPN15] = "FINAL/Zp15.root";
	files[names::ZPN30] = "FINAL/Zp30.root";
	files[names::ZPN20_BTAGUP] = "FINAL/Zp20.root";
	files[names::ZPN10_BTAGUP] = "FINAL/Zp10.root";
	files[names::ZPN12_BTAGUP] = "FINAL/Zp12.root";
	files[names::ZPN15_BTAGUP] = "FINAL/Zp15.root";
	files[names::ZPN30_BTAGUP] = "FINAL/Zp30.root";
	files[names::ZPN20_BTAGDN] = "FINAL/Zp20.root";
	files[names::ZPN10_BTAGDN] = "FINAL/Zp10.root";
	files[names::ZPN12_BTAGDN] = "FINAL/Zp12.root";
	files[names::ZPN15_BTAGDN] = "FINAL/Zp15.root";
	files[names::ZPN30_BTAGDN] = "FINAL/Zp30.root";
	files[names::ZPN20_PUUP] = "FINAL/Zp20.root";
	files[names::ZPN10_PUUP] = "FINAL/Zp10.root";
	files[names::ZPN12_PUUP] = "FINAL/Zp12.root";
	files[names::ZPN15_PUUP] = "FINAL/Zp15.root";
	files[names::ZPN30_PUUP] = "FINAL/Zp30.root";
	files[names::ZPN20_PUDN] = "FINAL/Zp20.root";
	files[names::ZPN10_PUDN] = "FINAL/Zp10.root";
	files[names::ZPN12_PUDN] = "FINAL/Zp12.root";
	files[names::ZPN15_PUDN] = "FINAL/Zp15.root";
	files[names::ZPN30_PUDN] = "FINAL/Zp30.root";
	files[names::ZPN20_SCALEUP] = "FINAL/Zp20_scaleup.root";
	files[names::ZPN10_SCALEUP] = "FINAL/Zp10_scaleup.root";
	files[names::ZPN12_SCALEUP] = "FINAL/Zp12_scaleup.root";
	files[names::ZPN15_SCALEUP] = "FINAL/Zp15_scaleup.root";
	files[names::ZPN30_SCALEUP] = "FINAL/Zp30_scaleup.root";
	files[names::ZPN20_PDFUP] = "FINAL/Zp20_pdf.root";
	files[names::ZPN10_PDFUP] = "FINAL/Zp10_pdf.root";
	files[names::ZPN12_PDFUP] = "FINAL/Zp12_pdf.root";
	files[names::ZPN15_PDFUP] = "FINAL/Zp15_pdf.root";
	files[names::ZPN30_PDFUP] = "FINAL/Zp30_pdf.root";
	files[names::ZPN20_PDFNOM] = "FINAL/Zp20_pdf.root";
	files[names::ZPN10_PDFNOM] = "FINAL/Zp10_pdf.root";
	files[names::ZPN12_PDFNOM] = "FINAL/Zp12_pdf.root";
	files[names::ZPN15_PDFNOM] = "FINAL/Zp15_pdf.root";
	files[names::ZPN30_PDFNOM] = "FINAL/Zp30_pdf.root";
	files[names::ZPN20_PDFDN] = "FINAL/Zp20_pdf.root";
	files[names::ZPN10_PDFDN] = "FINAL/Zp10_pdf.root";
	files[names::ZPN12_PDFDN] = "FINAL/Zp12_pdf.root";
	files[names::ZPN15_PDFDN] = "FINAL/Zp15_pdf.root";
	files[names::ZPN30_PDFDN] = "FINAL/Zp30_pdf.root";
	files[names::ZPN20_SCALEDN] = "FINAL/Zp20_scaledown.root";
	files[names::ZPN10_SCALEDN] = "FINAL/Zp10_scaledown.root";
	files[names::ZPN12_SCALEDN] = "FINAL/Zp12_scaledown.root";
	files[names::ZPN15_SCALEDN] = "FINAL/Zp15_scaledown.root";
	files[names::ZPN30_SCALEDN] = "FINAL/Zp30_scaledown.root";
	files[names::ZPN20_JERUP] = "FINAL/Zp20_ptsmearup.root";
	files[names::ZPN10_JERUP] = "FINAL/Zp10_ptsmearup.root";
	files[names::ZPN12_JERUP] = "FINAL/Zp12_ptsmearup.root";
	files[names::ZPN15_JERUP] = "FINAL/Zp15_ptsmearup.root";
	files[names::ZPN30_JERUP] = "FINAL/Zp30_ptsmearup.root";
	files[names::ZPN20_JERDN] = "FINAL/Zp20_ptsmeardown.root";
	files[names::ZPN10_JERDN] = "FINAL/Zp10_ptsmeardown.root";
	files[names::ZPN12_JERDN] = "FINAL/Zp12_ptsmeardown.root";
	files[names::ZPN15_JERDN] = "FINAL/Zp15_ptsmeardown.root";
	files[names::ZPN30_JERDN] = "FINAL/Zp30_ptsmeardown.root";
	files[names::ZPW20] = "FINAL/ZpW20.root";
	files[names::ZPW10] = "FINAL/ZpW10.root";
	files[names::ZPW12] = "FINAL/ZpW12.root";
	files[names::ZPW15] = "FINAL/ZpW15.root";
	files[names::ZPW30] = "FINAL/ZpW30.root";
	files[names::ZPW20_BTAGUP] = "FINAL/ZpW20.root";
	files[names::ZPW10_BTAGUP] = "FINAL/ZpW10.root";
	files[names::ZPW12_BTAGUP] = "FINAL/ZpW12.root";
	files[names::ZPW15_BTAGUP] = "FINAL/ZpW15.root";
	files[names::ZPW30_BTAGUP] = "FINAL/ZpW30.root";
	files[names::ZPW20_BTAGDN] = "FINAL/ZpW20.root";
	files[names::ZPW10_BTAGDN] = "FINAL/ZpW10.root";
	files[names::ZPW12_BTAGDN] = "FINAL/ZpW12.root";
	files[names::ZPW15_BTAGDN] = "FINAL/ZpW15.root";
	files[names::ZPW30_BTAGDN] = "FINAL/ZpW30.root";
	files[names::ZPW20_PUUP] = "FINAL/ZpW20.root";
	files[names::ZPW10_PUUP] = "FINAL/ZpW10.root";
	files[names::ZPW12_PUUP] = "FINAL/ZpW12.root";
	files[names::ZPW15_PUUP] = "FINAL/ZpW15.root";
	files[names::ZPW30_PUUP] = "FINAL/ZpW30.root";
	files[names::ZPW20_PUDN] = "FINAL/ZpW20.root";
	files[names::ZPW10_PUDN] = "FINAL/ZpW10.root";
	files[names::ZPW12_PUDN] = "FINAL/ZpW12.root";
	files[names::ZPW15_PUDN] = "FINAL/ZpW15.root";
	files[names::ZPW30_PUDN] = "FINAL/ZpW30.root";
	files[names::ZPW20_SCALEUP] = "FINAL/ZpW20_scaleup.root";
	files[names::ZPW10_SCALEUP] = "FINAL/ZpW10_scaleup.root";
	files[names::ZPW12_SCALEUP] = "FINAL/ZpW12_scaleup.root";
	files[names::ZPW15_SCALEUP] = "FINAL/ZpW15_scaleup.root";
	files[names::ZPW30_SCALEUP] = "FINAL/ZpW30_scaleup.root";
	files[names::ZPW20_SCALEDN] = "FINAL/ZpW20_scaledown.root";
	files[names::ZPW10_SCALEDN] = "FINAL/ZpW10_scaledown.root";
	files[names::ZPW12_SCALEDN] = "FINAL/ZpW12_scaledown.root";
	files[names::ZPW15_SCALEDN] = "FINAL/ZpW15_scaledown.root";
	files[names::ZPW30_SCALEDN] = "FINAL/ZpW30_scaledown.root";
	files[names::ZPW20_PDFUP] = "FINAL/ZpW20_pdf.root";
	files[names::ZPW10_PDFUP] = "FINAL/ZpW10_pdf.root";
	files[names::ZPW12_PDFUP] = "FINAL/ZpW12_pdf.root";
	files[names::ZPW15_PDFUP] = "FINAL/ZpW15_pdf.root";
	files[names::ZPW30_PDFUP] = "FINAL/ZpW30_pdf.root";
	files[names::ZPW20_PDFNOM] = "FINAL/ZpW20_pdf.root";
	files[names::ZPW10_PDFNOM] = "FINAL/ZpW10_pdf.root";
	files[names::ZPW12_PDFNOM] = "FINAL/ZpW12_pdf.root";
	files[names::ZPW15_PDFNOM] = "FINAL/ZpW15_pdf.root";
	files[names::ZPW30_PDFNOM] = "FINAL/ZpW30_pdf.root";
	files[names::ZPW20_PDFDN] = "FINAL/ZpW20_pdf.root";
	files[names::ZPW10_PDFDN] = "FINAL/ZpW10_pdf.root";
	files[names::ZPW12_PDFDN] = "FINAL/ZpW12_pdf.root";
	files[names::ZPW15_PDFDN] = "FINAL/ZpW15_pdf.root";
	files[names::ZPW30_PDFDN] = "FINAL/ZpW30_pdf.root";
	files[names::ZPW20_JERUP] = "FINAL/ZpW20_ptsmearup.root";
	files[names::ZPW10_JERUP] = "FINAL/ZpW10_ptsmearup.root";
	files[names::ZPW12_JERUP] = "FINAL/ZpW12_ptsmearup.root";
	files[names::ZPW15_JERUP] = "FINAL/ZpW15_ptsmearup.root";
	files[names::ZPW30_JERUP] = "FINAL/ZpW30_ptsmearup.root";
	files[names::ZPW20_JERDN] = "FINAL/ZpW20_ptsmeardown.root";
	files[names::ZPW10_JERDN] = "FINAL/ZpW10_ptsmeardown.root";
	files[names::ZPW12_JERDN] = "FINAL/ZpW12_ptsmeardown.root";
	files[names::ZPW15_JERDN] = "FINAL/ZpW15_ptsmeardown.root";
	files[names::ZPW30_JERDN] = "FINAL/ZpW30_ptsmeardown.root";
	files[names::RSG20] = "FINAL/RSG20.root";
	files[names::RSG10] = "FINAL/RSG10.root";
	files[names::RSG15] = "FINAL/RSG15.root";
	files[names::RSG30] = "FINAL/RSG30.root";
	files[names::RSG14] = "FINAL/RSG14.root";
	files[names::RSG18] = "FINAL/RSG18.root";
	files[names::RSG25] = "FINAL/RSG25.root";
	files[names::RSG20_BTAGUP] = "FINAL/RSG20.root";
	files[names::RSG10_BTAGUP] = "FINAL/RSG10.root";
	files[names::RSG15_BTAGUP] = "FINAL/RSG15.root";
	files[names::RSG30_BTAGUP] = "FINAL/RSG30.root";
	files[names::RSG14_BTAGUP] = "FINAL/RSG14.root";
	files[names::RSG18_BTAGUP] = "FINAL/RSG18.root";
	files[names::RSG25_BTAGUP] = "FINAL/RSG25.root";
	files[names::RSG20_BTAGDN] = "FINAL/RSG20.root";
	files[names::RSG10_BTAGDN] = "FINAL/RSG10.root";
	files[names::RSG15_BTAGDN] = "FINAL/RSG15.root";
	files[names::RSG30_BTAGDN] = "FINAL/RSG30.root";
	files[names::RSG14_BTAGDN] = "FINAL/RSG14.root";
	files[names::RSG18_BTAGDN] = "FINAL/RSG18.root";
	files[names::RSG25_BTAGDN] = "FINAL/RSG25.root";
	files[names::RSG20_PUUP] = "FINAL/RSG20.root";
	files[names::RSG10_PUUP] = "FINAL/RSG10.root";
	files[names::RSG15_PUUP] = "FINAL/RSG15.root";
	files[names::RSG30_PUUP] = "FINAL/RSG30.root";
	files[names::RSG14_PUUP] = "FINAL/RSG14.root";
	files[names::RSG18_PUUP] = "FINAL/RSG18.root";
	files[names::RSG25_PUUP] = "FINAL/RSG25.root";
	files[names::RSG20_PUDN] = "FINAL/RSG20.root";
	files[names::RSG10_PUDN] = "FINAL/RSG10.root";
	files[names::RSG15_PUDN] = "FINAL/RSG15.root";
	files[names::RSG30_PUDN] = "FINAL/RSG30.root";
	files[names::RSG14_PUDN] = "FINAL/RSG14.root";
	files[names::RSG18_PUDN] = "FINAL/RSG18.root";
	files[names::RSG25_PUDN] = "FINAL/RSG25.root";
	files[names::RSG20_SCALEUP] = "FINAL/RSG20_scaleup.root";
	files[names::RSG10_SCALEUP] = "FINAL/RSG10_scaleup.root";
	files[names::RSG15_SCALEUP] = "FINAL/RSG15_scaleup.root";
	files[names::RSG30_SCALEUP] = "FINAL/RSG30_scaleup.root";
	files[names::RSG14_SCALEUP] = "FINAL/RSG14_scaleup.root";
	files[names::RSG18_SCALEUP] = "FINAL/RSG18_scaleup.root";
	files[names::RSG25_SCALEUP] = "FINAL/RSG25_scaleup.root";
	files[names::RSG20_SCALEDN] = "FINAL/RSG20_scaledown.root";
	files[names::RSG10_SCALEDN] = "FINAL/RSG10_scaledown.root";
	files[names::RSG15_SCALEDN] = "FINAL/RSG15_scaledown.root";
	files[names::RSG30_SCALEDN] = "FINAL/RSG30_scaledown.root";
	files[names::RSG14_SCALEDN] = "FINAL/RSG14_scaledown.root";
	files[names::RSG18_SCALEDN] = "FINAL/RSG18_scaledown.root";
	files[names::RSG25_SCALEDN] = "FINAL/RSG25_scaledown.root";
	files[names::RSG20_PDFUP] = "FINAL/RSG20_pdf.root";
	files[names::RSG10_PDFUP] = "FINAL/RSG10_pdf.root";
	files[names::RSG15_PDFUP] = "FINAL/RSG15_pdf.root";
	files[names::RSG30_PDFUP] = "FINAL/RSG30_pdf.root";
	files[names::RSG14_PDFUP] = "FINAL/RSG14_pdf.root";
	files[names::RSG18_PDFUP] = "FINAL/RSG18_pdf.root";
	files[names::RSG25_PDFUP] = "FINAL/RSG25_pdf.root";
	files[names::RSG20_PDFNOM] = "FINAL/RSG20_pdf.root";
	files[names::RSG10_PDFNOM] = "FINAL/RSG10_pdf.root";
	files[names::RSG15_PDFNOM] = "FINAL/RSG15_pdf.root";
	files[names::RSG30_PDFNOM] = "FINAL/RSG30_pdf.root";
	files[names::RSG14_PDFNOM] = "FINAL/RSG14_pdf.root";
	files[names::RSG18_PDFNOM] = "FINAL/RSG18_pdf.root";
	files[names::RSG25_PDFNOM] = "FINAL/RSG25_pdf.root";
	files[names::RSG20_PDFDN] = "FINAL/RSG20_pdf.root";
	files[names::RSG10_PDFDN] = "FINAL/RSG10_pdf.root";
	files[names::RSG15_PDFDN] = "FINAL/RSG15_pdf.root";
	files[names::RSG30_PDFDN] = "FINAL/RSG30_pdf.root";
	files[names::RSG14_PDFDN] = "FINAL/RSG14_pdf.root";
	files[names::RSG18_PDFDN] = "FINAL/RSG18_pdf.root";
	files[names::RSG25_PDFDN] = "FINAL/RSG25_pdf.root";
	files[names::RSG20_JERUP] = "FINAL/RSG20_ptsmearup.root";
	files[names::RSG10_JERUP] = "FINAL/RSG10_ptsmearup.root";
	files[names::RSG15_JERUP] = "FINAL/RSG15_ptsmearup.root";
	files[names::RSG30_JERUP] = "FINAL/RSG30_ptsmearup.root";
	files[names::RSG14_JERUP] = "FINAL/RSG14_ptsmearup.root";
	files[names::RSG18_JERUP] = "FINAL/RSG18_ptsmearup.root";
	files[names::RSG25_JERUP] = "FINAL/RSG25_ptsmearup.root";
	files[names::RSG20_JERDN] = "FINAL/RSG20_ptsmeardown.root";
	files[names::RSG10_JERDN] = "FINAL/RSG10_ptsmeardown.root";
	files[names::RSG15_JERDN] = "FINAL/RSG15_ptsmeardown.root";
	files[names::RSG30_JERDN] = "FINAL/RSG30_ptsmeardown.root";
	files[names::RSG14_JERDN] = "FINAL/RSG14_ptsmeardown.root";
	files[names::RSG18_JERDN] = "FINAL/RSG18_ptsmeardown.root";
	files[names::RSG25_JERDN] = "FINAL/RSG25_ptsmeardown.root";



	ofstream event_list;
	event_list.open("data_events.txt");

	TFile *triggerFile = new TFile("trigger_weights.root");
	TH1F *triggerHist = (TH1F *) triggerFile->Get("triggerHist");

	TFile *pileupFile = new TFile("PileupSystWeights.root");
	TH1F *PileupUpH = (TH1F *) pileupFile->Get("upPileupW");
	TH1F *PileupDnH = (TH1F *) pileupFile->Get("dnPileupW");


	for (int proc = 0; proc < names::NUM_PROCS; proc++){

	TChain *currentChain = new TChain("treeVars");
	currentChain->Add(files[proc]); 
	int nEvents = currentChain->GetEntries();


	float mttMass, mttMassPred, mistagWt, deltaY, jet1Pt, jet2Pt, jet1tau32, jet1Eta, jet2Eta, jet2tau32, jet1BDisc, jet2BDisc, jet1SubjetMaxBDisc, jet2SubjetMaxBDisc, jet1bSF, jet1bSFErrUp, jet1bSFErrDn, jet2bSF, jet2bSFErrUp, jet2bSFErrDn, jetPtForMistag;
	float cutflow, index, ptReweight, pdfWeightUp, pdfWeightDn;
	
	int run, lumi, npv, event;











	currentChain->SetBranchAddress("mttMass", &mttMass);
	currentChain->SetBranchAddress("mttMassPred", &mttMassPred);
	currentChain->SetBranchAddress("cutflow", &cutflow);
	currentChain->SetBranchAddress("index", &index);
	currentChain->SetBranchAddress("mistagWt", &mistagWt);
	currentChain->SetBranchAddress("jet1Pt", &jet1Pt);
	currentChain->SetBranchAddress("jet2Pt", &jet2Pt);
	currentChain->SetBranchAddress("jet1tau32", &jet1tau32);
	currentChain->SetBranchAddress("jet2tau32", &jet2tau32);
	currentChain->SetBranchAddress("jet1SubjetMaxBDisc", &jet1SubjetMaxBDisc);
	currentChain->SetBranchAddress("jet1BDisc", &jet1BDisc);
	currentChain->SetBranchAddress("jet2SubjetMaxBDisc", &jet2SubjetMaxBDisc);
	currentChain->SetBranchAddress("jet2BDisc", &jet2BDisc);
	currentChain->SetBranchAddress("jetPtForMistag", &jetPtForMistag);
	currentChain->SetBranchAddress("ptReweight", &ptReweight);	
	currentChain->SetBranchAddress("deltaY", &deltaY);
	currentChain->SetBranchAddress("run", &run);
	currentChain->SetBranchAddress("event", &event);
	currentChain->SetBranchAddress("lumi", &lumi);
	currentChain->SetBranchAddress("npv", &npv);
	currentChain->SetBranchAddress("jet1bSF", &jet1bSF);
	currentChain->SetBranchAddress("jet1bSFErrUp", &jet1bSFErrUp);
	currentChain->SetBranchAddress("jet1bSFErrDn", &jet1bSFErrDn);
	currentChain->SetBranchAddress("jet2bSF", &jet2bSF);
	currentChain->SetBranchAddress("jet2bSFErrUp", &jet2bSFErrUp);
	currentChain->SetBranchAddress("jet2bSFErrDn", &jet2bSFErrDn);
	currentChain->SetBranchAddress("jet1Eta", &jet1Eta);	
	currentChain->SetBranchAddress("jet2Eta", &jet2Eta);	
	currentChain->SetBranchAddress("pdfWeightUp", &pdfWeightUp);
	currentChain->SetBranchAddress("pdfWeightDn", &pdfWeightDn);

	cout << "Processing sample " << labels[proc] << " with   " << nEvents << " events" << endl;


	int QCDcount = 0;
	for (int i = 0; i < nEvents; i++){

		currentChain->GetEntry(i);


		float upPU = PileupUpH->GetBinContent(npv);
		float dnPU = PileupDnH->GetBinContent(npv); 

		float tau32cut = 0.7;
		float csvCut = 0.679;

		bool btag0 = jet1SubjetMaxBDisc < csvCut && jet2SubjetMaxBDisc < csvCut && jet1tau32 < tau32cut && jet2tau32 < tau32cut;
		bool btag1 = (jet1SubjetMaxBDisc > csvCut || jet2SubjetMaxBDisc > csvCut) && (jet1SubjetMaxBDisc < csvCut || jet2SubjetMaxBDisc <  csvCut) && jet1tau32 < tau32cut && jet2tau32 < tau32cut;
		bool btag2 = jet1SubjetMaxBDisc > csvCut && jet2SubjetMaxBDisc > csvCut && jet1tau32 < tau32cut && jet2tau32 < tau32cut;
		bool passDY = abs(deltaY) < 1.0;
	
		bool passed = (btag0 || btag1 || btag2);
	
		float subSF1, subSF2 = 1.063;
		if (abs(jet1Eta) > 1.0) subSF1 = 0.906;
		if (abs(jet2Eta) > 1.0) subSF2 = 0.906;
		float triggerWt = triggerHist->GetBinContent( triggerHist->FindBin( jet1Pt + jet2Pt ) );
		float subSF = subSF1*subSF2;
		float bScore, tauScore = -999.9;
		subSF = 1.0;

		float w1 = jet1bSF;
		float w2 = jet2bSF;

		if (proc == names::TT7_BTAGUP || proc == names::TT10_BTAGUP || (proc >= names::ZPN10_BTAGUP && proc <= names::ZPN30_BTAGUP) || (proc >= names::ZPW10_BTAGUP && proc <= names::ZPW30_BTAGUP) || (proc >= names::RSG10_BTAGUP && proc <= names::RSG30_BTAGUP)){

			w1 = jet1bSFErrUp;
			w2 = jet1bSFErrUp;
		}
		
		if (proc == names::TT7_BTAGDN || proc == names::TT10_BTAGDN || (proc >= names::ZPN10_BTAGDN && proc <= names::ZPN30_BTAGDN) || (proc >= names::ZPW10_BTAGDN && proc <= names::ZPW30_BTAGDN) || (proc >= names::RSG10_BTAGDN && proc <= names::RSG30_BTAGDN)){
			w1 = jet1bSFErrDn;
			w2 = jet2bSFErrDn;
		}



		float bSF0 = 1;
		float bSF1 = w1*w2;
		float bSF2 = w1*w2;
		float bSF2to1 = (1-w1)*w2 + (1-w2)*w1;
		float bSF2to0 = (1-w1)*(1-w2);
		float bSF1to0 = (1-w1)*w2 + (1-w2)*w1;


		if (proc == names::ZPN20){

		cout << "migrations   " << bSF2to1 << " " << bSF2to0 << " " << bSF1to0 << endl;
		cout << "b SFs    " << bSF0 << " " << bSF1 << " " << bSF2 << endl;

		}

		
		if (proc == names::DATA || proc == names::ZPN20){
			subSF = 1.0;
			triggerWt = 1.0;
			bSF0 = 1.0;
			bSF1 = 1.0;
			bSF2 = 1.0;
			bSF2to1 = 0.0;
			bSF2to0 = 0.0;
			bSF1to0 = 0.0;
		}

		if (proc <= names::TT10 && cutflow == 4 && index == 0 && passed) {

			if (proc == names::DATA && (btag2) && jet1Eta > -1 && jet1Eta < 0.0) event_list << run << " " << lumi << " " << event << endl;

			ptReweight = 1.0;
			
			if (proc >= names::ZPN10_PDFUP && proc <= names::ZPN30_PDFUP) triggerWt = 1.0*pdfWeightUp;
			if (proc >= names::ZPN10_PDFDN && proc <= names::ZPN30_PDFDN) triggerWt = 1.0*pdfWeightDn;
			if (proc >= names::ZPN10_PDFNOM && proc <= names::ZPN30_PDFNOM) triggerWt = 1.0;
			if (proc >= names::ZPW10_PDFUP && proc <= names::ZPW30_PDFUP) triggerWt = 1.0*pdfWeightUp;
			if (proc >= names::ZPW10_PDFDN && proc <= names::ZPW30_PDFDN) triggerWt = 1.0*pdfWeightDn;
			if (proc >= names::ZPW10_PDFNOM && proc <= names::ZPW30_PDFNOM) triggerWt = 1.0;
			if (proc >= names::RSG10_PDFUP && proc <= names::RSG30_PDFUP) triggerWt = 1.0*pdfWeightUp;
			if (proc >= names::RSG10_PDFDN && proc <= names::RSG30_PDFDN) triggerWt = 1.0*pdfWeightDn;
			if (proc >= names::RSG10_PDFNOM && proc <= names::RSG30_PDFNOM) triggerWt = 1.0;
			if (proc >= names::ZPN10_PUUP && proc <= names::ZPN30_PUUP) ptReweight = upPU;
			if (proc >= names::ZPN10_PUDN && proc <= names::ZPN30_PUDN) ptReweight = dnPU;
			if (proc >= names::ZPW10_PUUP && proc <= names::ZPW30_PUUP) ptReweight = upPU;
			if (proc >= names::ZPW10_PUDN && proc <= names::ZPW30_PUDN) ptReweight = dnPU;
			if (proc >= names::RSG10_PUUP && proc <= names::RSG30_PUUP) ptReweight = upPU;
			if (proc >= names::RSG10_PUDN && proc <= names::RSG30_PUDN) ptReweight = dnPU;
			if (proc == names::TT7_PUUP || proc == names::TT10_PUUP) ptReweight = upPU; 
			if (proc == names::TT7_PUDN || proc == names::TT10_PUDN) ptReweight = dnPU; 

			if (proc >= names::TT7_SCALEUP && proc <= names::TT10 ) {

				if (proc == names::TT10_PDFUP || proc == names::TT7_PDFUP) triggerWt = 1.0*pdfWeightUp;
				if (proc == names::TT10_PDFDN || proc == names::TT7_PDFDN) triggerWt = 1.0*pdfWeightDn;
				if (proc == names::TT10_TOPPTUP || proc == names::TT7_TOPPTUP) { ptReweight = ptReweight + (abs(ptReweight - 1.0) / 2.0); cout << "UP " << ptReweight << endl; }
				if (proc == names::TT10_TOPPTDN || proc == names::TT7_TOPPTDN) { ptReweight = ptReweight - (abs(ptReweight - 1.0) / 2.0); cout << "DN " << ptReweight << endl; }





				histos[proc][6]->Fill(mttMass, ptReweight*subSF*triggerWt);
				if (btag0 && passDY) histos[proc][0]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF0);
				if (btag1 && passDY) {
					histos[proc][1]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF1);
					histos[proc][0]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF1to0);
				}
				if (btag2 && passDY) {
					histos[proc][2]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2);
					histos[proc][1]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2to1);
					histos[proc][0]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2to0);
				}
				if (btag0 && !passDY) histos[proc][3]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF0);
				if (btag1 && !passDY) {
					histos[proc][4]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF1);
					histos[proc][3]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF1to0);
				}
				if (btag2 && !passDY) {
					histos[proc][5]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2);
					histos[proc][4]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2to1);
					histos[proc][3]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2to0);
				}
			}
			else {
				histos[proc][6]->Fill(mttMass, ptReweight*subSF*triggerWt);
				if (btag0 && passDY) histos[proc][0]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF0);
				if (btag1 && passDY) {
					histos[proc][1]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF1);
					histos[proc][0]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF1to0);
				}
				if (btag2 && passDY) {
					histos[proc][2]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2);
					histos[proc][1]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2to1);
					histos[proc][0]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2to0);
				}
				if (btag0 && !passDY) histos[proc][3]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF0);
				if (btag1 && !passDY) {
					histos[proc][4]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF1);
					histos[proc][3]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF1to0);
				}
				if (btag2 && !passDY) {
					histos[proc][5]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2);
					histos[proc][4]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2to1);
					histos[proc][3]->Fill(mttMass, ptReweight*subSF*triggerWt*bSF2to0);
				}
			}
		}
		else if (proc > names::TT10 && index == 1 && passed) {

			if (proc == names::TT7_MIS || proc == names::TT10_MIS) {
				histos[proc][6]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt);
				if (btag0 && passDY) histos[proc][0]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF0);
				if (btag1 && passDY){
					 histos[proc][1]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF1);
					 histos[proc][0]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF1to0);
				}
				if (btag2 && passDY) {
					histos[proc][2]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF2);
					histos[proc][1]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF2to1);
					histos[proc][0]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF2to0);
				}
				if (btag0 && !passDY) histos[proc][3]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF0);
				if (btag1 && !passDY) {
					histos[proc][4]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF1);
					histos[proc][3]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF1to0);
				}
				if (btag2 && !passDY) {
					histos[proc][5]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF2);
					histos[proc][4]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF2to1);
					histos[proc][3]->Fill(mttMassPred, mistagWt*ptReweight*triggerWt*bSF2to0);
				}
			}
			else {
				histos[proc][6]->Fill(mttMassPred, mistagWt);
				if (btag0 && passDY) histos[proc][0]->Fill(mttMassPred, mistagWt);
				if (btag1 && passDY) histos[proc][1]->Fill(mttMassPred, mistagWt);
				if (btag2 && passDY) histos[proc][2]->Fill(mttMassPred, mistagWt);
				if (btag0 && !passDY) histos[proc][3]->Fill(mttMassPred, mistagWt);
				if (btag1 && !passDY) histos[proc][4]->Fill(mttMassPred, mistagWt);
				if (btag2 && !passDY) histos[proc][5]->Fill(mttMassPred, mistagWt);
			}	
			if (jet1Pt == jetPtForMistag) {
				bScore = jet1SubjetMaxBDisc;
				tauScore = jet1tau32;
			}
			if (jet2Pt == jetPtForMistag) {
				bScore = jet2SubjetMaxBDisc;
				tauScore = jet2tau32;
			}
			if (proc == names::QCD) {
				int mtt_bin = histos[names::QCD][3]->FindBin(mttMassPred) - 1;
				if (mtt_bin == 50) mtt_bin = 49;
				binFillHist->Fill(mtt_bin, mistagH->FindBin(ptMap(jetPtForMistag), bMap(bScore), tauMap(tauScore)));	
				arrayForErrors[6][ mtt_bin ][ int(floor(ptMap(jetPtForMistag)))][ int(floor(bMap(bScore)))][int(floor(tauMap(tauScore)))] += 1;		
				if (btag0 && passDY) arrayForErrors[0][ mtt_bin ][ int(floor(ptMap(jetPtForMistag)))][ int(floor(bMap(bScore)))][int(floor(tauMap(tauScore)))] += 1;		
				if (btag1 && passDY) arrayForErrors[1][ mtt_bin ][ int(floor(ptMap(jetPtForMistag)))][ int(floor(bMap(bScore)))][int(floor(tauMap(tauScore)))] += 1;		
				if (btag2 && passDY) arrayForErrors[2][ mtt_bin ][ int(floor(ptMap(jetPtForMistag)))][ int(floor(bMap(bScore)))][int(floor(tauMap(tauScore)))] += 1;		
				if (btag0 && !passDY) arrayForErrors[3][ mtt_bin ][ int(floor(ptMap(jetPtForMistag)))][ int(floor(bMap(bScore)))][int(floor(tauMap(tauScore)))] += 1;		
				if (btag1 && !passDY) arrayForErrors[4][ mtt_bin ][ int(floor(ptMap(jetPtForMistag)))][ int(floor(bMap(bScore)))][int(floor(tauMap(tauScore)))] += 1;		
				if (btag2 && !passDY) arrayForErrors[5][ mtt_bin ][ int(floor(ptMap(jetPtForMistag)))][ int(floor(bMap(bScore)))][int(floor(tauMap(tauScore)))] += 1;		
		
			}
		}

	
	}	

	delete currentChain;



	}


	cout << "Done with Filling Histograms" << endl;

	event_list.close();

	//Keep track of b-tag SF cascading
	float btagSF[names::NUM_PROCS][nTagCats];
	float wB = 0.972;
	for (int proc = 0; proc < names::NUM_PROCS; proc++){

			btagSF[proc][0] = 1.0;	
			btagSF[proc][1] = 1.0;	
			btagSF[proc][2] = 1.0;	
			btagSF[proc][3] = 1.0;	
			btagSF[proc][4] = 1.0;	
			btagSF[proc][5] = 1.0;	

	
		//btagSF[proc][2] = wB * wB;
		//btagSF[proc][5] = wB * wB;

		//btagSF[proc][1] = wB + histos[proc][2]->GetEntries() / histos[proc][1]->GetEntries() * 2 * wB * ( 1 - wB );
		//btagSF[proc][4] = wB + histos[proc][5]->GetEntries() / histos[proc][4]->GetEntries() * 2 * wB * ( 1 - wB );
	
		//btagSF[proc][0] = 1 + histos[proc][2]->GetEntries() / histos[proc][0]->GetEntries() * (1-wB) * (1-wB) + histos[proc][1]->GetEntries() / histos[proc][0]->GetEntries() * (1-wB);	
		//btagSF[proc][3] = 1 + histos[proc][5]->GetEntries() / histos[proc][3]->GetEntries() * (1-wB) * (1-wB) + histos[proc][4]->GetEntries() / histos[proc][3]->GetEntries() * (1-wB);	

		btagSF[proc][6] = 1.0;



	}	



	//Computing QCD errors

	for (int tag = 0; tag < nTagCats; tag++){

	txtfile << "=====================================" << tag << "============================" << endl;



	for (int bin0 = 0; bin0 < nBins; bin0++){


		float a1 = 0.0;
		float a2 = 0.0;
		float den = 0.0;



		for (int bin1 = 0; bin1 < mistagBins[0]; bin1++){

			for (int bin2 = 0; bin2 < mistagBins[1]; bin2++){ 

				for (int bin3 = 0; bin3 < mistagBins[2]; bin3++){

					a1 += pow(mistagH->GetBinContent(bin1, bin2, bin3) * sqrt(arrayForErrors[tag][bin0][bin1][bin2][bin3]), 2);
					a2 += pow(mistagH->GetBinError(bin1, bin2, bin3) * arrayForErrors[tag][bin0][bin1][bin2][bin3], 2);
					den += mistagH->GetBinContent(bin1, bin2, bin3) * arrayForErrors[tag][bin0][bin1][bin2][bin3]; 
					
					 txtfile << bin0 << " " << bin1 << " " << bin2 << " " << bin3 << " " << arrayForErrors[tag][bin0][bin1][bin2][bin3] << endl;

				}
			}
		}





		//Do calcluation
		double x = den;
      		double dx=0.0;
      		//if ( x > 0.0 ) dx = sqrt(a1*a1+a2*a2)/ den;
      		if ( x > 0.0 ) dx = sqrt(a1+a2)/ den;

		histos[names::QCD][tag]->SetBinError( bin0, dx * histos[names::QCD][tag]->GetBinContent(bin0));
		cout << bin0 << "  " << dx << endl;
	}

	}


	//make plot


	float lumi = 19700.;

	for (int tag = 0; tag < nTagCats; tag++){


	histos[names::QCD][tag]->SetMarkerSize(0);

	histos[names::TT7][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7][tag] * 0.074 / 3082812. );
	histos[names::TTLOW][tag]->Scale( 245.8 * lumi  * btagSF[names::TTLOW][tag] * (204. / 16.) / 21675970. );
	histos[names::TT10][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10][tag] * 0.014 / 1249111. );
	histos[names::TT7_SCALEUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_SCALEUP][tag] * 0.074 / 3082812. );
	histos[names::TT10_SCALEUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_SCALEUP][tag] * 0.014 / 1249111. );
	histos[names::TT7_SCALEDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_SCALEDN][tag] * 0.074 / 3082812. );
	histos[names::TT10_SCALEDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_SCALEDN][tag] * 0.014 / 1249111. );
	histos[names::TT7_PUUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_PUUP][tag] * 0.074 / 3082812. );
	histos[names::TT10_PUUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_PUUP][tag] * 0.014 / 1249111. );
	histos[names::TT7_PUDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_PUDN][tag] * 0.074 / 3082812. );
	histos[names::TT10_PUDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_PUDN][tag] * 0.014 / 1249111. );
	histos[names::TT7_JERUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_JERUP][tag] * 0.074 / 3082812. );
	histos[names::TT10_JERUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_JERUP][tag] * 0.014 / 1249111. );
	histos[names::TT7_JERDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_JERDN][tag] * 0.074 / 3082812. );
	histos[names::TT10_JERDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_JERDN][tag] * 0.014 / 1249111. );
	histos[names::TT7_BTAGUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_BTAGUP][tag] * 0.074 / 3082812. );
	histos[names::TT10_BTAGUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_BTAGUP][tag] * 0.014 / 1249111. );
	histos[names::TT7_BTAGDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_BTAGDN][tag] * 0.074 / 3082812. );
	histos[names::TT10_BTAGDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_BTAGDN][tag] * 0.014 / 1249111. );
	histos[names::TT7_Q2UP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_Q2UP][tag] * 0.069 / 2243672. );
	histos[names::TT10_Q2UP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_Q2UP][tag] * 0.013 / 1241650. );
	histos[names::TT7_TOPPTUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_TOPPTUP][tag] * 0.074 / 3082812. );
	histos[names::TT10_TOPPTUP][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_TOPPTUP][tag] * 0.014 / 1249111. );
	histos[names::TT7_Q2DN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_Q2DN][tag] * 0.078 / 2170074. );
	histos[names::TT10_Q2DN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_Q2DN][tag] * 0.016 / 1308090. );
	histos[names::TT7_TOPPTDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_TOPPTDN][tag] * 0.074 / 3082812. );
	histos[names::TT10_TOPPTDN][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_TOPPTDN][tag] * 0.014 / 1249111. );
	histos[names::TT7_MIS][tag]->Scale( 245.8 * lumi  * btagSF[names::TT7_MIS][tag] * 0.074 / 3082812. );
	histos[names::TT10_MIS][tag]->Scale( 245.8 * lumi  * btagSF[names::TT10_MIS][tag] * 0.014 / 1249111. );

	histos[names::TT7_PUUP][tag]->Scale( histos[names::TT7][tag]->Integral() / histos[names::TT7_PUUP][tag]->Integral() );
	histos[names::TT10_PUUP][tag]->Scale( histos[names::TT10][tag]->Integral() / histos[names::TT10_PUUP][tag]->Integral() );
	histos[names::TT7_PUDN][tag]->Scale( histos[names::TT7][tag]->Integral() / histos[names::TT7_PUDN][tag]->Integral() );
	histos[names::TT10_PUDN][tag]->Scale( histos[names::TT10][tag]->Integral() / histos[names::TT10_PUDN][tag]->Integral() );
	histos[names::ZPN10_PUUP][tag]->Scale( histos[names::ZPN10][tag]->Integral() / histos[names::ZPN10_PUUP][tag]->Integral() );
	histos[names::ZPN12_PUUP][tag]->Scale( histos[names::ZPN12][tag]->Integral() / histos[names::ZPN12_PUUP][tag]->Integral() );
	histos[names::ZPN15_PUUP][tag]->Scale( histos[names::ZPN15][tag]->Integral() / histos[names::ZPN15_PUUP][tag]->Integral() );
	histos[names::ZPN20_PUUP][tag]->Scale( histos[names::ZPN20][tag]->Integral() / histos[names::ZPN20_PUUP][tag]->Integral() );
	histos[names::ZPN30_PUUP][tag]->Scale( histos[names::ZPN30][tag]->Integral() / histos[names::ZPN30_PUUP][tag]->Integral() );
	histos[names::ZPN10_PUDN][tag]->Scale( histos[names::ZPN10][tag]->Integral() / histos[names::ZPN10_PUDN][tag]->Integral() );
	histos[names::ZPN12_PUDN][tag]->Scale( histos[names::ZPN12][tag]->Integral() / histos[names::ZPN12_PUDN][tag]->Integral() );
	histos[names::ZPN15_PUDN][tag]->Scale( histos[names::ZPN15][tag]->Integral() / histos[names::ZPN15_PUDN][tag]->Integral() );
	histos[names::ZPN20_PUDN][tag]->Scale( histos[names::ZPN20][tag]->Integral() / histos[names::ZPN20_PUDN][tag]->Integral() );
	histos[names::ZPN30_PUDN][tag]->Scale( histos[names::ZPN30][tag]->Integral() / histos[names::ZPN30_PUDN][tag]->Integral() );
	histos[names::ZPW10_PUUP][tag]->Scale( histos[names::ZPW10][tag]->Integral() / histos[names::ZPW10_PUUP][tag]->Integral() );
	histos[names::ZPW12_PUUP][tag]->Scale( histos[names::ZPW12][tag]->Integral() / histos[names::ZPW12_PUUP][tag]->Integral() );
	histos[names::ZPW15_PUUP][tag]->Scale( histos[names::ZPW15][tag]->Integral() / histos[names::ZPW15_PUUP][tag]->Integral() );
	histos[names::ZPW20_PUUP][tag]->Scale( histos[names::ZPW20][tag]->Integral() / histos[names::ZPW20_PUUP][tag]->Integral() );
	histos[names::ZPW30_PUUP][tag]->Scale( histos[names::ZPW30][tag]->Integral() / histos[names::ZPW30_PUUP][tag]->Integral() );
	histos[names::ZPW10_PUDN][tag]->Scale( histos[names::ZPW10][tag]->Integral() / histos[names::ZPW10_PUDN][tag]->Integral() );
	histos[names::ZPW12_PUDN][tag]->Scale( histos[names::ZPW12][tag]->Integral() / histos[names::ZPW12_PUDN][tag]->Integral() );
	histos[names::ZPW15_PUDN][tag]->Scale( histos[names::ZPW15][tag]->Integral() / histos[names::ZPW15_PUDN][tag]->Integral() );
	histos[names::ZPW20_PUDN][tag]->Scale( histos[names::ZPW20][tag]->Integral() / histos[names::ZPW20_PUDN][tag]->Integral() );
	histos[names::ZPW30_PUDN][tag]->Scale( histos[names::ZPW30][tag]->Integral() / histos[names::ZPW30_PUDN][tag]->Integral() );
	histos[names::RSG10_PUUP][tag]->Scale( histos[names::RSG10][tag]->Integral() / histos[names::RSG10_PUUP][tag]->Integral() );
	histos[names::RSG25_PUUP][tag]->Scale( histos[names::RSG25][tag]->Integral() / histos[names::RSG25_PUUP][tag]->Integral() );
	histos[names::RSG15_PUUP][tag]->Scale( histos[names::RSG15][tag]->Integral() / histos[names::RSG15_PUUP][tag]->Integral() );
	histos[names::RSG20_PUUP][tag]->Scale( histos[names::RSG20][tag]->Integral() / histos[names::RSG20_PUUP][tag]->Integral() );
	histos[names::RSG30_PUUP][tag]->Scale( histos[names::RSG30][tag]->Integral() / histos[names::RSG30_PUUP][tag]->Integral() );
	histos[names::RSG14_PUUP][tag]->Scale( histos[names::RSG14][tag]->Integral() / histos[names::RSG14_PUUP][tag]->Integral() );
	histos[names::RSG18_PUUP][tag]->Scale( histos[names::RSG18][tag]->Integral() / histos[names::RSG18_PUUP][tag]->Integral() );
	histos[names::RSG10_PUDN][tag]->Scale( histos[names::RSG10][tag]->Integral() / histos[names::RSG10_PUDN][tag]->Integral() );
	histos[names::RSG25_PUDN][tag]->Scale( histos[names::RSG25][tag]->Integral() / histos[names::RSG25_PUDN][tag]->Integral() );
	histos[names::RSG15_PUDN][tag]->Scale( histos[names::RSG15][tag]->Integral() / histos[names::RSG15_PUDN][tag]->Integral() );
	histos[names::RSG20_PUDN][tag]->Scale( histos[names::RSG20][tag]->Integral() / histos[names::RSG20_PUDN][tag]->Integral() );
	histos[names::RSG30_PUDN][tag]->Scale( histos[names::RSG30][tag]->Integral() / histos[names::RSG30_PUDN][tag]->Integral() );
	histos[names::RSG14_PUDN][tag]->Scale( histos[names::RSG14][tag]->Integral() / histos[names::RSG14_PUDN][tag]->Integral() );
	histos[names::RSG18_PUDN][tag]->Scale( histos[names::RSG18][tag]->Integral() / histos[names::RSG18_PUDN][tag]->Integral() );

	histos[names::ZPN10][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  / 101697. ); 
	histos[names::ZPN12][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15][tag]->Scale( 1. * lumi * btagSF[names::ZPN15][tag]  /  94371. ); 
	histos[names::ZPN20][tag]->Scale( 1. * lumi * btagSF[names::ZPN20][tag]  /  90778. ); 
	histos[names::ZPN30][tag]->Scale( 1. * lumi * btagSF[names::ZPN30][tag]  /  91209. ); 
	histos[names::ZPN10_MIS][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_MIS][tag]  / 101697. ); 
	histos[names::ZPN15_MIS][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_MIS][tag]  /  94371. ); 
	histos[names::ZPN20_MIS][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_MIS][tag]  /  90778. ); 
	histos[names::ZPN30_MIS][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_MIS][tag]  /  91209. ); 
	histos[names::ZPN10_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_SCALEUP][tag]  / 101697. ); 
	histos[names::ZPN12_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_SCALEUP][tag]  /  94371. ); 
	histos[names::ZPN20_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_SCALEUP][tag]  /  90778. ); 
	histos[names::ZPN30_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_SCALEUP][tag]  /  91209. ); 
	histos[names::ZPN10_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_SCALEDN][tag]  / 101697. ); 
	histos[names::ZPN12_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_SCALEDN][tag]  /  94371. ); 
	histos[names::ZPN20_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_SCALEDN][tag]  /  90778. ); 
	histos[names::ZPN30_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_SCALEDN][tag] /  91209. ); 
	histos[names::ZPN10_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_PUUP][tag]  / 101697. ); 
	histos[names::ZPN12_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_PUUP][tag]  /  94371. ); 
	histos[names::ZPN20_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_PUUP][tag]  /  90778. ); 
	histos[names::ZPN30_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_PUUP][tag]  /  91209. ); 
	histos[names::ZPN10_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_PUDN][tag]  / 101697. ); 
	histos[names::ZPN12_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_PUDN][tag]  /  94371. ); 
	histos[names::ZPN20_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_PUDN][tag]  /  90778. ); 
	histos[names::ZPN30_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_PUDN][tag] /  91209. ); 
	histos[names::ZPN10_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_BTAGUP][tag]  / 101697. ); 
	histos[names::ZPN12_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_BTAGUP][tag]  /  94371. ); 
	histos[names::ZPN20_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_BTAGUP][tag]  /  90778. ); 
	histos[names::ZPN30_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_BTAGUP][tag]  /  91209. ); 
	histos[names::ZPN10_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_BTAGDN][tag]  / 101697. ); 
	histos[names::ZPN12_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_BTAGDN][tag]  /  94371. ); 
	histos[names::ZPN20_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_BTAGDN][tag]  /  90778. ); 
	histos[names::ZPN30_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_BTAGDN][tag] /  91209. ); 
	histos[names::ZPN10_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_JERUP][tag]  / 101697. ); 
	histos[names::ZPN12_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_JERUP][tag]  /  94371. ); 
	histos[names::ZPN20_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_JERUP][tag]  /  90778. ); 
	histos[names::ZPN30_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_JERUP][tag]  /  91209. ); 
	histos[names::ZPN10_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN10_JERDN][tag]  / 101697. ); 
	histos[names::ZPN12_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN10][tag]  /  95793. ); 
	histos[names::ZPN15_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN15_JERDN][tag]  /  94371. ); 
	histos[names::ZPN20_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN20_JERDN][tag]  /  90778. ); 
	histos[names::ZPN30_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPN30_JERDN][tag]  /  91209. ); 
	histos[names::ZPW10][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 102618. ); 
	histos[names::ZPW12][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15][tag]->Scale( 1. * lumi * btagSF[names::ZPW15][tag]  /  96246. ); 
	histos[names::ZPW20][tag]->Scale( 1. * lumi * btagSF[names::ZPW20][tag]  /  93869. ); 
	histos[names::ZPW30][tag]->Scale( 1. * lumi * btagSF[names::ZPW30][tag]  /  97234. ); 
	histos[names::ZPW10_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW10_SCALEUP][tag]  / 102618. ); 
	histos[names::ZPW12_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW15_SCALEUP][tag]  /  96246. ); 
	histos[names::ZPW20_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW20_SCALEUP][tag]  /  93869. ); 
	histos[names::ZPW30_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW30_SCALEUP][tag]  /  97234. ); 
	histos[names::ZPW10_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW10_SCALEDN][tag]  / 102618. ); 
	histos[names::ZPW12_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW15_SCALEDN][tag]  /  96246. ); 
	histos[names::ZPW20_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW20_SCALEDN][tag]  /  93869. ); 
	histos[names::ZPW30_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW30_SCALEDN][tag]  /  97234. ); 
	histos[names::ZPW10_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW10_PUUP][tag]  / 102618. ); 
	histos[names::ZPW12_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW15_PUUP][tag]  /  96246. ); 
	histos[names::ZPW20_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW20_PUUP][tag]  /  93869. ); 
	histos[names::ZPW30_PUUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW30_PUUP][tag]  /  97234. ); 
	histos[names::ZPW10_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW10_PUDN][tag]  / 102618. ); 
	histos[names::ZPW12_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW15_PUDN][tag]  /  96246. ); 
	histos[names::ZPW20_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW20_PUDN][tag]  /  93869. ); 
	histos[names::ZPW30_PUDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW30_PUDN][tag]  /  97234. ); 
	histos[names::ZPW10_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW10_BTAGUP][tag]  / 102618. ); 
	histos[names::ZPW12_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW15_BTAGUP][tag]  /  96246. ); 
	histos[names::ZPW20_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW20_BTAGUP][tag]  /  93869. ); 
	histos[names::ZPW30_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW30_BTAGUP][tag]  /  97234. ); 
	histos[names::ZPW10_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW10_BTAGDN][tag]  / 102618. ); 
	histos[names::ZPW12_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW15_BTAGDN][tag]  /  96246. ); 
	histos[names::ZPW20_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW20_BTAGDN][tag]  /  93869. ); 
	histos[names::ZPW30_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW30_BTAGDN][tag]  /  97234. ); 
	histos[names::ZPW10_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW10_JERUP][tag]  / 102618. ); 
	histos[names::ZPW12_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW15_JERUP][tag]  /  96246. ); 
	histos[names::ZPW20_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW20_JERUP][tag]  /  93869. ); 
	histos[names::ZPW30_JERUP][tag]->Scale( 1. * lumi * btagSF[names::ZPW30_JERUP][tag]  /  97234. ); 
	histos[names::ZPW10_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW10_JERDN][tag]  / 102618. ); 
	histos[names::ZPW12_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW10][tag]  / 98848. ); 
	histos[names::ZPW15_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW15_JERDN][tag]  /  96246. ); 
	histos[names::ZPW20_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW20_JERDN][tag]  /  93869. ); 
	histos[names::ZPW30_JERDN][tag]->Scale( 1. * lumi * btagSF[names::ZPW30_JERDN][tag]  /  97234. ); 
	histos[names::RSG10][tag]->Scale( 1. * lumi * btagSF[names::RSG10][tag]  / 100000. ); 
	histos[names::RSG15][tag]->Scale( 1. * lumi * btagSF[names::RSG15][tag]  /  100000. ); 
	histos[names::RSG20][tag]->Scale( 1. * lumi * btagSF[names::RSG20][tag]  /  100000. ); 
	histos[names::RSG30][tag]->Scale( 1. * lumi * btagSF[names::RSG30][tag]  /  100000. ); 
	histos[names::RSG14][tag]->Scale( 1. * lumi * btagSF[names::RSG14][tag]  /  100000. ); 
	histos[names::RSG18][tag]->Scale( 1. * lumi * btagSF[names::RSG18][tag]  /  100000. ); 
	histos[names::RSG25][tag]->Scale( 1. * lumi * btagSF[names::RSG25][tag]  /  100000. ); 
	histos[names::RSG10_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::RSG10_SCALEUP][tag]  / 100000. ); 
	histos[names::RSG15_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::RSG15_SCALEUP][tag]  /  100000. ); 
	histos[names::RSG20_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::RSG20_SCALEUP][tag]  /  100000. ); 
	histos[names::RSG30_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::RSG30_SCALEUP][tag]  /  100000. ); 
	histos[names::RSG14_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::RSG14_SCALEUP][tag]  /  100000. ); 
	histos[names::RSG18_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::RSG18_SCALEUP][tag]  /  100000. ); 
	histos[names::RSG25_SCALEUP][tag]->Scale( 1. * lumi * btagSF[names::RSG25_SCALEUP][tag]  /  100000. ); 
	histos[names::RSG10_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::RSG10_SCALEDN][tag]  / 100000. ); 
	histos[names::RSG15_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::RSG15_SCALEDN][tag]  /  100000. ); 
	histos[names::RSG20_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::RSG20_SCALEDN][tag]  /  100000. ); 
	histos[names::RSG30_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::RSG30_SCALEDN][tag]  /  100000. ); 
	histos[names::RSG14_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::RSG14_SCALEDN][tag]  /  100000. ); 
	histos[names::RSG18_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::RSG18_SCALEDN][tag]  /  100000. ); 
	histos[names::RSG25_SCALEDN][tag]->Scale( 1. * lumi * btagSF[names::RSG25_SCALEDN][tag]  /  100000. ); 
	histos[names::RSG10_PUUP][tag]->Scale( 1. * lumi * btagSF[names::RSG10_PUUP][tag]  / 100000. ); 
	histos[names::RSG15_PUUP][tag]->Scale( 1. * lumi * btagSF[names::RSG15_PUUP][tag]  /  100000. ); 
	histos[names::RSG20_PUUP][tag]->Scale( 1. * lumi * btagSF[names::RSG20_PUUP][tag]  /  100000. ); 
	histos[names::RSG30_PUUP][tag]->Scale( 1. * lumi * btagSF[names::RSG30_PUUP][tag]  /  100000. ); 
	histos[names::RSG14_PUUP][tag]->Scale( 1. * lumi * btagSF[names::RSG14_PUUP][tag]  /  100000. ); 
	histos[names::RSG18_PUUP][tag]->Scale( 1. * lumi * btagSF[names::RSG18_PUUP][tag]  /  100000. ); 
	histos[names::RSG25_PUUP][tag]->Scale( 1. * lumi * btagSF[names::RSG25_PUUP][tag]  /  100000. ); 
	histos[names::RSG10_PUDN][tag]->Scale( 1. * lumi * btagSF[names::RSG10_PUDN][tag]  / 100000. ); 
	histos[names::RSG15_PUDN][tag]->Scale( 1. * lumi * btagSF[names::RSG15_PUDN][tag]  /  100000. ); 
	histos[names::RSG20_PUDN][tag]->Scale( 1. * lumi * btagSF[names::RSG20_PUDN][tag]  /  100000. ); 
	histos[names::RSG30_PUDN][tag]->Scale( 1. * lumi * btagSF[names::RSG30_PUDN][tag]  /  100000. ); 
	histos[names::RSG14_PUDN][tag]->Scale( 1. * lumi * btagSF[names::RSG14_PUDN][tag]  /  100000. ); 
	histos[names::RSG18_PUDN][tag]->Scale( 1. * lumi * btagSF[names::RSG18_PUDN][tag]  /  100000. ); 
	histos[names::RSG25_PUDN][tag]->Scale( 1. * lumi * btagSF[names::RSG25_PUDN][tag]  /  100000. ); 
	histos[names::RSG10_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::RSG10_BTAGUP][tag]  / 100000. ); 
	histos[names::RSG15_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::RSG15_BTAGUP][tag]  /  100000. ); 
	histos[names::RSG20_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::RSG20_BTAGUP][tag]  /  100000. ); 
	histos[names::RSG30_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::RSG30_BTAGUP][tag]  /  100000. ); 
	histos[names::RSG14_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::RSG14_BTAGUP][tag]  /  100000. ); 
	histos[names::RSG18_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::RSG18_BTAGUP][tag]  /  100000. ); 
	histos[names::RSG25_BTAGUP][tag]->Scale( 1. * lumi * btagSF[names::RSG25_BTAGUP][tag]  /  100000. ); 
	histos[names::RSG10_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::RSG10_BTAGDN][tag]  / 100000. ); 
	histos[names::RSG15_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::RSG15_BTAGDN][tag]  /  100000. ); 
	histos[names::RSG20_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::RSG20_BTAGDN][tag]  /  100000. ); 
	histos[names::RSG30_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::RSG30_BTAGDN][tag]  /  100000. ); 
	histos[names::RSG14_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::RSG14_BTAGDN][tag]  /  100000. ); 
	histos[names::RSG18_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::RSG18_BTAGDN][tag]  /  100000. ); 
	histos[names::RSG25_BTAGDN][tag]->Scale( 1. * lumi * btagSF[names::RSG25_BTAGDN][tag]  /  100000. ); 
	histos[names::RSG10_JERUP][tag]->Scale( 1. * lumi * btagSF[names::RSG10_JERUP][tag]  / 100000. ); 
	histos[names::RSG15_JERUP][tag]->Scale( 1. * lumi * btagSF[names::RSG15_JERUP][tag]  /  100000. ); 
	histos[names::RSG20_JERUP][tag]->Scale( 1. * lumi * btagSF[names::RSG20_JERUP][tag]  /  100000. ); 
	histos[names::RSG30_JERUP][tag]->Scale( 1. * lumi * btagSF[names::RSG30_JERUP][tag]  /  100000. ); 
	histos[names::RSG14_JERUP][tag]->Scale( 1. * lumi * btagSF[names::RSG14_JERUP][tag]  /  100000. ); 
	histos[names::RSG18_JERUP][tag]->Scale( 1. * lumi * btagSF[names::RSG18_JERUP][tag]  /  100000. ); 
	histos[names::RSG25_JERUP][tag]->Scale( 1. * lumi * btagSF[names::RSG25_JERUP][tag]  /  100000. ); 
	histos[names::RSG10_JERDN][tag]->Scale( 1. * lumi * btagSF[names::RSG10_JERDN][tag]  / 100000. ); 
	histos[names::RSG15_JERDN][tag]->Scale( 1. * lumi * btagSF[names::RSG15_JERDN][tag]  /  100000. ); 
	histos[names::RSG20_JERDN][tag]->Scale( 1. * lumi * btagSF[names::RSG20_JERDN][tag]  /  100000. ); 
	histos[names::RSG30_JERDN][tag]->Scale( 1. * lumi * btagSF[names::RSG30_JERDN][tag]  /  100000. ); 
	histos[names::RSG14_JERDN][tag]->Scale( 1. * lumi * btagSF[names::RSG14_JERDN][tag]  /  100000. ); 
	histos[names::RSG18_JERDN][tag]->Scale( 1. * lumi * btagSF[names::RSG18_JERDN][tag]  /  100000. ); 
	histos[names::RSG25_JERDN][tag]->Scale( 1. * lumi * btagSF[names::RSG25_JERDN][tag]  /  100000. ); 

	convertSyst( histos[names::TT10_PDFUP][tag], histos[names::TT10_PDFDN][tag], histos[names::TT10_PDFNOM][tag], histos[names::TT10][tag]);
	convertSyst( histos[names::TT7_PDFUP][tag], histos[names::TT7_PDFDN][tag], histos[names::TT7_PDFNOM][tag], histos[names::TT7][tag]);
	
	convertSyst( histos[names::ZPN10_PDFUP][tag], histos[names::ZPN10_PDFDN][tag], histos[names::ZPN10_PDFNOM][tag], histos[names::ZPN10][tag]);
	convertSyst( histos[names::ZPN12_PDFUP][tag], histos[names::ZPN12_PDFDN][tag], histos[names::ZPN12_PDFNOM][tag], histos[names::ZPN12][tag]);
	convertSyst( histos[names::ZPN15_PDFUP][tag], histos[names::ZPN15_PDFDN][tag], histos[names::ZPN15_PDFNOM][tag], histos[names::ZPN15][tag]);
	convertSyst( histos[names::ZPN20_PDFUP][tag], histos[names::ZPN20_PDFDN][tag], histos[names::ZPN20_PDFNOM][tag], histos[names::ZPN20][tag]);
	convertSyst( histos[names::ZPN30_PDFUP][tag], histos[names::ZPN30_PDFDN][tag], histos[names::ZPN30_PDFNOM][tag], histos[names::ZPN30][tag]);
	convertSyst( histos[names::ZPW10_PDFUP][tag], histos[names::ZPW10_PDFDN][tag], histos[names::ZPW10_PDFNOM][tag], histos[names::ZPW10][tag]);
	convertSyst( histos[names::ZPW12_PDFUP][tag], histos[names::ZPW12_PDFDN][tag], histos[names::ZPW12_PDFNOM][tag], histos[names::ZPW12][tag]);
	convertSyst( histos[names::ZPW15_PDFUP][tag], histos[names::ZPW15_PDFDN][tag], histos[names::ZPW15_PDFNOM][tag], histos[names::ZPW15][tag]);
	convertSyst( histos[names::ZPW20_PDFUP][tag], histos[names::ZPW20_PDFDN][tag], histos[names::ZPW20_PDFNOM][tag], histos[names::ZPW20][tag]);
	convertSyst( histos[names::ZPW30_PDFUP][tag], histos[names::ZPW30_PDFDN][tag], histos[names::ZPW30_PDFNOM][tag], histos[names::ZPW30][tag]);
	convertSyst( histos[names::RSG10_PDFUP][tag], histos[names::RSG10_PDFDN][tag], histos[names::RSG10_PDFNOM][tag], histos[names::RSG10][tag]);
	convertSyst( histos[names::RSG15_PDFUP][tag], histos[names::RSG15_PDFDN][tag], histos[names::RSG15_PDFNOM][tag], histos[names::RSG15][tag]);
	convertSyst( histos[names::RSG20_PDFUP][tag], histos[names::RSG20_PDFDN][tag], histos[names::RSG20_PDFNOM][tag], histos[names::RSG20][tag]);
	convertSyst( histos[names::RSG30_PDFUP][tag], histos[names::RSG30_PDFDN][tag], histos[names::RSG30_PDFNOM][tag], histos[names::RSG30][tag]);
	convertSyst( histos[names::RSG14_PDFUP][tag], histos[names::RSG14_PDFDN][tag], histos[names::RSG14_PDFNOM][tag], histos[names::RSG14][tag]);
	convertSyst( histos[names::RSG18_PDFUP][tag], histos[names::RSG18_PDFDN][tag], histos[names::RSG18_PDFNOM][tag], histos[names::RSG18][tag]);
	convertSyst( histos[names::RSG25_PDFUP][tag], histos[names::RSG25_PDFDN][tag], histos[names::RSG25_PDFNOM][tag], histos[names::RSG25][tag]);


	symmetrize( histos[names::TT10_Q2UP][tag], histos[names::TT10_Q2DN][tag], histos[names::TT10][tag] );
	symmetrize( histos[names::TT7_Q2UP][tag], histos[names::TT7_Q2DN][tag], histos[names::TT7][tag] );


	histos[names::TT7][tag]->Add(histos[names::TT10][tag]);
	//histos[names::TT7][tag]->Add(histos[names::TTLOW][tag]);
	histos[names::TT7_SCALEUP][tag]->Add(histos[names::TT10_SCALEUP][tag]);
	histos[names::TT7_SCALEDN][tag]->Add(histos[names::TT10_SCALEDN][tag]);
	histos[names::TT7_PUUP][tag]->Add(histos[names::TT10_PUUP][tag]);
	histos[names::TT7_PUDN][tag]->Add(histos[names::TT10_PUDN][tag]);
	histos[names::TT7_JERUP][tag]->Add(histos[names::TT10_JERUP][tag]);
	histos[names::TT7_JERDN][tag]->Add(histos[names::TT10_JERDN][tag]);
	histos[names::TT7_BTAGUP][tag]->Add(histos[names::TT10_BTAGUP][tag]);
	histos[names::TT7_BTAGDN][tag]->Add(histos[names::TT10_BTAGDN][tag]);
	histos[names::TT7_TOPPTUP][tag]->Add(histos[names::TT10_TOPPTUP][tag]);
	histos[names::TT7_TOPPTDN][tag]->Add(histos[names::TT10_TOPPTDN][tag]);
	histos[names::TT7_Q2DN][tag]->Add(histos[names::TT10_Q2DN][tag]);
	histos[names::TT7_Q2UP][tag]->Add(histos[names::TT10_Q2UP][tag]);
	histos[names::TT7_PDFDN][tag]->Add(histos[names::TT10_PDFDN][tag]);
	histos[names::TT7_PDFUP][tag]->Add(histos[names::TT10_PDFUP][tag]);





	histos[names::TT7][tag]->SetFillColor(kRed);
	histos[names::QCD][tag]->SetFillColor(kYellow);
	histos[names::QCD][tag]->Add(histos[names::TT7_MIS][tag], -1);
	histos[names::QCD][tag]->Add(histos[names::TT10_MIS][tag], -1);


	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,500);
	c1->Range(0,0,1,1);
	c1->Draw();

	TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.99,0.25);
	c1_1->Draw();
	TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.25,0.99,0.99);
	c1_2->Draw(); 

	c1->cd();
	c1_2->cd();
	c1_2->SetTopMargin(0.1);
	c1_2->SetBottomMargin(0.1);
	c1_2->SetRightMargin(0.05);
	c1_2->SetLeftMargin(0.1);
	c1_2->SetFillStyle(0);




	
	histos[names::DATA][tag]->SetMaximum(1.4 * histos[names::DATA][tag]->GetMaximum() );
	histos[names::DATA][tag]->SetMarkerStyle(20);
	histos[names::DATA][tag]->GetYaxis()->SetTitle("Events / (50 GeV)");
	histos[names::DATA][tag]->Draw("E");	
 
	THStack *stack = new THStack();
	stack->Add(histos[names::TT7][tag]);
	stack->Add(histos[names::QCD][tag]);
	stack->Draw("hist same");
	histos[names::DATA][tag]->SetMarkerStyle(21);
	histos[names::DATA][tag]->Draw("E same");

	histos[names::ZPN10][tag]->SetLineColor(kBlue);
	histos[names::ZPN20][tag]->SetLineColor(kGreen+1);
	histos[names::ZPN30][tag]->SetLineColor(kViolet);
	histos[names::ZPN10][tag]->SetLineWidth(2);
	histos[names::ZPN20][tag]->SetLineWidth(2);
	histos[names::ZPN30][tag]->SetLineWidth(2);

	histos[names::ZPN10][tag]->Draw("hist same");
	histos[names::ZPN20][tag]->Draw("hist same");
	histos[names::ZPN30][tag]->Draw("hist same");


	TH1F *totalH = new TH1F("totalH", "totalH", nBins, min, max);	
    	totalH = (TH1F *) histos[names::QCD][tag]->Clone("totalH");
	totalH->Add(histos[names::TT7][tag]);
	totalH->SetFillStyle(3004);
    	totalH->SetFillColor(kBlack);
    	totalH->Draw("same E2");

	TLatex *cmsLabel = new TLatex();
	cmsLabel->SetNDC();
	cmsLabel->DrawLatex(0.1,0.91, "CMS Preliminary, #sqrt{s} = 8 TeV, 19.7 fb^{-1}");

	if (tag == 0) cmsLabel->DrawLatex(0.7, 0.85, "0 b-tags, |#Deltay| < 1.0");
	if (tag == 1) cmsLabel->DrawLatex(0.7, 0.85, "1 b-tags, |#Deltay| < 1.0");
	if (tag == 2) cmsLabel->DrawLatex(0.7, 0.85, "2 b-tags, |#Deltay| < 1.0");
	if (tag == 3) cmsLabel->DrawLatex(0.7, 0.85, "0 b-tags, |#Deltay| > 1.0");
	if (tag == 4) cmsLabel->DrawLatex(0.7, 0.85, "1 b-tags, |#Deltay| > 1.0");
	if (tag == 5) cmsLabel->DrawLatex(0.7, 0.85, "2 b-tags, |#Deltay| > 1.0");
	if (tag == 6) cmsLabel->DrawLatex(0.7, 0.85, "All Signal Regions");



	gPad->RedrawAxis();
	
	TLegend *leg = new TLegend(0.7, 0.4, 0.94, 0.8);
	leg->AddEntry(histos[names::DATA][tag], "Data", "lp");
	leg->AddEntry(histos[names::QCD][tag], "NTMJ", "f");
	leg->AddEntry(histos[names::TT7][tag], "Top", "f");
	leg->AddEntry(histos[names::ZPN10][tag], "1 TeV Narrow Z'", "l");
	leg->AddEntry(histos[names::ZPN20][tag], "2 TeV Narrow Z'", "l");
	leg->AddEntry(histos[names::ZPN30][tag], "3 TeV Narrow Z'", "l");
	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->Draw("same");

	TH1F *ratioH = new TH1F();
	ratioH = (TH1F*) histos[names::DATA][tag]->Clone("ratio");
	ratioH->Sumw2();
	ratioH->Divide(ratioH, totalH, 1, 1, "B");




	c1_1->cd();
	c1_1->SetTopMargin(0.01);
	c1_1->SetBottomMargin(0.3);
	c1_1->SetRightMargin(0.05);
	c1_1->SetLeftMargin(0.1);
	c1_1->SetFillStyle(0);

	ratioH->GetYaxis()->SetRangeUser(0.,2.);
	ratioH->GetYaxis()->SetTitle("Data / BG Ratio");
	ratioH->GetYaxis()->SetTitleOffset(0.4);
	ratioH->GetYaxis()->SetTitleSize(0.11);
	ratioH->GetXaxis()->SetLabelSize(0.11);
	ratioH->GetXaxis()->SetTitleSize(0.11);
	ratioH->GetXaxis()->SetTitle( "t#bar{t} Invariant Mass [GeV]");
	ratioH->Draw("E");
	
	TF1 *line = new TF1("line", "1", min, max);
	line->SetLineColor(kBlack);
	line->Draw("same");
	
	gPad->RedrawAxis();






	double error = -999.9;
	histos[names::QCD][tag]->IntegralAndError(0, nBins+1, error);
	double ttError = -999.9;
	histos[names::TT7][tag]->IntegralAndError(0, nBins+1, ttError);

	double totError = sqrt(error*error + ttError*ttError);


	cout << tag << " B-TAG CATEGORY" << endl;
	cout << "==================================" << endl;
	cout << "QCD       :" << histos[names::QCD][tag]->Integral() << "  +/-  " << error << endl;
	cout << "ttbar     :" << histos[names::TT7][tag]->Integral() << " +/- " << ttError << endl;
	cout << "ttbar (low)" << histos[names::TTLOW][tag]->Integral() << endl;
	cout << "Total BG  :" << histos[names::QCD][tag]->Integral() + histos[names::TT7][tag]->Integral() << " +/- " << totError << endl;
	cout << "Data      :" << histos[names::DATA][tag]->Integral() << endl;




	




	c1->SaveAs("errors"+tagLabels[tag]+".pdf");
	c1->SaveAs("errors"+tagLabels[tag]+".png");
	c1_2->SetLogy(1);
	c1->SaveAs("errors"+tagLabels[tag]+"_log.pdf");
	c1->SaveAs("errors"+tagLabels[tag]+"_log.png");
	}


	outFile->cd();
	//binFillHist->Write();
	for (int tag = 0; tag < 6; tag++){



		TH1F *upQCDH = (TH1F *) histos[names::QCD][tag]->Clone("upQCDH");
		TH1F *dnQCDH = (TH1F *) histos[names::QCD][tag]->Clone("dnQCDH");

		for (int Nbin = 0; Nbin < histos[names::QCD][tag]->GetNbinsX(); Nbin++){

			if (histos[names::QCD][tag]->GetBinContent(Nbin) <= 0.0) histos[names::QCD][tag]->SetBinContent(Nbin, 0.0);
			upQCDH->SetBinContent( Nbin, histos[names::QCD][tag]->GetBinContent(Nbin) + histos[names::QCD][tag]->GetBinError(Nbin) ); 
			dnQCDH->SetBinContent( Nbin, histos[names::QCD][tag]->GetBinContent(Nbin) - histos[names::QCD][tag]->GetBinError(Nbin) ); 
			//histos[names::QCD][tag]->SetBinError(Nbin, sqrt(histos[names::QCD][tag]->GetBinContent(Nbin)));
			

		}



		cout << "TAG CATEGORY " << tag << endl;

		cout.precision(2);	
		cout << "NRW Z' 1000 :  "  << histos[names::ZPN10][tag]->GetEntries() / 101697. << endl;
		cout << "NRW Z' 1500 :  "  << histos[names::ZPN15][tag]->GetEntries() / 94371. << endl;
		cout << "NRW Z' 2000 :  "  << histos[names::ZPN20][tag]->GetEntries() / 90778. << endl;
		cout << "NRW Z' 3000 :  "  << histos[names::ZPN30][tag]->GetEntries() / 91209. << endl;
 	
		cout << "WDE Z' 1000 :  "  << histos[names::ZPW10][tag]->GetEntries() / 102618. << endl;
		cout << "WDE Z' 1500 :  "  << histos[names::ZPW15][tag]->GetEntries() / 96246. << endl;
		cout << "WDE Z' 2000 :  "  << histos[names::ZPW20][tag]->GetEntries() / 93869. << endl;
		cout << "WDE Z' 3000 :  "  << histos[names::ZPW30][tag]->GetEntries() / 97234. << endl;
		
		cout << "RSG    1000 :  "  << histos[names::RSG10][tag]->GetEntries() / 100000. << endl;
		cout << "RSG    1500 :  "  << histos[names::RSG15][tag]->GetEntries() / 100000. << endl;
		cout << "RSG    2000 :  "  << histos[names::RSG20][tag]->GetEntries() / 100000. << endl;
		cout << "RSG    3000 :  "  << histos[names::RSG30][tag]->GetEntries() / 100000. << endl;

		upQCDH->Write( Form("btag%d__qcd__misErr__plus", tag) );
		dnQCDH->Write( Form("btag%d__qcd__misErr__minus", tag) );

		histos[names::QCD][tag]->Write( Form("btag%d__qcd", tag) );
		histos[names::TT7][tag]->Write( Form("btag%d__ttbar", tag) );
		histos[names::TT7_SCALEUP][tag]->Write( Form("btag%d__ttbar__jec__plus", tag) );
		histos[names::TT7_SCALEDN][tag]->Write( Form("btag%d__ttbar__jec__minus", tag) );
		histos[names::TT7_PUUP][tag]->Write( Form("btag%d__ttbar__pileup__plus", tag) );
		histos[names::TT7_PUDN][tag]->Write( Form("btag%d__ttbar__pileup__minus", tag) );
		histos[names::TT7_JERUP][tag]->Write( Form("btag%d__ttbar__jer__plus", tag) );
		histos[names::TT7_JERDN][tag]->Write( Form("btag%d__ttbar__jer__minus", tag) );
		histos[names::TT7_BTAGUP][tag]->Write( Form("btag%d__ttbar__subjbtag__plus", tag) );
		histos[names::TT7_BTAGDN][tag]->Write( Form("btag%d__ttbar__subjbtag__minus", tag) );
		histos[names::TT7_Q2UP][tag]->Write( Form("btag%d__ttbar__q2__plus", tag) );
		histos[names::TT7_Q2DN][tag]->Write( Form("btag%d__ttbar__q2__minus", tag) );
		//histos[names::TT7_TOPPTUP][tag]->Write( Form("btag%d__ttbar__toppt__plus", tag) );
		//histos[names::TT7_TOPPTDN][tag]->Write( Form("btag%d__ttbar__toppt__minus", tag) );
		histos[names::TT7_PDFUP][tag]->Write( Form("btag%d__ttbar__pdf__plus", tag) );
		histos[names::TT7_PDFDN][tag]->Write( Form("btag%d__ttbar__pdf__minus", tag) );
		histos[names::DATA][tag]->Write( Form("btag%d__DATA", tag) );
		bool writeZPN = ZPNflag;
		if (writeZPN){
		histos[names::ZPN10][tag]->Write( Form("btag%d__Zprime1000", tag));
		histos[names::ZPN12][tag]->Write( Form("btag%d__Zprime1250", tag));
		histos[names::ZPN15][tag]->Write( Form("btag%d__Zprime1500", tag) );
		histos[names::ZPN20][tag]->Write( Form("btag%d__Zprime2000", tag) );
		histos[names::ZPN30][tag]->Write( Form("btag%d__Zprime3000", tag) );
		//histos[names::ZPN10_MIS][tag]->Write( Form("btag%d__Zprime1000_mis", tag));
		//histos[names::ZPN15_MIS][tag]->Write( Form("btag%d__Zprime1500_mis", tag) );
		//histos[names::ZPN20_MIS][tag]->Write( Form("btag%d__Zprime2000_mis", tag) );
		//histos[names::ZPN30_MIS][tag]->Write( Form("btag%d__Zprime3000_mis", tag) );
		histos[names::ZPN10_SCALEUP][tag]->Write( Form("btag%d__Zprime1000__jec__plus", tag));
		histos[names::ZPN12_SCALEUP][tag]->Write( Form("btag%d__Zprime1250__jec__plus", tag));
		histos[names::ZPN15_SCALEUP][tag]->Write( Form("btag%d__Zprime1500__jec__plus", tag) );
		histos[names::ZPN20_SCALEUP][tag]->Write( Form("btag%d__Zprime2000__jec__plus", tag) );
		histos[names::ZPN30_SCALEUP][tag]->Write( Form("btag%d__Zprime3000__jec__plus", tag) );
		histos[names::ZPN10_SCALEDN][tag]->Write( Form("btag%d__Zprime1000__jec__minus", tag));
		histos[names::ZPN12_SCALEDN][tag]->Write( Form("btag%d__Zprime1250__jec__minus", tag));
		histos[names::ZPN15_SCALEDN][tag]->Write( Form("btag%d__Zprime1500__jec__minus", tag) );
		histos[names::ZPN20_SCALEDN][tag]->Write( Form("btag%d__Zprime2000__jec__minus", tag) );
		histos[names::ZPN30_SCALEDN][tag]->Write( Form("btag%d__Zprime3000__jec__minus", tag) );
		histos[names::ZPN10_PUUP][tag]->Write( Form("btag%d__Zprime1000__pileup__plus", tag));
		histos[names::ZPN12_PUUP][tag]->Write( Form("btag%d__Zprime1250__pileup__plus", tag));
		histos[names::ZPN15_PUUP][tag]->Write( Form("btag%d__Zprime1500__pileup__plus", tag) );
		histos[names::ZPN20_PUUP][tag]->Write( Form("btag%d__Zprime2000__pileup__plus", tag) );
		histos[names::ZPN30_PUUP][tag]->Write( Form("btag%d__Zprime3000__pileup__plus", tag) );
		histos[names::ZPN10_PUDN][tag]->Write( Form("btag%d__Zprime1000__pileup__minus", tag));
		histos[names::ZPN12_PUDN][tag]->Write( Form("btag%d__Zprime1250__pileup__minus", tag));
		histos[names::ZPN15_PUDN][tag]->Write( Form("btag%d__Zprime1500__pileup__minus", tag) );
		histos[names::ZPN20_PUDN][tag]->Write( Form("btag%d__Zprime2000__pileup__minus", tag) );
		histos[names::ZPN30_PUDN][tag]->Write( Form("btag%d__Zprime3000__pileup__minus", tag) );
		histos[names::ZPN10_PDFUP][tag]->Write( Form("btag%d__Zprime1000__pdf__plus", tag));
		histos[names::ZPN12_PDFUP][tag]->Write( Form("btag%d__Zprime1250__pdf__plus", tag));
		histos[names::ZPN15_PDFUP][tag]->Write( Form("btag%d__Zprime1500__pdf__plus", tag) );
		histos[names::ZPN20_PDFUP][tag]->Write( Form("btag%d__Zprime2000__pdf__plus", tag) );
		histos[names::ZPN30_PDFUP][tag]->Write( Form("btag%d__Zprime3000__pdf__plus", tag) );
		histos[names::ZPN10_PDFDN][tag]->Write( Form("btag%d__Zprime1000__pdf__minus", tag));
		histos[names::ZPN12_PDFDN][tag]->Write( Form("btag%d__Zprime1250__pdf__minus", tag));
		histos[names::ZPN15_PDFDN][tag]->Write( Form("btag%d__Zprime1500__pdf__minus", tag) );
		histos[names::ZPN20_PDFDN][tag]->Write( Form("btag%d__Zprime2000__pdf__minus", tag) );
		histos[names::ZPN30_PDFDN][tag]->Write( Form("btag%d__Zprime3000__pdf__minus", tag) );
		histos[names::ZPN10_BTAGUP][tag]->Write( Form("btag%d__Zprime1000__subjbtag__plus", tag));
		histos[names::ZPN12_BTAGUP][tag]->Write( Form("btag%d__Zprime1250__subjbtag__plus", tag));
		histos[names::ZPN15_BTAGUP][tag]->Write( Form("btag%d__Zprime1500__subjbtag__plus", tag) );
		histos[names::ZPN20_BTAGUP][tag]->Write( Form("btag%d__Zprime2000__subjbtag__plus", tag) );
		histos[names::ZPN30_BTAGUP][tag]->Write( Form("btag%d__Zprime3000__subjbtag__plus", tag) );
		histos[names::ZPN10_BTAGDN][tag]->Write( Form("btag%d__Zprime1000__subjbtag__minus", tag));
		histos[names::ZPN12_BTAGDN][tag]->Write( Form("btag%d__Zprime1250__subjbtag__minus", tag));
		histos[names::ZPN15_BTAGDN][tag]->Write( Form("btag%d__Zprime1500__subjbtag__minus", tag) );
		histos[names::ZPN20_BTAGDN][tag]->Write( Form("btag%d__Zprime2000__subjbtag__minus", tag) );
		histos[names::ZPN30_BTAGDN][tag]->Write( Form("btag%d__Zprime3000__subjbtag__minus", tag) );
		histos[names::ZPN10_JERUP][tag]->Write( Form("btag%d__Zprime1000__jer__plus", tag));
		histos[names::ZPN12_JERUP][tag]->Write( Form("btag%d__Zprime1250__jer__plus", tag));
		histos[names::ZPN15_JERUP][tag]->Write( Form("btag%d__Zprime1500__jer__plus", tag) );
		histos[names::ZPN20_JERUP][tag]->Write( Form("btag%d__Zprime2000__jer__plus", tag) );
		histos[names::ZPN30_JERUP][tag]->Write( Form("btag%d__Zprime3000__jer__plus", tag) );
		histos[names::ZPN10_JERDN][tag]->Write( Form("btag%d__Zprime1000__jer__minus", tag));
		histos[names::ZPN12_JERDN][tag]->Write( Form("btag%d__Zprime1250__jer__minus", tag));
		histos[names::ZPN15_JERDN][tag]->Write( Form("btag%d__Zprime1500__jer__minus", tag) );
		histos[names::ZPN20_JERDN][tag]->Write( Form("btag%d__Zprime2000__jer__minus", tag) );
		histos[names::ZPN30_JERDN][tag]->Write( Form("btag%d__Zprime3000__jer__minus", tag) );
		}
		bool writeZPW = ZPWflag;
		if (writeZPW){
		histos[names::ZPW10][tag]->Write( Form("btag%d__ZprimeWide1000", tag));
		histos[names::ZPW12][tag]->Write( Form("btag%d__ZprimeWide1250", tag));
		histos[names::ZPW15][tag]->Write( Form("btag%d__ZprimeWide1500", tag) );
		histos[names::ZPW20][tag]->Write( Form("btag%d__ZprimeWide2000", tag) );
		histos[names::ZPW30][tag]->Write( Form("btag%d__ZprimeWide3000", tag) );
		histos[names::ZPW10_SCALEUP][tag]->Write( Form("btag%d__ZprimeWide1000__jec__plus", tag));
		histos[names::ZPW12_SCALEUP][tag]->Write( Form("btag%d__ZprimeWide1250__jec__plus", tag));
		histos[names::ZPW15_SCALEUP][tag]->Write( Form("btag%d__ZprimeWide1500__jec__plus", tag) );
		histos[names::ZPW20_SCALEUP][tag]->Write( Form("btag%d__ZprimeWide2000__jec__plus", tag) );
		histos[names::ZPW30_SCALEUP][tag]->Write( Form("btag%d__ZprimeWide3000__jec__plus", tag) );
		histos[names::ZPW10_SCALEDN][tag]->Write( Form("btag%d__ZprimeWide1000__jec__minus", tag));
		histos[names::ZPW12_SCALEDN][tag]->Write( Form("btag%d__ZprimeWide1250__jec__minus", tag));
		histos[names::ZPW15_SCALEDN][tag]->Write( Form("btag%d__ZprimeWide1500__jec__minus", tag) );
		histos[names::ZPW20_SCALEDN][tag]->Write( Form("btag%d__ZprimeWide2000__jec__minus", tag) );
		histos[names::ZPW30_SCALEDN][tag]->Write( Form("btag%d__ZprimeWide3000__jec__minus", tag) );
		histos[names::ZPW10_PUUP][tag]->Write( Form("btag%d__ZprimeWide1000__pileup__plus", tag));
		histos[names::ZPW12_PUUP][tag]->Write( Form("btag%d__ZprimeWide1250__pileup__plus", tag));
		histos[names::ZPW15_PUUP][tag]->Write( Form("btag%d__ZprimeWide1500__pileup__plus", tag) );
		histos[names::ZPW20_PUUP][tag]->Write( Form("btag%d__ZprimeWide2000__pileup__plus", tag) );
		histos[names::ZPW30_PUUP][tag]->Write( Form("btag%d__ZprimeWide3000__pileup__plus", tag) );
		histos[names::ZPW10_PUDN][tag]->Write( Form("btag%d__ZprimeWide1000__pileup__minus", tag));
		histos[names::ZPW12_PUDN][tag]->Write( Form("btag%d__ZprimeWide1250__pileup__minus", tag));
		histos[names::ZPW15_PUDN][tag]->Write( Form("btag%d__ZprimeWide1500__pileup__minus", tag) );
		histos[names::ZPW20_PUDN][tag]->Write( Form("btag%d__ZprimeWide2000__pileup__minus", tag) );
		histos[names::ZPW30_PUDN][tag]->Write( Form("btag%d__ZprimeWide3000__pileup__minus", tag) );
		histos[names::ZPW10_PDFUP][tag]->Write( Form("btag%d__ZprimeWide1000__pdf__plus", tag));
		histos[names::ZPW12_PDFUP][tag]->Write( Form("btag%d__ZprimeWide1250__pdf__plus", tag));
		histos[names::ZPW15_PDFUP][tag]->Write( Form("btag%d__ZprimeWide1500__pdf__plus", tag) );
		histos[names::ZPW20_PDFUP][tag]->Write( Form("btag%d__ZprimeWide2000__pdf__plus", tag) );
		histos[names::ZPW30_PDFUP][tag]->Write( Form("btag%d__ZprimeWide3000__pdf__plus", tag) );
		histos[names::ZPW10_PDFDN][tag]->Write( Form("btag%d__ZprimeWide1000__pdf__minus", tag));
		histos[names::ZPW12_PDFDN][tag]->Write( Form("btag%d__ZprimeWide1250__pdf__minus", tag));
		histos[names::ZPW15_PDFDN][tag]->Write( Form("btag%d__ZprimeWide1500__pdf__minus", tag) );
		histos[names::ZPW20_PDFDN][tag]->Write( Form("btag%d__ZprimeWide2000__pdf__minus", tag) );
		histos[names::ZPW30_PDFDN][tag]->Write( Form("btag%d__ZprimeWide3000__pdf__minus", tag) );
		histos[names::ZPW10_BTAGUP][tag]->Write( Form("btag%d__ZprimeWide1000__subjbtag__plus", tag));
		histos[names::ZPW12_BTAGUP][tag]->Write( Form("btag%d__ZprimeWide1250__subjbtag__plus", tag));
		histos[names::ZPW15_BTAGUP][tag]->Write( Form("btag%d__ZprimeWide1500__subjbtag__plus", tag) );
		histos[names::ZPW20_BTAGUP][tag]->Write( Form("btag%d__ZprimeWide2000__subjbtag__plus", tag) );
		histos[names::ZPW30_BTAGUP][tag]->Write( Form("btag%d__ZprimeWide3000__subjbtag__plus", tag) );
		histos[names::ZPW10_BTAGDN][tag]->Write( Form("btag%d__ZprimeWide1000__subjbtag__minus", tag));
		histos[names::ZPW12_BTAGDN][tag]->Write( Form("btag%d__ZprimeWide1250__subjbtag__minus", tag));
		histos[names::ZPW15_BTAGDN][tag]->Write( Form("btag%d__ZprimeWide1500__subjbtag__minus", tag) );
		histos[names::ZPW20_BTAGDN][tag]->Write( Form("btag%d__ZprimeWide2000__subjbtag__minus", tag) );
		histos[names::ZPW30_BTAGDN][tag]->Write( Form("btag%d__ZprimeWide3000__subjbtag__minus", tag) );
		histos[names::ZPW10_JERUP][tag]->Write( Form("btag%d__ZprimeWide1000__jer__plus", tag));
		histos[names::ZPW12_JERUP][tag]->Write( Form("btag%d__ZprimeWide1250__jer__plus", tag));
		histos[names::ZPW15_JERUP][tag]->Write( Form("btag%d__ZprimeWide1500__jer__plus", tag) );
		histos[names::ZPW20_JERUP][tag]->Write( Form("btag%d__ZprimeWide2000__jer__plus", tag) );
		histos[names::ZPW30_JERUP][tag]->Write( Form("btag%d__ZprimeWide3000__jer__plus", tag) );
		histos[names::ZPW10_JERDN][tag]->Write( Form("btag%d__ZprimeWide1000__jer__minus", tag));
		histos[names::ZPW12_JERDN][tag]->Write( Form("btag%d__ZprimeWide1250__jer__minus", tag));
		histos[names::ZPW15_JERDN][tag]->Write( Form("btag%d__ZprimeWide1500__jer__minus", tag) );
		histos[names::ZPW20_JERDN][tag]->Write( Form("btag%d__ZprimeWide2000__jer__minus", tag) );
		histos[names::ZPW30_JERDN][tag]->Write( Form("btag%d__ZprimeWide3000__jer__minus", tag) );
		}
		bool writeRSG = RSGflag;
		if (writeRSG){
		histos[names::RSG10][tag]->Write( Form("btag%d__RSgluon1000", tag));
		histos[names::RSG15][tag]->Write( Form("btag%d__RSgluon1500", tag) );
		histos[names::RSG20][tag]->Write( Form("btag%d__RSgluon2000", tag) );
		histos[names::RSG30][tag]->Write( Form("btag%d__RSgluon3000", tag) );
		histos[names::RSG14][tag]->Write( Form("btag%d__RSgluon1400", tag) );
		histos[names::RSG18][tag]->Write( Form("btag%d__RSgluon1800", tag) );
		histos[names::RSG25][tag]->Write( Form("btag%d__RSgluon2500", tag) );
		histos[names::RSG10_SCALEUP][tag]->Write( Form("btag%d__RSgluon1000__jec__plus", tag));
		histos[names::RSG15_SCALEUP][tag]->Write( Form("btag%d__RSgluon1500__jec__plus", tag) );
		histos[names::RSG20_SCALEUP][tag]->Write( Form("btag%d__RSgluon2000__jec__plus", tag) );
		histos[names::RSG30_SCALEUP][tag]->Write( Form("btag%d__RSgluon3000__jec__plus", tag) );
		histos[names::RSG14_SCALEUP][tag]->Write( Form("btag%d__RSgluon1400__jec__plus", tag) );
		histos[names::RSG18_SCALEUP][tag]->Write( Form("btag%d__RSgluon1800__jec__plus", tag) );
		histos[names::RSG25_SCALEUP][tag]->Write( Form("btag%d__RSgluon2500__jec__plus", tag) );
		histos[names::RSG10_SCALEDN][tag]->Write( Form("btag%d__RSgluon1000__jec__minus", tag));
		histos[names::RSG15_SCALEDN][tag]->Write( Form("btag%d__RSgluon1500__jec__minus", tag) );
		histos[names::RSG20_SCALEDN][tag]->Write( Form("btag%d__RSgluon2000__jec__minus", tag) );
		histos[names::RSG30_SCALEDN][tag]->Write( Form("btag%d__RSgluon3000__jec__minus", tag) );
		histos[names::RSG14_SCALEDN][tag]->Write( Form("btag%d__RSgluon1400__jec__minus", tag) );
		histos[names::RSG18_SCALEDN][tag]->Write( Form("btag%d__RSgluon1800__jec__minus", tag) );
		histos[names::RSG25_SCALEDN][tag]->Write( Form("btag%d__RSgluon2500__jec__minus", tag) );
		histos[names::RSG10_PUUP][tag]->Write( Form("btag%d__RSgluon1000__pileup__plus", tag));
		histos[names::RSG15_PUUP][tag]->Write( Form("btag%d__RSgluon1500__pileup__plus", tag) );
		histos[names::RSG20_PUUP][tag]->Write( Form("btag%d__RSgluon2000__pileup__plus", tag) );
		histos[names::RSG30_PUUP][tag]->Write( Form("btag%d__RSgluon3000__pileup__plus", tag) );
		histos[names::RSG14_PUUP][tag]->Write( Form("btag%d__RSgluon1400__pileup__plus", tag) );
		histos[names::RSG18_PUUP][tag]->Write( Form("btag%d__RSgluon1800__pileup__plus", tag) );
		histos[names::RSG25_PUUP][tag]->Write( Form("btag%d__RSgluon2500__pileup__plus", tag) );
		histos[names::RSG10_PUDN][tag]->Write( Form("btag%d__RSgluon1000__pileup__minus", tag));
		histos[names::RSG15_PUDN][tag]->Write( Form("btag%d__RSgluon1500__pileup__minus", tag) );
		histos[names::RSG20_PUDN][tag]->Write( Form("btag%d__RSgluon2000__pileup__minus", tag) );
		histos[names::RSG30_PUDN][tag]->Write( Form("btag%d__RSgluon3000__pileup__minus", tag) );
		histos[names::RSG14_PUDN][tag]->Write( Form("btag%d__RSgluon1400__pileup__minus", tag) );
		histos[names::RSG18_PUDN][tag]->Write( Form("btag%d__RSgluon1800__pileup__minus", tag) );
		histos[names::RSG25_PUDN][tag]->Write( Form("btag%d__RSgluon2500__pileup__minus", tag) );
		histos[names::RSG10_BTAGUP][tag]->Write( Form("btag%d__RSgluon1000__subjbtag__plus", tag));
		histos[names::RSG15_BTAGUP][tag]->Write( Form("btag%d__RSgluon1500__subjbtag__plus", tag) );
		histos[names::RSG20_BTAGUP][tag]->Write( Form("btag%d__RSgluon2000__subjbtag__plus", tag) );
		histos[names::RSG30_BTAGUP][tag]->Write( Form("btag%d__RSgluon3000__subjbtag__plus", tag) );
		histos[names::RSG14_BTAGUP][tag]->Write( Form("btag%d__RSgluon1400__subjbtag__plus", tag) );
		histos[names::RSG18_BTAGUP][tag]->Write( Form("btag%d__RSgluon1800__subjbtag__plus", tag) );
		histos[names::RSG25_BTAGUP][tag]->Write( Form("btag%d__RSgluon2500__subjbtag__plus", tag) );
		histos[names::RSG10_BTAGDN][tag]->Write( Form("btag%d__RSgluon1000__subjbtag__minus", tag));
		histos[names::RSG15_BTAGDN][tag]->Write( Form("btag%d__RSgluon1500__subjbtag__minus", tag) );
		histos[names::RSG20_BTAGDN][tag]->Write( Form("btag%d__RSgluon2000__subjbtag__minus", tag) );
		histos[names::RSG30_BTAGDN][tag]->Write( Form("btag%d__RSgluon3000__subjbtag__minus", tag) );
		histos[names::RSG14_BTAGDN][tag]->Write( Form("btag%d__RSgluon1400__subjbtag__minus", tag) );
		histos[names::RSG18_BTAGDN][tag]->Write( Form("btag%d__RSgluon1800__subjbtag__minus", tag) );
		histos[names::RSG25_BTAGDN][tag]->Write( Form("btag%d__RSgluon2500__subjbtag__minus", tag) );
		histos[names::RSG10_PDFUP][tag]->Write( Form("btag%d__RSgluon1000__pdf__plus", tag));
		histos[names::RSG15_PDFUP][tag]->Write( Form("btag%d__RSgluon1500__pdf__plus", tag) );
		histos[names::RSG20_PDFUP][tag]->Write( Form("btag%d__RSgluon2000__pdf__plus", tag) );
		histos[names::RSG30_PDFUP][tag]->Write( Form("btag%d__RSgluon3000__pdf__plus", tag) );
		histos[names::RSG14_PDFUP][tag]->Write( Form("btag%d__RSgluon1400__pdf__plus", tag) );
		histos[names::RSG18_PDFUP][tag]->Write( Form("btag%d__RSgluon1800__pdf__plus", tag) );
		histos[names::RSG25_PDFUP][tag]->Write( Form("btag%d__RSgluon2500__pdf__plus", tag) );
		histos[names::RSG10_PDFDN][tag]->Write( Form("btag%d__RSgluon1000__pdf__minus", tag));
		histos[names::RSG15_PDFDN][tag]->Write( Form("btag%d__RSgluon1500__pdf__minus", tag) );
		histos[names::RSG20_PDFDN][tag]->Write( Form("btag%d__RSgluon2000__pdf__minus", tag) );
		histos[names::RSG30_PDFDN][tag]->Write( Form("btag%d__RSgluon3000__pdf__minus", tag) );
		histos[names::RSG14_PDFDN][tag]->Write( Form("btag%d__RSgluon1400__pdf__minus", tag) );
		histos[names::RSG18_PDFDN][tag]->Write( Form("btag%d__RSgluon1800__pdf__minus", tag) );
		histos[names::RSG25_PDFDN][tag]->Write( Form("btag%d__RSgluon2500__pdf__minus", tag) );
		histos[names::RSG10_JERUP][tag]->Write( Form("btag%d__RSgluon1000__jer__plus", tag));
		histos[names::RSG15_JERUP][tag]->Write( Form("btag%d__RSgluon1500__jer__plus", tag) );
		histos[names::RSG20_JERUP][tag]->Write( Form("btag%d__RSgluon2000__jer__plus", tag) );
		histos[names::RSG30_JERUP][tag]->Write( Form("btag%d__RSgluon3000__jer__plus", tag) );
		histos[names::RSG14_JERUP][tag]->Write( Form("btag%d__RSgluon1400__jer__plus", tag) );
		histos[names::RSG18_JERUP][tag]->Write( Form("btag%d__RSgluon1800__jer__plus", tag) );
		histos[names::RSG25_JERUP][tag]->Write( Form("btag%d__RSgluon2500__jer__plus", tag) );
		histos[names::RSG10_JERDN][tag]->Write( Form("btag%d__RSgluon1000__jer__minus", tag));
		histos[names::RSG15_JERDN][tag]->Write( Form("btag%d__RSgluon1500__jer__minus", tag) );
		histos[names::RSG20_JERDN][tag]->Write( Form("btag%d__RSgluon2000__jer__minus", tag) );
		histos[names::RSG30_JERDN][tag]->Write( Form("btag%d__RSgluon3000__jer__minus", tag) );
		histos[names::RSG14_JERDN][tag]->Write( Form("btag%d__RSgluon1400__jer__minus", tag) );
		histos[names::RSG18_JERDN][tag]->Write( Form("btag%d__RSgluon1800__jer__minus", tag) );
		histos[names::RSG25_JERDN][tag]->Write( Form("btag%d__RSgluon2500__jer__minus", tag) );
		}





	}
	

	outFile->Close();

	return 0;


}


	
float ptMap(float pt){


	float out = -999.;

	if (pt < 405) out = 0.1;
	else if (pt < 410) out = 1.1;
	else if (pt < 415) out = 2.1;
	else if (pt < 420) out = 3.1;
	else if (pt < 425) out = 4.1;
	else if (pt < 430) out = 5.1;
	else if (pt < 435) out = 6.1;
	else if (pt < 440) out = 7.1;
	else if (pt < 445) out = 8.1;
	else if (pt < 450) out = 9.1;
	else if (pt < 460) out = 10.1;
	else if (pt < 470) out = 11.1;
	else if (pt < 480) out = 12.1;
	else if (pt < 490) out = 13.1;
	else if (pt < 500) out = 14.1;
	else if (pt < 525) out = 15.1;
	else if (pt < 550) out = 16.1;
	else if (pt < 575) out = 17.1;
	else if (pt < 600) out = 18.1;
	else if (pt < 650) out = 19.1;
	else if (pt < 700) out = 20.1;
	else if (pt < 750) out = 21.1;
	else if (pt < 800) out = 22.1;
	else if (pt < 900) out = 23.1;
	else if (pt < 1000) out = 24.1;
	else if (pt < 1100) out = 25.1;
	else if (pt < 1250) out = 26.1;
	else if (pt < 1500) out = 27.1;
	else out = 28.1;

	return out;

}

float bMap(float b){

	float out = -999.;

	if (b < 0.0) out = 0.1;
	else if (b < 0.244) out = 1.1;
	else if (b < 0.679) out = 2.1;
	else if (b < 10.00) out = 3.1;

	return out;

}


float tauMap(float t){

	float out = -999.;

	if (t < 0.4) out = 0.1;
	else if (t < 0.5) out = 1.1;
	else if (t < 0.6) out = 2.1;
	else if (t < 0.7) out = 3.1;
	else if (t < 0.8) out = 4.1;
	else if (t < 0.9) out = 5.1;
	else if (t < 1.0) out = 6.1;
	else if (t < 1.2) out = 7.1;
	else out = 8.1;

	return out;



} 


void convertSyst(TH1F *upH, TH1F *dnH, TH1F *nomH, TH1F *targetH){


	float wUP[upH->GetNbinsX()];
	float wDN[dnH->GetNbinsX()];
	for (int bin = 1; bin < upH->GetNbinsX() ; bin++){

		wUP[bin - 1] = 1.0;
		wDN[bin - 1] = 1.0;

		if (nomH->GetBinContent(bin) > 0.00){
			cout << upH->GetBinContent(bin) << "  " << nomH->GetBinContent(bin) << "  " << dnH->GetBinContent(bin) << endl;
			wUP[bin - 1] = upH->GetBinContent(bin) / nomH->GetBinContent(bin);
			wDN[bin - 1] = dnH->GetBinContent(bin) / nomH->GetBinContent(bin);
		}

		cout << " UP / DN " << wUP[bin - 1] << " " << wDN[bin - 1] << endl;


	}

	for (int bin = 1; bin < targetH->GetNbinsX(); bin++){
	
		upH->SetBinContent( bin, targetH->GetBinContent(bin) * wUP[bin - 1] );
		dnH->SetBinContent( bin, targetH->GetBinContent(bin) * wDN[bin - 1] );
	}



}

void symmetrize( TH1F *upH, TH1F *dnH, TH1F *nomH) {


	float upCount, dnCount, nomCount, delta;
	for (int bin = 1; bin < nomH->GetNbinsX(); bin++){

		upCount = upH->GetBinContent(bin);
		dnCount = dnH->GetBinContent(bin);
		nomCount = nomH->GetBinContent(bin);

		delta = max ( abs(upCount - nomCount), abs(dnCount - nomCount) );

		upH->SetBinContent( bin, nomCount + delta );
		dnH->SetBinContent( bin, nomCount - delta );

	}	
}

	

