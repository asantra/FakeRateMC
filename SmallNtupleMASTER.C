#define SmallNtupleQQQ_cxx
#include "SmallNtupleQQQ.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TMath.h>
#include <iostream>
#include <TSystem.h>
#include <TString.h>
#include <algorithm>
#include <TProfile.h>

using namespace std;

void photonEffectiveAreas(double _eta, double* _effA)
{
  double& effACH(_effA[0]);
  double& effANH(_effA[1]);
  double& effAPh(_effA[2]);
  double& effAW(_effA[3]);

  // Source: CutBasedPhotonID Run2 twiki
  if(_eta < 1.){
    effACH = 0.0234;
    effANH = 0.0053;
    effAPh = 0.078;
    effAW  = 0.075; // from 2012
  }
  else if(_eta > 1. && _eta < 1.479){
    effACH = 0.0189;
    effANH = 0.0103;
    effAPh = 0.0629;
    effAW  = 0.0617; // from 2012
  }
  else if(_eta > 1.479 && _eta < 2.0){
    effACH = 0.0171;
    effANH = 0.0057;
    effAPh = 0.0264;
    effAW  = 0.0617; // from 2012
  }
  else if(_eta > 2.0 && _eta < 2.2){
    effACH = 0.0129;
    effANH = 0.0070;
    effAPh = 0.0462;
    effAW  = 0.0617; // from 2012
  }
  else if(_eta > 2.2 && _eta < 2.3){
    effACH = 0.0110;
    effANH = 0.0152;
    effAPh = 0.0740;
    effAW  = 0.0617; // from 2012
  }
  else if(_eta > 2.3 && _eta < 2.4){
    effACH = 0.0074;
    effANH = 0.0232;
    effAPh = 0.0924;
    effAW  = 0.0617; // from 2012
  }
  else if(_eta > 2.4){
    effACH = 0.0035;
    effANH = 0.1709;
    effAPh = 0.1484;
    effAW  = 0.0617; // from 2012
  }
}

float dRCalc(float etaLead, float phiLead, float etaTrail, float phiTrail){
    
  float dphi = fabs(phiLead - phiTrail);
  if (dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  float deta = fabs(etaLead - etaTrail);
  float dR = sqrt(deta*deta + dphi*dphi);
  return dR;
    
}

float dPhiCalc(float phiLead, float phiTrail){
  float dphi = fabs(phiLead - phiTrail);
  if(dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  return dphi;
}

float findDiEMPt(float ELead, float EtaLead, float PhiLead, float ETrail, float EtaTrail, float PhiTrail){
  float theta1 = 2*atan(exp(-EtaLead));
  float theta2 = 2*atan(exp(-EtaTrail));
  float PX1 = ELead*sin(theta1)*cos(PhiLead);
  float PY1 = ELead*sin(theta1)*sin(PhiLead);
  float PX2 = ETrail*sin(theta2)*cos(PhiTrail);
  float PY2 = ETrail*sin(theta2)*sin(PhiTrail);
  float DiEMPt = sqrt((PX1+PX2)*(PX1+PX2)+(PY1+PY2)*(PY1+PY2));
  return DiEMPt;
}

float findMass(float ELead, float EtaLead, float PhiLead, float ETrail, float EtaTrail, float PhiTrail){
  float theta1 = 2*atan(exp(-EtaLead));
  float theta2 = 2*atan(exp(-EtaTrail));
  float PX1 = ELead*sin(theta1)*cos(PhiLead);
  float PY1 = ELead*sin(theta1)*sin(PhiLead);
  float PX2 = ETrail*sin(theta2)*cos(PhiTrail);
  float PY2 = ETrail*sin(theta2)*sin(PhiTrail);
  float PZ1 = ELead*cos(theta1);
  float PZ2 = ETrail*cos(theta2);
  float Mass = sqrt((ELead+ETrail)*(ELead+ETrail)-((PX1+PX2)*(PX1+PX2)+(PY1+PY2)*(PY1+PY2)+(PZ1+PZ2)*(PZ1+PZ2)));
  return Mass;
}

void Make1DPlot(TH1F *hist1D, const TString str1D, bool logY, const int xlow, const int xhigh){
  TCanvas *BBB = new TCanvas("BBB", "MET Response",600,450);
  BBB->cd();
  hist1D->GetXaxis()->SetRangeUser(xlow, xhigh); 
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.88); 
  gStyle->SetOptStat(111111);
  if(logY)BBB->SetLogy();  
  hist1D->Draw();
  TString name1, name2;
  name1 = str1D+".eps";
  name2 = str1D+".pdf";
  BBB->SaveAs(name1);
  BBB->SaveAs(name2);
  BBB->Delete();
}


void SmallNtupleQQQ::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SmallNtupleQQQ.C
//      Root > SmallNtupleQQQ t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  
  
   //TFile *f1 = new TFile("WJetFile74_Pt_RRR_WWW.root", "RECREATE"); 
   //TFile *f1 = new TFile("DYJetFile74_Pt_RRR_WWW.root", "RECREATE");
   //TFile *f1 = new TFile("DYtoEEFile74_Pt_RRR_WWW.root", "RECREATE");
   //TFile *f1 = new TFile("TTJetFile74_Pt_RRR_WWW.root", "RECREATE");
   //TFile *f1 = new TFile("GJetFile1_Pt_RRR_74_WWW.root", "RECREATE");
   //TFile *f1 = new TFile("GJetFile2_Pt_RRR_74_WWW.root", "RECREATE");
  
  TFile *f1 = new TFile("FILENAME", "RECREATE");
  
   /* Defining Histograms */
   TH1::SetDefaultSumw2();
   TH1F* h_DrAll(new TH1F("h_DrAll", "DeltaR between two objects, ; DeltaR;Events", 50, 0 ,1));
   TH1F* h_Dr_MCElectron_Photon(new TH1F("h_Dr_MCElectron_Photon", "DeltaR between MCElectron and photon, ; DeltaR;Events", 50, 0 ,1));
   TH1F* h_Photon_MET(new TH1F("h_Photon_MET", "PfMet for photons, ; PfMet;Events", 100, 0 ,500));
   TH1F* h_Fake_MET(new TH1F("h_Fake_MET", "PfMet for fakes, ; PfMet;Events", 100, 0 ,500));
   TH1F* h_Fake_Photon_Pt(new TH1F("h_Fake_Photon_Pt","p_{T}for fake photons;p_{T} (GeV);Event/5 GeV", 100, 0, 500));
   TH1F* h_Proxy_Electron_Pt(new TH1F("h_Proxy_Electron_Pt","p_{T} for proxy electrons;p_{T};Events/5 GeV", 100, 0, 500));
   TH1F* h_Proxy_Electron_Matched_Pt(new TH1F("h_Proxy_Electron_Matched_Pt","p_{T} for proxy electrons matched to gen-electron;p_{T};Events/5 GeV", 100, 0, 500));
   TH1F* h_Fake_Photon_Nvtx(new TH1F("h_Fake_Photon_Nvtx","Nvtx distribution for events with fake photons;nVtx;Events/2 ",20,0,40));
   TH1F* h_Fake_Photon_Ntrks(new TH1F("h_Fake_Photon_Ntrks","Ntrks distribution for events with fake photons;nTrks;Events/4 ",100,0,400));
   TH1F* h_Fake_Photon_JetMul(new TH1F("h_Fake_Photon_JetMul","JetMultiplicity distribution for events with fake photons;nJets;Events",30,0,30));
   TH1F* h_Electron_MET(new TH1F("h_Electron_MET", "PfMet for electrons, ; PfMet;Events", 100, 0 ,500));
   TH1F* h_SinglePhoton_MET(new TH1F("h_SinglePhoton_MET", "PfMet for single photon, ; PfMet;Events", 100, 0 ,500));
   TH1F* h_SingleFake_MET(new TH1F("h_SingleFake_MET", "PfMet for single fake, ; PfMet;Events", 100, 0 ,500));
   TH1F* h_SingleElectron_MET(new TH1F("h_SingleElectron_MET", "PfMet for single electron, ; PfMet;Events", 100, 0 ,500));
   TH1F* h_CutFlow_Electron(new TH1F("h_CutFlow_Electron", "CutFlow for electron, ;Cuts;Events",10,0,10));
   TH1F* h_CutFlow_Photon(new TH1F("h_CutFlow_Photon", "CutFlow for photon, ;Cuts;Events",10,0,10));
   TH1F* h_HoverE_Electron(new TH1F("h_HoverE_Electron", "HoverE distribution, Electron, ;H/E;Number",40,0,0.1));
   TH1F* h_HoverE_Photon(new TH1F("h_HoverE_Photon", "HoverE distribution, Photon, ;H/E;Number",40,0,0.1));
   TH1F* h_ShowerShape(new TH1F("h_ShowerShape", "Showershape distribution, ;showershape;Number",80,0,0.05));
   TH1F* h_Pt(new TH1F("h_Pt", "Pt distribution, ;p_{T};Number",100,0,100));
   TH1F* h_Dr_Electrons(new TH1F("h_Dr_Electrons", "#Delta R between mcElectron and recoElectron from PHYS14 selection, ;#Delta R;Number",100,0,3));
   TH1F* h_Dr_ElectronJets(new TH1F("h_Dr_ElectronJets", "#Delta R between recoElectron and jets from PHYS14 selection, ;#Delta R;Number",100,0,3));
   TH1F* h_Dr_MCElectronJets(new TH1F("h_Dr_MCElectronJets", "#Delta R between mcElectron and jets from PHYS14 selection, ;#Delta R;Number",100,0,3));
   TH1F* h_Ht(new TH1F("h_Ht", "Distribtuion of Ht;Ht;Events/4 GeV",1000,0,4000));
   TH1F* h_NTrks_All(new TH1F("h_NTrks_All", "Distribtuion of All Ntrks;NTrks;Events/4 GeV",50,0,200));
   
   TH1F* h_DeltaPhi_LeadPhotonMET(new TH1F("h_DeltaPhi_LeadPhotonMET", "Distribtuion of #Delta #phi between lead photon and MET;#Delta#phi;Events",25,0,3.5));
   TH1F* h_DeltaPhi_TrailPhotonMET(new TH1F("h_DeltaPhi_TrailPhotonMET", "Distribtuion of #Delta #phi between trail photon and MET;#Delta#phi;Events",25,0,3.5));
   TH1F* h_DeltaPhi_LeadElectronMET(new TH1F("h_DeltaPhi_LeadElectronMET", "Distribtuion of #Delta #phi between lead electron and MET;#Delta#phi;Events",25,0,3.5));
   TH1F* h_DeltaPhi_TrailElectronMET(new TH1F("h_DeltaPhi_TrailElectronMET", "Distribtuion of #Delta #phi between trail electron and MET;#Delta#phi;Events",25,0,3.5));
   TH1F* h_DeltaPhi_LeadFakeMET(new TH1F("h_DeltaPhi_LeadFakeMET", "Distribtuion of #Delta #phi between lead fake and MET;#Delta#phi;Events",25,0,3.5));
   TH1F* h_DeltaPhi_TrailFakeMET(new TH1F("h_DeltaPhi_TrailFakeMET", "Distribtuion of #Delta #phi between trail fake and MET;#Delta#phi;Events",25,0,3.5));
   
   
   float xbins[] = {0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,100,125,175,250,400};
   int binnum = sizeof(xbins)/sizeof(Float_t) - 1;
   
   TH1F* h_DoublePhoton_MET(new TH1F("h_DoublePhoton_MET","MET for double photon; MET (GeV); Events",binnum,xbins));
   TH1F* h_DoubleElectron_MET(new TH1F("h_DoubleElectron_MET","MET for double electron; MET (GeV); Events",binnum,xbins));
   TH1F* h_DoubleElectron_MET_Reweighted(new TH1F("h_DoubleElectron_MET_Reweighted","MET for double electron, diEMPt reweighted; MET (GeV); Events",binnum,xbins));
   TH1F* h_DoubleFake_MET(new TH1F("h_DoubleFake_MET","MET for double fake; MET (GeV); Events",binnum,xbins));
   TH1F* h_DoubleFake_MET_Reweighted(new TH1F("h_DoubleFake_MET_Reweighted","MET for double fake, diempt reweighted; MET (GeV); Events",binnum,xbins));
   TH1F* h_ElectronPhoton_MET(new TH1F("h_ElectronPhoton_MET","MET for electron-photon; MET (GeV); Events",binnum,xbins));
   TH1F* h_ElectronPhoton_InvMass(new TH1F("h_ElectronPhoton_InvMass","Invariant Mass for electron-photon; Invariant Mass (GeV); Events/4 GeV",50,0,200));
   TH1F* h_Scaled_ElectronPhoton_MET(new TH1F("h_Scaled_ElectronPhoton_MET","MET for electron-photon, scaled with fake rate; MET (GeV); Events",binnum,xbins));
   TH1F* h_DoubleElectron_DiEMPt(new TH1F("h_DoubleElectron_DiEMPt","DiEMPt for di-electron; DiEMPt (GeV); Events/4 GeV",50,0,200));
   TH1F* h_DoublePhoton_DiEMPt(new TH1F("h_DoublePhoton_DiEMPt","DiEMPt for di-photon; DiEMPt (GeV); Events/4 GeV",50,0,200));
   TH1F* h_DoubleFake_DiEMPt(new TH1F("h_DoubleFake_DiEMPt","DiEMPt for di-fake; DiEMPt (GeV); Events/4 GeV",50,0,200));
   TH1F* h_DoubleElectron_Rho(new TH1F("h_DoubleElectron_Rho","Rho for di-electron; Rho; Events/4 GeV",50,0,200));
   TH1F* h_DoublePhoton_Rho(new TH1F("h_DoublePhoton_Rho","Rho for di-photon; Rho; Events/4 GeV",50,0,200));
   TH1F* h_DoubleFake_Rho(new TH1F("h_DoubleFake_Rho","Rho for di-fake; Rho; Events/4 GeV",50,0,200));
   
   TH1F *h_DoubleElectron_MET_DiEMPtReweighted_ErrorPropagatedFromDiEMPt[1000];
    
   char *histNameee        = new char[50];
   char *histtitleee       = new char[50];

   for (int d=0;d<1000; ++d) {
     sprintf(histNameee, "DiEMPtReweightedee%d",d+1);
     sprintf(histtitleee,"Reweighted ee MET with DIEMPT %d",d+1);
     h_DoubleElectron_MET_DiEMPtReweighted_ErrorPropagatedFromDiEMPt[d]=new TH1F(histNameee,histtitleee,binnum, xbins);
   }
   
   TH1F *h_DoubleFake_MET_DiEMPtReweighted_ErrorPropagatedFromDiEMPt[1000];

   for (int d=0;d<1000; ++d) {
     sprintf(histNameee, "DiEMPtReweightedff%d",d+1);
     sprintf(histtitleee,"Reweighted ff MET with DIEMPT %d",d+1);
     h_DoubleFake_MET_DiEMPtReweighted_ErrorPropagatedFromDiEMPt[d]=new TH1F(histNameee,histtitleee,binnum, xbins);
   }
   
   TProfile* h_Dz_Electron_Ntrks(new TProfile("h_Dz_Electron_Ntrks", "Electron Dz vs Ntrks;Ntrks; Dz",40,0,200,0,10));
   TProfile* h_Dz_Electron_Ntrks_NoJet_AllSpace(new TProfile("h_Dz_Electron_Ntrks_NoJet_AllSpace", "Electron Dz vs Ntrks, No Jet in All Space;Ntrks; Dz",40,0,200,0,10));
   TProfile* h_Dz_Electron_Ntrks_NoJet_Barrel(new TProfile("h_Dz_Electron_Ntrks_NoJet_Barrel", "Electron Dz vs Ntrks, No Jet in Barrel;Ntrks; Dz",40,0,200,0,10));
   
   TH1F* h_Fake_Photon_Pt_AllPtBin(new TH1F("h_Fake_Photon_Pt_AllPtBin","Fake photon pt;Pt (GeV);Events",20,0,200));
   TH1F* h_Proxy_Electron_Matched_Pt_AllPtBin(new TH1F("h_Proxy_Electron_Matched_Pt_AllPtBin","Proxy Electron Pt; Pt (GeV); Events",20,0,200));
   TH1F* h_Fake_Photon_AllNtrksBin(new TH1F("h_Fake_Photon_AllNtrksBin","Fake photon Ntrks; Ntrks; Events",40,0,200));
   TH1F* h_Proxy_Electron_Matched_AllNtrksBin(new TH1F("h_Proxy_Electron_Matched_AllNtrksBin","Proxy electron Ntrks;Ntrks;Events",40,0,200));
   TH1F* h_Fake_Photon_AllNvtxBin(new TH1F("h_Fake_Photon_AllNvtxBin","Fake photon nVtx; Nvtx; Events",8,0,40));
   TH1F* h_Proxy_Electron_Matched_AllNvtxBin(new TH1F("h_Proxy_Electron_Matched_AllNvtxBin","proxy electron nVtx;Nvtx; Events",8,0,40));
   TH1F* h_Fake_Photon_AllJetMulBin(new TH1F("h_Fake_Photon_AllJetMulBin","Fake photon jetmul; JetMul; Events",15,0,30));
   TH1F* h_Proxy_Electron_Matched_AllJetMulBin(new TH1F("h_Proxy_Electron_Matched_AllJetMulBin","Proxy electron jetmul; JetMul; Events",15,0,30));
   
   vector<TH1F*> h_Fake_Photon_Pt_DiffPtBin, h_Proxy_Electron_Matched_Pt_DiffPtBin;
   vector<TH1F*> h_Fake_Photon_Pt_DiffNtrksBin, h_Proxy_Electron_Matched_Pt_DiffNtrksBin;
   vector<TH1F*> h_Fake_Photon_Pt_DiffNvtxBin, h_Proxy_Electron_Matched_Pt_DiffNvtxBin;
   vector<TH1F*> h_Fake_Photon_Pt_DiffJetMulBin, h_Proxy_Electron_Matched_Pt_DiffJetMulBin;
   
   
   TFile *Fnew = new TFile("DIEMPtRatioRRR.root","READ");
   TH1F *eephotondiempt = (TH1F*)Fnew->Get("eephotondiempt");
   TH1F *ffphotondiempt = (TH1F*)Fnew->Get("ffphotondiempt");
   
   TFile *FRandom = new TFile("RandomRatioErrorPropRRR.root","READ");
   TH1F *DiEMPtRatioGaus[1000];
   TH1F *DiEMPtRatioFFGaus[1000];
   
   // getting error propagation from diEMPt ratio
   for(int g=0; g<1000; ++g){
     char *nameGraph        = new char[40];

     sprintf(nameGraph,"DiEMPtee%d",g+1);
     DiEMPtRatioGaus[g] = (TH1F*)FRandom->Get(nameGraph);
     sprintf(nameGraph,"DiEMPtff%d",g+1);
     DiEMPtRatioFFGaus[g] = (TH1F*)FRandom->Get(nameGraph);
   }
   
   //float WEIGHT = 50100*1000*0.318*5/10017930.0;
   
   float fakerate = FFF;
   
   //TH1F* h_ShowerShape_Electron(new TH1F("h_ShowerShape_Electron", "Showershape distribution, electron;showershape;Events",40,0,0.1));
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   bool Debug = false;
   
   int PhotonNu(0), FakeNu(0), ElectronNu(0), totalphoton(0), totalmismatchphoton(0);
   int EventCutFlow[10] = {0};//NeutralCut(0), EventChargedCut(0), EventPhotonIsoCut(0), EventShowerShapeCut(0), EventHoverECut(0), EventAfterPtCut(0);
   int EventElectronDr(0), EventOnlyElectron(0), EventHEElectron(0), TotalFake(0);
   TString CutFlowLabels[9] = {"nocut", "neutral", "charged", "photoniso","showershape","Pt","PixelMatch","HoverE","DeltaR"};
   cout << "Total Events: " << nentries << endl;
   
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if(jentry %10000 == 0)cout << "processed: " << jentry << endl;
      //if(jentry > 200000)break;
      if(Debug)cout << "Just in event loop " << jentry << endl; 
      //if(jentry <  8810000)continue;
      
      int realphoton(0), fake(0), electron(0);
      int NeutralCut(0), ChargedCut(0), PhotonIsoCut(0), ShowerShapeCut(0), HoverECut(0), AfterPtCut(0), HoverECutElectron(0);
      int OnlyPhoton(0), OnlyElectron(0), NoCut(0);
      vector<int> realphoton_tracker, fake_tracker, electron_tracker, showershape_tracker, pt_tracker;
      realphoton_tracker.clear(); fake_tracker.clear(); electron_tracker.clear(); showershape_tracker.clear(); pt_tracker.clear();
      //if(MCPID->size()>0)cout <<"MC: " << MCPID->size() <<endl;
      /*if(jentry==17){
        cout << "nPho: " << nPho << endl;
        cout << "phoSize: " << phoEta->size() << endl;
      }*/ 
      int fakephoton_count(0);
      vector <float> PhotonEta, ElectronEta;
      vector <float> PhotonPhi, ElectronPhi;
      vector <float> GoodMuonEta, GoodMuonPhi;
      vector <float> FakePhotonPt, ElectronMatchedPtVect, DzVect;
      vector <float> NTrksVect;
      PhotonEta.clear(); ElectronEta.clear();
      GoodMuonEta.clear(); GoodMuonPhi.clear();
      PhotonPhi.clear(); ElectronPhi.clear();
      FakePhotonPt.clear(); ElectronMatchedPtVect.clear(); DzVect.clear();
      NTrksVect.clear();
      
      h_NTrks_All->Fill(NTrks, WEIGHT);
      
      
      // photon loop starts
      for(size_t i=0;i < PhoEta->size(); ++i){
        //cout << i << endl;
        if(Debug)cout << "Just in photon loop " << jentry << endl;
        float absEta = fabs(PhoEta->at(i));
        double effA[4];
        photonEffectiveAreas(absEta, effA);
        if((absEta > 1.479) || (PhoR9->at(i) > 1.0))continue;
	// containing newest cuts fron https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2
	
        float pt = fabs(PhoEt->at(i)); 
        bool neutral     = (PhoPFNeuIso->at(i) - Rho * effA[1] -  exp(0.0028*pt+0.5408)) < 7.23;
        bool charged     = (PhoPFChIso->at(i) - Rho * effA[0]) < 2.67;
        bool chargedfake = (PhoPFChIso->at(i) - Rho * effA[0]) > 2.67;
        bool chargedlim  = (PhoPFChIso->at(i) - Rho * effA[0]) < 15.;
        bool photoniso   = (PhoPFPhoIso->at(i) - Rho * effA[2] - 0.0014*pt) < 2.11;
        //printf("absEta: %f and effACh: %f effPh:%f and effN:%f\n", absEta,effA[0],effA[2],effA[1]);
        bool worstiso    = (PhoPFChWorstIso->at(i) -Rho*effA[3]) < 2.6;
        bool showercut   = PhoSigmaIEtaIEta->at(i) < 0.0107 ; 
        bool shower      = PhoSigmaIEtaIEta->at(i) < 0.014 ;
        bool pixelcut    = PhohasPixelSeed->at(i) == 0;
        //cout << phohasPixelSeed->at(i) << endl;
        bool HoverE      = PhoHoverE->at(i) < 0.028;
        bool SymPtCut    = pt > 30;
        
        if(Debug)cout << "photon selection ends " << jentry << endl;
        //if(jentry==535519)printf("yy pixelseed head:%d\n",photon_pixelseed->at(i));
        NoCut++;
        if(neutral){
          NeutralCut++;
          if(charged){
            ChargedCut++;
            if(photoniso){
              PhotonIsoCut++;
              showershape_tracker.push_back(i);
              //h_ShowerShape->Fill(phoSigmaIEtaIEta->at(i));
              if(showercut){
                ShowerShapeCut++;
                pt_tracker.push_back(i);
                if(SymPtCut){
                  AfterPtCut++;
                  if(pixelcut){
                    OnlyPhoton++;
                    h_HoverE_Photon->Fill(PhoHoverE->at(i));
                    if(HoverE){
                      HoverECut++;
                    }
                  }
                  if(!pixelcut){
                    OnlyElectron++;
                    h_HoverE_Electron->Fill(PhoHoverE->at(i));
                    if(HoverE){
                      HoverECutElectron++;
                    }
                  }
                }
              }
            } 
          }
        }
        //if(neutral && charged && photoniso && pixelcut && HoverE && SymPtCut )h_ShowerShape_Photon->Fill(phoSigmaIEtaIEta->at(i));
        //if(neutral && charged && photoniso && !pixelcut && HoverE && SymPtCut )h_ShowerShape_Electron->Fill(phoSigmaIEtaIEta->at(i));
        
        
        
        
        /////////////// photon //////////////////////////////////
        if(neutral && charged && photoniso && showercut && pixelcut && HoverE && SymPtCut ){//
          //cout << "photon inside loop: " << jentry << endl; 
	  // photon GenLevel Matching
	  bool GenMatchPhoton = false;
	  for(size_t m=0; m < MCPID->size(); ++m){
	    //if(abs(MCPID->at(m))==22)cout << "Photon MC Status: " << MCStatus->at(m) << endl;
	    if(!(abs(MCPID->at(m))==22))continue;
	    float drminMC(999);
	    float drPhotonMC = dRCalc(MCEta->at(m), MCPhi->at(m), PhoEta->at(i), PhoPhi->at(i));
	    if(!(drPhotonMC < 0.2 && drPhotonMC < drminMC))continue;
	    drminMC = drPhotonMC;
	    GenMatchPhoton = true;
	  }
	  if(true){
	    realphoton++;
	    totalphoton++;
	    PhotonEta.push_back(PhoEta->at(i));
	    PhotonPhi.push_back(PhoPhi->at(i));
	    realphoton_tracker.push_back(i);
	  }
	  float drmin(9999);
	  float FakePhoPt(0), Dz(0);
	  for(size_t m=0; m < MCPID->size(); ++m){
	    //cout << "=========== entry: " << jentry << "=============" << endl;
	    //cout << "Inside Photon PID: " << MCPID->at(m) << " and Status: " << MCStatus->at(m) << endl;
	    if(FFFF)continue; 
	    float drfake = dRCalc(MCEta->at(m), MCPhi->at(m), PhoEta->at(i), PhoPhi->at(i));
	    h_Dr_MCElectron_Photon->Fill(drfake);
	    if(!(drfake < 0.2 && drfake < drmin))continue;
	    drmin = drfake;
	    FakePhoPt = PhoEt->at(i);
	    float drminforele(9999);
	    for(size_t q=0; q < ElePt->size();++q){
	      float drMCandE = dRCalc(MCEta->at(m), MCPhi->at(m), EleEta->at(q), ElePhi->at(q));
	      if(!(drMCandE < 0.2 && drMCandE < drminforele))continue;
	      drminforele = drMCandE;
	      Dz = EleDz->at(q);
	    }
	
	  }
          if(FakePhoPt!=0){
            FakePhotonPt.push_back(FakePhoPt);
	    DzVect.push_back(fabs(Dz));
	    NTrksVect.push_back(NTrks);
	    h_Fake_Photon_Pt_AllPtBin->Fill(FakePhoPt, WEIGHT);
            h_Fake_Photon_AllNtrksBin->Fill(NTrks, WEIGHT);
            h_Fake_Photon_AllNvtxBin->Fill(NVtx, WEIGHT);
            h_Fake_Photon_Pt->Fill(FakePhoPt, WEIGHT);
	    h_Dz_Electron_Ntrks->Fill(NTrks, fabs(Dz), WEIGHT);
            fakephoton_count++;
          }
          if(PhohasPixelSeed->at(i) != PhoEleVeto->at(i))totalmismatchphoton++;
          if(Debug)cout << "doublephoton selection  " << jentry << endl;
        } //worstiso not used*/
        
        
        
        
        
        //////////////////// fake ////////////////////////////////
        if(neutral && !charged && chargedlim && photoniso && showercut && pixelcut  && HoverE && SymPtCut ){
          fake++;
          fake_tracker.push_back(i);
          if(Debug)cout << "doublefake selection  " << jentry << endl;
        }
        
        
        
        
        
        
        ////////////////// electron ///////////////////////////////////
        if(neutral && charged && photoniso && showercut && !pixelcut && HoverE && SymPtCut ){// 
           electron++;
	  //cout << "electron inside loop: " << jentry << endl;
	  ElectronEta.push_back(PhoEta->at(i));
	  ElectronPhi.push_back(PhoPhi->at(i));
	  electron_tracker.push_back(i);
          float drmin(9999);
          float ElectronMatchedPt(0);
          for(size_t m=0; m < MCPID->size(); ++m){
	    if(FFFF)continue;  
            float drele = dRCalc(MCEta->at(m), MCPhi->at(m), PhoEta->at(i), PhoPhi->at(i));
            if(!(drele < 0.2 && drele < drmin))continue;
            drmin = drele;
            ElectronMatchedPt = PhoEt->at(i);
          }
          if(ElectronMatchedPt!=0){
	   
	    
            ElectronMatchedPtVect.push_back(ElectronMatchedPt);
            h_Proxy_Electron_Matched_Pt->Fill(ElectronMatchedPt);
	    h_Proxy_Electron_Matched_Pt_AllPtBin->Fill(ElectronMatchedPt, WEIGHT);
            h_Proxy_Electron_Matched_AllNtrksBin->Fill(NTrks, WEIGHT);
            h_Proxy_Electron_Matched_AllNvtxBin->Fill(NVtx, WEIGHT);
          }
          h_Proxy_Electron_Pt->Fill(PhoEt->at(i));
          if(Debug)cout << "doubleelectron selection  " << jentry << endl;
        }
        if(Debug)cout << "doubleelectron selection ends  " << jentry << endl;
      }// photon loop ends
      
      
      
      
      
      
      //muon selection
      vector<int> GoodMuonTracker;
      GoodMuonTracker.clear();
      if(Debug)cout << "Muon selection for entry: " << jentry << endl;
      for(size_t w=0; w < MuPt->size();++w){
        if(!(MuIsGood->at(w)==1 && fabs(MuEta->at(w))<1.479 && MuPt->at(w)>5))continue;
        //cout << "max: " << MuEta->at(w) << endl;
        float combIso  = MuPFChIso->at(w) + TMath::Max(0.,MuPFNeuIso->at(w) + MuPFPhoIso->at(w) - 0.5*MuPFPUIso->at(w));
        if( combIso/MuPt->at(w) > 0.20 )continue;
        GoodMuonEta.push_back(MuEta->at(w));
        GoodMuonPhi.push_back(MuPhi->at(w));
        GoodMuonTracker.push_back(w);
      }
      
      
//       int MCENumber(0);
//       for(size_t m=0; m < MCPID->size(); ++m){
//         if(!(abs(MCPID->at(m))==11 && (abs(MCMomPID->at(m))==23  || abs(MCMomPID->at(m)==22))))continue;
// 	MCENumber++;
//       }
      
      
      
      //jet selection
      unsigned int jetNum(0), jetNumAll(0);
      float ht(0);
      for(size_t w=0; w < JetPt->size(); ++w){
        if(!(JetPt->at(w) > 30 && JetPFLooseId))continue;//&& ) 
        float drminpho(9999), drminele(9999), drminmu(9999);
        if(Debug)cout << "photon: " << PhotonEta.size() << " electron: " << ElectronEta.size() << " muon: " << GoodMuonEta.size() << " for entry: " << jentry << endl;
        for(size_t r=0; r < PhotonEta.size(); ++r){
          float drJetPho = dRCalc(JetEta->at(w), JetPhi->at(w), PhotonEta.at(r), PhotonPhi.at(r));
          if(!(drJetPho < drminpho))continue;
          drminpho = drJetPho;
        }
        for(size_t r=0; r < ElectronEta.size(); ++r){
          float drJetEle = dRCalc(JetEta->at(w), JetPhi->at(w), ElectronEta.at(r), ElectronPhi.at(r));
          if(!(drJetEle < drminele))continue;
          drminele = drJetEle;
        }
        for(size_t r=0; r < GoodMuonEta.size(); ++r){
          float drJetMu = dRCalc(JetEta->at(w), JetPhi->at(w), GoodMuonEta.at(r), GoodMuonPhi.at(r));
          if(!(drJetMu < drminmu))continue;
          drminmu = drJetMu;
        }
        if( PhotonEta.size() > 0 && drminpho < 0.3 )continue;
        if( ElectronEta.size() > 0 && drminele < 0.3 )continue;
        if( GoodMuonEta.size() > 0 && drminmu < 0.3 )continue;
        if(Debug)cout << "Inside jet selection loop for entry: " << jentry << endl;
        jetNumAll++;
	if(fabs(JetEta->at(w) > 2.4))continue;
	//if(MCENumber < 2)continue;
	jetNum++;
        ht += JetPt->at(w);
        //cout << "Ht: " << ht << endl;
      }
      h_Ht->Fill(ht, WEIGHT);
      for(size_t w=0; w < ElectronMatchedPtVect.size(); ++w){
        h_Proxy_Electron_Matched_AllJetMulBin->Fill(jetNum, WEIGHT);
      }
      for(size_t w=0; w < FakePhotonPt.size(); ++w){
        h_Fake_Photon_AllJetMulBin->Fill(jetNum, WEIGHT);
      }
      for(unsigned int s = 0; s < NTrksVect.size(); ++s){
	if(jetNumAll == 0)h_Dz_Electron_Ntrks_NoJet_AllSpace->Fill(NTrksVect.at(s), DzVect.at(s), WEIGHT);
	if(jetNum == 0) h_Dz_Electron_Ntrks_NoJet_Barrel->Fill(NTrksVect.at(s), DzVect.at(s), WEIGHT);
	
      }
      
      if(Debug)cout << "total jets: " << jetNum << " for entry:" << jentry << endl;
      
      
      
      if(fakephoton_count !=0){
        h_Fake_Photon_Nvtx->Fill(NVtx, WEIGHT);
        h_Fake_Photon_Ntrks->Fill(NTrks, WEIGHT);
        h_Fake_Photon_JetMul->Fill(jetNum, WEIGHT);
      }
      // if (Cut(ientry) < 0) continue;
      if(realphoton == 1)h_SinglePhoton_MET->Fill(PFMET);
      if(realphoton >=2){
        PhotonNu++;
        float DR = dRCalc(PhoEta->at(realphoton_tracker.at(0)), PhoPhi->at(realphoton_tracker.at(0)), PhoEta->at(realphoton_tracker.at(1)), PhoPhi->at(realphoton_tracker.at(1)));
	float DiEMPt = findDiEMPt(PhoE->at(realphoton_tracker.at(0)), PhoEta->at(realphoton_tracker.at(0)), PhoPhi->at(realphoton_tracker.at(0)), PhoE->at(realphoton_tracker.at(1)), PhoEta->at(realphoton_tracker.at(1)), PhoPhi->at(realphoton_tracker.at(1)));
        if(DR>0.5){
          h_DoublePhoton_DiEMPt->Fill(DiEMPt, WEIGHT);
	  h_DoublePhoton_Rho->Fill(Rho, WEIGHT);
	  h_DoublePhoton_MET->Fill(PFMET, WEIGHT);
	  if(PhoEt->at(realphoton_tracker.at(0)) > PhoEt->at(realphoton_tracker.at(1))){
	    h_DeltaPhi_LeadPhotonMET->Fill(dPhiCalc(PhoPhi->at(realphoton_tracker.at(0)), PFMETPhi));
	    h_DeltaPhi_TrailPhotonMET->Fill(dPhiCalc(PhoPhi->at(realphoton_tracker.at(1)), PFMETPhi));
	  }
	  else{
	    h_DeltaPhi_LeadPhotonMET->Fill(dPhiCalc(PhoPhi->at(realphoton_tracker.at(1)), PFMETPhi));
	    h_DeltaPhi_TrailPhotonMET->Fill(dPhiCalc(PhoPhi->at(realphoton_tracker.at(0)), PFMETPhi));
	  }
          EventCutFlow[8]++;
          if(Debug)cout << "inside doublephoton  " << jentry << endl;
          h_Photon_MET->Fill(PFMET);
        }
      }
      
      // e-gamma sample selection //
      if(electron==1 && realphoton==1){
	float DR = dRCalc(PhoEta->at(realphoton_tracker.at(0)), PhoPhi->at(realphoton_tracker.at(0)), PhoEta->at(electron_tracker.at(0)), PhoPhi->at(electron_tracker.at(0)));
	float Mass = findMass(PhoE->at(realphoton_tracker.at(0)), PhoEta->at(realphoton_tracker.at(0)), PhoPhi->at(realphoton_tracker.at(0)), PhoE->at(electron_tracker.at(0)), PhoEta->at(electron_tracker.at(0)), PhoPhi->at(electron_tracker.at(0)));
        if(DR > 0.5){
	  h_ElectronPhoton_MET->Fill(PFMET, WEIGHT);
	  h_ElectronPhoton_InvMass->Fill(Mass, WEIGHT);
	  h_Scaled_ElectronPhoton_MET->Fill(PFMET, fakerate/(1-fakerate)*WEIGHT);
	}
      }
      
      if(fake == 1)h_SingleFake_MET->Fill(PFMET);
      if(fake >=2){
        FakeNu++;
        float DR = dRCalc(PhoEta->at(fake_tracker.at(0)), PhoPhi->at(fake_tracker.at(0)), PhoEta->at(fake_tracker.at(1)), PhoPhi->at(fake_tracker.at(1)));
	float DiEMPt = findDiEMPt(PhoE->at(fake_tracker.at(0)), PhoEta->at(fake_tracker.at(0)), PhoPhi->at(fake_tracker.at(0)), PhoE->at(fake_tracker.at(1)), PhoEta->at(fake_tracker.at(1)), PhoPhi->at(fake_tracker.at(1)));
        if(DR>0.5){
	  if(PhoEt->at(fake_tracker.at(0)) > PhoEt->at(fake_tracker.at(1))){
	    h_DeltaPhi_LeadFakeMET->Fill(dPhiCalc(PhoPhi->at(fake_tracker.at(0)), PFMETPhi));
	    h_DeltaPhi_TrailFakeMET->Fill(dPhiCalc(PhoPhi->at(fake_tracker.at(1)), PFMETPhi));
	  }
	  else{
	    h_DeltaPhi_LeadFakeMET->Fill(dPhiCalc(PhoPhi->at(fake_tracker.at(1)), PFMETPhi));
	    h_DeltaPhi_TrailFakeMET->Fill(dPhiCalc(PhoPhi->at(fake_tracker.at(0)), PFMETPhi));
	  }
          
	  h_DoubleFake_MET->Fill(PFMET, WEIGHT);
	  h_DoubleFake_DiEMPt->Fill(DiEMPt, WEIGHT);
	  h_DoubleFake_Rho->Fill(Rho, WEIGHT);
	  int BinC=ffphotondiempt->FindBin(DiEMPt);
          float WeightC=ffphotondiempt->GetBinContent(BinC);
	  h_DoubleFake_MET_Reweighted->Fill(PFMET, WEIGHT*WeightC);
	  for( int s=0 ; s<1000; ++s){
            float Weightprop = DiEMPtRatioFFGaus[s]->GetBinContent(BinC);
            h_DoubleFake_MET_DiEMPtReweighted_ErrorPropagatedFromDiEMPt[s]->Fill(PFMET, Weightprop*WEIGHT);
          }
          TotalFake++;
          h_Fake_MET->Fill(PFMET);
          if(Debug)cout << "inside doublefake  " << jentry << endl;
        }
      }
      if(electron == 1)h_SingleElectron_MET->Fill(PFMET);
      if(electron >=2){
        ElectronNu++;
        float DR = dRCalc(PhoEta->at(electron_tracker.at(0)), PhoPhi->at(electron_tracker.at(0)), PhoEta->at(electron_tracker.at(1)), PhoPhi->at(electron_tracker.at(1)));
	float Mass = findMass(PhoE->at(electron_tracker.at(0)), PhoEta->at(electron_tracker.at(0)), PhoPhi->at(electron_tracker.at(0)), PhoE->at(electron_tracker.at(1)), PhoEta->at(electron_tracker.at(1)), PhoPhi->at(electron_tracker.at(1)));
	float DiEMPt = findDiEMPt(PhoE->at(electron_tracker.at(0)), PhoEta->at(electron_tracker.at(0)), PhoPhi->at(electron_tracker.at(0)), PhoE->at(electron_tracker.at(1)), PhoEta->at(electron_tracker.at(1)), PhoPhi->at(electron_tracker.at(1)));
        h_DrAll->Fill(DR);
        if(DR > 0.5 && Mass >=75 && Mass <=105){
          h_DoubleElectron_DiEMPt->Fill(DiEMPt, WEIGHT);
	  h_DoubleElectron_Rho->Fill(Rho, WEIGHT);
	  if(PhoEt->at(electron_tracker.at(0)) > PhoEt->at(electron_tracker.at(1))){
	    h_DeltaPhi_LeadElectronMET->Fill(dPhiCalc(PhoPhi->at(electron_tracker.at(0)), PFMETPhi));
	    h_DeltaPhi_TrailElectronMET->Fill(dPhiCalc(PhoPhi->at(electron_tracker.at(1)), PFMETPhi));
	  }
	  else{
	    h_DeltaPhi_LeadElectronMET->Fill(dPhiCalc(PhoPhi->at(electron_tracker.at(1)), PFMETPhi));
	    h_DeltaPhi_TrailElectronMET->Fill(dPhiCalc(PhoPhi->at(electron_tracker.at(0)), PFMETPhi));
	  }
	  int BinC=eephotondiempt->FindBin(DiEMPt);
          float WeightC=eephotondiempt->GetBinContent(BinC);
	  h_DoubleElectron_MET_Reweighted->Fill(PFMET, WEIGHT*WeightC);
	  h_DoubleElectron_MET->Fill(PFMET, WEIGHT);
	  for( int s=0 ; s<1000; ++s){
            float Weightprop = DiEMPtRatioGaus[s]->GetBinContent(BinC);
            h_DoubleElectron_MET_DiEMPtReweighted_ErrorPropagatedFromDiEMPt[s]->Fill(PFMET, Weightprop*WEIGHT);
          }
          for(int l = 0; l < electron; ++l){
            //cout << "Entry: " << jentry << " electron pixel seed: " << PhohasPixelSeed->at(electron_tracker.at(l)) << " conversion safe electron veto: " << PhoEleVeto->at(electron_tracker.at(l)) << endl;
            vector<float> MinDR, MinDRJet;
            //if(jentry==8813874)continue;
            MinDR.clear(); MinDRJet.clear();
            /// working with MC electron ////
            for(size_t m = 0; m < MCPID->size(); ++m){
              if(Debug)cout << "inside MCElectron  " << jentry << endl;
              if(Debug && jentry> 8813873){
                cout << "inside MCElectron  " << jentry << endl;
                cout << "PID(" << m << "): " << MCPID->at(m) << endl;
                cout << "MCEta(" << m << "): " << MCEta->at(m) << endl;
                cout << "MCPhi(" << m << "): " << MCPhi->at(m) << endl;
              }
              if(abs(MCPID->at(m))==11 && MCStatus->at(m)==1){
                if(Debug)cout << "MCElectron!!!" << endl;
                if(Debug)cout << MCStatus->at(m) << endl;
                float DRElectron = dRCalc(PhoEta->at(electron_tracker.at(l)), PhoPhi->at(electron_tracker.at(l)), MCEta->at(m), MCPhi->at(m));
                MinDR.push_back(DRElectron);
              }
            }
            float mindr = 0;
            if(MinDR.size()!=0) mindr = *std::min_element(MinDR.begin(),MinDR.end());
            if(Debug && jentry > 8813873)cout << "MinDr: " << mindr << endl;
            if(mindr!=0)h_Dr_Electrons->Fill(mindr);
          
            /// working with Jets ////
            for(size_t m = 0; m < JetEta->size(); ++m){
              if(Debug)cout << "inside Jet " << jentry << endl;
              if(Debug){
                cout << "inside Jet  " << jentry << endl;
                cout << "JetPt(" << m << "): " << JetPt->at(m) << endl;
                cout << "JetEta(" << m << "): " << JetEta->at(m) << endl;
                cout << "JetPhi(" << m << "): " << JetPhi->at(m) << endl;
              }
              if(Debug)cout << "Jets!!!" << endl;
              float DRJet = dRCalc(PhoEta->at(electron_tracker.at(l)), PhoPhi->at(electron_tracker.at(l)), JetEta->at(m), JetPhi->at(m));
              MinDRJet.push_back(DRJet);
            }
            float mindrjet = 0;
            if(MinDRJet.size()!=0) mindrjet = *std::min_element(MinDRJet.begin(),MinDRJet.end());
            if(Debug)cout << "MinDrJet: " << mindrjet << endl;
            if(mindrjet!=0)h_Dr_ElectronJets->Fill(mindrjet);
          }
          if(Debug)cout << "inside doubleelectron  " << jentry << endl;
          EventElectronDr++;
          h_Electron_MET->Fill(PFMET);
        }
      }
      
      
      
      
      
      
      // working with MC eletron and photon ///
      int MCE = 0, MCPhoton(0);
      vector<int> MCE_tracker;
      MCE_tracker.clear();
      for(size_t r=0; r < MCPID->size(); ++r){
        if(MCPID->at(r)==11){MCE++;MCE_tracker.push_back(r);}
      }
      
      //cout << MCE << endl;
      if(MCE >=2){
        float drMCE = dRCalc(MCEta->at(MCE_tracker.at(0)), MCPhi->at(MCE_tracker.at(0)), MCEta->at(MCE_tracker.at(1)), MCPhi->at(MCE_tracker.at(1))); 
        if(drMCE > 0.5){
          vector<float> MCEJet;
          MCEJet.clear();
          for(int q=0; q<MCE; ++q){
            float mindrmcjet(0); 
            for(size_t s=0; s < JetEta->size(); ++s){
              float drMCEJet = dRCalc(MCEta->at(MCE_tracker.at(q)), MCPhi->at(MCE_tracker.at(q)), JetEta->at(s), JetPhi->at(s));
              MCEJet.push_back(drMCEJet);
            }
            if(MCEJet.size()!=0)mindrmcjet = *std::min_element(MCEJet.begin(),MCEJet.end());
            if(mindrmcjet!=0)h_Dr_MCElectronJets->Fill(mindrmcjet); 
          }
        }
      }
      
      
      
      
      if(NoCut>=2)EventCutFlow[0]++;
      if(NeutralCut>=2)EventCutFlow[1]++;
      if(ChargedCut>=2)EventCutFlow[2]++;
      if(PhotonIsoCut>=2){
        EventCutFlow[3]++;
        for(int k = 0; k < PhotonIsoCut; ++k){
          h_ShowerShape->Fill(PhoSigmaIEtaIEta->at(showershape_tracker.at(k))); 
        }
      }
      if(ShowerShapeCut>=2){
        EventCutFlow[4]++;
        if(ShowerShapeCut == 2){
          for(int k = 0; k < ShowerShapeCut; ++k){
            h_Pt->Fill(PhoEt->at(pt_tracker.at(k))); 
          }
        }
      }
      if(AfterPtCut>=2)EventCutFlow[5]++;
      if(OnlyPhoton>=2)EventCutFlow[6]++;
      if(HoverECut>=2)EventCutFlow[7]++;
      if(OnlyElectron>=2)EventOnlyElectron++;
      if(HoverECutElectron>=2)EventHEElectron++;
      // if (Cut(ientry) < 0) continue;
   }// event loop ends
   
   
   
   
   
   cout << "Double Photon events: " << EventCutFlow[8] << endl;
   cout << "DOuble Fake events: " << TotalFake << endl;
   cout << "Double Electron events: " << EventElectronDr << endl;
   cout << "Number of photons: " << totalphoton << endl;
   cout << "number of photons with pixel seed and csev mismatch: " << totalmismatchphoton << endl;
   for(int i = 0; i < 9; ++i){
     h_CutFlow_Photon->SetBinContent(i+1, EventCutFlow[i]);
     h_CutFlow_Photon->GetXaxis()->SetBinLabel(i+1,CutFlowLabels[i]);
     if(i > 5)continue;
     h_CutFlow_Electron->SetBinContent(i+1, EventCutFlow[i]);
     h_CutFlow_Electron->GetXaxis()->SetBinLabel(i+1,CutFlowLabels[i]);
   }
   h_CutFlow_Electron->SetBinContent(7, EventOnlyElectron);
   h_CutFlow_Electron->GetXaxis()->SetBinLabel(7, CutFlowLabels[6]);
   h_CutFlow_Electron->SetBinContent(8, EventHEElectron);
   h_CutFlow_Electron->GetXaxis()->SetBinLabel(8, CutFlowLabels[7]);
   h_CutFlow_Electron->SetBinContent(9, EventElectronDr);
   h_CutFlow_Electron->GetXaxis()->SetBinLabel(9, CutFlowLabels[8]);
     
   gROOT->SetBatch(true);
  
   
   
   
   f1->Write();
}