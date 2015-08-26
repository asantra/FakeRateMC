//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 13 01:00:20 2015 by ROOT version 5.34/18
// from TTree tree/WJets SUSY tree
// found on file: /eos/uscms/store/user/asantra4/MyNtuples_Phys14/WJets_SmallNtupleQQQ.root
//////////////////////////////////////////////////////////

#ifndef SmallNtupleQQQ_h
#define SmallNtupleQQQ_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class SmallNtupleQQQ {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Lumis;
   Long64_t        Event;
   Int_t           NVtx;
   Int_t           NTrks;
   Float_t         Rho;
   vector<int>     *MCPID;
   vector<float>   *MCVtx_x;
   vector<float>   *MCVtx_y;
   vector<float>   *MCVtx_z;
   vector<float>   *MCPt;
   vector<float>   *MCMass;
   vector<float>   *MCEta;
   vector<float>   *MCPhi;
   vector<float>   *MCStatus;
   vector<int>     *MCGMomPID;
   vector<int>     *MCMomPID;
   vector<float>   *MCMomPt;
   vector<float>   *MCMomMass;
   vector<float>   *MCMomEta;
   vector<float>   *MCMomPhi;
   Int_t           NPUInfo;
   vector<int>     *NPU;
   vector<int>     *PUBX;
   vector<float>   *PUTrue;
   Float_t         PFMET;
   Float_t         PFMETPhi;
   Float_t         PFMETsumEt;
   Float_t         PFMETmEtSig;
   Float_t         PFMETSig;
   vector<float>   *ElePt;
   vector<float>   *EleEta;
   vector<float>   *ElePhi;
   vector<float>   *EleD0;
   vector<float>   *EleDz;
   vector<float>   *PhoE;
   vector<float>   *PhoEt;
   vector<float>   *PhoEta;
   vector<float>   *PhoPhi;
   vector<int>     *PhohasPixelSeed;
   vector<int>     *PhoEleVeto;
   vector<float>   *PhoR9;
   vector<float>   *PhoHoverE;
   vector<float>   *PhoSigmaIEtaIEta;
   vector<float>   *PhoSigmaIEtaIPhi;
   vector<float>   *PhoSigmaIPhiIPhi;
   //vector<float>   *PhoSigmaIEtaIEta_2012;
   //vector<float>   *PhoSigmaIEtaIPhi_2012;
   //vector<float>   *PhoSigmaIPhiIPhi_2012;
   vector<float>   *PhoPFChIso;
   vector<float>   *PhoPFPhoIso;
   vector<float>   *PhoPFNeuIso;
   vector<float>   *PhoPFChWorstIso;
   vector<float>   *MuPt;
   vector<float>   *MuEta;
   vector<float>   *MuPhi;
   vector<int>     *MuCharge;
   vector<int>     *MuType;
   vector<int>     *MuIsGood;
   vector<float>   *MuD0;
   vector<float>   *MuDz;
   vector<float>   *MuChi2NDF;
   vector<int>     *MuTrkLayers;
   vector<int>     *MuPixelLayers;
   vector<int>     *MuPixelHits;
   vector<int>     *MuMuonHits;
   vector<int>     *MuStations;
   vector<int>     *MuTrkQuality;
   vector<float>   *MuIsoTrk;
   vector<float>   *MuPFChIso;
   vector<float>   *MuPFPhoIso;
   vector<float>   *MuPFNeuIso;
   vector<float>   *MuPFPUIso;
   vector<float>   *JetPt;
   vector<float>   *JetEta;
   vector<float>   *JetPhi;
   vector<float>   *JetCHF;
   vector<float>   *JetNHF;
   vector<float>   *JetCEF;
   vector<float>   *JetNEF;
   vector<int>     *JetNCH;
   vector<float>   *JetHFHAE;
   vector<float>   *JetHFEME;
   vector<int>     *JetNConstituents;
   vector<int>     *JetPartonID;
   vector<bool>    *JetPFLooseId;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Lumis;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_NTrks;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_MCPID;   //!
   TBranch        *b_MCVtx_x;   //!
   TBranch        *b_MCVtx_y;   //!
   TBranch        *b_MCVtx_z;   //!
   TBranch        *b_MCPt;   //!
   TBranch        *b_MCMass;   //!
   TBranch        *b_MCEta;   //!
   TBranch        *b_MCPhi;   //!
   TBranch        *b_MCStatus; //!
   TBranch        *b_MCGMomPID;   //!
   TBranch        *b_MCMomPID;   //!
   TBranch        *b_MCMomPt;   //!
   TBranch        *b_MCMomMass;   //!
   TBranch        *b_MCMomEta;   //!
   TBranch        *b_MCMomPhi;   //!
   TBranch        *b_NPUInfo;   //!
   TBranch        *b_NPU;   //!
   TBranch        *b_PUBX;   //!
   TBranch        *b_PUTrue;   //!
   TBranch        *b_PFMET;   //!
   TBranch        *b_PFMETPhi;   //!
   TBranch        *b_PFMETsumEt;   //!
   TBranch        *b_PFMETmEtSig;   //!
   TBranch        *b_PFMETSig;   //!
   TBranch        *b_ElePt; //!
   TBranch        *b_EleEta; //!
   TBranch        *b_ElePhi; //!
   TBranch        *b_EleD0; //!
   TBranch        *b_EleDz; //!
   TBranch        *b_PhoE;   //!
   TBranch        *b_PhoEt;   //!
   TBranch        *b_PhoEta;   //!
   TBranch        *b_PhoPhi;   //!
   TBranch        *b_PhohasPixelSeed;   //!
   TBranch        *b_PhoEleVeto;   //!
   TBranch        *b_PhoR9;   //!
   TBranch        *b_PhoHoverE;   //!
   TBranch        *b_PhoSigmaIEtaIEta;   //!
   TBranch        *b_PhoSigmaIEtaIPhi;   //!
   TBranch        *b_PhoSigmaIPhiIPhi;   //!
   //TBranch        *b_PhoSigmaIEtaIEta_2012;   //!
   //TBranch        *b_PhoSigmaIEtaIPhi_2012;   //!
   //TBranch        *b_PhoSigmaIPhiIPhi_2012;   //!
   TBranch        *b_PhoPFChIso;   //!
   TBranch        *b_PhoPFPhoIso;   //!
   TBranch        *b_PhoPFNeuIso;   //!
   TBranch        *b_PhoPFChWorstIso;   //!
   TBranch        *b_MuPt; //!
   TBranch        *b_MuEta; //!
   TBranch        *b_MuPhi; //!
   TBranch        *b_MuCharge; //!
   TBranch        *b_MuType; //!
   TBranch        *b_MuIsGood; //!
   TBranch        *b_MuD0; //!
   TBranch        *b_MuDz; //!
   TBranch        *b_MuChi2NDF; //!
   TBranch        *b_MuTrkLayers; //!
   TBranch        *b_MuPixelLayers; //!
   TBranch        *b_MuPixelHits; //!
   TBranch        *b_MuMuonHits; //!
   TBranch        *b_MuStations; //!
   TBranch        *b_MuTrkQuality; //!
   TBranch        *b_MuIsoTrk; //!
   TBranch        *b_MuPFChIso;
   TBranch        *b_MuPFPhoIso;
   TBranch        *b_MuPFNeuIso;
   TBranch        *b_MuPFPUIso;
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetCHF;   //!
   TBranch        *b_JetNHF;   //!
   TBranch        *b_JetCEF;   //!
   TBranch        *b_JetNEF;   //!
   TBranch        *b_JetNCH;   //!
   TBranch        *b_JetHFHAE;   //!
   TBranch        *b_JetHFEME;   //!
   TBranch        *b_JetNConstituents;   //!
   TBranch        *b_JetPartonID;
   TBranch        *b_JetPFLooseId;

   SmallNtupleQQQ(TTree *tree=0);
   virtual ~SmallNtupleQQQ();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SmallNtupleQQQ_cxx
SmallNtupleQQQ::SmallNtupleQQQ(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TString fname = "PPP";
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
      if (!f || !f->IsOpen()) {
         f = new TFile(fname);
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

SmallNtupleQQQ::~SmallNtupleQQQ()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SmallNtupleQQQ::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SmallNtupleQQQ::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SmallNtupleQQQ::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   MCPID = 0;
   MCVtx_x = 0;
   MCVtx_y = 0;
   MCVtx_z = 0;
   MCPt = 0;
   MCMass = 0;
   MCEta = 0;
   MCPhi = 0;
   MCStatus = 0;
   MCGMomPID = 0;
   MCMomPID = 0;
   MCMomPt = 0;
   MCMomMass = 0;
   MCMomEta = 0;
   MCMomPhi = 0;
   NPU = 0;
   PUBX = 0;
   PUTrue = 0;
   ElePt = 0;
   EleEta = 0;
   ElePhi = 0;
   EleD0 = 0;
   EleDz = 0;
   PhoE = 0;
   PhoEt = 0;
   PhoEta = 0;
   PhoPhi = 0;
   PhohasPixelSeed = 0;
   PhoEleVeto = 0;
   PhoR9 = 0;
   PhoHoverE = 0;
   PhoSigmaIEtaIEta = 0;
   PhoSigmaIEtaIPhi = 0;
   PhoSigmaIPhiIPhi = 0;
   //PhoSigmaIEtaIEta_2012 = 0;
   //PhoSigmaIEtaIPhi_2012 = 0;
   //PhoSigmaIPhiIPhi_2012 = 0;
   PhoPFChIso = 0;
   PhoPFPhoIso = 0;
   PhoPFNeuIso = 0;
   PhoPFChWorstIso = 0;
   MuPt = 0;
   MuEta = 0;
   MuPhi = 0;
   MuCharge = 0;
   MuType = 0;
   MuIsGood = 0;
   MuD0 = 0;
   MuDz = 0;
   MuChi2NDF = 0;
   MuTrkLayers = 0;
   MuPixelLayers = 0;
   MuPixelHits = 0;
   MuMuonHits = 0;
   MuStations = 0;
   MuTrkQuality = 0;
   MuIsoTrk = 0;
   MuPFChIso = 0;
   MuPFPhoIso = 0;
   MuPFNeuIso = 0;
   MuPFPUIso = 0;
   JetPt = 0;
   JetEta = 0;
   JetPhi = 0;
   JetCHF = 0;
   JetNHF = 0;
   JetCEF = 0;
   JetNEF = 0;
   JetNCH = 0;
   JetHFHAE = 0;
   JetHFEME = 0;
   JetNConstituents = 0;
   JetPartonID = 0;
   JetPFLooseId = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Lumis", &Lumis, &b_Lumis);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
   fChain->SetBranchAddress("NTrks", &NTrks, &b_NTrks);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("MCPID", &MCPID, &b_MCPID);
   fChain->SetBranchAddress("MCVtx_x", &MCVtx_x, &b_MCVtx_x);
   fChain->SetBranchAddress("MCVtx_y", &MCVtx_y, &b_MCVtx_y);
   fChain->SetBranchAddress("MCVtx_z", &MCVtx_z, &b_MCVtx_z);
   fChain->SetBranchAddress("MCPt", &MCPt, &b_MCPt);
   fChain->SetBranchAddress("MCMass", &MCMass, &b_MCMass);
   fChain->SetBranchAddress("MCEta", &MCEta, &b_MCEta);
   fChain->SetBranchAddress("MCPhi", &MCPhi, &b_MCPhi);
   fChain->SetBranchAddress("MCStatus", &MCStatus, &b_MCStatus);
   fChain->SetBranchAddress("MCGMomPID", &MCGMomPID, &b_MCGMomPID);
   fChain->SetBranchAddress("MCMomPID", &MCMomPID, &b_MCMomPID);
   fChain->SetBranchAddress("MCMomPt", &MCMomPt, &b_MCMomPt);
   fChain->SetBranchAddress("MCMomMass", &MCMomMass, &b_MCMomMass);
   fChain->SetBranchAddress("MCMomEta", &MCMomEta, &b_MCMomEta);
   fChain->SetBranchAddress("MCMomPhi", &MCMomPhi, &b_MCMomPhi);
   fChain->SetBranchAddress("NPUInfo", &NPUInfo, &b_NPUInfo);
   fChain->SetBranchAddress("NPU", &NPU, &b_NPU);
   fChain->SetBranchAddress("PUBX", &PUBX, &b_PUBX);
   fChain->SetBranchAddress("PUTrue", &PUTrue, &b_PUTrue);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   fChain->SetBranchAddress("PFMETPhi", &PFMETPhi, &b_PFMETPhi);
   fChain->SetBranchAddress("PFMETsumEt", &PFMETsumEt, &b_PFMETsumEt);
   fChain->SetBranchAddress("PFMETmEtSig", &PFMETmEtSig, &b_PFMETmEtSig);
   fChain->SetBranchAddress("PFMETSig", &PFMETSig, &b_PFMETSig);
   fChain->SetBranchAddress("ElePt", &ElePt, &b_ElePt);
   fChain->SetBranchAddress("EleEta", &EleEta, &b_EleEta);
   fChain->SetBranchAddress("ElePhi", &ElePhi, &b_ElePhi);
   fChain->SetBranchAddress("EleD0", &EleD0, &b_EleD0);
   fChain->SetBranchAddress("EleDz", &EleDz, &b_EleDz);
   fChain->SetBranchAddress("PhoE", &PhoE, &b_PhoE);
   fChain->SetBranchAddress("PhoEt", &PhoEt, &b_PhoEt);
   fChain->SetBranchAddress("PhoEta", &PhoEta, &b_PhoEta);
   fChain->SetBranchAddress("PhoPhi", &PhoPhi, &b_PhoPhi);
   fChain->SetBranchAddress("PhohasPixelSeed", &PhohasPixelSeed, &b_PhohasPixelSeed);
   fChain->SetBranchAddress("PhoEleVeto", &PhoEleVeto, &b_PhoEleVeto);
   fChain->SetBranchAddress("PhoR9", &PhoR9, &b_PhoR9);
   fChain->SetBranchAddress("PhoHoverE", &PhoHoverE, &b_PhoHoverE);
   fChain->SetBranchAddress("PhoSigmaIEtaIEta", &PhoSigmaIEtaIEta, &b_PhoSigmaIEtaIEta);
   fChain->SetBranchAddress("PhoSigmaIEtaIPhi", &PhoSigmaIEtaIPhi, &b_PhoSigmaIEtaIPhi);
   fChain->SetBranchAddress("PhoSigmaIPhiIPhi", &PhoSigmaIPhiIPhi, &b_PhoSigmaIPhiIPhi);
   //fChain->SetBranchAddress("PhoSigmaIEtaIEta_2012", &PhoSigmaIEtaIEta_2012, &b_PhoSigmaIEtaIEta_2012);
   //fChain->SetBranchAddress("PhoSigmaIEtaIPhi_2012", &PhoSigmaIEtaIPhi_2012, &b_PhoSigmaIEtaIPhi_2012);
   //fChain->SetBranchAddress("PhoSigmaIPhiIPhi_2012", &PhoSigmaIPhiIPhi_2012, &b_PhoSigmaIPhiIPhi_2012);
   fChain->SetBranchAddress("PhoPFChIso", &PhoPFChIso, &b_PhoPFChIso);
   fChain->SetBranchAddress("PhoPFPhoIso", &PhoPFPhoIso, &b_PhoPFPhoIso);
   fChain->SetBranchAddress("PhoPFNeuIso", &PhoPFNeuIso, &b_PhoPFNeuIso);
   fChain->SetBranchAddress("PhoPFChWorstIso", &PhoPFChWorstIso, &b_PhoPFChWorstIso);
   fChain->SetBranchAddress("MuPt", &MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEta", &MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", &MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuCharge", &MuCharge, &b_MuCharge);
   fChain->SetBranchAddress("MuType", &MuType, &b_MuType);
   fChain->SetBranchAddress("MuIsGood", &MuIsGood, &b_MuIsGood);
   fChain->SetBranchAddress("MuD0", &MuD0, &b_MuD0);
   fChain->SetBranchAddress("MuDz", &MuDz, &b_MuDz);
   fChain->SetBranchAddress("MuChi2NDF", &MuChi2NDF, &b_MuChi2NDF);
   fChain->SetBranchAddress("MuTrkLayers", &MuTrkLayers,  &b_MuTrkLayers);
   fChain->SetBranchAddress("MuPixelLayers", &MuPixelLayers, &b_MuPixelLayers);
   fChain->SetBranchAddress("MuPixelHits", &MuPixelHits, &b_MuPixelHits);
   fChain->SetBranchAddress("MuMuonHits", &MuMuonHits, &b_MuMuonHits);
   fChain->SetBranchAddress("MuStations", &MuStations, &b_MuStations);
   fChain->SetBranchAddress("MuTrkQuality", &MuTrkQuality, &b_MuTrkQuality);
   fChain->SetBranchAddress("MuIsoTrk", &MuIsoTrk, &b_MuIsoTrk);
   fChain->SetBranchAddress("MuPFChIso", &MuPFChIso, &b_MuPFChIso);
   fChain->SetBranchAddress("MuPFPhoIso", &MuPFPhoIso, &b_MuPFPhoIso);
   fChain->SetBranchAddress("MuPFNeuIso", &MuPFNeuIso, &b_MuPFNeuIso);
   fChain->SetBranchAddress("MuPFPUIso", &MuPFPUIso, &b_MuPFPUIso);
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetCHF", &JetCHF, &b_JetCHF);
   fChain->SetBranchAddress("JetNHF", &JetNHF, &b_JetNHF);
   fChain->SetBranchAddress("JetCEF", &JetCEF, &b_JetCEF);
   fChain->SetBranchAddress("JetNEF", &JetNEF, &b_JetNEF);
   fChain->SetBranchAddress("JetNCH", &JetNCH, &b_JetNCH);
   fChain->SetBranchAddress("JetHFHAE", &JetHFHAE, &b_JetHFHAE);
   fChain->SetBranchAddress("JetHFEME", &JetHFEME, &b_JetHFEME);
   fChain->SetBranchAddress("JetNConstituents", &JetNConstituents, &b_JetNConstituents);
   fChain->SetBranchAddress("JetPartonID", &JetPartonID, &b_JetPartonID);
   fChain->SetBranchAddress("JetPFLooseId", & JetPFLooseId, &b_JetPFLooseId);
   Notify();
}

Bool_t SmallNtupleQQQ::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SmallNtupleQQQ::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SmallNtupleQQQ::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SmallNtupleQQQ_cxx
