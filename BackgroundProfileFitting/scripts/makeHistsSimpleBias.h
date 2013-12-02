//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec  1 09:34:15 2013 by ROOT version 5.34/03
// from TTree muTree/muTree
// found on file: root://eoscms//eos/cms/store/group/phys_higgs/meridian/BiasStudy/legacy_freeze_v2_massFac/cat0_mu0.0_mass110/StandardBiasStudyOut_cat0_job0.root
//////////////////////////////////////////////////////////

#ifndef makeHistsSimpleBias_h
#define makeHistsSimpleBias_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

//#include <boost/python.hpp>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class makeHistsSimpleBias {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //   boost::python::list plot;
   struct histo_properties
   {
     int nbins;
     float low;
     float high;
   };

   std::vector<string> plot_vec;
   std::vector<int> nbins_vec;
   std::vector<float> low_vec;
   std::vector<float> high_vec;

   int mycat;

   typedef std::map<string,histo_properties> histo_map;
   histo_map histos;

   string outFile;

   void FillHist(TH1F* histo, string plotType);
   TH1F* NewHist(TString name, string plotType);

   // Declaration of leaf types
   Int_t           mass;
   Float_t         mu;
   UInt_t          cat;
   Float_t         muTruth;
   Float_t         sigma_mu;
   Float_t         bkgTrue1fwhm;
   Float_t         bkgTrue2fwhm;
   Float_t         bkgSig1fwhm;
   Float_t         bkgSig2fwhm;
   Float_t         bkgErrSig1fwhm;
   Float_t         bkgErrSig2fwhm;
   Float_t         bkgErrNormSig1fwhm;
   Float_t         bkgErrNormSig2fwhm;
   Char_t          genFun[10];
   Char_t          fitFun[10];
   UInt_t          iToy;
   Int_t           fitStatus;

   // List of branches
   TBranch        *b_mass;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_cat;   //!
   TBranch        *b_muTruth;   //!
   TBranch        *b_sigma_mu;   //!
   TBranch        *b_bkgTrue1fwhm;   //!
   TBranch        *b_bkgTrue2fwhm;   //!
   TBranch        *b_bkgSig1fwhm;   //!
   TBranch        *b_bkgSig2fwhm;   //!
   TBranch        *b_bkgErrSig1fwhm;   //!
   TBranch        *b_bkgErrSig2fwhm;   //!
   TBranch        *b_bkgErrNormSig1fwhm;   //!
   TBranch        *b_bkgErrNormSig2fwhm;   //!
   TBranch        *b_genFun;   //!
   TBranch        *b_fitFun;   //!
   TBranch        *b_iToy;   //!
   TBranch        *b_fitStatus;   //!

   makeHistsSimpleBias(TTree *tree=0);
   virtual ~makeHistsSimpleBias();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef makeHistsSimpleBias_cxx
makeHistsSimpleBias::makeHistsSimpleBias(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_higgs/meridian/BiasStudy/legacy_freeze_v2_massFac/cat0_mu0.0_mass110/StandardBiasStudyOut_cat0_job0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_higgs/meridian/BiasStudy/legacy_freeze_v2_massFac/cat0_mu0.0_mass110/StandardBiasStudyOut_cat0_job0.root");
      }
      f->GetObject("muTree",tree);

   }
   Init(tree);
}

makeHistsSimpleBias::~makeHistsSimpleBias()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t makeHistsSimpleBias::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t makeHistsSimpleBias::LoadTree(Long64_t entry)
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

void makeHistsSimpleBias::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("cat", &cat, &b_cat);
   fChain->SetBranchAddress("muTruth", &muTruth, &b_muTruth);
   fChain->SetBranchAddress("sigma_mu", &sigma_mu, &b_sigma_mu);
   fChain->SetBranchAddress("bkgTrue1fwhm", &bkgTrue1fwhm, &b_bkgTrue1fwhm);
   fChain->SetBranchAddress("bkgTrue2fwhm", &bkgTrue2fwhm, &b_bkgTrue2fwhm);
   fChain->SetBranchAddress("bkgSig1fwhm", &bkgSig1fwhm, &b_bkgSig1fwhm);
   fChain->SetBranchAddress("bkgSig2fwhm", &bkgSig2fwhm, &b_bkgSig2fwhm);
   fChain->SetBranchAddress("bkgErrSig1fwhm", &bkgErrSig1fwhm, &b_bkgErrSig1fwhm);
   fChain->SetBranchAddress("bkgErrSig2fwhm", &bkgErrSig2fwhm, &b_bkgErrSig2fwhm);
   fChain->SetBranchAddress("bkgErrNormSig1fwhm", &bkgErrNormSig1fwhm, &b_bkgErrNormSig1fwhm);
   fChain->SetBranchAddress("bkgErrNormSig2fwhm", &bkgErrNormSig2fwhm, &b_bkgErrNormSig2fwhm);
   fChain->SetBranchAddress("genFun", genFun, &b_genFun);
   fChain->SetBranchAddress("fitFun", fitFun, &b_fitFun);
   fChain->SetBranchAddress("iToy", &iToy, &b_iToy);
   fChain->SetBranchAddress("fitStatus", &fitStatus, &b_fitStatus);
   Notify();
}

Bool_t makeHistsSimpleBias::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void makeHistsSimpleBias::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t makeHistsSimpleBias::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef makeHistsSimpleBias_cxx
