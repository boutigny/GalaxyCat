//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 15 15:01:30 2013 by ROOT version 5.34/10
// from TTree Risa/Risa
// found on file: RW_truth_190.root
//////////////////////////////////////////////////////////

#ifndef Risa_h
#define Risa_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Risa {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           id;
   Int_t           index;
   Float_t         ra;
   Float_t         dec;
   Float_t         z;
   Float_t         gamma1;
   Float_t         gamma2;
   Float_t         kappa;
   Float_t         size;
   Float_t         eps1;
   Float_t         eps2;
   Float_t         mag;
   Float_t         teps1;
   Float_t         teps2;
   Float_t         tra;
   Float_t         tdec;
   Float_t         mu;
   Float_t         tsize;

   // List of branches
   TBranch        *b_id;   //!
   TBranch        *b_index;   //!
   TBranch        *b_ra;   //!
   TBranch        *b_dec;   //!
   TBranch        *b_z;   //!
   TBranch        *b_gamma1;   //!
   TBranch        *b_gamma2;   //!
   TBranch        *b_kappa;   //!
   TBranch        *b_size;   //!
   TBranch        *b_eps1;   //!
   TBranch        *b_eps2;   //!
   TBranch        *b_mag;   //!
   TBranch        *b_teps1;   //!
   TBranch        *b_teps2;   //!
   TBranch        *b_tra;   //!
   TBranch        *b_tdec;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_tsize;   //!

   Risa(TTree *tree=0);
   virtual ~Risa();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
//   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Risa_cxx
Risa::Risa(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RW_truth_190.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RW_truth_190.root");
      }
      f->GetObject("Risa",tree);

   }
   Init(tree);
}

Risa::~Risa()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Risa::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Risa::LoadTree(Long64_t entry)
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

void Risa::Init(TTree *tree)
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

   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("index", &index, &b_index);
   fChain->SetBranchAddress("ra", &ra, &b_ra);
   fChain->SetBranchAddress("dec", &dec, &b_dec);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("gamma1", &gamma1, &b_gamma1);
   fChain->SetBranchAddress("gamma2", &gamma2, &b_gamma2);
   fChain->SetBranchAddress("kappa", &kappa, &b_kappa);
   fChain->SetBranchAddress("size", &size, &b_size);
   fChain->SetBranchAddress("eps1", &eps1, &b_eps1);
   fChain->SetBranchAddress("eps2", &eps2, &b_eps2);
   fChain->SetBranchAddress("mag", &mag, &b_mag);
   fChain->SetBranchAddress("teps1", &teps1, &b_teps1);
   fChain->SetBranchAddress("teps2", &teps2, &b_teps2);
   fChain->SetBranchAddress("tra", &tra, &b_tra);
   fChain->SetBranchAddress("tdec", &tdec, &b_tdec);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("tsize", &tsize, &b_tsize);
   Notify();
}

Bool_t Risa::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Risa::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Risa::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Risa_cxx
