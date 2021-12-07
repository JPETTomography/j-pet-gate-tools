#ifndef Hits_h
#define Hits_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

#include "Hit.h"

struct LoopResults {
  LoopResults() : real_time(0.),
                  multiplicities(vector<int>()),
                  counter_all_compton_hits(0),
                  counter_compton_hits_over_the_ETH0(0),
                  counter_compton_hits_over_the_ETH(0) {}
  double real_time;
  vector<int> multiplicities;
  int counter_all_compton_hits;
  int counter_compton_hits_over_the_ETH0;
  int counter_compton_hits_over_the_ETH;
};

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Hits {
public :

   std::string tree_name = std::string(getenv("GOJA_TREE_NAME"));

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           PDGEncoding;
   Double_t        time;
   Float_t         edep;
   Float_t         posX;
   Float_t         posY;
   Float_t         posZ;
   Int_t           rsectorID;
   Int_t           layerID;
   Int_t           nPhantomCompton;
   Int_t           nCrystalCompton;
   Int_t           nPhantomRayleigh;
   Int_t           nCrystalRayleigh;
   Float_t         sourcePosX;
   Float_t         sourcePosY;
   Float_t         sourcePosZ;
   Int_t           sourceID;
   Int_t           eventID;
   Int_t           runID;
   Float_t         axialPos;
   Float_t         rotationAngle;
   Int_t           volumeID[10];
   Char_t          processName[20];
   Char_t          comptVolName[5];
   Char_t          RayleighVolName[5];

   // List of branches
   TBranch        *b_PDGEncoding;   //!
   TBranch        *b_time;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_posX;   //!
   TBranch        *b_posY;   //!
   TBranch        *b_posZ;   //!
   TBranch        *b_rsectorID;
   TBranch        *b_layerID;   //!
   TBranch        *b_nPhantomCompton;   //!
   TBranch        *b_nCrystalCompton;   //!
   TBranch        *b_nPhantomRayleigh;   //!
   TBranch        *b_nCrystalRayleigh;   //!
   TBranch        *b_sourcePosX;   //!
   TBranch        *b_sourcePosY;   //!
   TBranch        *b_sourcePosZ;   //!
   TBranch        *b_sourceID;   //!
   TBranch        *b_eventID;   //!
   TBranch        *b_runID;   //!
   TBranch        *b_axialPos;   //!
   TBranch        *b_rotationAngle;   //!
   TBranch        *b_volumeID;   //!
   TBranch        *b_processName;   //!
   TBranch        *b_comptVolName;   //!
   TBranch        *b_RayleighVolName;   //!

   Hits(TTree *tree=0);
   virtual ~Hits();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual LoopResults Loop(bool singles);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Hits_cxx
Hits::Hits(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject(tree_name.c_str(),tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain(tree_name.c_str(),"");
      char* root_filename;
      root_filename = getenv("GOJA_ROOT_FILENAME");
      std::string root_filename_string = std::string(root_filename) + std::string("/") + tree_name;
      chain->Add(root_filename_string.c_str());

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

Hits::~Hits()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Hits::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Hits::LoadTree(Long64_t entry)
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

void Hits::Init(TTree *tree)
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

   if (tree_name == "Hits") {
     fChain->SetBranchAddress("edep", &edep, &b_edep);
     fChain->SetBranchAddress("posX", &posX, &b_posX);
     fChain->SetBranchAddress("posY", &posY, &b_posY);
     fChain->SetBranchAddress("posZ", &posZ, &b_posZ);
     fChain->SetBranchAddress("nPhantomCompton", &nPhantomCompton, &b_nPhantomCompton);
     fChain->SetBranchAddress("nCrystalCompton", &nCrystalCompton, &b_nCrystalCompton);
     fChain->SetBranchAddress("nPhantomRayleigh", &nPhantomRayleigh, &b_nPhantomRayleigh);
     fChain->SetBranchAddress("nCrystalRayleigh", &nCrystalRayleigh, &b_nCrystalRayleigh);
     fChain->SetBranchAddress("PDGEncoding", &PDGEncoding, &b_PDGEncoding);
     fChain->SetBranchAddress("processName", processName, &b_processName);
     fChain->SetBranchAddress("volumeID", volumeID, &b_volumeID);
   }
   else if (tree_name == "HESingles") {
     fChain->SetBranchAddress("energy", &edep, &b_edep);
     fChain->SetBranchAddress("globalPosX", &posX, &b_posX);
     fChain->SetBranchAddress("globalPosY", &posY, &b_posY);
     fChain->SetBranchAddress("globalPosZ", &posZ, &b_posZ);
     fChain->SetBranchAddress("comptonPhantom", &nPhantomCompton, &b_nPhantomCompton);
     fChain->SetBranchAddress("comptonCrystal", &nCrystalCompton, &b_nCrystalCompton);
     fChain->SetBranchAddress("RayleighPhantom", &nPhantomRayleigh, &b_nPhantomRayleigh);
     fChain->SetBranchAddress("RayleighCrystal", &nCrystalRayleigh, &b_nCrystalRayleigh);
   }
   fChain->SetBranchAddress("time", &time, &b_time);
   string systemType = string(getenv("GOJA_SYSTEM_TYPE"));
   if (systemType == "cylindricalPET")
     fChain->SetBranchAddress("rsectorID", &rsectorID, &b_rsectorID);
   fChain->SetBranchAddress("layerID", &layerID, &b_layerID);
   fChain->SetBranchAddress("sourcePosX", &sourcePosX, &b_sourcePosX);
   fChain->SetBranchAddress("sourcePosY", &sourcePosY, &b_sourcePosY);
   fChain->SetBranchAddress("sourcePosZ", &sourcePosZ, &b_sourcePosZ);
   fChain->SetBranchAddress("sourceID", &sourceID, &b_sourceID);
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("axialPos", &axialPos, &b_axialPos);
   fChain->SetBranchAddress("rotationAngle", &rotationAngle, &b_rotationAngle);
   fChain->SetBranchAddress("comptVolName", comptVolName, &b_comptVolName);
   fChain->SetBranchAddress("RayleighVolName", RayleighVolName, &b_RayleighVolName);
   Notify();
}

Bool_t Hits::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Hits::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

#endif // #ifdef Hits_cxx
