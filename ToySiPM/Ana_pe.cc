{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Tue May 11 17:04:00 2010 by ROOT version5.22/00a)
//   from TTree t1/t1
//   found on file: 10Watt_nolead.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();

   TH1F * H0 = new TH1F("H0","  L0 Pes ",1000,0.,1000.);
 
 
   TH1F * H1 = new TH1F("H1","  L1 Pes ",1000,0.,1000.);
 
 
   TH1F * H2 = new TH1F("H2","  L2 Pes ",1000,0.,1000.);
 
 
   TH1F * H3 = new TH1F("H3","  L3 Pes ",1000,0.,1000.);
 
 
   TH1F * H4 = new TH1F("H4","  L4 Pes ",1000,0.,1000.);
 
 
   TH1F * H5 = new TH1F("H5","  L5 Pes ",1000,0.,1000.);
 
 
   TH1F * H6 = new TH1F("H6","  L6 Pes ",1000,0.,1000.);
 
 
   TH1F * H7 = new TH1F("H7","  L7 Pes ",1000,0.,1000.);
 
 
   TH1F * H9 = new TH1F("H8","  L8 Pes ",1000,0.,1000.);
 
 
   TH1F * H9 = new TH1F("H9","  L9 Pes ",1000,0.,1000.);
 
 

   TFile outf("../temp.root","recreate");
   TTree *t1 = new TTree("t1","t1");







   t1 -> ReadFile("LayerPes500GeV.txt","slices[19]/D");



      //
   t1 -> Write();

   t1->Print();

//Declaration of leaves types

   Double_t           slices[19];
   Double_t weight;
   Double_t sum;
 

   // Set branch addresses.

   t1->SetBranchAddress("slices[19]",slices);



//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// t1->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Int_t icount = 0;

   Long64_t nentries = t1->GetEntries();

   Long64_t nbytes = 0;
   TCanvas * c1 = new TCanvas;

//Main loop over events in ntuple
//
   for (Long64_t i=0; i<nentries;i++) {
     nbytes += t1->GetEntry(i);
     
     if(i==0) {  cout << i  << endl; }
     
     
     for (Int_t ibin=0; ibin<18;ibin++) {
       
       H1 -> AddBinContent(ibin,weight);
       
     } // end of loop of bins (ibin)
     


   }// end of loop over events (i)


   outf.Close();

} // end of routine
