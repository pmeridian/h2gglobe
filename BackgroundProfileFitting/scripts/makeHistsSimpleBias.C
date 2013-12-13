#define makeHistsSimpleBias_cxx
#include "makeHistsSimpleBias.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void makeHistsSimpleBias::FillHist(TH1F* histo, string plotType)
{
  if (plotType=="mu")
    histo->Fill(mu);
  else if (plotType=="errMu")
    histo->Fill(sigma_mu);
  else if (plotType=="errBkg")
    histo->Fill(bkgErrSig1fwhm/bkgSig1fwhm);
  else if (plotType=="pullMu")
    histo->Fill(pull);
  else if (plotType=="pullBkg")
    histo->Fill( (bkgSig1fwhm-bkgTrue1fwhm)/bkgErrSig1fwhm);
}

TH1F* makeHistsSimpleBias::NewHist(TString name, string plotType)
{
  TH1F* histo=0;
  string myName(name.Data());
  histo=new TH1F(myName.c_str(),myName.c_str(),histos[plotType].nbins, histos[plotType].low, histos[plotType].high);
  std::cout << "[INFO]::Creating new histo " << myName << ":("<< histos[plotType].nbins << "," << histos[plotType].low << "," <<  histos[plotType].high << ")" << std::endl;
  return histo;
}

void makeHistsSimpleBias::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L makeHistsSimpleBias.C
//      Root > makeHistsSimpleBias t
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
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   std::map<TString,TH1F*> histos1D;

   for (int i=0;i<plot_vec.size();++i)
     {
       histo_properties thisHist;
       thisHist.nbins=nbins_vec[i];
       thisHist.low=low_vec[i];
       thisHist.high=high_vec[i];
       histos[plot_vec[i]]=thisHist;
     }

//    //Convert python list in a vector
//    boost::python::ssize_t len = boost::python::len(plot);
//    std::vector<string> plot_vec(len);
//    for(int i=0; i<len;i++)
//      plot_vec[i]=boost::python::extract<string>(plot[i]);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry%10000==0)
       std::cout << jentry << std::endl;
     if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //Check quality of the fit
//       if (sigma_mu<=0.1)
// 	continue;
//       if (fabs(fabs(mu)-10)<=0.1)
// 	continue;

//       if (sigma_mu>10)
// 	continue;
//       if (TMath::Abs(mu-muTruth)>8)
// 	continue;

      TString evtType=Form("cat%d_truth_%s_test_%s_muInj_%2.1f_muConstr_%3.2f_mass_%d",int(mycat),genFun,fitFun,muTruth,sigma_mu_constraint,mass);
      for (histo_map::const_iterator it=histos.begin();it!=histos.end();++it)
	{
	  TString key=evtType+Form("_%s",it->first.c_str());
	  if (histos1D.count(key)>0)
	      FillHist(&(*histos1D[key]),it->first);
	  else
	    {
	      histos1D[key]=NewHist(key,it->first);
	      FillHist(&(*histos1D[key]),it->first);
	    }
	}
      // if (Cut(ientry) < 0) continue;
   }

   TFile* f=TFile::Open(outFile.c_str(),"RECREATE");
   f->cd();
   for(std::map<TString,TH1F*>::const_iterator histo=histos1D.begin();histo!=histos1D.end();++histo)
     histo->second->Write();
   f->Close();

   //Clear Map
   for(std::map<TString,TH1F*>::const_iterator histo=histos1D.begin();histo!=histos1D.end();++histo)
     delete histo->second;
   
}
