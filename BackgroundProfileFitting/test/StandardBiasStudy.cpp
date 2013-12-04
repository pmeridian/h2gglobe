//Designed for individual categories

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"

#include "../interface/ProfileMultiplePdfs.h"
#include "../interface/PdfModelBuilder.h"

using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = boost::program_options;

void readDatFile(string datFileName, int cat, vector<pair<int,pair<string,string> > > &toysMap, vector<pair<int,pair<string,string> > > &testMap){
 
  ifstream datfile;
  datfile.open(datFileName.c_str());
  bool foundCat=false;
  if (datfile.fail()) return;
  while (datfile.good()){
    string line;
    getline(datfile,line);
    if (line=="\n" || line.substr(0,1)=="#" || line==" " || line.empty()) continue;
    if (line.substr(0,line.find("="))=="cat") {
      foundCat=true;
      if (atoi((line.substr(line.find("=")+1,string::npos)).c_str())!=cat) foundCat=false;
    }
    if (!foundCat) continue;
    string type = line.substr(0,line.find("="));
    string pdfs = line.substr(line.find("=")+1,string::npos);
    
    if (type=="truth" && starts_with(pdfs,"Hybrid")){
      vector<string> els;
      split(els,pdfs,boost::is_any_of(":"));
      string masses=els[1];
      string funcs=els[2];
      toysMap.push_back(pair<int,pair<string,string> >(-1,make_pair(masses,funcs)));
      continue;
    }
    if (type=="truth" && starts_with(pdfs,"KeysPdf")){
      vector<string> els;
      split(els,pdfs,boost::is_any_of(":"));
      string rho=els[1];
      string name=els[2];
      toysMap.push_back(pair<int,pair<string,string> > (-2,make_pair(rho,name)));
      continue;
    }
    if (type=="truth" && starts_with(pdfs,"File")){
      vector<string> els;
      split(els,pdfs,boost::is_any_of(":"));
      toysMap.push_back(pair<int,pair<string,string> > (-3,make_pair(els[2],els[0])));
      continue;
    }

    if (type=="truth" || type=="test" ){
      vector<string> els;
      split(els,pdfs,boost::is_any_of(":"));
      int order=atoi(els[1].c_str());
      string title=els[0];
      string name=els[2];
      if (type=="truth") {
        toysMap.push_back(pair<int,pair<string,string> >(order,make_pair(name,title)));
      }
      if (type=="test") {
        testMap.push_back(pair<int,pair<string,string> >(order,make_pair(name,title)));
      }
    }
  }
}

void printOptionsMap(vector<pair<int,pair<string,string> > > opts){
  for (vector<pair<int,pair<string,string> > >::iterator it=opts.begin(); it!=opts.end(); it++){
    cout << "\t" << it->first << " " << it->second.first << " " << it->second.second << endl;
  }
}

int main(int argc, char* argv[]){

  string bkgFileName;
  string sigFileName;
  string sigWSName;
  string bkgWSName;
  string outFileName;
  string datFileName;
  string outDir;
  int cat;
  int ntoys;
  int jobn;
  int seed;
  float mu_low;
  float mu_high;
  float expectSignal;
  int expectSignalMass;
  bool skipPlots=false;
  int verbosity;
  bool throwHybridToys=false;
  vector<float> switchMass;
  vector<string> switchFunc;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("sigfilename,s", po::value<string>(&sigFileName),                                          "Signal file name")
    ("bkgfilename,b", po::value<string>(&bkgFileName),                                          "Background file name")
    ("sigwsname", po::value<string>(&sigWSName)->default_value("cms_hgg_workspace"),            "Signal workspace name")
    ("bkgwsname", po::value<string>(&bkgWSName)->default_value("cms_hgg_workspace"),            "Background workspace name")
    ("outfilename,o", po::value<string>(&outFileName)->default_value("BiasStudyOut.root"),      "Output file name")
    ("datfile,d", po::value<string>(&datFileName)->default_value("config.dat"),                 "Name of datfile containing pdf info")
    ("outDir,D", po::value<string>(&outDir)->default_value("./"),                               "Name of out directory for plots")
    ("cat,c", po::value<int>(&cat),                                                             "Category")
    ("ntoys,t", po::value<int>(&ntoys)->default_value(0),                                       "Number of toys to run")
    ("jobn,j", po::value<int>(&jobn)->default_value(0),                                         "Job number")
    ("seed,r", po::value<int>(&seed)->default_value(0),                                         "Set random seed")
    ("mulow,L", po::value<float>(&mu_low)->default_value(-3.),                                  "Value of mu to start scan")
    ("muhigh,H", po::value<float>(&mu_high)->default_value(3.),                                 "Value of mu to end scan")
    ("expectSignal", po::value<float>(&expectSignal)->default_value(0.),                        "Inject signal into toy")
    ("expectSignalMass", po::value<int>(&expectSignalMass)->default_value(125),                 "Inject signal at this mass")
    ("skipPlots",                                                                               "Skip full profile and toy plots")                        
    ("verbosity,v", po::value<int>(&verbosity)->default_value(0),                               "Verbosity level")
  ;    
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("skipPlots")) skipPlots=true;
  if (expectSignalMass!=110 && expectSignalMass!=115 && expectSignalMass!=120 && expectSignalMass!=125 && expectSignalMass!=130 && expectSignalMass!=135 && expectSignalMass!=140 && expectSignalMass!=145 && expectSignalMass!=150){
    cerr << "ERROR - expectSignalMass has to be integer in range (110,150,5)" << endl;
    exit(1);
  }

  vector<pair<int,pair<string,string> > > toysMap;
  vector<pair<int,pair<string,string> > > testMap;
  readDatFile(datFileName,cat,toysMap,testMap);
  
  cout << "Toy vector.." << endl;
  printOptionsMap(toysMap);
  cout << "Test vector.." << endl;
  printOptionsMap(testMap);

  
  TStopwatch sw;
  sw.Start();
 
  if (verbosity<1) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }
  
  TFile *bkgFile = TFile::Open(bkgFileName.c_str());
  TFile *sigFile = TFile::Open(sigFileName.c_str());

  //RooWorkspace *bkgWS = (RooWorkspace*)bkgFile->Get("cms_hgg_workspace");
  RooWorkspace *bkgWS = (RooWorkspace*)bkgFile->Get(bkgWSName.c_str());
  RooWorkspace *sigWS = (RooWorkspace*)sigFile->Get(sigWSName.c_str());

  if (!bkgWS || !sigWS){
    cerr << "ERROR - one of signal or background workspace is NULL" << endl;
    cerr << " (looked for ) signal = " << sigWSName.c_str() << ", background = " << bkgWSName.c_str() <<endl;
    exit(1);
  }

  RooRealVar *mass = (RooRealVar*)bkgWS->var("CMS_hgg_mass");
  RooRealVar *mu = new RooRealVar("mu","mu",0.,mu_low,mu_high);

  TFile *outFile = new TFile(outFileName.c_str(),"RECREATE");

  //  int toyn;
  Float_t mu_;
  Float_t sigma_mu_;
  Float_t bkgSig1fwhm_;
  Float_t bkgSig2fwhm_;
  Float_t bkgTrue1fwhm_;
  Float_t bkgTrue2fwhm_;
  Float_t bkgErrSig1fwhm_;
  Float_t bkgErrSig2fwhm_;
  Float_t bkgErrNormSig1fwhm_;
  Float_t bkgErrNormSig2fwhm_;
  UInt_t iToy_;
  UInt_t iCat_;
  Int_t fitStatus_;
  char genName_[30];
  char fitName_[30];

  TTree *muTree = new TTree("muTree","muTree");
  muTree->Branch("mass", &expectSignalMass, "mass/I");
  muTree->Branch("mu", &mu_, "mu/F");
  muTree->Branch("cat", &iCat_, "cat/i");
  muTree->Branch("muTruth", &expectSignal, "muTruth/F");
  muTree->Branch("sigma_mu", &sigma_mu_, "sigma_mu/F");
  muTree->Branch("bkgTrue1fwhm", &bkgTrue1fwhm_, "bkgTrue1fwhm/F");
  muTree->Branch("bkgTrue2fwhm", &bkgTrue2fwhm_, "bkgTrue2fwhm/F");
  muTree->Branch("bkgSig1fwhm", &bkgSig1fwhm_, "bkgSig1fwhm/F");
  muTree->Branch("bkgSig2fwhm", &bkgSig2fwhm_, "bkgSig2fwhm/F");
  muTree->Branch("bkgErrSig1fwhm", &bkgErrSig1fwhm_, "bkgErrSig1fwhm/F");
  muTree->Branch("bkgErrSig2fwhm", &bkgErrSig2fwhm_, "bkgErrSig2fwhm/F");
  muTree->Branch("bkgErrNormSig1fwhm", &bkgErrNormSig1fwhm_, "bkgErrNormSig1fwhm/F");
  muTree->Branch("bkgErrNormSig2fwhm", &bkgErrNormSig2fwhm_, "bkgErrNormSig2fwhm/F");
  muTree->Branch("genFun", genName_, "genFun/C");
  muTree->Branch("fitFun", fitName_, "fitFun/C");
  muTree->Branch("iToy", &iToy_, "iToy/i");
  muTree->Branch("fitStatus",&fitStatus_, "fitStatus/I");

  
  //TH1F *muDistFab = new TH1F("muDistFab","muDistFab",int(20*(mu_high-mu_low)),mu_low,mu_high);
  //TH1F *muDistPaul = new TH1F("muDistPaul","muDistPaul",int(20*(mu_high-mu_low)),mu_low,mu_high);
  //TH1F *muDistChi2 = new TH1F("muDistChi2","muDistChi2",int(20*(mu_high-mu_low)),mu_low,mu_high);
  //TH1F *muDistAIC = new TH1F("muDistAIC","muDistAIC",int(20*(mu_high-mu_low)),mu_low,mu_high);
  
  mass->setBins(320);
  RooDataSet *data = (RooDataSet*)bkgWS->data(Form("data_mass_cat%d",cat));
  //RooDataSet *data = (RooDataSet*)bkgWS->data(Form("data_cat%d_7TeV",cat));
  RooDataHist *dataBinned = new RooDataHist(Form("roohist_data_mass_cat%d",cat),Form("roohist_data_mass_cat%d",cat),RooArgSet(*mass),*data);
  RooDataSet *sigMC = (RooDataSet*)sigWS->data(Form("sig_ggh_mass_m%d_cat%d",expectSignalMass,cat));
  RooDataSet *sigMC_vbf = (RooDataSet*)sigWS->data(Form("sig_vbf_mass_m%d_cat%d",expectSignalMass,cat));
  RooDataSet *sigMC_wh = (RooDataSet*)sigWS->data(Form("sig_wh_mass_m%d_cat%d",expectSignalMass,cat));
  RooDataSet *sigMC_zh = (RooDataSet*)sigWS->data(Form("sig_zh_mass_m%d_cat%d",expectSignalMass,cat));
  RooDataSet *sigMC_tth = (RooDataSet*)sigWS->data(Form("sig_tth_mass_m%d_cat%d",expectSignalMass,cat));
  std::cout << "Signal Model Building " << std::endl; 
  sigMC->Print(); sigMC_vbf->Print(); sigMC_wh->Print();sigMC_zh->Print();sigMC_tth->Print();
  sigMC->append(*sigMC_vbf);
  sigMC->append(*sigMC_wh);
  sigMC->append(*sigMC_zh);
  sigMC->append(*sigMC_tth);
  //RooExtendPdf *ggh_pdf = (RooExtendPdf*)sigWS->pdf(Form("sigpdfsmrel_cat%d_7TeV_ggh",cat));
  //RooExtendPdf *vbf_pdf = (RooExtendPdf*)sigWS->pdf(Form("sigpdfsmrel_cat%d_7TeV_vbf",cat));
  //RooExtendPdf *wzh_pdf = (RooExtendPdf*)sigWS->pdf(Form("sigpdfsmrel_cat%d_7TeV_wzh",cat));
  //RooExtendPdf *tth_pdf = (RooExtendPdf*)sigWS->pdf(Form("sigpdfsmrel_cat%d_7TeV_tth",cat));
  //RooAbsPdf *sigPdf = new RooAddPdf(Form("sigpdfsmrel_cat%d_7TeV",cat),Form("sigpdfsmrel_cat%d_7TeV",cat),RooArgList(*ggh_pdf,*vbf_pdf,*wzh_pdf,*tth_pdf));
  
  if (!dataBinned || !sigMC){
    cerr << "ERROR -- one of data or signal is NULL" << endl;
    exit(1);
  }
  
  // set of truth models to throw toys from
  PdfModelBuilder toysModel;
  toysModel.setObsVar(mass);
  toysModel.setSignalModifier(mu);
  // add truth pdfs from config datfile these need to be cached
  // to throw a toy from the SB fit make sure that the cache happens at makeSBPdfs
  for (vector<pair<int,pair<string,string> > >::iterator it=toysMap.begin(); it!=toysMap.end(); it++){
    if (it->first==-1) { // this is a hyrbid toy
      throwHybridToys=true;
      vector<string> temp;
      split(temp,it->second.first,boost::is_any_of(","));
      split(switchFunc,it->second.second,boost::is_any_of(","));
      for (unsigned int i=0; i<temp.size(); i++){
        switchMass.push_back(atof(temp[i].c_str()));
      }
      continue; 
    }
    if (it->first==-2) { // this is a keys pdf toy
      double rho = lexical_cast<double>(it->second.first);
      toysModel.setKeysPdfAttributes(data,rho);
      toysModel.addBkgPdf("KeysPdf",0,Form("truth_%s_cat%d",it->second.second.c_str(),cat),false);
      continue;
    }
    if (it->first==-3) { // this is read pdf from file
      toysModel.addBkgPdf(it->second.second,it->first,it->second.first,false);
      continue;
    }
    toysModel.addBkgPdf(it->second.second,it->first,Form("truth_%s_cat%d",it->second.first.c_str(),cat),false); 
  }
  toysModel.setSignalPdfFromMC(sigMC);
  //toysModel.setSignalPdf(sigPdf);
  toysModel.makeSBPdfs(true);
  map<string,RooAbsPdf*> toyBkgPdfs = toysModel.getBkgPdfs();
  map<string,RooAbsPdf*> toySBPdfs = toysModel.getSBPdfs();
  toysModel.setSeed(seed);

  // tests chosen model
  PdfModelBuilder testModel;
  testModel.setObsVar(mass);
  testModel.setSignalModifier(mu);
  // add pdfs from config datfile - should be no need to cache these
  for (vector<pair<int,pair<string,string> > >::iterator it=testMap.begin(); it!=testMap.end(); it++){
    testModel.addBkgPdf(it->second.second,it->first,Form("test_%s_cat%d",it->second.first.c_str(),cat),false); 
  }
  testModel.setSignalPdfFromMC(sigMC);
  //testModel.setSignalPdf(sigPdf);
  testModel.makeSBPdfs(false);
  map<string,RooAbsPdf*> testBkgPdfs = testModel.getBkgPdfs();
  map<string,RooAbsPdf*> testSBPdfs = testModel.getSBPdfs();


  if (!skipPlots) {
    system(Form("mkdir -p %s/plots/truthToData",outDir.c_str()));
    system(Form("mkdir -p %s/plots/envelopeNlls",outDir.c_str()));
    system(Form("mkdir -p %s/plots/toys",outDir.c_str()));
  }
  
  // throw toys - only need to fit data once as result will be cached
  cout << "------ FITTING TRUTH TO DATA ------" << endl;
  // sometimes useful to do best fit first to get reasonable starting value
  toysModel.setSignalModifierVal(0);
  toysModel.setSignalModifierConstant(true);
  toysModel.fitToData(dataBinned,false,false,false);
  // -----
  //  toysModel.setSignalModifierVal(expectSignal);
  //  toysModel.setSignalModifierConstant(true);
  toysModel.fitToData(dataBinned,false,true,true);
  //  toysModel.setSignalModifierVal(expectSignal);
  //  toysModel.setSignalModifierConstant(false);
  toysModel.saveWorkspace(outFile);
  if (!skipPlots) 
      toysModel.plotPdfsToData(dataBinned,80,Form("%s/plots/truthToData/datafit_mu%3.1f",outDir.c_str(),expectSignal),false);

  iCat_=cat;

  for (int toy=jobn*ntoys; toy<(jobn+1)*ntoys; toy++){
    cout << "---------------------------" << endl;
    cout << "--- RUNNING TOY " << toy << " / " << (jobn+1)*ntoys << " ----" << endl;
    cout << "---------------------------" << endl;

    iToy_=toy;
    // throw toy
    map<string,RooAbsData*> toys; 
    //     if (throwHybridToys) {
    //       toysModel.throwHybridToy(Form("truth_job%d_toy%d",jobn,toy),dataBinned->sumEntries(),switchMass,switchFunc,false,true,true,true);
    //       toys = toysModel.getHybridToyData();
    //       if (!skipPlots) toysModel.plotToysWithPdfs(Form("%s/plots/toys/job%d_toy%d",outDir.c_str(),jobn,toy),80,false);
    //       if (!skipPlots) toysModel.plotHybridToy(Form("%s/plots/toys/job%d_toy%d",outDir.c_str(),jobn,toy),80,switchMass,switchFunc,false);
    //     }
    //     else {
    toysModel.throwToy(Form("truth_job%d_toy%d",jobn,toy),dataBinned->sumEntries(),false,true,true,true,true,expectSignal);
    toys = toysModel.getToyData();
    if (!skipPlots) 
      toysModel.plotToysWithPdfs(Form("%s/plots/toys/job%d_toy%d",outDir.c_str(),jobn,toy),80,false);
    //     }
    for (map<string,RooAbsData*>::iterator it=toys.begin(); it!=toys.end(); it++){
       // ----- USEFUL DEBUG -----------
       //  --- this can be a useful check that the truth model values are being cached properly ---
       //toysModel.fitToData(it->second,true,false,true);
       //toysModel.plotPdfsToData(it->second,80,Form("%s/plots/toys/job%d_toy%d",outDir.c_str(),jobn,toy),true,"NONE");

      testModel.setSignalModifierVal(expectSignal);
      testModel.setSignalModifierConstant(false);
      testModel.fitToData(it->second,false,true,true,true);
      if (!skipPlots)
	testModel.plotPdfsToData(it->second,80,Form("%s/plots/toys/fit_%s",outDir.c_str(),it->first.c_str()),false,"",true);

      //      testModel.wsCache->Print("v");

      // -----  SAVE RESULTS IN THE TREE ----- 
      for (map<string,RooAbsPdf*>::iterator itf=testSBPdfs.begin(); itf!=testSBPdfs.end(); itf++){
	if (!testModel.wsCache->loadSnapshot(itf->first.c_str()))
	  {
	    cout << "CANNOT GET WS SNAPSHOT " << itf->first << endl;
	    continue;
	  }

	mu_= testModel.wsCache->var(mu->GetName())->getVal();
	sigma_mu_ = testModel.wsCache->var(mu->GetName())->getError();
	
	typedef vector< string > split_vector_type;
    	split_vector_type splitVec_gen; // #2: Search for tokens
    	split_vector_type splitVec_fit; // #2: Search for tokens
	split( splitVec_gen, it->first, boost::is_any_of("_") );
	split( splitVec_fit, itf->first, boost::is_any_of("_") );
	sprintf(genName_,"%s",splitVec_gen[2].c_str());	
	sprintf(fitName_,"%s",splitVec_fit[2].c_str());
	muTree->Fill();
	testModel.wsCache->Clear();
      }
      delete it->second;
    }

  }

  outFile->cd();
  muTree->Write();
  cout << "Done." << endl;
  cout << "Whole process took..." << endl;
  cout << "\t "; sw.Print();
 
  outFile->Close();

  return 0;
}
