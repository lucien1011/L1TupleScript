#include "L1Ntuple.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"

class trigTiming : public L1Ntuple {

 public :

  trigTiming() {}
  ~trigTiming() {}

  void Loop(const std::string & fname);

  void bookHistograms();
  void loadTrigMap();
  void printJetinfo(const float & );
  void printRunEventLumi(const int &);
  bool checkTriggerBit(const int & ,const int & );

private :

  TH1F *hNevts;
  TObjArray* Hlist;

  std::map <string,int> L1TriggerBitMap;
  std::map<string,int>::iterator trigbit_iter;

  std::map<string, TH1F*> m_HistNames;
  std::map<string, TH1F*>::iterator hid;

};

void trigTiming::Loop(const std::string & oname){

  // string MyTrigger="HLT_Any";
  string MyTrigger="HLT_L1SingleMu3p5_v1";
  // string MyTrigger="HLT_L1SingleEG5_v1";
  // string MyTrigger="HLT_ZeroBias_v1";

  loadTrigMap();
  bookHistograms();

  Int_t nevents = GetEntries();
  std::cout << "Running on " << nevents << " events." << std::endl;

  // nevents=20000;
  for (Long64_t event=0; event<nevents; event++)
    {     
      Long64_t ientry = LoadTree(event); if (ientry < 0) break;
      GetEntry(event);

      if (event%5000 == 0) {
	printRunEventLumi(event);
      }

      //if (event_->run == 247644 && (event_->lumi>100 && event_->lumi<120)) continue;


      bool accept=false;
      if (MyTrigger=="HLT_Any") {
	accept=true;
      }else{
	int ntrigs=event_->hlt.size();
	for (int itrig=0; itrig<ntrigs; itrig++){
	  if (event_->hlt[itrig] == MyTrigger) accept=true;
	}
      }
      if ( not accept) continue;

      hNevts->Fill(event_->lumi);

      Int_t nBx=gt_->tw1.size();

      for (trigbit_iter = L1TriggerBitMap.begin(); trigbit_iter != L1TriggerBitMap.end(); trigbit_iter++){
        string hname=trigbit_iter->first;
        int ibit=trigbit_iter->second;
        //std::cout << hname << "\t" << ibit << endl;
        for(Int_t ibx=0; ibx < nBx; ibx++){

          bool Fired = checkTriggerBit(ibit,ibx);
          if (Fired) {
            hid=m_HistNames.find(hname);
            if (hid==m_HistNames.end())
              std::cout << "%fillHist -- Could not find histogram with name: " << hname << std::endl;
            else
              hid->second->Fill(ibx); 
	  }
	}
      }// end of loop over trigger bits
 
    }// end of Loop over Events

  hNevts->Scale(1./23.3);

  string outhist=oname + "_" + MyTrigger + ".root";
  // outhist="/dev/null";
  std::cout << "Output written to: " << outhist << std::endl;
  TFile *of = new TFile(outhist.c_str(),"RECREATE");

  Hlist->Write();
  of->Close();

  return;
}

void trigTiming::bookHistograms(){

  Hlist = new TObjArray();

  hNevts = new TH1F("hNevts","Number of Events per Lumi Section",1050,0,1050.);   Hlist->Add(hNevts);

  for (trigbit_iter = L1TriggerBitMap.begin(); trigbit_iter != L1TriggerBitMap.end(); trigbit_iter++)
    {
      std::cout <<" booking histogram from trigger: "<< trigbit_iter->first << std::endl;
      string hname=trigbit_iter->first;
      string htitle="Rate vs  BX -- " + trigbit_iter->first;
      m_HistNames[hname]= new TH1F(hname.c_str(),htitle.c_str(),5,0.,5);
      m_HistNames[hname]->GetXaxis()->SetTitle("BX");
      Hlist->Add(m_HistNames[hname]);
    }
}

bool trigTiming::checkTriggerBit(const int & ibit,const int & ibx){

  bool Fired(false);
  if (ibit<64){
    Fired = (gt_->tw1[ibx]>>ibit)&1;
  }else
    Fired = (gt_->tw2[ibx]>>(ibit-64))&1;

  return Fired;
}

void trigTiming::loadTrigMap(){ // these bits can change! check these are correct for the run you are using

  L1TriggerBitMap["L1_ZeroBias"]  =0;

  L1TriggerBitMap["L1_MinimumBiasHF1_OR"] =4;
  L1TriggerBitMap["L1_MinimumBiasHF2_OR"] =5;
  L1TriggerBitMap["L1_MinimumBiasHF1_AND"]=6;
  L1TriggerBitMap["L1_MinimumBiasHF2_AND"]=7;

  L1TriggerBitMap["L1_SingleJet8_BptxAND"]=10;
  L1TriggerBitMap["L1_SingleJet12_BptxAND"]=11;
  L1TriggerBitMap["L1_SingleJet16"]=12;
  L1TriggerBitMap["L1_SingleJet20"]=13;
  L1TriggerBitMap["L1_SingleJet36"]=14;
  L1TriggerBitMap["L1_SingleJet68"]=15;
  L1TriggerBitMap["L1_SingleJet200"]=16;

  L1TriggerBitMap["L1_DoubleJet28"]=21;

  L1TriggerBitMap["L1_SingleEG2_BptxAND"] =30;
  L1TriggerBitMap["L1_SingleEG5"] =31;
  L1TriggerBitMap["L1_SingleEG20"]=32;

  L1TriggerBitMap["L1_ETT15"]=40;
  L1TriggerBitMap["L1_ETT40"]=41;
  L1TriggerBitMap["L1_ETT60"]=42;
  L1TriggerBitMap["L1_ETT90"]=43;
  L1TriggerBitMap["L1_ETT130"]=44;

  L1TriggerBitMap["L1_SingleMu3p5"]=51;
  L1TriggerBitMap["L1_SingleMuBeamHalo"]=53;
  L1TriggerBitMap["L1_SingleMuOpen_NotBptxOR"]=54;

  L1TriggerBitMap["L1_SingleMuOpen"]=98;
  L1TriggerBitMap["L1_DoubleMuOpen"]=99;
}

void trigTiming::printRunEventLumi(const int & count){
  std::cout << "Processed " << count << " events."
	    << "\tCurrent Run: " << event_->run 
	    << "\tEvent: "       << event_->event 
	    << "\tLumi: "        << event_->lumi 
	    << std::endl;
}

void trigTiming::printJetinfo(const float & minPt){


  for(Int_t ue=0; ue < gt_->Njet; ue++) {
    Int_t bx = gt_ -> Bxjet[ue];
    Float_t pt = gt_ -> Rankjet[ue];
    // std::cout << "Jet --  XXX " << pt << std::endl;
    if (pt>=minPt) {
      std::cout << "Jet --  Candidate pT: " << pt 
		<< "\tEta: " << gt_ -> Etajet[ue]
		<< "\tPhi: " << gt_ -> Phijet[ue]
		<< "\tRel BX: " << gt_ -> Bxjet[ue]
		<< std::endl;
    }
  }
  return;
}

void RunL1(std::string runnumber){


  //std::string L1NtupleFileName = "/eos/uscms/store/user/lpctrig/apana/Collisions2015/run"+runnumber+"/L1Tree_combined.root";
  std::string L1NtupleFileName = "L1Tree.root";
  std::string OutputFileName="triggerEff_"+runnumber+"_ExpressDS";

  trigTiming a;
  a.Open(L1NtupleFileName);
  a.Loop(OutputFileName);

 return;
}
