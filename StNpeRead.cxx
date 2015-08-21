#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <set>
#include "StNpeRead.h"
#include "StCuts.h"
//#include "StDmesonEvent.h"
#include "StTRIGGERS.h"
//#include "StNpeMaker/StDmesonEvent.h"
#include "prescales.h"
#include "StDmesonMaker/StDmesonTrack.h"

#include "StDmesonMaker/StElectronPair.h"
#include "StDmesonMaker/StKaonKaon.h"

#include "StLorentzVectorF.hh"

#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TRandom3.h"

//Float_t KAONMASS2 = 0.49367 * 0.49367;
//Float_t PHIMASS = 1.01945; // GeV/c/c

using namespace std;

ClassImp(StNpeRead)
//-----------------------------------------------------------------------------
StNpeRead::StNpeRead(const char* outName)
{
  numPtBins = 14; maxBufferSize = 10;
  writeXiaozhiHists = kFALSE; writeDataQA = kFALSE; writeRunQA = kFALSE; // Set flags for QA writing (to minimize file size when QA unnecessary)
  pi = 3.1415826;
  // Initialize the hadron weighting function (if can't find, exit)
  if(!readEff())
    exit(1);

  mOutputFile = new TFile(outName, "RECREATE");
  TH1F::SetDefaultSumw2();
  mPrescales=prescales::Instance();
  file_runNumber.open("/global/homes/z/zamiller/NPEhPhiAnalysis2015/StRoot/StNpeRead/runIDPicoSet.txt");
  //  outfile.open("/global/u1/x/xiao00/work/runIDPico.txt",ios::out|ios::app);  
  SetRunID_VPD();
}
void  StNpeRead::SetRunID_VPD( )
{

  //runID_VPDMB
  Int_t run_index;

  if(!file_runNumber.is_open()){
    cout<< "  Error open file ";
    exit(1);
  }
 
  if(file_runNumber.is_open())
    {
      // std::string line;                                                                                                        
      // start from beginning                                                  
                                               
      //  file_runNumber.clear() ;                                                                                                
      // file_runNumber.seekg(0, ios::beg) ;  

      runID_List.clear();
    
      // std::istringstream iss(line); 
      int runID;
      int runindex=0;
      // while (std::getline(file_runNumber, line))
      while(file_runNumber>>runID)
	{
	  // iss>>runID;  
	  // cout << runID << endl;
	  runID_List[runID]=runindex;
	  runindex++;
	}
      cout<<"The file has been opened."<<endl;
    }
  
  //  cout << runID_List.size() << endl;;

  for(map<int,int>::iterator iter=runID_List.begin();iter!=runID_List.end();iter++)
    {
      // std::cout<< iter->second<<"   \t"<<"  "<<iter->first<<std::endl;
    }
  
  
}
/*

  int runIndex=-1, runID=-999; 

  std::string line;
  // start from beginning
  file_runNumber.clear() ;
  file_runNumber.seekg(0, ios::beg) ;

  while (std::getline(file_runNumber, line)) {
    // while(!file_runNumber.eof()) {
    //file_runNumber>>runIndex>>runID;
    std::istringstream iss(line);
    iss >> runIndex >> runID;
    // std::cout<< runIndex<<"  "<<runID<<"  "<<run_ID<<std::endl;
    
    if(runID==run_ID) {
      run_index=runIndex;
      // std::cout<< runIndex<<"  "<<runID<<"  "<< run_ID<<std::endl;
      
    }
    
  }
  return run_index;
   
  }*/
//runID_List();

//-----------------------------------------------------------------------------
StNpeRead::~StNpeRead()
{
   /*  */
}

//-----------------------------------------------------------------------------
void StNpeRead::bookObjects()
{
  mNpeEvent = new StDmesonEvent;
 
  /// Z. Histos
  for(Int_t trg=0; trg<4; trg++)
    {
      mh1PtAllTracks[trg]      = new TH1F(Form("mh1PtAllTracks_%i",trg),"",400,0,20);
      mh1PtETracks[trg]        = new TH1F(Form("mh1PtETracks_%i",trg),"",400,0,20);
      mh2nSigmaEPt[trg]        = new TH2F(Form("mh2nSigmaEPt_%i",trg),"",600,-30,30,400,0,20);
      mh2nSigmaEPt_eID[trg]    = new TH2F(Form("mh2nSigmaEPt_eID_%i",trg),"",600,-30,30,400,0,20);
      mh2PoePt[trg]            = new TH2F(Form("mh2PoePt_%i",trg),"",600,-30,30,400,0,20);
      mh2nPhiPt[trg]           = new TH2F(Form("mh2nPhiPt_%i",trg),"",20,0,20,400,0,20);
      mh2nEtaPt[trg]           = new TH2F(Form("mh2nEtaPt_%i",trg),"",20,0,20,400,0,20);
      mh2PhiDistPt[trg]        = new TH2F(Form("mh2PhiDistPt_%i",trg),"",400,-10,10,400,0,20);
      mh2ZDistPt[trg]          = new TH2F(Form("mh2ZDistPt_%i",trg),"",400,-20,20,400,0,20);
      mh2PhiQPt[trg]           = new TH2F(Form("mh2PhiQPt_%i",trg),"",400,-10,10,800,-20,20);
      mh2TofPtAll[trg]         = new TH2F(Form("mh2TofPtAll_%i",trg),"",400,-20,20,400,0,20);
      mh2TofPtE[trg]           = new TH2F(Form("mh2TofPtE_%i",trg),"",400,-20,20,400,0,20);
      mh2InvMassPtAll[trg]     = new TH2F(Form("mh2InvMassPtAll_%i",trg),"",1000,0,10,1000,0,10);
      mh2InvMassPtE[trg]       = new TH2F(Form("mh2InvMassPtE_%i",trg),"",1000,0,10,1000,0,10);
      mh2InvMassPtLS[trg]      = new TH2F(Form("mh2InvMassPtLS_%i",trg),"",1000,0,10,1000,0,10);
      mh2InvMassPtUS[trg]      = new TH2F(Form("mh2InvMassPtUS_%i",trg),"",1000,0,10,1000,0,10);

      mh2nSigmaPionPt[trg]       = new TH2F(Form("mh2nSigmaPionPt_%i",trg),"",1000,-10,10,1000,0,10);

      mh3DelPhiIncl[trg]       = new TH3F(Form("mh3DelPhiIncl_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh3DelPhiPhotLS[trg]     = new TH3F(Form("mh3DelPhiPhotLS_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh3DelPhiPhotUS[trg]     = new TH3F(Form("mh3DelPhiPhotUS_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh3DelPhiPhotUSNP[trg]   = new TH3F(Form("mh3DelPhiPhotUSNP_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh3DelPhiPhotLSNP[trg]   = new TH3F(Form("mh3DelPhiPhotLSNP_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh3DelPhiPhotInclNP[trg] = new TH3F(Form("mh3DelPhiPhotInclNP_%i",trg),"",200,-10,10,200,0,20,40,0,20);

      mh3DelPhiInclWt[trg]     = new TH3F(Form("mh3DelPhiInclWt_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh3DelPhiPhotLSWt[trg]   = new TH3F(Form("mh3DelPhiPhotLSWt_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh3DelPhiPhotUSWt[trg]   = new TH3F(Form("mh3DelPhiPhotUSWt_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh3DelPhiInclWt[trg]->Sumw2();mh3DelPhiPhotLSWt[trg]->Sumw2();mh3DelPhiPhotUSWt[trg]->Sumw2();

      mh3DelPhiHadHad[trg]     = new TH3F(Form("mh3DelPhiHadHad_%i",trg),"",200,-10,10,200,0,20,40,0,20);
      mh1PtHadTracks[trg]      = new TH1F(Form("mh1PtHadTracks_%i",trg),"",400,0,20);

      mh3nTracksZdcx[trg]      = new TH3F(Form("mh3nTracksZdcx_%i",trg),"",200,0,20000,30,0,30,10,0,2);

    }

  /// Mixed Events
  mh3MixedDelPhi          = new TH3F("mh3MixedDelPhi","",400,-10,10,200,0,20,40,0,20);
  mh3MixedDelEta          = new TH3F("mh3MixedDelEta","",400,-10,10,200,0,20,40,0,20);
  mh3MixedEtaPhi          = new TH3F("mh3MixedEtaPhi","",400,-10,10,200,-5,5,40,0,20);
  

  /*
  Int_t bin1D[5]={1400,1400,1400,1400,1400}; Double_t xMin1D[5]={-7,-7,-7,-7,-7};  Double_t xMax1D[5]={7,7,7,7,7};
  Int_t binx2D[5]={700,1400,1400,1400,1400}; Double_t xMin2D[5]={-40,-7,-7,-7,-7}; Double_t xMax2D[5]={40,7,7,7,7};
  Int_t biny2D[5]={80,1400,1400,1400,1400};  Double_t yMin2D[5]={0,-7,-7,-7,-7};   Double_t yMax2D[5]={4,7,7,7,7};
  for(Int_t ii=0; ii<5; ii++)
    {
      testHist1D[ii] = new TH1F(Form("testHist1D_%i",ii),"",bin1D[ii],xMin1D[ii],xMax1D[ii]);
      testHist2D[ii] = new TH2F(Form("testHist2D_%i",ii),"",binx2D[ii],xMin2D[ii],xMax2D[ii],biny2D[ii],yMin2D[ii],yMax2D[ii]);
      }*/

  /// X. Histos
  /* if(1)// doesn't work right now -> if(writeXiaozhiHists)
   {
      HT0_HT2=new TH2F("HT0_HT2","",2,0,2,2,0,2);
      
      // ADC_nocut=new TH3F("ADC_nocut","",100,0,100,20,0,2000,20,0,2000);
      // ADC_cut=new TH3F("ADC_cut","",100,0,100,20,0,200,200,0,2000);
      mh1Vz_BBCMBTOF0=new TH1F(Form("mh1Vz_BBCMBTOF0"),"",200,-120,120);
      mh1Vz_BBCMB=new TH1F(Form("mh1Vz_BBCMB"),"",200,-120,120);
      
      mh1VPDVz_BBCMBTOF0=new TH1F(Form("mh1VPDVz_BBCMBTOF0"),"",200,-120,120);
      mh1VPDVz_BBCMB=new TH1F(Form("mh1VPDVz_BBCMB"),"",200,-120,120);
      
      mh1Vz_BBCMBTOF0PS=new TH1F(Form("mh1Vz_BBCMBTOF0PS"),"",200,-120,120);
      mh1Vz_BBCMBPS=new TH1F(Form("mh1Vz_BBCMBPS"),"",200,-120,120);
      mh1Vz_VPD= new TH1F("mh1Vz_VPD","",240,-60,60);
      mh1Vz_VPDPS= new TH1F("mh1Vz_VPDPS","",240,-60,60);
      mh1Vz_VPDVz= new TH1F("mh1Vz_VPDVz","",1000,-50,50);
   
      //   for(Int_t iTrg=0;iTrg<4;iTrg++ )
      //	{
	  
	  // Run QA
	  /*mh2TPCVzRunID[iTrg]= new TH2F(Form("mh2TPCVzRunIDTrg%i",iTrg),"",840,-0.5,839.5,120,-60,60);
	    mh2VPDVzRunID[iTrg]= new TH2F(Form("mh2VPDVzRunIDTrg%i",iTrg),"",840,-0.5,839.5,120,-60,60);
	    mh2TPCVz_VPDVz[iTrg]= new TH2F(Form("mh2TPCVz_VPDVzTrg%i",iTrg),"",240,-60,60,120,-60,60);
	    
	    mh2ZDCRunID[iTrg]=new TH2F(Form("mh2ZDCRunIDTrg%i",iTrg),"",840,-0.5,839.5,100,0,50);
	    mh2BBCRunID[iTrg]=new TH2F(Form("mh2BBCRunIDTrg%i",iTrg),"",840,-0.5,839.5,100,0,100);
	    mh2MultiPosRunID[iTrg]=new TH2F(Form("mh2MultiPosRunIDTrg%i",iTrg),"",840,-0.5,839.5,40,0,40);
	    mh2MultiNegRunID[iTrg]=new TH2F(Form("mh2MultiNegRunIDTrg%i",iTrg),"",840,-0.5,839.5,40,0,40);
	    mh2MultiRunID[iTrg]=new TH2F(Form("mh2MultiRunIDTrg%i",iTrg),"",840,-0.5,839.5,60,0,60);
	    
	    mh2gDCARunID[iTrg]=new TH2F(Form("mh2gDCARunIDTrg%i",iTrg),"",840,-0.5,839.5,100,-1.e-6,3-1.e-6);
	    mh2PtRunID[iTrg]=new TH2F(Form("mh2PtRunIDTrg%i",iTrg),"",840,-0.5,839.5,80,0,20);
	    mh2EtaRunID[iTrg]=new TH2F(Form("mh2EtaRunIDTrg%i",iTrg),"",840,-0.5,839.5,40,-1.5,1.5);
	    mh2PhiRunID[iTrg]=new TH2F(Form("mh2PhiRunIDTrg%i",iTrg),"",840,-0.5,839.5,100,-4,4);
	    mh2BetaRunID[iTrg]=new TH2F(Form("mh2BetaRunIDTrg%i",iTrg),"",840,-0.5,839.5,100,-0.5,5);
	    mh2Beta_RunID[iTrg]=new TH2F(Form("mh2Beta_RunIDTrg%i",iTrg),"",840,-0.5,839.5,100,-0.5,5);
	    mh2NhitFitRunID[iTrg]=new TH2F(Form("mh2NhitFitRunIDTrg%i",iTrg),"",840,-0.5,839.5,50,0,50);
	    mh2NhitDedxRunID[iTrg]=new TH2F(Form("mh2NhitDedxRunIDTrg%i",iTrg),"",840,-0.5,839.5,50,0,50);
	    mh2DedxRunID[iTrg]=new TH2F(Form("mh2DedxRunIDTrg%i",iTrg),"",840,-0.5,839.5,100,-1,30-1.e-6);
	    mh2nsigmaERunID[iTrg]=new TH2F(Form("mh2nsigmaERunIDTrg%i",iTrg),"",840,-0.5,839.5,100,-20,20);
	    
	    mh1PassTofMatchRunID[iTrg]=new TH1F(Form("mh1PassTofMatchRunIDTrg%i",iTrg),"",840,-0.5,839.5);
	    mh1WPassTofMatchRunID[iTrg]=new TH1F(Form("mh1WPassTofMatchRunIDTrg%i",iTrg),"",840,-0.5,839.5);
	    
	    
	    //.........................
	    
	    mh2PhiVsEta_Range1[iTrg]=new TH2F(Form("mh2PhiVsEta_Range1Trg%i",iTrg),"",100,-4,4,48,-1.2,1.2);
	    mh2PhiVsEta_Range2[iTrg]=new TH2F(Form("mh2PhiVsEta_Range2Trg%i",iTrg),"",100,-4,4,48,-1.2,1.2);
	    mh2PhiVsEta_Range3[iTrg]=new TH2F(Form("mh2PhiVsEta_Range3Trg%i",iTrg),"",100,-4,4,48,-1.2,1.2);
	  */
	  
	  //.........................
	  
	  //Tprofile
	  /* 
	     mh2TPCVzRunID_Profile[iTrg]= new TProfile (Form("mh2TPCVzRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2VPDVzRunID_Profile[iTrg]= new TProfile (Form("mh2VPDVzRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     // mh2TPCVz_VPDVz[iTrg]= new TProfile (Form("mh2TPCVz_VPDVzTrg%i",iTrg),"",240,-60,60,120,-60,60);
	     
	     mh2ZDCRunID_Profile[iTrg]=new TProfile (Form("mh2ZDCRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2BBCRunID_Profile[iTrg]=new TProfile (Form("mh2BBCRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2MultiPosRunID_Profile[iTrg]=new TProfile (Form("mh2MultiPosRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2MultiNegRunID_Profile[iTrg]=new TProfile (Form("mh2MultiNegRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2MultiRunID_Profile[iTrg]=new TProfile(Form("mh2MultiRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     
	     mh2gDCARunID_Profile[iTrg]=new TProfile(Form("mh2gDCARunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2PtRunID_Profile[iTrg]=new TProfile(Form("mh2PtRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2EtaRunID_Profile[iTrg]=new TProfile(Form("mh2EtaRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2PhiRunID_Profile[iTrg]=new TProfile(Form("mh2PhiRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2BetaRunID_Profile[iTrg]=new TProfile(Form("mh2BetaRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2Beta_RunID_Profile[iTrg]=new TProfile(Form("mh2Beta_RunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2NhitFitRunID_Profile[iTrg]=new TProfile(Form("mh2NhitFitRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2NhitDedxRunID_Profile[iTrg]=new TProfile(Form("mh2NhitDedxRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2DedxRunID_Profile[iTrg]=new TProfile(Form("mh2DedxRunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     mh2nsigmaERunID_Profile[iTrg]=new TProfile(Form("mh2nsigmaERunID_ProfileTrg%i",iTrg),"",840,-0.5,839.5);
	     
	     
	     mh2VXY[iTrg]= new TH2F(Form("mh2VXYTrg%i",iTrg),"",100,-0.5,0.5,100,-0.5,0.5);
	  
	  
	  
	  mh2DSMADC_Inclusive[iTrg]=new TH2F(Form("mh2DSMADC_InclusiveTrg%i",iTrg),"",100,0,100,200,0,20);
	  mh1TowerID_Inclusive[iTrg]=new TH1F(Form("mh1TowerID_InclusiveTrg%i",iTrg),"",4800,0.5,4800.5);
	  mh1TowerID_cut_Inclusive[iTrg]=new TH1F(Form("mh1TowerID_cut_InclusiveTrg%i",iTrg),"",4800,0.5,4800.5);
	  mh1TowerID_all_Inclusive[iTrg]=new TH1F(Form("mh1TowerID_all_InclusiveTrg%i",iTrg),"",4800,0.5,4800.5);
	  
	  mh1TowerID_Inclusive_noHOt[iTrg]=new TH1F(Form("mh1TowerID_Inclusive_noHOtTrg%i",iTrg),"",4800,0.5,4800.5);
	  mh1TowerID_cut_Inclusive_noHOt[iTrg]=new TH1F(Form("mh1TowerID_cut_Inclusive_noHOtTrg%i",iTrg),"",4800,0.5,4800.5);
	  mh1TowerID_all_Inclusive_noHOt[iTrg]=new TH1F(Form("mh1TowerID_all_Inclusive_noHOtTrg%i",iTrg),"",4800,0.5,4800.5);
	  
	  mh2DSMADC_PhotoUnlike[iTrg]=new TH2F(Form("mh2DSMADC_PhotoUnlikeTrg%i",iTrg),"",100,0,100,200,0,20);
	  mh2DSMADC_Photolike[iTrg]=new TH2F(Form("mh2DSMADC_PhotolikeTrg%i",iTrg),"",100,0,100,200,0,20);
	  
	  mh2ADC0_Inclusive[iTrg]=new TH2F(Form("mh2ADC0_InclusiveTrg%i",iTrg),"",1500,0,1500,200,0,20);
	  mh2ADC0_PhotoUnlike[iTrg]=new TH2F(Form("mh2ADC0_PhotoUnlikeTrg%i",iTrg),"",1500,0,1500,200,0,20);
	  mh2ADC0_Photolike[iTrg]=new TH2F(Form("mh2ADC0_PhotolikeTrg%i",iTrg),"",1500,0,1500,200,0,20);
	  
	  mh2ADC0_DSMadcInclusive[iTrg]=new TH2F(Form("mh2ADC0_DSMadcInclusiveTrg%i",iTrg),"",1500,0,1500,100,0,100);
	  mh2ADC0_DSMadcPhotoUnlike[iTrg]=new TH2F(Form("mh2ADC0_DSMadcPhotoUnlikeTrg%i",iTrg),"",1500,0,1500,100,0,100);
	  mh2ADC0_DSMadcPhotolike[iTrg]=new TH2F(Form("mh2ADC0_DSMadcPhotolikeTrg%i",iTrg),"",1500,0,1500,100,0,100);
	  
	  mh2DsmADCUnlike[iTrg]=new TH2F(Form("mh2DsmADCUnlikeTrg%i",iTrg),"",100,0,100,200,0,20);
	  mh2DsmADClike[iTrg]=new TH2F(Form("mh2DsmADClikeTrg%i",iTrg),"",100,0,100,200,0,20);
	  
	  mh2Btowadc0[iTrg]=new TH2F(Form("mh2Btowadc0Trg%i",iTrg),"",2000,0,2000,200,0,20);
	  mh2InvMassPoeUnlike[iTrg] = new TH2F(Form("mh2InvMassPoeUnlikeTrg%i",iTrg),"",40,0,0.4,200,0,20.);
	  mh2PoeUnlike[iTrg] = new TH2F(Form("mh2PoeUnlikeTrg%i",iTrg),"",100,0,2.5,200,0,20.);
	  mh2InvMassPoelike[iTrg] = new TH2F(Form("mh2InvMassPoelikeTrg%i",iTrg),"",40,0,0.4,200,0,20.);
	  mh2Poelike[iTrg] = new TH2F(Form("mh2PoelikeTrg%i",iTrg),"",100,0,2.5,200 ,0,20.);
	  mh2InvMassDzUnlike[iTrg] = new TH2F(Form("mh2InvMassDzUnlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
	  mh2DzUnlike[iTrg] = new TH2F(Form("mh2DzUnlikeTrg%i",iTrg),"",100,-5,5,200 ,0,20.);
	  mh2InvMassDzlike[iTrg] = new TH2F(Form("mh2InvMassDzlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
	  mh2Dzlike[iTrg] = new TH2F(Form("mh2DzlikeTrg%i",iTrg),"",100,-5,5,200 ,0,20.);
	  
	  
	  mh2InvMassDpUnlike[iTrg] = new TH2F(Form("mh2InvMassDpUnlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);  
	  mh2DpUnlike[iTrg] = new TH2F(Form("mh2DpUnlikeTrg%i",iTrg),"",50,-0.1,0.1,200 ,0,20.);
	  mh2InvMassDplike[iTrg] = new TH2F(Form("mh2InvMassDplikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
	  mh2Dplike[iTrg] = new TH2F(Form("mh2DplikeTrg%i",iTrg),"",50,-0.1,0.1,200 ,0,20.);
	  
	  mh2InvMassEMCUnlike[iTrg] = new TH2F(Form("mh2InvMassEMCUnlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);  
	  mh2InvMassEMClike[iTrg] = new TH2F(Form("mh2InvMassEMClikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
	  mh2NPhiUnlike[iTrg] = new TH2F(Form("mh2NPhiUnlikeTrg%i",iTrg),"",10,0,10,200 ,0,20.);
	  mh2NPhilike[iTrg] = new TH2F(Form("mh2NPhilikeTrg%i",iTrg),"",10,0,10,200 ,0,20.);
	  
	  mh2NEtaUnlike[iTrg] = new TH2F(Form("mh2NEtaUnlikeTrg%i",iTrg),"",10,0,10,200 ,0,20.);
	  
	  mh2NEtalike[iTrg] = new TH2F(Form("mh2NEtalikeTrg%i",iTrg),"",10,0,10,200 ,0,20.);
	  
	  mh2InvMassNEMCUnlike[iTrg] = new TH2F(Form("mh2InvMassNEMCUnlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);  
	  mh2InvMassNEMClike[iTrg] = new TH2F(Form("mh2InvMassNEMClikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
	  mh2InvMassNEEMCUnlike[iTrg] = new TH2F(Form("mh2InvMassNEEMCUnlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);  
	  mh2InvMassNEEMClike[iTrg] = new TH2F(Form("mh2InvMassNEEMClikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
	  mh2InvMassNPEMCUnlike[iTrg] = new TH2F(Form("mh2InvMassNPEMCUnlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);  
	  mh2InvMassNPEMClike[iTrg] = new TH2F(Form("mh2InvMassNPEMClikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
	  
	  mh3nSigPart_EMCUnlike[iTrg] = new TH3F(Form("mh3nSigPart_EMCUnlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  mh3nSigPart_EMClike[iTrg] = new TH3F(Form("mh3nSigPart_EMClikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  
	  mh3nSigPart_TREMCUnlike[iTrg] = new TH3F(Form("mh3nSigPart_TREMCUnlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  mh3nSigPart_TREMClike[iTrg] = new TH3F(Form("mh3nSigPart_TREMClikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  
	  mh3nSigPart_ADCTREMCUnlike[iTrg] = new TH3F(Form("mh3nSigPart_ADCTREMCUnlikeTrg%i",iTrg),"",140,-7,7,40 ,0,40.,40,0,0.4);
	  mh3nSigPart_ADCTREMClike[iTrg] = new TH3F(Form("mh3nSigPart_ADCTREMClikeTrg%i",iTrg),"",140,-7,7,40 ,0,40.,40,0,0.4);
	  
	  mh3nSigEUnlike[iTrg] = new TH3F(Form("mh3nSigEUnlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  mh3nSigPartUnlike[iTrg] = new TH3F(Form("mh3nSigPartUnlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  mh3nSigPartlike[iTrg] = new TH3F(Form("mh3nSigPartlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  mh3nSigElike[iTrg] = new TH3F(Form("mh3nSigElikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  //all cuts applied
	  mh2InvMassUnlike[iTrg]=new TH2F(Form("mh2InvMassUnlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
	  mh2InvMasslike[iTrg]=new TH2F(Form("mh2InvMasslikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20);
	  
	  mh3EMC_PartUnlike[iTrg]=new TH3F(Form("mh3EMC_PartUnlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
	  mh3EMC_Partlike[iTrg]=new TH3F(Form("mh3EMC_PartlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20,40,0,0.4);
	  mh3EMC_ADCPartUnlike[iTrg]=new TH3F(Form("mh3EMC_ADCPartUnlikeTrg%i",iTrg),"",140,-7,7,40 ,0,40.,40,0,0.4);
	  mh3EMC_ADCPartlike[iTrg]=new TH3F(Form("mh3EMC_ADCPartlikeTrg%i",iTrg),"",140,-7,7,40 ,0,40,40,0,0.4);     
	  mh3EtaPhiUnlike[iTrg]=new TH3F(Form("mh3EtaPhiUnlikeTrg%i",iTrg),"",20 ,-1,1,80,-4,4,200,0,20);
	  mh3EtaPhilike[iTrg]=new TH3F(Form("mh3EtaPhilikeTrg%i",iTrg),"",20 ,-1,1,80,-4,4,200,0,20);
	  
	  //for inclusive elctron
	  mh2nSigmaElec[iTrg] = new TH2F(Form("mh2nSigmaElecTrg%i",iTrg),"",200,-10,10,200 ,0,20.);
	  mh1electronPt[iTrg]=new TH1F(Form("mh1electronPtTrg%i",iTrg),"",200,0,20.);
	  
	  mh1Vz[iTrg] = new TH1F(Form("mh1VzTrg%i",iTrg),"",200,-120,120);
	  // mh1Vz_VPDVz= new TH1F("mh1Vz_VPDVz","",1000,-50,50);
	  mh1VzPS[iTrg] = new TH1F(Form("mh1VzPSTrg%i",iTrg),"",200,-120,120);     
	  //      mh1Vz_BBCMB=new TH1F(Form("mh1Vz_BBCMB"),"",200,-120,120);
	  mh2Vxy[iTrg]=new TH2F(Form("mh2VxyTrg%i",iTrg),"",200,0,0.3,200,0,0.1);
	  //for dataQA
	  mPhi_ptUnlike[iTrg]=new TH2F(Form("mPhi_ptUnlikeTrg%i",iTrg)," ",200,0,20,400,-4,4);
	  mEta_ptUnlike[iTrg]=new TH2F(Form("mEta_ptUnlikeTrg%i",iTrg)," ",200,0,20,300,-1.5,1.5);
	  mPhi_ptlike[iTrg]=new TH2F(Form("mPhi_ptlikeTrg%i",iTrg)," ",200,0,20,400,-4,4);
	  mEta_ptlike[iTrg]=new TH2F(Form("mEta_ptlikeTrg%i",iTrg)," ",200,0,20,300,-1.5,1.5);
	  
	  mHitFit_ptUnlike[iTrg] =new TH2F(Form("mHitFit_ptUnlikeTrg%i",iTrg)," ",200,0,20,60,0,60);
	  mHitFit_ptlike[iTrg] =new TH2F(Form("mHitFit_ptlikeTrg%i",iTrg)," ",200,0,20,60,0,60);
	  mNSMDEta_ptUnlike[iTrg]=new TH2F(Form("mNSMDEta_ptUnlikeTrg%i",iTrg)," ",200,0,20,10,0,10);
	  mNSMDEta_ptlike[iTrg]=new TH2F(Form("mNSMDEta_ptlikeTrg%i",iTrg)," ",200,0,20,10,0,10);
	  mNSMDPhi_ptUnlike[iTrg]=new TH2F(Form("mNSMDPhi_ptUnlikeTrg%i",iTrg)," ",200,0,20,10,0,10);
	  mNSMDPhi_ptlike[iTrg]=new TH2F(Form("mNSMDPhi_ptlikeTrg%i",iTrg)," ",200,0,20,10,0,10);
	  mHitsDedxUnlike[iTrg]=new TH2F(Form("mHitsDedxUnlikeTrg%i",iTrg)," ",200,0,20,60,0,60);
	  mHitsDedxlike[iTrg]=new TH2F(Form("mHitsDedxlikeTrg%i",iTrg)," ",200,0,20,60,0,60);
	  
	  mNTracklike[iTrg]=new TH1F(Form("mNTracklikeTrg%i",iTrg)," ",200,0,20);
	  mNTrackUnlike[iTrg]=new TH1F(Form("mNTrackUnlikeTrg%i",iTrg)," ",200,0,20);
	  
	  mNTrack_cutUnlike[iTrg]=new TH1F(Form("mNTrack_cutUnlikeTrg%i",iTrg)," ",200,0,20);
	  mNTrack_cutlike[iTrg]=new TH1F(Form("mNTrack_cutlikeTrg%i",iTrg)," ",200,0,20);     
	  
	  mNTrack_cut25Unlike[iTrg]=new TH1F(Form("mNTrack_cut25UnlikeTrg%i",iTrg)," ",200,0,20);
	  mNTrack_cut25like[iTrg]=new TH1F(Form("mNTrack_cut25likeTrg%i",iTrg)," ",200,0,20);
	  
	  mNTrack_cut20Unlike[iTrg]=new TH1F(Form("mNTrack_cut20UnlikeTrg%i",iTrg)," ",200,0,20);
	  mNTrack_cut20like[iTrg]=new TH1F(Form("mNTrack_cut20likeTrg%i",iTrg)," ",200,0,20);
	  
	  
	  mFitPos_ptlike[iTrg]=new TH2F(Form("mFitPos_ptlikeTrg%i",iTrg)," ",200,0,20,100,0,1);
	  mFitPos_ptUnlike[iTrg]=new TH2F(Form("mFitPos_ptUnlikeTrg%i",iTrg)," ",200,0,20,100,0,1);
	  
	  mgDcalike[iTrg]=new TH2F(Form("mgDcalikeTrg%i",iTrg),"",200,0,20,50,0,5);
	  mgDcaUnlike[iTrg]=new TH2F(Form("mgDcaUnlikeTrg%i",iTrg),"",200,0,20,50,0,5);
	  
	  mNsigElike[iTrg]=new TH2F(Form("mNsigElikeTrg%i",iTrg),"",200,0,20,50,-5,5);
	  mNsigEUnlike[iTrg]=new TH2F(Form("mNsigEUnlikeTrg%i",iTrg),"",200,0,20,50,-5,5);
	  
	  mDedxlike[iTrg]=new TH2F(Form("mDedxlikeTrg%i",iTrg),"",200,0,20,60,1,6);
	  mDedxUnlike[iTrg]=new TH2F(Form("mDedxUnlikeTrg%i",iTrg),"",200,0,20,60,1,6);
	  
	  mPoeUnlike[iTrg] = new TH2F(Form("mPoeUnlikeTrg%i",iTrg),"",100,0,2.5,200,0,20.);
	  mPoelike[iTrg] = new TH2F(Form("mPoelikeTrg%i",iTrg),"",100,0,2.5,200,0,20.);
	  //for dataQA
	}     
      //----------------------------------- low Pt----------------------------
      // low pt inclusive
      mh2nSigmaElec_VPD=new TH2F("mh2nSigmaElec_VPD","",400,-19.995,20.005,80 ,0,4.);
      mh3_nSigmaElec_Beta_VPD=new TH3F("mh3_nSigmaElec_Beta_VPD","",125,-10,15,80,-0.2,0.6,80,0,4.);
      mh2_nSigmaElec_BetaP_VPD=new TH2F("mh2_nSigmaElec_BetaP_VPD","",400,0,4.,600,0,3);
      mh2_nSigmaElec_BetaPt_VPD=new TH2F("mh2_nSigmaElec_BetaPt_VPD","",400,0,4.,600,0,3);
      
      mh2_Pion_nSigmaElec_VPD=new TH2F("mh2_Pion_nSigmaElec_VPD","",700,-34.995,35.005,80 ,0,4.);
      mh2_Kaon_nSigmaElec_VPD=new TH2F("mh2_Kaon_nSigmaElec_VPD","",700,-34.995,35.005,80 ,0,4.);
      mh2_Proton_nSigmaElec_VPD=new TH2F("mh2_Proton_nSigmaElec_VPD","",700,-34.995,35.005,80 ,0,4.);
      mh2_nsigamE_pt_VPD=new TH2F("mh2_nsigamE_pt_VPD","",700,-34.995,35.005,120 ,0,6.);
      mh2_nsigamE_p_VPD=new TH2F("mh2_nsigamE_p_VPD","",700,-34.995,35.005,120 ,0,6.);
      
      mh2_InvMass_KaonUnlike_VPD=new TH2F("mh2_InvMass_KaonUnlike_VPD","",100,0.98,1.08,80,0,4);
      mh2_NsigmaE_KaonUnlike_VPD=new TH2F("mh2_NsigmaE_KaonUnlike_VPD","",250,-10,15,80,0,4);
      mh2_InvMass_Kaonlike_VPD=new TH2F("mh2_InvMass_Kaonlike_VPD","",100,0.98,1.08,80,0,4);
      mh2_NsigmaE_Kaonlike_VPD=new TH2F("mh2_NsigmaE_Kaonlike_VPD","",250,-10,15,80,0,4);
      
      mh2_KaonnSigmaEle_eta_VPD=new TH2F("mh2_KaonnSigmaEle_eta_VPD","",80 ,0,4.,30,-1.5,1.5);
      mh2_ProtonnSigmaEle_eta_VPD=new TH2F("mh2_ProtonnSigmaEle_eta_VPD","",80,0,4.,30,-1.5,1.5);  
      mh2_PionnSigmaEle_eta_VPD=new TH2F("mh2_PionSigmaEle_eta_VPD","",80,0,4.,30,-1.5,1.5);  
      
      mh2_PionnSigmaElecDiff_VPD=new TH2F("mh2_PionnSigmaElecDiff_VPD","",80 ,0,4.,400,-20,20);
      mh2_MergrePionnSigmaElecDiff_VPD=new TH2F("mh2_MergePionnSigmaElecDiff_VPD","",80 ,0,4.,400,-20,20);
      mh2_KaonnSigmaElecDiff_VPD=new TH2F("mh2_KaonnSigmaElecDiff_VPD","",80 ,0,4.,400,-20,20);
      mh2_ProtonnSigmaElecDiff_VPD=new TH2F("mh2_ProtonnSigmaElecDiff_VPD","",80,0,4.,400,-20,20);
      
      mh2_InvMass_VPD=new TH2F("mh2_InvMass_VPD","",320,0,4,560,-0.2,1.2);
      mh2electronPtEta_VPD= new TH2F("mh2electronPtEta_VPD","",30,-1.5,1.5,80,0,4);
      mh1electronPt_VPD= new TH1F("mh1electronPt_VPD","",80,0,4);
      
      mh1electronPt_PSVPD= new TH1F("mh1electronPt_PSVPD","",80,0,4);
      mh1electronPt_MBPSVPD= new TH1F("mh1electronPt_MBPSVPD","",160,0,8);
      
      // low pt photonic no tofbeta
      mh2InvTofBetaUnlike_VPD=new TH2F("mh2InvTofBetaUnlike_VPD","",40,0,0.4,80,0,4);
      mh2InvTofBetalike_VPD=new TH2F("mh2InvTofBetalike_VPD","",40,0,0.4,80,0,4);
      
      mh3TofBetaUnlike_VPD=new TH3F("mh3TofBetaUnlike_VPD","",60,0.94,1.06,80,0,4,40,0,0.4);
      mh3TofBetalike_VPD=new TH3F("mh3TofBetalike_VPD","",60,0.94,1.06,80,0,4,40,0,0.4);
      
      mh3TofBeta_PartnerUnlike_VPD=new TH3F("mh3TofBeta_PartnerUnlike_VPD","",200,0.96,2,80,0,4,40,0,0.4);
      mh3TofBeta_Partnerlike_VPD=new TH3F("mh3TofBeta_Partnerlike_VPD","",200,0.96,2,80,0,4,40,0,0.4);
      
      mh3TofBeta_CutPartnerUnlike_VPD=new TH3F("mh3TofBeta_CutPartnerUnlike_VPD","",200,0.96,2,80,0,4,40,0,0.4);
      mh3TofBeta_CutPartnerlike_VPD=new TH3F("mh3TofBeta_CutPartnerlike_VPD","",200,0.96,2,80,0,4,40,0,0.4);
      
      // no  yloacal 
      mh2InvTofYlocalUnlike_VPD=new TH2F("mh2InvTofYlocalUnlike_VPD","",40,0,0.4,80,0,4);
      mh2InvTofYlocallike_VPD=new TH2F("mh2InvTofYlocallike_VPD","",40,0,0.4,80,0,4);
      
      mh2TofYlocalUnlike_VPD=new TH2F("mh2TofYlocalUnlike_VPD","",70,-3.5,3.5,80,0,4);
      mh2TofYlocallike_VPD=new TH2F("mh2TofYlocallike_VPD","",70,-3.5,3.5,80,0,4);
      // all the cuts applied 
      mh2InvUnlike_VPD=new TH2F("mh2InvUnlike_VPD","",40,0,0.4,80,0,4);
      mh2Invlike_VPD=new TH2F("mh2Invlike_VPD","",40,0,0.4,80,0,4);
      
      mh2InvUnlike_PSVPD=new TH2F("mh2InvUnlike_PSVPD","",40,0,0.4,80,0,4);
      mh2Invlike_PSVPD=new TH2F("mh2Invlike_PSVPD","",40,0,0.4,80,0,4);
      
      mh2InvPartnerUnlike_VPD=new TH2F("mh2InvPartnerUnlike_VPD","",40,0,0.4,80,0,4);
      mh2InvPartnerlike_VPD=new TH2F("mh2InvPartnerlike_VPD","",40,0,0.4,80,0,4);
      
      mh2CutPartnerUnlike_VPD=new TH2F("mh2CutPartnerUnlike_VPD","",40,0,0.4,80,0,4);
      mh2CutPartnerlike_VPD=new TH2F("mh2CutPartnerlike_VPD","",40,0,0.4,80,0,4);
      
      mh2nSigePartnerUnlike_VPD=new TH2F("mh2nSigePartnerUnlike_VPD","",100,-10,10,80,0,4);
      mh2nSigePartnerlike_VPD=new TH2F("mh2nSigePartnerlike_VPD","",100,-10,10,80,0,4);
      
      //for low pt Data   QA
      mPhi_ptUnlike_VPD=new TH2F("mPhi_ptUnlike_VPD","",40,-4,4,80,0,4);
      mPhi_ptlike_VPD=new TH2F("mPhi_ptlike_VPD","",40,-4,4,80,0,4);
      
      mEta_ptUnlike_VPD=new TH2F("mEta_ptUnlike_VPD","",30,-1.5,1.5,80,0,4);
      mEta_ptlike_VPD=new TH2F("mEta_ptlike_VPD","",30,-1.5,1.5,80,0,4);
      
      mTofBetaUnlike_VPD=new TH2F("mTofBetaUnlike_VPD","",60,0.94,1.06,80,0,4);
      mTofBetalike_VPD=new TH2F("mTofBetalike_VPD","",60,0.94,1.06,80,0,4);
      
      // Check the tof beta when apply the tofMachflag>0
      mTofBetaUnlike_match_VPD=new TH2F("mTofBetaUnlike_match_VPD","",200,-0.15,0.15,80,0,4);
      mTofBetalike_match_VPD=new TH2F("mTofBetalike_match_VPD","",100,-0.15,0.15,80,0,4);
      
      mTofYlocalUnlike_VPD=new TH2F("mTofYlocalUnlike_VPD","",70,-3.5,3.5,80,0,4);
      mTofYlocallike_VPD=new TH2F("mTofYlocallike_VPD","",70,-3.5,3.5,80,0,4);
      
      mHitFit_ptUnlike_VPD =new TH2F ("mHitFit_ptUnlike_VPD"," ",60,0,60,80,0,4);
      mHitFit_ptlike_VPD =new TH2F ("mHitFit_ptlike_VPD"," ",60,0,60,80,0,4);
      
      mHitsDedxUnlike_VPD=new TH2F ("mHitsDedxUnlike_VPD"," ",60,0,60,80,0,4);
      mHitsDedxlike_VPD=new TH2F ("mHitsDedxlike_VPD"," ",60,0,60,80,0,4);
      
      mgDcalike_VPD=new TH2F ("mgDcalike_VPD","",50,0,2.5,80,0,4);
      mgDcaUnlike_VPD=new TH2F ("mgDcaUnlike_VPD","",50,0,2.5,80,0,4);
      
      mNsigElike_VPD=new TH2F ("mNsigElike_VPD","",50,-5,5,80,0,4);
      mNsigEUnlike_VPD=new TH2F ("mNsigEUnlike_VPD","",50,-5,5,80,0,4);
      
      mDedxlike_VPD=new TH2F ("mDedxlike_VPD","",60,1,6,80,0,4);
      mDedxUnlike_VPD=new TH2F ("mDedxUnlike_VPD","",60,1,6,80,0,4);
      
      mNTrack_cutUnlike_VPD=new TH1F("mNTrack_cutUnlike_VPD","",80,0,4);
      mNTrack_cutlike_VPD=new TH1F("mNTrack_cutlike_VPD","",80,0,4);
      
      mFitPos_ptUnlike_VPD=new TH2F("mFitPos_ptUnlike_VPD","",50,0,1,80,0,4);
      mFitPos_ptlike_VPD=new TH2F("mFitPos_ptlike_VPD","",50,0,1,80,0,4);
      
      // for the purity study
      
      mh3nSigmaEUnlike_VPD=new TH3F ("mh3nSigmaEUnlike_VPD","",50,-5,5,80,0,4,8,0,0.4);
      mh3nSigmaElike_VPD=new TH3F ("mh3nSigmaElike_VPD","",50,-5,5,80,0,4,8,0,0.4);  
      }*/
}

//-----------------------------------------------------------------------------
void StNpeRead::writeObjects()
{
   mOutputFile->cd();

   /// Z. Histos
   for(Int_t trg=0; trg<4; trg++)
     {
       mh1PtAllTracks[trg]     -> Write();
       mh1PtETracks[trg]       -> Write();
       mh2nSigmaEPt[trg]       -> Write();
       mh2nSigmaEPt_eID[trg]   -> Write();
       mh2PoePt[trg]           -> Write();
       mh2nPhiPt[trg]          -> Write();
       mh2nEtaPt[trg]          -> Write();
       mh2PhiDistPt[trg]       -> Write();
       mh2ZDistPt[trg]         -> Write();
       mh2PhiQPt[trg]          -> Write();
       mh2TofPtAll[trg]        -> Write();
       mh2TofPtE[trg]          -> Write();
       mh2InvMassPtAll[trg]    -> Write();
       mh2InvMassPtE[trg]      -> Write();
       mh2InvMassPtUS[trg]     -> Write();
       mh2InvMassPtLS[trg]     -> Write();
       mh2nSigmaPionPt[trg]    -> Write();
       mh3DelPhiIncl[trg]      -> Write();
       mh3DelPhiPhotLS[trg]    -> Write();
       mh3DelPhiPhotUS[trg]    -> Write();
       mh3DelPhiPhotUSNP[trg]  -> Write();
       mh3DelPhiPhotLSNP[trg]  -> Write();
       mh3DelPhiPhotInclNP[trg]-> Write();
       mh3DelPhiInclWt[trg]    -> Write();
       mh3DelPhiPhotLSWt[trg]  -> Write();
       mh3DelPhiPhotUSWt[trg]  -> Write();
       mh3DelPhiHadHad[trg]    -> Write();
       mh1PtHadTracks[trg]     -> Write();

       // Pileup Histos
       mh3nTracksZdcx[trg]     -> Write();
       
       // Mixed Events
       mh3MixedDelPhi          -> Write();
       mh3MixedDelEta          -> Write();
       mh3MixedEtaPhi          -> Write();
     }
   
   /* for(Int_t ii=0; ii<5; ii++)
     {
       testHist1D[ii] -> Write();
       testHist2D[ii] -> Write();
       }*/

   /// X. Histos
   /* if(writeXiaozhiHists)
     {
       mh1Vz_BBCMB->Write();
       mh1Vz_BBCMBTOF0->Write();
       mh1VPDVz_BBCMBTOF0->Write();
       mh1VPDVz_BBCMB->Write();
       mh1Vz_BBCMBPS->Write();
       mh1Vz_BBCMBTOF0PS->Write();
       HT0_HT2->Write();
       mh1Vz_VPDVz->Write();
       mh1Vz_VPD->Write();
       mh1Vz_VPDPS->Write();
       
       
       for(Int_t iTrg=0;iTrg<4;iTrg++)
	 {    
	   // runQA
	   if(writeRunQA)
	     {
	       /*
	       mh2TPCVzRunID[iTrg]->Write();
	       mh2VPDVzRunID[iTrg]->Write();
	       mh2TPCVz_VPDVz[iTrg]->Write();
	       
	       mh2ZDCRunID[iTrg]->Write();
	       mh2BBCRunID[iTrg]->Write();
	       mh2MultiPosRunID[iTrg]->Write();
	       mh2MultiNegRunID[iTrg]->Write();
	       mh2MultiRunID[iTrg]->Write();
	       
	       mh2gDCARunID[iTrg]->Write();
	       mh2PtRunID[iTrg]->Write();
	       mh2EtaRunID[iTrg]->Write();
	       mh2PhiRunID[iTrg]->Write();
	       mh2BetaRunID[iTrg]->Write();
	       mh2Beta_RunID[iTrg]->Write();
	       mh2NhitFitRunID[iTrg]->Write();
	       mh2NhitDedxRunID[iTrg]->Write();
	       mh2DedxRunID[iTrg]->Write();
	       mh2nsigmaERunID[iTrg]->Write();
	       mh2PhiVsEta_Range1[iTrg]->Write();
	       mh2PhiVsEta_Range2[iTrg]->Write();
	       mh2PhiVsEta_Range3[iTrg]->Write();
	       mh1PassTofMatchRunID[iTrg]->Write();
	       mh1WPassTofMatchRunID[iTrg]->Write();
	       
	       mh2VXY[iTrg]->Write();
	       
	       // profile
	       mh2TPCVzRunID_Profile[iTrg]->Write();
	       mh2VPDVzRunID_Profile[iTrg]->Write();
	       	       
	       mh2ZDCRunID_Profile[iTrg]->Write();
	       mh2BBCRunID_Profile[iTrg]->Write();
	       mh2MultiPosRunID_Profile[iTrg]->Write();
	       mh2MultiNegRunID_Profile[iTrg]->Write();
	       mh2MultiRunID_Profile[iTrg]->Write();
	       
	       mh2gDCARunID_Profile[iTrg]->Write();
	       mh2PtRunID_Profile[iTrg]->Write();
	       mh2EtaRunID_Profile[iTrg]->Write();
	       mh2PhiRunID_Profile[iTrg]->Write();
	       mh2BetaRunID_Profile[iTrg]->Write();
	       mh2Beta_RunID_Profile[iTrg]->Write();
	       mh2NhitFitRunID_Profile[iTrg]->Write();
	       mh2NhitDedxRunID_Profile[iTrg]->Write();
	       mh2DedxRunID_Profile[iTrg]->Write();
	       mh2nsigmaERunID_Profile[iTrg]->Write();
	     }
	   
	   continue;
	   //for trigger sdudy
	   mh2DSMADC_Inclusive[iTrg]->Write();
	   mh2DSMADC_PhotoUnlike[iTrg]->Write();
	   mh2DSMADC_Photolike[iTrg]->Write();
	   
	   mh2ADC0_Inclusive[iTrg]->Write();
	   mh2ADC0_PhotoUnlike[iTrg]->Write();
	   mh2ADC0_Photolike[iTrg]->Write();
	   
	   
	   mh1TowerID_Inclusive[iTrg]->Write();
	   mh1TowerID_cut_Inclusive[iTrg]->Write();
	   mh1TowerID_all_Inclusive[iTrg]->Write();
	   
	   mh1TowerID_Inclusive_noHOt[iTrg]->Write();
	   mh1TowerID_cut_Inclusive_noHOt[iTrg]->Write();
	   mh1TowerID_all_Inclusive_noHOt[iTrg]->Write();
	   
	   mh2ADC0_DSMadcInclusive[iTrg]->Write();
	   mh2ADC0_DSMadcPhotoUnlike[iTrg]->Write();
	   mh2ADC0_DSMadcPhotolike[iTrg]->Write();
	   //
	   
	   mh2DsmADCUnlike[iTrg]->Write();
	   mh2DsmADClike[iTrg]->Write();
	   mh2Btowadc0[iTrg]->Write();
	   mh2InvMassPoeUnlike[iTrg]->Write();
	   mh2PoeUnlike[iTrg]->Write();
	   mh2InvMassPoelike[iTrg]->Write();  
	   mh2Poelike[iTrg]->Write();
	   mh2InvMassDzUnlike[iTrg]->Write();
	   mh2DzUnlike[iTrg]->Write();
	   mh2InvMassDzlike[iTrg]->Write();
	   mh2Dzlike[iTrg]->Write();
	   
	   mh2InvMassDpUnlike[iTrg]->Write();
	   mh2DpUnlike[iTrg]->Write();
	   mh2InvMassDplike[iTrg]->Write();
	   mh2Dplike[iTrg]->Write();
	   
	   mh2InvMassEMCUnlike[iTrg]->Write();
	   mh2InvMassEMClike[iTrg]->Write();
	   
	   mh2InvMassNEMCUnlike[iTrg]->Write();
	   mh2InvMassNEMClike[iTrg]->Write();
	   mh2InvMassNEEMCUnlike[iTrg]->Write();
	   mh2InvMassNEEMClike[iTrg]->Write();
	   mh2InvMassNPEMCUnlike[iTrg]->Write();
	   mh2InvMassNPEMClike[iTrg]->Write();
	   mh2NPhiUnlike[iTrg]->Write();
	   mh2NPhilike[iTrg]->Write();
	   mh2NEtaUnlike[iTrg]->Write();
	   mh2NEtalike[iTrg]->Write();
	   
	   mh3nSigPart_EMCUnlike[iTrg]->Write();
	   mh3nSigPart_EMClike[iTrg]->Write();
	   
	   mh3nSigPart_TREMCUnlike[iTrg]->Write();
	   mh3nSigPart_TREMClike[iTrg]->Write();
	   
	   mh3nSigPart_ADCTREMCUnlike[iTrg]->Write();
	   mh3nSigPart_ADCTREMClike[iTrg]->Write();
	   
	   mh3nSigEUnlike[iTrg]->Write();
	   mh3nSigPartUnlike[iTrg]->Write();
	   mh3nSigElike[iTrg]->Write();
	   mh3nSigPartlike[iTrg]->Write();
	   
	   mh3EMC_PartUnlike[iTrg]->Write();
	   mh3EMC_Partlike[iTrg]->Write();
	   mh3EMC_ADCPartUnlike[iTrg]->Write();
	   mh3EMC_ADCPartlike[iTrg]->Write();
	   mh3EtaPhiUnlike[iTrg]->Write();
	   mh3EtaPhilike[iTrg]->Write();     
	   mh2InvMasslike[iTrg]->Write();
	   mh2InvMassUnlike[iTrg]->Write();
	   mh2nSigmaElec[iTrg]->Write();
	   mh1electronPt[iTrg]->Write();
	   //  mh2Partner[Pt];
	   //for Event
	   mh1Vz[iTrg]->Write();
	   // mh1Vz_VPDVz->Write();
	   mh1VzPS[iTrg]->Write();     
	   //       mh1Vz_BBCMB->Write();
	   mh2Vxy[iTrg]->Write();
	   //dataQA
	   if(writeDataQA){
	     mPhi_ptUnlike[iTrg]->Write();
	     mEta_ptUnlike[iTrg]->Write();
	     mPhi_ptlike[iTrg]->Write();
	     mEta_ptlike[iTrg]->Write();
	     
	     mHitFit_ptUnlike[iTrg]->Write();
	     mNSMDEta_ptUnlike[iTrg]->Write();
	     mHitsDedxUnlike[iTrg]->Write();
	     mHitFit_ptlike[iTrg]->Write();
	     mNSMDEta_ptlike[iTrg]->Write();
	     mHitsDedxlike[iTrg]->Write();
	     
	     mNTrackUnlike[iTrg]->Write();
	     mNTrack_cutUnlike[iTrg]->Write();
	     mNTracklike[iTrg]->Write();
	     mNTrack_cutlike[iTrg]->Write();
	     mFitPos_ptlike[iTrg]->Write();
	     mgDcalike[iTrg]->Write();
	     mFitPos_ptUnlike[iTrg]->Write();
	     mgDcaUnlike[iTrg]->Write();
	     mNSMDPhi_ptUnlike[iTrg]->Write();
	     mNSMDPhi_ptlike[iTrg]->Write();
	     mNTrack_cut25Unlike[iTrg]->Write();
	     mNTrack_cut25like[iTrg]->Write();
	     mNTrack_cut20Unlike[iTrg]->Write();
	     mNTrack_cut20like[iTrg]->Write();
	     
	     mNsigElike[iTrg]->Write();
	     mNsigEUnlike[iTrg]->Write();
	     
	     mDedxlike[iTrg]->Write();
	     mDedxUnlike[iTrg]->Write();
	     mPoeUnlike[iTrg] ->Write();
	     mPoelike[iTrg] ->Write();
	   }
	 }
       
       //--------------------------------------------------low pt-------------------------
       mh2nSigmaElec_VPD->Write();
       mh3_nSigmaElec_Beta_VPD->Write();
       mh2_nSigmaElec_BetaP_VPD->Write();
       mh2_nSigmaElec_BetaPt_VPD->Write();
       
       mh2_PionnSigmaElecDiff_VPD->Write();
       mh2_KaonnSigmaElecDiff_VPD->Write();
       mh2_ProtonnSigmaElecDiff_VPD->Write();
       mh2_MergrePionnSigmaElecDiff_VPD->Write();
       
       
       mh2_nsigamE_pt_VPD->Write();
       mh2_nsigamE_p_VPD->Write();
       
       
       mh2_InvMass_KaonUnlike_VPD->Write();
       mh2_NsigmaE_KaonUnlike_VPD->Write();
       mh2_InvMass_Kaonlike_VPD->Write();
       mh2_NsigmaE_Kaonlike_VPD->Write();
       
       
       
       mh1electronPt_VPD->Write();
       mh1electronPt_PSVPD->Write();
       mh1electronPt_MBPSVPD->Write();
       mh2electronPtEta_VPD->Write();
       mh2_InvMass_VPD->Write();
       mh2_Pion_nSigmaElec_VPD->Write();
       mh2_Kaon_nSigmaElec_VPD->Write();
       mh2_Proton_nSigmaElec_VPD->Write();
       
       mh2_KaonnSigmaEle_eta_VPD->Write();
       mh2_ProtonnSigmaEle_eta_VPD->Write();
       mh2_PionnSigmaEle_eta_VPD->Write();
       // low pt photonic no tofbeta                                                                                                 
       mh2InvTofBetaUnlike_VPD->Write();
       mh2InvTofBetalike_VPD->Write();
              
       mh3TofBetaUnlike_VPD->Write();
       mh3TofBetalike_VPD->Write();
       
       mh3TofBeta_PartnerUnlike_VPD->Write();
       mh3TofBeta_Partnerlike_VPD->Write();
       
       mh3TofBeta_CutPartnerUnlike_VPD->Write();
       mh3TofBeta_CutPartnerlike_VPD->Write();
       
       
       // no  yloacal                                                                                                 
       
       mh2InvTofYlocalUnlike_VPD->Write();
       mh2InvTofYlocallike_VPD->Write();
       
       mh2InvPartnerUnlike_VPD->Write(); 
       mh2InvPartnerlike_VPD->Write();
       
       mh2TofYlocalUnlike_VPD->Write();
       mh2TofYlocallike_VPD->Write();
       // all the cuts applied                                                                                                       
       mh2CutPartnerUnlike_VPD->Write();
       mh2CutPartnerlike_VPD->Write();
       
       mh2InvUnlike_VPD->Write();
       mh2Invlike_VPD->Write();
       
       mh2InvUnlike_PSVPD->Write();
       mh2Invlike_PSVPD->Write();
       mh2nSigePartnerUnlike_VPD->Write();
       mh2nSigePartnerlike_VPD->Write();
       
       mh3nSigmaEUnlike_VPD->Write();
       mh3nSigmaElike_VPD->Write();
       
       // QA
       mPhi_ptUnlike_VPD->Write();
       mPhi_ptlike_VPD->Write();
       
       mEta_ptUnlike_VPD->Write();
       mEta_ptlike_VPD->Write();
       
       mTofBetaUnlike_VPD->Write();
       mTofBetalike_VPD->Write();
       mTofBetaUnlike_match_VPD->Write();
       mTofBetalike_match_VPD->Write();
       
       mTofYlocalUnlike_VPD->Write();
       mTofYlocallike_VPD->Write();
       
       mHitFit_ptUnlike_VPD->Write();
       mHitFit_ptlike_VPD->Write();
       
       mHitsDedxUnlike_VPD->Write();
       mHitsDedxlike_VPD->Write();
       
       mgDcalike_VPD->Write();
       mgDcaUnlike_VPD->Write();
       
       mNsigElike_VPD->Write();
       mNsigEUnlike_VPD->Write();
       
       mDedxlike_VPD->Write();
       mDedxUnlike_VPD->Write();
       
       mNTrack_cutUnlike_VPD->Write();
       mFitPos_ptUnlike_VPD->Write();
       
       mNTrack_cutlike_VPD->Write();
       mFitPos_ptlike_VPD->Write();
       }*/
   
   mOutputFile->Close();
}

//-----------------------------------------------------------------------------
void StNpeRead::read(TString fileName)
{  
  // Read Input File and Grab the Tree, setting Events to NpeEvent
  cout<<"input filename = "<<fileName.Data()<<endl;
  TFile* inFile = new TFile(fileName.Data(),"READ");
  TTree* tree = (TTree*)inFile->Get("T");
   
  tree->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  tree->SetBranchAddress("dEvent", &mNpeEvent);
   
  TClonesArray* aTracks = 0;
  TClonesArray* aPairs = 0;

  // events loop
  cout<<"Event"<<tree->GetEntriesFast()<<endl;
  for (UInt_t i = 0; i < tree->GetEntriesFast(); i++)
    {
      tree->GetEntry(i);
      
      //  Do_run_QA(mNpeEvent); 
      
      if(isBadRun(mNpeEvent)) continue;
      
      Int_t bTrg=999;
      if(!isGoodEvent(mNpeEvent,bTrg)) continue;
      TH1F::SetDefaultSumw2();
      // cout<< "pass Good Event"<<endl;
      
      // Do_run_QA(mNpeEvent);
      
      aPairs=mNpeEvent->electronPair();
      aTracks=mNpeEvent->tracks();
 
      // Sort based on trigger type (High Tower #, Beam-Beam Min Bias, Time of Flight, Vertex Position Min Bias)
      if(mNpeEvent->isHT0_BBCMB_TOF0())
	{ 
	  bTrg=0;
	  
	  Double_t  ps=mPrescales->GetPrescale(mNpeEvent->runId(),HT0BBCMBTOF0);
	  if(ps<0) continue;
	  zFillHists( bTrg, mNpeEvent, ps);
	  /* Fill_RunQA(bTrg,mNpeEvent);  
	  Fill_PhotonicE_hist ( bTrg, mNpeEvent , ps );
	  Fill_Inclusive_hist ( bTrg, mNpeEvent , ps );
	  Fill_event_hist(mNpeEvent,bTrg,ps);*/	
	}
      if(mNpeEvent->isHT1_BBCMB_TOF0())
	{
	  bTrg=1;
	  
	  Double_t  ps=mPrescales->GetPrescale(mNpeEvent->runId(),HT1BBCMBTOF0);
	  if(ps<0) continue;
	  zFillHists( bTrg, mNpeEvent, ps);
	  // Fill_RunQA(bTrg,mNpeEvent);  
	}
      if(mNpeEvent->isHT2_BBCMB())
	
	{
	  bTrg=2;
	  
	  Double_t  ps=mPrescales->GetPrescale(mNpeEvent->runId(),HT2BBCMB);
	  if(ps<0) continue;
	  zFillHists( bTrg, mNpeEvent, ps);

	  /* Fill_RunQA(bTrg,mNpeEvent);  
	  Fill_PhotonicE_hist ( bTrg, mNpeEvent , ps );
	  Fill_Inclusive_hist ( bTrg, mNpeEvent , ps );
	  Fill_event_hist(mNpeEvent,bTrg,ps);	*/
	}
      if(mNpeEvent->isVPDMB())
	{
	  bTrg=3;
	  
	  Double_t ps=mPrescales->GetPrescale(mNpeEvent->runId(),VPDMB);
	  if(ps<0) continue;
	  zFillHists( bTrg, mNpeEvent, ps);
	  /* Fill_RunQA(bTrg,mNpeEvent);  
	  Fill_PhotonicE_hist ( bTrg, mNpeEvent , ps );
	  Fill_Inclusive_hist ( bTrg, mNpeEvent , ps );
	  Fill_event_hist(mNpeEvent,bTrg,ps);
	  Fill_Kaon_Kaon(mNpeEvent,bTrg);*/
	} 
      // </ end trigger sort >
      
      if(i%1000==0) cout<< " Working on event: "<<i<<endl;
      
    } // ... end event loop
  
  inFile->Close();
  //zFillProjections();
}

/// ZACH for quick search

void StNpeRead::zFillHists (Int_t bTrg,StDmesonEvent * mNpeEvent ,Double_t ps ) // Zach Analysis Loops, only for HT triggers to match Dunkleberger Analysis (comparison)
{
  zFill_Inclusive(bTrg, mNpeEvent, ps);
  zFill_Photonic (bTrg, mNpeEvent, ps);
  //zFill_Pileup   (bTrg, mNpeEvent, ps); // Not Used
}

void StNpeRead::zFill_Inclusive (Int_t trg,StDmesonEvent * mNpeEvent ,Double_t ps ) // Fill histograms for events that are non-paired electrons                             
{
  Int_t pileupCounter = 0;
  
  TClonesArray* aTracks = 0;
  aTracks=mNpeEvent->tracks();
  TClonesArray* aPairs = 0;
  aPairs=mNpeEvent->electronPair();
  isAddedToBuffer = kFALSE;
  Double_t vz = mNpeEvent->primaryVertex().z();
  for(Int_t it=0;it<mNpeEvent->nTracks();it++)
    {
      Bool_t isInPair = kFALSE;
      pileupCounter = 0; // Clear for each event
      StDmesonTrack* trk = (StDmesonTrack*)aTracks->At(it);
      int  Run_ID=-1;
      map <int,int>::iterator iter=runID_List.find(mNpeEvent->runId());
      if(iter!=runID_List.end())
        Run_ID=iter->second;

      Float_t Phi = trk->gMom().phi();
      Float_t nSig = trk->nSigmaElectron();
      Float_t pT = trk->gMom().perp();
      Float_t q = trk->charge();
      Float_t beta = trk->btofBeta();
      Float_t p    = trk->gMom().mag();
      Float_t m_m  = p*p*(1/(beta*beta)-1);

      mh1PtAllTracks[trg] -> Fill(pT);
      mh2nSigmaEPt[trg]   -> Fill(nSig,pT);
      mh2TofPtAll[trg]    -> Fill(1/beta -1, pT);
      mh2InvMassPtAll[trg]-> Fill(m_m,pT);
      mh2nSigmaPionPt[trg]-> Fill(trk->nSigmaPion(),pT);
      if(isHotTower(trk,trg)) // If in a hot tower, don't do any other checks, skip track
	continue;
      if(pass_cut_GoodTrack(trk) && pass_cut_nsigmaE(trk) && pass_cut_Pt_Eta(trk) && pass_cut_ADC(trg,trk) &&
	 pass_cut_Match_EMC_Dz(trk) && pass_cut_Match_EMC_Dphi(trk) && pass_cut_poe(trk) && pass_cut_EMC(trk)) // Is track an electron and the ADC0 big enough to be the trigger?
	{
	  // if(!(trk->trgTowDsmAdc() < 0.1*trk->adc0())) // This matches the data analysis with embedding, according to Xiaozhi who did the embedding
	  //  continue;
	  // DEBUG cout << "!!There is NPE!!" << endl;
       
	  for(Int_t ip=0;ip<mNpeEvent->nElectronPair();ip++) // loop over pair branch to check if singleTrack also in Pair 
	    {
	      StElectronPair* pair = (StElectronPair*)aPairs->At(ip);
	      StDmesonTrack*  etrk = (StDmesonTrack*)aTracks->At(pair->electronId());
	      StDmesonTrack*  ptrk = (StDmesonTrack*)aTracks->At(pair->partnerId());
	      
	      if((trk == etrk || trk == ptrk) && pair->pairDca() < cuts::pairDCA && pair->m() < cuts::massDCA && etrk->charge()!=ptrk->charge()) // pair cuts + check for track id match
		{
		  // DEBUG cout << "Tracks are linked." << endl;
		  isInPair = kTRUE;
		}
	    }
	  // Mixed Event
	  if(!isAddedToBuffer) // Only compute it once per event, to avoid calculating with possibly correlated tracks 
	    {
	      computeMixedEvents(trk,vz);
	      isAddedToBuffer=kTRUE;
	    }

	  Float_t ePhi = Phi;
	  Float_t poe  = trk->gMom().mag()/trk->e0();
	  Float_t nPhi = trk->nPhi();
	  Float_t nEta = trk->nEta();
	  Float_t phiDist = trk->phiDist();
	  Float_t zDist   = trk->zDist();
	  Float_t epT  = pT;
	  Float_t eq   = q;
	  
	  mh1PtETracks[trg] -> Fill(epT);
	  mh2PhiDistPt[trg] -> Fill(phiDist,epT);
	  mh2ZDistPt[trg]   -> Fill(zDist,epT);
	  mh2nPhiPt[trg]    -> Fill(nPhi,epT);
	  mh2nEtaPt[trg]    -> Fill(nEta,epT);	  
	  mh2PoePt[trg]     -> Fill(poe,epT);
	  mh2nSigmaEPt_eID[trg] -> Fill(nSig,epT);
	  mh2TofPtE[trg]    -> Fill(1/beta -1, pT);
	  mh2InvMassPtE[trg]-> Fill(m_m,pT);

	  Int_t printCheck = 0;
	  for(Int_t ih = 0; ih < mNpeEvent->nTracks(); ih++) // Want to loop over all tracks looking for hads. Not going to double count, since there's only 1 NPE-e/evt on average (in events with NPE, which are rare)
	    {
	      StDmesonTrack* htrk = (StDmesonTrack*)aTracks->At(ih);
	      Float_t hpT   = htrk->pMom().perp();
	      
	      if(trk != htrk && pass_cut_hTrack(htrk)) // Is this track pass as hadron track quality AND not the same track
		{
	
		  // Make mixed events buffer
		  addToHadBuffer(htrk,vz);

		  Float_t hp    = htrk->pMom().mag();
		  Float_t hbeta = htrk->btofBeta();
		  Float_t hm_m  = p*p*(1/(beta*beta)-1);
		  Float_t hPhi  = htrk->pMom().phi();
		  Float_t hpoe  = htrk->pMom().mag()/trk->e0();
		  Float_t hq    = htrk->charge();
		  Float_t dPhi  = ePhi-hPhi;
		  Float_t hEta  = htrk->gMom().pseudoRapidity();
		  Float_t wt    = getHadronWt(hpT,hEta);
		  /* DEBUG if(printCheck < 20){
		    cout << "WEIGHT: " << wt << endl;
		    printCheck++;}*/
				   
		  if(dPhi > (3.*pi)/2.) dPhi = dPhi-2*pi;
		  if(dPhi < (-1*pi)/2.) dPhi = dPhi+2*pi;
		  mh2PhiQPt[trg]     -> Fill(hPhi,hq*hpT);
		  mh3DelPhiIncl[trg] -> Fill(dPhi,epT,hpT);
		  mh3DelPhiInclWt[trg] -> Fill(dPhi,epT,hpT,wt);
		  if(!isInPair)
		    mh3DelPhiPhotInclNP[trg] -> Fill(dPhi,epT,hpT);
		  
		}		  
	    }
	
	  Float_t hadptCuts[4]={0.2,0.5,1.0,1.5};
	  Float_t hpTCut = 0.2;
	  for(Int_t hC=0;hC<4;hC++) //  Study pileup as function of hpT cut
	    {
	      hpTCut = hadptCuts[hC];
	      pileupCounter = 0;
	      for(Int_t ih = 0; ih < mNpeEvent->nTracks(); ih++) 
		{
		  StDmesonTrack* htrk = (StDmesonTrack*)aTracks->At(ih);
		  Float_t hpT   = htrk->pMom().perp();
		  
		  if(trk != htrk && pass_cut_hTrack(htrk) && hpT > hpTCut)
		    {
		      pileupCounter++;
		    }
		}
	      mh3nTracksZdcx[trg]->Fill(mNpeEvent->ZDCx(),pileupCounter,hpTCut); // Fill on a per event basis
	    }
	}
      
      if(pass_cut_GoodTrack(trk) && pass_cut_nsigmaPi(trk) && pass_cut_Pt_Eta(trk)) // Is this track part of pure pion sample
        {
	  
	  Float_t phi = trk->pMom().phi();
	  Float_t pT  = trk->gMom().perp();
          mh1PtHadTracks[trg] -> Fill(pT);

	  for(Int_t ih = 0; ih < mNpeEvent->nTracks(); ih++) // Want to loop over all tracks looking for hads. Not going to double count, since there's only 1 NPE-e/evt on average \
	    (in events with NPE, which are rare)                                                                                                                                                 
	      {
		StDmesonTrack* htrk = (StDmesonTrack*)aTracks->At(ih);
		Float_t hpT   = htrk->pMom().perp();

		if(trk != htrk && pass_cut_hTrack(htrk)) // Is this track pass as hadron track quality AND not the same track              
		  {
		    Float_t hp    = htrk->pMom().mag();
		    Float_t hPhi  = htrk->pMom().phi();
		    Float_t dPhi  = phi-hPhi;
		    Float_t hEta  = htrk->gMom().pseudoRapidity();
		    Float_t wt    = getHadronWt(hpT,hEta);
		    /* DEBUG if(printCheck < 20){                                                                                                                                      
		       cout << "WEIGHT: " << wt << endl;                                                                                                                                
		       printCheck++;}*/

		    if(dPhi > (3.*pi)/2.) dPhi = dPhi-2*pi;
		    if(dPhi < (-1*pi)/2.) dPhi = dPhi+2*pi;
		    mh3DelPhiHadHad[trg] -> Fill(dPhi,pT,hpT);
		  }
	      }
	}
    }
}
  
void StNpeRead::zFill_Photonic (Int_t bTrg,StDmesonEvent * mNpeEvent ,Double_t ps )
{
  TClonesArray* aTracks = 0;
  TClonesArray* aPairs = 0;
  aPairs=mNpeEvent->electronPair();
  aTracks=mNpeEvent->tracks();

  for(Int_t ip=0;ip<mNpeEvent->nElectronPair();ip++)
    {
      StElectronPair* pair = (StElectronPair*)aPairs->At(ip);
      StDmesonTrack*  etrk = (StDmesonTrack*)aTracks->At(pair->electronId());
      StDmesonTrack*  ptrk = (StDmesonTrack*)aTracks->At(pair->partnerId());
    
      if(pair->pairDca() < cuts::pairDCA && pair->m() < cuts::massDCA) // Pair Conditions
	{ 
	  // if(pass_cut_Pt_Eta(etrk) && pass_cut_nsigmaE(etrk) && pass_cut_GoodTrack(etrk) && pass_cut_ADC(bTrg,etrk)) // Primary Track Conditions
	  // {
	  if(isHotTower(etrk,bTrg)) // If in a hot tower, don't do any other checks, skip track
	    continue;
	  if(pass_cut_GoodTrack(etrk) && pass_cut_nsigmaE(etrk) && pass_cut_Pt_Eta(etrk) && pass_cut_ADC(bTrg,etrk) &&
	     pass_cut_Match_EMC_Dz(etrk) && pass_cut_Match_EMC_Dphi(etrk) && pass_cut_poe(etrk) && pass_cut_EMC(etrk)) // Is track an electron and the ADC0 big enough to be the trigger?
	    {
	      // if(!(etrk->trgTowDsmAdc() < 0.1*etrk->adc0())) // This matches the data analysis with embedding, according to Xiaozhi who did the embedding
		//	continue;
	      
	      Float_t ep = etrk->gMom().mag();
	      Float_t ebeta = etrk->btofBeta();
	      Float_t em_m = ep*ep*(1/(ebeta*ebeta)-1);
	      Float_t ePhi = etrk->gMom().phi();
	      Float_t epoe = etrk->gMom().mag()/etrk->e0();
	      Float_t eq = etrk->charge();
	      Float_t epT = etrk->gMom().perp();
	      
	      Float_t pp = ptrk -> gMom().mag();
	      Float_t pbeta = ptrk -> btofBeta();
	      Float_t pm_m = pp*pp*(1/(pbeta*pbeta)-1);
	      Float_t pPhi = ptrk -> gMom().phi();
	      Float_t ppoe = ptrk -> gMom().mag()/ptrk->e0();
	      Float_t pq = ptrk -> charge();
	      
	      /// For Pair information sorting w/o Hadrons
	      if(eq == pq)
		{
		  mh2InvMassPtLS[bTrg]  -> Fill(pair->m(),epT);
		}
	      if(eq != pq)
		{
		  mh2InvMassPtUS[bTrg]  -> Fill(pair->m(),epT);
		}
	      
	      for(Int_t ih = ip; ih < mNpeEvent->nTracks(); ih++) // loop over all tracks in the event
		{
		  StDmesonTrack* htrk = (StDmesonTrack*)aTracks->At(ih);
		  Float_t hpT = htrk->pMom().perp();
		  if(etrk != htrk && pass_cut_hTrack(htrk)) // Is this track a hadron and not the same track or e-?
		    {
		      Float_t hPhi = htrk->pMom().phi();
		      Float_t dPhi  = ePhi-hPhi;
		      Float_t hEta  = htrk->gMom().pseudoRapidity();
		      Float_t wt = getHadronWt(hpT,hEta);
		      if(dPhi > (3.*pi)/2.) dPhi = dPhi-2*pi;
		      if(dPhi < -1*pi/2.) dPhi = dPhi+2*pi;
		      if(eq == pq)
			{
			  mh3DelPhiPhotLS[bTrg] -> Fill(dPhi,epT,hpT);
			  mh3DelPhiPhotLSWt[bTrg] -> Fill(dPhi,epT,hpT,wt);
			   if(ptrk != htrk)
			     mh3DelPhiPhotLSNP[bTrg] -> Fill(dPhi,epT,hpT);
			}
		      if(eq != pq)
			{
			  mh3DelPhiPhotUS[bTrg] -> Fill(dPhi,epT,hpT);
			  mh3DelPhiPhotUSWt[bTrg] -> Fill(dPhi,epT,hpT,wt);
			  if(ptrk != htrk)
			     mh3DelPhiPhotUSNP[bTrg] -> Fill(dPhi,epT,hpT);
			}
		    }
		}
	    }
	}
    }
}


// Fill Pileup Checking histos
void StNpeRead::zFill_Pileup (Int_t trg,StDmesonEvent * mNpeEvent ,Double_t ps ) // Fill histograms for events that are non-paired electrons
{
  // Get scalers of interest
  Float_t zdcx = mNpeEvent->ZDCx();
  Float_t bbcx = mNpeEvent->BBCx();
  // Get tracks info
  Float_t numTracks = mNpeEvent->nTracks();

  // Do etrack monitor filling by event when going through inclusive loop
}

/*void StNpeRead::Do_run_QA(StDmesonEvent * mNpeEvent) // Run Quality Assurance, compares values across all runIDs !!! UNUSED !!!
{
  //  cout<< " Do run QA"<<endl;
  return;
  if(!mNpeEvent) return;
  Float_t vz = mNpeEvent->primaryVertex().z();
  //  if (!(mNpeEvent->ranking() >0 )) continue;
  // if(fabs(vz) > cuts::vz) continue;
  // if(mNpeEvent->isVPDMB() && fabs(vz-mNpeEvent->vzVpd())>6) continue;
 
  Int_t bTrg=-1;
  if(mNpeEvent->isVPDMB())
    
    bTrg=0;
  
  //  Fill_RunQA(bTrg,RunID,trk,StDmesonEvent * mNpeEvent);
  if(mNpeEvent->isHT0_BBCMB_TOF0())
    bTrg=1;
  if(mNpeEvent->isHT1_BBCMB_TOF0())
    bTrg=2;
  if(mNpeEvent->isHT2_BBCMB())
    bTrg=3;
  if(bTrg==-1) return;
  
  int  Run_ID=-999;
  map <int,int>::iterator iter=runID_List.find(mNpeEvent->runId());
  if(iter!=runID_List.end())
    Run_ID=iter->second;

  
  // mh2ZDCRunID[bTrg]->Fill(Run_ID,mNpeEvent->ZDCx()/1000.);
  // mh2BBCRunID[bTrg]->Fill(Run_ID,mNpeEvent->BBCx()/1000.);
  // mh2MultiPosRunID[bTrg]->Fill(Run_ID,mNpeEvent->refMultPos());
  // mh2MultiNegRunID[bTrg]->Fill(Run_ID,mNpeEvent->refMultNeg());
  // mh2MultiRunID[bTrg]->Fill(Run_ID,mNpeEvent->refMult());
  // mh2TPCVzRunID[bTrg]->Fill(Run_ID,vz);
  // mh2VPDVzRunID[bTrg]->Fill(Run_ID,mNpeEvent->vzVpd());
  // mh2TPCVz_VPDVz[bTrg]->Fill(vz,mNpeEvent->vzVpd());

  /*mh2ZDCRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->ZDCx()/1000.);
  mh2BBCRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->BBCx()/1000.);
  mh2MultiPosRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->refMultPos());
  mh2MultiNegRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->refMultNeg());
  mh2MultiRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->refMult());
  mh2TPCVzRunID_Profile[bTrg]->Fill(Run_ID,vz);
  mh2VPDVzRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->vzVpd());
  mh2TPCVz_VPDVz[bTrg]->Fill(mNpeEvent->primaryVertex().z(),mNpeEvent->vzVpd());

  mh2VXY[bTrg]->Fill(mNpeEvent->primaryVertex().x(),mNpeEvent->primaryVertex().y());


  TClonesArray* aTracks = 0;
  aTracks=mNpeEvent->tracks();
  
  for(Int_t it=0;it<mNpeEvent->nTracks();it++)
    {
      StDmesonTrack* trk = (StDmesonTrack*)aTracks->At(it);
      
      int  Run_ID=-1;
      map <int,int>::iterator iter=runID_List.find(mNpeEvent->runId());
      if(iter!=runID_List.end())
        Run_ID=iter->second;
      if(!pass_loose_track_qaulity(trk,bTrg)) continue;

    // Run QA                                                                 
    /*  mh2gDCARunID_Profile[bTrg]->Fill(Run_ID,trk->dca());
      mh2PtRunID_Profile[bTrg]->Fill(Run_ID,trk->gMom().perp());
      mh2EtaRunID_Profile[bTrg]->Fill(Run_ID,trk->gMom().pseudoRapidity());
      mh2PhiRunID_Profile[bTrg]->Fill(Run_ID,trk->gMom().phi());
      mh2BetaRunID_Profile[bTrg]->Fill(Run_ID,1/trk->btofBeta());
      mh2DedxRunID_Profile[bTrg]->Fill(Run_ID,trk->dEdx());
      mh2nsigmaERunID_Profile[bTrg]->Fill(Run_ID,trk->nSigmaElectron());
      
      mh2Beta_RunID_Profile[bTrg]->Fill(Run_ID,trk->btofBeta());
      mh2NhitFitRunID_Profile[bTrg]->Fill(Run_ID,trk->nHitsFit());
      mh2NhitDedxRunID_Profile[bTrg]->Fill(Run_ID,trk->nHitsDedx());
    // Run QA                                                                 
      mh2gDCARunID[bTrg]->Fill(Run_ID,trk->dca());
      mh2PtRunID[bTrg]->Fill(Run_ID,trk->gMom().perp());
      mh2EtaRunID[bTrg]->Fill(Run_ID,trk->gMom().pseudoRapidity());
      mh2PhiRunID[bTrg]->Fill(Run_ID,trk->gMom().phi());
      mh2BetaRunID[bTrg]->Fill(Run_ID,1/trk->btofBeta());
      mh2DedxRunID[bTrg]->Fill(Run_ID,trk->dEdx());
      mh2nsigmaERunID[bTrg]->Fill(Run_ID,trk->nSigmaElectron());
      
      mh2Beta_RunID[bTrg]->Fill(Run_ID,trk->btofBeta());
      mh2NhitFitRunID[bTrg]->Fill(Run_ID,trk->nHitsFit());
      mh2NhitDedxRunID[bTrg]->Fill(Run_ID,trk->nHitsDedx());
      if(Run_ID<353)
	mh2PhiVsEta_Range1[bTrg]->Fill(trk->gMom().phi(),trk->gMom().pseudoRapidity());
      else if(353<=Run_ID && Run_ID<441)
	mh2PhiVsEta_Range2[bTrg]->Fill(trk->gMom().phi(),trk->gMom().pseudoRapidity());
      else
	mh2PhiVsEta_Range3[bTrg]->Fill(trk->gMom().phi(),trk->gMom().pseudoRapidity());
    
    }
    }*/

/*void StNpeRead::Fill_RunQA(Int_t bTrg,StDmesonEvent * mNpeEvent) // Run Quality (consistency) Check, Fill values vs RunID
{

  if(!mNpeEvent) return;

  int  Run_ID=-999;
  map <int,int>::iterator iter=runID_List.find(mNpeEvent->runId());
  if(iter!=runID_List.end())
    Run_ID=iter->second;

  // Fill Event based values
  /* mh2ZDCRunID[bTrg]->Fill(Run_ID,mNpeEvent->ZDCx()/1000.);
  mh2BBCRunID[bTrg]->Fill(Run_ID,mNpeEvent->BBCx()/1000.);
  mh2MultiPosRunID[bTrg]->Fill(Run_ID,mNpeEvent->refMultPos());
  mh2MultiNegRunID[bTrg]->Fill(Run_ID,mNpeEvent->refMultNeg());
  mh2MultiRunID[bTrg]->Fill(Run_ID,mNpeEvent->refMult());
  mh2TPCVzRunID[bTrg]->Fill(Run_ID,mNpeEvent->primaryVertex().z());
  mh2VPDVzRunID[bTrg]->Fill(Run_ID,mNpeEvent->vzVpd());
  mh2TPCVz_VPDVz[bTrg]->Fill(mNpeEvent->primaryVertex().z(),mNpeEvent->vzVpd());
  
  mh2ZDCRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->ZDCx()/1000.);
  mh2BBCRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->BBCx()/1000.);
  mh2MultiPosRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->refMultPos());
  mh2MultiNegRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->refMultNeg());
  mh2MultiRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->refMult());
  mh2TPCVzRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->primaryVertex().z());
  mh2VPDVzRunID_Profile[bTrg]->Fill(Run_ID,mNpeEvent->vzVpd());
  mh2TPCVz_VPDVz[bTrg]->Fill(mNpeEvent->primaryVertex().z(),mNpeEvent->vzVpd());

  mh2VXY[bTrg]->Fill(mNpeEvent->primaryVertex().x(),mNpeEvent->primaryVertex().y());

  // Loop over all tracks in the event
  TClonesArray* aTracks = 0;
  aTracks=mNpeEvent->tracks();
  for(Int_t it=0;it<mNpeEvent->nTracks();it++)
    {
      StDmesonTrack* trk = (StDmesonTrack*)aTracks->At(it);
   
      if(!(pass_cut_GoodTrack(trk) && 0.2<trk->gMom().perp() && fabs(trk->gMom().pseudoRapidity())<1)) continue;
      // DEBUG - if(trk->nHitsFit()<20) cout<< trk->nHitsFit()<< "  "<<endl;

      /* mh2gDCARunID_Profile[bTrg]->Fill(Run_ID,trk->dca());
      mh2PtRunID_Profile[bTrg]->Fill(Run_ID,trk->gMom().perp());
      mh2EtaRunID_Profile[bTrg]->Fill(Run_ID,trk->gMom().pseudoRapidity());
      mh2PhiRunID_Profile[bTrg]->Fill(Run_ID,trk->gMom().phi());
      mh2BetaRunID_Profile[bTrg]->Fill(Run_ID,1/trk->btofBeta());
      mh2DedxRunID_Profile[bTrg]->Fill(Run_ID,trk->dEdx());
      mh2nsigmaERunID_Profile[bTrg]->Fill(Run_ID,trk->nSigmaElectron());
      mh2Beta_RunID_Profile[bTrg]->Fill(Run_ID,trk->btofBeta());
      mh2NhitFitRunID_Profile[bTrg]->Fill(Run_ID,trk->nHitsFit());
      mh2NhitDedxRunID_Profile[bTrg]->Fill(Run_ID,trk->nHitsDedx());
      // Run QA                                                                 
      mh2gDCARunID[bTrg]->Fill(Run_ID,trk->dca());
      mh2PtRunID[bTrg]->Fill(Run_ID,trk->gMom().perp());
      mh2EtaRunID[bTrg]->Fill(Run_ID,trk->gMom().pseudoRapidity());
      mh2PhiRunID[bTrg]->Fill(Run_ID,trk->gMom().phi());
      mh2BetaRunID[bTrg]->Fill(Run_ID,1/trk->btofBeta());
      mh2DedxRunID[bTrg]->Fill(Run_ID,trk->dEdx());
      mh2nsigmaERunID[bTrg]->Fill(Run_ID,trk->nSigmaElectron());
      mh2Beta_RunID[bTrg]->Fill(Run_ID,trk->btofBeta());
      mh2NhitFitRunID[bTrg]->Fill(Run_ID,trk->nHitsFit());
      mh2NhitDedxRunID[bTrg]->Fill(Run_ID,trk->nHitsDedx());
      
      if(Run_ID<353)
	mh2PhiVsEta_Range1[bTrg]->Fill(trk->gMom().phi(),trk->gMom().pseudoRapidity());
      else if(353<=Run_ID && Run_ID<441)
	mh2PhiVsEta_Range2[bTrg]->Fill(trk->gMom().phi(),trk->gMom().pseudoRapidity());
      else
	mh2PhiVsEta_Range3[bTrg]->Fill(trk->gMom().phi(),trk->gMom().pseudoRapidity());
      

      if(pass_cut_Pt_Eta(trk) && 0<trk->pMom().perp()){
	mh1WPassTofMatchRunID[bTrg]->Fill(Run_ID);
	if(pass_cut_Tof_Match(trk)){
	  mh1PassTofMatchRunID[bTrg]->Fill(Run_ID);
	}
	}
      
    }
    }

void StNpeRead::Fill_Inclusive_hist (Int_t bTrg,StDmesonEvent * mNpeEvent ,Double_t ps ) // Fill histograms for events that are non-paired electrons
{
  TClonesArray* aTracks = 0;
  aTracks=mNpeEvent->tracks();
  
  for(Int_t it=0;it<mNpeEvent->nTracks();it++)
    {
      StDmesonTrack* trk = (StDmesonTrack*)aTracks->At(it);

      int  Run_ID=-1;
      map <int,int>::iterator iter=runID_List.find(mNpeEvent->runId());
      if(iter!=runID_List.end())
	Run_ID=iter->second;
      //      std::cout <<Run_ID<<endl;
      // // Run QA
      // mh2gDCARunID_VPD->Fill(Run_ID,trk->dca());
      // mh2PtRunID_VPD->Fill(Run_ID,trk->gMom().perp());
      // mh2EtaRunID_VPD->Fill(Run_ID,trk->gMom().pseudoRapidity());
      // mh2PhiRunID_VPD->Fill(Run_ID,trk->gMom().phi());
      // mh2BetaRunID_VPD->Fill(Run_ID,1/trk->btofBeta());
      // mh2DedxRunID_VPD->Fill(Run_ID,trk->dEdx());
      // mh2nsigmaERunID_VPD->Fill(Run_ID,trk->nSigmaElectron());

      // mh2Beta_RunID_VPD->Fill(Run_ID,trk->btofBeta());
      // mh2NhitFitRunID_VPD->Fill(Run_ID,trk->nHitsFit());
      // mh2NhitDedxRunID_VPD->Fill(Run_ID,trk->nHitsDedx());
      // // Run QA
      
      if(bTrg==3 && pass_cut_GoodTrack(trk) && pass_cut_Pt_Eta(trk) && pass_cut_Tof_Match(trk))
	{
	  Float_t M[3]={0.938,0.140,0.494};
	  Float_t p=trk->gMom().mag();
	  Float_t beta=trk->btofBeta();
	  Float_t m_m=p*p*(1/(beta*beta)-1);
	  //	  if(m_m>0.1)
	  //	  cout << " MMMM"<<m_m<<"   beta"<< beta<< " pp"<< p<<endl;
	
	  mh2_InvMass_VPD->Fill(trk->gMom().perp(),m_m);

	  if(m_m<0.9&&0.86<m_m && fabs(trk->nSigmaProton())<3)
	    {
	      mh2_Proton_nSigmaElec_VPD->Fill(trk->nSigmaElectron(),trk->gMom().perp());
	      mh2_ProtonnSigmaEle_eta_VPD->Fill(trk->gMom().perp(),trk->gMom().pseudoRapidity());
	    }
	  if(m_m<0.021&&0.016<m_m  && fabs(trk->nSigmaPion())<3)
	    {
	      mh2_Pion_nSigmaElec_VPD->Fill(trk->nSigmaElectron(),trk->gMom().perp());
	      mh2_PionnSigmaEle_eta_VPD->Fill(trk->gMom().perp(),trk->gMom().pseudoRapidity());
	    }
	  if(m_m<0.247&&0.237<m_m  && fabs(trk->nSigmaKaon())<3 )
	    {
	      mh2_Kaon_nSigmaElec_VPD->Fill(trk->nSigmaElectron(),trk->gMom().perp());
	      mh2_KaonnSigmaEle_eta_VPD->Fill(trk->gMom().perp(),trk->gMom().pseudoRapidity());
	    }
	}      
      
      if(bTrg==3 && pass_cut_GoodTrack(trk) && pass_cut_Pt_Eta(trk))
	{
	  mh2_PionnSigmaElecDiff_VPD->Fill(trk->gMom().perp(),trk->nSigmaElectron()-trk->nSigmaPion());
	  mh2_MergrePionnSigmaElecDiff_VPD->Fill(trk->gMom().perp(),trk->nSigmaElectron()-2*trk->nSigmaPion());
	  mh2_KaonnSigmaElecDiff_VPD->Fill(trk->gMom().perp(),trk->nSigmaElectron()-trk->nSigmaKaon());
	  mh2_ProtonnSigmaElecDiff_VPD->Fill(trk->gMom().perp(),trk->nSigmaElectron()-trk->nSigmaProton());
	}
      
      if(bTrg==3 && pass_cut_GoodTrack(trk) && pass_cut_Pt_Eta(trk))// eta cut debug
	{
	  float beta  =trk->btofBeta();
	  float pathl = beta*30.;
	  float pp    = trk->pMom().mag();
	  Double_t     dbeta = 1 - pathl/2.99792458e1*TMath::Sqrt(1-0.000511*0.000511/pp/pp);
	  mh3_nSigmaElec_Beta_VPD->Fill(trk->nSigmaElectron(),dbeta,trk->gMom().perp());
	  //if(pass_cut_nsigmaE(trk))
	  // {
	       //  float tof  = trk->btofBeta();
	  mh2_nSigmaElec_BetaPt_VPD->Fill(trk->gMom().perp(),1./beta);
	  mh2_nSigmaElec_BetaP_VPD->Fill(pp,1./beta);
	  // }
       	  if(pass_cut_Tof_Match(trk) && pass_cut_Tof(trk) )
	    {
	  // Fill low Pt histgram
	      mh2nSigmaElec_VPD->Fill(trk->nSigmaElectron(),trk->gMom().perp());
	      mh2_nsigamE_p_VPD->Fill(trk->nSigmaElectron(),trk->gMom().mag());
	      mh2_nsigamE_pt_VPD->Fill(trk->nSigmaElectron(),trk->gMom().perp());
	 	      
	      if( pass_cut_nsigmaE(trk))
		{
		  mh1electronPt_VPD->Fill(trk->gMom().perp());
		  mh1electronPt_PSVPD->Fill(trk->gMom().perp(),ps);
		  mh1electronPt_MBPSVPD->Fill(trk->gMom().perp(),ps);
		  mh2electronPtEta_VPD->Fill(trk->gMom().pseudoRapidity(),trk->gMom().perp());
		}
	    }
	}
      // Fill high pt

      if(!((bTrg==0    &&  trk->trgTowDsmAdc()>11 && trk->adc0()>180) || (bTrg==2  && trk->trgTowDsmAdc()>18 && trk->adc0()>300) || bTrg==3))  continue;
      //  if(!(bTrg==0 || bTrg==2  ))  continue;  
      if(pass_cut_GoodTrack(trk)&& pass_cut_Pt_Eta(trk) && pass_cut_Match_EMC_Dz(trk) && pass_cut_Match_EMC_Dphi(trk)  &&  pass_cut_poe(trk) &&  pass_cut_EMC(trk))	 
	{
	  mh2nSigmaElec[bTrg]->Fill(trk->nSigmaElectron(),trk->gMom().perp(),ps);
	  if( pass_cut_nsigmaE(trk))
	    {
	      mh1TowerID_all_Inclusive[bTrg]->Fill(trk->btowId());
	      if(!isHotTower(trk, bTrg))
		mh1TowerID_all_Inclusive_noHOt[bTrg]->Fill(trk->btowId());
	      if(trk->trgTowDsmAdc()>trk->adc0()*0.1)
		{
		  mh1TowerID_Inclusive[bTrg]->Fill(trk->btowId());
		  if(!isHotTower(trk, bTrg))	
		    mh1TowerID_Inclusive_noHOt[bTrg]->Fill(trk->btowId()); 
		}
	      if(trk->trgTowDsmAdc()<trk->adc0()*0.1)
		{
		  mh1TowerID_cut_Inclusive[bTrg]->Fill(trk->btowId());
		  if(!isHotTower(trk, bTrg))
		    {
		      mh1TowerID_cut_Inclusive_noHOt[bTrg]->Fill(trk->btowId());
		      mh1electronPt[bTrg]->Fill(trk->gMom().perp(),ps);
		      mh2DSMADC_Inclusive[bTrg]->Fill(trk->trgTowDsmAdc(),trk->gMom().perp());
		      mh2ADC0_Inclusive[bTrg]->Fill(trk->adc0(),trk->gMom().perp());
		      mh2ADC0_DSMadcInclusive[bTrg]->Fill(trk->adc0(),trk->trgTowDsmAdc());
		    }
		}
	    }
	}
    }                   
}
//-------------------Fill kaon kaon--
void StNpeRead::Fill_Kaon_Kaon(StDmesonEvent * mNpeEvent,Int_t bTrg)
{
  TClonesArray* aTracks = 0;
  TClonesArray* aPairs  = 0;
  aPairs=mNpeEvent->kaonKaon();
  aTracks=mNpeEvent->tracks();

  Float_t M[3]={0.938,0.140,0.494};
  
  for(Int_t ip=0;ip<mNpeEvent->nKaonKaon();ip++)
    {
      StKaonKaon* pair = (StKaonKaon*)aPairs->At(ip);
      StDmesonTrack* Kaon1Trk = (StDmesonTrack*)aTracks->At(pair->kaon1Id());
      StDmesonTrack* Kaon2Trk = (StDmesonTrack*)aTracks->At(pair->kaon2Id());
      
      Float_t p=Kaon1Trk->gMom().mag();
      Float_t beta=Kaon1Trk->btofBeta();
      Float_t m_m=p*p*(1/(beta*beta)-1);
      
      //   cout<< pair->m()<<endl;
      if(pass_cut_GoodTrack(Kaon1Trk)&& pass_cut_Pt_Eta(Kaon1Trk)&& fabs(Kaon1Trk->nSigmaKaon())<3 && fabs(m_m*m_m-M[2]*M[2]<0.06) )
	{
	  if(Kaon1Trk->charge()!=Kaon2Trk->charge())
	    {
	      mh2_InvMass_KaonUnlike_VPD->Fill(pair->m(),Kaon1Trk->gMom().perp());
	      
	    }
	  if(Kaon1Trk->charge()==Kaon2Trk->charge())
	    {
	      mh2_InvMass_Kaonlike_VPD->Fill(pair->m(),Kaon1Trk->gMom().perp());
	      
	    }
	}
      // --
      if(pass_cut_GoodTrack(Kaon1Trk)&& pass_cut_Pt_Eta(Kaon1Trk) && fabs(pair->m()-1.019)<0.004  && fabs(Kaon2Trk->nSigmaKaon())<3 && fabs(m_m*m_m-M[2]*M[2]<0.06))
	{
	  if(Kaon1Trk->charge()!=Kaon2Trk->charge())
	    {
	      mh2_NsigmaE_KaonUnlike_VPD->Fill(Kaon1Trk->nSigmaElectron(),Kaon1Trk->gMom().perp());
      
	    }
	  if(Kaon1Trk->charge()==Kaon2Trk->charge())
	    {
	      mh2_NsigmaE_Kaonlike_VPD->Fill(Kaon1Trk->nSigmaElectron(),Kaon1Trk->gMom().perp());
	      
	    }
	}
    }
}

//-----------------------------------
void StNpeRead::Fill_PhotonicE_hist (Int_t bTrg,StDmesonEvent * mNpeEvent ,Double_t ps )
{
  TClonesArray*   aTracks = 0;
  TClonesArray* aPairs=0;
  aPairs=mNpeEvent->electronPair();
  aTracks=mNpeEvent->tracks();

  for(Int_t ip=0;ip<mNpeEvent->nElectronPair();ip++)
    {
      StElectronPair* pair = (StElectronPair*)aPairs->At(ip);
      StDmesonTrack*  eTrk = (StDmesonTrack*)aTracks->At(pair->electronId());
      StDmesonTrack*  pTrk = (StDmesonTrack*)aTracks->At(pair->partnerId());
      // Fill low pt pair
      if(bTrg==3 && pair->m()<0.24 && pair->pairDca()< cuts::pairDCA && pass_cut_Pt_Eta(eTrk))
	Fill_pair_hist_low_pt(bTrg,pair,eTrk,pTrk,ps);
      if(isHotTower(eTrk, bTrg))continue;
      if(eTrk->trgTowDsmAdc()>eTrk->adc0()*0.1) continue;
      if(!((pair->m()<0.1) && pair->pairDca()<cuts::pairDCA  && pass_cut_Pt_Eta(eTrk)) )  continue;
      //---------------------------------------------------------------                                                                                  
      if((bTrg == 0  && eTrk->trgTowDsmAdc() > 11 && eTrk->adc0() > 180)  || (bTrg == 2  && eTrk->trgTowDsmAdc() > 18 && eTrk->adc0() > 300) || bTrg == 3)

      //  if(bTrg==0|| bTrg==2  )
	{
	  Fill_pair_hist_high_pt(bTrg,pair,eTrk,pTrk,ps);
	}
    }
}

void StNpeRead::Fill_pair_hist_low_pt(Int_t Trg,StElectronPair *pair,StDmesonTrack *eTrk,StDmesonTrack *pTrk, Double_t ps)
{
  // all cuts applied except tof eid beta 
  if( pass_cut_nsigmaE(eTrk) && pass_cut_Tof_Match(eTrk) && pass_cut_GoodTrack(eTrk) &&  pTrk->nSigmaElectron()<3 && -1.5<pTrk-> nSigmaElectron() )
    {
      if(eTrk->charge()!=pTrk->charge())
	{
	  mh2InvTofBetaUnlike_VPD->Fill(pair->m(),eTrk->gMom().perp());
	  mh3TofBetaUnlike_VPD->Fill(1/eTrk->btofBeta(),eTrk->gMom().perp(),pair->m());
	}   
      if(eTrk->charge()==pTrk->charge())   
	{
	  //
	  mh2InvTofBetalike_VPD->Fill(pair->m(),eTrk->gMom().perp());
	  mh3TofBetalike_VPD->Fill(1/eTrk->btofBeta(),eTrk->gMom().perp(),pair->m());
	}
    }
  // all cuts applied except tof match
  if(pass_cut_nsigmaE(eTrk) && pass_cut_Tof(eTrk) &&  pass_cut_GoodTrack(eTrk) &&  pTrk->nSigmaElectron()<3 && -3<pTrk-> nSigmaElectron() )
    {
      if(eTrk->charge()!=pTrk->charge())
        {
	  mh2InvTofYlocalUnlike_VPD->Fill(pair->m(),eTrk->gMom().perp());
	  mh2TofYlocalUnlike_VPD->Fill(eTrk->btofYLocal(),eTrk->gMom().perp());
	}
      if(eTrk->charge()==pTrk->charge())
        {
	  mh2InvTofYlocallike_VPD->Fill(pair->m(),eTrk->gMom().perp());
	  mh2TofYlocallike_VPD->Fill(eTrk->btofYLocal(),eTrk->gMom().perp());
	}
    }
  // all cuts applied except nsigmaE 
  if(pass_cut_Tof(eTrk) && pass_cut_Tof_Match(eTrk) && pass_cut_GoodTrack(eTrk) &&  pTrk->nSigmaElectron()<3 && -3<pTrk-> nSigmaElectron() )
    {
      if(eTrk->charge()!=pTrk->charge())
	{
	  mh3nSigmaEUnlike_VPD->Fill(eTrk->nSigmaElectron(),eTrk->gMom().perp(),pair->m());	
	}
      if(eTrk->charge()==pTrk->charge())
	{
	  mh3nSigmaElike_VPD->Fill(eTrk->nSigmaElectron(),eTrk->gMom().perp(),pair->m());    
	}
    }


  // all the cuts applied low pt
  if(pass_cut_nsigmaE(eTrk) && pass_cut_Tof(eTrk) && pass_cut_Tof_Match(eTrk) && pass_cut_GoodTrack(eTrk)  )
    {
      if(eTrk->charge()!=pTrk->charge())
        {
	  
	  mh2nSigePartnerUnlike_VPD->Fill(pTrk->nSigmaElectron(),eTrk->gMom().perp());
	  
	  if(pTrk->nSigmaElectron()<3 && -3<pTrk-> nSigmaElectron())
	    {
	      mh2InvUnlike_VPD->Fill(pair->m(),eTrk->gMom().perp());
	      
	      mh2InvUnlike_PSVPD->Fill(pair->m(),eTrk->gMom().perp(),ps);
	      
	      mh2InvPartnerUnlike_VPD->Fill(pair->m(),pTrk->gMom().perp());	       
	      
	      if(pass_cut_Tof_Match(pTrk) && 0<pTrk->pMom().perp())
		{
		  mh2CutPartnerUnlike_VPD->Fill(pair->m(),pTrk->gMom().perp());
		  
		  mh3TofBeta_CutPartnerUnlike_VPD->Fill(1/pTrk->btofBeta(),pTrk->gMom().perp(),pair->m());
		  //mTofBetaUnlike_match_VPD->Fill(btofBeta,pTrk->gMom().perp());
		}
	      mh3TofBeta_PartnerUnlike_VPD->Fill(1/pTrk->btofBeta(),pTrk->gMom().perp(),pair->m());
	      // QA
	      mTofBetaUnlike_VPD->Fill(1/eTrk->btofBeta(),eTrk->gMom().perp());
	      mTofYlocalUnlike_VPD->Fill(eTrk->btofYLocal(),eTrk->gMom().perp());
	      mHitFit_ptUnlike_VPD->Fill(eTrk->nHitsFit(),eTrk->gMom().perp());
	      mHitsDedxUnlike_VPD->Fill(eTrk->nHitsDedx(),eTrk->gMom().perp());
	      mNTrack_cutUnlike_VPD->Fill(eTrk->gMom().perp());
	      mFitPos_ptUnlike_VPD->Fill(eTrk->nHitsFit()/(Float_t)eTrk->nHitsMax(),eTrk->gMom().perp());
	      mgDcaUnlike_VPD->Fill(eTrk->dca(),eTrk->gMom().perp());
	      mPhi_ptUnlike_VPD->Fill(eTrk->gMom().phi(),eTrk->gMom().perp());
	      mEta_ptUnlike_VPD->Fill(eTrk->gMom().pseudoRapidity(),eTrk->gMom().perp());
	      mDedxUnlike_VPD->Fill(eTrk->dEdx(),eTrk->gMom().perp());
	      
	      mNsigEUnlike_VPD->Fill(eTrk->nSigmaElectron(),eTrk->gMom().perp());
	    }
    }
      if(eTrk->charge()==pTrk->charge())
	{
	  
	  mh2nSigePartnerlike_VPD->Fill(pTrk->nSigmaElectron(),eTrk->gMom().perp());
	  if(pTrk->nSigmaElectron()<3 && -3<pTrk-> nSigmaElectron())
	    {
	      mh2Invlike_VPD->Fill(pair->m(),eTrk->gMom().perp());
	      mh2Invlike_PSVPD->Fill(pair->m(),eTrk->gMom().perp(),ps);
	      mh2InvPartnerlike_VPD->Fill(pair->m(),pTrk->gMom().perp());
	      
	      if(pass_cut_Tof_Match(pTrk) && 0<pTrk->pMom().perp())
		{
		  mh2CutPartnerlike_VPD->Fill(pair->m(),pTrk->gMom().perp());
		  mh3TofBeta_CutPartnerlike_VPD->Fill(1/pTrk->btofBeta(),pTrk->gMom().perp(),pair->m());
		  // mTofBetalike_match_VPD->Fill(btofBeta,pTrk->gMom().perp());	 
		}
	      
	      mh3TofBeta_Partnerlike_VPD->Fill(1/pTrk->btofBeta(),pTrk->gMom().perp(),pair->m());
	      
	      mTofBetalike_VPD->Fill(1/eTrk->btofBeta(),eTrk->gMom().perp());
	      mTofYlocallike_VPD->Fill(eTrk->btofYLocal(),eTrk->gMom().perp());
	      mHitFit_ptlike_VPD->Fill(eTrk->nHitsFit(),eTrk->gMom().perp());
	      mHitsDedxlike_VPD->Fill(eTrk->nHitsDedx(),eTrk->gMom().perp());
	      mNTrack_cutlike_VPD->Fill(eTrk->gMom().perp());
	      mFitPos_ptlike_VPD->Fill(eTrk->nHitsFit()/(Float_t)eTrk->nHitsMax(),eTrk->gMom().perp());
	      mgDcalike_VPD->Fill(eTrk->dca(),eTrk->gMom().perp());
	      mPhi_ptlike_VPD->Fill(eTrk->gMom().phi(),eTrk->gMom().perp());
	      mEta_ptlike_VPD->Fill(eTrk->gMom().pseudoRapidity(),eTrk->gMom().perp());
	      mDedxlike_VPD->Fill(eTrk->dEdx(),eTrk->gMom().perp());
	      mNsigElike_VPD->Fill(eTrk->nSigmaElectron(),eTrk->gMom().perp());         
	    }
	}
    }
  
}// Fill_Pair_hist_low_pt

 void StNpeRead::Fill_pair_hist_high_pt(Int_t iTrg,StElectronPair * pair, StDmesonTrack* eTrk,StDmesonTrack * pTrk ,Double_t ps )    //Fill ht0
{
  if(pass_cut_nsigmaE(eTrk) && pass_cut_Match_EMC_Dz(eTrk) && pass_cut_Match_EMC_Dphi(eTrk) && pass_cut_EMC(eTrk)&& pass_cut_GoodTrack( eTrk))  //all cuts except poe
    {   Float_t poe=eTrk->gMom().mag()/eTrk->e0();
      
      if(eTrk->charge()!=pTrk->charge())
	{ 
	  //  cout<<"poe"<<poe<<endl;   
	  mh2InvMassPoeUnlike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2PoeUnlike[iTrg]->Fill(poe,eTrk->gMom().perp());   
	}
      if(eTrk->charge()==pTrk->charge()) 
	{
	  mh2InvMassPoelike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2Poelike[iTrg]->Fill(poe,eTrk->gMom().perp());
	}
    }
  if(pass_cut_nsigmaE(eTrk) && pass_cut_poe(eTrk) && pass_cut_Match_EMC_Dphi(eTrk) && pass_cut_GoodTrack( eTrk)&&  pass_cut_EMC(eTrk))  //all cuts except Zdist
    {   Float_t Dz=eTrk->zDist();
      if(eTrk->charge()!=pTrk->charge())
	{   // cout<<"Dz     "<<Dz<<endl;
	  mh2InvMassDzUnlike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2DzUnlike[iTrg]->Fill(Dz,eTrk->gMom().perp());   
	}
      if(eTrk->charge()==pTrk->charge()) 
	{
	  mh2InvMassDzlike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2Dzlike[iTrg]->Fill(Dz,eTrk->gMom().perp());
	}
    }
  
  if(pass_cut_nsigmaE(eTrk) && pass_cut_poe(eTrk) && pass_cut_Match_EMC_Dz(eTrk) && pass_cut_EMC(eTrk) && pass_cut_GoodTrack( eTrk))  //all cuts except phiDist 
    {   Float_t Dp=eTrk->phiDist();
      // cout<<"Dp    "<<Dp<<endl;
      if(eTrk->charge()!=pTrk->charge())
	{    
	  mh2InvMassDpUnlike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2DpUnlike[iTrg]->Fill(Dp,eTrk->gMom().perp());   
	}
      if(eTrk->charge()==pTrk->charge()) 
	{
	  mh2InvMassDplike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2Dplike[iTrg]->Fill(Dp,eTrk->gMom().perp());
	}
    }
  if(pass_cut_nsigmaE(eTrk) && pass_cut_poe(eTrk) && fabs( eTrk->zDist())<cuts::Dz_cut && fabs(eTrk->phiDist())<cuts::Dphi_cut && pass_cut_GoodTrack( eTrk) &&  1<eTrk->nEta() && 1<eTrk->nPhi())  //all cuts applied except neta nphi  
    { //  Float_t Dp=eTrk->phiDist();
      // cout<<"Dp    "<<Dp<<endl;
      if(eTrk->charge()!=pTrk->charge())
	{    
	  mh2InvMassNEMCUnlike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2NEtaUnlike[iTrg]->Fill(eTrk->nEta(),eTrk->gMom().perp());
	  mh2NPhiUnlike[iTrg]->Fill(eTrk->nPhi(),eTrk->gMom().perp());
	}
      if(eTrk->charge()==pTrk->charge()) 
	{
	  mh2InvMassNEMClike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2NEtalike[iTrg]->Fill(eTrk->nEta(),eTrk->gMom().perp());
	  mh2NPhilike[iTrg]->Fill(eTrk->nPhi(),eTrk->gMom().perp());
	}
    }
  if(pass_cut_nsigmaE(eTrk) && pass_cut_GoodTrack( eTrk));// && pair->m()<0.05)  // have no emc cut
  {
    Float_t DSMadc=eTrk->trgTowDsmAdc();//<pTrk->trgTowDsmAdc()?eTrk->trgTowDsmAdc():pTrk->trgTowDsmAdc();
    // cout<<"DSMadc"<<eTrk->trgTowDsmAdc()<<endl;
    Float_t nSigma_partner=pTrk->nSigmaElectron();
    if(eTrk->charge()!=pTrk->charge())
      {    
	mh2DSMADC_PhotoUnlike[iTrg]->Fill(DSMadc,eTrk->gMom().perp(),ps);
	
	mh2InvMassEMCUnlike[iTrg]->Fill(pair->m(),pTrk->gMom().perp(),ps);
	mh3nSigPart_TREMCUnlike[iTrg]->Fill(nSigma_partner,eTrk->gMom().perp(),pair->m());   
	mh3nSigPart_ADCTREMCUnlike[iTrg]->Fill(nSigma_partner,eTrk->trgTowDsmAdc(),pair->m()); 	
      }
    if(eTrk->charge()==pTrk->charge()) 
      {
	mh2DSMADC_Photolike[iTrg]->Fill(DSMadc,eTrk->gMom().perp(),ps);
	
	mh2InvMassEMClike[iTrg]->Fill(pair->m(),pTrk->gMom().perp(),ps);
	mh3nSigPart_TREMClike[iTrg]->Fill(nSigma_partner,eTrk->gMom().perp(),pair->m());
	mh3nSigPart_ADCTREMClike[iTrg]->Fill(nSigma_partner,eTrk->trgTowDsmAdc(),pair->m());
      }
  }
  if( pass_cut_poe(eTrk) && fabs(eTrk->zDist())<cuts::Dz_cut && fabs(eTrk->phiDist())<cuts::Dphi_cut && pass_cut_EMC(eTrk) && pass_cut_GoodTrack( eTrk))  // all cuts except nsigmaE
    {
      Float_t nSigmaE=eTrk->nSigmaElectron();
      // cout<<"nSigmaE"<<nSigmaE<<endl;
      Float_t nSigma_partner=pTrk->nSigmaElectron();
      if(eTrk->charge()!=pTrk->charge())
	{    
	  mh3nSigEUnlike[iTrg]->Fill(nSigmaE,eTrk->gMom().perp(),pair->m());
	  mh3nSigPartUnlike[iTrg]->Fill(nSigma_partner,eTrk->gMom().perp(),pair->m());   
	}
      if(eTrk->charge()==pTrk->charge()) 
	{
	  mh3nSigElike[iTrg]->Fill(nSigmaE,eTrk->gMom().perp(),pair->m());
	  mh3nSigPartlike[iTrg]->Fill(nSigma_partner,eTrk->gMom().perp(),pair->m());
	}
    }
  /*
    if( pass_cut_poe(eTrk) && fabs(eTrk->zDist())<cuts::Dz_cut && fabs(eTrk->phiDist())<cuts::Dphi_cut && pass_cut_nsigmaE(eTrk) &&  pass_cut_EMC(eTrk) && 25 <eTrk->nHitsFit() &&  15 <eTrk->nHitsDedx() &&  eTrk->dca()<1.5 && eTrk->firstTpcPointR()<73 )  //all cuts applied except good track quality cuts
    {
    if(eTrk->charge()!=pTrk->charge())
    mNTrack_cut25Unlike[iTrg]->Fill(eTrk->gMom().perp());
      if(eTrk->charge()==pTrk->charge())
      mNTrack_cut25like[iTrg]->Fill(eTrk->gMom().perp());
      
      }
  *//*
  if( pass_cut_poe(eTrk) && fabs(eTrk->zDist())<cuts::Dz_cut && fabs(eTrk->phiDist())<cuts::Dphi_cut && pass_cut_nsigmaE(eTrk) &&  pass_cut_EMC(eTrk) && 20 <eTrk->nHitsFit() &&  15 <eTrk->nHitsDedx() &&  eTrk->dca()<1.5 && eTrk->firstTpcPointR()<73)   //all cuts applied except good track quality cuts     
    {
      if(eTrk->charge()!=pTrk->charge())
	{	  mNTrack_cut20Unlike[iTrg]->Fill(eTrk->gMom().perp());
	  if(25 <eTrk->nHitsFit())
	    mNTrack_cut25Unlike[iTrg]->Fill(eTrk->gMom().perp());	 
	}
      if(eTrk->charge()==pTrk->charge())
	{
	  mNTrack_cut20like[iTrg]->Fill(eTrk->gMom().perp());
	  if(25 <eTrk->nHitsFit())
	    mNTrack_cut25like[iTrg]->Fill(eTrk->gMom().perp());
	  
	}
      
    }
  
  if( pass_cut_poe(eTrk) && fabs(eTrk->zDist())<cuts::Dz_cut && fabs(eTrk->phiDist())<cuts::Dphi_cut && pass_cut_nsigmaE(eTrk)     &&  pass_cut_EMC(eTrk) && pass_cut_GoodTrack( eTrk))  //all cuts applied
    if(pass_cut_nsigmaE(eTrk)&&pass_cut_GoodTrack(eTrk))  //all cuts applied
      
      { 
	Float_t nSigma_partner=pTrk->nSigmaElectron();
	// mh2DsmADC[iTrg]->Fill(eTrk->trgTowDsmAdc(),eTrk->gMom().perp());
	mh2Btowadc0[iTrg]->Fill(eTrk->adc0(),eTrk->gMom().perp());
	
	if(eTrk->charge()!=pTrk->charge())
	  { 
	    //	      mh3ADC_PhotoUnlike[iTrg]->Fill(eTrk->trgTowDsmAdc(),eTrk->adc0(),eTrk->gMom().perp());
	    mh2DsmADCUnlike[iTrg]->Fill(eTrk->trgTowDsmAdc(),eTrk->gMom().perp());
	    // mh2DSMADC_PhotoUnlike[iTrg]->Fill(eTrk->trgTowDsmAdc(),eTrk->gMom().perp());
	    mh2ADC0_PhotoUnlike[iTrg]->Fill(eTrk->adc0(),eTrk->gMom().perp());
	    mh2ADC0_DSMadcPhotoUnlike[iTrg]->Fill(eTrk->adc0(),eTrk->trgTowDsmAdc());
	    mh2InvMassUnlike[iTrg]->Fill(pair->m(),eTrk->gMom().perp(),ps);
	    //	mh2nSigPartUnlike[iTrg]->Fill(nSigma_partner,eTrk->gMom().perp());
	    mh3EMC_PartUnlike[iTrg]->Fill(nSigma_partner,eTrk->gMom().perp(),pair->m());
	    mh3EtaPhiUnlike[iTrg]->Fill(eTrk->gMom().pseudoRapidity(),eTrk->gMom().phi(),eTrk->gMom().perp());
	    mh3EMC_ADCPartUnlike[iTrg]->Fill(nSigma_partner,eTrk->trgTowDsmAdc(),pair->m());
	    //dataQA
	    ///、、   cout<<" eTrk->nHitsFit()"<<eTrk->nHitsFit()<<endl;
	    // cout<<"   eTrk->nHitsMax()"<<eTrk->nHitsMax()<<endl;
	    // cout<<"   eTrk->nHitsDedx()"<<eTrk->nHitsDedx()<<endl;
	    // cout<<"   eTrk->nHitsFit()/(Float_t)eTrk->nHitsMax()"<<eTrk->nHitsFit()/(Float_t)eTrk->nHitsMax()<<endl; 
	    
	    
	    mHitFit_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsFit());
	    mNSMDEta_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nEta());
	    
	    mHitsDedxUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsDedx());
	    //  mNTrackUnlike->[iTrg]->Fill(eTrk->gMom().perp());
	    mNTrack_cutUnlike[iTrg]->Fill(eTrk->gMom().perp());
	    mFitPos_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsFit()/(Float_t)eTrk->nHitsMax());
	    mgDcaUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->dca());
	    
	    mPhi_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->gMom().phi());
	    mEta_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->gMom().pseudoRapidity());
	    
	    mDedxUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->dEdx());
	    
	    
	    mPoeUnlike[iTrg]->Fill(eTrk->gMom().mag()/eTrk->e0(),eTrk->gMom().perp());
	  }
	if(eTrk->charge()==pTrk->charge()) 
	  {
	    // mh3ADC_Photolike[iTrg]->Fill(eTrk->trgTowDsmAdc(),eTrk->adc0(),eTrk->gMom().perp());
	    // mh2DSMADC_Photolike[iTrg]->Fill(eTrk->trgTowDsmAdc(),eTrk->gMom().perp());
	    //  mh2ADC0_Photolike[iTrg]->Fill(eTrk->adc0(),eTrk->gMom().perp());	    
	    mh2DsmADClike[iTrg]->Fill(eTrk->trgTowDsmAdc(),eTrk->gMom().perp());
	    mh2ADC0_DSMadcPhotolike[iTrg]->Fill(eTrk->adc0(),eTrk->trgTowDsmAdc());
	    mh2InvMasslike[iTrg]->Fill(pair->m(),eTrk->gMom().perp(),ps);
	    //  mh2nSigPartlike[iTrg]->Fill(nSigma_partner,eTrk->gMom().perp());
	    mh3EMC_Partlike[iTrg]->Fill(nSigma_partner,eTrk->gMom().perp(),pair->m());
	    mh3EtaPhilike[iTrg]->Fill(eTrk->gMom().pseudoRapidity(),eTrk->gMom().phi(),eTrk->gMom().perp());
	    mh3EMC_ADCPartlike[iTrg]->Fill(nSigma_partner,eTrk->trgTowDsmAdc(),pair->m());
	    mHitFit_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsFit());
	    
	    //	    mNSMDEta_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nEta());
	    mHitsDedxlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsDedx());
	    
	    mNTrack_cutlike[iTrg]->Fill(eTrk->gMom().perp());
	    mFitPos_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsFit()/(Float_t)eTrk->nHitsMax());
	    mgDcalike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->dca());	 
	    mPhi_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->gMom().phi());
	    mEta_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->gMom().pseudoRapidity());  
	    mDedxlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->dEdx()); 
	    mPoelike[iTrg]->Fill(eTrk->gMom().mag()/eTrk->e0(),eTrk->gMom().perp());
	  }
	
      }
      }*/

/*Bool_t StNpeRead:: isGoodEvent(StDmesonEvent * evt,Bool_t & bVPDMB,Bool_t & bHT1_BBCMB_TOF0 ,Bool_t& bHT2_BBCMB)
  
{
  if(!(evt->isVPDMB() || evt->isHT1_BBCMB_TOF0() || evt->isHT2_BBCMB() || evt->isBBCMB_TOF0() || evt->isBBCMB())) return kFALSE; 
  
  // reject bad runs
  if(badruns::is_bad_run(evt->runId(),-1)) return kFALSE;
  
  bVPDMB = kFALSE;
  bHT1_BBCMB_TOF0 = kFALSE;
  bHT2_BBCMB = kFALSE;
  
  Bool_t bPass = kTRUE;
  Float_t vz = evt->primaryVertex().z();
  if (!(evt->ranking() > 0)) bPass = kFALSE;
  if (fabs(vz) > cuts::vz) bPass = kFALSE;
  
  bVPDMB = evt->isVPDMB() && bPass && fabs(vz - evt->vzVpd()) < cuts::vzVpdVz;
  bHT1_BBCMB_TOF0 = evt->isHT1_BBCMB_TOF0() && bPass && !badruns::is_bad_run(evt->runId(),HT1BBCMBTOF0);
  bHT2_BBCMB = evt->isHT2_BBCMB() && bPass; // there is no other condition for this trigger
  
  bPass = bVPDMB || bHT1_BBCMB_TOF0 || bHT2_BBCMB;
  return bPass;
  }

//-----------------------------------------------------------------------------
Bool_t StNpeRead::isTrgTrack(StDmesonTrack* trk,UShort_t htTrg)
{
if(trk->trgTowDsmAdc()<=11) return kFALSE;

if(htTrg==0 && trk->trgTowDsmAdc()>11) return kTRUE;
if(htTrg==1 && trk->trgTowDsmAdc()>15) return kTRUE;
if(htTrg==2 && trk->trgTowDsmAdc()>18) return kTRUE;

   return kFALSE;
   }
*/

  /*void StNpeRead::Fill_event_hist( StDmesonEvent* Event,Int_t bTrg,Double_t ps)
{
  //  if(evt->isHT0_BBCMB_TOF0() || evt->isHT2_BBCMB())
  Float_t a=4;
  Float_t b=4; 
      if(Event->isHT0_BBCMB_TOF0()) a=0.5;
      else a=1.5;
      if(Event->isHT2_BBCMB()) b=0.5;
      else b=1.5;
       HT0_HT2->Fill(a,b);

       //      bVPDMB = evt->isVPDMB() && bPass && fabs(vz - evt->vzVpd()) 

  // Double_t ps=0; 
  //Event
  Float_t Vz=Event->primaryVertex().z();
  Float_t Vy=Event->primaryVertex().y();
  Float_t Vx=Event->primaryVertex().x();
 


   if(bTrg==0)
   ps=mPrescales->GetPrescale(Event->runId(),HT0BBCMBTOF0);

  if(bTrg==2)
   ps=mPrescales->GetPrescale(Event->runId(),HT2BBCMB);

  if(ps>0)
    {  
  // cout<<" btrg  "<<bTrg<<"ps   "<<ps<<endl;
  mh1Vz[bTrg]->Fill(Vz);
  mh1VzPS[bTrg]->Fill(Vz,ps);
 mh2Vxy[bTrg]->Fill(Vx,Vy);
    }
}*/


//eID cuts
Bool_t StNpeRead::pass_cut_nsigmaE(StDmesonTrack* trk)
{ 
  if(!trk) return kFALSE;
  Float_t nSigmaE=trk->nSigmaElectron();
  
      if(nSigmaE>cuts::nsigmae_low && nSigmaE<cuts::nsigmae_high)  
	return kTRUE;
  else return kFALSE;
}

//Select pure pion sample                                                                                                        
Bool_t StNpeRead::pass_cut_nsigmaPi(StDmesonTrack* trk)
{
  if(!trk) return kFALSE;
  Float_t nSigmaP=trk->nSigmaPion();
   
  if(nSigmaP>cuts::nsigmap_low && nSigmaP<cuts::nsigmap_high)
    return kTRUE;
  else return kFALSE;
 }
  
Bool_t StNpeRead::pass_cut_ADC(Int_t trg, StDmesonTrack* trk)
{
  // Different cuts for different trigs (not in StCuts, since not constant)
  Float_t adc0Cuts[4] = {180.,180.,300.,0.};
  Float_t dsmAdcCuts[4] = {11.,11.,18.,0.};
  if(trg == 3) // ADC cuts not for MB events
    return kTRUE;

  if(trk->adc0() > adc0Cuts[trg] && trk->trgTowDsmAdc() > dsmAdcCuts[trg])
    return kTRUE;
  else return kFALSE;
}

Bool_t StNpeRead::pass_cut_poe(StDmesonTrack * trk)
{
  Float_t P=trk->gMom().mag();
  Float_t E=trk->e0();
  //cout<<""
  if(P/E>cuts::poe_low && P/E<cuts::poe_high && cuts::Pt_cut<P)
    return kTRUE;
  else return kFALSE;
}

Bool_t StNpeRead::pass_cut_hTrack(StDmesonTrack * trk)
{
  Float_t Pt=trk->gMom().perp();
  Float_t eta=trk->gMom().pseudoRapidity();
  Int_t nhitDedx=trk->nHitsDedx();
  Int_t nhitsFit=trk->nHitsFit();
  Int_t nhitsMax=trk->nHitsMax();
  Float_t dca = trk->dca();
  //cout<<""
  if(Pt > cuts::hadPtMin && fabs(eta) < cuts::hadEta && nhitDedx >= cuts::hadDedx && nhitsFit >= cuts::hadHitsFit && dca < cuts::hadDCA && dca > 0)
    return kTRUE;
  else return kFALSE;
}

Bool_t StNpeRead::pass_loose_track_qaulity( StDmesonTrack *trk, Int_t bTrg)
{
  Float_t Pt=trk->gMom().perp();
  Float_t eTa=trk->gMom().pseudoRapidity();
  Int_t nhitDedx=trk->nHitsDedx();
  Int_t nhitsFit=trk->nHitsFit();
  Int_t nhitsMax=trk->nHitsMax();
  Float_t gDca=trk->dca();  

  if(gDca<3 && fabs(eTa)<1 && 0.2<Pt && fabs(trk->charge())==1 && (bTrg!=0 || ((bTrg==0) &&  trk->btofMatchFlag()>0)))
       return kTRUE;
        
     else return kFALSE;
}
Bool_t StNpeRead::pass_cut_GoodTrack(StDmesonTrack * trk)
{
  
  Float_t Pt=trk->gMom().perp();
  Float_t eTa=trk->gMom().pseudoRapidity();
  Int_t nhitDedx=trk->nHitsDedx();
  Int_t nhitsFit=trk->nHitsFit();
  Int_t nhitsMax=trk->nHitsMax();
  
  
  //  cout<<"pt    "<<Pt<<endl;
  // cout<<"nEta"<<nEta<<endl;
  // cout<<"nPh1"<<nPhi<<endl;
  
  // Float_t Dz=trk->zDist();
  //   Float_t Dphi=trk->phiDist();
  // cout<<"Dz  "<<Dz<<endl;
  // cout<<"DPhi"<<Dphi<<endl;
  Float_t gDca=trk->dca();  
  Float_t first_point= trk->firstTpcPointR();
  //cout<<"first"<<Ds_fst<<endl;
  if(gDca<cuts::gDca && cuts::nhitDedx<nhitDedx &&  nhitsFit>cuts::nhitsFit && (Float_t) nhitsFit/nhitsMax>cuts::nhit_Possmax  &&  first_point<cuts::first_point_cut)
    
    return kTRUE;
  else return kFALSE;
}

Bool_t StNpeRead::pass_cut_Pt_Eta(StDmesonTrack *trk)
{
  Float_t Pt=trk->gMom().perp();
  Float_t eTa=trk->gMom().pseudoRapidity();
  if(fabs(eTa)<cuts::eta_low && Pt>cuts::pt)
    return kTRUE;
  else return kFALSE;
}
Bool_t StNpeRead::pass_cut_EMC(StDmesonTrack *trk)
{
  Int_t nEta=trk->nEta();
  Int_t nPhi=trk->nPhi();
  if(nEta>cuts::NEta && nPhi>cuts::NPhi && cuts::Pt_cut<=trk->gMom().perp())
    return kTRUE;
  else return kFALSE;
}

Bool_t StNpeRead::is_EMC_Track(StDmesonTrack *trk)
{
  Int_t nEta=trk->nEta();
  Int_t nPhi=trk->nPhi();
  if(nEta > 0 && nPhi > 0 && cuts::Pt_cut<=trk->gMom().perp())
    return kTRUE;
  else return kFALSE;
}

Bool_t StNpeRead::pass_cut_Match_EMC_Dz(StDmesonTrack * trk){

  //  fabs(eTrk->zDist())<cuts::Dz_cut && fabs(eTrk->phiDist())<cuts::Dphi_cut
  //m/
  if(cuts::Dz_cut>fabs(trk->zDist()) && cuts::Pt_cut<=trk->gMom().perp())
    return kTRUE;
  else return kFALSE;
}
Bool_t StNpeRead::pass_cut_Match_EMC_Dphi(StDmesonTrack * trk){
  
  //
  if(cuts::Dphi_cut>fabs(trk->phiDist()) && cuts::Pt_cut<=trk->gMom().perp())
    return kTRUE;
  else return kFALSE;
}
Bool_t StNpeRead::pass_cut_Tof_Match( StDmesonTrack *trk){

  //
  Float_t tof_YTlocal=trk->btofYLocal();
  // trk->btofBeta()>0
  if(fabs(tof_YTlocal)<cuts::tof_Ylocal &&  trk->btofMatchFlag()>0 && trk->btofBeta()>0)
    return kTRUE;
  else  return kFALSE;
}
Bool_t StNpeRead:: pass_cut_Tof( StDmesonTrack *trk){
  //

  Float_t tofbeta=trk->btofBeta();
  if(fabs(1/tofbeta-1)<cuts::tofbeta)

  return kTRUE;

  else return kFALSE;
}

//-----------------------------------------------------------------------------

Bool_t StNpeRead::isHotTower(StDmesonTrack *trk,Int_t bTrg)
{

  // Int_t  HT0_Hot_towerlist[] ={32,52,115,246,268,276,294,386,510,562,682,750,773,894,987,994,1043,1064,1143,1233,1264,1285,1307,1487,1593,1710,1733,1823,1824,1851,1946,2022,2044,2064,2110,2146,2163,2165,2203,2291,2314,2522,2530,2634,2653,2835,2864,2866,2973,3006,3062,3533,3545,3727,3862,3949,4051,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767};

  // Int_t HT2_Hot_towerlist[]={32,52,115,268,276,294,295,510,562,682,750,987,994,1064,1143,1233,1264,1285,1307,1487,1593,1733,1824,1851,1946,2044,2163,2203,2634,2653,2835,2864,2866,2973,3006,3727,3862,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767};


  //ht0 hot tower
  //  if(bTrg==0)
 // cout<<sizeof(HT0_Hot_towerlist)/sizeof(HT0_Hot_towerlist[0])<<endl;
 // cout<<sizeof(HT2_Hot_towerlist)/sizeof(HT2_Hot_towerlist[0])<<endl;     

  //  Int_t  Hot_towerlist[] ={32,52,115,246,268,276,294,295,386,510,562,682,750,773,894,987,994,1043,1064,1143,1233,1264,1285,1307,1487,1593,1710,1733,1823,1824,1851,1946,2022,2044,2064,2110,2146,2163,2165,2203,2291,2314,2522,2530,2634,2653,2835,2864,2866,2973,3006,3062,3533,3545,3691,3727,3862,3949,4051,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767};


  Int_t  Hot_towerlist[]={22,30,31,114,251,275,308,509,552,555,681,682,691,749,772,1063,1263,1268,1284,1304,1306,1329,1486,1592,1709,1768,1770,1823,1882,1904,1909,1945,2022,2042,2043,2048,2051,2067,2145,2146,2162,2171,2190,2202,2272,2288,2290,2493,2504,2522,2529,2549,2706,2723,2712,2750,2863,2865,2868,2874,2952,3061,3007,3061,3092,3112,3154,3264,3166,3269,3271,3307,3326,3330,3331,3349,3365,3373,3532,3544,3626,3692,3821,3861,3932,3964,4130,4169,4226,4232,4249,4252,4262,4353,4430,4546,4749,4727,4766};


 for(Int_t i=0;i<sizeof(Hot_towerlist)/sizeof(Hot_towerlist[0]);i++ )
	{
	  if(Hot_towerlist[i]==trk->btowId())

	    return kTRUE;
	}
 return kFALSE;
}

Bool_t StNpeRead::isBadRun(StDmesonEvent * evt)
{
  Int_t bTrg=-999;
  if(evt->isVPDMB())
    { bTrg=3; if(badruns::is_bad_run(evt->runId(),bTrg)) 

	
      return kTRUE;}
  if(evt->isHT0_BBCMB_TOF0())
    { bTrg=0;if(badruns::is_bad_run(evt->runId(),bTrg)) return kTRUE;}
  if(evt->isHT1_BBCMB_TOF0())
    { bTrg=1;if(badruns::is_bad_run(evt->runId(),bTrg)) return kTRUE;}
  if(evt->isHT2_BBCMB())
    { bTrg=2;if(badruns::is_bad_run(evt->runId(),bTrg)) return kTRUE;}
  // if(badruns::is_bad_run(evt->runId(),bTrg)) return kTRUE;
  else return kFALSE;
  
}
Bool_t StNpeRead::isGoodEvent(StDmesonEvent* evt,Int_t bTrg)
{
  if(!evt) return kFALSE;
  // if(badruns::is_bad_run(evt->runId(),bTrg)) return kFALSE;
  //  outfile<<evt->runId()<<endl;
  // runIDPico.insert(evt->runId());
  // cuts
  Float_t vz = evt->primaryVertex().z();
  if (!(evt->ranking() >0 )) return kFALSE;
  if(fabs(vz) > cuts::vz) return kFALSE;
  if(evt->isVPDMB() && fabs(vz-evt->vzVpd())>6) return kFALSE;
  
  
    
   if(evt->isVPDMB())
     {
  
       Double_t ps_VPDMB=mPrescales->GetPrescale(evt->runId(),VPDMB);

       if(ps_VPDMB>0)
	 {		 
	   //mh1Vz_VPDVz->Fill(vz-evt->vzVpd());
	   //mh1Vz_VPD->Fill(vz);
	   //mh1Vz_VPDPS->Fill(vz,ps_VPDMB);
	 }
  
     }

   if(evt->isBBCMB() ||evt->isHT2_BBCMB())
     { 
  
       Double_t ps_BBCMB=mPrescales->GetPrescale(evt->runId(),BBCMB);
       Double_t ps_HT2BBCMB=mPrescales->GetPrescale(evt->runId(),HT2BBCMB);
       if( ps_BBCMB>0 && ps_HT2BBCMB>0)                
	 {
	   if(evt->isBBCMB())
	     {
	       //mh1Vz_BBCMB->Fill(vz);
	       //mh1VPDVz_BBCMB->Fill(evt->vzVpd());
	       //mh1Vz_BBCMBPS->Fill(vz,ps_BBCMB);
	     }	 
	 }    
  
     }


   if(evt->isBBCMB_TOF0() ||evt->isHT0_BBCMB_TOF0())
     {
  
       Double_t ps_HT0_BBCMB_TOF0=mPrescales->GetPrescale(evt->runId(),HT0BBCMBTOF0);
       Double_t ps_BBCMB_TOF0=mPrescales->GetPrescale(evt->runId(),BBCMBTOF0);
       if( ps_BBCMB_TOF0>0 && ps_HT0_BBCMB_TOF0>0)
	 {
	   if( evt->isBBCMB_TOF0())     
	     {  
	       //mh1Vz_BBCMBTOF0->Fill(vz);
	       //mh1VPDVz_BBCMBTOF0->Fill(evt->vzVpd());
	       //mh1Vz_BBCMBTOF0PS->Fill(vz,ps_BBCMB_TOF0);
	     }	

	 }
     }

   // if(evt->isHT1_BBCMB_TOF0())
   //   std::cout<<" HT1"<<std::endl;
   if(evt->isHT0_BBCMB_TOF0() || evt->isHT2_BBCMB() || evt->isVPDMB() || evt->isHT1_BBCMB_TOF0())
     // if(evt->isVPDMB())	
     // if( evt->isBBCMB_TOF0()||evt->isBBCMB())
     
     
     return kTRUE;
   else 
     return kFALSE; 
}

Double_t StNpeRead::getHadronWt(Double_t pt, Double_t eta){
	Int_t etabin = floor((eta+1)/0.1);
	if(etabin<0||etabin>=20) return 0.;
	fEff->SetParameters(effPars[etabin]);
	Double_t wt = fEff->Eval(pt);
	if(wt>0) return 1./wt;
	else return 0.;
}

 Int_t StNpeRead::readEff()
 {
   ifstream inf("/global/homes/z/zamiller/NPEhPhiAnalysis2015/collabCode/bingchu_comb_eff_fit.txt",std::ifstream::in);
   if(!inf.good()){ 
     cout<<"No input efficiency tables!"<<endl;
     return -1;
   }
   for(int i=0;i<20;i++){
     for(int j=0;j<3;j++) inf>>effPars[i][j];
   }
   inf.close();
   fEff = new TF1("fEff","[0]*exp(-pow([1]/x,[2]))",0.2,15);
   return 1;
 }

void StNpeRead::addToHadBuffer(StDmesonTrack *trk, Double_t vz)
 {
   Int_t vzbin = (Int_t)vz+35; // add 35 since there are 70 bins for -35,35. makes -35 = 0. 
   if(vzbin<0 || vzbin >70){
     cout << "VZ OUT OF RANGE" << endl;
     return;
   }
   //cout << "Vz: " << vzbin << " size: " << hadVec[vzbin].size() << endl;
   if(hadVec[vzbin].size() < maxBufferSize)
     hadVec[vzbin].push_back(trk); // Stores the event itself, not the pointer
   else
     {
       TRandom3* gRand = new TRandom3();
       Int_t eventPoint = (int) gRand->Uniform(0,maxBufferSize-1e-6);
       hadVec[vzbin][eventPoint] = trk; // Stores the event itself, not the pointer
       delete gRand;
     }
   //cout << "Vz(after add): " << vzbin << " size: " << hadVec[vzbin].size() << endl;
 }

 void StNpeRead::computeMixedEvents(StDmesonTrack* trk, Double_t vz)
 {
   //DEBUG  cout << "in compute" << endl;
   Float_t Phi = trk->gMom().phi();
   Float_t pT  = trk->gMom().perp();
   Float_t Eta = trk->gMom().pseudoRapidity();
   Int_t vzbin = (Int_t)vz+35;
   if(vzbin<0 || vzbin>70)
     return;
   
   if(hadVec[vzbin].size()<=0)
     return;
   
   cout << "at hadVec for" << endl;
   for(Int_t it=0; it < hadVec[vzbin].size(); it++)
     {
       StDmesonTrack* htrk = hadVec[vzbin][it]; 
       Float_t hpT   = htrk->pMom().perp();
       
       if(trk != htrk)
	 {
	   cout << "actually have mixed had track" << endl;
	   Float_t hPhi = htrk->pMom().phi();
	   Float_t hpT  = htrk->pMom().perp();
	   Float_t hEta = htrk->pMom().pseudoRapidity();
	   
	   Float_t dPhi = Phi-hPhi;
	   if(dPhi > (3.*pi)/2.) dPhi = dPhi-2*pi;
	   if(dPhi < -1*pi/2.)   dPhi = dPhi+2*pi;
	   Float_t dEta = Eta - hEta;
	   mh3MixedDelPhi -> Fill(dPhi, pT, hpT);
	   mh3MixedDelEta -> Fill(dEta, pT, hpT);
	   mh3MixedEtaPhi -> Fill(dPhi, dEta, pT);
	 }
     }
 }
