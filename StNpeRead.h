
#ifndef StNpeRead_h
#define StNpeRead_h

//#include "StNpeMaker/StDmesonEvent.h"
#include "StDmesonMaker/StDmesonEvent.h"

#include <exception>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <vector>
#include "TProfile.h"

class TString;
class TF1;
class TH1D;
class TH1F;
class TH2F;
class TH3F;
class TFile;
class TNtuple;
class StDmesonTrack;
class StDmesonEvent;
class prescales;
//class StDmesonEvent;
class StElectronPair;
// TH1F::SetDefaultSumw2();
class StNpeRead
{
 public:
  //  prescales *mPrescales;
  StNpeRead(const char *outName);
  virtual ~StNpeRead();
  // TH1F::SetDefaultSumw2();    
  void bookObjects();
  void read(TString fileName);
  void writeObjects();
  //    Bool_t isGoodEvent(StDmesonEvent* evt,Bool_t& bVPDMB,Bool_t& bHT1_BBCMB_TOF0,Bool_t& bHT2_BBCMB);
  
  
  Bool_t isTrgTrack(StDmesonTrack*,UShort_t htTrg);
  void SetRunID_VPD( );  
  Bool_t isBadRun(StDmesonEvent *);
  Bool_t isGoodEvent(StDmesonEvent*,Int_t  );
  Bool_t isHotTower(StDmesonTrack *,Int_t );  
  Bool_t isElectron(StDmesonTrack*);
  //Fill histogram
  void Fill_event_hist(StDmesonEvent *,Int_t,Double_t);
  
  void Fill_pair_hist_high_pt(Int_t,StElectronPair*,StDmesonTrack*,StDmesonTrack*,Double_t);
  void Fill_pair_hist_low_pt(Int_t,StElectronPair*,StDmesonTrack*,StDmesonTrack*,Double_t);
  void Fill_Inclusive_hist (Int_t ,StDmesonEvent * ,Double_t );
  void Fill_PhotonicE_hist (Int_t ,StDmesonEvent * ,Double_t );

  //// Added Z. Miller
  void zFillHists( Int_t, StDmesonEvent *, Double_t ); 
  void zFill_Inclusive( Int_t, StDmesonEvent*, Double_t );
  void zFill_Photonic( Int_t, StDmesonEvent*, Double_t );
  //void zFillProjections( ); // MOVED TO OFFLINE RECONSTRUCTION
  //// end add

  void Do_run_QA(StDmesonEvent * );   
  void Fill_RunQA(Int_t,StDmesonEvent * );
  void Fill_Kaon_Kaon(StDmesonEvent * ,Int_t);
  //  void Fill_Inclusive_hist_low_pt(Int_t ,StDmesonEvent * ,Double_t );
  // void Fill_PhotonicE_hist_low_pt(Int_t ,StDmesonEvent * ,Double_t ); 
//for eID
    Bool_t pass_cut_GoodTrack(StDmesonTrack *);
    Bool_t pass_loose_track_qaulity( StDmesonTrack *, Int_t );
    Bool_t pass_cut_Pt_Eta(StDmesonTrack *);
    Bool_t pass_cut_nsigmaE(StDmesonTrack* );
    Bool_t pass_cut_poe(StDmesonTrack *);
    Bool_t pass_cut_EMC(StDmesonTrack *);      
    Bool_t pass_cut_Match_EMC_Dphi(StDmesonTrack * );
    Bool_t pass_cut_Match_EMC_Dz(StDmesonTrack * );
    Bool_t pass_cut_Tof_Match( StDmesonTrack *);
    Bool_t pass_cut_Tof( StDmesonTrack *);
    //// Added Z. Miller
    Bool_t pass_cut_ADC( Int_t, StDmesonTrack *);
    Bool_t is_EMC_Track(StDmesonTrack *); // Just check nPhi, nEta
    Bool_t pass_cut_hTrack(StDmesonTrack *); // check for hadrons tracks
    Double_t getHadronWt(Double_t, Double_t);
    Int_t readEff(); // Read in Efficiencies from exterior files
    void addToHadBuffer(StDmesonTrack *, Double_t); // Create a buffer for hadrons
    void computeMixedEvents(StDmesonTrack *, Double_t);
    Bool_t pass_cut_nsigmaPi(StDmesonTrack *);
    Float_t correct_dPhi(Float_t); // controls wrap around, so I only have to change values in one place
    Double_t getTrgEff(Int_t, Double_t); // Get electron trigger efficiency from XB embedding
    //// end Add

  private:
    TFile* mOutputFile;
    StDmesonEvent* mNpeEvent;
    prescales *mPrescales;
    ifstream file_runNumber;
    ofstream  outfile;//("runID.txt",ios::out|ios::app);
       map <int,int> runID_List;
       set  <Int_t> runIDPico;
       TH3F * ADC_nocut;
       TH3F * ADC_cut;
       TH1F *	mh1Vz_VPDVz;
       TH1F *	mh1Vz_VPD;
       TH1F *  mh1VzVPD_VPD;
       TH1F *	mh1Vz_VPDPS;
       
	// Added Z. Miller
	//TH1F *  testHist1D[5];
	//TH2F *  testHist2D[5];
	TH1F *  mh1PtAllTracks[4];
	TH1F *  mh1PtETracks[4];
	TH2F *  mh2nSigmaEPt[4];
	TH2F *  mh2nSigmaEPt_eID[4];
	TH2F *  mh2PoePt[4];
	TH2F *  mh2nPhiPt[4];
	TH2F *  mh2nEtaPt[4];
	TH2F *  mh2PhiDistPt[4];
	TH2F *  mh2ZDistPt[4];
	TH2F *  mh2PhiQPt[4];
	TH2F *  mh2TofPtAll[4];
	TH2F *  mh2TofPtE[4];
	TH2F *  mh2InvMassPtAll[4];
	TH2F *  mh2InvMassPtE[4];
	TH2F *  mh2InvMassPtUS[4];
	TH2F *  mh2InvMassPtLS[4];
	TH2F *  mh2nSigmaPionPt[4];
	TH3F *  mh3nTracksZdcx[4][4];
	TH3F *  mh3nTracksZdcxUS[4][4];
	TH3F *  mh3nTracksZdcxLS[4][4];
	TH3F *  mh3nTracksZdcxHad[4][4];
        TH2F *  mh2PtEZdcxHad[4];
	TH2F *  mh2PtEZdcx[4];
	TH2F *  mh2PtEZdcxUS[4];
	TH2F *  mh2PtEZdcxLS[4];
	TH3F *  mh3DelPhiIncl[4];
	TH3F *  mh3DelPhiPhotLS[4];
	TH3F *  mh3DelPhiPhotUS[4];
	TH3F *  mh3DelPhiPhotUSNP[4];
	TH3F *  mh3DelPhiPhotLSNP[4];
	TH3F *  mh3DelPhiPhotInclNP[4];
	TH3F *  mh3DelPhiInclWt[4];
	TH3F *  mh3DelPhiPhotLSWt[4];
	TH3F *  mh3DelPhiPhotUSWt[4];
	TH3F *  mh3DelPhiHadHad[4];
	TH3F *  mh3DelPhiHadHadWt[4];
	TH1F *  mh1PtHadTracks[4];

	TH3F *  mh3MixedDelPhi;
	TH3F *  mh3MixedDelEta;
	TH3F *  mh3MixedEtaPhi;
	TH3F *  mh3MixedDelPhiWt;
        TH3F *  mh3MixedDelEtaWt;
        TH3F *  mh3MixedEtaPhiWt;

	// For comparison to 2.5-3.5 GeV bin prev analysis
	TH1F *  mh1delPhiIncl;
	TH1F *  mh1delPhiUS;
	TH1F *  mh1delPhiLS;
	TH1F *  mh1delPhiHad;
	TH1F *  mh1TrigCount;

	//Projections[ptbin][trig]
	/*TH1D *  projHPhi[14][4];
	TH1D *  projnSigmaE[14][4];
	TH1D *  projDelPhiIncl[14][4];
	TH1D *  projnSigmaE_eID[14][4];
	TH1D *  projDelPhiPhotLS[14][4];
	TH1D *  projDelPhiPhotUS[14][4];
	TH1D *  projInvMassUS[14][4];
	TH1D *  projInvMassLS[14][4];*/
	
	Bool_t writeRunQA; // flags for write control
	Bool_t writeDataQA;
	Bool_t writeXiaozhiHists;

	// Hadron weighting
	Double_t effPars[20][3];
	TF1 *fEff;

	Float_t pi; // just to make it a single definition
	Int_t numPtBins;
	Int_t maxBufferSize;
	Bool_t isAddedToBuffer;
	std::vector<Float_t> hadPhi[70];
	std::vector<Float_t> hadEta[70];
	std::vector<Float_t> hadPt[70];
	ClassDef(StNpeRead, 1)
	  };

#endif
