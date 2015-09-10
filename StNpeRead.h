
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
  void zFill_Pileup( Int_t, StDmesonEvent*, Double_t );
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
       
       /* TH2F *	mh2Vz_VPD; */
       /* TH2F *	mh2VPDVz_VPD; */
       /* TH2F *	mh2Vz_VPDVPD; */
       
       /*	TH2F *	mh2TPCVzRunID[4];
	TH2F *	mh2VPDVzRunID[4];
	TH2F *	mh2TPCVz_VPDVz[4];
	TH2F *	mh2ZDCRunID[4];
	TH2F *	mh2BBCRunID[4];
	TH2F *	mh2MultiPosRunID[4];
	TH2F *	mh2MultiNegRunID[4];
	TH2F*	mh2MultiRunID[4];
	TH2F *  mh2PtRunID[4];
	TH2F *	mh2EtaRunID[4];
	TH2F *	mh2PhiRunID[4];
	TH2F *	mh2BetaRunID[4];
	TH2F *	mh2DedxRunID[4];
	TH2F *	mh2Beta_RunID[4];
	TH2F *	mh2NhitFitRunID[4];
	TH2F *	mh2NhitDedxRunID[4];
	TH2F *	mh2nsigmaERunID[4];
	TH2F *	mh2gDCARunID[4];
	TH2F *	mh2PhiVsEta_Range1[4];
	TH2F *	mh2PhiVsEta_Range2[4];
	TH2F *	mh2PhiVsEta_Range3[4];
	TH1F *  mh1PassTofMatchRunID[4];
	TH1F *  mh1WPassTofMatchRunID[4];
	TH2F *	mh2VXY[4];
	// profile

	TProfile*	mh2TPCVzRunID_Profile[4];
	TProfile*	mh2VPDVzRunID_Profile[4];

	TProfile*	mh2ZDCRunID_Profile[4];
	TProfile*	mh2BBCRunID_Profile[4];
	TProfile*	mh2MultiPosRunID_Profile[4];
	TProfile*	mh2MultiNegRunID_Profile[4];
	TProfile *	mh2MultiRunID_Profile[4];
	TProfile*        mh2PtRunID_Profile[4];
	TProfile*	mh2EtaRunID_Profile[4];
	TProfile*	mh2PhiRunID_Profile[4];
	TProfile*	mh2BetaRunID_Profile[4];
	TProfile*	mh2DedxRunID_Profile[4];
	TProfile*	mh2Beta_RunID_Profile[4];
	TProfile*	mh2NhitFitRunID_Profile[4];
	TProfile*	mh2NhitDedxRunID_Profile[4];
	TProfile*	mh2nsigmaERunID_Profile[4];
	TProfile*	mh2gDCARunID_Profile[4];*/

	//
 
       /*	TH2F *	mh2DSMADC_Inclusive[4];
	TH2F*	mh2DSMADC_PhotoUnlike[4];
	TH2F *	mh2DSMADC_Photolike[4];

	TH2F *	mh2ADC0_Inclusive[4];
	TH2F *	mh2ADC0_PhotoUnlike[4];
	TH2F *	mh2ADC0_Photolike[4];

	TH1F*	mh1TowerID_Inclusive[4];
	TH1F*	mh1TowerID_cut_Inclusive[4];
	TH1F*	mh1TowerID_all_Inclusive[4];

	TH1F *	mh1TowerID_Inclusive_noHOt[4];
	TH1F*	mh1TowerID_cut_Inclusive_noHOt[4];
	TH1F*	mh1TowerID_all_Inclusive_noHOt[4];

	TH2F *  mh2ADC0_DSMadcInclusive[4];
        TH2F *  mh2ADC0_DSMadcPhotoUnlike[4];
        TH2F *  mh2ADC0_DSMadcPhotolike[4];

        TH2F* mh2InvMassPoeUnlike[4];
	TH2F* mh2PoeUnlike[4];
	TH2F* mh2InvMassPoelike[4];  
	TH2F* mh2Poelike[4];

	TH2F* mh2InvMassDzUnlike[4];
	TH2F* mh2DzUnlike[4];
	TH2F* mh2InvMassDzlike[4];  
	TH2F* mh2Dzlike[4];
    
	TH2F* mh2InvMassDpUnlike[4];
	TH2F* mh2DpUnlike[4];
	TH2F* mh2InvMassDplike[4];  
	TH2F* mh2Dplike[4];

	TH2F*	mh2InvMassEMCUnlike[4];
	TH2F* mh2InvMassEMClike[4];    
	TH2F*	mh2InvMassNEMCUnlike[4];
	TH2F* mh2InvMassNEMClike[4];    
	TH2F*	mh2InvMassNEEMCUnlike[4];
	TH2F* mh2InvMassNEEMClike[4];    
	TH2F*	mh2InvMassNPEMCUnlike[4];
	TH2F* mh2InvMassNPEMClike[4];    

	TH2F*	mh2NEtaUnlike[4];
	TH2F*	mh2NEtalike[4];
	TH2F*	mh2NPhiUnlike[4];
	TH2F*	mh2NPhilike[4];

	TH3F* mh3nSigPart_EMCUnlike[4];
	TH3F* mh3nSigPart_EMClike[4];

	TH3F* mh3nSigPart_TREMCUnlike[4];
	TH3F* mh3nSigPart_TREMClike[4];
	TH3F* mh3nSigPart_ADCTREMCUnlike[4];
	TH3F* mh3nSigPart_ADCTREMClike[4];

	TH3F* mh3nSigEUnlike[4];
	TH3F* mh3nSigPartUnlike[4];
	TH3F* mh3nSigElike[4];  
	TH3F* mh3nSigPartlike[4];


	//all cuts apllied
	TH2F*     mh2DsmADCUnlike[4];
	TH2F*     mh2DsmADClike[4];
	TH2F*  mh2Btowadc0[4];
	TH3F*	mh3EMC_PartUnlike[4];
	TH3F*	mh3EMC_Partlike[4];
	TH3F*	mh3EMC_ADCPartUnlike[4];
	TH3F*	mh3EMC_ADCPartlike[4];


	TH3F* mh3EtaPhiUnlike[4];
	TH3F* mh3EtaPhilike[4];	
	TH2F * mh2InvMasslike[4];
	TH2F * mh2InvMassUnlike[4];
	TH2F* mh2nSigmaElec[4];
	TH1F * mh1electronPt[4];

	//data QA%
	TH2F*	mPhi_ptUnlike[4];
	TH2F*	mEta_ptUnlike[4];
	TH2F*	mPhi_ptlike[4];
	TH2F*	mEta_ptlike[4];
	
	TH2F*	mHitFit_ptUnlike[4];
	TH2F*	mNSMDEta_ptUnlike[4];
	TH2F*	mHitsDedxUnlike[4];
	TH2F*	mHitFit_ptlike[4];
	TH2F* 	mNSMDEta_ptlike[4];
	TH2F* 	mHitsDedxlike[4];
	
	TH1F*	mNTrackUnlike[4];
	TH1F* 	mNTrack_cutUnlike[4];
	TH1F* 	mNTracklike[4];
	TH1F* 	mNTrack_cutlike[4];
	TH2F* 	mFitPos_ptlike[4];
	TH2F* 	mgDcalike[4];
	TH2F* 	mFitPos_ptUnlike[4];
	TH2F* 	mgDcaUnlike[4];
	
	TH2F* mNSMDPhi_ptUnlike[4];
	TH2F* mNSMDPhi_ptlike[4];
	TH1F* mNTrack_cut25Unlike[4];
	TH1F* mNTrack_cut25like[4];

	TH1F* mNTrack_cut20Unlike[4];
        TH1F* mNTrack_cut20like[4];
	
	TH2F* mNsigElike[4];
	TH2F* mNsigEUnlike[4];
	
	TH2F* mDedxlike[4];
	TH2F* mDedxUnlike[4];
	TH2F *mPoelike[4];
	TH2F *mPoeUnlike[4]; 	
	TH1F *mh1Vz[4];
	TH1F *mh1VzPS[4];
	TH2F* mh2Vxy[4]; 
	TH1F *mh1Vz_BBCMB;
	TH1F *mh1VPDVz_BBCMB;
	TH1F *mh1Vz_BBCMBTOF0;
	TH1F *mh1VPDVz_BBCMBTOF0;
	TH1F *mh1Vz_BBCMBPS;
	TH1F *mh1Vz_BBCMBTOF0PS;
	TH2F * HT0_HT2;	
	//-----------------------------low pt-----------------------------------------
	TH2F *	mh2nSigmaElec_VPD;
	TH2F *  mh2_PionnSigmaElecDiff_VPD;
	TH2F *  mh2_MergrePionnSigmaElecDiff_VPD;
	TH3F *	mh3_nSigmaElec_Beta_VPD;


	TH2F * mh2_nSigmaElec_BetaP_VPD;
	TH2F *  mh2_nSigmaElec_BetaPt_VPD;
	TH2F *  mh2_KaonnSigmaElecDiff_VPD;
	TH2F *  mh2_ProtonnSigmaElecDiff_VPD;
	TH2F *    mh2_nsigamE_pt_VPD;
	TH2F * mh2_nsigamE_p_VPD;

	TH1F *	mh1electronPt_VPD;
	TH1F *  mh1electronPt_PSVPD;
	TH1F *	mh1electronPt_MBPSVPD;
	TH2F *  mh2electronPtEta_VPD;
	TH2F *  mh2_InvMass_VPD;
	TH2F *	mh2_Pion_nSigmaElec_VPD;
	TH2F *	mh2_Kaon_nSigmaElec_VPD;
	TH2F *	mh2_Proton_nSigmaElec_VPD;

	TH2F *  mh2_KaonnSigmaEle_eta_VPD;
	TH2F *  mh2_ProtonnSigmaEle_eta_VPD;
	TH2F *  mh2_PionnSigmaEle_eta_VPD;

	TH2F * mh2_InvMass_KaonUnlike_VPD;
	TH2F * mh2_NsigmaE_KaonUnlike_VPD;
	TH2F * mh2_InvMass_Kaonlike_VPD;
	TH2F * mh2_NsigmaE_Kaonlike_VPD;
	
	TH2F * 	mh2InvTofBetaUnlike_VPD ;
	TH2F * 	mh2InvTofBetalike_VPD ;
	
	TH3F * 	mh3TofBetaUnlike_VPD ;
	TH3F * 	mh3TofBetalike_VPD;
	TH3F * mh3TofBeta_PartnerUnlike_VPD;
	TH3F *	mh3TofBeta_Partnerlike_VPD;

	TH3F *	mh3TofBeta_CutPartnerUnlike_VPD;
	TH3F *	mh3TofBeta_CutPartnerlike_VPD;
	// no  yloacal                                                                                                                   
	TH2F * 	mh2InvTofYlocalUnlike_VPD ;
	TH2F * 	mh2InvTofYlocallike_VPD ;
	
	TH2F * 	mh2TofYlocalUnlike_VPD ;
	TH2F * 	mh2TofYlocallike_VPD ;
	TH2F *	mTofBetaUnlike_match_VPD;
	TH2F *  mTofBetalike_match_VPD;
	// all the cuts applied                                                                                                          
	TH2F * 	mh2InvUnlike_VPD ;
	TH2F * 	mh2Invlike_VPD ;

	TH2F *  mh2InvUnlike_PSVPD ;
        TH2F *  mh2Invlike_PSVPD ;

	TH2F *	mh2InvPartnerUnlike_VPD;
	TH2F *  mh2InvPartnerlike_VPD;
	TH2F *	mh2CutPartnerUnlike_VPD;
	TH2F *	mh2CutPartnerlike_VPD;
	TH2F * 	mh2nSigePartnerUnlike_VPD ;
	TH2F * 	mh2nSigePartnerlike_VPD ;
	TH3F* mh3nSigmaEUnlike_VPD;
	TH3F* mh3nSigmaElike_VPD;
	//	QA
	TH2F *   mPhi_ptUnlike_VPD;
	TH2F * mPhi_ptlike_VPD;

	TH2F *   mEta_ptUnlike_VPD;
	TH2F *  mEta_ptlike_VPD;
	
	TH2F *  mTofBetaUnlike_VPD;
	TH2F *  mTofBetalike_VPD;
	
	TH2F *   mTofYlocalUnlike_VPD;
	TH2F *  mTofYlocallike_VPD;
	
	TH2F *  mHitFit_ptUnlike_VPD;
	TH2F *  mHitFit_ptlike_VPD;
	
	TH2F *  mHitsDedxUnlike_VPD;
	TH2F * mHitsDedxlike_VPD;
	
	TH2F *  mgDcalike_VPD;
	TH2F *  mgDcaUnlike_VPD;
	
	TH2F *  mNsigElike_VPD;
	TH2F *  mNsigEUnlike_VPD;
	
	TH2F *  mDedxlike_VPD;
	TH2F * mDedxUnlike_VPD;

	TH1F*	mNTrack_cutUnlike_VPD;
	TH1F*	mNTrack_cutlike_VPD;
	
	TH2F *	mFitPos_ptUnlike_VPD;
	TH2F *	mFitPos_ptlike_VPD;
       */
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
	TH3F *  mh3nTracksZdcx[4];
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
	TH1F *  mh1PtHadTracks[4];

	TH3F *  mh3MixedDelPhi;
	TH3F *  mh3MixedDelEta;
	TH3F *  mh3MixedEtaPhi;


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
