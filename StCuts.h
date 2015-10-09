
#include "StTRIGGERS.h"

#ifndef CUTS_H
#define CUTS_H

namespace badruns
{
  // marked as Junk but produced in picoDst// 13057052 13059087 VPDMB vz have strange point  13039122,13039126,13039131,13052068 not fire the trigger
  const int n_bad = 3;
  int bad_arr[n_bad] = {13069009,13054046,13049075};//,13039132,13039135,13039136,13052068,13057053,13060001};

// the following runs have no monitor bbcmbtof0
const int n_bad_ht1bbcmbtof0 = 13;
int bad_ht1bbcmbtof0_arr[n_bad_ht1bbcmbtof0] ={13067005,13067006,13067007,13067008,13067009,13067010,13067011,13067012,13067013,13067014,13067015,13067016,13067017};

const int n_bad_ht0bbcmbtof0 = 25;
int bad_ht0bbcmbtof0_arr[n_bad_ht0bbcmbtof0]={
  13044118,
  13047014,
  13047046,
  13049052,
  13052063,
  13053021,
  13054020,
  13054057,
  13058047,
  13058048,
  13061024,
  13061025,
  13061026,
  13064030,
  13066101,
  13066102,
  13066104,
  13067001,
  13069068,
  13070019,
  13070030,
  13071004,
  13071006,
  13071034,
  13071037};
const int n_bad_ht2bbcmb = 32;
int bad_ht2bbcmb_arr[n_bad_ht2bbcmb]={13041002,
  13042003,
  13044118,
  13047004,
  13047014,
  13047044,
  13047046,
  13047050,
  13049052,
  13053021,
  13054020,
  13054057,
  13058047,
  13058048,
  13061024,
  13061026,
  13064030,
  13066101,
  13066102,
  13066104,
  13067001,
  13070056,
  13070057,
  13070058,
  13070059,
  13070060,
  13070061,
  13071003,
  13071004,
  13071006,
  13071008,
  13071034};
const int n_bad_vpdmb = 17;
int bad_vpdmb_arr[n_bad_vpdmb]={13044118,
  13047014,
  13048031,
  13049052,
  13052063,
  13054057,
  13057053,
  13058047,
  13058048,
  13060001,
  13061024,
  13061026,
  13066101,
  13066102,
  13066104,
  13066109,
  13067001};
bool is_bad_run(int in,int trg) 
{
  int* arr = 0;
  int size = 0;

  if(trg == 1) {arr = bad_ht1bbcmbtof0_arr; size = n_bad_ht1bbcmbtof0;}

  else if(trg ==0) {arr = bad_ht0bbcmbtof0_arr; size = n_bad_ht0bbcmbtof0;}
  else if(trg ==2) {arr = bad_ht2bbcmb_arr; size = n_bad_ht2bbcmb;}
  else  if(trg ==3) {arr =bad_vpdmb_arr; size = n_bad_vpdmb;} 


  // {arr = bad_arr; size=n_bad;}

  //for(unsigned int i = 0; i<sizeof(arr)/sizeof(int); i++) // you cannot get the size of the array using sizeof(POINTER-TO-ARRAY)
  for(int i = 0; i<size; i++)
  {  
    if(in == arr[i]) { return true;}
  }   

  return false;
}

}

namespace cuts
{
  //event
  //	const Float_t vz = 100.0;
  const Float_t vzVpdVz = 6.0;

  //pions 
  const Float_t nSigmaPion = 2.0;
  const Float_t pionPt = 0.2;
  const Float_t pionEta = 1.0;

  //soft 
  const Float_t softNSigmaPion = 2.0;
  const Float_t softPionPt = 0.15;

  //TOF kaons 
  const Float_t kaonPt = 0.2;
  const Float_t tofYLocal = 1.8;

  // TPC kaons
  const Float_t nSigmaKaon = 2.0;
  const Float_t kaonEta = 1.0;



  //NPE cuts
  //event
  const Float_t vz = 35.;
  const Float_t diff_vzVpdVz = 6.0;
  // const Float_t vr_cut=1.0;
  //        const Float_t field_cut=0.1;
  const Float_t RefMult_cut=8.0;
  //tracking
  const Int_t nHitsFit = 20;
  const Float_t nHitsFitnHitsMax = 0.52;
  const Float_t nHitsDedx=15;

  const Float_t eta = 1.0;
  const Float_t eta_low =0.7;
  const Float_t gDca = 1.5;
  const Float_t pt=0.2;
  const Float_t Pt_cut=2.0;

  const Int_t nsigmae_low=-1;
  const Int_t nsigmae_high=3;
  const Int_t NEta =1;
  const Int_t NPhi=1;
  const Int_t nhitDedx=15;
  const Int_t nhitsFit=20;
  const Int_t nhit_Possmax=0.52;     
  //	const Float_t pt_cut =2.;
  const Float_t Dz_cut=3;

  const Float_t Dphi_cut =0.015;
  const Float_t poe_low=0.3;
  const Float_t poe_high=1.5;
  const Float_t first_point_cut =73;
  // tof cut
  const Float_t tofbeta=0.03;
  const Float_t tof_Ylocal=1.8;
  //for Pair;

  const Float_t massDCA=0.24;
  const Float_t pairDCA=1.0;

  // For Had ID (Bingchu cuts)
  const Float_t hadPtMin = 0.2;
  const Float_t hadEta = 1.0;
  const Float_t hadDedx = 10;
  const Float_t hadHitsFit = 15;
  const Float_t hadDCA = 1.0;

  // pion selection
  const Float_t nsigmap_low=-0.1;
  const Float_t nsigmap_high=0.1;

  //dEta Cut
  const Float_t dEtaCut = 0.5;
}
#endif
