#ifndef softDropGroomer_h
#define softDropGroomer_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include "TF1.h"

#include "../PU14/PU14.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "TMath.h"
#include "TH2.h"


#include "jetCollection.hh"
#include "jewelMatcher.hh"
//---------------------------------------------------------------
// Description
// This class runs SoftDrop grooming on a set of jets
// Author: M. Verweij
//---------------------------------------------------------------

class softDropGroomer {

private :
  double zcut_;
  double beta_;
  double r0_;
  int fReclusterAlgo;
  Float_t fqhat;                        //qhat
  Float_t fxlength;                     //medium length
  Bool_t bflagAdded;
  Int_t  fAdditionalTracks;
  Int_t kRecursiveAlgo;
  double eventWeight;


  

  std::vector<fastjet::PseudoJet> fjInputs_;   //ungroomed jets
  std::vector<fastjet::PseudoJet> fjOutputs_;  //groomed jets
  std::vector<double>             zg_;         //zg of groomed jets
  std::vector<int>                drBranches_; //dropped branches
  std::vector<double>             dr12_;       //distance between the two subjets
  std::vector<double>             sjmass_;       //subleading subjet mass
  std::vector<double>             logdr12_;     
  std::vector<double>             logzgdr12_;
  std::vector<double>             recur_jetpt_;
  std::vector<double>             recur_logdr12_;     
  std::vector<double>             recur_logzgdr12_;
  std::vector<int>                recur_n_;
  std::vector<double>             sjleadingtrack_;
  std::vector<double>             injectedtrack_z_;
  std::vector<double>             injectedtrack_theta_;
  std::vector<double>             injectedtrack_pt_;
  std::vector<double>             jetProfile0_;
  std::vector<double>             jetProfile1_;
  std::vector<double>             jetProfile2_;
  std::vector<double>             jetProfile3_;
  std::vector<double>             jetProfile4_;
  std::vector<double>             jetProfile5_;
  std::vector<double>             jetProfile6_;
  std::vector<double>             jetProfile7_;

  std::vector<std::vector<fastjet::PseudoJet>> constituents_;
  std::vector<std::vector<fastjet::PseudoJet>> constituents1_;
  std::vector<std::vector<fastjet::PseudoJet>> constituents2_;

public :
  softDropGroomer(double zcut = 0.1, double beta = 0., double r0 = 0.4);
  void setZcut(double c);
  void setBeta(double b);
  void setR0(double r);
  void setInputJets(std::vector<fastjet::PseudoJet> v);
  void setReclusteringAlgo(int r);
  void setMediumParameters(Float_t t, Float_t c);
  void setInjectTracks(Bool_t b, Int_t n);
  void setRecursiveAlgo(Int_t i);
  void setEventWeight(double e);
  void RecursiveParents(fastjet::PseudoJet fJet,Bool_t bflagAdded, Int_t fAdditionalTracks, Int_t ReclusterAlgo);
  std::vector<fastjet::PseudoJet> getGroomedJets() const;
  std::vector<double> getZgs() const;
  std::vector<int> getNDroppedSubjets() const;
  std::vector<double> getDR12() const;
  std::vector<double> getSubJetMass() const;
  std::vector<double> getLogDR12() const;
  std::vector<double> getLogZgDR12() const;
  std::vector<double> getRecur_JetPt() const;
  std::vector<double> getRecur_LogDR12() const;
  std::vector<double> getRecur_LogZgDR12() const;
  std::vector<double> getInjectedz() const;
  std::vector<double> getInjectedtheta() const;
  std::vector<double> getInjectedpt() const;
  std::vector<double> getJetProfile(Int_t i) const;

  std::vector<int> getRecur_N() const;
  std::vector<double> getSubJetLeadingTrackPt() const;

  std::vector<fastjet::PseudoJet> doGrooming(jetCollection &c);
  std::vector<fastjet::PseudoJet> doGrooming(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> doGrooming();

  TH2F* calculateProfile(jetCollection &c);
  TH2F* calculateProfile(std::vector<fastjet::PseudoJet> v);
  TH2F* calculateProfile();
  
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents() {return constituents_;}
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents1() {return constituents1_;}
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents2() {return constituents2_;}

  std::vector<fastjet::PseudoJet> doGroomingWithJewelSub(jetCollection &c, std::vector<fastjet::PseudoJet> particlesDummy);
  std::vector<fastjet::PseudoJet> doGroomingWithJewelSub(std::vector<fastjet::PseudoJet> v, std::vector<fastjet::PseudoJet> particlesDummy);
  std::vector<fastjet::PseudoJet> doGroomingWithJewelSub(std::vector<fastjet::PseudoJet> particlesDummy);
  double RelativePhi(Double_t mphi,Double_t vphi);

  Double_t leading_track_pt(fastjet::PseudoJet jet);
};

softDropGroomer::softDropGroomer(double zcut, double beta, double r0)
   : zcut_(zcut), beta_(beta), r0_(r0)
{
}

void softDropGroomer::setZcut(double c)
{
   zcut_ = c;
}

void softDropGroomer::setBeta(double b)
{
   beta_ = b;
}

void softDropGroomer::setR0(double r)
{
   r0_ = r;
}

void softDropGroomer::setInputJets(const std::vector<fastjet::PseudoJet> v)
{
   fjInputs_ = v;
}

void softDropGroomer::setReclusteringAlgo(int r)
{
  fReclusterAlgo = r;
}

void softDropGroomer::setMediumParameters(Float_t t, Float_t c)
{
  fqhat=t; fxlength=c;
}

void softDropGroomer::setInjectTracks(Bool_t b, Int_t n)
{
  bflagAdded = b;
  fAdditionalTracks = n;
}

void softDropGroomer::setRecursiveAlgo(Int_t i)
{
  kRecursiveAlgo = i;
}

void softDropGroomer::setEventWeight(double e)
{
  eventWeight = e;
}


std::vector<fastjet::PseudoJet> softDropGroomer::getGroomedJets() const
{
   return fjOutputs_;
}

std::vector<double> softDropGroomer::getZgs() const
{
   return zg_;
}

std::vector<int> softDropGroomer::getNDroppedSubjets() const
{
   return drBranches_;
}

std::vector<double> softDropGroomer::getDR12() const
{
   return dr12_;
}

std::vector<double> softDropGroomer::getSubJetMass() const
{
   return sjmass_;
}

std::vector<double> softDropGroomer::getSubJetLeadingTrackPt() const
{
   return sjleadingtrack_;
}

std::vector<double> softDropGroomer::getLogDR12() const
{
   return logdr12_;
}

std::vector<double> softDropGroomer::getLogZgDR12() const
{
   return logzgdr12_;
}

std::vector<double> softDropGroomer::getRecur_JetPt() const
{
   return recur_jetpt_;
}

std::vector<double> softDropGroomer::getRecur_LogDR12() const
{
   return recur_logdr12_;
}

std::vector<double> softDropGroomer::getRecur_LogZgDR12() const
{
   return recur_logzgdr12_;
}

std::vector<int> softDropGroomer::getRecur_N() const
{
   return recur_n_;
}

std::vector<double> softDropGroomer::getInjectedz() const
{
   return injectedtrack_z_;
}

std::vector<double> softDropGroomer::getInjectedtheta() const
{
   return injectedtrack_theta_;
}

std::vector<double> softDropGroomer::getInjectedpt() const
{
   return injectedtrack_pt_;
}

std::vector<double> softDropGroomer::getJetProfile(Int_t i) const
{
  if(i==0) return jetProfile0_;
  if(i==1) return jetProfile1_;
  if(i==2) return jetProfile2_;
  if(i==3) return jetProfile3_;
  if(i==4) return jetProfile4_;
  if(i==5) return jetProfile5_;
  if(i==6) return jetProfile6_;
  if(i==7) return jetProfile7_;

}

std::vector<fastjet::PseudoJet> softDropGroomer::doGrooming(jetCollection &c)
{
   return doGrooming(c.getJet());
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGrooming(std::vector<fastjet::PseudoJet> v)
{
   setInputJets(v);
   return doGrooming();
}

TH2F* softDropGroomer::calculateProfile(jetCollection &c)
{
  return calculateProfile(c.getJet());
}

TH2F* softDropGroomer::calculateProfile(std::vector<fastjet::PseudoJet> v)
{
   setInputJets(v);
   return calculateProfile();
}

TH2F* softDropGroomer::calculateProfile()
{
  TH2F *h2PtJetRho = new TH2F("h2PtJetRhotemp","",150,0,150,20,0,0.5);
  int nJets = fjInputs_.size();
  double min_track_pt = 1;
  double max_track_r  = 0.5; 
  for(int i = 0;i<nJets;i++) {
      
    std::vector<fastjet::PseudoJet> tracks_all = fjInputs_[i].constituents();
     std::vector<fastjet::PseudoJet> tracks, tracks_neutral;
    SelectorIsCharged().sift(tracks_all, tracks, tracks_neutral);
    Double_t jetPt = fjInputs_[i].perp();
    Double_t jetEta = fjInputs_[i].pseudorapidity();
    Double_t jetPhi = fjInputs_[i].phi();
    Double_t trackPt,trackEta,trackPhi,trackDelR,delEta,delPhi;
    for(int j=0;j<tracks.size();j++){
      trackEta = tracks[j].pseudorapidity();
      trackPhi = tracks[j].phi();
      trackPt = tracks[j].perp();
      delEta = trackEta-jetEta;
      delPhi = RelativePhi(trackPhi,jetPhi);
      trackDelR = TMath::Sqrt(delEta*delEta+delPhi*delPhi);
      if(trackPt > min_track_pt && trackDelR < max_track_r) {
	h2PtJetRho->Fill(jetPt,trackDelR,eventWeight*(trackPt/jetPt));
      }
    }

 }
  return h2PtJetRho;
}

Double_t softDropGroomer::leading_track_pt(fastjet::PseudoJet jet)
{
  std::vector< fastjet::PseudoJet > constituents = jet.constituents();
  fastjet::PseudoJet leadingtrack = constituents[0];
  Double_t track_pt;
  for(size_t i=0; i<constituents.size(); i++){
    track_pt=constituents[i].perp();
    if (track_pt > leadingtrack.perp()){
      leadingtrack = constituents[i];
    }
  }
  return leadingtrack.perp();
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGrooming()
{
   fjOutputs_.reserve(fjInputs_.size());
   zg_.reserve(fjInputs_.size());
   drBranches_.reserve(fjInputs_.size());
   dr12_.reserve(fjInputs_.size());
   sjmass_.reserve(fjInputs_.size());
   logdr12_.reserve(fjInputs_.size());
   logzgdr12_.reserve(fjInputs_.size());
   sjleadingtrack_.reserve(fjInputs_.size());
   jetProfile0_.reserve(fjInputs_.size());
   jetProfile1_.reserve(fjInputs_.size());
   jetProfile2_.reserve(fjInputs_.size());
   jetProfile3_.reserve(fjInputs_.size());
   jetProfile4_.reserve(fjInputs_.size());
   jetProfile5_.reserve(fjInputs_.size());
   jetProfile6_.reserve(fjInputs_.size());
   jetProfile7_.reserve(fjInputs_.size());

   for(fastjet::PseudoJet& jet : fjInputs_) {
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm,999.);
      fastjet::ClusterSequence cs(particles, jet_def);

      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());
      if(tempJets.size()<1) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(-1.);
         drBranches_.push_back(-1.);
         continue;
      }
      
      Double_t jetProfile0=0,jetProfile1=0,jetProfile2=0,jetProfile3=0,jetProfile4=0,jetProfile5=0,jetProfile6=0,jetProfile7=0;
      std::vector<fastjet::PseudoJet> tracks = tempJets[0].constituents();
      Double_t jetEta = tempJets[0].pseudorapidity();
      Double_t jetPhi = tempJets[0].phi();
      Double_t trackEta,trackPhi,trackDelR,delEta,delPhi;
      for(int i=0;i<tempJets[0].constituents().size();i++){
	trackEta = tracks[i].pseudorapidity();
	trackPhi = tracks[i].phi();
	delEta = trackEta-jetEta;
	delPhi = trackPhi-jetPhi;
	trackDelR = TMath::Sqrt(delEta*delEta+delPhi*delPhi);
	if(trackDelR<0.05) jetProfile0+=tracks[i].perp()/tempJets[0].perp();
	if(trackDelR>0.05 && trackDelR<0.1) jetProfile1+=tracks[i].perp()/tempJets[0].perp();
      	if(trackDelR>0.1 && trackDelR<0.15) jetProfile2+=tracks[i].perp()/tempJets[0].perp();
	if(trackDelR>0.15 && trackDelR<0.2) jetProfile3+=tracks[i].perp()/tempJets[0].perp();
	if(trackDelR>0.2 && trackDelR<0.25) jetProfile4+=tracks[i].perp()/tempJets[0].perp();
	if(trackDelR>0.25 && trackDelR<0.3) jetProfile5+=tracks[i].perp()/tempJets[0].perp();
	if(trackDelR>0.3 && trackDelR<0.35) jetProfile6+=tracks[i].perp()/tempJets[0].perp();
	if(trackDelR>0.35 && trackDelR<0.4) jetProfile7+=tracks[i].perp()/tempJets[0].perp();
      }
      jetProfile0_.push_back(jetProfile0/0.05);
      jetProfile1_.push_back(jetProfile1/0.05);
      jetProfile2_.push_back(jetProfile2/0.05);
      jetProfile3_.push_back(jetProfile3/0.05);
      jetProfile4_.push_back(jetProfile4/0.05);
      jetProfile5_.push_back(jetProfile5/0.05);
      jetProfile6_.push_back(jetProfile6/0.05);
      jetProfile7_.push_back(jetProfile7/0.05);

      
      fastjet::contrib::SoftDrop * sd = new fastjet::contrib::SoftDrop(beta_, zcut_, r0_ );
      sd->set_verbose_structure(true);
      fastjet::contrib::Recluster *reclusterer;
      if(fReclusterAlgo == 2) reclusterer = new fastjet::contrib::Recluster(fastjet::kt_algorithm,1,true);
      if(fReclusterAlgo == 1) reclusterer = new fastjet::contrib::Recluster(fastjet::antikt_algorithm,1,true);
      if(fReclusterAlgo == 0) reclusterer = new fastjet::contrib::Recluster(fastjet::cambridge_algorithm,1,true);  
      sd->set_reclustering(true,reclusterer);

      RecursiveParents(tempJets[0],bflagAdded,fAdditionalTracks,kRecursiveAlgo);
      
      fastjet::PseudoJet transformedJet = tempJets[0];
      if ( transformedJet == 0 ) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(-1.);
         drBranches_.push_back(0.);
         if(sd) { delete sd; sd = 0;}
         continue;
      }
      transformedJet = (*sd)(transformedJet);

      double zg = transformedJet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
      int ndrop = transformedJet.structure_of<fastjet::contrib::SoftDrop>().dropped_count();
      //std::vector<double> dropped_symmetry = transformedJet.structure_of<fastjet::contrib::SoftDrop>().dropped_symmetry();
            
      //get distance between the two subjets
      std::vector<fastjet::PseudoJet> subjets;
      if ( transformedJet.has_pieces() ) {
        subjets = transformedJet.pieces();
        fastjet::PseudoJet subjet1 = subjets[0];
        fastjet::PseudoJet subjet2 = subjets[1];
        dr12_.push_back(subjet1.delta_R(subjet2));
	if(subjet1.perp() > subjet2.perp()){
	  sjmass_.push_back(subjet2.m()/subjet2.perp());
	  sjleadingtrack_.push_back(leading_track_pt(subjet2));
	}
	else{
	  sjmass_.push_back(subjet1.m()/subjet1.perp());
	  sjleadingtrack_.push_back(leading_track_pt(subjet1));
	}

	if(subjet1.delta_R(subjet2)!=0){
	  logdr12_.push_back(TMath::Log(1/subjet1.delta_R(subjet2)));
	  logzgdr12_.push_back(TMath::Log(zg*subjet1.delta_R(subjet2)));
	}
	else{
	  logdr12_.push_back(0.);
	  logzgdr12_.push_back(0.);
	}
      } else {
        dr12_.push_back(-1.);
	sjmass_.push_back(-1.);
	sjleadingtrack_.push_back(-1.);
	logdr12_.push_back(0.);
	logzgdr12_.push_back(0.);
      }

      fjOutputs_.push_back( transformedJet ); //put CA reclusterd jet after softDrop into vector
      zg_.push_back(zg);
      drBranches_.push_back(ndrop);

      if(sd) { delete sd; sd = 0;}
   }
   return fjOutputs_;
}

//_________________________________________________________________________
void softDropGroomer::RecursiveParents(fastjet::PseudoJet fJet, Bool_t bflagAdded, Int_t fAdditionalTracks, Int_t ReclusterAlgo){
 
  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  double xflagalgo=0; 
  
  fInputVectors = fJet.constituents();
  if(bflagAdded && fJet.m()!=0){
    //fastjet::PseudoJet MyJet;
    fastjet::PseudoJet PseudoTracksCMS;

    //add tracks to the jet prior to the reclusterer in case of iterative mapping of splittings
    //MyJet.reset(fJet.four_mom());
    Double_t omegac=0.5*fqhat*fxlength*fxlength/0.2;
    Double_t thetac=TMath::Sqrt(12*0.2/(fqhat*TMath::Power(fxlength,3)));
    Double_t xQs=TMath::Sqrt(fqhat*fxlength);				
   

    for(Int_t i=0;i<fAdditionalTracks;i++){

      Double_t ppx,ppy,ppz,kTscale,lim2o,lim1o;
      Double_t lim2=xQs;   
      Double_t lim1=10000;
    
      //generation of kT according to 1/kT^4, with minimum QS=2 GeV and maximum ~sqrt(ptjet*T)	
      TF1 *fTf1Kt= new TF1("fTf1Kt","1/(x*x*x*x)",lim2,lim1);
      kTscale=fTf1Kt->GetRandom();
      //generation within the jet cone
    
      //generation of w according to 1/w, with minimum wc
      //omega needs to be larger than kT so to have well defined angles
      lim2o=kTscale;
      lim1o=kTscale/TMath::Sin(0.1);
      TF1 *fTf1Omega= new TF1("fTf1Omega","1/x",lim2o,lim1o);
      Double_t omega=fTf1Omega->GetRandom();
     
      Double_t sinpptheta=kTscale/omega;
      Double_t pptheta=TMath::ASin(sinpptheta);
      //cout<<"angle_omega_kt"<<pptheta<<" "<<omega<<" "<<kTscale<<endl;
      if(pptheta>0.4) continue;
     
      PseudoTracksCMS.reset(kTscale/TMath::Sqrt(2),kTscale/TMath::Sqrt(2),omega*TMath::Cos(pptheta),omega);
      //boost the particle in the rest frame of the jet to the lab frame
      //cout<<"jet pt = "<<fJet.perp()<<endl<<"jet E = "<<fJet.E()<<endl<<"jet M = "<<fJet.m()<<endl;
      fastjet::PseudoJet PseudoTracksLab=PseudoTracksCMS.boost(fJet);
      Double_t modJetpt = fJet.perp()+PseudoTracksLab.perp();
      Double_t nTracks = (fJet.constituents()).size();
      PseudoTracksLab.set_user_index(i+nTracks+100);											 
      fInputVectors.push_back(PseudoTracksLab);
      injectedtrack_z_.push_back(PseudoTracksLab.perp()/modJetpt);
      injectedtrack_theta_.push_back(pptheta);
      injectedtrack_pt_.push_back(PseudoTracksLab.perp());
      //in the frame of the jet 
      //xflagAdded=1;
    }

  }
  
  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
 
  if(ReclusterAlgo==0){
    jetalgo=fastjet::cambridge_algorithm;
  }
  if(ReclusterAlgo==1){
    jetalgo=fastjet::antikt_algorithm;
  } 
  if(ReclusterAlgo==2){ 
    jetalgo=fastjet::kt_algorithm ;
  }     
  fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 


  fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
  std::vector<fastjet::PseudoJet>   fOutputJets;
  fOutputJets.clear();
  fOutputJets=fClustSeqSA.inclusive_jets(0);
  double jet_pT = fOutputJets[0].perp();
  
  fastjet::PseudoJet jj;
  fastjet::PseudoJet j1;
  fastjet::PseudoJet j2;
  jj=fOutputJets[0];
  int n=0;
  while(jj.has_parents(j1,j2)){
    n++;
    if(j1.perp() < j2.perp()) swap(j1,j2);
    double delta_R=j1.delta_R(j2);
    double z=j2.perp()/(j1.perp()+j2.perp());
    recur_jetpt_.push_back(jet_pT);
    recur_logdr12_.push_back(TMath::Log(1.0/delta_R));
    recur_logzgdr12_.push_back(TMath::Log(z*delta_R));
    recur_n_.push_back(n);
    jj=j1;
  } 




  

 
  
  
  return;

  
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGroomingWithJewelSub(jetCollection &c, std::vector<fastjet::PseudoJet> particlesDummy)
{
  return doGroomingWithJewelSub(c.getJet(),particlesDummy);
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGroomingWithJewelSub(std::vector<fastjet::PseudoJet> v, std::vector<fastjet::PseudoJet> particlesDummy)
{
  setInputJets(v);
  return doGroomingWithJewelSub(particlesDummy);
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGroomingWithJewelSub(std::vector<fastjet::PseudoJet> particlesDummy)
{
   fjOutputs_.reserve(fjInputs_.size());
   zg_.reserve(fjInputs_.size());
   drBranches_.reserve(fjInputs_.size());
   dr12_.reserve(fjInputs_.size());
   constituents_.reserve(fjInputs_.size());
   constituents1_.reserve(fjInputs_.size());
   constituents2_.reserve(fjInputs_.size());
   
   for(fastjet::PseudoJet& jet : fjInputs_) {
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
      fastjet::ClusterSequence cs(particles, jet_def);

      //std::cout << "reclustered with CA" << std::endl;
      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());
      if(tempJets.size()<1) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(-1.);
         drBranches_.push_back(-1.);
         constituents1_.push_back(std::vector<fastjet::PseudoJet>());
         constituents2_.push_back(std::vector<fastjet::PseudoJet>());
         continue;
      }

      fastjet::PseudoJet CurrentJet = tempJets[0];
      RecursiveParents(tempJets[0],bflagAdded,fAdditionalTracks,kRecursiveAlgo);
      fastjet::PseudoJet Part1, Part2;
      fastjet::PseudoJet sj1, sj2;
      double zg = 0.;
      int ndrop = 0;

      // std::cout << "start grooming procedure" << std::endl;
      while(CurrentJet.has_parents(Part1, Part2)) {

        if (CurrentJet.pt2() <= 0) break;
        
        zg = 0.;

        double deltaRsq = Part1.squared_distance(Part2);
        double cut = zcut_ * std::pow(deltaRsq / r0_*r0_, 0.5*beta_);

        sj1 = GetCorrectedJet(Part1,particlesDummy);
        sj2 = GetCorrectedJet(Part2,particlesDummy);
        
        if(sj1.pt() + sj2.pt() > 0 && sj1.E()>0. && sj2.E()>0. && sj1.m()>0. && sj2.m()>0.)
          zg = min(sj1.pt(), sj2.pt()) / (sj1.pt() + sj2.pt());
        
        if(zg<cut) {
          if(sj1.pt() > sj2.pt())
            CurrentJet = Part1;
          else
            CurrentJet = Part2;
          zg = -0.;
          ++ndrop;
        } else {
          break;
        }
      }
      
      //build groomed jet
      fastjet::PseudoJet transformedJet;
      if(zg>0.)
        transformedJet = fastjet::PseudoJet(sj1.px()+sj2.px(),sj1.py()+sj2.py(),sj1.pz()+sj2.pz(),sj1.E()+sj2.E());

      if ( transformedJet == 0 ) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(0.);
	 dr12_.push_back(0.);
         drBranches_.push_back(-1.);
         constituents1_.push_back(std::vector<fastjet::PseudoJet>());
         constituents2_.push_back(std::vector<fastjet::PseudoJet>());
      } else {
        //storing unsubtracted constituents of the subjets
        if(Part1.has_constituents ()) constituents1_.push_back(Part1.constituents());
        if(Part2.has_constituents ()) constituents2_.push_back(Part2.constituents());
        
        //get distance between the two subjets
        //double deltaR = std::sqrt((sj1.eta() -  sj2.eta())*(sj1.eta() - sj2.eta())  (sj1.delta_phi_to(sj2))*(sj2.delta_phi_to(sj1)));
        double deltaR = sj1.delta_R(sj2);
        
        fjOutputs_.push_back( transformedJet ); //put CA reclusterd jet after softDrop into vector
        zg_.push_back(zg);
        dr12_.push_back(deltaR);
        drBranches_.push_back(ndrop);
      }
   } //jet loop
  return fjOutputs_;
}

Double_t softDropGroomer::RelativePhi(Double_t mphi,Double_t vphi){

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}

#endif
