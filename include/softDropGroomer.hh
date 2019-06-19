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
#include "fastjet/contrib/Recluster.hh"

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
  std::vector<double>             recur_erad_;
  std::vector<double>             recur_z_;
  std::vector<double>             sjleadingtrack_;
  std::vector<double>             injectedtrack_z_;
  std::vector<double>             injectedtrack_theta_;
  std::vector<double>             injectedtrack_pt_;

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
  void RecursiveParents(fastjet::PseudoJet fJet,Bool_t bflagAdded, Int_t fAdditionalTracks, Int_t RecursiveAlgo);
  std::vector<fastjet::PseudoJet> getGroomedJets() const;
  std::vector<double> getZgs() const;
  std::vector<int>    getNDroppedSubjets() const;
  std::vector<double> getDR12() const;
  std::vector<double> getSubJetMass() const;
  std::vector<double> getLogDR12() const;
  std::vector<double> getLogZgDR12() const;
  std::vector<double> getRecur_JetPt() const;
  std::vector<double> getRecur_LogDR12() const;
  std::vector<double> getRecur_LogZgDR12() const;
  std::vector<double> getRecur_z() const;
  std::vector<double> getRecur_Erad() const;
  std::vector<double> getInjectedz() const;
  std::vector<double> getInjectedtheta() const;
  std::vector<double> getInjectedpt() const;

  std::vector<int> getRecur_N() const;
  std::vector<double> getSubJetLeadingTrackPt() const;

  std::vector<fastjet::PseudoJet> doGrooming(jetCollection &c);
  std::vector<fastjet::PseudoJet> doGrooming(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> doGrooming();
  
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents() {return constituents_;}
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents1() {return constituents1_;}
  std::vector<std::vector<fastjet::PseudoJet>> getConstituents2() {return constituents2_;}

  std::vector<fastjet::PseudoJet> doGroomingWithJewelSub(jetCollection &c, std::vector<fastjet::PseudoJet> particlesDummy);
  std::vector<fastjet::PseudoJet> doGroomingWithJewelSub(std::vector<fastjet::PseudoJet> v, std::vector<fastjet::PseudoJet> particlesDummy);
  std::vector<fastjet::PseudoJet> doGroomingWithJewelSub(std::vector<fastjet::PseudoJet> particlesDummy);

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

std::vector<double> softDropGroomer::getRecur_z() const
{
   return recur_z_;
}

std::vector<double> softDropGroomer::getRecur_Erad() const
{
   return recur_erad_;
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

std::vector<fastjet::PseudoJet> softDropGroomer::doGrooming(jetCollection &c)
{
   return doGrooming(c.getJet());
}

std::vector<fastjet::PseudoJet> softDropGroomer::doGrooming(std::vector<fastjet::PseudoJet> v)
{
   setInputJets(v);
   return doGrooming();
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

   //NEW STUFF that is not properly treated
   // logdr12_.reserve(fjInputs_.size());
   // logzgdr12_.reserve(fjInputs_.size());
   // sjleadingtrack_.reserve(fjInputs_.size());
   //--
   
   for(fastjet::PseudoJet& jet : fjInputs_) {
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm,fastjet::JetDefinition::max_allowable_R);
      fastjet::ClusterSequence cs(particles, jet_def);

      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());
      if(tempJets.size()<1) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(-1.);
         drBranches_.push_back(-1.);
         dr12_.push_back(-1.);
         // sjmass_.push_back(-1.);
         continue;
      }
      
      fastjet::contrib::SoftDrop * sd = new fastjet::contrib::SoftDrop(beta_, zcut_, r0_ );
      sd->set_verbose_structure(true);

      // fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
      // if(fReclusterAlgo == 1)      jetalgo=fastjet::antikt_algorithm;
      // else if(fReclusterAlgo == 2) jetalgo=fastjet::kt_algorithm;
      // else if(fReclusterAlgo == 3) jetalgo=fastjet::genkt_algorithm;
      
      // fastjet::Recluster recluster(jetalgo,1,fastjet::Recluster::keep_only_hardest);
      // // if(fReclusterAlgo == 3) recluster = fastjet::Recluster(jetalgo,1,0.5);
      // sd->set_reclustering(true,&recluster);

      // RecursiveParents(tempJets[0],bflagAdded,fAdditionalTracks,kRecursiveAlgo); //do we have to do this here?
      
      fastjet::PseudoJet transformedJet = tempJets[0];
      if ( transformedJet == 0 ) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(-1.);
         drBranches_.push_back(0.);
         dr12_.push_back(0.);
         sjmass_.push_back(0.);
         if(sd) { delete sd; sd = 0;}
         //if(reclusterer) { delete reclusterer; sd = 0;}
         continue;
      }
      transformedJet = (*sd)(transformedJet);

      double zg = transformedJet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
      int ndrop = transformedJet.structure_of<fastjet::contrib::SoftDrop>().dropped_count();
      zg_.push_back(zg);
      drBranches_.push_back(ndrop);

      //get distance between the two subjets
      std::vector<fastjet::PseudoJet> subjets;
      if ( transformedJet.has_pieces() ) {
        subjets = transformedJet.pieces();
        fastjet::PseudoJet subjet1 = subjets[0];
        fastjet::PseudoJet subjet2 = subjets[1];
        dr12_.push_back(subjet1.delta_R(subjet2));

	// if(subjet1.perp() > subjet2.perp()){
	//   sjmass_.push_back(subjet2.m()/subjet2.perp());
	//   sjleadingtrack_.push_back(leading_track_pt(subjet2));
	// }
	// else{
	//   sjmass_.push_back(subjet1.m()/subjet1.perp());
	//   sjleadingtrack_.push_back(leading_track_pt(subjet1));
	// }

	// if(subjet1.delta_R(subjet2)!=0){
	//   logdr12_.push_back(TMath::Log(1/subjet1.delta_R(subjet2)));
	//   logzgdr12_.push_back(TMath::Log(zg*subjet1.delta_R(subjet2)));
	// }
	// else{
	//   logdr12_.push_back(0.);
	//   logzgdr12_.push_back(0.);
	// }
      } else {
        dr12_.push_back(-1.);
	// sjmass_.push_back(-1.);
	// sjleadingtrack_.push_back(-1.);
	// logdr12_.push_back(0.);
	// logzgdr12_.push_back(0.);
      }

      fjOutputs_.push_back( transformedJet ); //put CA reclustered jet after softDrop into vector

      if(sd) { delete sd; sd = 0;}
      //   if(reclusterer) { delete reclusterer; sd = 0;}
   }
   return fjOutputs_;
}

//_________________________________________________________________________
void softDropGroomer::RecursiveParents(fastjet::PseudoJet fJet, Bool_t bflagAdded, Int_t fAdditionalTracks, Int_t RecursiveAlgo){
 
  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  //double xflagalgo=0; 
  Float_t kpval=0;
  fInputVectors = fJet.constituents();
  if(bflagAdded && fJet.m()!=0){
    //fastjet::PseudoJet MyJet;
    fastjet::PseudoJet PseudoTracksCMS;

    //add tracks to the jet prior to the reclusterer in case of iterative mapping of splittings
    //MyJet.reset(fJet.four_mom());
    //Double_t omegac=0.5*fqhat*fxlength*fxlength/0.2;
    //Double_t thetac=TMath::Sqrt(12*0.2/(fqhat*TMath::Power(fxlength,3)));
    Double_t xQs=TMath::Sqrt(fqhat*fxlength);				
   

    for(Int_t i=0;i<fAdditionalTracks;i++){

      Double_t kTscale,lim2o,lim1o;
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
  
  fastjet::JetAlgorithm jetalgo(fastjet::genkt_algorithm);
  if(RecursiveAlgo==0){
    //CA
    kpval=0;
  }
  if(RecursiveAlgo==1){
    //kT
    kpval=-1;
  } 
  if(RecursiveAlgo==2){ 
    //antikt
    kpval=1;
  }
  if(RecursiveAlgo==3){ 
    //tform-ordered
    kpval=0.5;
  }
  
  fastjet::JetDefinition fJetDef(jetalgo, 1., kpval, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 

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
    recur_erad_.push_back(j1.e()+j2.e());
    recur_z_.push_back(z);
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

#endif
