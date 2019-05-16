//! JetTools workshop Bergen 2019
//! studying jet quenching in working group 2 

#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
//! root headers to produce histograms
#include "TH1D.h"
#include "TFile.h"

//! list of variables that we want to study -
//! leading hadron, dijet values, leading photon,
//! substructure - for the jets, we want to look at subjet given by the recursive study

//! list of observables -
//! delta phi - trigger object (leading hadron or photon), charged particles  
//! delta phi - trigger object and recoil leading jet 
//! delta phi - trigger object and recoil sub-leading jet 
//! delta phi - trigger object and all recoil jets normalized per trigger 
//! given trigger object, on the recoil side,
//! 4-momentum jet axis delta phi w/ WTA axis - redoing clustering for jet constituents
//! test that effect with clustering and jet finding or jet finding and clustering ()
//! 4-momentum jet axis vs leading jet subjet axis
//! WTA jet axis vs leading jet subjet axis [0.4 jet radii and 0.1 subjet]
//! leading and subleading subjets within the jet - store their axis information for delta phi
//! do the two subjet analyis
//! softdrop zg and rg
//! star dijet hardcore cuts - Aj for both hardcore and matched and delta phi between matched
//! all the other histograms are for full jets that are not entirely matched 
//! 



namespace Rivet {

  class JETTOOLS_COINCIDENCE : public Analysis {
  public:

    /// Constructor
    JETTOOLS_COINCIDENCE(string name = "JETTOOLS_COINCIDENCE")
      : Analysis(name)
    {
      setNeedsCrossSection(true);
      _mode=0;
      //! mode 0 will be for hadron+jet coincidence
      //! mode 1 will be for photon+jet events
      //! mode 2 will be for dijet coincidence  
    }


    //! http://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.111501
    //! Anti kT jet is reclustered using CA and min(pT1, pT2)/(pT1+pT2) > z_cut * (R12/Rj)^\beta
    //! for sub-jets that satisfy this condition, z_g = min(pT1,pT2)/(pT1+pT2)
    std::pair<double,double>  SoftDrop(const fastjet::PseudoJet &j, double jetR, double z_cut, double beta)
    {
      //! give the soft drop groomer a short name
      //! Use a symmetry cut z > z_cut R^beta
      //! By default, there is no mass-drop requirement
      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
      fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, jetR);
      PseudoJet sd_jet = sd(j);    
      if(sd_jet != 0){
	double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	return std::pair<double,double>(z,r);
      }else 
	return std::pair<double,double>(-1,-1);
    }
    

    void init() {

      //! jet selections  
      _pTCut = 10;
      _jetR = 0.4;
      
      //! parameters for splitting function
      z_cut=0.1;
      beta = 0;

      //! event normalizations
      _Ntriggers = 0.0;
      _Ndijets = 0.0;


      /// final state and jet projection
      FinalState fs(-5.0, 5.0, 0.2*GeV);
      addProjection(fs, "FS");
      FinalState HCfs(-5.0, 5.0, 2.0*GeV);
      addProjection(HCfs, "HCFS");
      VetoedFinalState vfs(fs);
      ChargedFinalState cfs(-5.0, 5.0, 0.2*GeV);
      addProjection(cfs, "CFS");

      if (_mode == 0){
	//! this is just hadron+jet
	LeadingParticlesFinalState hadronfs(FinalState(-1.0, 1.0, 10.0*GeV));
	hadronfs.addParticleId(111);
	hadronfs.addParticleId(211);
	hadronfs.addParticleId(-211);
	addProjection(hadronfs, "LeadingHadron");
	vfs.addVetoOnThisFinalState(hadronfs);
	addProjection(vfs, "VFS");
	// addProjection(FastJets(vfs, FastJets::ANTIKT, _jetR), "Jets");
	fout = new TFile("jettools_jewel_AuAu_cent05_wrecoil_200_hadron_jet_study_correlation_histograms.root", "RECREATE");
      } else if (_mode == 1) {
	//! this is photon+jet 
	LeadingParticlesFinalState photonfs(FinalState(-2.5, 2.5, 10.0*GeV));
	photonfs.addParticleId(PID::PHOTON);
	addProjection(photonfs, "LeadingPhoton");
	vfs.addVetoOnThisFinalState(photonfs);
	addProjection(vfs, "VFS");
	// addProjection(FastJets(vfs, FastJets::ANTIKT, _jetR), "Jets");
	fout = new TFile("jettools_jewel_AuAu_cent05_wrecoil_200_gamma_jet_study_correlation_histograms.root", "RECREATE");
      } else if (_mode == 2) {
	//! this is dijet+jet
	// addProjection(FastJets(HCfs, FastJets::ANTIKT, _jetR), "JetsHC");	  
	// addProjection(FastJets(fs, FastJets::ANTIKT, _jetR), "Jets");	  
	fout = new TFile("jettools_jewel_AuAu_cent05_wrecoil_200_dijet_study_correlation_histograms.root", "RECREATE");
      }

      
      //! declare all the histograms that you want to study -
      hDeltaPhi_Trigger_ChargedParticles = new TH1D("hDeltaPhi_Trigger_ChargedParticles", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_LeadingRecoilJet = new TH1D("hDeltaPhi_Trigger_LeadingRecoilJet", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis = new TH1D("hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis_largeR = new TH1D("hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis_largeR", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_SubLeadingRecoilJet = new TH1D("hDeltaPhi_Trigger_SubLeadingRecoilJet", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_AllRecoilJets = new TH1D("hDeltaPhi_Trigger_AllRecoilJets", "", 60, 0, 3.3);
      hDeltaR_RecoilJet_JetAxis_WTAAxis = new TH1D("hDeltaR_RecoilJet_JetAxis_WTAAxis", "", 60, 0, 0.4);
      hDeltaR_RecoilJet_JetAxis_WTAAxis_largeR = new TH1D("hDeltaR_RecoilJet_JetAxis_WTAAxis_largeR", "", 60, 0, 0.4);
      hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis = new TH1D("hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis", "", 60, 0, 0.4);
      hDeltaR_RecoilJet_GroomedJetAxis_WTAAxis_largeR = new TH1D("hDeltaR_RecoilJet_GroomedJetAxis_WTAAxis_largeR", "", 60, 0, 0.4);
      h_RecoilJet_TwoSubJet_Theta_R0p1 = new TH1D("h_RecoilJet_TwoSubJet_Theta_R0p1", "", 60, 0, 0.5);
      h_RecoilJet_TwoSubJet_Z_R0p1 = new TH1D("h_RecoilJet_TwoSubJet_Z_R0p1", "", 60, 0, 0.5);
      hSoftDrop_Zg_LeadingRecoilJet = new TH1D("hSoftDrop_Zg_LeadingRecoilJet", "", 60, 0, 0.5);
      hSoftDrop_Rg_LeadingRecoilJet = new TH1D("hSoftDrop_Rg_LeadingRecoilJet", "", 60, 0, 0.5);

      h_HardCore_Dijet_Aj = new TH1D("h_HardCore_Dijet_Aj", "", 20, 0, 1);
      h_HardCore_RecoilJetYield = new TH1D("h_HardCore_RecoilJetYield", "", 40, 0, 40);
      hDeltaPhi_HardCore_Dijet = new TH1D("hDeltaPhi_HardCore_Dijet", "", 60, 0, 3.3);
      h_HardCore_TriggerJet_Zg = new TH1D("h_HardCore_TriggerJet_Zg", "", 10, 0, 0.5);
      h_HardCore_RecoilJet_Zg = new TH1D("h_HardCore_RecoilJet_Zg", "", 10, 0, 0.5);
      h_HardCore_TriggerJet_Rg = new TH1D("h_HardCore_TriggerJet_Rg", "", 10, 0, 0.5);
      h_HardCore_RecoilJet_Rg = new TH1D("h_HardCore_RecoilJet_Rg", "", 10, 0, 0.5);
      h_HardCore_TriggerJet_Twosubjet_Z_R0p1 = new TH1D("h_HardCore_TriggerJet_Twosubjet_Z_R0p1", "", 10, 0, 0.5);
      h_HardCore_RecoilJet_Twosubjet_Z_R0p1 = new TH1D("h_HardCore_RecoilJet_Twosubjet_Z_R0p1", "", 10, 0, 0.5);
      h_HardCore_TriggerJet_Twosubjet_Theta_R0p1 = new TH1D("h_HardCore_TriggerJet_Twosubjet_Theta_R0p1", "", 10, 0, 0.5);
      h_HardCore_RecoilJet_Twosubjet_Theta_R0p1 = new TH1D("h_HardCore_RecoilJet_Twosubjet_Theta_R0p1", "", 10, 0, 0.5);
      h_HardCore_TriggerJet_deltaR_JetAxis_LeadSubJetAxis = new TH1D("h_HardCore_TriggerJet_deltaR_JetAxis_LeadSubJetAxis", "", 60, 0, 0.4);
      h_HardCore_RecoilJet_deltaR_JetAxis_LeadSubJetAxis = new TH1D("h_HardCore_RecoilJet_deltaR_JetAxis_LeadSubJetAxis", "", 60, 0, 0.4);

      h_Matched_Dijet_Aj = new TH1D("h_Matched_Dijet_Aj", "", 20, 0, 1);
      h_Matched_RecoilJetYield = new TH1D("h_Matched_RecoilJetYield", "", 40, 0, 40);
      hDeltaPhi_Matched_Dijet = new TH1D("hDeltaPhi_Matched_Dijet", "", 60, 0, 3.3);
      h_Matched_TriggerJet_Zg = new TH1D("h_Matched_TriggerJet_Zg", "", 10, 0, 0.5);
      h_Matched_RecoilJet_Zg = new TH1D("h_Matched_RecoilJet_Zg", "", 10, 0, 0.5);
      h_Matched_TriggerJet_Rg = new TH1D("h_Matched_TriggerJet_Rg", "", 10, 0, 0.5);
      h_Matched_RecoilJet_Rg = new TH1D("h_Matched_RecoilJet_Rg", "", 10, 0, 0.5);
      h_Matched_TriggerJet_Twosubjet_Z_R0p1 = new TH1D("h_Matched_TriggerJet_Twosubjet_Z_R0p1", "", 10, 0, 0.5);
      h_Matched_RecoilJet_Twosubjet_Z_R0p1 = new TH1D("h_Matched_RecoilJet_Twosubjet_Z_R0p1", "", 10, 0, 0.5);
      h_Matched_TriggerJet_Twosubjet_Theta_R0p1 = new TH1D("h_Matched_TriggerJet_Twosubjet_Theta_R0p1", "", 10, 0, 0.5);
      h_Matched_RecoilJet_Twosubjet_Theta_R0p1 = new TH1D("h_Matched_RecoilJet_Twosubjet_Theta_R0p1", "", 10, 0, 0.5);
      h_Matched_TriggerJet_deltaR_JetAxis_LeadSubJetAxis = new TH1D("h_Matched_TriggerJet_deltaR_JetAxis_LeadSubJetAxis", "", 60, 0, 0.4);
      h_Matched_RecoilJet_deltaR_JetAxis_LeadSubJetAxis = new TH1D("h_Matched_RecoilJet_deltaR_JetAxis_LeadSubJetAxis", "", 60, 0, 0.4);

      hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis = new TH1D("hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis", "", 60, 0, 0.6);
      hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis = new TH1D("hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis", "", 60, 0, 3.3);
    
      hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis_largeR = new TH1D("hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis_largeR", "", 60, 0, 0.4);
      hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis_largeR = new TH1D("hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis_largeR", "", 60, 0, 3.3);
      

      hTriggerPhi = new TH1D("hTriggerPhi", "", 100, -7, 7);
      hJetPhi = new TH1D("hJetPhi", "", 100, -7, 7);
      

    }//! init method

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      ParticleVector Trigger;
      const FinalState& chargedfs = applyProjection<FinalState>(event, "CFS");      
      const FinalState& finalstate = applyProjection<FinalState>(event, "FS");      
      
      //! event common fastjet objects - 
      fastjet::Selector select_eta  = fastjet::SelectorAbsEtaMax(1.0 - _jetR);
      fastjet::Selector select_pt_lo   = fastjet::SelectorPtMin(_pTCut);
      fastjet::Selector select_both = select_pt_lo && select_eta;
      fastjet::JetDefinition jetd(fastjet::antikt_algorithm, _jetR);
      //! get wta axis for the jet -
      fastjet::JetDefinition wta_jet_def(fastjet::antikt_algorithm, _jetR, fastjet::WTA_pt_scheme);
      //! get wta axis for the jet -
      fastjet::JetDefinition wta_larger_jet_def(fastjet::antikt_algorithm, 1.0, fastjet::WTA_pt_scheme);
      //! get subjets -
      fastjet::JetDefinition sub_jetd(fastjet::antikt_algorithm, 0.1);
      //! fill the softdrop variable
      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
      fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, _jetR);

      //! Now do the analysis for both triggered samples -
      if(_mode == 1 || _mode == 0){

	if (_mode == 0){
	  //! leading hadron -
	  Trigger = applyProjection<FinalState>(event, "LeadingHadron").particles();
	  if(Trigger.size() != 1) vetoEvent;
	} else if (_mode == 1){
	  //! leading photon 
	  Trigger = applyProjection<FinalState>(event, "LeadingPhoton").particles();
	  if (Trigger.size() != 1) vetoEvent;
	}

	Particle trigger = Trigger.at(0);
	
	if(trigger.momentum().pt() < 10)
	  vetoEvent;

	_Ntriggers+=weight;
	
	double trigger_phi = trigger.momentum().phi();
	hTriggerPhi->Fill(trigger_phi);
	
	//! get all the particles from the vetoed final state
	const FinalState& vetoedfs = applyProjection<FinalState>(event, "VFS");

	PseudoJets vetoedparticles;
	foreach(Particle part, vetoedfs.particles()){
	  FourMomentum particle = part.momentum();
	  double px, py, pz, E;
	  px = particle.px();
	  py = particle.py();
	  pz = particle.pz();
	  E  = particle.E();	  
	  PseudoJet consts(px,py,pz,E);
	  vetoedparticles.push_back(consts);
	}

	//! get subjets -
	fastjet::ClusterSequence cs_recoiljets(vetoedparticles, jetd);
	PseudoJets recoilJETS = sorted_by_pt(cs_recoiljets.inclusive_jets());
	PseudoJets rJets = select_both(recoilJETS);

	//! enforce the delta phi requirement to be > pi/2
	vector<fastjet::PseudoJet> recoilJets;
	foreach (const fastjet::PseudoJet& jet, rJets){
	  if(deltaPhi(trigger.momentum().phi(), jet.phi()) > PI/2)
	    recoilJets.push_back(jet);
	}
	
	const int _Njet = recoilJets.size();
	if(_Njet==0) vetoEvent;

	//! fill the histograms for the delta phi between trigger object and all charged particles
	/// now loop through all the charged particles 
	foreach (Particle part, chargedfs.particles()){
	  hDeltaPhi_Trigger_ChargedParticles->Fill(deltaPhi(trigger.momentum().phi(), part.momentum().phi()), weight);	  
	}
	fastjet::PseudoJet& leading_recoiljet = recoilJets[0];
	hJetPhi->Fill(leading_recoiljet.phi());
	hDeltaPhi_Trigger_LeadingRecoilJet->Fill(deltaPhi(trigger.momentum().phi(), leading_recoiljet.phi()), weight);

	if(_Njet > 1){
	  fastjet::PseudoJet& subleading_recoiljet = recoilJets[1];
	  hDeltaPhi_Trigger_SubLeadingRecoilJet->Fill(deltaPhi(trigger.momentum().phi(), subleading_recoiljet.phi()), weight);
	}

	foreach(const fastjet::PseudoJet& jet, recoilJets){
	  hDeltaPhi_Trigger_AllRecoilJets->Fill(deltaPhi(trigger.momentum().phi(), jet.phi()), weight);
	}
	
	PseudoJets lead_recoiljet_consts = leading_recoiljet.constituents();
	fastjet::ClusterSequence cs_wta_leadingrecoiljet(lead_recoiljet_consts, wta_jet_def);
	PseudoJets wta_lead_recoiljets = cs_wta_leadingrecoiljet.inclusive_jets();
	PseudoJet wta_lead_recoiljet = wta_lead_recoiljets.at(0);

	hDeltaR_RecoilJet_JetAxis_WTAAxis->Fill(leading_recoiljet.delta_R(wta_lead_recoiljet), weight);
	hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis->Fill(deltaPhi(trigger.momentum().phi(), wta_lead_recoiljet.phi()), weight);

	fastjet::ClusterSequence cs_wta_larger_leadingrecoiljet(lead_recoiljet_consts, wta_larger_jet_def);
	PseudoJets wta_larger_lead_recoiljets = cs_wta_larger_leadingrecoiljet.inclusive_jets();
	PseudoJet wta_larger_lead_recoiljet = wta_larger_lead_recoiljets.at(0);

	hDeltaR_RecoilJet_JetAxis_WTAAxis_largeR->Fill(leading_recoiljet.delta_R(wta_larger_lead_recoiljet), weight);
	hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis_largeR->Fill(deltaPhi(trigger.momentum().phi(), wta_larger_lead_recoiljet.phi()), weight);
	
	fastjet::ClusterSequence clust_seq_lead_recoil_subjet(lead_recoiljet_consts, sub_jetd);	
	PseudoJets subjets = sorted_by_pt(clust_seq_lead_recoil_subjet.inclusive_jets());

	hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis->Fill(leading_recoiljet.delta_R(subjets.at(0)), weight);

	if(subjets.size() >= 2){
	  //! three subjets - no grooming. just straight up two subjets 
	  double zsj = subjets.at(1).pt()/(subjets.at(0).pt()+subjets.at(1).pt());
	  double rsj = subjets.at(0).delta_R(subjets.at(1));
	  h_RecoilJet_TwoSubJet_Theta_R0p1->Fill(rsj, weight);
	  h_RecoilJet_TwoSubJet_Z_R0p1->Fill(zsj, weight);
	}

	PseudoJet sd_jet = sd(leading_recoiljet);
	double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();

	if(sd_jet!=0){
	  hDeltaR_RecoilJet_GroomedJetAxis_WTAAxis_largeR->Fill(deltaPhi(sd_jet.phi(), wta_lead_recoiljet.phi()), weight);	
	  hSoftDrop_Zg_LeadingRecoilJet->Fill(z, weight);
	  hSoftDrop_Rg_LeadingRecoilJet->Fill(r, weight);
	}


      }//! mode = 1 or 0 
      
      
      if (_mode == 2){
	
	
							 
	//! dijet selections - only run on dijet samples ofcourse 
	PseudoJets pJet_hardcore;
	//! first add the hard core objects in the event. 	      
	foreach ( const Particle& p, finalstate.particles()) {
	  if(p.pt() >= 2.0*GeV){
	    pJet_hardcore.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.E()));
	  }
	}

	
	//! Hard core Jets 	
	fastjet::ClusterSequence cs_hardcore(pJet_hardcore, jetd);
	PseudoJets Jets_hardcore = sorted_by_pt(cs_hardcore.inclusive_jets(_pTCut));
	PseudoJets HC_Jets = select_both(Jets_hardcore);
	if(HC_Jets.size() < 2){
	  vetoEvent;
	}

	PseudoJets hardcore_dijets;
	if(HC_Jets.at(0).pt() > 16 && HC_Jets.at(1).pt() > 8 && //! dijet pT cut 
	   deltaPhi(HC_Jets.at(0).phi(), HC_Jets.at(1).phi()) > PI/2 //! delta phi requirement
	   ){
	   hardcore_dijets.push_back(HC_Jets.at(0));
	   hardcore_dijets.push_back(HC_Jets.at(1));
	}

	if(hardcore_dijets.size() != 2)
	  vetoEvent;
	
	//! get the matched jet collection -
	PseudoJets pJet;
	//! first add the hard core objects in the event. 	      
	foreach ( const Particle& p, finalstate.particles()) {
	  if(p.pt() >= 0.2*GeV){
	    pJet.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.E()));
	  }
	}

	//! Hard core Jets 	
	fastjet::ClusterSequence cs_fullevent(pJet, jetd);
	PseudoJets Jets = sorted_by_pt(cs_fullevent.inclusive_jets(_pTCut));
	PseudoJets full_Jets = select_both(Jets);
	if(full_Jets.size() < 2){
	  vetoEvent;
	}

	PseudoJets matched_dijets;
	//! only select the dijet that matched the hardcore leading and subleading
	//! match the leading jet
	foreach(const PseudoJet& jet, full_Jets){
	  if(jet.delta_R(hardcore_dijets.at(0)) < _jetR)
	    matched_dijets.push_back(jet);
	}
	//! match the subleading jet 
	foreach(const PseudoJet& jet, full_Jets){
	  if(jet.delta_R(hardcore_dijets.at(1)) < _jetR)
	    matched_dijets.push_back(jet);
	}
	
	if(matched_dijets.size() != 2)
	  vetoEvent;

	_Ndijets+=weight;

	//! now you have both matched and hardore dijet collections - fill the histograms for the observables

	//! Aj - 
	double HC_Aj = (hardcore_dijets.at(0).pt() - hardcore_dijets.at(1).pt())/(hardcore_dijets.at(0).pt() + hardcore_dijets.at(1).pt());
	h_HardCore_Dijet_Aj->Fill(HC_Aj, weight);
	double MC_Aj = (matched_dijets.at(0).pt() - matched_dijets.at(1).pt())/(matched_dijets.at(0).pt() + matched_dijets.at(1).pt());
	h_Matched_Dijet_Aj->Fill(MC_Aj, weight);

	//! recoil jet yield -
	h_HardCore_RecoilJetYield->Fill(hardcore_dijets.at(1).pt(), weight);
	h_Matched_RecoilJetYield->Fill(matched_dijets.at(1).pt(), weight);

	//! correlations in jet phi and also look at the WTA axis for the recoil matched jets
	hDeltaPhi_HardCore_Dijet->Fill(deltaPhi(hardcore_dijets.at(0).phi(), hardcore_dijets.at(1).phi()), weight);
	hDeltaPhi_Matched_Dijet->Fill(deltaPhi(matched_dijets.at(0).phi(), matched_dijets.at(1).phi()), weight);

	//! get WTA axis for matched jet - 
	PseudoJets recoil_matchedjet_consts = matched_dijets.at(1).constituents();
	fastjet::ClusterSequence cs_wta_matchedrecoiljet(recoil_matchedjet_consts, wta_jet_def);
	PseudoJets wta_matched_recoiljets = cs_wta_matchedrecoiljet.inclusive_jets();
	PseudoJet wta_matched_recoiljet = wta_matched_recoiljets.at(0);

	hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis->Fill(matched_dijets.at(1).delta_R(wta_matched_recoiljet), weight);
	hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis->Fill(deltaPhi(matched_dijets.at(0).phi(), wta_matched_recoiljet.phi()), weight);

	fastjet::ClusterSequence cs_wta_larger_matchedrecoiljet(recoil_matchedjet_consts, wta_larger_jet_def);
	PseudoJets wta_larger_matched_recoiljets = cs_wta_larger_matchedrecoiljet.inclusive_jets();
	PseudoJet wta_larger_matched_recoiljet = wta_larger_matched_recoiljets.at(0);

	hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis_largeR->Fill(matched_dijets.at(1).delta_R(wta_larger_matched_recoiljet), weight);
	hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis_largeR->Fill(deltaPhi(matched_dijets.at(0).phi(), wta_larger_matched_recoiljet.phi()), weight);
	
	//! Do the softdrop zg/rg and two subjet theta and z analysis -
	//! loop over the hardcore dijet and get the values for two subjet z and theta
	vector<double> hardcore_zg;
	vector<double> hardcore_rg;
	vector<double> hardcore_twosubjet_z_r0p1;
	vector<double> hardcore_twosubjet_theta_r0p1;
	vector<double> deltaR_hardcorejetaxis_leadsubjetaxis;
	
	foreach(const PseudoJet& jet, hardcore_dijets){

	  PseudoJet sd_jet = sd(jet);
	  double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	  double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	  hardcore_zg.push_back(z);
	  hardcore_rg.push_back(r);
	  
	  PseudoJets jet_consts = jet.constituents();
	  fastjet::ClusterSequence clust_seq_subjet(jet_consts, sub_jetd);	
	  PseudoJets subjets = sorted_by_pt(clust_seq_subjet.inclusive_jets());
	  deltaR_hardcorejetaxis_leadsubjetaxis.push_back(jet.delta_R(subjets.at(0)));

	  if(subjets.size() >= 2){
	    //! three subjets - no grooming. just straight up two subjets 
	    double zsj = subjets.at(1).pt()/(subjets.at(0).pt()+subjets.at(1).pt());
	    double rsj = subjets.at(0).delta_R(subjets.at(1));
	    hardcore_twosubjet_z_r0p1.push_back(zsj);
	    hardcore_twosubjet_theta_r0p1.push_back(rsj);
	  }else{
	    hardcore_twosubjet_z_r0p1.push_back(-999);
	    hardcore_twosubjet_theta_r0p1.push_back(-999);
	  }
	}

	h_HardCore_TriggerJet_Zg->Fill(hardcore_zg.at(0), weight);
	h_HardCore_TriggerJet_Rg->Fill(hardcore_rg.at(0), weight);	
	h_HardCore_RecoilJet_Zg->Fill(hardcore_zg.at(1), weight);
	h_HardCore_RecoilJet_Rg->Fill(hardcore_rg.at(1), weight);	
	h_HardCore_TriggerJet_Twosubjet_Z_R0p1->Fill(hardcore_twosubjet_z_r0p1.at(0), weight);
	h_HardCore_RecoilJet_Twosubjet_Z_R0p1->Fill(hardcore_twosubjet_z_r0p1.at(1), weight);
	h_HardCore_TriggerJet_Twosubjet_Theta_R0p1->Fill(hardcore_twosubjet_theta_r0p1.at(0), weight);
	h_HardCore_RecoilJet_Twosubjet_Theta_R0p1->Fill(hardcore_twosubjet_theta_r0p1.at(1), weight);
	h_HardCore_TriggerJet_deltaR_JetAxis_LeadSubJetAxis->Fill(deltaR_hardcorejetaxis_leadsubjetaxis.at(0), weight);
	h_HardCore_RecoilJet_deltaR_JetAxis_LeadSubJetAxis->Fill(deltaR_hardcorejetaxis_leadsubjetaxis.at(1), weight);
	
	vector<double> matched_zg;
	vector<double> matched_rg;
	vector<double> matched_twosubjet_z_r0p1;
	vector<double> matched_twosubjet_theta_r0p1;
	vector<double> deltaR_matchedjetaxis_leadsubjetaxis;
	
	foreach(const PseudoJet& jet, matched_dijets){

	  PseudoJet sd_jet = sd(jet);
	  double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	  double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	  matched_zg.push_back(z);
	  matched_rg.push_back(r);
	  
	  PseudoJets jet_consts = jet.constituents();
	  fastjet::ClusterSequence clust_seq_subjet(jet_consts, sub_jetd);	
	  PseudoJets subjets = sorted_by_pt(clust_seq_subjet.inclusive_jets());
	  deltaR_matchedjetaxis_leadsubjetaxis.push_back(jet.delta_R(subjets.at(0)));

	  if(subjets.size() >= 2){
	    //! three subjets - no grooming. just straight up two subjets 
	    double zsj = subjets.at(1).pt()/(subjets.at(0).pt()+subjets.at(1).pt());
	    double rsj = subjets.at(0).delta_R(subjets.at(1));
	    matched_twosubjet_z_r0p1.push_back(zsj);
	    matched_twosubjet_theta_r0p1.push_back(rsj);
	  }else {
	    matched_twosubjet_z_r0p1.push_back(-999);
	    matched_twosubjet_theta_r0p1.push_back(-999);
	  }
	}

	
	h_Matched_TriggerJet_Zg->Fill(matched_zg.at(0), weight);
	h_Matched_RecoilJet_Zg->Fill(matched_zg.at(1), weight);
	h_Matched_TriggerJet_Rg->Fill(matched_rg.at(0), weight);
	h_Matched_RecoilJet_Rg->Fill(matched_rg.at(1), weight);	
	h_Matched_TriggerJet_Twosubjet_Z_R0p1->Fill(matched_twosubjet_z_r0p1.at(0), weight);
	h_Matched_RecoilJet_Twosubjet_Z_R0p1->Fill(matched_twosubjet_z_r0p1.at(1), weight);
	h_Matched_TriggerJet_Twosubjet_Theta_R0p1->Fill(matched_twosubjet_theta_r0p1.at(0), weight);
	h_Matched_RecoilJet_Twosubjet_Theta_R0p1->Fill(matched_twosubjet_theta_r0p1.at(1), weight);
	h_Matched_TriggerJet_deltaR_JetAxis_LeadSubJetAxis->Fill(deltaR_matchedjetaxis_leadsubjetaxis.at(0), weight);
	h_Matched_RecoilJet_deltaR_JetAxis_LeadSubJetAxis->Fill(deltaR_matchedjetaxis_leadsubjetaxis.at(1), weight);
	
	
      }//! mode = 2
      

    }//! analyze method 
    
    /// Normalise histograms etc., after the run
    void finalize() {

      fout->cd();

      //! normalize and write histograms 
      if(_mode==0 || _mode == 1){
	
	hDeltaPhi_Trigger_ChargedParticles->Scale(1./_Ntriggers);
	hDeltaPhi_Trigger_ChargedParticles->Write();
	
	hDeltaPhi_Trigger_LeadingRecoilJet->Scale(1./_Ntriggers);
	hDeltaPhi_Trigger_LeadingRecoilJet->Write();
	
	hDeltaPhi_Trigger_SubLeadingRecoilJet->Scale(1./_Ntriggers);
	hDeltaPhi_Trigger_SubLeadingRecoilJet->Write();
	
	hDeltaPhi_Trigger_AllRecoilJets->Scale(1./_Ntriggers);
	hDeltaPhi_Trigger_AllRecoilJets->Write();

	hDeltaR_RecoilJet_JetAxis_WTAAxis->Scale(1./_Ntriggers);
	hDeltaR_RecoilJet_JetAxis_WTAAxis->Write();
	
	hDeltaR_RecoilJet_JetAxis_WTAAxis_largeR->Scale(1./_Ntriggers);
	hDeltaR_RecoilJet_JetAxis_WTAAxis_largeR->Write();
	
	hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis->Scale(1./_Ntriggers);
	hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis->Write();
	
	h_RecoilJet_TwoSubJet_Theta_R0p1->Scale(1./_Ntriggers);
	h_RecoilJet_TwoSubJet_Theta_R0p1->Write();
	
	h_RecoilJet_TwoSubJet_Z_R0p1->Scale(1./_Ntriggers);
	h_RecoilJet_TwoSubJet_Z_R0p1->Write();
	
	hSoftDrop_Zg_LeadingRecoilJet->Scale(1./_Ntriggers);
	hSoftDrop_Zg_LeadingRecoilJet->Write();
	
	hSoftDrop_Rg_LeadingRecoilJet->Scale(1./_Ntriggers);
	hSoftDrop_Rg_LeadingRecoilJet->Write();
	
	hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis->Scale(1./_Ntriggers);
	hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis->Write();
	
	hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis_largeR->Scale(1./_Ntriggers);
	hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis_largeR->Write();
	
	hTriggerPhi->Scale(1./_Ntriggers);
	hTriggerPhi->Write();
	
	hJetPhi->Scale(1./_Ntriggers);
	hJetPhi->Write();
	
      }

      if(_mode == 2){

	h_HardCore_Dijet_Aj->Scale(1./_Ndijets);
	h_HardCore_RecoilJetYield->Scale(1./_Ndijets);
	hDeltaPhi_HardCore_Dijet->Scale(1./_Ndijets);
	h_HardCore_TriggerJet_Zg->Scale(1./_Ndijets);
	h_HardCore_RecoilJet_Zg->Scale(1./_Ndijets);
	h_HardCore_TriggerJet_Rg->Scale(1./_Ndijets);
	h_HardCore_RecoilJet_Rg->Scale(1./_Ndijets);
	h_HardCore_TriggerJet_Twosubjet_Z_R0p1->Scale(1./_Ndijets);
	h_HardCore_RecoilJet_Twosubjet_Z_R0p1->Scale(1./_Ndijets);
	h_HardCore_TriggerJet_Twosubjet_Theta_R0p1->Scale(1./_Ndijets);
	h_HardCore_RecoilJet_Twosubjet_Theta_R0p1->Scale(1./_Ndijets);
	h_HardCore_TriggerJet_deltaR_JetAxis_LeadSubJetAxis->Scale(1./_Ndijets);
	h_HardCore_RecoilJet_deltaR_JetAxis_LeadSubJetAxis->Scale(1./_Ndijets);
	h_Matched_Dijet_Aj->Scale(1./_Ndijets);
	h_Matched_RecoilJetYield->Scale(1./_Ndijets);
	hDeltaPhi_Matched_Dijet->Scale(1./_Ndijets);
	h_Matched_TriggerJet_Zg->Scale(1./_Ndijets);
	h_Matched_RecoilJet_Zg->Scale(1./_Ndijets);
	h_Matched_TriggerJet_Rg->Scale(1./_Ndijets);
	h_Matched_RecoilJet_Rg->Scale(1./_Ndijets);
	h_Matched_TriggerJet_Twosubjet_Z_R0p1->Scale(1./_Ndijets);
	h_Matched_RecoilJet_Twosubjet_Z_R0p1->Scale(1./_Ndijets);
	h_Matched_TriggerJet_Twosubjet_Theta_R0p1->Scale(1./_Ndijets);
	h_Matched_RecoilJet_Twosubjet_Theta_R0p1->Scale(1./_Ndijets);
	h_Matched_TriggerJet_deltaR_JetAxis_LeadSubJetAxis->Scale(1./_Ndijets);
	h_Matched_RecoilJet_deltaR_JetAxis_LeadSubJetAxis->Scale(1./_Ndijets);

	hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis->Scale(1./_Ndijets);
	hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis->Scale(1./_Ndijets);
    
	hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis_largeR->Scale(1./_Ndijets);
	hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis_largeR->Scale(1./_Ndijets);

	
	h_HardCore_Dijet_Aj->Write();
	h_HardCore_RecoilJetYield->Write();
	hDeltaPhi_HardCore_Dijet->Write();
	h_HardCore_TriggerJet_Zg->Write();
	h_HardCore_RecoilJet_Zg->Write();
	h_HardCore_TriggerJet_Rg->Write();
	h_HardCore_RecoilJet_Rg->Write();
	h_HardCore_TriggerJet_Twosubjet_Z_R0p1->Write();
	h_HardCore_RecoilJet_Twosubjet_Z_R0p1->Write();
	h_HardCore_TriggerJet_Twosubjet_Theta_R0p1->Write();
	h_HardCore_RecoilJet_Twosubjet_Theta_R0p1->Write();
	h_HardCore_TriggerJet_deltaR_JetAxis_LeadSubJetAxis->Write();
	h_HardCore_RecoilJet_deltaR_JetAxis_LeadSubJetAxis->Write();
	h_Matched_Dijet_Aj->Write();
	h_Matched_RecoilJetYield->Write();
	hDeltaPhi_Matched_Dijet->Write();
	h_Matched_TriggerJet_Zg->Write();
	h_Matched_RecoilJet_Zg->Write();
	h_Matched_TriggerJet_Rg->Write();
	h_Matched_RecoilJet_Rg->Write();
	h_Matched_TriggerJet_Twosubjet_Z_R0p1->Write();
	h_Matched_RecoilJet_Twosubjet_Z_R0p1->Write();
	h_Matched_TriggerJet_Twosubjet_Theta_R0p1->Write();
	h_Matched_RecoilJet_Twosubjet_Theta_R0p1->Write();
	h_Matched_TriggerJet_deltaR_JetAxis_LeadSubJetAxis->Write();
	h_Matched_RecoilJet_deltaR_JetAxis_LeadSubJetAxis->Write();

	hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis->Write();
	hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis->Write();
    
	hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis_largeR->Write();
	hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis_largeR->Write();
	
      }
      
      fout->Close();
      
    }//! finalize


  protected:

    size_t _mode;
    
    double _jetR;
    int RADIUS;
    double _pTCut;

    double z_cut;
    double beta;

    double _Ntriggers;
    double _Ndijets;
        
  private:

    TH1D * hDeltaPhi_Trigger_ChargedParticles;
    TH1D * hDeltaPhi_Trigger_LeadingRecoilJet;
    TH1D * hDeltaPhi_Trigger_SubLeadingRecoilJet;
    TH1D * hDeltaPhi_Trigger_AllRecoilJets;
    TH1D * hDeltaR_RecoilJet_JetAxis_WTAAxis;
    TH1D * hDeltaR_RecoilJet_JetAxis_WTAAxis_largeR;
    TH1D * hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis;
    TH1D * h_RecoilJet_TwoSubJet_Theta_R0p1;
    TH1D * h_RecoilJet_TwoSubJet_Z_R0p1;
    TH1D * hSoftDrop_Zg_LeadingRecoilJet;
    TH1D * hSoftDrop_Rg_LeadingRecoilJet;
    TH1D * hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis;
    TH1D * hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis_largeR;
    TH1D * hDeltaR_RecoilJet_GroomedJetAxis_WTAAxis_largeR;
    TH1D * hTriggerPhi;
    TH1D * hJetPhi;

    TH1D * h_HardCore_Dijet_Aj;
    TH1D * h_HardCore_RecoilJetYield;
    TH1D * hDeltaPhi_HardCore_Dijet;
    TH1D * h_HardCore_TriggerJet_Zg;
    TH1D * h_HardCore_RecoilJet_Zg;
    TH1D * h_HardCore_TriggerJet_Rg;
    TH1D * h_HardCore_RecoilJet_Rg;
    TH1D * h_HardCore_TriggerJet_Twosubjet_Z_R0p1;
    TH1D * h_HardCore_RecoilJet_Twosubjet_Z_R0p1;
    TH1D * h_HardCore_TriggerJet_Twosubjet_Theta_R0p1;
    TH1D * h_HardCore_RecoilJet_Twosubjet_Theta_R0p1;
    TH1D * h_HardCore_TriggerJet_deltaR_JetAxis_LeadSubJetAxis;
    TH1D * h_HardCore_RecoilJet_deltaR_JetAxis_LeadSubJetAxis;
    TH1D * h_Matched_Dijet_Aj;
    TH1D * h_Matched_RecoilJetYield;
    TH1D * hDeltaPhi_Matched_Dijet;
    TH1D * h_Matched_TriggerJet_Zg;
    TH1D * h_Matched_RecoilJet_Zg;
    TH1D * h_Matched_TriggerJet_Rg;
    TH1D * h_Matched_RecoilJet_Rg;
    TH1D * h_Matched_TriggerJet_Twosubjet_Z_R0p1;
    TH1D * h_Matched_RecoilJet_Twosubjet_Z_R0p1;
    TH1D * h_Matched_TriggerJet_Twosubjet_Theta_R0p1;
    TH1D * h_Matched_RecoilJet_Twosubjet_Theta_R0p1;
    TH1D * h_Matched_TriggerJet_deltaR_JetAxis_LeadSubJetAxis;
    TH1D * h_Matched_RecoilJet_deltaR_JetAxis_LeadSubJetAxis;

    TH1D * hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis;
    TH1D * hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis;
    
    TH1D * hDeltaR_MatchedRecoilJet_JetAxis_WTAAxis_largeR;
    TH1D * hDeltaPhi_MatchedRecoilJet_MatchedRecoilJet_WTAAxis_largeR;
    
    
    TFile * fout;
    
  };

  //! photon jet coincidence 
  class JETTOOLS_PHOTON_JET_COINCIDENCE : public JETTOOLS_COINCIDENCE {
  public:
    JETTOOLS_PHOTON_JET_COINCIDENCE()
      : JETTOOLS_COINCIDENCE("JETTOOLS_PHOTON_JET_COINCIDENCE")
    {
      _mode = 1;
    }
  };
  
  //! dijet coincidence 
  class JETTOOLS_DIJET_COINCIDENCE : public JETTOOLS_COINCIDENCE {
  public:
    JETTOOLS_DIJET_COINCIDENCE()
      : JETTOOLS_COINCIDENCE("JETTOOLS_DIJET_COINCIDENCE")
    {
      _mode = 2;
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JETTOOLS_COINCIDENCE);
  DECLARE_RIVET_PLUGIN(JETTOOLS_PHOTON_JET_COINCIDENCE);
  DECLARE_RIVET_PLUGIN(JETTOOLS_DIJET_COINCIDENCE);
  
  
}

