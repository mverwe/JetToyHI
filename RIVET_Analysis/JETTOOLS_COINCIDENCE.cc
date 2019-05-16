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

      /// final state and jet projection
      FinalState fs(-5.0, 5.0, 0.2*GeV);
      addProjection(fs, "FS");
      FinalState HCfs(-5.0, 5.0, 2.0*GeV);
      addProjection(HCfs, "HCFS");
      VetoedFinalState vfs(fs);
      ChargedFinalState cfs(-5.0, 5.0, 0.2*GeV);
      addProjection(cfs, "CFS");

      if (_mode == 0)
	{
	  //! this is just hadron+jet
	  LeadingParticlesFinalState hadronfs(FinalState(-1.0, 1.0, 10.0*GeV));
	  hadronfs.addParticleId(111);
	  hadronfs.addParticleId(211);
	  hadronfs.addParticleId(-211);
	  addProjection(hadronfs, "LeadingHadron");
	  vfs.addVetoOnThisFinalState(hadronfs);
	  addProjection(vfs, "VFS");
	  // addProjection(FastJets(vfs, FastJets::ANTIKT, _jetR), "Jets");
	}
      else if (_mode == 1)
	{
	  //! this is photon+jet 
	  LeadingParticlesFinalState photonfs(FinalState(-2.5, 2.5, 10.0*GeV));
	  photonfs.addParticleId(PID::PHOTON);
	  addProjection(photonfs, "LeadingPhoton");
	  vfs.addVetoOnThisFinalState(photonfs);
	  addProjection(vfs, "VFS");
	  // addProjection(FastJets(vfs, FastJets::ANTIKT, _jetR), "Jets");
	}
      else if (_mode == 2)
	{
	  //! this is dijet+jet
	  // addProjection(FastJets(HCfs, FastJets::ANTIKT, _jetR), "JetsHC");	  
	  // addProjection(FastJets(fs, FastJets::ANTIKT, _jetR), "Jets");	  
	}

      fout = new TFile("jettools_jewel_pp_200_dijet_gammaJet_study_correlation_histograms.root", "RECREATE");
      
      //! declare all the histograms that you want to study -
      hDeltaPhi_Trigger_ChargedParticles = new TH1D("hDeltaPhi_Trigger_ChargedParticles", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_LeadingRecoilJet = new TH1D("hDeltaPhi_Trigger_LeadingRecoilJet", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis = new TH1D("hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_SubLeadingRecoilJet = new TH1D("hDeltaPhi_Trigger_SubLeadingRecoilJet", "", 60, 0, 3.3);
      hDeltaPhi_Trigger_AllRecoilJets = new TH1D("hDeltaPhi_Trigger_AllRecoilJets", "", 60, 0, 3.3);
      hDeltaR_RecoilJet_JetAxis_WTAAxis = new TH1D("hDeltaR_RecoilJet_JetAxis_WTAAxis", "", 60, 0, 0.4);
      hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis = new TH1D("hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis", "", 60, 0, 0.4);
      h_RecoilJet_TwoSubJet_Theta_R0p1 = new TH1D("h_RecoilJet_TwoSubJet_Theta_R0p1", "", 60, 0, 0.5);
      h_RecoilJet_TwoSubJet_Z_R0p1 = new TH1D("h_RecoilJet_TwoSubJet_Z_R0p1", "", 60, 0, 0.5);
      hSoftDrop_Zg_LeadingRecoilJet = new TH1D("hSoftDrop_Zg_LeadingRecoilJet", "", 60, 0, 0.5);
      hSoftDrop_Rg_LeadingRecoilJet = new TH1D("hSoftDrop_Rg_LeadingRecoilJet", "", 60, 0, 0.5);

      hTriggerPhi = new TH1D("hTriggerPhi", "", 100, -7, 7);
      hJetPhi = new TH1D("hJetPhi", "", 100, -7, 7);
      

    }//! init method

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      ParticleVector Trigger;
      const FinalState& chargedfs = applyProjection<FinalState>(event, "CFS");      
      
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

      double trigger_phi = trigger.momentum().phi();
      hTriggerPhi->Fill(trigger_phi);

      //! Now do the analysis for both triggered samples -
      if(_mode == 1 || _mode == 0){

	fastjet::Selector select_eta  = fastjet::SelectorAbsEtaMax(1.0 - _jetR);
	fastjet::Selector select_pt_lo   = fastjet::SelectorPtMin(_pTCut);
	fastjet::Selector select_both = select_pt_lo && select_eta;

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
	fastjet::JetDefinition jetd(fastjet::antikt_algorithm, 0.1);
	fastjet::ClusterSequence cs_recoiljets(vetoedparticles, jetd);
	PseudoJets recoilJETS = sorted_by_pt(cs_recoiljets.inclusive_jets());
	PseudoJets recoilJets = select_both(recoilJETS); 
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

	//! get wta axis for the jet -
	JetDefinition wta_jet_def(fastjet::antikt_algorithm, _jetR, fastjet::WTA_pt_scheme);
	
	PseudoJets lead_recoiljet_consts = leading_recoiljet.constituents();
	fastjet::ClusterSequence cs_wta_leadingrecoiljet(lead_recoiljet_consts, wta_jet_def);
	PseudoJets wta_lead_recoiljets = cs_wta_leadingrecoiljet.inclusive_jets();
	PseudoJet wta_lead_recoiljet = wta_lead_recoiljets.at(0);

	hDeltaR_RecoilJet_JetAxis_WTAAxis->Fill(leading_recoiljet.delta_R(wta_lead_recoiljet), weight);
	hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis->Fill(deltaPhi(trigger.momentum().phi(), wta_lead_recoiljet.phi()), weight);

	//! get subjets -
	fastjet::JetDefinition sub_jetd(fastjet::antikt_algorithm, 0.1);
	fastjet::ClusterSequence clust_seq_lead_recoil_subjet(lead_recoiljet_consts, sub_jetd);\	
	PseudoJets subjets = sorted_by_pt(clust_seq_lead_recoil_subjet.inclusive_jets());

	hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis->Fill(leading_recoiljet.delta_R(subjets.at(0)), weight);

	if(subjets.size() >= 2){
	  //! three subjets - no grooming. just straight up two subjets 
	  double zsj = subjets.at(1).pt()/(subjets.at(0).pt()+subjets.at(1).pt());
	  double rsj = subjets.at(0).delta_R(subjets.at(1));
	  h_RecoilJet_TwoSubJet_Theta_R0p1->Fill(rsj, weight);
	  h_RecoilJet_TwoSubJet_Z_R0p1->Fill(zsj, weight);
	}

	//! fill the softdrop variable
	fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
	fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, _jetR);
	PseudoJet sd_jet = sd(leading_recoiljet);
	double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();

	hSoftDrop_Zg_LeadingRecoilJet->Fill(z, weight);
	hSoftDrop_Rg_LeadingRecoilJet->Fill(r, weight);


      }//! mode = 1 or 0 
      
      
      if (_mode == 2)
	{
	  // //! dijet selections 
	  // const FastJets& HCJets = applyProjection<FastJets>(event, "JetsHC");
	  // Cut cuts = Cuts::etaIn(-1.0, 1.0) & (Cuts::pT > 10.0*GeV);
	  // const Jets hcJets = HCJets.jetsByPt(cuts);
	  // const int _NjetHC = hcJets.size();
	  // if(_NjetHC < 1) vetoEvent;
	}
      

    }//! analyze method 

    /// Normalise histograms etc., after the run
    void finalize() {

      fout->cd();
      fout->Write();
      fout->Close();
      
    }//! finalize


  protected:

    size_t _mode;
    
    double _jetR;
    int RADIUS;
    double _pTCut;

    double z_cut;
    double beta;
        
  private:

    TH1D * hDeltaPhi_Trigger_ChargedParticles;
    TH1D * hDeltaPhi_Trigger_LeadingRecoilJet;
    TH1D * hDeltaPhi_Trigger_SubLeadingRecoilJet;
    TH1D * hDeltaPhi_Trigger_AllRecoilJets;
    TH1D * hDeltaR_RecoilJet_JetAxis_WTAAxis;
    TH1D * hDeltaR_RecoilJet_JetAxis_LeadSubJetAxis;
    TH1D * h_RecoilJet_TwoSubJet_Theta_R0p1;
    TH1D * h_RecoilJet_TwoSubJet_Z_R0p1;
    TH1D * hSoftDrop_Zg_LeadingRecoilJet;
    TH1D * hSoftDrop_Rg_LeadingRecoilJet;
    TH1D * hDeltaPhi_Trigger_LeadingRecoilJet_WTAAxis;
    TH1D * hTriggerPhi;
    TH1D * hJetPhi;

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
  class JETTOOLS_JET_JET_COINCIDENCE : public JETTOOLS_COINCIDENCE {
  public:
    JETTOOLS_JET_JET_COINCIDENCE()
      : JETTOOLS_COINCIDENCE("JETTOOLS_JET_JET_COINCIDENCE")
    {
      _mode = 2;
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JETTOOLS_COINCIDENCE);
  DECLARE_RIVET_PLUGIN(JETTOOLS_PHOTON_JET_COINCIDENCE);
  DECLARE_RIVET_PLUGIN(JETTOOLS_JET_JET_COINCIDENCE);
  
  
}

