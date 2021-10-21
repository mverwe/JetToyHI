#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"

#include "include/jetCollection.hh"
#include "include/softDropGroomer.hh"
#include "include/softDropCounter.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/csSubtractor.hh"
#include "include/csSubtractorFullEvent.hh"

using namespace std;
using namespace fastjet;

// This class runs time reclustering with background

// ./runSubjetFF -hard samples/PythiaEventsTune14PtHat120_10k.pu14 -pileup samples/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10

int main (int argc, char ** argv) {

  auto start_time = std::chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  std::cout << "will run on " << nEvent << " events" << std::endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trwSig("jetTreeSig");
 
  //Jet definition
  double R                   = 0.5;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 100.;//3.0;
  fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && iev < nEvent )
  {
    // increment event number    
    iev++;

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    std::vector<fastjet::PseudoJet> particlesMergedAll = mixer.particles();

    std::vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // extract hard partons that initiated the jets
    fastjet::Selector parton_selector = SelectorVertexNumber(-1);
    vector<PseudoJet> partons = parton_selector(particlesMergedAll);

    // extract hard partons from first splitting
    fastjet::Selector parton_selector_split = SelectorVertexNumber(-2);
    vector<PseudoJet> partonsFirstSplit = parton_selector_split(particlesMergedAll);

    // select final state particles from hard event only
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    // select final state particles from background event only
    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );

    //charged particles
    fastjet::Selector charged_selector = SelectorIsCharged();
    vector<PseudoJet> particlesSigCh = charged_selector(particlesSig);

    //std::cout << "#particles: " << particlesSig.size() << " of which charged: " << particlesSigCh.size() << std::endl;

 
    //---------------------------------------------------------------------------
    //   look at first splitting of hard partons
    //---------------------------------------------------------------------------
    std::vector<double> drsplit;
    std::vector<double> tfsplit;
    double hbarc = 0.19732697;
    double GeVtofm = 1./hbarc; //~5.068;
    int id = 0;
    for(int ip = 0; ip<partons.size(); ++ip) {
      PseudoJet p = partons[ip];
      PseudoJet d1 = partonsFirstSplit[id];
      PseudoJet d2 = partonsFirstSplit[id+1];
      double dr = d1.delta_R(d2);
      drsplit.push_back(dr);
      double z1 = max(d1.e(),d2.e())/p.e();
      double z2 = min(d1.e(),d2.e())/p.e();
      tfsplit.push_back(1./(2.*z1*z2*p.e()*GeVtofm*(1-fastjet::cos_theta(d1,d2))));
    
      id+=2;
    }
    
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    // run the clustering, extract the signal charged jets
    fastjet::ClusterSequenceArea csSigCh(particlesSigCh, jet_def, area_def);
    jetCollection jetCollectionSigCh(sorted_by_pt(jet_selector(csSigCh.inclusive_jets(10.))));
    
    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for signal jets
    //---------------------------------------------------------------------------
    
    softDropCounter sdcSig(0.0,0.0,R,0.0);
    sdcSig.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSig.run(jetCollectionSig);

    jetCollectionSig.addVector("sigJetRecur_jetpt",     sdcSig.getPts());
    jetCollectionSig.addVector("sigJetRecur_z",         sdcSig.getZgs());
    jetCollectionSig.addVector("sigJetRecur_dr12",      sdcSig.getDRs());
    jetCollectionSig.addVector("sigJetRecur_erad",      sdcSig.getErads());
    jetCollectionSig.addVector("sigJetRecur_logdr12",   sdcSig.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecur_logztheta", sdcSig.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecur_tf",        sdcSig.getTfs());
    jetCollectionSig.addVector("sigJetRecur_tfe",       sdcSig.getTfes());
    jetCollectionSig.addVector("sigJetRecur_nSD",       sdcSig.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecur_zSD",       sdcSig.calculateNSD(1.0));
    jetCollectionSig.addVector("sigJetRecur_kt",         sdcSig.getKts());
    
    //find closest parton for each jet
    std::vector<int> partonmatch;
    std::vector<double> partonmatchdr;
    std::vector<fastjet::PseudoJet> sigJets =  jetCollectionSig.getJet();
    for(fastjet::PseudoJet p : sigJets) {
      int ipmin = -1;
      double drmin = 999.;
      for(int ip = 0; ip<partons.size(); ++ip) {
        double dr = p.delta_R(partons[ip]);
        if(dr<drmin) {
          drmin = dr;
          ipmin = ip;
        }
      }
      partonmatch.push_back(ipmin);
      partonmatchdr.push_back(drmin);
    }
    jetCollectionSig.addVector("sigJetRecur_partonMatchID", partonmatch);
    jetCollectionSig.addVector("sigJetRecur_partonMatchDr", partonmatchdr);

    //---------------------------------------------------------------------------
    //   jet clustering of charged-particle signal jets
    //---------------------------------------------------------------------------
    softDropCounter sdcSigCh(0.,0.0,R,0.0);
    sdcSigCh.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSigCh.run(jetCollectionSigCh);

    
    jetCollectionSigCh.addVector("sigJetChRecur_jetpt",     sdcSigCh.getPts());
    jetCollectionSigCh.addVector("sigJetChRecur_z",         sdcSigCh.getZgs());
    jetCollectionSigCh.addVector("sigJetChRecur_dr12",      sdcSigCh.getDRs());
    jetCollectionSigCh.addVector("sigJetChRecur_erad",      sdcSigCh.getErads());
    jetCollectionSigCh.addVector("sigJetChRecur_logdr12",   sdcSigCh.getLog1DRs());
    jetCollectionSigCh.addVector("sigJetChRecur_logztheta", sdcSigCh.getLogzDRs());
    jetCollectionSigCh.addVector("sigJetChRecur_tf",        sdcSigCh.getTfs());
    jetCollectionSigCh.addVector("sigJetChRecur_tfe",       sdcSigCh.getTfes());
    jetCollectionSigCh.addVector("sigJetChRecur_nSD",       sdcSigCh.calculateNSD(0.0));
    jetCollectionSigCh.addVector("sigJetChRecur_zSD",       sdcSigCh.calculateNSD(1.0));
    jetCollectionSigCh.addVector("sigJetChRecur_kt",         sdcSigCh.getKts());

    //find closest parton for each charged jet
    std::vector<int> partonmatchCh;
    std::vector<double> partonmatchdrCh;
    std::vector<fastjet::PseudoJet> sigJetsCh =  jetCollectionSigCh.getJet();
    for(fastjet::PseudoJet p : sigJetsCh) {
      int ipmin = -1;
      double drmin = 999.;
      for(int ip = 0; ip<partons.size(); ++ip) {
        double dr = p.delta_R(partons[ip]);
        if(dr<drmin) {
          drmin = dr;
          ipmin = ip;
        }
      }
      partonmatchCh.push_back(ipmin);
      partonmatchdrCh.push_back(drmin);
    }
    jetCollectionSigCh.addVector("sigJetChRecur_partonMatchID", partonmatchCh);
    jetCollectionSigCh.addVector("sigJetChRecur_partonMatchDr", partonmatchdrCh);
    

    //---------------------------------------------------------------------------
    //   jet clustering of signal+background jets
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csRaw(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionRaw(sorted_by_pt(jet_selector(csRaw.inclusive_jets(10.))));

    //---------------------------------------------------------------------------
    //   background subtraction
    //---------------------------------------------------------------------------

    //run jet-by-jet constituent subtraction on mixed (hard+UE) event
    csSubtractor csSub(R, 0., -1, 0.005,ghostRapMax,2.5);
    csSub.setInputParticles(particlesMerged);
    jetCollection jetCollectionCS(csSub.doSubtraction());

    //Background densities used by constituent subtraction
    std::vector<double> rho;
    std::vector<double> rhom;
    rho.push_back(csSub.getRho());
    rhom.push_back(csSub.getRhoM());

    //match CS jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionCS);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();

    jmCS.reorderedToTag(jetCollectionCS);

    //---------------------------------------------------------------------------
    //  recluster subtracted jets
    //---------------------------------------------------------------------------
    softDropCounter sdcCS(0.0,0.0,R,0.0);
    sdcCS.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcCS.run(jetCollectionCS);

    jetCollectionCS.addVector("csJetRecur_jetpt",     sdcCS.getPts());
    jetCollectionCS.addVector("csJetRecur_z",         sdcCS.getZgs());
    jetCollectionCS.addVector("csJetRecur_dr12",      sdcCS.getDRs());
    jetCollectionCS.addVector("csJetRecur_erad",      sdcCS.getErads());
    jetCollectionCS.addVector("csJetRecur_logdr12",   sdcCS.getLog1DRs());
    jetCollectionCS.addVector("csJetRecur_logztheta", sdcCS.getLogzDRs());
    jetCollectionCS.addVector("csJetRecur_tf",        sdcCS.getTfs());
    jetCollectionCS.addVector("csJetRecur_tfe",       sdcCS.getTfes());
    jetCollectionCS.addVector("csJetRecur_nSD",       sdcCS.calculateNSD(0.0));
    jetCollectionCS.addVector("csJetRecur_zSD",       sdcCS.calculateNSD(1.0));
    jetCollectionCS.addVector("csJetRecur_kt",         sdcCS.getKts());

    //find closest parton for each jet
    std::vector<int> partonmatchCS;
    std::vector<double> partonmatchdrCS;
    std::vector<fastjet::PseudoJet> csJets =  jetCollectionCS.getJet();
    for(fastjet::PseudoJet p : csJets) {
      int ipmin = -1;
      double drmin = 999.;
      for(int ip = 0; ip<partons.size(); ++ip) {
        double dr = p.delta_R(partons[ip]);
        if(dr<drmin) {
          drmin = dr;
          ipmin = ip;
        }
      }
      partonmatchCS.push_back(ipmin);
      partonmatchdrCS.push_back(drmin);
    }
    jetCollectionCS.addVector("csJetRecur_partonMatchID", partonmatchCS);
    jetCollectionCS.addVector("csJetRecur_partonMatchDr", partonmatchdrCS);

    //---------------------------------------------------------------------------
    //   Full event constituent subtraction
    //---------------------------------------------------------------------------
    csSubtractorFullEvent csSubFull( 0., 0.25, 0.005, 2.5);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setRho(csSub.getRho());
    csSubFull.setRhom(csSub.getRhoM());
    csSubFull.setInputParticles(particlesMerged);
    
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtraction(), jet_def, area_def);
    jetCollection jetCollectionCSFull(sorted_by_pt(jet_selector(fullSig.inclusive_jets(10.)))); 

    //match CSFull jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCSFull);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();

    jmCSFull.reorderedToTag(jetCollectionCSFull);

    //---------------------------------------------------------------------------
    //  recluster subtracted (full event CS) jets
    //---------------------------------------------------------------------------
    softDropCounter sdcCSFull(0.0,0.0,R,0.0);
    sdcCSFull.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcCSFull.run(jetCollectionCSFull);

    jetCollectionCSFull.addVector("csFullJetRecur_jetpt",     sdcCSFull.getPts());
    jetCollectionCSFull.addVector("csFullJetRecur_z",         sdcCSFull.getZgs());
    jetCollectionCSFull.addVector("csFullJetRecur_dr12",      sdcCSFull.getDRs());
    jetCollectionCSFull.addVector("csFullJetRecur_erad",      sdcCSFull.getErads());
    jetCollectionCSFull.addVector("csFullJetRecur_logdr12",   sdcCSFull.getLog1DRs());
    jetCollectionCSFull.addVector("csFullJetRecur_logztheta", sdcCSFull.getLogzDRs());
    jetCollectionCSFull.addVector("csFullJetRecur_tf",        sdcCSFull.getTfs());
    jetCollectionCSFull.addVector("csFullJetRecur_tfe",       sdcCSFull.getTfes());
    jetCollectionCSFull.addVector("csFullJetRecur_nSD",       sdcCSFull.calculateNSD(0.0));
    jetCollectionCSFull.addVector("csFullJetRecur_zSD",       sdcCSFull.calculateNSD(1.0));
    jetCollectionCSFull.addVector("csFullJetRecur_kt",         sdcCSFull.getKts());

    //find closest parton for each jet
    std::vector<int> partonmatchCSFull;
    std::vector<double> partonmatchdrCSFull;
    std::vector<fastjet::PseudoJet> csFullJets =  jetCollectionCSFull.getJet();
    for(fastjet::PseudoJet p : csFullJets) {
      int ipmin = -1;
      double drmin = 999.;
      for(int ip = 0; ip<partons.size(); ++ip) {
        double dr = p.delta_R(partons[ip]);
        if(dr<drmin) {
          drmin = dr;
          ipmin = ip;
        }
      }
      partonmatchCSFull.push_back(ipmin);
      partonmatchdrCSFull.push_back(drmin);
    }
    jetCollectionCSFull.addVector("csFullJetRecur_partonMatchID", partonmatchCSFull);
    jetCollectionCSFull.addVector("csFullJetRecur_partonMatchDr", partonmatchdrCSFull);

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("eventWeight",   eventWeight);
    trwSig.addCollection("csRho",         rho);
    trwSig.addCollection("csRhom",        rhom);

    trwSig.addPartonCollection("partons",       partons);
    trwSig.addPartonCollection("partonsFirstSplit",       partonsFirstSplit);
    trwSig.addDoubleCollection("drSplit", drsplit);
    trwSig.addDoubleCollection("tfSplit", tfsplit);
    
    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetCh",      jetCollectionSigCh);

    trwSig.addCollection("csJet",         jetCollectionCS);
    trwSig.addCollection("csFullJet",         jetCollectionCSFull);
    
    trwSig.fillTree();  //signal jets
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile("JetToyHIResultSubjetFF.root","RECREATE");
  trwSig.getTree()->Write();
  
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
