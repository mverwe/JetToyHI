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
#include "include/csSubtractor.hh"
#include "include/csSubtractorFullEvent.hh"
#include "include/skSubtractor.hh"
#include "include/softDropGroomer.hh"
#include "include/softDropCounter.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/randomCones.hh"
#include "include/Angularity.hh"
#include "include/jetProfile.hh"

using namespace std;
using namespace fastjet;

// ./runJetProfile -hard PythiaEventsTune14PtHat120.pu14  -nev 1

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
  treeWriter trwEmb("jetTreeEmb");

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 3.0;
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

    std::vector<fastjet::PseudoJet> particlesMerged = mixer.particles();

    std::vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // cluster hard event only
    std::vector<fastjet::PseudoJet> particlesBkg, particlesSig;
    SelectorIsHard().sift(particlesMerged, particlesSig, particlesBkg); // this sifts the full event into two vectors of PseudoJet, one for the hard event, one for the underlying event
    
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    //---------------------------------------------------------------------------
    //   background subtraction
    //---------------------------------------------------------------------------
    
    //run jet-by-jet constituent subtraction on mixed (hard+UE) event
    std::vector<jetCollection> jetCollectionCSs;
    std::vector<double> rho;
    std::vector<double> rhom; 
    csSubtractor csSub(R, 1., -1, 0.005,ghostRapMax,jetRapMax);
    csSub.setInputParticles(particlesMerged);
    jetCollection jetCollectionCS(csSub.doSubtraction());
    //Background densities used by constituent subtraction
    rho.push_back(csSub.getRho());
    rhom.push_back(csSub.getRhoM());

    //---------------------------------------------------------------------------
    //   Calculate jet radial profile for signal jets
    //---------------------------------------------------------------------------
    jetProfile jetProf(jetCollectionSig, 1.);
    const int np = 8;
    double min_delR[np] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35};
    double max_delR[np] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4};
    std::vector<double> vmin;
    std::vector<double> vmax;
    for(int j = 0; j<np; ++j) {
      vmin.push_back(min_delR[j]);
      vmax.push_back(max_delR[j]);
    }
    jetProf.setBoundariesMin(vmin);
    jetProf.setBoundariesMax(vmax);
    jetProf.calculateProfile();
    jetCollectionSig.addVector(Form("sigJetProfile"), jetProf.getJetProfiles());


    //---------------------------------------------------------------------------
    //   Soft Drop for the signal jets
    //---------------------------------------------------------------------------
    
    //Using soft drop grooming class without actually doing grooming, kT ordering
    softDropGroomer sdgSig(0.0, 0.0, R);
    sdgSig.setReclusteringAlgo(2);//0 = CA 1 = AKT 2 = KT
 
    std::vector<fastjet::PseudoJet> groomedJets_Sig = sdgSig.doGrooming(jetCollectionSig);
    jetCollection jetCollectionSigSD(groomedJets_Sig);
    jetCollectionSigSD.addVector("sigJetSDzg",    sdgSig.getZgs());
    jetCollectionSigSD.addVector("sigJetSDndrop", sdgSig.getNDroppedSubjets());
    jetCollectionSigSD.addVector("sigJetSDdr12",  sdgSig.getDR12());
    jetCollectionSigSD.addVector("sigJetSDlogdr12",  sdgSig.getLogDR12());
    jetCollectionSigSD.addVector("sigJetSDlogztheta",  sdgSig.getLogZgDR12());
    jetCollectionSigSD.addVector("sigJetSDmass",  sdgSig.getSubJetMass());
    jetCollectionSigSD.addVector("sigJetSDleadingtrack_pt",  sdgSig.getSubJetLeadingTrackPt());

    //Using soft drop grooming class with classical grooming (zcut=0.1, beta=0.0), CA ordering
    softDropGroomer sdgSigZ01B00(0.1, 0.0, R);
    sdgSigZ01B00.setReclusteringAlgo(0);//0 = CA 1 = AKT 2 = KT
 
    std::vector<fastjet::PseudoJet> groomedJets_SigZ01B00 = sdgSigZ01B00.doGrooming(jetCollectionSig);
    jetCollection jetCollectionSigSDZ01B00(groomedJets_SigZ01B00);
    jetCollectionSigSDZ01B00.addVector("sigJetSDZ01B00zg",    sdgSigZ01B00.getZgs());
    jetCollectionSigSDZ01B00.addVector("sigJetSDZ01B00ndrop", sdgSigZ01B00.getNDroppedSubjets());
    jetCollectionSigSDZ01B00.addVector("sigJetSDZ01B00dr12",  sdgSigZ01B00.getDR12());
    jetCollectionSigSDZ01B00.addVector("sigJetSDZ01B00logdr12",  sdgSigZ01B00.getLogDR12());
    jetCollectionSigSDZ01B00.addVector("sigJetSDZ01B00logztheta",  sdgSigZ01B00.getLogZgDR12());
    jetCollectionSigSDZ01B00.addVector("sigJetSDZ01B00mass",  sdgSigZ01B00.getSubJetMass());
    jetCollectionSigSDZ01B00.addVector("sigJetSDZ01B00leadingtrack_pt",  sdgSigZ01B00.getSubJetLeadingTrackPt());

    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for signal jets
    //---------------------------------------------------------------------------
    
    softDropCounter sdcSig(0.0,0.0,R,0.0);
    sdcSig.setRecursiveAlgo(2);//0 = CA 1 = AKT 2 = KT
    sdcSig.run(jetCollectionSig);

    jetCollectionSigSD.addVector("sigJetRecur_jetpt",     sdcSig.getPts());
    jetCollectionSigSD.addVector("sigJetRecur_z",         sdcSig.getZgs());
    jetCollectionSigSD.addVector("sigJetRecur_dr12",      sdcSig.getDRs());
    jetCollectionSigSD.addVector("sigJetRecur_erad",      sdcSig.getErads());
    jetCollectionSigSD.addVector("sigJetRecur_logdr12",   sdcSig.getLog1DRs());
    jetCollectionSigSD.addVector("sigJetRecur_logztheta", sdcSig.getLogzDRs());


    //---------------------------------------------------------------------------
    //   Calculate jet radial profile for embedded jets
    //---------------------------------------------------------------------------
    jetProfile jetProfCS(jetCollectionCS, 1.);
    jetProfCS.setBoundariesMin(vmin);
    jetProfCS.setBoundariesMax(vmax);
    jetProfCS.calculateProfile();
    jetCollectionCS.addVector(Form("csJetProfile"), jetProfCS.getJetProfiles());

    //---------------------------------------------------------------------------
    //   Soft Drop for the embedded CS jets
    //---------------------------------------------------------------------------
    
    //Using soft drop grooming class without actually doing grooming
    softDropGroomer sdgCS(0.0, 0.0, R);
    sdgCS.setReclusteringAlgo(2);//0 = CA 1 = AKT 2 = KT
 
    std::vector<fastjet::PseudoJet> groomedJets_CS = sdgCS.doGrooming(jetCollectionCS);
    jetCollection jetCollectionCSSD(groomedJets_CS);
    jetCollectionCSSD.addVector("csJetSDzg",    sdgCS.getZgs());
    jetCollectionCSSD.addVector("csJetSDndrop", sdgCS.getNDroppedSubjets());
    jetCollectionCSSD.addVector("csJetSDdr12",  sdgCS.getDR12());
    jetCollectionCSSD.addVector("csJetSDlogdr12",  sdgCS.getLogDR12());
    jetCollectionCSSD.addVector("csJetSDlogztheta",  sdgCS.getLogZgDR12());
    jetCollectionCSSD.addVector("csJetSDmass",  sdgCS.getSubJetMass());
    jetCollectionCSSD.addVector("csJetSDleadingtrack_pt",  sdgCS.getSubJetLeadingTrackPt());

    //Using soft drop grooming class with classical grooming (zcut=0.1, beta=0.0), CA ordering
    softDropGroomer sdgCSZ01B00(0.1, 0.0, R);
    sdgCSZ01B00.setReclusteringAlgo(0);//0 = CA 1 = AKT 2 = KT
 
    std::vector<fastjet::PseudoJet> groomedJets_CSZ01B00 = sdgCSZ01B00.doGrooming(jetCollectionCS);
    jetCollection jetCollectionCSSDZ01B00(groomedJets_CSZ01B00);
    jetCollectionCSSDZ01B00.addVector("csJetSDZ01B00zg",    sdgCSZ01B00.getZgs());
    jetCollectionCSSDZ01B00.addVector("csJetSDZ01B00ndrop", sdgCSZ01B00.getNDroppedSubjets());
    jetCollectionCSSDZ01B00.addVector("csJetSDZ01B00dr12",  sdgCSZ01B00.getDR12());
    jetCollectionCSSDZ01B00.addVector("csJetSDZ01B00logdr12",  sdgCSZ01B00.getLogDR12());
    jetCollectionCSSDZ01B00.addVector("csJetSDZ01B00logztheta",  sdgCSZ01B00.getLogZgDR12());
    jetCollectionCSSDZ01B00.addVector("csJetSDZ01B00mass",  sdgCSZ01B00.getSubJetMass());
    jetCollectionCSSDZ01B00.addVector("csJetSDZ01B00leadingtrack_pt",  sdgCSZ01B00.getSubJetLeadingTrackPt());


    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for embedded CS jets
    //---------------------------------------------------------------------------
    
    softDropCounter sdcCS(0.0,0.0,R,0.0);
    sdcCS.setRecursiveAlgo(2);//0 = CA 1 = AKT 2 = KT
    sdcCS.run(jetCollectionCS);

    jetCollectionCSSD.addVector("csJetRecur_jetpt",     sdcCS.getPts());
    jetCollectionCSSD.addVector("csJetRecur_z",         sdcCS.getZgs());
    jetCollectionCSSD.addVector("csJetRecur_dr12",      sdcCS.getDRs());
    jetCollectionCSSD.addVector("csJetRecur_erad",      sdcCS.getErads());
    jetCollectionCSSD.addVector("csJetRecur_logdr12",   sdcCS.getLog1DRs());
    jetCollectionCSSD.addVector("csJetRecur_logztheta", sdcCS.getLogzDRs());
  
        
    /* Question from MV: do uoi want matching between embedded and signal jets here?
    //match the CS jets to signal jets
    for(int ics = 0; ics<ncs; ++ics) {
      jetMatcher jmCS(R);
      jmCS.setBaseJets(jetCollectionCSs[ics]);
      jmCS.setTagJets(jetCollectionSig);
      jmCS.matchJets();

      jmCS.reorderedToTag(jetCollectionCSs[ics]);
      jmCS.reorderedToTag(jetCollectionCSSDs[ics]);
    }
  */

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("eventWeight",   eventWeight);
    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetSD",      jetCollectionSigSD);
    trwSig.addCollection("sigJetSDZ01B00",jetCollectionSigSDZ01B00);

    trwEmb.addCollection("csRho",         rho);
    trwEmb.addCollection("csRhom",        rhom);
    trwEmb.addCollection("eventWeight",   eventWeight);
    trwEmb.addCollection("csJet",         jetCollectionCS);
    trwEmb.addCollection("csJetSD",       jetCollectionCSSD);
    trwEmb.addCollection("csJetSDZ01B00", jetCollectionCSSDZ01B00);
        

    trwSig.fillTree();  //signal jets
    trwEmb.fillTree();  //embedded jets
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile("JetToyHIResultJetProfile.root","RECREATE");
  trwSig.getTree()->Write();
  trwEmb.getTree()->Write();

  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
