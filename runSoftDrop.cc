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

using namespace std;
using namespace fastjet;

// This class runs soft drop for signal jets (no background subtraction)

// ./runSoftDrop -hard PythiaEventsTune14PtHat120.pu14  -nev 1

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
    //   Soft Drop for the signal jets
    //---------------------------------------------------------------------------
    
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


    softDropGroomer sdgSigZ01B00KT(0.1, 0.0, R);
    sdgSigZ01B00KT.setReclusteringAlgo(2);//0 = CA 1 = AKT 2 = KT
 
    std::vector<fastjet::PseudoJet> groomedJets_SigZ01B00KT = sdgSigZ01B00KT.doGrooming(jetCollectionSig);
    jetCollection jetCollectionSigSDZ01B00KT(groomedJets_SigZ01B00KT);
    jetCollectionSigSDZ01B00KT.addVector("sigJetSDZ01B00KTzg",    sdgSigZ01B00KT.getZgs());
    jetCollectionSigSDZ01B00KT.addVector("sigJetSDZ01B00KTndrop", sdgSigZ01B00KT.getNDroppedSubjets());
    jetCollectionSigSDZ01B00KT.addVector("sigJetSDZ01B00KTdr12",  sdgSigZ01B00KT.getDR12());
    jetCollectionSigSDZ01B00KT.addVector("sigJetSDZ01B00KTlogdr12",  sdgSigZ01B00KT.getLogDR12());
    jetCollectionSigSDZ01B00KT.addVector("sigJetSDZ01B00KTlogztheta",  sdgSigZ01B00KT.getLogZgDR12());
    jetCollectionSigSDZ01B00KT.addVector("sigJetSDZ01B00KTmass",  sdgSigZ01B00KT.getSubJetMass());
    jetCollectionSigSDZ01B00KT.addVector("sigJetSDZ01B00KTleadingtrack_pt",  sdgSigZ01B00KT.getSubJetLeadingTrackPt());
    
    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for signal jets
    //---------------------------------------------------------------------------
    
    softDropCounter sdcSig(0.0,0.0,R,0.0);
    sdcSig.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT
    sdcSig.run(jetCollectionSig);

    jetCollectionSig.addVector("sigJetRecur_jetpt",     sdcSig.getPts());
    jetCollectionSig.addVector("sigJetRecur_z",         sdcSig.getZgs());
    jetCollectionSig.addVector("sigJetRecur_dr12",      sdcSig.getDRs());
    jetCollectionSig.addVector("sigJetRecur_erad",      sdcSig.getErads());
    jetCollectionSig.addVector("sigJetRecur_logdr12",   sdcSig.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecur_logztheta", sdcSig.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecur_nSD",       sdcSig.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecur_zSD",       sdcSig.calculateNSD(1.0));

    softDropCounter sdcSigKT(0.0,0.0,R,0.0);
    sdcSigKT.setRecursiveAlgo(2);//0 = CA 1 = AKT 2 = KT
    sdcSigKT.run(jetCollectionSig);

    jetCollectionSig.addVector("sigJetRecurKT_jetpt",     sdcSigKT.getPts());
    jetCollectionSig.addVector("sigJetRecurKT_z",         sdcSigKT.getZgs());
    jetCollectionSig.addVector("sigJetRecurKT_dr12",      sdcSigKT.getDRs());
    jetCollectionSig.addVector("sigJetRecurKT_erad",      sdcSigKT.getErads());
    jetCollectionSig.addVector("sigJetRecurKT_logdr12",   sdcSigKT.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecurKT_logztheta", sdcSigKT.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecurKT_nSD",       sdcSigKT.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurKT_zSD",       sdcSigKT.calculateNSD(1.0));

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("eventWeight",   eventWeight);
    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetSDZ01B00",jetCollectionSigSDZ01B00);
    trwSig.addCollection("sigJetSDZ01B00KT",jetCollectionSigSDZ01B00KT);
    
    trwSig.fillTree();  //signal jets
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile("JetToyHIResultSoftDrop.root","RECREATE");
  trwSig.getTree()->Write();
  
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
