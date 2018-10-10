#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"

#include "include/jetCollection.hh"
#include "include/sharedLayerSubtractor.hh"

#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"

using namespace std;
using namespace fastjet;

// ./runSharedLayerSubtraction -hard  /eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -pileup  /eos/project/j/jetquenching/JetWorkshop2017/samples/thermal/Mult7000/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10

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
  treeWriter trw("jetTree");

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

  //Angularity width(1.,1.,R);
  //Angularity pTD(0.,2.,R);
    
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

    // run the clustering, extract the unsubtracted jets
    ClusterSequenceArea csMerged(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionMerged(sorted_by_pt(jet_selector(csMerged.inclusive_jets())));
    
    //---------------------------------------------------------------------------
    //   background subtraction
    //---------------------------------------------------------------------------
    
    //run jet-by-jet constituent subtraction on mixed (hard+UE) event
    sharedLayerSubtractor sharedLayerSub(R,0.005,ghostRapMax,jetRapMax);
    sharedLayerSub.setInputParticles(particlesMerged);
    jetCollection jetCollectionSL(sharedLayerSub.doSubtraction());

    std::vector<double> rho;
    rho.push_back(sharedLayerSub.getRho());
    std::vector<double> rhoSigma;
    rhoSigma.push_back(sharedLayerSub.getRhoSigma());

    std::vector<double> pTDBkg;
    pTDBkg.push_back(sharedLayerSub.getPTDBkg());
    std::vector<double> pTDBkgSigma;
    pTDBkgSigma.push_back(sharedLayerSub.getPTDBkgSigma());
    
    //match the subtracted jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionSL);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();
    jmCS.reorderedToTag(jetCollectionSL);
  

    //match the unsubtracted jets to signal jets
    jetMatcher jmUnSub(R);
    jmUnSub.setBaseJets(jetCollectionMerged);
    jmUnSub.setTagJets(jetCollectionSig);
    jmUnSub.matchJets();
    jmUnSub.reorderedToTag(jetCollectionMerged);


    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trw.addCollection("sigJet",        jetCollectionSig, true);
    trw.addCollection("slJet",         jetCollectionSL,  true);
    
    trw.addCollection("rho",           rho);
    trw.addCollection("rhoSigma",      rhoSigma);

    trw.addCollection("pTDBkg",        pTDBkg);
    trw.addCollection("pTDBkgSigma",   pTDBkgSigma);

    trw.addCollection("unsubJet",      jetCollectionMerged, true);
    
    trw.addCollection("eventWeight",   eventWeight);
        
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResultSharedLayers.root","RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  std::cout << "Check JetToyHIResultSharedLayers.root for results" << std::endl;

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
