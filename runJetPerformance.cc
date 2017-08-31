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
#include "include/csSubtractor.hh"
#include "include/csSubtractorFullEvent.hh"
#include "include/skSubtractor.hh"
#include "include/softDropGroomer.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/randomCones.hh"
#include "include/Angularity.hh"

using namespace std;
using namespace fastjet;

//  /eos/project/j/jetquenching/JetWorkshop2017/samples/thermal/Mult7000/ThermalEventsMult7000PtAv1.20_0.pu14
// /eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14

// ./runJetPerformance -hard  /eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -pileup  /eos/project/j/jetquenching/JetWorkshop2017/samples/thermal/Mult7000/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10

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

  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);
    
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
    csSubtractor csSub(R, 1., -1, 0.005,ghostRapMax,jetRapMax);
    csSub.setInputParticles(particlesMerged);
    jetCollection jetCollectionCS(csSub.doSubtraction());

    //Background densities used by constituent subtraction
    std::vector<double> rho;    rho.push_back(csSub.getRho());
    std::vector<double> rhom;   rhom.push_back(csSub.getRhoM());

    //run full event constituent subtraction on mixed (hard+UE) event
    //csSubtractorFullEvent csSubGlobal(1., R, 0.005,ghostRapMax);
    //csSubGlobal.setRho(rho[0]);
    //csSubGlobal.setRhom(rhom[0]);
    //csSubGlobal.setInputParticles(particlesMerged);
    //std::vector<fastjet::PseudoJet> csEvent = csSubGlobal.doSubtraction();

    //cluster jets from constituent subtracted event
    //fastjet::ClusterSequenceArea csGlobal(csEvent, jet_def, area_def);
    //jetCollection jetCollectionCSGlobal(sorted_by_pt(jet_selector(csGlobal.inclusive_jets())));

    //run soft killer on mixed event
    const int nsk = 7;
    double gridsize[nsk] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5};
    //skSubtractor skSub[nsk];
    std::vector<std::vector<fastjet::PseudoJet>> skEvents;// skEvents.reserve(nsk);
    std::vector<std::vector<double>> skPtThresholds;      // skPtThresholds.reserve(nsk);
    std::vector<jetCollection> jetCollectionSKs;          // jetCollectionSKs.reserve(nsk);
    
    for(int is = 0; is<nsk; ++is) {
      skSubtractor skSub(gridsize[is], 3.0);
      skSub.setInputParticles(particlesMerged);
      std::vector<fastjet::PseudoJet> skEvent = skSub.doSubtraction();
      std::vector<double> skPtThreshold;
      skPtThreshold.push_back(skSub.getPtThreshold()); //SoftKiller pT threshold

      //cluster jets for soft killed event
      fastjet::ClusterSequenceArea csSK(skEvent, jet_def, area_def);
      jetCollection jetCollectionSK(sorted_by_pt(jet_selector(csSK.inclusive_jets())));
      
      skEvents.push_back(skEvent);
      skPtThresholds.push_back(skPtThreshold);
      jetCollectionSKs.push_back(jetCollectionSK);
    }

    //-------------------------------------------------------------
    //        jet matching
    //-------------------------------------------------------------
    
    //match the CS jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionCS);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();

    jmCS.reorderedToTag(jetCollectionCS);

    //match the SK jets to signal jets
    for(int is = 0; is<nsk; ++is) {
      jetMatcher jmSK(R);
      jmSK.setBaseJets(jetCollectionSKs[is]);
      jmSK.setTagJets(jetCollectionSig);
      jmSK.matchJets();
    
      jmSK.reorderedToTag(jetCollectionSKs[is]);
    }

    //match the jets from full-event CS subtraction to signal jets
    //jetMatcher jmCSGlobal(R);
    //jmCSGlobal.setBaseJets(jetCollectionCSGlobal);
    //jmCSGlobal.setTagJets(jetCollectionSig);
    //jmCSGlobal.matchJets();
    
    //jmCSGlobal.reorderedToTag(jetCollectionCSGlobal);
    
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("csJet",         jetCollectionCS);
    for(int is = 0; is<nsk; ++is) {
      trw.addCollection(Form("skJet%d",(int)(gridsize[is]*100.)),         jetCollectionSKs[is]);
      trw.addCollection(Form("skPtThreshold%d",(int)(gridsize[is]*100.)), skPtThresholds[is]);
    }
    //trw.addCollection("skJet",         jetCollectionSK);
    //trw.addCollection("csGlobJet",     jetCollectionCSGlobal);

    trw.addCollection("csRho",         rho);
    trw.addCollection("csRhom",        rhom);
    trw.addCollection("eventWeight",   eventWeight);

    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResultJetPerformance.root","RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runJetPerformance: " << time_in_seconds << std::endl;
}
