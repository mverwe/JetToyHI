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

//./runSDGenVarious -hard  /eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -nev 10

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
    
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    
    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------
    
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta00Z01.addVector("zgSigSDBeta00Z01",    sdgSigBeta00Z01.getZgs());
    jetCollectionSigSDBeta00Z01.addVector("ndropSigSDBeta00Z01", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01.addVector("dr12SigSDBeta00Z01",  sdgSigBeta00Z01.getDR12());

    softDropGroomer sdgSigBeta00Z02(0.2, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z02(sdgSigBeta00Z02.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta00Z02.addVector("zgSigSDBeta00Z02",    sdgSigBeta00Z02.getZgs());
    jetCollectionSigSDBeta00Z02.addVector("ndropSigSDBeta00Z02", sdgSigBeta00Z02.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z02.addVector("dr12SigSDBeta00Z02",  sdgSigBeta00Z02.getDR12());
    
    softDropGroomer sdgSigBeta15Z05(0.5, 1.5, R);
    jetCollection jetCollectionSigSDBeta15Z05(sdgSigBeta15Z05.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta15Z05.addVector("zgSigSDBeta15Z05",    sdgSigBeta15Z05.getZgs());
    jetCollectionSigSDBeta15Z05.addVector("ndropSigSDBeta15Z05", sdgSigBeta15Z05.getNDroppedSubjets());
    jetCollectionSigSDBeta15Z05.addVector("dr12SigSDBeta15Z05",  sdgSigBeta15Z05.getDR12());

    softDropGroomer sdgSigBetam1Z01(0.1, -1., R);
    jetCollection jetCollectionSigSDBetam1Z01(sdgSigBetam1Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBetam1Z01.addVector("zgSigSDBetam1Z01",    sdgSigBetam1Z01.getZgs());
    jetCollectionSigSDBetam1Z01.addVector("ndropSigSDBetam1Z01", sdgSigBetam1Z01.getNDroppedSubjets());
    jetCollectionSigSDBetam1Z01.addVector("dr12SigSDBetam1Z01",  sdgSigBetam1Z01.getDR12());

    softDropGroomer sdgSigBetam1Z02(0.2, -1., R);
    jetCollection jetCollectionSigSDBetam1Z02(sdgSigBetam1Z02.doGrooming(jetCollectionSig));
    jetCollectionSigSDBetam1Z02.addVector("zgSigSDBetam1Z02",    sdgSigBetam1Z02.getZgs());
    jetCollectionSigSDBetam1Z02.addVector("ndropSigSDBetam1Z02", sdgSigBetam1Z02.getNDroppedSubjets());
    jetCollectionSigSDBetam1Z02.addVector("dr12SigSDBetam1Z02",  sdgSigBetam1Z02.getDR12());

    softDropGroomer sdgSigBetam2Z01(0.1, -2., R);
    jetCollection jetCollectionSigSDBetam2Z01(sdgSigBetam2Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBetam2Z01.addVector("zgSigSDBetam2Z01",    sdgSigBetam2Z01.getZgs());
    jetCollectionSigSDBetam2Z01.addVector("ndropSigSDBetam2Z01", sdgSigBetam2Z01.getNDroppedSubjets());
    jetCollectionSigSDBetam2Z01.addVector("dr12SigSDBetam2Z01",  sdgSigBetam2Z01.getDR12());

    softDropGroomer sdgSigBetam2Z005(0.05, -2., R);
    jetCollection jetCollectionSigSDBetam2Z005(sdgSigBetam2Z005.doGrooming(jetCollectionSig));
    jetCollectionSigSDBetam2Z005.addVector("zgSigSDBetam2Z005",    sdgSigBetam2Z005.getZgs());
    jetCollectionSigSDBetam2Z005.addVector("ndropSigSDBetam2Z005", sdgSigBetam2Z005.getNDroppedSubjets());
    jetCollectionSigSDBetam2Z005.addVector("dr12SigSDBetam2Z005",  sdgSigBetam2Z005.getDR12());

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);

    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("sigJetSDBeta00Z01",      jetCollectionSigSDBeta00Z01);
    trw.addCollection("sigJetSDBeta00Z02",      jetCollectionSigSDBeta00Z02);
    trw.addCollection("sigJetSDBeta15Z05",      jetCollectionSigSDBeta15Z05);
    trw.addCollection("sigJetSDBetam1Z01",      jetCollectionSigSDBetam1Z01);
    trw.addCollection("sigJetSDBetam1Z02",      jetCollectionSigSDBetam1Z02);
    trw.addCollection("sigJetSDBetam2Z01",      jetCollectionSigSDBetam2Z01);
    trw.addCollection("sigJetSDBetam2Z005",     jetCollectionSigSDBetam2Z005);
    
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResultSDGen.root","RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
