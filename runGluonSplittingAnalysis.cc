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

#include "include/extraInfo.hh"
#include "include/jetCollection.hh"
#include "include/softDropGroomer.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/Angularity.hh"

using namespace std;
using namespace fastjet;

// ./runGluonSplittingAnalysis -hard /Users/mverweij/mnt/eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -nev 10

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  cout << "will run on " << nEvent << " events" << endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 3.0;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);

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

    vector<PseudoJet> particlesMergedAll = mixer.particles();

    vector<double> eventWeight;
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

    //---------------------------------------------------------------------------
    //   look at first splitting of hard partons
    //---------------------------------------------------------------------------
    std::vector<double> drsplit;
    std::vector<double> drmom1;
    std::vector<double> drmom2;
    int id = 0;
    for(int ip = 0; ip<partons.size(); ++ip) {
      PseudoJet p = partons[ip];
      PseudoJet d1 = partonsFirstSplit[id];
      PseudoJet d2 = partonsFirstSplit[id+1];
      double dr = d1.delta_R(d2);
      drsplit.push_back(dr);
      
      drmom1.push_back(p.delta_R(d1));
      drmom2.push_back(p.delta_R(d2));
      id+=2;
    }
    

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    //calculate some angularities
    vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);
    
    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta00Z01.addVector("zgSigSDBeta00Z01",    sdgSigBeta00Z01.getZgs());
    jetCollectionSigSDBeta00Z01.addVector("ndropSigSDBeta00Z01", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01.addVector("dr12SigSDBeta00Z01",  sdgSigBeta00Z01.getDR12());
  

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);
    trw.addPartonCollection("partons",       partons);
    trw.addPartonCollection("partonsFirstSplit",       partonsFirstSplit);
    trw.addDoubleCollection("drSplit", drsplit);
    trw.addDoubleCollection("drmom1", drmom1);
    trw.addDoubleCollection("drmom2", drmom2);
    
    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("sigJetSDBeta00Z01",      jetCollectionSigSDBeta00Z01);
    
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIResultGluonSplittingnalysis.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
