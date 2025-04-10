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
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"

#include "include/pythiaEvent.hh"
#include "include/thermalAlice.hh"

#include "include/csSubFullEventIterative.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");
  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIThermalAlice.root").c_str(), "RECREATE");

  double R = cmdline.value<double>("-R",0.2);

  // pthat for pythia jet generation
  double       ptHat = cmdline.value<double>("-pthat",60);

  cout << "will run on " << nEvent << " events" << endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double ghostRapMax         = 1.0;
  double ghost_area          = 0.5;
  int    active_area_repeats = 1;
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);

  double jetRapMax = 0.9-R;
  Selector jet_selector = SelectorAbsEtaMax(jetRapMax);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  //EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  //PYTHIA generator
  pythiaEvent pyt(ptHat, 14, -0.95, 0.95, false, false, 0);

  //THERMAL generator
  thermalAlice thrmAlice;

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( iev < nEvent )
  {
    // increment event number    
    iev++;

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    //---------------------------------------------------------------------------
    //   Creating PYTHIA and thermal background
    //---------------------------------------------------------------------------
    // Pythia
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();

    // Thermal
    vector<PseudoJet> particlesBkg = thrmAlice.createThermalEventAlice();

    // Mix them together
    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );

    //---------------------------------------------------------------------------
    //   jet clustering of Pythia
    //---------------------------------------------------------------------------
    JetDefinition jet_def(antikt_algorithm, R);

    fastjet::ClusterSequenceArea sigTruth(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig_Truth(sorted_by_pt(jet_selector(sigTruth.inclusive_jets(10.))));

    trw.addCollection("sigJet_Truth",        jetCollectionSig_Truth);
    //---------------------------------------------------------------------------
    //   jet clustering of Pythia+Thermal
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea sigEmbedded(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionSig_Embedded(sorted_by_pt(jet_selector(sigEmbedded.inclusive_jets(10.))));

    trw.addCollection("sigJet_Embedded",        jetCollectionSig_Embedded);
    //---------------------------------------------------------------------------
    //   Subtracting background
    //---------------------------------------------------------------------------
    csSubFullEventIterative csSub( {0.0} , {0.1}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSub.setInputParticles(particlesMerged);
    csSub.setMaxEta(1.0);
    fastjet::ClusterSequenceArea csEmbedded(csSub.Subtract(), jet_def, area_def);
    jetCollection jetCollectionSig_Subtracted(sorted_by_pt(jet_selector(csEmbedded.inclusive_jets(1.))));  

    trw.addCollection("sigJet_Subtracted",        jetCollectionSig_Subtracted);
    //---------------------------------------------------------------------------
    //   match Pythia jets to Pythia+Thermal jets
    //---------------------------------------------------------------------------
    jetMatcher jmCSFull(0.2);
    jmCSFull.setBaseJets(jetCollectionSig_Subtracted);
    jmCSFull.setTagJets(jetCollectionSig_Truth);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionSig_Subtracted);    

    trw.addCollection("sigJet_Matched",        jetCollectionSig_Subtracted);

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  fout->cd();
  TTree *trOut = trw.getTree();
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
