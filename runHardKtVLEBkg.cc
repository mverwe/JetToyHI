#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

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
#include "include/Angularity.hh"
#include "include/dyGroomer.hh"

using namespace std;
using namespace fastjet;

// This class runs large kT VLE analysis
// ./runHardKtVLEBkg -hard samples/PythiaEventsTune14PtHat120_10k.pu14 -pileup samples/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10

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
  fastjet::JetDefinition jet_def_ca(cambridge_algorithm, 999.);
  
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
    //std::cout << "iev: " << iev << std::endl;
    //if(iev<6) { iev++; continue;}
    //std::cout << "iev: " << iev << std::endl;
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

    // select final state particles from hard event only
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    // select final state particles from background event only
    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );

     
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    //---------------------------------------------------------------------------
    //   Dynamical grooming
    //---------------------------------------------------------------------------
    //std::cout << "do dynamical grooming signal jets" << std::endl;
    dyGroomer dygKTDSig(1);
    jetCollection jetCollectionSigDYKTD(dygKTDSig.doGrooming(jetCollectionSig));
    jetCollectionSigDYKTD.addVector("kappaSigDYKTD",    dygKTDSig.getKappas());
    jetCollectionSigDYKTD.addVector("rgSigDYKTD",       dygKTDSig.getDR12());   
    jetCollectionSigDYKTD.addVector("zgSigDYKTD",       dygKTDSig.getZgs());   
    jetCollectionSigDYKTD.addVector("ktgSigDYKTD",       dygKTDSig.getKts());

    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for signal jets
    //---------------------------------------------------------------------------
    softDropCounter sdcSig(0.1,0.0,R,0.0);
    sdcSig.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSig.run(jetCollectionSig);
    jetCollectionSig.addVector("sigJetRecur_nSD",       sdcSig.calculateNSD(0.0));

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
    //   Dynamical grooming for jet-by-jet subtraction
    //---------------------------------------------------------------------------
    dyGroomer dygKTDCS(1);
    jetCollection jetCollectionCSDYKTD(dygKTDCS.doGrooming(jetCollectionCS));
    jetCollectionCSDYKTD.addVector("kappaCSDYKTD",    dygKTDCS.getKappas());
    jetCollectionCSDYKTD.addVector("rgCSDYKTD",       dygKTDCS.getDR12());   
    jetCollectionCSDYKTD.addVector("zgCSDYKTD",       dygKTDCS.getZgs());   
    jetCollectionCSDYKTD.addVector("ktgCSDYKTD",      dygKTDCS.getKts());

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
    //   Dynamical grooming for full event CS subtraction
    //---------------------------------------------------------------------------
    dyGroomer dygKTDCSFull(1);
    jetCollection jetCollectionCSFullDYKTD(dygKTDCSFull.doGrooming(jetCollectionCSFull));
    jetCollectionCSFullDYKTD.addVector("kappaCSFullDYKTD",    dygKTDCSFull.getKappas());
    jetCollectionCSFullDYKTD.addVector("rgCSFullDYKTD",       dygKTDCSFull.getDR12());   
    jetCollectionCSFullDYKTD.addVector("zgCSFullDYKTD",       dygKTDCSFull.getZgs());   
    jetCollectionCSFullDYKTD.addVector("ktgCSFullDYKTD",      dygKTDCSFull.getKts());

    
    
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported
    trwSig.addCollection("eventWeight",   eventWeight);
  
    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetDYkTD",   jetCollectionSigDYKTD);
    
    trwSig.addCollection("csJet", jetCollectionCS);
    trwSig.addCollection("csJetDYkTD", jetCollectionCSDYKTD);

    trwSig.addCollection("csFullJet", jetCollectionCSFull);
    trwSig.addCollection("csFullJetDYkTD", jetCollectionCSFullDYKTD);
    
    trwSig.fillTree();  //signal jets
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile("JetToyHIResultHardKtVLEBkg.root","RECREATE");
  trwSig.getTree()->Write();
  
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
