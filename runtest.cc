#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "PU14/CmdLine.hh"

#include "include/ProgressBar.h"

#include "include/jetCollection.hh"
#include "include/thermalEvent.hh"
#include "include/pythiaEvent.hh"
#include "include/csSubtractor.hh"
#include "include/skSubtractor.hh"
#include "include/softDropGroomer.hh"
#include "include/softDropCounter.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/randomCones.hh"

using namespace std;
using namespace fastjet;

//to run do ./runtest -nev 100 #or whatever other number of events you want
int main (int argc, char ** argv)
{
  auto start_time = std::chrono::steady_clock::now();

  CmdLine cmdline(argc,argv);
  // inputs read from command line
  unsigned int nEvent = cmdline.value<unsigned int>("-nev",10);  // first argument: command line option; second argument: default value
  
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  // // Number of events, generated and listed ones.
  // unsigned int nEvent    = 10000;

  //to write info to root tree
  treeWriter trw("jetTree");

  //event generators
  thermalEvent thrm(12000, 0.7, -3.0, 3.0);
  pythiaEvent pyt(120., 14, -3.0, 3.0);

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 3.0;
  fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  for(unsigned int ie = 0; ie < nEvent; ie++)
  {
    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);

    //---------------------------------------------------------------------------
    //   produce event
    //---------------------------------------------------------------------------
    
    //create thermal event
    std::vector<fastjet::PseudoJet> particlesBkg = thrm.createThermalEvent();
    
    //create pythia event
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();

    //Embed pythia event into thermal event
    std::vector<fastjet::PseudoJet> particlesMerged;
    particlesMerged.reserve(particlesSig.size() + particlesBkg.size()); // preallocate memory
    particlesMerged.insert(particlesMerged.end(), particlesSig.begin(), particlesSig.end());
    particlesMerged.insert(particlesMerged.end(), particlesBkg.begin(), particlesBkg.end());

    // std::cout << " nBkg: " << particlesBkg.size() << " nSig: " << particlesSig.size() << " nMerged: " << particlesMerged.size() << std::endl;

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the jets
    fastjet::ClusterSequenceArea csMerged(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionMerged(sorted_by_pt(jet_selector(csMerged.inclusive_jets())));

    // fastjet::ClusterSequenceArea csBkg(particlesBkg, jet_def, area_def);
    // jetCollection jetCollectionBkg(sorted_by_pt(csBkg.inclusive_jets()));

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    //---------------------------------------------------------------------------
    //   background subtraction
    //---------------------------------------------------------------------------
    
    //run constituent subtraction on hybrid/embedded/merged event
    csSubtractor csSub(R, 1., -1, 0.005,ghostRapMax,jetRapMax);
    csSub.setInputParticles(particlesMerged);
    //csSub.setInputJets(jetCollectionMerged.getJet());
    jetCollection jetCollectionCS(csSub.doSubtraction());
        
    std::vector<double> rho;    rho.push_back(csSub.getRho());
    std::vector<double> rhom;   rhom.push_back(csSub.getRhoM());

    // randomCones rc(4,R,2.3,rho[0]);
    // rc.setInputParticles(particlesMerged);
    // jetCollection jetCollectionRC(rc.run());

    //run soft killer on hybrid/embedded/merged event
    skSubtractor skSub(0.4, 3.0);
    skSub.setInputParticles(particlesMerged);
    std::vector<fastjet::PseudoJet> skEvent = skSub.doSubtraction();
    std::vector<double> skPtThreshold;
    skPtThreshold.push_back(skSub.getPtThreshold());

    fastjet::ClusterSequenceArea csSK(skEvent, jet_def, area_def);
    jetCollection jetCollectionSK(sorted_by_pt(jet_selector(csSK.inclusive_jets())));

    //std::cout << "njets CS: " << jetsCS.size() << " SK: " << jetsSK.size() << std::endl;
    
    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------
    
    //SoftDrop grooming classic for CS jets
    softDropGroomer sdgCS(0.1, 0.0, R);
    jetCollection jetCollectionCSSD(sdgCS.doGrooming(jetCollectionCS));
    jetCollectionCSSD.addVector("zgCSSD",    sdgCS.getZgs());
    jetCollectionCSSD.addVector("ndropCSSD", sdgCS.getNDroppedBranches());

    //SoftDrop emission counting for CS jets
    softDropCounter sdcCS(0.1, 0.0, 0.4, 0.1);
    sdcCS.run(jetCollectionCS);
    jetCollectionCS.addVector("nCSSD", sdcCS.calculateNSD(0.0));
    jetCollectionCS.addVector("zCSSD", sdcCS.calculateNSD(1.0));

    //SoftDrop grooming classic for signal jets
    softDropGroomer sdgSig(0.1, 0.0, R);
    jetCollection jetCollectionSigSD(sdgSig.doGrooming(jetCollectionSig));
    jetCollectionSigSD.addVector("zgSigSD",    sdgSig.getZgs());
    jetCollectionSigSD.addVector("ndropSigSD", sdgSig.getNDroppedBranches());

    //SoftDrop emission counting for signal jets
    softDropCounter sdcSig(0.1, 0.0, 0.4, 0.1);
    sdcSig.run(jetCollectionSig);
    jetCollectionSig.addVector("nSigSD", sdcSig.calculateNSD(0.0));
    jetCollectionSig.addVector("zSigSD", sdcSig.calculateNSD(1.0));

    //match the CS jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionCS);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();

    jmCS.reorderedToTag(jetCollectionCS);
    jmCS.reorderedToTag(jetCollectionCSSD);

    //match the SK jets to signal jets
    jetMatcher jmSK(R);
    jmSK.setBaseJets(jetCollectionSK);
    jmSK.setTagJets(jetCollectionSig);
    jmSK.matchJets();
    
    jmSK.reorderedToTag(jetCollectionSK);

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

    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("csJet",         jetCollectionCS);
    trw.addCollection("sigJetSD",      jetCollectionSigSD);
    trw.addCollection("csJetSD",       jetCollectionCSSD);
    trw.addCollection("skJet",         jetCollectionSK);
    //trw.addCollection("randomCones",   jetCollectionRC);

    trw.addCollection("csRho",         rho);
    trw.addCollection("csRhom",        rhom);
    trw.addCollection("skPtThreshold", skPtThreshold);

    trw.addCollection("unsubJet",      jetCollectionMerged);
        
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResult.root","RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runtest: " << time_in_seconds << std::endl;
  
}
