#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include <iostream>

using namespace fastjet;
using namespace std;

//root stuff
#include "TFile.h"
#include "TTree.h"

#include "thermalEvent.hh"
#include "pythiaEvent.hh"
#include "csSubtractor.hh"
#include "softDropGroomer.hh"
#include "softDropCounter.hh"
#include "treeWriter.hh"
#include "jetMatcher.hh"

int main () {

  // Number of events, generated and listed ones.
  unsigned int nEvent    = 10;

  //to write info to root tree
  treeWriter trw;
  trw.setTreeName("jetTree");

  //event generators
  thermalEvent thrm;
  pythiaEvent pyt;

  //Jet definition
  double R = 0.4;
  double ghost_etamax = 6.0;
  double ghost_area    = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghost_etamax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  unsigned int entryDiv = ((unsigned int)(nEvent/20));
  if(entryDiv==0) entryDiv = 1;
  for(unsigned int ie = 0; ie<nEvent; ++ie) {

    if(ie%entryDiv == 0) std::cout << "Event #" << ie << std::endl;

    //---------------------------------------------------------------------------
    //   produce event
    //---------------------------------------------------------------------------
    
    //create thermal event
    std::vector<fastjet::PseudoJet> particlesBkg = thrm.createThermalEvent();
    
    //create pythia event
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();

    //Embed pythia event into thermal event
    std::vector<fastjet::PseudoJet> particlesMerged;
    particlesMerged.reserve( particlesSig.size() + particlesBkg.size() ); // preallocate memory
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );
    particlesMerged.insert( particlesMerged.end(), particlesBkg.begin(), particlesBkg.end() );

    // std::cout << " nBkg: " << particlesBkg.size() << " nSig: " << particlesSig.size() << " nMerged: " << particlesMerged.size() << std::endl;

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the jets
    fastjet::ClusterSequenceArea csMerged(particlesMerged, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsMerged = sorted_by_pt(csMerged.inclusive_jets());

    // fastjet::ClusterSequenceArea csBkg(particlesBkg, jet_def, area_def);
    // std::vector<fastjet::PseudoJet> jetsBkg = sorted_by_pt(csBkg.inclusive_jets());

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsSig = sorted_by_pt(csSig.inclusive_jets(10.));

    //---------------------------------------------------------------------------
    //   background subtraction
    //---------------------------------------------------------------------------
    
    //run constituent subtraction on hybrid/embedded/merged event
    csSubtractor csSub;
    csSub.setInputParticles(particlesMerged);
    std::vector<fastjet::PseudoJet> jetsCS = csSub.doSubtraction();
    std::vector<double> rho;
    rho.push_back(csSub.getRho());
    std::vector<double> rhom;
    rhom.push_back(csSub.getRhoM());

    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------
    
    //SoftDrop grooming classic for CS jets
    softDropGroomer sdgCS;
    sdgCS.setZcut(0.1);
    sdgCS.setBeta(0.);
    sdgCS.setR0(R);
    sdgCS.setInputJets(jetsCS);
    std::vector<fastjet::PseudoJet> jetsCSSD = sdgCS.doGrooming();
    std::vector<double> zgCSSD = sdgCS.getZgs();
    std::vector<int> ndropCSSD = sdgCS.getNDroppedBranches();

    //SoftDrop emission counting for CS jets
    softDropCounter sdcCS(0.1, 0.0, 0.4, 0.1);
    sdcCS.run(jetsCS);
    std::vector<double> nCSSD = sdcCS.calculateNSD(0.0);

    //SoftDrop grooming classic for signal jets
    softDropGroomer sdgSig;
    sdgSig.setZcut(0.1);
    sdgSig.setBeta(0.);
    sdgSig.setR0(R);
    sdgSig.setInputJets(jetsSig);
    std::vector<fastjet::PseudoJet> jetsSigSD = sdgSig.doGrooming();
    std::vector<double> zgSigSD = sdgSig.getZgs();
    std::vector<int> ndropSigSD = sdgSig.getNDroppedBranches();

    //SoftDrop emission counting for signal jets
    softDropCounter sdcSig(0.1, 0.0, 0.4, 0.1);
    sdcSig.run(jetsSig);
    std::vector<double> nSigSD = sdcSig.calculateNSD(0.0);

    //match the CS jets to signal jets
    jetMatcher jmCS;
    jmCS.setMaxDist(R);
    jmCS.setBaseJets(jetsCS);
    jmCS.setTagJets(jetsSig);
    jmCS.matchJets();
    std::vector<fastjet::PseudoJet> jetsCSMatch = jmCS.getBaseJetsOrderedToTag();
    
    //reorder groomed observables for CS jets to gen jet order using CS jet matching indices
    //Note: it will break down if your grooming procedure throws away jets. In that case run matching for your collection
    std::vector<fastjet::PseudoJet> jetsCSSDMatch = jmCS.reorderedToTag(jetsCSSD);
    std::vector<double> zgCSSDMatch = jmCS.reorderedToTag(zgCSSD);
    std::vector<int> ndropCSSDMatch = jmCS.reorderedToTag(ndropCSSD);

    std::vector<double> nCSSDMatch = jmCS.reorderedToTag(nCSSD);

    //std::cout << "nCS: " << jetsCS.size() << " nCSSD: " << jetsCSSD.size() << " nCSSD: " << nCSSD.size() << std::endl;

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'double', 'int', and 'fastjet::PseudoJet' are supported
    trw.addJetCollection("sigJet",jetsSig);
    trw.addJetCollection("csJet",jetsCSMatch);
    trw.addJetCollection("sigJetSD",jetsSigSD);
    trw.addDoubleCollection("zgGenSD",zgSigSD);
    trw.addIntCollection("ndropGenSD",ndropSigSD);
    trw.addDoubleCollection("nSigSD",nSigSD);
    trw.addJetCollection("csJetSD",jetsCSSDMatch);
    trw.addDoubleCollection("zgCSSD",zgCSSDMatch);
    trw.addIntCollection("ndropCSSD",ndropCSSDMatch);
    trw.addDoubleCollection("nCSSD",nCSSDMatch);
    trw.addDoubleCollection("rho",rho);
    trw.addDoubleCollection("rhom",rhom);
    
    trw.fillTree();
    
  }//event loop

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResult.root","RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();
}
