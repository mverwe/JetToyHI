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
#include "treeWriter.hh"
#include "jetMatcher.hh"

int main () {

  // Number of events, generated and listed ones.
  unsigned int nEvent    = 100;

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

  //the stuff you want to store in the tree
  std::vector<double> zgSigSD;

  unsigned int entryDiv = ((unsigned int)(nEvent/20));
  for(unsigned int ie = 0; ie<nEvent; ++ie) {

    if(ie%entryDiv == 0) std::cout << "Event #" << ie << std::endl;
  
    std::vector<fastjet::PseudoJet> particlesBkg = thrm.createThermalEvent();

    //create pythia event
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();

    //Embed pythia event into thermal event
    std::vector<fastjet::PseudoJet> particlesMerged;
    particlesMerged.reserve( particlesSig.size() + particlesBkg.size() ); // preallocate memory
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );
    particlesMerged.insert( particlesMerged.end(), particlesBkg.begin(), particlesBkg.end() );

    // std::cout << " nBkg: " << particlesBkg.size() << " nSig: " << particlesSig.size() << " nMerged: " << particlesMerged.size() << std::endl;
    
    // run the clustering, extract the jets
    fastjet::ClusterSequenceArea csMerged(particlesMerged, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsMerged = sorted_by_pt(csMerged.inclusive_jets());

    fastjet::ClusterSequenceArea csBkg(particlesBkg, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsBkg = sorted_by_pt(csBkg.inclusive_jets());

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsSig = sorted_by_pt(csSig.inclusive_jets(10.));
  
    // print out some infos
    //    cout << "Clustering with " << jet_def.description() << endl;

    // print the jets
    // cout <<   "Bkg        pt y phi" << endl;
    // for(fastjet::PseudoJet& jet : jetsBkg) {
    //   vector<PseudoJet> constituents = jet.constituents();
    //   if(jet.pt()>100.)
    //     std::cout << "jet " << jet.pt() << " " 
    //               << jet.rap() << " " << jet.phi()
    //               << " area: " << jet.area() << " nconst: " << constituents.size() << std::endl;
    // }

    // cout <<   "Sig        pt y phi" << endl;
    // for(fastjet::PseudoJet& jet : jetsSig) {
    //   vector<PseudoJet> constituents = jet.constituents();
    //   if(jet.pt()>10.)
    //     std::cout << "jet " << jet.pt() << " " 
    //               << jet.rap() << " " << jet.phi()
    //               << " area: " << jet.area() << " nconst: " << constituents.size() << std::endl;
    // }

    // std::cout << "#jetsMerged: " << jetsMerged.size() << std::endl;
    // cout <<   "Merged        pt y phi" << endl;
    // for(fastjet::PseudoJet& jet : jetsMerged) {
    //   vector<PseudoJet> constituents = jet.constituents();
    //   if(jet.pt()>100.)
    //     std::cout << "jet " << jet.pt() << " " 
    //               << jet.rap() << " " << jet.phi()
    //               << " area: " << jet.area() << " nconst: " << constituents.size() << std::endl;
    // }

    csSubtractor csSub;
    csSub.setInputParticles(particlesMerged);
    std::vector<fastjet::PseudoJet> jetsCS = csSub.doSubtraction();

    // std::cout << "subtracted #jetsCS: " << jetsCS.size() << std::endl;
    // cout <<   "CS        pt y phi" << endl;
    // for(fastjet::PseudoJet& jet : jetsCS) {
    //   if(jet.pt()>20.)
    //     std::cout << "jet " << jet.pt() << " " 
    //               << jet.rap() << " " << jet.phi() << std::endl;
    // }

    //SoftDrop grooming classic
    softDropGroomer sdgCS;
    sdgCS.setZcut(0.1);
    sdgCS.setBeta(0.);
    sdgCS.setR0(R);
    sdgCS.setInputJets(jetsCS);
    std::vector<fastjet::PseudoJet> jetsCSSD = sdgCS.doGrooming();
    std::vector<double> zgCSSD = sdgCS.getZgs();
    std::vector<int> ndropCSSD = sdgCS.getNDroppedBranches();

    // cout <<   "CS SD      pt y phi zg ndrop" << endl;
    // for(unsigned int ij = 0; ij<jetsCSSD.size(); ++ij) {
    //   if(jetsCSSD[ij].pt()>20.)
    //     std::cout << "jet " << jetsCSSD[ij].pt() << " " 
    //               << jetsCSSD[ij].rap() << " " << jetsCSSD[ij].phi() << " "
    //               << zgCSSD[ij] << " " << ndropCSSD[ij] << std::endl;

    // }

    softDropGroomer sdgSig;
    sdgSig.setZcut(0.1);
    sdgSig.setBeta(0.);
    sdgSig.setR0(R);
    sdgSig.setInputJets(jetsSig);
    std::vector<fastjet::PseudoJet> jetsSigSD = sdgSig.doGrooming();
    //std::vector<double> zgSigSD = sdgSig.getZgs();
    zgSigSD = sdgSig.getZgs();
    std::vector<int> ndropSigSD = sdgSig.getNDroppedBranches();

    // cout <<   "Sig SD      pt y phi zg ndrop" << endl;
    // for(unsigned int ij = 0; ij<jetsSigSD.size(); ++ij) {
    //   if(jetsSigSD[ij].pt()>10.)
    //     std::cout << "jet " << jetsSigSD[ij].pt() << " " 
    //               << jetsSigSD[ij].rap() << " " << jetsSigSD[ij].phi() << " "
    //               << zgSigSD[ij] << " " << ndropSigSD[ij] << std::endl;

    // }

    jetMatcher jmCS;
    jmCS.setMaxDist(R);
    jmCS.setBaseJets(jetsCS);
    jmCS.setTagJets(jetsSig);
    jmCS.matchJets();
    std::vector<fastjet::PseudoJet> jetsCSMatch = jmCS.getBaseJetsOrderedToTag();

    jetMatcher jmCSSD;
    jmCSSD.setMaxDist(R);
    jmCSSD.setBaseJets(jetsCSSD);
    jmCSSD.setTagJets(jetsSig);
    jmCSSD.matchJets();
    std::vector<fastjet::PseudoJet> jetsCSSDMatch = jmCSSD.getBaseJetsOrderedToTag();
    std::vector<double> zgCSSDMatch = jmCSSD.reorderedToTag(zgCSSD);
    std::vector<int> ndropCSSDMatch = jmCSSD.reorderedToTag(ndropCSSD);

    jetMatcher jmSigSD;
    jmSigSD.setMaxDist(R);
    jmSigSD.setBaseJets(jetsSigSD);
    jmSigSD.setTagJets(jetsSig);
    jmSigSD.matchJets();
    std::vector<fastjet::PseudoJet> jetsSigSDMatch = jmSigSD.getBaseJetsOrderedToTag();
    std::vector<double> zgSigSDMatch = jmSigSD.reorderedToTag(zgSigSD);
    std::vector<int> ndropSigSDMatch = jmSigSD.reorderedToTag(ndropSigSD);
    
    trw.addJetCollection("sigJet",jetsSig);
    trw.addJetCollection("csJet",jetsCSMatch);
    trw.addJetCollection("sigJetSD",jetsSigSDMatch);
    trw.addDoubleCollection("zgGenSD",zgSigSDMatch);
    trw.addIntCollection("ndropGenSD",ndropSigSDMatch);
    trw.addJetCollection("csJetSD",jetsCSSDMatch);
    trw.addDoubleCollection("zgCSSD",zgCSSDMatch);
    trw.addIntCollection("ndropCSSD",ndropCSSDMatch);

    // trw.addIntCollection("ndropGen",ndropSigSD);

    trw.fillTree();
    
  }//event loop

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResult.root","RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();
}
