#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/ConstituentSubtractor.hh"

#include <iostream>

using namespace fastjet;
using namespace std;

#include "thermalEvent.hh"
#include "pythiaEvent.hh"
#include "csSubtractor.hh"

int main () {

  // Number of events, generated and listed ones.
  unsigned int nEvent    = 1;

  for(unsigned int ie = 0; ie<nEvent; ++ie) {
  
    thermalEvent thrm;
    std::vector<fastjet::PseudoJet> particlesBkg = thrm.createThermalEvent();

    //create pythia event
    pythiaEvent pyt;
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();

    //Embed pythia event into thermal event
    std::vector<fastjet::PseudoJet> particlesMerged;
    particlesMerged.reserve( particlesSig.size() + particlesBkg.size() ); // preallocate memory
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );
    particlesMerged.insert( particlesMerged.end(), particlesBkg.begin(), particlesBkg.end() );

    std::cout << " nBkg: " << particlesBkg.size() << " nSig: " << particlesSig.size() << " nMerged: " << particlesMerged.size() << std::endl;
        
    // choose a jet definition
    double R = 0.4;

    double ghost_etamax = 6.0;
    double ghost_area    = 0.005;
    int    active_area_repeats = 1;

    fastjet::GhostedAreaSpec ghost_spec(ghost_etamax, active_area_repeats, ghost_area);
    fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
    JetDefinition jet_def(antikt_algorithm, R);
  
    // run the clustering, extract the jets
    fastjet::ClusterSequenceArea csMerged(particlesMerged, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsMerged = sorted_by_pt(csMerged.inclusive_jets());

    fastjet::ClusterSequenceArea csBkg(particlesBkg, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsBkg = sorted_by_pt(csBkg.inclusive_jets());

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jetsSig = sorted_by_pt(csSig.inclusive_jets());
  
    // print out some infos
    cout << "Clustering with " << jet_def.description() << endl;

    // print the jets
    cout <<   "Bkg        pt y phi" << endl;
    for(fastjet::PseudoJet& jet : jetsBkg) {
      vector<PseudoJet> constituents = jet.constituents();
      if(jet.pt()>100.)
        std::cout << "jet " << jet.pt() << " " 
                  << jet.rap() << " " << jet.phi()
                  << " area: " << jet.area() << " nconst: " << constituents.size() << std::endl;
    }

    cout <<   "Sig        pt y phi" << endl;
    for(fastjet::PseudoJet& jet : jetsSig) {
      vector<PseudoJet> constituents = jet.constituents();
      if(jet.pt()>10.)
        std::cout << "jet " << jet.pt() << " " 
                  << jet.rap() << " " << jet.phi()
                  << " area: " << jet.area() << " nconst: " << constituents.size() << std::endl;
    }

    std::cout << "#jetsMerged: " << jetsMerged.size() << std::endl;
    cout <<   "Merged        pt y phi" << endl;
    for(fastjet::PseudoJet& jet : jetsMerged) {
      vector<PseudoJet> constituents = jet.constituents();
      if(jet.pt()>100.)
        std::cout << "jet " << jet.pt() << " " 
                  << jet.rap() << " " << jet.phi()
                  << " area: " << jet.area() << " nconst: " << constituents.size() << std::endl;
    }

    csSubtractor csSub;
    csSub.setInputParticles(particlesMerged);
    std::vector<fastjet::PseudoJet> jetsCS = csSub.doSubtraction();

    std::cout << "subtracted #jetsCS: " << jetsCS.size() << std::endl;
    cout <<   "CS        pt y phi" << endl;
    for(fastjet::PseudoJet& jet : jetsCS) {
      //vector<PseudoJet> constituents = jet.constituents();
      if(jet.pt()>20.)
        std::cout << "jet " << jet.pt() << " " 
                  << jet.rap() << " " << jet.phi() << std::endl; //
      //                << " area: " << jet.area() << " nconst: " << std::endl; // constituents.size() << std::endl;

    }
  }//event loop
}
