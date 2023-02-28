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
// ./runHardKtVLEJewel -hard /Users/mverweij/wrk/jetsPbPb/HardKtVLE/samples/jewel/pu14/simple-10_0.pu14 -nev 5000

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
    fastjet::Selector dummy_selector = SelectorVertexNumber(-1);
    vector<PseudoJet> particlesDummy = dummy_selector(particlesMergedAll);

    // select final state particles from hard event only
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    // select final state particles from background event only
    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );

    //remove ghosts from jewel dummies
    for(int i = 0; i < (int)particlesDummy.size(); i++) {
      if(particlesDummy[i].perp() < 1e-5 && fabs(particlesDummy[i].pz()) > 2000) {
        particlesDummy.erase(particlesDummy.begin() + i);
        i = i - 1;
      }
    }

    //---------------------------------------------------------------------------
    //   subtract medium response in full event
    //---------------------------------------------------------------------------  
    fastjet::contrib::ConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);  // distance in eta-phi plane
    subtractor.set_max_distance(0.5);  // free parameter for the maximal allowed distance between particle i and ghost k
    subtractor.set_alpha(0.);  // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the package alpha was multiplied by two but in newer versions this is not the case anymore
    //subtractor.set_scale_fourmomentum();
    subtractor.set_remove_all_zero_pt_particles(true);

    std::vector<fastjet::PseudoJet> subtracted_particles = subtractor.do_subtraction(particlesSig, particlesDummy);
    
     
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(subtracted_particles, jet_def, area_def);
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
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported
    trwSig.addCollection("eventWeight",   eventWeight);
  
    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetDYkTD",   jetCollectionSigDYKTD);
    
    trwSig.fillTree();  //signal jets
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile("JetToyHIResultHardKtVLEJewel.root","RECREATE");
  trwSig.getTree()->Write();
  
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
