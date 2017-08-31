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

// ./runCSVariations -hard  /eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -pileup  /eos/project/j/jetquenching/JetWorkshop2017/samples/thermal/Mult7000/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10

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
    const int ncs = 5;
    double alpha[ncs] = {-2.,-1.,0.,1.,2.};
    std::vector<jetCollection> jetCollectionCSs;
    std::vector<double> rho;
    std::vector<double> rhom; 
    for(int ics = 0; ics<ncs; ++ics) {
      csSubtractor csSub(R, alpha[ics], -1, 0.005,ghostRapMax,jetRapMax);
      csSub.setInputParticles(particlesMerged);
      jetCollection jetCollectionCS(csSub.doSubtraction());

      if(ics==3) {
        //Background densities used by constituent subtraction
        rho.push_back(csSub.getRho());
        rhom.push_back(csSub.getRhoM());
      }
      jetCollectionCSs.push_back(jetCollectionCS);
    }


    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------
    

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSig(0.1, 0.0, R);
    jetCollection jetCollectionSigSD(sdgSig.doGrooming(jetCollectionSig));
    jetCollectionSigSD.addVector("zgSigSD",    sdgSig.getZgs());
    jetCollectionSigSD.addVector("ndropSigSD", sdgSig.getNDroppedSubjets());
    jetCollectionSigSD.addVector("dr12SigSD",  sdgSig.getDR12());

    //SoftDrop grooming classic for CS jets (zcut=0.1, beta=0)
    std::vector<jetCollection> jetCollectionCSSDs;
    
    for(int ics = 0; ics<ncs; ++ics) {
      softDropGroomer sdgCS(0.1, 0.0, R);
      jetCollection jetCollectionCSSD(sdgCS.doGrooming(jetCollectionCSs[ics]));
      jetCollectionCSSD.addVector(Form("csJetSD_%dzg",ics),    sdgCS.getZgs());
      jetCollectionCSSD.addVector(Form("csJetSD_%dndrop",ics), sdgCS.getNDroppedSubjets());
      jetCollectionCSSD.addVector(Form("csJetSD_%ddr12",ics),  sdgCS.getDR12());

      jetCollectionCSSDs.push_back(jetCollectionCSSD);
      //std::cout << "n CS jets: " << jetCollectionCSs[ics].getJet().size() <<
      //  "  n CSSD jets: " << jetCollectionCSSDs[ics].getJet().size() << std::endl;
    }
    
    
    //match the CS jets to signal jets
    for(int ics = 0; ics<ncs; ++ics) {
      jetMatcher jmCS(R);
      jmCS.setBaseJets(jetCollectionCSs[ics]);
      jmCS.setTagJets(jetCollectionSig);
      jmCS.matchJets();

      jmCS.reorderedToTag(jetCollectionCSs[ics]);
      jmCS.reorderedToTag(jetCollectionCSSDs[ics]);
    }

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("sigJetSD",      jetCollectionSigSD);

    for(int ics = 0; ics<ncs; ++ics) {
      trw.addCollection(Form("csJet_%d",ics),         jetCollectionCSs[ics]);
      trw.addCollection(Form("csJetSD_%d",ics),       jetCollectionCSSDs[ics]);
    }
    trw.addCollection("csRho",         rho);
    trw.addCollection("csRhom",        rhom);

    trw.addCollection("eventWeight",   eventWeight);
        
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResultCSVariations.root","RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
