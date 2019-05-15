#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2.h"

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
  treeWriter trw2("recursiveTreeSig");
  treeWriter trw3("recursiveTreeEmb");
  treeWriter trw4("recursiveTreeToySig");
  treeWriter trw5("recursiveTreeToyEmb");
  treeWriter trw6("injectedTree");

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  TH2F *h2PtJetRho = new TH2F("h2PtJetRho","",150,0,150,20,0,0.5);
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
    /*
    // cluster hard event only
    std::vector<fastjet::PseudoJet> particlesBkg, particlesSig;
    std::vector<fastjet::PseudoJet> particlesNCh, particlesCh;   
    //   SelectorIsCharged().sift(particlesMerged, particlesCh, particlesNCh); // this sifts the full event into two vectors of PseudoJet, one for the hard event, one for the underlying event
    //SelectorIsHard().sift(particlesCh,particlesSig,particlesBkg);
    SelectorIsHard().sift(particlesMerged,particlesSig,particlesBkg);
    */
    vector<PseudoJet> particlesDummy, particlesReal;
    vector<PseudoJet> particlesBkg, particlesSig;
    SelectorVertexNumber(-1).sift(particlesMerged, particlesDummy, particlesReal);
    SelectorVertexNumber(0).sift(particlesReal, particlesSig, particlesBkg);

    for(int i = 0; i < (int)particlesDummy.size(); i++) {
      if(particlesDummy[i].perp() < 1e-5 && fabs(particlesDummy[i].pz()) > 2000) {
        particlesDummy.erase(particlesDummy.begin() + i);
        i = i - 1;
      }
    }
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));
    jetCollection jetCollectionSigJewel(GetCorrectedJets(jetCollectionSig.getJet(), particlesDummy));


    
    //---------------------------------------------------------------------------
    //   background subtraction
    //---------------------------------------------------------------------------
    
    //run jet-by-jet constituent subtraction on mixed (hard+UE) event
    const int ncs = 1;
    double alpha[ncs] = {1.};
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
    softDropGroomer sdgSig_Toy(0.1, 0.0, R);
    sdgSig.setEventWeight(mixer.hard_weight());
    sdgSig.setReclusteringAlgo(0);//0 = CA 1 = AKT 2 = KT
    sdgSig.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT
    sdgSig.setInjectTracks(kFALSE,0);
    sdgSig_Toy.setReclusteringAlgo(0);
    sdgSig_Toy.setRecursiveAlgo(0);
    sdgSig_Toy.setInjectTracks(kTRUE,1);
    sdgSig_Toy.setMediumParameters(1,3);
    std::vector<fastjet::PseudoJet> groomedJets_Sig = sdgSig.doGroomingWithJewelSub(jetCollectionSig, particlesDummy);
    std::vector<fastjet::PseudoJet> groomedJets_Sig_Toy = sdgSig_Toy.doGroomingWithJewelSub(jetCollectionSig, particlesDummy);
    TH2F *temph2PtJetRho = sdgSig.calculateProfile(jetCollectionSig);
    h2PtJetRho->Add(temph2PtJetRho,1);
    jetCollection jetCollectionSigSD(groomedJets_Sig);
    jetCollectionSigSD.addVector("zgSigSD",    sdgSig.getZgs());
    jetCollectionSigSD.addVector("ndropSigSD", sdgSig.getNDroppedSubjets());
    jetCollectionSigSD.addVector("dr12SigSD",  sdgSig.getDR12());
    jetCollectionSigSD.addVector("logdr12",  sdgSig.getLogDR12());
    jetCollectionSigSD.addVector("logztheta",  sdgSig.getLogZgDR12());
    jetCollectionSigSD.addVector("sjMass",  sdgSig.getSubJetMass());
    jetCollectionSigSD.addVector("leadingtrack_pt",  sdgSig.getSubJetLeadingTrackPt());
    for(int j=0;j<8;j++){
      jetCollectionSigSD.addVector(Form("jetProfile%d",j), sdgSig.getJetProfile(j));
    }

    jetCollection jetCollectionSigSD_Recur(groomedJets_Sig);
    jetCollectionSigSD_Recur.addVector("recur_jetpt",  sdgSig.getRecur_JetPt());
    jetCollectionSigSD_Recur.addVector("recur_logdr12",  sdgSig.getRecur_LogDR12());
    jetCollectionSigSD_Recur.addVector("recur_logztheta",  sdgSig.getRecur_LogZgDR12());
    jetCollectionSigSD_Recur.addVector("recur_n",  sdgSig.getRecur_N());


    jetCollection jetCollectionSigSD_Recur_Toy(groomedJets_Sig_Toy);
    jetCollectionSigSD_Recur_Toy.addVector("recur_jetpt",  sdgSig_Toy.getRecur_JetPt());
    jetCollectionSigSD_Recur_Toy.addVector("recur_logdr12",  sdgSig_Toy.getRecur_LogDR12());
    jetCollectionSigSD_Recur_Toy.addVector("recur_logztheta",  sdgSig_Toy.getRecur_LogZgDR12());
    jetCollectionSigSD_Recur_Toy.addVector("recur_n",  sdgSig_Toy.getRecur_N());

    jetCollection jetCollectionSigSD_Recur_Injected(groomedJets_Sig_Toy);
    jetCollectionSigSD_Recur_Injected.addVector("injected_pt",  sdgSig_Toy.getInjectedpt());
    jetCollectionSigSD_Recur_Injected.addVector("injected_z",  sdgSig_Toy.getInjectedz());
    jetCollectionSigSD_Recur_Injected.addVector("injected_theta",  sdgSig_Toy.getInjectedtheta());
    


    //SoftDrop grooming classic for CS jets (zcut=0.1, beta=0)
    std::vector<jetCollection> jetCollectionCSSDs;
    std::vector<jetCollection> jetCollectionCSSDs_Recur;
    std::vector<jetCollection> jetCollectionCSSDs_Recur_Toy;

    for(int ics = 0; ics<ncs; ++ics) {
      softDropGroomer sdgCS(0.1, 0.0, R);
      softDropGroomer sdgCS_Toy(0.1, 0.0, R);
      sdgCS.setReclusteringAlgo(0);//0 = CA 1 = AKT 2 = KT
      sdgCS.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT
      sdgCS.setInjectTracks(kFALSE,0);
      sdgCS_Toy.setReclusteringAlgo(0);
      sdgCS_Toy.setRecursiveAlgo(0);
      sdgCS_Toy.setInjectTracks(kTRUE,1);
      sdgCS_Toy.setMediumParameters(1,3);
      std::vector<fastjet::PseudoJet> groomedJets_CS = sdgCS.doGroomingWithJewelSub(jetCollectionCSs[ics], particlesDummy);
      std::vector<fastjet::PseudoJet> groomedJets_CS_Toy = sdgCS_Toy.doGroomingWithJewelSub(jetCollectionCSs[ics], particlesDummy);
      jetCollection jetCollectionCSSD(groomedJets_CS);
      jetCollectionCSSD.addVector(Form("csJetSD_%dzg",ics),    sdgCS.getZgs());
      jetCollectionCSSD.addVector(Form("csJetSD_%dndrop",ics), sdgCS.getNDroppedSubjets());
      jetCollectionCSSD.addVector(Form("csJetSD_%ddr12",ics),  sdgCS.getDR12());
      jetCollectionCSSD.addVector(Form("log_%ddr12",ics),  sdgCS.getLogDR12());
      jetCollectionCSSD.addVector(Form("log_%dztheta",ics),  sdgCS.getLogZgDR12());
      jetCollectionCSSD.addVector(Form("sjMass_CS%d",ics),  sdgCS.getSubJetMass());
      jetCollectionCSSD.addVector(Form("leadingtrack_pt_CS%d",ics),  sdgCS.getSubJetLeadingTrackPt());
      for(int j=0;j<8;j++){
      jetCollectionCSSD.addVector(Form("jetProfile%d",j), sdgCS.getJetProfile(j));
      }
      jetCollection jetCollectionCSSD_Recur(groomedJets_CS);
      jetCollectionCSSD_Recur.addVector(Form("recur%d_jetpt",ics),  sdgCS.getRecur_JetPt());
      jetCollectionCSSD_Recur.addVector(Form("recur%d_logdr12",ics),  sdgCS.getRecur_LogDR12());
      jetCollectionCSSD_Recur.addVector(Form("recur%d_logztheta",ics),  sdgCS.getRecur_LogZgDR12());
      jetCollectionCSSD_Recur.addVector(Form("recur%d_n",ics),  sdgCS.getRecur_N());

      jetCollection jetCollectionCSSD_Recur_Toy(groomedJets_CS_Toy);
      jetCollectionCSSD_Recur_Toy.addVector(Form("recur_Toy%d_jetpt",ics),  sdgCS_Toy.getRecur_JetPt());
      jetCollectionCSSD_Recur_Toy.addVector(Form("recur_Toy%d_logdr12",ics),  sdgCS_Toy.getRecur_LogDR12());
      jetCollectionCSSD_Recur_Toy.addVector(Form("recur_Toy%d_logztheta",ics),  sdgCS_Toy.getRecur_LogZgDR12());
      jetCollectionCSSD_Recur_Toy.addVector(Form("recur_Toy%d_n",ics),  sdgCS_Toy.getRecur_N());
    
      jetCollectionCSSDs.push_back(jetCollectionCSSD);
      jetCollectionCSSDs_Recur.push_back(jetCollectionCSSD_Recur);
      jetCollectionCSSDs_Recur_Toy.push_back(jetCollectionCSSD_Recur_Toy);

      //std::cout << "n CS jets: " << jetCollectionCSs[ics].getJet().size() <<
      //  "  n CSSD jets: " << jetCollectionCSSDs[ics].getJet().size() << std::endl;
    }
    
    /* 
    //match the CS jets to signal jets
    for(int ics = 0; ics<ncs; ++ics) {
      jetMatcher jmCS(R);
      jmCS.setBaseJets(jetCollectionCSs[ics]);
      jmCS.setTagJets(jetCollectionSig);
      jmCS.matchJets();

      jmCS.reorderedToTag(jetCollectionCSs[ics]);
      jmCS.reorderedToTag(jetCollectionCSSDs[ics]);

     // jetMatcher jmCS_Recur(R);
     // jmCS_Recur.setBaseJets(jetCollectionCSSDs_Recur[ics]);
     // jmCS_Recur.setTagJets(jetCollectionCSSDs_Recur_Toy[ics]);
     // jmCS_Recur.setTagJets(jetCollectionSigSD_Recur);
     // jmCS_Recur.setTagJets(jetCollectionSigSD_Recur_Toy);
     // jmCS_Recur.matchJets();

      jmCS_Recur.reorderedToTag(jetCollectionCSSDs_Recur[ics]);
    }
    */
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("sigJetSD",      jetCollectionSigSD);
    trw2.addRecurCollection("sigJetSD_Recur",       jetCollectionSigSD_Recur);
    trw4.addRecurCollection("sigJetSD_Recur_Toy",       jetCollectionSigSD_Recur_Toy);
    trw6.addRecurCollection("injectedTracks",       jetCollectionSigSD_Recur_Injected);
    for(int ics = 0; ics<ncs; ++ics) {
      trw.addCollection(Form("csJet_%d",ics),         jetCollectionCSs[ics]);
      trw.addCollection(Form("csJetSD_%d",ics),       jetCollectionCSSDs[ics]);
      trw3.addRecurCollection(Form("csJetSD_%d_Recur",ics),       jetCollectionCSSDs_Recur[ics]);
      trw5.addRecurCollection(Form("csJetSD_%d_Recur_Toy",ics),       jetCollectionCSSDs_Recur_Toy[ics]);
    }
    trw.addCollection("csRho",         rho);
    trw.addCollection("csRhom",        rhom);

    trw.addCollection("eventWeight",   eventWeight);
        
    trw.fillTree();
    trw2.fillTree();
    trw3.fillTree();
    trw4.fillTree();
    trw5.fillTree();
    trw6.fillTree();


  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();
  TTree *trOut2 = trw2.getTree();
  TTree *trOut3 = trw3.getTree();
  TTree *trOut4 = trw4.getTree();
  TTree *trOut5 = trw5.getTree();
  TTree *trOut6 = trw6.getTree();

  TFile *fout = new TFile("JetToyHIResultCSVariations.root","RECREATE");
  trOut->Write();
  trOut2->Write();
  trOut3->Write();
  trOut4->Write();
  trOut5->Write();
  trOut6->Write();

  h2PtJetRho->Write();

  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
