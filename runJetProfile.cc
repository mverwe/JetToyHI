#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

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

void FillPtHistogramm(jetCollection jColl, TH1D* hPt){
  vector<fastjet::PseudoJet> jSig=jColl.getJet();
  for(int ijet=0;ijet<jSig.size();ijet++){
     fastjet::PseudoJet jet=jSig[ijet];
     hPt->Fill(jet.pt());
  }
}

void DesignRecurColl(jetCollection *jColl, softDropGroomer sDG, int ncs=0){
    jColl->addVector(Form("recur%i_jetpt",ncs),  sDG.getRecur_JetPt());
    jColl->addVector(Form("recur%i_logdr12",ncs),  sDG.getRecur_LogDR12());
    jColl->addVector(Form("recur%i_logztheta",ncs),  sDG.getRecur_LogZgDR12());
    jColl->addVector(Form("recur%i_n",ncs),  sDG.getRecur_N());
    jColl->addVector(Form("recur%i_z",ncs),  sDG.getRecur_z());
    jColl->addVector(Form("recur%i_erad",ncs),  sDG.getRecur_Erad());
}

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

  TH1D* hJetSigPt=new TH1D("hJetSigPt","hJetSigPt",600,0,300);
  TH1D* hJetEmbPt=new TH1D("hJetEmbPt","hJetEmbPt",600,0,300);
  TH1D* hJetToySigPt=new TH1D("hJetToySigPt","hJetToySigPt",600,0,300);
  TH1D* hJetToyEmbPt=new TH1D("hJetToyEmbPt","hJetToyEmbPt",600,0,300);
  TH1D* hJetInjectedPt=new TH1D("hJetInjectedPt","hJetInjectedPt",600,0,300);


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
    softDropGroomer sdgSig(0.0, 0.0, R);
    softDropGroomer sdgSig_Toy(0.0, 0.0, R);
    sdgSig.setReclusteringAlgo(2);//0 = CA 1 = AKT 2 = KT
    sdgSig.setRecursiveAlgo(2);//0 = CA 1 = AKT 2 = KT
    sdgSig.setInjectTracks(kFALSE,0);
    sdgSig_Toy.setReclusteringAlgo(2);
    sdgSig_Toy.setRecursiveAlgo(2);
    sdgSig_Toy.setInjectTracks(kTRUE,1);
    sdgSig_Toy.setMediumParameters(1,3);
    std::vector<fastjet::PseudoJet> groomedJets_Sig = sdgSig.doGrooming(jetCollectionSig);
    std::vector<fastjet::PseudoJet> groomedJets_Sig_Toy = sdgSig_Toy.doGrooming(jetCollectionSig);
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
    DesignRecurColl(&jetCollectionSigSD_Recur,sdgSig);
    jetCollection jetCollectionSigSD_Recur_Toy(groomedJets_Sig_Toy);
    DesignRecurColl(&jetCollectionSigSD_Recur_Toy,sdgSig_Toy);

    jetCollection jetCollectionSigSD_Recur_Injected(groomedJets_Sig_Toy);
    jetCollectionSigSD_Recur_Injected.addVector("injected_pt",  sdgSig_Toy.getInjectedpt());
    jetCollectionSigSD_Recur_Injected.addVector("injected_z",  sdgSig_Toy.getInjectedz());
    jetCollectionSigSD_Recur_Injected.addVector("injected_theta",  sdgSig_Toy.getInjectedtheta());

    //SoftDrop grooming classic for CS jets (zcut=0.1, beta=0)
    std::vector<jetCollection> jetCollectionCSSDs;
    std::vector<jetCollection> jetCollectionCSSDs_Recur;
    std::vector<jetCollection> jetCollectionCSSDs_Recur_Toy;

    for(int ics = 0; ics<ncs; ++ics) {
      softDropGroomer sdgCS(0.0, 0.0, R);
      softDropGroomer sdgCS_Toy(0.0, 0.0, R);
      sdgCS.setReclusteringAlgo(2);//0 = CA 1 = AKT 2 = KT
      sdgCS.setRecursiveAlgo(2);//0 = CA 1 = AKT 2 = KT
      sdgCS.setInjectTracks(kFALSE,0);
      sdgCS_Toy.setReclusteringAlgo(2);
      sdgCS_Toy.setRecursiveAlgo(2);
      sdgCS_Toy.setInjectTracks(kTRUE,1);
      sdgCS_Toy.setMediumParameters(1,3);
      std::vector<fastjet::PseudoJet> groomedJets_CS = sdgCS.doGrooming(jetCollectionCSs[ics]);
      std::vector<fastjet::PseudoJet> groomedJets_CS_Toy = sdgCS_Toy.doGrooming(jetCollectionCSs[ics]);
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
      DesignRecurColl(&jetCollectionCSSD_Recur,sdgCS,ics);
      jetCollection jetCollectionCSSD_Recur_Toy(groomedJets_CS_Toy);
      DesignRecurColl(&jetCollectionCSSD_Recur_Toy,sdgCS_Toy,ics);
    
      jetCollectionCSSDs.push_back(jetCollectionCSSD);
      jetCollectionCSSDs_Recur.push_back(jetCollectionCSSD_Recur);
      jetCollectionCSSDs_Recur_Toy.push_back(jetCollectionCSSD_Recur_Toy);

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

      /*jetMatcher jmCS_Recur(R);
      jmCS_Recur.setBaseJets(jetCollectionCSSDs_Recur[ics]);
      jmCS_Recur.setTagJets(jetCollectionCSSDs_Recur_Toy[ics]);
      jmCS_Recur.setTagJets(jetCollectionSigSD_Recur);
      jmCS_Recur.setTagJets(jetCollectionSigSD_Recur_Toy);
      jmCS_Recur.matchJets();

      jmCS_Recur.reorderedToTag(jetCollectionCSSDs_Recur[ics]);*/
    }

    //---------------------------------------------------------------------------
    //Generate jet-pt histograms
    //

    FillPtHistogramm(jetCollectionSigSD_Recur,hJetSigPt);
    FillPtHistogramm(jetCollectionCSSDs_Recur[0],hJetEmbPt);
    FillPtHistogramm(jetCollectionSigSD_Recur_Toy,hJetToySigPt);
    FillPtHistogramm(jetCollectionCSSDs_Recur_Toy[0],hJetToyEmbPt);
    FillPtHistogramm(jetCollectionSigSD_Recur_Injected,hJetInjectedPt);

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
        
    trw.fillTree(); //jetTree
    trw2.fillTree();    //TreeSig
    trw3.fillTree();    //TreeEmb
    trw4.fillTree();    //TreeToySig
    trw5.fillTree();    //TreeToyEmb
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

  hJetSigPt->Write();
  hJetEmbPt->Write();
  hJetToySigPt->Write();
  hJetToyEmbPt->Write();
  hJetInjectedPt->Write();

  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
