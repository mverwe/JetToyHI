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
#include "include/softDropGroomer.hh"
#include "include/softDropCounter.hh"
#include "include/treeWriter.hh"
#include "include/JetCharge.hh"

#include "Pythia8/Pythia.h"
#include "include/extraInfo.hh"

#include "fastjet/contrib/QCDAwarePlugin.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib::QCDAwarePlugin;

// This class runs time reclustering with background

// ./runJetCharge -hard PythiaEventsTune14PtHat120.pu14  -nev 1

int main (int argc, char ** argv) {

  auto start_time = std::chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");
  double       ptHat = cmdline.value<double>("-pthat",40);//120.;
  unsigned int tune  = cmdline.value<int>("-tune",14);
  int          showerType = cmdline.value<int>("-shower",1);
  
  std::cout << "will run on " << nEvent << " events" << std::endl;
  std::cout << "showerType: " << showerType << std::endl;
  
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trwSig("jetTreeSig");

  bool runFull = true;
  bool runCharged = false;//true;
 
  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

  AntiKtMeasure *akt04dm = new AntiKtMeasure(0.4);
  QCDAwarePlugin *qcdawareakt04 = new QCDAwarePlugin(akt04dm);
  std::cout << "using " << qcdawareakt04->description() << std::endl;
  
  double jetRapMax = 100.;//3.0;
  fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  JetCharge jetCharge(0.5, -1.);
  JetCharge jetChargePtMin1(0.5, 1.);
  JetCharge jetChargePtMin2(0.5, 2.);

  //Pythia settings
  Pythia8::Pythia pythia;
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("HardQCD:all = on");
  pythia.readString(Form("PhaseSpace:pTHatMin = %.1f",ptHat));
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.readString(Form("PartonShowers:Model = %d",showerType)); //1: default 2: VINCIA  3: DIRE
  if(showerType==1) pythia.readString(Form("Tune:pp = %d",tune));
  pythia.init();
  
  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);
  
  // loop over events
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  for(unsigned int ie = 0; ie < nEvent; ie++) {
    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);
  
    //--------------------------------------------------------
    //     CREATE PYTHIA EVENT
    //--------------------------------------------------------
    pythia.next(); //generate next event

    std::vector<PseudoJet> particlesSig;
    std::vector<PseudoJet> finalPartons;
    std::vector<PseudoJet> partons;
    
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].isFinal()) { //all final state particles
        fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
        //p.set_user_info(new extraInfo(pythia.event[i].id(), 0));
        p.set_user_info(new PU14(pythia.event[i].id(), particlesSig.size(), 0));
        particlesSig.push_back(p);
      }

      if(pythia.event[i].isFinalPartonLevel()) {
        fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
        //p.set_user_info(new extraInfo(pythia.event[i].id(), -1));
        p.set_user_info(new PU14(pythia.event[i].id(), finalPartons.size(), -1));
        p.set_user_index(pythia.event[i].id());
        //std::cout << "pt: " << p.perp() << "  pdg: " << pythia.event[i].id() << std::endl;
        finalPartons.push_back(p);
      }

      if (pythia.event[i].status()==-23) {
        fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
        //p.set_user_info(new extraInfo(pythia.event[i].id(), -1));
        p.set_user_info(new PU14(pythia.event[i].id(), partons.size(), -1));
        partons.push_back(p);
      }

      // //find last partons in the parton shower
      // if(pythia.event[i].isParton()) { //true for a gluon, a quark or antiquark up to the b (but excluding top), and a diquark or antidiquark consisting of quarks up to the b.
      //   std::vector<int> daughters = pythia.event[i].daughterList();
      //   bool acc = true;
      //   for(int id : daughters) {
      //     if(pythia.event[id].isParton()) acc = false;
      //   }
      //   if(acc) {
      //     fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      //     p.set_user_info(new extraInfo(pythia.event[i].id(), 0));
      //     finalPartons.push_back(p);
      //   }
      // }
    }//pythia event loop
    // std::cout << "finished generating pythia event" << std::endl;


    std::vector<double> eventWeight;
    eventWeight.push_back(pythia.info.weight());
    trwSig.addCollection("eventWeight",   eventWeight);
    
    //charged particles
    fastjet::Selector charged_selector = SelectorIsCharged();
    vector<PseudoJet> particlesSigCh = charged_selector(particlesSig);

    //std::cout << "#particles: " << particlesSig.size() << " of which charged: " << particlesSigCh.size() << std::endl;

    //---------------------------------------------------------------------------
    //   clustering of final state partons
    //---------------------------------------------------------------------------
    fastjet::ClusterSequence qcdawareakt04cs(finalPartons, qcdawareakt04);
    const vector<PseudoJet> akt04PartonJets = sorted_by_pt(qcdawareakt04cs.inclusive_jets());
        
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------

    if(runFull) {
      // run the clustering, extract the signal jets
      fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
      jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

      //find closest parton for each jet
      std::vector<int> partonmatch;
      std::vector<double> partonmatchdr;
      std::vector<double> partonmatchpdg;
      std::vector<fastjet::PseudoJet> sigJets =  jetCollectionSig.getJet();
      for(fastjet::PseudoJet p : sigJets) {
        int ipmin = -1;
        double drmin = 999.;
        int pdgmin = -999;
        for(int ip = 0; ip<akt04PartonJets.size(); ++ip) {
          double dr = p.delta_R(akt04PartonJets[ip]);
          if(dr<drmin) {
            drmin = dr;
            ipmin = ip;
            pdgmin = akt04PartonJets[ip].user_index();
          }
        }
        partonmatch.push_back(ipmin);
        partonmatchdr.push_back(drmin);
        partonmatchpdg.push_back(pdgmin);
      }
      jetCollectionSig.addVector("sigJet_partonMatchID", partonmatch);
      jetCollectionSig.addVector("sigJet_partonMatchDr", partonmatchdr);
      jetCollectionSig.addVector("sigJet_partonMatchPDG", partonmatchpdg);
  
      //Calculate jet charge
      vector<double> jetChargeSig; jetChargeSig.reserve(jetCollectionSig.getJet().size());
      vector<double> jetChargePtmin1Sig; jetChargePtmin1Sig.reserve(jetCollectionSig.getJet().size());
      vector<double> jetChargePtmin2Sig; jetChargePtmin2Sig.reserve(jetCollectionSig.getJet().size());
      for(PseudoJet jet : jetCollectionSig.getJet()) {
        jetChargeSig.push_back(jetCharge.result(jet));
        jetChargePtmin1Sig.push_back(jetChargePtMin1.result(jet));
        jetChargePtmin2Sig.push_back(jetChargePtMin2.result(jet));
      }
      jetCollectionSig.addVector("sigJetCharge", jetChargeSig);
      jetCollectionSig.addVector("sigJetChargePtmin1", jetChargePtmin1Sig);
      jetCollectionSig.addVector("sigJetChargePtmin2", jetChargePtmin2Sig);

      trwSig.addCollection("sigJet",        jetCollectionSig);
    }

    if(runCharged) {
      // run the clustering, extract the signal charged jets
      fastjet::ClusterSequenceArea csSigCh(particlesSigCh, jet_def, area_def);
      jetCollection jetCollectionSigCh(sorted_by_pt(jet_selector(csSigCh.inclusive_jets(10.))));
    
      //find closest parton for each charged jet
      std::vector<int> partonmatchCh;
      std::vector<double> partonmatchdrCh;
      std::vector<double> partonmatchpdgCh;
      std::vector<fastjet::PseudoJet> sigJetsCh =  jetCollectionSigCh.getJet();
      for(fastjet::PseudoJet p : sigJetsCh) {
        int ipmin = -1;
        double drmin = 999.;
        int pdgmin = -999;
        for(int ip = 0; ip<akt04PartonJets.size(); ++ip) {
          double dr = p.delta_R(akt04PartonJets[ip]);
          if(dr<drmin) {
            drmin = dr;
            ipmin = ip;
            pdgmin = akt04PartonJets[ip].user_index();
          }
        }
        partonmatchCh.push_back(ipmin);
        partonmatchdrCh.push_back(drmin);
        partonmatchpdgCh.push_back(pdgmin);
      }
      jetCollectionSigCh.addVector("sigJetCh_partonMatchID", partonmatchCh);
      jetCollectionSigCh.addVector("sigJetCh_partonMatchDr", partonmatchdrCh);
      jetCollectionSigCh.addVector("sigJetCh_partonMatchPDG", partonmatchpdgCh);
      
      //Calculate jet charge
      vector<double> jetChargeSigCh; jetChargeSigCh.reserve(jetCollectionSigCh.getJet().size());
      for(PseudoJet jet : jetCollectionSigCh.getJet()) {
        jetChargeSigCh.push_back(jetCharge.result(jet));
      }
      jetCollectionSigCh.addVector("sigJetChCharge", jetChargeSigCh);
    
      //---------------------------------------------------------------------------
      //   write tree
      //---------------------------------------------------------------------------
    
      //Give variable we want to write out to treeWriter.
      //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported
      trwSig.addCollection("sigJetCh",      jetCollectionSigCh);
    }
    
    trwSig.addPartonCollection("partons",       partons);
    //trwSig.addPartonCollection("finalPartons",  finalPartons);
    trwSig.fillTree();
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TString strOut = Form("JetToyHIResultJetChargeQCDAware_PYTHIA_PTHAT%.0f.root",ptHat);
  if(showerType==2) strOut = Form("JetToyHIResultJetChargeQCDAware_VINCIA_PTHAT%.0f.root",ptHat);
  if(showerType==3) strOut = Form("JetToyHIResultJetChargeQCDAware_DIRE_PTHAT%.0f.root",ptHat);
  TFile *fout = new TFile(strOut.Data(),"RECREATE");
  trwSig.getTree()->Write();
  
  fout->Write();
  fout->Close();

  delete akt04dm;
  delete qcdawareakt04;
  
  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
