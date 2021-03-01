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

using namespace std;
using namespace fastjet;

// This class runs time reclustering with background

// ./runTimeClus -hard PythiaEventsTune14PtHat120.pu14  -nev 1

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
  double R                   = 0.5;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

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

    // extract hard partons from first splitting
    fastjet::Selector parton_selector_split = SelectorVertexNumber(-2);
    vector<PseudoJet> partonsFirstSplit = parton_selector_split(particlesMergedAll);

    // select final state particles from hard event only
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    // select final state particles from background event only
    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );

    //charged particles
    fastjet::Selector charged_selector = SelectorIsCharged();
    vector<PseudoJet> particlesSigCh = charged_selector(particlesSig);

    //std::cout << "#particles: " << particlesSig.size() << " of which charged: " << particlesSigCh.size() << std::endl;

 
    //---------------------------------------------------------------------------
    //   look at first splitting of hard partons
    //---------------------------------------------------------------------------
    std::vector<double> drsplit;
    std::vector<double> tfsplit;
    double hbarc = 0.19732697;
    double GeVtofm = 1./hbarc; //~5.068;
    int id = 0;
    for(int ip = 0; ip<partons.size(); ++ip) {
      PseudoJet p = partons[ip];
      PseudoJet d1 = partonsFirstSplit[id];
      PseudoJet d2 = partonsFirstSplit[id+1];
      double dr = d1.delta_R(d2);
      drsplit.push_back(dr);
      double z1 = max(d1.e(),d2.e())/p.e();
      double z2 = min(d1.e(),d2.e())/p.e();
      tfsplit.push_back(1./(2.*z1*z2*p.e()*GeVtofm*(1-fastjet::cos_theta(d1,d2))));
    
      id+=2;
    }
    
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    // run the clustering, extract the signal charged jets
    fastjet::ClusterSequenceArea csSigCh(particlesSigCh, jet_def, area_def);
    jetCollection jetCollectionSigCh(sorted_by_pt(jet_selector(csSigCh.inclusive_jets(10.))));
    
    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for signal jets
    //---------------------------------------------------------------------------
    
    softDropCounter sdcSig(0.0,0.0,R,0.0);
    sdcSig.setRecursiveAlgo(3);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSig.run(jetCollectionSig);

    jetCollectionSig.addVector("sigJetRecur_jetpt",     sdcSig.getPts());
    jetCollectionSig.addVector("sigJetRecur_z",         sdcSig.getZgs());
    jetCollectionSig.addVector("sigJetRecur_dr12",      sdcSig.getDRs());
    jetCollectionSig.addVector("sigJetRecur_erad",      sdcSig.getErads());
    jetCollectionSig.addVector("sigJetRecur_logdr12",   sdcSig.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecur_logztheta", sdcSig.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecur_tf",        sdcSig.getTfs());
    jetCollectionSig.addVector("sigJetRecur_tfe",       sdcSig.getTfes());
    jetCollectionSig.addVector("sigJetRecur_nSD",       sdcSig.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecur_zSD",       sdcSig.calculateNSD(1.0));

    //find closest parton for each jet
    std::vector<int> partonmatch;
    std::vector<double> partonmatchdr;
    std::vector<fastjet::PseudoJet> sigJets =  jetCollectionSig.getJet();
    for(fastjet::PseudoJet p : sigJets) {
      int ipmin = -1;
      double drmin = 999.;
      for(int ip = 0; ip<partons.size(); ++ip) {
        double dr = p.delta_R(partons[ip]);
        if(dr<drmin) {
          drmin = dr;
          ipmin = ip;
        }
      }
      partonmatch.push_back(ipmin);
      partonmatchdr.push_back(drmin);
    }
    jetCollectionSig.addVector("sigJetRecur_partonMatchID", partonmatch);
    jetCollectionSig.addVector("sigJetRecur_partonMatchDr", partonmatchdr);

    softDropCounter sdcSigzcut(0.1,0.0,R,0.0);
    sdcSigzcut.setRecursiveAlgo(3);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSigzcut.run(jetCollectionSig);

    
    jetCollectionSig.addVector("sigJetRecurZcut_jetpt",     sdcSigzcut.getPts());
    jetCollectionSig.addVector("sigJetRecurZcut_z",         sdcSigzcut.getZgs());
    jetCollectionSig.addVector("sigJetRecurZcut_dr12",      sdcSigzcut.getDRs());
    jetCollectionSig.addVector("sigJetRecurZcut_erad",      sdcSigzcut.getErads());
    jetCollectionSig.addVector("sigJetRecurZcut_logdr12",   sdcSigzcut.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecurZcut_logztheta", sdcSigzcut.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecurZcut_tf",        sdcSigzcut.getTfs());
    jetCollectionSig.addVector("sigJetRecurZcut_tfe",       sdcSigzcut.getTfes());
    jetCollectionSig.addVector("sigJetRecurZcut_nSD",       sdcSigzcut.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurZcut_zSD",       sdcSigzcut.calculateNSD(1.0));


    //charged jets
    softDropCounter sdcSigzcutCh(0.1,0.0,R,0.0);
    sdcSigzcutCh.setRecursiveAlgo(3);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSigzcutCh.run(jetCollectionSigCh);

    
    jetCollectionSigCh.addVector("sigJetChRecurZcut_jetpt",     sdcSigzcutCh.getPts());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_z",         sdcSigzcutCh.getZgs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_dr12",      sdcSigzcutCh.getDRs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_erad",      sdcSigzcutCh.getErads());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_logdr12",   sdcSigzcutCh.getLog1DRs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_logztheta", sdcSigzcutCh.getLogzDRs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_tf",        sdcSigzcutCh.getTfs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_tfe",       sdcSigzcutCh.getTfes());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_nSD",       sdcSigzcutCh.calculateNSD(0.0));
    jetCollectionSigCh.addVector("sigJetChRecurZcut_zSD",       sdcSigzcutCh.calculateNSD(1.0));

    //find closest parton for each charged jet
    std::vector<int> partonmatchCh;
    std::vector<double> partonmatchdrCh;
    std::vector<fastjet::PseudoJet> sigJetsCh =  jetCollectionSigCh.getJet();
    for(fastjet::PseudoJet p : sigJetsCh) {
      int ipmin = -1;
      double drmin = 999.;
      for(int ip = 0; ip<partons.size(); ++ip) {
        double dr = p.delta_R(partons[ip]);
        if(dr<drmin) {
          drmin = dr;
          ipmin = ip;
        }
      }
      partonmatchCh.push_back(ipmin);
      partonmatchdrCh.push_back(drmin);
    }
    jetCollectionSigCh.addVector("sigJetChRecur_partonMatchID", partonmatchCh);
    jetCollectionSigCh.addVector("sigJetChRecur_partonMatchDr", partonmatchdrCh);
    

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("eventWeight",   eventWeight);

    trwSig.addPartonCollection("partons",       partons);
    trwSig.addPartonCollection("partonsFirstSplit",       partonsFirstSplit);
    trwSig.addDoubleCollection("drSplit", drsplit);
    trwSig.addDoubleCollection("tfSplit", tfsplit);
    
    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetCh",      jetCollectionSigCh);
    
    
    trwSig.fillTree();  //signal jets
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile("JetToyHIResultTimeClus.root","RECREATE");
  trwSig.getTree()->Write();
  
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
