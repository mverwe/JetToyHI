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

// This class runs time reclustering with background

// ./runMiguelsVariablesBkg -hard samples/PythiaEventsTune14PtHat120_10k.pu14 -pileup samples/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10
// ./runMiguelsVariablesBkg -hard samples/nikhef/PythiaEventsTune14PtHat120_0.pu14 -pileup samples/nikhef/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10
//https://indico.cern.ch/event/974749/contributions/4104667/attachments/2147059/3619205/Strong2020-JetObservables.pdf
//(SD_tau2, kappaTD), (tau2,SD_pt), (SD_rz, kappaTD), (SD_ptD, kappaktD)

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

  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);
  fastjet::contrib::Nsubjettiness  nSub2_beta1(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(1.));

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
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    //calculate some angularities
    //std::cout << "calc angularities signal jets" << std::endl;
    vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    vector<double> tau2Sig;  tau2Sig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
      tau2Sig.push_back(nSub2_beta1(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);
    jetCollectionSig.addVector("tau2Sig", tau2Sig);

    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------
    //std::cout << "groom signal jets" << std::endl;
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta00Z01.addVector("zgSigSDBeta00Z01",    sdgSigBeta00Z01.getZgs());
    jetCollectionSigSDBeta00Z01.addVector("ndropSigSDBeta00Z01", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01.addVector("dr12SigSDBeta00Z01",  sdgSigBeta00Z01.getDR12());

    //calculate some angularities
    //std::cout << "calc angularities groomed jets" << std::endl;
    vector<double> widthSigSD; widthSigSD.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> pTDSigSD;   pTDSigSD.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> tau2SigSD;  tau2SigSD.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    //need to get list of constituents of groomed jets
    std::vector<std::vector<fastjet::PseudoJet>> listOfConstituents = sdgSigBeta00Z01.getConstituents();
    int ij = 0;
    for(PseudoJet jet : jetCollectionSigSDBeta00Z01.getJet()) {
      fastjet::ClusterSequenceArea csSigSD(listOfConstituents[ij], jet_def_ca, area_def);
      std::vector<fastjet::PseudoJet> tempjets = csSigSD.inclusive_jets();
      //std::cout << "pt jet: " << jet.pt() << "  reclustered: " << tempjets[0].pt() << std::endl;
      widthSigSD.push_back(width.result(tempjets[0]));
      pTDSigSD.push_back(pTD.result(tempjets[0]));
      tau2SigSD.push_back(nSub2_beta1(tempjets[0]));
      ++ij;
    }
    jetCollectionSigSDBeta00Z01.addVector("widthSigSD", widthSigSD);
    jetCollectionSigSDBeta00Z01.addVector("pTDSigSD", pTDSigSD);
    jetCollectionSigSDBeta00Z01.addVector("tau2SigSD", tau2SigSD);

    //---------------------------------------------------------------------------
    //   Dynamical grooming
    //---------------------------------------------------------------------------
    //std::cout << "do dynamical grooming signal jets" << std::endl;
    dyGroomer dygTDSig(2);
    jetCollection jetCollectionSigDYTD(dygTDSig.doGrooming(jetCollectionSig));
    jetCollectionSigDYTD.addVector("kappaSigDYTD",    dygTDSig.getKappas());
    
    dyGroomer dygKTDSig(1);
    jetCollection jetCollectionSigDYKTD(dygKTDSig.doGrooming(jetCollectionSig));
    jetCollectionSigDYKTD.addVector("kappaSigDYKTD",    dygKTDSig.getKappas());
   
    //---------------------------------------------------------------------------
    // run the clustering, extract the signal charged jets
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csSigCh(particlesSigCh, jet_def, area_def);
    jetCollection jetCollectionSigCh(sorted_by_pt(jet_selector(csSigCh.inclusive_jets(10.))));
    //calculate some angularities
    vector<double> widthSigCh; widthSigCh.reserve(jetCollectionSigCh.getJet().size());
    vector<double> pTDSigCh;   pTDSigCh.reserve(jetCollectionSigCh.getJet().size());
    for(PseudoJet jet : jetCollectionSigCh.getJet()) {
      widthSigCh.push_back(width.result(jet));
      pTDSigCh.push_back(pTD.result(jet));
    }
    jetCollectionSig.addVector("widthSigCh", widthSigCh);
    jetCollectionSig.addVector("pTDSigCh", pTDSigCh);

   
    //std::cout << "start with bkg" << std::endl;
    

    //---------------------------------------------------------------------------
    //   jet clustering of signal+background jets
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csRaw(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionRaw(sorted_by_pt(jet_selector(csRaw.inclusive_jets(10.))));

    //---------------------------------------------------------------------------
    //   background subtraction
    //---------------------------------------------------------------------------

    //run jet-by-jet constituent subtraction on mixed (hard+UE) event
    csSubtractor csSub(R, 0., -1, 0.005,ghostRapMax,2.5);
    csSub.setInputParticles(particlesMerged);
    jetCollection jetCollectionCS(csSub.doSubtraction());

    //Background densities used by constituent subtraction
    std::vector<double> rho;
    std::vector<double> rhom;
    rho.push_back(csSub.getRho());
    rhom.push_back(csSub.getRhoM());

    //match CS jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionCS);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();

    jmCS.reorderedToTag(jetCollectionCS);

    //---------------------------------------------------------------------------
    //   Groom the CS jets
    //---------------------------------------------------------------------------
    // std::cout << "groom CS jets " << std::endl;
    // //SoftDrop grooming classic for CS jets (zcut=0.1, beta=0)
    // softDropGroomer sdgCSBeta00Z01(0.1, 0.0, R);
    // jetCollection jetCollectionCSSDBeta00Z01(sdgCSBeta00Z01.doGrooming(jetCollectionCS));
    // jetCollectionCSSDBeta00Z01.addVector("zgCSSDBeta00Z01",    sdgCSBeta00Z01.getZgs());
    // jetCollectionCSSDBeta00Z01.addVector("ndropCSSDBeta00Z01", sdgCSBeta00Z01.getNDroppedSubjets());
    // jetCollectionCSSDBeta00Z01.addVector("dr12CSSDBeta00Z01",  sdgCSBeta00Z01.getDR12());

    // std::cout << "calc angularities groomed CS jets" << std::endl;
    // //calculate some angularities
    // vector<double> widthCSSD; widthCSSD.reserve(jetCollectionCSSDBeta00Z01.getJet().size());
    // vector<double> pTDCSSD;   pTDCSSD.reserve(jetCollectionCSSDBeta00Z01.getJet().size());
    // vector<double> tau2CSSD;  tau2CSSD.reserve(jetCollectionCSSDBeta00Z01.getJet().size());

    // int ijc = 0;
    // std::vector<std::vector<fastjet::PseudoJet>> listOfConstituentsSDCS = sdgCSBeta00Z01.getConstituents();
    // for(PseudoJet jet : jetCollectionCSSDBeta00Z01.getJet()) {
    //   if(jet.pt() > 50.) {
    //     fastjet::ClusterSequenceArea csCSSD(listOfConstituentsSDCS[ijc], jet_def_ca, area_def);
    //     std::vector<fastjet::PseudoJet> tempjets = csCSSD.inclusive_jets();
    //     //std::cout << "CS pt jet: " << jet.pt() << "  reclustered: " << tempjets[0].pt() << std::endl;
    //     widthCSSD.push_back(width.result(tempjets[0]));
    //     pTDCSSD.push_back(pTD.result(tempjets[0]));
    //     tau2CSSD.push_back(nSub2_beta1(tempjets[0]));
    //   } else {
    //     widthCSSD.push_back(-999.);
    //     pTDCSSD.push_back(-999.);
    //     tau2CSSD.push_back(-999.);
    //   }
    //   ++ijc;
    // }
    // jetCollectionCSSDBeta00Z01.addVector("widthCSSD", widthSigCSSD);
    // jetCollectionCSSDBeta00Z01.addVector("pTDCSSD", pTDSigCSSD);
    // jetCollectionCSSDBeta00Z01.addVector("tau2CSSD", tau2CSSD);

    // //---------------------------------------------------------------------------
    // //   Dynamical grooming CS jets
    // //---------------------------------------------------------------------------
    // std::cout << "dyn grooming CS" << std::endl;
    // dyGroomer dygTDCS(2);
    // jetCollection jetCollectionCSDYTD(dygTDCS.doGrooming(jetCollectionCS));
    // jetCollectionCSDYTD.addVector("kappaCSDYTD",    dygTDCS.getKappas());

    // dyGroomer dygKTDCS(1);
    // jetCollection jetCollectionCSDYKTD(dygKTDCS.doGrooming(jetCollectionCS));
    // jetCollectionCSDYKTD.addVector("kappaCSDYKTD",    dygKTDCS.getKappas());


    //---------------------------------------------------------------------------
    //   Full event constituent subtraction
    //---------------------------------------------------------------------------
    //std::cout << "do full event CS" << std::endl;
    csSubtractorFullEvent csSubFull( 0., 0.25, 0.005, 2.5);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setRho(csSub.getRho());
    csSubFull.setRhom(csSub.getRhoM());
    csSubFull.setInputParticles(particlesMerged);
    
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtraction(), jet_def, area_def);

    std::vector<fastjet::PseudoJet> csFullJets = sorted_by_pt(jet_selector(fullSig.inclusive_jets(50.)));

    std::vector<fastjet::PseudoJet> csFullJetsClean;
    for(fastjet::PseudoJet jet : csFullJets) {
      if(jet.has_constituents())
        csFullJetsClean.push_back(jet);
    }
    
    jetCollection jetCollectionCSFull(csFullJetsClean);

    //match CSFull jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCSFull);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();

    jmCSFull.reorderedToTag(jetCollectionCSFull);

    //---------------------------------------------------------------------------
    //   Groom the CSFull jets
    //---------------------------------------------------------------------------
    //std::cout << "SD grooming CSFull" << std::endl;
    //SoftDrop grooming classic for CSFull jets (zcut=0.1, beta=0)
    softDropGroomer sdgCSFullBeta00Z01(0.1, 0.0, R);
    std::vector<fastjet::PseudoJet> groomedCSFullJets = sdgCSFullBeta00Z01.doGrooming(jetCollectionCSFull);
    //std::cout << "done SD grooming CSFull" << std::endl;
    jetCollection jetCollectionCSFullSDBeta00Z01(groomedCSFullJets);
    //std::cout << "SD CSFull jets in jetcoll" << std::endl;
    jetCollectionCSFullSDBeta00Z01.addVector("zgCSFullSDBeta00Z01",    sdgCSFullBeta00Z01.getZgs());
    jetCollectionCSFullSDBeta00Z01.addVector("ndropCSFullSDBeta00Z01", sdgCSFullBeta00Z01.getNDroppedSubjets());
    jetCollectionCSFullSDBeta00Z01.addVector("dr12CSFullSDBeta00Z01",  sdgCSFullBeta00Z01.getDR12());
    //std::cout << "size ungroomed CSFull jets: " << jetCollectionCSFull.getJet().size() << std::endl;
    //std::cout << "size groomed  CSFull jets: " << jetCollectionCSFullSDBeta00Z01.getJet().size() << std::endl;
    
    //calculate some angularities
    //std::cout << "angularities SD grooming CSFull" << std::endl;
    vector<double> widthCSFullSD; widthCSFullSD.reserve(jetCollectionCSFullSDBeta00Z01.getJet().size());
    vector<double> pTDCSFullSD;   pTDCSFullSD.reserve(jetCollectionCSFullSDBeta00Z01.getJet().size());
    vector<double> tau2CSFullSD;  tau2CSFullSD.reserve(jetCollectionCSFullSDBeta00Z01.getJet().size());

    int ijcf = 0;
    std::vector<std::vector<fastjet::PseudoJet>> listOfConstituentsSDCSFull = sdgCSFullBeta00Z01.getConstituents();
    for(PseudoJet jet : jetCollectionCSFullSDBeta00Z01.getJet()) {
      //std::cout << "CSF pt jet: " << jet.pt() << std::endl;
      if(jet.pt() > 50. || listOfConstituentsSDCSFull[ijcf].size()>10000) {
        //std::cout << "n constituents: " << listOfConstituentsSDCSFull[ijcf].size() << std::endl;
        fastjet::ClusterSequenceArea csCSFSD(listOfConstituentsSDCSFull[ijcf], jet_def_ca, area_def);
        std::vector<fastjet::PseudoJet> tempjets = csCSFSD.inclusive_jets();
        //std::cout << "CSF pt jet: " << jet.pt() << "  reclustered: " << tempjets[0].pt() << std::endl;
        widthCSFullSD.push_back(width.result(tempjets[0]));
        pTDCSFullSD.push_back(pTD.result(tempjets[0]));
        tau2CSFullSD.push_back(nSub2_beta1(tempjets[0]));
        //std::cout << "stored angularities " << std::endl;
      } else {
        widthCSFullSD.push_back(-999.);
        pTDCSFullSD.push_back(-999.);
        tau2CSFullSD.push_back(-999.);
      }
      ijcf++;
    }
    //std::cout << "add angulaity vectors" << std::endl;
    jetCollectionCSFullSDBeta00Z01.addVector("widthCSFullSD", widthCSFullSD);
    jetCollectionCSFullSDBeta00Z01.addVector("pTDCSFullSD", pTDCSFullSD);
    jetCollectionCSFullSDBeta00Z01.addVector("tau2CSFullSD", tau2CSFullSD);
    //std::cout << "add angulaity vectors  --p-- done" << std::endl;
    
    //---------------------------------------------------------------------------
    //   Dynamical grooming CSFull jets
    //---------------------------------------------------------------------------
    //std::cout << "dyn grooming CSFull" << std::endl;
    dyGroomer dygTDCSFull(2);
    jetCollection jetCollectionCSFullDYTD(dygTDCSFull.doGrooming(jetCollectionCSFull));
    jetCollectionCSFullDYTD.addVector("kappaCSFullDYTD",    dygTDCSFull.getKappas());

    dyGroomer dygKTDCSFull(1);
    jetCollection jetCollectionCSFullDYKTD(dygKTDCSFull.doGrooming(jetCollectionCSFull));
    jetCollectionCSFullDYKTD.addVector("kappaCSFullDYKTD",    dygKTDCSFull.getKappas());
    
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("eventWeight",   eventWeight);
    trwSig.addCollection("csRho",         rho);
    trwSig.addCollection("csRhom",        rhom);

    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetCh",      jetCollectionSigCh);

    trwSig.addCollection("sigJetSDBeta00Z01",  jetCollectionSigSDBeta00Z01);
    trwSig.addCollection("sigJetDYTD",         jetCollectionSigDYTD);
    trwSig.addCollection("sigJetDYkTD",        jetCollectionSigDYKTD);
    
    trwSig.addCollection("csJet",         jetCollectionCS);
    trwSig.addCollection("csFullJet",     jetCollectionCSFull);

    //trwSig.addCollection("csJetSDBeta00Z01",  jetCollectionCSSDBeta00Z01);
    //trwSig.addCollection("csJetDYTD",     jetCollectionCSDYTD);
    //trwSig.addCollection("csJetDYKTD",     jetCollectionCSDYKTD);
    trwSig.addCollection("csFullJetSDBeta00Z01",  jetCollectionCSFullSDBeta00Z01);
    trwSig.addCollection("csFullJetDYTD", jetCollectionCSFullDYTD);
    trwSig.addCollection("csFullJetDYKTD", jetCollectionCSFullDYKTD);
    
    trwSig.fillTree();  //signal jets
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile("JetToyHIResultMiguelsVariablesBkg.root","RECREATE");
  trwSig.getTree()->Write();
  
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
