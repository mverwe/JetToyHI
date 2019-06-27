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
#include "include/jewelMatcher.hh"
#include "include/gridSubtractor.hh"

using namespace std;
using namespace fastjet;

//./runJewelSub -hard  /eos/project/j/jetquenching/JetWorkshop2017/samples/jewel/DiJet/RecoilOn_0_10/Jewel_0_T_0.pu14 -nev 10

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  cout << "will run on " << nEvent << " events" << endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 3.0;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);

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

    vector<PseudoJet> particlesMerged = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // cluster hard event only
    vector<PseudoJet> particlesDummy, particlesReal;
    vector<PseudoJet> particlesBkg, particlesSig;
    SelectorVertexNumber(-1).sift(particlesMerged, particlesDummy, particlesReal);
    SelectorVertexNumber(0).sift(particlesReal, particlesSig, particlesBkg);

    for(int i = 0; i < (int)particlesDummy.size(); i++)
    {
       if(particlesDummy[i].perp() < 1e-5 && fabs(particlesDummy[i].pz()) > 2000)
       {
          particlesDummy.erase(particlesDummy.begin() + i);
          i = i - 1;
       }
    }

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));
    jetCollection jetCollectionSigJewel(GetCorrectedJets(jetCollectionSig.getJet(), particlesDummy));

    //---------------------------------------------------------------------------
    //   Groom the jets using JEWEL 4MomSub for recoils
    //---------------------------------------------------------------------------
    
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSig(0.1, 0.0, R);
    jetCollection jetCollectionSigSD(sdgSig.doGrooming(jetCollectionSig));
    jetCollectionSigSD.addVector("sigJetSDZg",    sdgSig.getZgs());
    jetCollectionSigSD.addVector("sigJetSDndrop", sdgSig.getNDroppedSubjets());
    jetCollectionSigSD.addVector("sigJetSDdr12",  sdgSig.getDR12());
    
    jetCollection jetCollectionSigSDJewel(GetCorrectedJets(sdgSig.getConstituents(), particlesDummy));
    vector<pair<PseudoJet, PseudoJet>> SigSDJewel = GetCorrectedSubJets(sdgSig.getConstituents1(), sdgSig.getConstituents2(), particlesDummy);
    jetCollectionSigSDJewel.addVector("sigJetSDJewelzg",   CalculateZG(SigSDJewel));
    jetCollectionSigSDJewel.addVector("sigJetSDJeweldr12", CalculateDR(SigSDJewel));

    //---------------------------------------------------------------------------
    //   Subtract recoils using JEWEL grid subtraction
    //---------------------------------------------------------------------------
    gridSubtractor gridSub(-5.,-M_PI,5.,M_PI,R,0.05);
    fastjet::ClusterSequenceArea csSub(gridSub.doGridSub1(jetCollectionSig, particlesDummy), jet_def, area_def);
    jetCollection jetCollectionSigGridSub(sorted_by_pt(jet_selector(csSub.inclusive_jets(10.))));

    //---------------------------------------------------------------------------
    //   Groom de grid subtracted jets
    //---------------------------------------------------------------------------
    softDropGroomer sdgSigGridSub(0.1, 0.0, R);
    jetCollection jetCollectionSigGridSubSD(sdgSigGridSub.doGrooming(jetCollectionSigGridSub));
    jetCollectionSigGridSubSD.addVector("sigJetGridSubSDZg",    sdgSigGridSub.getZgs());
    jetCollectionSigGridSubSD.addVector("sigJetGridSubSDndrop", sdgSigGridSub.getNDroppedSubjets());
    jetCollectionSigGridSubSD.addVector("sigJetGridSubSDdr12",  sdgSigGridSub.getDR12());
    
   
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);
    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("sigJetJewel",   jetCollectionSigJewel);
    trw.addCollection("sigJetSD",      jetCollectionSigSD);
    trw.addCollection("sigJetSDJewel", jetCollectionSigSDJewel);
    trw.addCollection("sigJetJewelGridSub", jetCollectionSigGridSub);
    trw.addCollection("sigJetJewelGridSubSD", jetCollectionSigGridSubSD);
    
            
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIResulJewelSub.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
