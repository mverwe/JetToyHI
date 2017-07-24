#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

// #include "include/jetCollection.hh"
#include "include/thermalEvent.hh"
// #include "include/pythiaEvent.hh"
// #include "include/csSubtractor.hh"
// #include "include/skSubtractor.hh"
// #include "include/softDropGroomer.hh"
// #include "include/softDropCounter.hh"
// #include "include/treeWriter.hh"
// #include "include/jetMatcher.hh"
// #include "include/randomCones.hh"

using namespace std;
using namespace fastjet;

int main ()
{
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);
  
  // Number of events, generated and listed ones.
  unsigned int nEvent    = 100;

  //event generators
  unsigned int mult = 12000;
  double       ptAve = 0.7;
  thermalEvent thrm(mult,ptAve, -3.0, 3.0);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  //output text file
  ofstream fout;
  const char *dir = getenv("PWD");//"/eos/user/m/mverweij/JetWorkshop2017/samples/";
  TString outFileName = Form("%s/ThermalEventsMult%dPtAv%.2f.pu14",dir,mult,ptAve);
  
  fout.open(outFileName.Data());
  
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  for(unsigned int ie = 0; ie < nEvent; ie++) {
    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);

    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);
    
    //---------------------------------------------------------------------------
    //   produce event
    //---------------------------------------------------------------------------

    fout << "# event " << ie << "\n";
    
    //create thermal event
    std::vector<fastjet::PseudoJet> particlesBkg = thrm.createThermalEvent();

    for(fastjet::PseudoJet p : particlesBkg) {
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << 211 << " " << 1 << "\n";
      //fout << p.pt() << " " << p.rap() << " " << p.phi() << " " << p.m() << " " << 211 << " " << 1 << "\n"; 
    }
    fout << "end\n";
  }

  fout.close();
    
}
