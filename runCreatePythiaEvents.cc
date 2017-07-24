#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "include/pythiaEvent.hh"
#include "include/extraInfo.hh"

using namespace std;
using namespace fastjet;

int main ()
{
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);
  
  // Number of events, generated and listed ones.
  unsigned int nEvent    = 100;

  //event generator settings
  double       ptHat = 120.;
  unsigned int tune  = 14;
  pythiaEvent pyt(ptHat, tune, -3.0, 3.0);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  //output text file
  ofstream fout;
  const char *dir = getenv("PWD");//"/eos/user/m/mverweij/JetWorkshop2017/samples/";
  TString outFileName = Form("%s/PythiaEventsTune%dPtHat%.0f.pu14",dir,tune,ptHat);
  
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

    //create pythia event
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();
   
    for(fastjet::PseudoJet p : particlesSig) {
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << pdgid << " " << vtx << "\n";
      //fout << p.pt() << " " << p.rap() << " " << p.phi() << " " << p.m() << " " << pdgid << " " << vtx << "\n"; 
    }
    fout << "end\n";
  }

  fout.close();
    
}
