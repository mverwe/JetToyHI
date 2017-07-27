#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "include/thermalEvent.hh"

#include "PU14/CmdLine.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv)
{
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  CmdLine cmdline(argc,argv);
  // inputs read from command line
  unsigned int nEvent = cmdline.value<unsigned int>("-nev",1);  // first argument: command line option; second argument: default value
  
  //event generator settings
  //double       ptAve = cmdline.value<double>("-ptAve",1.2);
  //unsigned int mult  = cmdline.value<unsigned int>("-mult",7000);

  int centBin = cmdline.value<int>("-ncent",0);  // first argument: command line option; second argument: default value
  if(centBin>3) {
    std::cout << "provided centBin too large (centBin=" << centBin << ")" << std::endl;
    return -1;
  }
  unsigned int multArr[4] = {7000,4500,1700,100};
  double ptAveArr[4] = {1.2,1.0,0.9,0.85}; //for central: 1.2. peripheral: 1.

  unsigned int mult = multArr[centBin];
  double ptAve = ptAveArr[centBin];
  
  std::cout << "generating " << nEvent << " events with ptAve = " << ptAve << " and multiplicity = " << mult << std::endl;
  
  //unsigned int mult = 12000;
  //double       ptAve = 0.7;
  thermalEvent thrm(mult,ptAve, -3.0, 3.0, 0.5);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  //output text file
  ofstream fout;
  const char *dir = getenv("PWD");//"/eos/user/m/mverweij/JetWorkshop2017/samples/";
  int jobId = cmdline.value<int>("-jobId",0); 
  TString outFileName = Form("%s/ThermalEventsMult%dPtAv%.2f_%d.pu14",dir,mult,ptAve,jobId);
  
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
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << pdgid << " " << vtx << "\n";
    }
    fout << "end\n";
  }

  fout.close();

  std::cout <<  std::endl;
    
}
