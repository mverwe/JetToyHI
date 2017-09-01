#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include <fstream>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "include/pythiaEvent.hh"
#include "include/extraInfo.hh"

#include "PU14/CmdLine.hh"

using namespace std;
using namespace fastjet;

int main(int argc, char *argv[])
{
   if(argc != 2)
   {
      cerr << "Usage: " << argv[0] << " RootFileName" << endl;
      return -1;
   }

   TFile File(argv[1]);

   TTree *T = (TTree *)File.Get("Tdata");

   vector<double> *px = NULL;   T->SetBranchAddress("px", &px);
   vector<double> *py = NULL;   T->SetBranchAddress("py", &py);
   vector<double> *pz = NULL;   T->SetBranchAddress("pz", &pz);
   vector<double> *m = NULL;    T->SetBranchAddress("m", &m);
   vector<int> *id = NULL;      T->SetBranchAddress("id", &id);
   double w;                    T->SetBranchAddress("xSec", &w);

   ProgressBar Bar(cerr, T->GetEntries());
   Bar.SetStyle(3);

   //output text file
   std::ofstream fout;
   TString outFileName = Form("%s",argv[1]);
   outFileName.ReplaceAll("pthat120","pu14/pthat120");
   outFileName.ReplaceAll(".root",".pu14");
   std::cout << "outFileName: " << outFileName << std::endl;

   fout.open(outFileName.Data());
   
   int EntryCount = T->GetEntries();
   for(int iE = 0; iE < EntryCount; iE++)
   {
      T->GetEntry(iE);

      Bar.Update(iE);
      Bar.PrintWithMod(EntryCount / 300);

      //cout << px << " " << id << endl;

      fout << "# event " << iE << "\n";
      fout << "weight " << w <<  "\n";
      for(int i = 0; i < (int)id->size(); i++)
        fout << (*px)[i] << " " << (*py)[i] << " " << (*pz)[i] << " " << (*m)[i]
             << " " << (*id)[i] << " " << 0 <<  "\n";
      fout << "end" << "\n";
   }

   fout.close();

   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   File.Close();


   return 0;
}




