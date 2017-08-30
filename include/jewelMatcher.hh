#ifndef jewelMatcher_h
#define jewelMatcher_h

#include <iostream>

#include "fastjet/PseudoJet.hh"

fastjet::PseudoJet GetCorrection(std::vector<fastjet::PseudoJet> &Constituents, std::vector<fastjet::PsuedoJet> ThermalParticles);

fastjet::PseudoJet GetCorrection(std::vector<fastjet::PseudoJet> &Constituents, std::vector<fastjet::PsuedoJet> ThermalParticles)
{
   fastjet::PseudoJet Correction;

   for(fastjet::PseudoJet p : Constituents)
   {
      if(p.E() > 0.0001)   // definitely not a dummy - this should speed things up a lot
         continue;

      for(fastjet::PseudoJet &j : ThermalParticles)
      {
         if(p.squared_distance(j) > 1e-10)
            continue;

         Correction = Correction + j;
         j.reset(0, 0, 0, 0);
      }
   }

   return Correction;
}

#endif
