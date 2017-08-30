#ifndef jewelMatcher_h
#define jewelMatcher_h

#include <iostream>
#include <vector>

#include "fastjet/PseudoJet.hh"

fastjet::PseudoJet GetCorrection(std::vector<fastjet::PseudoJet> Constituents, std::vector<fastjet::PseudoJet> ThermalParticles);
fastjet::PseudoJet GetJetCorrection(fastjet::PseudoJet Jet, std::vector<fastjet::PseudoJet> ThermalParticles);
fastjet::PseudoJet GetCorrectedJet(fastjet::PseudoJet Jet, std::vector<fastjet::PseudoJet> ThermalParticles);
std::vector<fastjet::PseudoJet> GetCorrectedJets(std::vector<fastjet::PseudoJet> Jets, std::vector<fastjet::PseudoJet> ThermalParticles);

fastjet::PseudoJet GetCorrection(std::vector<fastjet::PseudoJet> Constituents, std::vector<fastjet::PseudoJet> ThermalParticles)
{
   fastjet::PseudoJet Correction;

   for(fastjet::PseudoJet p : Constituents)
   {
      if(p.E() > 0.01)   // definitely not a dummy - this should speed things up a lot
         continue;

      for(fastjet::PseudoJet &j : ThermalParticles)
      {
         if(p.squared_distance(j) > 1e-8)
            continue;

         Correction = Correction + j;
         j.reset(0, 0, 0, 0);
      }
   }

   return Correction;
}

fastjet::PseudoJet GetJetCorrection(fastjet::PseudoJet Jet, std::vector<fastjet::PseudoJet> ThermalParticles)
{
   fastjet::PseudoJet Correction = GetCorrection(Jet.constituents(), ThermalParticles);
   return Correction;
}

fastjet::PseudoJet GetCorrectedJet(fastjet::PseudoJet Jet, std::vector<fastjet::PseudoJet> ThermalParticles)
{
   fastjet::PseudoJet Correction = GetCorrection(Jet.constituents(), ThermalParticles);
   return (Jet - Correction);
}

std::vector<fastjet::PseudoJet> GetCorrectedJets(std::vector<fastjet::PseudoJet> Jets, std::vector<fastjet::PseudoJet> ThermalParticles)
{
   std::vector<fastjet::PseudoJet> Result;
   for(fastjet::PseudoJet j : Jets)
      Result.push_back(GetCorrectedJet(j, ThermalParticles));
   return Result;
}

#endif
