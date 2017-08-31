#ifndef jewelMatcher_h
#define jewelMatcher_h

#include <iostream>
#include <vector>
#include <map>

#include "fastjet/PseudoJet.hh"

fastjet::PseudoJet GetCorrection(std::vector<fastjet::PseudoJet> Constituents, std::vector<fastjet::PseudoJet> ThermalParticles);
fastjet::PseudoJet GetJetCorrection(fastjet::PseudoJet Jet, std::vector<fastjet::PseudoJet> ThermalParticles);
fastjet::PseudoJet GetCorrectedJet(fastjet::PseudoJet Jet, std::vector<fastjet::PseudoJet> ThermalParticles);
std::vector<fastjet::PseudoJet> GetCorrectedJets(std::vector<fastjet::PseudoJet> Jets, std::vector<fastjet::PseudoJet> ThermalParticles);
std::vector<fastjet::PseudoJet> GetCorrectedJets(std::vector<std::vector<fastjet::PseudoJet>> Jets, std::vector<fastjet::PseudoJet> ThermalParticles);
std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> GetCorrectedSubJets(std::vector<fastjet::PseudoJet> Jets, std::vector<fastjet::PseudoJet> ThermalParticles);
std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> GetCorrectedSubJets(std::vector<std::vector<fastjet::PseudoJet>> Jets1, std::vector<std::vector<fastjet::PseudoJet>> Jets2, std::vector<fastjet::PseudoJet> ThermalParticles);
std::vector<double> CalculateDR(std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> SubJets);
std::vector<double> CalculateZG(std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> SubJets);

fastjet::PseudoJet GetCorrection(std::vector<fastjet::PseudoJet> Constituents, std::vector<fastjet::PseudoJet> ThermalParticles)
{
   fastjet::PseudoJet Correction;

   for(fastjet::PseudoJet p : Constituents)
   {
   

     if(p.E() > 0.01)   // definitely not a dummy - this should speed things up a lot
       continue;

      for(fastjet::PseudoJet j : ThermalParticles)
      {
         double deltaR = std::sqrt((p.eta() -  j.eta())*(p.eta() - j.eta()) + (p.delta_phi_to(j))*(p.delta_phi_to(j)));
         
         if(deltaR > 1e-5) //1e-8)
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
   for(fastjet::PseudoJet &j : Jets)
      Result.push_back(GetCorrectedJet(j, ThermalParticles));
   return Result;
}

std::vector<fastjet::PseudoJet> GetCorrectedJets(std::vector<std::vector<fastjet::PseudoJet>> Jets, std::vector<fastjet::PseudoJet> ThermalParticles)
{
   std::vector<fastjet::PseudoJet> Result;
   for(std::vector<fastjet::PseudoJet> j : Jets)
      Result.push_back(join(j) - GetCorrection(j, ThermalParticles));
   return Result;
}

std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> GetCorrectedSubJets(std::vector<fastjet::PseudoJet> Jets, std::vector<fastjet::PseudoJet> ThermalParticles)
{
   std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> Result;
   for(fastjet::PseudoJet j : Jets)
   {
      fastjet::PseudoJet j1, j2;
      if(j.has_parents(j1, j2))
      {
         j1 = GetCorrectedJet(j1, ThermalParticles);
         j2 = GetCorrectedJet(j2, ThermalParticles);
      }
      else
      {
         j1.reset(0, 0, 0, 0);
         j2.reset(0, 0, 0, 0);
      }
      Result.emplace_back(j1, j2);
   }
   return Result;
}

std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> GetCorrectedSubJets(std::vector<std::vector<fastjet::PseudoJet>> Jets1, std::vector<std::vector<fastjet::PseudoJet>> Jets2, std::vector<fastjet::PseudoJet> ThermalParticles)
{
   std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> Result;
   for(int i = 0; i < (int)Jets1.size(); i++)
   {
      fastjet::PseudoJet j1 = join(Jets1[i]) - GetCorrection(Jets1[i], ThermalParticles);
      fastjet::PseudoJet j2 = join(Jets2[i]) - GetCorrection(Jets2[i], ThermalParticles);
      Result.emplace_back(j1, j2);
   }
   return Result;
}

std::vector<double> CalculateDR(std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> SubJets)
{
   std::vector<double> Result;
   for(std::pair<fastjet::PseudoJet, fastjet::PseudoJet> p : SubJets)
      Result.push_back(std::sqrt(p.first.squared_distance(p.second)));
   return Result;
}

std::vector<double> CalculateZG(std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> SubJets)
{
   std::vector<double> Result;
   for(std::pair<fastjet::PseudoJet, fastjet::PseudoJet> p : SubJets)
   {
      double pt1 = p.first.perp();
      double pt2 = p.second.perp();
      double zg = -1;
      if(pt1 + pt2 > 0)
         zg = min(pt1, pt2) / (pt1 + pt2);
      Result.push_back(zg);
   }
   return Result;
}

#endif
