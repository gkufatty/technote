// check_final_states.C
// Answers: for events passing eventLevelCut with >=1 true primary neutron,
// what are the most common true final states (by primary particle content)?
//
// Motivation: the neutron/pion topology cuts in NeutronMultiplicityCuts.h
// (kOneMuon1NeutronOnly, kOneMuonNNeutronOnly, kOneMuon1pNNeutronsOnly,
// kOneMuonNpNNeutronsOnly) all reject any event with a primary proton
// (hasOtherPrimary=true the moment pdg==2212 is seen). Since numu(bar) CC
// interactions on nuclear targets very commonly produce a primary proton,
// this macro checks directly, rather than assuming, what fraction of
// >=1-neutron events those four cuts actually cover -- and what the
// proton-inclusive final states (uncovered by any of the four) look like.
//
// Each event's final state is summarized as counts of:
//   mu (|pdg|==13), p (pdg==2212), n (pdg==2112),
//   pi (|pdg|==211 or pdg==111), other (anything else)
// e.g. "mu1_p1_n1_pi0_other0"
//
// Run (needs --numubarccinc so the linker pulls in the compiled
// ana::xsec::numubarcc::* cuts from NDAna/numubarcc_inc -- see the same
// flag on make_mc_spectra.C's own Run: line):
//   cafe -bq --numubarccinc check_final_states.C

#ifdef __CINT__
void check_final_states()
{
  std::cout << "Sorry, you must run in compiled mode" << std::endl;
}
#else

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "3FlavorAna/Cuts/NumuCuts2024.h"

#include "StandardRecord/Proxy/SRProxy.h"
#include "NDAna/numubarcc_inc/NumubarCCIncCuts.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>

using namespace ana;

namespace
{
  std::string FinalStateKey(const caf::SRProxy* sr)
  {
    int nMu = 0, nP = 0, nN = 0, nPi = 0, nOther = 0;
    for (const auto& prim : sr->mc.nu[0].prim) {
      const int pdg = prim.pdg;
      if      (std::abs(pdg) == 13)                          ++nMu;
      else if (pdg == 2212)                                  ++nP;
      else if (pdg == 2112)                                  ++nN;
      else if (std::abs(pdg) == 211 || pdg == 111)            ++nPi;
      else                                                    ++nOther;
    }
    return "mu" + std::to_string(nMu) + "_p" + std::to_string(nP)
         + "_n" + std::to_string(nN) + "_pi" + std::to_string(nPi)
         + "_other" + std::to_string(nOther);
  }

  // Runs one beam, tallying final-state signatures for events with
  // eventLevelCut && >=1 primary neutron. Uses a side-effect Var (the
  // standard CAFAna trick for per-event diagnostics that don't fit a
  // histogram axis) to fill `counts` as the loader runs.
  void CountFinalStates(const std::string& defname,
                         const ana::Cut& eventLevelCut,
                         std::map<std::string, long>& counts,
                         long& nTotal, long& nWithNeutron)
  {
    SpectrumLoader loader(defname);

    const ana::Cut hasPrimaryNeutron([](const caf::SRProxy* sr) {
      if (!sr || sr->mc.nnu == 0) return false;
      for (const auto& prim : sr->mc.nu[0].prim)
        if (prim.pdg == 2112) return true;
      return false;
    });

    const ana::Var kTally([&](const caf::SRProxy* sr) -> double {
      ++nTotal;
      if (sr->mc.nnu == 0) return 0.0;
      bool hasNeutron = false;
      for (const auto& prim : sr->mc.nu[0].prim)
        if (prim.pdg == 2112) { hasNeutron = true; break; }
      if (hasNeutron) {
        ++nWithNeutron;
        ++counts[FinalStateKey(sr)];
      }
      return 0.0;
    });

    HistAxis dummyAxis("dummy", Binning::Simple(1, 0, 1), kTally);
    Spectrum sTally(loader, dummyAxis, eventLevelCut, kNoShift, kUnweighted);

    loader.Go();
  }

  void PrintTopFinalStates(const std::string& beam,
                            const std::map<std::string, long>& counts,
                            long nWithNeutron, int topN = 15)
  {
    std::vector<std::pair<std::string, long>> sorted(counts.begin(), counts.end());
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    std::cout << "\n=== " << beam << ": most common final states "
              << "(events with >=1 true primary neutron, N = " << nWithNeutron << ") ===\n";
    std::cout << std::left << std::setw(28) << "final state"
              << std::right << std::setw(10) << "count"
              << std::setw(10) << "%\n";

    long shown = 0;
    for (int i = 0; i < topN && i < (int)sorted.size(); ++i) {
      const double pct = 100.0 * sorted[i].second / nWithNeutron;
      std::cout << std::left << std::setw(28) << sorted[i].first
                << std::right << std::setw(10) << sorted[i].second
                << std::setw(9) << std::fixed << std::setprecision(2) << pct << "%\n";
      shown += sorted[i].second;
    }
    std::cout << std::left << std::setw(28) << "(all others)"
              << std::right << std::setw(10) << (nWithNeutron - shown)
              << std::setw(9) << std::fixed << std::setprecision(2)
              << 100.0 * (nWithNeutron - shown) / nWithNeutron << "%\n";

    // Cross-check against the four existing topology cuts in
    // NeutronMultiplicityCuts.h, which all require zero primary protons.
    long noProtonNoOther = 0;
    for (const auto& kv : counts) {
      // key format: mu<M>_p<P>_n<N>_pi<PI>_other<O>
      if (kv.first.find("_p0_") != std::string::npos &&
          kv.first.find("_other0") != std::string::npos)
        noProtonNoOther += kv.second;
    }
    std::cout << "-- events with zero primary protons AND zero 'other' primaries "
              << "(i.e. covered by kOneMuon*NeutronOnly/kOneMuon*pnN cuts, "
              << "modulo their nMu==1 requirement): "
              << noProtonNoOther << " / " << nWithNeutron << " = "
              << std::fixed << std::setprecision(2)
              << 100.0 * noProtonNoOther / nWithNeutron << "%\n";
  }

  // Parses a "mu<M>_p<P>_n<N>_pi<PI>_other<O>" key back into its five counts.
  void ParseFinalStateKey(const std::string& key, int& nMu, int& nP, int& nN, int& nPi, int& nOther)
  {
    sscanf(key.c_str(), "mu%d_p%d_n%d_pi%d_other%d", &nMu, &nP, &nN, &nPi, &nOther);
  }

  std::string BinLabel(int n, int cap, const char* prefix)
  {
    return n >= cap ? (std::string(prefix) + std::to_string(cap) + "+")
                     : (std::string(prefix) + std::to_string(n));
  }
y
  // Groups the full final-state breakdown (not just the top N) by proton
  // count x pion count, so the "(all others)" tail from PrintTopFinalStates
  // becomes a small, readable table instead of hundreds of <1.5% rows.
  // Binned by nProtons/nPions only -- nMu and nOther aren't split out here,
  // so every event with >=1 primary neutron lands in exactly one cell and
  // the table sums to exactly nWithNeutron. The nProtons==0 row (all pion
  // bins) is exactly the union of the four kOneMuon*Only/kOneMuon*pnN cuts,
  // modulo their nMu==1 requirement -- same check as noProtonNoOther above,
  // now visible as a single row rather than a separate cross-check line.
  void PrintFSIGrouping(const std::string& beam,
                         const std::map<std::string, long>& counts,
                         long nWithNeutron)
  {
    const int kProtonCap = 5; // bins: 0,1,2,3,4,5+
    const int kPionCap   = 3; // bins: 0,1,2,3+

    std::map<std::pair<int,int>, long> table; // (protonBin, pionBin) -> count
    std::vector<long> rowTotal(kProtonCap + 1, 0);
    std::vector<long> colTotal(kPionCap + 1, 0);

    for (const auto& kv : counts) {
      int nMu, nP, nN, nPi, nOther;
      ParseFinalStateKey(kv.first, nMu, nP, nN, nPi, nOther);
      const int pBin = std::min(nP,  kProtonCap);
      const int piBin = std::min(nPi, kPionCap);
      table[{pBin, piBin}] += kv.second;
      rowTotal[pBin]  += kv.second;
      colTotal[piBin] += kv.second;
    }

    std::cout << "\n=== " << beam << ": FSI grouping by N(proton) x N(pion) "
              << "(events with >=1 true primary neutron, N = " << nWithNeutron << ") ===\n";

    std::cout << std::left << std::setw(10) << "protons";
    for (int piBin = 0; piBin <= kPionCap; ++piBin)
      std::cout << std::right << std::setw(14) << BinLabel(piBin, kPionCap, "pi=");
    std::cout << std::right << std::setw(16) << "row total" << "\n";

    for (int pBin = 0; pBin <= kProtonCap; ++pBin) {
      std::cout << std::left << std::setw(10) << BinLabel(pBin, kProtonCap, "");
      for (int piBin = 0; piBin <= kPionCap; ++piBin) {
        const long c = table.count({pBin, piBin}) ? table.at({pBin, piBin}) : 0;
        const double pct = 100.0 * c / nWithNeutron;
        std::cout << std::right << std::setw(9) << c
                  << " (" << std::fixed << std::setprecision(1) << std::setw(2) << pct << "%)";
      }
      const double rowPct = 100.0 * rowTotal[pBin] / nWithNeutron;
      std::cout << std::right << std::setw(10) << rowTotal[pBin]
                << " (" << std::fixed << std::setprecision(1) << rowPct << "%)\n";
    }

    std::cout << std::left << std::setw(10) << "col total";
    for (int piBin = 0; piBin <= kPionCap; ++piBin) {
      const double colPct = 100.0 * colTotal[piBin] / nWithNeutron;
      std::cout << std::right << std::setw(9) << colTotal[piBin]
                << " (" << std::fixed << std::setprecision(1) << colPct << "%)";
    }
    std::cout << "\n";
  }
}

void check_final_states()
{
  const ana::Cut eventLevelCut =
      ana::xsec::numubarcc::kQualityCut &&
      ana::xsec::numubarcc::kContainmentCut &&
      ana::xsec::numubarcc::kRecoVtxNumuFiducialCut &&
      ana::xsec::numubarcc::kMuonIDCut;

  
  const std::string defFHC= "prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_filematchedSystematics_nominal_v1";
  const std::string defRHC = "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_ndphysics_systs_filematch_v1";

  std::map<std::string, long> countsFHC, countsRHC;
  long nTotalFHC = 0, nWithNeutronFHC = 0;
  long nTotalRHC = 0, nWithNeutronRHC = 0;

  CountFinalStates(defFHC, eventLevelCut, countsFHC, nTotalFHC, nWithNeutronFHC);
  CountFinalStates(defRHC, eventLevelCut, countsRHC, nTotalRHC, nWithNeutronRHC);

  PrintTopFinalStates("FHC", countsFHC, nWithNeutronFHC);
  PrintTopFinalStates("RHC", countsRHC, nWithNeutronRHC);

  PrintFSIGrouping("FHC", countsFHC, nWithNeutronFHC);
  PrintFSIGrouping("RHC", countsRHC, nWithNeutronRHC);
}

#endif
