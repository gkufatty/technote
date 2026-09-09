// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs_multiverse_filtered.C
//
// GENIE multiverse (every listed systematic thrown simultaneously, each
// drawn from an independent unit Gaussian per universe, universe 0 =
// nominal -- see GenieMultiverseSpectra ctor), restricted to the
// systematics that actually moved a diagnostic in the category-by-category
// study (make_neutron_spectra_systs_diagnostics_{high,medium,low}.C +
// plot_neutron_spectra_systs_diagnostics_ranking.C -- all three tiers now
// checked, LowPriority included). Modeled directly on
// ../../../spectrums/make_neutron_spectra_systs_multiverse_filtered.C, but
// with a whitelist re-derived from that category study instead of a single
// ntrue_per_event ranking.
//
// Threshold raised from the first cut of this list (chi2_shape>~1,
// "Moderate", 39 systs) to chi2_shape>5 ("Strong") in at least ONE
// diagnostic, across ALL 77 systematics in all three tiers -- a global
// re-ranking, not just a per-category one; see the per-systematic
// max-chi2-across-every-diagnostic table behind this cut for the full
// ranking. This keeps the ensemble to the systematics with an
// unambiguous, clearly-nonzero effect and drops the "Moderate but not
// clearly more than noise" middle ground, which contributed little to the
// multiverse's covariance while doubling its cost.
//
// Whitelist (26 of the 77 systematics across FSI+MEC+RES+QE+DIS+LowPriority),
// by category, with what was dropped and why -- see
// neutron_spectra_systs_diagnostics_*_metrics.csv for the actual
// per-systematic chi2 values behind this:
//   FSI (5/5 kept)  -- all five clear chi2>5 somewhere (KE, prongs, or the
//     N(primary)-N(visible) confusion diagnostic); hNFSI_MFP_2024 (179.5)
//     and FormZone2024 (77.3) dominate, the three FateFrac EV knobs trail
//     but still clear the bar (5.5-12.2).
//   MEC (6/6 kept)  -- all six clear chi2>5, cleanly split by nu/antinu
//     beam in q0 and/or final-state topology.
//   RES (6/10 kept) -- dropped RESvpvnRatioNuXSecSyst (max chi2 2.26, was
//     kept at the looser threshold), RESvpvnRatioNubarXSecSyst (0.30),
//     Theta_Delta2Npi (0.98), and MvNCRES (0.83) -- all four stayed well
//     below Strong despite Table 1 flagging the first two as
//     highest-priority for setting the target-nucleon p-vs-n ratio.
//   QE (3/7 kept)   -- kept only ZNormCCQE, RPAShapesupp2020, and
//     ZExpAxialFFSyst2020_EV1; dropped RPAShapeenh2020 (3.65),
//     ZExpAxialFFSyst2020_EV{2,3,4} (4.0, 1.2, 0.005). Unlike the
//     shape-only filtered multiverse in spectrums/, ZNormCCQE IS kept
//     here: that earlier list dropped it as "just a rate effect", but
//     this study's QE diagnostics explicitly target rate (Table 2: QE
//     sets the clean-neutron population size / overall CC rate), and
//     ZNormCCQE's ~16% rate swing on the clean-CCQE sample and ~6-8% swing
//     on the inclusive rate is exactly the effect being captured -- and it
//     clears chi2>5 outright (104.2) regardless.
//   DIS (6/40 kept) -- kept only the Q0/Q1-hadronic-energy knobs and the
//     four largest CC-multipion knobs (DISvnCC1pi_2020, DISvnCC2pi_2020,
//     DISvpCC2pi_2020, DISvpCC3pi_2020); the other ~9 soft-pion knobs that
//     cleared chi2>~1 at the looser threshold (DISvbarpCC1pi_2020,
//     DISvbarpCC2pi_2020, DISvbarnCC2pi_2020, DISvpNC1-3pi_2020,
//     DISvnNC2-3pi_2020, DISvnCC3pi_2020) didn't clear chi2>5, nor did any
//     of the remaining ~25 soft-pion/hadronization knobs.
//   LowPriority (0/9 kept) -- confirmed negligible as the technote expects
//     (no nucleon knockout / out-of-selection channel), with one
//     near-miss worth flagging: RDecBR1eta (Delta radiative eta-decay
//     branching) reaches chi2 3.64 (FHC, true-multiplicity sanity check)
//     -- real, but still below the chi2>5 bar, so it stays out.
//
// Produces multiverse spectra for 5 axes, FHC and RHC:
//   ntrue_per_event      -- true neutron prong count (post-FSI)
//   ngenie_per_event      -- N primary GENIE neutrons (pre-FSI)
//   nsel_per_event        -- reco NeutronicLikeStrict prong count
//   hadvise_q0            -- HadVisE [MeV]: the single most broadly
//                            sensitive diagnostic in the study (MEC/RES/QE);
//                            restricted to events with >=1 primary GENIE
//                            neutron AND a primary GENIE muon matching the
//                            beam (neutronMuonCut)
//   final_state_topology  -- 8-bin proton/neutron/pion breakdown: the other
//                            most sensitive diagnostic (MEC/RES/DIS); same
//                            neutronMuonCut restriction as hadvise_q0
// all thrown from the SAME set of universes (not split by category), so
// the resulting spread captures cross-category correlations -- that's the
// whole point of a multiverse over a per-syst +-1sigma scan.

#ifdef __CINT__
void make_neutron_spectra_systs_multiverse_filtered()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "../../../cuts.C"

// NeutronMultSystLists.cxx (included below) defines all five category
// getters unconditionally, so all six syst-family headers are needed
// regardless of which categories the whitelist below actually keeps.
#include "CAFAna/Systs/FSISysts.h"
#include "CAFAna/Systs/MECSysts.h"
#include "CAFAna/Systs/RESSysts.h"
#include "CAFAna/Systs/RPASysts.h"
#include "CAFAna/Systs/DISSysts.h"
#include "CAFAna/Systs/XSecSysts.h"  // GetGenieKnobSyst(), rwgt::
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/XSec/GenieMultiverseSyst.h"

#include "TFile.h"
#include "TDirectory.h"
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// cafana weights includes
#include "CAFAna/Weights/GenieWeights.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"
using namespace ana;

namespace
{
#include "NeutronMultSystLists.cxx"
#include "NeutronMultSystDiagnostics.cxx"
// kRelevantSysts (see header comment above for the per-category derivation
// and what was dropped) is shared with plot_relevant_systs_ranking.C so the
// two can't drift apart.
#include "RelevantSystsWhitelist.cxx"

  std::vector<const ana::ISyst*> FilterRelevant(const std::vector<const ana::ISyst*>& all)
  {
    std::vector<const ana::ISyst*> out;
    std::vector<std::string> dropped;
    for (const ana::ISyst* s : all) {
      if (kRelevantSysts.count(s->ShortName())) out.push_back(s);
      else dropped.push_back(s->ShortName());
    }
    std::cout << "Filtered " << all.size() << " -> " << out.size()
              << " systematics (" << dropped.size() << " dropped as negligible).\n";
    // Sanity check: warn (don't silently swallow) if the whitelist references
    // a name the category getters didn't actually provide -- usually means
    // the syst lists changed since this whitelist was derived.
    if (out.size() != kRelevantSysts.size()) {
      std::cout << "[WARN] " << kRelevantSysts.size() << " names were whitelisted but only "
                << out.size() << " were found -- the syst lists may have changed; "
                << "re-derive kRelevantSysts from a fresh diagnostics ranking.\n";
    }
    return out;
  }

  void ProcessBeam(const std::string& defname, const std::string& beam,
                    const std::vector<const ana::ISyst*>& systs,
                    unsigned int nUniverses, TFile* fOut)
  {
    // hadvise_q0 and final_state_topology are restricted to events with at
    // least one primary (pre-FSI) GENIE neutron AND a genuine primary GENIE
    // muon for this beam (mu- for FHC, mu+ for RHC) -- see
    // MakeHasPrimaryNeutronCut()/MakeHasPrimaryMuonCut() in cuts.C.
    const ana::Cut neutronMuonCut =
      kNeutronMultBaseCut && MakeHasPrimaryNeutronCut() && MakeHasPrimaryMuonCut(beam);

    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    const HistAxis axisTrue(  "N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar());
    const HistAxis axisGENIE( "N GENIE neutrons per event",      kBinNSel, MakeNGENIENeutronsVar());
    const HistAxis axisNSel(  "N_{NLS} per event",                kBinNSel, MakeNLSCountVar(beam));
    const HistAxis axisQ0(    "HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV);
    const HistAxis axisTopo(  "Final-state topology bin", kBinFinalStateTopology, kFinalStateTopologyBinVar);

    // universe 0 = nominal (shift=kNoShift); universes 1..N-1 each get every
    // *relevant* syst in `systs` shifted by an independent Gaus(0,1) draw,
    // thrown together across all five categories in the same universes.
    GenieMultiverseSpectra multiverseTrue( nUniverses, loader, axisTrue,  kNeutronMultBaseCut, kNoShift, kUnweighted, systs);
    GenieMultiverseSpectra multiverseGENIE(nUniverses, loader, axisGENIE, kNeutronMultBaseCut, kNoShift, kUnweighted, systs);
    GenieMultiverseSpectra multiverseNSel( nUniverses, loader, axisNSel,  kNeutronMultBaseCut, kNoShift, kUnweighted, systs);
    GenieMultiverseSpectra multiverseQ0(   nUniverses, loader, axisQ0,    neutronMuonCut, kNoShift, kUnweighted, systs);
    GenieMultiverseSpectra multiverseTopo( nUniverses, loader, axisTopo,  neutronMuonCut, kNoShift, kUnweighted, systs);

    std::cout << "[" << beam << "] " << nUniverses
              << " GENIE multiverse universes registered (over " << systs.size()
              << " relevant systematics) for ntrue_per_event + ngenie_per_event + "
              << "nsel_per_event + hadvise_q0 + final_state_topology; "
              << "running loader.Go() ...\n";

    loader.Go();

    TDirectory* beamDir = GetOrMkdir(fOut, beam);
    TDirectory* topDir  = GetOrMkdir(beamDir, "relevant_multiverse");

    multiverseTrue .SaveTo(topDir, "ntrue_per_event_multiverse");
    multiverseGENIE.SaveTo(topDir, "ngenie_per_event_multiverse");
    multiverseNSel .SaveTo(topDir, "nsel_per_event_multiverse");
    multiverseQ0   .SaveTo(topDir, "hadvise_q0_multiverse");
    multiverseTopo .SaveTo(topDir, "final_state_topology_multiverse");

    std::cout << beam << "/relevant_multiverse done.\n";
  }
}

void make_neutron_spectra_systs_multiverse_filtered()
{
  // All three tiers (FSI+MEC+RES+QE+DIS+LowPriority) -- LowPriority is now
  // included in the filtering pass (confirmed negligible, see header
  // comment) rather than skipped outright.
  const std::vector<const ana::ISyst*> allSysts = getNeutronMultSysts_All();
  std::cout << "Loaded " << allSysts.size()
            << " GENIE xsec systematics from all three priority tiers "
               "(FSI+MEC+RES+QE+DIS+LowPriority)\n";

  const std::vector<const ana::ISyst*> systs = FilterRelevant(allSysts);

  const unsigned int nUniverses = 50;  // universe 0 = nominal; lower this for a quicker test

  TFile* fOut = new TFile("make_neutron_spectra_systs_multiverse_filtered.root", "RECREATE");

  ProcessBeam(kNeutronMultDefFHC, "FHC", systs, nUniverses, fOut);
  ProcessBeam(kNeutronMultDefRHC, "RHC", systs, nUniverses, fOut);

  fOut->Write();
  fOut->Close();
  delete fOut;
  std::cout << "Saved to make_neutron_spectra_systs_multiverse_filtered.root\n";
}

#endif
