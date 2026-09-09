// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs_multiverse_filtered.C
//
// q0lt100MeV variant of ../genie/make_neutron_spectra_systs_multiverse_filtered.C:
// GENIE multiverse (every listed systematic thrown simultaneously, each
// drawn from an independent unit Gaussian per universe, universe 0 =
// nominal), restricted to HadVisE < 100 MeV (Q0Lt100MeVCommon.cxx) and to
// the systematics whitelisted in Q0Lt100MeVRelevantSystsWhitelist.cxx --
// re-derived from this folder's own diagnostics study
// (make_neutron_spectra_systs_diagnostics_{high,medium,low}.C +
// plot_neutron_spectra_systs_diagnostics_ranking.C), NOT the same 26
// systematics as ../genie/'s inclusive multiverse.
//
// Same chi2_shape > 5 ("Strong") threshold as ../genie/, globally re-ranked
// across all 77 systematics in getNeutronMultSysts_All(). See
// Q0Lt100MeVRelevantSystsWhitelist.cxx for the full per-category
// derivation (11 of 77 kept, down from 26 inclusively) and
// neutron_spectra_systs_diagnostics_*_metrics.csv for the per-systematic
// chi2 values behind it. Headline result: EVERY FSI systematic (the
// dominant category inclusively -- hNFSI_MFP_2024 alone reached chi2=179.5)
// drops below threshold once HadVisE is cut below 100 MeV, since FSI
// rescattering/absorption mostly acts by adding hadronic activity, i.e. by
// pushing events to HIGHER HadVisE -- exactly the phase space this cut
// removes. MEC/RES/QE/DIS each keep roughly half their inclusive members;
// none newly appear that weren't already relevant inclusively.
//
// Produces multiverse spectra for the SAME 5 axes as ../genie/, FHC and RHC,
// all restricted to HadVisE < 100 MeV:
//   ntrue_per_event      -- true neutron prong count (post-FSI)
//   ngenie_per_event      -- N primary GENIE neutrons (pre-FSI)
//   nsel_per_event        -- reco NeutronicLikeStrict prong count
//   hadvise_q0            -- HadVisE [MeV], binned 0-100 MeV (kBinHadVisMeV100
//                            -- binning it to 300 MeV while cutting at 100
//                            would leave the upper bins permanently empty);
//                            restricted to events with >=1 primary GENIE
//                            neutron AND a primary GENIE muon matching the
//                            beam (neutronMuonCut), same as ../genie/
//   final_state_topology  -- 8-bin proton/neutron/pion breakdown, same
//                            neutronMuonCut restriction as hadvise_q0
// all thrown from the SAME set of universes (not split by category), so the
// resulting spread captures cross-category correlations.

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
#include "../genie/NeutronMultSystLists.cxx"
#include "../genie/NeutronMultSystDiagnostics.cxx"
#include "Q0Lt100MeVCommon.cxx"
// kRelevantSysts (see header comment above and this file for the
// per-category derivation) is shared with any ranking plot in this folder
// so the two can't drift apart.
#include "Q0Lt100MeVRelevantSystsWhitelist.cxx"

  std::vector<const ana::ISyst*> FilterRelevant(const std::vector<const ana::ISyst*>& all)
  {
    std::vector<const ana::ISyst*> out;
    std::vector<std::string> dropped;
    for (const ana::ISyst* s : all) {
      if (kRelevantSysts.count(s->ShortName())) out.push_back(s);
      else dropped.push_back(s->ShortName());
    }
    std::cout << "Filtered " << all.size() << " -> " << out.size()
              << " systematics (" << dropped.size() << " dropped as negligible "
              << "for HadVisE < 100 MeV).\n";
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
    // muon for this beam (mu- for FHC, mu+ for RHC), on top of the same
    // HadVisE < 100 MeV cut every axis here uses -- see
    // MakeHasPrimaryNeutronCut()/MakeHasPrimaryMuonCut() in cuts.C and
    // kQ0Lt100MeVCut in Q0Lt100MeVCommon.cxx.
    const ana::Cut neutronMuonCut =
      kQ0Lt100MeVCut && MakeHasPrimaryNeutronCut() && MakeHasPrimaryMuonCut(beam);

    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    const HistAxis axisTrue(  "N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar());
    const HistAxis axisGENIE( "N GENIE neutrons per event",      kBinNSel, MakeNGENIENeutronsVar());
    const HistAxis axisNSel(  "N_{NLS} per event",                kBinNSel, MakeNLSCountVar(beam));
    const HistAxis axisQ0(    "HadVisE [MeV]", kBinHadVisMeV100, kNumuHadVisMeV);
    const HistAxis axisTopo(  "Final-state topology bin", kBinFinalStateTopology, kFinalStateTopologyBinVar);

    // universe 0 = nominal (shift=kNoShift); universes 1..N-1 each get every
    // *relevant* syst in `systs` shifted by an independent Gaus(0,1) draw,
    // thrown together across all categories in the same universes.
    GenieMultiverseSpectra multiverseTrue( nUniverses, loader, axisTrue,  kQ0Lt100MeVCut, kNoShift, kUnweighted, systs);
    GenieMultiverseSpectra multiverseGENIE(nUniverses, loader, axisGENIE, kQ0Lt100MeVCut, kNoShift, kUnweighted, systs);
    GenieMultiverseSpectra multiverseNSel( nUniverses, loader, axisNSel,  kQ0Lt100MeVCut, kNoShift, kUnweighted, systs);
    GenieMultiverseSpectra multiverseQ0(   nUniverses, loader, axisQ0,    neutronMuonCut, kNoShift, kUnweighted, systs);
    GenieMultiverseSpectra multiverseTopo( nUniverses, loader, axisTopo,  neutronMuonCut, kNoShift, kUnweighted, systs);

    std::cout << "[" << beam << "] " << nUniverses
              << " GENIE multiverse universes registered (over " << systs.size()
              << " relevant systematics, HadVisE < 100 MeV) for ntrue_per_event + "
              << "ngenie_per_event + nsel_per_event + hadvise_q0 + "
              << "final_state_topology; running loader.Go() ...\n";

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
  // All three tiers (FSI+MEC+RES+QE+DIS+LowPriority) -- FilterRelevant()
  // below drops FSI and LowPriority to nothing for this HadVisE < 100 MeV
  // phase space, but they're still loaded here so the filter can confirm
  // that (see the [WARN] sanity check).
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
