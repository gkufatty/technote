// GENIE_vs_GiBUU_common.h
//
// Shared between make_GENIE_vs_GiBUU_histos.C (reads the raw "Convenient
// files" and fills histograms) and plot_GENIE_vs_GiBUU.C (reads the ROOT
// file the former writes and makes the PDF) -- split out of the original
// single-file GENIE_vs_GiBUU_simplified.C so the plotting half can be
// iterated on without re-running over the raw data every time.
//
// Constants, the tune list, and the GenHistos/MultVsXProfiles struct
// layouts live here since both scripts need to agree on them: the "make"
// script uses them to build+fill+name histograms, the "plot" script uses
// the exact same names to Get() them back out, and the exact same struct
// layouts so drawComparePad()'s `TH1D* GenHistos::*member` pointers-to-
// member (defined in the plot script) type-check against a GenHistos this
// header defines once, not two independently-maintained copies.
#pragma once

#include "TH1D.h"
#include "TProfile.h"
#include "TString.h"
#include <string>
#include <vector>

// ── kinematic index convention (see GENIE_comparison.C) ─────────────────────
const int kIdxE = 0;
const double kMinParticleE = 0.005;  // GeV -- KE floor for the multiplicity counts
const double kEavailCut    = 0.300;  // GeV -- E_avail selection threshold

// binning for the neutron kinetic-energy-distribution plot (MeV, KE = E - m;
// no kMinParticleE floor -- every final-state neutron is filled regardless
// of whether it'd pass the 5 MeV cut used for the multiplicity counts)
const int    kNNeutronKEBins = 20;
const double kNeutronKELo    = 0.0;
const double kNeutronKEHi    = 200.0;

// binning for the E_avail distribution plot (MeV) -- the selection cut
// itself (kEavailCut above) stays in GeV; only this display histogram's
// fill value and axis are converted to MeV.
const int    kNEavailBins = 10;
const double kEavailLo    = 0.0;
const double kEavailHi    = 50.0;

// ── particle masses (GeV) ────────────────────────────────────────────────────
const double kMassProton   = 0.938272;
const double kMassNeutron  = 0.939565;
const double kMassPiC      = 0.139570;
const double kMassPi0      = 0.134977;

// species colors for the avg-multiplicity-vs-Enu plot (color = species,
// line style = tune -- see drawMultiplicityVsXPad() in the plot script)
const int kColorProton    = kRed+1;
const int kColorNeutron   = kBlue+1;
const int kColorPiCharged = kOrange+7;
const int kColorPi0       = kGreen+2;

// Enu (true incoming neutrino energy) binning for the average-FS-
// multiplicity-vs-Enu profile plot -- same variable-width edges used for
// the E_nu histogram in GENIE_comparison.C / colin_numubarCCinc_xsec.C.
const std::vector<double> kEnuEdges =
  {0.0, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0, 2.50, 3.0, 4.0};
const int kNEnuBins = (int)kEnuEdges.size() - 1;

// multiplicity histogram bins {0, 1, 2, "3+"} -- used by makeCountHist() in
// the make script; the plot script only needs it for Get()-back consistency
const int kNCountBins = 4;

// ── generator configurations ───────────────────────────────────────────────
const std::string kBase     = "/exp/nova/data/groups/nd/convenient/ConvenientOutputs_NOvA";
const std::string kFluxHead = "DeriveFlux50MeVwidth_ppfx";
const std::string kOutDir   = "G18_10j_AR23_GiBUU2025";

struct GenConfig {
  std::string label;
  std::string fhcBase;  // .../GENIE/<version_tune>/FHCDeriveFlux50MeVwidth_ppfx
  std::string rhcBase;
  int color;
  // per-flavor subdirectory names -- GENIE Convenient files use
  // "numu(bar)_only", GiBUU's use "numu(bar)CC_only"
  std::string numuSubdir    = "numu_only";
  std::string numubarSubdir = "numubar_only";
};

const std::vector<GenConfig> kGens = {
  {"GENIE hN18",
   kBase+"/GENIE/v3.00.06_G18_10j_00_000/FHC"+kFluxHead,
   kBase+"/GENIE/v3.00.06_G18_10j_00_000/RHC"+kFluxHead,
   kCyan+2},
  {"GENIE hA18",
   kBase+"/GENIE/v3.04.00_AR23_20i_00_000/FHC"+kFluxHead,
   kBase+"/GENIE/v3.04.00_AR23_20i_00_000/RHC"+kFluxHead,
   kRed+1},
  {"GiBUU R2025",
   kBase+"/GiBUU/R2025_03_Defaultparams2025/FHC"+kFluxHead,
   kBase+"/GiBUU/R2025_03_Defaultparams2025/RHC"+kFluxHead,
   kAzure-2, "numuCC_only", "numubarCC_only"},
};

// kGens indices, named -- used to pick out GiBUU R2025 vs GENIE hN18 only
// for the avg-multiplicity-vs-Enu plot (unlike the count/KE/E_avail plots,
// which overlay all three tunes).
enum { kIdxHN18 = 0, kIdxHA18 = 1, kIdxGiBUU = 2 };

// Four-species TProfile set for "average FS multiplicity vs. Enu" (proton,
// neutron, pi_charged, pi0 -- no "other", unlike GENIE_comparison.C's
// five-species version). Unrestricted by the E_avail cut (unlike every
// other histogram here), and not unit-area normalized -- each point is a
// genuine average multiplicity via TProfile, not a shape.
struct MultVsXProfiles {
  TProfile* proton    = nullptr;
  TProfile* neutron   = nullptr;
  TProfile* piCharged = nullptr;
  TProfile* pi0       = nullptr;
};

// ── per-generator results ──────────────────────────────────────────────────
struct GenHistos {
  TH1D* hNeutron       = nullptr;  // N_neutrons (KE >= 5 MeV)
  TH1D* hProton        = nullptr;  // N_protons
  TH1D* hPiCharged     = nullptr;  // N_pi+ + N_pi-
  TH1D* hPiZero        = nullptr;  // N_pi0
  TH1D* hNeutronKE     = nullptr;  // neutron count vs. kinetic energy, no KE floor
  TH1D* hEavail        = nullptr;  // event count vs. E_avail, no E_avail cut
  MultVsXProfiles multVsEnu;       // avg. FS multiplicity vs. Enu, no E_avail cut
};
inline bool haveData(const GenHistos& H) { return H.hNeutron != nullptr; }

// Object names as written by make_GENIE_vs_GiBUU_histos.C / read back by
// plot_GENIE_vs_GiBUU.C -- centralized here so the two scripts can't drift
// out of sync with each other.
inline TString NameNeutron(const std::string& hc, int idx)   { return Form("nfsn_%s_%d",   hc.c_str(), idx); }
inline TString NameProton(const std::string& hc, int idx)    { return Form("nfsp_%s_%d",   hc.c_str(), idx); }
inline TString NamePiCharged(const std::string& hc, int idx) { return Form("nfspic_%s_%d", hc.c_str(), idx); }
inline TString NamePiZero(const std::string& hc, int idx)    { return Form("nfspi0_%s_%d", hc.c_str(), idx); }
inline TString NameNeutronKE(const std::string& hc, int idx) { return Form("nfsnKE_%s_%d", hc.c_str(), idx); }
inline TString NameEavail(const std::string& hc, int idx)    { return Form("hEavail_%s_%d", hc.c_str(), idx); }
inline TString NameMultVsEnu(const char* species, const std::string& hc, int idx) {
  return Form("p%sVsEnu_%s_%d", species, hc.c_str(), idx);
}
