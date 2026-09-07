// GENIE_vs_GiBUU_simplified.C
//
// Simplified 3-way generator comparison: GENIE G18_10j_00_000, GENIE
// AR23_20i_00_000, and GiBUU R2025_03_Defaultparams2025. Trimmed down from
// GENIE_comparison.C -- same "Convenient files" data source, selection, and
// E_avail convention, and the same core multiplicity/energy plots, but:
//   - only these three tunes (not all 16 in GENIE_comparison.C's kGens)
//   - no FSI-topology grid, no kinematic profiles (Enu/Emu/MuPt/CosThetaMu),
//     no cut_report.txt / topology report, no "comparison group" abstraction
//     (there's only one three-way comparison here, so no need for it)
//   - dropped plot types: neutron_count_allEavail, neutron_ke,
//     neutron_vs_muon_costheta (kept: the four multiplicity counts, neutron
//     energy, and E_avail)
//
// Selection (identical to GENIE_comparison.C): flagCC, PDGnu == +14 (FHC,
// numu) / -14 (RHC, numubar), n_FS_Muons > 0, E_avail < 300 MeV. E_avail
// follows Colin Weber's addtoEavail() convention: summed true kinetic
// energy of protons and pi+/pi-, full energy of electrons/gammas/pi0/other
// (i.e. everything except muons, neutrons, and neutrinos).
//
// Per passing event, four multiplicities are histogrammed, each counting
// only final-state particles with kinetic energy (E - m) >= 5 MeV:
//   N_neutrons, N_protons, N_pi_charged (pi+ + pi-), N_pi0
// Also, every final-state neutron's kinetic energy KE = E - m (MeV, no
// kMinParticleE floor) is filled into a neutron-count-vs-KE histogram, and
// every CC/flavor/muon-selected event's E_avail (MeV; before the E_avail
// cut itself, so the cut's location is visible on the distribution) is
// filled into an E_avail histogram. Every plotted histogram is normalized
// to unit area, so the three tunes are compared on shape only (their
// GenScaleFactor/EventWeight conventions and MC exposures aren't guaranteed
// to share an absolute normalization).
//
// Also, for GiBUU R2025 and GENIE hN18 ONLY (not GENIE hA18) -- per
// CC/flavor/muon-selected event, unlike everything above NOT restricted by
// the E_avail cut, since this is meant to show the full multiplicity-vs-
// kinematics trend rather than just the low-hadronic-activity subsample --
// average final-state multiplicity (proton, neutron, pi_charged, pi0; same
// KE >= 5 MeV floor) is profiled (TProfile, so each point is a genuine
// average, not unit-area normalized) vs. Enu (true incoming neutrino
// energy), using the same variable-width bins as GENIE_comparison.C's Enu
// histogram.
//
// Output: a single multi-page PDF, ./G18_10j_AR23_GiBUU2025_comparison.pdf,
// one page per quantity (neutron/proton/pi_charged/pi0 count, neutron KE,
// E_avail, avg. multiplicity vs. Enu). Each page is a canvas with
// two side-by-side pads (FHC left, RHC right). The count/KE/E_avail
// pages overlay all three tunes by color; the multiplicity-vs-Enu page
// overlays only GiBUU R2025 (solid) vs GENIE hN18 (dashed), species by
// color. A pad is left blank with a "No <beam> data" note if a tune has no
// files on disk for that beam.
//
// Run: root -l -b -q GENIE_vs_GiBUU_simplified.C

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TStyle.h"
#include "TROOT.h"
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

// ── kinematic index convention (see GENIE_comparison.C) ─────────────────────
const int kIdxE = 0;
const double kMinParticleE = 0.005;  // GeV -- KE floor for the multiplicity counts
const double kEavailCut    = 0.300;  // GeV -- E_avail selection threshold

// binning for the neutron kinetic-energy-distribution plot (MeV, KE = E - m;
// no kMinParticleE floor -- every final-state neutron is filled regardless
// of whether it'd pass the 5 MeV cut used for the multiplicity counts)
const int    kNNeutronKEBins = 50;
const double kNeutronKELo    = 0.0;
const double kNeutronKEHi    = 500.0;

// binning for the E_avail distribution plot (MeV) -- the selection cut
// itself (kEavailCut above) stays in GeV; only this display histogram's
// fill value and axis are converted to MeV.
const int    kNEavailBins = 60;
const double kEavailLo    = 0.0;
const double kEavailHi    = 600.0;

// ── particle masses (GeV) ────────────────────────────────────────────────────
const double kMassProton   = 0.938272;
const double kMassNeutron  = 0.939565;
const double kMassPiC      = 0.139570;
const double kMassPi0      = 0.134977;

// species colors for the avg-multiplicity-vs-Enu plot (color = species,
// line style = tune -- see drawMultiplicityVsXPad())
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
// for the avg-multiplicity-vs-Enu plot (unlike the count/energy/E_avail
// plots above, which overlay all three tunes).
enum { kIdxHN18 = 0, kIdxHA18 = 1, kIdxGiBUU = 2 };

// ── helpers (same conventions as GENIE_comparison.C) ────────────────────────
std::vector<std::string> listRootFiles(const std::string& dir) {
  std::vector<std::string> paths;
  if (!fs::exists(dir)) return paths;
  for (const auto& e : fs::directory_iterator(dir))
    if (e.is_regular_file() && e.path().extension() == ".root")
      paths.push_back(e.path().string());
  return paths;
}

// summed kinetic energy of a species' particle list (index 0 of each inner
// vector is energy, per the FS_* convention)
double sumKE(const TTreeReaderArray<std::vector<float>>& arr, int n, double mass) {
  double sum = 0.0;
  for (int i = 0; i < n && i < (int)arr.GetSize(); ++i) {
    if (arr[i].empty()) continue;
    sum += std::max(0.f, arr[i][0] - (float)mass);
  }
  return sum;
}

// counts particles with kinetic energy (E - m) >= kMinParticleE
int countAboveKECut(const TTreeReaderArray<std::vector<float>>& arr, int n, double mass) {
  int count = 0;
  for (int i = 0; i < n && i < (int)arr.GetSize(); ++i) {
    if (arr[i].empty()) continue;
    if (arr[i][0] - mass >= kMinParticleE) ++count;
  }
  return count;
}

// scale a histogram so its bin contents sum to 1 -- every plotted histogram
// goes through this, so the three tunes are compared on shape only.
void normalizeToUnitArea(TH1D* h) {
  // Integral(0, nbins+1) includes the underflow/overflow bins -- without
  // this, hEavail (range truncated to 600 MeV) and hNeutronKE (truncated to
  // 500 MeV) would only sum the visible slice, silently rescaling it to
  // 100% regardless of how many events actually fall outside the plotted
  // range, rather than reflecting what fraction of all events that slice
  // really represents.
  const double integral = h->Integral(0, h->GetNbinsX() + 1);
  if (integral > 0) h->Scale(1.0 / integral);
}

// multiplicity histogram with bins {0, 1, 2, "3+"}
const int kNCountBins = 4;
TH1D* makeCountHist(const char* name, const char* xtitle) {
  TH1D* h = new TH1D(name, (std::string(";") + xtitle + ";Fraction of events").c_str(),
                      kNCountBins, -0.5, kNCountBins - 0.5);
  h->GetXaxis()->SetBinLabel(1, "0");
  h->GetXaxis()->SetBinLabel(2, "1");
  h->GetXaxis()->SetBinLabel(3, "2");
  h->GetXaxis()->SetBinLabel(4, "3+");
  return h;
}
void fillCount(TH1D* h, int n, double w) { h->Fill(std::min(n, kNCountBins - 1), w); }

// Four-species TProfile set for "average FS multiplicity vs. Enu" (proton,
// neutron, pi_charged, pi0 -- "other" dropped, unlike GENIE_comparison.C's
// five-species version, per request). Unrestricted by the E_avail cut
// (unlike every other plot in this file), and not unit-area normalized --
// each point is a genuine average multiplicity via TProfile, not a shape.
struct MultVsXProfiles {
  TProfile* proton    = nullptr;
  TProfile* neutron   = nullptr;
  TProfile* piCharged = nullptr;
  TProfile* pi0       = nullptr;
};

MultVsXProfiles makeMultVsXProfilesVar(const char* tag, const std::string& hc, int idx,
                                        const std::vector<double>& edges) {
  const int nb = (int)edges.size() - 1;
  MultVsXProfiles P;
  auto make = [&](const char* species) {
    return new TProfile(Form("p%s%s_%s_%d", species, tag, hc.c_str(), idx),
      ";E_{#nu} [GeV];Avg. multiplicity per neutrino", nb, edges.data());
  };
  P.proton    = make("Proton");
  P.neutron   = make("Neutron");
  P.piCharged = make("PiCharged");
  P.pi0       = make("Pi0");
  return P;
}

void fillMultVsX(MultVsXProfiles& P, double x, int nProton, int nNeutron,
                  int nPiCharged, int nPi0, double w) {
  P.proton->Fill(x, nProton, w);
  P.neutron->Fill(x, nNeutron, w);
  P.piCharged->Fill(x, nPiCharged, w);
  P.pi0->Fill(x, nPi0, w);
}

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
bool haveData(const GenHistos& H) { return H.hNeutron != nullptr; }

GenHistos processGen(const GenConfig& cfg, const std::string& hc, int idx) {
  const bool isFHC = (hc == "FHC");
  const std::string& base = isFHC ? cfg.fhcBase : cfg.rhcBase;
  const std::string dir = base + "/" + (isFHC ? cfg.numuSubdir : cfg.numubarSubdir);
  const int pdgTarget = isFHC ? 14 : -14;

  GenHistos H;
  auto paths = listRootFiles(dir);
  if (paths.empty()) {
    std::cerr << "[WARN] No files under " << dir << "\n";
    return H;
  }

  bool newFmt, newSF;
  {
    TFile* f0 = TFile::Open(paths[0].c_str(), "READ");
    TTree* t0 = f0 ? (TTree*)f0->Get("generator_data") : nullptr;
    newFmt = t0 && t0->GetBranch("PDGnu");
    newSF  = t0 && t0->GetBranch("GenScaleFactor");
    if (f0) { f0->Close(); delete f0; }
  }
  std::cout << "  " << hc << " (" << paths.size() << " files, "
            << (newFmt ? "new" : "old") << " flavor fmt)\n";
  const char* pdgBr = newFmt ? "PDGnu"          : "NuIn.Flavor";
  const char* sfBr  = newSF  ? "GenScaleFactor" : "fScaleFactor";
  const char* enuBr = newFmt ? "Enu"            : "NuIn.E_in";
  double      sfMul = newSF  ? 1.0              : 0.1;

  H.hNeutron   = makeCountHist(Form("nfsn_%s_%d",   hc.c_str(), idx), "N_{neutrons} (KE #geq 5 MeV)");
  H.hProton    = makeCountHist(Form("nfsp_%s_%d",   hc.c_str(), idx), "N_{protons} (KE #geq 5 MeV)");
  H.hPiCharged = makeCountHist(Form("nfspic_%s_%d", hc.c_str(), idx), "N_{#pi^{#pm}} (KE #geq 5 MeV)");
  H.hPiZero    = makeCountHist(Form("nfspi0_%s_%d", hc.c_str(), idx), "N_{#pi^{0}} (KE #geq 5 MeV)");
  H.hNeutronKE = new TH1D(Form("nfsnKE_%s_%d", hc.c_str(), idx),
    ";KE_{n} [MeV];Neutrons / bin [a.u.]", kNNeutronKEBins, kNeutronKELo, kNeutronKEHi);
  H.hNeutronKE->Sumw2();
  H.hEavail = new TH1D(Form("hEavail_%s_%d", hc.c_str(), idx),
    ";E_{avail} [MeV];Events / bin [a.u.]", kNEavailBins, kEavailLo, kEavailHi);
  H.hEavail->Sumw2();
  H.multVsEnu = makeMultVsXProfilesVar("VsEnu", hc, idx, kEnuEdges);

  for (const auto& path : paths) {
    TFile* file = TFile::Open(path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "[WARN] Could not open " << path << "\n";
      continue;
    }

    TTreeReader reader("generator_data", file);
    TTreeReaderValue<bool> flagCC(reader, "flagCC");
    TTreeReaderValue<int>  pdgnu (reader, pdgBr);
    TTreeReaderValue<int>  namu  (reader, "n_FS_Muons");
    TTreeReaderValue<int>  nneut (reader, "n_FS_Neutrons");
    TTreeReaderArray<std::vector<float>> aneut(reader, "FS_Neutrons");
    TTreeReaderValue<int>  nprot (reader, "n_FS_Protons");
    TTreeReaderValue<int>  npip  (reader, "n_FS_PiPs");
    TTreeReaderValue<int>  npim  (reader, "n_FS_PiMs");
    TTreeReaderValue<int>  npi0  (reader, "n_FS_Pi0s");
    TTreeReaderValue<double> gensf(reader, sfBr);
    TTreeReaderValue<double> evwgt(reader, "EventWeight");
    TTreeReaderValue<float>  enu  (reader, enuBr);

    // Eavail inputs: everything except muons, neutrons, and neutrinos
    TTreeReaderValue<int> nElectrons(reader, "n_FS_Electrons");
    TTreeReaderValue<int> nGammas   (reader, "n_FS_Gammas");
    TTreeReaderValue<int> nOthers   (reader, "n_FS_Others");
    TTreeReaderArray<std::vector<float>> electrons(reader, "FS_Electrons");
    TTreeReaderArray<std::vector<float>> gammas   (reader, "FS_Gammas");
    TTreeReaderArray<std::vector<float>> protons  (reader, "FS_Protons");
    TTreeReaderArray<std::vector<float>> pips     (reader, "FS_PiPs");
    TTreeReaderArray<std::vector<float>> pims     (reader, "FS_PiMs");
    TTreeReaderArray<std::vector<float>> pi0s     (reader, "FS_Pi0s");
    TTreeReaderArray<std::vector<float>> others   (reader, "FS_Others");

    while (reader.Next()) {
      if (!*flagCC) continue;              // CC only
      if (*pdgnu != pdgTarget) continue;    // numu (FHC) or numubar (RHC)
      if (*namu <= 0) continue;             // >=1 FS (anti)muon

      double eavail = sumKE(electrons, *nElectrons, 0.0)
                    + sumKE(gammas,    *nGammas,    0.0)
                    + sumKE(protons,   *nprot,      kMassProton)
                    + sumKE(pips,      *npip,       kMassPiC)
                    + sumKE(pims,      *npim,       kMassPiC)
                    + sumKE(pi0s,      *npi0,       0.0)
                    + sumKE(others,    *nOthers,    0.0);
      double w = (*gensf * sfMul) * *evwgt;

      // Every CC/flavor/muon-selected event, unrestricted by the E_avail
      // cut itself, so the cut's location is visible on the distribution.
      // eavail itself stays in GeV (compared against kEavailCut below) --
      // only the fill value here is converted to MeV for display.
      H.hEavail->Fill(eavail * 1000.0, w);

      // Computed for every CC/flavor/muon-selected event (not gated by the
      // E_avail cut below), since multVsEnu -- unlike everything else in
      // this function -- is meant to show the full multiplicity-vs-Enu
      // trend, not just the low-hadronic-activity subsample.
      int nGoodNeutrons = countAboveKECut(aneut, *nneut, kMassNeutron);
      int nGoodProtons  = countAboveKECut(protons, *nprot, kMassProton);
      int nGoodPiC      = countAboveKECut(pips, *npip, kMassPiC)
                        + countAboveKECut(pims, *npim, kMassPiC);
      int nGoodPi0      = countAboveKECut(pi0s, *npi0, kMassPi0);
      fillMultVsX(H.multVsEnu, *enu, nGoodProtons, nGoodNeutrons, nGoodPiC, nGoodPi0, w);

      if (eavail >= kEavailCut) continue;   // this pass's E_avail selection

      for (int i = 0; i < *nneut && i < (int)aneut.GetSize(); ++i) {
        const auto& p = aneut[i];
        if ((int)p.size() <= kIdxE) continue;
        // no kMinParticleE floor here -- every final-state neutron's KE is
        // filled, regardless of the 5 MeV cut used above. p[kIdxE] is total
        // energy E (GeV); KE = E - m, converted to MeV for display.
        H.hNeutronKE->Fill((p[kIdxE] - kMassNeutron) * 1000.0, w);
      }
      fillCount(H.hNeutron,   nGoodNeutrons, w);
      fillCount(H.hProton,    nGoodProtons,  w);
      fillCount(H.hPiCharged, nGoodPiC,      w);
      fillCount(H.hPiZero,    nGoodPi0,      w);
    }
    file->Close();
    delete file;
  }

  normalizeToUnitArea(H.hNeutron);
  normalizeToUnitArea(H.hProton);
  normalizeToUnitArea(H.hPiCharged);
  normalizeToUnitArea(H.hPiZero);
  normalizeToUnitArea(H.hNeutronKE);
  normalizeToUnitArea(H.hEavail);
  return H;
}

// ── drawing ─────────────────────────────────────────────────────────────────
void drawComparePad(const std::vector<TH1D*>& hists,
                     const std::vector<int>& colors,
                     const std::vector<std::string>& labels,
                     const char* header,
                     const char* drawOpt = "HIST") {
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);

  double ymax = 0;
  for (auto* h : hists) ymax = std::max(ymax, h->GetMaximum());
  std::string opt0 = drawOpt, optN = std::string(drawOpt) + " SAME";
  for (size_t i = 0; i < hists.size(); ++i) {
    hists[i]->SetLineColor(colors[i]);
    hists[i]->SetLineWidth(2);
    hists[i]->SetMarkerColor(colors[i]);
    hists[i]->SetMarkerStyle(20);
    hists[i]->SetMarkerSize(0.8);
    hists[i]->GetYaxis()->SetRangeUser(0, ymax * 1.35);
    hists[i]->Draw(i == 0 ? opt0.c_str() : optN.c_str());
  }

  TLegend* leg = new TLegend(0.48, 0.68, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);
  for (size_t i = 0; i < hists.size(); ++i)
    leg->AddEntry(hists[i], labels[i].c_str(), "l");
  leg->Draw();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextSize(0.034);
  tex.DrawLatex(0.13, 0.93, header);
}

void drawNoDataPad(const char* beamLabel) {
  TLatex tex;
  tex.SetNDC();
  tex.SetTextSize(0.05);
  tex.DrawLatex(0.3, 0.5, Form("No %s data", beamLabel));
}

// ── avg. multiplicity vs. Enu: GiBUU R2025 (solid) vs GENIE hN18 (dashed)
// only -- color = species, line style = tune. Same layout as
// GENIE_comparison.C's drawMultiplicityVsXPad(), just four species instead
// of five (no "other").
void drawMultiplicityVsXPad(const MultVsXProfiles& A, const MultVsXProfiles& B,
                             const std::vector<int>& colors,
                             const std::vector<std::string>& speciesLabels,
                             const std::string& labelA, const std::string& labelB,
                             const char* header) {
  std::vector<TProfile*> profA = {A.proton, A.neutron, A.piCharged, A.pi0};
  std::vector<TProfile*> profB = {B.proton, B.neutron, B.piCharged, B.pi0};

  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);

  double ymax = 0;
  for (auto* p : profA) ymax = std::max(ymax, p->GetMaximum());
  for (auto* p : profB) ymax = std::max(ymax, p->GetMaximum());

  for (size_t i = 0; i < profA.size(); ++i) {
    profA[i]->SetLineColor(colors[i]);
    profA[i]->SetLineWidth(2);
    profA[i]->SetLineStyle(1);  // solid
    profA[i]->GetYaxis()->SetRangeUser(0, ymax * 1.2);
    profA[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
  }
  for (size_t i = 0; i < profB.size(); ++i) {
    profB[i]->SetLineColor(colors[i]);
    profB[i]->SetLineWidth(2);
    profB[i]->SetLineStyle(2);  // dashed
    profB[i]->Draw("HIST SAME");
  }

  TLegend* leg = new TLegend(0.16, 0.63, 0.55, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  for (size_t i = 0; i < profA.size(); ++i)
    leg->AddEntry(profA[i], speciesLabels[i].c_str(), "l");
  leg->Draw();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextSize(0.034);
  tex.DrawLatex(0.13, 0.93, header);
  tex.SetTextSize(0.028);
  tex.DrawLatex(0.16, 0.58, Form("solid: %s    dashed: %s", labelA.c_str(), labelB.c_str()));
}

void GENIE_vs_GiBUU_simplified() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLabelSize(0.05, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gROOT->ForceStyle();

  const std::string hcs[] = {"FHC", "RHC"};
  const std::string hcTex[] = {
    "#nu_{#mu} CC, N_{FS #mu} > 0, E_{avail} < 300 MeV, FHC",
    "#bar{#nu}_{#mu} CC, N_{FS #bar{#mu}} > 0, E_{avail} < 300 MeV, RHC"
  };

  // ── Phase 1: process every tune, both beams, before any drawing ─────────
  std::vector<GenHistos>     allH[2];
  std::vector<int>           colors[2];
  std::vector<std::string>   labels[2];
  std::vector<size_t>        genIdx[2];   // kGens index of each allH[hci] entry
  std::map<size_t, GenHistos*> byIdx[2];  // kGens index -> its GenHistos, for
                                           // picking out GiBUU R2025/GENIE hN18
                                           // specifically in Phase 3 below

  for (int hci = 0; hci < 2; ++hci) {
    std::cout << "\n=== " << hcs[hci] << " ===\n";
    for (size_t i = 0; i < kGens.size(); ++i) {
      std::cout << kGens[i].label << ":\n";
      GenHistos H = processGen(kGens[i], hcs[hci], (int)i);
      if (!haveData(H)) continue;
      allH[hci].push_back(H);
      colors[hci].push_back(kGens[i].color);
      labels[hci].push_back(kGens[i].label);
      genIdx[hci].push_back(i);
    }
    if (allH[hci].empty())
      std::cerr << "[WARN] No data for " << hcs[hci]
                << " -- its pad will be left blank in every canvas below\n";
    // Built only now, after every push_back for this beam is done, so the
    // pointers below can't be invalidated by a later vector reallocation.
    for (size_t k = 0; k < allH[hci].size(); ++k)
      byIdx[hci][genIdx[hci][k]] = &allH[hci][k];
  }

  // ── Phase 2: one two-pad canvas per quantity ─────────────────────────────
  struct Quantity { TH1D* GenHistos::*member; const char* tag; const char* drawOpt; };
  const std::vector<Quantity> quantities = {
    {&GenHistos::hNeutron,       "neutron_count",     "HIST"},
    {&GenHistos::hProton,        "proton_count",      "HIST"},
    {&GenHistos::hPiCharged,     "pi_charged_count",  "HIST"},
    {&GenHistos::hPiZero,        "pi0_count",         "HIST"},
    {&GenHistos::hNeutronKE,     "neutron_ke",        "E"},
    {&GenHistos::hEavail,        "eavail",            "HIST"},
  };

  const TString outpdf = kOutDir + "_comparison.pdf";
  bool isFirst = true;
  for (const auto& q : quantities) {
    auto* c = new TCanvas(Form("c_%s", q.tag), Form("c_%s", q.tag), 1600, 600);
    c->Divide(2, 1);
    for (int hci = 0; hci < 2; ++hci) {
      c->cd(hci + 1);
      if (allH[hci].empty()) { drawNoDataPad(hcs[hci].c_str()); continue; }
      std::vector<TH1D*> hs;
      for (auto& H : allH[hci]) hs.push_back(H.*q.member);
      drawComparePad(hs, colors[hci], labels[hci], hcTex[hci].c_str(), q.drawOpt);
    }
    if (isFirst) { c->Print(outpdf+"("); isFirst = false; } else c->Print(outpdf);
    delete c;
  }

  // ── Phase 3: avg. multiplicity vs. Enu -- GiBUU R2025 vs GENIE hN18 only ──
  {
    const std::vector<int> multColors =
      {kColorProton, kColorNeutron, kColorPiCharged, kColorPi0};
    const std::vector<std::string> multLabels =
      {"proton", "neutron", "#pi^{#pm}", "#pi^{0}"};

    auto* cMult = new TCanvas("c_multiplicity_vs_Enu", "c_multiplicity_vs_Enu", 1600, 600);
    cMult->Divide(2, 1);
    for (int hci = 0; hci < 2; ++hci) {
      cMult->cd(hci + 1);
      if (!byIdx[hci].count(kIdxGiBUU) || !byIdx[hci].count(kIdxHN18)) {
        drawNoDataPad(hcs[hci].c_str());
        continue;
      }
      GenHistos* Ha = byIdx[hci].at(kIdxGiBUU);
      GenHistos* Hb = byIdx[hci].at(kIdxHN18);
      const std::string hdr = hcTex[hci] + " -- avg. FS multiplicity vs. E_{#nu} (no E_{avail} cut)";
      drawMultiplicityVsXPad(Ha->multVsEnu, Hb->multVsEnu, multColors, multLabels,
                              kGens[kIdxGiBUU].label, kGens[kIdxHN18].label, hdr.c_str());
    }
    if (isFirst) { cMult->Print(outpdf+"("); isFirst = false; } else cMult->Print(outpdf);
    delete cMult;
  }

  if (!isFirst) {
    // Close the multi-page PDF (ROOT convention: an empty last c->Print(pdf+")")).
    TCanvas dummy;
    dummy.Print(outpdf + ")");
  }

  std::cout << "\nDone. Saved to " << outpdf << "\n";
}
