// make_GENIE_vs_GiBUU_histos.C
//
// Reads the raw "Convenient files" for GENIE hN18 (G18_10j_00_000), GENIE
// hA18 (AR23_20i_00_000), and GiBUU R2025 (R2025_03_Defaultparams2025), and
// writes every histogram/profile plot_GENIE_vs_GiBUU.C needs to a single
// ROOT file. Split out of the original single-file
// GENIE_vs_GiBUU_simplified.C specifically so that script (styling, binning
// tweaks, which tunes go on which page, etc.) can be iterated on without
// re-running over the raw data every time -- only re-run this one when the
// selection/fill logic itself changes.
//
// Selection: flagCC, PDGnu == +14 (FHC, numu) / -14 (RHC, numubar),
// n_FS_Muons > 0, E_avail < 300 MeV. E_avail follows Colin Weber's
// addtoEavail() convention: summed true kinetic energy of protons and
// pi+/pi-, full energy of electrons/gammas/pi0/other (i.e. everything
// except muons, neutrons, and neutrinos).
//
// Per passing event, four multiplicities are histogrammed, each counting
// only final-state particles with kinetic energy (E - m) >= 5 MeV:
//   N_neutrons, N_protons, N_pi_charged (pi+ + pi-), N_pi0
// Also, every final-state neutron's kinetic energy KE = E - m (MeV, no
// kMinParticleE floor) is filled into a neutron-count-vs-KE histogram, and
// every CC/flavor/muon-selected event's E_avail (MeV; before the E_avail
// cut itself, so the cut's location is visible on the distribution) is
// filled into an E_avail histogram. Both are unit-area normalized, like
// the four multiplicity counts, so the three tunes are compared on shape
// only (their GenScaleFactor/EventWeight conventions and MC exposures
// aren't guaranteed to share an absolute normalization).
//
// Also, per CC/flavor/muon-selected event -- NOT restricted by the E_avail
// cut, since this is meant to show the full multiplicity-vs-kinematics
// trend rather than just the low-hadronic-activity subsample -- average
// final-state multiplicity (proton, neutron, pi_charged, pi0; same KE >= 5
// MeV floor) is profiled (TProfile, not unit-area normalized -- each point
// is a genuine average) vs. Enu (true incoming neutrino energy).
//
// Output: GENIE_vs_GiBUU_histos.root, containing every GenHistos member for
// every (beam, tune) pair that had files on disk, named via the
// NameNeutron()/NameProton()/etc. helpers in GENIE_vs_GiBUU_common.h (so
// plot_GENIE_vs_GiBUU.C can Get() them back by the exact same names).
// (beam, tune) pairs with no files on disk are simply absent from the
// output file -- plot_GENIE_vs_GiBUU.C treats a missing object the same
// way this script treats an empty file listing.
//
// Run: cafe -bq -s <N> make_GENIE_vs_GiBUU_histos.C
//   (or: root -l -b -q make_GENIE_vs_GiBUU_histos.C -- no cafe-specific
//   flags are needed, this reads flat Convenient-file trees directly, same
//   as GENIE_comparison.C)

#include "GENIE_vs_GiBUU_common.h"

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <iostream>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

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

// scale a histogram so its bin contents sum to 1 -- every histogram written
// out (except the multVsEnu profiles) goes through this, so the three tunes
// are compared on shape only.
void normalizeToUnitArea(TH1D* h) {
  // Integral(0, nbins+1) includes the underflow/overflow bins -- without
  // this, hNeutronKE (range truncated to 200 MeV) and hEavail (range
  // truncated to 50 MeV) would only sum the visible slice, silently
  // rescaling it to 100% regardless of how many events actually fall
  // outside the plotted range, rather than reflecting what fraction of all
  // events that slice really represents.
  const double integral = h->Integral(0, h->GetNbinsX() + 1);
  if (integral > 0) h->Scale(1.0 / integral);
}

// multiplicity histogram with bins {0, 1, 2, "3+"}
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

MultVsXProfiles makeMultVsXProfilesVar(const std::string& hc, int idx,
                                        const std::vector<double>& edges) {
  const int nb = (int)edges.size() - 1;
  MultVsXProfiles P;
  auto make = [&](const char* species) {
    return new TProfile(NameMultVsEnu(species, hc, idx),
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

  H.hNeutron   = makeCountHist(NameNeutron(hc, idx),   "N_{neutrons} (KE #geq 5 MeV)");
  H.hProton    = makeCountHist(NameProton(hc, idx),    "N_{protons} (KE #geq 5 MeV)");
  H.hPiCharged = makeCountHist(NamePiCharged(hc, idx), "N_{#pi^{#pm}} (KE #geq 5 MeV)");
  H.hPiZero    = makeCountHist(NamePiZero(hc, idx),    "N_{#pi^{0}} (KE #geq 5 MeV)");
  H.hNeutronKE = new TH1D(NameNeutronKE(hc, idx),
    ";KE_{n} [MeV];Neutrons / bin [a.u.]", kNNeutronKEBins, kNeutronKELo, kNeutronKEHi);
  H.hNeutronKE->Sumw2();
  H.hEavail = new TH1D(NameEavail(hc, idx),
    ";E_{avail} [MeV];Events / bin [a.u.]", kNEavailBins, kEavailLo, kEavailHi);
  H.hEavail->Sumw2();
  H.multVsEnu = makeMultVsXProfilesVar(hc, idx, kEnuEdges);

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

void make_GENIE_vs_GiBUU_histos() {
  const char* outroot = "GENIE_vs_GiBUU_histos.root";
  TFile* fOut = new TFile(outroot, "RECREATE");

  const std::string hcs[] = {"FHC", "RHC"};
  for (const auto& hc : hcs) {
    std::cout << "\n=== " << hc << " ===\n";
    for (size_t i = 0; i < kGens.size(); ++i) {
      std::cout << kGens[i].label << ":\n";
      GenHistos H = processGen(kGens[i], hc, (int)i);
      if (!haveData(H)) continue;  // no files on disk -- nothing written,
                                    // plot_GENIE_vs_GiBUU.C's Get() will
                                    // just come back null for this pair
      fOut->cd();
      H.hNeutron->Write();
      H.hProton->Write();
      H.hPiCharged->Write();
      H.hPiZero->Write();
      H.hNeutronKE->Write();
      H.hEavail->Write();
      H.multVsEnu.proton->Write();
      H.multVsEnu.neutron->Write();
      H.multVsEnu.piCharged->Write();
      H.multVsEnu.pi0->Write();
    }
  }

  fOut->Close();
  std::cout << "\nDone. Saved to " << outroot << "\n";
}
