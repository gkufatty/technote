// Run:
// cafe -bq plot_relevant_systs_ranking.C
//
// q0lt100MeV variant of ../genie/plot_relevant_systs_ranking.C: reads every
// neutron_spectra_systs_diagnostics_*_metrics.csv IN THIS FOLDER (produced
// by this folder's make_neutron_spectra_systs_diagnostics_{high,medium,low}.C
// + plot_neutron_spectra_systs_diagnostics_ranking.C, all restricted to
// HadVisE < 100 MeV) and, for the 11 systematics whitelisted in
// Q0Lt100MeVRelevantSystsWhitelist.cxx, draws a single ranked bar chart of
// each systematic's largest shape chi2 across every diagnostic and beam it
// was evaluated in -- same technique as ../genie/'s copy, different (smaller,
// HadVisE < 100 MeV) whitelist. Companion to
// make_neutron_spectra_systs_multiverse_filtered.C's whitelist derivation
// and to relevant_systs_table_q0lt100MeV.tex's category-comparison table.
//
// This is a plain CSV-driven plot (no CAFAna Spectrum/SpectrumLoader
// involved), so it only needs ROOT, not the cafana weights/cuts includes
// the other scripts in this folder need.
//
// Produces neutron_spectra_relevant_systs_ranking.pdf and a companion CSV
// (neutron_spectra_relevant_systs_ranking.csv) with, per whitelisted
// systematic: category, max shape chi2, and which (diagnostic, beam) gave
// that max -- the same numbers behind relevant_systs_table_q0lt100MeV.tex.

#ifdef __CINT__
void plot_relevant_systs_ranking()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
#include "Q0Lt100MeVRelevantSystsWhitelist.cxx"

  // {diagnostic label (matches the *_metrics.csv filename stem, minus the
  // "neutron_spectra_systs_diagnostics_" prefix and "_metrics.csv" suffix),
  // physics category}. Same 14 diagnostics as ../genie/'s copy.
  const std::vector<std::pair<std::string, std::string>> kDiagnosticFiles = {
    {"FSI_ke",        "FSI"},
    {"FSI_prongs",    "FSI"},
    {"FSI_confusion", "FSI"},
    {"MEC_q0",        "MEC"},
    {"MEC_mult",      "MEC"},
    {"MEC_topo",      "MEC"},
    {"RES_q0",        "RES"},
    {"RES_mult",      "RES"},
    {"RES_topo",      "RES"},
    {"QE_clean",      "QE"},
    {"QE_incl",       "QE"},
    {"DIS_mult",      "DIS"},
    {"DIS_topo",      "DIS"},
    {"LowPriority_mult", "LowPriority"},
  };

  struct Best { double chi2 = -1.0; std::string diagnostic, beam; };

  std::vector<std::string> SplitCSVLine(const std::string& line)
  {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) out.push_back(field);
    return out;
  }

  // Reads one metrics CSV (beam,syst,chi2_shape_up,chi2_shape_down,...) and
  // updates `best[syst]` with the largest |chi2_shape| seen, tagging which
  // diagnostic/beam produced it.
  void ScanCSV(const std::string& path, const std::string& diagnosticLabel,
               std::map<std::string, Best>& best)
  {
    std::ifstream in(path);
    if (!in.is_open()) {
      std::cerr << "[WARN] could not open " << path << " -- skipping\n";
      return;
    }
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
      if (line.empty()) continue;
      const std::vector<std::string> f = SplitCSVLine(line);
      if (f.size() < 4) continue;
      const std::string& beam = f[0];
      const std::string& syst = f[1];
      const double chi2Up   = std::fabs(std::stod(f[2]));
      const double chi2Down = std::fabs(std::stod(f[3]));
      const double chi2 = std::max(chi2Up, chi2Down);

      Best& b = best[syst];
      if (chi2 > b.chi2) { b.chi2 = chi2; b.diagnostic = diagnosticLabel; b.beam = beam; }
    }
  }

  int ColorForCategory(const std::string& category)
  {
    if (category == "FSI") return kRed + 1;
    if (category == "MEC") return kAzure + 2;
    if (category == "RES") return kGreen + 2;
    if (category == "QE")  return kOrange + 1;
    if (category == "DIS") return kMagenta + 1;
    return kGray + 1; // LowPriority, shouldn't appear (none whitelisted)
  }

  // ── NOvA watermark, matching nu_interactions_plot_topology_only.C's
  // DrawWatermark so every technote plot carries the same identification.
  // No per-beam POT label here (unlike the other plotters in this folder)
  // since this ranking mixes FHC and RHC together by construction. ───────
  void DrawWatermark()
  {
    TLatex t; t.SetNDC(); t.SetTextSize(0.031); t.SetTextColor(kGray+1);
    t.DrawLatex(0.2, 0.755, "NOvA ND Simulation");
    t.DrawLatex(0.2, 0.730, "Work In Progress");
  }
}

void plot_relevant_systs_ranking()
{
  gStyle->SetOptStat(0);

  std::map<std::string, Best> best;
  for (size_t i = 0; i < kDiagnosticFiles.size(); ++i) {
    const std::string& label = kDiagnosticFiles[i].first;
    const std::string path = "neutron_spectra_systs_diagnostics_" + label + "_metrics.csv";
    ScanCSV(path, label, best);
  }

  // Keep only the whitelisted (relevant) systematics, sorted by descending
  // max chi2 -- this IS the ranking.
  struct Row { std::string syst, category, diagnostic, beam; double chi2; };
  std::vector<Row> rows;
  for (const auto& kv : best) {
    if (!kRelevantSysts.count(kv.first)) continue;
    const std::string category = kRelevantSystCategory.count(kv.first)
                                    ? kRelevantSystCategory.at(kv.first) : "?";
    rows.push_back({kv.first, category, kv.second.diagnostic, kv.second.beam, kv.second.chi2});
  }
  std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.chi2 > b.chi2; });

  std::cout << "Ranked " << rows.size() << " / " << kRelevantSysts.size()
            << " whitelisted systematics (largest shape chi2 first, HadVisE < 100 MeV):\n";
  for (const Row& r : rows) {
    std::cout << "  " << r.syst << " [" << r.category << "] chi2=" << r.chi2
              << " (" << r.diagnostic << ", " << r.beam << ")\n";
  }
  if (rows.size() != kRelevantSysts.size()) {
    std::cout << "[WARN] " << kRelevantSysts.size() << " systematics are whitelisted but only "
              << rows.size() << " were found in the metrics CSVs -- some CSVs may be missing "
              << "or the whitelist may be out of date.\n";
  }

  // ── bar chart, color-coded by category ───────────────────────────────
  const int n = (int)rows.size();
  TCanvas c("c", "c", std::max(1200, 32 * n), 750);
  c.SetBottomMargin(0.32);
  c.SetTopMargin(0.08);
  c.SetLeftMargin(0.09);
  c.SetLogy();
  c.SetGrid();

  TH1D h("h", "", n, 0, n);
  for (int i = 0; i < n; ++i) {
    h.SetBinContent(i + 1, rows[i].chi2);
    h.GetXaxis()->SetBinLabel(i + 1, rows[i].syst.c_str());
  }
  h.SetTitle("Relevant GENIE systematics, ranked by largest shape #chi^{2} "
             "across all diagnostics (FHC+RHC), HadVisE < 100 MeV");
  h.GetYaxis()->SetTitle("max #chi^{2}_{shape} (any diagnostic, either beam)");
  h.GetXaxis()->LabelsOption("v");
  h.GetXaxis()->SetLabelSize(0.030);
  h.SetMinimum(0.5);
  h.SetMaximum(h.GetMaximum() * 3.0);

  // One histogram per category so each set of bars gets its own fill
  // color; bins outside a category's own rows stay at zero (invisible).
  std::map<std::string, std::unique_ptr<TH1D>> byCategory;
  for (const auto& catPair : kDiagnosticFiles) {
    const std::string& category = catPair.second;
    if (byCategory.count(category)) continue;
    auto hc = std::make_unique<TH1D>(("h_" + category).c_str(), "", n, 0, n);
    hc->SetFillColor(ColorForCategory(category));
    hc->SetLineColor(ColorForCategory(category));
    byCategory[category] = std::move(hc);
  }
  for (int i = 0; i < n; ++i)
    byCategory[rows[i].category]->SetBinContent(i + 1, rows[i].chi2);

  h.Draw("axis");
  for (auto& kv : byCategory) kv.second->Draw("hist same");

  TLine thresh(0, 5.0, n, 5.0);
  thresh.SetLineColor(kBlack);
  thresh.SetLineStyle(2);
  thresh.SetLineWidth(2);
  thresh.Draw();

  TLegend leg(0.80, 0.70, 0.98, 0.92);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  for (const std::string& category : {"FSI", "MEC", "RES", "QE", "DIS"}) {
    if (byCategory.count(category))
      leg.AddEntry(byCategory[category].get(), category.c_str(), "f");
  }
  leg.AddEntry(&thresh, "#chi^{2}=5 cut", "l");
  leg.Draw();

  DrawWatermark();

  const std::string pdfName = "neutron_spectra_relevant_systs_ranking.pdf";
  c.Print(pdfName.c_str());
  std::cout << "Saved ranking plot to " << pdfName << "\n";

  // ── companion CSV (same numbers, machine-readable) ───────────────────
  const std::string csvName = "neutron_spectra_relevant_systs_ranking.csv";
  std::ofstream csv(csvName);
  csv << "rank,syst,category,max_chi2_shape,diagnostic,beam\n";
  for (int i = 0; i < n; ++i) {
    csv << (i + 1) << ',' << rows[i].syst << ',' << rows[i].category << ','
        << rows[i].chi2 << ',' << rows[i].diagnostic << ',' << rows[i].beam << '\n';
  }
  csv.close();
  std::cout << "Saved ranking table to " << csvName << "\n";
}

#endif
