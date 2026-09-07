// plot_GENIE_vs_GiBUU.C
//
// Reads GENIE_vs_GiBUU_histos.root (output of make_GENIE_vs_GiBUU_histos.C)
// and produces the comparison PDF. Split out of the original single-file
// GENIE_vs_GiBUU_simplified.C specifically so this half -- styling, binning
// display tweaks, which tunes go on which page, legend text, etc. -- can be
// iterated on freely without re-running over the raw "Convenient files"
// every time; only re-run make_GENIE_vs_GiBUU_histos.C when the underlying
// selection/fill logic changes.
//
// Styled to match nu_interactions_plot_topology_only.C: same NOvA
// watermark, same top-right beam label convention (just "FHC"/"RHC", not a
// full selection description), same forced y-axis scientific notation.
//
// Output: a single multi-page PDF, ./G18_10j_AR23_GiBUU2025_comparison.pdf,
// one page per quantity (neutron/proton/pi_charged/pi0 count, neutron KE,
// E_avail, avg. multiplicity vs. Enu). Each page is a canvas with two
// side-by-side pads (FHC left, RHC right). The count/KE/E_avail pages
// overlay all three tunes by color; the multiplicity-vs-Enu page overlays
// only GiBUU R2025 (solid) vs GENIE hN18 (dashed), species by color. A pad
// is left blank with a "No <beam> data" note if a tune has no histograms
// in the input file for that beam (i.e. had no files on disk when
// make_GENIE_vs_GiBUU_histos.C ran).
//
// Run: root -l -b -q 'plot_GENIE_vs_GiBUU.C("GENIE_vs_GiBUU_histos.root")'

#include "GENIE_vs_GiBUU_common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <algorithm>
#include <map>

// ── loading ─────────────────────────────────────────────────────────────────
// Clone+detach a TH1D/TProfile by name; nullptr (with no error) if absent --
// same convention make_GENIE_vs_GiBUU_histos.C uses for a (beam,tune) pair
// with no files on disk, so a missing object here just means that pair
// wasn't written, not that something went wrong.
TH1D* LoadH1(TFile* f, const TString& name) {
  TH1D* h = dynamic_cast<TH1D*>(f->Get(name));
  if (!h) return nullptr;
  h = (TH1D*)h->Clone();
  h->SetDirectory(nullptr);
  return h;
}
TProfile* LoadProfile(TFile* f, const TString& name) {
  TProfile* p = dynamic_cast<TProfile*>(f->Get(name));
  if (!p) return nullptr;
  p = (TProfile*)p->Clone();
  p->SetDirectory(nullptr);
  return p;
}

GenHistos LoadGenHistos(TFile* f, const std::string& hc, int idx) {
  GenHistos H;
  H.hNeutron   = LoadH1(f, NameNeutron(hc, idx));
  H.hProton    = LoadH1(f, NameProton(hc, idx));
  H.hPiCharged = LoadH1(f, NamePiCharged(hc, idx));
  H.hPiZero    = LoadH1(f, NamePiZero(hc, idx));
  H.hNeutronKE = LoadH1(f, NameNeutronKE(hc, idx));
  H.hEavail    = LoadH1(f, NameEavail(hc, idx));
  H.multVsEnu.proton    = LoadProfile(f, NameMultVsEnu("Proton",    hc, idx));
  H.multVsEnu.neutron   = LoadProfile(f, NameMultVsEnu("Neutron",   hc, idx));
  H.multVsEnu.piCharged = LoadProfile(f, NameMultVsEnu("PiCharged", hc, idx));
  H.multVsEnu.pi0       = LoadProfile(f, NameMultVsEnu("Pi0",       hc, idx));
  return H;
}

// ── NOvA watermark ─────────────────────────────────────────────────────────
void DrawWatermark() {
  TLatex t; t.SetNDC(); t.SetTextSize(0.031); t.SetTextColor(kGray+1);
  t.DrawLatex(0.2, 0.805, "NOvA ND Simulation");
  t.DrawLatex(0.2, 0.780, "Work In Progress");
}

// Top-right beam label -- just "FHC"/"RHC", right-aligned. Same convention
// (and the same y=0.92 position) as nu_interactions_plot_topology_only.C's
// DrawBeamLabel(), simplified further since there's no POT string to split
// out here.
void DrawBeamLabel(const std::string& beam) {
  if (beam.empty()) return;
  TLatex t; t.SetNDC(); t.SetTextSize(0.042); t.SetTextFont(62);
  t.SetTextAlign(31);
  t.DrawLatex(0.93, 0.92, beam.c_str());
}

// ── drawing ─────────────────────────────────────────────────────────────────
void drawComparePad(const std::vector<TH1D*>& hists,
                     const std::vector<int>& colors,
                     const std::vector<std::string>& labels,
                     const std::string& beam,
                     const char* drawOpt = "HIST") {
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.12);

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
  // Force scientific notation on the y-axis rather than plain digits (same
  // as nu_interactions_plot_topology_only.C) -- TGaxis::SetMaxDigits caps
  // how many digits an axis label can show before switching to "N x10^k".
  TGaxis::SetMaxDigits(3);
  hists[0]->GetYaxis()->SetNoExponent(kFALSE);

  TLegend* leg = new TLegend(0.48, 0.68, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);
  for (size_t i = 0; i < hists.size(); ++i)
    leg->AddEntry(hists[i], labels[i].c_str(), "l");
  leg->Draw();

  DrawBeamLabel(beam);
  DrawWatermark();
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
                             const std::string& beam) {
  std::vector<TProfile*> profA = {A.proton, A.neutron, A.piCharged, A.pi0};
  std::vector<TProfile*> profB = {B.proton, B.neutron, B.piCharged, B.pi0};

  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.12);

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
  TGaxis::SetMaxDigits(3);
  profA[0]->GetYaxis()->SetNoExponent(kFALSE);

  TLegend* leg = new TLegend(0.55, 0.63, 0.94, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  for (size_t i = 0; i < profA.size(); ++i)
    leg->AddEntry(profA[i], speciesLabels[i].c_str(), "l");
  leg->Draw();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextSize(0.028);
  tex.SetTextAlign(31);
  tex.DrawLatex(0.80, 0.58, Form("solid: %s    dashed: %s", labelA.c_str(), labelB.c_str()));

  DrawBeamLabel(beam);
  DrawWatermark();
}

void plot_GENIE_vs_GiBUU(const char* infile = "GENIE_vs_GiBUU_histos.root") {
  TFile* f = TFile::Open(infile);
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << infile << "\n"; return; }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLabelSize(0.05, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gROOT->ForceStyle();

  const std::string hcs[] = {"FHC", "RHC"};

  // ── Phase 1: load every tune, both beams, before any drawing ────────────
  std::vector<GenHistos>     allH[2];
  std::vector<int>           colors[2];
  std::vector<std::string>   labels[2];
  std::vector<size_t>        genIdx[2];   // kGens index of each allH[hci] entry
  std::map<size_t, GenHistos*> byIdx[2];  // kGens index -> its GenHistos, for
                                           // picking out GiBUU R2025/GENIE hN18
                                           // specifically in Phase 3 below

  for (int hci = 0; hci < 2; ++hci) {
    for (size_t i = 0; i < kGens.size(); ++i) {
      GenHistos H = LoadGenHistos(f, hcs[hci], (int)i);
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
      drawComparePad(hs, colors[hci], labels[hci], hcs[hci], q.drawOpt);
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
      drawMultiplicityVsXPad(Ha->multVsEnu, Hb->multVsEnu, multColors, multLabels,
                              kGens[kIdxGiBUU].label, kGens[kIdxHN18].label, hcs[hci]);
    }
    if (isFirst) { cMult->Print(outpdf+"("); isFirst = false; } else cMult->Print(outpdf);
    delete cMult;
  }

  if (!isFirst) {
    // Close the multi-page PDF (ROOT convention: an empty last c->Print(pdf+")")).
    TCanvas dummy;
    dummy.Print(outpdf + ")");
  }

  f->Close();
  std::cout << "\nDone. Saved to " << outpdf << "\n";
}
