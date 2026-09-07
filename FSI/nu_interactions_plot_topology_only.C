// nu_interactions_plot_topology_only.C
// Plots for nu_interactions_spectra_topology_only.root (output of
// nu_interactions_spectra_topology_only.C):
//   - visible hadronic energy (ana::kNumuHadVisMeV, MeV) spectrum split by
//     FSI grouping (see kFSI_* in cuts.C), inclusive spectrum drawn as a
//     gray background, one page per beam, plotted up to 350 MeV (FHC) /
//     200 MeV (RHC). The category list is DIFFERENT per beam -- each only
//     carries categories that are actually populated for that beam (per
//     check_final_states.C's proton x pion table and top-final-state
//     list; see CategoriesForBeam() below and the matching comment in
//     nu_interactions_spectra_topology_only.C):
//       FHC: 1N1pi, 1p+0pi, >=2p+0pi, 1p+Npi, >=2p+Npi
//       RHC: 1N, >=2N, 1N1pi, 1p+0pi, >=2p+0pi, 1p+Npi
//     NOTE: neither beam's category list exhausts every >=1-neutron event
//     -- 0-proton events with >=2 pions (both beams), or (FHC only)
//     0-proton events with 0/1 pion, still fall outside everything listed
//     for that beam, so "Total" will visibly undershoot the gray
//     inclusive background.
// One multi-page PDF, hadVisE_topology.pdf (FHC then RHC), following the
// same open/close-bracket convention as plot_mc_spectra.C's per-sample PDFs.
//
// Styled to match plot_mc_spectra.C: same POT constants and beam label,
// same NOvA watermark, same LoadH1/StyleHist helpers.
//
// Run:
//   root -l -b -q 'nu_interactions_plot_topology_only.C("nu_interactions_spectra_topology_only.root")'

#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TStyle.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

const double kFHCPOT = 14.2283e20;
const double kRHCPOT = 12.5003e20;

namespace
{
  // ── NOvA watermark ─────────────────────────────────────────────────────────
  void DrawWatermark()
  {
    TLatex t; t.SetNDC(); t.SetTextSize(0.031); t.SetTextColor(kGray+1);
    t.DrawLatex(0.2, 0.805, "NOvA ND Simulation");
    t.DrawLatex(0.2, 0.780, "Work In Progress");
  }

  // Whole label ("Prod 5.1 - RHC POT 15e20"), right-aligned as one piece.
  // Previously split into a left piece (near the y-axis) and a right
  // piece, but the left piece sat in the same top-left corner as the
  // y-axis's auto scientific-notation exponent (e.g. "x10^3") -- keeping
  // everything on the right avoids that entirely, so this can sit lower
  // than the split version needed to.
  void DrawBeamLabel(const std::string& text)
  {
    if (text.empty()) return;
    TLatex t; t.SetNDC(); t.SetTextSize(0.042); t.SetTextFont(62);
    t.SetTextAlign(31);
    t.DrawLatex(0.93, 0.92, text.c_str());
  }

  // Fixed top-right label for each beam -- the only text DrawBeamLabel
  // should ever be given (no sample name, no other variation).
  std::string BeamPOTLabel(const std::string& beam)
  {
    if (beam == "FHC") return "Prod 5.1 - FHC POT 14e20";
    if (beam == "RHC") return "Prod 5.1 - RHC POT 15e20";
    return "";
  }

  // ── Load a 1D Spectrum saved via Spectrum::SaveTo(dir, key). Its "hist" key
  // is a plain TH1D, rescaled from its stored POT to targetPOT. ─────────────
  TH1D* LoadH1(TDirectory* dir, const char* key, double targetPOT)
  {
    if (!dir) { std::cerr << "Null dir for key " << key << "\n"; return nullptr; }
    TDirectory* sub = dynamic_cast<TDirectory*>(dir->Get(key));
    if (!sub) { std::cerr << "Cannot find spectrum: " << key << "\n"; return nullptr; }
    TH1D* h = dynamic_cast<TH1D*>(sub->Get("hist"));
    if (!h) { std::cerr << "No hist in spectrum: " << key << "\n"; return nullptr; }
    h = (TH1D*)h->Clone();
    h->SetDirectory(nullptr);
    auto* pp = dynamic_cast<TH1*>(sub->Get("pot"));
    const double stored = (pp && pp->GetBinContent(1) > 0) ? pp->GetBinContent(1) : 1.0;
    if (targetPOT > 0) h->Scale(targetPOT / stored);
    return h;
  }

  void StyleHist(TH1D* h, Color_t col, Style_t lstyle, Width_t lw, bool fill)
  {
    if (!h) return;
    h->SetLineColor(col);
    h->SetLineStyle(lstyle);
    h->SetLineWidth(lw);
    h->SetFillStyle(fill ? 1001 : 0);
    if (fill) h->SetFillColorAlpha(col, 0.15);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetLabelSize(0.042);
    h->SetStats(0);
  }

  // One FSI grouping category: the TrueHadVisE_<key>_<beam> spectrum key
  // suffix, its legend label, and its line color.
  struct Category
  {
    std::string key;
    std::string label;
    Color_t     color;
  };

  // Different category list per beam -- see the file-header comment for
  // why each entry is included/excluded for that beam.
  std::vector<Category> CategoriesForBeam(const std::string& beam)
  {
    std::vector<Category> cats;
    if (beam == "RHC") {
      cats.push_back({"1N",  "1 Neutron",        kBlue});
      cats.push_back({"g2N", "#geq2 Neutrons",   kAzure+2});
    }
    cats.push_back({"1N1pi", "0p, 1 Neutron, 1 Pion",       kCyan+2});
    cats.push_back({"1p0pi", "1p, N Neutrons, 0 Pions",     kRed});
    cats.push_back({"2p0pi", "#geq2p, N Neutrons, 0 Pions", kOrange+7});
    cats.push_back({"1pNpi", "1p, N Neutrons, N Pions",     kGreen+2});
    if (beam == "FHC") {
      cats.push_back({"2pNpi", "#geq2p, N Neutrons, N Pions", kMagenta});
    }
    return cats;
  }

  // ── Visible hadronic energy split by FSI grouping, inclusive drawn as a
  // gray background. xMax > 0 clips the x-axis; 0 means full range. `hists`
  // holds one loaded histogram per CategoriesForBeam() entry, same order.
  // See kFSI_* in cuts.C and the proton x pion table in
  // check_final_states.C that motivated the split (proton count, not pion
  // count, dominates the events the old proton-free-only cuts missed). ───
  void DrawTopologyOverlay(TH1D* hIncl,
                            const std::vector<Category>& cats,
                            const std::vector<TH1D*>& hists,
                            const std::string& beamPOT,
                            const TString& canvName,
                            const TString& pdf, bool& isFirst,
                            double xMax = 0.0)
  {
    StyleHist(hIncl, kGray, 1, 1, true);
    for (size_t i = 0; i < cats.size(); ++i)
      StyleHist(hists[i], cats[i].color, 1, 2, false);

    TH1D* hTot = (TH1D*)hists[0]->Clone("hTot_"+canvName);
    hTot->SetDirectory(nullptr);
    for (size_t i = 1; i < hists.size(); ++i)
      hTot->Add(hists[i]);
    StyleHist(hTot, kBlack, 1, 2, false);

    if (xMax > 0) {
      hIncl->GetXaxis()->SetRangeUser(0, xMax);
      hTot  ->GetXaxis()->SetRangeUser(0, xMax);
    }
    hIncl->GetXaxis()->SetTitle("Visible Hadronic Energy [MeV]");
    hIncl->GetYaxis()->SetTitle("Events");
    // Force scientific notation (e.g. "3 x10^3") on the y-axis rather than
    // plain digits -- TGaxis::SetMaxDigits caps how many digits an axis
    // label can show before ROOT switches to the exponent form; 3 forces
    // it for any value >=1000, which these event counts always are.
    TGaxis::SetMaxDigits(3);
    hIncl->GetYaxis()->SetNoExponent(kFALSE);

    TCanvas* c = new TCanvas(canvName, "", 800, 700);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.14); c->SetTopMargin(0.12);

    double ymax = std::max(hIncl->GetMaximum(), hTot->GetMaximum());
    hIncl->SetMaximum(ymax * 1.30); hIncl->SetMinimum(0);
    hIncl->SetTitle("");
    hIncl->Draw("HIST");
    hTot->Draw("HIST SAME");
    for (TH1D* h : hists) h->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.42, 0.50, 0.93, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->AddEntry(hTot, "Total", "l");
    for (size_t i = 0; i < cats.size(); ++i)
      leg->AddEntry(hists[i], cats[i].label.c_str(), "l");
    leg->AddEntry(hIncl, "Inclusive", "f");
    leg->Draw();

    DrawBeamLabel(beamPOT);
    DrawWatermark();

    if (isFirst) { c->Print(pdf+"("); isFirst = false; } else c->Print(pdf);

    delete leg; delete hTot; delete c;
  }

  void ProcessBeam(TFile* f, const std::string& beam, double pot,
                    const TString& pdf, bool& isFirst)
  {
    TH1D* hIncl = LoadH1(f, ("TrueHadVisE_"+beam).c_str(), pot);

    const std::vector<Category> cats = CategoriesForBeam(beam);
    std::vector<TH1D*> hists;
    hists.reserve(cats.size());
    bool ok = (hIncl != nullptr);
    for (const auto& cat : cats) {
      TH1D* h = LoadH1(f, ("TrueHadVisE_"+cat.key+"_"+beam).c_str(), pot);
      hists.push_back(h);
      ok = ok && (h != nullptr);
    }

    if (!ok) {
      std::cerr << "Missing spectra for beam " << beam << "\n";
      delete hIncl;
      for (TH1D* h : hists) delete h;
      return;
    }

    // FHC binned to 350 MeV, RHC to 200 MeV (see nu_interactions_spectra_
    // topology_only.C's bins_HadVisE_FHC/RHC) -- RHC's FSI activity sits
    // at lower HadVisE (dominated by proton-free, pion-free CCQE-like
    // events per check_final_states.C's proton x pion table).
    const double xMax = (beam == "FHC") ? 350.0 : 200.0;
    DrawTopologyOverlay(hIncl, cats, hists, BeamPOTLabel(beam),
                        "c_topology_"+beam, pdf, isFirst, xMax);

    delete hIncl;
    for (TH1D* h : hists) delete h;
  }
}

void nu_interactions_plot_topology_only(const char* infile = "nu_interactions_spectra_topology_only.root",
                                         const char* outpdf = "hadVisE_topology.pdf")
{
  TFile* f = TFile::Open(infile);
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << infile << "\n"; return; }

  bool isFirst = true;
  ProcessBeam(f, "FHC", kFHCPOT, outpdf, isFirst);
  ProcessBeam(f, "RHC", kRHCPOT, outpdf, isFirst);

  if (!isFirst) {
    // Close the multi-page PDF (ROOT convention: an empty last c->Print(pdf+")")).
    TCanvas dummy;
    dummy.Print(TString(outpdf) + ")");
  }

  f->Close();
  std::cout << "Saved to " << outpdf << "\n";
}
