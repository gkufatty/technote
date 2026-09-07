// plot_mc_spectra.C
// Plots for make_mc_spectra.root (output of make_mc_spectra.C):
//   - GENIE GENIE neutron KE vs. visible ("left a prong") neutron KE,
//     0-250 MeV, with a visible/primary ratio panel, plus a shape-normalized
//     (each scaled to its own 100%) version of the same comparison
//   - printed to stdout: percent of GENIE neutrons that are visible
//     (integral ratio, per beam/sample)
//   - true neutron KE vs. matched prong's true energy (2D correlation)
//   - visible neutron KE split by how many prongs that neutron produced
//     (1, 2, >=3 -- prongs grouped by matching parent truth 4-momentum)
//   - prong-multiplicity count histogram (how many neutrons produce
//     1 / 2 / 3+ prongs)
//   - GENIE vs. visible neutron count per event, with a ratio panel
//   - confusion matrix: N GENIE neutrons vs. N visible neutrons
//     (% of all events per cell)
// One multi-page PDF per sample (both beams together), e.g.
// make_mc_spectra_plots_sample1_q0Lo.pdf, make_mc_spectra_plots_sample3_full.pdf
// -- kept separate so pages from different samples are never mixed together.
// Plus one more PDF, make_mc_spectra_plots_prongs_per_neutron_overlay.pdf,
// a single page split into two side-by-side panels (FHC left, RHC right),
// each overlaying "prongs per visible neutron" for the full sample against
// the low-hadronic-energy sample.
//
// Run:
//   root -l -b -q 'plot_mc_spectra.C("make_mc_spectra.root")'

#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TParameter.h"
#include "TColor.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <utility>
#include <vector>

const double kFHCPOT = 14.2283e20;
const double kRHCPOT = 12.5003e20;


namespace
{
  // ── NOvA watermark ─────────────────────────────────────────────────────────
  void DrawWatermark()
  {
    TLatex t; t.SetNDC(); t.SetTextSize(0.031); t.SetTextColor(kGray+1);
    t.DrawLatex(0.2, 0.855, "NOvA ND Simulation");
    t.DrawLatex(0.2, 0.830, "Work In Progress");
  }

  // Splits "Prod 5.1 - FHC POT 14e20" into a left piece ("Prod 5.1 - FHC")
  // and a right piece ("POT 14e20"), drawn at the same height, with more
  // clearance above the plot frame than a single right-aligned line.
  void DrawBeamLabel(const std::string& text)
  {
    if (text.empty()) return;
    const double y = 0.92;
    const size_t potPos = text.find(" POT ");
    if (potPos == std::string::npos) {
      TLatex t; t.SetNDC(); t.SetTextSize(0.045); t.SetTextFont(62);
      t.SetTextAlign(31);
      t.DrawLatex(0.93, y, text.c_str());
      return;
    }
    const std::string left  = text.substr(0, potPos);
    const std::string right = text.substr(potPos + 1);

    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.042); tl.SetTextFont(62);
    tl.SetTextAlign(11);
    tl.DrawLatex(0.15, y, left.c_str());

    TLatex tr; tr.SetNDC(); tr.SetTextSize(0.042); tr.SetTextFont(62);
    tr.SetTextAlign(31);
    tr.DrawLatex(0.93, y, right.c_str());
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
  // is a plain TH1D (CAFAna flattens even multi-axis spectra into a product-
  // binned 1D histogram there, which is why 2D spectra need LoadH2 instead). ─
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

  // ── Load a raw TH2 written directly under "h2_<key>" (see make_mc_spectra.C:
  // ProcessBeam, which exports s2_KE_vs_prongE this way, since Spectrum::SaveTo
  // can't preserve a genuine 2D histogram). POT stored alongside as
  // "h2_<key>_pot". ──────────────────────────────────────────────────────────
  TH2* LoadH2(TDirectory* dir, const char* key, double targetPOT)
  {
    if (!dir) { std::cerr << "Null dir for key " << key << "\n"; return nullptr; }
    const std::string h2name = std::string("h2_") + key;
    TH2* h = dynamic_cast<TH2*>(dir->Get(h2name.c_str()));
    if (!h) { std::cerr << "No TH2 object: " << h2name << "\n"; return nullptr; }
    h = (TH2*)h->Clone();
    h->SetDirectory(nullptr);
    auto* pp = dynamic_cast<TParameter<double>*>(dir->Get((h2name+"_pot").c_str()));
    const double stored = (pp && pp->GetVal() > 0) ? pp->GetVal() : 1.0;
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

  // ── Primary-vs-visible overlay (KE or total E) with a visible/primary
  // ratio panel ───────────────────────────────────────────────────────────
  void DrawKEOverlay(TH1D* hPrim, TH1D* hVis,
                      const char* title, const std::string& beamSample,
                      double xMax, const TString& canvName,
                      const TString& pdf, bool& isFirst,
                      const char* xTitle = "Neutron KE [MeV]",
                      const char* yTitle = "Neutrons")
  {
    StyleHist(hPrim, TColor::GetColor(105, 70,166), 2, 3, false);
    StyleHist(hVis,  TColor::GetColor(100,149,237), 1, 2, true );
    if (hPrim) { hPrim->GetXaxis()->SetTitle(xTitle); hPrim->GetYaxis()->SetTitle(yTitle); }
    if (hVis)  { hVis ->GetXaxis()->SetTitle(xTitle); hVis ->GetYaxis()->SetTitle(yTitle); }
    if (hPrim && xMax > 0) hPrim->GetXaxis()->SetRangeUser(0, xMax);
    if (hVis  && xMax > 0) hVis ->GetXaxis()->SetRangeUser(0, xMax);

    TCanvas* c = new TCanvas(canvName, "", 800, 800);
    TPad* pTop = new TPad("pT_"+canvName, "", 0, 0.30, 1, 1);
    pTop->SetBottomMargin(0.02); pTop->SetLeftMargin(0.15); pTop->SetTopMargin(0.10);
    pTop->Draw(); pTop->cd();

    double ymax = 0;
    if (hPrim) ymax = std::max(ymax, hPrim->GetMaximum());
    if (hVis)  ymax = std::max(ymax, hVis ->GetMaximum());
    if (hPrim) { hPrim->SetMaximum(ymax*1.45); hPrim->SetMinimum(0);
                 hPrim->GetXaxis()->SetLabelSize(0); hPrim->Draw("HIST"); }
    if (hVis)  hVis->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.52, 0.68, 0.93, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.044);
    if (hPrim) leg->AddEntry(hPrim, "All GENIE neutrons",       "l");
    if (hVis)  leg->AddEntry(hVis,  "Visible Neutrons",     "lf");
    leg->Draw();
    DrawBeamLabel(beamSample);
    pTop->Update();

    c->cd();
    TPad* pBot = new TPad("pB_"+canvName, "", 0, 0, 1, 0.30);
    pBot->SetTopMargin(0.02); pBot->SetBottomMargin(0.30);
    pBot->SetLeftMargin(0.15); pBot->Draw(); pBot->cd();
    TH1D* hRatio = nullptr;
    if (hVis && hPrim) {
      hRatio = (TH1D*)hVis->Clone("hRatio_"+canvName);
      hRatio->SetDirectory(nullptr);
      hRatio->Divide(hPrim);
      hRatio->SetLineColor(kBlack); hRatio->SetLineWidth(2); hRatio->SetFillStyle(0);
      hRatio->GetYaxis()->SetTitle("Visible / Primary");
      hRatio->GetYaxis()->SetTitleSize(0.11); hRatio->GetYaxis()->SetTitleOffset(0.55);
      hRatio->GetYaxis()->SetLabelSize(0.10); hRatio->GetYaxis()->SetNdivisions(504);
      hRatio->GetXaxis()->SetTitle(xTitle);
      hRatio->GetXaxis()->SetTitleSize(0.13); hRatio->GetXaxis()->SetLabelSize(0.10);
      if (xMax > 0) hRatio->GetXaxis()->SetRangeUser(0, xMax);
      hRatio->SetMinimum(0); hRatio->SetMaximum(1.1);
      hRatio->Draw("HIST");
      TLine* unity = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                               hRatio->GetXaxis()->GetXmax(), 1.0);
      unity->SetLineColor(kGray+1); unity->SetLineStyle(2); unity->SetLineWidth(2);
      unity->Draw();
    }
    pBot->Update();
    c->cd(); DrawWatermark();
    if (isFirst) { c->Print(pdf+"("); isFirst = false; } else c->Print(pdf);
    delete hRatio; delete c;
  }

  // ── Visible-neutron KE split by how many prongs that neutron produced ──────
  void DrawKEByProngMultiplicity(TH1D* h1, TH1D* h2, TH1D* h3plus,
                                  const std::string& beamSample,
                                  const TString& canvName,
                                  const TString& pdf, bool& isFirst)
  {
    StyleHist(h1,     TColor::GetColor(100,149,237), 1, 3, false);
    StyleHist(h2,     TColor::GetColor(119,198,110), 1, 3, false);
    StyleHist(h3plus, TColor::GetColor(229,115,115), 1, 3, false);

    // Each curve rescaled to unit area (counts / total), so the three
    // multiplicity classes compare KE *shape* rather than absolute rate.
    auto normalizeShape = [](TH1D* h) {
      if (!h) return;
      const double total = h->Integral(0, h->GetNbinsX()+1);
      if (total > 0) h->Scale(1.0 / total);
    };
    normalizeShape(h1); normalizeShape(h2); normalizeShape(h3plus);

    TCanvas* c = new TCanvas(canvName, "", 800, 650);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.14);

    double ymax = 0;
    if (h1)     ymax = std::max(ymax, h1->GetMaximum());
    if (h2)     ymax = std::max(ymax, h2->GetMaximum());
    if (h3plus) ymax = std::max(ymax, h3plus->GetMaximum());

    TH1D* hFirst = h1 ? h1 : (h2 ? h2 : h3plus);
    if (hFirst) {
      hFirst->GetXaxis()->SetTitle("Visible neutron KE [MeV]");
      hFirst->GetYaxis()->SetTitle("Counts (normalized)");
      hFirst->SetMaximum(ymax * 1.45); hFirst->SetMinimum(0);
      hFirst->Draw("HIST");
    }
    if (h1     && h1     != hFirst) h1    ->Draw("HIST SAME");
    if (h2     && h2     != hFirst) h2    ->Draw("HIST SAME");
    if (h3plus && h3plus != hFirst) h3plus->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.55, 0.65, 0.93, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
    if (h1)     leg->AddEntry(h1,     "1 prong",     "l");
    if (h2)     leg->AddEntry(h2,     "2 prongs",    "l");
    if (h3plus) leg->AddEntry(h3plus, "#geq3 prongs","l");
    leg->Draw();

    DrawBeamLabel(beamSample);
    DrawWatermark();

    if (isFirst) { c->Print(pdf+"("); isFirst = false; } else c->Print(pdf);
    delete c;
  }

  // ── How many prongs each visible neutron produces (count histogram) ────────
  void DrawProngMultiplicity(TH1D* h, const std::string& beamSample,
                              const TString& canvName,
                              const TString& pdf, bool& isFirst)
  {
    if (!h) return;
    StyleHist(h, TColor::GetColor(100,149,237), 1, 2, true);
    const double total = h->Integral(0, h->GetNbinsX()+1);
    if (total > 0) h->Scale(100.0 / total);
    h->GetXaxis()->SetTitle("Prongs per Visible Neutron");
    h->GetYaxis()->SetTitle("Percentage of neutrons");
    h->GetXaxis()->SetBinLabel(1, "1"); h->GetXaxis()->SetBinLabel(2, "2");
    h->GetXaxis()->SetBinLabel(3, "3+");

    TCanvas* c = new TCanvas(canvName, "", 800, 650);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.14);
    h->SetMaximum(h->GetMaximum() * 1.30); h->SetMinimum(0);
    h->Draw("HIST");

    DrawBeamLabel(beamSample);
    DrawWatermark();

    if (isFirst) { c->Print(pdf+"("); isFirst = false; } else c->Print(pdf);
    delete c;
  }

  // ── Prongs-per-neutron overlay: full sample vs. low-hadronic-energy region ─
  // Drawn into a caller-supplied pad, so two beams can share one canvas.
  void DrawProngMultiplicitySampleOverlayPad(TPad* pad, TH1D* hFull, TH1D* hLowE,
                                              const std::string& beamPOT)
  {
    pad->cd();
    pad->SetLeftMargin(0.16); pad->SetBottomMargin(0.14); pad->SetRightMargin(0.04);

    // Pastel line palette (ColorBrewer "Pastel1"-style), solid, no fill.
    const Int_t colFull = TColor::GetColor(179, 205, 227); // pastel blue
    const Int_t colLowE = TColor::GetColor(251, 180, 174); // pastel pink

    auto normalizeShape = [](TH1D* h) {
      if (!h) return;
      const double total = h->Integral(0, h->GetNbinsX()+1);
      if (total > 0) h->Scale(1.0 / total);
    };
    normalizeShape(hFull); normalizeShape(hLowE);

    auto stylePastel = [](TH1D* h, Int_t col) {
      if (!h) return;
      h->SetLineColor(col);
      h->SetLineStyle(1); // solid
      h->SetLineWidth(3);
      h->SetFillStyle(0);
      h->SetStats(0);
      h->GetXaxis()->SetTitleSize(0.050); h->GetYaxis()->SetTitleSize(0.050);
      h->GetXaxis()->SetLabelSize(0.042); h->GetYaxis()->SetLabelSize(0.042);
    };
    stylePastel(hFull, colFull);
    stylePastel(hLowE, colLowE);

    auto labelProngBins = [](TH1D* h) {
      if (!h) return;
      h->GetXaxis()->SetTitle("Prongs per Visible Neutron");
      h->GetYaxis()->SetTitle("Counts (normalized)");
      h->GetXaxis()->SetBinLabel(1, "1"); h->GetXaxis()->SetBinLabel(2, "2");
      h->GetXaxis()->SetBinLabel(3, "3+");
    };
    labelProngBins(hFull);
    labelProngBins(hLowE);

    double ymax = 0;
    if (hFull) ymax = std::max(ymax, hFull->GetMaximum());
    if (hLowE) ymax = std::max(ymax, hLowE->GetMaximum());
    TH1D* hFirst = hFull ? hFull : hLowE;
    if (hFirst) { hFirst->SetMaximum(ymax * 1.35); hFirst->SetMinimum(0); hFirst->Draw("HIST"); }
    if (hLowE && hLowE != hFirst) hLowE->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.42, 0.68, 0.95, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
    if (hFull) leg->AddEntry(hFull, "Full sample",   "l");
    if (hLowE) leg->AddEntry(hLowE, "VisE<150 MeV",  "l");
    leg->Draw();

    DrawBeamLabel(beamPOT);
    DrawWatermark();
    pad->Update();
  }

  // One PDF, one page total: a single canvas split into two side-by-side
  // pads, FHC on the left and RHC on the right, each overlaying "prongs per
  // visible neutron" for the full sample against the low-hadronic-energy
  // sample.
  void PlotProngMultiplicityOverlay(TFile* f, const TString& outpdf)
  {
    TCanvas* c = new TCanvas("c_prongs_overlay_fhc_rhc", "", 1400, 650);
    c->Divide(2, 1, 0.006, 0.006);

    // Held until after c->Print() below: TH1::Draw() only registers the
    // pointer as a pad primitive, it doesn't copy the data, so deleting
    // these before Print() (as this loop used to, right after each beam's
    // draw call) left both pads pointing at freed memory by the time the
    // canvas was actually painted/printed -- which is why this PDF came
    // out blank.
    TH1D* toDelete[4] = {nullptr, nullptr, nullptr, nullptr};

    const char* beams[2] = {"FHC", "RHC"};
    for (int i = 0; i < 2; ++i) {
      const char* beam = beams[i];
      const std::string pathFull = std::string(beam) + "/sample3_full";
      const std::string pathLowE = std::string(beam) + "/sample1_q0Lo";
      TDirectory* dirFull = dynamic_cast<TDirectory*>(f->Get(pathFull.c_str()));
      TDirectory* dirLowE = dynamic_cast<TDirectory*>(f->Get(pathLowE.c_str()));
      if (!dirFull || !dirLowE) {
        std::cerr << "Missing directory for " << beam << " prongs-per-neutron overlay\n";
        continue;
      }
      const double pot = (i == 0) ? kFHCPOT : kRHCPOT;
      TH1D* hFull = LoadH1(dirFull, "prongs_per_visible_neutron", pot);
      TH1D* hLowE = LoadH1(dirLowE, "prongs_per_visible_neutron", pot);
      toDelete[2*i] = hFull; toDelete[2*i+1] = hLowE;

      TPad* pad = (TPad*)c->cd(i + 1);
      DrawProngMultiplicitySampleOverlayPad(pad, hFull, hLowE, BeamPOTLabel(beam));
    }

    c->Print(outpdf); // single page -- no open/close bracket dance needed
    for (TH1D* h : toDelete) delete h;
    delete c;
    std::cout << "Saved to " << outpdf << "\n";
  }

  // ── 2D: true neutron KE vs. matched prong's proton kinetic energy ──────────
  void DrawKEvsProngE(TH2* h2, const std::string& beamSample,
                       const TString& canvName,
                       const TString& pdf, bool& isFirst)
  {
    if (!h2) return;
    TCanvas* c = new TCanvas(canvName, "", 800, 700);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.14); c->SetRightMargin(0.17);
    h2->SetStats(0);
    h2->SetTitle("");
    h2->GetXaxis()->SetTitleSize(0.045); h2->GetXaxis()->SetLabelSize(0.038);
    h2->GetYaxis()->SetTitleSize(0.045); h2->GetYaxis()->SetLabelSize(0.038);
    h2->GetZaxis()->SetLabelSize(0.032);
    h2->Draw("COLZ");

    DrawBeamLabel(beamSample);
    DrawWatermark();

    if (isFirst) { c->Print(pdf+"("); isFirst = false; } else c->Print(pdf);
    delete c;
  }

  // ── Confusion matrix: N GENIE neutrons vs. N visible neutrons ────────────
  // Globally normalized: every cell divided by the total event count, so a
  // cell reads as "% of ALL events with exactly this many GENIE neutrons
  // AND exactly this many visible neutrons" (the whole matrix sums to 100%).
  void DrawConfusionMatrix(TH2* h2, const std::string& beamSample,
                            const TString& canvName,
                            const TString& pdf, bool& isFirst)
  {
    if (!h2) return;

    const int nx = h2->GetNbinsX();
    const int ny = h2->GetNbinsY();
    double total = 0;
    for (int ix = 1; ix <= nx; ++ix)
      for (int iy = 1; iy <= ny; ++iy)
        total += h2->GetBinContent(ix, iy);
    if (total > 0) h2->Scale(100.0 / total);

    double zmax = 0;
    for (int ix = 1; ix <= nx; ++ix)
      for (int iy = 1; iy <= ny; ++iy)
        zmax = std::max(zmax, h2->GetBinContent(ix, iy));

    TCanvas* c = new TCanvas(canvName, "", 800, 700);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.14); c->SetRightMargin(0.17);
    h2->SetStats(0);
    h2->SetTitle("");
    h2->GetXaxis()->SetTitle("N GENIE neutrons");
    h2->GetYaxis()->SetTitle("N visible neutrons");
    h2->GetXaxis()->SetRangeUser(0.5, 5.5);  // N GENIE >= 1 by construction now
    h2->GetYaxis()->SetRangeUser(-0.5, 5.5);
    h2->GetXaxis()->SetTitleSize(0.045); h2->GetXaxis()->SetLabelSize(0.038);
    h2->GetYaxis()->SetTitleSize(0.045); h2->GetYaxis()->SetLabelSize(0.038);
    h2->GetZaxis()->SetLabelSize(0.032);
    h2->SetMinimum(0); h2->SetMaximum(zmax > 0 ? zmax * 1.1 : 100);
    h2->SetMarkerSize(1.6);
    gStyle->SetPaintTextFormat("4.1f"); // one decimal place in each cell
    h2->Draw("COLZ TEXT45");

    DrawBeamLabel(beamSample);
    DrawWatermark();

    if (isFirst) { c->Print(pdf+"("); isFirst = false; } else c->Print(pdf);
    delete c;
  }

  void ProcessSample(TFile* f, const std::string& beam, const std::string& sample,
                      const TString& pdf, bool& isFirst,
                      std::map<std::pair<std::string,std::string>, double>& visibilityTable)
  {
    const std::string path = beam + "/" + sample;
    TDirectory* dir = dynamic_cast<TDirectory*>(f->Get(path.c_str()));
    if (!dir) { std::cerr << "Missing directory: " << path << "\n"; return; }

    const double pot = (beam == "FHC") ? kFHCPOT : kRHCPOT;
    const std::string label   = beam + " -- " + sample; // stdout logging only
    const std::string beamPOT = BeamPOTLabel(beam);      // fixed top-right plot label

    TH1D* hPrimAll  = LoadH1(dir, "ke_primary_neutrons", pot);
    TH1D* hVisAll   = LoadH1(dir, "ke_visible_neutrons",  pot);
    TH2*  h2KEvsE   = LoadH2(dir, "ke_vs_prongE_visible",  pot);

    DrawKEOverlay(hPrimAll, hVisAll, "Prod 5.1 POT", beamPOT,
                  0.0, ("c_ke_"+beam+"_"+sample).c_str(), pdf, isFirst);

    // Percent of GENIE neutrons that are visible (integral ratio, not
    // POT-scale-dependent -- both histograms carry the same scale factor).
    const double totalPrim = hPrimAll ? hPrimAll->Integral(0, hPrimAll->GetNbinsX()+1) : 0.0;
    const double totalVis  = hVisAll  ? hVisAll ->Integral(0, hVisAll ->GetNbinsX()+1) : 0.0;
    const double pctVisible = (totalPrim > 0) ? 100.0 * totalVis / totalPrim : 0.0;
    std::cout << label << ": visible/primary = " << pctVisible << "%  ("
               << totalVis << " / " << totalPrim << ")\n";
    visibilityTable[{beam, sample}] = pctVisible;

    // Shape comparison: each histogram scaled to its own 100%, so the
    // overlay compares KE *shape* rather than absolute rate.
    TH1D* hPrimShape = nullptr;
    TH1D* hVisShape  = nullptr;
    if (hPrimAll && totalPrim > 0) {
      hPrimShape = (TH1D*)hPrimAll->Clone(("hPrimShape_"+beam+"_"+sample).c_str());
      hPrimShape->SetDirectory(nullptr);
      hPrimShape->Scale(100.0 / totalPrim);
    }
    if (hVisAll && totalVis > 0) {
      hVisShape = (TH1D*)hVisAll->Clone(("hVisShape_"+beam+"_"+sample).c_str());
      hVisShape->SetDirectory(nullptr);
      hVisShape->Scale(100.0 / totalVis);
    }
    const TString shapeTitle = TString::Format(
      "Primary vs. visible neutron KE shape (visible = %.1f%% of primary)", pctVisible);
    DrawKEOverlay(hPrimShape, hVisShape, shapeTitle.Data(), beamPOT,
                  0.0, ("c_ke_shape_"+beam+"_"+sample).c_str(), pdf, isFirst,
                  "Neutron KE [MeV]", "Percentage of neutrons");
    delete hPrimShape; delete hVisShape;

    DrawKEvsProngE(h2KEvsE, beamPOT, ("c_ke_vs_prongE_"+beam+"_"+sample).c_str(), pdf, isFirst);

    delete hPrimAll; delete hVisAll; delete h2KEvsE;

    // ── Visible neutrons split by prong multiplicity ────────────────────────
    TH1D* h1prong = LoadH1(dir, "ke_visible_neutrons_1prong",     pot);
    TH1D* h2prong = LoadH1(dir, "ke_visible_neutrons_2prong",     pot);
    TH1D* h3plus  = LoadH1(dir, "ke_visible_neutrons_3plusprong", pot);
    TH1D* hNProng = LoadH1(dir, "prongs_per_visible_neutron",     pot);

    DrawKEByProngMultiplicity(h1prong, h2prong, h3plus, beamPOT,
                              ("c_ke_multiplicity_"+beam+"_"+sample).c_str(), pdf, isFirst);
    DrawProngMultiplicity(hNProng, beamPOT,
                          ("c_nprong_"+beam+"_"+sample).c_str(), pdf, isFirst);

    delete h1prong; delete h2prong; delete h3plus; delete hNProng;

    // ── Neutron count per event: GENIE (GENIE) vs. visible ────────────────
    TH1D* hNPrimPerEvent = LoadH1(dir, "n_primary_neutrons_per_event", pot);
    TH1D* hNVisPerEvent  = LoadH1(dir, "n_visible_neutrons_per_event", pot);

    DrawKEOverlay(hNPrimPerEvent, hNVisPerEvent, "Primary vs. visible neutron count per event", beamPOT,
                  0.0, ("c_ncount_"+beam+"_"+sample).c_str(), pdf, isFirst,
                  "Neutrons per event");

    delete hNPrimPerEvent; delete hNVisPerEvent;

    // ── Confusion matrix: N GENIE neutrons vs. N visible neutrons ─────────
    TH2* h2Confusion = LoadH2(dir, "confusion_neutron_count", pot);
    DrawConfusionMatrix(h2Confusion, beamPOT,
                        ("c_confusion_"+beam+"_"+sample).c_str(), pdf, isFirst);
    delete h2Confusion;
  }
}

void plot_mc_spectra(const char* infile = "make_mc_spectra_bck.root",
                      const char* outpdf_prefix = "make_mc_spectra_plots")
{
  TFile* f = TFile::Open(infile);
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << infile << "\n"; return; }

  const std::vector<std::string> beams   = {"FHC", "RHC"};
  const std::vector<std::string> samples = {"sample1_q0Lo", "sample3_full"};

  // One PDF per sample (all beams for that sample together), so pages from
  // different samples never end up mixed in the same file.
  std::map<std::pair<std::string,std::string>, double> visibilityTable;
  for (const auto& sample : samples) {
    const TString pdf = TString::Format("%s_%s.pdf", outpdf_prefix, sample.c_str());
    bool isFirst = true;

    for (const auto& beam : beams)
      ProcessSample(f, beam, sample, pdf, isFirst, visibilityTable);

    if (!isFirst) {
      // Close the multi-page PDF (ROOT convention: an empty last c->Print(pdf+")")).
      TCanvas dummy;
      dummy.Print(pdf + ")");
    }
    std::cout << "Saved to " << pdf << "\n";
  }

  // Cross-sample overlay: full sample vs. low-hadronic region, one page
  // per beam, kept in its own PDF since it isn't part of either per-sample
  // file above.
  const TString overlayPdf = TString::Format("%s_prongs_per_neutron_overlay.pdf", outpdf_prefix);
  PlotProngMultiplicityOverlay(f, overlayPdf);

  f->Close();

  // ── Neutron visibility summary table ──────────────────────────────────────
  // "All Interactions" = sample3_full (no q0 cut).
  // "Low Hadronic Interactions" = sample1_q0Lo -- NOTE: as currently defined
  // in make_mc_spectra.C this cut is q0 < 0.150 GeV (150 MeV), not the
  // 100 MeV named in this table's header. Fix the cut upstream, or relabel
  // this header, to make the two agree.
  auto getPct = [&](const char* b, const char* s) -> double {
    auto it = visibilityTable.find({b, s});
    return (it != visibilityTable.end()) ? it->second : 0.0;
  };
  std::cout << "\nNeutron Visibility\t\tAll Interactions\tLow Hadronic Interactions (<100 MeV)\n";
  for (const char* beam : {"RHC", "FHC"}) {
    std::cout << beam << "\t\t\t"
               << TString::Format("%.2f%%", getPct(beam, "sample3_full")) << "\t\t\t"
               << TString::Format("%.2f%%", getPct(beam, "sample1_q0Lo")) << "\n";
  }
}
