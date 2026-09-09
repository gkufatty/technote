// Run:
// cafe -bq plot_neutron_spectra_systs_multiverse_filtered.C
//
// Reads make_neutron_spectra_systs_multiverse_filtered.root (produced by
// this folder's make_neutron_spectra_systs_multiverse_filtered.C -- the
// 26-systematic, chi2>5-filtered GENIE multiverse) and draws, for each beam
// and each of the five spectra it saved (ntrue_per_event, ngenie_per_event,
// nsel_per_event, hadvise_q0, final_state_topology), three pages. Modeled
// directly on ../../../spectrums/plot_neutron_spectra_systs_multiverse_filtered.C
// -- same technique, just reading this folder's directory layout
// (beam/relevant_multiverse/<objName> instead of
// beam/topo4_All_q0lt100MeV/<objName>) and five spectra instead of three.
//
//  Page A (shape): nominal (universe 0) overlaid with a +/-RMS envelope
//    computed per bin across all thrown universes relative to nominal
//    (RMS_b = sqrt(mean over universes u>=1 of (u_b - nominal_b)^2)),
//    each RESCALED to nominal's total first so only shape differences are
//    visible. Bottom panel shows % difference from nominal, bin by bin.
//    Title reports the shape-only chi2 and the raw chi2 (rate+shape
//    together), the raw rate change, and the shift in mean.
//
//  Page B (spaghetti, raw): every universe as a thin gray line, nominal in
//    bold black on top. Shows the RAW (rate+shape) spread, so a systematic
//    with a large rate effect (e.g. ZNormCCQE) can dominate the envelope.
//
//  Page C (spaghetti, shape only): the same universes, each individually
//    rescaled to nominal's total first -- isolates genuine bin-to-bin shape
//    variation across the whole ensemble, decluttered from rate effects.
//
// A companion CSV (neutron_spectra_systs_multiverse_filtered_metrics.csv)
// records chi2_shape/chi2_raw/rate_change/dmean for every (beam, spectrum).

#ifdef __CINT__
void plot_neutron_spectra_systs_multiverse_filtered()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/XSec/GenieMultiverseSyst.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace ana;

namespace
{
  // Rescales `h` in place to have the same total as `nominal`, and returns
  // {chi2, deltaMean} of the (now shape-only) comparison.
  std::pair<double, double> RescaleToShapeOnly(const TH1D* nominal, TH1D* h)
  {
    const int nb = nominal->GetNbinsX();
    double totNom = 0.0, totH = 0.0;
    for (int b = 1; b <= nb; ++b) { totNom += nominal->GetBinContent(b); totH += h->GetBinContent(b); }
    if (totH > 0.0) h->Scale(totNom / totH);

    double chi2 = 0.0, meanNom = 0.0, meanH = 0.0;
    for (int b = 1; b <= nb; ++b) {
      const double n = nominal->GetBinContent(b);
      const double s = h->GetBinContent(b);
      const double x = nominal->GetBinCenter(b);
      if (n > 0.0) chi2 += (s - n) * (s - n) / n;
      meanNom += x * n;
      meanH   += x * s;
    }
    const double deltaMean = (totNom > 0.0) ? (meanH - meanNom) / totNom : 0.0;
    return {chi2, deltaMean};
  }

  struct RawMetrics { double chi2 = 0.0; double rateChangePct = 0.0; };

  RawMetrics ComputeRawMetrics(const TH1D* nominal, const TH1D* h)
  {
    RawMetrics out;
    const int nb = nominal->GetNbinsX();
    double totNom = 0.0, totH = 0.0;
    for (int b = 1; b <= nb; ++b) {
      const double n = nominal->GetBinContent(b);
      const double s = h->GetBinContent(b);
      if (n > 0.0) out.chi2 += (s - n) * (s - n) / n;
      totNom += n;
      totH   += s;
    }
    out.rateChangePct = (totNom > 0.0) ? 100.0 * (totH - totNom) / totNom : 0.0;
    return out;
  }

  // Per-bin RMS across universes[1..end], relative to nominal (universe 0):
  // RMS_b = sqrt( mean_u (universe_u,b - nominal_b)^2 ). Returns
  // {upBand = nominal + RMS, downBand = nominal - RMS} (clamped at 0, since
  // a raw count can't go negative).
  std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>
  BuildRMSEnvelope(const TH1D* nominal, const std::vector<std::unique_ptr<TH1D>>& hists)
  {
    std::unique_ptr<TH1D> up(  static_cast<TH1D*>(nominal->Clone("upBand")));
    std::unique_ptr<TH1D> down(static_cast<TH1D*>(nominal->Clone("downBand")));
    const int nb = nominal->GetNbinsX();
    const int nThrows = int(hists.size()) - 1;  // exclude universe 0 (nominal itself)
    for (int b = 1; b <= nb; ++b) {
      const double n = nominal->GetBinContent(b);
      double sumsq = 0.0;
      for (size_t u = 1; u < hists.size(); ++u) {
        const double d = hists[u]->GetBinContent(b) - n;
        sumsq += d * d;
      }
      const double rms = (nThrows > 0) ? std::sqrt(sumsq / nThrows) : 0.0;
      up->SetBinContent(b,   n + rms);
      down->SetBinContent(b, std::max(0.0, n - rms));
    }
    return {std::move(up), std::move(down)};
  }

  struct Metrics
  {
    double chi2ShapeUp = 0.0, chi2ShapeDown = 0.0;
    double chi2RawUp   = 0.0, chi2RawDown   = 0.0;
    double dMeanUp      = 0.0, dMeanDown      = 0.0;
    double rateChangeUp = 0.0, rateChangeDown = 0.0;
  };

  void WriteCSVRow(std::ofstream& csv, const std::string& beam, const std::string& spec, const Metrics& m)
  {
    csv << beam << ',' << spec << ','
        << m.chi2ShapeUp << ',' << m.chi2ShapeDown << ','
        << m.chi2RawUp   << ',' << m.chi2RawDown   << ','
        << m.dMeanUp      << ',' << m.dMeanDown      << ','
        << m.rateChangeUp << ',' << m.rateChangeDown << '\n';
  }

  // For final_state_topology_multiverse, the bin index alone is opaque --
  // see NeutronMultSystDiagnostics.cxx's kFinalStateTopologyBinVar.
  void ApplyBinLabels(TH1D* h, const std::vector<std::string>& labels)
  {
    for (size_t i = 0; i < labels.size() && (int)i < h->GetNbinsX(); ++i)
      h->GetXaxis()->SetBinLabel((int)i + 1, labels[i].c_str());
  }

  // ── NOvA watermark and beam/POT label, matching
  // nu_interactions_plot_topology_only.C's DrawWatermark/DrawBeamLabel so
  // every technote plot carries the same identification. ─────────────────
  void DrawWatermark()
  {
    TLatex t; t.SetNDC(); t.SetTextSize(0.031); t.SetTextColor(kGray+1);
    t.DrawLatex(0.2, 0.805, "NOvA ND Simulation");
    t.DrawLatex(0.2, 0.780, "Work In Progress");
  }

  void DrawBeamLabel(const std::string& text)
  {
    if (text.empty()) return;
    TLatex t; t.SetNDC(); t.SetTextSize(0.042); t.SetTextFont(62);
    t.SetTextAlign(31);
    t.DrawLatex(0.93, 0.92, text.c_str());
  }

  std::string BeamPOTLabel(const std::string& beam, double pot)
  {
    return TString::Format("Prod 5.1 - %s POT %.0fe20", beam.c_str(), pot / 1e20).Data();
  }

  void DrawShapePage(TCanvas& c, const std::string& pdfName, const std::string& beam,
                      const std::string& specLabel, double pot,
                      TH1D* nomSrc, TH1D* upRawSrc, TH1D* downRawSrc, Metrics& m,
                      const std::vector<std::string>& binLabels)
  {
    // Clone before rescaling -- RescaleToShapeOnly mutates in place, and the
    // raw metrics below must see the un-rescaled counts.
    std::unique_ptr<TH1D> up(  static_cast<TH1D*>(upRawSrc->Clone()));
    std::unique_ptr<TH1D> down(static_cast<TH1D*>(downRawSrc->Clone()));

    const RawMetrics rawUp   = ComputeRawMetrics(nomSrc, upRawSrc);
    const RawMetrics rawDown = ComputeRawMetrics(nomSrc, downRawSrc);
    m.chi2RawUp   = rawUp.chi2;    m.rateChangeUp   = rawUp.rateChangePct;
    m.chi2RawDown = rawDown.chi2;  m.rateChangeDown = rawDown.rateChangePct;

    std::tie(m.chi2ShapeUp,   m.dMeanUp)   = RescaleToShapeOnly(nomSrc, up.get());
    std::tie(m.chi2ShapeDown, m.dMeanDown) = RescaleToShapeOnly(nomSrc, down.get());

    std::unique_ptr<TH1D> nom(static_cast<TH1D*>(nomSrc->Clone()));
    ApplyBinLabels(nom.get(), binLabels);

    // Bottom-panel ratio: % difference from nominal, bin by bin.
    std::unique_ptr<TH1D> ratioUp(static_cast<TH1D*>(nom->Clone()));
    std::unique_ptr<TH1D> ratioDown(static_cast<TH1D*>(nom->Clone()));
    for (int b = 1; b <= nom->GetNbinsX(); ++b) {
      const double n = nom->GetBinContent(b);
      ratioUp->SetBinContent(b,   n > 0.0 ? 100.0 * (up->GetBinContent(b)   - n) / n : 0.0);
      ratioDown->SetBinContent(b, n > 0.0 ? 100.0 * (down->GetBinContent(b) - n) / n : 0.0);
    }

    c.Clear();
    c.cd();
    TPad pad1("pad1", "", 0.0, 0.32, 1.0, 1.0);
    TPad pad2("pad2", "", 0.0, 0.0,  1.0, 0.32);
    pad1.SetBottomMargin(0.02);
    pad1.SetTopMargin(0.16);
    pad2.SetTopMargin(0.02);
    pad2.SetBottomMargin(0.30);
    pad1.Draw();
    pad2.Draw();

    pad1.cd();
    nom->SetLineColor(kBlack); nom->SetLineWidth(2); nom->SetLineStyle(1);
    up->SetLineColor(kRed);    up->SetLineWidth(2);  up->SetLineStyle(2);
    down->SetLineColor(kBlue); down->SetLineWidth(2); down->SetLineStyle(2);

    const double ymax = std::max({nom->GetMaximum(), up->GetMaximum(), down->GetMaximum()});
    nom->SetMaximum(ymax * 1.3);
    nom->SetMinimum(0.0);
    nom->SetTitle(TString::Format(
      "#splitline{%s: %s -- GENIE multiverse (26 relevant systs, shape only)  "
      "#chi^{2}_{shape,+RMS}=%.2f  #chi^{2}_{shape,-RMS}=%.2f  "
      "#Delta#LT x#GT_{+RMS}=%+.4f  #Delta#LT x#GT_{-RMS}=%+.4f}"
      "{#chi^{2}_{raw,+RMS}=%.2f  #chi^{2}_{raw,-RMS}=%.2f   "
      "rate: %+.2f%% / %+.2f%%}",
      beam.c_str(), specLabel.c_str(), m.chi2ShapeUp, m.chi2ShapeDown, m.dMeanUp, m.dMeanDown,
      m.chi2RawUp, m.chi2RawDown, m.rateChangeUp, m.rateChangeDown));
    nom->GetXaxis()->SetLabelSize(0.0);

    nom->Draw("hist");
    up->Draw("hist same");
    down->Draw("hist same");

    TLegend leg(0.56, 0.68, 0.88, 0.86);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(nom.get(), "nominal (universe 0)",  "l");
    leg.AddEntry(up.get(),  "nominal+RMS (shape only)", "l");
    leg.AddEntry(down.get(),"nominal-RMS (shape only)", "l");
    leg.Draw();

    DrawBeamLabel(BeamPOTLabel(beam, pot));
    DrawWatermark();

    pad2.cd();
    double axMax = 1.0;
    for (int b = 1; b <= nom->GetNbinsX(); ++b)
      axMax = std::max({axMax, std::abs(ratioUp->GetBinContent(b)), std::abs(ratioDown->GetBinContent(b))});

    ratioUp->SetTitle("");
    ratioUp->SetMaximum(axMax * 1.3);
    ratioUp->SetMinimum(-axMax * 1.3);
    ratioUp->GetYaxis()->SetTitle("% diff from nominal");
    ratioUp->GetYaxis()->SetTitleSize(0.10);
    ratioUp->GetYaxis()->SetTitleOffset(0.45);
    ratioUp->GetYaxis()->SetLabelSize(0.09);
    ratioUp->GetYaxis()->SetNdivisions(505);
    ratioUp->GetXaxis()->SetTitle(specLabel.c_str());
    ratioUp->GetXaxis()->SetLabelSize(0.10);
    ratioUp->GetXaxis()->SetTitleSize(0.10);
    ratioUp->SetLineColor(kRed);     ratioUp->SetLineWidth(2);
    ratioUp->SetMarkerColor(kRed);   ratioUp->SetMarkerStyle(20);
    ratioDown->SetLineColor(kBlue);  ratioDown->SetLineWidth(2);
    ratioDown->SetMarkerColor(kBlue); ratioDown->SetMarkerStyle(21);

    ratioUp->Draw("PL");
    ratioDown->Draw("PL same");

    TLine zero(ratioUp->GetXaxis()->GetXmin(), 0.0, ratioUp->GetXaxis()->GetXmax(), 0.0);
    zero.SetLineColor(kBlack);
    zero.SetLineStyle(2);
    zero.Draw();

    c.cd();
    c.Print(pdfName.c_str());
  }

  void DrawSpaghettiPage(TCanvas& c, const std::string& pdfName, const std::string& beam,
                        const std::string& specLabel, double pot,
                        const std::vector<std::unique_ptr<TH1D>>& hists, bool shapeOnly,
                        const std::vector<std::string>& binLabels)
  {
    c.Clear();
    c.cd();
    c.SetLeftMargin(0.15); c.SetBottomMargin(0.14); c.SetTopMargin(0.12);

    TH1D* nom = hists[0].get();
    ApplyBinLabels(nom, binLabels);
    double ymax = 0.0;
    for (auto& h : hists) ymax = std::max(ymax, h->GetMaximum());

    nom->SetTitle(TString::Format(
      "%s: %s GENIE multiverse (26 relevant systs%s, %zu universes)",
      beam.c_str(), specLabel.c_str(),
      shapeOnly ? ", shape only -- each universe rescaled to nominal's total" : ", raw",
      hists.size()));
    nom->SetMaximum(ymax * 1.3);
    nom->SetMinimum(0.0);
    nom->SetLineColor(kBlack);
    nom->SetLineWidth(1);
    nom->Draw("hist");

    for (size_t u = 1; u < hists.size(); ++u) {
      hists[u]->SetLineColor(kGray + 1);
      hists[u]->SetLineWidth(1);
      hists[u]->Draw("hist same");
    }

    // redraw nominal on top, bold, so it isn't buried under the spaghetti
    nom->SetLineWidth(3);
    nom->Draw("hist same");

    TLegend leg(0.55, 0.78, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(nom, "nominal", "l");
    if (hists.size() > 1)
      leg.AddEntry(hists[1].get(), shapeOnly ? "GENIE universes (shape only)" : "GENIE universes (raw)", "l");
    leg.Draw();

    DrawBeamLabel(BeamPOTLabel(beam, pot));
    DrawWatermark();

    c.Print(pdfName.c_str());
  }

  void PlotSpectrum(TFile* fIn, TCanvas& c, const std::string& pdfName,
                     const std::string& beam, const std::string& objName,
                     const std::string& specLabel, std::ofstream& csv,
                     const std::vector<std::string>& binLabels)
  {
    TDirectory* topDir = fIn->GetDirectory((beam + "/relevant_multiverse").c_str());
    if (!topDir) {
      std::cerr << "[WARN] no directory for beam " << beam << "\n";
      return;
    }

    std::unique_ptr<GenieMultiverseSpectra> multiverse =
      GenieMultiverseSpectra::LoadFrom(topDir, objName.c_str());
    if (!multiverse) {
      std::cerr << "[WARN] could not load " << objName << " for " << beam << "\n";
      return;
    }

    const double pot = multiverse->Nominal()->POT();
    const std::vector<const Spectrum*> universes = multiverse->AllUniverses();
    std::cout << "[" << beam << "] " << objName << ": plotting " << universes.size() << " universes...\n";

    std::vector<std::unique_ptr<TH1D>> hists;
    hists.reserve(universes.size());
    for (const Spectrum* u : universes) hists.emplace_back(u->ToTH1(pot));

    TH1D* nom = hists[0].get();
    std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>> band = BuildRMSEnvelope(nom, hists);
    TH1D* upRaw = band.first.get();
    TH1D* downRaw = band.second.get();

    Metrics m;
    DrawShapePage(c, pdfName, beam, specLabel, pot, nom, upRaw, downRaw, m, binLabels);
    WriteCSVRow(csv, beam, specLabel, m);

    DrawSpaghettiPage(c, pdfName, beam, specLabel, pot, hists, /*shapeOnly=*/false, binLabels);

    // Shape-only spaghetti: clone every universe and rescale each one
    // individually to nominal's total, so the spread shown here is genuine
    // bin-to-bin shape variation, not rate riding along with it.
    std::vector<std::unique_ptr<TH1D>> histsShape;
    histsShape.reserve(hists.size());
    for (const auto& h : hists) histsShape.emplace_back(static_cast<TH1D*>(h->Clone()));
    for (auto& h : histsShape) RescaleToShapeOnly(nom, h.get());

    DrawSpaghettiPage(c, pdfName, beam, specLabel, pot, histsShape, /*shapeOnly=*/true, binLabels);
  }
}

void plot_neutron_spectra_systs_multiverse_filtered()
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.032);
  TH1::AddDirectory(kFALSE);
  // Force scientific notation on event-count y-axes, same as
  // nu_interactions_plot_topology_only.C.
  TGaxis::SetMaxDigits(3);

  TFile* fIn = TFile::Open("make_neutron_spectra_systs_multiverse_filtered.root", "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERR] could not open make_neutron_spectra_systs_multiverse_filtered.root -- "
                 "run make_neutron_spectra_systs_multiverse_filtered.C first\n";
    return;
  }

  const std::string pdfName = "neutron_spectra_systs_multiverse_filtered_plots.pdf";
  const std::string csvName = "neutron_spectra_systs_multiverse_filtered_metrics.csv";
  TCanvas c("c", "c", 800, 700);

  // Same short labels used in plot_neutron_spectra_systs_diagnostics_ranking.C
  // and plot_relevant_systs_ranking.C for the categorical topology axis.
  const std::vector<std::string> kTopologyBinLabels = {
    "0p 1n", "0p #geq2n", "0p 1n 1#pi", "1p 0#pi",
    "#geq2p 0#pi", "1p N#pi", "#geq2p N#pi", "other/0n",
  };
  const std::vector<std::string> kNoLabels;

  std::ofstream csv(csvName);
  csv << "beam,spectrum,"
         "chi2_shape_up,chi2_shape_down,"
         "chi2_raw_up,chi2_raw_down,"
         "dmean_up,dmean_down,"
         "rate_change_pct_up,rate_change_pct_down\n";

  c.Print((pdfName + "[").c_str());
  for (const std::string& beam : {"FHC", "RHC"}) {
    PlotSpectrum(fIn, c, pdfName, beam, "ntrue_per_event_multiverse",     "N true neutron prongs",   csv, kNoLabels);
    PlotSpectrum(fIn, c, pdfName, beam, "ngenie_per_event_multiverse",    "N GENIE neutrons",        csv, kNoLabels);
    PlotSpectrum(fIn, c, pdfName, beam, "nsel_per_event_multiverse",      "N_{NLS} per event",       csv, kNoLabels);
    PlotSpectrum(fIn, c, pdfName, beam, "hadvise_q0_multiverse",          "HadVisE [MeV]",           csv, kNoLabels);
    PlotSpectrum(fIn, c, pdfName, beam, "final_state_topology_multiverse","Final-state topology bin",csv, kTopologyBinLabels);
  }
  c.Print((pdfName + "]").c_str());

  csv.close();
  fIn->Close();
  std::cout << "Saved shape + spaghetti plots to " << pdfName << "\n";
  std::cout << "Saved metrics table to " << csvName << "\n";
}

#endif
