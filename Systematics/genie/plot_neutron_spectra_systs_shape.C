// Run:
// cafe -bq plot_neutron_spectra_systs_shape.C
//
// Reads make_neutron_spectra_systs.root and, for each of the two variables
// it contains (ntrue_per_event / N true neutron prongs, ngenie_per_event /
// N GENIE neutrons -- see make_neutron_spectra_systs.C) and each beam and
// each GENIE systematic, draws two pages:
//
//  Page A (shape): nominal ntrue_per_event overlaid with its +1 sigma and
//    -1 sigma shifted versions, each shifted histogram RESCALED to
//    nominal's total first so only shape differences are visible, plus a
//    bottom panel of % difference from nominal, bin by bin. Title reports
//    both the shape-only chi2 (sum (s_i-n_i)^2/n_i after rescaling to
//    nominal's total -- isolates shape distortion) and the raw chi2 (same
//    formula, no rescaling -- captures shape AND rate together), the
//    resulting shift in mean multiplicity <N>, the raw rate change, and
//    the shift in variance of the (shape-only) distribution.
//
//  Page B (diagnostics): cumulative distribution (CDF) of nominal vs the
//    shape-only +/-1 sigma shifts, with the KS-type max-CDF-distance
//    reported per shift; and a bottom panel of the per-bin "curvature"
//    (up + down - 2*nominal, on RAW un-rescaled counts) -- a linearity/
//    asymmetry check: zero everywhere means the +1/-1 sigma response is a
//    perfect mirror image of nominal, nonzero flags a bin where the
//    response is asymmetric or nonlinear.
//
// Pages are ordered largest-shape-chi2-first (same convention as
// plot_neutron_spectra_systs_chi2ranking.C). Each variable gets its own
// multi-page PDF and companion CSV: ntrue_per_event keeps the original
// filenames (neutron_spectra_systs_shape_plots.pdf /
// neutron_spectra_systs_metrics.csv) unchanged, since that exact CSV is
// parsed by neutron_spectra.py, neutron_spectra_extended.C, and read by
// make_neutron_spectra_systs_multiverse_filtered.C /
// make_neutron_spectra_mode_comparison.C for their top-N systematic lists.
// ngenie_per_event writes to new, separate files
// (neutron_spectra_systs_shape_plots_ngenie.pdf /
// neutron_spectra_systs_metrics_ngenie.csv) so none of those consumers are
// affected.

#ifdef __CINT__
void plot_neutron_spectra_systs_shape()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "CAFAna/Core/Spectrum.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
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
  std::vector<std::string> ListSubdirNames(TDirectory* dir)
  {
    std::vector<std::string> names;
    if (!dir) return names;
    TIter next(dir->GetListOfKeys());
    while (TKey* key = dynamic_cast<TKey*>(next())) {
      if (std::string(key->GetClassName()) == "TDirectoryFile")
        names.push_back(key->GetName());
    }
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
  }

  TH1D* LoadHist(TDirectory* topDir, const std::string& subpath, const std::string& objName,
                 double* potOut = nullptr)
  {
    TDirectory* parentDir = topDir->GetDirectory(subpath.c_str());
    if (!parentDir) return nullptr;
    std::unique_ptr<Spectrum> spec = Spectrum::LoadFrom(parentDir, objName.c_str());
    const double pot = spec->POT();
    if (potOut) *potOut = pot;
    return spec->ToTH1(pot);
  }

  // ── NOvA watermark and beam/POT label, matching
  // nu_interactions_plot_topology_only.C's DrawWatermark/DrawBeamLabel so
  // every technote plot carries the same identification. ─────────────────
  void DrawWatermark()
  {
    TLatex t; t.SetNDC(); t.SetTextSize(0.031); t.SetTextColor(kGray+1);
    t.DrawLatex(0.2, 0.755, "NOvA ND Simulation");
    t.DrawLatex(0.2, 0.730, "Work In Progress");
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

  // Same chi2 formula as RescaleToShapeOnly, but with NO rescaling first --
  // this bakes in both the rate (normalization) change and the shape
  // change, unlike the shape-only chi2 above.
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

  // Population variance of the multiplicity distribution `h` (treated as a
  // discrete distribution over its bin centers).
  double Variance(const TH1D* h)
  {
    const int nb = h->GetNbinsX();
    double tot = 0.0, mean = 0.0;
    for (int b = 1; b <= nb; ++b) { const double c = h->GetBinContent(b); tot += c; mean += c * h->GetBinCenter(b); }
    if (tot <= 0.0) return 0.0;
    mean /= tot;
    double var = 0.0;
    for (int b = 1; b <= nb; ++b) {
      const double c  = h->GetBinContent(b);
      const double dx = h->GetBinCenter(b) - mean;
      var += c * dx * dx;
    }
    return var / tot;
  }

  // Max absolute difference between the cumulative distributions of `h`
  // and `nominal` (each normalized to its own total), in percentage
  // points -- a KS-type distance for the shape difference.
  double KSDistance(const TH1D* nominal, const TH1D* h)
  {
    const int nb = nominal->GetNbinsX();
    double totNom = 0.0, totH = 0.0;
    for (int b = 1; b <= nb; ++b) { totNom += nominal->GetBinContent(b); totH += h->GetBinContent(b); }
    if (totNom <= 0.0 || totH <= 0.0) return 0.0;
    double cumNom = 0.0, cumH = 0.0, maxDiff = 0.0;
    for (int b = 1; b <= nb; ++b) {
      cumNom += nominal->GetBinContent(b) / totNom;
      cumH   += h->GetBinContent(b)       / totH;
      maxDiff = std::max(maxDiff, std::abs(cumH - cumNom));
    }
    return 100.0 * maxDiff;
  }

  // Linearity/asymmetry check on RAW (un-rescaled) counts: if the +1/-1
  // sigma response were a perfect linear mirror image of nominal, then
  // up_b + down_b - 2*nom_b = 0 in every bin. Returns the sum of |curvature|
  // over all bins, as a % of the nominal total.
  double AsymmetryMetric(const TH1D* nominal, const TH1D* upRaw, const TH1D* downRaw)
  {
    const int nb = nominal->GetNbinsX();
    double totNom = 0.0, sumAbsCurv = 0.0;
    for (int b = 1; b <= nb; ++b) {
      const double n = nominal->GetBinContent(b);
      const double u = upRaw->GetBinContent(b);
      const double d = downRaw->GetBinContent(b);
      sumAbsCurv += std::abs(u + d - 2.0 * n);
      totNom += n;
    }
    return (totNom > 0.0) ? 100.0 * sumAbsCurv / totNom : 0.0;
  }

  std::unique_ptr<TH1D> BuildCDF(const TH1D* h, const char* name)
  {
    std::unique_ptr<TH1D> cdf(static_cast<TH1D*>(h->Clone(name)));
    const int nb = h->GetNbinsX();
    double tot = 0.0;
    for (int b = 1; b <= nb; ++b) tot += h->GetBinContent(b);
    double cum = 0.0;
    for (int b = 1; b <= nb; ++b) {
      cum += h->GetBinContent(b);
      cdf->SetBinContent(b, tot > 0.0 ? cum / tot : 0.0);
    }
    return cdf;
  }

  // Per-bin curvature (%), on RAW un-rescaled counts: 100*(up+down-2*nom)/nom.
  std::unique_ptr<TH1D> BuildCurvature(const TH1D* nominal, const TH1D* upRaw,
                                       const TH1D* downRaw, const char* name)
  {
    std::unique_ptr<TH1D> curv(static_cast<TH1D*>(nominal->Clone(name)));
    const int nb = nominal->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
      const double n = nominal->GetBinContent(b);
      const double u = upRaw->GetBinContent(b);
      const double d = downRaw->GetBinContent(b);
      curv->SetBinContent(b, n > 0.0 ? 100.0 * (u + d - 2.0 * n) / n : 0.0);
    }
    return curv;
  }

  struct RankedSyst
  {
    std::string name;
    double chi2ShapeUp = 0.0, chi2ShapeDown = 0.0;   // shape-only (rescaled to nominal total)
    double chi2RawUp   = 0.0, chi2RawDown   = 0.0;    // raw (rate + shape together)
    double dMeanUp      = 0.0, dMeanDown      = 0.0;   // shift in <N> (shape only)
    double rateChangeUp = 0.0, rateChangeDown = 0.0;   // % change in total events (raw)
    double dVarUp        = 0.0, dVarDown        = 0.0; // shift in variance of N (shape only)
    double ksUp           = 0.0, ksDown           = 0.0; // KS-type max CDF distance, % (shape only)
    double asymmetry      = 0.0;                         // linearity/asymmetry check, % (raw)
  };

  std::vector<RankedSyst> RankSysts(TDirectory* topDir, const TH1D* nominal, const std::string& objName)
  {
    std::vector<RankedSyst> out;
    TDirectory* systsDir = topDir->GetDirectory("systs");
    const double varNom = Variance(nominal);

    for (const std::string& name : ListSubdirNames(systsDir)) {
      std::unique_ptr<TH1D> upRaw(  LoadHist(systsDir, name + "/p1sigma", objName));
      std::unique_ptr<TH1D> downRaw(LoadHist(systsDir, name + "/m1sigma", objName));
      if (!upRaw || !downRaw) continue;

      RankedSyst r;
      r.name = name;

      const RawMetrics rawUp   = ComputeRawMetrics(nominal, upRaw.get());
      const RawMetrics rawDown = ComputeRawMetrics(nominal, downRaw.get());
      r.chi2RawUp     = rawUp.chi2;           r.rateChangeUp   = rawUp.rateChangePct;
      r.chi2RawDown   = rawDown.chi2;         r.rateChangeDown = rawDown.rateChangePct;
      r.asymmetry = AsymmetryMetric(nominal, upRaw.get(), downRaw.get());

      // Clone before rescaling -- RescaleToShapeOnly mutates in place, and
      // the raw metrics above must see the un-rescaled counts.
      std::unique_ptr<TH1D> up(static_cast<TH1D*>(upRaw->Clone()));
      std::unique_ptr<TH1D> down(static_cast<TH1D*>(downRaw->Clone()));
      std::tie(r.chi2ShapeUp,   r.dMeanUp)   = RescaleToShapeOnly(nominal, up.get());
      std::tie(r.chi2ShapeDown, r.dMeanDown) = RescaleToShapeOnly(nominal, down.get());

      r.dVarUp   = Variance(up.get())   - varNom;
      r.dVarDown = Variance(down.get()) - varNom;
      r.ksUp   = KSDistance(nominal, up.get());
      r.ksDown = KSDistance(nominal, down.get());

      out.push_back(r);
    }
    std::sort(out.begin(), out.end(), [](const RankedSyst& a, const RankedSyst& b) {
      return std::max(a.chi2ShapeUp, a.chi2ShapeDown) > std::max(b.chi2ShapeUp, b.chi2ShapeDown);
    });
    return out;
  }

  void WriteCSVRow(std::ofstream& csv, const std::string& beam, const RankedSyst& r)
  {
    csv << beam << ',' << r.name << ','
        << r.chi2ShapeUp << ',' << r.chi2ShapeDown << ','
        << r.chi2RawUp   << ',' << r.chi2RawDown   << ','
        << r.dMeanUp      << ',' << r.dMeanDown      << ','
        << r.rateChangeUp << ',' << r.rateChangeDown << ','
        << r.dVarUp        << ',' << r.dVarDown        << ','
        << r.ksUp           << ',' << r.ksDown           << ','
        << r.asymmetry << '\n';
  }

  void DrawShapePage(TCanvas& c, const std::string& pdfName,
                      const std::string& beam, double pot, const RankedSyst& r,
                      TH1D* nomSrc, TDirectory* systsDir,
                      const std::string& objName, const std::string& varLabel)
  {
    std::unique_ptr<TH1D> up(  LoadHist(systsDir, r.name + "/p1sigma", objName));
    std::unique_ptr<TH1D> down(LoadHist(systsDir, r.name + "/m1sigma", objName));
    if (!up || !down) return;
    RescaleToShapeOnly(nomSrc, up.get());
    RescaleToShapeOnly(nomSrc, down.get());

    // nomSrc is reused across every page for this beam; clone it so we
    // never mutate/rescale the shared object, and so a stale SetMaximum()
    // from a previous page can't leak into this one's axis range.
    std::unique_ptr<TH1D> nom(static_cast<TH1D*>(nomSrc->Clone()));

    // Bottom-panel ratio: % difference from nominal, bin by bin. This is
    // the actual signal -- on the linear absolute-count scale above, bin 0
    // (tens of thousands of events) dwarfs bins 1-4 (hundreds), so a real
    // 5-6% shift in a small bin is invisible there. Down here it isn't.
    std::unique_ptr<TH1D> ratioUp(static_cast<TH1D*>(nom->Clone()));
    std::unique_ptr<TH1D> ratioDown(static_cast<TH1D*>(nom->Clone()));
    for (int b = 1; b <= nom->GetNbinsX(); ++b) {
      const double n = nom->GetBinContent(b);
      ratioUp->SetBinContent(b,   n > 0.0 ? 100.0 * (up->GetBinContent(b)   - n) / n : 0.0);
      ratioDown->SetBinContent(b, n > 0.0 ? 100.0 * (down->GetBinContent(b) - n) / n : 0.0);
    }

    c.Clear();
    c.cd();  // pad1/pad2 below are locals destroyed at the end of every call;
              // make sure they attach to c, not a stale gPad from the last one
    TPad pad1("pad1", "", 0.0, 0.32, 1.0, 1.0);
    TPad pad2("pad2", "", 0.0, 0.0,  1.0, 0.32);
    pad1.SetBottomMargin(0.02);
    pad1.SetTopMargin(0.16);  // two-line title needs more headroom than the default
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
      "#splitline{%s: %s (shape only)  "
      "#chi^{2}_{shape,+1#sigma}=%.2f  #chi^{2}_{shape,-1#sigma}=%.2f  "
      "#Delta#LT N#GT_{+1#sigma}=%+.4f  #Delta#LT N#GT_{-1#sigma}=%+.4f}"
      "{#chi^{2}_{raw,+1#sigma}=%.2f  #chi^{2}_{raw,-1#sigma}=%.2f   "
      "rate: %+.2f%% / %+.2f%%   #Delta#sigma^{2}_{N}: %+.4f / %+.4f}",
      beam.c_str(), r.name.c_str(), r.chi2ShapeUp, r.chi2ShapeDown, r.dMeanUp, r.dMeanDown,
      r.chi2RawUp, r.chi2RawDown, r.rateChangeUp, r.rateChangeDown, r.dVarUp, r.dVarDown));
    nom->GetXaxis()->SetLabelSize(0.0);  // bottom panel carries the bin labels

    nom->Draw("hist");
    up->Draw("hist same");
    down->Draw("hist same");

    TLegend leg(0.60, 0.68, 0.88, 0.86);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(nom.get(), "nominal",              "l");
    leg.AddEntry(up.get(),  "+1#sigma (shape only)", "l");
    leg.AddEntry(down.get(),"-1#sigma (shape only)", "l");
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
    ratioUp->GetXaxis()->SetTitle(varLabel.c_str());
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

  void DrawDiagnosticsPage(TCanvas& c, const std::string& pdfName,
                            const std::string& beam, double pot, const RankedSyst& r,
                            TH1D* nomSrc, TDirectory* systsDir,
                            const std::string& objName, const std::string& varLabel)
  {
    std::unique_ptr<TH1D> upRaw(  LoadHist(systsDir, r.name + "/p1sigma", objName));
    std::unique_ptr<TH1D> downRaw(LoadHist(systsDir, r.name + "/m1sigma", objName));
    if (!upRaw || !downRaw) return;

    std::unique_ptr<TH1D> nom(static_cast<TH1D*>(nomSrc->Clone()));

    // Shape-only clones for the CDF, so it isolates shape exactly like
    // page A -- an overall rate change would otherwise dominate the CDF
    // gap from bin 0 onward and hide the shape effect underneath it.
    std::unique_ptr<TH1D> upShape(static_cast<TH1D*>(upRaw->Clone()));
    std::unique_ptr<TH1D> downShape(static_cast<TH1D*>(downRaw->Clone()));
    RescaleToShapeOnly(nom.get(), upShape.get());
    RescaleToShapeOnly(nom.get(), downShape.get());

    std::unique_ptr<TH1D> cdfNom(  BuildCDF(nom.get(),      "cdfNom"));
    std::unique_ptr<TH1D> cdfUp(   BuildCDF(upShape.get(),  "cdfUp"));
    std::unique_ptr<TH1D> cdfDown( BuildCDF(downShape.get(),"cdfDown"));

    // Curvature uses RAW counts on purpose (see AsymmetryMetric comment):
    // this is a check on the response itself, not on the shape-only view.
    std::unique_ptr<TH1D> curv(BuildCurvature(nom.get(), upRaw.get(), downRaw.get(), "curv"));

    c.Clear();
    c.cd();
    TPad pad1("pad1d", "", 0.0, 0.32, 1.0, 1.0);
    TPad pad2("pad2d", "", 0.0, 0.0,  1.0, 0.32);
    pad1.SetBottomMargin(0.02);
    pad1.SetTopMargin(0.16);
    pad2.SetTopMargin(0.02);
    pad2.SetBottomMargin(0.30);
    pad1.Draw();
    pad2.Draw();

    pad1.cd();
    cdfNom->SetLineColor(kBlack); cdfNom->SetLineWidth(2);
    cdfUp->SetLineColor(kRed);    cdfUp->SetLineWidth(2);  cdfUp->SetLineStyle(2);
    cdfDown->SetLineColor(kBlue); cdfDown->SetLineWidth(2); cdfDown->SetLineStyle(2);
    cdfNom->SetMaximum(1.15);
    cdfNom->SetMinimum(0.0);
    cdfNom->GetXaxis()->SetLabelSize(0.0);
    cdfNom->GetYaxis()->SetTitle("cumulative fraction");
    cdfNom->SetTitle(TString::Format(
      "#splitline{%s: %s -- cumulative distribution (shape only)}"
      "{KS_{+1#sigma}=%.2f%%   KS_{-1#sigma}=%.2f%%   "
      "linearity/asymmetry (raw curvature) = %.2f%%}",
      beam.c_str(), r.name.c_str(), r.ksUp, r.ksDown, r.asymmetry));

    cdfNom->Draw("hist");
    cdfUp->Draw("hist same");
    cdfDown->Draw("hist same");

    TLegend leg(0.60, 0.18, 0.88, 0.36);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(cdfNom.get(), "nominal",              "l");
    leg.AddEntry(cdfUp.get(),  "+1#sigma (shape only)", "l");
    leg.AddEntry(cdfDown.get(),"-1#sigma (shape only)", "l");
    leg.Draw();

    DrawBeamLabel(BeamPOTLabel(beam, pot));
    DrawWatermark();

    pad2.cd();
    double axMax = 1.0;
    for (int b = 1; b <= curv->GetNbinsX(); ++b)
      axMax = std::max(axMax, std::abs(curv->GetBinContent(b)));

    curv->SetTitle("");
    curv->SetMaximum(axMax * 1.3);
    curv->SetMinimum(-axMax * 1.3);
    curv->GetYaxis()->SetTitle("curvature (raw): (up+down-2#timesnom)/nom [%]");
    curv->GetYaxis()->SetTitleSize(0.09);
    curv->GetYaxis()->SetTitleOffset(0.5);
    curv->GetYaxis()->SetLabelSize(0.09);
    curv->GetYaxis()->SetNdivisions(505);
    curv->GetXaxis()->SetTitle(varLabel.c_str());
    curv->GetXaxis()->SetLabelSize(0.10);
    curv->GetXaxis()->SetTitleSize(0.10);
    curv->SetLineColor(kMagenta + 2);
    curv->SetLineWidth(2);
    curv->SetMarkerColor(kMagenta + 2);
    curv->SetMarkerStyle(20);
    curv->Draw("PL");

    TLine zero(curv->GetXaxis()->GetXmin(), 0.0, curv->GetXaxis()->GetXmax(), 0.0);
    zero.SetLineColor(kBlack);
    zero.SetLineStyle(2);
    zero.Draw();

    c.cd();
    c.Print(pdfName.c_str());
  }

  void PlotBeam(TFile* fIn, TCanvas& c, const std::string& pdfName,
                const std::string& beam, std::ofstream& csv,
                const std::string& objName, const std::string& varLabel)
  {
    // NOTE: make_neutron_spectra_systs.C writes this cut region under
    // "HadVisE_100MeV" (see ProcessBeam there), not "topo4_All_q0lt100MeV".
    TDirectory* topDir = fIn->GetDirectory((beam + "/HadVisE_100MeV").c_str());
    if (!topDir) {
      std::cerr << "[WARN] no directory for beam " << beam << "\n";
      return;
    }

    double pot = 0.0;
    std::unique_ptr<TH1D> nominal(LoadHist(topDir, "nominal", objName, &pot));
    if (!nominal) return;

    const std::vector<RankedSyst> ranked = RankSysts(topDir, nominal.get(), objName);
    std::cout << "[" << beam << "] " << objName << ": plotting " << ranked.size()
              << " systematics, largest shape chi2 first...\n";

    TDirectory* systsDir = topDir->GetDirectory("systs");
    for (const RankedSyst& r : ranked) {
      DrawShapePage(c, pdfName, beam, pot, r, nominal.get(), systsDir, objName, varLabel);
      DrawDiagnosticsPage(c, pdfName, beam, pot, r, nominal.get(), systsDir, objName, varLabel);
      WriteCSVRow(csv, beam, r);
    }
  }
}

void plot_neutron_spectra_systs_shape()
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.032);
  gStyle->SetTitleY(0.97);  // nudge the splitline chi2/rate title down from the pad's top edge  // two-line splitline titles need headroom, not size
  TH1::AddDirectory(kFALSE);
  // Force scientific notation on event-count y-axes, same as
  // nu_interactions_plot_topology_only.C.
  TGaxis::SetMaxDigits(3);

  TFile* fIn = TFile::Open("make_neutron_spectra_systs.root", "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERR] could not open make_neutron_spectra_systs.root\n";
    return;
  }

  TCanvas c("c", "c", 800, 700);

  // {object name in the ROOT file, axis label, output PDF, output CSV}.
  // ntrue_per_event keeps its original filenames untouched -- see header
  // comment for why. ngenie_per_event is purely additive, new filenames.
  struct VarSpec { std::string objName, varLabel, pdfName, csvName; };
  const std::vector<VarSpec> variables = {
    {"ntrue_per_event",  "N true neutron prongs",
     "neutron_spectra_systs_shape_plots.pdf",        "neutron_spectra_systs_metrics.csv"},
    {"ngenie_per_event", "N GENIE neutrons",
     "neutron_spectra_systs_shape_plots_ngenie.pdf", "neutron_spectra_systs_metrics_ngenie.csv"},
  };

  for (const VarSpec& v : variables) {
    std::ofstream csv(v.csvName);
    csv << "beam,syst,"
           "chi2_shape_up,chi2_shape_down,"
           "chi2_raw_up,chi2_raw_down,"
           "dmean_up,dmean_down,"
           "rate_change_pct_up,rate_change_pct_down,"
           "dvar_up,dvar_down,"
           "ks_pct_up,ks_pct_down,"
           "asymmetry_pct\n";

    c.Print((v.pdfName + "[").c_str());
    PlotBeam(fIn, c, v.pdfName, "FHC", csv, v.objName, v.varLabel);
    PlotBeam(fIn, c, v.pdfName, "RHC", csv, v.objName, v.varLabel);
    c.Print((v.pdfName + "]").c_str());

    csv.close();
    std::cout << "Saved shape-only + diagnostics plots to " << v.pdfName << "\n";
    std::cout << "Saved per-syst metrics table to " << v.csvName << "\n";
  }

  fIn->Close();
}

#endif
