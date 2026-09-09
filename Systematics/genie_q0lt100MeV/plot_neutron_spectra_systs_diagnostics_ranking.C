// Run:
// cafe -bq plot_neutron_spectra_systs_diagnostics_ranking.C
//
// q0lt100MeV variant of ../genie/plot_neutron_spectra_systs_diagnostics_ranking.C.
// Unchanged from that copy -- this script only reads whatever ROOT files
// and beam/category/nominal/systs structure it finds in the current
// directory, so it doesn't need to know that every spectrum here was built
// with an extra HadVisE < 100 MeV cut. Reads
// make_neutron_spectra_systs_diagnostics_{high,medium,low}.root (produced
// by this folder's own make_neutron_spectra_systs_diagnostics_{high,medium,
// low}.C) and answers, for each of the technote's five physics categories
// (FSI, MEC, RES, QE, DIS) and its recommended diagnostic variable(s), does
// shifting the category's systematics by +-1 sigma actually move that
// diagnostic within this restricted region -- and does the answer differ
// from the fully-inclusive study in ../genie/?
//
// Same metric machinery as plot_neutron_spectra_systs_shape.C (shape-only
// chi2, raw chi2, rate change, KS distance, variance shift, curvature/
// asymmetry -- see that file's header for definitions), generalized to
// loop over (category, diagnostic) pairs instead of a single fixed
// directory. For every pair: one multi-page PDF (largest-shape-chi2-first,
// same convention as the existing ranking/shape scripts) and one CSV of
// per-systematic metrics, plus a stdout summary ranking every systematic
// within its category by shape chi2 -- the direct "is this relevant?"
// readout.

#ifdef __CINT__
void plot_neutron_spectra_systs_diagnostics_ranking()
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
#include <iomanip>
#include <iostream>
#include <map>
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

  // ── NOvA watermark, matching nu_interactions_plot_topology_only.C's
  // DrawWatermark so every technote plot carries the same identification.
  // The beam/POT label used to be a separate corner label (DrawBeamLabel)
  // but is now folded directly into each page's title via BeamPOTLabel(). ──
  void DrawWatermark()
  {
    TLatex t; t.SetNDC(); t.SetTextSize(0.031); t.SetTextColor(kGray+1);
    t.DrawLatex(0.2, 0.755, "NOvA ND Simulation");
    t.DrawLatex(0.2, 0.730, "Work In Progress");
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

  struct RankedSyst
  {
    std::string name;
    double chi2ShapeUp = 0.0, chi2ShapeDown = 0.0;
    double chi2RawUp   = 0.0, chi2RawDown   = 0.0;
    double dMeanUp      = 0.0, dMeanDown      = 0.0;
    double rateChangeUp = 0.0, rateChangeDown = 0.0;
    double dVarUp        = 0.0, dVarDown        = 0.0;
    double ksUp           = 0.0, ksDown           = 0.0;
    double asymmetry      = 0.0;
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

  // For categorical diagnostics (final-state topology, prongs-per-neutron)
  // the bin index alone is opaque -- replace the numeric ticks with short
  // text labels. No-op if `labels` is empty (the default, for diagnostics
  // whose numeric bin value is already self-explanatory).
  void ApplyBinLabels(TH1D* h, const std::vector<std::string>& labels)
  {
    for (size_t i = 0; i < labels.size() && (int)i < h->GetNbinsX(); ++i)
      h->GetXaxis()->SetBinLabel((int)i + 1, labels[i].c_str());
  }

  void DrawShapePage(TCanvas& c, const std::string& pdfName,
                      const std::string& beam, double pot, const RankedSyst& r,
                      TH1D* nomSrc, TDirectory* systsDir,
                      const std::string& objName, const std::string& varLabel,
                      const std::string& category,
                      const std::vector<std::string>& binLabels)
  {
    std::unique_ptr<TH1D> up(  LoadHist(systsDir, r.name + "/p1sigma", objName));
    std::unique_ptr<TH1D> down(LoadHist(systsDir, r.name + "/m1sigma", objName));
    if (!up || !down) return;
    RescaleToShapeOnly(nomSrc, up.get());
    RescaleToShapeOnly(nomSrc, down.get());

    std::unique_ptr<TH1D> nom(static_cast<TH1D*>(nomSrc->Clone()));
    ApplyBinLabels(nom.get(), binLabels);

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
    // Full chi2/rate/mean/variance breakdown is still in the companion CSV
    // (WriteCSVRow); the beam/POT identification that used to be a separate
    // corner label (DrawBeamLabel) is now folded into this single title.
    nom->SetTitle(TString::Format("%s - %s",
      r.name.c_str(), BeamPOTLabel(beam, pot).c_str()));
    nom->GetXaxis()->SetLabelSize(0.0);

    nom->Draw("hist");
    up->Draw("hist same");
    down->Draw("hist same");

    TLegend leg(0.60, 0.63, 0.88, 0.81);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(nom.get(), "nominal",              "l");
    leg.AddEntry(up.get(),  "+1#sigma (rescaled)", "l");
    leg.AddEntry(down.get(),"-1#sigma (rescaled)", "l");
    leg.Draw();

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

  // objName/varLabel/pdfName/csvName for one (category, diagnostic) pair.
  // binLabels is optional -- only set for categorical diagnostics whose
  // bin index isn't self-explanatory (see ApplyBinLabels).
  struct VarSpec
  {
    std::string category, objName, varLabel, pdfName, csvName;
    std::vector<std::string> binLabels;
  };

  // Each priority tier is produced by its own script/job (see
  // make_neutron_spectra_systs_diagnostics_{high,medium,low}.C) into its
  // own ROOT file, so each category's spectra live in a different file.
  std::string RootFileForCategory(const std::string& category)
  {
    if (category == "FSI" || category == "MEC" || category == "RES")
      return "make_neutron_spectra_systs_diagnostics_high.root";
    if (category == "QE" || category == "DIS")
      return "make_neutron_spectra_systs_diagnostics_medium.root";
    return "make_neutron_spectra_systs_diagnostics_low.root"; // LowPriority
  }

  void PlotBeam(TFile* fIn, TCanvas& c, const std::string& pdfName,
                const std::string& beam, std::ofstream& csv, const VarSpec& v)
  {
    TDirectory* topDir = fIn->GetDirectory((beam + "/" + v.category).c_str());
    if (!topDir) {
      std::cerr << "[WARN] no directory for " << beam << "/" << v.category << "\n";
      return;
    }

    double pot = 0.0;
    std::unique_ptr<TH1D> nominal(LoadHist(topDir, "nominal", v.objName, &pot));
    if (!nominal) {
      std::cerr << "[WARN] no nominal/" << v.objName << " under " << beam << "/" << v.category << "\n";
      return;
    }

    const std::vector<RankedSyst> ranked = RankSysts(topDir, nominal.get(), v.objName);
    std::cout << "\n=== " << beam << " / " << v.category << " / " << v.objName
              << " -- ranked " << ranked.size()
              << " systematics, largest shape chi2 first ===\n";
    std::cout << std::left << std::setw(28) << "systematic"
              << std::right << std::setw(14) << "chi2_shape(+)"
              << std::setw(14) << "chi2_shape(-)"
              << std::setw(14) << "rate%(+)"
              << std::setw(14) << "rate%(-)" << "\n";
    for (const RankedSyst& r : ranked) {
      std::cout << std::left << std::setw(28) << r.name
                << std::right << std::setw(14) << r.chi2ShapeUp
                << std::setw(14) << r.chi2ShapeDown
                << std::setw(14) << r.rateChangeUp
                << std::setw(14) << r.rateChangeDown << "\n";
    }

    TDirectory* systsDir = topDir->GetDirectory("systs");
    for (const RankedSyst& r : ranked) {
      DrawShapePage(c, pdfName, beam, pot, r, nominal.get(), systsDir, v.objName, v.varLabel, v.category, v.binLabels);
      WriteCSVRow(csv, beam, r);
    }
  }
}

void plot_neutron_spectra_systs_diagnostics_ranking()
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.032);
  gStyle->SetTitleY(0.97);  // nudge the splitline chi2/rate title down from the pad's top edge
  TH1::AddDirectory(kFALSE);
  // Force scientific notation on event-count y-axes, same as
  // nu_interactions_plot_topology_only.C.
  TGaxis::SetMaxDigits(3);

  // Each tier's production script writes its own ROOT file (see
  // RootFileForCategory) so any subset of {high,medium,low} can be run
  // (and re-run) independently -- open whichever are present; a missing
  // tier just means that tier's categories are skipped below with a
  // warning, not a hard failure.
  const std::vector<std::string> kTierFiles = {
    "make_neutron_spectra_systs_diagnostics_high.root",
    "make_neutron_spectra_systs_diagnostics_medium.root",
    "make_neutron_spectra_systs_diagnostics_low.root",
  };
  std::map<std::string, TFile*> tierFiles;
  for (const std::string& fname : kTierFiles) {
    TFile* f = TFile::Open(fname.c_str(), "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "[WARN] could not open " << fname << " -- skipping its categories\n";
      continue;
    }
    tierFiles[fname] = f;
  }
  if (tierFiles.empty()) {
    std::cerr << "[ERR] none of the tier ROOT files were found -- run at least one of "
                 "make_neutron_spectra_systs_diagnostics_{high,medium,low}.C first\n";
    return;
  }

  TCanvas c("c", "c", 800, 700);

  // Short axis labels for the two categorical diagnostics -- see
  // ../genie/NeutronMultSystDiagnostics.cxx's kFinalStateTopologyBinVar for
  // the exact nP/nN/nPi definition of each bin (mirrors cuts.C's kFSI_*).
  const std::vector<std::string> kTopologyBinLabels = {
    "0p 1n", "0p #geq2n", "0p 1n 1#pi", "1p 0#pi",
    "#geq2p 0#pi", "1p N#pi", "#geq2p N#pi", "other/0n",
  };
  const std::vector<std::string> kProngsPerNeutronBinLabels = {"1", "2", "3+"};

  const std::string kOutPrefix = "neutron_spectra_systs_diagnostics_";
  const std::vector<VarSpec> variables = {
    {"FSI", "neutron_ke_visible",      "Visible neutron KE [MeV]",
     kOutPrefix + "FSI_ke_plots.pdf",          kOutPrefix + "FSI_ke_metrics.csv"},
    {"FSI", "prongs_per_neutron",      "Prongs per visible neutron",
     kOutPrefix + "FSI_prongs_plots.pdf",      kOutPrefix + "FSI_prongs_metrics.csv",
     kProngsPerNeutronBinLabels},
    {"FSI", "n_primary_minus_visible", "N primary - N visible neutrons",
     kOutPrefix + "FSI_confusion_plots.pdf",   kOutPrefix + "FSI_confusion_metrics.csv"},

    {"MEC", "hadvise_q0",           "HadVisE [MeV]",
     kOutPrefix + "MEC_q0_plots.pdf",     kOutPrefix + "MEC_q0_metrics.csv"},
    {"MEC", "ntrue_neutron_mult",   "N true neutron prongs",
     kOutPrefix + "MEC_mult_plots.pdf",   kOutPrefix + "MEC_mult_metrics.csv"},
    {"MEC", "final_state_topology", "Final-state topology bin",
     kOutPrefix + "MEC_topo_plots.pdf",   kOutPrefix + "MEC_topo_metrics.csv",
     kTopologyBinLabels},

    {"RES", "hadvise_q0",           "HadVisE [MeV]",
     kOutPrefix + "RES_q0_plots.pdf",     kOutPrefix + "RES_q0_metrics.csv"},
    {"RES", "ntrue_neutron_mult",   "N true neutron prongs",
     kOutPrefix + "RES_mult_plots.pdf",   kOutPrefix + "RES_mult_metrics.csv"},
    {"RES", "final_state_topology", "Final-state topology bin",
     kOutPrefix + "RES_topo_plots.pdf",   kOutPrefix + "RES_topo_metrics.csv",
     kTopologyBinLabels},

    {"QE", "hadvise_cleanccqe", "HadVisE [MeV] (0p, 0#pi)",
     kOutPrefix + "QE_clean_plots.pdf", kOutPrefix + "QE_clean_metrics.csv"},
    {"QE", "hadvise_inclusive", "HadVisE [MeV] (inclusive)",
     kOutPrefix + "QE_incl_plots.pdf",  kOutPrefix + "QE_incl_metrics.csv"},

    {"DIS", "ntrue_neutron_mult_tail", "N true neutron prongs",
     kOutPrefix + "DIS_mult_plots.pdf", kOutPrefix + "DIS_mult_metrics.csv"},
    {"DIS", "final_state_topology",    "Final-state topology bin",
     kOutPrefix + "DIS_topo_plots.pdf", kOutPrefix + "DIS_topo_metrics.csv",
     kTopologyBinLabels},

    {"LowPriority", "ntrue_neutron_mult_sanity", "N true neutron prongs",
     kOutPrefix + "LowPriority_mult_plots.pdf", kOutPrefix + "LowPriority_mult_metrics.csv"},
  };

  const std::vector<std::string> beams = {"FHC", "RHC"};

  for (const VarSpec& v : variables) {
    const std::string tierFile = RootFileForCategory(v.category);
    auto it = tierFiles.find(tierFile);
    if (it == tierFiles.end()) {
      std::cerr << "[WARN] skipping " << v.category << "/" << v.objName
                << " -- " << tierFile << " not open\n";
      continue;
    }
    TFile* fIn = it->second;

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
    for (const std::string& beam : beams)
      PlotBeam(fIn, c, v.pdfName, beam, csv, v);
    c.Print((v.pdfName + "]").c_str());

    csv.close();
    std::cout << "Saved " << v.category << "/" << v.objName << " plots to " << v.pdfName
              << " and metrics to " << v.csvName << "\n";
  }

  for (auto& kv : tierFiles) kv.second->Close();
}

#endif
