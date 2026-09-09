// Run:
// cafe -bq plot_neutron_spectra_systs_ranking.C
//
// Reads make_neutron_spectra_systs.root (produced by
// make_neutron_spectra_systs.C) and, for each beam, ranks every GENIE
// systematic in getAllXsecNuTruthSysts_2024() by its impact on
// ntrue_per_event: the largest |% difference from nominal| in any single
// bin, evaluated at the +1 sigma and -1 sigma shifts. Draws a two-sided
// bar chart (impact bars above/below zero for +1 sigma / -1 sigma),
// sorted largest-impact-first, one page per beam, saved into a PDF.

#ifdef __CINT__
void plot_neutron_spectra_systs_ranking()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "CAFAna/Core/Spectrum.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
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

  TH1D* LoadNTrueHist(TDirectory* topDir, const std::string& subpath)
  {
    TDirectory* parentDir = topDir->GetDirectory(subpath.c_str());
    if (!parentDir) return nullptr;
    std::unique_ptr<Spectrum> spec = Spectrum::LoadFrom(parentDir, "ntrue_per_event");
    return spec->ToTH1(spec->POT());
  }

  struct Impact
  {
    std::string name;
    double up   = 0.0;  // % diff at the most-affected bin, +1 sigma shift
    double down = 0.0;  // % diff at that same bin,          -1 sigma shift
  };

  std::vector<Impact> RankSysts(TDirectory* topDir)
  {
    std::vector<Impact> out;
    std::unique_ptr<TH1D> nominal(LoadNTrueHist(topDir, "nominal"));
    if (!nominal) return out;

    TDirectory* systsDir = topDir->GetDirectory("systs");
    for (const std::string& name : ListSubdirNames(systsDir)) {
      std::unique_ptr<TH1D> up(  LoadNTrueHist(systsDir, name + "/p1sigma"));
      std::unique_ptr<TH1D> down(LoadNTrueHist(systsDir, name + "/m1sigma"));
      if (!up || !down) continue;

      double bestMag = -1.0, bestUp = 0.0, bestDown = 0.0;
      for (int b = 1; b <= nominal->GetNbinsX(); ++b) {
        const double n = nominal->GetBinContent(b);
        if (n == 0.0) continue;
        const double du = 100.0 * (up->GetBinContent(b)   - n) / n;
        const double dd = 100.0 * (down->GetBinContent(b) - n) / n;
        const double mag = std::max(std::abs(du), std::abs(dd));
        if (mag > bestMag) { bestMag = mag; bestUp = du; bestDown = dd; }
      }
      if (bestMag < 0.0) continue;
      out.push_back({name, bestUp, bestDown});
    }

    std::sort(out.begin(), out.end(), [](const Impact& a, const Impact& b) {
      return std::max(std::abs(a.up), std::abs(a.down)) > std::max(std::abs(b.up), std::abs(b.down));
    });
    return out;
  }

  void DrawRanking(TCanvas& c, const std::string& pdfName,
                    const std::string& beam, const std::vector<Impact>& impacts)
  {
    const int n = (int)impacts.size();
    if (n == 0) return;

    c.Clear();
    c.cd();
    c.SetBottomMargin(0.30);
    c.SetTopMargin(0.10);
    c.SetGrid();

    TH1D hUp("hUp", "", n, 0, n);
    TH1D hDown("hDown", "", n, 0, n);
    double axMax = 1.0;
    for (int i = 0; i < n; ++i) {
      hUp.SetBinContent(i + 1, impacts[i].up);
      hDown.SetBinContent(i + 1, impacts[i].down);
      hUp.GetXaxis()->SetBinLabel(i + 1, impacts[i].name.c_str());
      axMax = std::max({axMax, std::abs(impacts[i].up), std::abs(impacts[i].down)});
    }

    hUp.SetTitle(TString::Format("%s: ntrue_per_event GENIE systematic ranking "
                                  "(max |%% diff from nominal| in any bin)", beam.c_str()));
    hUp.GetYaxis()->SetTitle("% difference from nominal");
    hUp.GetXaxis()->LabelsOption("v");
    hUp.GetXaxis()->SetLabelSize(0.028);
    hUp.SetFillColor(kRed);
    hUp.SetFillStyle(1001);
    hUp.SetLineColor(kRed);
    hUp.SetBarWidth(0.8);
    hUp.SetBarOffset(0.1);
    hUp.SetMaximum(axMax * 1.15);
    hUp.SetMinimum(-axMax * 1.15);

    hDown.SetFillColor(kBlue);
    hDown.SetFillStyle(1001);
    hDown.SetLineColor(kBlue);
    hDown.SetBarWidth(0.8);
    hDown.SetBarOffset(0.1);

    hUp.Draw("B");
    hDown.Draw("B same");

    TLine zero(0, 0, n, 0);
    zero.SetLineColor(kBlack);
    zero.Draw();

    c.Print(pdfName.c_str());
  }
}

void plot_neutron_spectra_systs_ranking()
{
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  TFile* fIn = TFile::Open("make_neutron_spectra_systs.root", "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERR] could not open make_neutron_spectra_systs.root\n";
    return;
  }

  // Rank both beams first so the canvas can be sized once, up front, for
  // the largest of the two -- resizing an already-open multi-page PDF
  // mid-stream is an avoidable risk.
  const std::vector<std::string> beams = {"FHC", "RHC"};
  std::vector<std::vector<Impact>> allImpacts;
  int maxN = 0;
  for (const std::string& beam : beams) {
    TDirectory* topDir = fIn->GetDirectory((beam + "/topo4_All_q0lt100MeV").c_str());
    std::vector<Impact> impacts = topDir ? RankSysts(topDir) : std::vector<Impact>();
    std::cout << "[" << beam << "] ranked " << impacts.size() << " systematics\n";
    maxN = std::max(maxN, (int)impacts.size());
    allImpacts.push_back(std::move(impacts));
  }

  const std::string pdfName = "neutron_spectra_systs_ranking.pdf";
  TCanvas c("c", "c", std::max(1200, 22 * maxN), 700);

  c.Print((pdfName + "[").c_str());
  for (size_t i = 0; i < beams.size(); ++i)
    DrawRanking(c, pdfName, beams[i], allImpacts[i]);
  c.Print((pdfName + "]").c_str());

  fIn->Close();
  std::cout << "Saved ranking plot to " << pdfName << "\n";
}

#endif
