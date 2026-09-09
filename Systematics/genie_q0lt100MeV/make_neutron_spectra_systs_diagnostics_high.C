// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs_diagnostics_high.C
//
// q0lt100MeV variant of ../genie/make_neutron_spectra_systs_diagnostics_high.C:
// same high-priority tier (FSI + MEC + RES) and same diagnostics, but every
// one of them is now restricted to HadVisE < 100 MeV (see
// Q0Lt100MeVCommon.cxx) instead of the fully inclusive sample -- the
// low-hadronic-energy region the original spectrums/
// make_neutron_spectra_systs.C studied for ntrue_per_event alone. Answers:
// do the same systematics that matter inclusively still matter once you
// restrict to this region, and does their relative ranking change?
//
// Diagnostics (see ../genie/'s header comment for the full rationale):
//   FSI -- visible-neutron KE spectrum, prongs-per-visible-neutron, and
//          N(primary) - N(visible) neutrons, all within HadVisE < 100 MeV.
//   MEC -- HadVisE (q0, now binned 0-100 MeV), true neutron multiplicity,
//          final-state topology, all within HadVisE < 100 MeV.
//   RES -- same three diagnostics as MEC, RES systs instead.
//
// Read the output with plot_neutron_spectra_systs_diagnostics_ranking.C
// (this folder's copy, unchanged from ../genie/'s -- it's cut-agnostic).

#ifdef __CINT__
void make_neutron_spectra_systs_diagnostics_high()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "../../../cuts.C"

// NeutronMultSystLists.cxx (included below) defines all five category
// getters unconditionally, regardless of which ones this tier actually
// calls, so all six syst-family headers are needed here too -- not just
// FSI/MEC/RES's.
#include "CAFAna/Systs/FSISysts.h"
#include "CAFAna/Systs/MECSysts.h"
#include "CAFAna/Systs/RESSysts.h"
#include "CAFAna/Systs/RPASysts.h"
#include "CAFAna/Systs/DISSysts.h"
#include "CAFAna/Systs/XSecSysts.h"  // GetGenieKnobSyst(), rwgt::
#include "CAFAna/Core/SystShifts.h"

#include "TFile.h"
#include "TDirectory.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

// cafana weights includes
#include "CAFAna/Weights/GenieWeights.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"
using namespace ana;

namespace
{
#include "../genie/NeutronMultSystLists.cxx"
#include "../genie/NeutronMultSystDiagnostics.cxx"
#include "Q0Lt100MeVCommon.cxx"

  void ProcessBeam(const std::string& defname, const std::string& beam, TFile* fOut)
  {
    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    const std::vector<const ISyst*> fsiSysts = getNeutronMultSysts_FSI();
    const std::vector<const ISyst*> mecSysts = getNeutronMultSysts_MEC();
    const std::vector<const ISyst*> resSysts = getNeutronMultSysts_RES();

    const ana::Cut topologyCut = kQ0Lt100MeVCut && MakeHasPrimaryNeutronCut();

    // ── FSI: KE spectrum, prong multiplicity, confusion (collapsed) ────────
    std::unique_ptr<Spectrum> fsiNomKE, fsiNomProngs, fsiNomNMinusV;
    std::map<std::string, ShiftPair> fsiShKE, fsiShProngs, fsiShNMinusV;
    FillShiftedSpectra(loader,
      MultiVarHistAxis("Visible neutron KE [MeV]", kBinNeutronKE, MakeVisibleNeutronKEDedupMV()),
      topologyCut, fsiSysts, fsiNomKE, fsiShKE);
    FillShiftedSpectra(loader,
      MultiVarHistAxis("Prongs per visible neutron", kBinProngsPerNeutron, MakeProngsPerVisibleNeutronMV()),
      topologyCut, fsiSysts, fsiNomProngs, fsiShProngs);
    FillShiftedSpectra(loader,
      HistAxis("N primary - N visible neutrons", kBinNPrimMinusVis, kNPrimMinusVisibleVar),
      topologyCut, fsiSysts, fsiNomNMinusV, fsiShNMinusV);

    // ── MEC: q0 (0-100 MeV), true neutron multiplicity, final-state topology ─
    std::unique_ptr<Spectrum> mecNomQ0, mecNomMult, mecNomTopo;
    std::map<std::string, ShiftPair> mecShQ0, mecShMult, mecShTopo;
    FillShiftedSpectra(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV100, kNumuHadVisMeV),
      kQ0Lt100MeVCut, mecSysts, mecNomQ0, mecShQ0);
    FillShiftedSpectra(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      kQ0Lt100MeVCut, mecSysts, mecNomMult, mecShMult);
    FillShiftedSpectra(loader,
      HistAxis("Final-state topology bin", kBinFinalStateTopology, kFinalStateTopologyBinVar),
      topologyCut, mecSysts, mecNomTopo, mecShTopo);

    // ── RES: same three diagnostics as MEC, RES systs instead ──────────────
    std::unique_ptr<Spectrum> resNomQ0, resNomMult, resNomTopo;
    std::map<std::string, ShiftPair> resShQ0, resShMult, resShTopo;
    FillShiftedSpectra(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV100, kNumuHadVisMeV),
      kQ0Lt100MeVCut, resSysts, resNomQ0, resShQ0);
    FillShiftedSpectra(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      kQ0Lt100MeVCut, resSysts, resNomMult, resShMult);
    FillShiftedSpectra(loader,
      HistAxis("Final-state topology bin", kBinFinalStateTopology, kFinalStateTopologyBinVar),
      topologyCut, resSysts, resNomTopo, resShTopo);

    std::cout << "[" << beam << "] registered " << fsiSysts.size() << " FSI, "
              << mecSysts.size() << " MEC, " << resSysts.size()
              << " RES systematics (HadVisE < 100 MeV); running loader.Go() ...\n";

    loader.Go();

    TDirectory* beamDir = GetOrMkdir(fOut, beam);

    TDirectory* fsiDir = GetOrMkdir(beamDir, "FSI");
    SaveDiagnostic(fsiDir, "neutron_ke_visible",      fsiNomKE,      fsiShKE,      fsiSysts);
    SaveDiagnostic(fsiDir, "prongs_per_neutron",      fsiNomProngs,  fsiShProngs,  fsiSysts);
    SaveDiagnostic(fsiDir, "n_primary_minus_visible", fsiNomNMinusV, fsiShNMinusV, fsiSysts);

    TDirectory* mecDir = GetOrMkdir(beamDir, "MEC");
    SaveDiagnostic(mecDir, "hadvise_q0",           mecNomQ0,   mecShQ0,   mecSysts);
    SaveDiagnostic(mecDir, "ntrue_neutron_mult",   mecNomMult, mecShMult, mecSysts);
    SaveDiagnostic(mecDir, "final_state_topology", mecNomTopo, mecShTopo, mecSysts);

    TDirectory* resDir = GetOrMkdir(beamDir, "RES");
    SaveDiagnostic(resDir, "hadvise_q0",           resNomQ0,   resShQ0,   resSysts);
    SaveDiagnostic(resDir, "ntrue_neutron_mult",   resNomMult, resShMult, resSysts);
    SaveDiagnostic(resDir, "final_state_topology", resNomTopo, resShTopo, resSysts);

    std::cout << beam << " done.\n";
  }
}

void make_neutron_spectra_systs_diagnostics_high()
{
  TFile* fOut = new TFile("make_neutron_spectra_systs_diagnostics_high.root", "RECREATE");

  ProcessBeam(kNeutronMultDefFHC, "FHC", fOut);
  ProcessBeam(kNeutronMultDefRHC, "RHC", fOut);

  fOut->Write();
  fOut->Close();
  delete fOut;
  std::cout << "Saved to make_neutron_spectra_systs_diagnostics_high.root\n";
}

#endif
