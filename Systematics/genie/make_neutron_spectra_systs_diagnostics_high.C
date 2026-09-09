// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs_diagnostics_high.C
//
// High-priority tier only (FSI + MEC + RES -- the systematics the technote
// says directly set the neutron/proton split or the FSI cascade, per
// NeutronMultSystLists.cxx's getNeutronMultSysts_{FSI,MEC,RES}()). Split
// out from the medium/low tiers into its own job so it isn't held up by
// the ~38-knob DIS/soft-pion medium tier: this tier is the one the
// technote calls "highest priority", so it's the one worth iterating on
// fastest.
//
// For each category, produces +-1 sigma shifted spectra (FHC and RHC) of
// the diagnostic(s) the technote's "Diagnostic Variables by Systematic
// Category" table recommends:
//   FSI -- visible-neutron KE spectrum, prongs-per-visible-neutron, and
//          N(primary) - N(visible) neutrons (confusion matrix collapsed to
//          one axis). FSI acts on propagation/absorption, not production,
//          so the true multiplicity histogram alone would likely miss it.
//   MEC -- HadVisE (q0), true neutron multiplicity, final-state topology.
//   RES -- same three diagnostics as MEC (RES reshapes q0 and the
//          target-nucleon/topology split the same way MEC's knobs do).
//
// Read the output with plot_neutron_spectra_systs_diagnostics_ranking.C,
// which computes per-systematic chi2/rate/KS metrics for every
// (category, diagnostic) pair -- the actual "is this relevant?" answer.

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
#include "NeutronMultSystLists.cxx"
#include "NeutronMultSystDiagnostics.cxx"

  void ProcessBeam(const std::string& defname, const std::string& beam, TFile* fOut)
  {
    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    const std::vector<const ISyst*> fsiSysts = getNeutronMultSysts_FSI();
    const std::vector<const ISyst*> mecSysts = getNeutronMultSysts_MEC();
    const std::vector<const ISyst*> resSysts = getNeutronMultSysts_RES();

    const ana::Cut topologyCut = kNeutronMultBaseCut && MakeHasPrimaryNeutronCut();

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

    // ── MEC: q0, true neutron multiplicity, final-state topology ───────────
    std::unique_ptr<Spectrum> mecNomQ0, mecNomMult, mecNomTopo;
    std::map<std::string, ShiftPair> mecShQ0, mecShMult, mecShTopo;
    FillShiftedSpectra(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      kNeutronMultBaseCut, mecSysts, mecNomQ0, mecShQ0);
    FillShiftedSpectra(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      kNeutronMultBaseCut, mecSysts, mecNomMult, mecShMult);
    FillShiftedSpectra(loader,
      HistAxis("Final-state topology bin", kBinFinalStateTopology, kFinalStateTopologyBinVar),
      topologyCut, mecSysts, mecNomTopo, mecShTopo);

    // ── RES: same three diagnostics as MEC, RES systs instead ──────────────
    std::unique_ptr<Spectrum> resNomQ0, resNomMult, resNomTopo;
    std::map<std::string, ShiftPair> resShQ0, resShMult, resShTopo;
    FillShiftedSpectra(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      kNeutronMultBaseCut, resSysts, resNomQ0, resShQ0);
    FillShiftedSpectra(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      kNeutronMultBaseCut, resSysts, resNomMult, resShMult);
    FillShiftedSpectra(loader,
      HistAxis("Final-state topology bin", kBinFinalStateTopology, kFinalStateTopologyBinVar),
      topologyCut, resSysts, resNomTopo, resShTopo);

    std::cout << "[" << beam << "] registered " << fsiSysts.size() << " FSI, "
              << mecSysts.size() << " MEC, " << resSysts.size()
              << " RES systematics; running loader.Go() ...\n";

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
