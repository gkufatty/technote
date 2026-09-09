// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs_diagnostics_medium.C
//
// q0lt100MeV variant of ../genie/make_neutron_spectra_systs_diagnostics_medium.C:
// same medium-priority tier (QE + DIS) and same diagnostics, but every one
// of them is now restricted to HadVisE < 100 MeV (see Q0Lt100MeVCommon.cxx)
// instead of the fully inclusive sample.
//
// Diagnostics (see ../genie/'s header comment for the full rationale):
//   QE  -- HadVisE restricted to the proton-free/pion-free ("clean CCQE-
//          like") topology, and inclusive HadVisE (now binned 0-100 MeV),
//          both additionally restricted to HadVisE < 100 MeV.
//   DIS -- true neutron multiplicity (fine tail binning) and final-state
//          topology, both restricted to HadVisE < 100 MeV.
//
// Read the output with plot_neutron_spectra_systs_diagnostics_ranking.C
// (this folder's copy, unchanged from ../genie/'s -- it's cut-agnostic).

#ifdef __CINT__
void make_neutron_spectra_systs_diagnostics_medium()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "../../../cuts.C"

// NeutronMultSystLists.cxx (included below) defines all five category
// getters unconditionally, regardless of which ones this tier actually
// calls, so all six syst-family headers are needed here too -- not just
// QE/DIS's.
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

    const std::vector<const ISyst*> qeSysts  = getNeutronMultSysts_QE();
    const std::vector<const ISyst*> disSysts = getNeutronMultSysts_DIS();

    const ana::Cut topologyCut  = kQ0Lt100MeVCut && MakeHasPrimaryNeutronCut();
    const ana::Cut cleanCCQECut = kQ0Lt100MeVCut && (kFSI_1Neutron || kFSI_2PlusNeutrons);

    // ── QE: HadVisE for the clean-CCQE-like topology, and inclusive ────────
    std::unique_ptr<Spectrum> qeNomClean, qeNomIncl;
    std::map<std::string, ShiftPair> qeShClean, qeShIncl;
    FillShiftedSpectra(loader,
      HistAxis("HadVisE [MeV] (0p, 0#pi)", kBinHadVisMeV100, kNumuHadVisMeV),
      cleanCCQECut, qeSysts, qeNomClean, qeShClean);
    FillShiftedSpectra(loader,
      HistAxis("HadVisE [MeV] (inclusive)", kBinHadVisMeV100, kNumuHadVisMeV),
      kQ0Lt100MeVCut, qeSysts, qeNomIncl, qeShIncl);

    // ── DIS: true neutron multiplicity tail, final-state topology ──────────
    std::unique_ptr<Spectrum> disNomMult, disNomTopo;
    std::map<std::string, ShiftPair> disShMult, disShTopo;
    FillShiftedSpectra(loader,
      HistAxis("N true neutron prongs per event", kBinTruePrim, MakeTruePrimCountVar()),
      kQ0Lt100MeVCut, disSysts, disNomMult, disShMult);
    FillShiftedSpectra(loader,
      HistAxis("Final-state topology bin", kBinFinalStateTopology, kFinalStateTopologyBinVar),
      topologyCut, disSysts, disNomTopo, disShTopo);

    std::cout << "[" << beam << "] registered " << qeSysts.size() << " QE, "
              << disSysts.size() << " DIS systematics (HadVisE < 100 MeV); "
              << "running loader.Go() ...\n";

    loader.Go();

    TDirectory* beamDir = GetOrMkdir(fOut, beam);

    TDirectory* qeDir = GetOrMkdir(beamDir, "QE");
    SaveDiagnostic(qeDir, "hadvise_cleanccqe", qeNomClean, qeShClean, qeSysts);
    SaveDiagnostic(qeDir, "hadvise_inclusive", qeNomIncl,  qeShIncl,  qeSysts);

    TDirectory* disDir = GetOrMkdir(beamDir, "DIS");
    SaveDiagnostic(disDir, "ntrue_neutron_mult_tail", disNomMult, disShMult, disSysts);
    SaveDiagnostic(disDir, "final_state_topology",    disNomTopo, disShTopo, disSysts);

    std::cout << beam << " done.\n";
  }
}

void make_neutron_spectra_systs_diagnostics_medium()
{
  TFile* fOut = new TFile("make_neutron_spectra_systs_diagnostics_medium.root", "RECREATE");

  ProcessBeam(kNeutronMultDefFHC, "FHC", fOut);
  ProcessBeam(kNeutronMultDefRHC, "RHC", fOut);

  fOut->Write();
  fOut->Close();
  delete fOut;
  std::cout << "Saved to make_neutron_spectra_systs_diagnostics_medium.root\n";
}

#endif
