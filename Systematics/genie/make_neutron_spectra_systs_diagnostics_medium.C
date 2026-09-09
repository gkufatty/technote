// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs_diagnostics_medium.C
//
// Medium-priority tier only (QE + DIS -- systematics that shift QE/DIS mode
// composition, with only an indirect effect on the neutron sample, per
// NeutronMultSystLists.cxx's getNeutronMultSysts_{QE,DIS}()). Split out on
// its own because DIS alone is ~38 knobs (32 soft-pion + Q0/Q1 + AGKY/
// Bodek-Yang) -- by far the most expensive tier to reweight -- so it
// shouldn't gate the high-priority (FSI/MEC/RES) job's turnaround.
//
// Diagnostics (see the technote's "Diagnostic Variables by Systematic
// Category" table):
//   QE  -- HadVisE restricted to the proton-free/pion-free ("clean CCQE-
//          like") topology, and inclusive HadVisE (whose integral is the
//          overall CC rate) -- CCQE sets the clean-neutron population.
//   DIS -- true neutron multiplicity (fine tail binning, since DIS mainly
//          affects the high-multiplicity tail) and the final-state
//          topology table.
//
// Read the output with plot_neutron_spectra_systs_diagnostics_ranking.C.

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
#include "NeutronMultSystLists.cxx"
#include "NeutronMultSystDiagnostics.cxx"

  void ProcessBeam(const std::string& defname, const std::string& beam, TFile* fOut)
  {
    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    const std::vector<const ISyst*> qeSysts  = getNeutronMultSysts_QE();
    const std::vector<const ISyst*> disSysts = getNeutronMultSysts_DIS();

    const ana::Cut topologyCut  = kNeutronMultBaseCut && MakeHasPrimaryNeutronCut();
    const ana::Cut cleanCCQECut = kNeutronMultBaseCut && (kFSI_1Neutron || kFSI_2PlusNeutrons);

    // ── QE: HadVisE for the clean-CCQE-like topology, and inclusive ────────
    std::unique_ptr<Spectrum> qeNomClean, qeNomIncl;
    std::map<std::string, ShiftPair> qeShClean, qeShIncl;
    FillShiftedSpectra(loader,
      HistAxis("HadVisE [MeV] (0p, 0#pi)", kBinHadVisMeV300, kNumuHadVisMeV),
      cleanCCQECut, qeSysts, qeNomClean, qeShClean);
    FillShiftedSpectra(loader,
      HistAxis("HadVisE [MeV] (inclusive)", kBinHadVisMeV300, kNumuHadVisMeV),
      kNeutronMultBaseCut, qeSysts, qeNomIncl, qeShIncl);

    // ── DIS: true neutron multiplicity tail, final-state topology ──────────
    std::unique_ptr<Spectrum> disNomMult, disNomTopo;
    std::map<std::string, ShiftPair> disShMult, disShTopo;
    FillShiftedSpectra(loader,
      HistAxis("N true neutron prongs per event", kBinTruePrim, MakeTruePrimCountVar()),
      kNeutronMultBaseCut, disSysts, disNomMult, disShMult);
    FillShiftedSpectra(loader,
      HistAxis("Final-state topology bin", kBinFinalStateTopology, kFinalStateTopologyBinVar),
      topologyCut, disSysts, disNomTopo, disShTopo);

    std::cout << "[" << beam << "] registered " << qeSysts.size() << " QE, "
              << disSysts.size() << " DIS systematics; running loader.Go() ...\n";

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
