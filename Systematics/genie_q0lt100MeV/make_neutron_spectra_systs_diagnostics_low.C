// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs_diagnostics_low.C
//
// q0lt100MeV variant of ../genie/make_neutron_spectra_systs_diagnostics_low.C:
// same LowPriority tier (COH, NC elastic, radiative/2nd-class-current
// knobs) and same single true-multiplicity sanity check, but restricted to
// HadVisE < 100 MeV (see Q0Lt100MeVCommon.cxx) instead of the fully
// inclusive sample.
//
// Read the output with plot_neutron_spectra_systs_diagnostics_ranking.C
// (this folder's copy, unchanged from ../genie/'s -- it's cut-agnostic).

#ifdef __CINT__
void make_neutron_spectra_systs_diagnostics_low()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "../../../cuts.C"

// NeutronMultSystLists.cxx (included below) defines all five category
// getters unconditionally, regardless of which ones this tier actually
// calls, so all six syst-family headers are needed here too -- not just
// LowPriority's.
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

    const std::vector<const ISyst*> lowSysts = getNeutronMultSysts_LowPriority();

    // ── single true-multiplicity sanity check, HadVisE < 100 MeV ───────────
    std::unique_ptr<Spectrum> lowNomMult;
    std::map<std::string, ShiftPair> lowShMult;
    FillShiftedSpectra(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      kQ0Lt100MeVCut, lowSysts, lowNomMult, lowShMult);

    std::cout << "[" << beam << "] registered " << lowSysts.size()
              << " LowPriority systematics (HadVisE < 100 MeV); running loader.Go() ...\n";

    loader.Go();

    TDirectory* beamDir = GetOrMkdir(fOut, beam);
    TDirectory* lowDir  = GetOrMkdir(beamDir, "LowPriority");
    SaveDiagnostic(lowDir, "ntrue_neutron_mult_sanity", lowNomMult, lowShMult, lowSysts);

    std::cout << beam << " done.\n";
  }
}

void make_neutron_spectra_systs_diagnostics_low()
{
  TFile* fOut = new TFile("make_neutron_spectra_systs_diagnostics_low.root", "RECREATE");

  ProcessBeam(kNeutronMultDefFHC, "FHC", fOut);
  ProcessBeam(kNeutronMultDefRHC, "RHC", fOut);

  fOut->Write();
  fOut->Close();
  delete fOut;
  std::cout << "Saved to make_neutron_spectra_systs_diagnostics_low.root\n";
}

#endif
