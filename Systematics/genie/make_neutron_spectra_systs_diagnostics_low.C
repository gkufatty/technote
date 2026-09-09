// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs_diagnostics_low.C
//
// Low-priority tier only (COH, NC elastic, radiative/2nd-class-current
// knobs -- per NeutronMultSystLists.cxx's getNeutronMultSysts_LowPriority()).
// The technote expects these to be negligible for a numu(bar) CC neutron-
// multiplicity measurement (no nucleon knockout, or an out-of-selection
// channel), so this is a single true-multiplicity sanity check, not a
// dedicated study -- cheapest of the three tiers, kept separate so it can
// run (and be re-run) independently of the higher tiers.
//
// Read the output with plot_neutron_spectra_systs_diagnostics_ranking.C.

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
#include "NeutronMultSystLists.cxx"
#include "NeutronMultSystDiagnostics.cxx"

  void ProcessBeam(const std::string& defname, const std::string& beam, TFile* fOut)
  {
    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    const std::vector<const ISyst*> lowSysts = getNeutronMultSysts_LowPriority();

    // ── single true-multiplicity sanity check ──────────────────────────────
    std::unique_ptr<Spectrum> lowNomMult;
    std::map<std::string, ShiftPair> lowShMult;
    FillShiftedSpectra(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      kNeutronMultBaseCut, lowSysts, lowNomMult, lowShMult);

    std::cout << "[" << beam << "] registered " << lowSysts.size()
              << " LowPriority systematics; running loader.Go() ...\n";

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
