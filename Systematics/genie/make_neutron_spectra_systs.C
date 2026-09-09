// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_systs.C
//
// Produces GENIE cross-section systematic-shifted (+-1 sigma) spectra for
// make_neutron_spectra.C's ntrue_per_event histogram (MakeTruePrimCountVar,
// true neutron prong count) AND ngenie_per_event (MakeNGENIENeutronsVar, N
// GENIE neutrons per event), for FHC and RHC, topo4_All (inclusive
// topology), with a q0 < 100 MeV cut (kNumuHadVisE(sr) < 100 MeV). Shifts
// come from every systematic in ana::getAllXsecNuTruthSysts_2024()
// (novasoft/CAFAna/Systs/XSecSystLists.cxx). Both variables are saved as
// separate objects in the same "nominal"/"systs/<name>/p1sigma"/"m1sigma"
// directories, so ntrue_per_event stays exactly where every downstream
// plot_neutron_spectra_systs_*.C script and neutron_spectra.py/
// neutron_spectra_extended.C already expect it -- ngenie_per_event is
// purely additive.

#ifdef __CINT__
void make_neutron_spectra_systs()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "../cuts.C"

#include "CAFAna/Systs/XSecSystLists.h"
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
  std::unique_ptr<Spectrum> MakeNTruePerEventSpectrum(SpectrumLoader& loader,
                                                       const ana::Cut& fullCut,
                                                       const ana::SystShifts& shift)
  {
    return std::make_unique<Spectrum>(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      fullCut, shift, kUnweighted);
  }

  std::unique_ptr<Spectrum> MakeNGeniePerEventSpectrum(SpectrumLoader& loader,
                                                        const ana::Cut& fullCut,
                                                        const ana::SystShifts& shift)
  {
    return std::make_unique<Spectrum>(loader,
      HistAxis("N GENIE neutrons per event", kBinNSel, MakeNGENIENeutronsVar()),
      fullCut, shift, kUnweighted);
  }

  void ProcessBeam(const std::string&                  defname,
                    const std::string&                  beam,
                    const ana::Cut&                     baseCut,
                    const std::vector<const ana::ISyst*>& systs,
                    TFile*                               fOut)
  {
    // All topologies in the range HadVisE < 100 MeV (kNumuHadVisE(sr) < 100 MeV)
    const ana::Cut fullCut = baseCut && MakeVisHadEMeVCut(0.0, 100.0);

    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    std::unique_ptr<Spectrum> nominalTrue  = MakeNTruePerEventSpectrum(loader, fullCut, kNoShift);
    std::unique_ptr<Spectrum> nominalGenie = MakeNGeniePerEventSpectrum(loader, fullCut, kNoShift);

    struct ShiftPair { std::unique_ptr<Spectrum> p1sigma, m1sigma; };
    std::map<std::string, ShiftPair> shiftedTrue, shiftedGenie;

    for (const ana::ISyst* syst : systs) {
      const std::string& name = syst->ShortName();

      ShiftPair& spTrue = shiftedTrue[name];
      spTrue.p1sigma = MakeNTruePerEventSpectrum(loader, fullCut, ana::SystShifts(syst, +1.0));
      spTrue.m1sigma = MakeNTruePerEventSpectrum(loader, fullCut, ana::SystShifts(syst, -1.0));

      ShiftPair& spGenie = shiftedGenie[name];
      spGenie.p1sigma = MakeNGeniePerEventSpectrum(loader, fullCut, ana::SystShifts(syst, +1.0));
      spGenie.m1sigma = MakeNGeniePerEventSpectrum(loader, fullCut, ana::SystShifts(syst, -1.0));
    }

    std::cout << "[" << beam << "] " << systs.size()
              << " systematics registered (+-1 sigma each) for ntrue_per_event + "
              << "ngenie_per_event; running loader.Go() ...\n";

    loader.Go();

    TDirectory* beamDir = fOut->GetDirectory(beam.c_str());
    if (!beamDir) beamDir = fOut->mkdir(beam.c_str());
    TDirectory* VisEDir = beamDir->mkdir("HadVisE_100MeV");

    TDirectory* nomDir = VisEDir->mkdir("nominal");
    nominalTrue ->SaveTo(nomDir, "ntrue_per_event");
    nominalGenie->SaveTo(nomDir, "ngenie_per_event");

    TDirectory* systsDir = VisEDir->mkdir("systs");
    for (const ana::ISyst* syst : systs) {
      const std::string& name = syst->ShortName();
      TDirectory* sDir  = systsDir->mkdir(name.c_str());
      TDirectory* p1Dir = sDir->mkdir("p1sigma");
      TDirectory* m1Dir = sDir->mkdir("m1sigma");

      shiftedTrue[name] .p1sigma->SaveTo(p1Dir, "ntrue_per_event");
      shiftedTrue[name] .m1sigma->SaveTo(m1Dir, "ntrue_per_event");
      shiftedGenie[name].p1sigma->SaveTo(p1Dir, "ngenie_per_event");
      shiftedGenie[name].m1sigma->SaveTo(m1Dir, "ngenie_per_event");
    }

    std::cout << beam << "/HadVisE_100MeV done.\n";
  }
}

void make_neutron_spectra_systs()
{
  const ana::Cut baseCut =
    ana::xsec::numubarcc::kQualityCut             &&
    ana::xsec::numubarcc::kContainmentCut         &&
    ana::xsec::numubarcc::kRecoVtxNumuFiducialCut &&
    ana::xsec::numubarcc::kMuonIDCut;

  const std::string defFHC =
    "prod_sumdecaf_R20-11-25-prod5.1reco.t_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_menate_v1_numu2020";
  const std::string defRHC =
    "prod_sumdecaf_R20-11-25-prod5.1reco.t_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_menate_v1_numu2020";

  const std::vector<const ana::ISyst*> systs = ana::getAllXsecNuTruthSysts_2024();
  std::cout << "Loaded " << systs.size()
            << " GENIE xsec systematics from getAllXsecNuTruthSysts_2024()\n";

  TFile* fOut = new TFile("make_neutron_spectra_systs.root", "RECREATE");

  ProcessBeam(defFHC, "FHC", baseCut, systs, fOut);
  ProcessBeam(defRHC, "RHC", baseCut, systs, fOut);

  fOut->Write();
  fOut->Close();
  delete fOut;
  std::cout << "Saved to make_neutron_spectra_systs.root\n";
}

#endif
