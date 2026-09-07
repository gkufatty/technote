// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 nu_interactions_spectra_topology_only.C

#ifdef __CINT__
void nu_interactions_spectra_topology_only()
{
  std::cout << "Sorry, you must run in compiled mode" << std::endl;
}
#else

// ana::kNumuHadVisE and the six kFSI_* grouping cuts live in selection_vars.h
// / cuts.C -- same single include make_mc_spectra.C uses, not the ad hoc,
// unresolvable "NeutronMultiplicityVars.h"/"NeutronMultiplicityCuts.h"
// this script used to reach for.
#include "../../cuts.C"

#include "TFile.h"
#include <string>

// cafana weights includes -- same weight as make_mc_spectra.C
#include "CAFAna/Weights/GenieWeights.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"

using namespace ana;

// Produces the spectra needed by the FSI grouping overlay in
// nu_interactions_plot_topology_only.C, using a different category list
// per beam so each only carries categories that are actually populated
// (per check_final_states.C's proton x pion table and top-final-state
// list):
//   FHC: TrueHadVisE_{1N1pi,1p0pi,2p0pi,1pNpi,2pNpi}_FHC
//     -- kFSI_1Neutron/kFSI_2PlusNeutrons (0p, 0pi) are dropped for FHC:
//        combined they're only ~1.9% of FHC's >=1-neutron sample. FHC is
//        RES/DIS-driven, not CCQE-driven, so a clean proton-free,
//        pion-free final state is rare there.
//   RHC: TrueHadVisE_{1N,g2N,1N1pi,1p0pi,2p0pi,1pNpi}_RHC
//     -- kFSI_2PlusProtons_NPion (>=2p, >=1pi) is dropped for RHC: it's
//        the highest-multiplicity category and needs more hadronic energy
//        than RHC's other categories to produce, so within RHC's 0-100
//        MeV display range (tuned for the CCQE-like 1N/g2N peak) it was
//        mostly flat tail, not a meaningful curve.
// Both beams get kFSI_1Neutron_1Pion (0p, 1n, 1pi) on its own: it's FHC's
// single most common final state (9.49%) and RHC's second most common
// (10.51%), so it's relevant either way, unlike the plain 0p+0pi case.
//
// Plus the inclusive TrueHadVisE_{FHC,RHC} background, all vs. visible
// hadronic energy (ana::kNumuHadVisMeV) rather than true neutrino energy.
//
// NOTE: neither beam's list sums to N(>=1 neutron) -- 0-proton events with
// >=2 pions, or (for FHC) 0-proton+0/1-pion, still fall outside everything
// listed for that beam; see cuts.C's kFSI_* comments for exact shares.
void nu_interactions_spectra_topology_only()
{
  // Same weight as make_mc_spectra.C
  const ana::Weight Wei = ana::kPPFXFluxCVWgt * ana::kXSecCVWgt2024;

  // MakeHasPrimaryNeutronCut() (cuts.C) restores the ">=1 primary neutron"
  // scoping that used to be baked into the old kTrueNeutrinoE MultiVar's
  // own "hasNeutron" gate. ana::kNumuHadVisE has no such gate -- it's a
  // plain Var that fires on every event -- so without this, "Inclusive"
  // would include 0-neutron events, no longer matching what the six
  // kFSI_* categories (all implicitly or explicitly n>=1) are subsets of.
  //
  // kCutPrimMuon (cuts.C, |pdg|==13) additionally requires >=1 true primary
  // muon/antimuon: kMuonIDCut above is a RECO-level muon-ID cut and does
  // NOT guarantee this at truth level -- check_final_states.C found ~2% of
  // events with a true primary neutron have zero true primary muons
  // (e.g. mu0_p0_n1_pi1_other1). Since none of the six kFSI_* cuts check
  // muon count themselves, this ensures every grouping (Inclusive and all
  // six) is a subset of {true muon>0} x {true neutron>0}.
  const ana::Cut eventLevelCut =
      ana::xsec::numubarcc::kQualityCut &&
      ana::xsec::numubarcc::kContainmentCut &&
      ana::xsec::numubarcc::kRecoVtxNumuFiducialCut &&
      ana::xsec::numubarcc::kMuonIDCut &&
      kCutPrimMuon &&
      MakeHasPrimaryNeutronCut();

  // Same samples as make_mc_spectra.C
  const std::string defFHC =
    "prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_filematchedSystematics_nominal_v1";
  const std::string defRHC =
    "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_ndphysics_systs_filematch_v1";

  SpectrumLoader loaderFHC(defFHC);
  SpectrumLoader loaderRHC(defRHC);

  // Same spill cut as make_mc_spectra.C
  loaderFHC.SetSpillCut(kStandardSpillCuts);
  loaderRHC.SetSpillCut(kStandardSpillCuts);

  // Same bin width (10 MeV) for both beams, different upper edge: RHC's
  // FSI activity is concentrated at lower HadVisE than FHC's (see the
  // proton x pion table in check_final_states.C -- RHC is dominated by
  // proton-free, pion-free CCQE-like events), so 100 MeV covers RHC's
  // populated range while FHC uses 350 MeV.
  const Binning bins_HadVisE_FHC = Binning::Simple(35, 0, 350); // MeV
  const Binning bins_HadVisE_RHC = Binning::Simple(10, 0, 100); // MeV

  // ana::kNumuHadVisMeV (selection_vars.h) = kNumuHadVisE(sr)*1000.0, i.e.
  // the same standard CAFAna visible-hadronic-energy Var used elsewhere in
  // this technote (see cuts.C's MakeVisHadEMeVCut), just in MeV. It's a
  // plain Var, not a MultiVar, so this uses HistAxis rather than
  // MultiVarHistAxis.
  HistAxis TrueHadVisEAxis_FHC(
      "Visible Hadronic Energy [MeV]", bins_HadVisE_FHC, kNumuHadVisMeV
  );
  HistAxis TrueHadVisEAxis_RHC(
      "Visible Hadronic Energy [MeV]", bins_HadVisE_RHC, kNumuHadVisMeV
  );

  // =========================================================
  // FHC
  // =========================================================
  Spectrum TrueHadVisE_FHC       (loaderFHC, TrueHadVisEAxis_FHC, eventLevelCut,                            kNoShift, Wei);
  Spectrum TrueHadVisE_1N1pi_FHC (loaderFHC, TrueHadVisEAxis_FHC, eventLevelCut && kFSI_1Neutron_1Pion,     kNoShift, Wei);
  Spectrum TrueHadVisE_1p0pi_FHC (loaderFHC, TrueHadVisEAxis_FHC, eventLevelCut && kFSI_1Proton_0Pion,      kNoShift, Wei);
  Spectrum TrueHadVisE_2p0pi_FHC (loaderFHC, TrueHadVisEAxis_FHC, eventLevelCut && kFSI_2PlusProtons_0Pion, kNoShift, Wei);
  Spectrum TrueHadVisE_1pNpi_FHC (loaderFHC, TrueHadVisEAxis_FHC, eventLevelCut && kFSI_1Proton_NPion,      kNoShift, Wei);
  Spectrum TrueHadVisE_2pNpi_FHC (loaderFHC, TrueHadVisEAxis_FHC, eventLevelCut && kFSI_2PlusProtons_NPion, kNoShift, Wei);

  // =========================================================
  // RHC
  // =========================================================
  // No 2+p_Npi_RHC: that category (>=2 protons AND >=1 pion) is the
  // highest-multiplicity of the six and needs more hadronic energy than
  // the other RHC categories to produce -- within RHC's 0-100 MeV display
  // range (tuned for the CCQE-like 1N/g2N peak) it's mostly just tail, not
  // a meaningful curve. Dropped rather than shown near-zero across the
  // whole range.
  Spectrum TrueHadVisE_RHC       (loaderRHC, TrueHadVisEAxis_RHC, eventLevelCut,                            kNoShift, Wei);
  Spectrum TrueHadVisE_1N_RHC    (loaderRHC, TrueHadVisEAxis_RHC, eventLevelCut && kFSI_1Neutron,           kNoShift, Wei);
  Spectrum TrueHadVisE_g2N_RHC   (loaderRHC, TrueHadVisEAxis_RHC, eventLevelCut && kFSI_2PlusNeutrons,      kNoShift, Wei);
  Spectrum TrueHadVisE_1N1pi_RHC (loaderRHC, TrueHadVisEAxis_RHC, eventLevelCut && kFSI_1Neutron_1Pion,     kNoShift, Wei);
  Spectrum TrueHadVisE_1p0pi_RHC (loaderRHC, TrueHadVisEAxis_RHC, eventLevelCut && kFSI_1Proton_0Pion,      kNoShift, Wei);
  Spectrum TrueHadVisE_2p0pi_RHC (loaderRHC, TrueHadVisEAxis_RHC, eventLevelCut && kFSI_2PlusProtons_0Pion, kNoShift, Wei);
  Spectrum TrueHadVisE_1pNpi_RHC (loaderRHC, TrueHadVisEAxis_RHC, eventLevelCut && kFSI_1Proton_NPion,      kNoShift, Wei);

  loaderFHC.Go();
  loaderRHC.Go();

  TFile* fOut = new TFile("nu_interactions_spectra_topology_only.root", "RECREATE");

  TrueHadVisE_FHC.SaveTo(fOut,       "TrueHadVisE_FHC");
  TrueHadVisE_1N1pi_FHC.SaveTo(fOut, "TrueHadVisE_1N1pi_FHC");
  TrueHadVisE_1p0pi_FHC.SaveTo(fOut, "TrueHadVisE_1p0pi_FHC");
  TrueHadVisE_2p0pi_FHC.SaveTo(fOut, "TrueHadVisE_2p0pi_FHC");
  TrueHadVisE_1pNpi_FHC.SaveTo(fOut, "TrueHadVisE_1pNpi_FHC");
  TrueHadVisE_2pNpi_FHC.SaveTo(fOut, "TrueHadVisE_2pNpi_FHC");

  TrueHadVisE_RHC.SaveTo(fOut,       "TrueHadVisE_RHC");
  TrueHadVisE_1N_RHC.SaveTo(fOut,    "TrueHadVisE_1N_RHC");
  TrueHadVisE_g2N_RHC.SaveTo(fOut,   "TrueHadVisE_g2N_RHC");
  TrueHadVisE_1N1pi_RHC.SaveTo(fOut, "TrueHadVisE_1N1pi_RHC");
  TrueHadVisE_1p0pi_RHC.SaveTo(fOut, "TrueHadVisE_1p0pi_RHC");
  TrueHadVisE_2p0pi_RHC.SaveTo(fOut, "TrueHadVisE_2p0pi_RHC");
  TrueHadVisE_1pNpi_RHC.SaveTo(fOut, "TrueHadVisE_1pNpi_RHC");

  fOut->Close();
}

#endif
