// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_mc_spectra.C

#ifdef __CINT__
void make_mc_spectra()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "../../cuts.C"

#include "TFile.h"
#include "TH2.h"
#include "TParameter.h"
#include <string>
#include <vector>

// cafana weights includes
#include "CAFAna/Weights/GenieWeights.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"
using namespace ana;

namespace
{
  // ── per-q0sample-per-beam helper ─────────────────────────────────────────
  void ProcessBeam(const std::string& defname,
                   const std::string& beam,
                   const ana::Cut&    baseCut,
                   const ana::Cut&    baseCutNoPt,
                   const ana::Cut&    topoCut,
                   const std::string& topoName,
                   TFile*             fOut)
  {
    const ana::Cut fullCut = baseCut && topoCut;
    const ana::Weight Wei = ana::kPPFXFluxCVWgt * ana::kXSecCVWgt2024;

    //Compare the spectrum of energy of the neutrons: both genie neutrons and the visible neutrons. 

    // Kinetic energy [MeV] of a prong's true parent particle (E - m_n).
    ExtractFn fNeutronKE = [](const caf::SRProngProxy* p){
      return (double)(p->truth.motherp.E * 1000.0 - 939.565);
    };
    ExtractFn fNeutronE = [](const caf::SRProngProxy* p){
      return (double)(p->truth.motherp.E * 1000.0);
    };
    // Kinetic energy [MeV] of the induced prong, assuming it's a proton (E - m_p).
    ExtractFn fProtonKE = [](const caf::SRProngProxy* p){ return (double)(p->truth.p.E * 1000.0 - 938.272); };

    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    const ana::Binning kBinNeutronKE    = ana::Binning::Simple(40, 0.0, 100.0);  // MeV, full range
    // Coarser than kBinNeutronKE: the 2D KE-vs-protonKE correlation is
    // sparsest exactly where it's most populated (low KE), so finer bins
    // there just spread MC statistics thinner and produce visual noise.
    const ana::Binning kBinKE = ana::Binning::Simple(20, 0.0, 80.0);
    const ana::Binning kBinE    = ana::Binning::Simple(20, 0.0, 80.0);

    // ── all primary GENIE neutrons (pre-FSI truth list) ────────────────────
    Spectrum sKE_prim_all(loader,
      MultiVarHistAxis("Primary Neutrons KE [MeV]", kBinNeutronKE, MakeGENIENeutronKineticEnergyMV()),
      fullCut, kNoShift, Wei);

    // ── visible primary neutrons (primary neutron ancestor left a prong) ───
    // One entry per unique neutron (deduplicated by parent 4-momentum), not
    // per prong -- a neutron that leaves 2+ prongs must only count once.
    Spectrum sKE_prim_visible(loader,
      MultiVarHistAxis("Visible Neutrons KE [MeV]", kBinNeutronKE, MakeVisibleNeutronKEDedupMV()),
      fullCut, kNoShift, Wei);

    // ── visible neutrons: true KE vs. the matched prong's proton KE ─────────
    Spectrum s2_KE_vs_prongE(loader,
      MultiVarHistAxis("Primary Neutron KE [MeV]", kBinKE, MakeTrueNeutronVar(fNeutronKE),
                       "Proton Kinetic Energy [MeV]",      kBinE,     MakeTrueNeutronVar(fProtonKE)),
      fullCut, kNoShift, Wei);

    // ── visible neutrons split by how many prongs they produced ─────────────
    // Prongs sharing the exact same parent (truth.motherpdg==2112) truth
    // 4-momentum came from the same physical neutron. Each such neutron's
    // true KE goes into exactly one of these three, by prong multiplicity.
    Spectrum sKE_visible_1prong(loader,
      MultiVarHistAxis("Visible Neutron KE [MeV] (1 prong)", kBinNeutronKE,
                       MakeVisibleNeutronKEByProngCountMV(1)),
      fullCut, kNoShift, Wei);
    Spectrum sKE_visible_2prong(loader,
      MultiVarHistAxis("Visible Neutron KE [MeV] (2 prongs)", kBinNeutronKE,
                       MakeVisibleNeutronKEByProngCountMV(2)),
      fullCut, kNoShift, Wei);
    Spectrum sKE_visible_3plusprong(loader,
      MultiVarHistAxis("Visible Neutron KE [MeV] (#geq3 prongs)", kBinNeutronKE,
                       MakeVisibleNeutronKEByProngCountMV(3)),
      fullCut, kNoShift, Wei);

    // ── how many prongs each visible neutron produces (one entry/neutron) ───
    const ana::Binning kBinProngsPerNeutron = ana::Binning::Simple(3, 0.5, 3.5); // 1, 2, 3+
    Spectrum sProngsPerNeutron(loader,
      MultiVarHistAxis("Prongs per visible neutron", kBinProngsPerNeutron,
                       MakeProngsPerVisibleNeutronMV()),
      fullCut, kNoShift, Wei);

    // ── neutron count per event: primary (GENIE) vs. visible ────────────────
    const ana::Binning kBinNeutronCount = ana::Binning::Simple(8, -0.5, 7.5);
    Spectrum sNPrimNeutronsPerEvent(loader,
      HistAxis("N neutrons per event", kBinNeutronCount, MakeNGENIENeutronsVar()),
      fullCut, kNoShift, Wei);
    Spectrum sNVisibleNeutronsPerEvent(loader,
      HistAxis("N neutrons per event", kBinNeutronCount, MakeVisibleNeutronCountVar()),
      fullCut, kNoShift, Wei);

    // ── confusion matrix: N primary neutrons vs. N visible neutrons ─────────
    // One entry per event: (N GENIE primary neutrons, N unique visible
    // neutrons). Restricted to events with >=1 primary neutron, since
    // events with zero primary neutrons (and necessarily zero visible ones)
    // would otherwise dominate the (0,0) cell and dilute the comparison.
    const ana::Cut fullCutHasPrimNeutron = fullCut && MakeHasPrimaryNeutronCut();
    Spectrum s2_confusionNeutronCount(loader,
      HistAxis("N primary neutrons", kBinNeutronCount, MakeNGENIENeutronsVar(),
               "N visible neutrons", kBinNeutronCount, MakeVisibleNeutronCountVar()),
      fullCutHasPrimNeutron, kNoShift, Wei);

    loader.Go();

    TDirectory* beamDir = fOut->GetDirectory(beam.c_str());
    if (!beamDir) beamDir = fOut->mkdir(beam.c_str());
    TDirectory* dir = beamDir->mkdir(topoName.c_str());

    sKE_prim_all.SaveTo(dir, "ke_primary_neutrons");
    sKE_prim_visible.SaveTo(dir, "ke_visible_neutrons");
    sKE_visible_1prong    .SaveTo(dir, "ke_visible_neutrons_1prong");
    sKE_visible_2prong    .SaveTo(dir, "ke_visible_neutrons_2prong");
    sKE_visible_3plusprong.SaveTo(dir, "ke_visible_neutrons_3plusprong");
    sProngsPerNeutron     .SaveTo(dir, "prongs_per_visible_neutron");
    sNPrimNeutronsPerEvent   .SaveTo(dir, "n_primary_neutrons_per_event");
    sNVisibleNeutronsPerEvent.SaveTo(dir, "n_visible_neutrons_per_event");

    // Spectrum::SaveTo() always flattens multi-axis histograms into a 1D
    // "hist" (product binning), so the genuine 2D spectrum must also be
    // exported as a raw TH2 for plotting.
    const double potKEvsE = sKE_prim_visible.POT();
    dir->cd();
    { TParameter<double> potParam("h2_ke_vs_prongE_visible_pot", potKEvsE); potParam.Write(); }
    TH2* h2KEvsProngE = s2_KE_vs_prongE.ToTH2(potKEvsE);
    h2KEvsProngE->SetName("h2_ke_vs_prongE_visible");
    dir->cd();
    h2KEvsProngE->Write();
    delete h2KEvsProngE;

    const double potConfusion = sNPrimNeutronsPerEvent.POT();
    dir->cd();
    { TParameter<double> potParam("h2_confusion_neutron_count_pot", potConfusion); potParam.Write(); }
    TH2* h2Confusion = s2_confusionNeutronCount.ToTH2(potConfusion);
    h2Confusion->SetName("h2_confusion_neutron_count");
    dir->cd();
    h2Confusion->Write();
    delete h2Confusion;

    std::cout << beam << "/" << topoName << " done.\n";
  }
}

void make_mc_spectra()
{
  const ana::Weight Wei = ana::kPPFXFluxCVWgt * ana::kXSecCVWgt2024;
  const ana::Cut baseCutNoPt =
    ana::xsec::numubarcc::kQualityCut             &&
    ana::xsec::numubarcc::kContainmentCut         &&
    ana::xsec::numubarcc::kRecoVtxNumuFiducialCut &&
    ana::xsec::numubarcc::kMuonIDCut;

  const ana::Cut baseCut = baseCutNoPt;

  const std::string defFHC= "prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_filematchedSystematics_nominal_v1";
  const std::string defRHC = "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_ndphysics_systs_filematch_v1";

  //const std::string defname_FHC = "prod_sumdecaf_R20-11-25-prod5.1reco.t_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_menate_v1_numu2020";
  //const std::string defname_RHC = "prod_sumdecaf_R20-11-25-prod5.1reco.t_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_menate_v1_numu2020";


  // Two q0 samples: low (q0 < 0.3 GeV) and high (q0 >= 0.3 GeV),
  // both using the inclusive neutron topology.
  using CutFn = std::function<ana::Cut(const std::string&)>;
  const std::vector<std::pair<std::string, CutFn>> samples = {
    {"sample1_q0Lo", [](const std::string& b){ return kTopoInclusive(b) && MakeQ0Cut(0.0,  0.150);   }},
    {"sample3_full", [](const std::string& b){ return kTopoInclusive(b);                             }},
  };

  TFile* fOut = new TFile("make_mc_spectra.root", "RECREATE");

  for (const auto& sample : samples) {
    const std::string& sampleName = sample.first;
    const CutFn&       cutFn      = sample.second;
    ProcessBeam(defFHC, "FHC", baseCut, baseCutNoPt, cutFn("FHC"), sampleName, fOut);
    ProcessBeam(defRHC, "RHC", baseCut, baseCutNoPt, cutFn("RHC"), sampleName, fOut);
  }

  fOut->Write();
  fOut->Close();
  delete fOut;
  std::cout << "Saved to make_mc_spectra.root\n";
}

#endif
