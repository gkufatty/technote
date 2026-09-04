// Run:
// cafe -bq --numubarccinc --neutron_multiplicity -s 50 make_neutron_spectra_q0.C

#ifdef __CINT__
void make_mc_spectra()
{
  std::cout << "Run in compiled mode (ACLiC)." << std::endl;
}
#else

#include "../cuts.C"

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
    const ana::Cut fullCut      = baseCut     && topoCut;
    const ana::Cut fullCutNoPt  = baseCutNoPt && topoCut;

    // variable extractors
    ExtractFn fCalE  = [](const caf::SRProngProxy* p){ return (double)p->calE       * 1000.0; };  // GeV→MeV
    ExtractFn fLen   = [](const caf::SRProngProxy* p){ return (double)p->len; };
    ExtractFn fVisE  = [](const caf::SRProngProxy* p){ return (double)p->truth.visE * 1000.0; };  // GeV→MeV
    ExtractFn fMom   = [](const caf::SRProngProxy* p){ return TrueMom(p); };
    ExtractFn fNhits = [](const caf::SRProngProxy* p){ return (double)p->nhit; };

    SpectrumLoader loader(defname);
    loader.SetSpillCut(kStandardSpillCuts);

    // ── selected (all NeutronicLikeStrict prongs) ──────────────────────────
    Spectrum sCalE_sel(loader,
      MultiVarHistAxis("calE [MeV]",       kBinCalE, MakeNLSVar(beam, kCatSelected, fCalE)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sLen_sel(loader,
      MultiVarHistAxis("len [cm]",         kBinLen,  MakeNLSVar(beam, kCatSelected, fLen)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sVisE_sel(loader,
      MultiVarHistAxis("true visE [MeV]",  kBinVisE, MakeNLSVar(beam, kCatSelected, fVisE)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sMom_sel(loader,
      MultiVarHistAxis("true |p| [GeV/c]", kBinMom,  MakeNLSVar(beam, kCatSelected, fMom)),
      fullCut, kNoShift, kUnweighted);

    // ── signal (NeutronicLikeStrict + motherpdg == neutron) ────────────────
    Spectrum sCalE_sig(loader,
      MultiVarHistAxis("calE [MeV]",       kBinCalE, MakeNLSVar(beam, kCatSignal, fCalE)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sLen_sig(loader,
      MultiVarHistAxis("len [cm]",         kBinLen,  MakeNLSVar(beam, kCatSignal, fLen)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sVisE_sig(loader,
      MultiVarHistAxis("true visE [MeV]",  kBinVisE, MakeNLSVar(beam, kCatSignal, fVisE)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sMom_sig(loader,
      MultiVarHistAxis("true |p| [GeV/c]", kBinMom,  MakeNLSVar(beam, kCatSignal, fMom)),
      fullCut, kNoShift, kUnweighted);

    // ── background (NeutronicLikeStrict + motherpdg != neutron) ───────────
    Spectrum sCalE_bkg(loader,
      MultiVarHistAxis("calE [MeV]",       kBinCalE, MakeNLSVar(beam, kCatBackground, fCalE)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sLen_bkg(loader,
      MultiVarHistAxis("len [cm]",         kBinLen,  MakeNLSVar(beam, kCatBackground, fLen)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sVisE_bkg(loader,
      MultiVarHistAxis("true visE [MeV]",  kBinVisE, MakeNLSVar(beam, kCatBackground, fVisE)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sMom_bkg(loader,
      MultiVarHistAxis("true |p| [GeV/c]", kBinMom,  MakeNLSVar(beam, kCatBackground, fMom)),
      fullCut, kNoShift, kUnweighted);

    // ── event-level NLS count ─────────────────────────────────────────────────
    Spectrum sNSelPerEvent(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNTruePerEvent(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNGENIEPerEvent(loader,
      HistAxis("N GENIE neutrons per event", kBinNSel, MakeNGENIENeutronsVar()),
      fullCut, kNoShift, kUnweighted);

    // ── event-level |q| and Q² distributions ────────────────────────────────
    const ana::Cut kCutQ0Lo = MakeQ0Cut(0.0,  0.10);
    const ana::Cut kCutQ0Hi = MakeQ0Cut(0.10, 9999.0);

    Spectrum sQmag_sel_dist(loader,
      HistAxis("|q| [GeV/c]",      kBinQmag, kRecoQmag),
      fullCut,                                   kNoShift, kUnweighted);
    Spectrum sQmag_true(loader,
      HistAxis("|q| [GeV/c]",      kBinQmag, kRecoQmag),
      fullCut && MakeTrueNeutronPresentCut(),    kNoShift, kUnweighted);
    Spectrum sQ2_sel(loader,
      HistAxis("Q^{2} [GeV^{2}]",  kBinQ2,   kRecoQ2),
      fullCut,                                   kNoShift, kUnweighted);
    Spectrum sQ2_true(loader,
      HistAxis("Q^{2} [GeV^{2}]",  kBinQ2,   kRecoQ2),
      fullCut && MakeTrueNeutronPresentCut(),    kNoShift, kUnweighted);

    // ── q0 (hadronic energy) 1D ───────────────────────────────────────────────
    Spectrum sQ0(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut, kNoShift, kUnweighted);

    // ── HadVisE split by true primary particle species ───────────────────────
    Spectrum sQ0_prim_muon(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutPrimMuon,    kNoShift, kUnweighted);
    Spectrum sQ0_prim_proton(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutPrimProton,  kNoShift, kUnweighted);
    Spectrum sQ0_prim_neutron(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutPrimNeutron, kNoShift, kUnweighted);
    Spectrum sQ0_prim_pion(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutPrimPion,    kNoShift, kUnweighted);
    Spectrum sQ0_prim_other(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutPrimOther,   kNoShift, kUnweighted);

    // ── HadVisE split by prong truth.motherpdg species (all prongs) ─────────
    Spectrum sHadVisE_prongmother_muon(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutProngMotherMuon,    kNoShift, kUnweighted);
    Spectrum sHadVisE_prongmother_proton(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutProngMotherProton,  kNoShift, kUnweighted);
    Spectrum sHadVisE_prongmother_neutron(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutProngMotherNeutron, kNoShift, kUnweighted);
    Spectrum sHadVisE_prongmother_pion(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutProngMotherPion,    kNoShift, kUnweighted);
    Spectrum sHadVisE_prongmother_other(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutProngMotherOther,   kNoShift, kUnweighted);

    // ── HadVisE split by NeutronicLike prong truth.motherpdg species ─────────
    const ana::Cut kCutNLMotherMuon    = MakeNLProngMotherPDGCut(beam, {13});
    const ana::Cut kCutNLMotherProton  = MakeNLProngMotherPDGCut(beam, {2212});
    const ana::Cut kCutNLMotherNeutron = MakeNLProngMotherPDGCut(beam, {2112});
    const ana::Cut kCutNLMotherPion    = MakeNLProngMotherPDGCut(beam, {211, 111});
    const ana::Cut kCutNLMotherOther   = MakeNLProngMotherOtherCut(beam);

    Spectrum sHadVisE_nlmother_muon(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutNLMotherMuon,    kNoShift, kUnweighted);
    Spectrum sHadVisE_nlmother_proton(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutNLMotherProton,  kNoShift, kUnweighted);
    Spectrum sHadVisE_nlmother_neutron(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutNLMotherNeutron, kNoShift, kUnweighted);
    Spectrum sHadVisE_nlmother_pion(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutNLMotherPion,    kNoShift, kUnweighted);
    Spectrum sHadVisE_nlmother_other(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV300, kNumuHadVisMeV),
      fullCut && kCutNLMotherOther,   kNoShift, kUnweighted);

    // ── NL prong motherpdg composition vs HadVisE (per-prong fill) ───────────
    // Each NeutronicLike prong fills the event's HadVisE once.
    // Divide species by total per bin in the plot to get composition fractions.
    Spectrum sNLcomp_total(loader,
      MultiVarHistAxis("HadVisE [MeV]", kBinHadVisMeV300, MakeNLProngTotalHadVisMV(beam)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNLcomp_muon(loader,
      MultiVarHistAxis("HadVisE [MeV]", kBinHadVisMeV300, MakeNLProngMotherHadVisMV(beam, {13})),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNLcomp_proton(loader,
      MultiVarHistAxis("HadVisE [MeV]", kBinHadVisMeV300, MakeNLProngMotherHadVisMV(beam, {2212})),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNLcomp_neutron(loader,
      MultiVarHistAxis("HadVisE [MeV]", kBinHadVisMeV300, MakeNLProngMotherHadVisMV(beam, {2112})),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNLcomp_pion(loader,
      MultiVarHistAxis("HadVisE [MeV]", kBinHadVisMeV300, MakeNLProngMotherHadVisMV(beam, {211, 111})),
      fullCut, kNoShift, kUnweighted);

    // ── 2D: q0 vs |q| (all events + per mode) ────────────────────────────────
    Spectrum s2_q0_qmag(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag,
               "HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut,              kNoShift, kUnweighted);
    Spectrum s2_q0_qmag_QE(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag,
               "HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut && kCutQE,    kNoShift, kUnweighted);
    Spectrum s2_q0_qmag_Res(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag,
               "HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut && kCutRes,   kNoShift, kUnweighted);
    Spectrum s2_q0_qmag_DIS(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag,
               "HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut && kCutDIS,   kNoShift, kUnweighted);
    Spectrum s2_q0_qmag_MEC(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag,
               "HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut && kCutMEC,   kNoShift, kUnweighted);

    // ── NLS per event in HadVisE and q0 slices ───────────────────────────────
    const ana::Cut kCutHadVisE_0_100   = MakeQ0Cut(0.0,   0.100);
    const ana::Cut kCutHadVisE_100_200 = MakeQ0Cut(0.100, 0.200);
    const ana::Cut kCutHadVisE_200_400 = MakeQ0Cut(0.200, 0.400);
    const ana::Cut kCutHadVisE_400inf  = MakeQ0Cut(0.400, 9999.0);

    Spectrum sNSelPerEvent_hadVisE_0_100(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCut && kCutHadVisE_0_100,   kNoShift, kUnweighted);
    Spectrum sNSelPerEvent_hadVisE_100_200(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCut && kCutHadVisE_100_200, kNoShift, kUnweighted);
    Spectrum sNSelPerEvent_hadVisE_200_400(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCut && kCutHadVisE_200_400, kNoShift, kUnweighted);
    Spectrum sNSelPerEvent_hadVisE_400inf(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCut && kCutHadVisE_400inf,  kNoShift, kUnweighted);

    Spectrum sNSelPerEvent_q0Lo(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCut && kCutQ0Lo, kNoShift, kUnweighted);
    Spectrum sNTruePerEvent_q0Lo(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      fullCut && kCutQ0Lo, kNoShift, kUnweighted);
    Spectrum sNGENIEPerEvent_q0Lo(loader,
      HistAxis("N GENIE neutrons per event", kBinNSel, MakeNGENIENeutronsVar()),
      fullCut && kCutQ0Lo, kNoShift, kUnweighted);
    Spectrum sNSelPerEvent_q0Hi(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCut && kCutQ0Hi, kNoShift, kUnweighted);
    Spectrum sNTruePerEvent_q0Hi(loader,
      HistAxis("N true neutron prongs per event", kBinNSel, MakeTruePrimCountVar()),
      fullCut && kCutQ0Hi, kNoShift, kUnweighted);
    Spectrum sNGENIEPerEvent_q0Hi(loader,
      HistAxis("N GENIE neutrons per event", kBinNSel, MakeNGENIENeutronsVar()),
      fullCut && kCutQ0Hi, kNoShift, kUnweighted);

    // ── |q| by interaction mode (one entry per event) ─────────────────────────
    Spectrum sQmag_sel(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag),
      fullCut,               kNoShift, kUnweighted);
    Spectrum sQmag_QE(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag),
      fullCut && kCutQE,     kNoShift, kUnweighted);
    Spectrum sQmag_Res(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag),
      fullCut && kCutRes,    kNoShift, kUnweighted);
    Spectrum sQmag_DIS(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag),
      fullCut && kCutDIS,    kNoShift, kUnweighted);
    Spectrum sQmag_MEC(loader,
      HistAxis("|q| [GeV/c]", kBinQmag, kRecoQmag),
      fullCut && kCutMEC,    kNoShift, kUnweighted);

    // ── q0 (hadronic energy) by interaction mode ──────────────────────────────
    Spectrum sQ0_QE(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut && kCutQE,     kNoShift, kUnweighted);
    Spectrum sQ0_Res(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut && kCutRes,    kNoShift, kUnweighted);
    Spectrum sQ0_DIS(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut && kCutDIS,    kNoShift, kUnweighted);
    Spectrum sQ0_MEC(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV),
      fullCut && kCutMEC,    kNoShift, kUnweighted);

    // ── NLS count per event in muon pT slices ────────────────────────────────
    Spectrum sNSel_pt0(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCutNoPt && MakeMuonPtCut(kMuPtEdges[0], kMuPtEdges[1]), kNoShift, kUnweighted);
    Spectrum sNSel_pt1(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCutNoPt && MakeMuonPtCut(kMuPtEdges[1], kMuPtEdges[2]), kNoShift, kUnweighted);
    Spectrum sNSel_pt2(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCutNoPt && MakeMuonPtCut(kMuPtEdges[2], kMuPtEdges[3]), kNoShift, kUnweighted);
    Spectrum sNSel_pt3(loader,
      HistAxis("N_{NLS} per event", kBinNSel, MakeNLSCountVar(beam)),
      fullCutNoPt && MakeMuonPtCut(kMuPtEdges[3], kMuPtEdges[4]), kNoShift, kUnweighted);

    // ── nhits 1D for selected, signal, background (NLS) ──────────────────────
    Spectrum sNhits_sel(loader,
      MultiVarHistAxis("nhits", kBinNhits, MakeNLSVar(beam, kCatSelected,   fNhits)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNhits_sig(loader,
      MultiVarHistAxis("nhits", kBinNhits, MakeNLSVar(beam, kCatSignal,     fNhits)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNhits_bkg(loader,
      MultiVarHistAxis("nhits", kBinNhits, MakeNLSVar(beam, kCatBackground, fNhits)),
      fullCut, kNoShift, kUnweighted);

    // ── true-neutron prongs ───────────────────────────────────────────────────
    Spectrum sCalE_true(loader,
      MultiVarHistAxis("calE [MeV]",       kBinCalE,  MakeTrueNeutronVar(fCalE)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sLen_true(loader,
      MultiVarHistAxis("len [cm]",         kBinLen,   MakeTrueNeutronVar(fLen)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sVisE_true(loader,
      MultiVarHistAxis("true visE [MeV]",  kBinVisE,  MakeTrueNeutronVar(fVisE)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sMom_true(loader,
      MultiVarHistAxis("true |p| [GeV/c]", kBinMom,   MakeTrueNeutronVar(fMom)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sNhits_true(loader,
      MultiVarHistAxis("nhits",            kBinNhits, MakeTrueNeutronVar(fNhits)),
      fullCut, kNoShift, kUnweighted);

    // ── vertex displacement 1D ───────────────────────────────────────────────
    Spectrum sDispl_sel(loader,
      MultiVarHistAxis("displ [cm]", kBinDispl, MakeNLSDisplVar(beam, kCatSelected)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sDispl_sig(loader,
      MultiVarHistAxis("displ [cm]", kBinDispl, MakeNLSDisplVar(beam, kCatSignal)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sDispl_true(loader,
      MultiVarHistAxis("displ [cm]", kBinDispl, MakeTrueNeutronDisplVar()),
      fullCut, kNoShift, kUnweighted);

    // ── 2D: displacement vs nhits ─────────────────────────────────────────────
    Spectrum s2_displ_nhits_sel(loader,
      MultiVarHistAxis("displ [cm]", kBinDispl, MakeNLSDisplVar(beam, kCatSelected),
                       "nhits",       kBinNhits, MakeNLSVar(beam, kCatSelected, fNhits)),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_displ_nhits_q0Lo(loader,
      MultiVarHistAxis("displ [cm]", kBinDispl, MakeNLSDisplVar(beam, kCatSelected),
                       "nhits",       kBinNhits, MakeNLSVar(beam, kCatSelected, fNhits)),
      fullCut && kCutQ0Lo, kNoShift, kUnweighted);
    Spectrum s2_displ_nhits_q0Hi(loader,
      MultiVarHistAxis("displ [cm]", kBinDispl, MakeNLSDisplVar(beam, kCatSelected),
                       "nhits",       kBinNhits, MakeNLSVar(beam, kCatSelected, fNhits)),
      fullCut && kCutQ0Hi, kNoShift, kUnweighted);

    // ── mu_frac 1D ────────────────────────────────────────────────────────────
    Spectrum sMuFrac_sel(loader,
      MultiVarHistAxis("mu_frac", kBinMuFrac, MakeNLSMuFracVar(beam, kCatSelected)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sMuFrac_sig(loader,
      MultiVarHistAxis("mu_frac", kBinMuFrac, MakeNLSMuFracVar(beam, kCatSignal)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sMuFrac_bkg(loader,
      MultiVarHistAxis("mu_frac", kBinMuFrac, MakeNLSMuFracVar(beam, kCatBackground)),
      fullCut, kNoShift, kUnweighted);
    Spectrum sMuFrac_true(loader,
      MultiVarHistAxis("mu_frac", kBinMuFrac, MakeTrueNeutronMuFracVar()),
      fullCut, kNoShift, kUnweighted);

    // ── 2D: mu_frac vs nhits ──────────────────────────────────────────────────
    Spectrum s2_muFrac_nhits_sel(loader,
      MultiVarHistAxis("mu_frac", kBinMuFrac, MakeNLSMuFracVar(beam, kCatSelected),
                       "nhits",   kBinNhits,  MakeNLSVar(beam, kCatSelected,   fNhits)),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_muFrac_nhits_sig(loader,
      MultiVarHistAxis("mu_frac", kBinMuFrac, MakeNLSMuFracVar(beam, kCatSignal),
                       "nhits",   kBinNhits,  MakeNLSVar(beam, kCatSignal,     fNhits)),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_muFrac_nhits_bkg(loader,
      MultiVarHistAxis("mu_frac", kBinMuFrac, MakeNLSMuFracVar(beam, kCatBackground),
                       "nhits",   kBinNhits,  MakeNLSVar(beam, kCatBackground, fNhits)),
      fullCut, kNoShift, kUnweighted);

    // ── 2D: nhits vs pT ───────────────────────────────────────────────────────
    Spectrum s2_pT_nhits_sig(loader,
      MultiVarHistAxis("true p_{T} [GeV/c]", kBinMom,   MakeNLSVar(beam, kCatSignal,     fMom),
                       "nhits",               kBinNhits, MakeNLSVar(beam, kCatSignal,     fNhits)),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_pT_nhits_bkg(loader,
      MultiVarHistAxis("true p_{T} [GeV/c]", kBinMom,   MakeNLSVar(beam, kCatBackground, fMom),
                       "nhits",               kBinNhits, MakeNLSVar(beam, kCatBackground, fNhits)),
      fullCut, kNoShift, kUnweighted);


    // ── 2D: N reco signal vs N true primary ───────────────────────────────────
    Spectrum s2_recoSig_truePrim(loader,
      MultiVarHistAxis("N reco signal",  kBinRecoSig,  MakeRecoSigCountMV(beam),
                       "N true primary", kBinTruePrim, MakeTruePrimCountMV()),
      fullCut, kNoShift, kUnweighted);

    const ana::Cut fullCutNTrueGt0 = fullCut && MakeTrueNeutronPresentCut();
    Spectrum s2_recoSig_truePrim_nTgt0(loader,
      MultiVarHistAxis("N reco signal",  kBinRecoSig,  MakeRecoSigCountMV(beam),
                       "N true primary", kBinTruePrim, MakeTruePrimCountMV()),
      fullCutNTrueGt0, kNoShift, kUnweighted);

    // ── 2D: N reco vs N true primary ─────────────────────────────────────────────
    Spectrum s2_recoAll_truePrim(loader,
      MultiVarHistAxis("N reco",         kBinNSel,    MakeNLSCountMV(beam),
                       "N true primary", kBinTruePrim, MakeTruePrimCountMV()),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_recoAll_truePrim_nTgt0(loader,
      MultiVarHistAxis("N reco",         kBinNSel,    MakeNLSCountMV(beam),
                       "N true primary", kBinTruePrim, MakeTruePrimCountMV()),
      fullCutNTrueGt0, kNoShift, kUnweighted);

    // ── 2D: N NLS / N true per event vs muon pT (uniform binning for profiles) ─
    Spectrum s2_nLSCount_muPt(loader,
      MultiVarHistAxis("muon p_{T} [GeV/c]", kBinMom,  MakeMuonTruePtMV(),
                       "N NLS per event",     kBinNSel, MakeNLSCountMV(beam)),
      fullCutNoPt, kNoShift, kUnweighted);
    Spectrum s2_nTrue_vs_muPt(loader,
      HistAxis("muon p_{T} [GeV/c]", kBinMom,  MakeMuonTruePtVar(),
               "N true per event",   kBinNSel, MakeTruePrimCountVar()),
      fullCutNoPt, kNoShift, kUnweighted);

    // ── 2D: N NLS / N true / N signal NLS per event vs muon pT ──────────────
    Spectrum s2_nSigNLS_vs_muPt(loader,
      HistAxis("muon p_{T} [GeV/c]",  kBinMom,  MakeMuonTruePtVar(),
               "N signal NLS per event", kBinNSel, MakeSignalNLSCountVar(beam)),
      fullCutNoPt, kNoShift, kUnweighted);

    // ── 2D: N NLS / N true / N signal NLS per event vs |q| ──────────────────
    Spectrum s2_nLS_vs_qmag(loader,
      HistAxis("|q| [GeV/c]",       kBinQmag, kRecoQmag,
               "N_{NLS} per event",  kBinNSel, MakeNLSCountVar(beam)),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_nTrue_vs_qmag(loader,
      HistAxis("|q| [GeV/c]",       kBinQmag, kRecoQmag,
               "N true per event",   kBinNSel, MakeTruePrimCountVar()),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_nSigNLS_vs_qmag(loader,
      HistAxis("|q| [GeV/c]",            kBinQmag, kRecoQmag,
               "N signal NLS per event",  kBinNSel, MakeSignalNLSCountVar(beam)),
      fullCut, kNoShift, kUnweighted);

    // ── 2D: muon pT × N neutrons (fraction-vs-pT profiles) ───────────────────
    Spectrum s2_muPt_nRecoNLS(loader,
      HistAxis("muon p_{T} [GeV/c]", kBinMuPtFine,  MakeMuonTruePtVar(),
               "N NLS per event",    kBinNNeutrons, MakeNLSCountUncappedVar(beam)),
      fullCutNoPt, kNoShift, kUnweighted);
    Spectrum s2_muPt_nTruePrim(loader,
      HistAxis("muon p_{T} [GeV/c]", kBinMuPtFine,  MakeMuonTruePtVar(),
               "N true primary",     kBinNNeutrons, MakeTruePrimCountVar()),
      fullCutNoPt, kNoShift, kUnweighted);

    // ── 2D: |q| × N neutrons (fraction-vs-q profiles, same structure as muPt) ──
    Spectrum s2_q_nRecoNLS(loader,
      HistAxis("|q| [GeV/c]",     kBinQmag,      kRecoQmag,
               "N NLS per event", kBinNNeutrons, MakeNLSCountUncappedVar(beam)),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_q_nTruePrim(loader,
      HistAxis("|q| [GeV/c]",    kBinQmag,      kRecoQmag,
               "N true primary",  kBinNNeutrons, MakeTruePrimCountVar()),
      fullCut, kNoShift, kUnweighted);

    // ── 2D: N NLS / N true / N signal NLS per event vs q0 ────────────────────
    Spectrum s2_nLS_vs_q0(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV,
               "N_{NLS} per event",      kBinNSel, MakeNLSCountVar(beam)),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_nTrue_vs_q0(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV,
               "N true per event",       kBinNSel, MakeTruePrimCountVar()),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_nSigNLS_vs_q0(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV,
               "N signal NLS per event", kBinNSel, MakeSignalNLSCountVar(beam)),
      fullCut, kNoShift, kUnweighted);

    // ── 2D: q0 × N neutrons (fraction-vs-q0 profiles) ────────────────────────
    Spectrum s2_q0_nRecoNLS(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV,
               "N NLS per event",    kBinNNeutrons, MakeNLSCountUncappedVar(beam)),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_q0_nTruePrim(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV,
               "N true primary",     kBinNNeutrons, MakeTruePrimCountVar()),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_q0_nGENIEneutrons(loader,
      HistAxis("HadVisE [MeV]", kBinHadVisMeV, kNumuHadVisMeV,
               "N GENIE neutrons",   kBinNNeutrons, MakeNGENIENeutronsVar()),
      fullCut, kNoShift, kUnweighted);

    // ── 2D: N GENIE neutrons / N true prong neutrons vs HadVisE (for profiles) ─
    Spectrum s2_nGENIE_vs_q0(loader,
      HistAxis("HadVisE [MeV]",             kBinHadVisMeV, kNumuHadVisMeV,
               "N GENIE neutrons per event", kBinNSel,      MakeNGENIENeutronsVar()),
      fullCut, kNoShift, kUnweighted);
    Spectrum s2_nTrue_prongN_vs_q0(loader,
      HistAxis("HadVisE [MeV]",                  kBinHadVisMeV, kNumuHadVisMeV,
               "N true neutron prongs per event", kBinNSel,      MakeTruePrimCountVar()),
      fullCut, kNoShift, kUnweighted);

    loader.Go();

    TDirectory* beamDir = fOut->GetDirectory(beam.c_str());
    if (!beamDir) beamDir = fOut->mkdir(beam.c_str());
    TDirectory* dir = beamDir->mkdir(topoName.c_str());

    sCalE_sel.SaveTo(dir, "calE_selected");
    sLen_sel .SaveTo(dir, "len_selected");
    sVisE_sel.SaveTo(dir, "visE_selected");
    sMom_sel .SaveTo(dir, "mom_selected");

    sCalE_sig.SaveTo(dir, "calE_signal");
    sLen_sig .SaveTo(dir, "len_signal");
    sVisE_sig.SaveTo(dir, "visE_signal");
    sMom_sig .SaveTo(dir, "mom_signal");

    sCalE_bkg.SaveTo(dir, "calE_background");
    sLen_bkg .SaveTo(dir, "len_background");
    sVisE_bkg.SaveTo(dir, "visE_background");
    sMom_bkg .SaveTo(dir, "mom_background");

    sNSelPerEvent  .SaveTo(dir, "nsel_per_event");
    sNTruePerEvent .SaveTo(dir, "ntrue_per_event");
    sNGENIEPerEvent.SaveTo(dir, "ngenie_per_event");

    sQmag_sel_dist.SaveTo(dir, "qmag_dist_selected");
    sQmag_true    .SaveTo(dir, "qmag_dist_true");
    sQ2_sel       .SaveTo(dir, "q2_selected");
    sQ2_true      .SaveTo(dir, "q2_true");

    sQ0          .SaveTo(dir, "q0");
    sQ0_prim_muon   .SaveTo(dir, "q0_prim_muon");
    sQ0_prim_proton .SaveTo(dir, "q0_prim_proton");
    sQ0_prim_neutron.SaveTo(dir, "q0_prim_neutron");
    sQ0_prim_pion   .SaveTo(dir, "q0_prim_pion");
    sQ0_prim_other  .SaveTo(dir, "q0_prim_other");
    sHadVisE_prongmother_muon   .SaveTo(dir, "hadvisE_prongmother_muon");
    sHadVisE_prongmother_proton .SaveTo(dir, "hadvisE_prongmother_proton");
    sHadVisE_prongmother_neutron.SaveTo(dir, "hadvisE_prongmother_neutron");
    sHadVisE_prongmother_pion   .SaveTo(dir, "hadvisE_prongmother_pion");
    sHadVisE_prongmother_other  .SaveTo(dir, "hadvisE_prongmother_other");
    sHadVisE_nlmother_muon   .SaveTo(dir, "hadvisE_nlmother_muon");
    sHadVisE_nlmother_proton .SaveTo(dir, "hadvisE_nlmother_proton");
    sHadVisE_nlmother_neutron.SaveTo(dir, "hadvisE_nlmother_neutron");
    sHadVisE_nlmother_pion   .SaveTo(dir, "hadvisE_nlmother_pion");
    sHadVisE_nlmother_other  .SaveTo(dir, "hadvisE_nlmother_other");
    sNLcomp_total  .SaveTo(dir, "nlcomp_total");
    sNLcomp_muon   .SaveTo(dir, "nlcomp_muon");
    sNLcomp_proton .SaveTo(dir, "nlcomp_proton");
    sNLcomp_neutron.SaveTo(dir, "nlcomp_neutron");
    sNLcomp_pion   .SaveTo(dir, "nlcomp_pion");
    s2_q0_qmag   .SaveTo(dir, "q0_qmag");
    s2_q0_qmag_QE .SaveTo(dir, "q0_qmag_QE");
    s2_q0_qmag_Res.SaveTo(dir, "q0_qmag_Res");
    s2_q0_qmag_DIS.SaveTo(dir, "q0_qmag_DIS");
    s2_q0_qmag_MEC.SaveTo(dir, "q0_qmag_MEC");

    sNSelPerEvent_hadVisE_0_100  .SaveTo(dir, "nsel_hadVisE_0_100");
    sNSelPerEvent_hadVisE_100_200.SaveTo(dir, "nsel_hadVisE_100_200");
    sNSelPerEvent_hadVisE_200_400.SaveTo(dir, "nsel_hadVisE_200_400");
    sNSelPerEvent_hadVisE_400inf .SaveTo(dir, "nsel_hadVisE_400inf");

    sNSelPerEvent_q0Lo  .SaveTo(dir, "nsel_q0Lo");
    sNTruePerEvent_q0Lo .SaveTo(dir, "ntrue_q0Lo");
    sNGENIEPerEvent_q0Lo.SaveTo(dir, "ngenie_q0Lo");
    sNSelPerEvent_q0Hi  .SaveTo(dir, "nsel_q0Hi");
    sNTruePerEvent_q0Hi .SaveTo(dir, "ntrue_q0Hi");
    sNGENIEPerEvent_q0Hi.SaveTo(dir, "ngenie_q0Hi");

    sQmag_sel.SaveTo(dir, "qmag_selected");
    sQmag_QE .SaveTo(dir, "qmag_QE");
    sQmag_Res.SaveTo(dir, "qmag_Res");
    sQmag_DIS.SaveTo(dir, "qmag_DIS");
    sQmag_MEC.SaveTo(dir, "qmag_MEC");

    sQ0_QE .SaveTo(dir, "q0_QE");
    sQ0_Res.SaveTo(dir, "q0_Res");
    sQ0_DIS.SaveTo(dir, "q0_DIS");
    sQ0_MEC.SaveTo(dir, "q0_MEC");

    sNSel_pt0.SaveTo(dir, "nsel_muPt_bin0");
    sNSel_pt1.SaveTo(dir, "nsel_muPt_bin1");
    sNSel_pt2.SaveTo(dir, "nsel_muPt_bin2");
    sNSel_pt3.SaveTo(dir, "nsel_muPt_bin3");

    sNhits_sel.SaveTo(dir, "nhits_selected");
    sNhits_sig.SaveTo(dir, "nhits_signal");
    sNhits_bkg.SaveTo(dir, "nhits_background");

    sCalE_true .SaveTo(dir, "calE_trueneutron");
    sLen_true  .SaveTo(dir, "len_trueneutron");
    sVisE_true .SaveTo(dir, "visE_trueneutron");
    sMom_true  .SaveTo(dir, "mom_trueneutron");
    sNhits_true.SaveTo(dir, "nhits_trueneutron");

    sDispl_sel .SaveTo(dir, "displ_selected");
    sDispl_sig .SaveTo(dir, "displ_signal");
    sDispl_true.SaveTo(dir, "displ_trueneutron");

    sMuFrac_sel .SaveTo(dir, "muFrac_selected");
    sMuFrac_sig .SaveTo(dir, "muFrac_signal");
    sMuFrac_bkg .SaveTo(dir, "muFrac_background");
    sMuFrac_true.SaveTo(dir, "muFrac_trueneutron");

    s2_muFrac_nhits_sel.SaveTo(dir, "muFrac_nhits_selected");
    s2_muFrac_nhits_sig.SaveTo(dir, "muFrac_nhits_signal");
    s2_muFrac_nhits_bkg.SaveTo(dir, "muFrac_nhits_background");

    s2_pT_nhits_sig.SaveTo(dir, "pT_nhits_signal");
    s2_pT_nhits_bkg.SaveTo(dir, "pT_nhits_background");

    s2_recoSig_truePrim      .SaveTo(dir, "recoSig_truePrim");
    s2_recoSig_truePrim_nTgt0.SaveTo(dir, "recoSig_truePrim_nTgt0");
    s2_recoAll_truePrim      .SaveTo(dir, "recoAll_truePrim");
    s2_recoAll_truePrim_nTgt0.SaveTo(dir, "recoAll_truePrim_nTgt0");
    s2_nLSCount_muPt         .SaveTo(dir, "nLSCount_muPt");
    s2_muPt_nRecoNLS         .SaveTo(dir, "muPt_nRecoNLS");
    s2_muPt_nTruePrim        .SaveTo(dir, "muPt_nTruePrim");
    s2_q_nRecoNLS            .SaveTo(dir, "q_nRecoNLS");
    s2_q_nTruePrim           .SaveTo(dir, "q_nTruePrim");
    s2_nLS_vs_q0             .SaveTo(dir, "nLS_vs_q0");
    s2_nTrue_vs_q0           .SaveTo(dir, "nTrue_vs_q0");
    s2_nSigNLS_vs_q0         .SaveTo(dir, "nSigNLS_vs_q0");
    s2_q0_nRecoNLS           .SaveTo(dir, "q0_nRecoNLS");
    s2_q0_nTruePrim          .SaveTo(dir, "q0_nTruePrim");
    s2_q0_nGENIEneutrons     .SaveTo(dir, "q0_nGENIEneutrons");
    s2_nGENIE_vs_q0          .SaveTo(dir, "nGENIE_vs_q0");
    s2_nTrue_prongN_vs_q0    .SaveTo(dir, "nTrueProngN_vs_q0");

    // Raw TH2 saves — use the dataset's actual POT for normalization
    const double actualPOT = sMom_sel.POT();
    dir->cd();
    { TParameter<double> potParam("h2_pot", actualPOT); potParam.Write(); }

    auto saveH2 = [&](const ana::Spectrum& s, const char* name) {
      TH2* h = s.ToTH2(actualPOT);
      h->SetName(name); dir->cd(); h->Write(); delete h;
    };

    saveH2(s2_pT_nhits_sig,           "h2_pT_nhits_signal");
    saveH2(s2_pT_nhits_bkg,           "h2_pT_nhits_background");
    saveH2(s2_displ_nhits_sel,        "h2_displ_nhits_selected");
    saveH2(s2_displ_nhits_q0Lo,       "h2_displ_nhits_q0Lo");
    saveH2(s2_displ_nhits_q0Hi,       "h2_displ_nhits_q0Hi");
    saveH2(s2_muFrac_nhits_sel,       "h2_muFrac_nhits_selected");
    saveH2(s2_muFrac_nhits_sig,       "h2_muFrac_nhits_signal");
    saveH2(s2_muFrac_nhits_bkg,       "h2_muFrac_nhits_background");
    saveH2(s2_recoSig_truePrim,       "h2_recoSig_truePrim");
    saveH2(s2_recoSig_truePrim_nTgt0, "h2_recoSig_truePrim_nTgt0");
    saveH2(s2_recoAll_truePrim,       "h2_recoAll_truePrim");
    saveH2(s2_recoAll_truePrim_nTgt0, "h2_recoAll_truePrim_nTgt0");
    saveH2(s2_nLSCount_muPt,          "h2_nLSCount_muPt");
    saveH2(s2_nTrue_vs_muPt,          "h2_nTrue_vs_muPt");
    saveH2(s2_nSigNLS_vs_muPt,        "h2_nSigNLS_vs_muPt");
    saveH2(s2_nLS_vs_qmag,            "h2_nLS_vs_qmag");
    saveH2(s2_nTrue_vs_qmag,          "h2_nTrue_vs_qmag");
    saveH2(s2_nSigNLS_vs_qmag,        "h2_nSigNLS_vs_qmag");
    saveH2(s2_muPt_nRecoNLS,          "h2_muPt_nRecoNLS");
    saveH2(s2_muPt_nTruePrim,         "h2_muPt_nTruePrim");
    saveH2(s2_q_nRecoNLS,             "h2_q_nRecoNLS");
    saveH2(s2_q_nTruePrim,            "h2_q_nTruePrim");
    saveH2(s2_nLS_vs_q0,             "h2_nLS_vs_q0");
    saveH2(s2_nTrue_vs_q0,           "h2_nTrue_vs_q0");
    saveH2(s2_nSigNLS_vs_q0,         "h2_nSigNLS_vs_q0");
    saveH2(s2_q0_nRecoNLS,           "h2_q0_nRecoNLS");
    saveH2(s2_q0_nTruePrim,          "h2_q0_nTruePrim");
    saveH2(s2_q0_nGENIEneutrons,     "h2_q0_nGENIEneutrons");
    saveH2(s2_nGENIE_vs_q0,          "h2_nGENIE_vs_q0");
    saveH2(s2_nTrue_prongN_vs_q0,    "h2_nTrueProngN_vs_q0");
    saveH2(s2_q0_qmag,                "h2_q0_qmag");
    saveH2(s2_q0_qmag_QE,             "h2_q0_qmag_QE");
    saveH2(s2_q0_qmag_Res,            "h2_q0_qmag_Res");
    saveH2(s2_q0_qmag_DIS,            "h2_q0_qmag_DIS");
    saveH2(s2_q0_qmag_MEC,            "h2_q0_qmag_MEC");

    std::cout << beam << "/" << topoName << " done.\n";
  }
}

void make_neutron_spectra_q0()
{
  const ana::Weight Wei = ana::kPPFXFluxCVWgt * ana::kXSecCVWgt2024;
  const ana::Cut baseCutNoPt =
    ana::xsec::numubarcc::kQualityCut             &&
    ana::xsec::numubarcc::kContainmentCut         &&
    ana::xsec::numubarcc::kRecoVtxNumuFiducialCut &&
    ana::xsec::numubarcc::kMuonIDCut;

  const ana::Cut baseCut = baseCutNoPt && MakeMuonPtCut(0.0, 0.8);

  const std::string defFHC =
    "prod_sumdecaf_R20-11-25-prod5.1reco.t_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_menate_v1_numu2020";
  const std::string defRHC =
    "prod_sumdecaf_R20-11-25-prod5.1reco.t_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_menate_v1_numu2020";

  // Two q0 samples: low (q0 < 0.3 GeV) and high (q0 >= 0.3 GeV),
  // both using the inclusive neutron topology.
  using CutFn = std::function<ana::Cut(const std::string&)>;
  const std::vector<std::pair<std::string, CutFn>> samples = {
    {"sample1_q0Lo", [](const std::string& b){ return kTopoInclusive(b) && MakeQ0Cut(0.0,  0.10);   }},
    {"sample2_q0Hi", [](const std::string& b){ return kTopoInclusive(b) && MakeQ0Cut(0.10, 9999.0); }},
    {"sample3_full", [](const std::string& b){ return kTopoInclusive(b);                             }},
  };

  TFile* fOut = new TFile("make_neutron_spectra_q0.root", "RECREATE");

  for (const auto& sample : samples) {
    const std::string& sampleName = sample.first;
    const CutFn&       cutFn      = sample.second;
    ProcessBeam(defFHC, "FHC", baseCut, baseCutNoPt, cutFn("FHC"), sampleName, fOut);
    ProcessBeam(defRHC, "RHC", baseCut, baseCutNoPt, cutFn("RHC"), sampleName, fOut);
  }

  fOut->Write();
  fOut->Close();
  delete fOut;
  std::cout << "Saved to make_neutron_spectra_q0.root\n";
}

#endif
