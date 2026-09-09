// -------------------------------------------------------------------------------
// Shared plumbing for the make_neutron_spectra_systs_diagnostics_{high,medium,
// low}.C scripts: the diagnostic Vars/Binnings not already in
// selection_vars.h/cuts.C, plus the generic "nominal + per-syst +-1 sigma"
// spectrum builder and its ROOT-file writer. #include this (after cuts.C
// and NeutronMultSystLists.cxx) inside each tier script's own anonymous
// namespace -- kept as one file so the three tier scripts can't drift out
// of sync on how a diagnostic is defined or saved.

  // How many of the event's primary (pre-FSI) GENIE neutrons never left a
  // visible prong (0 = every primary neutron is visible; >0 = that many
  // were absorbed/invisible). Same physics content as the 2D confusion
  // matrix in make_mc_spectra.C, collapsed to one axis so it can be
  // shifted/ranked by the same 1D machinery as every other diagnostic here.
  const ana::Var kNPrimMinusVisibleVar([](const caf::SRProxy* sr) -> double {
    if (!sr || !sr->vtx.elastic.IsValid || sr->mc.nnu == 0) return -1.0;
    int nPrim = 0;
    for (const auto& prim : sr->mc.nu[0].prim)
      if (prim.pdg == 2112) ++nPrim;
    const int nVis = (int)GroupVisibleNeutronsByProngCount(sr).nProngs.size();
    return (double)(nPrim - nVis);
  });

  // Final-state topology bin, mirroring cuts.C's kFSI_* grouping cuts
  // exactly (same nProton/nNeutron/nPion thresholds) so this single
  // categorical histogram reproduces check_final_states.C's proton x pion
  // table as one fillable axis:
  //   0: kFSI_1Neutron            (0p, 1n,  0pi)
  //   1: kFSI_2PlusNeutrons       (0p, >=2n,0pi)
  //   2: kFSI_1Neutron_1Pion      (0p, 1n,  1pi)
  //   3: kFSI_1Proton_0Pion       (1p, >=1n,0pi)
  //   4: kFSI_2PlusProtons_0Pion  (>=2p,>=1n,0pi)
  //   5: kFSI_1Proton_NPion       (1p, >=1n,>=1pi)
  //   6: kFSI_2PlusProtons_NPion  (>=2p,>=1n,>=1pi)
  //   7: no primary neutron, OR the uncovered gap (0p with >=2 pions, or
  //      >=2 neutrons plus any pion) -- see cuts.C's kFSI_* comments.
  const ana::Var kFinalStateTopologyBinVar([](const caf::SRProxy* sr) -> double {
    if (!sr || sr->mc.nnu == 0) return 7.0;
    int nP = 0, nN = 0, nPi = 0;
    for (const auto& prim : sr->mc.nu[0].prim) {
      const int pdg = prim.pdg;
      if      (pdg == 2212)                                  ++nP;
      else if (pdg == 2112)                                  ++nN;
      else if (std::abs(pdg) == 211 || std::abs(pdg) == 111) ++nPi;
    }
    if (nN < 1)                       return 7.0;
    if (nP == 0 && nPi == 0)          return nN == 1 ? 0.0 : 1.0;
    if (nP == 0 && nPi == 1 && nN==1) return 2.0;
    if (nP == 1 && nPi == 0)          return 3.0;
    if (nP >= 2 && nPi == 0)          return 4.0;
    if (nP == 1 && nPi >= 1)          return 5.0;
    if (nP >= 2 && nPi >= 1)          return 6.0;
    return 7.0;
  });

  const ana::Binning kBinNeutronKE           = ana::Binning::Simple(40, 0.0, 100.0);  // MeV
  const ana::Binning kBinProngsPerNeutron    = ana::Binning::Simple(3, 0.5, 3.5);     // 1, 2, 3+
  const ana::Binning kBinNPrimMinusVis       = ana::Binning::Simple(8, -0.5, 7.5);
  const ana::Binning kBinFinalStateTopology  = ana::Binning::Simple(8, -0.5, 7.5);

  // ── generic nominal + (+-1 sigma per syst) spectrum builder ─────────────
  struct ShiftPair { std::unique_ptr<Spectrum> p1sigma, m1sigma; };

  template <class AxisT>
  void FillShiftedSpectra(SpectrumLoader& loader, const AxisT& axis, const ana::Cut& cut,
                           const std::vector<const ISyst*>& systs,
                           std::unique_ptr<Spectrum>& nominal,
                           std::map<std::string, ShiftPair>& shifted)
  {
    nominal = std::make_unique<Spectrum>(loader, axis, cut, kNoShift, kUnweighted);
    for (const ISyst* syst : systs) {
      ShiftPair& sp = shifted[syst->ShortName()];
      sp.p1sigma = std::make_unique<Spectrum>(loader, axis, cut, ana::SystShifts(syst, +1.0), kUnweighted);
      sp.m1sigma = std::make_unique<Spectrum>(loader, axis, cut, ana::SystShifts(syst, -1.0), kUnweighted);
    }
  }

  // Saves one diagnostic's nominal + per-syst +-1 sigma spectra into
  // catDir/nominal/<objName>, catDir/systs/<syst>/p1sigma/<objName>,
  // catDir/systs/<syst>/m1sigma/<objName> -- multiple diagnostics for the
  // same category share the same per-syst directories, one object per
  // diagnostic, same convention as ntrue_per_event/ngenie_per_event coexist
  // in make_neutron_spectra_systs.C.
  void SaveDiagnostic(TDirectory* catDir, const std::string& objName,
                       std::unique_ptr<Spectrum>& nominal,
                       std::map<std::string, ShiftPair>& shifted,
                       const std::vector<const ISyst*>& systs)
  {
    TDirectory* nomDir = catDir->GetDirectory("nominal");
    if (!nomDir) nomDir = catDir->mkdir("nominal");
    nominal->SaveTo(nomDir, objName.c_str());

    TDirectory* systsDir = catDir->GetDirectory("systs");
    if (!systsDir) systsDir = catDir->mkdir("systs");
    for (const ISyst* syst : systs) {
      const std::string& name = syst->ShortName();
      TDirectory* sDir = systsDir->GetDirectory(name.c_str());
      if (!sDir) sDir = systsDir->mkdir(name.c_str());
      TDirectory* p1Dir = sDir->GetDirectory("p1sigma");
      if (!p1Dir) p1Dir = sDir->mkdir("p1sigma");
      TDirectory* m1Dir = sDir->GetDirectory("m1sigma");
      if (!m1Dir) m1Dir = sDir->mkdir("m1sigma");

      shifted[name].p1sigma->SaveTo(p1Dir, objName.c_str());
      shifted[name].m1sigma->SaveTo(m1Dir, objName.c_str());
    }
  }

  TDirectory* GetOrMkdir(TDirectory* parent, const std::string& name)
  {
    TDirectory* d = parent->GetDirectory(name.c_str());
    return d ? d : parent->mkdir(name.c_str());
  }

  // Same base event selection and GENIE-truth-matched samples every tier
  // script uses.
  const ana::Cut kNeutronMultBaseCut =
    ana::xsec::numubarcc::kQualityCut             &&
    ana::xsec::numubarcc::kContainmentCut         &&
    ana::xsec::numubarcc::kRecoVtxNumuFiducialCut &&
    ana::xsec::numubarcc::kMuonIDCut;

  const std::string kNeutronMultDefFHC =
    "prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_filematchedSystematics_nominal_v1";
  const std::string kNeutronMultDefRHC =
    "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_ndphysics_systs_filematch_v1";
