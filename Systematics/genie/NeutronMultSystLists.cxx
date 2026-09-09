// -------------------------------------------------------------------------------
// Systematics relevant to neutron-multiplicity analysis, ranked by priority.
// Subset of getAllXsecNuTruthSysts_2024(); see technote Section X for the
// physics justification behind each tier.
//
// Split by physics category (FSI/MEC/RES/QE/DIS), matching the rows of the
// technote's "Systematics Relevant to Neutron Multiplicity" table, so that
// make_neutron_spectra_systs_diagnostics.C can pair each category with the
// diagnostic variable(s) that actually track its effect (see the technote's
// "Diagnostic Variables by Systematic Category" table). The *Priority()
// tiers below are unchanged in meaning and content -- just now composed
// from the category lists instead of listing systs directly.

  /// FSI: mean free path / rescattering / absorption -- the dominant
  /// mechanism for producing >1 visible neutron via secondary scattering.
  std::vector<const ISyst*> getNeutronMultSysts_FSI()
  {
    return std::vector<const ISyst*>
    {
      &khNFSISyst2024_EV1,
      &khNFSISyst2024_EV2,
      &khNFSISyst2024_EV3,
      &khNFSISyst2024_MFP,
      &kFormZone_2024,
    };
  }

  /// MEC: sets n vs. p split for 2p2h (init-state fraction) and reshapes
  /// the 2p2h q0 spectrum, which correlates with neutron visibility.
  std::vector<const ISyst*> getNeutronMultSysts_MEC()
  {
    return std::vector<const ISyst*>
    {
      &kMECInitStateNPFracSyst2020Nu,
      &kMECInitStateNPFracSyst2020AntiNu,
      &kMECShapeSyst2024Nu,
      &kMECShapeSyst2024AntiNu,
      &kMECEnuShapeSyst2020Nu,
      &kMECEnuShapeSyst2020AntiNu,
    };
  }

  /// RES: sets target nucleon (p vs n), the Delta-decay n/p partition, and
  /// how much of the sample sits in the Delta-dominated regime vs. higher
  /// mass resonances.
  std::vector<const ISyst*> getNeutronMultSysts_RES()
  {
    return std::vector<const ISyst*>
    {
      &kRESvpvnRatioNuXSecSyst,
      &kRESvpvnRatioNubarXSecSyst,
      GetGenieKnobSyst(rwgt::fReweightTheta_Delta2Npi),   ///< Delta -> N + pi angular distribution
      &kRESLowQ2SuppressionSyst2020,
      &kRESDeltaScaleSyst,
      &kRESOtherScaleSyst,
      GetGenieKnobSyst(rwgt::fReweightMaCCRES),
      GetGenieKnobSyst(rwgt::fReweightMvCCRES),
      GetGenieKnobSyst(rwgt::fReweightMaNCRES),
      GetGenieKnobSyst(rwgt::fReweightMvNCRES),
    };
  }

  /// QE: CCQE rate/shape -- the dominant proton-free/pion-free channel in
  /// RHC, so it directly sets the "clean neutron" population size.
  std::vector<const ISyst*> getNeutronMultSysts_QE()
  {
    return std::vector<const ISyst*>
    {
      &kRPACCQEEnhSyst2020,
      &kRPACCQESuppSyst2020,
      &kZExpEV1Syst2020,
      &kZExpEV2Syst2020,
      &kZExpEV3Syst2020,
      &kZExpEV4Syst2020,
      GetGenieKnobSyst(rwgt::fReweightZNormCCQE),
    };
  }

  /// DIS/hadronization: soft n-pion (32 knobs) + AGKY/Bodek-Yang
  /// hadronization -- mainly affects higher-multiplicity, higher-W final
  /// states, so its effect is concentrated in the multiplicity tail.
  std::vector<const ISyst*> getNeutronMultSysts_DIS()
  {
    return std::vector<const ISyst*>
    {
      &kDISvpCC0pi_2020,    &kDISvpCC1pi_2020,    &kDISvpCC2pi_2020,    &kDISvpCC3pi_2020,
      &kDISvpNC0pi_2020,    &kDISvpNC1pi_2020,    &kDISvpNC2pi_2020,    &kDISvpNC3pi_2020,
      &kDISvnCC0pi_2020,    &kDISvnCC1pi_2020,    &kDISvnCC2pi_2020,    &kDISvnCC3pi_2020,
      &kDISvnNC0pi_2020,    &kDISvnNC1pi_2020,    &kDISvnNC2pi_2020,    &kDISvnNC3pi_2020,
      &kDISvbarpCC0pi_2020, &kDISvbarpCC1pi_2020, &kDISvbarpCC2pi_2020, &kDISvbarpCC3pi_2020,
      &kDISvbarpNC0pi_2020, &kDISvbarpNC1pi_2020, &kDISvbarpNC2pi_2020, &kDISvbarpNC3pi_2020,
      &kDISvbarnCC0pi_2020, &kDISvbarnCC1pi_2020, &kDISvbarnCC2pi_2020, &kDISvbarnCC3pi_2020,
      &kDISvbarnNC0pi_2020, &kDISvbarnNC1pi_2020, &kDISvbarnNC2pi_2020, &kDISvbarnNC3pi_2020,
      &kDISNuHadronQ1Syst,
      &kDISNuBarHadronQ0Syst,
      GetGenieKnobSyst(rwgt::fReweightAGKY_xF1pi),
      GetGenieKnobSyst(rwgt::fReweightAGKY_pT1pi),
      GetGenieKnobSyst(rwgt::fReweightAhtBY),
      GetGenieKnobSyst(rwgt::fReweightBhtBY),
      GetGenieKnobSyst(rwgt::fReweightCV1uBY),
      GetGenieKnobSyst(rwgt::fReweightCV2uBY),
    };
  }

  /// Highest/high priority: directly set the neutron/proton split, drive the
  /// FSI cascade that produces extra visible neutrons, or reshape the q0
  /// (hadronic energy transfer) spectrum that neutron visibility tracks.
  /// = FSI + MEC + RES.
  std::vector<const ISyst*> getNeutronMultSysts_HighPriority()
  {
    std::vector<const ISyst*> systs = getNeutronMultSysts_FSI();
    auto mec = getNeutronMultSysts_MEC();
    auto res = getNeutronMultSysts_RES();
    systs.insert(systs.end(), mec.begin(), mec.end());
    systs.insert(systs.end(), res.begin(), res.end());
    return systs;
  }

  /// Medium priority: shift QE/DIS mode composition, with an indirect
  /// effect on the neutron sample. = QE + DIS.
  std::vector<const ISyst*> getNeutronMultSysts_MediumPriority()
  {
    std::vector<const ISyst*> systs = getNeutronMultSysts_QE();
    auto dis = getNeutronMultSysts_DIS();
    systs.insert(systs.end(), dis.begin(), dis.end());
    return systs;
  }

  /// Low priority: no nucleon knockout (COH), out-of-selection channel
  /// (NC elastic, nu_e-specific radiative/2nd-class-current knobs), or too
  /// small a phase space (rare Delta radiative/eta-decay branching) to
  /// matter for a numu/numubar CC neutron-multiplicity measurement.
  std::vector<const ISyst*> getNeutronMultSysts_LowPriority()
  {
    return std::vector<const ISyst*>
    {
      GetGenieKnobSyst(rwgt::fReweightMaNCEL),
      GetGenieKnobSyst(rwgt::fReweightEtaNCEL),
      &kCOHCCScaleSyst2018,
      &kCOHNCScaleSyst2018,
      &kRadCorrNue,
      &kRadCorrNuebar,
      &k2ndClassCurrs,
      GetGenieKnobSyst(rwgt::fReweightBR1gamma),
      GetGenieKnobSyst(rwgt::fReweightBR1eta),
    };
  }

  /// Convenience: all three tiers concatenated (equivalent to the full
  /// getAllXsecNuTruthSysts_2024() list, re-ordered by priority).
  std::vector<const ISyst*> getNeutronMultSysts_All()
  {
    std::vector<const ISyst*> systs = getNeutronMultSysts_HighPriority();
    auto med = getNeutronMultSysts_MediumPriority();
    auto low = getNeutronMultSysts_LowPriority();
    systs.insert(systs.end(), med.begin(), med.end());
    systs.insert(systs.end(), low.begin(), low.end());
    return systs;
  }
