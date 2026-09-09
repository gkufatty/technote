// -------------------------------------------------------------------------------
// Systematics retained after the category-by-category diagnostics study
// (make_neutron_spectra_systs_diagnostics_{high,medium,low}.C +
// plot_neutron_spectra_systs_diagnostics_ranking.C), all three tiers
// checked. Kept: chi2_shape > 5 ("Strong") in at least one diagnostic,
// globally re-ranked across all 77 systematics in
// getNeutronMultSysts_All() -- not just a per-category cut. See
// make_neutron_spectra_systs_multiverse_filtered.C's header comment for the
// full per-category derivation (what was dropped and why) and
// neutron_spectra_systs_diagnostics_*_metrics.csv for the actual
// per-systematic chi2 values behind this.
//
// #include this (after NeutronMultSystLists.cxx) inside the including
// file's own anonymous namespace. Shared by
// make_neutron_spectra_systs_multiverse_filtered.C and
// plot_relevant_systs_ranking.C so the two can't drift apart on which 26
// systematics count as "relevant".

  const std::set<std::string> kRelevantSysts = {
    // FSI (5/5)
    "hNFSI_MFP_2024", "FormZone2024",
    "hNFSI_FateFracEV1_2024", "hNFSI_FateFracEV2_2024", "hNFSI_FateFracEV3_2024",
    // MEC (6/6)
    "MECInitStateNPFrac2020Nu", "MECInitStateNPFrac2020AntiNu",
    "MECShape2024Nu", "MECShape2024AntiNu",
    "MECEnuShape2020Nu", "MECEnuShape2020AntiNu",
    // RES (6/10 -- dropped RESvpvnRatioNuXSecSyst, RESvpvnRatioNubarXSecSyst,
    // Theta_Delta2Npi, MvNCRES)
    "MaCCRES", "RESDeltaScaleSyst", "MvCCRES", "RESOtherScaleSyst",
    "MaNCRES", "LowQ2RESSupp2020",
    // QE (3/7 -- dropped RPAShapeenh2020, ZExpAxialFFSyst2020_EV{2,3,4})
    "ZNormCCQE", "RPAShapesupp2020", "ZExpAxialFFSyst2020_EV1",
    // DIS (6/40)
    "DISNuHadronQ1Syst", "DISNuBarHadronQ0Syst",
    "DISvnCC1pi_2020", "DISvnCC2pi_2020",
    "DISvpCC2pi_2020", "DISvpCC3pi_2020",
    // LowPriority (0/9 -- all confirmed negligible)
  };

  // Category tag per whitelisted systematic, for plotting/tables (color
  // coding, grouping) -- kept in sync with kRelevantSysts above by
  // construction (same grouping, just as a lookup instead of a comment).
  const std::map<std::string, std::string> kRelevantSystCategory = {
    {"hNFSI_MFP_2024", "FSI"}, {"FormZone2024", "FSI"},
    {"hNFSI_FateFracEV1_2024", "FSI"}, {"hNFSI_FateFracEV2_2024", "FSI"},
    {"hNFSI_FateFracEV3_2024", "FSI"},

    {"MECInitStateNPFrac2020Nu", "MEC"}, {"MECInitStateNPFrac2020AntiNu", "MEC"},
    {"MECShape2024Nu", "MEC"}, {"MECShape2024AntiNu", "MEC"},
    {"MECEnuShape2020Nu", "MEC"}, {"MECEnuShape2020AntiNu", "MEC"},

    {"MaCCRES", "RES"}, {"RESDeltaScaleSyst", "RES"}, {"MvCCRES", "RES"},
    {"RESOtherScaleSyst", "RES"}, {"MaNCRES", "RES"}, {"LowQ2RESSupp2020", "RES"},

    {"ZNormCCQE", "QE"}, {"RPAShapesupp2020", "QE"}, {"ZExpAxialFFSyst2020_EV1", "QE"},

    {"DISNuHadronQ1Syst", "DIS"}, {"DISNuBarHadronQ0Syst", "DIS"},
    {"DISvnCC1pi_2020", "DIS"}, {"DISvnCC2pi_2020", "DIS"},
    {"DISvpCC2pi_2020", "DIS"}, {"DISvpCC3pi_2020", "DIS"},
  };
