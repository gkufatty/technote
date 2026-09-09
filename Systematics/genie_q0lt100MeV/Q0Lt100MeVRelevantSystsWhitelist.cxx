// -------------------------------------------------------------------------------
// q0lt100MeV variant of ../genie/RelevantSystsWhitelist.cxx: systematics
// retained after re-running the SAME category-by-category diagnostics study
// (make_neutron_spectra_systs_diagnostics_{high,medium,low}.C +
// plot_neutron_spectra_systs_diagnostics_ranking.C, this folder's copies)
// with every diagnostic restricted to HadVisE < 100 MeV (Q0Lt100MeVCommon.cxx)
// instead of the fully inclusive sample. Same chi2_shape > 5 ("Strong")
// threshold, globally re-ranked across all 77 systematics in
// getNeutronMultSysts_All() -- not just a per-category cut.
//
// Comparison against ../genie/RelevantSystsWhitelist.cxx's 26-syst whitelist
// (max chi2_shape across the same 14 diagnostics, inclusive sample vs. this
// HadVisE < 100 MeV sample; see neutron_spectra_systs_diagnostics_*_metrics.csv
// in both folders for the full per-systematic numbers):
//
//   FSI (0/5 kept, down from 5/5) -- EVERY FSI knob drops below chi2>5 once
//     HadVisE is cut below 100 MeV: hNFSI_MFP_2024 179.5 -> 2.00,
//     FormZone2024 77.3 -> 0.17, the three FateFrac EV knobs 5.5-12.2 -> all
//     <1. Physically expected: these knobs move MFP/rescattering/absorption,
//     which mostly acts by adding secondary-scattering hadronic activity --
//     i.e. by pushing events to HIGHER HadVisE. Restricting to the clean,
//     low-q0 region removes almost exactly the phase space where FSI
//     matters, so it is NOT a relevant systematic for a q0<100 MeV analysis.
//   MEC (4/6 kept, down from 6/6) -- kept MECInitStateNPFrac2020{Nu,AntiNu}
//     (15.8, 12.5) and MECShape2024{Nu,AntiNu} (30.1, 34.6 -- actually LARGER
//     than inclusive's 66.4/19.5 split, since MEC's q0 dependence is
//     concentrated at low q0); dropped MECEnuShape2020Nu (14.0 -> 3.0) and
//     MECEnuShape2020AntiNu (9.2 -> 4.8, just below threshold).
//   RES (3/10 kept, down from 6/10) -- kept MaCCRES (58.5 -> 7.1),
//     RESDeltaScaleSyst (55.2 -> 19.8), LowQ2RESSupp2020 (5.8 -> 8.7, grows
//     under the cut); dropped MvCCRES (31.9 -> 2.4), RESOtherScaleSyst
//     (14.5 -> 0.2), MaNCRES (9.6 -> 2.1).
//   QE (2/3 kept, down from 3/7) -- kept ZNormCCQE (104.2 -> 13.3) and
//     RPAShapesupp2020 (29.5 -> 11.6); dropped ZExpAxialFFSyst2020_EV1
//     (8.7 -> 0.76) -- the low-q0 cut removes the tail sensitivity that let
//     it clear the bar inclusively.
//   DIS (2/6 kept, down from 6/40) -- kept DISNuHadronQ1Syst (184.2 -> 18.7)
//     and DISvnCC1pi_2020 (166.6 -> 19.4); dropped DISvnCC2pi_2020
//     (16.6 -> 0.0), DISvpCC2pi_2020 (17.0 -> 0.71), DISvpCC3pi_2020
//     (6.3 -> 0.0), DISNuBarHadronQ0Syst (14.3 -> 1.3) -- DIS multipion
//     effects live mostly at higher hadronic energy, outside this cut.
//   LowPriority (0/9 kept) -- unchanged, still negligible.
//
// Net: 11 of 77 systematics (down from 26), and NO systematic newly clears
// the bar here that didn't already clear it inclusively -- restricting to
// HadVisE < 100 MeV only ever suppresses significance in this study, never
// reveals a new relevant systematic. FSI drops out entirely; MEC/RES/QE/DIS
// each lose roughly half their inclusive members.
//
// #include this (after ../genie/NeutronMultSystLists.cxx and
// Q0Lt100MeVCommon.cxx) inside the including file's own anonymous namespace.
// Shared by make_neutron_spectra_systs_multiverse_filtered.C and any ranking
// plot in this folder so they can't drift apart on which 11 systematics
// count as "relevant" for the q0<100 MeV phase space.

  const std::set<std::string> kRelevantSysts = {
    // MEC (4/6)
    "MECInitStateNPFrac2020Nu", "MECInitStateNPFrac2020AntiNu",
    "MECShape2024Nu", "MECShape2024AntiNu",
    // RES (3/10)
    "MaCCRES", "RESDeltaScaleSyst", "LowQ2RESSupp2020",
    // QE (2/7)
    "ZNormCCQE", "RPAShapesupp2020",
    // DIS (2/40)
    "DISNuHadronQ1Syst", "DISvnCC1pi_2020",
    // FSI (0/5) and LowPriority (0/9): none clear chi2>5 in this phase space.
  };

  // Category tag per whitelisted systematic, for plotting/tables (color
  // coding, grouping) -- kept in sync with kRelevantSysts above by
  // construction (same grouping, just as a lookup instead of a comment).
  const std::map<std::string, std::string> kRelevantSystCategory = {
    {"MECInitStateNPFrac2020Nu", "MEC"}, {"MECInitStateNPFrac2020AntiNu", "MEC"},
    {"MECShape2024Nu", "MEC"}, {"MECShape2024AntiNu", "MEC"},

    {"MaCCRES", "RES"}, {"RESDeltaScaleSyst", "RES"}, {"LowQ2RESSupp2020", "RES"},

    {"ZNormCCQE", "QE"}, {"RPAShapesupp2020", "QE"},

    {"DISNuHadronQ1Syst", "DIS"}, {"DISvnCC1pi_2020", "DIS"},
  };
