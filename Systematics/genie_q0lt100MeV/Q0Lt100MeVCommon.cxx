// -------------------------------------------------------------------------------
// Shared restriction for the q0lt100MeV variant of the category diagnostics
// study: every diagnostic here is evaluated within HadVisE < 100 MeV
// (ana::kNumuHadVisE(sr) < 100 MeV), the low-hadronic-energy region the
// original spectrums/make_neutron_spectra_systs.C studied for
// ntrue_per_event alone. This folder asks the same "which systematics are
// relevant, and for which diagnostic" question as ../genie/, but restricted
// to that region instead of the fully inclusive sample.
//
// #include this (after ../genie/NeutronMultSystLists.cxx and
// ../genie/NeutronMultSystDiagnostics.cxx, which this depends on for
// kNeutronMultBaseCut/MakeVisHadEMeVCut) inside each tier script's own
// anonymous namespace.

  // HadVisE < 100 MeV on top of the same base event selection every tier
  // script uses. Every diagnostic in this folder is built on this cut (or
  // this cut && a topology/species restriction) instead of the plain
  // kNeutronMultBaseCut used in ../genie/.
  const ana::Cut kQ0Lt100MeVCut = kNeutronMultBaseCut && MakeVisHadEMeVCut(0.0, 100.0);

  // HadVisE itself is one of the diagnostics (hadvise_q0 / hadvise_cleanccqe
  // / hadvise_inclusive) -- binning it 0-300 MeV (../genie/'s
  // kBinHadVisMeV300) while also cutting HadVisE < 100 MeV would leave the
  // 100-300 MeV bins permanently empty. Same 10 MeV/bin width as
  // kBinHadVisMeV300, just over the restricted range.
  const ana::Binning kBinHadVisMeV100 = ana::Binning::Simple(10, 0.0, 100.0);
