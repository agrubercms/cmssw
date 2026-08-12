from __future__ import absolute_import, division, print_function
import ROOT
from DataFormats.FWLite import Events
from modules.EvtData import EvtData, add_product
from modules.GenTools import load_fwlitelibs
load_fwlitelibs()

import glob
import sys
from array import array

# Get input and output files from command line or use defaults
if len(sys.argv) > 2:
    files = [sys.argv[1]]  # Use the provided input file
    output_filename = sys.argv[2]  # Use provided output filename
elif len(sys.argv) > 1:
    files = [sys.argv[1]]  # Use the provided input file
    output_filename = "histograms.root"  # Use default output filename
else:
    # Default file patterns
    #files1 = glob.glob('/eos/user/a/agruber/samples/HLT_Upgrade_L1filter/VBF_ParT_4global/Phase2_L1P2GT_HLT*.root')
    files1 = glob.glob('/eos/user/a/agruber/samples/HLT_Upgrade_L1filter/ParT_unfiltered_16_1_1/Phase2_L1P2GT_HLT_*.root')
    #files1 = glob.glob('/afs/cern.ch/user/a/agruber/public/forOsman/Phase2_L1P2GT_HLT_*.root')
    #files1 = glob.glob('/eos/user/a/agruber/samples/HLT_Upgrade_L1filter/VBF_ParT_300k/Phase2_L1P2GT_HLT_*.root')
    #files1 = glob.glob('/eos/user/a/agruber/samples/HLT_Upgrade_L1filter/ParT_splitSamples/Phase2_L1P2GT_HLT*.root')
    #files1 = glob.glob('/data/agruber/tauPOG/HLT/l1seed/CMSSW_15_1_0_pre4/src/29727.0_QQToHToTauTau_14TeV+Run4D110/step2.root')
    #files1 = glob.glob('/afs/cern.ch/work/a/agruber/private/HLT_Phase2/CMSSW_15_0_0_pre1/src/Phase2_L1P2GT_HLT*.root')
    #files2 = glob.glob('/eos/user/a/agruber/samples/HLT_Upgrade_L1filter/VBF_HTauTau_14_2_0_pre3_l1TauFilterNoPtCut_ext/Phase2_L1P2GT_HLT*.root')
    files = files1 #+ files2
    output_filename = "histograms.root"

ROOT.gSystem.Load('libDataFormatsL1Trigger.so')
events = Events(files)

# Define trigger paths and their corresponding filters
trigger_paths = [
    #["L1P2GT_DoubleNNTau52", "hltL1P2GTTau"],
    ["HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1", "hltHpsDoublePFTau40TrackPt1MediumChargedIsolation"],
    ["HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1", "hltHpsDoublePFTau35MediumDitauWPDeepTau"],
    ["HLT_DoublePNetTauh", "hltDoublePFJets30PNetTauhTagMediumWPL2DoubleTau"],
    ["HLT_DoubleMediumPFPuppiParTTauh30_eta2p1", "hltDoublePFJets30ParTTauhTagMediumWPL2DoubleTau"],
]

# Define the binning configuration in one place
pt_bins = array('d', [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 120, 140, 160, 200])
histogram_config = {
    "pt": {"bins": pt_bins},
    "eta": {"bins": 25, "min": -2.5, "max": 2.5},
    "n": {"bins": 21, "min": -0.5, "max": 20.5}
}

histogram2d_config = {
    "pt": {"bins": pt_bins}
}

histogram_categories = [
    "pt_all_gen", "pt_matched_gen", "eta_all_gen", "eta_matched_gen", "pt_all_filterobj", "eta_all_filterobj",
    "pt_pass_gen", "eta_pass_gen", "pt_leading_all_gen", "pt_subleading_all_gen", "eta_leading_all_gen", "eta_subleading_all_gen",
    "pt_leading_pass_gen", "pt_subleading_pass_gen", "eta_leading_pass_gen", "eta_subleading_pass_gen",
    "pt_matched_filterobj", "eta_matched_filterobj", "pt_leading_fake_filterobj", "pt_subleading_fake_filterobj", "eta_leading_fake_filterobj", "eta_subleading_fake_filterobj",
    "pt_leading_matched_gen", "pt_subleading_matched_gen", "eta_leading_matched_gen", "eta_subleading_matched_gen",
    "pt_leading_l1matched_gen", "pt_subleading_l1matched_gen", "eta_leading_l1matched_gen", "eta_subleading_l1matched_gen",
    "pt_leading_matched_filterobj", "pt_subleading_matched_filterobj", "eta_leading_matched_filterobj", "eta_subleading_matched_filterobj",
    "pt_fake_filterobj", "eta_fake_filterobj", "pt_leading_all_filterobj", "pt_subleading_all_filterobj", "eta_leading_all_filterobj", "eta_subleading_all_filterobj",
    "pt_leading_barrel_matched_filterobj", "pt_leading_endcap_matched_filterobj", "eta_subleading_barrel_matched_filterobj", "eta_subleading_endcap_matched_filterobj",
    "eta_leading_barrel_matched_filterobj", "eta_leading_endcap_matched_filterobj", "pt_subleading_barrel_matched_filterobj", "pt_subleading_endcap_matched_filterobj",
    "pt_leading_barrel_all_gen", "pt_leading_endcap_all_gen", "eta_subleading_barrel_all_gen", "eta_subleading_endcap_all_gen",
    "eta_leading_barrel_all_gen", "eta_leading_endcap_all_gen", "pt_subleading_barrel_all_gen", "pt_subleading_endcap_all_gen",
    "pt_leading_filterobj", "pt_subleading_filterobj", "eta_leading_filterobj", "eta_subleading_filterobj",
    "pt_leading_barrel_l1matched_filterobj", "pt_leading_endcap_l1matched_filterobj", "eta_subleading_barrel_l1matched_filterobj", "eta_subleading_endcap_l1matched_filterobj",
    "eta_leading_barrel_l1matched_filterobj", "eta_leading_endcap_l1matched_filterobj", "pt_subleading_barrel_l1matched_filterobj", "pt_subleading_endcap_l1matched_filterobj",
    "n_pfjets", "pt_all_pfjet", "eta_all_pfjet", "pt_leading_pfjet", "eta_leading_pfjet", "pt_subleading_pfjet", "eta_subleading_pfjet",
    "pt_pfjet_matched", "eta_pfjet_matched",

    # NEW: gen-level, requiring BOTH gen taus to be matched
    "pt_leading_matched_gen_both",
    "eta_leading_matched_gen_both",
    "pt_subleading_matched_gen_both",
    "eta_subleading_matched_gen_both",
]

# Add new histogram categories for 2D histograms
histogram_categories_2d = [
    "pt_gen_vs_reco_matched_leading",
    "pt_gen_vs_reco_matched_subleading",
    "pt_gen_vs_reco_all_leading",
    "pt_gen_vs_reco_all_subleading"
]

# Define the maximum number of reconstructed taus to consider
max_reco_taus = 5  # Adjust this number as needed

histograms = {}

# Initialize histograms
for path, _ in trigger_paths:
    histograms[path] = {}
    # Initialize 1D histograms
    for category in histogram_categories:
        if category == "n_pfjets":
            histograms[path][category] = ROOT.TH1D(
                f"{category}_{path}",
                f"{category.replace('_', ' ').title()} for {path}",
                histogram_config["n"]["bins"], histogram_config["n"]["min"], histogram_config["n"]["max"]
            )
        elif category.startswith("pt") and "pfjet" in category:
            histograms[path][category] = ROOT.TH1D(
                f"{category}_{path}",
                f"{category.replace('_', ' ').title()} for {path}",
                len(pt_bins)-1, pt_bins
            )
        elif category.startswith("eta") and "pfjet" in category:
            histograms[path][category] = ROOT.TH1D(
                f"{category}_{path}",
                f"{category.replace('_', ' ').title()} for {path}",
                histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"]
            )
        elif "pt" in category:
            histograms[path][category] = ROOT.TH1D(
                f"{category}_{path}",
                f"{category.replace('_', ' ').title()} for {path}",
                len(pt_bins)-1,
                pt_bins
            )
        elif "eta" in category:
            histograms[path][category] = ROOT.TH1D(
                f"{category}_{path}", f"{category.replace('_', ' ').title()} for {path}",
                histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"]
            )
    # Initialize 2D histograms
    for category in histogram_categories_2d:
        histograms[path][category] = ROOT.TH2D(
            f"{category}_{path}",
            f"{category.replace('_', ' ').title()} for {path}",
            100, 0, pt_bins[-1],   # x-axis: 100 bins from 0 to pt_bins[-1]
            100, 0, pt_bins[-1]    # y-axis: 100 bins from 0 to pt_bins[-1]
        )
    histograms[path]["pt_reco_tau_overflow_barrel"] = ROOT.TH1D(
    f"pt_reco_tau_overflow_barrel_{path}",
    f"Pt of reconstructed taus beyond tau{max_reco_taus-1} for {path} in barrel",
    len(pt_bins)-1,
    pt_bins
    )
    histograms[path]["pt_reco_tau_overflow_endcap"] = ROOT.TH1D(
        f"pt_reco_tau_overflow_endcap_{path}",
        f"Pt of reconstructed taus beyond tau{max_reco_taus-1} for {path} in endcap",
        len(pt_bins)-1,
        pt_bins
    )
    for i in range(max_reco_taus):
        histograms[path][f"pt_reco_tau{i}_barrel"] = ROOT.TH1D(
            f"pt_reco_tau{i}_barrel_{path}",
            f"Pt of reconstructed tau{i} for {path} in barrel",
            len(pt_bins)-1,
            pt_bins
        )
        histograms[path][f"pt_reco_tau{i}_endcap"] = ROOT.TH1D(
            f"pt_reco_tau{i}_endcap_{path}",
            f"Pt of reconstructed tau{i} for {path} in endcap",
            len(pt_bins)-1,
            pt_bins
        )
    histograms[path]["eta_reco_tau_overflow"] = ROOT.TH1D(
        f"eta_reco_tau_overflow_{path}",
        f"Eta of reconstructed taus beyond tau{max_reco_taus-1} for {path}",
        histogram_config["eta"]["bins"],
        histogram_config["eta"]["min"],
        histogram_config["eta"]["max"]
    )
    for i in range(max_reco_taus):
        histograms[path][f"eta_reco_tau{i}"] = ROOT.TH1D(
            f"eta_reco_tau{i}_{path}",
            f"Eta of reconstructed tau{i} for {path}",
            histogram_config["eta"]["bins"],
            histogram_config["eta"]["min"],
            histogram_config["eta"]["max"]
        )
    histograms[path]["num_reco_taus"] = ROOT.TH1D(
        f"num_reco_taus_{path}",
        f"Number of reconstructed taus for {path}",
        max_reco_taus + 1, 0, max_reco_taus + 1
    )
    # Add new PFJet matched histograms split for leading and subleading taus
    histograms[path]["pt_leading_pfjet_matched"] = ROOT.TH1D(
         f"pt_leading_pfjet_matched_{path}",
         f"Leading tau PFJet Matched Pt for {path}",
         len(pt_bins)-1, pt_bins)
    histograms[path]["eta_leading_pfjet_matched"] = ROOT.TH1D(
         f"eta_leading_pfjet_matched_{path}",
         f"Leading tau PFJet Matched Eta for {path}",
         histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"])
    histograms[path]["pt_subleading_pfjet_matched"] = ROOT.TH1D(
         f"pt_subleading_pfjet_matched_{path}",
         f"Subleading tau PFJet Matched Pt for {path}",
         len(pt_bins)-1, pt_bins)
    histograms[path]["eta_subleading_pfjet_matched"] = ROOT.TH1D(
         f"eta_subleading_pfjet_matched_{path}",
         f"Subleading tau PFJet Matched Eta for {path}",
         histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"])

histograms_single_cuts = {}
for path, _ in trigger_paths:
    histograms_single_cuts[path] = {}
    for category in histogram_categories:
        # Initialize 1D histograms
        if category == "n_pfjets":
            histograms_single_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withSingleCuts",
                f"{category.replace('_', ' ').title()} for {path} with one tau passing cuts",
                histogram_config["n"]["bins"], histogram_config["n"]["min"], histogram_config["n"]["max"]
            )
        elif category.startswith("pt") and "pfjet" in category:
            histograms_single_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withSingleCuts",
                f"{category.replace('_', ' ').title()} for {path} with one tau passing cuts",
                len(pt_bins)-1, pt_bins
            )
        elif category.startswith("eta") and "pfjet" in category:
            histograms_single_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withSingleCuts",
                f"{category.replace('_', ' ').title()} for {path} with one tau passing cuts",
                histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"]
            )
        elif "pt" in category:
            histograms_single_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withSingleCuts",
                f"{category.replace('_', ' ').title()} for {path} with one tau passing cuts",
                len(pt_bins)-1,
                pt_bins
            )
        elif "eta" in category:
            histograms_single_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withSingleCuts", f"{category.replace('_', ' ').title()} for {path} with one tau passing cuts",
                histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"]
            )
    # Initialize 2D histograms
    for category in histogram_categories_2d:
        histograms_single_cuts[path][category] = ROOT.TH2D(
            f"{category}_{path}_withSingleCuts",
            f"{category.replace('_', ' ').title()} for {path} with one tau passing cuts",
            100, 0, pt_bins[-1],
            100, 0, pt_bins[-1]
        )

histograms_all_cuts = {}
for path, _ in trigger_paths:
    histograms_all_cuts[path] = {}
    for category in histogram_categories:
        # Initialize 1D histograms
        if category == "n_pfjets":
            histograms_all_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withAllCuts",
                f"{category.replace('_', ' ').title()} for {path} with both taus passing pt>40, |eta|<2.1",
                histogram_config["n"]["bins"], histogram_config["n"]["min"], histogram_config["n"]["max"]
            )
        elif category.startswith("pt") and "pfjet" in category:
            histograms_all_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withAllCuts",
                f"{category.replace('_', ' ').title()} for {path} with both taus passing pt>40, |eta|<2.1",
                len(pt_bins)-1, pt_bins
            )
        elif category.startswith("eta") and "pfjet" in category:
            histograms_all_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withAllCuts",
                f"{category.replace('_', ' ').title()} for {path} with both taus passing pt>40, |eta|<2.1",
                histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"]
            )
        elif "pt" in category:
            histograms_all_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withAllCuts",
                f"{category.replace('_', ' ').title()} for {path} with both taus passing pt>40, |eta|<2.1",
                len(pt_bins)-1,
                pt_bins
            )
        elif "eta" in category:
            histograms_all_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withAllCuts", f"{category.replace('_', ' ').title()} for {path} with both taus passing pt>40, |eta|<2.1",
                histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"]
            )
    # Initialize 2D histograms
    for category in histogram_categories_2d:
        histograms_all_cuts[path][category] = ROOT.TH2D(
            f"{category}_{path}_withAllCuts",
            f"{category.replace('_', ' ').title()} for {path} with both taus passing pt>40, |eta|<2.1",
            100, 0, pt_bins[-1],
            100, 0, pt_bins[-1]
        )

# Initialize cutflow counters dynamically from trigger_paths so renaming
# or adding/removing paths does not require manual key updates.
cutflow = {
    'TotalEvents': 0,
    'TwoGenTaus': 0,
    'BothTausPassCuts': 0,
    'filterobj': 0,
}

for path_name, _ in trigger_paths:
    cutflow[path_name] = 0
    cutflow[f'sub_matched_{path_name}'] = 0
    cutflow[f'lead_matched_to_sub_{path_name}'] = 0
    cutflow[f'sub_matched_to_lead_{path_name}'] = 0

products = []
add_product(products, "trig_sum", "trigger::TriggerEvent", "hltTriggerSummaryAOD")
add_product(products, "trig_res", "edm::TriggerResults", "TriggerResults")
add_product(products, "genparts", "std::vector<reco::GenParticle>", "genParticles")
add_product(products, "taugenjets", "std::vector<reco::GenJet>", "tauGenJets")
add_product(products, "genVisTaus", "std::vector<reco::GenParticle>", "genVisTaus")
add_product(products, "pfjets", "std::vector<reco::PFJet>", "hltPFJetForBtag")

evtdata = EvtData(products, verbose=False)
trigger_names = None
hadron_counter = 0
deltaR2_threshold = 0.1

def deltaR2(obj1, obj2):
    return ROOT.reco.deltaR2(obj1.eta(), obj1.phi(), obj2.eta(), obj2.phi())


for event_nr, event in enumerate(events):
    matched_objects = [[None, None] for _ in trigger_paths]

    cutflow['TotalEvents'] += 1
    if event_nr % 100 == 0:
        print(f"Processing event {event_nr}")

    evtdata.get_handles(event)
    trigger_names = event.object().triggerNames(evtdata.get("trig_res"))
    trigger_event = evtdata.get("trig_sum")


    # Retrieve PFJets safely: default to empty list if None
    pf_jets = evtdata.get("pfjets") or []
    sorted_pfjets = sorted(pf_jets, key=lambda jet: jet.pt(), reverse=True)
    # Compute matched PFJets using a ΔR threshold (e.g. 0.4→ΔR² < 0.16)
    matched_pfjets = []
    """for jet in sorted_pfjets:
        for tau in evtdata.get("genVisTaus"):
            if ROOT.reco.deltaR2(jet, tau) < 0.16:  # 0.4² = 0.16
                matched_pfjets.append(jet)
                break"""

    # Select hadronic taus from genVisTaus
    hadronic_taus = []
    """for idx, particle in enumerate(evtdata.get("genparts")):
        if abs(particle.pdgId()) == 15:
            hadronic_taus.append(particle)
            hadron_counter += 1"""
    genVisTaus = evtdata.get("genVisTaus")
    for idx, particle in enumerate(genVisTaus):
        hadronic_taus.append(particle)
        hadron_counter += 1

    # Warn if more than two gen visible taus are present

    if len(hadronic_taus) >= 2:
        cutflow['TwoGenTaus'] += 1
    else:
        continue  # Skip events with fewer than 2 hadronic taus

    lead_filterobj = [None for _ in trigger_paths]
    sub_filterobj = [None for _ in trigger_paths]

    hadronic_taus.sort(key=lambda x: x.pt(), reverse=True)
    leading_tau = hadronic_taus[0] if len(hadronic_taus) > 0 else None
    subleading_tau = hadronic_taus[1] if len(hadronic_taus) > 1 else None

    # NEW: Fill pt_all_gen and eta_all_gen histograms with both taus (guarded)
    for path, _ in trigger_paths:
        if leading_tau is not None:
            histograms[path]["pt_all_gen"].Fill(leading_tau.pt())
            histograms[path]["eta_all_gen"].Fill(leading_tau.eta())
        if subleading_tau is not None:
            histograms[path]["pt_all_gen"].Fill(subleading_tau.pt())
            histograms[path]["eta_all_gen"].Fill(subleading_tau.eta())

    # Check if both taus pass the pt and eta cuts (guarded)
    lead_pass_pt  = (leading_tau is not None) and (leading_tau.pt() > 40)
    lead_pass_eta = (leading_tau is not None) and (abs(leading_tau.eta()) < 2.1)
    sub_pass_pt   = (subleading_tau is not None) and (subleading_tau.pt() > 40)
    sub_pass_eta  = (subleading_tau is not None) and (abs(subleading_tau.eta()) < 2.1)

    lead_cut_both = lead_pass_pt and lead_pass_eta
    sub_cut_both  = sub_pass_pt and sub_pass_eta
    both_taus_pass_cuts = lead_cut_both and sub_cut_both

    if both_taus_pass_cuts:
        cutflow['BothTausPassCuts'] += 1

    for path_index, (path, filter_name) in enumerate(trigger_paths):
        # Retrieve the trigger index
        trigger_index = None
        for i in range(trigger_names.size()):
            if path in str(trigger_names.triggerName(i)):
                trigger_index = i
                #print(f"Found trigger path {path} at index {trigger_index}")
                break

        if trigger_index is None:
            print(f"Warning: No loose match for trigger path {path} in event {event_nr}")
            continue
            
        histograms[path]["pt_leading_all_gen"].Fill(leading_tau.pt())
        histograms[path]["eta_leading_all_gen"].Fill(leading_tau.eta())
        histograms[path]["pt_subleading_all_gen"].Fill(subleading_tau.pt())
        histograms[path]["eta_subleading_all_gen"].Fill(subleading_tau.eta())
        if lead_pass_pt:
            histograms_single_cuts[path]["eta_leading_all_gen"].Fill(leading_tau.eta())
        if lead_pass_eta:
            histograms_single_cuts[path]["pt_leading_all_gen"].Fill(leading_tau.pt())
        if sub_pass_pt:
            histograms_single_cuts[path]["eta_subleading_all_gen"].Fill(subleading_tau.eta())
        if sub_pass_eta:
            histograms_single_cuts[path]["pt_subleading_all_gen"].Fill(subleading_tau.pt())
        if lead_cut_both and sub_pass_pt:
            histograms_all_cuts[path]["eta_subleading_all_gen"].Fill(subleading_tau.eta())
            if abs(subleading_tau.eta()) < 1.5:
                histograms_all_cuts[path]["eta_subleading_barrel_all_gen"].Fill(subleading_tau.eta())
            else:
                histograms_all_cuts[path]["eta_subleading_endcap_all_gen"].Fill(subleading_tau.eta())
        if lead_cut_both and sub_pass_eta:
            histograms_all_cuts[path]["pt_subleading_all_gen"].Fill(subleading_tau.pt())
            if abs(subleading_tau.pt()) < 1.5:
                histograms_all_cuts[path]["pt_subleading_barrel_all_gen"].Fill(subleading_tau.pt())
            else:
                histograms_all_cuts[path]["pt_subleading_endcap_all_gen"].Fill(subleading_tau.pt())
        if sub_cut_both and lead_pass_pt:
            histograms_all_cuts[path]["eta_leading_all_gen"].Fill(leading_tau.eta())
            if abs(leading_tau.eta()) < 1.5:
                histograms_all_cuts[path]["eta_leading_barrel_all_gen"].Fill(leading_tau.eta())
            else:
                histograms_all_cuts[path]["eta_leading_endcap_all_gen"].Fill(leading_tau.eta())
        if sub_cut_both and lead_pass_eta:
            histograms_all_cuts[path]["pt_leading_all_gen"].Fill(leading_tau.pt())
            if abs(leading_tau.pt()) < 1.5:
                histograms_all_cuts[path]["pt_leading_barrel_all_gen"].Fill(leading_tau.pt())
            else:
                histograms_all_cuts[path]["pt_leading_endcap_all_gen"].Fill(leading_tau.pt())
        if abs(leading_tau.eta()) < 1.5:
            histograms[path]["pt_leading_barrel_all_gen"].Fill(leading_tau.pt())
            histograms[path]["eta_leading_barrel_all_gen"].Fill(leading_tau.eta())
        else:
            histograms[path]["pt_leading_endcap_all_gen"].Fill(leading_tau.pt())
            histograms[path]["eta_leading_endcap_all_gen"].Fill(leading_tau.eta())
        if abs(subleading_tau.eta()) < 1.5:
            histograms[path]["pt_subleading_barrel_all_gen"].Fill(subleading_tau.pt())
            histograms[path]["eta_subleading_barrel_all_gen"].Fill(subleading_tau.eta())
        else:
            histograms[path]["pt_subleading_endcap_all_gen"].Fill(subleading_tau.pt())
            histograms[path]["eta_subleading_endcap_all_gen"].Fill(subleading_tau.eta())

        
        # Fill PFJet histograms using the same cuts as for taus: pt > 40 and |eta| < 2.1
        histograms[path]["n_pfjets"].Fill(len(sorted_pfjets))
        for i, jet in enumerate(sorted_pfjets):
            if jet.pt() > 40 and abs(jet.eta()) < 2.1:  # apply tau-like cuts on jets
                histograms[path]["pt_all_pfjet"].Fill(jet.pt())
                histograms[path]["eta_all_pfjet"].Fill(jet.eta())
                if i == 0:
                    histograms[path]["pt_leading_pfjet"].Fill(jet.pt())
                    histograms[path]["eta_leading_pfjet"].Fill(jet.eta())
                elif i == 1:
                    histograms[path]["pt_subleading_pfjet"].Fill(jet.pt())
                    histograms[path]["eta_subleading_pfjet"].Fill(jet.eta())

        # Instead of looping over all matched PFJets, fill only once per tau:
        if leading_tau is not None:
            if any(ROOT.reco.deltaR2(jet, leading_tau) < 0.16 for jet in sorted_pfjets):
                histograms[path]["pt_leading_pfjet_matched"].Fill(leading_tau.pt())
                histograms[path]["eta_leading_pfjet_matched"].Fill(leading_tau.eta())
        if subleading_tau is not None:
            if any(ROOT.reco.deltaR2(jet, subleading_tau) < 0.16 for jet in sorted_pfjets):
                histograms[path]["pt_subleading_pfjet_matched"].Fill(subleading_tau.pt())
                histograms[path]["eta_subleading_pfjet_matched"].Fill(subleading_tau.eta())

        trigger_accept = evtdata.get("trig_res").accept(trigger_index)

        if trigger_accept:
            cutflow[path] += 1
            histograms[path]["pt_leading_pass_gen"].Fill(leading_tau.pt())
            histograms[path]["eta_leading_pass_gen"].Fill(leading_tau.eta())
            histograms[path]["pt_subleading_pass_gen"].Fill(subleading_tau.pt())
            histograms[path]["eta_subleading_pass_gen"].Fill(subleading_tau.eta())
            if lead_pass_eta:
                histograms_single_cuts[path]["pt_leading_pass_gen"].Fill(leading_tau.pt())
            if lead_pass_pt:
                histograms_single_cuts[path]["eta_leading_pass_gen"].Fill(leading_tau.eta())
            if sub_pass_eta:
                histograms_single_cuts[path]["pt_subleading_pass_gen"].Fill(subleading_tau.pt())
            if sub_pass_pt:
                histograms_single_cuts[path]["eta_subleading_pass_gen"].Fill(subleading_tau.eta())
            if lead_pass_pt and lead_pass_eta and sub_pass_pt:
                histograms_all_cuts[path]["eta_subleading_pass_gen"].Fill(subleading_tau.eta())
            if lead_pass_pt and lead_pass_eta and sub_pass_eta:
                histograms_all_cuts[path]["pt_subleading_pass_gen"].Fill(subleading_tau.pt())
            if lead_pass_pt and sub_pass_pt and sub_pass_eta:
                histograms_all_cuts[path]["eta_leading_pass_gen"].Fill(leading_tau.eta())
            if lead_pass_eta and sub_pass_pt and sub_pass_eta:
                histograms_all_cuts[path]["pt_leading_pass_gen"].Fill(leading_tau.pt())


        # Construct the InputTag with filter name and process name
        filter_tag = ROOT.edm.InputTag(filter_name, "", "HLTX")
        # Get index of the filter using filterIndex method with InputTag
        filter_index = trigger_event.filterIndex(filter_tag)
        if filter_index == trigger_event.sizeFilters():
            continue  # Skip to the next path if filter is not there
        #print(f"Filter {filter_name} found at index {filter_index}")
        # Get trigger objects passing the filter
        filter_keys = trigger_event.filterKeys(filter_index)
        all_trigger_objects = trigger_event.getObjects()
        selected_objects = []
        for key in filter_keys:
            obj = all_trigger_objects[key]
            selected_objects.append(obj)
        lead_matched = False
        sub_matched = False
        selected_objects.sort(key=lambda obj: obj.pt(), reverse=True)
        obj_lead = None
        obj_sub = None
        unmatched_objs = []
        matched_trigger_objects = set()
        #print (f"number of taus in filter {filter_name}: {len(selected_objects)}")
        #if filter_name == "hltL1P2GTTau":
        #    print(f"filter objects found: {len(selected_objects)}")

        for i, tau in enumerate(selected_objects):
            if i < max_reco_taus:
                if abs(tau.eta()) < 1.5:
                    histograms[path][f"pt_reco_tau{i}_barrel"].Fill(tau.pt())
                else:
                    histograms[path][f"pt_reco_tau{i}_endcap"].Fill(tau.pt())
                histograms[path][f"eta_reco_tau{i}"].Fill(tau.eta())
            else:
                if abs(tau.eta()) < 1.5:
                    histograms[path]["pt_reco_tau_overflow_barrel"].Fill(tau.pt())
                else:
                    histograms[path]["pt_reco_tau_overflow_endcap"].Fill(tau.pt())
                histograms[path]["eta_reco_tau_overflow"].Fill(tau.eta())
        histograms[path]["num_reco_taus"].Fill(len(selected_objects))

        for obj in selected_objects:
            # Fill histograms for all filter objects (only when trigger fired)
            if trigger_accept:
                histograms[path]["pt_all_filterobj"].Fill(obj.pt())
                histograms[path]["eta_all_filterobj"].Fill(obj.eta())
                if obj.pt() > 40:
                    histograms_single_cuts[path]["eta_all_filterobj"].Fill(obj.eta())
                if abs(obj.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_all_filterobj"].Fill(obj.pt())

            # Match to ANY genVisTau (not just the first two)
            matched_any_gen = any(
                ROOT.reco.deltaR2(tau.eta(), tau.phi(), obj.eta(), obj.phi()) < 0.01
                for tau in hadronic_taus
            )
            if matched_any_gen:
                # Optional: fill generic matched filter-object histos
                histograms[path]["pt_matched_filterobj"].Fill(obj.pt())
                histograms[path]["eta_matched_filterobj"].Fill(obj.eta())

            # Keep leading/subleading-specific matching for dedicated histos (guarded)
            if (leading_tau is not None and
                ROOT.reco.deltaR2(leading_tau.eta(), leading_tau.phi(), obj.eta(), obj.phi()) < 0.01 and not lead_matched):
                lead_matched = True
                histograms[path]["pt_leading_matched_filterobj"].Fill(obj.pt())
                histograms[path]["eta_leading_matched_filterobj"].Fill(obj.eta())
                if obj.pt() > 40:
                    histograms_single_cuts[path]["eta_leading_matched_filterobj"].Fill(obj.eta())
                if abs(obj.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_leading_matched_filterobj"].Fill(obj.pt())
                if sub_cut_both and obj.pt() > 40:
                    histograms_all_cuts[path]["eta_leading_matched_filterobj"].Fill(obj.eta())
                if sub_cut_both and abs(obj.eta()) < 2.1:
                    histograms_all_cuts[path]["pt_leading_matched_filterobj"].Fill(obj.pt())
                obj_lead = obj
                lead_filterobj[path_index] = obj_lead
                matched_objects[path_index][0] = obj_lead

            elif (subleading_tau is not None and
                  ROOT.reco.deltaR2(subleading_tau.eta(), subleading_tau.phi(), obj.eta(), obj.phi()) < 0.01 and not sub_matched):
                sub_matched = True
                histograms[path]["pt_subleading_matched_filterobj"].Fill(obj.pt())
                histograms[path]["eta_subleading_matched_filterobj"].Fill(obj.eta())
                if obj.pt() > 40:
                    histograms_single_cuts[path]["eta_subleading_matched_filterobj"].Fill(obj.eta())    
                if abs(obj.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_subleading_matched_filterobj"].Fill(obj.pt())
                if lead_cut_both and obj.pt() > 40:
                    histograms_all_cuts[path]["eta_subleading_matched_filterobj"].Fill(obj.eta())
                if lead_cut_both and abs(obj.eta()) < 2.1:
                    histograms_all_cuts[path]["pt_subleading_matched_filterobj"].Fill(obj.pt())
                obj_sub = obj
                sub_filterobj[path_index] = obj_sub
                matched_objects[path_index][1] = obj_sub
                cutflow[f'sub_matched_{path}'] += 1

            # Only count as fake if it did NOT match ANY genVisTau
            if not matched_any_gen:
                unmatched_objs.append(obj)

        """obj_lead = selected_objects[0] if len(selected_objects) > 0 else None
        obj_sub = selected_objects[1] if len(selected_objects) > 1 else None

        if obj_lead is not None:
            histograms[path]["pt_leading_all_filterobj"].Fill(obj_lead.pt())
            histograms[path]["eta_leading_all_filterobj"].Fill(obj_lead.eta())
            if obj_lead.pt() > 40:
                histograms_single_cuts[path]["eta_leading_all_filterobj"].Fill(obj_lead.eta())
            if abs(obj_lead.eta()) < 2.1:
                histograms_single_cuts[path]["pt_leading_all_filterobj"].Fill(obj_lead.pt())
            if ROOT.reco.deltaR2(leading_tau, obj_lead) < deltaR2_threshold:
                lead_matched = True
                matched_trigger_objects.add(obj_lead)
                lead_filterobj[path_index] = obj_lead
                histograms[path]["pt_leading_matched_filterobj"].Fill(obj_lead.pt())
                histograms[path]["eta_leading_matched_filterobj"].Fill(obj_lead.eta())
                if obj_lead.pt() > 40:
                    histograms_single_cuts[path]["eta_leading_matched_filterobj"].Fill(obj_lead.eta())
                if abs(obj_lead.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_leading_matched_filterobj"].Fill(obj_lead.pt())
                if sub_cut_both and obj_lead.pt() > 40:
                    histograms_all_cuts[path]["eta_leading_matched_filterobj"].Fill(obj_lead.eta())
                if sub_cut_both and abs(obj_lead.eta()) < 2.1:
                    histograms_all_cuts[path]["pt_leading_matched_filterobj"].Fill(obj_lead.pt())
            else: 
                unmatched_objs.append(obj_lead)
        if obj_sub is not None:
            histograms[path]["pt_subleading_all_filterobj"].Fill(obj_sub.pt())
            histograms[path]["eta_subleading_all_filterobj"].Fill(obj_sub.eta())
            if obj_sub.pt() > 40:
                histograms_single_cuts[path]["eta_subleading_all_filterobj"].Fill(obj_sub.eta())
            if abs(obj_sub.eta()) < 2.1:
                histograms_single_cuts[path]["pt_subleading_all_filterobj"].Fill(obj_sub.pt())
            if ROOT.reco.deltaR2(subleading_tau, obj_sub) < deltaR2_threshold:
                sub_matched = True
                matched_trigger_objects.add(obj_sub)
                sub_filterobj[path_index] = obj_sub
                histograms[path]["pt_subleading_matched_filterobj"].Fill(obj_sub.pt())
                histograms[path]["eta_subleading_matched_filterobj"].Fill(obj_sub.eta())
                if obj_sub.pt() > 40:
                    histograms_single_cuts[path]["eta_subleading_matched_filterobj"].Fill(obj_sub.eta())
                if abs(obj_sub.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_subleading_matched_filterobj"].Fill(obj_sub.pt())
                if lead_cut_both and obj_sub.pt() > 40:
                    histograms_all_cuts[path]["eta_subleading_matched_filterobj"].Fill(obj_sub.eta())
                if lead_cut_both and abs(obj_sub.eta()) < 2.1:
                    histograms_all_cuts[path]["pt_subleading_matched_filterobj"].Fill(obj_sub.pt())
                cutflow[f'sub_matched_{path}'] += 1
            else:
                unmatched_objs.append(obj_sub)
        if len(selected_objects) > 2:
            unmatched_objs.append(selected_objects[2])
        # Matching leading tau
        for obj in selected_objects:
            if obj in matched_trigger_objects:
                continue  # Skip if the trigger object is already matched
            if deltaR2(leading_tau, obj) < 0.1:
                lead_matched = True
                obj_lead = obj
                matched_trigger_objects.add(obj)
                lead_filterobj[path_index] = obj_lead
                # Fill histograms for leading matched object
                histograms[path]["pt_leading_matched_filterobj"].Fill(obj_lead.pt())
                histograms[path]["eta_leading_matched_filterobj"].Fill(obj_lead.eta())
                if obj_lead.pt() > 40:
                    histograms_single_cuts[path]["eta_leading_matched_filterobj"].Fill(obj_lead.eta())
                if abs(obj_lead.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_leading_matched_filterobj"].Fill(obj_lead.pt())
                break  # Stop after finding the closest match

        # Matching subleading tau
        for obj in selected_objects:
            if obj in matched_trigger_objects:
                continue  # Skip if the trigger object is already matched
            if deltaR2(subleading_tau, obj) < 0.1:
                sub_matched = True
                obj_sub = obj
                matched_trigger_objects.add(obj)
                sub_filterobj[path_index] = obj_sub
                # Fill histograms for subleading matched object
                histograms[path]["pt_subleading_matched_filterobj"].Fill(obj_sub.pt())
                histograms[path]["eta_subleading_matched_filterobj"].Fill(obj_sub.eta())
                if obj_sub.pt() > 40:
                    histograms_single_cuts[path]["eta_subleading_matched_filterobj"].Fill(obj_sub.eta())
                if abs(obj_sub.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_subleading_matched_filterobj"].Fill(obj_sub.pt())
                break  # Stop after finding the closest match"""

        # Collect unmatched trigger objects
        # Fill histograms for unmatched (fake) trigger objects ONLY if the trigger fired
        if trigger_accept:
            for fake in unmatched_objs:
                histograms[path]["pt_fake_filterobj"].Fill(fake.pt())
                histograms[path]["eta_fake_filterobj"].Fill(fake.eta())
                if fake.pt() > 40:
                    histograms_single_cuts[path]["eta_fake_filterobj"].Fill(fake.eta())
                if abs(fake.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_fake_filterobj"].Fill(fake.pt())

            unmatched_objs.sort(key=lambda obj: obj.pt(), reverse=True)
            if len(unmatched_objs) > 0:
                histograms[path]["pt_leading_fake_filterobj"].Fill(unmatched_objs[0].pt())
                histograms[path]["eta_leading_fake_filterobj"].Fill(unmatched_objs[0].eta())
                if unmatched_objs[0].pt() > 40:
                    histograms_single_cuts[path]["eta_leading_fake_filterobj"].Fill(unmatched_objs[0].eta())
                if abs(unmatched_objs[0].eta()) < 2.1:
                    histograms_single_cuts[path]["pt_leading_fake_filterobj"].Fill(unmatched_objs[0].pt())

            if len(unmatched_objs) > 1:
                cutflow['filterobj'] += 1
                histograms[path]["pt_subleading_fake_filterobj"].Fill(unmatched_objs[1].pt())
                histograms[path]["eta_subleading_fake_filterobj"].Fill(unmatched_objs[1].eta())
                if unmatched_objs[1].pt() > 40:
                    histograms_single_cuts[path]["eta_subleading_fake_filterobj"].Fill(unmatched_objs[1].eta())
                if abs(unmatched_objs[1].eta()) < 2.1:
                    histograms_single_cuts[path]["pt_subleading_fake_filterobj"].Fill(unmatched_objs[1].pt())
        # ...existing code...

        # Fill histograms for all filter objects if trigger is accepted
        if trigger_accept:
            if obj_lead is not None:
                histograms[path]["pt_leading_all_filterobj"].Fill(obj_lead.pt())
                histograms[path]["eta_leading_all_filterobj"].Fill(obj_lead.eta())
                if obj_lead.pt() > 40:
                    histograms_single_cuts[path]["eta_leading_all_filterobj"].Fill(obj_lead.eta())
                if abs(obj_lead.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_leading_all_filterobj"].Fill(obj_lead.pt())
            if obj_sub is not None:
                histograms[path]["pt_subleading_all_filterobj"].Fill(obj_sub.pt())
                histograms[path]["eta_subleading_all_filterobj"].Fill(obj_sub.eta())
                if obj_sub.pt() > 40:
                    histograms_single_cuts[path]["eta_subleading_all_filterobj"].Fill(obj_sub.eta())
                if abs(obj_sub.eta()) < 2.1:
                    histograms_single_cuts[path]["pt_subleading_all_filterobj"].Fill(obj_sub.pt())

        # Additional cuts and histogram filling
        if obj_lead is not None and obj_sub is not None:
            if obj_lead.pt() > 40 and obj_sub.pt() > 40:
                histograms_all_cuts[path]["eta_leading_matched_filterobj"].Fill(obj_lead.eta())
                histograms_all_cuts[path]["eta_subleading_matched_filterobj"].Fill(obj_sub.eta())
            if abs(obj_lead.eta()) < 2.1 and abs(obj_sub.eta()) < 2.1:
                histograms_all_cuts[path]["pt_leading_matched_filterobj"].Fill(obj_lead.pt())
                histograms_all_cuts[path]["pt_subleading_matched_filterobj"].Fill(obj_sub.pt())
            # Update cutflow counts (guarded)
            if (leading_tau is not None) and (deltaR2(leading_tau, obj_sub) < 0.1):
                cutflow[f'lead_matched_to_sub_{path}'] += 1
            if (subleading_tau is not None) and (deltaR2(subleading_tau, obj_lead) < 0.1):
                cutflow[f'sub_matched_to_lead_{path}'] += 1

        # Fill 2D histograms for leading matched tau (guarded)
        if lead_matched and obj_lead is not None and leading_tau is not None:
            histograms[path]["pt_gen_vs_reco_matched_leading"].Fill(leading_tau.pt(), obj_lead.pt())
            if lead_pass_pt or lead_pass_eta:
                histograms_single_cuts[path]["pt_gen_vs_reco_matched_leading"].Fill(leading_tau.pt(), obj_lead.pt())
            if both_taus_pass_cuts:
                histograms_all_cuts[path]["pt_gen_vs_reco_matched_leading"].Fill(leading_tau.pt(), obj_lead.pt())

        # Fill 2D histograms for subleading matched tau (guarded)
        if sub_matched and obj_sub is not None and subleading_tau is not None:
            histograms[path]["pt_gen_vs_reco_matched_subleading"].Fill(subleading_tau.pt(), obj_sub.pt())
            if sub_pass_pt or sub_pass_eta:
                histograms_single_cuts[path]["pt_gen_vs_reco_matched_subleading"].Fill(subleading_tau.pt(), obj_sub.pt())
            if both_taus_pass_cuts:
                histograms_all_cuts[path]["pt_gen_vs_reco_matched_subleading"].Fill(subleading_tau.pt(), obj_sub.pt())

        # Fill histograms for matched gen taus (guarded)
        if lead_matched and trigger_accept and (leading_tau is not None):
            histograms[path]["pt_leading_matched_gen"].Fill(leading_tau.pt())
            histograms[path]["eta_leading_matched_gen"].Fill(leading_tau.eta())
            if obj_lead is not None:
                if abs(obj_lead.eta()) < 1.5:
                    histograms[path]["pt_leading_barrel_matched_filterobj"].Fill(obj_lead.pt())
                    histograms[path]["eta_leading_barrel_matched_filterobj"].Fill(obj_lead.eta())
                else:
                    histograms[path]["pt_leading_endcap_matched_filterobj"].Fill(obj_lead.pt())
                    histograms[path]["eta_leading_endcap_matched_filterobj"].Fill(obj_lead.eta())

        if lead_matched and sub_matched and trigger_accept and (subleading_tau is not None):
            histograms[path]["pt_subleading_matched_gen"].Fill(subleading_tau.pt())
            histograms[path]["eta_subleading_matched_gen"].Fill(subleading_tau.eta())
            if obj_sub is not None:
                if abs(obj_sub.eta()) < 1.5:
                    histograms[path]["pt_subleading_barrel_matched_filterobj"].Fill(obj_sub.pt())
                    histograms[path]["eta_subleading_barrel_matched_filterobj"].Fill(obj_sub.eta())
                else:
                    histograms[path]["pt_subleading_endcap_matched_filterobj"].Fill(obj_sub.pt())
                    histograms[path]["eta_subleading_endcap_matched_filterobj"].Fill(obj_sub.eta())
            
        if matched_objects[0][0] is not None and matched_objects[path_index][0] is not None:
            if trigger_accept:  # Add trigger accept check
                for l1_tau in matched_objects[0]:
                    if l1_tau is not None:
                        if ROOT.reco.deltaR2(l1_tau.eta(), l1_tau.phi(), matched_objects[path_index][0].eta(), matched_objects[path_index][0].phi()) < 0.1:
                            histograms[path]["pt_leading_l1matched_gen"].Fill(leading_tau.pt())
                            histograms[path]["eta_leading_l1matched_gen"].Fill(leading_tau.eta())
                            if obj_lead is not None:
                                if abs(obj_lead.eta()) < 1.5:
                                    histograms[path]["pt_leading_barrel_l1matched_filterobj"].Fill(obj_lead.pt())
                                    histograms[path]["eta_leading_barrel_l1matched_filterobj"].Fill(obj_lead.eta())
                                else:
                                    histograms[path]["pt_leading_endcap_l1matched_filterobj"].Fill(obj_lead.pt())
                                    histograms[path]["eta_leading_endcap_l1matched_filterobj"].Fill(obj_lead.eta())
                            if lead_pass_pt:
                                histograms_single_cuts[path]["eta_leading_l1matched_gen"].Fill(leading_tau.eta())
                            if lead_pass_eta:
                                histograms_single_cuts[path]["pt_leading_l1matched_gen"].Fill(leading_tau.pt())
                            if sub_cut_both and lead_pass_pt:  # Corrected condition order
                                histograms_all_cuts[path]["eta_leading_l1matched_gen"].Fill(leading_tau.eta())
                            if sub_cut_both and lead_pass_eta:  # Corrected condition order
                                histograms_all_cuts[path]["pt_leading_l1matched_gen"].Fill(leading_tau.pt())
                       
        if matched_objects[0][1] is not None and matched_objects[path_index][1] is not None:
            if trigger_accept:  # Add trigger accept check
                for l1_tau in matched_objects[0]:
                    if l1_tau is not None:
                        if ROOT.reco.deltaR2(l1_tau.eta(), l1_tau.phi(), matched_objects[path_index][1].eta(), matched_objects[path_index][1].phi()) < 0.1:
                            histograms[path]["pt_subleading_l1matched_gen"].Fill(subleading_tau.pt())
                            histograms[path]["eta_subleading_l1matched_gen"].Fill(subleading_tau.eta())
                            if obj_sub is not None:
                                if abs(obj_sub.eta()) < 1.5:
                                    histograms[path]["pt_subleading_barrel_l1matched_filterobj"].Fill(obj_sub.pt())
                                    histograms[path]["eta_subleading_barrel_l1matched_filterobj"].Fill(obj_sub.eta())
                                else:
                                    histograms[path]["pt_subleading_endcap_l1matched_filterobj"].Fill(obj_sub.pt())
                                    histograms[path]["eta_subleading_endcap_l1matched_filterobj"].Fill(obj_sub.eta())
                            if sub_pass_pt: 
                                histograms_single_cuts[path]["eta_subleading_l1matched_gen"].Fill(subleading_tau.eta())
                            if sub_pass_eta:
                                histograms_single_cuts[path]["pt_subleading_l1matched_gen"].Fill(subleading_tau.pt())
                            if lead_cut_both and sub_pass_pt:
                                histograms_all_cuts[path]["eta_subleading_l1matched_gen"].Fill(subleading_tau.eta())
                            if lead_cut_both and sub_pass_eta:    
                                histograms_all_cuts[path]["pt_subleading_l1matched_gen"].Fill(subleading_tau.pt())
                         
        if lead_matched and trigger_accept and lead_pass_pt:
            histograms_single_cuts[path]["eta_leading_matched_gen"].Fill(leading_tau.eta())
        if lead_matched and trigger_accept and lead_pass_eta:
            histograms_single_cuts[path]["pt_leading_matched_gen"].Fill(leading_tau.pt())
        if sub_matched and trigger_accept and sub_pass_pt:
            histograms_single_cuts[path]["eta_subleading_matched_gen"].Fill(subleading_tau.eta())
        if sub_matched and trigger_accept and sub_pass_eta:
            histograms_single_cuts[path]["pt_subleading_matched_gen"].Fill(subleading_tau.pt())

        if lead_matched and sub_matched and trigger_accept and lead_cut_both and sub_pass_pt:
            histograms_all_cuts[path]["eta_subleading_matched_gen"].Fill(subleading_tau.eta())
            if obj_sub is not None:
                if abs(obj_sub.eta()) < 1.5:
                    histograms_all_cuts[path]["eta_subleading_barrel_matched_filterobj"].Fill(obj_sub.eta())
                else:
                    histograms_all_cuts[path]["eta_subleading_endcap_matched_filterobj"].Fill(obj_sub.eta())
        if lead_matched and sub_matched and trigger_accept and lead_cut_both and sub_pass_eta:
            histograms_all_cuts[path]["pt_subleading_matched_gen"].Fill(subleading_tau.pt())
            if obj_sub is not None:
                if abs(obj_sub.eta()) < 1.5:
                    histograms_all_cuts[path]["pt_subleading_barrel_matched_filterobj"].Fill(obj_sub.pt())
                else:
                    histograms_all_cuts[path]["pt_subleading_endcap_matched_filterobj"].Fill(obj_sub.pt())
        if lead_matched and sub_matched and trigger_accept and sub_cut_both and lead_pass_pt:
            histograms_all_cuts[path]["eta_leading_matched_gen"].Fill(leading_tau.eta())
            if obj_lead is not None:
                if abs(obj_lead.eta()) < 1.5:
                    histograms_all_cuts[path]["eta_leading_barrel_matched_filterobj"].Fill(obj_lead.eta())
                else:
                    histograms_all_cuts[path]["eta_leading_endcap_matched_filterobj"].Fill(obj_lead.eta())
        if lead_matched and sub_matched and trigger_accept and sub_cut_both and lead_pass_eta:
            histograms_all_cuts[path]["pt_leading_matched_gen"].Fill(leading_tau.pt())
            if obj_lead is not None:
                if abs(obj_lead.eta()) < 1.5:
                    histograms_all_cuts[path]["pt_leading_barrel_matched_filterobj"].Fill(obj_lead.pt())
                else:
                    histograms_all_cuts[path]["pt_leading_endcap_matched_filterobj"].Fill(obj_lead.pt())

        # NEW: gen-only histos, requiring BOTH gen taus to be matched
        if lead_matched and sub_matched and trigger_accept:
            if leading_tau is not None:
                histograms[path]["pt_leading_matched_gen_both"].Fill(leading_tau.pt())
                histograms[path]["eta_leading_matched_gen_both"].Fill(leading_tau.eta())
            if subleading_tau is not None:
                histograms[path]["pt_subleading_matched_gen_both"].Fill(subleading_tau.pt())
                histograms[path]["eta_subleading_matched_gen_both"].Fill(subleading_tau.eta())

            if lead_pass_pt:
                histograms_single_cuts[path]["eta_leading_matched_gen_both"].Fill(leading_tau.eta())
            if lead_pass_eta:
                histograms_single_cuts[path]["pt_leading_matched_gen_both"].Fill(leading_tau.pt())
            if sub_pass_pt:
                histograms_single_cuts[path]["eta_subleading_matched_gen_both"].Fill(subleading_tau.eta())
            if sub_pass_eta:
                histograms_single_cuts[path]["pt_subleading_matched_gen_both"].Fill(subleading_tau.pt())

            if lead_pass_pt and sub_cut_both:
                histograms_all_cuts[path]["eta_leading_matched_gen_both"].Fill(leading_tau.eta())
            if lead_pass_eta and sub_cut_both:
                histograms_all_cuts[path]["pt_leading_matched_gen_both"].Fill(leading_tau.pt())
            if sub_pass_pt and lead_cut_both:
                histograms_all_cuts[path]["eta_subleading_matched_gen_both"].Fill(subleading_tau.eta())
            if sub_pass_eta and lead_cut_both:
                histograms_all_cuts[path]["pt_subleading_matched_gen_both"].Fill(subleading_tau.pt())


output_file = ROOT.TFile(output_filename, "RECREATE")

# Print cutflow summary
print("Cutflow Summary:")
for key, value in cutflow.items():
    print(f"{key}: {value}")

print(f"Histograms have been saved to {output_filename}.")
print(f"Total events with at least two hadronic taus: {hadron_counter}")
for path_hists in histograms.values():
    for hist in path_hists.values():
        hist.Write()
for path_hists in histograms_single_cuts.values():
    for hist in path_hists.values():
        hist.Write()
for path_hists in histograms_all_cuts.values():
    for hist in path_hists.values():
        hist.Write()
output_file.Close()
