from __future__ import absolute_import, division, print_function
import ROOT
from DataFormats.FWLite import Events
from modules.EvtData import EvtData, add_product
from modules.GenTools import load_fwlitelibs
from modules.tau_selection import is_hadronic_tau
load_fwlitelibs()

import glob

#files = glob.glob('/eos/user/a/agruber/samples/HLT_Upgrade_L1filter/VBF_HTauTau_14_1_0_pre7_TriggerOnly/Phase2_L1P2GT_HLT*.root')
files = glob.glob('/eos/user/a/agruber/samples/HLT_Upgrade_L1filter/VBF_HTauTau_14_2_0_pre2_TrigOnly/Phase2_L1P2GT_HLT*.root')
ROOT.gSystem.Load('libDataFormatsL1Trigger.so')
events = Events(files)
#events = Events("Phase2_L1P2GT_HLT.root")

# Define trigger paths and their corresponding filters
trigger_paths = [
    ["L1P2GT_DoubleNNTau52", "hltL1P2GTCandTau"],
    ["HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1", "hltHpsDoublePFTau40TrackPt1MediumChargedIsolation"],
    ["HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1", "hltHpsDoublePFTau35MediumDitauWPDeepTau"]
]

# Define the binning configuration in one place
histogram_config = {
    "pt": {"bins": 50, "min": 0, "max": 200},
    "eta": {"bins": 25, "min": -2.5, "max": 2.5}
}

histogram_categories = [
    "pt_all_gen", "pt_matched_gen", "eta_all_gen", "eta_matched_gen", "pt_all_filterobj", "eta_all_filterobj",
    "pt_pass_gen", "eta_pass_gen",
    "pt_matched_filterobj", "eta_matched_filterobj", "pt_leading_all_gen", "pt_subleading_all_gen", "eta_leading_all_gen", "eta_subleading_all_gen",
    "pt_leading_matched_gen", "pt_subleading_matched_gen", "eta_leading_matched_gen", "eta_subleading_matched_gen",
    "pt_leading_matched_filterobj", "pt_subleading_matched_filterobj", "eta_leading_matched_filterobj", "eta_subleading_matched_filterobj",
    "pt_fake_filterobj", "eta_fake_filterobj"  
]

histograms = {}
for path, _ in trigger_paths:
    histograms[path] = {}
    for category in histogram_categories:
        if "pt" in category:
            histograms[path][category] = ROOT.TH1D(
                f"{category}_{path}", f"{category.replace('_', ' ').title()} for {path}",
                histogram_config["pt"]["bins"], histogram_config["pt"]["min"], histogram_config["pt"]["max"]
            )
        elif "eta" in category:
            histograms[path][category] = ROOT.TH1D(
                f"{category}_{path}", f"{category.replace('_', ' ').title()} for {path}",
                histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"]
            )

histograms_with_cuts = {}
for path, _ in trigger_paths:
    histograms_with_cuts[path] = {}
    for category in histogram_categories:
        if "pt" in category:
            histograms_with_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withCuts", f"{category.replace('_', ' ').title()} for {path} with both taus passing pt>35, |eta|<2.1",
                histogram_config["pt"]["bins"], histogram_config["pt"]["min"], histogram_config["pt"]["max"]
            )
        elif "eta" in category:
            histograms_with_cuts[path][category] = ROOT.TH1D(
                f"{category}_{path}_withCuts", f"{category.replace('_', ' ').title()} for {path} with both taus passing pt>35, |eta|<2.1",
                histogram_config["eta"]["bins"], histogram_config["eta"]["min"], histogram_config["eta"]["max"]
            )
# Initialize cutflow counters
cutflow = {
    'TotalEvents': 0,
    'TwoGenTaus': 0,
    'BothTausPassCuts': 0,
    'L1P2GT_DoubleNNTau52': 0,
    'HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1': 0,
    'HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1': 0
}

products = []
add_product(products, "trig_sum", "trigger::TriggerEvent", "hltTriggerSummaryAOD")
add_product(products, "trig_res", "edm::TriggerResults", "TriggerResults")
add_product(products, "genparts", "std::vector<reco::GenParticle>", "genParticles")

evtdata = EvtData(products, verbose=False)
trigger_names = None
hadron_counter = 0

for event_nr, event in enumerate(events):
    cutflow['TotalEvents'] += 1
    if event_nr % 100 == 0:
        print(f"Processing event {event_nr}")

    evtdata.get_handles(event)

    trigger_names = event.object().triggerNames(evtdata.get("trig_res"))

    trigger_event = evtdata.get("trig_sum")
    gen_particles = evtdata.get("genparts")

    # Select hadronic taus from genParticles
    hadronic_taus = []
    for idx, particle in enumerate(gen_particles):
        if abs(particle.pdgId()) == 15 and particle.statusFlags().isPrompt() and particle.statusFlags().isLastCopy():
            if is_hadronic_tau(particle):
                hadronic_taus.append(particle)
                hadron_counter += 1

    hadronic_taus.sort(key=lambda x: x.pt(), reverse=True)

    leading_tau = hadronic_taus[0] if len(hadronic_taus) > 0 else None
    subleading_tau = hadronic_taus[1] if len(hadronic_taus) > 1 else None

    if len(hadronic_taus) == 2:
        cutflow['TwoGenTaus'] += 1
    else:
        continue  # Skip events with fewer than 2 hadronic taus

    # Check if both taus pass the pt and eta cuts
    passes_pt_lead = leading_tau.pt() > 35
    passes_eta_lead = abs(leading_tau.eta()) < 2.1
    passes_pt_sub = subleading_tau.pt() > 35
    passes_eta_sub = abs(subleading_tau.eta()) < 2.1

    passes_pt_eta_lead = passes_pt_lead and passes_eta_lead
    passes_pt_eta_sub = passes_pt_sub and passes_eta_sub
    both_taus_pass_cuts = passes_pt_eta_lead and passes_pt_eta_sub

    if both_taus_pass_cuts:
        cutflow['BothTausPassCuts'] += 1

    for path, filter_name in trigger_paths:
        # Fill histograms before applying the filter for all hadronic taus
        event_counter = 0
        for tau in hadronic_taus:
            histograms[path]["pt_all_gen"].Fill(tau.pt())
            histograms[path]["eta_all_gen"].Fill(tau.eta())

        if both_taus_pass_cuts:
            for tau in hadronic_taus:
                histograms_with_cuts[path]["pt_all_gen"].Fill(tau.pt())
                histograms_with_cuts[path]["eta_all_gen"].Fill(tau.eta())

        if leading_tau:
            # Fill leading tau histograms
            histograms[path]["pt_leading_all_gen"].Fill(leading_tau.pt())
            histograms[path]["eta_leading_all_gen"].Fill(leading_tau.eta())

            if both_taus_pass_cuts:
                histograms_with_cuts[path]["pt_leading_all_gen"].Fill(leading_tau.pt())
                histograms_with_cuts[path]["eta_leading_all_gen"].Fill(leading_tau.eta())

        # retrieving all trigger names and converting them to str() is necessary
        # because triggerNames are returned as basic_string_view<char,char_traits<char>>
        trigger_index = None
        for i in range(trigger_names.size()):
            if path in str(trigger_names.triggerName(i)):
                trigger_index = i
                break

        # If no match was found, `trigger_index` remains None
        if trigger_index is None:
            print(f"Warning: No loose match for trigger path {path} in event {event_nr}")
            continue

        trigger_accept = evtdata.get("trig_res").accept(trigger_index)

        if trigger_accept:
            cutflow[path] += 1
            for tau in hadronic_taus:
                histograms[path]["pt_pass_gen"].Fill(tau.pt())
                histograms[path]["eta_pass_gen"].Fill(tau.eta())

            if both_taus_pass_cuts:
                for tau in hadronic_taus:
                    histograms_with_cuts[path]["pt_pass_gen"].Fill(tau.pt())
                    histograms_with_cuts[path]["eta_pass_gen"].Fill(tau.eta())

        # Construct the InputTag with filter name and process name
        filter_tag = ROOT.edm.InputTag(filter_name, "", "HLTX")

        # Get index of the filter using filterIndex method with InputTag
        filter_index = trigger_event.filterIndex(filter_tag)
        if filter_index == trigger_event.sizeFilters():
            # print(f"Filter {filter_name} not found in this event.")
            continue  # Skip to the next path if filter is not there

        # Get trigger objects passing the filter
        filter_keys = trigger_event.filterKeys(filter_index)
        all_trigger_objects = trigger_event.getObjects()
        selected_objects = []
        for key in filter_keys:
            obj = all_trigger_objects[key]
            selected_objects.append(obj)
        for obj in selected_objects:
            # Fill histograms for all accepted objects per filter
            histograms[path]["pt_all_filterobj"].Fill(obj.pt())
            histograms[path]["eta_all_filterobj"].Fill(obj.eta())

            # Check if trigger object is matched to any hadronic tau
            matched_to_tau = False
            for tau in hadronic_taus:
                dR2 = ROOT.reco.deltaR2(tau.eta(), tau.phi(), obj.eta(), obj.phi())
                if dR2 < 0.01:
                    matched_to_tau = True
                    break
            if not matched_to_tau:
                # This is a fake trigger object
                histograms[path]["pt_fake_filterobj"].Fill(obj.pt())
                histograms[path]["eta_fake_filterobj"].Fill(obj.eta())

            if both_taus_pass_cuts:
                histograms_with_cuts[path]["pt_all_filterobj"].Fill(obj.pt())
                histograms_with_cuts[path]["eta_all_filterobj"].Fill(obj.eta())

                if not matched_to_tau:
                    # Fake object in events where both taus pass cuts
                    histograms_with_cuts[path]["pt_fake_filterobj"].Fill(obj.pt())
                    histograms_with_cuts[path]["eta_fake_filterobj"].Fill(obj.eta())

        if leading_tau:
            # Matching for Leading Tau
            matched_leading = False
            matched_obj_lead = None
            for obj in selected_objects:
                dR2 = ROOT.reco.deltaR2(leading_tau.eta(), leading_tau.phi(), obj.eta(), obj.phi())
                if dR2 < 0.01:
                    matched_leading = True
                    matched_obj_lead = obj
                    break

            if matched_leading:
                histograms[path]["pt_matched_gen"].Fill(leading_tau.pt())
                histograms[path]["eta_matched_gen"].Fill(leading_tau.eta())
                histograms[path]["pt_leading_matched_gen"].Fill(leading_tau.pt())
                histograms[path]["eta_leading_matched_gen"].Fill(leading_tau.eta())
                # Fill matched filter object histograms
                histograms[path]["pt_leading_matched_filterobj"].Fill(matched_obj_lead.pt())
                histograms[path]["eta_leading_matched_filterobj"].Fill(matched_obj_lead.eta())

                # Fill with cuts only if both taus pass cuts
                if both_taus_pass_cuts:
                    histograms_with_cuts[path]["pt_matched_gen"].Fill(leading_tau.pt())
                    histograms_with_cuts[path]["eta_matched_gen"].Fill(leading_tau.eta())
                    histograms_with_cuts[path]["pt_leading_matched_gen"].Fill(leading_tau.pt())
                    histograms_with_cuts[path]["eta_leading_matched_gen"].Fill(leading_tau.eta())
                    histograms_with_cuts[path]["pt_leading_matched_filterobj"].Fill(matched_obj_lead.pt())
                    histograms_with_cuts[path]["eta_leading_matched_filterobj"].Fill(matched_obj_lead.eta())

            # Process Subleading Tau Only If Leading Tau Is Matched and Subleading Tau exists
            if matched_leading and subleading_tau:
                # Fill subleading tau histograms
                histograms[path]["pt_subleading_all_gen"].Fill(subleading_tau.pt())
                histograms[path]["eta_subleading_all_gen"].Fill(subleading_tau.eta())

                # Apply cuts
                if passes_pt_sub:
                    histograms_with_cuts[path]["eta_subleading_all_gen"].Fill(subleading_tau.eta())
                if passes_eta_sub:
                    histograms_with_cuts[path]["pt_subleading_all_gen"].Fill(subleading_tau.pt())

                # Matching for Subleading Tau
                matched_subleading = False
                matched_obj_sub = None
                for obj in selected_objects:
                    dR2 = ROOT.reco.deltaR2(subleading_tau.eta(), subleading_tau.phi(), obj.eta(), obj.phi())
                    if dR2 < 0.01:
                        matched_subleading = True
                        matched_obj_sub = obj
                        break

                if matched_subleading:
                    histograms[path]["pt_matched_gen"].Fill(subleading_tau.pt())
                    histograms[path]["eta_matched_gen"].Fill(subleading_tau.eta())
                    histograms[path]["pt_subleading_matched_gen"].Fill(subleading_tau.pt())
                    histograms[path]["eta_subleading_matched_gen"].Fill(subleading_tau.eta())
                    # Fill matched filter object histograms
                    histograms[path]["pt_subleading_matched_filterobj"].Fill(matched_obj_sub.pt())
                    histograms[path]["eta_subleading_matched_filterobj"].Fill(matched_obj_sub.eta())

                    # Fill with cuts only if both taus pass cuts
                    if both_taus_pass_cuts:
                        histograms_with_cuts[path]["pt_matched_gen"].Fill(subleading_tau.pt())
                        histograms_with_cuts[path]["eta_matched_gen"].Fill(subleading_tau.eta())
                        histograms_with_cuts[path]["pt_subleading_matched_gen"].Fill(subleading_tau.pt())
                        histograms_with_cuts[path]["eta_subleading_matched_gen"].Fill(subleading_tau.eta())
                        histograms_with_cuts[path]["pt_subleading_matched_filterobj"].Fill(matched_obj_sub.pt())
                        histograms_with_cuts[path]["eta_subleading_matched_filterobj"].Fill(matched_obj_sub.eta())

# Save histograms to a ROOT file
output_file = ROOT.TFile("histograms.root", "RECREATE")
for path_hists in histograms.values():
    for hist in path_hists.values():
        hist.Write()
for path_hists in histograms_with_cuts.values():
    for hist in path_hists.values():
        hist.Write()
output_file.Close()

# Print cutflow summary
print("Cutflow Summary:")
for key, value in cutflow.items():
    print(f"{key}: {value}")

print("Histograms have been saved to histograms.root.")
print(f"Total events with at least two hadronic taus: {hadron_counter}")
