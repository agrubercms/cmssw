#!/usr/bin/env python3
# plot_efficiencies.py

import sys
import os
import ROOT
from ROOT import (
    TFile,
    TCanvas,
    TLegend,
    gStyle,
    kBlue,
    kRed,
    kBlack,
    kGreen,
    kMagenta,
    kCyan,
    kYellow,
    TEfficiency,
    gPad
)

# Define the display mapping globally
display_map = {
    #"HLT_DoublePNetTauh": "hltParT"  # Use this name in plots
}

def load_histogram(file, hist_name):
    """Load a histogram by name"""
    hist = file.Get(hist_name)
    if not hist or hist.GetEntries() == 0:
        return None
    hist.SetDirectory(0)
    hist.SetStats(False)
    return hist

def create_efficiency(numerator, denominator, color):
    if not TEfficiency.CheckConsistency(numerator, denominator):
        return None

    efficiency = TEfficiency(numerator, denominator)
    efficiency.SetLineColor(color)
    efficiency.SetLineWidth(2)
    efficiency.SetMarkerColor(color)
    efficiency.SetMarkerStyle(20)
    return efficiency

def setup_canvas(title):
    canvas = TCanvas(title, title, 800, 600)
    canvas.SetGridx()
    canvas.SetGridy()
    return canvas

def draw_efficiency_frame(canvas, var, y_label="Efficiency", y_min=0, y_max=1.20, x_max=200):
    """Draw a frame with proper axis labels before overlaying TEfficiency objects."""
    if var == "pt":
        frame = canvas.DrawFrame(0, y_min, x_max, y_max)
        frame.GetXaxis().SetTitle("p_{T} [GeV]")
    elif var == "eta":
        frame = canvas.DrawFrame(-2.5, y_min, 2.5, y_max)
        frame.GetXaxis().SetTitle("#eta")
    else:
        frame = canvas.DrawFrame(0, y_min, x_max, y_max)
    frame.GetYaxis().SetTitle(y_label)
    return frame

def plot_efficiencies(file, trigger_paths, variables, color_map, plot_types, output_dir, cut_levels):
    for var in variables:
        for plot_type in plot_types:
            for cut_level in cut_levels:
                plot_name = plot_type["name"]
                numerator_suffix = plot_type["numerator_suffix"]
                denominator_suffix = plot_type["denominator_suffix"]
                canvas_title = f"{plot_name}_efficiency_{var}{cut_level}"
                canvas = setup_canvas(canvas_title)
                frame = draw_efficiency_frame(canvas, var)
                
                legend = TLegend(0.13, 0.7, 0.3, 0.9)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.03)

                efficiencies = []
                for path in trigger_paths:
                    numerator = load_histogram(file, f"{var}_{numerator_suffix}_{path}{cut_level}")
                    denominator = load_histogram(file, f"{var}_{denominator_suffix}_{path}{cut_level}")

                    if not numerator or not denominator:
                        continue

                    efficiency = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
                    if not efficiency:
                        continue

                    efficiencies.append(efficiency)
                    efficiency.Draw("P SAME")

                    legend.AddEntry(efficiency, display_map.get(path, path), "lep")

                # ---- New code: overlay PFJet efficiency as a fourth entry ----
                """rep = trigger_paths[2]
                if var == "pt":
                    pfjet_num_name = f"pt_pfjet_matched_{rep}{cut_level}"
                    pfjet_den_name = f"pt_all_gen_{rep}{cut_level}"
                elif var == "eta":
                    pfjet_num_name = f"eta_pfjet_matched_{rep}{cut_level}"
                    pfjet_den_name = f"eta_all_gen_{rep}{cut_level}"
                pfjet_numer = load_histogram(file, pfjet_num_name)
                pfjet_denom = load_histogram(file, pfjet_den_name)
                if pfjet_numer and pfjet_denom:
                    pfjet_eff = create_efficiency(pfjet_numer, pfjet_denom, kMagenta)
                    if pfjet_eff:
                        pfjet_eff.SetLineStyle(2)  # dashed style for differentiation
                        pfjet_eff.SetMarkerStyle(24)
                        pfjet_eff.Draw("P SAME")
                        gPad.Update()  # make sure the jet curve is drawn
                        efficiencies.append(pfjet_eff)
                        legend.AddEntry(pfjet_eff, "PFJet efficiency", "lep")
                # ---- End New code ----"""

                if efficiencies:
                    legend.Draw()
                    os.makedirs(output_dir, exist_ok=True)
                    output_filename = f"{plot_name}_efficiency_{var}{cut_level}.png"
                    output_path_png = os.path.join(output_dir, output_filename)
                    canvas.SaveAs(output_path_png)
                    canvas.Close()
                    print(f"Saved plot: {output_path_png}")

def plot_distributions(file, output_dir, cut_levels, trigger_paths, hist_types):
    for path in trigger_paths:
        if path.startswith("L1"):
            continue
        for cut_level in cut_levels:
            for hist_type in hist_types:
                suffix = hist_type["suffix"]
                label = hist_type["label"]
                canvas_title = f"{path}_{suffix}_pt_distribution{cut_level}"
                canvas = setup_canvas(canvas_title)

                legend = TLegend(0.7, 0.7, 0.9, 0.9)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.03)

                histograms_plotted = False

                for tau, color in zip(["subleading", "leading"], [kBlue, kRed]):
                    hist_name = f"pt_{tau}_{suffix}_{path}{cut_level}"
                    histogram = load_histogram(file, hist_name)
                    if not histogram:
                        print(f"Histogram '{hist_name}' not found. Skipping.")
                        continue

                    histogram.SetLineColor(color)
                    histogram.SetLineWidth(2)
                    histogram.SetMarkerColor(color)
                    histogram.SetMarkerStyle(20)
                    gPad.Update()

                    draw_option = "HIST" if not histograms_plotted else "HIST SAME"
                    histogram.Draw(draw_option)
                    histograms_plotted = True
                    gPad.Update()
                    histogram.GetXaxis().SetTitle("p_{T} [GeV]")
                    histogram.GetYaxis().SetTitle("Events")

                    legend.AddEntry(histogram, f"{tau.capitalize()} Tau", "l")

                if histograms_plotted:
                    legend.Draw()
                    os.makedirs(output_dir, exist_ok=True)
                    output_filename = f"z_{path}_{suffix}_pt_distribution{cut_level}.png"
                    output_path_png = os.path.join(output_dir, output_filename)
                    canvas.SaveAs(output_path_png)
                    canvas.Close()
                    print(f"Saved plot: {output_path_png}")
                else:
                    canvas.Close()
                    print(f"No histograms were plotted for {canvas_title}")

def plot_2d_histograms(file, output_dir, trigger_paths, cut_levels):
    # Updated list: added pt_gen_vs_reco_matched_barrel and pt_gen_vs_reco_matched_endcap
    histogram_names_2d = [
        "pt_gen_vs_reco_matched_leading",
        "pt_gen_vs_reco_matched_subleading",
        "pt_gen_vs_reco_matched_barrel",
        "pt_gen_vs_reco_matched_endcap"
    ]

    for path in trigger_paths:
        if path.startswith("L1"):
            continue
            
        for hist_name in histogram_names_2d:
            for cut_level in cut_levels:
                full_hist_name = f"{hist_name}_{path}{cut_level}"
                histogram = load_histogram(file, full_hist_name)
                if not histogram:
                    continue

                canvas_title = f"{hist_name}_{path}{cut_level}"
                canvas = setup_canvas(canvas_title)

                histogram.SetTitle(f"{hist_name.replace('_', ' ').title()} for {path}")
                histogram.GetXaxis().SetTitle("Gen p_{T} [GeV]")
                histogram.GetYaxis().SetTitle("Reco p_{T} [GeV]")

                histogram.Draw("COLZ")
                line = ROOT.TLine(histogram.GetXaxis().GetXmin(), histogram.GetXaxis().GetXmin(),
                                histogram.GetXaxis().GetXmax(), histogram.GetXaxis().GetXmax())
                line.SetLineColor(ROOT.kRed)
                line.SetLineStyle(1)
                line.Draw("SAME")
                os.makedirs(output_dir, exist_ok=True)
                output_filename = f"{hist_name}_{path}{cut_level}.png"
                output_path_png = os.path.join(output_dir, output_filename)
                canvas.SaveAs(output_path_png)
                canvas.Close()
                print(f"Saved 2D histogram plot: {output_path_png}")

def plot_fake_efficiencies(file, output_dir, cut_levels, trigger_paths, variables, color_map):
    for var in variables:
        for cut_level in cut_levels:
            canvas_title = f"fake_efficiency_{var}{cut_level}"
            canvas = setup_canvas(canvas_title)
            legend = TLegend(0.13, 0.7, 0.3, 0.9)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextSize(0.03)
            frame = draw_efficiency_frame(canvas, var, y_label="Fake Rate", y_max=1.0)
            efficiencies = []
            for path in trigger_paths:
                num_name = f"{var}_fake_filterobj_{path}{cut_level}"
                den_name = f"{var}_all_filterobj_{path}{cut_level}"
                numerator = load_histogram(file, num_name)
                denominator = load_histogram(file, den_name)
                if not numerator or not denominator:
                    continue
                color = color_map.get(path, kMagenta)
                efficiency = create_efficiency(numerator, denominator, color)
                if not efficiency:
                    continue
                efficiencies.append(efficiency)
                efficiency.Draw("P SAME")
                legend.AddEntry(efficiency, display_map.get(path, path), "lep")
            if efficiencies:
                legend.Draw()
                os.makedirs(output_dir, exist_ok=True)
                output_filename = f"fake_efficiency_{var}{cut_level}.png"
                output_path = os.path.join(output_dir, output_filename)
                canvas.SaveAs(output_path)
                canvas.Close()
                print(f"Saved fake efficiency plot: {output_path}")
            else:
                canvas.Close()

def plot_pfjet_distributions(file, output_dir, cut_levels, trigger_paths, pfjet_types):
    from ROOT import kBlack
    for path in trigger_paths:
        for cut_level in cut_levels:
            for pfjet in pfjet_types:
                hist_name = f"{pfjet['suffix']}_{path}{cut_level}"
                histogram = load_histogram(file, hist_name)
                if not histogram:
                    print(f"Histogram '{hist_name}' not found. Skipping.")
                    continue
                canvas_title = f"{path}_{pfjet['suffix']}_distribution{cut_level}"
                canvas = setup_canvas(canvas_title)
                histogram.SetLineColor(kBlack)
                histogram.SetLineWidth(2)
                histogram.Draw("HIST")
                histogram.GetXaxis().SetTitle(pfjet['xlabel'])
                histogram.GetYaxis().SetTitle(pfjet['ylabel'])
                os.makedirs(output_dir, exist_ok=True)
                output_filename = f"{pfjet['suffix']}_distribution_{path}{cut_level}.png"
                output_path = os.path.join(output_dir, output_filename)
                canvas.SaveAs(output_path)
                canvas.Close()
                print(f"Saved PFJet plot: {output_path}")

def plot_pfjet_efficiencies(file, output_dir, cut_levels, trigger_paths, color_map, display_map):
    from ROOT import kBlack
    for cut_level in cut_levels:
        if cut_level != "":
            continue
        canvas = setup_canvas(f"PFJet_pt_efficiency{cut_level}")
        frame = draw_efficiency_frame(canvas, "pt")
        legend = TLegend(0.15, 0.7, 0.35, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        efficiencies = []
        for path in trigger_paths:
            num_name = f"pt_pfjet_matched_{path}{cut_level}"
            den_name = f"pt_all_gen_{path}{cut_level}"
            numerator = load_histogram(file, num_name)
            denominator = load_histogram(file, den_name)
            if not numerator or not denominator:
                continue
            eff = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
            if not eff:
                continue
            eff.Draw("P SAME")
            efficiencies.append(eff)
            legend.AddEntry(eff, display_map.get(path, path), "lep")
        if efficiencies:
            #legend.Draw()
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"pfjet_pt_efficiency{cut_level}.png")
            canvas.SaveAs(output_path)
            canvas.Close()
            print(f"Saved PFJet efficiency plot: {output_path}")
        else:
            canvas.Close()

def plot_pfjet_eta_efficiencies(file, output_dir, cut_levels, trigger_paths, color_map, display_map):
    from ROOT import kBlack
    for cut_level in cut_levels:
        if cut_level != "":
            continue
        canvas = setup_canvas(f"PFJet_eta_efficiency{cut_level}")
        frame = draw_efficiency_frame(canvas, "eta")
        legend = TLegend(0.15, 0.7, 0.35, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        efficiencies = []
        for path in trigger_paths:
            num_name = f"eta_pfjet_matched_{path}{cut_level}"
            den_name = f"eta_all_gen_{path}{cut_level}"
            numerator = load_histogram(file, num_name)
            denominator = load_histogram(file, den_name)
            if not numerator or not denominator:
                continue
            eff = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
            if not eff:
                continue
            eff.Draw("P SAME")
            efficiencies.append(eff)
            legend.AddEntry(eff, display_map.get(path, path), "lep")
        if efficiencies:
            #legend.Draw()
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"pfjet_eta_efficiency{cut_level}.png")
            canvas.SaveAs(output_path)
            canvas.Close()
            print(f"Saved PFJet eta efficiency plot: {output_path}")
        else:
            canvas.Close()

def plot_pfjet_leading_pt_efficiency(file, output_dir, trigger_paths, color_map, display_map):
    from ROOT import kBlack
    canvas = setup_canvas("PFJet_leading_pt_efficiency")
    frame = draw_efficiency_frame(canvas, "pt")
    legend = TLegend(0.15, 0.7, 0.35, 0.9)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.03)
    efficiencies = []
    for path in trigger_paths:
        num_name = f"pt_leading_pfjet_matched_{path}"
        den_name = f"pt_leading_all_gen_{path}"
        numerator = load_histogram(file, num_name)
        denominator = load_histogram(file, den_name)
        if not numerator or not denominator:
            continue
        eff = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
        if not eff:
            continue
        eff.Draw("P SAME")
        efficiencies.append(eff)
        legend.AddEntry(eff, display_map.get(path, path), "lep")
    if efficiencies:
        #legend.Draw()
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "PFJet_leading_pt_efficiency.png")
        canvas.SaveAs(output_path)
        canvas.Close()
        print(f"Saved PFJet leading pt efficiency plot: {output_path}")
    else:
        canvas.Close()

def plot_pfjet_subleading_pt_efficiency(file, output_dir, trigger_paths, color_map, display_map):
    from ROOT import kBlack
    canvas = setup_canvas("PFJet_subleading_pt_efficiency")
    frame = draw_efficiency_frame(canvas, "pt")
    legend = TLegend(0.15, 0.7, 0.35, 0.9)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.03)
    efficiencies = []
    for path in trigger_paths:
        num_name = f"pt_subleading_pfjet_matched_{path}"
        den_name = f"pt_subleading_all_gen_{path}"
        numerator = load_histogram(file, num_name)
        denominator = load_histogram(file, den_name)
        if not numerator or not denominator:
            continue
        eff = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
        if not eff:
            continue
        eff.Draw("P SAME")
        efficiencies.append(eff)
        legend.AddEntry(eff, display_map.get(path, path), "lep")
    if efficiencies:
        #legend.Draw()
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "PFJet_subleading_pt_efficiency.png")
        canvas.SaveAs(output_path)
        canvas.Close()
        print(f"Saved PFJet subleading pt efficiency plot: {output_path}")
    else:
        canvas.Close()

def plot_pfjet_leading_eta_efficiency(file, output_dir, trigger_paths, color_map, display_map):
    from ROOT import kBlack
    canvas = setup_canvas("PFJet_leading_eta_efficiency")
    frame = draw_efficiency_frame(canvas, "eta")
    legend = TLegend(0.15, 0.7, 0.35, 0.9)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.03)
    efficiencies = []
    for path in trigger_paths:
        num_name = f"eta_leading_pfjet_matched_{path}"
        den_name = f"eta_leading_all_gen_{path}"
        numerator = load_histogram(file, num_name)
        denominator = load_histogram(file, den_name)
        if not numerator or not denominator:
            continue
        eff = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
        if not eff:
            continue
        eff.Draw("P SAME")
        efficiencies.append(eff)
        legend.AddEntry(eff, display_map.get(path, path), "lep")
    if efficiencies:
        #legend.Draw()
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "PFJet_leading_eta_efficiency.png")
        canvas.SaveAs(output_path)
        canvas.Close()
        print(f"Saved PFJet leading eta efficiency plot: {output_path}")
    else:
        canvas.Close()

def plot_pfjet_subleading_eta_efficiency(file, output_dir, trigger_paths, color_map, display_map):
    from ROOT import kBlack
    canvas = setup_canvas("PFJet_subleading_eta_efficiency")
    frame = draw_efficiency_frame(canvas, "eta")
    legend = TLegend(0.15, 0.7, 0.35, 0.9)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.03)
    efficiencies = []
    for path in trigger_paths:
        num_name = f"eta_subleading_pfjet_matched_{path}"
        den_name = f"eta_subleading_all_gen_{path}"
        numerator = load_histogram(file, num_name)
        denominator = load_histogram(file, den_name)
        if not numerator or not denominator:
            continue
        eff = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
        if not eff:
            continue
        eff.Draw("P SAME")
        efficiencies.append(eff)
        legend.AddEntry(eff, display_map.get(path, path), "lep")
    if efficiencies:
        #legend.Draw()
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "PFJet_subleading_eta_efficiency.png")
        canvas.SaveAs(output_path)
        canvas.Close()
        print(f"Saved PFJet subleading eta efficiency plot: {output_path}")
    else:
        canvas.Close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_efficiencies.py <histograms.root>")
        sys.exit(1)

    root_file_path = sys.argv[1]
    ROOT.gROOT.SetBatch(True)

    if not os.path.isfile(root_file_path):
        print(f"Error: File '{root_file_path}' does not exist.")
        sys.exit(1)

    trigger_paths = [
        "L1P2GT_DoubleNNTau52",
        "HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1",
        "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1",
        #"HLT_DoublePNetTauh",
        "HLT_DoubleMediumPFPuppiParTTauh30_eta2p1",
    ]

    variables = ["pt", "eta"]

    color_map = {
        "L1P2GT_DoubleNNTau52": kBlue,
        "HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1": kRed,
        "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1": kGreen+2,
        #"HLT_DoublePNetTauh": kBlack,
        "HLT_DoubleMediumPFPuppiParTTauh30_eta2p1": kYellow+2,
    }

    plot_types = [
        {
            "name": "matched_gen_leading",
            "numerator_suffix": "leading_matched_gen",
            "denominator_suffix": "leading_all_gen"
        },
        {
            "name": "matched_gen_subleading",
            "numerator_suffix": "subleading_matched_gen",
            "denominator_suffix": "subleading_all_gen"
        },
        {
            "name": "l1matched_gen_leading",
            "numerator_suffix": "leading_l1matched_gen",
            "denominator_suffix": "leading_all_gen"
        },
        {
            "name": "l1matched_gen_subleading",
            "numerator_suffix": "subleading_l1matched_gen",
            "denominator_suffix": "subleading_all_gen"
        },
        {
            "name": "pass_gen_leading",
            "numerator_suffix": "leading_pass_gen",
            "denominator_suffix": "leading_all_gen"
        },
        {
            "name": "pass_gen_subleading",
            "numerator_suffix": "subleading_pass_gen",
            "denominator_suffix": "subleading_all_gen"
        },

        # NEW: efficiencies where BOTH gen taus are matched
        {
            "name": "matched_gen_both_leading",
            "numerator_suffix": "leading_matched_gen_both",
            "denominator_suffix": "leading_all_gen"
        },
        {
            "name": "matched_gen_both_subleading",
            "numerator_suffix": "subleading_matched_gen_both",
            "denominator_suffix": "subleading_all_gen"
        },
    ]

    hist_types = [
        {
            "suffix": "all_gen",
            "label": "All Gen Taus"
        },
        {
            "suffix": "matched_gen",
            "label": "Matched Gen Taus"
        },
        {
            "suffix": "l1matched_gen",
            "label": "L1 Matched Gen Taus"
        },
        {
            "suffix": "all_filterobj",
            "label": "All Filter Objects"
        },
        {
            "suffix": "matched_filterobj",
            "label": "Matched Filter Objects"
        },
        # NEW: distributions for events where BOTH gen taus are matched
        {
            "suffix": "matched_gen_both",
            "label": "Matched Gen Taus (both)"
        },
    ]

    cut_levels = ['', '_withSingleCuts', '_withAllCuts']

    file = ROOT.TFile.Open(root_file_path, "READ")
    if not file or file.IsZombie():
        print(f"Error: Cannot open ROOT file '{root_file_path}'.")
        sys.exit(1)

    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)

    output_dir = "efficiency_plots"
    os.makedirs(output_dir, exist_ok=True)

    plot_efficiencies(file, trigger_paths, variables, color_map, plot_types, output_dir, cut_levels)
    plot_distributions(file, output_dir, cut_levels, trigger_paths, hist_types)
    plot_2d_histograms(file, output_dir, trigger_paths, cut_levels)
    plot_fake_efficiencies(file, output_dir, cut_levels, trigger_paths, variables, color_map)

    pfjet_types = [
         {"suffix": "n_pfjets", "xlabel": "Number of PFJets", "ylabel": "Events"},
         {"suffix": "pt_all_pfjet", "xlabel": "p_{T} [GeV]", "ylabel": "Events"},
         {"suffix": "eta_all_pfjet", "xlabel": "#eta", "ylabel": "Events"},
         {"suffix": "pt_leading_pfjet", "xlabel": "p_{T} [GeV]", "ylabel": "Events"},
         {"suffix": "eta_leading_pfjet", "xlabel": "#eta", "ylabel": "Events"},
         {"suffix": "pt_subleading_pfjet", "xlabel": "p_{T} [GeV]", "ylabel": "Events"},
         {"suffix": "eta_subleading_pfjet", "xlabel": "#eta", "ylabel": "Events"}
    ]
    plot_pfjet_distributions(file, output_dir, cut_levels, trigger_paths, pfjet_types)

    plot_pfjet_efficiencies(file, output_dir, cut_levels, trigger_paths, color_map, display_map)
    plot_pfjet_eta_efficiencies(file, output_dir, cut_levels, trigger_paths, color_map, display_map)

    plot_pfjet_leading_pt_efficiency(file, output_dir, trigger_paths, color_map, display_map)
    plot_pfjet_subleading_pt_efficiency(file, output_dir, trigger_paths, color_map, display_map)
    plot_pfjet_leading_eta_efficiency(file, output_dir, trigger_paths, color_map, display_map)
    plot_pfjet_subleading_eta_efficiency(file, output_dir, trigger_paths, color_map, display_map)

    file.Close()

    print("All plots have been saved in the 'efficiency_plots' directory.")

if __name__ == "__main__":
    main()
