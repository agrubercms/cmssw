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
    TEfficiency,
    gPad
)

# Define the display mapping globally
display_map = {
    "HLT_DoublePNetTauh": "HLT_DoubleParTtauh"  # Use this name in plots
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

def plot_efficiencies(file, trigger_paths, variables, color_map, plot_types, output_dir, cut_levels):
    for var in variables:
        for plot_type in plot_types:
            for cut_level in cut_levels:
                plot_name = plot_type["name"]
                numerator_suffix = plot_type["numerator_suffix"]
                denominator_suffix = plot_type["denominator_suffix"]
                canvas_title = f"{plot_name}_efficiency_{var}{cut_level}"
                canvas = setup_canvas(canvas_title)
                
                legend = TLegend(0.13, 0.7, 0.3, 0.9)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.03)

                efficiencies = []
                for path in trigger_paths:
                    # Update histogram paths to include cut level
                    numerator = load_histogram(file, f"{var}_{numerator_suffix}_{path}{cut_level}")
                    denominator = load_histogram(file, f"{var}_{denominator_suffix}_{path}{cut_level}")

                    if not numerator or not denominator:
                        continue

                    efficiency = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
                    if not efficiency:
                        continue

                    efficiencies.append(efficiency)
                    gPad.Update()
                    # Determine the draw option
                    draw_option = "AP" if len(efficiencies) == 1 else "P SAME"
                    efficiency.Draw(draw_option)

                    # Update the pad to ensure the graph is painted
                    gPad.Update()

                    # Retrieve the painted graph
                    graph = efficiency.GetPaintedGraph()

                    # Set the axis titles and range on the graph
                    if graph:
                        if var == "pt":
                            graph.GetXaxis().SetTitle("p_{T} [GeV]")
                        elif var == "eta":
                            graph.GetXaxis().SetTitle("#eta")
                        graph.GetYaxis().SetTitle("Efficiency")
                        graph.SetMinimum(0)
                        graph.SetMaximum(1)

                    # Add to legend
                    legend.AddEntry(efficiency, display_map.get(path, path), "lep")

                if efficiencies:
                    legend.Draw()
                    # Ensure the output directory exists
                    os.makedirs(output_dir, exist_ok=True)
                    output_filename = f"{plot_name}_efficiency_{var}{cut_level}.png"
                    output_path_png = os.path.join(output_dir, output_filename)
                    canvas.SaveAs(output_path_png)
                    canvas.Close()
                    print(f"Saved plot: {output_path_png}")

def plot_distributions(file, output_dir, cut_levels, trigger_paths, hist_types):
    for path in trigger_paths:
        # Skip L1 histograms until we unpack them correctly
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

                    # Set histogram style
                    histogram.SetLineColor(color)
                    histogram.SetLineWidth(2)
                    histogram.SetMarkerColor(color)
                    histogram.SetMarkerStyle(20)
                    gPad.Update()

                    # Determine the draw option
                    draw_option = "HIST" if not histograms_plotted else "HIST SAME"
                    histogram.Draw(draw_option)
                    histograms_plotted = True
                    gPad.Update()
                    # Set axis titles
                    histogram.GetXaxis().SetTitle("p_{T} [GeV]")
                    histogram.GetYaxis().SetTitle("Events")

                    # Add to legend
                    legend.AddEntry(histogram, f"{tau.capitalize()} Tau", "l")

                if histograms_plotted:
                    legend.Draw()
                    # Ensure the output directory exists
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
    # Define the names of the 2D histograms
    histogram_names_2d = [
        "pt_gen_vs_reco_matched_leading",
        "pt_gen_vs_reco_matched_subleading",
    ]

    for path in trigger_paths:
        # Skip L1 histograms until we unpack them correctly
        if path.startswith("L1"):
            continue
            
        for hist_name in histogram_names_2d:
            for cut_level in cut_levels:
                # Construct the full histogram name
                full_hist_name = f"{hist_name}_{path}{cut_level}"
                histogram = load_histogram(file, full_hist_name)
                if not histogram:
                    continue

                canvas_title = f"{hist_name}_{path}{cut_level}"
                canvas = setup_canvas(canvas_title)

                # Set histogram titles
                histogram.SetTitle(f"{hist_name.replace('_', ' ').title()} for {path}")
                histogram.GetXaxis().SetTitle("Gen p_{T} [GeV]")
                histogram.GetYaxis().SetTitle("Reco p_{T} [GeV]")

                # Draw the 2D histogram
                histogram.Draw("COLZ")
                line = ROOT.TLine(histogram.GetXaxis().GetXmin(), histogram.GetXaxis().GetXmin(),
                                histogram.GetXaxis().GetXmax(), histogram.GetXaxis().GetXmax())
                line.SetLineColor(ROOT.kRed)
                line.SetLineStyle(1)
                line.Draw("SAME")
                # Ensure the output directory exists
                os.makedirs(output_dir, exist_ok=True)
                output_filename = f"{hist_name}_{path}{cut_level}.png"
                output_path_png = os.path.join(output_dir, output_filename)
                canvas.SaveAs(output_path_png)
                canvas.Close()
                print(f"Saved 2D histogram plot: {output_path_png}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_efficiencies.py <histograms.root>")
        sys.exit(1)

    root_file_path = sys.argv[1]
    ROOT.gROOT.SetBatch(True)

    if not os.path.isfile(root_file_path):
        print(f"Error: File '{root_file_path}' does not exist.")
        sys.exit(1)

    # Define trigger paths as used for file access (no renaming here)
    trigger_paths = [
        "L1P2GT_DoubleNNTau52",
        "HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1",
        "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1",
    ]

    variables = ["pt", "eta"]

    # Define color mapping for trigger paths
    color_map = {
        "L1P2GT_DoubleNNTau52": ROOT.kBlue,
        "HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1": ROOT.kRed,
        "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1": ROOT.kGreen+2,
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
        }
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
        }
    ]

    cut_levels = ['', '_withSingleCuts', '_withAllCuts']

    # Open the ROOT file and create plots
    file = ROOT.TFile.Open(root_file_path, "READ")
    if not file or file.IsZombie():
        print(f"Error: Cannot open ROOT file '{root_file_path}'.")
        sys.exit(1)

    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)

    output_dir = "efficiency_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Plot efficiencies
    plot_efficiencies(file, trigger_paths, variables, color_map, plot_types, output_dir, cut_levels)

    # Plot distributions
    plot_distributions(file, output_dir, cut_levels, trigger_paths, hist_types)

    # Plot 2D histograms
    plot_2d_histograms(file, output_dir, trigger_paths, cut_levels)

    file.Close()

    print("All plots have been saved in the 'efficiency_plots' directory.")

if __name__ == "__main__":
    main()
