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
    TEfficiency,
    gPad
)

def load_histogram(file, hist_path):

    hist = file.Get(hist_path)
    if not hist:
        return None
    if hist.GetEntries() == 0:
        return None
    hist.SetDirectory(0)  # Detach from file to prevent issues when file is closed
    hist.SetStats(False)  # Disable the histogram's statistics box
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

def plot_efficiencies(file, trigger_paths, variables, color_map, plot_types, output_dir):

    for var in variables:
        for plot_type in plot_types:
            for cut in ["", "_withCuts"]:
                plot_name = plot_type["name"]
                numerator_suffix = plot_type["numerator_suffix"]
                denominator_suffix = plot_type["denominator_suffix"]
                canvas_title = f"{plot_name}_efficiency_{var}{cut}"
                canvas = setup_canvas(canvas_title)
                if any("subleading" in name for name in plot_name.split("_")):
                    legend = TLegend(0.35, 0.13, 0.9, 0.3)
                else:
                    legend = TLegend(0.13, 0.7, 0.3, 0.9)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.03)

                efficiencies = []
                for path in trigger_paths:
                    numerator_name = f"{var}_{numerator_suffix}_{path}{cut}"
                    denominator_name = f"{var}_{denominator_suffix}_{path}{cut}"

                    numerator = load_histogram(file, numerator_name)
                    denominator = load_histogram(file, denominator_name)

                    if not numerator or not denominator:
                        print(f"Histograms not found for path '{path}' with numerator '{numerator_name}' and denominator '{denominator_name}'. Skipping.")
                        continue

                    efficiency = create_efficiency(numerator, denominator, color_map.get(path, kBlack))
                    if not efficiency:
                        print(f"Failed to create TEfficiency for path '{path}' with numerator '{numerator_name}' and denominator '{denominator_name}'. Skipping.")
                        print(f"numerator bins: {numerator.GetNbinsX()}, denominator bins: {denominator.GetNbinsX()}")
                        continue

                    # Append to the efficiencies list
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
                    legend.AddEntry(efficiency, path, "lep")

                if efficiencies:
                    legend.Draw()
                    output_path_png = os.path.join(output_dir, f"{plot_name}_efficiency_{var}{cut}.png")
                    canvas.SaveAs(output_path_png)
                    canvas.Close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_efficiencies.py <histograms.root>")
        sys.exit(1)

    root_file_path = sys.argv[1]

    if not os.path.isfile(root_file_path):
        print(f"Error: File '{root_file_path}' does not exist.")
        sys.exit(1)

    # Define trigger paths as in the histogram generation script
    trigger_paths = [
        "L1P2GT_DoubleNNTau52",
        "HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1",
        "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1"
    ]

    # Define variables to plot
    variables = ["pt", "eta"]

    # Define color mapping for trigger paths
    color_map = {
        "L1P2GT_DoubleNNTau52": ROOT.kBlue,
        "HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1": ROOT.kRed,
        "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1": ROOT.kGreen+2
    }

    # Define plot types
    plot_types = [
        {
            "name": "matched_gen",
            "numerator_suffix": "matched_gen",
            "denominator_suffix": "all_gen"
        },
        {
            "name": "pass_gen",
            "numerator_suffix": "pass_gen",
            "denominator_suffix": "all_gen"
        },
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
            "name": "fake_rate_filterobj",
            "numerator_suffix": "fake_filterobj",
            "denominator_suffix": "all_filterobj"
        }
    ]

    # Open the ROOT file
    file = ROOT.TFile.Open(root_file_path, "READ")
    if not file or file.IsZombie():
        print(f"Error: Cannot open ROOT file '{root_file_path}'.")
        sys.exit(1)

    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)  # Disable statistics box globally

    # Prepare output directory
    output_dir = "efficiency_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Plot efficiencies
    plot_efficiencies(file, trigger_paths, variables, color_map, plot_types, output_dir)

    # Close the ROOT file
    file.Close()

    print("All efficiency plots have been saved in the 'efficiency_plots' directory.")

if __name__ == "__main__":
    main()