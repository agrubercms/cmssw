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
    """
    Load a histogram from a ROOT file.

    Parameters:
    - file: TFile object
    - hist_path: Full path to the histogram within the ROOT file

    Returns:
    - TH1 object if valid, else None
    """
    hist = file.Get(hist_path)
    if not hist:
        return None
    if hist.GetEntries() == 0:
        return None
    hist.SetDirectory(0)  # Detach from file to prevent issues when file is closed
    hist.SetStats(False)  # Disable the histogram's statistics box
    return hist

def create_efficiency(numerator, denominator, color):
    """
    Create a TEfficiency object from numerator and denominator histograms.

    Parameters:
    - numerator: TH1 object representing passed events
    - denominator: TH1 object representing total events
    - color: ROOT color code for styling

    Returns:
    - TEfficiency object if successful, else None
    """
    if not TEfficiency.CheckConsistency(numerator, denominator):
        return None

    efficiency = TEfficiency(numerator, denominator)
    efficiency.SetLineColor(color)
    efficiency.SetLineWidth(2)
    efficiency.SetMarkerColor(color)
    efficiency.SetMarkerStyle(20)
    return efficiency

def setup_canvas(title):
    """
    Setup a ROOT canvas with grid and return it.

    Parameters:
    - title: String title for the canvas

    Returns:
    - TCanvas object
    """
    canvas = TCanvas(title, title, 800, 600)
    canvas.SetGridx()
    canvas.SetGridy()
    return canvas

def plot_efficiencies(file, trigger_paths, variables, color_map, plot_types, output_dir):
    """
    Plot various types of efficiencies based on the specified plot types.

    Parameters:
    - file: TFile object
    - trigger_paths: List of trigger path names
    - variables: List of variables to plot ('pt' and 'eta')
    - color_map: Dictionary mapping trigger paths to ROOT color codes
    - plot_types: List of dictionaries specifying plot types
    - output_dir: Directory to save the plots
    """
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
                            graph.GetXaxis().SetTitle("Gen p_{T}^{#tau} [GeV]")
                        elif var == "eta":
                            graph.GetXaxis().SetTitle("Gen #eta")
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
        # Add more plot types here if needed, e.g., leading and subleading
        {
            "name": "matched_gen_leading",
            "numerator_suffix": "leading_matched_gen",
            "denominator_suffix": "leading_all_gen"
        },
        {
            "name": "matched_gen_subleading",
            "numerator_suffix": "subleading_matched_gen",
            "denominator_suffix": "subleading_all_gen"
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


"""#!/usr/bin/env python3
# plot_efficiencies_extended.py

import sys
import os
import ROOT
import argparse
import logging
from ROOT import (
    TFile,
    TCanvas,
    TLegend,
    gStyle,
    kBlue,
    kRed,
    kGreen,
    kBlack,
    TEfficiency,
)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot trigger efficiencies from a ROOT file, including leading and subleading objects.")
    parser.add_argument("root_file", help="Path to the ROOT file containing histograms.")
    parser.add_argument("-o", "--output", default="efficiency_plots", help="Output directory for plots.")
    parser.add_argument("-v", "--variables", nargs="+", default=["pt", "eta", "pt_leading", "pt_subleading", "eta_leading", "eta_subleading"], help="Variables to plot.")
    parser.add_argument("--list", action="store_true", help="List all histograms in the ROOT file.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")
    return parser.parse_args()

def setup_logging(verbose=False):
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')

def load_and_check_histogram(file, hist_path):
    hist = file.Get(hist_path)
    if not hist:
        logging.error(f"Histogram '{hist_path}' not found.")
        return None
    if hist.GetEntries() == 0:
        logging.warning(f"Histogram '{hist_path}' is empty.")
    if not isinstance(hist, ROOT.TH1):
        logging.error(f"Object '{hist_path}' is not a histogram.")
        return None
    hist.SetDirectory(0)  # Detach from file to prevent issues when file is closed
    hist.SetStats(False)  # Disable the histogram's statistics box
    return hist

def create_efficiency_histogram(numerator, denominator, color, path, var):
    if not numerator or not denominator:
        logging.error(f"Invalid histograms provided for TEfficiency: path='{path}', var='{var}'.")
        return None

    # Check consistency
    if not TEfficiency.CheckConsistency(numerator, denominator):
        logging.error(f"Numerator and Denominator histograms are not consistent for path '{path}', var='{var}'.")
        return None

    # Create TEfficiency object
    efficiency = TEfficiency(numerator, denominator)
    if not isinstance(efficiency, ROOT.TEfficiency):
        logging.error(f"TEfficiency object is invalid for path '{path}', var='{var}'.")
        return None

    efficiency.SetLineColor(color)
    efficiency.SetLineWidth(2)
    efficiency.SetMarkerColor(color)
    efficiency.SetMarkerStyle(20)
    return efficiency

def setup_canvas(title):
    canvas = TCanvas(title, title, 800, 600)
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.15)
    return canvas

def list_all_histograms(file):
    keys = file.GetListOfKeys()
    for key in keys:
        obj = key.ReadObj()
        if isinstance(obj, ROOT.TH1):
            print(f" - {obj.GetName()}: {obj.GetTitle()}")

def get_axis_titles(var):
    if var.startswith("pt"):
        x_title = "Gen p_{T}^{#tau} [GeV]"
    elif var.startswith("eta"):
        x_title = "Gen #eta"
    else:
        x_title = var  # Fallback title

    if "leading" in var:
        x_title += " (Leading)"
    elif "subleading" in var:
        x_title += " (Subleading)"
    return x_title

def plot_efficiencies(file, variables, trigger_paths, color_map, output_dir, suffix=''):

    for var in variables:
        logging.debug(f"Processing variable: '{var}' with suffix '{suffix}'")
        canvas_title = f"efficiency_{var}{suffix}"
        canvas = setup_canvas(canvas_title)
        legend = TLegend(0.7, 0.7, 0.9, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)

        first_efficiency = True
        efficiencies = []  # List to hold TEfficiency objects
        for path in trigger_paths:
            logging.debug(f"Processing trigger path: '{path}' for variable '{var}' with suffix '{suffix}'")

            # Define histogram names based on your naming convention
            # Assuming 'matched' refers to matched_gen in numerator
            numerator_name = f"{var}_matched_gen_{path}{suffix}"
            denominator_name = f"{var}_all_gen_{path}{suffix}"

            # Load histograms
            numerator = load_and_check_histogram(file, numerator_name)
            denominator = load_and_check_histogram(file, denominator_name)

            if not numerator or not denominator:
                logging.warning(f"Skipping efficiency for path '{path}' due to missing histograms with suffix '{suffix}'.")
                continue

            # Create efficiency object
            efficiency = create_efficiency_histogram(
                numerator, denominator, color_map.get(path, kBlack), path, var
            )
            if not efficiency:
                logging.warning(f"Failed to create TEfficiency for path '{path}', variable '{var}' with suffix '{suffix}'. Skipping.")
                continue

            # Keep reference to the efficiency object
            efficiencies.append(efficiency)

            # Draw efficiency
            try:
                if first_efficiency:
                    efficiency.Draw("AP")  # A: Axis, P: Points
                    painted_graph = efficiency.GetPaintedGraph()
                    if painted_graph:
                        print(f"painted_graph: {painted_graph}")
                        x_title = get_axis_titles(var)
                        painted_graph.GetXaxis().SetTitle(x_title)
                        painted_graph.GetYaxis().SetTitle("Efficiency")
                        # Optionally, set axis ranges
                        painted_graph.GetYaxis().SetRangeUser(0, 1.1)
                    first_efficiency = False
                else:
                    efficiency.Draw("P SAME")  # Draw on the same canvas
                    painted_graph = efficiency.GetPaintedGraph()
                    if painted_graph:
                        painted_graph.GetYaxis().SetRangeUser(0, 1.1)
            except Exception as e:
                logging.error(f"Exception occurred while drawing TEfficiency for path '{path}', var '{var}' with suffix '{suffix}': {e}")
                continue

            # Add to legend
            legend.AddEntry(efficiency, path, "lep")

        if not first_efficiency:
            # Draw legend
            logging.debug(f"Drawing legend for variable '{var}' with suffix '{suffix}'.")
            legend.Draw()

            var_dir = output_dir
            os.makedirs(var_dir, exist_ok=True)

            output_filename = f"trigger_efficiency_{var}{suffix}.png"
            output_path_png = os.path.join(var_dir, output_filename)
            try:
                canvas.SaveAs(output_path_png)
                logging.info(f"Saved plot: {output_path_png}")
            except Exception as e:
                logging.error(f"Failed to save canvas for variable '{var}' with suffix '{suffix}': {e}")
        else:
            logging.info(f"No efficiencies were plotted for variable '{var}' with suffix '{suffix}'.")

        canvas.Close()
                # --- Plotting Filter Object Efficiencies ---
        canvas = setup_canvas(f"filterobj_efficiency_{var}")
        legend = TLegend(0.7, 0.7, 0.9, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)

        efficiencies = []


def main():
    args = parse_arguments()
    setup_logging(args.verbose)

    root_file_path = args.root_file

    if not os.path.isfile(root_file_path):
        logging.error(f"File '{root_file_path}' does not exist.")
        sys.exit(1)

    # Define trigger paths as in the histogram generation script
    trigger_paths = [
        "L1P2GT_DoubleNNTau52",
        "HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1",
        "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1"
    ]

    # Define variables to plot (including leading and subleading)
    variables = args.variables

    # Define color mapping for trigger paths
    color_map = {
        "L1P2GT_DoubleNNTau52": kBlue,
        "HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1": kRed,
        "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1": kGreen+2
    }

    # Open the ROOT file
    logging.debug(f"Opening ROOT file: '{root_file_path}'")
    file = TFile.Open(root_file_path, "READ")
    if not file or file.IsZombie():
        logging.error(f"Cannot open ROOT file '{root_file_path}'.")
        sys.exit(1)

    if args.list:
        logging.info("Listing all histograms in the ROOT file:")
        list_all_histograms(file)
        file.Close()
        sys.exit(0)

    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)  
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    # First round of plotting without cuts
    plot_efficiencies(
        file=file,
        variables=variables,
        trigger_paths=trigger_paths,
        color_map=color_map,
        output_dir=output_dir,
        suffix=''  # No suffix for the first round
    )

    # Second round of plotting with '_withCuts' suffix
    plot_efficiencies(
        file=file,
        variables=variables,
        trigger_paths=trigger_paths,
        color_map=color_map,
        output_dir=output_dir,
        suffix='_withCuts'  # Suffix for the second round
    )

    file.Close()

if __name__ == "__main__":
    main()"""