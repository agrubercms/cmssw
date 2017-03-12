#!/usr/bin/env python

import argparse
import os
import sys
import shutil
import subprocess

def readInput():
	parser = argparse.ArgumentParser(description="Perform the skimming for all samples in a folder")
	parser.add_argument("path",
	                    help="Input folder.")
	parser.add_argument("-c", "--config",
	                    default="$CMSSW_BASE/src/Validation/RecoTau/Tools/GetRecoTauVFromDQM_MC_cff.py",
	                    #default="$CMSSW_BASE/src/Validation/RecoTau/Tools/GetRecoTauVFromDQM_MC_cff_2.py",
	                    help="CMSSW configuration. [default: %(default)s]")
	parser.add_argument("-o", "--output-dir",
	                    default="$CMSSW_BASE/src/Validation/RecoTau/data/relval_plots",
	                    help="Output directory. [default: %(default)s]")

	args = parser.parse_args()
	arg.output_dir = os.path.expandvars(arg.output_dir)
	
	return args


def main():
	args = readInput()

	input_list = os.listdir(args.path)

	for sample in input_list:
		sample_name_new = sample.split("__RelVal")[-1]

		postfix = ""
		if "ZTT" in sample_name_new or "Zprime" in sample_name_new:
			postfix = "ZTT"
		if "ZMM" in sample_name_new:
			postfix = "ZMM"
		if "ZEE" in sample_name_new:
			postfix = "ZEE"
		if "TTbar" in sample_name_new or "QCD" in sample_name_new:
			postfix = "QCD"

		exe = "python " + args.config + " " + args.path + "/" + sample + " " + args.path + "/"  + sample_name_new + " " + postfix
		print exe
		p = subprocess.Popen(exe,
		stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE, stdin = subprocess.PIPE)
		(out, err) = p.communicate()
		print out
		print err


if __name__ == "__main__":
	main()
