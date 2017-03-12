#!/usr/bin/env python

import argparse
import os
import sys
import shutil
import subprocess

import Validation.RecoTau.tools as tools


def readInput():
	parser = argparse.ArgumentParser(description="Plot all standard validation distributions")
	parser.add_argument("paths", nargs=2,
	                    help="Two paths to compare.")
	parser.add_argument("--old-dm", action="store_true", default=False,
	                    help="Use old DecayModes instead of new ones. [default: %(default)s]")
	parser.add_argument("-o", "--output-dir", default="$CMSSW_BASE/src/Validation/RecoTau/data/relval_plots",
	                    help="Output directory. [default: %(default)s]")
	
	args = parser.parse_args()
	args.output_dir = os.path.expandvars(args.output_dir)
	
	return args


def main():
	args = readInput()

	exe = "python "+os.path.expandvars("$CMSSW_BASE/src/Validation/RecoTau/Tools/MultipleCompare.py")

	#exe + '-T ' + testsample  + "-R " + referencesample + "-t" +testtag+ " -r " + referencetag + distribution + " --logScaleY --maxLogY=100000 --maxYR=2.0 --rebin=10 --maxXaxis=80 -o " + outputfile
	if args.old_dm:
		discriminator_list = [
				["DecayModeFindingOldDMsEff", "/MVA6*ElectronRejectionEff", "_AntiEle"],
				["DecayModeFindingOldDMsEff","*CombinedIsolationDBSumPtCorr3HitsEff", "_Comb_Iso"],
				["DecayModeFindingOldDMsEff", "*IsolationMVArun2v1DBoldDMwLTEff", "_MVA_Iso"],
				["DecayModeFindingOldDMsEff", "/*MuonRejection3Eff", "_AntiMu"],
		]
	else:
		discriminator_list = [
				["DecayModeFindingNewDMsEff", "/MVA6*ElectronRejectionEff", "_AntiEle"],
				["DecayModeFindingNewDMsEff","*CombinedIsolationDBSumPtCorr3HitsEff", "_Comb_Iso"],
				["DecayModeFindingNewDMsEff", "*IsolationMVArun2v1DBnewDMwLTEff", "_MVA_Iso"],
				["DecayModeFindingNewDMsEff", "/*MuonRejection3Eff", "_AntiMu"],
		]

	type_list = ["eff", "fake"]

	kinematic_list = ["pt", "eta", "phi", "pileup"]

	distribution_list = []
	for kinematic in kinematic_list:
		for discriminator_tuple in discriminator_list:
			temp_list = []
			for discriminator in discriminator_tuple[:-1]:
				distribution = discriminator + kinematic
				temp_list.append(distribution)
			temp_list.append(discriminator_tuple[-1])
			distribution_list.append(temp_list)

	input_list_t = os.listdir(os.path.dirname(args.paths[0]))
	input_list_r = os.listdir(os.path.dirname(args.paths[1]))
	#print input_list_t
	#print input_list_r

	release_t = args.paths[0][args.paths[0].rfind("_CMSSW_")+1:].replace(".root", "")
	release_r = args.paths[1][args.paths[1].rfind("_CMSSW_")+1:].replace(".root", "")
	#print release_t
	#print release_r
	
	output_dir = os.path.join(args.output_dir, release_t+"_vs_"+release_r)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	input_output_tuple_list = []
	for sample in input_list_t:
		output_name = sample.split(".ro")[0]
		identifier = sample.split("_13_")[0]
		sampleMatch = ""
		for sample2 in input_list_r:
			if identifier in sample2:
				sampleMatch = sample2
		input_output_tuple_list.append([output_name, sample, sampleMatch])

	for type in type_list:
		for tuple in input_output_tuple_list:
			for distribution_tuple in distribution_list:

				type_postfix = ""
				rebin = ""
				maxXaxis = ""
				maxLogY = " --maxLogY=1000 "
				if "eff" in type:
					maxLogY = " --maxLogY=100 "
					if not "ZTT" in str(tuple[0]):
						continue
				else:
					type_postfix = "--fakeRate"
					rebin = " --rebin=10 "
					if "_AntiEle" in distribution_tuple[-1]and not "ZEE" in tuple[0]:
						continue
					if "_AntiMu" in distribution_tuple[-1] and not "ZMM" in tuple[0]:
						continue
					if "_Iso" in distribution_tuple[-1] and not ("TTbar" in tuple[0] or "QCD" in tuple[0]):
						continue

				distribution_str = " "
				for distribution in distribution_tuple[:-1]:
					distribution_str +=  " '" + distribution + "' "
				kin_str = ""
				if "eta" in distribution_str:
					kin_str = "_eta"
				elif "phi" in distribution_str:
					kin_str = "_phi"
					maxXaxis = " --maxXaxis=200 --minXaxis=-200 "
				elif "pileup" in distribution_str:
					kin_str = "_pileup"
					maxXaxis = " --maxXaxis=80  "
				else:
					kin_str = "_pt"
					maxXaxis = " --maxXaxis=200  "


				command = exe + " -T '" + os.path.join(os.path.dirname(args.paths[0]), tuple[1]) + "' -R '" + os.path.join(os.path.dirname(args.paths[1]), tuple[2]) + \
						"' -t '" + release_t + "' -r '" + release_r + "'" + str(distribution_str) + rebin + maxXaxis + maxLogY + \
						" --logScaleY --maxYR=2.0 " + \
						" -o " + os.path.join(output_dir,  type + distribution_tuple[-1] + kin_str + ".png") + " " + type_postfix
				print command
				(out, err) = tools.call_command(command)
				print command

# --logScaleY --maxLogY=100000 --maxYR=2.0 --rebin=10 --maxXaxis=80 " +



if __name__ == "__main__":
	main()
