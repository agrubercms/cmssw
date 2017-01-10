#!/usr/bin/env python

import os
import sys
import optparse
import shutil
import subprocess

def readInput():
	parser = optparse.OptionParser(description='Perform the skimming for all samples in a folder',
	                               usage='usage: %prog [options] folder ')
	parser.add_option('-c', '--config', action='store', metavar='config', dest="config",
	                  default="$CMSSW_BASE/src/Validation/RecoTau/Tools/GetRecoTauVFromDQM_MC_cff_2.py",
	                  help='Additional filter on the filenames. [default: %default]') dest="output_dir",
	parser.add_option('-o', '--output_dir', action='store', metavar='output_dir',
	                  default="/disk1/knutzen/TauPOG/RelVal/samples/",
	                  help='Output directory. [default: %default]')

	(options, args) = parser.parse_args()

	return options, args[0]



def main():
	options, path = readInput()

	input_list = os.listdir(path)

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

		exe = "python " + options.config + " " + path + "/" + sample + " " + path + "/"  + sample_name_new + " " + postfix
		print exe
		p = subprocess.Popen(exe,
		stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE, stdin = subprocess.PIPE)
		(out, err) = p.communicate()
		print out
		print err


if __name__ == '__main__':
	main()
