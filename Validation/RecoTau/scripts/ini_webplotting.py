#!/usr/bin/env python

import argparse
import os
import Validation.RecoTau.webplotting as webplotting


def readInput():
	parser = argparse.ArgumentParser(description="Plot all standard validation distributions")
	parser.add_argument("-o", "--output-dir", default="/afs/cern.ch/work/d/dmroy/CMSSW_9_0_0_pre1/src/Validation/RecoTau/compare_ZTT",
	                    help="Folder to create index.html for")
	
	args = parser.parse_args()
	args.output_dir = os.path.expandvars(args.output_dir)
	
	return args


def main():
	args = readInput()

	webplotting.webplotting(input_dir=args.output_dir, recursive=True)


if __name__ == "__main__":
	main()
