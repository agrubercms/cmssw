#!/usr/bin/env python

import argparse
import os
import sys
import shutil

import Validation.RecoTau.tools as tools


def readInput():
	parser = argparse.ArgumentParser(description="Get the names of the input files from web.")
	parser.add_argument("identifier",
	                    help="Identifier (e.g. CMSSW_8_1_0_pre5)")
	parser.add_argument("-f", "--filter", default="",
	                    help="Additional filter on the filenames. [default: %(default)s]")
	parser.add_argument("-t", "--tag", default="",
	                    help="Additional name-tag for the folder. [default: %(default)s]")
	parser.add_argument("-v", "--veto", default="",
	                    help="Additional veto on the filenames. [default: %(default)s]")
	parser.add_argument("-o", "--output-dir", default="$CMSSW_BASE/src/Validation/RecoTau/data/relval_inputs",
	                    help="Output directory. [default: %(default)s]")
	parser.add_argument("-n", "--n-processes", type=int, default=1,
	                    help="Number of parallel processes. [default: %(default)s]")
	
	args = parser.parse_args()
	
	if not "CMSSW" in args.identifier:
		args.identifier = "CMSSW_"+args.identifier
	
	args.output_dir = os.path.expandvars(args.output_dir)
	
	return args

def copyFiles(exe, filelist, identifier, identifier_folder, args):
	outDir_parent = str(args.output_dir)
	outDir = os.path.join(outDir_parent, identifier + "_" + str(args.tag))
	if os.path.isdir(outDir):
		print ""
		input_str_outDir = raw_input(outDir + " already exists. Do you want to delete it? ([n] y): ")
		if input_str_outDir == "y":
			shutil.rmtree(outDir, ignore_errors = True)
	
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	commands = []
	for filename in filelist:
		outputFile = os.path.join(outDir, filename)
		if not os.path.exists(outputFile):
			commands.append("{exe}{filename} -o {outputFile}".format(exe=exe, filename=filename, outputFile=outputFile))
	
	tools.parallelize(tools.call_command, commands, n_processes=args.n_processes)

def main():
	args = readInput()
	
	converted_key = "userkey.pem"
	tools.convert_certificate_key_rsa("~/.globus/userkey.pem", converted_key)

	cmssw_version_raw_list = args.identifier.split("CMSSW_")[-1].split("_")

	if len(cmssw_version_raw_list) > 2:
		num0 = cmssw_version_raw_list[0]
		num1 = cmssw_version_raw_list[1]
		num2 = cmssw_version_raw_list[2]
	#elif len(cmssw_version_raw_list) == 4:
	#	num0 = cmssw_version_raw_list[0]
	#	num1 = cmssw_version_raw_list[1]
	#	num2 = cmssw_version_raw_list[2]
	else:
		print "Error: Cannot parse identifier."
		print "The identifier has to match the pattern: CMSSW_Z_Y_X(_A)"
		print "E.g.: CMSSW_8_1_0_pre5 or CMSSW_8_1_0 "
		sys.exit(1)

	identifier_folder = "CMSSW_" + num0 + "_" + num1 + "_x/"
	print identifier_folder

	exe = ("curl --anyauth --cert-type PEM --cert ~/.globus/usercert.pem --key {user_key} --key-type PEM -k " +
		   "https://cmsweb.cern.ch/dqm/relval/data/browse/ROOT/RelVal/{identifier_folder}").format(user_key=converted_key, identifier_folder=identifier_folder)

	print exe

	(curl_out, curl_err) = tools.call_command(exe)

	print curl_err

	filename_list = []

	#print curl_out

	output_list = curl_out.split("<tr><td><a href='")
	for line in output_list:
		if ".root'>" not in line or args.identifier+'-' not in line or str(args.filter) not in line or (str(args.veto) in line and str(args.veto)) \
			or not ("ZEE" in line or "ZMM" in line or "ZTT" in line or "QCD" in line or "TTbar" in line or "Zprime" in line or "TenTau" in line):
			continue
		filename = line.split(".root'>")[-1].split("</a></td><td>")[0]
		print filename
		filename_list.append(filename)

	input_files_decision = raw_input("Do you want to copy the above stated files to your local machine? (n [y]): ")
	if input_files_decision == "n":
	   sys.exit(0)

	copyFiles(exe, filename_list, args.identifier, identifier_folder, args)
	
	if os.path.exists(converted_key):
		os.remove(converted_key)

if __name__ == "__main__":
	main()
