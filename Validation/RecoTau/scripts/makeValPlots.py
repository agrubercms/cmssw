#!/usr/bin/env python

import argparse
import os
import sys
import shutil
import subprocess

import Validation.RecoTau.tools as tools
import Validation.RecoTau.webplotting as webplotting


def readInput():
	parser = argparse.ArgumentParser(description="Plot all standard validation distributions")
	parser.add_argument("paths", nargs=2,
	                    help="Two paths to compare.")
	parser.add_argument("-o", "--output-dir", default="$CMSSW_BASE/src/Validation/RecoTau/data/relval_plots",
	                    help="Output directory. [default: %(default)s]")
	parser.add_argument("-t", "--tag", default="",
	                    help="Possible tag for output dir (will create one additional sub directory). [default: %(default)s]")
	parser.add_argument("-n", "--n-processes", type=int, default=1,
	                    help="Number of parallel processes. [default: %(default)s]")
	
	args = parser.parse_args()
	args.output_dir = os.path.expandvars(args.output_dir)
	
	return args


def main():
	args = readInput()

	exe = "python "+os.path.expandvars("$CMSSW_BASE/src/Validation/RecoTau/Tools/MultipleCompare.py")
	
	parameterset = {
		'ZTT':[['AntiEle_pt_eff',"DecayModeFindingNewDMsEffpt /MVA6*ElectronRejectionEffpt --minXaxis=0.001 --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_pt_eff',"DecayModeFindingNewDMsEffpt *CombinedIsolationDBSumPtCorr3HitsEffpt --minXaxis=0.001 --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_pt_eff',"DecayModeFindingNewDMsEffpt *IsolationMVArun2v1DBnewDMwLTEffpt --minXaxis=0.001 --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_pt_eff',"DecayModeFindingNewDMsEffpt /*MuonRejection3Effpt --minXaxis=0.001 --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_eta_eff',"DecayModeFindingNewDMsEffeta /MVA6*ElectronRejectionEffeta --minXaxis=-3 --maxXaxis=3 --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_eta_eff',"DecayModeFindingNewDMsEffeta *CombinedIsolationDBSumPtCorr3HitsEffeta --minXaxis=-3 --maxXaxis=3 --maxLogY=1.4  --logScaleY  "],
			   ['MVA_Iso_eta_eff',"DecayModeFindingNewDMsEffeta *IsolationMVArun2v1DBnewDMwLTEffeta --minXaxis=-3 --maxXaxis=3 --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_eta_eff',"DecayModeFindingNewDMsEffeta /*MuonRejection3Effeta --minXaxis=-3 --maxXaxis=3 --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_phi_eff',"DecayModeFindingNewDMsEffphi /MVA6*ElectronRejectionEffphi  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_phi_eff',"DecayModeFindingNewDMsEffphi *CombinedIsolationDBSumPtCorr3HitsEffphi  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_phi_eff',"DecayModeFindingNewDMsEffphi *IsolationMVArun2v1DBnewDMwLTEffphi  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_phi_eff',"DecayModeFindingNewDMsEffphi /*MuonRejection3Effphi  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_pileup_eff',"DecayModeFindingNewDMsEffpileup /MVA6*ElectronRejectionEffpileup --minXaxis=0.001 --maxXaxis=80   --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_pileup_eff',"DecayModeFindingNewDMsEffpileup *CombinedIsolationDBSumPtCorr3HitsEffpileup --minXaxis=0.001 --maxXaxis=80   --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_pileup_eff',"DecayModeFindingNewDMsEffpileup *IsolationMVArun2v1DBnewDMwLTEffpileup --minXaxis=0.001 --maxXaxis=80   --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_pileup_eff',"DecayModeFindingNewDMsEffpileup /*MuonRejection3Effpileup --minXaxis=0.001 --maxXaxis=80   --maxLogY=1.4  --logScaleY "]],

		'TTbar':[['Comb_Iso_pt_fake',"DecayModeFindingNewDMsEffpt *CombinedIsolationDBSumPtCorr3HitsEffpt  --rebin=10 --minXaxis=0.001 --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['MVA_Iso_pt_fake',"DecayModeFindingNewDMsEffpt *IsolationMVArun2v1DBnewDMwLTEffpt  --rebin=10 --minXaxis=0.001 --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['Comb_Iso_eta_fake',"DecayModeFindingNewDMsEffeta *CombinedIsolationDBSumPtCorr3HitsEffeta  --rebin=10 --minXaxis=-3 --maxXaxis=3 --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['MVA_Iso_eta_fake',"DecayModeFindingNewDMsEffeta *IsolationMVArun2v1DBnewDMwLTEffeta  --rebin=10 --minXaxis=-3 --maxXaxis=3 --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['Comb_Iso_phi_fake',"DecayModeFindingNewDMsEffphi *CombinedIsolationDBSumPtCorr3HitsEffphi  --rebin=10  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['MVA_Iso_phi_fake',"DecayModeFindingNewDMsEffphi *IsolationMVArun2v1DBnewDMwLTEffphi  --rebin=10  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['Comb_Iso_pileup_fake',"DecayModeFindingNewDMsEffpileup *CombinedIsolationDBSumPtCorr3HitsEffpileup  --rebin=10 --minXaxis=0.001 --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['MVA_Iso_pileup_fake',"DecayModeFindingNewDMsEffpileup *IsolationMVArun2v1DBnewDMwLTEffpileup  --rebin=10 --minXaxis=0.001 --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"]],

		'ZEE':[['AntiEle_pt_fake',"DecayModeFindingNewDMsEffpt /MVA6*ElectronRejectionEffpt  --rebin=10 --minXaxis=0.001 --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_eta_fake',"DecayModeFindingNewDMsEffeta /MVA6*ElectronRejectionEffeta  --rebin=10 --minXaxis=-3 --maxXaxis=3 --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_phi_fake',"DecayModeFindingNewDMsEffphi /MVA6*ElectronRejectionEffphi  --rebin=10  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_pileup_fake',"DecayModeFindingNewDMsEffpileup /MVA6*ElectronRejectionEffpileup  --rebin=10 --minXaxis=0.001 --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"]],

		'ZMM':[['AntiMu_pt_fake',"DecayModeFindingNewDMsEffpt /*MuonRejection3Effpt  --rebin=10 --minXaxis=0.001 --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiMu_eta_fake',"DecayModeFindingNewDMsEffeta /*MuonRejection3Effeta  --rebin=10 --minXaxis=-3 --maxXaxis=3 --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiMu_phi_fake',"DecayModeFindingNewDMsEffphi /*MuonRejection3Effphi  --rebin=10  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiMu_pileup_fake',"DecayModeFindingNewDMsEffpileup /*MuonRejection3Effpileup  --rebin=10 --minXaxis=0.001 --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"]],

#====================

		'SingleMuon':[['AntiEle_pt_eff',"DecayModeFindingNewDMsEffpt /MVA6*ElectronRejectionEffpt  --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['Comb_Iso_pt_eff',"DecayModeFindingNewDMsEffpt *CombinedIsolationDBSumPtCorr3HitsEffpt  --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_pt_eff',"DecayModeFindingNewDMsEffpt *IsolationMVArun2v1DBnewDMwLTEffpt  --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_pt_eff',"DecayModeFindingNewDMsEffpt /*MuonRejection3Effpt  --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_eta_eff',"DecayModeFindingNewDMsEffeta /MVA6*ElectronRejectionEffeta  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['Comb_Iso_eta_eff',"DecayModeFindingNewDMsEffeta *CombinedIsolationDBSumPtCorr3HitsEffeta  --maxLogY=1.4  --logScaleY  "],
			   ['MVA_Iso_eta_eff',"DecayModeFindingNewDMsEffeta *IsolationMVArun2v1DBnewDMwLTEffeta  --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_eta_eff',"DecayModeFindingNewDMsEffeta /*MuonRejection3Effeta  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_phi_eff',"DecayModeFindingNewDMsEffphi /MVA6*ElectronRejectionEffphi  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['Comb_Iso_phi_eff',"DecayModeFindingNewDMsEffphi *CombinedIsolationDBSumPtCorr3HitsEffphi  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_phi_eff',"DecayModeFindingNewDMsEffphi *IsolationMVArun2v1DBnewDMwLTEffphi  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_phi_eff',"DecayModeFindingNewDMsEffphi /*MuonRejection3Effphi  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_pileup_eff',"DecayModeFindingNewDMsEffpileup /MVA6*ElectronRejectionEffpileup  --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['Comb_Iso_pileup_eff',"DecayModeFindingNewDMsEffpileup *CombinedIsolationDBSumPtCorr3HitsEffpileup  --maxXaxis=80   --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_pileup_eff',"DecayModeFindingNewDMsEffpileup *IsolationMVArun2v1DBnewDMwLTEffpileup  --maxXaxis=80   --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_pileup_eff',"DecayModeFindingNewDMsEffpileup /*MuonRejection3Effpileup  --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"]],
 	}

	target = args.paths[0][args.paths[0].rfind("CMSSW"):].rstrip(".root")
	reference = args.paths[1][args.paths[1].rfind("CMSSW"):].rstrip(".root")

	target_run = args.paths[0][args.paths[0].rfind("_R0")+2:args.paths[0].find("__")]
	reference_run = args.paths[1][args.paths[1].rfind("_R0")+2:args.paths[1].find("__")]

	if target_run == "000000001" and reference_run == "000000001": #Doing FullSim or FastSim
	    del parameterset['SingleMuon']
	else: #Doing Data
	    del parameterset['ZTT']
	    del parameterset['TTbar']
	    del parameterset['ZEE']
	    del parameterset['ZMM']
	print parameterset

	target_version = "0001"
	for file in os.listdir(os.path.dirname(args.paths[0])):
	    version = args.paths[0][args.paths[0].rfind("DQM_V")+5:args.paths[0].find("_R0")]
	    if int(version) > int(target_version): target_version = version
	target_ver = int(target_version)
	reference_version = "0001"
	for file in os.listdir(os.path.dirname(args.paths[0])):
	    version = args.paths[1][args.paths[1].rfind("DQM_V")+5:args.paths[1].find("_R0")]
	    if int(version) > int(reference_version): reference_version = version
	reference_ver = int(reference_version)

	out_folder = target + '_vs_' + reference
	out_folder = out_folder.replace("CMSSW_","").replace("__DQMIO","").replace('_mcRun2_asymptotic_', "").replace('v1', '').replace('v0', '')
	out_folder = os.path.join(args.output_dir, args.tag, out_folder)
	if not os.path.exists(out_folder):
		os.makedirs(out_folder)
	
	commands = []
	for tar_ver in xrange(1,target_ver+1,1):
		for ref_ver in xrange(1,reference_ver+1,1):
			for evtt, parameters in parameterset.iteritems():
				if evtt != 'SingleMuon': evtt = 'RelVal' + evtt + "_13"
				for params in parameters:
					command = "python $CMSSW_BASE/src/Validation/RecoTau/Tools/MultipleCompare.py -T {target} -R {reference} -t {tfolder} -r {rfolder} {param2} -o {output}".format(
							target=os.path.join(os.path.dirname(args.paths[0]), "DQM_V000{tar_ver}_R{target_run}__{evtt}__{target}.root".format(tar_ver=tar_ver, target_run=target_run, evtt=evtt, target=target)),
							reference=os.path.join(os.path.dirname(args.paths[1]), "DQM_V000{ref_ver}_R{reference_run}__{evtt}__{reference}.root".format(ref_ver=ref_ver, reference_run=reference_run, evtt=evtt, reference=reference)),
							tfolder=target,
							rfolder=reference,
							param2=params[1],
							output=os.path.join(out_folder, "{evtt}__{target}_{param1}.png".format(evtt=evtt, target=target, param1=params[0]))
					)
					commands.append(command)
	
	tools.parallelize(tools.call_command, commands, n_processes=args.n_processes)
	webplotting.webplotting(input_dir=os.path.join(args.output_dir, args.tag), recursive=True)


if __name__ == "__main__":
	main()
