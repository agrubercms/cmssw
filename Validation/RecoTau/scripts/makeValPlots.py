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
                ########### ZTT old
		'ZTT':[['AntiEle_pt_eff',"DecayModeFindingNewDMsEffpt /MVA6*ElectronRejectionEffpt  --rebin=10 --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_pt_eff',"DecayModeFindingNewDMsEffpt *CombinedIsolationDBSumPtCorr3HitsEffpt  --rebin=10 --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_pt_eff',"DecayModeFindingNewDMsEffpt *IsolationMVArun2v1DBnewDMwLTEffpt  --rebin=10 --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_pt_eff',"DecayModeFindingNewDMsEffpt /*MuonRejection3Effpt  --rebin=10 --maxXaxis=200   --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_eta_eff',"DecayModeFindingNewDMsEffeta /MVA6*ElectronRejectionEffeta --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_eta_eff',"DecayModeFindingNewDMsEffeta *CombinedIsolationDBSumPtCorr3HitsEffeta --rebin=10 --maxLogY=1.4  --logScaleY  "],
			   ['MVA_Iso_eta_eff',"DecayModeFindingNewDMsEffeta *IsolationMVArun2v1DBnewDMwLTEffeta --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_eta_eff',"DecayModeFindingNewDMsEffeta /*MuonRejection3Effeta --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_phi_eff',"DecayModeFindingNewDMsEffphi /MVA6*ElectronRejectionEffphi  --maxXaxis=200 --minXaxis=-200 --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_phi_eff',"DecayModeFindingNewDMsEffphi *CombinedIsolationDBSumPtCorr3HitsEffphi  --maxXaxis=200 --minXaxis=-200 --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_phi_eff',"DecayModeFindingNewDMsEffphi *IsolationMVArun2v1DBnewDMwLTEffphi  --maxXaxis=200 --minXaxis=-200 --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_phi_eff',"DecayModeFindingNewDMsEffphi /*MuonRejection3Effphi  --maxXaxis=200 --minXaxis=-200 --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_pileup_eff',"DecayModeFindingNewDMsEffpileup /MVA6*ElectronRejectionEffpileup  --maxXaxis=80  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_pileup_eff',"DecayModeFindingNewDMsEffpileup *CombinedIsolationDBSumPtCorr3HitsEffpileup  --maxXaxis=80  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_pileup_eff',"DecayModeFindingNewDMsEffpileup *IsolationMVArun2v1DBnewDMwLTEffpileup  --maxXaxis=80  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_pileup_eff',"DecayModeFindingNewDMsEffpileup /*MuonRejection3Effpileup  --maxXaxis=80  --rebin=10 --maxLogY=1.4  --logScaleY "]],

                ########### TenTau old
		'TenTau':[['AntiEle_pt_eff',"DecayModeFindingNewDMsEffpt /MVA6*ElectronRejectionEffpt  --maxXaxis=200  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_pt_eff',"DecayModeFindingNewDMsEffpt *CombinedIsolationDBSumPtCorr3HitsEffpt  --maxXaxis=200  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_pt_eff',"DecayModeFindingNewDMsEffpt *IsolationMVArun2v1DBnewDMwLTEffpt  --maxXaxis=200  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_pt_eff',"DecayModeFindingNewDMsEffpt /*MuonRejection3Effpt  --maxXaxis=200  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_eta_eff',"DecayModeFindingNewDMsEffeta /MVA6*ElectronRejectionEffeta --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_eta_eff',"DecayModeFindingNewDMsEffeta *CombinedIsolationDBSumPtCorr3HitsEffeta --rebin=10 --maxLogY=1.4  --logScaleY  "],
			   ['MVA_Iso_eta_eff',"DecayModeFindingNewDMsEffeta *IsolationMVArun2v1DBnewDMwLTEffeta --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_eta_eff',"DecayModeFindingNewDMsEffeta /*MuonRejection3Effeta --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_phi_eff',"DecayModeFindingNewDMsEffphi /MVA6*ElectronRejectionEffphi  --maxXaxis=200 --minXaxis=-200 --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_phi_eff',"DecayModeFindingNewDMsEffphi *CombinedIsolationDBSumPtCorr3HitsEffphi  --maxXaxis=200 --minXaxis=-200 --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_phi_eff',"DecayModeFindingNewDMsEffphi *IsolationMVArun2v1DBnewDMwLTEffphi  --maxXaxis=200 --minXaxis=-200 --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_phi_eff',"DecayModeFindingNewDMsEffphi /*MuonRejection3Effphi  --maxXaxis=200 --minXaxis=-200 --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiEle_pileup_eff',"DecayModeFindingNewDMsEffpileup /MVA6*ElectronRejectionEffpileup  --maxXaxis=80  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['Comb_Iso_pileup_eff',"DecayModeFindingNewDMsEffpileup *CombinedIsolationDBSumPtCorr3HitsEffpileup  --maxXaxis=80  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['MVA_Iso_pileup_eff',"DecayModeFindingNewDMsEffpileup *IsolationMVArun2v1DBnewDMwLTEffpileup  --maxXaxis=80  --rebin=10 --maxLogY=1.4  --logScaleY "],
			   ['AntiMu_pileup_eff',"DecayModeFindingNewDMsEffpileup /*MuonRejection3Effpileup  --maxXaxis=80  --rebin=10 --maxLogY=1.4  --logScaleY "]],

                ########### TTbar old
		'TTbar':[['Comb_Iso_pt_fake',"DecayModeFindingNewDMsEffpt *CombinedIsolationDBSumPtCorr3HitsEffpt  --rebin=10  --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['MVA_Iso_pt_fake',"DecayModeFindingNewDMsEffpt *IsolationMVArun2v1DBnewDMwLTEffpt  --rebin=10  --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['Comb_Iso_eta_fake',"DecayModeFindingNewDMsEffeta *CombinedIsolationDBSumPtCorr3HitsEffeta  --rebin=10  --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['MVA_Iso_eta_fake',"DecayModeFindingNewDMsEffeta *IsolationMVArun2v1DBnewDMwLTEffeta  --rebin=10  --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['Comb_Iso_phi_fake',"DecayModeFindingNewDMsEffphi *CombinedIsolationDBSumPtCorr3HitsEffphi  --rebin=10  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['MVA_Iso_phi_fake',"DecayModeFindingNewDMsEffphi *IsolationMVArun2v1DBnewDMwLTEffphi  --rebin=10  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['Comb_Iso_pileup_fake',"DecayModeFindingNewDMsEffpileup *CombinedIsolationDBSumPtCorr3HitsEffpileup  --rebin=10  --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"],
				 ['MVA_Iso_pileup_fake',"DecayModeFindingNewDMsEffpileup *IsolationMVArun2v1DBnewDMwLTEffpileup  --rebin=10  --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"]],

                ########### ZEE old
		'ZEE':[['AntiEle_pt_fake',"DecayModeFindingNewDMsEffpt /MVA6*ElectronRejectionEffpt  --rebin=10  --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_eta_fake',"DecayModeFindingNewDMsEffeta /MVA6*ElectronRejectionEffeta  --rebin=10  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_phi_fake',"DecayModeFindingNewDMsEffphi /MVA6*ElectronRejectionEffphi  --rebin=10  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiEle_pileup_fake',"DecayModeFindingNewDMsEffpileup /MVA6*ElectronRejectionEffpileup  --rebin=10  --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"]],

                ########### ZMM old
		'ZMM':[['AntiMu_pt_fake',"DecayModeFindingNewDMsEffpt /*MuonRejection3Effpt  --rebin=10  --maxXaxis=200   --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiMu_eta_fake',"DecayModeFindingNewDMsEffeta /*MuonRejection3Effeta  --rebin=10  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiMu_phi_fake',"DecayModeFindingNewDMsEffphi /*MuonRejection3Effphi  --rebin=10  --maxXaxis=200 --minXaxis=-200  --maxLogY=1.4  --logScaleY --fakeRate"],
			   ['AntiMu_pileup_fake',"DecayModeFindingNewDMsEffpileup /*MuonRejection3Effpileup  --rebin=10  --maxXaxis=80   --maxLogY=1.4  --logScaleY --fakeRate"]]
 	}

                ########### new for all samples
	for evtt, parameters in parameterset.iteritems():
		parameters.append(['pT',"tau_pt  --rebin=2  --normalize  --maxXaxis=200   --maxLogY=1.4  --logScaleY --events"])
		parameters.append(['eta',"tau_eta  --rebin=10  --normalize  --maxLogY=1.4  --logScaleY --events"])
		parameters.append(['phi',"tau_phi  --rebin=10  --normalize  --maxXaxis=200   --maxLogY=1.4  --logScaleY --events"])
		if evtt not in ['ZMM', 'ZEE']: parameters.append(['DeepTauVSjet',"tau_byDeepTau2017v2p1VSjetraw  --rebin=10  --normalize  --maxXaxis=200   --maxLogY=1.4  --logScaleY --events"])
		if evtt not in ['TTbar', 'ZEE']: parameters.append(['DeepTauVSmu',"tau_byDeepTau2017v2p1VSmuraw  --rebin=10  --normalize  --maxXaxis=200   --maxLogY=1.4  --logScaleY --events"])
		if evtt not in ['ZMM', 'TTbar']: parameters.append(['DeepTauVSe',"tau_byDeepTau2017v2p1VSeraw  --rebin=10  --normalize  --maxXaxis=200   --maxLogY=1.4  --logScaleY --events"])

	target = args.paths[0][args.paths[0].rfind("CMSSW"):].rstrip(".root")
	reference = args.paths[1][args.paths[1].rfind("CMSSW"):].rstrip(".root")

	out_folder = target + '_vs_' + reference
	out_folder = out_folder.replace("CMSSW_","").replace("__DQMIO","").replace('_mcRun2_asymptotic_', "").replace('v1', '').replace('v0', '')
	out_folder = os.path.join(args.output_dir, args.tag, out_folder)
	if not os.path.exists(out_folder):
		os.makedirs(out_folder)
	
	commands = []
	#filenameFormat = "DQM_V000([0-9])_R00000000([0-9])__RelVal(ZTT|TenTau|TTbar|ZEE|ZMM)_(13|13_UP18|15_500)__"
	#for filet in os.listdir(os.path.dirname(args.paths[0])):
	#	pattern_tar = re.match(filenameFormat+target+".root", os.path.dirname(args.paths[0])+"/"+filet)
	#	if pattern_tar == None: continue
	#	for filer in os.listdir(os.path.dirname(args.paths[1])):
	#		pattern_ref = re.match(filenameFormat+reference+".root", os.path.dirname(args.paths[1])+"/"+filer)
	#		if pattern_ref == None: continue
	#		if pattern_tar.group(3)!=pattern_ref.group(3) and pattern_tar.group(4)!=pattern_ref.group(4): continue
	#		break
	#	if pattern_ref == None: continue

	for numi in xrange(1,2,1):#xrange(1,ver.nr.+1,1) #For target
		for numj in range(1,2,1):#For reference
			for evtt, parameters in parameterset.iteritems():
				
				
				for params in parameters:
					TEV = "13"
					if "UP18" in args.paths[0]: TEV = "13_UP18"
					if "14TeV__CMSSW" in args.paths[0]: 
						TEV = "14TeV"
					if "14__CMSSW" in args.paths[0]: 
						TEV = "14"
						if evtt=="TTbar" :
							TEV="14TeV"
					if evtt == "TenTau": TEV = "15_500"

					thistarget = os.path.join(os.path.dirname(args.paths[0]), "DQM_V000{numi}_R000000001__RelVal{evtt}_{tev}__{target}.root".format(numi=numi, evtt=evtt, tev=TEV, target=target))
					thisreference = os.path.join(os.path.dirname(args.paths[1]), "DQM_V000{numj}_R000000001__RelVal{evtt}_{tev}__{reference}.root".format(numj=numj, evtt=evtt, tev=TEV, reference=reference))
					
					command = "python $CMSSW_BASE/src/Validation/RecoTau/Tools/MultipleCompare.py -T {target} -R {reference} -t {tfolder} -r {rfolder} {param2} -o {output}".format(
							target=thistarget,
							reference=thisreference,
							tfolder=target,
							rfolder=reference,
							param2=params[1],
							output=os.path.join(out_folder, "{evtt}_13__{target}_{param1}.png".format(evtt=evtt, target=target, param1=params[0]))
					)
					commands.append(command)
	
	
	tools.parallelize(tools.call_command, commands, n_processes=args.n_processes)
	# webplotting.webplotting(input_dir=os.path.join(args.output_dir, args.tag), recursive=True)


if __name__ == "__main__":
	main()
