#!/usr/bin/env python

import os
import sys
import optparse
import shutil

import Validation.tau_validation_tools.tools as tools


def readInput():
    parser = optparse.OptionParser( description='Get the names of the input files from web.', \
                                    usage='usage: %prog [options] identifier (e.g. CMSSW_8_1_0_pre5) ' )
    parser.add_option( '-f', '--filter', action='store', metavar='filter', default="", dest="filter",
                       help='Additional filter on the filenames. [default: %default]' )
    parser.add_option( '-t', '--tag', action='store', metavar='tag', default="", dest="tag",
                       help='Additional name-tag for the folder. [default: %default]' )
    parser.add_option( '-v', '--veto', action='store', metavar='veto', default="", dest="veto",
                       help='Additional veto on the filenames. [default: %default]' )
    parser.add_option( '-o', '--output_dir', action='store', metavar='output_dir', default="/disk1/knutzen/TauPOG/RelVal/samples/", dest="output_dir",
                       help='Output directory. [default: %default]' )

    ( options, args ) = parser.parse_args()

    if len( args ) > 1:
        print "Error: Please state exactly 1 identifier."
        sys.exit( 1 )

    return options, args[ 0 ]

def copyFiles(exe, filelist, identifier, identifier_folder, options):
    outDir_parent = str( options.output_dir )
    outDir = outDir_parent + "_" + identifier + "_" + str( options.tag )
    if os.path.isdir( outDir ):
        print ""
        input_str_outDir = raw_input( outDir + " already exists. It will be deleted. Do you want to continue? ([n] y): " )
        if input_str_outDir == "y":
           pass
        else:
           sys.exit( 0 )

    shutil.rmtree( outDir, ignore_errors = True )
    os.makedirs(outDir)

    for filename in filelist:
        exe_temp = exe + filename + " -o " + outDir + "/" + filename
        print exe_temp
        tools.call_command(exe_temp)

def main():
    converted_key = "userkey.pem"
    tools.convert_certificate_key_rsa("~/.globus/userkey.pem", converted_key)

    options, identifier = readInput()

    if not "CMSSW" in identifier:
        print "The identifier has to contain the string 'CMSSW'"
        sys.exit( 1 )
    cmssw_version_raw_list = identifier.split( "CMSSW_" )[ -1 ].split( "_" )

    if len( cmssw_version_raw_list ) == 3:
        num0 = cmssw_version_raw_list[ 0 ]
        num1 = cmssw_version_raw_list[ 1 ]
        num2 = cmssw_version_raw_list[ 2 ]
    elif len( cmssw_version_raw_list ) == 4:
        num0 = cmssw_version_raw_list[ 0 ]
        num1 = cmssw_version_raw_list[ 1 ]
        num2 = cmssw_version_raw_list[ 2 ]
    else:
        print "Error: Cannot parse identifier."
        print "The identifier has to match the pattern: CMSSW_Z_Y_X(_A)"
        print "E.g.: CMSSW_8_1_0_pre5 or CMSSW_8_1_0 "
        sys.exit( 1 )

    identifier_folder = "CMSSW_" + num0 + "_" + num1 + "_x/"
    print identifier_folder

    exe = ("curl --anyauth --cert-type PEM --cert ~/.globus/usercert.pem --key {user_key} --key-type PEM -k " +
           "https://cmsweb.cern.ch/dqm/relval/data/browse/ROOT/RelVal/{identifier_folder}").format(user_key=converted_key, identifier_folder=identifier_folder)

    print exe

    ( curl_out, curl_err ) = tools.call_command(exe)

    print curl_err

    filename_list = []

    #print curl_out

    output_list = curl_out.split( "<tr><td><a href='" )
    for line in output_list:
        #if ".root'>" not in line or identifier not in line \
        if ".root'>" not in line or identifier not in line or str( options.filter ) not in line or ( str( options.veto ) in line and str( options.veto ) ) \
            or not ( "ZEE" in line or "ZMM" in line or "ZTT" in line or "QCD" in line or "TTbar" in line or "Zprime" in line ):
            continue
        filename = line.split( ".root'>" )[ -1 ].split( "</a></td><td>" )[ 0 ]
        print filename
        filename_list.append( filename )

    input_files_decision = raw_input( "Do you want to copy the above stated files to your local machine? ([n] y): " )
    if input_files_decision == "y":
       pass
    else:
       sys.exit( 0 )

    copyFiles(exe, filename_list, identifier, identifier_folder, options)
    
    if os.path.exists(converted_key):
        os.remove(converted_key)

if __name__ == '__main__':
    main()
