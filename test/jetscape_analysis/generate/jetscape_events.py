#!/usr/bin/env python3

"""
  Class to launch the generation of JETSCAPE events over a set of pt-hat bins,
  and optionally over any set of additional JETSCAPE parameter values.
  If additional parameter values are specified, all combinations of the
  specified parameter values will be generated.

  Run from inside the JETSCAPE docker container with:
    python generate_jetscape_events.py -c /home/jetscape-user/JETSCAPE-analysis/config/jetscapeAnalysisConfig.yaml -o /my/outputdir -j /path/to/JETSCAPE

.. codeauthor:: James Mulligan <james.mulligan@berkeley.edu>, UC Berkeley
"""

from __future__ import print_function

import argparse
import fileinput
import os
import shutil
import subprocess
import sys
import yaml
import itertools
import re

# Base class
sys.path.append('../..')
from jetscape_analysis.base import common_base

################################################################
class GenerateJetscapeEvents(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file="", output_dir="", jetscape_dir="", **kwargs):
        super(GenerateJetscapeEvents, self).__init__(**kwargs)
        self.config_file = config_file
        self.output_dir = output_dir
        self.jetscape_dir = jetscape_dir

        # Create output dir
        if not self.output_dir.endswith("/"):
            self.output_dir = self.output_dir + "/"
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.initialize_config()

        print(self)

    # ---------------------------------------------------------------
    # Initialize config file into class members
    # ---------------------------------------------------------------
    def initialize_config(self):

        # Read config file
        with open(self.config_file, "r") as stream:
            config = yaml.safe_load(stream)

        self.debug_level = config["debug_level"]

        self.xml_user_file = config["xml_user_file"]
        self.xml_main_file = config["xml_main_file"]
                
        self.parameter_scan_dict = config['parameter_scan']
        self.pt_hat_bins = self.parameter_scan_dict['pt_hat_bins']['values']

    # ---------------------------------------------------------------
    # Main processing function
    # ---------------------------------------------------------------
    def generate_jetscape_events(self):

        # Store list of parameter labels
        parameter_labels = [self.parameter_scan_dict[key]['label'] for key in self.parameter_scan_dict]

        # Create list of all combinations of parameters
        parameter_values = [self.parameter_scan_dict[key]['values'] for key in self.parameter_scan_dict]
        parameter_combinations = list(itertools.product(*parameter_values))
        
        # Remove that last pt-hat bin edge
        n_combinations_per_pthat = int(len(parameter_combinations)/len(self.pt_hat_bins))
        parameter_combinations = parameter_combinations[:-n_combinations_per_pthat]

        # Loop through all parameter combinations
        for index, parameter_combination in enumerate(parameter_combinations):
        
            pt_hat_bin = int(index / n_combinations_per_pthat)
            if pt_hat_bin < len(self.pt_hat_bins) - 1:
                pt_hat_min = self.pt_hat_bins[pt_hat_bin]
                pt_hat_max = self.pt_hat_bins[pt_hat_bin + 1]
            else:
                continue
            if index % n_combinations_per_pthat == 0:
                print('Generating pt-hat: {} - {} ...'.format(pt_hat_min, pt_hat_max))

            # Create label for output directory
            dir_label = ''
            for index, value in enumerate(parameter_combination):
                if index == 0:
                    dir_label += str(pt_hat_bin)
                    continue
                dir_label += '_'
                dir_label += parameter_labels[index]
                dir_label += str(value)
            if len(parameter_combination) > 1:
                print('    Generating {}'.format(dir_label))
                
            # Create outputDir for each bin
            output_dir_bin = '{}{}'.format(self.output_dir, dir_label)
            if not output_dir_bin.endswith("/"):
                output_dir_bin = output_dir_bin + "/"
            if not os.path.exists(output_dir_bin):
                os.makedirs(output_dir_bin)

            # Copy XML files to pt-hat bin directory
            xml_main_file_copy = "{}{}".format(output_dir_bin, "jetscape_main.xml")
            cmd = "rsync {} {}".format(self.xml_main_file, xml_main_file_copy)
            os.system(cmd)

            xml_user_file_copy = "{}{}".format(output_dir_bin, "jetscape_user.xml")
            cmd = "rsync {} {}".format(self.xml_user_file, xml_user_file_copy)
            os.system(cmd)
                    
            # Set parameters in the Jetscape User XML configuration
            for index, value in enumerate(parameter_combination):
            
                parameter_label = parameter_labels[index]
            
                # Set pt-hat
                if parameter_label == 'pt_hat_bins':
                
                    for line in fileinput.input(xml_user_file_copy, inplace=True):
                
                        if 'pTHatMin' in line:
                            print(re.sub(r'{0}>\w*</{0}'.format('pTHatMin'), '{0}>{1}</{0}'.format('pTHatMin', pt_hat_min), line), end='')
                        elif 'pTHatMax' in line:
                            print(re.sub(r'{0}>\w*</{0}'.format('pTHatMax'), '{0}>{1}</{0}'.format('pTHatMax', pt_hat_max), line), end='')
                        else:
                            print(line, end='')
                      
                # Set other parameters
                else:
                
                    for line in fileinput.input(xml_user_file_copy, inplace=True):
                    
                        if parameter_label in line:
                            print(re.sub(r'{0}>\w*</{0}'.format(parameter_label), '{0}>{1}</{0}'.format(parameter_label, value), line), end='')
                        else:
                            print(line, end='')

            # cd into bin directory in order to write Jetscape output there
            os.chdir(output_dir_bin)

            # Run Jetscape executable
            logfile_name = os.path.join(output_dir_bin, "log_{}.txt".format(dir_label))
            with open(logfile_name, "w") as logfile:
                cmd = '{}/build/runJetscape jetscape_user.xml jetscape_main.xml'.format(self.jetscape_dir)
                subprocess.run(cmd, check=True, shell=True, stdout=logfile)
                
            os.chdir(self.output_dir)

##################################################################
if __name__ == "__main__":
    # Define arguments
    parser = argparse.ArgumentParser(description="Generate JETSCAPE events")
    parser.add_argument(
        "-c",
        "--configFile",
        action="store",
        type=str,
        metavar="configFile",
        default="/home/jetscape-user/JETSCAPE/test/jetscape_analysis/config/example.yaml",
        help="Path of config file for analysis",
    )
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        type=str,
        metavar="outputDir",
        default="/home/jetscape-user/JETSCAPE/test/TestOutput",
        help="Output directory for output to be written to",
    )
    parser.add_argument(
        "-j",
        "--jetscapeDir",
        action="store",
        type=str,
        metavar="jetscapeDir",
        default="/home/jetscape-user/JETSCAPE",
        help="Location of JETSCAPE repository",
    )

    # Parse the arguments
    args = parser.parse_args()

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File "{0}" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    analysis = GenerateJetscapeEvents(config_file=args.configFile, output_dir=args.outputDir, jetscape_dir=args.jetscapeDir)
    analysis.generate_jetscape_events()
