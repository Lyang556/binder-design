#!/usr/bin/env python3

# Script for designing protein binder
# Using PatchDock > ProteinMPNN > Rosetta
# Sample miniproteins in this repository were originated from the supplementary data in Cao et al., (2022). Nature. 605, 551-560.

import os
import time
import argparse

from tools import dock
from tools import mpnn
from tools import rosetta_relax
from tools import postproc


parser = argparse.ArgumentParser (description='Script for designing protein binder')

# dock
parser.add_argument('--patchdock', required=False, type=str, default='/opt/tools/PatchDock/', help='path to PatchDock')
parser.add_argument('--nstruct', required=False, type=int, default=1, help='nstruct argument for PatchDock')
parser.add_argument('--template', '-t', required=False, type=str, default='./templates/my.pdb', help='template file')

# lists
parser.add_argument('--binderlist', '-l', required=False, type=str, default='./binders.list', help='list of binders for docking')
parser.add_argument('--reslist', '-resl', required=False, type=str, default='./res.list', help='list of residue in receptor for PatchDock')

# dir. setup
parser.add_argument('--binderdir', '-bd',required=False, type=str, default='./binders/', help='path to binders')
parser.add_argument('--dockdir', '-dd', required=False, type=str, default='./01_patchdock', help='Path to prepared input files for docking')
parser.add_argument('--designdir', '-ded', required=False, type=str, default='./02_seq_design_files', help='Path to fasta files of designed sequences')
parser.add_argument('--packerdir', '-pd', required=False, type=str, default='./03_packed_files', help='Path to packed PDB files')
parser.add_argument('--minimdir', '-md', required=False, type=str, default='./04_minimized_files', help='Path to minimized files')
parser.add_argument('--dst', required=False, type=str, default='./08_results_pdb', help='Directory to final results')

# mpnn
# parser.add_argument('--mpnnpath', required=False, type=str, default='', help='path to ProteinMPNN')
# parser.add_argument('--mpnnconda', required=False, type=str, default='proteinmpnn', help='Name of conda env for ProteinMPNN')

# etc
parser.add_argument('--np', required=False, type=str, default='4', help='num cores')
parser.add_argument('--ops', required=False, type=str, default='linux', help='mac or linux')

# postproc
parser.add_argument('--csvname', required=False, type=str, default='./first_score.csv', help='Name of converted csv file')
parser.add_argument('--score', required=False, type=str, default='./final_score.txt', help='Scorefile to be converted to csv format')

# fn. sele
parser.add_argument('--fn', required=False, type=str, default='', help='A function name to use')



args = parser.parse_args()


class RunScript:

    def __init__ (self, args):

        self.args = args


    def fxn_select (self):

        start_init = time.time()

        dock_prep = dock.PatchDock(args)
        docking = dock.Parallel(args)
        seq_design_fn = mpnn.ProteinMPNN(self.args)
        seq_design = mpnn.RunProteinMPNN(self.args)
        rosetta_packer = rosetta_relax.Parallel(args)
        rosetta_minimize = rosetta_relax.RosettaRun(args)
        post_proc = postproc.PostProc(args)


        fn_dict = {
            'listgen':dock_prep.listgen,
            'patchdock':docking.patchdock_pararun,
            'seq_design':seq_design.mpnn_run,
            'packer':rosetta_packer.packer_pararun, 
            'minimize':rosetta_minimize.rosetta_minim,
            'txt2csv':post_proc.txt2csv,
            'pdbsorting':post_proc.pdbsorting
        }

        if args.fn in fn_dict:
            if args.fn == 'packer':
                # parsing mpnn_list that is required for rosetta_packer
                mpnn_list = seq_design_fn.mpnn_list()
            
                fn_ = fn_dict[args.fn]
                fn_(mpnn_list)

            elif args.fn == 'txt2csv':
                fn_ = fn_dict[args.fn]
                fn_(self.args.score)

            else:
                fn_ = fn_dict[args.fn]
                fn_()

            print ('Total elapsed time: ', time.time() - start_init, 'sec')
            quit()

        else:
            fxns = '''
            1) listgen (binderdir, binderlist)
            2) patchdock (patchdock, template, reslist, binderlist, dockdir, np)
            3) seq_design (dockdir, designdir)
            4) packer (ops, designdir, dockdir, packerdir, np)
            5) minimize (ops, packerdir, minimdir, np)
            6) txt2csv (score, csvname)
            7) pdbsorting (csvname, minimdir, designdir, dst)
            '''

            print ('Please select/enter the one of functions below and enter it.')
            print (fxns)
            quit()


    def auto_run (self):

        start_init = time.time ()

        dock_prep = dock.PatchDock(self.args)
        docking = dock.Parallel(self.args)
        seq_design = mpnn.RunProteinMPNN(self.args)
        rosetta_packer = rosetta_relax.Parallel(self.args)
        rosetta_minimize = rosetta_relax.RosettaRun(self.args)
        post_proc = postproc.PostProc(self.args)
        
        # 1st run
        print ('*** 1st run ***')
        dock_prep.listgen()
        docking.patchdock_pararun()
        mpnn_list = seq_design.mpnn_run()
        rosetta_packer.packer_pararun(mpnn_list)
        score_path = rosetta_minimize.rosetta_minim()
        post_proc.txt2csv(score_path)

        # new paths
        self.args.designdir = '05_seq_design_files_2'
        self.args.dockdir = '04_minimized_files'
        self.args.packerdir = '06_packed_files_2'
        self.args.minimdir = '07_minimized_files_2'
        self.args.csvname = './second_score.csv'

        # 2nd run
        print ('*** 2nd run ***')
        mpnn_list = seq_design.mpnn_run()
        rosetta_packer.packer_pararun(mpnn_list)
        score_path = rosetta_minimize.rosetta_minim()
        post_proc.txt2csv(score_path)
        post_proc.pdbsorting()

        print ('Total elapsed time: ', time.time() - start_init, 'sec')
        quit()



if __name__ == '__main__':

    run_script = RunScript(args)

    if args.fn == '':
        pass
    else:
        run_script.fxn_select()

    run_script.auto_run()