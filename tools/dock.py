#!/usr/bin/env python3

import os
import time
import multiprocessing
import subprocess
import shutil

class PatchDock:

    def __init__ (self, args):

        self.args = args

    def listgen (self): 
        
        binderdir = os.path.abspath (self.args.binderdir)
        binderdir_path = os.listdir(binderdir)
        miniproteins = [file for file in binderdir_path if file.endswith('.pdb')]

        binderlist = os.path.abspath (self.args.binderlist)

        if os.path.isfile(binderlist):
            print ('Pre-existed list file will be deleted, and new list file will be created')
            time.sleep(3.9)
            os.remove(binderlist)
        else:
            pass

        listfile = open (f'{binderlist}', 'a')

        for i in miniproteins:
            listfile.write(f'{binderdir}/{i}' + '\n')


    def dockprep (self) -> list:

        patchdock = os.path.abspath(self.args.patchdock)
        template = os.path.abspath(self.args.template)
        reslist = os.path.abspath(self.args.reslist)
        binderlist = os.path.abspath(self.args.binderlist)
        dockdir = os.path.abspath(self.args.dockdir)

        try:
            os.mkdir(dockdir)
        except FileExistsError:
            pass

        ls = open(binderlist)

        for i in ls:
            j = i.split('/')[-1].replace('.pdb','').rstrip()
            binder_protein = os.path.abspath(i).rstrip()
            with open (f'{dockdir}/{j}.params', 'w') as param:
                param.write(f'receptorPdb {template}' + '\n')
                param.write(f'ligandPdb {binder_protein}' + '\n')
                param.write(f'protLib {patchdock}/chem.lib' + '\n')
                param.write(f'log-file {j}.log' + '\n')
                param.write(f'log-level 0' + '\n')
                param.write(f'receptorSeg 10.0 20.0 1.5 1 0 1 0' + '\n')
                param.write(f'ligandSeg 10.0 20.0 1.5 1 0 1 0' + '\n')
                param.write(f'scoreParams 0.3 -5.0 0.5 0.0 0.0 1500 -8 -4 0 1 0' + '\n')
                param.write(f'desolvationParams 500.0 1.0' + '\n')
                param.write(f'clusterParams 0.1 4 2.0 3.0' + '\n')
                param.write(f'baseParams 4.0 13.0 2' + '\n')
                param.write(f'matchingParams 1.5 1.5 0.4 0.5 0.9' + '\n')
                param.write(f'matchAlgorithm 1' + '\n')
                param.write(f'receptorGrid 0.5 6.0 6.0' + '\n')
                param.write(f'ligandGrid 0.5 6.0 6.0' + '\n')
                param.write(f'receptorMs 10.0 1.8' + '\n')
                param.write(f'ligandMs 10.0 1.8' + '\n')
                param.write(f'receptorActiveSite {reslist}' + '\n')

                param.close()
        
        cmd_ls = []

        for i in os.listdir(f'{dockdir}'):
            outs = i.replace('.params','.out')
            cmd = [f'{patchdock}/patch_dock.Linux', f'{dockdir}/{i}', f'{dockdir}/{outs}']
            cmd_ls.append(cmd)

        return cmd_ls


    def pararun_func (self, cmd):

        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


    def gen_pdb (self) -> list:

        patchdock = os.path.abspath(self.args.patchdock)
        nstruct = self.args.nstruct

        dockdir = os.path.abspath(self.args.dockdir)
        dockdir_name = os.path.basename(dockdir)

        out_ls = [file for file in os.listdir(dockdir) if file.endswith('.out')]

        all_ls = [file for file in os.listdir(dockdir)]

        cmd_ls = []
        for i in out_ls:
            cmd = ['perl', f'{patchdock}/transOutput.pl', str(i), '1', str(nstruct)]
            cmd_ls.append(cmd)

        return cmd_ls, all_ls, dockdir_name


class Parallel:

    def __init__ (self, args):

        self.args = args


    def patchdock_pararun (self):

        cmd_ls = PatchDock(self.args).dockprep()

        if self.args.np == 'all':
            nc = os.cpu_count() - 1
        else:
            nc = self.args.np

        dockdir = os.path.abspath(self.args.dockdir)
        runpath = os.path.dirname(dockdir)

        print ('[PatchDock] Docking ...')
        os.chdir(dockdir)        
        p = multiprocessing.Pool(processes=int(nc))
        p.map(PatchDock(self.args).pararun_func, cmd_ls)
        p.close()
        os.chdir(runpath)
        
        print ('[PatchDock] PDB generating ...')
        cmd_ls, all_ls, dockdir_name = PatchDock(self.args).gen_pdb()
        os.chdir(dockdir)
        p = multiprocessing.Pool(processes=int(nc))
        p.map(PatchDock(self.args).pararun_func, cmd_ls)
        p.close()
        os.chdir(runpath)

        result_dir = os.path.join (dockdir, dockdir_name + '_result')

        os.mkdir (result_dir)
        for i in all_ls:
            shutil.move(f'{dockdir}/{i}', result_dir) 



