#!/usr/bin/env python3

import multiprocessing
import os
import subprocess



class RosettaRun:

    def __init__ (self, args):

        self.args = args

    def rosetta_packer (self, mpnn_seq):

        designdir = os.path.abspath(self.args.designdir)
        dockdir = os.path.abspath(self.args.dockdir)
        packerdir = os.path.abspath(self.args.packerdir)

        ops = self.args.ops.lower()

        if ops == 'mac':
            fixbb = 'fixbb.mpi.macosclangrelease'
        else:
            fixbb = 'fixbb.mpi.linuxgccrelease'

        try:
            os.mkdir(packerdir)
        except FileExistsError:
            pass


        # Prep. resfile
        with open (f'{designdir}/seqs/{mpnn_seq}', 'r') as fasta:
            seqs = fasta.readlines()[3].rstrip()
            seqs = list(seqs)
            mpnn_seq_f_name = mpnn_seq.replace('.fa','')
            with open (f'{mpnn_seq_f_name}_resfile.txt', 'w') as resfile:
                resfile.write ('NATRO' + '\n')
                resfile.write ('start' + '\n')            
                for num, seq in enumerate (seqs):
                    resfile.writelines (str(num + 1) + f' A PIKAA {seq}' + '\n')
                
                resfile.close()
            fasta.close()

        mpnn_seq_to_pdb = mpnn_seq.replace ('.fa', '.pdb')

        cmd = [fixbb, '-in:file:s', f'{dockdir}/{mpnn_seq_to_pdb}', '-in:file:fullatom', '-resfile', f'{mpnn_seq_f_name}_resfile.txt', '-out:path:all', packerdir, '-out:suffix', '_packer']
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        os.remove (f'{mpnn_seq_f_name}_resfile.txt')


    def rosetta_packer_postproc (self):

        packerdir = os.path.abspath(self.args.packerdir)
        packed_files = os.listdir(packerdir)
        packed_list_path = [file for file in packed_files if file.endswith('.pdb')]

        lsfile = open (f'{packerdir}/packer_pdb.list', 'a')

        for i in packed_list_path:
            lsfile.write (f'{packerdir}/{i}' + '\n')

        lsfile.close()


    def rosetta_minim (self) -> str:
        
        packerdir = os.path.abspath(self.args.packerdir)
        minimdir = os.path.abspath(self.args.minimdir)

        xmlStr = f"""
    <ROSETTASCRIPTS>

        <SCOREFXNS>
            <ScoreFunction name="beta_nov16" weights="beta_nov16" />
            
        </SCOREFXNS>

        <RESIDUE_SELECTORS>
        </RESIDUE_SELECTORS>

        <TASKOPERATIONS>
        </TASKOPERATIONS>

        <FILTERS>
        </FILTERS>

        <MOVERS>
            <AtomTree name="atomtree" docking_ft="true" />

            <MinMover name='minimize' jump='1' scorefxn='beta_nov16' chi='true' bb='false' bondangle='false' bondlength='false' cartesian='0' tolerance='0.001' />

            <InterfaceAnalyzerMover name="int_ana" packstat="true" interface_sc="true" interface="A_B" pack_separated="true" />
        </MOVERS>

        <PROTOCOLS>
        
            <Add mover="atomtree" />
            <Add mover="minimize" />
            <Add mover="int_ana" />
            
        </PROTOCOLS>

    </ROSETTASCRIPTS>
    """

        try:
            os.mkdir(minimdir)
        except FileExistsError:
            pass

        f = open(f'{minimdir}/minimize.xml', 'w')
        f.write(xmlStr)
        f.close()

        ops = self.args.ops.lower()

        if ops == 'mac':
            rosetta_scripts = 'rosetta_scripts.mpi.macosclangrelease'
        else:
            rosetta_scripts = 'rosetta_scripts.mpi.linuxgccrelease'

        num_core = self.args.np

        # Run command
        cmd = \
            ['mpirun', 
            '--use-hwthread-cpus', 
            rosetta_scripts,
            '-in:file:l',
            f'{packerdir}/packer_pdb.list',
            '-out:path:all',
            minimdir,
            '-out:suffix',
            '_minimize',
            '-parser:protocol',
            f'{minimdir}/minimize.xml',
            '-out:file:scorefile', 
            'score_minim.txt',
            '-corrections::beta_nov16',
            'true'
            ]

        print ('[Rosetta] Minimize running ...')
        if num_core == 'all':
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        else:
            cmd[1] = '-np'
            cmd.insert(2, str(num_core))
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        os.remove(f'{minimdir}/minimize.xml')

        score_path = os.path.join(minimdir, 'score_minim.txt')

        return score_path



class Parallel:

    def __init__ (self, args):

        self.args = args

    def packer_pararun(self, mpnn_list):

        if self.args.np == 'all':
            nc = os.cpu_count() - 1
        else:
            nc = self.args.np

        print ('[Rosetta] Packer running ...')
        p = multiprocessing.Pool(processes=int(nc))
        p.map(RosettaRun(self.args).rosetta_packer, mpnn_list)
        p.close()

        RosettaRun(self.args).rosetta_packer_postproc()