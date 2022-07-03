import os
import argparse

parser = argparse.ArgumentParser (description='Script for designing protein binder')

parser.add_argument('--binderdir', required=False, type=str, default='./binders/', help='path to binders')
parser.add_argument('--ops', required=False, type=str, default='mac', help='path to binders')
parser.add_argument('--np', required=False, type=str, default='4', help='num cores')

args = parser.parse_args()


def polyV (binder_dir=args.binderdir, ops=args.ops, np=args.np):

    if np == 'all':
        np = (os.cpu_count()) -1
    else:
        pass

    polyV_outputdir = 'polyV_binders'

    if polyV_outputdir in os.listdir('./'):
        pass

    else:
        os.mkdir(polyV_outputdir)

    polyV_flag="""<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="sfxn" weights="ref2015" />

    </SCOREFXNS>
    <TASKOPERATIONS>

         <RestrictAbsentCanonicalAAS name="polyV" keep_aas="V" />

    </TASKOPERATIONS>
    <MOVERS>
        <PackRotamersMover name="PackRotamers" scorefxn="sfxn" task_operations="polyV" /> 
    </MOVERS>
    <PROTOCOLS>

        <Add mover="PackRotamers" />

    </PROTOCOLS>
    <OUTPUT />
</ROSETTASCRIPTS>
"""

    with open (f'./{polyV_outputdir}/polyV.xml', 'a') as file:
        file.writelines(polyV_flag)
        file.close()

    print ('-----------------------------')
    print ('polyV-ize running ...')

    os.system(f"find {binder_dir} -name '*.pdb' > {binder_dir}/binders.list")

    if ops == 'mac':
        mpicmd = 'rosetta_scripts.mpi.macosclangrelease'
    elif ops == 'linux':
        mpicmd = 'rosetta_scripts.mpi.linuxgccrelease'

    os.system(f'mpirun -np {np} {mpicmd} -parser:protocol {polyV_outputdir}/polyV.xml -beta_nov16 -l {binder_dir}/binders.list -mute protocols.rosetta_scripts.ParsedProtocol.REPORT -mute protocols.rosetta_scripts.ParsedProtocol.REPORT -out:path:all ./{polyV_outputdir} > polyV_log.log 2>&1')
    os.system(f"find $(pwd)/{polyV_outputdir} -name '*.pdb' > {polyV_outputdir}/polyV_scaffolds.list")
    os.system(f'rm -rf ./{polyV_outputdir}/polyV.xml')


if __name__ == '__main__':
    polyV (binder_dir=args.binderdir, ops=args.ops, np=args.np)