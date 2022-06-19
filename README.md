# Protein binder design 

Script to design protein binder using Rosetta and ProteinMPNN.

## Requirements

```
Rosetta
ProteinMPNN (in conda environment)

Pytorch
Pandas
Numpy
python-dateutil
PyMOL
```

- - -

## *How to use*

<br>


### 1. Env setup

#### conda

```
conda create -n binderdesign python=3.7.0

conda install -c pytorch pytorch=1.11.0
conda install -c anaconda numpy pandas python-dateutil
conda install -c schrodinger pymol

```

#### ProteinMPNN

```
git clone https://github.com/dauparas/ProteinMPNN

```

and then move the ProteinMPNN directory to somewhere (eg. /opt/tools/ProteinMPNN/)


### 1. Template preparation

'Template' is a PDB file containing a target and a binder. 

Considering how rosetta docking works, once you have selected the surface on your target for the binding of binder, a binder should be properly positioned around that surface. 

In the template PDB file, binder should be set to chain A and target to chain B, Since the PyMOL could properly visualize your structure regardless of the order of chain arrangement, you should make sure that the chain order is sorted from A to B even in the PDB file opened in the text editor such as vim.

<br>

### 2. Binder library

You can find the library of miniproteins (binders) in the supplementary data of *Cao_2022_Nature* paper.

<br>

### 3. Working directory setting

```

    Working_dir / 

        template / 
            template.pdb

        binders / 
            binder_1.pdb
            binder_2.pdb
            binder_3.pdb
            ...

        protein_binder_design.py


```

<br>

### 4. run the script

flags

```
parser.add_argument('--binderdir', required=False, type=str, default='./miniproteins/', help='path to miniproteins')
parser.add_argument('--minilist', required=False, type=str, default='./miniproteins.list', help='list file' )
parser.add_argument('--template', required=False, type=str, default='./templates/my.pdb', help='template file')
parser.add_argument('--inpdir', required=False, type=str, default='./01_inputs', help='saving path to prepared input files')
parser.add_argument('--dockdir', required=False, type=str, default='./02_docking_files', help='saving path to docking files')
parser.add_argument('--designdir', required=False, type=str, default='./03_seq_design_files', help='saving path to designed sequences')
parser.add_argument('--packerdir', required=False, type=str, default='./04_packed_files', help='saving path to packed PDB files')
parser.add_argument('--minimdir', required=False, type=str, default='./05_minimized_files', help='saving path to minimized files')


parser.add_argument('--fn', required=False, type=str, default='', help='type a function name to use')
parser.add_argument('--design_after_design', required=False, type=str, default='true', help='Run design twice')


parser.add_argument('--mpnnpath', required=False, type=str, default='/opt/utilities/ProteinMPNN/', help='path to ProteinMPNN')
parser.add_argument('--mpnnconda', required=False, type=str, default='proteinmpnn', help='Name of conda env for ProteinMPNN')
parser.add_argument('--prodigyrun', required=False, type=str, default='false', help='Run PRODIGY or not (true/false)')
parser.add_argument('--prodigyconda', required=False, type=str, default='prodigy', help='Name of conda env for PRODIGY')


parser.add_argument('--np', required=False, type=str, default='4', help='num cores')
parser.add_argument('--ops', required=False, type=str, default='mac', help='mac or linux')
parser.add_argument('--nstruct', required=False, type=int, default=3, help='nstruct argument for rosetta')


parser.add_argument('--score', required=False, type=str, default='./final_score.txt', help='Scorefile to be converted to csv format')
parser.add_argument('--csvname', required=False, type=str, default='./final_score.csv', help='Name of converted csv file')
parser.add_argument('--dist', required=False, type=str, default='./09_results_pdb', help='directory to final results')
parser.add_argument('--postsele', required=False, type=str, default='false', help='')
```

<br>
<br>

```
python3 protein_binder_design.py \
    --binderdir ./binders/
    --template ./template/template.pdb
    --mpnnpath /opt/utilities/ProteinMPNN/
    --mpnnconda binderdesign
    --np 8
    --ops mac
```

The script automatically run the rosetta docking,

and `seq. design` > `rosetta packing` > `rosetta minimize` will be then repeated twice.

You can find 9 new folders in your working directory, and several csv files containing scores of rosetta minimize and ProteinMPNN.


<br>

### etc

If you want to use only one function, not entire process,

try to use `--fn` flag.

```
# check the function list and require arguments

python3 protein_binder_design.py --fn x

        listgen (binderdir, minilist)
        inputprep (template, minilist, inpdir)
        docking (inpdir, dockdir, np, ops, nstruct)
        design (mpnnpath, mpnnconda, designdir, dockdir)
        packer (designdir, dockdir, packerdir, np, ops)
        minim_intana (packerdir, minimdir, np, ops)
        txt2csv (score, csvname, designdir)
        pdbsorting (csvname, minimdir, dist, postsele)
        ((prodigy (minimdir, csvname, prodigyconda)))



# example

python3 protein_binder_design.py \
    --fn txt2csv 
    --score ./final_score.txt 
    --csvname ./final_score_csv.csv 
    --designdir ./03_seq_design_files/
```


You can utilize it when you want to change docking protocol that you want (such as Rifdock?).





