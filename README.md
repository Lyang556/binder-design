# Protein binder design 

Script to design protein binder using Rosetta and ProteinMPNN.

## Requirements

```
Required

    Rosetta
    PatchDock

    Conda env
        * name: binderdesign
            - Numpy
            - Pandas
            - PyMOL

        * name: ProteinMPNN
            - ProteinMPNN

        * name: dan
            - DeepAccNet
```

- - -

## *How to use*

<br>


### 1. Env setup

#### conda envs

```
conda create -n binderdesign python=3.7.0
conda install -c anaconda numpy pandas
conda install -c schrodinger pymol
```

for ProteinMPNN and DeepAccNet (dan), please find the related github page.

*[ProteinMPNN](https://github.com/dauparas/ProteinMPNN)*
<br>

*[DeepAccNet](https://github.com/hiranumn/DeepAccNet)*

<br>

#### ProteinMPNN

```
git clone https://github.com/dauparas/ProteinMPNN
```

and then move the ProteinMPNN directory to somewhere (eg. /opt/tools/ProteinMPNN/)

<br>

#### Rosetta

*https://www.rosettacommons.org*

It is required to compile mpi version of rosetta especially in the case of rosetta_script and fixbb.

In addition, the name of binary files (rosetta_script and fixbb) have to be   

`rosetta_scripts.mpi.macosclangrelease` and `fixbb.mpi.macosclangrelease` in macOS,

and `rosetta_scripts.mpi.linuxgccrelease` and `fixbb.mpi.linuxgccrelease` in linux.

<br>

### 2. Template preparation

'Template' is a PDB file of target whose chain set to 'B'.

<br>

### 3. Binder library

You can find the library of miniproteins (binders) in the supplementary data of *Cao_2021_Nature* paper.

<br>

### 4. binding region preparation

You have to create a file containing list of residue (first column) and chain ID (second column) as below.

Since the chain of template PDB file was assigned as B, the second column have to be set as B.

```
10 B
11 B
12 B
20 B
21 B
22 B
186 B
187 B
188 B
189 B
190 B
191 B
192 B
193 B
```

You can use the script in `prep_res_list.py` in `etc` directory to create a list file automatically.

```
python3 prep_res_list.py --residues 10-12+20-22+186-193 --chain B
```

<br>

### 5. Working directory setting

```

    Working_dir / 

        template / 
            template.pdb
            res.list

        binders / 
            binder_1.pdb
            binder_2.pdb
            binder_3.pdb
            ...

        protein_binder_design.py


```

<br>

### 6. run the script

*flags*

```
usage: protein_binder_design.py [-h] [--template TEMPLATE]
                                [--binderdir BINDERDIR] [--np NP] [--ops OPS]
                                [--nstruct NSTRUCT] [--minilist MINILIST]
                                [--reslist RESLIST] [--danpath DANPATH]
                                [--mpnnpath MPNNPATH]
                                [--patchdockpath PATCHDOCKPATH]
                                [--prodigyconda PRODIGYCONDA]
                                [--danconda DANCONDA] [--mpnnconda MPNNCONDA]
                                [--paramdir PARAMDIR] [--dockdir DOCKDIR]
                                [--designdir DESIGNDIR]
                                [--packerdir PACKERDIR] [--minimdir MINIMDIR]
                                [--dandir DANDIR] [--danresult DANRESULT]
                                [--sorting SORTING]
                                [--design_after_design DESIGN_AFTER_DESIGN]
                                [--score SCORE] [--csvname CSVNAME]
                                [--dst DST] [--prodigyrun PRODIGYRUN]
                                [--fn FN] [--postsele POSTSELE]
```

<br>
<br>

```
options:
  -h, --help            show this help message and exit
  --template TEMPLATE   template file
  --binderdir BINDERDIR
                        path to binders
  --np NP               num cores
  --ops OPS             mac or linux
  --nstruct NSTRUCT     nstruct argument for PatchDock
  --minilist MINILIST   list file
  --reslist RESLIST     residue list file for PatchDock
  --danpath DANPATH     Path to DeepAccNet
  --mpnnpath MPNNPATH   path to ProteinMPNN
  --patchdockpath PATCHDOCKPATH
                        path to PatchDock
  --prodigyconda PRODIGYCONDA
                        Name of conda env for PRODIGY
  --danconda DANCONDA   Name of conda env for DeepAccNet
  --mpnnconda MPNNCONDA
                        Name of conda env for ProteinMPNN
  --paramdir PARAMDIR   Path to prepared input files for docking
  --dockdir DOCKDIR     Path to result files of docking
  --designdir DESIGNDIR
                        Path to fasta files of designed sequences
  --packerdir PACKERDIR
                        Path to packed PDB files
  --minimdir MINIMDIR   Path to minimized files
  --dandir DANDIR       Path to binder files for DeepAccNet
  --danresult DANRESULT
                        Name of result csv file of DeepAccNet
  --sorting SORTING     Sort by ...
  --design_after_design DESIGN_AFTER_DESIGN
                        Run design twice
  --score SCORE         Scorefile to be converted to csv format
  --csvname CSVNAME     Name of converted csv file
  --dst DST             Directory to final results
  --prodigyrun PRODIGYRUN
                        Run PRODIGY or not (true/false)
  --fn FN               Type a function name to use
  --postsele POSTSELE   Post-selection (after finished) (true/false)

```

<br>
<br>

```
python3 protein_binder_design.py \
    --binderdir ./binders/
    --template ./template/template.pdb
    --reslist ./template/res.list
    --mpnnpath /opt/tools/ProteinMPNN/
    --patchdockpath /opt/tools/PatchDock/
    --danpath /opt/tools/DeepAccNet/
    --np all
    --ops linux
    --nstruct 1
```

The script automatically run the `PatchDock`

and `seq. design` > `rosetta packing` > `rosetta minimize` will be then repeated twice.

You can find 10 new folders in your working directory, and several csv files containing scores.


<br>

### etc

If you want to use only one function, not entire process,

try to use `--fn` flag.

```
# check the function list and require arguments

python3 protein_binder_design.py --fn x

        listgen (binderdir, minilist)
        patch_dock (template, minilist, paramdir, patchdockpath, reslist)
        gen_pdb(paramdir, dockdir, patchdockpath, nstruct)
        design (mpnnpath, mpnnconda, designdir, dockdir)
        packer (designdir, dockdir, packerdir, np, ops)
        minim_intana (packerdir, minimdir, np, ops)
        txt2csv (score, csvname, designdir)
        pdbsorting (csvname, minimdir, dst, postsele, sorting, designdir)
        dan (minimdir, dandir, danpath, np, danconda)
        ((prodigy (minimdir, csvname, prodigyconda)))
        ((score_dan_merge (csvname, danresult)))


# example

python3 protein_binder_design.py \
    --fn txt2csv 
    --score ./final_score.txt 
    --csvname ./final_score_csv.csv 
    --designdir ./03_seq_design_files/
```

You can utilize it when you want to change docking protocol that you want (such as Rifdock / zdock ...).