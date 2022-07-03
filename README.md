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

### 4. Working directory setting

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

*flags*

```
usage: protein_binder_design_ver2.py [-h] [--template TEMPLATE]
                                     [--binderdir BINDERDIR] [--np NP]
                                     [--ops OPS] [--nstruct NSTRUCT]
                                     [--minilist MINILIST] [--reslist RESLIST]
                                     [--danpath DANPATH] [--mpnnpath MPNNPATH]
                                     [--patchdockpath PATCHDOCKPATH]
                                     [--prodigyconda PRODIGYCONDA]
                                     [--danconda DANCONDA]
                                     [--mpnnconda MPNNCONDA] [--inpdir INPDIR]
                                     [--paramdir PARAMDIR] [--dockdir DOCKDIR]
                                     [--designdir DESIGNDIR]
                                     [--packerdir PACKERDIR]
                                     [--minimdir MINIMDIR] [--dandir DANDIR]
                                     [--danresult DANRESULT]
                                     [--sorting SORTING]
                                     [--design_after_design DESIGN_AFTER_DESIGN]
                                     [--score SCORE] [--csvname CSVNAME]
                                     [--dist DIST] [--prodigyrun PRODIGYRUN]
                                     [--fn FN] [--postsele POSTSELE]

Script for designing protein binder

optional arguments:
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
  --inpdir INPDIR       saving path to prepared input files for docking
  --paramdir PARAMDIR   saving path to prepared input files for docking
  --dockdir DOCKDIR     saving path to docking files
  --designdir DESIGNDIR
                        saving path to designed sequences
  --packerdir PACKERDIR
                        saving path to packed PDB files
  --minimdir MINIMDIR   saving path to minimized files
  --dandir DANDIR       saving path to binder files for DeepAccNet
  --danresult DANRESULT
                        Name of result csv file of DeepAccNet
  --sorting SORTING     Sort by ...
  --design_after_design DESIGN_AFTER_DESIGN
                        Run design twice
  --score SCORE         Scorefile to be converted to csv format
  --csvname CSVNAME     Name of converted csv file
  --dist DIST           directory to final results
  --prodigyrun PRODIGYRUN
                        Run PRODIGY or not (true/false)
  --fn FN               type a function name to use
  --postsele POSTSELE

```

<br>
<br>

```
python3 protein_binder_design.py \
    --binderdir ./binders/
    --template ./template/template.pdb
    --mpnnpath /opt/tools/ProteinMPNN/
    --patchdockpath /opt/tools/PatchDock/
    --danpath /opt/tools/DeepAccNet/
    --np 8
    --ops mac
    --nstruct 2
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
        pdbsorting (csvname, minimdir, dist, postsele, sorting, designdir)
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