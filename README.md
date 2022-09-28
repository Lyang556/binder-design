# Protein binder design 

Script to design protein binder using Rosetta and ProteinMPNN.

## Requirements

```
Required

    Rosetta
    PatchDock
    Conda
        - PyTorch
        - Numpy
        - Pandas
```

- - -

## *How to use*

<br>


### 1. Env setup

#### conda envs - dependencies

```
Python>=3.0
PyTorch
Numpy
Pandas
```

<br>

for ProteinMPNN, please find the related github page.

*[ProteinMPNN](https://github.com/dauparas/ProteinMPNN)*
<br>


#### ProteinMPNN

```
git clone https://github.com/dauparas/ProteinMPNN
```

The modified version of ProteinMPNN is already included in this repository (tools directory).

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

        tools /

        template / 
            template.pdb
            res_chain.list

        binders / 
            binder_1.pdb
            binder_2.pdb
            binder_3.pdb
            ...

        binder_design.py


```

<br>

### 6. run the script

*flags*

```
usage: binder_design.py [-h] [--patchdock PATCHDOCK] [--nstruct NSTRUCT] [--template TEMPLATE]
                        [--binderlist BINDERLIST] [--reslist RESLIST] [--binderdir BINDERDIR]
                        [--dockdir DOCKDIR] [--designdir DESIGNDIR] [--packerdir PACKERDIR]
                        [--minimdir MINIMDIR] [--dst DST] [--np NP] [--ops OPS]
                        [--csvname CSVNAME] [--score SCORE] [--fn FN]
```

<br>
<br>

```
options:
  -h, --help            show this help message and exit
  --patchdock PATCHDOCK
                        path to PatchDock
  --nstruct NSTRUCT     nstruct argument for PatchDock
  --template TEMPLATE, -t TEMPLATE
                        template file
  --binderlist BINDERLIST, -l BINDERLIST
                        list of binders for docking
  --reslist RESLIST, -resl RESLIST
                        list of residue in receptor for PatchDock
  --binderdir BINDERDIR, -bd BINDERDIR
                        path to binders
  --dockdir DOCKDIR, -dd DOCKDIR
                        Path to prepared input files for docking
  --designdir DESIGNDIR, -ded DESIGNDIR
                        Path to fasta files of designed sequences
  --packerdir PACKERDIR, -pd PACKERDIR
                        Path to packed PDB files
  --minimdir MINIMDIR, -md MINIMDIR
                        Path to minimized files
  --dst DST             Directory to final results
  --np NP               num cores
  --ops OPS             mac or linux
  --csvname CSVNAME     Name of converted csv file
  --score SCORE         Scorefile to be converted to csv format
  --fn FN               A function name to use
```

<br>
<br>

If you want to run automatically, there is no need to specify all of the directories.

Mostly, directories were specified as default values.

```
python3 binder_design.py \
    --template ./template/template.pdb
    --binderdir ./binders/
    --reslist ./template/res_chain.list
    --patchdock /opt/tools/PatchDock/
    --np 4
    --ops linux
    --nstruct 1
```

The script automatically run the `PatchDock`

and `seq. design` > `rosetta packing` > `rosetta minimize` will be then repeated twice.

You can find 10 new folders in your working directory, and several csv files containing scores.

When running, You will see the message below.

```
*** 1st run ***
[PatchDock] Docking ...
[PatchDock] PDB generating ...
[ProteinMPNN] Seq. design running ...
[Rosetta] Packer running ...
[Rosetta] Minimize running ...
*** 2nd run ***
[ProteinMPNN] Seq. design running ...
[Rosetta] Packer running ...
[Rosetta] Minimize running ...
Total elapsed time:  113.91192865371704 sec
```

<br>

### etc

If you want to use only one function, not entire process,

try to use `--fn` flag.

```
# check the function list and require arguments

python3 binder_design.py --fn x

    Please select/enter the one of functions below and enter it.

    1) listgen (binderdir, minilist)
    2) patchdock (patchdock, template, reslist, minilist, dockdir, np)
    3) seq_design (dockdir, designdir)
    4) packer (ops, designdir, dockdir, packerdir, np)
    5) minimize (ops, packerdir, minimdir, np)
    6) txt2csv (score, csvname)
    7) pdbsorting (csvname, minimdir, designdir, dst)

```


```
# example

python3 binder_design.py \
    --fn txt2csv 
    --score ./final_score.txt 
    --csvname ./final_score_csv.csv 
```