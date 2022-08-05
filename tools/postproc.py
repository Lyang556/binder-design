#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
import numpy as np
import shutil

class PostProc:

    def __init__ (self, args):

        self.args = args

    def txt2csv (self, score):

        csvname = self.args.csvname
        
        pd.read_csv(score, delimiter=r'\s+').to_csv(csvname)

        with open (csvname, 'r') as file:
            ls = file.readlines()
            with open (csvname, 'w') as newfile:
                for i, j in enumerate (ls):
                    if i == 0:
                        pass
                        # remove first line that is an unnecessary thing
                    else:
                        newfile.writelines(j)

        df = pd.read_csv(csvname)

        df = df[['description','dSASA_int','dSASA_hphobic','dSASA_polar','hbonds_int','total_score','sc_value','complex_normalized','dG_cross','dG_cross/dSASAx100','dG_separated','dG_separated/dSASAx100','delta_unsatHbonds','dslf_fa13','fa_atr','fa_dun_dev','fa_dun_rot','fa_dun_semi','fa_elec','fa_intra_atr_xover4','fa_intra_elec','fa_intra_rep_xover4','fa_intra_sol_xover4','fa_rep','fa_sol','hbond_E_fraction','hbond_bb_sc','hbond_lr_bb','hbond_sc','hbond_sr_bb','hxl_tors','lk_ball','lk_ball_bridge','lk_ball_bridge_uncpl','lk_ball_iso','nres_all','nres_int','omega','p_aa_pp','packstat','per_residue_energy_int','pro_close','rama_prepro','ref','side1_normalized','side1_score','side2_normalized','side2_score']]
        df[df.dSASA_int > 0].reset_index(drop=True).to_csv(csvname)


    def pdbsorting (self):

        csvname = os.path.abspath(self.args.csvname)
        designdir = os.path.abspath(self.args.designdir)
        # designdir_seqdir = os.path.join(designdir, 'seqs')
     
        df_orig = pd.read_csv(csvname)
        condition = (df_orig.sc_value > 0.4) & (df_orig.dSASA_int > 900) & (df_orig.hbonds_int > 3)

        df = df_orig[condition].sort_values(by=['dG_separated'], ascending=True).reset_index(drop=True)
        col_list = list(df.columns)
        col_list[0] = 'original_number'
        df.columns = col_list

        # set the list of pdb files for copying.
        pdblist = df['description'].values.tolist()[0:200]

        csvname = csvname.replace('.csv', '')
        dst = os.path.abspath(self.args.dst)
        minimdir = os.path.abspath(self.args.minimdir)
        minimdir_name = os.path.basename(minimdir)

        pdbpath = 'rank_' + minimdir_name
        pdbpath = os.path.join(dst, pdbpath)

        try:
            os.mkdir(dst)
        except FileExistsError:
            pass

        try:
            os.mkdir(pdbpath)
        except FileExistsError:
            pass

        if 'fastas' in os.listdir(pdbpath) and os.path.isdir(os.path.join(pdbpath, 'fastas')):
            pass
        else:
            os.mkdir(f'{pdbpath}/fastas')

        df.loc[0:200, :].to_csv(f'{dst}/final_score_top200.csv')

        # copy pdb files
        fa_ls = []
        for i, j in enumerate (pdblist):
            i += 1
            m = j[:-26]
            fa_ls.append(m)
            shutil.copyfile(f'{minimdir}/{j}.pdb', f'{pdbpath}/rank_{i}_in_{minimdir_name}_{j}.pdb')

        # copy fasta files
        i = 0
        for m in fa_ls:
            i += 1
            try:
                shutil.copyfile(f'{designdir}/seqs/{m}.fa', f'{pdbpath}/fastas/rank_{i}_in_{minimdir_name}_{m}.fa')
            except FileNotFoundError:
                pass

            