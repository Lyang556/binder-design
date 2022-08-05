
import os
import random
import string
import argparse
import glob

import time


parser = argparse.ArgumentParser (description='Script for designing protein binder')

parser.add_argument('--targetseq', required=False, type=str, default='', help='Target sequence')
parser.add_argument('--fastas', required=False, type=str, default='./03_seq_design_files', help='Path to designed sequences using ProteinMPNN')
parser.add_argument('--dst', required=False, type=str, default='./for_colab', help='')
parser.add_argument('--merge', action='store_true', help='merge fasta files')


args = parser.parse_args()



def fasta_parsing (fastas=args.fastas, targetseq=args.targetseq, dst=args.dst):

    fastas = os.path.abspath(fastas)
    dst = os.path.abspath(dst)

    list_path = os.listdir(fastas)
    mpnn_list = [file for file in list_path if file.endswith('.fa')]

    targetseq = targetseq.strip()
    targetseq = targetseq.replace('\n', '')
    targetseq = targetseq.replace(' ', '')

    if os.path.isdir(dst):
        pass
    else:
        os.mkdir (dst)

    
    fa = "".join([random.choice(string.ascii_letters) for _ in range(10)])

    for i in mpnn_list:
        j = i.replace('.fa', '')
        with open (f'{dst}/{j}.fasta', 'w') as fa:
            with open (f'{fastas}/{i}', 'r') as mpnn:
                seq = mpnn.readlines()[3].rstrip()
                fa.writelines(f'>{j}' + '\n')
                fa.writelines(f'{seq}:{targetseq}')

    if args.merge == True:
        fa_ls = glob.glob(f'{dst}/*.fasta')

        with open (f'{dst}/merged.fasta', 'a') as merged_fa:
            for i in fa_ls:
                with open (i, 'r') as fa_file:
                    fa_lines = fa_file.readlines()
                    merged_fa.writelines(fa_lines)
                    merged_fa.writelines('\n')


if __name__ == '__main__':

    start = time.time ()

    fasta_parsing (fastas=args.fastas, targetseq=args.targetseq, dst=args.dst)

    print ('Working done !')
    print ('Elapsed time: ', time.time() - start, 'sec')
    quit()        