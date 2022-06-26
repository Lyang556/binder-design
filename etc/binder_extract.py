
import os
import random
import string
import argparse

import time

from pymol import cmd


parser = argparse.ArgumentParser (description='Script for designing protein binder')

parser.add_argument('--pdbs', required=False, type=str, default='./03_seq_design_files', help='Path to designed sequences using ProteinMPNN')
parser.add_argument('--dst', required=False, type=str, default='./for_colab', help='')

args = parser.parse_args()



def binder_extract (pdbs=args.pdbs, dst=args.dst):

    list_path = os.listdir(f'./{pdbs}/')
    pdb_list = [file for file in list_path if file.endswith('.pdb')]

    dst = dst.replace('\n', '')
    dst = dst.replace('./', '')

    

    if dst in os.listdir('./'):
        pass
    else:
        os.mkdir(dst)
    
    
    for i in pdb_list:
        binder = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        binder = 'binder'
        cmd.load(f'./{pdbs}/{i}')
        cmd.create(binder, 'chain A')
        cmd.save(f'./{dst}/binder_{i}', f'{binder}')
        j = i.replace('.pdb','')
        cmd.remove(j)
        cmd.remove(binder)



    



if __name__ == '__main__':

    start = time.time ()

    binder_extract (pdbs=args.pdbs, dst=args.dst)

    print ('Working done !')
    print ('Elapsed time: ', time.time() - start, 'sec')
    quit()        