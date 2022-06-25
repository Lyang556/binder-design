
import os
import random
import string
import argparse

import time


parser = argparse.ArgumentParser (description='Script for designing protein binder')

parser.add_argument('--targetseq', required=False, type=str, default='', help='Target sequence')
parser.add_argument('--fastas', required=False, type=str, default='./03_seq_design_files', help='Path to designed sequences using ProteinMPNN')
parser.add_argument('--dst', required=False, type=str, default='./for_colab', help='')

args = parser.parse_args()



def fasta_parsing (fastas=args.fastas, targetseq=args.targetseq, dst=args.dst):

    list_path = os.listdir(f'./{fastas}/')
    mpnn_list = [file for file in list_path if file.endswith('.fa')]

    targetseq = targetseq.strip()
    targetseq = targetseq.replace('\n', '')
    targetseq = targetseq.replace(' ', '')

    fastas = fastas.replace('\n', '')
    fastas = fastas.replace('./', '')

    dst = dst.replace('\n', '')
    dst = dst.replace('./', '')

    if dst in os.listdir('./'):
        pass

    else:
        os.mkdir(dst)
    
    fa = "".join([random.choice(string.ascii_letters) for _ in range(10)])

    for i in mpnn_list:
        j = i.replace('.fa', '')
        with open (f'./{dst}/{j}.fasta', 'w') as fa:
            with open (f'./{fastas}/{i}', 'r') as mpnn:
                seq = mpnn.readlines()[3].rstrip()
                fa.writelines(f'>{j}' + '\n')
                fa.writelines(f'{seq}:{targetseq}')



if __name__ == '__main__':

    start = time.time ()

    fasta_parsing (fastas=args.fastas, targetseq=args.targetseq, dst=args.dst)

    print ('Working done !')
    print ('Elapsed time: ', time.time() - start, 'sec')
    quit()        