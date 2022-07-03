import os
import argparse


parser = argparse.ArgumentParser (description='Script for designing protein binder')

parser.add_argument('--residues', required=True, type=str, default='1-10', help='path to binders')
parser.add_argument('--chain', required=False, type=str, default='', help='path to binders')

args = parser.parse_args()



residues = args.residues

ls = residues.split('+')

for i in ls:
    i = i.split('-')
    # print(i)
    with open ('res.list', 'a') as file:
        if len(i) == 1:
            a = int(i[0])
            file.write(str(a) + '\n')

        else:
            a = int(i[0])
            b = int(i[1]) + 1
            for i in range (a, b):
                file.write(str(i) + '\n')     


if not args.chain == '':
    with open ('res_chain.list', 'a') as file:
        with open ('res.list', 'r') as file2:
            for i in file2:
                file.write(i.strip() + ' ' + args.chain + '\n')

print ('done')