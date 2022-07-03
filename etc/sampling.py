import os
import argparse
import random

parser = argparse.ArgumentParser (description='Script for designing protein binder')

parser.add_argument('--dir', required=False, type=str, default='./binders/', help='path to binders')
parser.add_argument('--ratio', required=False, type=int, default=0.1, help='sampling ratio')

args = parser.parse_args()



path = args.dir

file_list = os.listdir(path)
pdb_file_list = [file for file in file_list if file.endswith('.pdb')]

ls = ['HHH_bc', 'HHH_ems', 'HHHH_bc', 'HEEHE_lc', 'HHHH_ems', 'ferr_ems', 'ferr_hh', 'F4_ems']

def classify(x):
    x = [file for file in pdb_file_list if file[0:-10] == str(x)]
    return x



newpath = 'partial_' + str(path.lstrip('./'))

os.mkdir(newpath)

for i in ls:
    i = classify(i)
    newls = random.sample(i, round(len(i) * args.ratio))

    if len(i) * args.ratio < 1:
        f = ''.join(s for s in i)
        os.system(f'cp {path}/{f} {newpath}')
    else:
        for i in newls:
            os.system(f'cp {path}/{i} {newpath}')






