import os
import argparse
import shutil

import pandas as pd

parser = argparse.ArgumentParser (description='test')

parser.add_argument('--csv', required=False, type=str, default='./', help='')
parser.add_argument('--pdbdir', required=False, type=str, default='./', help='')
parser.add_argument('--sort', required=False, type=str, default='hbonds_int', help='')
parser.add_argument('-a', '--ascending', dest='ascending', action='store_true')


args = parser.parse_args()


if args.ascending == True:
    df = pd.read_csv(args.csv).sort_values(by=args.sort, ascending=True)
else:
    df = pd.read_csv(args.csv).sort_values(by=args.sort, ascending=False)

df.reset_index().to_csv(args.csv.replace('.csv', '') + '_sorted.csv')

pdblist = df['description'].values.tolist()[0:200]

sortlist = df[args.sort].values.tolist()[0:200]

dict_ = {key:value for key, value in zip(pdblist, sortlist)}

dir = os.path.abspath(args.pdbdir)
dst = os.path.join(os.getcwd(), 'pdb_sorting')

try:
    os.mkdir(dst)
except FileExistsError:
    print ('Please remove the directory: ', dst)
    exit()

num = 0
for key, value in dict_.items():
    key = key + '.pdb'
    value = str(value)
    num += 1
    shutil.copyfile(f'{dir}/{key}', f'{dst}/{num}_{value}_{key}')


print ('done')

