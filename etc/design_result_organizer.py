import os
import pandas as pd


mpnn_list = os.listdir('./')
fasta_list = [file for file in mpnn_list if file.endswith('.fa') ]


df1 = pd.DataFrame(columns=['description', 'score', 'recovery'])

for i in fasta_list:
    with open (i, 'r') as file:

        fa_header = file.readlines()[2]

        recov = fa_header.rstrip().split(',')[-1].split('=')[-1]
        score = fa_header.rstrip().split(',')[-2].split('=')[-1]

        df2 = pd.DataFrame([[i, score, recov]], columns=['description', 'score', 'recovery'])

        df1 = pd.concat([df1, df2], ignore_index=True)

        # df1 = df1.append(pd.DataFrame([[i, score, recov]], columns=['name', 'score', 'recovery']), ignore_index=True)

df1.to_csv('design_result_merged.csv')

