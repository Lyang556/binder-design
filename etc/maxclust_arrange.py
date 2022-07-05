#!/usr/bin/env python3

import os
import pandas as pd


os.system('/opt/tools/maxcluster64bit -l li.list > log.log')

with open ('log.log', 'r') as file:
    ls = file.readlines()


ls_num  = ls.index('INFO  : Item     Cluster\n')


sorted_ls = ls[ls_num + 1:-1]

df = pd.DataFrame(columns=['item', 'cluster', 'description'])

for i in sorted_ls:
    j = i.split(' ')
    info_clust = [x for x in j if x != '' and x != ':' and x != 'INFO']
    info_clust[2] = info_clust[2].rstrip()
    dict = {}
    dict['item'] = info_clust[0]
    dict['cluster'] = info_clust[1]
    dict['description'] = info_clust[2]

    df2 = pd.DataFrame([dict])

    df = pd.concat([df, df2], ignore_index=True)

df.to_csv('result.csv')