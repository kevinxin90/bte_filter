"""
Parser file for API
"""

import os
import pandas as pd
import numpy as np

def load_json(data_folder):

    # returns ngd and ngd*
    def get_ngd(doc):
        N = 2.7e+7 * 20 # from PubMed home page there are 27 million articles; avg 20 MeSH terms per article
        fx = df_freq[df_freq['DUI'] == doc['MESH1']]['Total'].item()
        fy = df_freq[df_freq['DUI'] == doc['MESH2']]['Total'].item()
        fxy = doc['overall_freq']['total']
        fxy_s = doc['starred_freq']['total']

        ngd = (max(np.log(fx),np.log(fy)) - np.log(fxy)) / (np.log(N) - min(np.log(fx),np.log(fy)))
        ngd_s = (max(np.log(fx),np.log(fy)) - np.log(fxy_s)) / (np.log(N) - min(np.log(fx),np.log(fy)))
        return ngd, ngd_s

    co_occur = os.path.join(data_folder, 'summary_CoOccurs_2020.txt')
    freq = os.path.join(data_folder, 'MH_freq_counts_2020.txt')
    assert os.path.exists(co_occur)
    assert os.path.exists(freq)

    df = pd.read_csv(co_occur, sep='|')
    df_freq = pd.read_csv(freq, sep='|')

    timeFrames = ['MBD', 'MED']
    for time in timeFrames:
        loc = 'Freq_' + time
        loc_s = 'StarFreq_' + time
        df[loc] = df[df['TimeFrame(SAB)'] == time]['Freq']
        df[loc_s] = df[df['TimeFrame(SAB)'] == time]['StarFreq(COF)']
    df.fillna(0, inplace=True)

    cols = ['Freq', 'StarFreq(COF)', 'Freq_MBD', 'StarFreq_MBD', 'Freq_MED', 'StarFreq_MED']
    df = df.groupby(['DUI1', 'DUI2', 'CUI1', 'CUI2'])[cols].sum().astype(int).reset_index()

    for i,row in df.iterrows():
        doc = {'_id' : '-'.join([row[0], row[1]]),
              'MESH1' : row[0], 'MESH2' : row[1], 'UMLS1' : row[2], 'UMLS2' : row[3],
              'combo' : ['-'.join([i,j]) for i in [row[0],row[2]] for j in [row[1],row[3]]],
              'overall_freq' : {'total' : row[4], 'most_current_10' : row[6], 'most_current_5' : row[8]},
              'starred_freq' : {'total' : row[5], 'most_current_10' : row[7], 'most_current_5' : row[9]}
              }
        ngd = get_ngd(doc)
        doc['ngd_overall'] = ngd[0]
        doc['ngd_starred'] = ngd[1]
        yield doc
