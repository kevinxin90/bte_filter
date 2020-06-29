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

        ngd = N if fxy == 0 else ((max(np.log(fx),np.log(fy)) - np.log(fxy)) /\
                                (np.log(N) - min(np.log(fx),np.log(fy))))
        ngd_s = N if fxy_s == 0 else ((max(np.log(fx),np.log(fy)) - np.log(fxy_s)) /\
                                (np.log(N) - min(np.log(fx),np.log(fy))))
        return ngd, ngd_s

    co_occur = os.path.join(data_folder, 'summary_CoOccurs_2020.txt')
    freq = os.path.join(data_folder, 'MH_freq_counts_2020.txt')
    assert os.path.exists(co_occur)
    assert os.path.exists(freq)

    df_freq = pd.read_csv(freq, sep='|')
    chunks = pd.read_csv(co_occur, sep='|', chunksize=10**7)

    df_list = []
    for chunk in chunks:
        timeFrames = ['MBD', 'MED']
        for time in timeFrames:
            loc = 'Freq_' + time
            loc_s = 'StarFreq_' + time
            chunk[loc] = chunk[chunk['TimeFrame(SAB)'] == time]['Freq']
            chunk[loc_s] = chunk[chunk['TimeFrame(SAB)'] == time]['StarFreq(COF)']
        chunk.fillna(0, inplace=True)

        id_cols = ['DUI1', 'DUI2', 'CUI1', 'CUI2']
        cols = ['Freq', 'StarFreq(COF)', 'Freq_MBD', 'StarFreq_MBD', 'Freq_MED', 'StarFreq_MED']
        df_list.append(chunk.groupby(id_cols)[cols].sum().astype(int).reset_index())
    df = pd.concat(df_list).groupby(id_cols)[cols].sum().reset_index()

    for i,row in df.iterrows():
        doc = {'_id' : '-'.join([row[0], row[1]]),
              'MESH1' : row[0], 'MESH2' : row[1], 'UMLS1' : row[2], 'UMLS2' : row[3],
              'combo' : ['-'.join([i,j]) for i in [row[0],row[2]] for j in [row[1],row[3]]] + ['-'.join([i,j]) for i in [row[1],row[3]] for j in [row[0],row[2]]],
              'overall_freq' : {'total' : row[4], 'most_current_10' : row[6], 'most_current_5' : row[8]},
              'starred_freq' : {'total' : row[5], 'most_current_10' : row[7], 'most_current_5' : row[9]}}
        ngd = get_ngd(doc)
        doc['ngd_overall'] = ngd[0]
        doc['ngd_starred'] = ngd[1]
        yield doc
