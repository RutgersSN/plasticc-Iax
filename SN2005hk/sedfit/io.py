#input and output
import os
import glob
import numpy as np
import pandas as pd
from astropy.table import Table,vstack,hstack,join,Column
from astropy.io import ascii,fits

def read_data(datafile=None,snlist=None,sample=None,datatype='unknown'):
    if datafile == None:
        # print datatype
        if datatype == 'lc':
            datafile = 'data/allsne.txt'
        elif datatype == 'host':
            datafile = 'data/host/allsne_host.txt'
        elif datatype == 'meta':
            datafile = 'data/meta/allsne_sninfo.txt'
        elif datatype == 'snstat':
            datafile = 'data/snstats.txt'
        else:
            print "Reading in unknown data type from: ",datafile
    else:
        datafile = datafile
    print "Reading in ", datatype ," data from", datafile
    data = ascii.read(datafile)
    if snlist != None:
        snnames = pd.Series(data['Name'])
        data = data[np.array(snnames.isin(snlist))]
    if sample != None:
        if datatype == 'lc':
            sample_names = pd.Series(data['Survey'])
        else:
            sample_names = pd.Series(data['Sample'])
        data = data[np.array(sample_names.isin(sample))]
    return data

