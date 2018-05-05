##Read in data from different sources and write to uniform format
import os
import glob
import numpy as np
import pandas as pd
from astropy.table import Table,vstack,hstack,join,Column
from astropy.io import ascii,fits
from astropy import units as u
from astropy.coordinates import SkyCoord,ICRS

def read_lc(sample,rootfolder=None):
    data = Table()
    data1 = Table()
    data2 = Table()

    # if sample == 'CSP_DR2':
    #     folder = 'CSP/CSP_Photometry_DR2'
    #     filefolder = os.path.join(rootfolder,folder)
    #     files = glob.glob(filefolder+"/*.dat")
    #     for f in files:
    #         # print f
    #         data0 = ascii.read(f)
    #         for i,col in enumerate(data0.meta['comments'][-1].split()):
    #             # print i, col
    #             if i == 0:
    #                 data0.rename_column(data0.colnames[i],col)
    #             elif i == 1 or i%2 != 0:
    #                 data0.rename_column(data0.colnames[i],col)
    #                 data0.rename_column(data0.colnames[i+1],col+'_err')
    #         data0['Name'] = f[len(filefolder)+1:f.find('opt')]
    #         data1 = vstack([data1, data0], join_type='outer')
    #     for name in set(data1['Name'].tolist()):
    #         data3 = Table()
    #         for filt in ['u','B','V','g','r','i','Y','J','H','K']:
    #             idx = data1['Name'] == name
    #             data2 = data1[idx]['MJD', filt, filt+'_err']
    #             data2.rename_column(filt,'Mag')
    #             data2.rename_column(filt+'_err','MagErr')
    #             data2['Filter'] = filt
    #             data2 = data2[np.all([~data2['Mag'].mask,data2['Mag']<99.],axis=0)].copy()
    #             data3 = vstack([data3,data2],join_type='outer')
    #         data3['Name'] = name
    #         data = vstack([data,data3],join_type='outer')

    if sample == 'CSP':
        folder = 'CSP/'
        filefolder = os.path.join(rootfolder,folder)
        f11 = filefolder+"/CSP3_opt_ajaa8df0t7_mrt.txt"
        f12 = filefolder+"/CSP3_opt2_ajaa8df0t8_ascii.txt"
        files2 = glob.glob(filefolder+"/CSP3_nir*.txt")

        t1 = Table()
        t11 = ascii.read(f11)
        newcolnames = ['Name','JD','u','u_err','g','g_err','r','r_err','i','i_err',
                       'B','B_err','V','V_err']
        for i,col in enumerate(t11.colnames):
            t11.rename_column(col,newcolnames[i])
        t11['Instrument'] = 'Swope'
        t1 = vstack([t1,t11])
        # print data1

        t12 = Table(names=newcolnames,
                  dtype=['S20']+['f8']*(len(newcolnames)-1))
        with open(f12) as f:
            for line in f:           
                if line[0] == '#':
                    continue
                if line[:2] == 'SN':
                    snname = line.strip().replace(' ','')
                    # print snname
                else:
                    row = []
                    for i,col in enumerate(line.split()):
                        if i == 0:
                            col = col.replace(',','')
                            row.append(float(col))
                        else:
                            c = col.split('(')
                            row.append(float(c[0]))
                            x = '0.'+c[1].split(')')[0]
                            row.append(float(x))
                    t12.add_row([snname]+row)
        t12['Instrument'] = 'duPont'
        t1 = vstack([t1,t12])
        # print t

        t2 = Table()
        for f in files2:
            t21 = ascii.read(f)
            newcolnames = ['Name','JD','Y','Y_err','J','J_err','H','H_err']
            for i,col in enumerate(t21.colnames):
                t21.rename_column(col,newcolnames[i])
            if 't10' in f:
                t21['Instrument'] = 'Swope1'
            if 't11' in f:
                t21['Instrument'] = 'Swope2'
            if 't12' in f:
                t21['Instrument'] = 'duPont'
            t2 = vstack([t2,t21])   
        data1 = join(t1,t2,join_type='outer')

        # print data1

        for name in set(data1['Name'].tolist()):
            data3 = Table()
            for filt in ['u','B','V','g','r','i','Y','J','H','K']:
                idx = data1['Name'] == name
                if not filt in data1[idx].colnames:
                    continue
                data2 = data1[idx]['JD', filt, filt+'_err','Instrument']
                data2.rename_column(filt,'Mag')
                data2.rename_column(filt+'_err','MagErr')
                data2['Filter'] = filt
                data2 = data2[np.all([~data2['Mag'].mask,data2['Mag']<99.,data2['Mag']!=0.],axis=0)].copy()
                data3 = vstack([data3,data2],join_type='outer')
            data3['Name'] = name
            data = vstack([data,data3],join_type='outer')
            data['MJD'] = data['JD'] - 2400000.5

    if sample == 'Swift':
        folder = 'Swift/individual_downloads'
        f_ia = os.path.join(rootfolder,folder,'Ialist.txt')
        ia = Table.read(f_ia,format = 'ascii')
        filefolder = os.path.join(rootfolder,folder)
        files = glob.glob(filefolder+"/*.dat")
        for f in files:
            # print f
            data0 = ascii.read(f)
            for ch in data0.meta['comments']:
                if ch[:4] == 'Filt':
                    colnames = ch.split()            
                    for i,col in enumerate(colnames):
                        # print i, col
                        data0.rename_column(data0.colnames[i],col.split('[')[0])
            df = data0.to_pandas()
            df = df.replace(['NULL','**********'],np.nan)
            df[df.columns[1:]] = df[df.columns[1:]].astype(float)
            data0 = Table.from_pandas(df)    
            data0['Name'] = f[len(filefolder)+1:f.find('_uvot')] 
            if data0['Name'][0] in ia['Name']:
                data = vstack([data, data0], join_type='outer')
            # print data0

    if sample in ['CfA','CfA1','CfA2']:
        if sample == 'CfA':
            folder = 'CfA/cfalc_allsnIa'
        elif sample == 'CfA1':
            folder = 'CfA/cfalc_R99snIa'
        elif sample == 'CfA2':
            folder = 'CfA/cfalc_J06snIa'
        filefolder = os.path.join(rootfolder,folder)
        files = glob.glob(filefolder+"/*.dat")
        for f in files:
            # print f
            data0 = ascii.read(f,format='no_header',delimiter='\s')
            for i,col in enumerate(data0.meta['comments'][-1].split()):
                if i == 0:
                    data0.rename_column(data0.colnames[i],col)
                elif i == 1 or i%2 != 0:
                    data0.rename_column(data0.colnames[i],col)
                    data0.rename_column(data0.colnames[i+1],col+'_err')
            for i,row in enumerate(data0):
                for j,d in enumerate(row):
                    if isinstance(d,str) and d[0]=='(':
                        data0[i][data0.colnames[j]] = float(d.translate(None,"()"))
            df = data0.to_pandas()
            df = df.astype(float)
            data0 = Table.from_pandas(df)    
            data0['Name'] = f[len(filefolder)+1:f.find('_',len(filefolder))] 
            # print data0
            # print data1     
            data1 = vstack([data1, data0], join_type='outer')

        for name in set(data1['Name'].tolist()):
            # print name
            data3 = Table()
            for filt in ['U','B','V','R','I','J','H','K']:
                idx = data1['Name'] == name
                # print data1.colnames
                list1 = data1[idx].colnames
                if not filt in list1:
                    # print filt
                    continue
                list2 = pd.Series(['HJD','JD-2400000','JD-2440000','JD-2450000','Julian_Day'])
                cols = list(list2[np.array(list2.isin(list1))])
                # print cols
                data2 = data1[idx][cols+[filt]+[filt+'_err']]
                data2.rename_column(filt,'Mag')
                data2.rename_column(filt+'_err','MagErr')
                data2['Filter'] = filt
                data2 = data2[np.all([~data2['Mag'].mask,data2['Mag']<99.],axis=0)].copy()
                data3 = vstack([data3,data2],join_type='outer')
            data3['Name'] = name
            # print data3
            data = vstack([data,data3],join_type='outer')

        data['MJD'] = 0.0
        if 'HJD' in data.colnames:
            data['MJD'][~data['HJD'].mask] = data['HJD'][~data['HJD'].mask] - 2400000.5
        if 'JD-2400000' in data.colnames:
            data['MJD'][~data['JD-2400000'].mask] = data['JD-2400000'][~data['JD-2400000'].mask] - 0.5 #+ 2400000. - 2400000.5
        if 'JD-2440000' in data.colnames:
            data['MJD'][~data['JD-2440000'].mask] = data['JD-2440000'][~data['JD-2440000'].mask] + 39999.5 #+ 2440000. - 2400000.5
        if 'JD-2450000' in data.colnames:
            data['MJD'][~data['JD-2450000'].mask] = data['JD-2450000'][~data['JD-2450000'].mask] + 49999.5 #+ 2450000. - 2400000.5
        if 'Julian_Day' in data.colnames:
            data['MJD'][~data['Julian_Day'].mask] = data['Julian_Day'][~data['Julian_Day'].mask] - 2400000.5

    if sample == 'CfA3':
        folder = 'CfA/'
        datafile = 'cfa3lightcurves.naturalsystem.txt'
        file = os.path.join(rootfolder,folder,datafile)
        data = Table(names=['Name','Fltnum','MJD','Mag','MagErr'],
                     dtype=['S10','i4','f8','f8','f8'])
        with open(file) as f:
            for line in f:           
                if line[0] == '#':
                    continue
                if line[:2] == 'sn':
                    snname = line.rstrip().replace('sn0','sn200')
                    # print snname
                else:
                    # print line.split()
                    data.add_row([snname]+line.split())
                # print data
        fltnum = [1,2,3,4,5,13,14]
        flt = ['U','B','V','R','I',"r'","i'"]
        data['Filter'] = 'xx'
        for i, filt in enumerate(flt):
            data['Filter'][data['Fltnum']==fltnum[i]] = filt
        data['Instrument'] = 'Keplercam'
        for name in set(data['Name']):
            if np.max(data['Fltnum'][data['Name'] == name]) < 10.:
                data['Instrument'][data['Name'] == name] = '4Shooter' 

    if sample == 'CfA4':
        folder = 'CfA/'
        datafile = 'cfa4.lc.natsystem.sort.ascii'
        f = os.path.join(rootfolder,folder,datafile)
        data = ascii.read(f,comment='###')
        for i,col in enumerate(data.meta['comments'][-1].split()):
            # print i, col
            data.rename_column(data.colnames[i],col)
        data.rename_column('SN','Name')
        data.rename_column('sigma_{total}','MagErr')

    if sample == 'CfAIR2':
        folder = 'CfAIR2/'
        datafile = 'cfair2.txt'
        f = os.path.join(rootfolder,folder,datafile)
        data = ascii.read(f)
        data.rename_column('SN','Name')
        data.rename_column('mag','Mag')
        data.rename_column('e_mag','MagErr')      
        data.rename_column('Band','Filter')      

    if sample == 'Pan-STARRS':
        folder = 'Pan-STARRS'
        datafile = 'rest2014.txt'
        f = os.path.join(rootfolder,folder,datafile)
        data = ascii.read(f)
        data.rename_column('SN','Name')  
        data.rename_column('mag','Mag')
        data.rename_column('e_mag','MagErr')    
        data.rename_column('Filt','Filter')

    if sample == 'SNLS':
        folder = 'SNLS/'
        datafile = 'guy2010.txt'
        f = os.path.join(rootfolder,folder,datafile)
        data = ascii.read(f)
        newcolnames = ['Name','Filter','MJD','Flux','FluxErr']
        for i,col in enumerate(data.colnames):
            data.rename_column(col,newcolnames[i])
        data['Mag'] = -2.5*np.log10(data['Flux']) + 30.
        data['MagErr'] = np.absolute(1.086/data['Flux']*data['FluxErr'])
        data = data[data['Mag']>0.]

    if sample == 'SDSS':
        # SNTYPE=118,120 are confirmed SNe Ia (Sako et al. 2014)
        folder = 'SDSS/SDSS_dataRelease-snana/SDSS_allCandidates'
        data_header = 'SDSS_allCandidates+BOSS_HEAD.FITS'
        data_phot = 'SDSS_allCandidates+BOSS_PHOT.FITS'
        header = os.path.join(rootfolder,folder,data_header)
        phot = os.path.join(rootfolder,folder,data_phot)
        f1 = fits.open(header)
        f2 = fits.open(phot)
        metadata = Table(f1[1].data)
        lcdata = Table(f2[1].data)
        for meta in metadata:
            if meta['SNTYPE'] == 118 or meta['SNTYPE'] == 120:
                start = meta['PTROBS_MIN']-1
                end = meta['PTROBS_MAX']
                data0 = lcdata[start:end]
                data0['Name'] = 'SDSS' + meta['SNID']
                data = vstack([data, data0], join_type='outer')
        data.rename_column('MAG','Mag')
        data.rename_column('MAGERR','MagErr')    
        data.rename_column('FLT','Filter')

    if sample == 'Hamuy96':
        folder = 'Hamuy96'
        datafile = 'hamuy1996.txt'
        f = os.path.join(rootfolder,folder,datafile)
        data1 = ascii.read(f)
        data1.rename_column('SN','Name')
        for name in set(data1['Name'].tolist()):
            # print name
            data3 = Table()
            for filt in ['B','V','R','I']:
                idx = data1['Name'] == name
                # print data1.colnames
                data2 = data1[idx]['HJD', filt+'mag', 'e_'+filt+'mag']
                data2.rename_column(filt+'mag','Mag')
                data2.rename_column('e_'+filt+'mag','MagErr')
                data2['Filter'] = filt
                data2 = data2[~data2['Mag'].mask].copy()
                data3 = vstack([data3,data2],join_type='outer')
            data3['Name'] = name
            data = vstack([data,data3],join_type='outer')
        data['MJD'] = data['HJD'] - 2400000.5

    if sample == 'Foundation':
        folder = 'Foundation/Foundation_DS17'
        filefolder = os.path.join(rootfolder,folder)
        files = glob.glob(filefolder+"/found*.txt")
        for f in files:   
            # print f
            data0 = ascii.read(f,data_start=14,data_end=-1,header_start=13)
            # print data0
            data0['Name'] = f[f.find('found')+5:f.find('.txt')] 
            data = vstack([data, data0], join_type='outer')
        data.rename_column('MAG','Mag')
        data.rename_column('MAGERR','MagErr')
        data.rename_column('FLT','Filter')
        data.remove_column('VARLIST:')  

    return data

def read_sn_info(sample,rootfolder=None,metafolder='metadata'):
    if rootfolder == None:
        rootfolder = '/home/mi/Desktop/project-data'

    filefolder = os.path.join(rootfolder,metafolder)

    meta = Table()
    t = Table()

    # if sample == 'CSP_DR2':
    #     files = glob.glob(filefolder+"/CSP*.txt")
    #     for f in files:
    #         print f
    #         t0 = ascii.read(f,format='no_header',delimiter='\t')
    #         newcolnames = ['Name','ra(2000)','dec(2000)','Host','Morphology','z_helio','Discovery','Discovery Group']
    #         for i,col in enumerate(t0.colnames):
    #             t0.rename_column(col,newcolnames[i])
    #         # print t0
    #         meta = vstack([meta,t0],join_type='outer')
    #     c = SkyCoord(meta['ra(2000)'],meta['dec(2000)'],unit=(u.hourangle,u.degree))
    #     meta['RA'] = c.ra.degree
    #     meta['DEC'] = c.dec.degree

    if sample == 'CSP':
        # files = glob.glob(filefolder+"/CSP3*.txt")
        f1 = filefolder+"/CSP3_ajaa8df0t1_ascii.txt"
        f2 = filefolder+"/CSP3_ajaa8df0t2_ascii.txt"
        for f in [f1,f2]:
            print f
            t0 = ascii.read(f,format='no_header',delimiter='\t')
            if f == f1:
                newcolnames = ['Name','ra(2000)','dec(2000)','Host','Morphology','z_helio(str)','Discovery','Discovery Group']
            elif f == f2:
                newcolnames = ['Name','delta_m15(str)','Type','Type2-Wang','Type3-Branch','N_opt','N_ir','T_sp']
            for i,col in enumerate(t0.colnames):
                t0.rename_column(col,newcolnames[i])
            # print t0
            if len(meta) == 0:
                meta = vstack([meta,t0],join_type='outer')
            else:
                meta = join(meta,t0,join_type='left')
        c = SkyCoord(meta['ra(2000)'],meta['dec(2000)'],unit=(u.hourangle,u.degree))
        meta['RA'] = c.ra.degree
        meta['DEC'] = c.dec.degree
        meta['z_helio'] = [float(x.split("^")[0]) for x in meta['z_helio(str)']]

    if sample == 'Swift':
        f= filefolder+"/Swift_onlinelist.txt"
        meta = ascii.read(f)
        meta['z_helio'] = meta['redshift']  ##need to verify this...

    if sample in ['CfA','CfA1','CfA2']:

        files1 = glob.glob(filefolder+"/CfA_Riess*.txt")
        files2 = glob.glob(filefolder+"/CfA_Jha*.txt")
        if sample in ['CfA','CfA1']:
            for f in files1:
                print f
                t0 = ascii.read(f)
                # print t0
                if len(meta) == 0:
                    meta = vstack([meta,t0])
                else:
                    meta = join(meta,t0,join_type='outer')
            meta.rename_column('Type','Morphology')
            meta.rename_column('Galaxy','Host')
            meta['ra(B1950)'] = [x.replace('~','') for x in meta['ra']]
            meta['dec(B1950)'] = [x.replace('~','') for x in meta['dec']]
            c1 = SkyCoord(meta['ra(B1950)'],meta['dec(B1950)'],
                          unit=(u.hourangle,u.degree),equinox='B1950',frame='fk4')
            c2 = SkyCoord(c1.ra,c1.dec)
            c = c1.transform_to(c2)
            meta['RA'] = c.ra.degree
            meta['DEC'] = c.dec.degree
            # print meta['Name','ra(B1950)','RA','dec(B1950)','DEC']

        if sample in ['CfA','CfA2']:
            for f in files2:
                print f
                t0 = ascii.read(f)
                meta = vstack([meta,t0],join_type='outer')

        meta['z_helio'] = 0.
        if 'cz_helio' in meta.colnames:
            meta['z_helio'][meta['cz_helio'].mask == False] = meta[meta['cz_helio'].mask == False]['cz_helio']/300000.
        if 'log_cz' in meta.colnames:
            meta['z_helio'][meta['log_cz'].mask == False] = np.power(10.,meta[meta['log_cz'].mask == False]['log_cz'])/300000.


    if sample == 'CfA3':
        files = glob.glob(filefolder+"/CfA3*.txt")
        for f in files:
            print f
            t0 = ascii.read(f,delimiter='\t')
            # print t0
            if len(meta) == 0:
                meta = vstack([meta,t0])
            else:
                meta = join(meta,t0,join_type='inner')  
        meta.rename_column('Galaxy','Host')
        c = SkyCoord(meta['Position'],unit=(u.hourangle,u.degree))
        meta['RA'] = c.ra.degree
        meta['DEC'] = c.dec.degree

    if sample == 'CfA4':
        files = glob.glob(filefolder+"/CfA4*.txt")
        for f in files:
            print f
            t0 = ascii.read(f,delimiter='\t')
            # print t0
            if len(meta) == 0:
                meta = vstack([meta,t0])
            else:
                meta = join(meta,t0,join_type='outer')  
        meta.rename_column('Galaxy','Host')
        meta.rename_column('E(B - V)','ebv')    
        c = SkyCoord(meta['Position'],unit=(u.hourangle,u.degree))
        meta['RA'] = c.ra.degree
        meta['DEC'] = c.dec.degree  

    if sample == 'CfAIR2':
        files = glob.glob(filefolder+"/CfAIR2*.txt")
        for f in files:
            print f
            t0 = ascii.read(f,comment='#',delimiter='\s')
            # print t0
            if len(meta) == 0:
                meta = vstack([meta,t0])
            else:
                meta = join(meta,t0,join_type='outer')  
        meta.rename_column('ra','ra(deg)')
        meta.rename_column('dec','dec(deg)')
        meta['RA'] = meta['ra(deg)']
        meta['DEC'] = meta['dec(deg)']


    if sample == 'SNLS':
        f = filefolder+'/SNLS_Sullivan2010_t1.txt'
        meta = ascii.read(f)
        meta['logM_val'] = [x.split('+')[0] for x in meta['logM']] 
        meta['logM_p'] = [x.split('+')[1].split('-')[0] for x in meta['logM']] 
        meta['logM_m'] = [x.split('+')[1].split('-')[1] for x in meta['logM']] 
        idx = meta['logSFR'] != '<-3.00'
        meta['logSFR_val'] = Column(dtype='S10',length=len(meta))
        meta['logSFR_p'] = Column(dtype='S10',length=len(meta))
        meta['logSFR_m'] = Column(dtype='S10',length=len(meta))
        meta['logSFR_val'] = 'N/A'
        meta['logSFR_p'] = 'N/A'
        meta['logSFR_m'] = 'N/A'
        meta['logSFR_val'][idx] = [x.split('+')[0] for x in meta[idx]['logSFR']] 
        meta['logSFR_p'][idx] = [x.split('+')[1].split('-')[0] for x in meta[idx]['logSFR']]
        meta['logSFR_m'][idx] = [x.split('+')[1].split('-')[1] for x in meta[idx]['logSFR']] 
        meta.rename_column('z_cmb','z_CMB')

        f2 = filefolder+'/SNLS_Guy_t9.txt'
        meta2 = ascii.read(f2)

        meta2.rename_column('z','z_helio')
        meta2.rename_column('e_z','dz_helio')

        meta = join(meta,meta2,join_type='outer')
        c = SkyCoord(meta['RAJ2000'],meta['DEJ2000'],unit=(u.hourangle,u.degree))
        meta['RA'] = c.ra.degree
        meta['DEC'] = c.dec.degree  

    if sample == 'SDSS':
        # SNTYPE=118,120 are confirmed SNe Ia (Sako et al. 2014)
        folder = 'SDSS/SDSS_dataRelease-snana/SDSS_allCandidates'
        data_header = 'SDSS_allCandidates+BOSS_HEAD.FITS'
        header = os.path.join(rootfolder,folder,data_header)
        f1 = fits.open(header)
        metadata = Table(f1[1].data)
        # print len(metadata)
        # print np.sum(metadata['SNTYPE']==118)
        # print np.sum(metadata['SNTYPE']==120)
        metadata = metadata[np.any([metadata['SNTYPE']==118,metadata['SNTYPE']==120],axis=0)]
        metadata['Name'] = ['SDSS' + x.strip() for x in metadata['SNID']]

        f = filefolder+'/SDSS_Sako2014_t1_free_format.txt'
        meta = ascii.read(f)
        meta['Name'] = ['SDSS' + str(x) for x in meta['CID']]
        meta['Host'] =  meta['objIDHost']
        # print len(meta['Name'][0]),len(metadata['Name'][0])
        # print len(meta),len(metadata)
        meta = join(meta,metadata['Name','SNTYPE'],join_type='inner')
        meta['z_helio'] = [float(x.strip()) for x in meta['zspecHelio']]
        meta['z_CMB'] = [float(x.strip()) for x in meta['zCMB']]
        # print meta

    if sample == 'Hamuy96':
        f = filefolder+'/Hamuy96_t1.csv'
        meta = ascii.read(f,format='csv')
        meta.rename_column('SN','Name')
        c = SkyCoord(meta['Ra(2000)'],meta['Dec(2000)'],unit=(u.hourangle,u.degree))
        meta['RA'] = c.ra.degree
        meta['DEC'] = c.dec.degree  

    if sample == 'Foundation':
        f = filefolder+'/Foundation_spreadsheet.csv'
        meta = ascii.read(f,format='csv')
        meta.rename_column('SN Name','Name')
        c = SkyCoord(meta['SN RA'],meta['Dec (2000)'],unit=(u.hourangle,u.degree))
        meta['RA'] = c.ra.degree
        meta['DEC'] = c.dec.degree  
        meta['z_helio'] = meta['Heliocentric Redshift']
        meta['dz_hel'] = meta['Uncertainty']

    if sample == 'Pan-STARRS':
        files = glob.glob(filefolder+"/PanSTARRS*.txt")
        for f in files:
            print f
            t0 = ascii.read(f,comment='#',delimiter='\t')
            # print t0
            if len(meta) == 0:
                meta = vstack([meta,t0])
            else:
                meta = join(meta,t0,join_type='outer') 
        meta.rename_column('c','c(err)')
        meta['z_hel'].fill_value = '-99.(-99.)' 
        meta['z_cmb'].fill_value = '-99.(-99.)' 
        meta['E(B - V)'].fill_value = '-99.(-99.)' 
        # print meta['z_hel','z_cmb']
        meta['z_helio'] = [np.float(x.split('(')[0]) for x in meta.filled()['z_hel']]
        meta['z_CMB'] = [np.float(x.split('(')[0]) for x in meta.filled()['z_cmb']]
        meta['ebv'] = [np.float(x.split('(')[0]) for x in meta.filled()['E(B - V)']]
        meta['debv'] = [np.float(x.split('(')[1].split(')')[0]) for x in meta.filled()['E(B - V)']]
        c = SkyCoord(meta['R.A. (J2000)'],meta['Decl. (J2000)'],unit=(u.hourangle,u.degree))
        meta['RA'] = c.ra.degree
        meta['DEC'] = c.dec.degree          

    if sample == 'Uddin':
        f1 = '/home/mi/Desktop/project-data/metadata/Combination_Uddin_2017_t8.fits'
        hdu1 = fits.open(f1)
        # print hdu[1].header
        t1 = Table(hdu1[1].data)
        print t1.colnames
        f2 = '/home/mi/Desktop/project-data/metadata/Combination_Uddin_2017_t7.fits'
        hdu2 = fits.open(f2)
        t2 = Table(hdu2[1].data)
        print t2.colnames
        meta = join(t1,t2,join_type='inner')
        # print meta
        meta.rename_column('name','Name')
        for m in meta:
            if m['Name'][0:4] == 'SNLS':
                m['Name'] = m['Name'][4:]
        # print len(meta)
        meta = meta[meta['type']=='s']
        # print len(meta)


    if sample == 'Neill':
        f1 = filefolder+"/Combine_lowz_Neill2009t1.txt"
        f2 = filefolder+"/Combine_lowz_Neill2009t2.txt"
        t1 = ascii.read(f1)
        t2 = ascii.read(f2)
        meta = join(t1[t1.colnames[0:-2]],t2,join_type='inner')  
        meta.rename_column('SN','Name')


    if sample == 'Wolf':
        f = filefolder+'/SDSS_Wolf2016_t4.txt'
        meta = ascii.read(f)
        cols = ['CID','IAUName','Type','z_helio','Color','e_Color',
                'Stretch','e_Stretch','HR','e_HR','DR8-OBJID','BPT',
                'LogMass','e_LogMass','LogMet','e_LogMet','LogSSFR','e_LogSSFR',
                'GFF','Source']
        for i,col in enumerate(meta.colnames):
            # print col,cols[i]
            meta.rename_column(col,cols[i])
        meta['Name'] = ['SDSS' + str(x) for x in meta['CID']]
        # print len(meta)
        meta = meta[meta['Type']=='S']
        # print len(meta)

    return meta

def write_lc(samples = None, datafolder = '/home/mi/Desktop/project-data',outprefix='allsne'):
    if samples == None:
        samples = ['CSP','Swift','CfA1','CfA2','CfA3','CfA4','CfAIR2','Pan-STARRS','SDSS','SNLS','Hamuy96','Foundation']    
    num = 0
    data_all = Table()
    for s in samples:
        print "Reading in sample:", s
        data = read_lc(s,datafolder)
        print "Columns:", data.colnames
        print "Number of SNe in the sample:", len(set(data['Name'].tolist()))
        num = num + len(set(data['Name'].tolist()))
        data['Survey'] = s
        data['Name_upper'] = map(str.upper,data['Name'])
        data_all = vstack([data_all,data])
        f = os.path.join('/home/mi/Desktop/project-sn-fitter/project/data',s+'.txt')
        print "Writing out to file:",f
        data.write(f,format='ascii',comment=False,overwrite=True)
    print "Total num of SNe in all samples:", len(set(data_all['Name_upper'].tolist()))
    f_all = os.path.join('/home/mi/Desktop/project-sn-fitter/project/data',outprefix+'.txt')
    print "Writing all SNe to: ", f_all
    data_all.write(f_all,format='ascii',comment=False,overwrite=True)

    return #data_all



def write_sn_info(samples = None, datafolder = None, outprefix='allsne'):
    if samples == None:
        samples = ['CSP','Swift','CfA1','CfA2','CfA3','CfA4','CfAIR2','Pan-STARRS','SNLS','SDSS','Foundation','Hamuy96']
    if datafolder == None:
        datafolder = '/home/mi/Desktop/project-data'  
    num = 0
    data_all = Table()
    for s in samples:
        print "Reading in info for sample:", s
        data = read_sn_info(s,datafolder)
        print "Columns:", data.colnames
        print "Number of SNe in the sample:", len(set(data['Name'].tolist()))
        num = num + len(set(data['Name'].tolist()))
        data['Sample'] = s
        data['Name_upper'] = map(str.upper,data['Name'])
        data_all = vstack([data_all,data])
        f = os.path.join('/home/mi/Desktop/project-sn-fitter/project/data/meta/',s+'_sninfo.txt')
        print "Writing out to file:",f
        data.write(f,format='ascii',comment=False,overwrite=True)
    print "Total num of SNe in all samples:", len(set(data_all['Name_upper'].tolist()))
    f_all = os.path.join('/home/mi/Desktop/project-sn-fitter/project/data/meta',outprefix+'_sninfo.txt')
    print "Writing all SNe to: ", f_all
    data_all.write(f_all,format='ascii',comment=False,overwrite=True)

    return #data_all



def write_host_info(samples = None, datafolder = None, outprefix='allsne'):
    if samples == None:
        samples = ['Pan-STARRS','SNLS','SDSS','Wolf','Neill','Uddin']  
    if datafolder == None:
        datafolder = '/home/mi/Desktop/project-data'  
    num = 0
    data_all = Table()
    for s in samples:
        print "Reading in info for sample:", s
        data = read_sn_info(s,datafolder)
        print "Columns:", data.colnames
        print "Number of SNe in the sample:", len(set(data['Name'].tolist()))
        num = num + len(set(data['Name'].tolist()))
        data['Sample'] = s
        data['Name_upper'] = map(str.upper,data['Name'])
        data_all = vstack([data_all,data])
        f = os.path.join('/home/mi/Desktop/project-sn-fitter/project/data/host/',s+'_host.txt')
        print "Writing out to file:",f
        data.write(f,format='ascii',comment=False,overwrite=True)
    print "Total num of SNe in all samples:", len(set(data_all['Name_upper'].tolist()))
    f_all = os.path.join('/home/mi/Desktop/project-sn-fitter/project/data/host',outprefix+'_host.txt')
    print "Writing all SNe to: ", f_all
    data_all.write(f_all,format='ascii',comment=False,overwrite=True)

    return #data_all


def rename_sn(data,rootfolder = '/home/mi/Desktop/project-sn-fitter/project/data/', 
              flist='rename.list',outfile='allsne.txt'):
    f = os.path.join(rootfolder,flist)
    rename_list = Table.read(f,format='ascii')
    data['Name_alt'] = data['Name']
    for i,name in enumerate(rename_list['Name2']):
        idx = data['Name'] == name
        data['Name'][idx]= [rename_list['Name1'][i]]*np.sum(idx)
    f_out = os.path.join('/home/mi/Desktop/project-sn-fitter/project/data',outfile)
    data['Name_upper'] = map(str.upper,data['Name'])
    print "Total num of SNe after rename:", len(set(data['Name_upper'].tolist()))
    print "Writing SNe to: ", f_out
    data.write(f_out,format='ascii',comment=False,overwrite=True)

    return #data