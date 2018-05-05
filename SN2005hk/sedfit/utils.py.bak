import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table,vstack,unique,Row
import sfdmap
import sncosmo
import astropy.units as units
import errno
import os
import pymc.gp as gp
import random 
from scipy.optimize import minimize

def plot_sed(sedmodels, day=[-10.,0.,20.], wave=[3000.,7000.], 
             npoints=50,color=None,labels=None,scale=1.,z=0.,restframe=True):
    if color == None:
        color = plt.cm.jet(np.linspace(0,1,len(sedmodels)))
    for i,sed in enumerate(sedmodels):
        ts = day
        w = np.linspace(wave[0],wave[1],npoints)
        for t in ts:
            label = 't={0:.2f}'.format(t)
            if labels:
                label = label + '(' + labels[i] +')'
            else:
                label = label + '(' + str(i) + ')'
            if restframe:
                z = 0.
            plt.plot(w*(1.+z),sed.flux(t,w)/scale,label=label,color=color[i],alpha=1)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid(True)

def register_bands_for_sncosmo():
    #register bands in sncosmo
    #SDSS bands
    ffolder = '/home/mi/Desktop/project-sn-fitter/project/data/filters/'
    for b in ['u','g','r','i','z']:
        band = sncosmo.get_bandpass('sdss'+b)
        band.name = 'SDSS-'+b
        sncosmo.register(band,force=True)
        band.name = 'SDSS-'+b.upper()
        sncosmo.register(band,force=True)
    #CSP bands
    for b in ['u','g','r','i','B','V']:   
        if b == 'V':
            b = 'V9844'
        band = sncosmo.get_bandpass('csp'+b.lower())
        band.name = 'CSP-'+b[0:1]
        sncosmo.register(band,force=True)
        band.name = 'CSP-'+b[0:1]+'s'
        sncosmo.register(band,force=True)
        band.name = 'CSP-'+b[0:1]+'d'
        sncosmo.register(band,force=True)
    for b in ['Y','J','H']:
        for inst in ['Swope','DuPont']:
            band = sncosmo.get_bandpass('csp'+b.lower()+inst[0:1].lower())
            band.name = 'CSP-'+b+inst[0:1]
            sncosmo.register(band,force=True)
    #Standard bands(CfA1,CfA2,Hamuy96)
    for b in ['U','B','V','R','I']:
        band = sncosmo.get_bandpass('standard::'+b.lower())
        band.name = 'CfA1-'+b
        sncosmo.register(band,force=True)
        band.name = 'CfA2-'+b
        sncosmo.register(band,force=True)
        band.name = 'Hamuy96-'+b
        sncosmo.register(band,force=True)
        band.name = 'Standard-'+b
        sncosmo.register(band,force=True)

    #CfA3, CfA4 - Keplercam, 4Shooter
    for b in ['U','B','V','R','I']:
        inst = '4shooter2'
        if b in ["U","u'"]:
            b1 = 'Us'
            band = sncosmo.get_bandpass(inst+"::"+b1.lower())
        else:
            band = sncosmo.get_bandpass(inst+"::"+b.lower())
        band.name = 'CfA3-'+b[0:1].upper()+inst[0:1]
        sncosmo.register(band,force=True)
    for b in ['U','B','V',"r'","i'","u'"]:
        inst = 'keplercam'
        if b in ["U","u'"]:
            b1 = 'Us'
            band = sncosmo.get_bandpass(inst+"::"+b1.lower())
        else:
            band = sncosmo.get_bandpass(inst+"::"+b[0:1].lower())
        band.name = 'CfA3-'+b+inst[0:1]        
        sncosmo.register(band,force=True)
        band.name = 'CfA4-'+b
        sncosmo.register(band,force=True)
    #SNLS - Megacam, !!test only, focal radius = 0., focal radii need to be calculated sn by sn
    for b in ['g','r','i','z']:
        inst = 'megacampsf'
        band = sncosmo.get_bandpass(inst+"::"+b, 8.)
        band.name = 'SNLS-'+b                
        sncosmo.register(band,force=True)
    #CfAIR2
    for b in ['H','J','K_s']:
        band = sncosmo.read_bandpass(ffolder+'CfAIR2/'+b+'_RSR.txt',wave_unit=units.um)
        band.name = 'CfAIR2-'+b
        sncosmo.register(band,force=True)
    #Swift
    for b in ['U','B','V','UVW1','UVW2','UVM2']:
        band = sncosmo.read_bandpass(ffolder+'Swift/Swift_UVOT.'+b+'.dat')
        band.name = 'Swift-'+b
        sncosmo.register(band,force=True)
    #Pan-STARRS, Foundation
    bands = Table.read(ffolder+'/Pan-STARRS/Pan-STARRS_apj425122t3_mrt.txt',format='ascii.no_header')
    cols = ['wave','open','g','r','i','z','y','w','Aero','Ray','Mol']
    for i,col in enumerate(cols):
        bands.rename_column('col'+str(i+1),col)
    for b in ['g','r','i','z']:
        band = sncosmo.Bandpass(bands['wave'],bands[b],wave_unit=units.nm)
        band.name = 'Foundation-'+b
        sncosmo.register(band,force=True)
        band.name = 'Pan-STARRS-'+b
        sncosmo.register(band,force=True)


def get_refmag(survey,band): ## !!need to double check all the values here are correct
    # print survey,band
    if survey == 'SNLS': #Guy2010
        magsys = 'bd17'
        tmref = Table([['u','g','r','i','z'],
                      [9.7688,9.6906,9.2183,8.9142,8.7736]],
                      names=('band', 'mref'))
        mref = tmref['mref'][tmref['band']==band][0]
        return mref,magsys

    if survey == 'CSP': #http://csp.obs.carnegiescience.edu/data/filters
        if band[0:1] in ['u','g','r','i']:
            magsys = 'bd17'
            tmref = Table([['u','g','r','i'],
                          [10.518,9.644,9.352,9.25]],
                          names=('band', 'mref'))
            mref = tmref['mref'][tmref['band']==band[0:1]][0]
        elif band[0:1] in ['B','V']:
            magsys = 'vega'
            tmref = Table([['B','V'],
                          [0.03,0.0096]],
                          names=('band', 'mref'))
            mref = tmref['mref'][tmref['band']==band[0:1]][0]
        else:
            magsys = 'vega'
            mref = 0.
        return mref,magsys

    if survey in ['Pan-STARRS','Foundation','Pan']:
        return 0.,'ab'

    if survey == 'SDSS': #Sako 2014
        magsys = 'ab'
        tmref = Table([['u','g','r','i','z'],
                      [0.0679,-0.0203,-0.0049,-0.0178,-0.0102]],
                      names=('band', 'mref'))
        mref = tmref['mref'][tmref['band']==band.lower()][0]
        return mref,magsys

    if survey == 'CfA3':
        magsys = 'bd17'
        if band in ['U4','B4','V4','R4','I4']:
            tmref = Table([['U4','B4','V4','R4','I4'],
                          [9.693,9.8744,9.4789,9.1554,8.8506]],
                          names=('band', 'mref'))
            mref = tmref['mref'][tmref['band']==band][0]        
        elif band in ['UK','BK','VK',"r'K","i'K"]:
            tmref = Table([['UK','BK','VK',"r'K","i'K"],
                          [9.6922,9.8803,9.4722,9.3524,9.2542]],
                          names=('band', 'mref'))
            mref = tmref['mref'][tmref['band']==band][0]               
        return mref,magsys

    if survey == 'CfA4':
        magsys = 'bd17'
        tmref = Table([['U','B','V',"r'","i'"],
                      [9.6922,9.8803,9.4722,9.3524,9.2542]],
                      names=('band', 'mref'))
        mref = tmref['mref'][tmref['band']==band][0]        
        return mref,magsys

    if survey in ['CfA1','CfA2','Hamuy96','Standard']:
        magsys = 'bd17'
        tmref = Table([['U','B','V',"R","I"],
                      [9.724,9.907,9.464,9.166,8.846]],
                      names=('band', 'mref'))
        mref = tmref['mref'][tmref['band']==band][0]        
        return mref,magsys

    if survey == 'Swift':
        return 0.,'vega'

    else:
        return 0.,'ab'

def salt2_params(lcdata,metadata,snlist=None,modelcov=True,register_bands=False,
                 write_to_file=True,outfile='data/salt2params.dat',usebands='all'):
    if register_bands:
        register_bands_for_sncosmo()

    fit_result = Table()
    dustmap = sfdmap.SFDMap("/home/mi/sfddata-master")
    model = sncosmo.Model(source='salt2',
                          effects=[sncosmo.F99Dust()],
                          effect_names=['mw'],
                          effect_frames=['obs'])
    if snlist == None:
        snlist = list(set(lcdata['Name_upper']))
        snlist = np.sort(np.array(snlist))

    f = open(outfile,'w')
    for i,sn in enumerate(snlist):
        try:
            zp = 27.5
            lc = lcdata[lcdata['Name_upper'] == sn]
            meta = metadata[metadata['Name_upper'] == sn]

            lc['flux'] = np.power(10.,-0.4*(lc['Mag']-zp))
            lc['flux_err'] = np.absolute(0.921*lc['flux']*lc['MagErr'])
            lc['zp'] = zp
            lc['zpsys'] = 'ab'.ljust(10)
            a1 = lc['Filter']
            c = ("-" * len(lcdata))
            a2 = lc['Survey']
            lc['Instrument'].fill_value = ""
            inst = [x[0:1] for x in lc['Instrument'].filled()]
            lc['Filter'] = map(''.join, zip(a2, c, a1,inst))
            lc = lc[lc['flux']>0.]
            z = meta[0]['z_helio']
            mwebv = dustmap.ebv(meta[0]['RA'],meta[0]['DEC'])
            print "fitting ", sn, "z=", z, "mwebv=", mwebv
            model.set(z=z,mwebv=mwebv)
            surveys = meta[0]['Surveys'].split("::")
            for s in surveys:
                print 'Survey = ', s
                lc0 = lc[lc['Survey'] == s]

                if usebands != 'all':
                    lc0 = lc0[np.array([x.split('-')[-1][0:1] in list(usebands) for x in lc0['Filter']])]
                    if len(lc0) == 0:
                        print "No data in selected bands:",usebands
                        continue

                for b in set(lc0['Filter']):
                    # mref,magsys = get_refmag(s,b.split('-')[1])
                    mref,magsys = get_refmag(s,b[b.find(s)+len(s):].split('-')[-1])
                    lc0['zp'][lc0['Filter']==b] = lc0['zp'][lc0['Filter']==b]-mref
                    lc0['zpsys'][lc0['Filter']==b] = magsys

                if s == 'SNLS':
                    xf = meta[0]['xf']
                    yf = meta[0]['yf']
                    rf = np.sqrt(xf**2 + yf**2)
                    print "calculate rf and re-register SNLS filters: rf= ",rf
                    for b in ['g','r','i','z']:
                        inst = 'megacampsf'
                        band = sncosmo.get_bandpass(inst+"::"+b, rf)
                        band.name = 'SNLS-'+b                
                        sncosmo.register(band,force=True)

                result, fitted_model = sncosmo.fit_lc(lc0,model,['t0','x0','x1','c'],
                                                  modelcov=modelcov,
                                                  bounds={'x0':(0,1.),
                                                          'x1':(-5.,5.),
                                                          'c':(-3.,3.)})
                # print result
                res = []
                res.append(sncosmo.flatten_result(result))
                res = Table(res)
                res['Name'] = sn.ljust(30)
                res['Survey'] = s.ljust(20)
                res['ErrorType'] = 'None'.ljust(50)
                res['ErrorMesg'] = 'None'.ljust(50)
                # print res
                if i == 0:
                    res.write(f,format='ascii')
                else:
                    res.write(f,format='ascii.no_header')
                f.flush()
                os.fsync(f.fileno())
                fit_result = vstack([fit_result,res],join_type='outer')
            # print fit_result
        except Exception as e: 
            print 'An error occured:'
            print type(e)
            print e
            if len(fit_result) > 0:
                res = list([-99.9]*(len(fit_result.colnames)-4))+list([sn,s,type(e),e.message])
                fit_result.add_row(res)           
                Table(fit_result[-1]).write(f,format='ascii.no_header')
    f.close()

    return fit_result



class GPResults(object):
    def __init__(self, meanlc, meanlc_cov, peakmjds, filters):
        self.meanlc = meanlc
        self.meanlc_cov = meanlc_cov
        self.peakmjds = peakmjds
        self.filters = filters

def Constant(x,c):

    return 0.*x + c


def fit_gp(sn,diff_degree=3.,scale=30.,plot=True):
    
    if plot:
        figure = plt.figure(figsize=(15,9))

    meanlc = []
    meanlc_cov = []
    peakmjds = []
    filters = []

    for i,f in enumerate(set(sn['Filter'])):
#         print "Filter:", f
        
        if plot:
            plt.subplot(len(set(sn['Filter']))//6+1,6,i+1)
        
        sn_f = sn[sn['Filter']==f]
        
        y = np.array(sn_f['Mag'])
        yerr = np.array(sn_f['MagErr'])
        x = np.array(sn_f['MJD'])

#         idx = yerr < 5.*np.median(yerr)
        idx = y>0    
        
        if np.sum(idx) == 0:
            continue

        y = y[idx]
        yerr = yerr[idx]
        x = x[idx]

        mean = np.median(y) 
        amp = np.std(y-mean)

        M = gp.Mean(Constant,c=mean)
        C = gp.Covariance(gp.matern.euclidean, diff_degree=diff_degree,
                          amp=amp, scale=scale)

        gp.observe(M, C, obs_mesh=x, obs_vals=y, obs_V=np.power(yerr,2))

        pred_list = []
        for i in range(3):
            pred = gp.Realization(M,C)
            pred_list.append(pred)

        xa = np.atleast_1d(np.linspace(x.min(),x.max(),100))
        ya_err = np.sqrt(C(xa))

#         gp.plot_envelope(M, C, mesh=xa)
        
        if plot:
            plt.errorbar(x, y, yerr, fmt='b.')
            plt.plot(xa,M(xa),'r-')
            plt.fill_between(xa, M(xa)-ya_err, M(xa)+ya_err,color='red',alpha=0.5)

            for p in pred_list:
                ya = p(xa)
                plt.plot(xa, ya, lw =0.2, color='g')
            plt.title(f)
            plt.ylim(plt.ylim()[::-1])
            plt.tight_layout()
               
        meanlc.append(M)
        meanlc_cov.append(C)
                
        ##estimate peak in each band
        peakres = minimize(M,x[np.argmin(y)])
        pmjd = peakres.x[0]
        if pmjd < x.min():
            print "No data before the peak for band", f
            pmjd = -99.
        peakmjds.append(pmjd)
        
        filters.append(f)

    return GPResults(meanlc,meanlc_cov,peakmjds,filters)


## estimate peakmjd in restframe B band
def get_b_band_peakmjd(filters,peakmjds,redshift):

    restfilter_names = ['Bessell' + x for x in ['b','v']]    
    offset = [0., 1.]

    peakmjd = -99.
    for i,fname in enumerate(restfilter_names):

        restf = sncosmo.get_bandpass(fname)
        print "Effective wavelength of {}: {}".format(restf.name,restf.wave_eff)
        
        efflam = []
        
        for f in filters:
            b = sncosmo.get_bandpass(f)
            efflam.append(b.wave_eff)

        efflam = np.array(efflam)

        fid = np.argmin(np.abs(efflam/(1.+redshift)-restf.wave_eff))
        
        print "wavelength difference:",efflam[fid]/(1.+redshift)-restf.wave_eff
        # print efflam/(1.+redshift)

        if np.abs(efflam[fid]/(1.+redshift)-restf.wave_eff) > 1000.:
            continue
        else:
            if peakmjds[fid] == -99.:
                print "peak not available in ",filters[fid]
                continue
            else:
                peakmjd = peakmjds[fid] + offset[i]
                break

    return peakmjd

