### Python version of spetrum mangling code by Saurabh

import numpy as np
import sncosmo
from .utils import get_refmag,plot_sed
import matplotlib.pyplot as plt
import copy
from operator import add
from scipy.interpolate import UnivariateSpline

class MangleResults(object):
  def __init__(self, insed, outsed, phase, minwave, maxwave, 
               scale, origdiff, filtnames):
     self.insed = insed
     self.outsed = outsed
     self.phase = phase
     self.minwave = minwave
     self.maxwave = maxwave
     self.scale = scale
     self.origdiff = origdiff
     self.filtnames = filtnames

def manglespec(phot,insed=None,redshift=0.,taperscale=3.,phaserange=[-50,100],imax=30,
               insed_is_rest=True,plot=True,spline=False,**kwargs):

    ## photometry should be astropy table with these columns: Mag,MagErr,Phase,Filter

    filternames = np.unique(phot['Filter'])

    ## insed needs to be sncosmo.Model()
    if not insed:
        insed = sncosmo.Model(source='hsiao')

    if insed_is_rest:
        insed.set(z=redshift)

    filters = []

    for f in filternames:
        band = sncosmo.get_bandpass(f)
        band.name = f
        filters.append(band)

    efflam = np.array([f.wave_eff for f in filters])
    minlam = np.array([f.minwave() for f in filters])
    maxlam = np.array([f.maxwave() for f in filters])

    sortid = np.argsort(efflam)

    filters_sorted = [filters[i] for i in sortid]

    print "Filters: ",[f.name for f in filters_sorted]
    efflam_sorted = efflam[sortid]
    print "efflam = ", efflam_sorted
    minlam_sorted = minlam[sortid]
    maxlam_sorted = maxlam[sortid]

    for f in filters_sorted:
        if f.name.split('-')[0] == 'SNLS':
            if 'rf' in kwargs:  
                rf = kwargs['rf']   
                print "re-register SNLS filters: rf= ", rf
                b = f.name.split('-')[-1]
                inst = 'megacampsf'
                band = sncosmo.get_bandpass(inst+"::"+b, rf)
                band.name = 'SNLS-'+b                
                sncosmo.register(band,force=True)
            else:
                print "rf needed for SNLS filters."

    phot.sort(keys='Phase')
    
    # for p in np.unique(phot['Phase']):
    if len(np.unique(phot['Phase']))>1:
        raise Exception("Only one phase can be processed at a time.")

    p = phot['Phase'][0]

    if p<phaserange[0] or p>phaserange[1]:
        raise Exception("Phase {0:.2f} out of range:".format(p),phaserange)
        # print "Phase {0:.2f} out of range:".format(p),phaserange
        # continue

    print "Phase =", p

    idx = phot['Phase'] == p
    phot_at_p = phot[idx]
    phot_at_p['sortid'] = 0
    for i,f in enumerate(filters_sorted):
        idx = phot_at_p['Filter'] == f.name
        phot_at_p['sortid'][idx] = i

    phot_at_p.sort(keys='sortid')

    print "Phot_at_p:"
    print phot_at_p['Phase','Filter','Mag','MagErr','sortid']

    modelmag0 = []
    for f in filters_sorted:
        mref,magsys = get_refmag(f.name.split('-')[0],f.name.split('-')[-1])
        # magsys = 'bd17'
        m = insed.bandmag(f,magsys,p) + mref  ##need zp? what's sncomso built in?
        modelmag0.append(m)
    # modelmag0 = np.array(modelmag0)+np.array(zp)

    print "modelmag0",np.array(modelmag0)
    diff0 = phot_at_p['Mag'] - modelmag0 
    print "phot_at_p['Mag'] - modelmag0:",np.array(diff0)

    # inwave = np.linspace(minlam.min(),maxlam.max(),200)
    inwave = np.linspace(insed.minwave(),insed.maxwave(),300)
    influx = insed.flux(p,inwave)

    ## change the normwave pivot and scale according to available photometry points
    wl = np.log10(efflam_sorted[0])
    wr = np.log10(efflam_sorted[-1])
    normpivot = np.power(10,np.mean([wl,wr]))
    wlmin = np.log10(minlam_sorted[0])
    wrmax = np.log10(maxlam_sorted[-1])
    normscale = 3./(wrmax-wlmin)

    effx = normwave(efflam_sorted/(1.+redshift),pivot=normpivot,scale=normscale)
    # print "effx:",effx
    vid = np.argmin(np.abs(effx))

    scale = np.power(10.,-0.4*(diff0[vid]))

    wave = inwave.copy()
    flux = np.atleast_2d(influx*scale)

    origflux = flux.copy()

    nrepeat = 5
    # print np.repeat(flux,n,axis=0).T

    phase = np.linspace(-50.,50.,nrepeat)

    newsource = sncosmo.TimeSeriesSource(phase,wave,np.repeat(flux,nrepeat,axis=0))
    newsed_at_p = sncosmo.Model(newsource)

    # plot_sed([newsed_at_p],day=[p], wave=[minlam.min(),maxlam.max()])
    # plt.show()

    modelmag = []
    for f in filters_sorted:
        mref,magsys = get_refmag(f.name.split('-')[0],f.name.split('-')[-1])
        # magsys = 'bd17'
        m = newsed_at_p.bandmag(f,magsys,p) + mref  ##need zp?
        modelmag.append(m)
    # modelmag = np.array(modelmag)+np.array(zp)

    print "modelmag:", np.array(modelmag)
    diff = phot_at_p['Mag'] - modelmag
    origdiff = diff.copy()
    print "phot_at_p['Mag'] - modelmag:",np.array(diff)
    print "diff/err:",(np.abs(diff/phot_at_p['MagErr'])).max()

    xx = normwave(wave/(1.+redshift),pivot=normpivot,scale=normscale)
    deg = len(filters_sorted)-1

    ## Start iteration

    for i in range(0,imax):
        print "iter =", i
        if (np.abs(diff/phot_at_p['MagErr'])).max()<1.:
            break
        else:        
            y = -0.4*diff*np.exp(taperscale*effx**2)

            ## polyfit/spline

            # polyfit
            if spline:
                spl = UnivariateSpline(effx,y)
                yy = spl(xx)
            else:
                polypars = np.polyfit(effx,y,deg)
                pfit = np.poly1d(polypars)
                yy = pfit(xx)

            ## some plots to check
            # plt.plot(effx,y,'o')
            # plt.plot(xx,yy)
            # plt.show()

            # plt.plot(xx,np.exp(taperscale*xx**2))
            # plt.show()

            # plt.plot(xx,np.exp(-taperscale*xx**2))
            # plt.show()

            lrfit = yy*np.exp(-taperscale*xx**2)
            newscale = np.power(10,lrfit)

            # plt.plot(wave,newscale,label='newscale')
            # plt.legend()
            # plt.show()

            flux = newsed_at_p.flux(p,wave)
            flux = np.atleast_2d(flux*newscale)

            newsource = sncosmo.TimeSeriesSource(phase,wave,np.repeat(flux,nrepeat,axis=0))
            newsed_at_p = sncosmo.Model(newsource)

            # plot_sed([newsed_at_p],day=[p], wave=[minlam.min(),maxlam.max()])
            # plt.show()

            modelmag = []
            for f in filters_sorted:
                
                mref,magsys = get_refmag(f.name.split('-')[0],f.name.split('-')[-1])
                # magsys = 'bd17'
                m = newsed_at_p.bandmag(f,magsys,p) + mref ##need zp?
                modelmag.append(m)
            # modelmag = np.array(modelmag)+np.array(zp)

            diff = phot_at_p['Mag'] - modelmag 

            print "modelmag:", np.array(modelmag)
            print "phot_at_p['Mag'] - modelmag:",np.array(diff)
            print "diff/err:",(np.abs(diff/phot_at_p['MagErr'])).max()

            ##plot outsed in each iteration                
            # plt.plot(wave,influx,label='insed')
            # # plt.plot(wave,origflux.ravel(),label='insed*scale')
            # flux = newsed_at_p.flux(p,wave)
            # plt.plot(wave,flux.ravel(),label='outsed')
            # plt.legend()
            # plt.show()

    ##fluxes are not scaled back to the right unit
    outflux = newsed_at_p.flux(p,wave)
    outflux = np.atleast_2d(outflux)
    outwave = wave/(1.+redshift)
    minwave = np.min(outwave)
    maxwave = np.max(outwave)

    influx = np.atleast_2d(influx)
    insource = sncosmo.TimeSeriesSource(phase,outwave,np.repeat(influx,nrepeat,axis=0))
    insed_plot = sncosmo.Model(insource)

    insource_scaled = sncosmo.TimeSeriesSource(phase,outwave,np.repeat(origflux,nrepeat,axis=0))
    insed_scaled = sncosmo.Model(insource_scaled)

    outsource = sncosmo.TimeSeriesSource(phase,outwave,np.repeat(outflux,nrepeat,axis=0))
    outsed = sncosmo.Model(outsource)

    if plot:
        fig, ax1 = plt.subplots()

        # ##plot sed in rest frame
        # print minlam.min(),maxlam.max()
        # plot_sed([insed_scaled,outsed],day=[p], 
        #           wave=[minlam.min()/(1.+redshift),maxlam.max()/(1.+redshift)],
        #           npoints=200,labels=['insed(scaled)','outsed'])

        ##plot sed in obs frame
        plot_sed([insed_scaled,outsed],day=[p], 
                  wave=[minwave,maxwave],
                  npoints=200,labels=['insed(scaled)','outsed'],z=redshift,restframe=False)

        ax2 = ax1.twinx()
        for f in filters_sorted:
            plotband = sncosmo.get_bandpass(f)
            bandwave = np.linspace(plotband.minwave(),plotband.maxwave(),50)
            ax2.plot(bandwave,plotband(bandwave),'--',color='grey',lw=0.8)

        # plt.show()

    filtnames = [f.name for f in filters_sorted] 

    return MangleResults(insed_scaled,outsed,p,minwave,maxwave,scale,origdiff,filtnames)


def normwave(wave,pivot=5500.,scale=3.):

    return np.log10(wave/pivot)*scale
