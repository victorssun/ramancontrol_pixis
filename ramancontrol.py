#!/usr/bin/env python
# Only takes Pixis images
# slit width: each division is 10 um, with numerical are 1 mm. min(0.00), max(3.00)m

from scipy import *
import pylab, pvcam, serial, sys, datetime, pickle
from serial import Serial
from scipy.interpolate import interp1d
import numpy as np
import time
import matplotlib.pylab as plt


'''
select_monomove()
get_monomove()
get_wlrange()
initalize_pixis()
class raman()
    close_pixis()
    change_binning()
    change_integration()
    change_averages()
    change_monomove()
    acquire_pixis()
    acquire_img()
    acquire_spec()
    acquire_plot()
    save_data()
    load_bg()
    plot_img()
    plot_spec()
    preview_live()
'''

def select_monomove(wl_centre=650, no_grating=0):
    # STILL REQUIRES INITILIZE BASH SCRIPT IF TURNED OFF
    # input centre wavelength and grating
    # 1: 1200 gr/mm (bl 700), 2: 600 gr/mm (bl 1000), 3: 300 gr/mm (bl 1000)
    if wl_centre == 0:
        input_monomove = raw_input(' input no_grating & wl_centre (ex. 690 3): ')
        wl_centre, no_grating = input_monomove.split(' ')
    else:
        no_grating = str(no_grating)
        wl_centre = str(wl_centre)
    # setup
    sp = Serial('/dev/ttySP', baudrate=9600, timeout=1) # Serial('/dev/ttyUSB0' ...)
    sp.readlines()
    # switch grating/wl
    sp.write('%s GRATING\r' %no_grating)
    sp.readlines()
    time.sleep(15)
    # switch centre wavelength
    sp.write('%s GOTO\r' %wl_centre)
    sp.readlines()
    # query new settings
    sp.write('?GRATING\r')
    gn = sp.readlines()
    gn = gn[-1].split()[0]
    sp.write('?NM\r')
    centre = sp.readlines()
    centre = int(float(centre[0].split()[0]) + 0.5)
    sp.close()
    print(' centre: ' + str(centre) + ', grating: ' + str(gn))

def get_monomove():
    # see where the monochromator is, to record info
    # gn, grating, centre, wl, wn information
    sp = Serial('/dev/ttySP', baudrate=9600, timeout=1)
    sp.readlines()
    sp.write('?GRATING\r')
    gn = sp.readlines()
    gn = gn[-1].split()[0]
    sp.write('?NM\r')
    centre = sp.readlines()
    centre = int(float(centre[0].split()[0]) + 0.5)
    sp.close()
    if gn == '1':
        grating = '1200 gr/mm'
        calFile = '/mnt/cluster-victor/lin_cam/1200grDispersion.dat'
    elif gn == '2':
        grating = '600 gr/mm'
        calFile = '/mnt/cluster-victor/lin_cam/600grDispersion.dat'
    elif gn == '3':
        grating = '300 gr/mm'
        calFile = '/mnt/cluster-victor/lin_cam/300grDispersion.dat'
    else:
        print('unknown grating\n')
        sys.exit(1)
    print(' grating: %s,  centre: %s' %(grating, centre)) # just to make sure it's correct
    # based on mono setting, figure out the pixel wavelengths
    wl_c, wl_low, wl_high = loadtxt(calFile, unpack=True)
    f = interp1d(wl_c, wl_low)
    x_low = f(centre)
    f = interp1d(wl_c, wl_high)
    x_high = f(centre)
    wl = linspace(x_low, x_high, 1340) # full range of wavelengths to plot
    # determine Stokes Raman shift in cm-1
    freq = 10000./(wl/1000.)
    stokes = 10000./.6328 - freq # full range of stokes shift to plot
    return gn, grating, centre, wl, stokes

def get_wlrange(wl_centre, no_grating):
    # get stokes shift range before actually moving monochromator
    from scipy.interpolate import interp1d
    gn = str(no_grating)
    centre = int(wl_centre)
    if gn == '1':
        calFile = '/mnt/cluster-victor/lin_cam/1200grDispersion.dat'
    if gn == '2':
        calFile = '/mnt/cluster-victor/lin_cam/600grDispersion.dat'
    if gn == '3':
        calFile = '/mnt/cluster-victor/lin_cam/300grDispersion.dat'
    else:
        print(' unknown grating')
    wl_c, wl_low, wl_high = loadtxt(calFile, unpack=True)
    f = interp1d(wl_c, wl_low)
    x_low = f(centre)
    f = interp1d(wl_c, wl_high)
    x_high = f(centre)
    print(' %d to %d nm' %(x_low, x_high))
    freq = 10000./(x_low/1000.)
    stokes_low = 10000./.6328 - freq
    freq = 10000./(x_high/1000.)
    stokes_high = 10000./.6328 - freq
    print(' %d to %d cm-1 shift' %(stokes_low, stokes_high))

def initialize_pixis(integration=20, averages=1):
    # set up pixis, set exposure time
    a = pvcam.model.DigitalCamera # do you need to make this into a variable?
    pixis = pvcam.PVCam(a, 0, 0)
    pixis.exposureTime.value = integration/1000. 
    #print pixis.exposureTime.value
    #pixis.resolutionFitter([100,200])
    # set up monochromator information
    print(' integration: %d, averages: %d' %(integration, averages))
    return pixis, a, integration, averages # , gn, grating, centre, wl, stokes, averages

class raman:
    # encloses pixis into raman class
    # includes: close_pixis, binning, acquire_img, acquire_spec, save_data, load_bg
    def __init__(self, direct, integration, averages, binning=400, bg_subtract=False):
        self.direct = direct # where to load/save data
        self.gn, self.grating, self.centre, self.wl, self.stokes = get_monomove()
        self.pixis, self.a_before, self.integration, self.averages = initialize_pixis(integration, averages)
        self.change_binning(binning) # does this work?
        self.settings = [self.gn, self.grating, self.centre, self.integration, self.averages, self.stokes] # ex. [3, 1200, 633, 500, 3, [range]]
        self.bg_subtract = bg_subtract

    def close_pixis(self):
        self.pixis.terminate()
        print(' pixis closed')
        return
    
    def change_binning(self, binning_value):
        self.pixis._binning = (1, binning_value) # bins x (never use), bins y (should either be 1 or 400, what happens when between the two?)
        self.binvalue = binning_value
        '''
        if self.binvalue == 400:
            print(' Raman spectrum selected.')
        elif self.binvalue == 1:
            print(' Raman image selected.')
        else:
            print(' not 400 or 1.')
        '''
        return 

    def change_integration(self, integration_value):
        self.pixis.exposureTime.value = integration_value/1000. # will this work?
        self.integration = integration_value
        print(' integration time: %d ms' %integration_value)
        return 
    
    def change_averages(self, averages_value):
        self.averages = averages_value
        print( ' averages: %d' %averages_value)
        return 

    def change_monomove(self, wl_centre, no_grating):
        # change monochromator
        select_monomove(wl_centre, no_grating)
        self.gn, self.grating, self.centre, self.wl, self.stokes = get_monomove()
        return

    def acquire_pixis(self):
        # NOT USING ATM
        # acquires pixis image/spectrum based on binning. Averaging only works for spectra
        #pixis._image_rect = (0, 1339, 0, 399)
        if self.averages == 1:
            im = self.pixis.data.get()
            if self.binvalue == 400: # complete y binning, you can just grab directly from within, so don't have to do it in later stages
                im = np.array(im[0]) # data in the first item in a one-item list(ince all being binned to one vector), converted from odemis to numpy array
            return im
        else:
            total = zeros(1340)
            for x in range(self.averages):
                im = self.pixis.data.get()
                total += im[0]
                total = np.array(total)
            return total/self.averages
    
    def acquire_img(self, save_img=False):
        # acquire img
        # can't do averages/bg or try converting to array then doing averages/bg?
        if self.binvalue == 1:
            im = self.pixis.data.get()
            self.data_img = np.array(im)
            if save_img == True:
                comment = raw_input(' saving img, add comment?')
                self.save_data(comment, is_img=True)
            elif isinstance(save_spec, basestring):
                self.save_data(save_spec)
        else:
            print(' error: binning is not 1.')
        return 

    def acquire_spec(self, save_spec=False, is_bg=False):
        # take pixis spectrum, binning must be 400
        # also can acquire bg
        if self.binvalue == 400:
            if self.averages == 1:
                im = self.pixis.data.get()[0] # does this work?
            else:
                total = np.array(zeros(1340))
                for x in range(self.averages):
                    im_temp = self.pixis.data.get()[0]
                    total += np.array(im_temp)
                    im = total / self.averages #np or no np array?
            im = np.array(im)
            if is_bg == False:
                self.data = im # save variable into self
                if save_spec == True:
                    comment = raw_input(' saving data, add comment? ')
                    self.save_data(comment) # default would either be 1 or blank/space
                elif isinstance(save_spec, basestring):
                    self.save_data(save_spec)
            elif is_bg == True:
                self.bg = im # save variable into self
                self.bg_subtract = True
                if save_spec == True:
                    comment = raw_input(' saving bg, add comment [atleast add bg]? ')
                    self.save_data(comment)
                elif isinstance(save_spec, basestring):
                    self.save_data(save_spec)
        else:
            print(' error: binning is not 400.')
        return 

    def acquire_plot(self, save_spec=False, is_bg=False):
        # acquire and plot img or spectrum
        if self.binvalue == 400:
            self.acquire_spec(save_spec, is_bg)
            self.plot_spec()
        elif self.binvalue == 1:
            self.acquire_img(save_spec)
            self.plot_img()

    def filenaming_raman(self, comment):
        # filename with raman parameters
        filename = '%sd%s_g%s_%snm_%dms_%davgs_%s.pkl' %(self.direct,
                datetime.datetime.now().strftime('%Y%m%d-%H%M%S'),
                self.gn, self.centre, self.integration, self.averages, comment)
        return filename

    def save_data(self, comment, is_img=False):
        # save pixis data and xrange, with spectrometer parameters, possible comment as pickle
        filename = self.filenaming_raman(comment)
        '''
        filename = '%sd%s_g%s_%snm_%dms_%davgs_%s.pkl' %(self.direct, 
                datetime.datetime.now().strftime('%Y%m%d-%H%M%S'), 
                self.gn, self.centre, self.integration, self.averages, comment)
        '''
        # pickle version
        if is_img == False:
            pickle.dump([[self.gn, self.centre, self.integration, self.averages, self.wl], self.stokes, self.data], open(filename, 'wb')) # [spectrometer parameters + wl], stokes shift, data 
        elif is_img == True:
            pickle.dump([[self.gn, self.centre, self.integration, self.averages, self.wl], self.stokes, self.data_img], open(filename, 'wb')) # [spectrometer parameters + wl], stokes shift, data 
        print(' %s data saved.' %filename)
        return

    def load_bg(self, picklename):
        # import bg spectrum from pickle file
        # ex. d20180915-093025_g3_633nm_500ms_3avgs_comment.pkl
        filename = self.direct + picklename
        settings_bg, stokes_bg, data_bg = pickle.load(open(filename), 'rb')
       # check if bg is the correct, there should be no errors
        if settings_bg[2] != self.integration:
            print(' error: integration times different. Intensity will differ.\n')
        if settings_bg[3] != self.averages:
            print(' error: different averages. Noise will differ. \n')
        if settings_bg[0] != self.gn:
            print(' error: different grating. Stokes range will differ.\n')
        if settings_bg[1] != self.centre:
            print(' error: different centre. Stokes range will differ.\n')
        self.bg = data_bg
        self.bg_subtract = True
        return 
    
    def plot_img(self, fig=''):
        # plot image usually companied with spectrum
        if self.binvalue != 1:
            print(' possible error: bin value is not currently 1')
        if fig == '':
            fig = plt.figure() # not required if using preview_live()
        ax1 = fig.add_subplot(211)
        ax1.set_yticks(range(0, 401, 100))
        ax1.set_xlabel('wavelength / nm')
        ax1.set_ylabel('Y pixel number')
        ax1.imshow(self.data_img, extent=[min(self.wl), max(self.wl), 1, 400], aspect='auto')
        plt.show()
        return ax1 # is this required?

    def plot_spec(self, fig='', manual_yrange=False, dual=False):
        # plot spectrum
        if self.binvalue != 400:
            print(' possible error: bin value is not currently 400')
        if fig == '':
            fig = plt.figure() # not required if using preview_live()
        if dual == True:
            ax2 = fig.add_subplot(212) # if want both image and spectrum
        elif dual == False:
            ax2 = fig.add_subplot(212) # change this to 111 with nice graph later?
        ax2.set_xlabel('Stokes shift, $\Delta \omega$ / cm$^{-1}$')
        ax2.set_ylabel('counts')
        ax2.set_xlim(min(self.stokes), max(self.stokes))
        if manual_yrange == True:
            ylim_low = raw_input(' ylim low?: ')
            ylim_low = int(ylim_low)
            ylim_high = raw_input(' ylim high?: ')
            ylim_high = int(ylim_high)
        elif manual_yrange == False:
            ylim_low = min(self.data)
            ylim_high = max(self.data)
        else:
            ylim_low = manual_yrange[0]
            ylim_high = manual_yrange[1]
        ax2.set_ylim([ylim_low, ylim_high])
        if self.bg_subtract == False:
            ax2.plot(self.stokes, self.data)
        elif self.bg_subtract == True:
            ax2.plot(self.stokes, self.data - self.bg) # does bg subtract work?
        ax2.set_title('%s, %d nm, %d ms, %d avgs, %d max' %(self.grating, self.centre, self.integration, self.averages, max(self.data)))
        plt.show()
        return ax2

    def preview_live(self, dual=False, manual_yrange=False):
        # preview img+spectrum, or just spectrum. No save
        # preview sometimes gives core dumped sometimes
        pylab.ion()
        pylab.rc('font', size=9)
        #fig = pylab.figure(figsize=(3.25, 8))
        fig = pylab.figure()
        fig.set_tight_layout(True)
        try:
            while True:
                #a = time.time()
                if dual == True:
                    # acquire + plot image
                    pylab.clf() # clear fig?
                    self.change_binning(1)
                    self.acquire_img(save_img=False)
                    ax1 = self.plot_img(fig=fig)
                    self.change_binning(400)
                # acquire + plot spectrum
                self.acquire_spec(save_spec=False, is_bg=False)
                if self.bg_subtract == True:
                    self.data2 = self.data - self.bg
                
                ax2 = self.plot_spec(fig=fig, dual=dual, manual_yrange=manual_yrange)

                ax2.figure.canvas.draw()
                ax2.clear()

                #b = time.time()
                #print b - a

        except KeyboardInterrupt:
            if dual == True:
                ax1 = self.plot_img(fig=fig)
            ax2.clear()
            ax2 = self.plot_spec(fig=fig, dual=dual, manual_yrange=False)
            ax2.figure.canvas.draw()
            savefig = raw_input('\n save fig w/ descr. (ensure no overwrite): ')
            if savefig != '':
                # saves if input something
                filename = '%sd%s_g%s_%snm_%dms_%davgs_%s.png' %(self.direct, 
                        datetime.datetime.now().strftime('%Y%m%d-%H%M%S'), 
                        self.gn, self.centre, self.integration, self.averages, savefig)
                pylab.savefig(filename)
                print(' saved figure.')
                comment = raw_input(' saving spectrum, add comment? ')
                self.save_data(comment)
                comment = raw_input(' saving img, add comment? ')
                self.save_data(comment, is_img=True)
        return


        


### MAIN CODE ###
# init
direct = '/mnt/cluster-victor/lin_cam/new_data/'
binning = 400 # for spectrum analysis

# initial parameters
if raw_input(' change grating, centre wl (y/n)? ') == 'y':
    select_monomove(0) # only do if settings are different
if raw_input(' change integration, averages (y/n)? ') == 'y':
    input1 = raw_input(' input integration averages (ex. 500 1): ')
    integration, averages = input1.split(' ')
    integration = int(integration)
    averages = int(averages)
else:
    integration = 500
    averages = 1

# placeholder for testing
if raw_input(' initialize raman class (y/n)? ') == 'y':
    raman = raman(direct, integration, averages, binning)


# popular functions
#raman.acquire_plot() # along with acuire_spec(), plot_spec()
#raman.preview_live(dual=False, manual_yrange=False)
#raman.close_pixis() # always have at the end





