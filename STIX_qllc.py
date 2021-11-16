import astropy
from astropy.io import fits
import glob
import pandas as pd
from astropy.table import Table, vstack, hstack
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time.core import Time, TimeDelta
import numpy as np

en_band_index = 4 # index of which stix energy band you are interested in

import get_stix_bkg as bkg
this_date = '2021-07-17T00:00:00' # Please, use this format: 'YYYY-MMDDThh:mm:ss'
bkg_spectrum = bkg.stix_bkg(this_date,qllc=True)

eta = 2.5 * u.us
tau = 12.5 * u.us

def _lfrac(trigger_rate):
    nin = trigger_rate / (1 - (trigger_rate * (tau+eta)))
    return np.exp(-eta * nin) / (1 + tau * nin)

def stix_counts_timeseries_nobkg(en_band_index):
	plt.plot_date(data_concat['time'].plot_date, data_concat['counts'][:,en_band_index],'-', label=f' STIX {energies[en_band_index]["e_low"]}-{energies[en_band_index]["e_high"]} (without bkg)')
	plt.legend()
	plt.ylabel('Counts') #units of counts depends when you call the function.. could be Counts or Counts $s^{-1}keV^{-1}$ so change label accordingly
	plt.xlabel('Time')
	plt.show()

#stix_int_timeseries is for integrating counts over certain time period and taking the mean, the choice of i depends of the int_factor and length of array data_concat
def stix_int_timeseries(int_factor, en_band_index,i):
	for i in range(len(df_energies['Electron_Bins_Text'])):
		stix_counts_int = np.asarray([[np.mean(data_concat['counts'][k*int_factor:(k+1)*int_factor,j]) for j in range(5)] for k in range(i)])
		time_int = data_concat['time'][::int_factor]
		plt.plot_date(time_int.plot_date,stix_counts_int[:,en_band_index], '-', label = f' STIX {energies[en_band_index]["e_low"]}-{energies[en_band_index]["e_high"]} keV (w/o bkg)')
		plt.yscale('log')
		plt.gcf().autofmt_xdate() 
		plt.ylabel('Flux \n(STIX: Counts $s^{-1}keV^{-1}$')
		plt.xlabel('Time')
		plt.legend()
		plt.savefig(f'path_to_directory/{energies[en_band_index]["e_low"]}-{energies[en_band_index]["e_high"]}_{df_energies["Electron_Bins_Text"][i][0]}.png')
		plt.close()


hdulist = fits.open('solo_L1_stix-ql-lightcurve_20210507_V01.fits')
hdulist2 = fits.open('solo_L1_stix-ql-lightcurve_20210508_V01.fits')
hdulist3 = fits.open('solo_L1_stix-ql-lightcurve_20210509_V01.fits') 

header = hdulist[0].header
control = Table(hdulist[1].data)
data = Table(hdulist[2].data)
energies = Table(hdulist[4].data)

header2 = hdulist2[0].header
control2 = Table(hdulist2[1].data)
data2 = Table(hdulist2[2].data)

header3 = hdulist3[0].header
control3 = Table(hdulist3[1].data)
data3 = Table(hdulist3[2].data)

data2['time'] = data2['time']+data['time'][-1]
data3['time'] = data3['time']+data2['time'][-1]

data_concat = vstack([data,data2,data3],join_type = 'outer') #(lines 46-63) these steps are for combining/stacking several ql lightcurves, not necessary if just looking at one

energy_delta = energies['e_high'] - energies['e_low'] << u.keV 

live_time = _lfrac(data_concat['triggers'].reshape(-1)/(16*data_concat['timedel']*u.s))

bkg_array = [[bkg_spectrum[i] for i in range(5)] for j in range(len(data_concat))]

data_concat['counts'] = data_concat['counts'] - bkg_array # subtract_stix_bkg() using Andrea's file get_stix_bkg

data_concat['time'] = Time(header['date_obs']) + TimeDelta(data_concat['time'] * u.s)

data_concat['counts'] = data_concat['counts'] / (live_time.reshape(-1, 1) * energy_delta) #to divide by energy and time in units of Counts $s^{-1}keV^{-1}$


