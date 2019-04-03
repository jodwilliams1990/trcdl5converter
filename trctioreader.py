"""
Executable for processing the R1 waveforms, and storing the reduced parameters
into a HDF5 file, openable as a `pandas.DataFrame`.
"""
import sys
sys.path.append('../')
sys.path.append('../../')
import numpy as np
import h5py
import target_io2
import read_lecroy_trc as trc
import glob
import os
import matplotlib.pyplot as plt
import argparse
from argparse import ArgumentDefaultsHelpFormatter as Formatter
import pandas as pd
from tqdm import tqdm
import json
from CHECLabPy.core.io import ReaderR1, DL1Writer
#from CHECLabPy.core.factory import WaveformReducerFactory
#from CHECLabPy.utils.waveform import BaselineSubtractor
from inspect import isabstract
from CHECLabPy.core.base_reducer import WaveformReducer
from CHECLabPy.core.spectrum_fitter import SpectrumFitter
from CHECLabPy.core.base_reducer import WaveformReducer
from scipy import interpolate
from scipy.ndimage import correlate1d
import operator
path =r'/users/chec/Desktop/TRC to TIO converter'
path2=r'/Users/chec/Software/CHECLabPy/CHECLabPy/data/checs_reference_pulse.txt'
fileprefix='C2new0_0_RC_-highTrans00000'
test=0
if test==1:
    samplelength=128
else:
    samplelength=1252
'''
for filename in glob.glob(os.path.join(path,fileprefix + '*.trc')):
    all_time, all_amplitude, d = trc.readTrc(filename)

print(len(all_time))
plt.plot(all_time, all_amplitude)
plt.show()
'''

def starter(self, n_pixels, n_samples, plot=False, reference_pulse_path='', **kwargs):
    super().__init__(n_pixels, n_samples, plot, **kwargs)
    ref = self.load_reference_pulse(reference_pulse_path)
    self.reference_pulse, self.y_1pe = ref
    self.cc = None

#@staticmethod
def load_reference_pulse(path2):
    file = np.loadtxt(path2)
    print("Loaded reference pulse: {}".format(path2))
    time_slice = 1E-9
    refx = file[:, 0]
    refy = file[:, 1]
    f = interpolate.interp1d(refx, refy, kind=3)
    max_sample = int(refx[-1] / time_slice)
    x = np.linspace(0, max_sample * time_slice, max_sample + 1)
    y = f(x)
    #plt.plot(y)
    #plt.show()
    pad = y.size - 2 * np.argmax(y)        # Put pulse in center so result peak time matches with input peak
    if pad > 0:
        y = np.pad(y, (pad, 0), mode='constant')
    else:
        y = np.pad(y, (0, -pad), mode='constant')
    y_1pe = y / np.trapz(y)                # Create 1p.e. pulse shape
    #print('0', max(y), np.trapz(y))
    #print('1', max(y_1pe), np.trapz(y_1pe))
    y = y / correlate1d(y_1pe, y).max()    # Make maximum of cc result == 1
    #plt.plot(y_1pe, label='y_1pe')
    #plt.plot(y, label='y')
    #plt.legend(loc='best')
    #plt.show()
    return {'y':y, 'y_1pe':y_1pe}
'''
def get_reference_pulse_at_t(self, t):
    ref_pad = np.pad(self.reference_pulse, self.n_samples, 'constant')
    ref_t_start = ref_pad.size // 2
    ref_t_end = ref_t_start + self.n_samples
    if t > self.n_samples:
        raise IndexError
    start = ref_t_start - t
    end = ref_t_end - t
    return ref_pad[start:end]

def _apply_cc(self, waveforms):
    cc = correlate1d(wfp, y_1pe)
    return cc

def _set_t_event(self, waveforms):
    self.cc = self._apply_cc(waveforms)
    super()._set_t_event(self.cc)
'''
def _get_charge(wft, wfp, y, y_1pe): #, waveforms):

    cext=1
    cc = correlate1d(wfp, y)  # _apply_cc(self, waveforms)
    if cext==1:     #JASONS
        #_set_t_event(self, waveforms)
        ind=np.argmax(wfp)
        charge = cc[ind]        #cc[:, self.t_event]   # <----- SORT THIS       THIS IS THE NUMBER OF PE/mVns
        cc_height = charge * y_1pe.max()    #HOW MANY PE x 1PE
    elif cext==2:
        charge=np.trapz(wfp)
        c2=np.trapz(y_1pe)
        cc_height = charge * y_1pe.max()  # /c2             #get_pulse_height(self,charge)
    elif cext==3:
        ind=np.argmax(cc)
        print(cc[ind])
        plt.plot(y_1pe[0:127]*cc)
        plt.plot(y_1pe, label='1pe')
        plt.plot(cc, label='cc')
        plt.legend(loc='best')
        plt.show()

    indj, valuea = max(enumerate(charge*y_1pe), key=operator.itemgetter(1))
    interpolate = np.zeros((5, 2))
    for s in range(0, 5):
        interpolate[s, 0] = indj + s - 2  # index[v] is the time of peak measured on line 52
        interpolate[s, 1] = charge*y_1pe[int(interpolate[s, 0])]  # voltage at time index[v] = pulse height
    other_x = np.linspace(indj - 2, indj + 2, 500)
    z = np.polyfit(interpolate[:, 0], interpolate[:, 1], 10)
    f = np.poly1d(z)
    b1 = f(other_x)
    tpulse, peak = max(enumerate(b1), key=operator.itemgetter(1))
    tpulse=other_x[tpulse]

    #ADD OTHER PARAMETERS HERE

    return {'charge':charge, 'cc_height':cc_height, 'tpulse':tpulse}

def main():
    for filename in glob.glob(os.path.join(path,fileprefix + '*.trc')):
        all_time, all_amplitude, d = trc.readTrc(filename)
        #print(d)
        if test==1:
            file = np.loadtxt(path2)
            print("Loaded reference pulse: {}".format(path2))
            all_time = file[:, 0]
            all_amplitude = file[:, 1]
            plt.plot(all_amplitude)
            plt.show()
        n_events = len(all_time)/samplelength
        n_pixels = 1
        n_cells = (3*1000/50)**2
        n_samples = 1
        df = pd.DataFrame()
        input_path = fileprefix
        output_path = fileprefix + "_dl1.h5"
        #print(output_path)
        #with DL1Writer(output_path, n_events*n_pixels, args.monitor) as writer:
        #    t_cpu = 0
        #    start_time = 0
        print(n_events)

        for wf in range (0,int(n_events)):
            print(wf)
            wft=all_time[(samplelength*wf):((samplelength*wf)+(samplelength-1))]
            wfp=all_amplitude[(samplelength*wf):((samplelength*wf)+(samplelength-1))]
            wfp=wfp-wfp[0]    #BASELINE SUBTRACTION
            iev = n_events

            lrp=load_reference_pulse(path2)
            y = lrp['y']
            y_1pe = lrp['y_1pe']
            #print(len(y),len(y_1pe))
            #_get_charge(wft, wfp, y_1pe)
            params = _get_charge(wft,wfp,2*y, y_1pe) #reducer.process(waveforms_bs)
            df_ev = pd.DataFrame(dict(
                wf=wf,
                iev=iev,
                pixel=1,
                cc_height=params['cc_height'],
                charge=params['charge'],
                tpulse=params['tpulse'],
                samplelength=samplelength,
                n_cells=n_cells,
                date_generated=pd.datetime.now(),
                input_path=input_path,
            ), index=[0])
            print(df)
            df=df.append(df_ev, ignore_index=True)
    df.to_hdf(output_path,key='df',mode='w')

if __name__ == '__main__':
    main()
    reader=pd.read_hdf(path+'/C2new0_0_RC_-highTrans00000_dl1.h5')
    wf=reader['wf']
    pixel=reader['pixel']
    charge=reader['charge']
    amp_pulse=reader['cc_height']
    tpulse=reader['tpulse']
    samplelength=['samplelength']
    print(charge, amp_pulse, tpulse)