# Paul Yollin

# This file contains a class and assorted functions used in the analysis of MEG recordings to investigate the
# presence of eddy current.

#   Required packages
#       -mne
#       -numpy
#       -matplotlib

#   Classes
#       -ReClass(file_name) (please note you must specify data_path locations and calibration_file location)

#   Functions
#       -find_norm(data)
#       -find_rms(data)
#       -find_drop(data)
#       -plot_norms(data,tmin,tmax)

import mne
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt


class RecClass:

    # This class loads in a raw file specific by the data path and file name, then applies SSS filtering using
    # the specific path to the calibration file. It then pulls out the magnetometer raw and filtered data along
    # with the gradiometer raw and filtered data. It stores these arrays in a dictionary called 'data'.

    # Path to the .fif file of the recording, and path to the SSS fine calibration file (machine specific!!)
    data_path = '/Users/cyoll/OneDrive/Documents/I-Labs/Research/Eddy_Current/eddy_current3/'
    calibration_file = '/Users/cyoll/OneDrive/Documents/I-Labs/Research/Eddy_Current/Analysis/sss_cal.dat'

    # initializing the object
    def __init__(self, file_name, start, stop):
        # reading in raw file and also fixing coil types
        self.raw = (mne.io.read_raw_fif(self.data_path + str(file_name),
                                        allow_maxshield=True,
                                        verbose='warning')).fix_mag_coil_types()
        # marking bad channels
        self.raw.info['bads'] = ['MEG1433', 'MEG1842', 'MEG1743']
        # Applying SSS without using the calibration file
        self.rawSSS = mne.preprocessing.maxwell_filter(self.raw,
                                                       coord_frame='meg',
                                                       verbose='warning')
        # Applying SSS using the calibration file
        self.rawSSSc = mne.preprocessing.maxwell_filter(self.raw,
                                                        coord_frame='meg',
                                                        verbose='warning',
                                                        calibration=self.calibration_file)
        # creates data arrays holding the recordings for easier analysis
        self.data = {'r_mdat': self.raw.get_data(start=start*1000, stop=stop*1000, picks='mag'),
                     'r_gdat': self.raw.get_data(start=start*1000, stop=stop*1000, picks='grad'),
                     # 'f_mdat': self.rawSSS.get_data(start=start*1000, stop=stop*1000, picks='mag'),
                     # 'f_gdat': self.rawSSS.get_data(start=start*1000, stop=stop*1000, picks='grad'),
                     'fc_mdat': self.rawSSSc.get_data(start=start*1000, stop=stop*1000, picks='mag'),
                     'fc_gdat': self.rawSSSc.get_data(start=start*1000, stop=stop*1000, picks='grad'),
                     'r_whole': self.raw.get_data(start=start*1000, stop=stop*1000),
                     'f_whole': self.rawSSS.get_data(start=start*1000, stop=stop*1000),
                     'fc_whole': self.rawSSSc.get_data(start=start*1000, stop=stop*1000)}
        # setting the bad gradiometer channels from the raw gradiometere array to zero
        self.data['r_gdat'][104,:] = 0
        self.data['r_gdat'][130, :] = 0
        self.data['r_gdat'][139, :] = 0


def find_norm(data):



    # This function takes in the data dictionary created from a RecClass object and returns a dictionary called
    # norms. These norms are the signal vector norms at a given point in time, so if your recording is 10,000 samples
    # then the norm vectors will be of length 10,000.

    # creating empty dictionary to hold norm vectors
    norms = {}

    # for loop going through and getting the norm vectors for all the data arrays in the RecClass.data object.
    # r_mdat means raw magnetometer data, f_gdat means filtered gradiometer data, and
    # fc_mdat means filtered using fine calibration file magnetometer data
    # for i in ['r_mdat', 'r_gdat', 'f_mdat', 'f_gdat', 'fc_mdat', 'fc_gdat']:
    #     norms[i] = LA.norm(data[i], axis=0)
    # del i

    for i in ['r_mdat', 'r_gdat', 'fc_mdat', 'fc_gdat']:
        norms[i] = LA.norm(data[i], axis=0)
    del i

    return norms


def find_rms(data):

    # This function takes data from RecClass.data and creates a dictionary containing rms values for each norm
    # vector created in the find_norm() function.

    # getting norms using find_norm() function
    norms = find_norm(data)

    # creating empty dictionary to hold rms values
    rms = {}

    # for loop going through and getting the rms values
    # rms is calculated by dotting norm vector with itself then taking the sqrt of that
    # for i in ['r_mdat', 'r_gdat', 'f_mdat', 'f_gdat', 'fc_mdat', 'fc_gdat']:
    #     rms[i] = np.sqrt(np.dot(norms[i], norms[i]))
    # del i

    for i in ['r_mdat', 'r_gdat', 'fc_mdat', 'fc_gdat']:
        rms[i] = np.sqrt(np.dot(norms[i], norms[i]))
    del i

    return rms


def find_drop(data):

    # this function takes the rms values calculated in the find_rms function and then uses those
    # to calculate the reduction factor. This is calculated by taking the raw rms over the filtered rms
    # doing this for both filtered using a fine calibration file and filtered without a fine
    # calibration file.

    #   running the find_rms() function defined above
    rms = find_rms(data)

    # creating an empty dictionary to hold calculated drops
    drop = {}

    # for loop going through and finding the drop for both magnetometers and gradiometers
    # using a fine calibration file and not using a fine calibration file
    for i in ['m', 'g']:

        #   defining name to be used to index into the rms dictionary
        # name = i + 'dat'

        #   calculating raw rms over filtered rms without calibration file
        # drop[name] = rms['r_' + i + 'dat'] / rms['f_' + i + 'dat']

        #   defining name to index into rms dict
        name_c = i + 'cdat'

        #   calculating raw rms over filtered rms with calibration file
        drop[name_c] = rms['r_' + i + 'dat'] / rms['fc_' + i + 'dat']

    return drop


def plot_norms(data, tmin, tmax):

    # sorry this plot function is messy, will document it better soon...

    norms = find_norm(data)
    fs = 1 / 1000
    start = int(tmin / fs)
    stop = int(tmax / fs)
    rec_len = len(norms['r_mdat'])
    xs = np.arange(0, rec_len) * fs

    fig, ax = plt.subplots(2, 1, figsize=(13, 9))

    for j, i in enumerate(['m', 'g']):
        name = i+'dat'
        line = ax[j].plot(
            xs[start:stop], norms['r_' + name][start:stop], 'k',
            #xs[start:stop], norms['f_' + name][start:stop], 'g',
            xs[start:stop], norms['fc_' + name][start:stop], 'b',
            linewidth=1,
        )
        ax[j].set(ylabel='Signal Norm')
        ax[j].grid()
        line[0].set_label('Raw Signal')
        line[1].set_label('Filtered w/ cal file Signal')

        #line[2].set_label('Filtered w/ cal file Signal')
        ax[j].legend(loc='upper right')
    ax[0].set(title='Raw vs filtered signal vector norms over time')
    ax[1].set(xlabel='Time (s)')


def ave_dat(data):

    rec_len = int(data.shape[1] / 1000)
    ns = np.arange(1, rec_len)

    int_mat = data[:, 0:1000]
    for n in ns:
        mat = int_mat + data[:, n * 1000:(n + 1) * 1000]
        int_mat = mat
    del n

    ave_dat = np.divide(mat, rec_len)

    return ave_dat


def ave_raw(data, info):

    a_dat = ave_dat(data)

    ave_raw = mne.io.RawArray(a_dat, info)

    return ave_raw


def ave_norm(data, info):

    ave_norms = {}

    rec_len = int(data['r_mdat'].shape[1] / 1000)
    ns = np.arange(1, rec_len)

    a_raw = ave_raw(data['r_whole'], info)
    a_sss = mne.preprocessing.maxwell_filter(a_raw, coord_frame='meg', verbose='warning')

    data['fc_mdat'] = a_sss.get_data(picks='mag')
    data['fc_gdat'] = a_sss.get_data(picks='grad')

    for j in ['fc_mdat', 'fc_gdat']:
        ave_norms[j] = LA.norm(data[j], axis=0)
    del j

    for i in ['r_mdat', 'r_gdat']:

        int_mat = data[i][:, 0:1000]
        for n in ns:
            mat = int_mat + data[i][:, n * 1000:(n + 1) * 1000]
            int_mat = mat
        del n

        a_dat = np.divide(mat, rec_len)
        ave_norms[i] = LA.norm(a_dat, axis=0)

    return ave_norms

def ave_drop(data, info):

    ave_drops = {}

    a_norms = ave_norm(data, info)

    rms = {}

    for i in ['r_mdat', 'r_gdat', 'fc_mdat', 'fc_gdat']:
        rms[i] = np.sqrt(np.dot(a_norms[i], a_norms[i]))
    del i

    drop = {}

    for i in ['m', 'g']:
        name_c = i + 'cdat'
        drop[name_c] = rms['r_' + i + 'dat'] / rms['fc_' + i + 'dat']
    del i

    return drop


def plot_ave_norms(data, info):

    # sorry this plot function is messy, will document it better soon...

    ave_norms = ave_norm(data, info)
    fs = 1 / 1000
    rec_len = len(ave_norms['r_mdat'])
    xs = np.arange(0, rec_len) * fs

    fig, ax = plt.subplots(2, 1, figsize=(13, 9))

    for j, i in enumerate(['m', 'g']):
        name = i+'dat'
        line = ax[j].plot(
            xs, ave_norms['r_' + name], 'k',
            xs, ave_norms['fc_' + name], 'b',
            linewidth=1,
        )
        ax[j].set(ylabel='Averaged Signal Norm')
        ax[j].grid()
        line[0].set_label('Raw Signal')
        line[1].set_label('Filtered w/ cal file Signal')
        ax[j].legend(loc='upper right')
    ax[0].set(title='Raw vs filtered signal norms after averaging')
    ax[1].set(xlabel='Time (s)')
    del j, i


def plot_ave_vs_raw(data, info):

    norms = find_norm(data)
    ave_norms = ave_norm(data, info)

    fig, ax = plt.subplots(2, 1, figsize=(13, 9))

    fs = 1/1000
    xs = np.arange(0, 1000) * fs

    for j, i in enumerate(['r_mdat', 'r_gdat']):
        line = ax[j].plot(
            xs, norms[i][0:1000], 'k',
            xs, ave_norms[i], 'b',
            linewidth=1,
        )
        ax[j].set(ylabel='Averaged Signal Norm')
        ax[j].grid()
        line[0].set_label('Raw Signal')
        line[1].set_label('Averaged Raw Signal')
        ax[j].legend(loc='upper right')
    ax[0].set(title='Raw vs Averaged Signal Norms')
    ax[1].set(xlabel='Time (s)')
    del j, i