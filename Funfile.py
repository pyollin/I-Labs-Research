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
    def __init__(self, file_name):
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
        self.data = {'r_mdat': self.raw.get_data(picks='mag'),
                     'r_gdat': self.raw.get_data(picks='grad'),
                     'f_mdat': self.rawSSS.get_data(picks='mag'),
                     'f_gdat': self.rawSSS.get_data(picks='grad'),
                     'fc_mdat': self.rawSSSc.get_data(picks='mag'),
                     'fc_gdat': self.rawSSSc.get_data(picks='grad')}
        # deleting the bad gradiometer channels from the raw gradiometere array
        self.data['r_gdat'] = np.delete(self.data['r_gdat'], [104, 130, 139], axis=0)


def find_norm(data):

    # This function takes in the data dictionary created from a RecClass object and returns a dictionary called
    # norms. These norms are the signal vector norms at a given point in time, so if your recording is 10,000 samples
    # then the norm vectors will be of length 10,000.

    # creating empty dictionary to hold norm vectors
    norms = {}

    # for loop going through and getting the norm vectors for all the data arrays in the RecClass.data object.
    # r_mdat means raw magnetometer data, f_gdat means filtered gradiometer data, and
    # fc_mdat means filtered using fine calibration file magnetometer data
    for i in ['r_mdat', 'r_gdat', 'f_mdat', 'f_gdat', 'fc_mdat', 'fc_gdat']:
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
    for i in ['r_mdat', 'r_gdat', 'f_mdat', 'f_gdat', 'fc_mdat', 'fc_gdat']:
        rms[i] = np.sqrt(np.dot(norms[i], norms[i]))

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
        name = i + 'dat'

        #   calculating raw rms over filtered rms without calibration file
        drop[name] = rms['r_' + i + 'dat'] / rms['f_' + i + 'dat']

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
            xs[start:stop], norms['f_' + name][start:stop], 'g',
            xs[start:stop], norms['fc_' + name][start:stop], 'b',
            linewidth=1,
        )
        ax[j].set(ylabel='Signal Norm')
        ax[j].grid()
        line[0].set_label('Raw Signal')
        line[1].set_label('Filtered w/o cal file Signal')
        line[2].set_label('Filtered w/ cal file Signal')
        ax[j].legend(loc='upper right')
    ax[0].set(title='Raw vs filtered signal vector norms over time')
    ax[1].set(xlabel='Time (s)')
