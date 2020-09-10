
# Paul Yollin


# Attempting to using mne/preprocessing/compute_fine_calibration() to generate a fine calibration file
# to apply to the eddy current recordings.


# Things that need to be changed to be able to replicate on different machine
#       data_path
#       .fif file name
#       write location for generated fine calibration file


#   importing required python packages
import mne
from mne.preprocessing._fine_cal import compute_fine_calibration, write_fine_calibration


#   path to .fif files, update before attempting to run script on a separate machine
#       something else to note, I changed the name of my .fif files for easier code
#       originally they were '5Hz_10v_erm_raw.fif', now they are just '5hz_raw.fif'
data_path = '/Users/pjyol/OneDrive/Documents/I-Labs/eddy_current3/'     # update before running!
fif_ext = '_raw.fif'        # originally '_10v_erm_raw.fif' for the eddy current recordings
write_path = '/Users/pjyol/OneDrive/Documents/I-Labs/'



# 5Hz recording


rec = '5hz'     # Defining variable for the recording


# reading the raw fif file, fix_mag_coil_types specifies T3 magnetometer type
raw = mne.io.read_raw_fif(
    fname = data_path + rec + fif_ext,
    allow_maxshield=True
).fix_mag_coil_types()


#   running compute_fine_calibration()
cal_file = compute_fine_calibration(
    raw=raw,            # raw from above
    n_imbalance=3,
    t_window=10         # making the time window for the segments 10 seconds
)


# making a .dat file from the cal_file variable above
#   note you need to change write path to where you want to save the .dat file
write_fine_calibration(
    fname = write_path + rec + '_cal_file.dat',
    calibration=cal_file[0]         # [0] is there since cal_file is a tuple variable
)


# deleting the raw file and cal_file to conserve memory
#   deleting rec to be used again in other recordings.
del raw, cal_file, rec


# 70hz recording


rec = '70hz'

raw = mne.io.read_raw_fif(
    fname=data_path+rec+fif_ext,
    allow_maxshield=True
).fix_mag_coil_types()

cal_file = compute_fine_calibration(
    raw=raw,
    n_imbalance=3,
    t_window=10
)

write_fine_calibration(
    fname=write_path+rec+'_cal_file.dat',
    calibration=cal_file[0]
)

del raw, cal_file, rec


# 105hz recording


rec = '105hz'

raw = mne.io.read_raw_fif(
    fname=data_path+rec+fif_ext,
    allow_maxshield=True
).fix_mag_coil_types()

cal_file = compute_fine_calibration(
    raw=raw,
    n_imbalance=3,
    t_window=10
)

write_fine_calibration(
    fname=write_path+rec+'_cal_file.dat',
    calibration=cal_file[0]
)

del raw, cal_file, rec


# empty room recording


rec = 'erm'

raw = mne.io.read_raw_fif(
    fname=data_path+rec+fif_ext,
    allow_maxshield=True
).fix_mag_coil_types()

cal_file = compute_fine_calibration(
    raw=raw,
    n_imbalance=3,
    t_window=10
)

write_fine_calibration(
    fname=write_path+rec+'_cal_file.dat',
    calibration=cal_file[0]
)

del raw, cal_file, rec

