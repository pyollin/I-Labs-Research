#%%
# Paul Yollin
# 07/09/2020

# creating a false calibration file to then test reducation factors


import Funfile as Fun
import numpy as np


params = {
    'f_ext': 'Hz_10v_erm_raw.fif',
    '5': {'start': 80, 'stop': 120},
    '20': {'start': 0, 'stop': 40},
    '70': {'start': 70, 'stop': 110},
    '105': {'start': 65, 'stop': 105}
}

rec= '5'
file_name = rec + params['f_ext']

diag_cal_mat = np.diag(np.append(np.linspace(-0.5, 0.5, 306), np.ones(27)))

dat = Fun.RecClass_Order_Cal(file_name=file_name,
                             start=params[rec]['start'],
                             stop=params[rec]['stop'],
                             internal=2,
                             external=2,
                             diag_cal_mat=diag_cal_mat)

drop = Fun.find_drop_order(dat.data)