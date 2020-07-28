#%%
# Paul Yollin
# started:      07/09/2020
# last update:  07/27/2020

# creating a false calibration file to then test reduction factors


import Funfile as Fun
import numpy as np
import matplotlib.pyplot as plt


params = {
    'f_ext': 'Hz_10v_erm_raw.fif',
    '5': {'start': 80, 'stop': 120},
    '20': {'start': 0, 'stop': 40},
    '70': {'start': 70, 'stop': 110},
    '105': {'start': 65, 'stop': 105}
}

xs = np.arange(1, 8)

rec = '5'
file_name = rec + params['f_ext']

#   Creating vector w/ first 306 elements randomly distributed (1,2], last 27 = 1
diagonal = np.append(1+np.random.rand(306), np.ones(27))

#   Creating a dictionary holding calibration error matrix and an identity matrix
Mats = {
    'Identity': np.diag(np.ones(333)),
    'Cal_error': np.diag(diagonal)
}

#   creating figure and axes to plot vector names 'diagonal'
fig, ax = plt.subplots(1, 1, figsize=(13,7))

#   plotting vector, adding gridlines, setting titles and axes labels
ax.plot(diagonal)
ax.grid()
ax.set(title='Diagonal values of generated calibration errors matrix',
       xlabel='Diagonal entry',
       ylabel='Value')

#%%

#   Creating dictionary of empty dictionaries to hold calculated reduction factors
drop = {
    'Identity' : {},
    'Cal_error' : {}
}

#   for loop going through and calculating reduction factors
for e in ['Identity', 'Cal_error']:
    for i in xs:
        dat = Fun.RecClass_Order_Cal(file_name=file_name,
                                     start=params[rec]['start'],
                                     stop=params[rec]['stop'],
                                     internal=2,
                                     external=i,
                                     diag_cal_mat=Mats[e])
        drop[str(e)][str(i)] = Fun.find_drop_order(dat.data)
del i, e

#   using plot_drop_order function found in FunFile to plot reduction factors
Fun.plot_drop_order(drop)
