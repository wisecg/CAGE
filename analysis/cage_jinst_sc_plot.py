"""
Code originally from function "plot_run_stability" in sql_to_pandas.py
Extracted to this script for .pkl file creation and styling changes if reviewers ask.
C. Wiseman
19 Sep 2025
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime as dt
from datetime import timedelta
from matplotlib import gridspec
import pickle

# From Grace running her 1460 keV stability code!
run = 66
df_file = f'./data/run{run}_SCDB.h5'
df_file_1460 = './data/run66_1460.h5'

with pd.HDFStore(df_file_1460, 'r') as f:
    print("Keys:", f.keys())

baseline = pd.read_hdf(df_file, key='/cage_baseline')
tophat_temp = pd.read_hdf(df_file, key='/cage_topHat_temp')
cp_temp = pd.read_hdf(df_file, key='/cage_coldPlate_temp')
pressure = pd.read_hdf(df_file, key='/cage_pressure')

df_1460 = pd.read_hdf(df_file_1460, key = '/data')

time_1460 = df_1460['time']
val_1460 = df_1460['val']

# normalize the uncalibrated 1460 keV value by the max value in this time period
# to make deviations in plot more understandable to non-familar persons
max_1460 = np.amax(val_1460)
norm_1460 = np.divide(val_1460, max_1460)

# timestamps in the dataframe from Grace were off by 7 hours. Correct that here
times_new = []
for time in time_1460:
    new_time1 = dt.strptime(time, '%Y-%m-%d %H:%M:%S')
    new_time = new_time1 + timedelta(hours=7)
    times_new.append(new_time)

# Make figure with 5 panels
fig = plt.figure(figsize=(8, 8))
fig.suptitle(f'System Stability for run {run}')
gs = gridspec.GridSpec(5, 1) #, height_ratios=[2, 1]



# the first subplot
ax0 = plt.subplot(gs[0])
line0, = ax0.plot(pressure["timestamp"],
                pressure["cage_pressure"], '-b')
ax0.set_ylabel('Pressure \n(hPa)', fontsize=10)
ax0.set_ylim(4.9e-8, 5.4e-8)
ax0.xaxis.set_visible(False)


# the second subplot
# shared axis X
ax1 = plt.subplot(gs[1], sharex = ax0)
line1, = ax1.plot(cp_temp["timestamp"], cp_temp["cage_coldPlate_temp"] + 273, "-r")
plt.setp(ax1.get_xticklabels(), visible=False)
# remove last tick label for the subplot
yticks = ax1.yaxis.get_major_ticks()
# yticks[-1].label1.set_visible(False)
ax1.set_ylabel('Cold Plate \nTemperature \n(K)', fontsize=10)
ax1.set_ylim(89.8, 90.4)
ax1.xaxis.set_visible(False)



# the third subplot
# shared axis X
ax2 = plt.subplot(gs[2], sharex = ax0)
line1, = ax2.plot(tophat_temp["timestamp"], tophat_temp["cage_topHat_temp"] + 273, "-g")
plt.setp(ax2.get_xticklabels(), visible=False)
# remove last tick label for the subplot
yticks2 = ax2.yaxis.get_major_ticks()
# yticks2[-1].label1.set_visible(False)
ax2.set_ylabel('IR shield \nTemperature \n(K)', fontsize=10)
ax2.set_ylim(96.9, 97.7)
ax2.xaxis.set_visible(False)

# the fourth subplot
# shared axis X
ax3 = plt.subplot(gs[3], sharex = ax0)
line1, = ax3.plot(baseline["timestamp"], baseline["cage_baseline"], "-c")
plt.setp(ax3.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks3 = ax3.yaxis.get_major_ticks()
# yticks3[0].label1.set_visible(False)
# yticks3[-1].label1.set_visible(False)
ax3.set_ylabel('Baseline Value \n(ADC)', fontsize=10)
ax3.set_ylim(-0.685, -0.665)
ax3.xaxis.set_visible(False)


# the fifth subplot
# shared axis X
ax4 = plt.subplot(gs[4], sharex = ax0)
line1, = ax4.plot(times_new, norm_1460, "-b")
# remove last tick label
yticks4 = ax4.yaxis.get_major_ticks()
print(yticks4[-1])
# yticks4[-1].set_visible(False)
ax4.set_ylabel('1460 keV position \n(uncal; norm)', fontsize=10)
ax4.set_ylim(0.997, 1.0005)
ax4.set_xlabel('Timestamp (UTC) \n[mm-dd hr]')

plt.tight_layout()

plt.subplots_adjust(hspace=.025)
plt.savefig('./plots/scm_test.png', dpi=200)
# plt.savefig('./plots/scm_test.pdf', dpi=200)

# jason requested a pickle file
# with open('./plots/scm_test.pkl', 'wb') as f:
#     pickle.dump(fig, f)

# plt.show()