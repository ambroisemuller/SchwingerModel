import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys;
import matplotlib as mpl
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 12

args = sys.argv

folder = args[1]
data_folder = folder+"observables/"
plot_folder = folder+"plots/"
ntherm = int(args[2]) if len(args) > 2 else 0

idxT = folder.find("T")
idx_ = folder.find("_")
idxL = folder.find("L")
idx__ = folder.find("/")
T = int(folder[idxT+1:idx_])   
L = int(folder[idxL+1:idx__-2]) 

batch_size = 1
log_plot = False

filename = f"{data_folder}V1V1.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != (T) * (L): 
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T, L) 

new_shape = (reshaped_array.shape[0] // batch_size, batch_size) + reshaped_array.shape[1:]
if reshaped_array.shape[0] % batch_size != 0:
    reshaped_array = reshaped_array[:-(reshaped_array.shape[0] % batch_size)]
reshaped_array = np.reshape(reshaped_array, new_shape)
reshaped_array = np.mean(reshaped_array, axis=1)

jackknife_means = reshaped_array.copy()[int(np.ceil(ntherm*1.0/batch_size)):, :, :]
for i in range(jackknife_means.shape[0]):
    jackknife_means[i,:,:] = (np.sum(jackknife_means[:i,:,:], axis=0) + np.sum(jackknife_means[i+1:,:,:], axis=0))/(jackknife_means.shape[0]-1)
mean_jackknife_mean = np.mean(jackknife_means, axis=0)
summed_over_0 = np.sum(jackknife_means, axis=1)
summed_over_0_avg = np.mean(summed_over_0, axis=0)
summed_over_0_std = np.std(summed_over_0, axis=0)*np.sqrt((jackknife_means.shape[0]-1)/jackknife_means.shape[0])
summed_over_1 = np.sum(jackknife_means, axis=2)
summed_over_1_avg = np.mean(summed_over_1, axis=0)
summed_over_1_std = np.std(summed_over_1, axis=0)*np.sqrt((jackknife_means.shape[0]-1)/jackknife_means.shape[0])

fig = plt.figure(figsize=(12, 6))
gs = GridSpec(10, 18, fig)
ax0 = plt.subplot(gs[:,:8])
ax1 = plt.subplot(gs[:4,10:])
ax2 = plt.subplot(gs[5:,10:])

if log_plot:
    heatmap = ax0.imshow(np.log(mean_jackknife_mean.transpose()), cmap='viridis', aspect='equal')
    contours = ax0.contour(np.log(mean_jackknife_mean.transpose()), colors='black', linewidths=0.5)
else:
    heatmap = ax0.imshow((mean_jackknife_mean.transpose())[1:,1:], cmap='viridis', aspect='equal')
    contours = ax0.contour((mean_jackknife_mean.transpose())[1:,1:], colors='black', linewidths=0.5)
ax0.clabel(contours, inline=True, fontsize=8)
# ax0.set_title()
cbar = fig.colorbar(heatmap, ax=ax0, fraction=0.05, pad=0.04)
cbar.set_label("(log scale)" if log_plot else "(linear scale)")
ax0.set_xlabel('t')
ax0.set_ylabel('x')

x_over_T = np.arange(L)/T
G_tilde = summed_over_0_avg*T
G_tilde_err = summed_over_0_std*T

ax1.plot(x_over_T, 4*(1/(np.sinh(2*np.pi*x_over_T)) + 1/(np.sinh(2*np.pi*(L/T - x_over_T)))), color='blue', linestyle='-')
ax1.errorbar(x_over_T, G_tilde, 3*G_tilde_err, color='black', fmt='.', capsize=2)
# ax1.set_yscale('log')
ax1.set_xlabel('x')

ax2.errorbar(np.arange(T)[1:], (summed_over_1_avg)[1:], 3*summed_over_1_std[1:], color='black', fmt='.-', capsize=2)
if log_plot:
    ax2.set_yscale('log')
ax2.set_xlabel('t')

# ax2.set_ylim(-0.5, 0.5)

fig.savefig(f"{plot_folder}V1V1_test.png")
plt.close(fig)