import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys;

args = sys.argv

# SETUP

folder = args[1]
data_folder = folder+"observables/"
plot_folder = folder+"plots/"
ntherm = int(args[2]) if len(args) > 2 else 0

idxT = folder.find("T")
idx_ = folder.find("_")
idxL = folder.find("L")
idx__ = folder.find("/")
T = int(folder[idxT+1:idx_])      # remove -1
L = int(folder[idxL+1:idx__-2])   # remove -1

# V0V0

filename = f"{data_folder}V0V0.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != (T) * (L): # remove +1
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T, L) # remove +1
# reshaped_array = reshaped_array[ntherm:,1:,1:] # remove line
jackknife_means_00 = reshaped_array.copy()[ntherm:, :, :]
for i in range(jackknife_means_00.shape[0]):
    jackknife_means_00[i,:,:] = (np.sum(jackknife_means_00[:i,:,:], axis=0) + np.sum(jackknife_means_00[i+1:,:,:], axis=0))/(jackknife_means_00.shape[0]-1)
mean_jackknife_mean_00 = np.mean(jackknife_means_00, axis=0)

# V0V1

filename = f"{data_folder}V0V1.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != (T) * (L): # remove +1
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T, L) # remove +1
# reshaped_array = reshaped_array[ntherm:,1:,1:] # remove line
jackknife_means_01 = reshaped_array.copy()[ntherm:, :, :]
for i in range(jackknife_means_01.shape[0]):
    jackknife_means_01[i,:,:] = (np.sum(jackknife_means_01[:i,:,:], axis=0) + np.sum(jackknife_means_01[i+1:,:,:], axis=0))/(jackknife_means_01.shape[0]-1)
mean_jackknife_mean_01 = np.mean(jackknife_means_01, axis=0)

# V1V0

filename = f"{data_folder}V1V0.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != (T) * (L): # remove +1
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T, L) # remove +1
# reshaped_array = reshaped_array[ntherm:,1:,1:] # remove line
jackknife_means_10 = reshaped_array.copy()[ntherm:, :, :]
for i in range(jackknife_means_10.shape[0]):
    jackknife_means_10[i,:,:] = (np.sum(jackknife_means_10[:i,:,:], axis=0) + np.sum(jackknife_means_10[i+1:,:,:], axis=0))/(jackknife_means_10.shape[0]-1)
mean_jackknife_mean_10 = np.mean(jackknife_means_10, axis=0)

# V1V1

filename = f"{data_folder}V1V1.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != (T) * (L): # remove +1
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T, L) # remove +1
# reshaped_array = reshaped_array[ntherm:,1:,1:] # remove line
jackknife_means_11 = reshaped_array.copy()[ntherm:, :, :]
for i in range(jackknife_means_11.shape[0]):
    jackknife_means_11[i,:,:] = (np.sum(jackknife_means_11[:i,:,:], axis=0) + np.sum(jackknife_means_11[i+1:,:,:], axis=0))/(jackknife_means_11.shape[0]-1)
mean_jackknife_mean_11 = np.mean(jackknife_means_11, axis=0)

# SUM OF 4 COMPONENTS

jackknife_means_tot = jackknife_means_00 + jackknife_means_11 - jackknife_means_10 - jackknife_means_01
mean_jackknife_mean_tot = np.mean(jackknife_means_tot, axis=0)

# FFT

x = np.arange(L)
t = np.arange(T)
X_mg, T_mg = np.meshgrid(x, t)

q_arr = np.arange(1, 8, 0.2)
Re_HVP_avg = []
Re_HVP_std = []
Im_HVP_avg = []
Im_HVP_std = []

for q in q_arr:
    exp_term_q = np.exp(1j*q*(X_mg+T_mg)/L)
    hvp_samples_q = np.array([(1/(2*q*q*L*T))*np.sum(np.sum(exp_term_q*jackknife_means_tot[i,:,:], axis=1), axis=0) for i in range(jackknife_means_tot.shape[0])])
    exp_term_2q = np.exp(1j*q*(X_mg+T_mg)/L)
    hvp_samples_2q = np.array([(1/(8*q*q*L*T))*np.sum(np.sum(exp_term_2q*jackknife_means_tot[i,:,:], axis=1), axis=0) for i in range(jackknife_means_tot.shape[0])])
    hvp_samples = hvp_samples_2q - hvp_samples_q
    factor = jackknife_means_tot.shape[0]/(jackknife_means_tot.shape[0]+1)
    hvp_re, hvp_im = np.real(hvp_samples), np.imag(hvp_samples)
    hvp_re_avg, hvp_re_std, hvp_im_avg, hvp_im_std = np.mean(hvp_re), np.std(hvp_re)*factor, np.mean(hvp_im), np.std(hvp_im)*factor
    Re_HVP_avg.append(hvp_re_avg)
    Re_HVP_std.append(hvp_re_std)
    Im_HVP_avg.append(hvp_im_avg)
    Im_HVP_std.append(hvp_im_std)

print(Re_HVP_avg)
print(Re_HVP_std)
print(Im_HVP_avg)
print(Im_HVP_std)

# PLOT

fig = plt.figure(figsize=(12, 6))

ax1 = fig.add_subplot(231) # , projection='3d'
# ax1.plot_wireframe(X_mg[1:-1,1:-1], T_mg[1:-1,1:-1], np.log(mean_jackknife_mean_00)[1:-1,1:-1], color='black')
# max_range = max(x[-1]-x[0], t[-1]+t[0])*0.5
# mid_x = (x[-1]+x[0]) * 0.5
# mid_t = (t[-1]+t[0]) * 0.5
# ax1.set_xlim(mid_x - max_range, mid_x + max_range)
# ax1.set_ylim(mid_t - max_range, mid_t + max_range)
heatmap = ax1.imshow(np.log(mean_jackknife_mean_00.transpose()), cmap='viridis', aspect='equal')
contours = ax1.contour(np.log(mean_jackknife_mean_00.transpose()), colors='black', linewidths=0.5)
ax1.clabel(contours, inline=True, fontsize=8)
ax1.set_title(r"$\langle J_0(t, x) J_0(0,0) \rangle$")
cbar = fig.colorbar(heatmap, ax=ax1, fraction=0.05, pad=0.04)
cbar.set_label("(log scale)")
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$x$')

ax2 = fig.add_subplot(232)
heatmap = ax2.imshow(np.log(mean_jackknife_mean_01.transpose()), cmap='viridis', aspect='equal')
contours = ax2.contour(np.log(mean_jackknife_mean_01.transpose()), colors='black', linewidths=0.5)
ax2.clabel(contours, inline=True, fontsize=8)
ax2.set_title(r"$\langle J_0(t, x) J_1(0,0) \rangle$")
cbar = fig.colorbar(heatmap, ax=ax2, fraction=0.05, pad=0.04)
cbar.set_label("(log scale)")
ax2.set_xlabel(r'$t$')
ax2.set_ylabel(r'$x$')

ax4 = fig.add_subplot(234)
heatmap = ax4.imshow(np.log(mean_jackknife_mean_10.transpose()), cmap='viridis', aspect='equal')
contours = ax4.contour(np.log(mean_jackknife_mean_10.transpose()), colors='black', linewidths=0.5)
ax4.clabel(contours, inline=True, fontsize=8)
ax4.set_title(r"$\langle J_1(t, x) J_0(0,0) \rangle$")
cbar = fig.colorbar(heatmap, ax=ax4, fraction=0.05, pad=0.04)
cbar.set_label("(log scale)")
ax4.set_xlabel(r'$t$')
ax4.set_ylabel(r'$x$')

ax5 = fig.add_subplot(235)
heatmap = ax5.imshow(np.log(mean_jackknife_mean_11.transpose()), cmap='viridis', aspect='equal')
contours = ax5.contour(np.log(mean_jackknife_mean_11.transpose()), colors='black', linewidths=0.5)
ax5.clabel(contours, inline=True, fontsize=8)
ax5.set_title(r"$\langle J_1(t, x) J_1(0,0) \rangle$")
cbar = fig.colorbar(heatmap, ax=ax5, fraction=0.05, pad=0.04)
cbar.set_label("(log scale)")
ax5.set_xlabel(r'$t$')
ax5.set_ylabel(r'$x$')

ax3 = fig.add_subplot(233)
heatmap = ax3.imshow(np.log(mean_jackknife_mean_tot.transpose()), cmap='viridis', aspect='equal')
contours = ax3.contour(np.log(mean_jackknife_mean_tot.transpose()), colors='black', linewidths=0.5)
ax3.clabel(contours, inline=True, fontsize=8)
ax3.set_title(r"$\sum_\mu \sum_\nu \delta_{\mu\nu} \langle J_\mu(t, x) J_\nu(0,0) \rangle $")
cbar = fig.colorbar(heatmap, ax=ax3, fraction=0.05, pad=0.04)
cbar.set_label("(log scale)")
ax3.set_xlabel(r'$t$')
ax3.set_ylabel(r'$x$')

ax6 = fig.add_subplot(236)
ax6.errorbar(q_arr, Re_HVP_avg, Re_HVP_std)
ax6.errorbar(q_arr, Im_HVP_avg, Im_HVP_std)
# ax6.set_yscale('log')

fig.tight_layout()    
fig.savefig(f"{plot_folder}HVP.png")

