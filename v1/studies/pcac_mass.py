import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import fsolve
import sys

args = sys.argv

folder = args[1]
data_folder = folder+"observables/"
plot_folder = folder+"plots/"
ntherm = int(args[2]) if len(args) > 2 else 0

idxT = folder.find("T")
idx_ = folder.find("_")
idxL = folder.find("L")
idx__ = folder.find("/")
T = int(folder[idxT+1:idx_])-2      # remove -1
L = int(folder[idxL+1:idx__-2])-2   # remove -1

sum_axis = 2
dim_sum = L if sum_axis == 1 else T
other = T if sum_axis == 1 else L


filename = f"{data_folder}PP.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != (T+2) * (L+2): # remove +1
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T+2, L+2) # remove +1
reshaped_array = reshaped_array[:,1:-1,1:-1] # remove line

#  generate ensemble of PP means using jackknife
PP_jackknife_means = reshaped_array.copy()[ntherm:, :, :]
for i in range(PP_jackknife_means.shape[0]):
    PP_jackknife_means[i,:,:] = (np.sum(PP_jackknife_means[:i,:,:], axis=0) + np.sum(PP_jackknife_means[i+1:,:,:], axis=0))/(PP_jackknife_means.shape[0]-1)
PP_summed_jackknife_means = np.sum(PP_jackknife_means, axis=sum_axis)
PP_summed_jackknife_mean_avg = np.mean(PP_summed_jackknife_means, axis=0)
PP_summed_jackknife_mean_std = np.std(PP_summed_jackknife_means, axis=0)*np.sqrt((PP_jackknife_means.shape[0]-1)/PP_jackknife_means.shape[0])

filename = f"{data_folder}PA1.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != (T+2) * (L+2): # remove +1
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T+2, L+2) # remove +1
reshaped_array = reshaped_array[:,1:-1,1:-1] # remove line

#  generate ensemble of AP means using jackknife
AP_jackknife_means = reshaped_array.copy()[ntherm:, :, :]
for i in range(AP_jackknife_means.shape[0]):
    AP_jackknife_means[i,:,:] = (np.sum(AP_jackknife_means[:i,:,:], axis=0) + np.sum(AP_jackknife_means[i+1:,:,:], axis=0))/(AP_jackknife_means.shape[0]-1)
AP_summed_jackknife_means = np.sum(AP_jackknife_means, axis=sum_axis)
np.random.shuffle(AP_jackknife_means)
AP_summed_jackknife_mean_avg = np.mean(AP_summed_jackknife_means, axis=0)
AP_summed_jackknife_mean_std = np.std(AP_summed_jackknife_means, axis=0)*np.sqrt((AP_jackknife_means.shape[0]-1)/AP_jackknife_means.shape[0])

# derivative of means
d_AP_summed_jackknife_means = np.zeros_like(AP_summed_jackknife_means)
d_AP_summed_jackknife_means[:,1:-1] = 0.5*(AP_summed_jackknife_means[:,2:] - AP_summed_jackknife_means[:,:-2])
d_AP_summed_jackknife_means[:,0] = 0.5*(AP_summed_jackknife_means[:,1] - AP_summed_jackknife_means[:,-1])
d_AP_summed_jackknife_means[:,-1] = 0.5*(AP_summed_jackknife_means[:,0] - AP_summed_jackknife_means[:,-2])
# d_AP_summed_jackknife_means[:,:-1] = AP_summed_jackknife_means[:,1:] - AP_summed_jackknife_means[:,:-1]
# d_AP_summed_jackknife_means[:,-1] = AP_summed_jackknife_means[:,0] - AP_summed_jackknife_means[:,-1]
d_AP_summed_jackknife_mean_avg = np.mean(d_AP_summed_jackknife_means, axis=0)
d_AP_summed_jackknife_mean_std = np.std(d_AP_summed_jackknife_means, axis=0)*np.sqrt((d_AP_summed_jackknife_means.shape[0]-1)/d_AP_summed_jackknife_means.shape[0])

# quotient
quotient_samples = 0.5*np.abs(d_AP_summed_jackknife_means)/PP_summed_jackknife_means
quotient_avg = np.mean(quotient_samples, axis=0)
quotient_std = np.std(quotient_samples, axis=0)*np.sqrt((quotient_samples.shape[0]-1)/quotient_samples.shape[0])



fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
ax1.errorbar(np.arange(dim_sum), AP_summed_jackknife_mean_avg, 3*AP_summed_jackknife_mean_std, fmt='.-', color='black', capsize=2)
ax1.set_yscale('log')
ax1.set_title(r"$\sum_t$ < $P(t, x) \ A_1(0, 0)$ >")
ax2.errorbar(np.arange(dim_sum), (d_AP_summed_jackknife_mean_avg), 3*d_AP_summed_jackknife_mean_std, fmt='.-', color='black', capsize=2)
# ax2.set_yscale('log')
ax2.set_title(r"$\sum_t$ $\partial_t$ < $P(t, x) \ A_1(0, 0)$ >")
ax3.errorbar(np.arange(dim_sum), np.abs(d_AP_summed_jackknife_mean_avg), 3*d_AP_summed_jackknife_mean_std, fmt='.-', color='black', capsize=2)
ax3.set_yscale('log')
ax3.set_title(r"|$\partial_x$ $\sum_t$ < $P(t, x) \ A_1(0, 0)$ >|")
ax4.errorbar(np.arange(dim_sum), quotient_avg, 3*quotient_std, fmt='.-', color='black', capsize=2)
ax4.plot(np.arange(dim_sum), np.ones(dim_sum)*(1.1/dim_sum), '--', color='blue')
ax4.set_title(r" |$\partial_x$ $\sum_t$ < $P(t, x) \ A_1(0, 0)$ >| / (2 $\sum_t$ < $P(t, x) \ P(0, 0)$ >)")
fig.savefig(f'{plot_folder}pcac_mass.png')








# summed_over_x_PP = np.sum(reshaped_array, axis=2)
# average_PP = np.mean(summed_over_x_PP[ntherm:, :], axis=0)

# filename = f"{data_folder}PA0.csv"
# df = pd.read_csv(filename)
# time = df.iloc[:, 0].values
# number_of_time_values = len(time)
# df = df.drop(df.columns[0], axis=1)
# if df.shape[1] != (T+1) * (L+1): # remove +1
#     raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
# reshaped_array = df.values.reshape(number_of_time_values, T+1, L+1) # remove +1
# reshaped_array = reshaped_array[:,1:,1:] # remove line
# summed_over_x = np.sum(reshaped_array, axis=2)
# average_summed_over_x = np.mean(summed_over_x[ntherm:, :], axis=0)
# # print(average_summed_over_x)
# std_summed_over_x = np.std(summed_over_x[ntherm:, :], axis=0)/np.sqrt(summed_over_x.shape[0]-ntherm)
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
# ax1.errorbar(np.arange(T), average_summed_over_x, 3*std_summed_over_x, fmt='.-', color='black', capsize=2)
# ax1.set_yscale('log')
# ax1.set_title(r"< $P(t) \ A_0(0)$ >")
# derivative = summed_over_x.copy()
# derivative[:, 1:-1] = 0.5*(summed_over_x[:, 2:] - summed_over_x[:, :-2])
# derivative[:, 0] = 0.5*(summed_over_x[:, 1] - summed_over_x[:, -1])
# derivative[:, -1] = 0.5*(summed_over_x[:, 0] - summed_over_x[:, -2])
# derivative_avg = np.mean(derivative[ntherm:], axis=0)
# derivative_std = np.std(derivative[ntherm:], axis=0)/np.sqrt(derivative.shape[0]-ntherm)
# ax2.errorbar(np.arange(T), (derivative_avg), 3*derivative_std, fmt='.-', color='black', capsize=2)
# # ax2.set_yscale('log')
# ax2.set_title(r"$\partial_t$ < $P(t) \ A_0(0)$ >")

# ax3.errorbar(np.arange(T), np.abs(derivative_avg), 3*derivative_std, fmt='.-', color='black', capsize=2)
# ax3.set_yscale('log')
# ax3.set_title(r"|$\partial_t$ < $P(t) \ A_0(0)$ >|")

# quotient = np.abs(derivative)/(2*summed_over_x_PP)
# quotient_avg = np.mean(quotient[ntherm:], axis=0)
# quotient_std = np.std(quotient[ntherm:], axis=0)/np.sqrt(quotient.shape[0]-ntherm)
# ax4.errorbar(np.arange(T), quotient_avg, 3*quotient_std, fmt='.-', color='black', capsize=2)
# ax4.plot(np.arange(T), np.ones(T)*(1.0/L), '--', color='blue')
# ax4.set_title(r"|$\partial_t$ < $P(t) \ A_0(0)$ >| / (2 < $P(t) \ P(0)$ >)")

# fig.savefig(f'{plot_folder}pcac_mass.png')