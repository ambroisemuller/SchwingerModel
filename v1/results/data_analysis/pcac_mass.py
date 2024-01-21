import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import fsolve
import sys

args = sys.argv

T = int(args[1])
L = int(args[2])
ntherm = int(args[3]) if len(args) > 3 else 0

filename = f"../observables/PP.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != T * L:
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T, L)
summed_over_x = np.sum(reshaped_array, axis=2)
average_PP = np.mean(summed_over_x[ntherm:, :], axis=0)


filename = f"../observables/PA0.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != T * L:
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T, L)
summed_over_x = np.sum(reshaped_array, axis=2)
average_summed_over_x = np.mean(summed_over_x[ntherm:, :], axis=0)
# print(average_summed_over_x)
std_summed_over_x = np.std(summed_over_x[ntherm:, :], axis=0)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
ax1.errorbar(np.arange(T), average_summed_over_x, std_summed_over_x, fmt='.-', color='black', capsize=2)
ax1.set_yscale('log')
ax1.set_title(r"< $P(t) \ A_0(0)$ >")
derivative = summed_over_x.copy()
derivative[:, 1:-1] = 0.5*(summed_over_x[:, 2:] - summed_over_x[:, :-2])
derivative[:, 0] = 0.5*(summed_over_x[:, 1] - summed_over_x[:, -1])
derivative[:, -1] = 0.5*(summed_over_x[:, 0] - summed_over_x[:, -2])
derivative_avg = np.mean(derivative[ntherm:], axis=0)
derivative_std = np.std(derivative[ntherm:], axis=0)
ax2.errorbar(np.arange(T), derivative_avg, derivative_std, fmt='.-', color='black', capsize=2)
# ax2.set_yscale('log')
ax2.set_title(r"$\partial_t$ < $P(t) \ A_0(0)$ >")

ax3.plot(derivative_avg + derivative_avg[::-1], '.-', color='black')

quotient = (derivative_avg + derivative_avg[::-1])/average_PP
ax4.plot(quotient, '.-', color='black')

fig.savefig('../plots/pcac_mass.png')