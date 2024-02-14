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
summed_over_x_PP = np.sum(reshaped_array, axis=2)
average_PP = np.mean(summed_over_x_PP[ntherm:, :], axis=0)


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
std_summed_over_x = np.std(summed_over_x[ntherm:, :], axis=0)/np.sqrt(summed_over_x.shape[0]-ntherm)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
ax1.errorbar(np.arange(T), average_summed_over_x, 3*std_summed_over_x, fmt='.-', color='black', capsize=2)
ax1.set_yscale('log')
ax1.set_title(r"< $P(t) \ A_0(0)$ >")
derivative = summed_over_x.copy()
derivative[:, 1:-1] = 0.5*(summed_over_x[:, 2:] - summed_over_x[:, :-2])
derivative[:, 0] = 0.5*(summed_over_x[:, 1] - summed_over_x[:, -1])
derivative[:, -1] = 0.5*(summed_over_x[:, 0] - summed_over_x[:, -2])
derivative_avg = np.mean(derivative[ntherm:], axis=0)
derivative_std = np.std(derivative[ntherm:], axis=0)/np.sqrt(derivative.shape[0]-ntherm)
ax2.errorbar(np.arange(T), (derivative_avg), 3*derivative_std, fmt='.-', color='black', capsize=2)
# ax2.set_yscale('log')
ax2.set_title(r"$\partial_t$ < $P(t) \ A_0(0)$ >")

ax3.errorbar(np.arange(T), np.abs(derivative_avg), 3*derivative_std, fmt='.-', color='black', capsize=2)
ax3.set_yscale('log')
ax3.set_title(r"|$\partial_t$ < $P(t) \ A_0(0)$ >|")

quotient = np.abs(derivative)/(2*summed_over_x_PP)
quotient_avg = np.mean(quotient[ntherm:], axis=0)
quotient_std = np.std(quotient[ntherm:], axis=0)/np.sqrt(quotient.shape[0]-ntherm)
ax4.errorbar(np.arange(T), quotient_avg, 3*quotient_std, fmt='.-', color='black', capsize=2)
ax4.plot(np.arange(T), np.ones(T)*(0.995/32), '--', color='blue')
ax4.set_title(r"|$\partial_t$ < $P(t) \ A_0(0)$ >| / (2 < $P(t) \ P(0)$ >)")

fig.savefig('../plots/pcac_mass.png')