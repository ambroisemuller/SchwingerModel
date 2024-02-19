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
T = int(folder[idxT+1:idx_])-1      # remove -1
L = int(folder[idxL+1:idx__-2])-1   # remove -1

if __name__=="__main__":
     # find pion mass
    filename = f"{data_folder}PP.csv"
    df = pd.read_csv(filename)
    time = df.iloc[:, 0].values
    number_of_time_values = len(time)
    df = df.drop(df.columns[0], axis=1)
    if df.shape[1] != (T+1) * (L+1): # remove +1
        raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
    reshaped_array = df.values.reshape(number_of_time_values, T+1, L+1) # remove +1 
    reshaped_array = reshaped_array[:,1:,1:] # remove line
    summed_over_x = np.sum(reshaped_array, axis=2)
    average_summed_over_x = np.mean(summed_over_x[ntherm:, :], axis=0)
    # print(average_summed_over_x)
    std_summed_over_x = np.std(summed_over_x[ntherm:, :], axis=0)/np.sqrt(summed_over_x.shape[0]-ntherm)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
    ax1.errorbar(np.arange(T), average_summed_over_x, 3*std_summed_over_x, fmt='.-', color='black', capsize=2)
    ax1.set_yscale('log')
    ax1.set_title(r'$< $P(t)  P(0)$ >$')

    # LOG

    meff = summed_over_x.copy()
    meff[:, :-1] /= summed_over_x[:, 1:]
    meff[:, -1] /= summed_over_x[:, 0]
    meff = np.log(meff)
    meff_avg = np.mean(meff[ntherm:, :], axis=0)
    meff_std = np.std(meff[ntherm:, :], axis=0)/np.sqrt(meff.shape[0]-ntherm)
    ax2.errorbar(np.arange(T)+0.5, meff_avg, 3*meff_std, fmt='.-', color='black', capsize=2)
    ax2.plot(np.arange(0, T//2), np.ones(T//2)*4.7/L, '--', color='blue')
    ax2.plot(np.arange(T//2, T), -np.ones(T-T//2)*4.7/L, '--', color='blue')
    ax2.set_title(r'effective mass (log)')

    # COSH solution
    
    meff = summed_over_x.copy()
    m0 = 1
    def equation_to_solve(m_eff, n_t, N_T, c_ratio):
        return c_ratio * np.cosh(m_eff * ((n_t + 1)//N_T - N_T/2)) - np.cosh(m_eff * (n_t - N_T/2))
    for sample in range(number_of_time_values):
        C_pp = summed_over_x[sample, :]
        for i in range(T):
            # print(sample)
            n_t_val = i
            N_T_val = T
            R_val = C_pp[i]/C_pp[(i+1)//T]
            m_eff_guess = 2
            m_eff_solution = fsolve(equation_to_solve, m_eff_guess, args=(n_t_val, N_T_val, R_val))
            meff[sample][i] = m_eff_solution

    meff_avg = np.mean(meff[ntherm:, :], axis=0)
    meff_std = np.std(meff[ntherm:, :], axis=0)/np.sqrt(meff.shape[0]-ntherm)
    ax3.errorbar(np.arange(T)[1:]-0.5, meff_avg[1:], 3*meff_std[1:], fmt='.-', color='black', capsize=2)
    # ax3.set_ylim(0, 0.5)
    ax3.plot(np.arange(T), np.ones(T)*4.7/L, '--', color='blue')
    ax3.set_title(r'effective mass (cosh solution)')


    fig.savefig(f'{plot_folder}pion_mass.png')