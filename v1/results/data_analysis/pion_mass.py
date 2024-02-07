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

if __name__=="__main__":
     # find pion mass
    filename = f"../observables/PP.csv"
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
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.errorbar(np.arange(T), average_summed_over_x, 3*std_summed_over_x, fmt='.-', color='black', capsize=2)
    ax1.set_yscale('log')
    meff = summed_over_x.copy()
    # meff[:, :-1] /= summed_over_x[:, 1:]
    # meff[:, -1] /= summed_over_x[:, 0]
    # meff = np.log(meff)

    m0 = 0.7
    def equation_to_solve(m_eff, n_t, N_T, c_ratio):
        return c_ratio * np.cosh(m_eff * ((n_t + 1)//N_T - N_T/2)) - np.cosh(m_eff * (n_t - N_T/2))
    for sample in range(number_of_time_values):
        C_pp = summed_over_x[sample, :]
        for i in range(T):
            print(sample)
            n_t_val = i
            N_T_val = T
            R_val = C_pp[i]/C_pp[(i+1)//T]
            m_eff_guess = 2
            m_eff_solution = fsolve(equation_to_solve, m_eff_guess, args=(n_t_val, N_T_val, R_val))
            meff[sample][i] = m_eff_solution

    meff_avg = np.mean(meff[ntherm:, :], axis=0)
    meff_std = np.std(meff[ntherm:, :], axis=0)/np.sqrt(meff.shape[0]-ntherm)
    ax2.errorbar(np.arange(T)+0.5, meff_avg, 3*meff_std, fmt='.-', color='black', capsize=2)
    fig.savefig('../plots/pion_mass.png')