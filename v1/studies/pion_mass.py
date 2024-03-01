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
T = int(folder[idxT+1:idx_])-0      # remove -1
L = int(folder[idxL+1:idx__-2])-0   # remove -1

sum_axis = 2
dim_sum = L if sum_axis == 1 else T
other = T if sum_axis == 1 else L

if __name__=="__main__":
     # find pion mass
    filename = f"{data_folder}PP.csv"
    df = pd.read_csv(filename)
    time = df.iloc[:, 0].values
    number_of_time_values = len(time)
    df = df.drop(df.columns[0], axis=1)
    if df.shape[1] != (T+0) * (L+0): # remove +1
        raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
    reshaped_array = df.values.reshape(number_of_time_values, T+0, L+0) # remove +1
    reshaped_array = reshaped_array[:,0:,0:] # remove line

    #  generate ensemble of means using jackknife
    jackknife_means = reshaped_array.copy()[ntherm:, :, :]
    for i in range(jackknife_means.shape[0]):
        jackknife_means[i,:,:] = (np.sum(jackknife_means[:i,:,:], axis=0) + np.sum(jackknife_means[i+1:,:,:], axis=0))/(jackknife_means.shape[0]-1)
    summed_jackknife_means = np.sum(jackknife_means, axis=sum_axis)
    summed_jackknife_mean_avg = np.mean(summed_jackknife_means, axis=0)
    summed_jackknife_mean_std = np.std(summed_jackknife_means, axis=0)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
    ax1.errorbar(np.arange(dim_sum), summed_jackknife_mean_avg, 3*summed_jackknife_mean_std, fmt='.-', color='black', capsize=2)
    ax1.set_yscale('log')
    ax1.set_title(r'$\langle P(t)  P(0) \rangle$')

    # LOG

    meff = summed_jackknife_means.copy()
    meff[:, :-1] /= summed_jackknife_means[:, 1:]
    meff[:, -1] /= summed_jackknife_means[:, 0]
    meff = np.log(meff)
    meff_avg = np.mean(meff, axis=0)
    meff_std = np.std(meff, axis=0)*np.sqrt((meff.shape[0]-1)/meff.shape[0])
    ax2.errorbar(np.arange(dim_sum)+0.5, other*np.abs(meff_avg), other*3*meff_std, fmt='.-', color='black', capsize=2)
    # ax2.plot(np.arange(0, dim_sum//2), np.ones(dim_sum//2)*4.7/other, '--', color='blue')
    # ax2.plot(np.arange(dim_sum//2, dim_sum), -np.ones(dim_sum-dim_sum//2)*4.7/other, '--', color='blue')
    ax2.plot(np.arange(dim_sum)+0.5, np.ones(dim_sum)*4.7, '--', color='blue')
    ax2.set_title(r'$m_{\pi} \dot L = L \times \left| \log \frac{\langle P(t)  P(0) \rangle}{\langle P(t+a)  P(0) \rangle} \right|$')

    # COSH solution
    
    meff = summed_jackknife_means.copy()
    m0 = 1
    def equation_to_solve(m_eff, n_t, N_T, c_ratio):
        return c_ratio * np.cosh(m_eff * ((n_t + 1)//N_T - N_T/2)) - np.cosh(m_eff * (n_t - N_T/2))
    for sample in range(summed_jackknife_means.shape[0]):
        C_pp = summed_jackknife_means[sample, :]
        for i in range(dim_sum):
            n_t_val = i
            N_T_val = dim_sum
            R_val = C_pp[i]/C_pp[(i+1)//dim_sum]
            m_eff_guess = 2
            m_eff_solution = fsolve(equation_to_solve, m_eff_guess, args=(n_t_val, N_T_val, R_val))
            meff[sample][i] = m_eff_solution
    meff[:,-1] = meff[:, 1]
    meff_avg = np.mean(meff, axis=0)
    meff_std = np.std(meff, axis=0)*np.sqrt((meff.shape[0]-1)/meff.shape[0])
    ax3.errorbar(np.arange(dim_sum)[1:]-0.5, other*meff_avg[1:], other*3*meff_std[1:], fmt='.-', color='black', capsize=2)
    ax3.plot(np.arange(dim_sum), np.ones(dim_sum)*4.7, '--', color='blue')
    ax3.set_title(r'effective mass (cosh solution)')

    print(meff_avg)
    print(meff_std)

    fig.savefig(f'{plot_folder}pion_mass.png')