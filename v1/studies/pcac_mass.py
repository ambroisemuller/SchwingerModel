import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import fsolve
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl
# mpl.rcParams['font.family'] = 'STIXGeneral'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex=True)
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.size'] = 12

batch_size = 10

def plot_pcac_mass(T, L, data_folder, plot_folder, ntherm, sum_axis, plot=False):

    dim_sum = L if sum_axis == 1 else T
    other = T if sum_axis == 1 else L

    filename = f"{data_folder}PP.csv"
    df = pd.read_csv(filename)
    time = df.iloc[:, 0].values
    number_of_time_values = len(time)
    df = df.drop(df.columns[0], axis=1)
    if df.shape[1] != T*L:
        raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
    reshaped_array = df.values.reshape(number_of_time_values, T, L)

    # batch batch_size consecutive values
    new_shape = (reshaped_array.shape[0] // batch_size, batch_size) + reshaped_array.shape[1:]
    if reshaped_array.shape[0] % batch_size != 0:
        reshaped_array = reshaped_array[:-(reshaped_array.shape[0] % batch_size)]
    reshaped_array = np.reshape(reshaped_array, new_shape)
    reshaped_array = np.mean(reshaped_array, axis=1)

    #  generate ensemble of PP means using jackknife
    PP_jackknife_means = reshaped_array.copy()[int(np.ceil(ntherm*1.0/batch_size)):, :, :]
    for i in range(PP_jackknife_means.shape[0]):
        PP_jackknife_means[i,:,:] = (np.sum(PP_jackknife_means[:i,:,:], axis=0) + np.sum(PP_jackknife_means[i+1:,:,:], axis=0))/(PP_jackknife_means.shape[0]-1)
    PP_summed_jackknife_means = np.sum(PP_jackknife_means, axis=sum_axis)
    PP_summed_jackknife_mean_avg = np.mean(PP_summed_jackknife_means, axis=0)
    PP_summed_jackknife_mean_std = np.std(PP_summed_jackknife_means, axis=0)*np.sqrt((PP_jackknife_means.shape[0]-1)/PP_jackknife_means.shape[0])

    filename = f"{data_folder}PA{'1' if sum_axis==1 else '0'}.csv"
    df = pd.read_csv(filename)
    time = df.iloc[:, 0].values
    number_of_time_values = len(time)
    df = df.drop(df.columns[0], axis=1)
    if df.shape[1] != T*L:
        raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
    reshaped_array1 = df.values.reshape(number_of_time_values, T, L) 

    filename = f"{data_folder}A{'1' if sum_axis==1 else '0'}P.csv"
    df = pd.read_csv(filename)
    time = df.iloc[:, 0].values
    number_of_time_values = len(time)
    df = df.drop(df.columns[0], axis=1)
    if df.shape[1] != T*L:
        raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
    reshaped_array2 = df.values.reshape(number_of_time_values, T, L) 

    min_dim = min(reshaped_array1.shape[0], reshaped_array2.shape[0])

    # average
    reshaped_array = 0.5*(reshaped_array1[:min_dim, :, :] - reshaped_array2[:min_dim, :, :])

    # batch batch_size consecutive values
    new_shape = (reshaped_array.shape[0] // batch_size, batch_size) + reshaped_array.shape[1:]
    if reshaped_array.shape[0] % batch_size != 0:
        reshaped_array = reshaped_array[:-(reshaped_array.shape[0] % batch_size)]
    reshaped_array = np.reshape(reshaped_array, new_shape)
    reshaped_array = np.mean(reshaped_array, axis=1)

    #  generate ensemble of AP means using jackknife
    AP_jackknife_means = reshaped_array.copy()[int(np.ceil(ntherm*1.0/batch_size)):, :, :]
    for i in range(AP_jackknife_means.shape[0]):
        AP_jackknife_means[i,:,:] = (np.sum(AP_jackknife_means[:i,:,:], axis=0) + np.sum(AP_jackknife_means[i+1:,:,:], axis=0))/(AP_jackknife_means.shape[0]-1)
    AP_summed_jackknife_means = np.sum(AP_jackknife_means, axis=sum_axis)
    AP_summed_jackknife_mean_avg = np.mean(AP_summed_jackknife_means, axis=0)
    AP_summed_jackknife_mean_std = np.std(AP_summed_jackknife_means, axis=0)*np.sqrt((AP_jackknife_means.shape[0]-1)/AP_jackknife_means.shape[0])

    # derivative of means
    d_AP_summed_jackknife_means = np.zeros_like(AP_summed_jackknife_means)
    d_AP_summed_jackknife_means[:,1:-1] = 0.5*(AP_summed_jackknife_means[:,2:] - AP_summed_jackknife_means[:,:-2])
    d_AP_summed_jackknife_means[:,0] = 0.5*(AP_summed_jackknife_means[:,1] - AP_summed_jackknife_means[:,-1])
    d_AP_summed_jackknife_means[:,-1] = 0.5*(AP_summed_jackknife_means[:,0] - AP_summed_jackknife_means[:,-2])
    # d_AP_summed_jackknife_means = np.abs(d_AP_summed_jackknife_means)
    d_AP_summed_jackknife_mean_avg = np.mean(d_AP_summed_jackknife_means, axis=0)
    d_AP_summed_jackknife_mean_std = np.std(d_AP_summed_jackknife_means, axis=0)*np.sqrt((d_AP_summed_jackknife_means.shape[0]-1)/d_AP_summed_jackknife_means.shape[0])

    n_samples = min(d_AP_summed_jackknife_means.shape[0], PP_summed_jackknife_means.shape[0])
    n_sites = d_AP_summed_jackknife_means.shape[1]

    # quotient
    np.random.shuffle(d_AP_summed_jackknife_means)
    min_dim = min(d_AP_summed_jackknife_means.shape[0], PP_summed_jackknife_means.shape[0])
    quotient_samples = 0.5*(d_AP_summed_jackknife_means[:min_dim])/PP_summed_jackknife_means[:min_dim]
    quotient_avg = np.mean(quotient_samples, axis=0)*L
    quotient_std = np.std(quotient_samples, axis=0)*np.sqrt((quotient_samples.shape[0]-1)/quotient_samples.shape[0])*L*3

    plateau_indices = np.arange(int(np.ceil(0.3*dim_sum)), int(np.floor(0.7*dim_sum)+1))
    x_plateau = np.arange(dim_sum)[plateau_indices]
    y_plateau = quotient_avg[plateau_indices]
    yerr_plateau = quotient_std[plateau_indices]
    # jackknife_means = np.array([np.mean(np.delete(y_plateau, i)) for i in range(len(y_plateau))])
    # jackknife_estimate = np.mean(jackknife_means)
    # n = len(y_plateau)
    # jackknife_error = np.sqrt((n-1)/n * np.sum((jackknife_means - jackknife_estimate)**2))*3
    weights = 1 / yerr_plateau**2
    weighted_mean = np.sum(weights * y_plateau) / np.sum(weights)
    error_weighted_mean = np.sqrt(1 / np.sum(weights))*3

    if np.isnan(weighted_mean) or np.isnan(error_weighted_mean):
        print("NaN detected")
        weighted_mean = np.mean(y_plateau)
        error_weighted_mean = np.std(y_plateau)*3/np.sqrt(len(plateau_indices))

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
    ax1.errorbar(np.arange(dim_sum)[1:], AP_summed_jackknife_mean_avg[1:], 3*AP_summed_jackknife_mean_std[1:], fmt='.-', color='black', capsize=2)
    ax1.set_title(r"$\sum_x \langle P(t, x) \ A_0(0, 0) - A_0(t, x) \ P(0, 0) \rangle$")
    ax2.errorbar(np.arange(dim_sum)[1:], PP_summed_jackknife_mean_avg[1:], 3*PP_summed_jackknife_mean_std[1:], fmt='.-', color='black', capsize=2)
    ax2.set_yscale('log')
    ax2.set_title(r"$\sum_x \langle P(t, x) \ P(0, 0) \rangle $")
    ax3.errorbar(np.arange(dim_sum)[2:-1], (d_AP_summed_jackknife_mean_avg)[2:-1], 3*d_AP_summed_jackknife_mean_std[2:-1], fmt='.-', color='black', capsize=2)
    ax3.set_yscale('log')
    ax3.set_title(r" $\partial_t \sum_x \langle P(t, x) \ A_0(0, 0) - A_0(t, x) \ P(0, 0) \rangle $")
    ax4.errorbar(np.arange(dim_sum)[2:-1], quotient_avg[2:-1], quotient_std[2:-1], fmt='.-', color='black', capsize=2)
    ax4.plot(np.arange(dim_sum)[2:-1], np.ones(dim_sum)[2:-1]*(1.0), '--', color='blue')
    ax4.plot(np.arange(dim_sum)[2:-1], np.zeros(dim_sum)[2:-1]*(1.0), '--', color='gray')
    ax4.plot(plateau_indices, weighted_mean*np.ones_like(plateau_indices), color='r', linestyle='--', label=r'$m_{PCAC} = %.3f \pm %.3f$' % (weighted_mean, error_weighted_mean))
    ax4.fill_between(plateau_indices, weighted_mean - error_weighted_mean, weighted_mean + error_weighted_mean, color='red', alpha=0.3)
    ax4.set_title(r" $m_{PCAC} \cdot L = L \times \frac{1}{2} \frac{\partial_t \sum_x \langle P(t, x) \ A_0(0, 0) - A_0(t, x) \ P(0, 0) \rangle}{\sum_x \langle P(t, x) \ P(0, 0) \rangle}$")
    ax4.legend()
    fig.tight_layout()
    if plot: fig.savefig(f'{plot_folder}pcac_mass.png')
    plt.close(fig)

    return weighted_mean, error_weighted_mean

if __name__=="__main__":
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

    sum_axis = 1

    plot_pcac_mass(T, L, data_folder, plot_folder, ntherm, sum_axis, plot=True)