import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import fsolve
import sys
import matplotlib as mpl
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 12

batch_size = 10

def plot_pion_mass(T, L, data_folder, plot_folder, ntherm, sum_axis, plot=False):

     # find pion mass
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

    #  generate ensemble of means using jackknife
    jackknife_means = reshaped_array.copy()[int(np.ceil(ntherm*1.0/batch_size)):, :, :]
    for i in range(jackknife_means.shape[0]):
        jackknife_means[i,:,:] = (np.sum(jackknife_means[:i,:,:], axis=0) + np.sum(jackknife_means[i+1:,:,:], axis=0))/(jackknife_means.shape[0]-1)
    summed_jackknife_means = np.sum(jackknife_means, axis=sum_axis)
    summed_jackknife_mean_avg = np.mean(summed_jackknife_means, axis=0)
    summed_jackknife_mean_std = np.std(summed_jackknife_means, axis=0)

    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.errorbar(np.arange(dim_sum), summed_jackknife_mean_avg, 3*summed_jackknife_mean_std, fmt='.-', color='black', capsize=2)
    ax1.set_yscale('log')
    ax1.set_title(r'$\langle P(x)  P(0) \rangle$')

    # LOG

    # meff = summed_jackknife_means.copy()
    # meff[:, :-1] /= summed_jackknife_means[:, 1:]
    # meff[:, -1] /= summed_jackknife_means[:, 0]
    # meff = np.log(meff)
    # meff_avg = np.mean(meff, axis=0)
    # meff_std = np.std(meff, axis=0)*np.sqrt((meff.shape[0]-1)/meff.shape[0])

    # ax2.errorbar(np.arange(dim_sum)+0.5, L*np.abs(meff_avg), L*3*meff_std, fmt='.-', color='black', capsize=2)
    # # ax2.plot(np.arange(0, dim_sum//2), np.ones(dim_sum//2)*4.7/other, '--', color='blue')
    # # ax2.plot(np.arange(dim_sum//2, dim_sum), -np.ones(dim_sum-dim_sum//2)*4.7/other, '--', color='blue')
    # ax2.plot(np.arange(dim_sum)+0.5, np.ones(dim_sum)*4.7, '--', color='blue')
    # ax2.set_title(r'$m_{\pi} \dot L = L \times \left| \log \frac{\langle P(t)  P(0) \rangle}{\langle P(t+a)  P(0) \rangle} \right|$')

    # COSH solution
    
    meff = np.zeros((summed_jackknife_means.shape[0], summed_jackknife_means.shape[1]))
    m0 = 1

    for sample in range(summed_jackknife_means.shape[0]):
        for i in range(summed_jackknife_means.shape[1]):
            C_i = summed_jackknife_means[sample, i%dim_sum]
            C_ip1 = summed_jackknife_means[sample, (i+1)%dim_sum]
            def eq(m_):
                return C_i/C_ip1 - np.cosh(m_*((i-(dim_sum)/2)))/np.cosh(m_*((i+1-(dim_sum)/2)))
            meff[sample, i] = np.abs(fsolve(eq, m0))

    meff_avg = np.mean(meff, axis=0)*L
    meff_std = np.std(meff, axis=0)*np.sqrt((meff.shape[0]-1)/meff.shape[0])*L*3

    ax3.errorbar(np.arange(dim_sum), meff_avg, meff_std, fmt='.-', color='black', capsize=2)
    ax3.plot(np.arange(dim_sum), np.ones(dim_sum)*4.7, '--', color='blue')

    plateau_region1 = np.arange(int(np.floor(0.3*dim_sum)), int(np.floor(0.5*dim_sum)+1))
    plateau_region2 = np.arange(int(np.floor(0.5*dim_sum)), int(np.floor(0.7*dim_sum)+1))
    plateau_indices = np.concatenate((plateau_region1, plateau_region2))
    x_plateau = np.arange(dim_sum)[plateau_indices]
    y_plateau = meff_avg[plateau_indices]
    yerr_plateau = meff_std[plateau_indices]
    jackknife_means = np.array([np.mean(np.delete(y_plateau, i)) for i in range(len(y_plateau))])
    jackknife_estimate = np.mean(jackknife_means)
    n = len(y_plateau)
    jackknife_error = np.sqrt((n-1)/n * np.sum((jackknife_means - jackknife_estimate)**2))*3
    
    ax3.plot(plateau_region1, jackknife_estimate*np.ones_like(plateau_region1), color='r', linestyle='--', label=r'$m_{eff} = %.3f \pm %.3f$' % (jackknife_estimate, jackknife_error))
    ax3.plot(plateau_region2, jackknife_estimate*np.ones_like(plateau_region2), color='r', linestyle='--')
    ax3.fill_between(plateau_region1, jackknife_estimate - jackknife_error, jackknife_estimate + jackknife_error, color='red', alpha=0.3)
    ax3.fill_between(plateau_region2, jackknife_estimate - jackknife_error, jackknife_estimate + jackknife_error, color='red', alpha=0.3)
    ax3.set_title(r'effective mass (cosh solution)')
    ax3.legend()

    if plot: fig.savefig(f'{plot_folder}pion_mass.png')
    plt.close(fig)

    return jackknife_estimate, jackknife_error


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

    plot_pion_mass(T, L, data_folder, plot_folder, ntherm, sum_axis, plot=True)

