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

def plot_correlator(corr, title, log_plot):
    filename = f"{data_folder}{corr}.csv"
    df = pd.read_csv(filename)
    time = df.iloc[:, 0].values
    number_of_time_values = len(time)
    df = df.drop(df.columns[0], axis=1)
    if df.shape[1] != (T) * (L): # remove +1
        raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
    reshaped_array = df.values.reshape(number_of_time_values, T, L) # remove +1
    # reshaped_array = reshaped_array[ntherm:,1:,1:] # remove line

    # batch batch_size consecutive values
    new_shape = (reshaped_array.shape[0] // batch_size, batch_size) + reshaped_array.shape[1:]
    if reshaped_array.shape[0] % batch_size != 0:
        reshaped_array = reshaped_array[:-(reshaped_array.shape[0] % batch_size)]
    reshaped_array = np.reshape(reshaped_array, new_shape)
    reshaped_array = np.mean(reshaped_array, axis=1)

    # generate ensemble of means using jackknife
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
        heatmap = ax0.imshow(np.log(mean_jackknife_mean.transpose())[1:,1:], cmap='viridis', aspect='equal')
        contours = ax0.contour(np.log(mean_jackknife_mean.transpose())[1:,1:], colors='black', linewidths=0.5)
    else:
        heatmap = ax0.imshow((mean_jackknife_mean.transpose())[1:,1:], cmap='viridis', aspect='equal')
        contours = ax0.contour((mean_jackknife_mean.transpose())[1:,1:], colors='black', linewidths=0.5)
    ax0.clabel(contours, inline=True, fontsize=8)
    ax0.set_title(title)
    cbar = fig.colorbar(heatmap, ax=ax0, fraction=0.05, pad=0.04)
    cbar.set_label("(log scale)" if log_plot else "(linear scale)")
    ax0.set_xlabel('t')
    ax0.set_ylabel('x')
    # if corr == 'V0V0' or corr == 'V1V1':
    #     temp = float(L)/float(T)
    #     # ax1.set_xlim(0, 5)
    #     # ax1.set_ylim(-1e-2,1e-1)
    #     ax1.errorbar(np.arange(L), np.arange(L), 3*np.arange(L), color='black', fmt='.-', capsize=2)
    #     # ax1.plot(np.arange(L), np.ones(L)/(2*np.pi)-(np.pi/3)*np.arange(L)**2*temp/L**2, '--', color='blue')
    #     # ax1.set_yscale('log')
    # else:
    ax1.errorbar(np.arange(L)[1:], (summed_over_0_avg)[1:], 3*summed_over_0_std[1:], color='black', fmt='.-', capsize=2)
    # ax1.errorbar(np.arange(L)+1, summed_over_0_avg[::-1], 3*summed_over_0_std[::-1], color='black', fmt='.:', capsize=2)
    if log_plot:
        ax1.set_yscale('log')
    ax1.set_xlabel('x')
    ax2.errorbar(np.arange(T)[1:], (summed_over_1_avg)[1:], 3*summed_over_1_std[1:], color='black', fmt='.-', capsize=2)
    # ax2.errorbar(np.arange(T)+1, summed_over_1_avg[::-1], 3*summed_over_1_std[::-1], color='black', fmt='.:', capsize=2)
    if log_plot:
        ax2.set_yscale('log')
    ax2.set_xlabel('t')
    
    fig.savefig(f"{plot_folder}{corr}.png")
    plt.close(fig)





    # average_values = np.mean(reshaped_array[ntherm:, :, :], axis=0)

    # fig = plt.figure(figsize=(12, 6))
    # gs = GridSpec(10, 18, fig)
    # ax0 = plt.subplot(gs[:,:8])
    # ax1 = plt.subplot(gs[:4,10:])
    # ax2 = plt.subplot(gs[5:,10:])

    # if log_plot:
    #     heatmap = ax0.imshow(np.log(average_values.transpose()), cmap='viridis', aspect='equal')
    #     contours = ax0.contour(np.log(average_values.transpose()), colors='black', linewidths=0.5)
    # else:
    #     heatmap = ax0.imshow((average_values.transpose()), cmap='viridis', aspect='equal')
    #     contours = ax0.contour((average_values.transpose()), colors='black', linewidths=0.5)
    # ax0.clabel(contours, inline=True, fontsize=8)
    # ax0.set_title(title)
    # cbar = fig.colorbar(heatmap, ax=ax0, fraction=0.05, pad=0.04)
    # cbar.set_label("(log scale)" if log_plot else "(linear scale)")
    # ax0.set_xlabel('t')
    # ax0.set_ylabel('x')

    # summed_values_0 = np.sum(reshaped_array, axis=1)
    # average_summed_values_0 = np.mean(summed_values_0[ntherm:, :], axis=0)
    # std_summed_values_0 = np.std(summed_values_0[ntherm:, :], axis=0)/np.sqrt(summed_values_0.shape[0]-ntherm)
    # ax1.errorbar(np.arange(L), average_summed_values_0, std_summed_values_0, color='black', fmt='.-', capsize=2)
    # if log_plot:
    #     ax1.set_yscale('log')
    # ax1.set_xlabel('x')

    # summed_values_1 = np.sum(reshaped_array, axis=2)
    # average_summed_values_1 = np.mean(summed_values_1[ntherm:, :], axis=0)
    # std_summed_values_1 = np.std(summed_values_1[ntherm:, :], axis=0)/np.sqrt(summed_values_1.shape[0]-ntherm)
    # ax2.errorbar(np.arange(T), average_summed_values_1, std_summed_values_1, color='black', fmt='.-', capsize=2)
    # if log_plot:
    #     ax2.set_yscale('log')
    # ax2.set_xlabel('t')
    
    # fig.savefig(f"{plot_folder}{corr}.png")

if __name__=="__main__":
    plot_correlator("PP", r"$\langle \ P (t, x) \ P (0, 0) \ \rangle$", False)
    plot_correlator("V0V0", r"$\langle \ J_0 (t, x) \ J_0 (0, 0) \ \rangle$", False)
    plot_correlator("V0V1", r"$\langle \ J_0 (t, x) \ J_1 (0, 0) \ \rangle$", False)
    plot_correlator("V1V0", r"$\langle \ J_1 (t, x) \ J_0 (0, 0) \ \rangle$", False)
    plot_correlator("V1V1", r"$\langle \ J_1 (t, x) \ J_1 (0, 0) \ \rangle$", False)
    plot_correlator("A0P", r"$\langle \ A_0 (t, x) \ P (0, 0) \ \rangle$", False)
    plot_correlator("A1P", r"$\langle \ A_1 (t, x) \ P (0, 0) \ \rangle$", False)
    plot_correlator("PA0", r"$\langle \ P (t, x) \ A_0 (0, 0) \ \rangle$", False)
    plot_correlator("PA1", r"$\langle \ P (t, x) \ A_1 (0, 0) \ \rangle$", False)

    # print("re-run and exclude thermalization trajectories in averages with\npython plot_corr.py %d %d ${ N_THERM }" % (T, L))