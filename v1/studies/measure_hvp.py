import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib.gridspec import GridSpec
import sys
from matplotlib import rc
import matplotlib as mpl
# mpl.rcParams['font.family'] = 'STIXGeneral'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
mpl.rcParams['mathtext.fontset'] = 'cm'
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

subitem = folder[idx__+1:]
idxk = subitem.find("k")
idx_ = subitem.find("_")
idxb = subitem.find("b")
idx__ = subitem.find("/")
kappa = float(subitem[idxk+1:idx_])
beta = float(subitem[idxb+1:idx__-2])

batch_size = 1
corr = 'V0V0'

filename = f"{data_folder}{corr}.csv"
df = pd.read_csv(filename)
time = df.iloc[:, 0].values
number_of_time_values = len(time)
df = df.drop(df.columns[0], axis=1)
if df.shape[1] != (T) * (L): 
    raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
reshaped_array = df.values.reshape(number_of_time_values, T, L)

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
mean_jackknife_mean = np.mean(jackknife_means, axis=0)[1:,1:]

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 6))

heatmap = ax1.imshow(np.log(mean_jackknife_mean.transpose()), cmap='viridis', aspect='equal')
contours = ax1.contour(np.log(mean_jackknife_mean.transpose()), colors='black', linewidths=0.5)
ax1.clabel(contours, inline=True, fontsize=8)
cbar = fig.colorbar(heatmap, ax=ax1, fraction=0.05, pad=0.04)
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$x$')

q_list = np.array([4*2*np.pi/L])
print(q_list)

xx, yy = np.meshgrid(np.arange(1, T), np.arange(1, L)) # maybe swap
print(mean_jackknife_mean.shape)
print(xx.shape)

PI = np.array([(1/q)*np.sum(mean_jackknife_mean*np.exp(1j*(q*xx+q*yy))) for q in q_list])
ax2.plot(q_list, np.real(PI))
print(PI)

fig.tight_layout()
fig.savefig(f"{plot_folder}hvp.png")
plt.close(fig)

