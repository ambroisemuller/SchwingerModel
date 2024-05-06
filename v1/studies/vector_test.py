import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
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

subitem = folder[idx__+1:]
idxk = subitem.find("k")
idx_ = subitem.find("_")
idxb = subitem.find("b")
idx__ = subitem.find("/")
kappa = float(subitem[idxk+1:idx_])
beta = float(subitem[idxb+1:idx__-2])

batch_size = 1
log_plot = False
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
mean_jackknife_mean = np.mean(jackknife_means, axis=0)
summed_over_0 = np.sum(jackknife_means, axis=1)
summed_over_0_avg = np.mean(summed_over_0, axis=0)
summed_over_0_std = np.std(summed_over_0, axis=0)*np.sqrt((jackknife_means.shape[0]-1)/jackknife_means.shape[0])
summed_over_1 = np.sum(jackknife_means, axis=2)
summed_over_1_avg = np.mean(summed_over_1, axis=0)
summed_over_1_std = np.std(summed_over_1, axis=0)*np.sqrt((jackknife_means.shape[0]-1)/jackknife_means.shape[0])

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 4))

m = 0.5*((1/kappa) - 4)

ax1.scatter(m*np.arange(L)[1:L//2], (summed_over_0_avg)[1:L//2]*m*(np.arange(L)[1:L//2]), color='black', marker='s')
if log_plot:
    ax1.set_yscale('log')
ax1.set_xlabel('x')

ax2.errorbar(np.arange(T)[1:], (summed_over_1_avg)[1:], 3*summed_over_1_std[1:], color='black', fmt='.-', capsize=2)
if log_plot:
    ax2.set_yscale('log')
ax2.set_xlabel('t')

def integrand(k, xm):
    return xm*(1/(2*np.pi))*np.exp(-2*xm*np.sqrt(k**2+1))*(1+(1-k**2)/(1+k**2))
def compute_integral(xm):
    result, error = quad(integrand, 0, np.inf, args=(xm,))
    return 2*result

scale = 1/(4*np.pi)
xm_values = m*np.linspace(0, L//2., 256)
integral_values = np.array([compute_integral(xm) for xm in xm_values])
ax1.plot(xm_values, scale*(integral_values), label=f"m={m}")

ax1.set_xlim(0, 3)

print("xm_values = [" + ', '.join([str(x) for x in xm_values]) + ']')
print("integral_values = [" + ', '.join([str(x) for x in integral_values]) + ']')
print("xm_lattice = [" + ', '.join([str(x) for x in m*np.arange(L)[1:L//2]]) + ']') 
print("xG00_lattice = [" + ', '.join([str(x) for x in (summed_over_0_avg)[1:L//2]*m*np.arange(L)[1:L//2]]) + ']') 

fig.tight_layout()
fig.savefig(f"{plot_folder}vector.png")
plt.close(fig)

print("m0*L = ", m*L)


