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
mean_jackknife_mean = np.mean(jackknife_means, axis=0)
summed_over_0 = np.sum(jackknife_means, axis=1)
summed_over_0_avg = np.mean(summed_over_0, axis=0)
summed_over_0_std = np.std(summed_over_0, axis=0)*np.sqrt((jackknife_means.shape[0]-1)/jackknife_means.shape[0])
summed_over_1 = np.sum(jackknife_means, axis=2)
summed_over_1_avg = np.mean(summed_over_1, axis=0)
summed_over_1_std = np.std(summed_over_1, axis=0)*np.sqrt((jackknife_means.shape[0]-1)/jackknife_means.shape[0])

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 3.3))

m = 0.5*((1/kappa) - 4) # ma, mass in inverse lattice spacing

ax1.plot(m*np.arange(L)[1:L//2], (summed_over_0_avg)[1:L//2]*(np.arange(L)[1:L//2]), color='black', marker='s', label=r'lattice measurement ($L_0 = %d$, $L = %d$)'%(T, L))
ax1.set_xlabel(r"$x \cdot m = (x/a) \cdot (m\cdot a)$")
ax1.set_ylabel(r"$(x/a)\cdot\sum_t \langle J_0(t, x) J_0(0, 0) \rangle$")

ax2.plot(m*np.arange(T)[1:T//2], (summed_over_1_avg)[1:T//2]/m, color='black', marker='s', label=r'lattice measurement ($L_0 = %d$, $L = %d$)'%(T, L))
ax2.set_xlabel(r"$t \cdot m = (t/a) \cdot (m\cdot a)$")
ax2.set_ylabel(r"$\frac{1}{ma} \cdot\sum_x \langle J_0(t, x) J_0(0, 0) \rangle$")

def integrand(k, xm):
    return 1/(np.pi)*np.exp(-2*xm*np.sqrt(1+(k/xm)**2))*(1+(xm**2-k**2)/(xm**2+k**2))
def compute_integral(xm):
    result, error = quad(integrand, 0, np.inf, args=(xm,))
    return 2*result

scale = 1
xm_values = m*np.linspace(0, L//2., 256)
integral_values = np.array([compute_integral(xm) for xm in xm_values])
ax1.plot(xm_values, scale*(integral_values), label=f"continuum infinite volume result")
tm_values = m*np.linspace(0, T//2., 256)
integral_values = np.array([compute_integral(tm) for tm in tm_values])
ax2.plot(tm_values, np.zeros_like(tm_values), label=f"continuum infinite volume result")

ax1.legend()
ax2.legend()

# ax1.set_xlim(-0.1, 4.1)

# print("xm_values = [" + ', '.join([str(x) for x in xm_values]) + ']')
# print("integral_values = [" + ', '.join([str(x) for x in integral_values]) + ']')
# print("xm_lattice = [" + ', '.join([str(x) for x in m*np.arange(L)[1:L//2]]) + ']') 
# print("xG00_lattice = [" + ', '.join([str(x) for x in (summed_over_0_avg)[1:L//2]*m*np.arange(L)[1:L//2]]) + ']') 

fig.tight_layout()
fig.savefig(f"{plot_folder}vector.png")
plt.close(fig)

print("m0*L = ", m*L)


