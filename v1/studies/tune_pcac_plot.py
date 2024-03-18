import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import fsolve
import os
import subprocess
import sys
import matplotlib as mpl
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 12
from pcac_mass import plot_pcac_mass

kappas = []
masses = []
errors = []

args = sys.argv
t_folder = args[1]

ntherm = int(args[2]) if len(args) > 2 else 0

base_path = os.getcwd()
t_dir_path = os.path.join(base_path, t_folder)
for subitem in os.listdir(t_dir_path):
    if os.path.isdir(os.path.join(t_dir_path, subitem)) and subitem.startswith('k'):

        subitem += "/"
        folder = os.path.join(t_folder, subitem)
        
        data_folder = folder+"observables/"
        plot_folder = folder+"plots/"

        idxT = t_folder.find("T")
        idx_ = t_folder.find("_")
        idxL = t_folder.find("L")
        idx__ = t_folder.find("/")
        T = int(t_folder[idxT+1:idx_])
        L = int(t_folder[idxL+1:idx__-2])

        idxk = subitem.find("k")
        idx_ = subitem.find("_")
        idxb = subitem.find("b")
        idx__ = subitem.find("/")
        kappa = float(subitem[idxk+1:idx_])
        beta = float(subitem[idxb+1:idx__-2])

        sum_axis = 2

        m_pcac, err = plot_pcac_mass(T, L, data_folder, plot_folder, ntherm, sum_axis, True)

        print(f"kappa = {kappa}, beta = {beta} -> m_pcac = {m_pcac} +/- {err}")

        kappas.append(kappa)
        masses.append(m_pcac)
        errors.append(err)

        m, c = np.polyfit(kappas, masses, 1)
        k_range = np.linspace(min(kappas), max(kappas))

        k_critical = (1-c)/m

        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.errorbar(kappas, masses, errors, fmt='d', capsize=3, color='black', label='Extracted masses')
        ax.plot(k_range, m*k_range+c, '--', label='Linear fit')
        ax.plot(k_range, np.ones_like(k_range), '--', color='black', label='Target physical mass')
        ax.plot([k_critical, k_critical], [min(masses)-errors[np.argmin(masses)], 1], '--', color='gray', label=r'$\kappa^* = $%.4f'%(k_critical))
        # ax.annotate(, (k_critical+0.0001, min(masses)))
        ax.set_xlabel(r'$\kappa$')
        ax.set_ylabel(r'$m_{PCAC}$')
        ax.legend()
        fig.tight_layout()
        fig.savefig(t_folder+'mass_tuning.png')
        plt.close(fig)
                