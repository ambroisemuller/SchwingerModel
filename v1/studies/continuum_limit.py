import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 12
from pcac_mass import plot_pcac_mass
from pion_mass import plot_pion_mass

args = sys.argv

temp = float(args[1])
ntherm = int(args[2]) if len(args) > 2 else 0

# Set the base_path to the current working directory
base_path = os.getcwd()

if temp == 8:
    directories = [

    ]
elif temp == 6:
    directories = [
        
    ]
elif temp == 5:
    directories = [
        'T3_L15_A/k0.2706_b1.758'
    ]
elif temp == 4:
    directories = [
        'T4_L16_A/k0.2677_b2.000', 
        'T5_L20_A/k0.2603_b3.125', 
        'T6_L24_A/k0.2564_b4.500',
        'T9_L36_A/k0.2519_b10.125'
    ]
elif temp == 3:
    directories = [
        'T6_L18_A/k0.2634_b2.531',
        'T8_L24_A/k0.2564_b4.500'
    ]
elif temp == 2.5:
    directories = [
        'T6_L15_A/k0.2706_b1.758'
    ]
elif temp == 2:
    directories = [
        'T10_L20_A/k0.2603_b3.125', 
        'T12_L24_A/k0.2564_b4.500'
        ]
elif temp == 1.5:
    directories = [
        'T12_L18_A/k0.2634_b2.531'
    ]
elif temp == 1:
    directories = [
        'T16_L16_A/k0.2680_b2.000',
        'T18_L18_A/k0.2634_b2.531'
    ]
elif temp == 0.5:
    directories = [
        
    ]

sum_axis_pcac = 1
sum_axis_pion = 1

a_arr = []
plaq_arr, plaq_err = [], []
cond_arr, cond_err = [], []
pcac_arr, pcac_err = [], []
pion_arr, pion_err = [], [] 

for folder in directories:
    print(f"- analyzing {folder}")

    idxT = folder.find("T")
    idx_ = folder.find("_")
    idxL = folder.find("L")
    idx__ = folder.find("/")
    T = int(folder[idxT+1:idx_])
    L = int(folder[idxL+1:idx__-2])

    a_arr.append((1/L)**2)
    
    data_folder = folder+"/observables/"
    plot_folder = folder+"/plots/"
    
    # plaquette
    df_plaq = pd.read_csv(f"{data_folder}plaq.csv")
    data_plaq = {col: df_plaq[col].values for col in df_plaq.columns}
    t = data_plaq['time']
    plaq = data_plaq['plaq']
    avg, std = np.mean(plaq[ntherm:]), np.std(plaq[ntherm:])*3/np.sqrt(len(plaq[ntherm:]))
    plaq_arr.append(avg)
    plaq_err.append(std)

    # condensate
    try:
        df_cond = pd.read_csv(f"{data_folder}condensate.csv")
        data_cond = {col: df_cond[col].values for col in df_cond.columns}
        t = data_cond['time']
        cond = data_cond['condensate']
        avg, std = np.mean(cond[ntherm:]), np.std(cond[ntherm:])*3/np.sqrt(len(cond[ntherm:]))
        cond_arr.append(avg)
        cond_err.append(std)
    except:
        pass

    # pcac mass (ideally constant)
    pcac_mass, pcac_mass_err = plot_pcac_mass(T, L, data_folder, plot_folder, ntherm, sum_axis_pcac, plot=True)
    pcac_arr.append(pcac_mass)
    pcac_err.append(pcac_mass_err)

    # pion mass
    pion_mass, pion_mass_err = plot_pion_mass(T, L, data_folder, plot_folder, ntherm, sum_axis_pion, plot=True)
    pion_arr.append(pion_mass) # *T/L here T is L0, so this is actually m/temp
    pion_err.append(pion_mass_err)



    fig, ((ax_plaq, ax_cond), (ax_pcac, ax_pion)) = plt.subplots(2, 2, figsize=(12, 8))

    ax_plaq.errorbar(a_arr, plaq_arr, plaq_err, fmt='x', capsize=3, color='black')
    ax_plaq.set_xlabel(r"Lattice spacing ($(a/L)^2$)")
    ax_plaq.set_ylabel(r"Plaquette")

    try:
        ax_cond.errorbar(a_arr, cond_arr, cond_err, fmt='x', capsize=3, color='black')
        ax_cond.set_xlabel(r"Lattice spacing ($(a/L)^2$)")
        ax_cond.set_ylabel(r"Chiral condensate $\langle \bar\psi \psi \rangle$")
    except:
        pass

    ax_pcac.errorbar(a_arr, pcac_arr, pcac_err, fmt='x', capsize=3, color='black')
    ax_pcac.plot([min(a_arr), max(a_arr)], [1, 1], '--', color='gray')
    ax_pcac.set_xlabel(r"Lattice spacing $(a/L)^2$")
    ax_pcac.set_ylabel(r"PCAC mass ($m_{PCAC} \cdot L$)")

    ax_pion.errorbar(a_arr, pion_arr, pion_err, fmt='x', capsize=3, color='black')
    ax_pion.set_xlabel(r"Lattice spacing ($(a/L)^2$)")
    ax_pion.set_ylabel(r"Pion mass ($m_\pi L$)")

    fig.tight_layout()
    fig.savefig(f'continuum_limit_temp{temp}.png')
    plt.close(fig)