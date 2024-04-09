import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys;
import matplotlib as mpl
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 12

if __name__=="__main__":

    args = sys.argv

    folder = args[1]
    data_folder = folder+"observables/"
    plot_folder = folder+"plots/"
    ntherm = int(args[2]) if len(args) > 2 else 0
    
    df_dH = pd.read_csv(f"{data_folder}dH.csv")
    data_dH = {col: df_dH[col].values for col in df_dH.columns}
    df_acc = pd.read_csv(f"{data_folder}acc.csv")
    data_acc = {col: df_acc[col].values for col in df_acc.columns}
    df_plaq = pd.read_csv(f"{data_folder}plaq.csv")
    data_plaq = {col: df_plaq[col].values for col in df_plaq.columns}
    df_qtop = pd.read_csv(f"{data_folder}qtop.csv")
    data_qtop = {col: df_qtop[col].values for col in df_qtop.columns}
    try:
        df_condensate = pd.read_csv(f"{data_folder}condensate.csv")
        data_condensate = {col: df_condensate[col].values for col in df_condensate.columns}
    except:
        pass

    t = data_dH['time']
    min_length = t.shape[0]
    dH = data_dH['dH']
    min_length = min(min_length, dH.shape[0])
    t = data_acc['time']
    min_length = min(min_length, t.shape[0])
    acc = data_acc['acc']
    min_length = min(min_length, acc.shape[0])
    t = data_plaq['time']
    min_length = min(min_length, t.shape[0])
    plaq = data_plaq['plaq']
    min_length = min(min_length, plaq.shape[0])
    t = data_qtop['time']
    min_length = min(min_length, t.shape[0])
    qtop = data_qtop['qtop']
    min_length = min(min_length, qtop.shape[0])
    try:
        condensate = data_condensate['condensate']
        min_length = min(min_length, condensate.shape[0])
    except:
        pass
    t, dH, acc, plaq, qtop, condensate = t[:min_length], dH[:min_length], acc[:min_length], plaq[:min_length], qtop[:min_length], condensate[:min_length]

    fig, ((ax_dH, ax_acc), (ax_plaq, ax_qtop)) = plt.subplots(2, 2, figsize=(12, 6))
    
    ax_dH.plot(t, dH, '.-', label='HMC measurements')
    ax_dH.plot(t[ntherm:], [np.mean(dH[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(dH[ntherm:]), np.std(dH[ntherm:])))
    ax_dH.set_xlabel(r"Simulation time $\tau$")
    ax_dH.set_ylabel(r"Energy change $dH$")
    ax_dH.legend()

    ax_acc.plot(t, acc, '.-', label='HMC measurements')
    ax_acc.plot(t[ntherm:], [np.mean(acc[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(acc[ntherm:]), np.std(acc[ntherm:])))
    ax_acc.set_xlabel(r"Simulation time $\tau$")
    ax_acc.set_ylabel(r"Acceptance")
    ax_acc.legend()

    ax_plaq.plot(t, plaq, '.-', label='HMC measurements')
    ax_plaq.plot(t[ntherm:], [np.mean(plaq[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(plaq[ntherm:]), np.std(plaq[ntherm:])))
    ax_plaq.set_xlabel(r"Simulation time $\tau$")
    ax_plaq.set_ylabel(r"Plaquette")
    ax_plaq.legend()

    ax_qtop.plot(t, qtop, '.-', label='HMC measurements')
    ax_qtop.plot(t[ntherm:], [np.mean(qtop[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(qtop[ntherm:]), np.std(qtop[ntherm:])))
    ax_qtop.set_xlabel(r"Simulation time $\tau$")
    ax_qtop.set_ylabel(r"Topological charge")
    ax_qtop.legend()

    fig.tight_layout()

    plot_file = f"{plot_folder}general.png"
    fig.savefig(plot_file)
    # print("General observable plots saved to " + plot_file)
    plt.close(fig)

    try:
        fig_cond, ax_cond = plt.subplots(1, 1, figsize=(12, 3))
        ax_cond.plot(t, condensate, '.-', label='Lattice measurements')
        ax_cond.plot(t[ntherm:], [np.mean(condensate[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(condensate[ntherm:]), np.std(condensate[ntherm:])))
        ax_cond.set_xlabel(r"Simulation time $\tau$")
        ax_cond.set_ylabel(r"Chiral fermion condensate")
        ax_cond.legend()
        fig_cond.tight_layout()
        fig_cond.savefig(f"{plot_folder}condensate.png")
        plt.close(fig_cond)
    except:
        pass
    