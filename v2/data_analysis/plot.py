import pandas as pd
import matplotlib.pyplot as plt

if __name__=="__main__":

    df_dH = pd.read_csv("../results/observables/dH.csv")
    data_dH = {col: df_dH[col].values for col in df_dH.columns}
    df_acc = pd.read_csv("../results/observables/acc.csv")
    data_acc = {col: df_acc[col].values for col in df_acc.columns}
    df_plaq = pd.read_csv("../results/observables/plaq.csv")
    data_plaq = {col: df_plaq[col].values for col in df_plaq.columns}
    df_qtop = pd.read_csv("../results/observables/qtop.csv")
    data_qtop = {col: df_qtop[col].values for col in df_qtop.columns}

    t = data_dH['t']
    dH = data_dH['dH']
    acc = data_acc['acc']
    plaq = data_plaq['plaq']
    qtop = data_qtop['qtop']

    fig, ((ax_dH, ax_acc), (ax_plaq, ax_qtop)) = plt.subplots(2, 2, figsize=(12, 6))
    ax_dH.plot(t, dH, '.-')
    ax_dH.set_xlabel(r"Simulation time $\tau$")
    ax_dH.set_ylabel(r"Energy change $dH$")
    ax_acc.plot(t, acc, '.-')
    ax_acc.set_xlabel(r"Simulation time $\tau$")
    ax_acc.set_ylabel(r"Acceptance")
    ax_plaq.plot(t, plaq, '.-')
    ax_plaq.set_xlabel(r"Simulation time $\tau$")
    ax_plaq.set_ylabel(r"Plaquette")
    ax_qtop.plot(t, qtop, '.-')
    ax_qtop.set_xlabel(r"Simulation time $\tau$")
    ax_qtop.set_ylabel(r"Topological charge")

    fig.tight_layout()

    plot_file = "../results/plots/general.png"
    fig.savefig(plot_file)
    print("Plot saved to "+plot_file)